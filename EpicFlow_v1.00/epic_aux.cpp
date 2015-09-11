#include "epic_aux.h"
#include "omp.h"
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <string.h>
#include <assert.h>

// include sgels from blas to solve overdetermined linear system
extern "C" {
    #define integer int
    #define real float
    extern int sgels_(char *trans, integer *m, integer *n, integer *nrhs, 
                        real *a, integer *lda, real *b, integer *ldb, 
                        real *work, integer *lwork, integer *info);
}

using namespace std;

const float INF = 1.0f/0.0f;


/*************************** NN-search on a graph ********************************/

/* auxiliary structures for nearest-neighbor search on a graph */
typedef vector<int> list_node_t;

struct node_dist_t {
  int node;
  float dis;
  node_dist_t(){}
  node_dist_t(int i, float f):node(i),dis(f){}
};
typedef vector<node_dist_t> list_node_dist_t;
typedef node_dist_t current_t;
template<typename T>
struct smallest_on_top {
  bool operator() (const T& a, const T& b) const {
    return a.dis > b.dis; }
};

/* Find nearest neighbors in a weighted directed graph (The matrix must be symmetric !) */
static int find_nn_graph_arr( csr_matrix* graph, int seed, int nmax, int* best, float* dist ) {

    assert(nmax>0);
    assert( graph->nr==graph->nc && (graph->indptr[graph->nr]%2==0) );
    const int* indptr = graph->indptr;

    // init done to INF
    float* done = NEWA(float,graph->nr);
    memset(done,0x7F,graph->nr*sizeof(float));
  
    // explore nodes in order of increasing distances
    priority_queue<current_t,vector<current_t>,smallest_on_top<current_t> > stack;
    stack.emplace(seed,0);
    done[seed] = 0;  // mark as done
  
    int n=0;
    while(stack.size()) {
        current_t cur = stack.top();
        stack.pop();
        if(cur.dis > done[cur.node]) continue;
        
        // insert result
        best[n] = cur.node;
        dist[n] = cur.dis;
        n++;
        if( n>= nmax ) break;
        
        // find nearest neighbors
        for(int i=indptr[cur.node]; i<indptr[cur.node+1]; i++) {
            int neigh = graph->indices[i];
            float newd = cur.dis + graph->data[i];
            if( newd>=done[neigh] ) continue;
            stack.emplace(neigh, newd); // add only if it makes sense
            done[neigh] = newd;
        }
    }
  
    free(done);

    // in case we do not get enough results
    memset(best+n,0xFF,(nmax-n)*sizeof(int));
    memset(dist+n,0x7F,(nmax-n)*sizeof(float));
    return n;
}


/******************* DISTANCE TRANSFORM **************/

static float arg_sweep( const float_image* cost, float_image* res, int_image* labels, const char x, const char y ) {
  int i, j;
  const int tx = res->tx, ty = res->ty;
  float* A = res->pixels;
  int* L = labels->pixels;
  const float* Cost = cost->pixels;
  
  const int bx = x>0 ? 0 : tx-1;
  const int by = y>0 ? 0 : ty-1;
  const int ex = x>0 ? tx : -1;
  const int ey = y>0 ? ty : -1;
  
  float t0, t1, t2, C, max_diff = 0.0;
  int l0, l1, l2;
  for(j=by; j!=ey; j+=y)
    for(i=bx; i!=ex; i+=x) {
      if(j==by) {
        t1 = INF;
        l1 = -1;
      } else {
        t1 = A[i + (j-y)*tx];
        l1 = L[i + (j-y)*tx];
      }
      if(i==bx) {
        t2 = INF;
        l2 = -1;
      } else {
        t2 = A[i-x + j*tx];
        l2 = L[i-x + j*tx];
      }
      float dt12 = fabs(t1-t2);
      C = Cost[i + j*tx];
      
      if( dt12 > C ) {  // handle degenerate case
        if( t1 < t2 ) {
          t0 = t1 + C;
          l0 = l1;
        } else {
          t0 = t2 + C;
          l0 = l2;
        }
      } else {
        t0 = 0.5*(t1 + t2 + sqrtf(2*C*C - dt12*dt12));
        l0 = (t1<t2) ? l1 : l2;
      }
      
      if( t0 < A[i + j*tx] )  {
        max_diff = MAX(max_diff, A[i + j*tx] - t0);
        A[i + j*tx] = t0;
        L[i + j*tx] = l0;
      }
  }
  
  return max_diff;
}

void set_default_dt_params( dt_params_t* params ) {
  params->max_iter = 40;
  params->min_change = 1;
}

#define DEFAULT_dt_params(dt_params)  \
  dt_params_t tmp_dt_params;\
  if(!dt_params){set_default_dt_params(&tmp_dt_params);dt_params=&tmp_dt_params;}

/* Compute distance map from a given seeds (in res) and a cost map.
  if labels!=NULL:  labels are propagated along with distance map 
                    (for each pixel, we remember the closest seed)*/
static float weighted_distance_transform( const float_image* cost, const dt_params_t* dt_params, 
                                   float_image* res, int_image* labels ) {
  assert( cost && res );
  ASSERT_SAME_SIZE(cost,res);
  if(labels)  ASSERT_SAME_SIZE(res,labels);
  DEFAULT_dt_params(dt_params)
  assert( dt_params->min_change >= 0 );
  assert(labels);

  const char x[4] = {-1,1,1,-1};
  const char y[4] = {1,1,-1,-1};
  int i = 0, end_iter = 4;
  float change = -1;
  while(++i <= end_iter) {
    change = arg_sweep(cost, res, labels, x[i%4], y[i%4]);
    if( change > dt_params->min_change )
      end_iter = MIN(dt_params->max_iter, i+3); // finish the turn
  }
  return change;
}



/****************** BUILD NEIGHBORHOOD GRAPH *************/

/* structure for the border between two regions in the assignment map */
struct border_t {
  float accu; /* cost as the minimum over the border */
  int nb; /* number of pixels */
  
  border_t():accu(0),nb(0){}
  
  void add( float v) {
    if(!nb || accu>v) accu=v;
    nb++;
  }
  
  float get_val() const {
     return accu;
  }
};

// a border between two int is represented as a long
static inline long key( long i, long j ) {
  if(j>i) SWAP(i,j,long);   // always i<j 
  return i + (j<<32);
}
static inline int first(long i) {return int(i);}
static inline int second(long i) {return int(i >> 32);}

template<typename Ti>
struct Tint_float_t {
  Ti i;
  float f;
  Tint_float_t(){}
  Tint_float_t(Ti i, float f):i(i),f(f){}
};
typedef Tint_float_t<int> int_float_t;
typedef Tint_float_t<long> long_float_t;
static int cmp_long_float ( const void* a, const void* b ) {
  long diff = ((long_float_t*)a)->i - ((long_float_t*)b)->i;
  return (diff>0) - (diff<0);
}


/* Find neighboring labels and store their relation in a sparse matrix */
static void ngh_labels_to_spmat( int ns, int_image* labels, float_image* dmap,  csr_matrix* csr ) {
  ASSERT_SAME_SIZE(labels, dmap);
  const int tx = labels->tx;
  const int ty = labels->ty;
  const int* lab = labels->pixels;
  const float* dis = dmap->pixels;
  
  typedef unordered_map<long,border_t> ngh_t;
  ngh_t ngh;
  for(int j=1; j<ty; j++)
    for(int i=1; i<tx; i++) {
      int l0 = lab[i+j*tx];
      int l1 = lab[i-1+j*tx];
      int l2 = lab[i+(j-1)*tx];
      if( l0!=l1 ) {
        long k = key(l0,l1);
        float d = dis[i+j*tx] + dis[i-1+j*tx];
        ngh[k].add(d);
      }
      if( l0!=l2 ) {
        long k = key(l0,l2);
        float d = dis[i+j*tx] + dis[i+(j-1)*tx];
        ngh[k].add(d);
      }
    }
  
  // convert result into a sparse graph
  vector<long_float_t> sorted;
  for(ngh_t::iterator it=ngh.begin(); it!=ngh.end(); ++it) {
    const float cost = it->second.get_val();
    sorted.emplace_back(it->first,cost);
    long symkey = (it->first>>32) + (it->first<<32);
    sorted.emplace_back(symkey,cost);  // add symmetric
  }
  // sort by (row,col) = (second,first)
  qsort(sorted.data(), sorted.size(), sizeof(long_float_t), cmp_long_float);
  
  csr->nr = ns;
  csr->nc = ns;
  assert(!csr->indptr);
  csr->indptr = NEWA(int, ns+1);
  csr->indices = NEWA(int, sorted.size());
  csr->data = NEWA(float, sorted.size());
  int n=0,r=0;
  csr->indptr[0] = 0;
  while( n<(signed)sorted.size() ) {
    int row = second(sorted[n].i);
    //assert(r<=row);
    while(r<row)  csr->indptr[++r] = n;
    csr->indices[n] = first(sorted[n].i);
    csr->data[n] = sorted[n].f;
    n++;
  }
  assert(r<ns);
  // finish row marker
  while(r<ns)
    csr->indptr[++r] = n;
}

/* Compute the neighboring matrix between seeds as well as the closest seed for each pixel */
static void distance_transform_and_graph( const int_image* seeds,  const float_image* cost, dt_params_t* dt_params,
                                int_image* labels, float_image* dmap, csr_matrix* ngh, int n_thread ) {
  const int tx = cost->tx;
  const int ty = cost->ty;
  assert(seeds->tx==2);
  const int ns = seeds->ty;
  
  
  // compute distance transform for all seeds altogether
  if(labels->pixels)
    assert(labels->tx==tx && labels->ty==ty);
  else
    *labels = {NEWA(int,tx*ty),tx,ty};
  if(dmap->pixels)
    assert(dmap->tx==tx && dmap->ty==ty);
  else
    *dmap = {NEWA(float,tx*ty),tx,ty};
  memset(dmap->pixels,0x7F,tx*ty*sizeof(float));
  for(int i=0; i<ns; i++) {
    int p = seeds->pixels[2*i] + seeds->pixels[2*i+1]*tx;
    labels->pixels[p] = i;
    dmap->pixels[p] = cost->pixels[p];
  }
  weighted_distance_transform( cost, dt_params, dmap, labels );
  
  
  // compute distances between neighboring seeds
  ngh_labels_to_spmat( ns, labels, dmap, ngh );
  
}



/* Compute the closest seeds using a geodesic distance for a subset of points, and the assignment of each pixel to the closest seeds
    best:    output containing the closest seeds for each query point
    dist:    output containing the distances to the closest seeds
    labels:  output containing the assignment of each pixel to the closest seed
    seeds:   2D positions of the seeds
    cost:    cost of going throw a pixel (ie that defines the geodesic distance)
    dt_params: distance transform parameters (NULL for default parameters)
    pixels:  2D positions of the query points
*/
void dist_trf_nnfield_subset( int_image* best, float_image* dist, int_image *labels,
                               const int_image* seeds,
                               const float_image* cost, dt_params_t* dt_params,
                               const int_image* pixels, const int n_thread ) {
  ASSERT_SAME_SIZE(best,dist);
  const int tx = cost->tx;
  const int npix = pixels->ty;
  assert(best->ty==npix && dist->ty==npix);
  const int ns = seeds->ty;
  const int nn = best->tx;
  assert(dist->tx==nn);
  
  float_image dmap = {0};
  int_image  nnf = {NEWA(int,ns*nn),nn,ns};
  float_image dis = {NEWA(float,ns*nn),nn,ns};

  // compute distance transform and build graph
  csr_matrix ngh = {0};
  distance_transform_and_graph( seeds,  cost, dt_params, labels, &dmap, &ngh, n_thread );
  
  // compute nearest neighbors using the graph
  #ifdef USE_OPENMP
  #pragma omp parallel for num_threads(n_thread)
  #endif
  for(int n=0; n<ns; n++)
    find_nn_graph_arr( &ngh, n, nn, nnf.pixels+n*nn, dis.pixels+n*nn );
  free(ngh.indptr);
  free(ngh.indices);
  free(ngh.data);  
  
  // for each pixel, just look at closest seed and copy
  #ifdef USE_OPENMP
  #pragma omp parallel for num_threads(n_thread)
  #endif
  for(int _i=0; _i<npix; _i++) {
    const int i = pixels->pixels[2*_i] + pixels->pixels[2*_i+1]*tx;
    // this pixel is at distance <d> from seed <s>
    int s = labels->pixels[i];
    float d = dmap.pixels[i];
    memcpy(best->pixels+_i*nn,nnf.pixels+s*nn,nn*sizeof(nn));
    for(int j=0; j<nn; j++)
      dist->pixels[_i*nn+j] = d+dis.pixels[s*nn+j];  // distance plus constant
  }
  
  free(nnf.pixels);
  free(dis.pixels);
  free(dmap.pixels);
}


/******************************** FITTING/APPLYING INTERPOLATION **********************/

/* fit nadaraya watson for a set of seeds 
    res:  (output) a nseeds*2 array containing the estimated displacement for each seed 
    nnf:   array of size nseeds*nn with index of the nn-closest seed
    dis:   as nnf but containing the distances to the corresponding seed
    vects: 2D vector of the matches
*/
void fit_nadarayawatson(float_image *res, const int_image *nnf, const float_image *dis, const float_image *vects, const int n_thread){
    const int nn = nnf->tx;
    #ifdef USE_OPENMP
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for(int i = 0 ; i<vects->ty ; i++){
        float u = 0.0f, v = 0.0f, s = 0.0f;
        for(int j = i*nn ; j<((i+1)*nn) ; j++){
            const float d = dis->pixels[j];
            const int jj = nnf->pixels[j];
            u += d*vects->pixels[2*jj];
            v += d*vects->pixels[2*jj+1];
            s += d;
	    }
        res->pixels[2*i  ] = u/s;
        res->pixels[2*i+1] = v/s;
    }
}

/* apply nadaraya watson interpolation
    newvects:    output containing the flow vector for each pixel 
    seedsvects:  input containing the estimated flow vector for each seed
    labels:      closest seed for each pixel
*/
void apply_nadarayawatson(float_image *newvects, const float_image *seedsvects, const int_image *labels, const int n_thread){
    #ifdef USE_OPENMP
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for( int n=0 ; n<labels->tx*labels->ty ; n++ ){
        const int s = labels->pixels[n];
        newvects->pixels[2*n  ] = seedsvects->pixels[2*s];
        newvects->pixels[2*n+1] = seedsvects->pixels[2*s+1];
    }
}

/* fit locally-weighted affine interpolation
    res: (output) a nseeds*6 array containing the estimated affine model for each seed
    nnf:   array of size nseeds*nn with index of the nn-closest seed
    dis:   as nnf but containing the distances to the corresponding seed
    seeds: original point of matches
    vects: 2D vector of the matches
*/
// in addition to the n points, we add 4 matches (close to the current math and with a low weight) in order to ensure non-colinearity
#define AFFINE_ADDPOINT(m,v,x,y,wx,wy,c) {m[0] = m[8] = (x)*(c); m[1] = m[9] = (y)*(c); m[4] = m[11] = (c); v[0] = ((x)+(wx))*(c); v[1] = ((y)+(wy))*(c); m+=12; v+=2;}
void fit_localaffine(float_image *res, const int_image *nnf, const float_image *dis, const int_image *seeds, const float_image *vects){
    const int nn = nnf->tx;
    {
        float *mat = NEWA(float, 12*(nn+4));
        float *vec = NEWA(float, 2*(nn+4));
        for(int i = 0 ; i<nnf->ty ; i++){
            memset(mat, 0, sizeof(float)*12*(nn+4));  
            float coefi=0.0f, *m = mat, *v = vec, *dist = dis->pixels + i*nn;
            int *ind = nnf->pixels + i*nn;
            for( int j = 0 ; j<nn ; j++){
                const int s = ind[j];
                float coef = dist[j];
                if(s==i) {coefi = 0.01f*coef; coef *= 0.96f;}
                AFFINE_ADDPOINT(m,v,seeds->pixels[2*s],seeds->pixels[2*s+1],vects->pixels[2*s],vects->pixels[2*s+1],coef);
            }
            float xi = (float) seeds->pixels[2*i], yi = (float) seeds->pixels[2*i+1], ui = vects->pixels[2*i], vi = vects->pixels[2*i+1];
            AFFINE_ADDPOINT(m,v,xi+0.1f,yi,ui,vi,coefi);
            AFFINE_ADDPOINT(m,v,xi,yi+0.1f,ui,vi,coefi);
            AFFINE_ADDPOINT(m,v,xi-0.1f,yi,ui,vi,coefi);
            AFFINE_ADDPOINT(m,v,xi,yi-0.1f,ui,vi,coefi);
            {
               integer info;
               integer mm=6, nrhs=1, lda=6, lwork=-1;
               integer n=2*(nn+4),ldb=n;
               float work_sz;
               sgels_((char*) "Transposed", &mm, &n, &nrhs, mat, &lda, 
                   vec, &ldb, &work_sz, &lwork, &info); 
               float *work;    
               lwork=(int)work_sz;
               work=NEWA(float,lwork);
               sgels_((char*) "Transposed", &mm, &n, &nrhs, mat, &lda, 
                   vec, &ldb, work, &lwork, &info); 
               free(work);
               assert(info>=0); /* there is always a result for coherent input */
            } 
            float *aff = &res->pixels[6*i];
            aff[0]=vec[0]; aff[1]=vec[1]; aff[2]=vec[4];
            aff[3]=vec[2]; aff[4]=vec[3]; aff[5]=vec[5];
        }
        free(mat);
        free(vec); 
    }
}

/* apply locally-weighted affine interpolation 
    newvects:    output containing the flow vector for each pixel
    seedaffine:  esimated affine transformation for each seed
    labels:      closest seed for each pixel
*/
void apply_localaffine(float_image* newvects, const float_image* seedsaffine, const int_image *labels, const int n_thread){
    #ifdef USE_OPENMP
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for( int j = 0 ; j<labels->ty ; j++){
        for(int i=0 ; i<labels->tx ; i++){
            const int n = j*labels->tx+i;
            const int s = labels->pixels[n];
            float *m = &seedsaffine->pixels[6*s];
            newvects->pixels[2*n  ] = m[0]*i + m[1]*j + m[2] -i;
            newvects->pixels[2*n+1] = m[3]*i + m[4]*j + m[5] -j;
	    }
    }
}
