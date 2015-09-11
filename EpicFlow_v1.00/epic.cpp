#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "epic.h"
#include "image.h"
#include "array_types.h"
#include "epic_aux.h"

#include "omp.h"

#include "io.h"

/* create a copy of input matches with 4 columns, with all points inside the image area*/
static float_image rectify_corres(const float_image* matches, const int w1, const int h1, const int w2, const int h2, const int n_thread){
    float_image res = empty_image(float, 4, matches->ty);
    int i;
    #if defined(USE_OPENMP)
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for(i=0 ; i<matches->ty ; i++){
        res.pixels[4*i  ] = MAX(0,MIN(matches->pixels[i*matches->tx  ],w1-1));
        res.pixels[4*i+1] = MAX(0,MIN(matches->pixels[i*matches->tx+1],h1-1));
        res.pixels[4*i+2] = MAX(0,MIN(matches->pixels[i*matches->tx+2],w2-1));
        res.pixels[4*i+3] = MAX(0,MIN(matches->pixels[i*matches->tx+3],h2-1));
    }
    return res;
}

/* given a set of matches, return the set of points in the first image where a match exists */
static int_image matches_to_seeds(const float_image *matches, const int n_thread){
    int_image res = empty_image(int, 2, matches->ty);
    int i;
    #if defined(USE_OPENMP)
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for(i=0;i<matches->ty;i++){
        res.pixels[2*i] = (int) matches->pixels[4*i];
        res.pixels[2*i+1] = (int) matches->pixels[4*i+1];
    }
    return res;
}

/* given a set of matches, return the set of vecotrs of the matches*/
static float_image matches_to_vects(const float_image *matches, const int n_thread){
    float_image res = empty_image(float, 2, matches->ty);
    int i;
    #if defined(USE_OPENMP)
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for(i=0;i<matches->ty;i++){
        res.pixels[2*i] = matches->pixels[4*i+2]-matches->pixels[4*i];
        res.pixels[2*i+1] = matches->pixels[4*i+3]-matches->pixels[4*i+1];
    }
    return res;
}

/* remove matches coming from a pixel with a low saliency */
static void apply_saliency_threshold(float_image *matches, const color_image_t *im, const float saliency_threshold){
    image_t *s = saliency(im, 0.8f, 1.0f);
    float_image tmp = empty_image(float, 4, matches->ty);
    int i, ii=0;
    for(i=0 ; i<matches->ty ; i++){
        if( s->data[ (int) (matches->pixels[i*4+1]*s->stride+matches->pixels[i*4]) ] >= saliency_threshold ){
            memcpy( &tmp.pixels[ii*4], &matches->pixels[i*4], sizeof(float)*4 );
            ii += 1;
        }
    }   
    image_delete(s);
    REALLOC(tmp.pixels, float, ii*4);
    matches->ty = ii;
    free(matches->pixels);
    matches->pixels = tmp.pixels;
}

/* remove matches where the nadaraya-watson estimation is too different from the input match */
static void prefiltering( float_image *matches, const float_image *edges, const int nn, const float threshold, const float coef_kernel, const int n_thread){
    const float th2 = threshold*threshold;
    const int nns = MIN(nn+1, matches->ty); // nn closest plus itself
    if( nns != nn+1 ) fprintf(stderr, "Warning: not enough matches for prefiltering\n");
    int_image seeds = matches_to_seeds(matches, n_thread);
    float_image vects = matches_to_vects(matches, n_thread);
    
    // compute closest matches
    int_image nnf = empty_image( int, nns, matches->ty);
    float_image dis = empty_image( float, nns, matches->ty);
    int_image labels = empty_image( int, edges->tx, edges->ty);
    dist_trf_nnfield_subset( &nnf, &dis, &labels, &seeds, edges,  NULL, &seeds, n_thread);
    
    // apply kernel to the distance
    #if defined(USE_OPENMP)
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for(int i=0 ; i<dis.tx*dis.ty ; i++){
        dis.pixels[i] = expf(-coef_kernel*dis.pixels[i])+1e-08;
    }
    
    // compute nadaraya-watson estimation
    float_image seedsvects = empty_image(float, 2, matches->ty);
    fit_nadarayawatson(&seedsvects, &nnf, &dis, &vects, n_thread);
       
    // remove matches if necessary
    float_image tmp = empty_image(float, 4, matches->ty);
    int ii=0;
    for( int i=0 ; i<matches->ty ; i++ ){
        if( pow2(seedsvects.pixels[2*i]-vects.pixels[2*i]) + pow2(seedsvects.pixels[2*i+1]-vects.pixels[2*i+1])<th2 ){
            memcpy( &tmp.pixels[ii*4], &matches->pixels[i*4], sizeof(float)*4 );
            ii += 1;            
        }   
    }
    REALLOC(tmp.pixels, float, ii*4);
    matches->ty = ii;
    free(matches->pixels);
    matches->pixels = tmp.pixels;

    // free memory
    free(seeds.pixels);
    free(vects.pixels);
    free(nnf.pixels);
    free(dis.pixels);
    free(labels.pixels);
    free(seedsvects.pixels);
}


/* set params to default value */
void epic_params_default(epic_params_t* params){
    strcpy(params->method, "LA");
    params->saliency_th = 0.045f;
    params->pref_nn = 25;
    params->pref_th = 5.0f;
    params->nn = 100;
    params->coef_kernel = 0.8f;
    params->euc = 0.001f;
    params->verbose = 0;
}

/* main function for edge-preserving interpolation of correspondences
    flowx                  x-component of the flow (output)
    flowy                  y-component of the flow (output)
    input_matches          input matches with each line a match and the first four columns containing x1 y1 x2 y2 
    im                     first image (in lab colorspace)
    edges                  edges cost (can be modified)
    params                 parameters
    n_thread               number of threads
*/
void epic(image_t *flowx, image_t *flowy, const color_image_t *im, const float_image *input_matches, float_image* edges, const epic_params_t* params, const int n_thread){   

    // copy matches and correct them if necessary
    float_image matches = rectify_corres(input_matches, im->width, im->height, im->width, im->height, n_thread);
    if( params->verbose ) printf("%d input matches\n", matches.ty);
    
        
    // eventually add a constant to edges cost
    if(params->euc){
        int i;
        #if defined(USE_OPENMP)
        #pragma omp parallel for num_threads(n_thread)
        #endif
        for(i=0 ; i<edges->tx*edges->ty ; i++){
            edges->pixels[i] += params->euc;
        }
    }

    // saliency filter
    if(params->saliency_th){
        apply_saliency_threshold( &matches, im, params->saliency_th);
        if( params->verbose ) printf("Saliency filtering, remaining %d matches\n", matches.ty);
    }
    // consistency filter
    if(params->pref_nn){
        prefiltering( &matches, edges, params->pref_nn, params->pref_th, params->coef_kernel, n_thread);
        if( params->verbose ) printf("Consistenct filter, remaining %d matches\n", matches.ty);
    } 
      
    // prepare variables
    const int nns = MIN(params->nn, matches.ty);
    if( nns < params->nn ) fprintf(stderr, "Warning: not enough matches for interpolating\n");
    if( params->verbose ) printf("Computing %d nearest neighbors for each match\n", nns);
    int_image seeds = matches_to_seeds(&matches, n_thread);
    float_image vects = matches_to_vects(&matches, n_thread);
    
    // compute nearest matches for each seed
    int_image nnf = empty_image( int, nns, matches.ty);
    float_image dis = empty_image( float, nns, matches.ty);
    int_image labels = empty_image( int, edges->tx, edges->ty);
    dist_trf_nnfield_subset( &nnf, &dis, &labels, &seeds, edges, NULL, &seeds, n_thread);  
           
    // apply kernel to the distance
    #if defined(USE_OPENMP)
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for(int i=0 ; i<dis.tx*dis.ty ; i++){
        dis.pixels[i] = expf(-params->coef_kernel*dis.pixels[i])+1e-08;
    }
    
    // interpolation
    if( params->verbose ) printf("Interpolation of matches using %s\n", params->method);
    float_image newvects = empty_image( float, 2, im->width*im->height);
    if( !strcmp( params->method, "LA") ){
        float_image seedsaffine = empty_image(float, 6, vects.ty);
        fit_localaffine(&seedsaffine, &nnf, &dis, &seeds, &vects);
        apply_localaffine(&newvects, &seedsaffine, &labels, n_thread);
        free(seedsaffine.pixels);
    } else if ( !strcmp( params->method, "NW") ){
        float_image seedsvects = empty_image(float, 2, vects.ty);
        fit_nadarayawatson(&seedsvects, &nnf, &dis, &vects, n_thread);
        apply_nadarayawatson(&newvects, &seedsvects, &labels, n_thread);        
        free(seedsvects.pixels);
    } else {
        fprintf(stderr, "method %s not recognized\n", params->method);
        exit(EXIT_FAILURE);
    }

    // copy result to the output
    #if defined(USE_OPENMP)
    #pragma omp parallel for num_threads(n_thread)
    #endif
    for(int i=0 ; i<im->height ; i++){
        for( int j=0 ; j<im->width ; j++){
            flowx->data[i*im->stride+j] = newvects.pixels[2*(i*im->width+j)];
            flowy->data[i*im->stride+j] = newvects.pixels[2*(i*im->width+j)+1];
        }
    }  
    
    // free memory
    free(seeds.pixels);
    free(vects.pixels);
    free(nnf.pixels);
    free(dis.pixels);
    free(labels.pixels);
    free(newvects.pixels);
    free(matches.pixels);
}

