#ifndef ___ARRAY_TYPES_H___
#define ___ARRAY_TYPES_H___

typedef unsigned char UBYTE;
typedef unsigned int UINT;

#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(sizeof(type),(n))
#define NEW(type) NEWA(type,1)
#define REALLOC(ptr,type,n) ptr = (type*)realloc(ptr, sizeof(type)*(n))

#define MIN(a,b)  (((a)<(b)) ? (a) : (b))
#define MAX(a,b)  (((a)>(b)) ? (a) : (b))
#define SWAP(a,b,type)  {type _t = a; a = b; b = _t;}
#define between(min,val,max)  (min<=val && val<=max)

static inline float pow2( float f ) {
  return f*f;
}


/*const float NaN = 0.0/0.0;
const int INT_MIN = 0x80000000;
const int INT_MAX = 0x7FFFFFFF;*/

/************************
* Tuples
*/

#define DEFINE_TUPLE(type,nval)  \
    typedef type type ## _Tuple ## nval;

DEFINE_TUPLE(int,2)
DEFINE_TUPLE(int,3)
DEFINE_TUPLE(int,4)
DEFINE_TUPLE(float,2)
DEFINE_TUPLE(float,3)
DEFINE_TUPLE(float,4)
DEFINE_TUPLE(float,8)


/************************
* 1D Array
*/

#define DEFINE_ARRAY(type)  \
    typedef struct {  \
      type* pixels; \
      int tx; \
    } type##_array;

DEFINE_ARRAY(UBYTE)
DEFINE_ARRAY(int)
DEFINE_ARRAY(UINT)
DEFINE_ARRAY(float)

#define ASSERT_ARRAY_ZEROS(arr) {int size=arr->tx; assert((arr->pixels[0]==0 && arr->pixels[size/2]==0 && arr->pixels[size-1]==0) || !"error: matrix " #arr "is supposed to be zeros");}


/************************
* 2D Image
*/

#define DEFINE_IMG(type)    \
    typedef struct { \
      type* pixels;\
      int tx,ty;\
    } type##_image;

DEFINE_IMG(UBYTE)
DEFINE_IMG(int)
DEFINE_IMG(UINT)
DEFINE_IMG(float)

#define ASSERT_SAME_SIZE  ASSERT_SAME_IMG_SIZE
#define ASSERT_IMG_SIZE  ASSERT_SAME_IMG_SIZE
#define ASSERT_SAME_IMG_SIZE(im1,im2)  if(im1 && im2)  assert(im1->tx==im2->tx && im1->ty==im2->ty);

#define ASSERT_IMAGE_ZEROS
#define ASSERT_IMG_ZEROS(img) {int size=img->tx*img->ty; assert((img->pixels[0]==0 && img->pixels[size/2]==0 && img->pixels[size-1]==0) || !"error: matrix " #img "is supposed to be zeros");}


/************************
* RGB Image
*/

#define DEFINE_IMG3(type) \
    typedef struct {  \
      type* pixels;  \
      int tx,ty;  \
    } type##_image3;

DEFINE_IMG3(UBYTE)
DEFINE_IMG3(float)


/************************
* 3D Image = Cube (Z coordinates are contiguous)
*/

#define DEFINE_CUBE(type) \
    typedef struct {  \
      type* pixels;  \
      int tx,ty,tz;  \
    } type##_cube;

DEFINE_CUBE(UBYTE)
DEFINE_CUBE(int)
DEFINE_CUBE(UINT)
DEFINE_CUBE(float)

#define ASSERT_CUBE_ZEROS(img) {int size=img->tx*img->ty*img->tz; assert((img->pixels[0]==0 && img->pixels[size/2]==0 && img->pixels[size-1]==0) || !"error: matrix " #img "is supposed to be zeros");}


/************************
* 3D Image = concatenation of XY layers
*/

#define DEFINE_LAYERS(type) \
    typedef struct {  \
      type* pixels;  \
      int tx,ty,tz;  \
    } type##_layers;

DEFINE_LAYERS(UBYTE)
DEFINE_LAYERS(int)
DEFINE_LAYERS(UINT)
DEFINE_LAYERS(float)


#define ASSERT_SAME_LAYERS_SIZE(im1,im2)  (if(im1 && im2)  assert(im1->tx==im2->tx && im1->ty==im2->ty && im1->tz==im2->tz));

#define ASSERT_LAYERS_ZEROS ASSERT_CUBE_ZEROS


/************************
* Sparse COO matrix
*/

typedef struct  {
  int* indices; //  indices of columns
  int* indptr;  //  indices of indices of rows
  float* data;  //  row i contains values data[indptr[i]:indptr[i+1]] at columns indices[indptr[i]:indptr[i+1]]
  int nr;
  int nc;
} csr_matrix;


/************************
* Sparse COO matrix
*/

typedef struct  {
  int* row;
  int* col;
  float* data;
  int nr,nc,n_elem;
} coo_matrix;

#define empty_array(type,tx)        ((type##_array){NEWA(type,(tx)),tx})
#define empty_image(type,tx,ty)     ((type##_image){NEWA(type,(tx)*(ty)),tx,ty})
#define empty_cube(type,tx,ty,tz)   ((type##_cube  ){NEWA(type,(tx)*(ty)*(tz)),tx,ty,tz})
#define empty_layers(type,tx,ty,tz) ((type##_layers){NEWA(type,(tx)*(ty)*(tz)),tx,ty,tz})

#define zeros_array(type,tx)        ((type##_array){NEWAC(type,(tx)),tx})
#define zeros_image(type,tx,ty)     ((type##_image){NEWAC(type,(tx)*(ty)),tx,ty})
#define zeros_cube(type,tx,ty,tz)   ((type##_cube  ){NEWAC(type,(tx)*(ty)*(tz)),tx,ty,tz})
#define zeros_layers(type,tx,ty,tz) ((type##_layers){NEWAC(type,(tx)*(ty)*(tz)),tx,ty,tz})

#define free_image(img) {free(img->pixels); free(img); img=NULL;}
#define free_array(img) {free(img->pixels); free(img); img=NULL;}


#endif























