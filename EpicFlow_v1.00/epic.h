/* This file contains the main functions for the sparse-to-dense interpolation */
#include "array_types.h"
#include "image.h"

/* parameter of epic */
typedef struct epic_params_s {
    char method[20];     // method for interpolation: la (locally-weighted affine) or nw (nadaraya-watson)
    float saliency_th;   // matches coming from pixels with a saliency below this threshold are removed before interpolation
    int pref_nn;         // number of neighbors for consistent checking
    float pref_th;       // threshold for the first prefiltering step
    int nn;              // number of neighbors to consider for the interpolation
    float coef_kernel;   // coefficient in the sigmoid of the interpolation kernel
    float euc;           // constant added to the edge cost
    int verbose;         // verbose mode
} epic_params_t;

/* set params to default value */
void epic_params_default(epic_params_t* params);


/* main function for edge-preserving interpolation of correspondences
    flowx                  x-component of the flow (output)
    flowy                  y-component of the flow (output)
    im                     first image (in lab colorspace)
    input_matches          input matches with each line a match and the first four columns containing x1 y1 x2 y2 
    edges                  edges cost (can be modified)
    params                 parameters
    n_thread               number of threads
*/
void epic(image_t *flowx, image_t *flowy, const color_image_t *im, const float_image *input_matches, float_image* edges, const epic_params_t* params, const int n_thread);

