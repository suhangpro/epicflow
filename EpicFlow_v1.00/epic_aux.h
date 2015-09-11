/* This file contains auxiliary function for the interpolation: looking fpr the nearest neighbors with a geodesic distance, fitting/applying the interpolation */
#include "array_types.h"

/* structure for distance transform parameters */
typedef struct {
  int max_iter;
  float min_change;
} dt_params_t;


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
                                const int_image* seeds, const float_image* cost, dt_params_t* dt_params,
                                const int_image* pixels, const int n_thread );
                       


/* fit nadaraya watson for a set of seeds 
    res:  (output) a nseeds*2 array containing the estimated displacement for each seed 
    nnf:   array of size nseeds*nn with index of the nn-closest seed
    dis:   as nnf but containing the distances to the corresponding seed
    vects: 2D vector of the matches
*/
void fit_nadarayawatson(float_image *res, const int_image *nnf, const float_image *dis, const float_image *vects, const int n_thread);

/* apply nadaraya watson interpolation
    newvects:    output containing the flow vector for each pixel 
    seedsvects:  input containing the estimated flow vector for each seed
    labels:      closest seed for each pixel
*/
void apply_nadarayawatson(float_image *newvects, const float_image *seedsvects, const int_image *labels, const int n_thread);


/* fit locally-weighted affine interpolation
    res: (output) a nseeds*6 array containing the estimated affine model for each seed
    nnf:   array of size nseeds*nn with index of the nn-closest seed
    dis:   as nnf but containing the distances to the corresponding seed
    seeds: original point of matches
    vects: 2D vector of the matches
*/
void fit_localaffine(float_image *res, const int_image *nnf, const float_image *dis, const int_image *seeds, const float_image *vects);


/* apply locally-weighted affine interpolation 
    newvects:    output containing the flow vector for each pixel
    seedaffine:  esimated affine transformation for each seed
    labels:      closest seed for each pixel
*/
void apply_localaffine(float_image *newvects, const float_image *seedsaffine, const int_image *labels, const int n_thread);

