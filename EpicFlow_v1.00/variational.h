#ifdef __cplusplus
extern "C" {
#endif

#ifndef __VARIATIONAL_H_
#define __VARIATIONAL_H_


#include <stdio.h>
#include <stdlib.h>

#include "image.h"
#include "array_types.h"

typedef struct variational_params_s {
  float alpha;             // smoothness weight
  float gamma;             // gradient constancy assumption weight
  float delta;             // color constancy assumption weight
  float sigma;             // presmoothing of the images
  int niter_outer;         // number of outer fixed point iterations
  int niter_inner;         // number of inner fixed point iterations
  int niter_solver;        // number of solver iterations 
  float sor_omega;         // omega parameter of sor method
} variational_params_t;

/* set flow parameters to default */
void variational_params_default(variational_params_t *params);

/* Compute a refinement of the optical flow (wx and wy are modified) between im1 and im2 */
void variational(image_t *wx, image_t *wy, const color_image_t *im1, const color_image_t *im2, variational_params_t *params);

#endif

#ifdef __cplusplus
}
#endif
