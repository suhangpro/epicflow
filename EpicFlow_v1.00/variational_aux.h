#ifdef __cplusplus
extern "C" {
#endif

#ifndef __VARIATIONAL_AUX_H_
#define __VARIATIONAL_AUX_H_

#include <stdlib.h>
#include "image.h"

/* warp a color image according to a flow. src is the input image, wx and wy, the input flow. dst is the warped image and mask contains 0 or 1 if the pixels goes outside/inside image boundaries */
void image_warp(color_image_t *dst, image_t *mask, const color_image_t *src, const image_t *wx, const image_t *wy);

/* compute image first and second order spatio-temporal derivatives of a color image */
void get_derivatives(const color_image_t *im1, const color_image_t *im2, const convolution_t *deriv, color_image_t *dx, color_image_t *dy, color_image_t *dt, color_image_t *dxx, color_image_t *dxy, color_image_t *dyy, color_image_t *dxt, color_image_t *dyt);

/* compute the smoothness term */
void compute_smoothness(image_t *dst_horiz, image_t *dst_vert, const image_t *uu, const image_t *vv, const image_t *dpsis_weight, const convolution_t *deriv_flow, const float half_alpha);

/* sub the laplacian (smoothness term) to the right-hand term */
void sub_laplacian(image_t *dst, const image_t *src, const image_t *weight_horiz, const image_t *weight_vert);

/* compute local smoothness weight as a sigmoid on image gradient*/
image_t* compute_dpsis_weight(color_image_t *im, float coef, const convolution_t *deriv);

/* compute the dataterm and the matching term
   a11 a12 a22 represents the 2x2 diagonal matrix, b1 and b2 the right hand side
   other (color) images are input */
void compute_data_and_match(image_t *a11, image_t *a12, image_t *a22, image_t *b1, image_t *b2, image_t *mask, image_t *du, image_t *dv, color_image_t *Ix, color_image_t *Iy, color_image_t *Iz, color_image_t *Ixx, color_image_t *Ixy, color_image_t *Iyy, color_image_t *Ixz, color_image_t *Iyz, const float half_delta_over3, const float half_gamma_over3);

#endif

#ifdef __cplusplus
}
#endif
