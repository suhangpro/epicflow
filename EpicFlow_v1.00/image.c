#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "image.h"

#include <xmmintrin.h>
typedef __v4sf v4sf;

/********** Create/Delete **********/

/* allocate a new image of size width x height */
image_t *image_new(const int width, const int height){
    image_t *image = (image_t*) malloc(sizeof(image_t));
    if(image == NULL){
        fprintf(stderr, "Error: image_new() - not enough memory !\n");
        exit(1);
    }
    image->width = width;
    image->height = height;  
    image->stride = ( (width+3) / 4 ) * 4;
    image->data = (float*) memalign(16, image->stride*height*sizeof(float));
    if(image->data == NULL){
        fprintf(stderr, "Error: image_new() - not enough memory !\n");
        exit(1);
    }
    return image;
}

/* allocate a new image and copy the content from src */
image_t *image_cpy(const image_t *src){
    image_t *dst = image_new(src->width, src->height);
    memcpy(dst->data, src->data, src->stride*src->height*sizeof(float));
    return dst;
}

/* set all pixels values to zeros */
void image_erase(image_t *image){
    memset(image->data, 0, image->stride*image->height*sizeof(float));
}


/* multiply an image by a scalar */
void image_mul_scalar(image_t *image, const float scalar){
    int i;
    v4sf* imp = (v4sf*) image->data;
    const v4sf scalarp = {scalar,scalar,scalar,scalar};
    for( i=0 ; i<image->stride/4*image->height ; i++){
        (*imp) *= scalarp;
        imp+=1;
    }
}

/* free memory of an image */
void image_delete(image_t *image){
    if(image == NULL){
        //fprintf(stderr, "Warning: Delete image --> Ignore action (image not allocated)\n");
    }else{
    free(image->data);
    free(image);
    }
}


/* allocate a new color image of size width x height */
color_image_t *color_image_new(const int width, const int height){
    color_image_t *image = (color_image_t*) malloc(sizeof(color_image_t));
    if(image == NULL){
        fprintf(stderr, "Error: color_image_new() - not enough memory !\n");
        exit(1);
    }
    image->width = width;
    image->height = height;  
    image->stride = ( (width+3) / 4 ) * 4;
    image->c1 = (float*) memalign(16, 3*image->stride*height*sizeof(float));
    if(image->c1 == NULL){
        fprintf(stderr, "Error: color_image_new() - not enough memory !\n");
        exit(1);
    }
    image->c2 =  image->c1+image->stride*height;
    image->c3 =  image->c2+image->stride*height;
    return image;
}

/* allocate a new color image and copy the content from src */
color_image_t *color_image_cpy(const color_image_t *src){
    color_image_t *dst = color_image_new(src->width, src->height);
    memcpy(dst->c1, src->c1, 3*src->stride*src->height*sizeof(float));
    return dst;
}

/* set all pixels values to zeros */
void color_image_erase(color_image_t *image){
    memset(image->c1, 0, 3*image->stride*image->height*sizeof(float));
}

/* free memory of a color image */
void color_image_delete(color_image_t *image){
    if(image){
        free(image->c1); // c2 and c3 was allocated at the same moment
        free(image);
    }
}


/************ Convolution ******/

/* return half coefficient of a gaussian filter
Details:
- return a float* containing the coefficient from middle to border of the filter, so starting by 0,
- it so contains half of the coefficient.
- sigma is the standard deviation.
- filter_order is an output where the size of the output array is stored */
float *gaussian_filter(const float sigma, int *filter_order){
    if(sigma == 0.0f){
        fprintf(stderr, "gaussian_filter() error: sigma is zeros\n");
        exit(1);
    }
    if(!filter_order){
        fprintf(stderr, "gaussian_filter() error: filter_order is null\n");
        exit(1);
    }
    // computer the filter order as 1 + 2* floor(3*sigma)
    *filter_order = floor(3*sigma)+1; 
    if ( *filter_order == 0 )
        *filter_order = 1; 
    // compute coefficients
    float *data = (float*) malloc(sizeof(float) * (2*(*filter_order)+1));
    if(data == NULL ){
        fprintf(stderr, "gaussian_filter() error: not enough memory\n");
        exit(1);
    }
    const float alpha = 1.0f/(2.0f*sigma*sigma);
    float sum = 0.0f;
    int i;
    for(i=-(*filter_order) ; i<=*filter_order ; i++){
        data[i+(*filter_order)] = exp(-i*i*alpha);
        sum += data[i+(*filter_order)];
    }
    for(i=-(*filter_order) ; i<=*filter_order ; i++){
        data[i+(*filter_order)] /= sum;
    }
    // fill the output
    float *data2 = (float*) malloc(sizeof(float)*(*filter_order+1));
    if(data2 == NULL ){
        fprintf(stderr, "gaussian_filter() error: not enough memory\n");
        exit(1);
    }
    memcpy(data2, &data[*filter_order], sizeof(float)*(*filter_order)+sizeof(float));
    free(data);
    return data2;
}

/* given half of the coef, compute the full coefficients and the accumulated coefficients */
static void convolve_extract_coeffs(const int order, const float *half_coeffs, float *coeffs, float *coeffs_accu, const int even){
    int i;
    float accu = 0.0;
    if(even){
        for(i = 0 ; i <= order; i++){
	        coeffs[order - i] = coeffs[order + i] = half_coeffs[i];
        }
        for(i = 0 ; i <= order; i++){
	        accu += coeffs[i];
	        coeffs_accu[2 * order - i] = coeffs_accu[i] = accu;
        }
    }else{
        for(i = 0; i <= order; i++){
	        coeffs[order - i] = +half_coeffs[i];
	        coeffs[order + i] = -half_coeffs[i];
        }
        for(i = 0 ; i <= order; i++){
            accu += coeffs[i];
	        coeffs_accu[i] = accu;
	        coeffs_accu[2 * order - i]= -accu;
        }
    }
}

/* create a convolution structure with a given order, half_coeffs, symmetric or anti-symmetric according to even parameter */
convolution_t *convolution_new(const int order, const float *half_coeffs, const int even){
    convolution_t *conv = (convolution_t *) malloc(sizeof(convolution_t));
    if(conv == NULL){
        fprintf(stderr, "Error: convolution_new() - not enough memory !\n");
        exit(1);
    }
    conv->order = order;
    conv->coeffs = (float *) malloc((2 * order + 1) * sizeof(float));
    if(conv->coeffs == NULL){
        fprintf(stderr, "Error: convolution_new() - not enough memory !\n");
        free(conv);
        exit(1);
    }
    conv->coeffs_accu = (float *) malloc((2 * order + 1) * sizeof(float));
    if(conv->coeffs_accu == NULL){
        fprintf(stderr, "Error: convolution_new() - not enough memory !\n");
        free(conv->coeffs);
        free(conv);
        exit(1);
    }
    convolve_extract_coeffs(order, half_coeffs, conv->coeffs,conv->coeffs_accu, even);
    return conv;
}

static void convolve_vert_fast_3(image_t *dst, const image_t *src, const convolution_t *conv){
    const int iterline = (src->stride>>2)+1;
    const float *coeff = conv->coeffs;
    //const float *coeff_accu = conv->coeffs_accu;
    v4sf *srcp = (v4sf*) src->data, *dstp = (v4sf*) dst->data;
    v4sf *srcp_p1 = (v4sf*) (src->data+src->stride);
    int i;
    for(i=iterline ; --i ; ){ // first line
        *dstp = (coeff[0]+coeff[1])*(*srcp) + coeff[2]*(*srcp_p1);
        dstp+=1; srcp+=1; srcp_p1+=1;
    }
    v4sf* srcp_m1 = (v4sf*) src->data; 
    for(i=src->height-1 ; --i ; ){ // others line
        int j;
        for(j=iterline ; --j ; ){
            *dstp = coeff[0]*(*srcp_m1) + coeff[1]*(*srcp) + coeff[2]*(*srcp_p1);
            dstp+=1; srcp_m1+=1; srcp+=1; srcp_p1+=1;
        }
    }       
    for(i=iterline ; --i ; ){ // last line
        *dstp = coeff[0]*(*srcp_m1) + (coeff[1]+coeff[2])*(*srcp);
        dstp+=1; srcp_m1+=1; srcp+=1; 
    }  
}

static void convolve_vert_fast_5(image_t *dst, const image_t *src, const convolution_t *conv){
    const int iterline = (src->stride>>2)+1;
    const float *coeff = conv->coeffs;
    //const float *coeff_accu = conv->coeffs_accu;
    v4sf *srcp = (v4sf*) src->data, *dstp = (v4sf*) dst->data;
    v4sf *srcp_p1 = (v4sf*) (src->data+src->stride);
    v4sf *srcp_p2 = (v4sf*) (src->data+2*src->stride);
    int i;
    for(i=iterline ; --i ; ){ // first line
        *dstp = (coeff[0]+coeff[1]+coeff[2])*(*srcp) + coeff[3]*(*srcp_p1) + coeff[4]*(*srcp_p2);
        dstp+=1; srcp+=1; srcp_p1+=1; srcp_p2+=1;
    }
    v4sf* srcp_m1 = (v4sf*) src->data;
    for(i=iterline ; --i ; ){ // second line
        *dstp = (coeff[0]+coeff[1])*(*srcp_m1) + coeff[2]*(*srcp) + coeff[3]*(*srcp_p1) + coeff[4]*(*srcp_p2);
        dstp+=1; srcp_m1+=1; srcp+=1; srcp_p1+=1; srcp_p2+=1;
    }   
    v4sf* srcp_m2 = (v4sf*) src->data;
    for(i=src->height-3 ; --i ; ){ // others line
        int j;
        for(j=iterline ; --j ; ){
            *dstp = coeff[0]*(*srcp_m2) + coeff[1]*(*srcp_m1) + coeff[2]*(*srcp) + coeff[3]*(*srcp_p1) + coeff[4]*(*srcp_p2);
            dstp+=1; srcp_m2+=1;srcp_m1+=1; srcp+=1; srcp_p1+=1; srcp_p2+=1;
        }
    }    
    for(i=iterline ; --i ; ){ // second to last line
        *dstp = coeff[0]*(*srcp_m2) + coeff[1]*(*srcp_m1) + coeff[2]*(*srcp) + (coeff[3]+coeff[4])*(*srcp_p1);
        dstp+=1; srcp_m2+=1;srcp_m1+=1; srcp+=1; srcp_p1+=1;
    }          
    for(i=iterline ; --i ; ){ // last line
        *dstp = coeff[0]*(*srcp_m2) + coeff[1]*(*srcp_m1) + (coeff[2]+coeff[3]+coeff[4])*(*srcp);
        dstp+=1; srcp_m2+=1;srcp_m1+=1; srcp+=1; 
    }  
}

static void convolve_horiz_fast_3(image_t *dst, const image_t *src, const convolution_t *conv){
    const int stride_minus_1 = src->stride-1;
    const int iterline = (src->stride>>2);
    const float *coeff = conv->coeffs;
    v4sf *srcp = (v4sf*) src->data, *dstp = (v4sf*) dst->data;
    // create shifted version of src
    float *src_p1 = (float*) malloc(sizeof(float)*src->stride),
        *src_m1 = (float*) malloc(sizeof(float)*src->stride);
    int j;
    for(j=0;j<src->height;j++){
        int i;
        float *srcptr = (float*) srcp;
        const float right_coef = srcptr[src->width-1];
        for(i=src->width;i<src->stride;i++)
            srcptr[i] = right_coef;
        src_m1[0] = srcptr[0];
        memcpy(src_m1+1, srcptr , sizeof(float)*stride_minus_1);
        src_p1[stride_minus_1] = right_coef;
        memcpy(src_p1, srcptr+1, sizeof(float)*stride_minus_1);
        v4sf *srcp_p1 = (v4sf*) src_p1, *srcp_m1 = (v4sf*) src_m1;
        
        for(i=0;i<iterline;i++){
            *dstp = coeff[0]*(*srcp_m1) + coeff[1]*(*srcp) + coeff[2]*(*srcp_p1);
            dstp+=1; srcp_m1+=1; srcp+=1; srcp_p1+=1;
        }
    }
    free(src_p1);
    free(src_m1);
}

static void convolve_horiz_fast_5(image_t *dst, const image_t *src, const convolution_t *conv){
    const int stride_minus_1 = src->stride-1;
    const int stride_minus_2 = src->stride-2;
    const int iterline = (src->stride>>2);
    const float *coeff = conv->coeffs;
    v4sf *srcp = (v4sf*) src->data, *dstp = (v4sf*) dst->data;
    float *src_p1 = (float*) malloc(sizeof(float)*src->stride*4);
    float *src_p2 = src_p1+src->stride;
    float *src_m1 = src_p2+src->stride;
    float *src_m2 = src_m1+src->stride;
    int j;
    for(j=0;j<src->height;j++){
        int i;
        float *srcptr = (float*) srcp;
        const float right_coef = srcptr[src->width-1];
        for(i=src->width;i<src->stride;i++)
            srcptr[i] = right_coef;
        src_m1[0] = srcptr[0];
        memcpy(src_m1+1, srcptr , sizeof(float)*stride_minus_1);
        src_m2[0] = srcptr[0];
        src_m2[1] = srcptr[0];
        memcpy(src_m2+2, srcptr , sizeof(float)*stride_minus_2);
        src_p1[stride_minus_1] = right_coef;
        memcpy(src_p1, srcptr+1, sizeof(float)*stride_minus_1);
        src_p2[stride_minus_1] = right_coef;
        src_p2[stride_minus_2] = right_coef;
        memcpy(src_p2, srcptr+2, sizeof(float)*stride_minus_2);
                
        v4sf *srcp_p1 = (v4sf*) src_p1, *srcp_p2 = (v4sf*) src_p2, *srcp_m1 = (v4sf*) src_m1, *srcp_m2 = (v4sf*) src_m2;
        
        for(i=0;i<iterline;i++){
            *dstp = coeff[0]*(*srcp_m2) + coeff[1]*(*srcp_m1) + coeff[2]*(*srcp) + coeff[3]*(*srcp_p1) + coeff[4]*(*srcp_p2);
            dstp+=1; srcp_m2 +=1; srcp_m1+=1; srcp+=1; srcp_p1+=1; srcp_p2+=1;
        }
    }
    free(src_p1);
}

/* perform an horizontal convolution of an image */
void convolve_horiz(image_t *dest, const image_t *src, const convolution_t *conv){
    if(conv->order==1){
        convolve_horiz_fast_3(dest,src,conv);
        return;
    }else if(conv->order==2){
        convolve_horiz_fast_5(dest,src,conv);
        return;    
    }
    float *in = src->data;
    float * out = dest->data;
    int i, j, ii;
    float *o = out;
    int i0 = -conv->order;
    int i1 = +conv->order;
    float *coeff = conv->coeffs + conv->order;
    float *coeff_accu = conv->coeffs_accu + conv->order;
    for(j = 0; j < src->height; j++){
        const float *al = in + j * src->stride;
        const float *f0 = coeff + i0;
        float sum;
        for(i = 0; i < -i0; i++){
	        sum=coeff_accu[-i - 1] * al[0];
	        for(ii = i1 + i; ii >= 0; ii--){
	            sum += coeff[ii - i] * al[ii];
            }
	        *o++ = sum;
        }
        for(; i < src->width - i1; i++){
	        sum = 0;
	        for(ii = i1 - i0; ii >= 0; ii--){
	            sum += f0[ii] * al[ii];
            }
	        al++;
	        *o++ = sum;
        }
        for(; i < src->width; i++){
	        sum = coeff_accu[src->width - i] * al[src->width - i0 - 1 - i];
	        for(ii = src->width - i0 - 1 - i; ii >= 0; ii--){
	            sum += f0[ii] * al[ii];
            }
	        al++;
	        *o++ = sum;
        }
        for(i = 0; i < src->stride - src->width; i++){
	        o++;
        }
    }
}

/* perform a vertical convolution of an image */
void convolve_vert(image_t *dest, const image_t *src, const convolution_t *conv){
    if(conv->order==1){
        convolve_vert_fast_3(dest,src,conv);
        return;
    }else if(conv->order==2){
        convolve_vert_fast_5(dest,src,conv);
        return;    
    }
    float *in = src->data;
    float *out = dest->data;
    int i0 = -conv->order;
    int i1 = +conv->order;
    float *coeff = conv->coeffs + conv->order;
    float *coeff_accu = conv->coeffs_accu + conv->order;
    int i, j, ii;
    float *o = out;
    const float *alast = in + src->stride * (src->height - 1);
    const float *f0 = coeff + i0;
    for(i = 0; i < -i0; i++){
        float fa = coeff_accu[-i - 1];
        const float *al = in + i * src->stride;
        for(j = 0; j < src->width; j++){
	        float sum = fa * in[j];
	        for(ii = -i; ii <= i1; ii++){
	            sum += coeff[ii] * al[j + ii * src->stride];
            }
	        *o++ = sum;
        }
        for(j = 0; j < src->stride - src->width; j++) 
	    {
	        o++;
        }
    }
    for(; i < src->height - i1; i++){
        const float *al = in + (i + i0) * src->stride;
        for(j = 0; j < src->width; j++){
	        float sum = 0;
	        const float *al2 = al;
	        for(ii = 0; ii <= i1 - i0; ii++){
	            sum += f0[ii] * al2[0];
	            al2 += src->stride;
            }
	        *o++ = sum;
	        al++;
        }
        for(j = 0; j < src->stride - src->width; j++){
	        o++;
        }
    }
    for(;i < src->height; i++){
        float fa = coeff_accu[src->height - i];
        const float *al = in + i * src->stride;
        for(j = 0; j < src->width; j++){
	        float sum = fa * alast[j];
	        for(ii = i0; ii <= src->height - 1 - i; ii++){
	            sum += coeff[ii] * al[j + ii * src->stride];
            }
	        *o++ = sum;
        }
        for(j = 0; j < src->stride - src->width; j++){
	        o++;
        }
    }
}

/* free memory of a convolution structure */
void convolution_delete(convolution_t *conv){
    if(conv)
    {
        free(conv->coeffs);
        free(conv->coeffs_accu);
        free(conv);
    }
}

/* perform horizontal and/or vertical convolution to a color image */
void color_image_convolve_hv(color_image_t *dst, const color_image_t *src, const convolution_t *horiz_conv, const convolution_t *vert_conv){
    const int width = src->width, height = src->height, stride = src->stride;
    // separate channels of images
    image_t src_red = {width,height,stride,src->c1}, src_green = {width,height,stride,src->c2}, src_blue = {width,height,stride,src->c3}, 
            dst_red = {width,height,stride,dst->c1}, dst_green = {width,height,stride,dst->c2}, dst_blue = {width,height,stride,dst->c3};
    // horizontal and vertical
    if(horiz_conv != NULL && vert_conv != NULL){
        float *tmp_data = malloc(sizeof(float)*stride*height);
        if(tmp_data == NULL){
	        fprintf(stderr,"error color_image_convolve_hv(): not enough memory\n");
	        exit(1);
        }  
        image_t tmp = {width,height,stride,tmp_data};   
        // perform convolution for each channel
        convolve_horiz(&tmp,&src_red,horiz_conv); 
        convolve_vert(&dst_red,&tmp,vert_conv); 
        convolve_horiz(&tmp,&src_green,horiz_conv);
        convolve_vert(&dst_green,&tmp,vert_conv); 
        convolve_horiz(&tmp,&src_blue,horiz_conv); 
        convolve_vert(&dst_blue,&tmp,vert_conv);
        free(tmp_data);
    }else if(horiz_conv != NULL && vert_conv == NULL){ // only horizontal
        convolve_horiz(&dst_red,&src_red,horiz_conv);
        convolve_horiz(&dst_green,&src_green,horiz_conv);
        convolve_horiz(&dst_blue,&src_blue,horiz_conv);
    }else if(vert_conv != NULL && horiz_conv == NULL){ // only vertical
        convolve_vert(&dst_red,&src_red,vert_conv);
        convolve_vert(&dst_green,&src_green,vert_conv);
        convolve_vert(&dst_blue,&src_blue,vert_conv);
    }
}

/************ Others **********/

float pow2( float f ) {return f*f;}
/* return a new image in lab color space */
color_image_t *rgb_to_lab(const color_image_t *im){

    color_image_t *res = color_image_new(im->width, im->height);
    const int npix = im->stride*im->height;

    const float T=0.008856;
    const float color_attenuation = 1.5f;
    int i;
    for(i=0 ; i<npix ; i++){
        const float r = im->c1[i]/255.f;
        const float g = im->c2[i]/255.f;
        const float b = im->c3[i]/255.f;
        float X=0.412453 * r + 0.357580 * g + 0.180423 * b;
        float Y=0.212671 * r + 0.715160 * g + 0.072169 * b;
        float Z=0.019334 * r + 0.119193 * g + 0.950227 * b;
        X/=0.950456;
        Z/=1.088754;
        float Y3 = pow(Y,1./3);
        float fX = X>T ? pow(X,1./3) : 7.787 * X + 16/116.;
        float fY = Y>T ? Y3 : 7.787 * Y + 16/116.;
        float fZ = Z>T ? pow(Z,1./3) : 7.787 * Z + 16/116.;
        float L = Y>T ? 116 * Y3 - 16.0 : 903.3 * Y;
        float A = 500 * (fX - fY);
        float B = 200 * (fY - fZ);
        // correct L*a*b*: dark area or light area have less reliable colors
        float correct_lab = exp(-color_attenuation*pow2(pow2(L/100) - 0.6)); 
        res->c1[i] = L;
        res->c2[i] = A*correct_lab;
        res->c3[i] = B*correct_lab;
    }
    return res;

}   

/* compute the saliency of a given image */
image_t* saliency(const color_image_t *im, float sigma_image, float sigma_matrix ){
    int width = im->width, height = im->height, filter_size;
    
    // smooth image
    color_image_t *sim = color_image_new(width, height);
    float *presmooth_filter = gaussian_filter(sigma_image, &filter_size);
    convolution_t *presmoothing = convolution_new(filter_size, presmooth_filter, 1);
    color_image_convolve_hv(sim, im, presmoothing, presmoothing);
    convolution_delete(presmoothing);
    free(presmooth_filter);
  
    // compute derivatives
    float deriv_filter[2] = {0.0f, -0.5f};
    convolution_t *deriv = convolution_new(1, deriv_filter, 0);
    color_image_t *imx = color_image_new(width, height), *imy = color_image_new(width, height);
    color_image_convolve_hv(imx, sim, deriv, NULL);
    color_image_convolve_hv(imy, sim, NULL, deriv);
    convolution_delete(deriv);
    
    // compute autocorrelation matrix
    image_t *imxx = image_new(width, height), *imxy = image_new(width, height), *imyy = image_new(width, height);
    v4sf *imx1p = (v4sf*) imx->c1, *imx2p = (v4sf*) imx->c2, *imx3p = (v4sf*) imx->c3, *imy1p = (v4sf*) imy->c1, *imy2p = (v4sf*) imy->c2, *imy3p = (v4sf*) imy->c3, 
        *imxxp = (v4sf*) imxx->data, *imxyp = (v4sf*) imxy->data, *imyyp = (v4sf*) imyy->data;
    int i;
    for(i = 0 ; i<height*im->stride/4 ; i++){
        *imxxp = (*imx1p)*(*imx1p) + (*imx2p)*(*imx2p) + (*imx3p)*(*imx3p);
        *imxyp = (*imx1p)*(*imy1p) + (*imx2p)*(*imy2p) + (*imx3p)*(*imy3p);
        *imyyp = (*imy1p)*(*imy1p) + (*imy2p)*(*imy2p) + (*imy3p)*(*imy3p);
        imxxp+=1; imxyp+=1; imyyp+=1;
        imx1p+=1; imx2p+=1; imx3p+=1;
        imy1p+=1; imy2p+=1; imy3p+=1;
    }
    
    // integrate autocorrelation matrix
    float *smooth_filter = gaussian_filter(sigma_matrix, &filter_size);
    convolution_t *smoothing = convolution_new(filter_size, smooth_filter, 1);
    image_t *tmp = image_new(width, height);
    convolve_horiz(tmp, imxx, smoothing);
    convolve_vert(imxx, tmp, smoothing);
    convolve_horiz(tmp, imxy, smoothing);
    convolve_vert(imxy, tmp, smoothing);    
    convolve_horiz(tmp, imyy, smoothing);
    convolve_vert(imyy, tmp, smoothing);    
    convolution_delete(smoothing);
    free(smooth_filter);
    
    // compute smallest eigenvalue
    v4sf vzeros = {0.0f,0.0f,0.0f,0.0f};
    v4sf vhalf = {0.5f,0.5f,0.5f,0.5f};
    v4sf *tmpp = (v4sf*) tmp->data;
    imxxp = (v4sf*) imxx->data; imxyp = (v4sf*) imxy->data; imyyp = (v4sf*) imyy->data;
    for(i = 0 ; i<height*im->stride/4 ; i++){
        (*tmpp) = vhalf*( (*imxxp)+(*imyyp) ) ;
        (*tmpp) = __builtin_ia32_sqrtps(__builtin_ia32_maxps(vzeros, (*tmpp) - __builtin_ia32_sqrtps(__builtin_ia32_maxps(vzeros, (*tmpp)*(*tmpp) + (*imxyp)*(*imxyp) - (*imxxp)*(*imyyp) ) )));
        tmpp+=1; imxyp+=1; imxxp+=1; imyyp+=1;
    }
    
    image_delete(imxx); image_delete(imxy); image_delete(imyy);
    color_image_delete(imx); color_image_delete(imy);
    color_image_delete(sim);
    
    return tmp;
}
