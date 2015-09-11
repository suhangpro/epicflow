#include <stdlib.h>
#include <stdio.h>
#include <jpeglib.h>
#include <assert.h>
#include <setjmp.h>
#include <png.h>
#include "io.h"

/*********** EDGES and MATCHES ***********/

/* read edges from a binary file containing width*height float32 values */
float_image read_edges(const char *filename, const int width, const int height){
    float_image res = empty_image(float, width, height);
    FILE *fid = fopen(filename, "rb");
    assert(fread(res.pixels, sizeof(float), width*height, fid)==width*height);
    fclose(fid);
    return res;   
}

/* read matches, stored as x1 y1 x2 y2 per line (other values on the same is not taking into account */
float_image read_matches(const char *filename){
    float_image res = empty_image(float, 4, 100000);
    FILE *fid = fopen(filename, "r");
    int nmatch = 0;
    float x1, x2, y1, y2;
	while(!feof(fid) && fscanf(fid, "%f %f %f %f%*[^\n]", &x1, &y1, &x2, &y2)==4){
	    res.pixels[4*nmatch  ] = x1;
	    res.pixels[4*nmatch+1] = y1;
	    res.pixels[4*nmatch+2] = x2;
	    res.pixels[4*nmatch+3] = y2;	    
	    nmatch++;
    }
    res.pixels = (float*) realloc(res.pixels, sizeof(float)*4*nmatch);
    res.ty = nmatch;
    fclose(fid);
    return res;
}


/******* FLOW ********/

/* read a flow file and returns a pointer with two images containing the flow along x and y axis */
image_t** readFlowFile(const char* filename){
    FILE *fid = fopen(filename, "rb");
    if (fid == 0){
        fprintf(stderr,"readFlow() error: could not open file  %s\n",filename);
        exit(1);
    }
    float help;
    fread(&help,sizeof(float),1,fid);
    int aXSize,aYSize;
    fread(&aXSize,sizeof(int),1,fid);
    fread(&aYSize,sizeof(int),1,fid);
    image_t** flow = (image_t**) malloc(sizeof(image_t*)*2);
    flow[0] = image_new(aXSize, aYSize);
    flow[1] = image_new(aXSize, aYSize);
    int x,y;
    for (y = 0; y < aYSize; y++)
        for (x = 0; x < aXSize ; x++) {
            fread(&(flow[0]->data[y*flow[0]->stride+x]),sizeof(float),1,fid);
            fread(&(flow[1]->data[y*flow[0]->stride+x]),sizeof(float),1,fid);
    }
    fclose(fid);
    return flow;
}

/* write a flow to a file */
void writeFlowFile(const char *filename, const image_t *flowx, const image_t *flowy){
    FILE *stream = fopen(filename, "wb");
    if (stream == 0){
        fprintf(stderr, "Error while opening %s\n",filename);
        exit(1);
    }
    const float help=202021.25;
    fwrite(&help,sizeof(float),1,stream);
    const int aXSize = flowx->width, aYSize = flowx->height;
    fwrite(&aXSize,sizeof(int),1,stream);
    fwrite(&aYSize,sizeof(int),1,stream);
    int y,x;
    for (y = 0; y < aYSize ; y++)
        for (x = 0; x < aXSize ; x++) {
	        fwrite(&flowx->data[y*flowx->stride+x],sizeof(float),1,stream);
	        fwrite(&flowy->data[y*flowy->stride+x],sizeof(float),1,stream);
        }
    fclose(stream);
}

/********************* IMAGE ***********************/

// PPM

typedef struct{
    int magic;
    int width;
    int height;
    int pixmax;
} ppm_hdr_t;

static void get_magic(FILE *fp, ppm_hdr_t *ppm_hdr){
    char str[1024];
    fgets(str, 1024, fp);
    if(str[0] == 'P' && (str[1] <= '6' || str[1] >= '1')){
        ppm_hdr->magic = str[1] - '0';
    }
}

static int skip_comment(FILE *fp){
    char c;
    do{
        c = (char) fgetc(fp);
    } 
	while (c == ' ' || c == '\t' || c == '\n');
    if(c == '#'){
        do {
            c = (char) fgetc(fp);

        } while(c != 0x0A);
        return 1;
    }else{
        ungetc(c, fp);
    }
    return 0;
}

/*----------------------------------------------------------------------------*/

static void skip_comments(FILE *fp){
    while(skip_comment(fp));
}

/*----------------------------------------------------------------------------*/

static int get_image_size(FILE *fp, ppm_hdr_t *ppm_hdr){
    skip_comments(fp);
    if(fscanf(fp, "%d %d", &ppm_hdr->width, &ppm_hdr->height) != 2){
        fprintf(stderr, "Warning: PGM --> File currupted\n");
        return 0;
    }
    return 1;
}

/*----------------------------------------------------------------------------*/

static int get_pixmax(FILE *fp, ppm_hdr_t *ppm_hdr){
    skip_comments(fp);
    ppm_hdr->pixmax = 1;
    if(ppm_hdr->magic == 2 || ppm_hdr->magic == 3 || ppm_hdr->magic == 5 || ppm_hdr->magic == 6){
        if(fscanf(fp, "%d", &ppm_hdr->pixmax) != 1){
            fprintf(stderr, "Warning: PGM --> pixmax not valid\n");
            return 0;
        }
    }
    fgetc(fp);
    return 1;
}

/*----------------------------------------------------------------------------*/

static int get_ppm_hdr(FILE *fp, ppm_hdr_t *ppm_hdr){
    get_magic(fp, ppm_hdr);
    if(!get_image_size(fp, ppm_hdr)){
        return 0;
    }
    if(!get_pixmax(fp, ppm_hdr)){
        return 0;
    }
    return 1;
}

static void raw_read_color(FILE *fp, color_image_t *image){
    int j;
    for( j=0 ; j<image->height ; j++){
        int o = j*image->stride, i;
        for( i=0 ; i<image->width ; i++,o++ ){
            image->c1[o] = (float) fgetc(fp);
            image->c2[o] = (float) fgetc(fp);
            image->c3[o] = (float) fgetc(fp);
        }
    }
}

color_image_t *color_image_pnm_load(FILE *fp){
    color_image_t *image = NULL;
    ppm_hdr_t ppm_hdr;
    if(!get_ppm_hdr(fp, &ppm_hdr)) 	{
        return NULL;
    }
    switch(ppm_hdr.magic)    {
        case 1: /* PBM ASCII */
        case 2: /* PGM ASCII */
        case 3: /* PPM ASCII */
        case 4: /* PBM RAW */
        case 5: /* PGM RAW */
            fprintf(stderr, "color_image_pnm_load: only PPM raw with maxval 255 supported\n");            
            break;
        case 6: /* PPM RAW */
            image = color_image_new(ppm_hdr.width, ppm_hdr.height);
            raw_read_color(fp, image);
            break;
    }
    return image;
}

// JPG

color_image_t *color_image_jpeg_load(FILE *fp){
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    JSAMPARRAY buffer;
    int row_stride;
    int index = 0;
    color_image_t *image = NULL;
    float *r_p, *g_p, *b_p;
    JSAMPROW buffer_p;
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, fp);
    jpeg_read_header(&cinfo, TRUE);
    cinfo.out_color_space = JCS_RGB;
    cinfo.quantize_colors = FALSE;
    image = color_image_new(cinfo.image_width, cinfo.image_height);
    if(image == NULL){
        return NULL;
    }
    jpeg_start_decompress(&cinfo);
    row_stride = cinfo.output_width * cinfo.output_components;
    buffer = (*cinfo.mem->alloc_sarray)
        ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

    r_p = image->c1;
    g_p = image->c2;
    b_p = image->c3;

    const int incr_line = image->stride-image->width;

    while (cinfo.output_scanline < cinfo.output_height){
        jpeg_read_scanlines(&cinfo, buffer, 1);
        buffer_p = buffer[0];
        index = cinfo.output_width;
        while(index--){
            *r_p++ = (float) *buffer_p++;
            *g_p++ = (float) *buffer_p++;
            *b_p++ = (float) *buffer_p++;
        }
        r_p += incr_line; g_p += incr_line; b_p += incr_line;
    }
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    return image;
}

// PNG

color_image_t *color_image_png_load( FILE* fp, const char* file_name ){
    // read the header
    png_byte header[8];
    fread(header, 1, 8, fp);
    
    if (png_sig_cmp(header, 0, 8)){
        fprintf(stderr, "error: %s is not a PNG.\n", file_name);
        fclose(fp);
        return 0;
    }
    
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr){
        fprintf(stderr, "error: png_create_read_struct returned 0.\n");
        fclose(fp);
        return 0;
    }
    
    // create png info struct
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr){
        fprintf(stderr, "error: png_create_info_struct returned 0.\n");
        png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
        fclose(fp);
        return 0;
    }
    
    // create png info struct
    png_infop end_info = png_create_info_struct(png_ptr);
    if (!end_info){
        fprintf(stderr, "error: png_create_info_struct returned 0.\n");
        png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);
        fclose(fp);
        return 0;
    }

    // the code in this if statement gets called if libpng encounters an error
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "error from libpng\n");
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
        fclose(fp);
        return 0;
    }

    // init png reading
    png_init_io(png_ptr, fp);

    // let libpng know you already read the first 8 bytes
    png_set_sig_bytes(png_ptr, 8);

    // read all the info up to the image data
    png_read_info(png_ptr, info_ptr);

    // variables to pass to get info
    int bit_depth, color_type;
    png_uint_32 temp_width, temp_height;

    // get info about png
    png_get_IHDR(png_ptr, info_ptr, &temp_width, &temp_height, &bit_depth, &color_type,
        NULL, NULL, NULL);

    // Update the png info struct.
    png_read_update_info(png_ptr, info_ptr);

    // Row size in bytes.
    int rowbytes = png_get_rowbytes(png_ptr, info_ptr);

    // Allocate the image_data as a big block, to be given to opengl
    png_byte * image_data;
    image_data = (png_byte*) malloc(sizeof(png_byte)*rowbytes*temp_height);
    assert(image_data!=NULL);

    // row_pointers is for pointing to image_data for reading the png with libpng
    png_bytep * row_pointers = (png_bytep*) malloc(sizeof(png_bytep)*temp_height);
    assert(row_pointers!=NULL);

    // set the individual row_pointers to point at the correct offsets of image_data
    unsigned int i;
    for (i = 0; i <temp_height; i++)
        row_pointers[i] = image_data + i * rowbytes;

    // read the png into image_data through row_pointers
    png_read_image(png_ptr, row_pointers);
    
    // copy into color image
    color_image_t* image = color_image_new(temp_width,temp_height);
    if( color_type==0 ) {
      assert((unsigned)rowbytes == temp_width || !"error: not a proper gray png image");
      for(i=0; i<temp_height; i++){
        uint j;
        for(j=0; j<temp_width; j++)
            image->c1[i*image->stride+j] = image->c2[i*image->stride+j] = image->c3[i*image->stride+j] = image_data[i*image->width+j];
      }
    } else if( color_type == 2 ) {
      assert((unsigned)rowbytes == 3*temp_width || !"error: not a proper color png image");
      for(i=0; i<temp_height; i++) {
        uint j;
        for(j=0; j<temp_width; j++){
          image->c1[i*image->stride+j] = image_data[3*(i*image->width+j)+0];
          image->c2[i*image->stride+j] = image_data[3*(i*image->width+j)+1];
          image->c3[i*image->stride+j] = image_data[3*(i*image->width+j)+2];
        }
      }
    } else
      assert(!"error: unknown PNG color type" );
    
    // clean up
    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    free(row_pointers);
    free(image_data);
    return image;
}

// GENERAL LOAD

/* load a color image from a file */
color_image_t *color_image_load(const char *fname){
    FILE *fp;
    char magic[2];
    unsigned short *magic_short = (unsigned short *) magic;
    color_image_t *image = NULL;
    if((fp = fopen(fname, "rb")) == NULL){
        fprintf(stderr, "Error in color_image_load() - can not open file `%s' !\n", fname);
        exit(1);
    }
    fread(magic, sizeof(char), 2, fp);
    rewind(fp);
    if(magic_short[0] == 0xd8ff){
        image = color_image_jpeg_load(fp);
    } else if(magic[0]=='P' && (magic[1]=='6' || magic[1]=='5')){ /* PPM raw */
        image = color_image_pnm_load(fp);
    } else if( magic[0]==-119 && magic[1]=='P' ) {
      image = color_image_png_load( fp, fname );
    } else{
        fprintf(stderr, "Error in color_image_load(%s) - image format not supported, can only read jpg or ppm\n",fname);
        exit(1);
    }
    fclose(fp);
    return image;
}
