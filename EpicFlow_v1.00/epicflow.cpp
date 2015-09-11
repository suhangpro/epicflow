#include <stdlib.h>
#include <string.h>
#include "epic.h"
#include "image.h"
#include "io.h"
#include "variational.h"


/* show usage information */
void usage(){
    printf("usage:\n");
    printf("    ./epicflow image1 image2 edges matches outputfile [options]\n");
    printf("Compute EpicFlow between two images using given matches and edges and store it into a .flo file\n");
    printf("Images must be in PPM, JPG or PNG format.\n");
    printf("Edges are read as width*height float32 values in a binary file\n");
    printf("Matches are read from a text file, each match in a different line, each line starting with 4 numbers corresponding to x1 y1 x2 y2\n");
    printf("\n");
    printf("options:\n"); 
    printf("    -h, -help                                                print this message\n");
    printf("  interpolation parameters\n");
    printf("    -nw                                                      use Nadaraya-Watson instead of LA interpolator in the interpolation\n");
    printf("    -p, -prefnn             <int>(25)                        number of neighbors for consisteny checking in the interpolation\n");
    printf("    -n, -nn                 <int>(100)                       number of neighnors for the interpolation\n");
    printf("    -k                      <float>(0.8)                     coefficient of the sigmoid of the Gaussian kernel used in the interpolation\n");
    printf("  energy minimization parameters\n");
    printf("    -i, -iter               <int>(5)                         number of iterations for the energy minimization\n");
    printf("    -a, -alpha              <float>(1.0)                     weight of smoothness term\n");
    printf("    -g, -gamma              <float>(3.0)                     weight of gradient constancy assumption\n");
    printf("    -d, -delta              <float>(2.0)                     weight of color constancy assumption\n");
    printf("    -s, -sigma              <float>(0.8)                     standard deviation of Gaussian presmoothing kernel\n");
    printf("  predefined parameters\n");
    printf("    -sintel                                                  set the parameters to the one optimized on (a subset of) the MPI-Sintel dataset\n");
    printf("    -middlebury                                              set the parameters to the one optimized on the Middlebury dataset\n");
    printf("    -kitti                                                   set the parameters to the one optimized on the KITTI dataset\n");
    printf("\n");
}


int main(int argc, char **argv){
    if( argc<6){
        if(argc>1) fprintf(stderr,"Error, not enough arguments\n");
        usage();
        exit(1);
    }

    // read arguments
    color_image_t *im1 = color_image_load(argv[1]);
    color_image_t *im2 = color_image_load(argv[2]);
    float_image edges = read_edges(argv[3], im1->width, im1->height);
    float_image matches = read_matches(argv[4]);
    const char *outputfile = argv[5];

    // prepare variables
    epic_params_t epic_params;
    epic_params_default(&epic_params);
    variational_params_t flow_params;
    variational_params_default(&flow_params);
    image_t *wx = image_new(im1->width, im1->height), *wy = image_new(im1->width, im1->height);
    
    // read optional arguments 
    #define isarg(key)  !strcmp(a,key)
    int current_arg = 6;
    while(current_arg < argc ){
        const char* a = argv[current_arg++];
        if( isarg("-h") || isarg("-help") ) 
            usage();
        else if( isarg("-nw") ) 
            strcpy(epic_params.method, "NW");
        else if( isarg("-p") || isarg("-prefnn") ) 
            epic_params.pref_nn = atoi(argv[current_arg++]);
        else if( isarg("-n") || isarg("-nn") ) 
            epic_params.nn = atoi(argv[current_arg++]); 
        else if( isarg("-k") ) 
            epic_params.coef_kernel = atof(argv[current_arg++]);
        else if( isarg("-i") || isarg("-iter") ) 
            flow_params.niter_outer = atoi(argv[current_arg++]); 
        else if( isarg("-a") || isarg("-alpha") ) 
            flow_params.alpha= atof(argv[current_arg++]);  
        else if( isarg("-g") || isarg("-gamma") ) 
            flow_params.gamma= atof(argv[current_arg++]);                                  
        else if( isarg("-d") || isarg("-delta") ) 
            flow_params.delta= atof(argv[current_arg++]);  
        else if( isarg("-s") || isarg("-sigma") ) 
            flow_params.sigma= atof(argv[current_arg++]); 
        else if( isarg("-sintel") ){ 
            epic_params.pref_nn= 25; 
            epic_params.nn= 160; 
            epic_params.coef_kernel = 1.1f; 
            flow_params.niter_outer = 5;
            flow_params.alpha = 1.0f;
            flow_params.gamma = 0.72f;
            flow_params.delta = 0.0f;
            flow_params.sigma = 1.1f;            
        }  
        else if( isarg("-kitti") ){ 
            epic_params.pref_nn= 25; 
            epic_params.nn= 160; 
            epic_params.coef_kernel = 1.1f;
            flow_params.niter_outer = 2;
            flow_params.alpha = 1.0f;
            flow_params.gamma = 0.77f;
            flow_params.delta = 0.0f;
            flow_params.sigma = 1.7f; 
        }
        else if( isarg("-middlebury") ){ 
            epic_params.pref_nn= 15; 
            epic_params.nn= 65; 
            epic_params.coef_kernel = 0.2f;       
            flow_params.niter_outer = 25;
            flow_params.alpha = 1.0f;
            flow_params.gamma = 0.72f;
            flow_params.delta = 0.0f;
            flow_params.sigma = 1.1f;  
        }
        else{
            fprintf(stderr, "unknown argument %s", a);
            usage();
            exit(1);
        }   
    }
    
    // compute interpolation and energy minimization
    color_image_t *imlab = rgb_to_lab(im1);
    epic(wx, wy, imlab, &matches, &edges, &epic_params, 1);
    // energy minimization
    variational(wx, wy, im1, im2, &flow_params);
    // write output file and free memory
    writeFlowFile(outputfile, wx, wy);
    
    color_image_delete(im1);
    color_image_delete(imlab);
    color_image_delete(im2);
    free(matches.pixels);
    free(edges.pixels);
    image_delete(wx);
    image_delete(wy);

    return 0;
}
