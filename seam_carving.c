#include "seamcarving.h"
#include "c_img.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MIN(x,y) (((x) < (y)) ? (x) : (y))


void calc_energy(struct rgb_img *im, struct rgb_img **grad){

    //Allocating memory to store seam energies
    create_img(grad, im->height, im->width);

    //Initializing
    int rx;
    int gx;
    int bx;
    int ry;
    int gy;
    int by;

    int dx2;
    int dy2;

    int pixel_energy;
    int gradient_energy;

    for(int y = 0; y<im->height; y++){
        for(int x = 0; x<im->width; x++){
            
            //LEFT OF IMAGE
            if(x == 0){
                rx = get_pixel(im, y, (im->width)-1, 0) - get_pixel(im, y, x+1, 0);
                gx = get_pixel(im, y, (im->width)-1, 1) - get_pixel(im, y, x+1, 1);
                bx = get_pixel(im, y, (im->width)-1, 2) - get_pixel(im, y, x+1, 2);
            //RIGHT OF IMAGE
            }else if(x == (im->width)-1){
                rx = get_pixel(im, y, x-1, 0) - get_pixel(im, y, 0, 0);
                gx = get_pixel(im, y, x-1, 1) - get_pixel(im, y, 0, 1);
                bx = get_pixel(im, y, x-1, 2) - get_pixel(im, y, 0, 2);
            }else{
                rx = get_pixel(im, y, x-1, 0) - get_pixel(im, y, x+1, 0);
                gx = get_pixel(im, y, x-1, 1) - get_pixel(im, y, x+1, 1);
                bx = get_pixel(im, y, x-1, 2) - get_pixel(im, y, x+1, 2);
            }

            //TOP OF IMAGE
            if(y==0){
                ry = get_pixel(im, (im->height)-1, x, 0) - get_pixel(im, y+1, x, 0);
                gy = get_pixel(im, (im->height)-1, x, 1) - get_pixel(im, y+1, x, 1);
                by = get_pixel(im, (im->height)-1, x, 2) - get_pixel(im, y+1, x, 2);
            //BOTTOM OF IMAGE
            }else if(y== (im->height)-1){
                ry = get_pixel(im, y-1, x, 0) - get_pixel(im, 0, x, 0);
                gy = get_pixel(im, y-1, x, 1) - get_pixel(im, 0, x, 1);
                by = get_pixel(im, y-1, x, 2) - get_pixel(im, 0, x, 2);
            }else{
                ry = get_pixel(im, y-1, x, 0) - get_pixel(im, y+1, x, 0);
                gy = get_pixel(im, y-1, x, 1) - get_pixel(im, y+1, x, 1);
                by = get_pixel(im, y-1, x, 2) - get_pixel(im, y+1, x, 2);
            }

            dx2 = (rx*rx) + (gx*gx) + (bx*bx);
            dy2 = (ry*ry) + (gy*gy) + (by*by);

            pixel_energy = sqrt(dx2 + dy2);

            gradient_energy = (uint8_t)(pixel_energy/10);

            set_pixel(*grad, y, x, gradient_energy, gradient_energy, gradient_energy);
        }
        }
}

void dynamic_seam(struct rgb_img *grad, double **best_arr){

    *best_arr = malloc((grad->height) * (grad->width) * sizeof(double) + 1);

    int width = grad->width;
    int height = grad->height;


    for(int i = 0; i<width; i++){
        (*best_arr)[i] = grad->raster[3*i];
    }

    for(int y = 1; y<height; y++){
        for(int x = 0; x<width; x++){
            if(x == 0){
                (*best_arr)[y*width+x] = (double)grad->raster[3 * (y*(grad->width) + x)] + (double)MIN((*best_arr)[(y-1)*width+x], (*best_arr)[(y-1)*width+(x+1)]);
            }else if(x == (width-1)){
                (*best_arr)[y*width+x] = (double)grad->raster[3 * (y*(grad->width) + x)] + (double)MIN((*best_arr)[(y-1)*width+x], (*best_arr)[(y-1)*width+(x-1)]);
            }else{
                double min = MIN((*best_arr)[(y-1)*width+x], (*best_arr)[(y-1)*width+(x-1)]);
                (*best_arr)[y*width+x] = (double)grad->raster[3 * (y*(grad->width) + x)] + (double)MIN(min, (*best_arr)[(y-1)*width+(x+1)]);   
            }
        }
    }

}

void recover_path(double *best, int height, int width, int **path){
    double min = 10000000000;
    int y = height;
    int min_ind = -1;

    *path = malloc(y*sizeof(int)+1);

    for (int i = 0; i<width; i++){
        if (best[(y-1)*width+i] < min){
            min = best[(y-1)*width+i];
            min_ind = i;
        }
    }

    (*path)[height-1] = min_ind;
    for (int i = y-1; i> 0; i--){
        //Find the min of the three/two values that are above it
        
        if (min_ind == 0){
            min = MIN(best[(i-1)*width], best[((i-1)*width) +1]); 
            if (min == best[(i-1)*width]){
                min_ind = 0;
            }else{
                min_ind = 1;
            }
        }else if (min_ind == width -1){
            min = MIN(best[(i-1)*width + min_ind], best[(i-1)*width + min_ind -1]);
            if (min == best[(i-1)*width +min_ind]){
                min_ind = width - 1;
            }else{
                min_ind = width - 2;
            }
        }else{
            double min_temp = MIN(best[(i-1)*width +min_ind], best[(i-1)*width + (min_ind - 1)]);
            min = MIN(best[(i-1)*width + (min_ind +1)], min_temp);
            if (min == best[(i-1)*width + min_ind -1]){
                min_ind = min_ind - 1;
            }if (min == best[(i-1)*width +(min_ind)]){
                min_ind = min_ind;
            }else{
                min_ind = min_ind +1;
            }
        }
        (*path)[i-1] = min_ind;
    }
}

void remove_seam(struct rgb_img *src, struct rgb_img **dest, int *path){

    create_img(dest, src->height, src->width-1);

    //Initialize
    int r;
    int g;
    int b;

    for (int y = 0; y<src->height;y++){
        for(int x = 0; x<src->width; x++){
            if(x == path[y]){
                continue;
            }else{
                r = get_pixel(src, y, x, 0);
                g = get_pixel(src, y, x, 1);
                b = get_pixel(src, y, x, 2); 
                if (x < path[y]){
                    set_pixel(*dest, y, x, r, g, b);
                }else{
                    set_pixel(*dest, y, x-1, r, g, b);
                }
            }
        }
    }

}
