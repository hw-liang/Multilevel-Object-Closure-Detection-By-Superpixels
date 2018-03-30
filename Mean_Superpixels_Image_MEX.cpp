#include <stdio.h>
#include <math.h>
#include "mex.h"

double max(double num1, double num2) {
    return num1 > num2 ? num1: num2;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
		const mxArray *prhs[])
{
	double *seg_image, *image;
	double *mean_superpixel;
    double *count_pixels_per_sup;
    
	int num_nodes = 0;
	
	seg_image = mxGetPr(prhs[0]);
	image = mxGetPr(prhs[1]);

    
    int M = (int)mxGetM(prhs[0]);
    int N = (int)mxGetN(prhs[0]);

    for (int y = 0; y < N; y++) {
        for (int x = 0; x < M; x++) {
            num_nodes = (int)(max((double)num_nodes, *(seg_image + y * M + x)));
        }
    }
    
//     mexPrintf("Number of nodes: %d \n", num_nodes);

    plhs[0] = mxCreateDoubleMatrix(num_nodes, 1, mxREAL);
    mean_superpixel = mxGetPr(plhs[0]);
    count_pixels_per_sup = new double[num_nodes];
    
    for (int i=0; i<num_nodes; i++) {
        count_pixels_per_sup[i] = 0;
        mean_superpixel[i] = 0;
    }
    
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < M; x++) {
            int cur_sup = (int)(*(seg_image + y * M + x)) - 1;
            
            if (cur_sup >= 0) {
                mean_superpixel[cur_sup] += *(image + y * M + x);
                count_pixels_per_sup[cur_sup]++;
            }
        }
    }
    
    for (int i=0; i<num_nodes; i++) {
        mean_superpixel[i] /= count_pixels_per_sup[i];
    }
    
    delete count_pixels_per_sup;
}
