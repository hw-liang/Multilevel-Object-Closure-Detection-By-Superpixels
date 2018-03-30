#include <stdio.h>
#include <math.h>
#include "mex.h"

double max(double num1, double num2) {
    return num1 > num2 ? num1: num2;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
		const mxArray *prhs[])
{
	double *seg_image, *sup_value;
	double *sup_value_image;
    
	int num_nodes = 0;
	
	seg_image = mxGetPr(prhs[0]);  // superpixel labels for MxN
	sup_value = mxGetPr(prhs[1]);  // selected superpixel labels

    
    int M = (int)mxGetM(prhs[0]);  // number of rows
    int N = (int)mxGetN(prhs[0]);  // number of columns

    // get the maximum number in the nodes (highest number of superpixel)
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < M; x++) {
            num_nodes = (int)(max((double)num_nodes, *(seg_image + y * M + x)));
        }
    }
    
//     mexPrintf("Number of nodes: %d \n", num_nodes);

    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);  // create a 2d array
    sup_value_image = mxGetPr(plhs[0]);
    
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < M; x++) {
            // for each pixel location
            int cur_sup = (int)(*(seg_image + y * M + x)) - 1; // current superpixel label - 1
            *(sup_value_image + y * M + x) = 0;  // make the current pixel label as 0
            if (cur_sup >= 0) {
                // If cur_sup is non-negative, then assign the corresponding selected label value.
                *(sup_value_image + y * M + x) = sup_value[cur_sup];
            }
        }
    }
}
