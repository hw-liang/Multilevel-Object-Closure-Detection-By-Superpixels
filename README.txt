General info:
-------------

This code computes multiple closure hypotheses in an image as described in :
"Optimal Contour Closure by Superpixel Grouping", Alex Levinshtein, Cristian Sminchisescu, Sven Dickinson. ECCV 2010.

When using this code please reference the above paper.

Our implementation included code from other sources:
1) Pb and globalPb from University of Berkeley
2) Normalized cuts
3) Contour fragment code from Andrew Stein
4) Parametric Maxflow from Vladimir Kolmogorov et al
5) Matlab helper functions from Peter Kovesi
6) MatlabBGL library for graph operations from David Gleich

We acknowledge the authors of the above software.
Our code, as well as the code included above is covered under GPL (see LICENSE)


Quick startup guide:
--------------------

1) Run the command 'DefinePath' from the Closure directory in Matlab
2) Run make.m (if you haven't done so before on this system)
3) Use the function 'ClosureMain' to compute closure (see header for usage)

The code was tested on Matlab 7.8 on Linux 32bit, Linux 64bit, and Windows.

Note:

1) For computational efficiency, we save the preprocessing results of the various stages (Pb, globalPb, superpixels). 
However, this can cause unexpexted behaviour. For example, if superpixel were computed with Pb, but then one desides to use 
globalPb, the superpixels won't be recomputed since the old superpixels are already stored.

2) The by default the code uses Pb and not globalPb to compute contours.
globalPb can be used by setting the use_gpb = true flag. We provide a compiled version of globalPb for Matlab 7.8
for 32 and 64 bit Linux. For a different platform the user will need to compile the globalPb code manually and put it under the globalPb
directory. The code can be downloaded from : http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/
 
Example:
--------
Run the following from matlab after completing the setup steps above:

>> ClosureMain('airplane.png', '.', 100, 10, 0.05, 'turbo');

This will run closure detection on the image 'airplane.png' in the current directory and store the results in the
current directory. This example uses 100 superpixels, returns at most 10 solutions, and uses an edge threshold (T_e)
of 0.05. The turbopixels method will be used to generate superpixels.

This will generate airplane_solution_???.jpg files, where each one contains a binary mask for one closure hypothesis.
The solutions are ordered by increasing closure cost.
We also generate airplane_multiplesolutions.jpg file the visualizes all the solutions in one image.
