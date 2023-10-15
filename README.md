# Cimplicity
Marching Cross Algorithm - (2dWaveletMesh)

TODO: speed up neighbor search. Currently uses brute force search, many ways to improve. 

Idea 1: give the search function an initial index as a guess, neigbors are lilkely not too far from each other. 

Idea 2: cache coordinate / index locations in a const array for fast repeated query. 


The purpose of this algorithm is to adapt simulation meshes for problems with sparse fractal data. It would be very poor for evaluating harmonic functions with regularly spaced data structures.

	Given a regular square grid that has size (N+1) where N is a power of 2:
	(ex 5x5) -> 	* * * * *
            	* * * * *
            	* * * * *
            	* * * * *
            	* * * * *

We wish to determine if the points the grid can be approximated with 
a plane drawn between the corner points. 

    ->  * - - - - - *
        |           |
        |           |
        |           |
        |           |
        * _ _ _ _ _ *
 
We want to know this recursively, 
with the smallest plane we will attempt to draw being:

    (3x3) ->    * * *    * - *
                * * * -> |   |
                * * *    * - *

To begin, evaluate the points (x) lying on lines connected to the corners

	->    * x *
	      x * x
	      * x *
    
	The evaluation is simple, the points () are called the 'ring' because
	 we can imagine them as a mathematical ring where each point is distance 4 
	    from itself. The points themselves are 'points' and each calculation compares the point to an interpolated point:


        
        - > p_interp = p_k + 1/2(r_k+1 - r_k)

    where k+1 is the next point in the ring. The difference
    ABS(p_interp - p_k) is summed for all 4 points, and compared to a threshold value. 

If the points in the 3x3 grid can be well approximated by a plane, we discard them (at least initially, we'll get back to that...)

Returning to the larger grid, the points remaining that have not been evaluated are:

    ->  *   *   * 
          *   *   
        *   *   *
          *   *  
        *   *   *

We now want to know if that center point (x) in each of the former 3x3 grids lies on a cross drawn through the outermost corners:

    - > *   *   *
          x   x
        *   *   *
          x   x
        *   *   *

We can perform this with a series of diagonal 'zigs' and 'zags' through the grid.
(In a larger grid, it is conceptually eazier to loop through the entire grid row by row)

    - > *   *   *       *   *   *
          \   /   . . .   \   /
        *   *   *       *   *   *
          /   \   . . .   /   \
        *   *   *       *   *   *
          .   .
          .   .
          .   .
        *   *   *
          \   /
        *   *   *
          /   \
        *   *   *

The analysis is the same as before, but rather than a ring, we have a staggering sequence of points shifting through a 1D array of 2 points [0,1] where the first point becomes the last with each step:

        (1)    (3) 
          \   /   \  ... 
           (2)    (4)

As we can see, for the first step, we have points (1) and (2) between which to approximate the point (x), and the next step we have (2) and (3). we keep point 2, shift it to the first position and pop in a new point to draw the next line. The evaluation is the same as before, but only a single point is evaluated for each 3x3 grid (see above) so the decision to keep/discard is point by point.

After these steps, we achieve the final grid!

        * * * * *                          *    *    *
        * * * * *
        * * * * *  - (marching cross) - >  *    *    *
        * * * * *
        * * * * *                          *    *    *

This process can be repeated, by starting with the new 3x3 grid (that was 5x5) for arbitrary initial grid sizes, as long as it is size length (N+1) where N is a power of 2. 

The astute will observe however, that much repeated calculation is performed.

For the 5x5 grid, there are 4 3x3 grids, and all 4 share a side. So the interpolation for each side is performed twice! This adds much calculation for very large grids. Furthermore, suppose
we have 2 grids A and B side by side,
 
        * * * * *
        * A x B *
        * * * * *

If B is refined, and A is not, we will end up with

        * * *   *
        * A   B   
        * * *   *

This is not good, as the decision by A to keep the point x was overriden by B
This may not matter strictly for compression, but when we get to drawing geomtric cells around our points, this sort of undefined behavior will increase the number of computations required to define the mesh geometry for the next iteration.

There is a very simple solution to this problem which exploits the symmetry of the oddnumbered geometry we are using. 

We can look at the index of the point A and B, and determine if we are at the top left corner, top edge, left side, or body of the grid. We will give each of these cases a number:

        Case (0) - > Top left Corner X * * * 
                                     *
                                     *
        
        Case (1) -> Left Edge   * * *
                                X * *
                                * * *

        Case (2) -> Top edge   * X *
                               * * *
                               * * *

        Case (3) -> otherwise  * * *
                               * X *
                               * * *

(In these illustrations, points are 3x3 grids, not grid points)

If we are in Case 0, the evaluation is as described above.

If we are in Case 1, we only need to evaluate the lower 2 rows:
        
            . . .
        ->  x * x
            * x *

	As the top row was already evaluated by the prior loop iteration. 

If we are in Case 2, we only need to evaluate the right 2 columns:

        -> . x *
           . * x
           . x *

And in case 3, the most common case, we only have 2 evaluations:

        -> . . *
           . * x
           * x *

This method not only reduces the number of calculations dramatically 
for large grids where N^2 >> 2N, but it also guarantews much more predictable behavior 
(we have directionality for our algorithm).

One more step however is required to avoid the 'disputed node' problem. for a grid 
of size 17x17 (the example chosen for this code) we divide the the grid into 64 zones.
Each zone is given one bit of a 64 bit number to represent if it has been kept, or discarded,
and we simply lookup the these logical values for the neighboring cells to see if a side is
needed before dicarding it. This does increase the number of total cells in the grid, 
but in our grid, we care about capturing sparse fractal data (As set out at the begining) 
There's more to the guts of this algorithm, but that should be enough for you to get your bearings.

Good luck! 