/*  -------------------------------------------------------------------------

    Author: Daniel Hoven 
    Date: 10/14/2023

    email: Daniel.Hoven@Honeywell.com

    (This project is in no way affiliated with Honeywell inc.)

    The purpose of this algorithm is to adapt simulation meshes for problems with sparse fractal
    data. It would be very poor for evaluating harmonic functions with regularly spaced
    data structures. 

    This module is responsible for Generating, Adapting, conmpressing, 
    and decompressing data stored in the MeshHandle_t data type. It makes
    use of the Marching Crosses algorithm described below. 

    Given a regular square grid that has size (N+1) where N is a power of 2:

    (ex 5x5) -> * * * * *
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
    
    The evaluation is simple, the points (*) are called the 'ring' because
    we can imagine them as a mathematical ring where each point is distance 4 
    from itself. The points themselves are 'points' and each calculation compares 
    the point to an interpolated point:

        - > p_interp = p_k + 1/2(r_k+1 - r_k)

    where k+1 is the next point in the ring. The difference
    ABS(p_interp - p_k) is summed for all 4 points, and compared to a threshold value. 

    If the points in the 3x3 grid can be well approximated by a plane, we discard them 
    (at least initially, we'll get back to that...)

    Returning to the larger grid, the points remaining that have not been evaluated are:

    ->  *   *   * 
          *   *   
        *   *   *
          *   *  
        *   *   *

    We now want to know if that center point (x) in each of the former 3x3 grids lies on a cross
    drawn through the outermost corners (*). 

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

        The analysis is the same as before, but rather than a ring, we have a staggering sequence
        of points shifting through a 1D array of 2 points [0,1] where the first point becomes the 
        last with each step:

        (1)    (3) 
          \   /   \  ... 
           (2)    (4)

    As we can see, for the first step, we have points (1) and (2) between which to approximate the 
    point (x), and the next step we have (2) and (3). we keep point 2, shift it to the first position
    and pop in a new point to draw the next line. The evaluation is the same as before, but only a single
    point is evaluated for each 3x3 grid (see above) so the decision to keep/discard is point by point.

        After these steps, we achieve the final grid!

        * * * * *                          *    *    *
        * * * * *
        * * * * *  - (marching cross) - >  *    *    *
        * * * * *
        * * * * *                          *    *    *

        This process can be repeated, by starting with the new 3x3 grid (that was 5x5) for arbitrary initial
        grid sizes, as long as it is size length (N+1) where N is a power of 2. 

        The astute will observe however, that much repeated calculation is performed.

        For the 5x5 grid, there are 4 3x3 grids, and all 4 share a side. So the interpolation for 
        each side is performed twice! This adds much calculation for very large grids. Furthermore, suppose
        we have 2 grids A and B side by side,
         
        * * * * *
        * A x B *
        * * * * *

        If B is refined, and A is not, we will end up with

        * * *   *
        * A   B   
        * * *   *

    This is not good, as the decision by A to keep the point x was overriden by B
    This may not matter strictly for compression, but when we get to drawing geomtric cells
    around our points, this sort of undefined behavior will increase the number of computations
    required to define the mesh geometry for the next iteration.

    There is a very simple solution to this problem which exploits the symmetry of the odd
    numbered geometry we are using. 

    We can look at the index of the point A and B, and determine if we are at the top left corner,
    top edge, left side, or body of the grid. We will give each of these cases a number

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

        (In these illustrations, '*'s are 3x3 grids, not actual points)
        
        If we are in Case 0, the evaluation is as described above.

        if we are in Case 1, we only need to evaluate the lower 2 rows:
        
            . . .
        ->  x * x
            * x *

        As the top row was already evaluated by the prior loop iteration. 

        if we are in Case 2, we only need to evaluate the right 2 columns:

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
*/


#include "header.h"
#include "DynamicMesh.h"

// Public functions

bool GenerateMesh(MeshHandle_t handle){

    if (Indexer_Create(&(handle->Indexer),LENGTH) && InitData(&(handle->DataField),LENGTH)){
        PopulateData(handle->DataField,LENGTH);
        SmoothData2D(handle->DataField,LENGTH,10);

        // set bitfields
        handle->ZeroLevel = UINT64_MAX;
        handle->N1Level    = UINT16_MAX;
        handle->N2Level   = 15u;
        handle->N3Level   = 1u;

        for (int i = 0; i < LENGTH; ++i)
        {
            handle->Indexer[i].data_ptr = (void*)(&(handle->DataField[i]));
            DataHandle_t p_data;
        }



        return true;
    } else {
        return false;
    }
}

void DestroyMesh(MeshHandle_t handle){
    Indexer_Destroy(handle->Indexer);
    Data_Destroy(handle->DataField);
}

void SetMeshThreshold(double p_threshold, MeshHandle_t handle){
    handle->threshold = p_threshold;
}

void getMeshThreshold(double* p_threshold, MeshHandle_t handle) {
    *p_threshold = (double)handle->threshold;
}

void AdaptMesh(MeshHandle_t handle){
    RefineMesh(handle);
    TransposeMesh(handle);
}

bool GetScalarByCoordinate(int x, int y, MeshHandle_t handle, double* data){
    
    IndexHandle_t p_index;
    DataHandle_t p_data;
    double data_out;
    Indexer_GetIndexByCoordinate(x,y,(handle->Indexer),&p_index);
    bool b_IsValidData = false, b_IsValidPointer = false;

    b_IsValidPointer = GetData(p_index->data_ptr,&p_data);
    b_IsValidData    = GetScalar(p_data,&data_out);

    if (b_IsValidData && b_IsValidPointer){
        *data = data_out;
        return true;
    }
    return false;
}

// Grid Logic tools

bool ZeroLevelIsEmpty(int x, int y, MeshHandle_t handle){
    int ind = x + (y*8);

    if ((ind<64) && !(ind<0) && !(x<0) && !(y<0) && (x<8) && (y<8)) {
        uint64_t mask = (uint64_t)(1ull << ind);
        return !((handle->ZeroLevel & mask) == mask);
    }
    return true;
}

bool N1LevelIsEmpty(int x, int y, MeshHandle_t handle){
    int ind = x + (y*4);
    if ((ind<16) && !(ind<0) && !(x<0) && !(y<0) && (x<4) && (y<4)){
        uint16_t mask = (uint16_t)(1 << ind);
        return !((handle->N1Level & mask) == mask);
    }
    return true;
}

bool N2LevelIsEmpty(int x, int y, MeshHandle_t handle){
    int ind = x + (y*2);
    if ((ind < 4)&& !(ind < 0) && !(x<0) && !(y<0) && (x<2) && (y<2)){
        uint8_t mask = (uint8_t)(1 << ind);
        return !((handle->N2Level & mask) == mask);
    }
    return true;
}

bool ZeroLevelCrossIsEmpty(int x, int y, MeshHandle_t handle){
    bool p_bool = false;
    p_bool |= (x > 0) ? !ZeroLevelIsEmpty(x-1,y,handle) : false;
    p_bool |= (x < 7) ? !ZeroLevelIsEmpty(x+1,y,handle) : false;
    p_bool |= (y > 0) ? !ZeroLevelIsEmpty(x,y-1,handle) : false;
    p_bool |= (y < 7) ? !ZeroLevelIsEmpty(x,y+1,handle) : false;
    p_bool |= !ZeroLevelIsEmpty(x,y,handle);

    return !p_bool;
}

bool N1LevelCrossIsEmpty(int x, int y, MeshHandle_t handle){
    bool p_bool = false;

    p_bool |= (x > 0) ? !N1LevelIsEmpty(x-1,y,handle) : false;
    p_bool |= (x < 3) ? !N1LevelIsEmpty(x+1,y,handle) : false;
    p_bool |= (y > 0) ? !N1LevelIsEmpty(x,y-1,handle) : false;
    p_bool |= (y < 3) ? !N1LevelIsEmpty(x,y+1,handle) : false;
    p_bool |= !N1LevelIsEmpty(x,y,handle);

    return !p_bool;
}

bool N2LevelCrossIsEmpty(int x, int y, MeshHandle_t handle){
    bool p_bool = false;

    p_bool |= (x == 0) ? !N2LevelIsEmpty(x+1,y,handle) : false;
    p_bool |= (x == 1) ? !N2LevelIsEmpty(x-1,y,handle) : false;
    p_bool |= (y == 1) ? !N2LevelIsEmpty(x,y-1,handle) : false;
    p_bool |= (y == 0) ? !N2LevelIsEmpty(x,y+1,handle) : false;
    p_bool |= !N2LevelIsEmpty(x,y,handle);

    return !p_bool;
}

bool ZeroLevelSquareIsEmpty(int x, int y, MeshHandle_t handle){
    bool p_bool = false;
    for (int i = -1; i < 3; ++i){
        for (int j = -1; j < 3; ++j){
            p_bool |= !ZeroLevelIsEmpty(2*x+j,2*y+i,handle);
        }
    }   return !p_bool;
}

bool N1LevelSquareIsEmpty(int x, int y, MeshHandle_t handle){
    bool p_bool = false;
    for (int i = -1; i < 3; ++i){
        for (int j = -1; j < 3; ++j){
            p_bool |= !N1LevelIsEmpty(2*x+j,2*y+i,handle);
        }
    }   return !p_bool;
}


/*/-------------------------------------------------------
    MESH REFINEMENT (SERIOUSLY DO NOT TOUCH)
/ /---------------------------------------------------- */

static void MarchingCrossTypeZero(MeshHandle_t handle,char*buff){

    /* ---------------------------------
        This is the initial marching cross
        algorithm. It requires the Indexer
        field of the mesh to be coalesced
    */

    // get Threshold
    double threshold = handle->threshold;

    // get data 
    double data_arr[LENGTH];
    int count = 0;
    while(count < (LENGTH)){
        DataHandle_t data_temp;
        data_temp = (DataHandle_t)(handle->Indexer[count].data_ptr);
        data_arr[count] = data_temp->data;
        count++;
    }

    double deltaSum = 0;

    for (int i = 0; i < 8; ++i){
        for (int j = 0; j < 8; ++j){
             // center point
            int X0 = 1 + (j * 2);
            int Y0 = 1 + (i * 2);

            // ------- lookup ring scalars ---------

            double r[5];
            r[0] = data_arr[GetIndex(X0+1,Y0-1)]; // check if r0 is in bounds
            r[1] = data_arr[GetIndex(X0-1,Y0-1)];
            r[2] = data_arr[GetIndex(X0-1,Y0+1)];
            r[3] = data_arr[GetIndex(X0+1,Y0+1)];
            

            // ------- lookup pivot scalars ________
            double p[4];
            p[0] = data_arr[GetIndex(X0  ,Y0+1)];
            p[1] = data_arr[GetIndex(X0-1,Y0  )];
            p[2] = data_arr[GetIndex(X0  ,Y0-1)];
            p[3] = data_arr[GetIndex(X0+1,Y0  )];

            // compute cross

            // calculate marching cross configuration. 
            int Xcheck = j;
            int Yckeck = i;

            int cross = ((Xcheck == 0) && (Yckeck == 0)) ? 0 : 0;
                cross = ((Xcheck == 0) && (Yckeck != 0)) ? 1 : cross;
                cross = ((Xcheck != 0) && (Yckeck == 0)) ? 2 : cross;
                cross = ((Xcheck != 0) && (Yckeck != 0)) ? 3 : cross;

            switch (cross)
            {
                case 0: // compute marching cross for case 0
                    
                    // ------- calculate pivot deltas _______
                    
                    for (int i = 0; i < 4; ++i)
                    {
                        int r1 = i;
                        int r2 = (i+1) % 4;

                        deltaSum += fabs((((r[r2] - r[r1]) / 2) + r[r1]) - p[i]);
                    }

                    if (deltaSum > threshold){

                        buff[GetIndex(X0-1,0)]  = '*';
                        buff[GetIndex(X0  ,0)]  = '*';
                        buff[GetIndex(X0+1,0)]  = '*';
                        buff[GetIndex(X0-1,1)]  = '*';
                        buff[GetIndex(X0  ,1)]  = '*';
                        buff[GetIndex(X0+1,1)]  = '*';
                        buff[GetIndex(X0-1,2)]  = '*';
                        buff[GetIndex(X0  ,2)]  = '*';
                        buff[GetIndex(X0+1,2)]  = '*';
                        KeepZoneZeroLevel(j,i,handle);
                    } else {

                        buff[GetIndex(X0-1,0)]  = '*';
                        buff[GetIndex(X0  ,0)]  = ' ';
                        buff[GetIndex(X0+1,0)]  = '*';
                        buff[GetIndex(X0-1,1)]  = ' ';
                        buff[GetIndex(X0  ,1)]  = '*';
                        buff[GetIndex(X0+1,1)]  = ' ';
                        buff[GetIndex(X0-1,2)]  = '*';
                        buff[GetIndex(X0  ,2)]  = ' ';
                        buff[GetIndex(X0+1,2)]  = '*';
                        DeleteZoneZeroLevel(j,i,handle);
                    } break;

                case 1: // compute marching cross for case 1


                    // ------- calculate pivot deltas _______
                    
                    deltaSum  = fabs((((r[1] - r[0]) / 2) + r[0]) - p[0]);
                    deltaSum += fabs((((r[3] - r[2]) / 2) + r[2]) - p[2]);
                    deltaSum += fabs((((r[0] - r[3]) / 2) + r[3]) - p[3]);

                    if (deltaSum > threshold){

                        buff[GetIndex(X0-1,Y0  )]  = '*';
                        buff[GetIndex(X0-1,Y0+1)]  = '*';
                        buff[GetIndex(X0  ,Y0  )]  = '*';
                        buff[GetIndex(X0  ,Y0+1)]  = '*';
                        buff[GetIndex(X0+1,Y0  )]  = '*';
                        buff[GetIndex(X0+1,Y0+1)]  = '*';
                        buff[GetIndex(X0  ,Y0-1)]  = '*';
                        KeepZoneZeroLevel(j,i,handle);
                    } else {

                        buff[GetIndex(X0-1,Y0  )]  = ' ';
                        buff[GetIndex(X0-1,Y0+1)]  = '*';
                        buff[GetIndex(X0  ,Y0  )]  = '*';
                        buff[GetIndex(X0  ,Y0+1)]  = ' ';
                        buff[GetIndex(X0+1,Y0  )]  = ' ';
                        buff[GetIndex(X0+1,Y0+1)]  = '*';
                        DeleteZoneZeroLevel(j,i,handle);
                    } break;

                case 2: // compute marching cross for case 2

                    // ------- calculate pivot deltas _______
                    
                    deltaSum  = fabs((((r[2] - r[1]) / 2) + r[1]) - p[1]);
                    deltaSum += fabs((((r[3] - r[2]) / 2) + r[2]) - p[2]);
                    deltaSum += fabs((((r[0] - r[3]) / 2) + r[3]) - p[3]);


                    if (deltaSum > threshold){

                        buff[GetIndex(X0  ,0)]  = '*';
                        buff[GetIndex(X0+1,0)]  = '*';
                        buff[GetIndex(X0  ,1)]  = '*';
                        buff[GetIndex(X0+1,1)]  = '*';
                        buff[GetIndex(X0  ,2)]  = '*';
                        buff[GetIndex(X0+1,2)]  = '*';
                        buff[GetIndex(X0-1,Y0)] = '*';
                        KeepZoneZeroLevel(j,i,handle);
                    } else {

                        buff[GetIndex(X0  ,0)]  = ' ';
                        buff[GetIndex(X0+1,0)]  = '*';
                        buff[GetIndex(X0  ,1)]  = '*';
                        buff[GetIndex(X0+1,1)]  = ' ';
                        buff[GetIndex(X0  ,2)]  = ' ';
                        buff[GetIndex(X0+1,2)]  = '*';
                        DeleteZoneZeroLevel(j,i,handle);
                    } break;

                case 3: // compute marching cross for case 3

                    // ------- calculate pivot deltas _______

                    deltaSum  = fabs((((r[3] - r[2]) / 2) + r[2]) - p[2]);
                    deltaSum += fabs((((r[0] - r[3]) / 2) + r[3]) - p[3]);

                    if (deltaSum > threshold){


                        buff[GetIndex(X0  ,Y0  )]  = '*';
                        buff[GetIndex(X0+1,Y0  )]  = '*';
                        buff[GetIndex(X0  ,Y0+1)]  = '*';
                        buff[GetIndex(X0+1,Y0+1)]  = '*';
                        buff[GetIndex(X0  ,Y0-1)]  = '*';
                        buff[GetIndex(X0-1,Y0  )]  = '*';
                        KeepZoneZeroLevel(j,i,handle);
                    } else {

                        buff[GetIndex(X0  ,Y0  )]  = '*';
                        buff[GetIndex(X0+1,Y0  )]  = ' ';
                        buff[GetIndex(X0  ,Y0+1)]  = ' ';
                        buff[GetIndex(X0+1,Y0+1)]  = '*';
                        DeleteZoneZeroLevel(j,i,handle);
                    } break;
            }
        }
    }
}

static void MarchingCrossTypeN(MeshHandle_t handle, int Nlevel, char* buff){

    // get Threshold
    double threshold = handle->threshold;

    // get data 
    double data_arr[LENGTH];
    int count = 0;
    double deltaSum = 0;

    int mult = (int)(1 << (unsigned)Nlevel);

    while(count < (17*17)){
        DataHandle_t data_temp;
        data_temp = (DataHandle_t)(handle->Indexer[count].data_ptr);
        data_arr[count] = data_temp->data;
        
        //printf("Index[%i] = %f\n",count, data[count]);
        count++;
    }

    for (int i = 0; i < (8/mult); ++i){
        for (int j = 0; j < (8/mult); ++j){
             // center point
            int X0 = mult + (j * 2*mult);
            int Y0 = mult + (i * 2*mult);

            // ------- lookup ring scalars ---------

            double r[5];
            r[0] = data_arr[GetIndex(X0+mult,Y0-mult)]; // check if r0 is in bounds
            r[1] = data_arr[GetIndex(X0-mult,Y0-mult)];
            r[2] = data_arr[GetIndex(X0-mult,Y0+mult)];
            r[3] = data_arr[GetIndex(X0+mult,Y0+mult)];
            

            // ------- lookup pivot scalars ________
            double p[4];
            p[0] = data_arr[GetIndex(X0  ,Y0+mult)];
            p[1] = data_arr[GetIndex(X0-mult,Y0  )];
            p[2] = data_arr[GetIndex(X0  ,Y0-mult)];
            p[3] = data_arr[GetIndex(X0+mult,Y0  )];

            // compute cross

            // calculate marching cross configuration. 
            int Xcheck = j;
            int Yckeck = i;

            int cross = ((Xcheck == 0) && (Yckeck == 0)) ? 0 : 0;
                cross = ((Xcheck == 0) && (Yckeck != 0)) ? 1 : cross;
                cross = ((Xcheck != 0) && (Yckeck == 0)) ? 2 : cross;
                cross = ((Xcheck != 0) && (Yckeck != 0)) ? 3 : cross;

            switch (cross) {

            case 0:
                for (int i = 0; i < 4; ++i){
                    int r1 = i;
                    int r2 = (i+1) % 4;

                    deltaSum += fabs((((r[r2] - r[r1]) / 2) + r[r1]) - p[i]);
                }
                
                if(deltaSum > threshold){
                    switch(Nlevel){
                    case 1:
                        KeepZoneN1Level(j,i,handle);
                        break;
                    case 2:
                        KeepZoneN2Level(j,i,handle);
                    }
                    

                } else {
                    bool p_bool;
                    switch (Nlevel){
                    case 1:
                        p_bool = ZeroLevelSquareIsEmpty(j,i,handle);
                        if (p_bool) DeleteZoneN1Level(j,i,handle);
                        break;
                    case 2:
                        p_bool = N1LevelSquareIsEmpty(j,i,handle);
                        if (p_bool) DeleteZoneN2Level(j,i,handle);
                        
                    }

                    if (p_bool){
                        buff[GetIndex(mult,0)] = ' ';
                        buff[GetIndex(mult,2*mult)] = ' ';
                        buff[GetIndex(0,mult)] = ' ';
                        buff[GetIndex(2*mult,mult)] = ' ';
                    }
                } break;

            case 1:

                 // ------- calculate pivot deltas _______
                    
                deltaSum  = fabs((((r[1] - r[0]) / 2) + r[0]) - p[0]);
                deltaSum += fabs((((r[3] - r[2]) / 2) + r[2]) - p[2]);
                deltaSum += fabs((((r[0] - r[3]) / 2) + r[3]) - p[3]);
                
                if(deltaSum > threshold){

                    switch(Nlevel){
                        case 1:
                            KeepZoneN1Level(j,i,handle);
                    break;
                        case 2:
                            KeepZoneN2Level(j,i,handle);
                    break;
                    }
                    buff[GetIndex(2,i*4)] = '*';
                } else {
                    bool p_bool;
                    switch (Nlevel){
                    case 1:
                        p_bool = ZeroLevelSquareIsEmpty(j,i,handle);
                        if (p_bool) DeleteZoneN1Level(j,i,handle);
                        break;
                    case 2:
                        p_bool = N1LevelSquareIsEmpty(j,i,handle);
                        if (p_bool) DeleteZoneN2Level(j,i,handle);
                        
                    }
                    
                    if (p_bool){
                        DeleteZoneN1Level(j,i,handle);
                        buff[GetIndex(0,Y0)] = ' ';
                        buff[GetIndex(2*mult,Y0)] = ' ';
                        buff[GetIndex(mult,Y0+mult)] = ' ';
                    } else {
                        buff[GetIndex(mult,Y0-mult)] = '*';
                    }
                } break;

            case 2:

                 // ------- calculate pivot deltas _______
                    
                deltaSum  = fabs((((r[2] - r[1]) / 2) + r[1]) - p[1]);
                deltaSum += fabs((((r[3] - r[2]) / 2) + r[2]) - p[2]);
                deltaSum += fabs((((r[0] - r[3]) / 2) + r[3]) - p[3]);

                
                if(deltaSum > threshold){
                    switch(Nlevel){
                        case 1:
                            KeepZoneN1Level(j,i,handle);
                    break;
                        case 2:
                            KeepZoneN2Level(j,i,handle);
                    break;
                    }
                    buff[GetIndex(X0-2,2)] = '*';

                } else {
                    bool p_bool;
                    switch (Nlevel){
                    case 1:
                        p_bool = ZeroLevelSquareIsEmpty(j,i,handle);
                        if (p_bool) DeleteZoneN1Level(j,i,handle);
                        break;
                    case 2:
                        p_bool = N1LevelSquareIsEmpty(j,i,handle);
                        if (p_bool) DeleteZoneN2Level(j,i,handle);
                        
                    }
                    
                    if (p_bool){
                        DeleteZoneN1Level(j,i,handle);
                        buff[GetIndex(X0,0)] = ' ';
                        buff[GetIndex(X0,2*mult)] = ' ';
                        buff[GetIndex(X0+mult,mult)] = ' ';
                       
                       
                    } else {
                        buff[GetIndex(X0-mult,mult)] = '*';
                    }
                } break;

            case 3:

                // ------- calculate pivot deltas _______

                deltaSum  = fabs((((r[3] - r[2]) / 2) + r[2]) - p[2]);
                deltaSum += fabs((((r[0] - r[3]) / 2) + r[3]) - p[3]);
                
                if(deltaSum > threshold){
                    switch(Nlevel){
                        case 1:
                            KeepZoneN1Level(j,i,handle);
                    break;
                        case 2:
                            KeepZoneN2Level(j,i,handle);
                    break;
                    }
                    buff[GetIndex(X0  ,Y0-2)] = '*';
                    buff[GetIndex(X0-2,Y0  )] = '*';
                } else {
                    bool p_bool;
                    switch (Nlevel){
                    case 1:
                        p_bool = ZeroLevelSquareIsEmpty(j,i,handle);
                        if (p_bool) DeleteZoneN1Level(j,i,handle);
                        break;
                    case 2:
                        p_bool = N1LevelSquareIsEmpty(j,i,handle);
                        if (p_bool) DeleteZoneN2Level(j,i,handle);
                        
                    }
                    
                    if (p_bool){
                        DeleteZoneN1Level(j,i,handle);
                        buff[GetIndex(X0  ,Y0+mult)] = ' ';

                        buff[GetIndex(X0+mult,Y0  )] = ' ';

                       
                    } else {
                        buff[GetIndex(X0  ,Y0-mult)] = '*';
                        buff[GetIndex(X0-mult,Y0)] = '*';
                    }
                } break;
            }
        }
    }
}

static void DiagonalMarch(MeshHandle_t handle, int Nlevel, char* buff){
    // quick input check
    if ((Nlevel > 2) || (Nlevel < 0)) return;

    // calculate local values
    int mult = (int)(1 << (unsigned)Nlevel);
    printf("mult = %i\n",mult);
    int bound = 8 / (mult);
    int c = -mult;

    // get data 
    double data_arr[LENGTH];
    int count = 0;

    while(count < (LENGTH)){
        DataHandle_t data_temp;
        data_temp = (DataHandle_t)(handle->Indexer[count].data_ptr);
        data_arr[count] = data_temp->data;
        count++;
    }

    // local point array
    double points[2];
    points[1] = data_arr[0];

    for (int j = 0; j < bound; ++j){
        for (int i = 0; i < bound; ++i) {
            int X0 = (i * (2*mult)) + mult;
            int Y0 = (j * (2*mult)) + mult;
            
            c *= -1;
            int Y_step = Y0 + c;
            int X_step = X0 + mult;


            double center_point = data_arr[GetIndex(X0,Y0)];
            double new_point = data_arr[GetIndex(X_step,Y_step)];

            points[0] = points[1];
            points[1] = new_point;

            double delta = fabs((((points[1] - points[0]) / 2) + points[0]) - center_point);

            if (delta > handle->threshold){
                buff[GetIndex(X0,Y0)] = '*';

            } else {
                bool p_bool;
                switch (Nlevel) {
                case 0:
                    p_bool = ZeroLevelCrossIsEmpty(i,j,handle);
                    break;
                case 1:
                    p_bool = N1LevelCrossIsEmpty(i,j,handle);
                    break;
                case 2:
                    p_bool = N2LevelCrossIsEmpty(i,j,handle);
                }

                if (p_bool){
                    buff[GetIndex(X0,Y0)] = ' ';
                } else {
                    buff[GetIndex(X0,Y0)] = '*';
                }
                

            }
        }
        if (j % 2 == 0){
            points[0] = data_arr[GetIndex(0,(j*(2*mult))+(4*mult))];
            c = mult;
        } else {
            points[0] = data_arr[GetIndex(0,((j-1)*(2*mult))+(4*mult))];
            c = -mult;
        }
    }
}

static void RefineMesh(MeshHandle_t handle){
    char* buff = malloc(LENGTH * sizeof (*buff));

    MarchingCrossTypeZero(handle,buff);
    DiagonalMarch(handle,0,buff);
    MarchingCrossTypeN(handle,1,buff);
    DiagonalMarch(handle,1,buff);
    MarchingCrossTypeN(handle,2,buff);
    DiagonalMarch(handle,2,buff);


    // print mesh
    int count = 0;
    for (int j = 0; j < 17; ++j){
        for (int i = 0; i < 17; ++i){
            printf("%c ",buff[i + j*17]);
            if (buff[i + j*17] == '*') count++;
        }
        printf("\n");
    }
    printf("\nCompression ratio: %.2f%%\n\n",(double)(((double)count*100.0f) / (double)(LENGTH)));
}


/*/-------------------------------------------------------
    DATA RECOALESCION (BE CAREFUL!!)
/ /---------------------------------------------------- */
static void TransposeMesh(MeshHandle_t handle){};


/*/-------------------------------------------------------
    ZONE DELETION (BE CAREFUL!!)
/ /---------------------------------------------------- */

static void KeepZoneZeroLevel(int x, int y, MeshHandle_t handle){
    if ((x < 8) && (y < 8) && (x >= 0) && (y >= 0)){
        int ind = x + 8*y;
        handle->ZeroLevel |= ((uint64_t)1 << ind);
    }
}

static void DeleteZoneZeroLevel(int x, int y, MeshHandle_t handle){
    if ((x < 8) && (y < 8) && (x >= 0) && (y >= 0)){
        int ind = x + 8*y;
        handle->ZeroLevel &= ~((uint64_t)1 << ind);
    }
}

static void KeepZoneN1Level(int x, int y, MeshHandle_t handle){
    if ((x < 4) && (y < 4) && (x >= 0) && (y >= 0)){
        int ind = x + 4*y;
        handle->N1Level |= ((uint16_t)1 << ind);
    }
}

static void DeleteZoneN1Level(int x, int y, MeshHandle_t handle){
    if ((x < 4) && (y < 4) && (x >= 0) && (y >= 0)){
        int ind = x + 4*y;
        handle->N1Level &= ~((uint16_t)1 << ind);
    }
}

static void KeepZoneN2Level(int x, int y, MeshHandle_t handle){
    if ((x < 2) && (y < 2) && (x >= 0) && (y >= 0)){
        int ind = x + 2*y;
        handle->N2Level |= ((uint16_t)1 << ind);
    }
}

static void DeleteZoneN2Level(int x, int y, MeshHandle_t handle){
    if ((x < 2) && (y < 2) && (x >= 0) && (y >= 0)){
        int ind = x + 2*y;
        handle->N2Level &= ~((uint16_t)1 << ind);
    }
}