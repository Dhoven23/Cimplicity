#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool InitData(DataHandle_t *handle, unsigned length){
    DataHandle_t p_handle;

    p_handle = malloc(length * (sizeof *p_handle));
    
    if (p_handle != NULL){
        printf("%lu bytes allocated from %p to %p\n",length * (sizeof *p_handle), p_handle, (p_handle+length));
        *handle = p_handle;
        return true;
    } else {
        return false;
    }
}


void FreeData(DataHandle_t handle){

    free(handle);
}

bool PopulateData(DataHandle_t handle, unsigned length){
    // If memory is allocated, populate with random values
    if (handle != NULL){
        for (int i = 0; i < length; ++i){
            handle[i].data = ((double)rand() / (double)RAND_MAX);
        }

        return true;

    } else {

        return false;
    }   
}

void PrintData(DataHandle_t handle, unsigned length){
    if (handle != NULL){
        for (int i = 0; i < length; ++i)
        {
            printf("data[%i] = ,%f,\n",i,handle[i].data);
        }
    }
}

void SmoothDataInMemory(DataHandle_t handle, unsigned length, unsigned iterations){
    for (unsigned i = 0; i < iterations; ++i)
    {
        for (unsigned j = 1; j < (length-1); ++j)
        {
            handle[j].data *= .9;
            handle[j].data += .1 * (handle[j-1].data + handle[j+1].data) / 2;
        }
    }
}

void SmoothDataLocal(DataHandle_t handle, unsigned length, unsigned iterations){
    double arr[length];
    for (int i = 0; i < length; ++i)
    {
        arr[i] = handle[i].data;
    }

    for (unsigned i = 0; i < iterations; ++i)
    {
        for (unsigned j = 1; j < (length-1); ++j)
        {
            arr[j] *= .9;
            arr[j] += .1 * (arr[j-1] + arr[j+1]) / 2;
        }
    }

    for (int i = 0; i < length; ++i)
    {
        handle[i].data = arr[i];
    }

}
void SmoothData2D(DataHandle_t handle, unsigned length, unsigned iterations){
    double arr[length];
    for (int i = 0; i < length; ++i)
    {
        arr[i] = handle[i].data;
    }
    for (unsigned i = 0; i < iterations; ++i)
    {
        for (unsigned j = 0; j < length; ++j)
        {
            int X = j % 17;
            int Y = (int)(floor(j / 17));

            if ((X != 0 && Y != 0) && (X < 16 && Y < 16)){
                arr[j] *= .9;
                arr[j] += .1 * (arr[j-1] + arr[j+1] + arr[j-17] + arr[j+17]) / 4;
            }
        }
    }
    for (int i = 0; i < length; ++i)
    {
        handle[i].data = arr[i];
    }
}

 void DerivativeEuler(DataHandle_t* out_handle, DataHandle_t handle, unsigned length){
    ((*out_handle))->data = 0.0f;
    for (unsigned i = 1; i < length; ++i){

        ((*out_handle)+i)->data = (double)(handle[i].data - handle[i-1].data) * (double)length;
    }
 }

 void DerivativeTrapezoidal(DataHandle_t* out_handle, DataHandle_t handle, unsigned length){
    ((*out_handle))->data = 0.0f;
    for (unsigned i = 1; i < (length - 1); ++i)
    {
        ((*out_handle)+i)->data = (handle[i+1].data - handle[i-1].data) * .5f * (double)length;
    }
 }

 int GetIndex(int X, int Y){
    return (X + Y*17);
 }

 void DeleteGridNodes2D(DataHandle_t data_handle, GridNodeHandle_t grid_handle, unsigned length, double threshold){

    // Local storage memory
    double data_arr[length];
    //int grid_coords[length][2];
    //uint32_t grid_neighbors[length];
    //void* address_arr[length];

    for (int i = 0; i < length; ++i){
        data_arr[i]         = data_handle[i].data;
        //grid_coords[i][0]   = grid_handle[i].coordinate[0];
        //grid_coords[i][1]   = grid_handle[i].coordinate[1];
        //grid_neighbors[i]   = grid_handle[i].neighbors;
        //address_arr[i]      = grid_handle[i].DataPointer;
    }

    char* buff;
    buff = malloc(17*17* (sizeof *buff));

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
           

            // center scalar
            
            r[4]  = data_arr[ GetIndex(X0,Y0) ]; 

            // ------- lookup point scalars --------

                    
            // pivot scalars
            double p[4];
            p[0] = data_arr[GetIndex(X0  ,Y0+1)];
            p[1] = data_arr[GetIndex(X0-1,Y0  )];
            p[2] = data_arr[GetIndex(X0  ,Y0-1)];
            p[3] = data_arr[GetIndex(X0+1,Y0  )];


            // calculate marching cross configuration. 
            int Xcheck = j;
            int Yckeck = i;

            int cross = ((Xcheck == 0) && (Yckeck == 0)) ? 0 : 0;
                cross = ((Xcheck == 0) && (Yckeck != 0)) ? 1 : cross;
                cross = ((Xcheck != 0) && (Yckeck == 0)) ? 2 : cross;
                cross = ((Xcheck != 0) && (Yckeck != 0)) ? 3 : cross;

            // Evaluate marching cross derivatives

            double deltaSum = 0;



            // Printing buffer

        
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

                        buff[X0-1     ]  = '*';
                        buff[X0       ]  = '*';
                        buff[X0+1     ]  = '*';
                        buff[X0-1 + 17]  = '*';
                        buff[X0   + 17]  = '*';
                        buff[X0+1 + 17]  = '*';
                        buff[X0-1 + 34]  = '*';
                        buff[X0   + 34]  = '*';
                        buff[X0+1 + 34]  = '*';

                    } else {

                        buff[X0-1     ]  = '*';
                        buff[X0       ]  = ' ';
                        buff[X0+1     ]  = '*';
                        buff[X0-1 + 17]  = ' ';
                        buff[X0   + 17]  = '*';
                        buff[X0+1 + 17]  = ' ';
                        buff[X0-1 + 34]  = '*';
                        buff[X0   + 34]  = ' ';
                        buff[X0+1 + 34]  = '*';
                    }
                    printf(" %.2f ",deltaSum);
                    
                    break;

                case 1: // compute marching cross for case 1


                    // ------- calculate pivot deltas _______
                    
                    deltaSum  = fabs((((r[1] - r[0]) / 2) + r[0]) - p[0]);
                    deltaSum += fabs((((r[3] - r[2]) / 2) + r[2]) - p[2]);
                    deltaSum += fabs((((r[0] - r[3]) / 2) + r[3]) - p[3]);

                    if (deltaSum > threshold){

                        buff[X0-1 + 17*2*i + 17]  = '*';
                        buff[X0   + 17*2*i + 17]  = '*';
                        buff[X0+1 + 17*2*i + 17]  = '*';
                        buff[X0-1 + 17*2*i + 34]  = '*';
                        buff[X0   + 17*2*i + 34]  = '*';
                        buff[X0+1 + 17*2*i + 34]  = '*';


                    } else {

                        buff[X0-1 + 17*2*i + 17]  = ' ';
                        buff[X0   + 17*2*i + 17]  = '*';
                        buff[X0+1 + 17*2*i + 17]  = ' ';
                        buff[X0-1 + 17*2*i + 34]  = '*';
                        buff[X0   + 17*2*i + 34]  = ' ';
                        buff[X0+1 + 17*2*i + 34]  = '*';

                    }
                    printf(" %.2f ",deltaSum);
                    
                    
                    
                    break;

                case 2: // compute marching cross for case 2

                    // ------- calculate pivot deltas _______
                    
                    deltaSum  = fabs((((r[2] - r[1]) / 2) + r[1]) - p[1]);
                    deltaSum += fabs((((r[3] - r[2]) / 2) + r[2]) - p[2]);
                    deltaSum += fabs((((r[0] - r[3]) / 2) + r[3]) - p[3]);


                    if (deltaSum > threshold){


                        buff[X0       ]  = '*';
                        buff[X0+1     ]  = '*';
                        buff[X0   + 17]  = '*';
                        buff[X0+1 + 17]  = '*';
                        buff[X0   + 34]  = '*';
                        buff[X0+1 + 34]  = '*';

                    } else {


                        buff[X0       ]  = ' ';
                        buff[X0+1     ]  = '*';
                        buff[X0   + 17]  = '*';
                        buff[X0+1 + 17]  = ' ';
                        buff[X0   + 34]  = ' ';
                        buff[X0+1 + 34]  = '*';
                    }
                    printf(" %.2f ",deltaSum);

                    break;

                case 3: // compute marching cross for case 3

                    // ------- calculate pivot deltas _______

                    deltaSum  = fabs((((r[3] - r[2]) / 2) + r[2]) - p[2]);
                    deltaSum += fabs((((r[0] - r[3]) / 2) + r[3]) - p[3]);

                    if (deltaSum > threshold){


                        buff[X0   + 17*i*2]  = '*';
                        buff[X0   + 17+ 17*i*2]  = '*';
                        buff[X0+1 + 17+ 17*i*2]  = '*';
                        buff[X0   + 34+ 17*i*2]  = '*';
                        buff[X0+1 + 34+ 17*i*2]  = '*';
                        buff[X0-1 + 17+ 17*i*2]  = '*';



                    } else {

                        buff[X0   + 17+ 17*i*2]  = '*';
                        buff[X0+1 + 17+ 17*i*2]  = ' ';
                        buff[X0   + 34+ 17*i*2]  = ' ';
                        buff[X0+1 + 34+ 17*i*2]  = '*';

                    }
                    printf(" %.2f ",deltaSum);

                    break;
            }
        }  
        printf("\n");     
    }



        for (int j = 0; j < 17; ++j){
            for (int i = 0; i < 17; ++i){
                printf("%c ",buff[i + j*17]);
            }
            printf("\n");
        }
    free(buff);
 }

 

