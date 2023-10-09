#include "header.h"
#include "DynamicMesh.h"


// private functions 

// -------- TODO : IMPLEMENT

void MarchingCrossTypeZero(MeshHandle_t handle){

    /* ---------------------------------
        This is the initial marching cross
        algorithm. It requires the Indexer
        field of the mesh to be coalesced
    */

    // get Threshold
    double threshold = handle->threshold;

    // get data 
    double data[17*17];
    int count = 0;

    while(count < (17*17)){
        DataHandle_t data_temp;
        data_temp = (DataHandle_t)(handle->Indexer[count].data_ptr);
        data[count] = data_temp->data;
        
        printf("Index[%i] = %f\n",count, data[count]);
        count++;
    }
    // get ring scalars
    // get points
    // compute cross
    // evaluate marching cross for case
}

// static void MarchingCrossTypeN(MeshHandle_t handle, int N);

// static void CompressMesh(MeshHandle_t handle);

// static void DecompressMesh(MeshHandle_t handle);

// static void RefineMesh(MeshHandle_t handle);

// static void TransposeMesh(MeshHandle_t handle);


// Public functions

void AdaptMesh(MeshHandle_t handle);

bool GenerateMesh(MeshHandle_t handle){

    if (Indexer_Create(&(handle->Indexer),LENGTH) && InitData(&(handle->DataField),LENGTH)){
        PopulateData(handle->DataField,LENGTH);
        SmoothData2D(handle->DataField,LENGTH,350);

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

bool ZeroLevelIsEmpty(int x, int y, MeshHandle_t handle){
    int ind = x + (y*8);
    if (ind < 64){
        uint64_t mask = (uint64_t)(1ull << ind);
        return ((handle->ZeroLevel & mask) == mask);
    }
    return false;
}

bool N1LevelIsEmpty(int x, int y, MeshHandle_t handle){
    int ind = x + (y*4);
    if (ind < 16){
        uint16_t mask = (uint16_t)(1 << ind);
        return ((handle->N1Level & mask) == mask);
    }
    return false;
}

bool N2LevelIsEmpty(int x, int y, MeshHandle_t handle){
    int ind = x + (y*2);
    if (ind < 4){
        uint8_t mask = (uint8_t)(1 << ind);
        return ((handle->N2Level & mask) == mask);
    }
    return false;
}

// Toy function for tweaking algo (does not modify mesh)
void DeleteGridNodes2D(MeshHandle_t handle){

    double threshold = handle->threshold;
    
    double data_arr[LENGTH];


    for (int i = 0; i < LENGTH; ++i){
        data_arr[i]         = handle->DataField[i].data;

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

