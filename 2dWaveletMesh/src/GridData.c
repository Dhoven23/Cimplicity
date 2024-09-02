#include <header.h>
#include "GridData.h"

// ---------------------------------------------------
/*  This Class handles the direct manipulation of data 
    in the 2D Grid. Use along with Indexer class to work 
    with dynamic data arrays. 
*/


bool InitData(DataHandle_t *handle, unsigned length){
    DataHandle_t p_handle;

    p_handle = malloc(length * (sizeof *p_handle));
    
    if (p_handle != NULL){
        printf("DATA: %lu bytes allocated from %p to %p\n",length * (sizeof *p_handle), p_handle, (p_handle+length));
        *handle = p_handle;
        return true;
    } else {
        return false;
    }
}


void Data_Destroy(DataHandle_t handle){

    free(handle);
}

bool PopulateData(DataHandle_t handle, unsigned length){
    // If memory is allocated, populate with random values
    if (handle != NULL){
        for (int i = 0; i < length; ++i){
            handle[i].data = .15*((double)rand() / (double)RAND_MAX);
            handle[i].boundary[0] = 0.0f;
            handle[i].boundary[1] = 0.0f;
            handle[i].boundary[2] = 0.0f;
            handle[i].boundary[3] = 0.0f;
        }
        handle[0].data = 0;
        return true;

    } else {

        return false;
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

int GetIndex(int X, int Y){
    return (X + Y*17);
}

bool GetData(void* data,DataHandle_t* p_data){
    if ((data != NULL) && (p_data != NULL)){
        *p_data = (DataHandle_t)(data);
        return true;
    }
    return false;
}

bool getBounds(float* bounds, DataHandle_t p_data){
    if((bounds != NULL) && (p_data != NULL)){
        for(int i = 0; i < 4; i++){
            *(bounds+i) = p_data->boundary[i];
        }
        return true;
    }
    return false;
}

bool setBounds(DataHandle_t handle, float b0, float b1, float b2, float b3){
    if (handle != NULL){
        handle->boundary[0] = b0;
        handle->boundary[1] = b1;
        handle->boundary[2] = b2;
        handle->boundary[3] = b3;



        return true;
    }
    return false;
}

bool GetScalar(DataHandle_t data, double* value){
    if(data != NULL){
        *value = data->data;
        return true;
    } else {
        return false;
    }
    return false;
}





