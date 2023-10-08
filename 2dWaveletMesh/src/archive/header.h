#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define ANS 0


struct priv_GridNode {

    int coordinate[2];
    uint32_t neighbors;
    void* DataPointer;
};

struct priv_Data{
    double data;
};

typedef struct priv_GridNode* GridNodeHandle_t;

typedef struct priv_Data* DataHandle_t;

bool InitData(DataHandle_t *handle, unsigned length);

bool InitGrid(GridNodeHandle_t *handle, unsigned length);

void FreeData(DataHandle_t handle);

bool PopulateData(DataHandle_t handle, unsigned length);

void PrintData(DataHandle_t handle, unsigned length);

void SmoothDataInMemory(DataHandle_t handle, unsigned length, unsigned iterations);

void SmoothDataLocal(DataHandle_t handle, unsigned length, unsigned iterations);

void SmoothData2D(DataHandle_t handle, unsigned length, unsigned iterations);

void DerivativeEuler(DataHandle_t* out_handle, DataHandle_t handle, unsigned length);

void DerivativeTrapezoidal(DataHandle_t* out_handle, DataHandle_t handle, unsigned length);

void DeleteGridNodes2D(DataHandle_t data_handle, GridNodeHandle_t grid_handle, unsigned length, double threshold);

int GetIndex(int X, int Y);


