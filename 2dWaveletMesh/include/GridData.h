#include <header.h>

#define GRIDDATA 1

struct priv_Data{
    double data;
};

typedef struct priv_Data* DataHandle_t;

bool InitData(DataHandle_t *handle, unsigned length);

void Data_Destroy(DataHandle_t handle);

bool PopulateData(DataHandle_t handle, unsigned length);

void PrintData(DataHandle_t handle, unsigned length);

void SmoothData2D(DataHandle_t handle, unsigned length, unsigned iterations);

int GetIndex(int X, int Y);

bool GetData(void* data,DataHandle_t* p_data);

bool GetScalar(DataHandle_t data, double* value);