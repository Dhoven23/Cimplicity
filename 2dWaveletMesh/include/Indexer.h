#define INDEXER 1

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <time.h>



#ifndef GRIDDATA
	#include "GridData.h"
#endif


struct priv_IndexNode {
	int coordinate[2];
	void* data_ptr;
	void* neighbors[4];
	bool b_IsInterp;
};


typedef struct priv_IndexNode* IndexHandle_t;

bool Indexer_Create(IndexHandle_t* handle, unsigned length);

bool Indexer_Destroy(IndexHandle_t handle);

void Indexer_Print(IndexHandle_t handle, unsigned length);

bool Indexer_GetNeighbor(IndexHandle_t* ret_handle, IndexHandle_t handle, int loc);

bool Indexer_SetNeighbor(int loc, int x, int y, IndexHandle_t handle, IndexHandle_t index);

bool Indexer_GetCoordinates(int* x, int* y, IndexHandle_t handle, int index);

bool priv_CoordinateSearch(int x, int y, IndexHandle_t p_handle, IndexHandle_t* p_found_handle);

bool Indexer_GetIndexByCoordinate(int x, int y, IndexHandle_t handle, IndexHandle_t* found_handle);

bool Indexer_GetDataByCoordinate(int x, int y, IndexHandle_t handle, DataHandle_t* data);

