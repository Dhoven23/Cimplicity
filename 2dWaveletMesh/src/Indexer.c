#include <header.h>
#include "Indexer.h"

// Cache for coordinates-by-index
static int Coords[LENGTH];

// Constructor
bool Indexer_Create(IndexHandle_t* handle, unsigned length){
	if (handle == NULL)
	{
		return false;
	}

	IndexHandle_t p_handle;
	p_handle = malloc(length * sizeof(*p_handle));

	for (int i = 0; i < LENGTH; ++i){
		Coords[i] = -1;
	}

	for (int i = 0; i < length; ++i)
	{
		int x = i % 17;
		int y = floor(i/17);
		
		(p_handle+i)->coordinate[0] = x;
		(p_handle+i)->coordinate[1] = y;
		(p_handle+i)->neighbors[0]  = (x < 16) ? (void*)(p_handle+i+1) : NULL;
		(p_handle+i)->neighbors[1]	= (y > 0) ? (void*)(p_handle+i-17) : NULL;
		(p_handle+i)->neighbors[2]	= (x > 0) ? (void*)(p_handle+i- 1) : NULL;
		(p_handle+i)->neighbors[3]	= (y < 16) ? (void*)(p_handle+i+17) : NULL;
		(p_handle+i)->data_ptr		= NULL;
		(p_handle+i)->b_IsInterp 	= false;	
		Coords[i] = x + y*17;
	}

	*handle = p_handle;
	printf("INDEX: %lu bytes allocated from %p to %p\n",length * (sizeof *p_handle), p_handle, (p_handle+length));
	return true;

}

void Indexer_SetCoordinateCache(int CacheIndex,int varIndex){
	Coords[CacheIndex] = varIndex;
}

void Indexer_GetCoordinateCache(int CacheIndex,int* p_varIndex){
	*p_varIndex = Coords[CacheIndex];
}

bool Indexer_SetInterpTrue(IndexHandle_t handle){
	if (handle != NULL){
		handle->b_IsInterp = true;
		return true;
	}
	return false;
}

bool Indexer_SetInterpFalse(IndexHandle_t handle){
	if (handle != NULL){
		handle->b_IsInterp = false;
		return true;
	}
	return false;
}

// Lookup neighbor
bool Indexer_GetNeighbor(IndexHandle_t* ret_handle, IndexHandle_t handle, int loc){
	if ((loc < 4)){
		IndexHandle_t p_handle = NULL;
		p_handle = (handle->neighbors[loc] != NULL) ? (IndexHandle_t)(handle->neighbors[loc]) : NULL;
	
		if (p_handle != NULL){
			*ret_handle = p_handle;
			printf("Neightbor Coords: (%i,%i)\n",p_handle->coordinate[0],p_handle->coordinate[1]);
			return true;
		} else {
			return false;
		}
	} else {
		return false;
	}
}

bool Indexer_SetNeighbor(int loc, int x, int y, IndexHandle_t handle, IndexHandle_t index){
	if ((loc < 4) && (loc > -1)){
		IndexHandle_t p_handle;
		if (Indexer_GetIndexByCoordinate(1, x, y, index, &p_handle)){
			handle->neighbors[loc] = p_handle;
			return true;
		}
	}
	return false;
}

// get XY coordinates of grid node (or at offset of index if index != 0)
bool Indexer_GetCoordinates(int* x, int* y, IndexHandle_t handle, int index){
	if ((handle+index)->data_ptr != NULL)
	{
		*x = (handle+index)->coordinate[0];
		*y = (handle+index)->coordinate[1];
		return true;
	}
	return false;
}

// free pointer
bool Indexer_Destroy(IndexHandle_t handle){
	if(handle != NULL){
		free(handle);
		return true;
	} else {
		return false;
	}
}

// Diagnostic print
void Indexer_Print(IndexHandle_t handle, unsigned length){
	for (int i = 0; i < length; ++i)
	{
		printf("%i: Index[%i][%i] -> %i \n",i, 
			(handle+i)->coordinate[0], 
			(handle+i)->coordinate[1],
			(int)(handle+i)->b_IsInterp
			);

	}
}

// get grid node by coordinates (if it exists)
bool Indexer_GetIndexByCoordinate(bool b_useCache,int x, int y, IndexHandle_t handle, IndexHandle_t* found_handle){

	bool b_found = false;
	int max_depth = 17;
	int tries = 0;
	int direction = 0, init_direction = 0;
	IndexHandle_t p_handle;

	if ((x > 16) || (y > 16) || (x < 0) || (y < 0)) {return false;}

	b_found = priv_CoordinateSearch(b_useCache, x,y,handle,&p_handle);
	if (b_found){
		*found_handle = p_handle;
		return true;
	} else {
		return false;
	}
}


static bool priv_CoordinateSearch(bool b_useCache, int x, int y, IndexHandle_t p_handle, IndexHandle_t* p_found_handle){

	IndexHandle_t temp_handle;

	bool b_found = false;
	int count = 0;

	int offset = 0;
	if ((x == 0) && (y == 0)){
		*p_found_handle = p_handle;
		return true;
	}

	// Search Cache
	int ind = x + 17*y;
	for (int i = 0; i < LENGTH; ++i){
		int p_ind = Coords[i];
		
		if (p_ind == ind){
			*p_found_handle = (p_handle+i);
			return true;
		} else if ((p_ind == -1) && b_useCache){
			return false;
		}
	}

	// if we missed the cache, brute force... 
	while(true){

		temp_handle = &(p_handle[count++]);
		int priv_x, priv_y;

		Indexer_GetCoordinates(&priv_x,&priv_y,temp_handle,0);

		if ((priv_y == y) && (priv_x == x)){
			if (temp_handle->data_ptr != NULL) {
				b_found = true;
			} else {
				b_found = false;
			}
		}

		if (b_found) {

			*p_found_handle = temp_handle;
			return true;
		} else if (count > 17*17) {
			return false;
		}
	
	}
	return false;
}
