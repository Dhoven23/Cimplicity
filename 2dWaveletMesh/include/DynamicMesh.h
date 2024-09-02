#include <GridData.h>
#include <Indexer.h>

#define DM_printIndex(handle) Indexer_Print(handle->Indexer,handle->length-1)

// Private Data type
struct priv_MeshHandle{
	IndexHandle_t Indexer;
	DataHandle_t  DataField;
	uint64_t ZeroLevel;
	uint16_t N1Level	;
	uint8_t  N2Level : 4;
	uint8_t  N3Level : 1;
	double threshold;
	unsigned length;
};

// Public data type
typedef struct priv_MeshHandle* MeshHandle_t;



// Public functions

void AdaptMesh(MeshHandle_t handle);

bool GenerateMesh(MeshHandle_t handle);

void DestroyMesh(MeshHandle_t handle);

void SetMeshThreshold(double p_threshold, MeshHandle_t handle);

void getMeshThreshold(double* p_threshold, MeshHandle_t handle);

bool GetDataByCoordinate(int x, int y, IndexHandle_t handle, DataHandle_t* data);

bool ZeroLevelCrossIsEmpty(int x, int y, MeshHandle_t handle);

bool N1LevelCrossIsEmpty(int x, int y, MeshHandle_t handle);

bool ZeroLevelSquareIsEmpty(int x, int y, MeshHandle_t handle);

bool N1LevelSquareIsEmpty(int x, int y, MeshHandle_t handle);

bool N2LevelCrossIsEmpty(int x, int y, MeshHandle_t handle);

bool ZeroLevelIsEmpty(int x, int y, MeshHandle_t handle);

bool N1LevelIsEmpty(int x, int y, MeshHandle_t handle);

bool N2LevelIsEmpty(int x, int y, MeshHandle_t handle);

bool N3LevelIsEmpty(int x, int y, MeshHandle_t handle);

