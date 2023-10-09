#include <GridData.h>
#include <Indexer.h>

#define DM_printIndex(handle,length) Indexer_Print(handle,length)

// Private Data type
struct priv_MeshHandle{
	IndexHandle_t Indexer;
	DataHandle_t  DataField;
	uint64_t ZeroLevel;
	uint16_t N1Level	;
	uint8_t  N2Level : 4;
	uint8_t  N3Level : 1;
	double threshold;
};

// Public data type
typedef struct priv_MeshHandle* MeshHandle_t;

// private functions 

void MarchingCrossTypeZero(MeshHandle_t handle);

static void MarchingCrossTypeN(MeshHandle_t handle, int N);

static void CompressMesh(MeshHandle_t handle);

static void DecompressMesh(MeshHandle_t handle);

static void RefineMesh(MeshHandle_t handle);

static void TransposeMesh(MeshHandle_t handle);

// toy functuons

void DeleteGridNodes2D(MeshHandle_t handle);

// Public functions

void AdaptMesh(MeshHandle_t handle);

bool GenerateMesh(MeshHandle_t handle);

void DestroyMesh(MeshHandle_t handle);

void SetMeshThreshold(double p_threshold, MeshHandle_t handle);

void getMeshThreshold(double* p_threshold, MeshHandle_t handle);

bool GetDataByCoordinate(int x, int y, IndexHandle_t handle, DataHandle_t* data);

bool ZeroLevelIsEmpty(int x, int y, MeshHandle_t handle);

bool N1LevelIsEmpty(int x, int y, MeshHandle_t handle);

bool N2LevelIsEmpty(int x, int y, MeshHandle_t handle);