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

// private functions
static void MarchingCrossTypeZero(MeshHandle_t handle, char* buff);

static void MarchingCrossTypeN(MeshHandle_t handle, int N, char* buff);

static void DiagonalMarch(MeshHandle_t handle, int Nlevel, char* buff);

static void CompressMesh(MeshHandle_t handle);

static void DecompressMesh(MeshHandle_t handle);

static void RefineMesh(MeshHandle_t handle, char* buff);

static void TransposeMesh(MeshHandle_t handle, char* buff);

static void ResetNeighbors(MeshHandle_t handle);

// Grid Logic
static void KeepZoneN2Level(int x, int y, MeshHandle_t handle);

static void DeleteZoneN2Level(int x, int y, MeshHandle_t handle);


static void KeepZoneZeroLevel(int x, int y, MeshHandle_t handle);

static void DeleteZoneZeroLevel(int x, int y, MeshHandle_t handle);

static void KeepZoneN1Level(int x, int y, MeshHandle_t handle);

static void DeleteZoneN1Level(int x, int y, MeshHandle_t handle);

