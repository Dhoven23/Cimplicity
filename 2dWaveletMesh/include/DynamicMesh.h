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
	// Tiling support: adjacent tile handles [right=0, up=1, left=2, down=3]
	// and this tile's origin in global node-space coordinates.
	struct priv_MeshHandle* tile_neighbors[4];
	int tile_origin[2];
};

// Public data type
typedef struct priv_MeshHandle* MeshHandle_t;



// Public functions

void AdaptMesh(MeshHandle_t handle);

bool GenerateMesh(MeshHandle_t handle);

void DestroyMesh(MeshHandle_t handle);

void GetGridLines(MeshHandle_t handle, int* numVertices, int* numIndices, float* vertices, int* indices);

void getGridTriangles(MeshHandle_t handle, int* numVertices, int* numIndices, float* vertices, int* indices);

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

// Tiling: register an adjacent tile in direction dir (0=right,1=up,2=left,3=down).
// Wires the shared border IndexNode neighbor pointers bidirectionally so that nodes
// at tile seams see their cross-tile neighbors after AdaptMesh().
// Call after both tiles have been adapted, and again whenever either is re-adapted.
bool Mesh_ConnectTile(MeshHandle_t handle, int direction, MeshHandle_t neighbor);

