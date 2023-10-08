#include "header.h"

#include "DynamicMesh.h"


int main()
{
    MeshHandle_t Mesh;

    clock_t begin = clock();

    if (GenerateMesh(Mesh)){
    	fprintf(stderr,"SUCCESS, Mesh generation Complete...\n");


    	// Put your program Here ...
    	
    	SetMeshThreshold(0.1f,Mesh);

    	DeleteGridNodes2D(Mesh);

    	// End your program Here ...
    	DestroyMesh(Mesh);

    } else {

    	fprintf(stderr, "ERROR, mesh not generated\n");
    }

    clock_t end = clock();

    double duration = ((double)end - (double)begin) / CLOCKS_PER_SEC;
    duration *= 1000;
    printf("Time Elapsed: %.3f ms\n",duration);
    return 0;
}
