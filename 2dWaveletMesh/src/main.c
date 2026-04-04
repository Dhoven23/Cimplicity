#include "header.h"

#include "DynamicMesh.h"


int main()
{
    // MeshHandle_t is a pointer type; allocate the backing struct first.
    MeshHandle_t Mesh = malloc(sizeof(*Mesh));
    if (Mesh == NULL){
        fprintf(stderr, "ERROR, malloc failed\n");
        return 1;
    }

    clock_t begin = clock();

    if (GenerateMesh(Mesh)){
    	fprintf(stderr,"SUCCESS, Mesh generation Complete...\n");


    	// Put your program Here ...

    	SetMeshThreshold(0.1f,Mesh);

    	AdaptMesh(Mesh);

        //DM_printIndex(Mesh);

        // for (int i = 0; i < 8; ++i){
        //     for (int j = 0; j < 8; ++j){
        //         printf(" %i ",(int)ZeroLevelIsEmpty(j,i,Mesh));
        //     } printf("\n");
        // }
        // printf("\n");
        // for (int i = 0; i < 4; ++i){
        //     for (int j = 0; j < 4; ++j){
        //         printf(" %i ",(int)!N1LevelIsEmpty(j,i,Mesh));
        //     } printf("\n");
        // }
    	// End your program Here ...
    	DestroyMesh(Mesh);

    } else {

    	fprintf(stderr, "ERROR, mesh not generated\n");
    }

    free(Mesh);

    clock_t end = clock();

    double duration = ((double)end - (double)begin) / CLOCKS_PER_SEC;
    duration *= 1000;
    printf("Time Elapsed: %.3f ms\n",duration);
    return 0;
}
