#ifndef _MESH_H
#define _MESH_H

#include <complex>
#include <armadillo>

# ifdef __cplusplus
#  include <cstdio>
# else
#  include <stdio.h>
# endif

#ifdef __cplusplus
# include <cstdlib>
#else
# include <stdlib.h>
#endif

void* meshMalloc(size_t size);
void* meshFree(void *ptr);

typedef struct Material_{
  char *name;
  struct Material_ *next; // linked-list next pointer
  int type;
  union{
    double s[2];
    double abcde[10];
  } eps;
} Material;

typedef struct Layer_{
	char *name;       // name of layer
	double thickness; // thickness of layer
	char *material;   // name of background material
	//Pattern pattern;  // See pattern.h
	char *copy;       // See below.
	struct Layer_ *next; // linked-list next pointer
} Layer;


typedef struct Solution_{
  int n_G;
  double *omega;
  double *transmissionFactor;
  char *outputName;
} Solution;

typedef struct Simulation_{
  int n_G;
  Material *material;
  Layer *layer;
  Solution *solution;
  //Options options;
} Simulation;

#ifdef __cplusplus
extern "C"{
#endif

void layerInit(Layer *L, const char *name, double thickness, const char *material, const char* copy);
void layerDestroy(Layer *L);

void materialInit(Material *M, const char *name, const double eps[2]);
void materialInitTensor(Material *M, const char *name, const double abcde[10]);
void materialDestroy(Material *M);

void simulationInit(Simulation *S);
void simulationDestroy(Simulation *S);
void simulationClone(Simulation *S, Simulation *T);

void simulationDestroySolution(Simulation *S);

Material* simulationAddMaterial(Simulation *S);
Layer* simulationAddLayer(Simulation *S);

int simulationSetNumG(Simulation *S, int n);

Layer* simulationGetLayerByName(const Simulation *S, const char *name, int *index);
Layer* simulationGetLayerByIndex(const Simulation *S, int i);

Material* simulationGetMaterialByName(const Simulation *S, const char *name, int *index);
Material* simulationGetMaterialByIndex(const Simulation *S, int i);

int simulationAddLayeredPatternCircle(Simulation *S, Layer *layer, int material, const double center[2], double radius);
int simulationAddLayeredPatternRectangle(Simulation *S, Layer *layer, int material, const double center[2], double halfwidths[2]);

int simulationRemoveLayerPatterns(Simulation *S, Layer *layer);
int simulationChangeLayerThickness(Simulation *S, Layer *layer, const double *thickness);

int simulationInitSolution(Simulation *S);
double simulationGetPoyntingFlux(Simulation *S, int *sourceLayers, int targetLayer);

#ifdef __cplusplus
}
#endif

#endif
