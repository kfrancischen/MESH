#include "Mesh"
#include <complex>
#include <armadillo>
#include <cstring>
#include <cstdlib>
#include <numalloc.h>

void* meshMalloc(size_t size){
  void* ret = malloc_aligned(size, 16);  //allocated memory aligned, c++11
  return ret;
}

void meshFree(void *ptr){
  free_aligned(ptr);
}

void layerInit(Layer *L,
  const *name,
  double thickness,
  const char *material,
  const char *copy){
    L->name = strdup(name);
    L->thickness = thickness;
    L->material = strdup(material);
    L->copy = strdup(copy);
  }

void layerDestroy(Layer *L){
  free(L->name);
  L->name = NULL;
  free(L->material);
  L->material = NULL;
  free(L->copy);
  L->copy = NULL

  return;
}

void materialInit(Material *M,
  const char *name,
  const double eps[2]){
    M->name = strdup(name);
    M->eps.s[0] = eps[0];
    M->eps.s[1] = eps[1];

    return;
}

void materialDestroy(Material *M){
  free(M->name);
  M->name = NULL;

  return;
}


void simulationInit(Simulation *S){
  S->n_G = 0;
  S->solution = NULL;
  S->material = NULL;
  S->layer = NULL;

  return;
}

void simulationDestroy(Simulation *S){
  if(S->solution != NULL){
    simulationDestroySolution(s);
  }

  while(S->layer != NULL){
    Layer *l = S->layer;
    S->layer = S->layer->next;
    layerDestroy(l);
    meshFree(l);
  }

  while(S->material != NULL){
    Material *m = S->material;
    S->material = S->material->next;
    materialDestroy(m);
    meshFree(m);
  }

  return;
}

void simulationDestroySolution(Simulation *S){
  Solution *sol = S->solution;
  if(sol == NULL){
    return;
  }
  if(sol->n_G != NULL){
    sol->n_G = NULL;
  }
  if(sol->omega != NULL){
    sol->omega = NULL;
  }
  if(sol->transmissionFactor != NULL){
    sol->transmissionFactor = NULL;
  }
  if(sol->outputName != NULL){
    sol->outputName = NULL
  }
  return;
}

Material* simulationAddMaterial(Simulation *S){

}

Layer* simulationAddLayer(Simulation *S){

}

int simulationSetNumG(Simulation *S, int n){

}

Layer* simulationGetLayerByName(const Simulation *S,
  const char *name,
  int *index){

}

Layer* simulationGetLayerByIndex(const Simulation *S, int i){

}

Material* simulationGetMaterialByName(const Simulation *S,
  const char *name,
  int *index){

}

Material* simulationGetMaterialByIndex(const Simulation *S, int i){

}

int simulationAddLayeredPatternCircle(Simulation *S,
  Layer *layer,
  int material,
  const double center[2],
  double radius){

}

int simulationAddLayeredPatternRectangle(Simulation *S,
  Layer *layer,
  int material,
  const double center[2],
  double halfwidths[2]){

}

int simulationRemoveLayerPatterns(Simulation *S, Layer *layer){

}

int simulationChangeLayerThickness(Simulation *S, Layer *layer){

}

int simulationInitSolution(Simulation *S){

}

double simulationGetPoyntingFlux(Simulation *S,
  int *sourceLayers,
  int targetLayer){

}
