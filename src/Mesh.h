#ifndef _MESH_H
#define _MESH_H

#include <complex>
#include <armadillo>
namespace MESH{

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
  double *omega;
  double *transmissionFactor;
  char *outputName;
} Solution;

typedef struct Simulation_{
  int G_x, G_y;
  Material *material;
  Layer *layer;
  Solution *solution;
  //Options options;
} Simulation;

void layerInit(Layer *L, const char *name, double thickness, const char *material, const char* copy);
void layerDestroy(Layer *L);
void materialInit(Material *M, const char *name, const double eps[2]);
void materialDestry(Material *M);

void simulationInit(Simulation *S);
void simulationDestroy(Simulation *S);



};

#endif
