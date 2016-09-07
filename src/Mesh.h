/* Copyright (C) 2016-2018, Stanford University
 * This file is part of MESH
 * Written by Kaifeng Chen (kfchen@stanford.edu)
 *
 * MESH is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * MESH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef _MESH_H
#define _MESH_H

#ifdef __cplusplus
#include <complex>
extern "C" {
#endif
#include "pattern/pattern.h"
#ifdef __cplusplus
}
#endif


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

/*============================================================
* Helper functions for mallocating and freeing memory
==============================================================*/
void* meshMalloc(size_t size);
void* meshFree(void *ptr);

/*============================================================
* material type
==============================================================*/
enum Type {SCALAR ,TENSOR};

/*============================================================
* material structure
==============================================================*/
typedef struct Material_{
  char *name;
  struct Material_ *next; // linked-list next pointer
  Type type;
  union{
    double s[2];
    double abcde[10];
  } eps;
} Material;

/*============================================================
* layer structure
==============================================================*/
typedef struct Layer_{
	char *name;       // name of layer
	double thickness; // thickness of layer
	char *material;   // name of background material
	Pattern pattern;  // See pattern.h
	char *copy;       // See below.
	struct Layer_ *next; // linked-list next pointer
  bool isSource;
  bool isTarget;
} Layer;
// If a layer is a copy, then `copy' is the name of the layer that should
// be copied, and `material' and `pattern' are inherited, and so they can
// be arbitrary. For non-copy layers, copy should be NULL.

/*============================================================
* solution structure
==============================================================*/
typedef struct Solution_{
  int **G;
  double *omega;
  double *transmissionFactor;
  char *outputName;
  double *kx, *ky;
} Solution;

/*============================================================
* simulation structure
==============================================================*/
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

/*============================================================
* functions to initiate, destroy, and modify a layer
==============================================================*/
void layerInit(Layer *L, const char *name, double thickness, const char *material, const char* copy);
void layerDestroy(Layer *L);
void setSource(Layer *L);
void setTarget(Layer *L);

/*============================================================
* functions to initiate, and destroy a material
==============================================================*/
void materialInit(Material *M, const char *name, const double eps[2]);
void materialInitTensor(Material *M, const char *name, const double abcde[10]);
void materialDestroy(Material *M);

/*============================================================
* functions to initiate, clone and destroy a simulation
==============================================================*/
void simulationInit(Simulation *S);
void simulationDestroy(Simulation *S);
void simulationClone(Simulation *S, Simulation *T);
void simulationDestroySolution(Simulation *S);

/*============================================================
* functions to add material, or layer, or setting n_G for the simulation
==============================================================*/
Material* simulationAddMaterial(Simulation *S);
Layer* simulationAddLayer(Simulation *S);
int simulationSetNumG(Simulation *S, int n);

/*============================================================
* functions to get n_G, layer and materials from a simulation
==============================================================*/
int simulationGetNumG(Simulation *S, int **G);
Layer* simulationGetLayerByName(const Simulation *S, const char *name, int *index);
Layer* simulationGetLayerByIndex(const Simulation *S, int i);
Material* simulationGetMaterialByName(const Simulation *S, const char *name, int *index);
Material* simulationGetMaterialByIndex(const Simulation *S, int i);

/*============================================================
* functions to add patterns to the layer in the simulation
==============================================================*/
int simulationAddLayeredPatternCircle(Simulation *S, Layer *layer, int material, const double center[2], double radius);
int simulationAddLayeredPatternRectangle(Simulation *S, Layer *layer, int material, const double center[2], double halfwidths[2]);

/*============================================================
* functions to change layers in the simulation
==============================================================*/
int simulationRemoveLayerPatterns(Simulation *S, Layer *layer);
int simulationChangeLayerThickness(Simulation *S, Layer *layer, const double *thickness);

/*============================================================
* functions to solve for the heat transfer
==============================================================*/
int simulationInitSolution(Simulation *S);
double simulationGetPoyntingFlux(Simulation *S);

#ifdef __cplusplus
}
#endif

#endif
