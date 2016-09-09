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

#include "Mesh.h"
#include <complex>
#include <armadillo>
#include <cstring>
#include <cstdlib>
#include <numalloc.h>


/*============================================================
* Helper functions for mallocating and freeing memory
==============================================================*/
void* meshMalloc(size_t size){
  void* ret = malloc_aligned(size, 16);  //allocated memory aligned, c++11
  return ret;
}

void meshFree(void *ptr){
  free_aligned(ptr);
  return;
}

/*============================================================
* functions to initiate, destroy, and modify a layer
==============================================================*/
void layerInit(Layer *L,
  const *name,
  double thickness,
  const char *material,
  const char *copy){
    L->name = strdup(name);
    L->thickness = thickness;
    L->pattern.nshapes = 0
    L->pattern.shapes = NULL;
    L->pattern.parent = NULL;
    L->material = strdup(material);
    L->copy = strdup(copy);
    L->isSource = false;
    L->isTarget = false;
    return;
  }

void layerDestroy(Layer *L){
  free(L->name);
  L->name = NULL;
  free(L->material);
  L->material = NULL;
  free(L->copy);
  L->copy = NULL
  free(L->pattern.shapes);
  L->pattern.shapes = NULL;
  free(l->pattern.parent);
  L->pattern.parent = NULL;

  return;
}

void setSource(Layer *L){
  L->isSource = true;
  return;
}

void setTarget(Layer *L){
  L->isTarget = true;
  return;
}

/*============================================================
* functions to initiate, and destroy a material
==============================================================*/

void materialInit(Material *M,
  const char *name,
  const double eps[2]){
    M->name = strdup(name);
    M->eps.s[0] = eps[0];
    M->eps.s[1] = eps[1];
    M->type = SCALAR;
    return;
}

void materialInitTensor(Material *M,
  const char *name,
  const double abcde[10]){
    // TODO
    return;
}

void materialDestroy(Material *M){
  free(M->name);
  M->name = NULL;

  return;
}

/*============================================================
* functions to initiate, clone and destroy a simulation
==============================================================*/
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

void simulationClone(Simulation *S, Simulation *T){
  if(S == NULL || T == NULL){
    return;
  }
  std::memcpy(T, S, sizeof(Simulation));
  Layer *l = S->layer;
  Layer **l2 = &T->layer;

  while(l != NULL){
		*l2 = (Layer*)meshMalloc(sizeof(Layer));
		layerInit(*l2, l->name, l->thickness, l->material, l->copy);

		// Copy pattern
		(*l2)->pattern.nshapes = l->pattern.nshapes;
		(*l2)->pattern.shapes = (shape*)malloc(sizeof(shape)*l->pattern.nshapes);
		std::memcpy((*l2)->pattern.shapes, l->pattern.shapes, sizeof(shape)*l->pattern.nshapes);
		(*l2)->pattern.parent = NULL;

		(*l2)->next = NULL;
		l2 = &((*l2)->next);
		l = l->next;
	}

	Material *m = S->material;
	Material **m2 = &T->material;
	while(m != NULL){
		*m2 = (Material*)meshMalloc(sizeof(Material));
		if(m->type == SCALAR){
			materialInit(*m2, m->name, m->eps.s);
		}else{
			materialInitTensor(*m2, m->name, m->eps.abcde);
		}

		(*m2)->next = NULL;
		m2 = &((*m2)->next);
		m = m->next;
	}
  return;
}

void simulationDestroySolution(Simulation *S){
  Solution *sol = S->solution;
  if(sol == NULL){
    return;
  }
  if(sol->G != NULL){
    meshFree(sol->G);
    sol->G = NULL;
  }
  if(sol->omega != NULL){
    meshFree(sol->omega);
    sol->omega = NULL;
  }
  if(sol->transmissionFactor != NULL){
    meshFree(sol->transmissionFactor);
    sol->transmissionFactor = NULL;
  }
  if(sol->outputName != NULL){
    meshFree(sol->outputName);
    sol->outputName = NULL
  }
  if(sol->kx != NULL){
    meshFree(sol->kx);
    sol->kx = NULL;
  }
  if(sol->ky != NULL){
    meshFree(sol->ky);
    sol->ky = NULL;
  }

  return;
}

/*============================================================
* functions to add material, or layer, or setting n_G for the simulation
==============================================================*/

Material* simulationAddMaterial(Simulation *S){
  Material *m;
  if(S->material == NULL){
    S->material = (Material *)meshMalloc(sizeof(Material));
    S->material->next = NULL;
    m = S->material;
  }
  else{
    m = S->material;
    while(m->next != NULL){
      m = m->next;
    }
    m->next = (Material *)meshMalloc(sizeof(Material));
    m = m->next;
    m->next = NULL;
  }
  materialInit(m, NULL, NULL);
  return m;
}

Layer* simulationAddLayer(Simulation *S){
  Layer *l;
  if(S->layer == NULL){
    S->layer = (Layer *)meshMalloc(sizeof(Layer));
    S->layer->next = NULL;
    l = S->layer;
  }
  else{
    l = S->layer;
    while(l->next != NULL){
      l = l->next;
    }
    l->next = (Layer *)meshMalloc(sizeof(Layer));
    l = l->next;
    l->next = NULL;
  }
  layerInit(l, NULL, 0, NULL, NULL);
  return l;
}

int simulationSetNumG(Simulation *S, int n){
  if(S == NULL){
    return -1;
  }
  simulationDestroySolution(S);
  S->n_G = n;
  return 0;
}

/*============================================================
* functions to get n_G, layer and materials from a simulation
==============================================================*/

int simulationGetNumG(Simulation *S, int **G){
  int result = 0;
  if(S == NULL){
    return -1;
  }
  result = S->n_G;
  if(G != NULL){
    if(S->solution != NULL){
      *G = S->solution->G;
    }
  }
  return result;
}

Layer* simulationGetLayerByName(const Simulation *S,
  const char *name,
  int *index){
  Layer* l = S->layer;
  int i = 0;
  while(l != NULL){
    if(strcmp(l->name, name) == 0){
      if(index != NULL){
        *index = i;
      }
      return l;
    }
    l = l->next;
    ++i;
  }
  return NULL;
}

Layer* simulationGetLayerByIndex(const Simulation *S, int i){
  Layer *l = S->Layer;
  while(l != NULL){
    if(i == 0){
      return l;
    }
    --i;
    l = l->next;
  }
  return l;
}

Material* simulationGetMaterialByName(const Simulation *S,
  const char *name,
  int *index){
  Material *m = S->material;
  int i = 0;
  while(M != NULL){
    if(strcmp(M->name, name) == 0){
      if(index != NULL){
        *index = i;
      }
      return M;
    }
    M = M->next;
    ++i;
  }
  return NULL;
}

Material* simulationGetMaterialByIndex(const Simulation *S, int i){
  Material *m = S->material;
  while(m != NULL){
    if(i == 0){
      return m;
    }
    --i;
    m = m->next;
  }
  return NULL;
}

/*============================================================
* functions to add patterns to the layer in the simulation
==============================================================*/

int simulationAddLayeredPatternCircle(Simulation *S,
  Layer *layer,
  int material,
  const double center[2],
  double radius){
    int ret = 0;
    if(S == NULL){
      ret = -1;
    }
    if(layer == NULL){
      ret = -2;
    }
    if(material < 0){
      ret = -3;
    }
    if(center == NULL){
      ret = -4;
    }
    if(radius < 0){
      ret = -5;
    }
    if(ret < 0){
      return ret;
    }
    simulationDestroySolution(S);
    int n = layer->patter.nshapes++;
    layer->pattern.nshapes = (shape*)realloc(layer->pattern.shapes, sizeof(shape)*layer->pattern.nshapes);
    if(layer->pattern.shapes == NULL){
      return 1;
    }
    shape *sh = &layer->pattern.nshapes[n];
    sh->type = CIRCLE;
    sh->center[0] = center[0];
    sh->center[1] = center[1];
    sh->vtab.circle.radius = radius;
    sh->tag = material;
    return 0;
}

int simulationAddLayeredPatternRectangle(Simulation *S,
  Layer *layer,
  int material,
  const double center[2],
  double angle,
  double halfwidths[2]){
    int ret = 0;
    if(S == NULL){
      ret = -1;
    }
    if(layer == NULL){
      ret = -2;
    }
    if(material < 0){
      ret = -3;
    }
    if(center == NULL){
      ret = -4;
    }
    if(halfwidths == NULL){
      ret = -6;
    }
    if(ret != 0){
      return;
    }
    simulationDestroySolution(S);
    int n = layer->patter.nshapes++;
    layer->pattern.nshapes = (shape*)realloc(layer->pattern.shapes, sizeof(shape)*layer->pattern.nshapes);
    if(layer->pattern.shapes == NULL){
      return 1;
    }
    shape *sh = &layer->pattern.nshapes[n];
    sh->type = RECTANGLE;
    sh->center[0] = center[0];
    sh->center[1] = center[1];
    sh->angle = angle;
    sh->vtab.rectange.halfwidths[0] = halfwidths[0];
    sh->vtab.rectange.halfwidths[1] = halfwidths[1];
    sh->tag = material;

    return 0;
}

/*============================================================
* functions to change layers in the simulation
==============================================================*/

int simulationRemoveLayerPatterns(Simulation *S, Layer *layer){

}

int simulationChangeLayerThickness(Simulation *S, Layer *layer){

}

/*============================================================
* functions to solve for the heat transfer
==============================================================*/

int simulationInitSolution(Simulation *S){

}

double simulationGetPoyntingFlux(Simulation *S,
  int *sourceLayers,
  int targetLayer){

}
