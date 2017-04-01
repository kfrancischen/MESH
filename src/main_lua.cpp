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
#include "setup.h"
#include <unistd.h>
extern "C"
{
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}
#include "luawrapper/luawrapper.hpp"
#include "luawrapper/luawrapperutil.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#include "luawrapper/lua_mpi.c"
#include "luawrapper/buffer_mpi.c"
#endif

using namespace MESH;

// LuaWrapper knows about primitive types like ints and floats, but it doesn't
// know about things like std::strings or other more complicated types.
// Sometimes, rather than register the type with LuaWrapper, it's easier to
// be able to convert it to and from Lua's primitive types, like strings or
// tables.
//
// To do this, you must write luaU_check, luaU_to and luaU_push functions for
// your type. You don't always need all three, it depends on if you're pushing
// objects to Lua, getting objects from Lua, or both.

// This example uses std::string, but if you have other custom string types it
// should be easy to write versions of those functions too

template<>
struct luaU_Impl<std::string>
{
    static std::string luaU_check(lua_State* L, int index)
    {
        return std::string(luaL_checkstring(L, index));
    }

    static std::string luaU_to(lua_State* L, int index )
    {
        return std::string(lua_tostring(L, index));
    }

    static void luaU_push(lua_State* L, const std::string& val)
    {
        lua_pushstring(L, val.c_str());
    }
};
/*======================================================*/
// wrappaer for the main class
/*=======================================================*/
// this function wraps setPeriodicity(const double p1, const double p2 = 0)
// @how to use:
// SetPeriodicity(p1) or
// Setperiodicity(p1, p2)
int MESH_SetPeriodicity(lua_State* L){
  int n = lua_gettop(L);
	if(n != 2 && n != 3){
		return luaL_error(L, "expecting 1 or 2 arguments");
	}
	Simulation* s = luaW_check<Simulation>(L, 1);
	double p1 = luaU_check<double>(L, 2);
	if(n == 2){
		s->setPeriodicity(p1);
	}
  else{
  	double p2 = luaU_check<double>(L, 3);
  	s->setPeriodicity(p1, p2);
  }
	return 1;
}
// this function wraps addMaterial(const std::string name, const std::string infile)
// @how to use:
// AddMaterial(material name, input file)
int MESH_AddMaterial(lua_State *L){
	Simulation* s = luaW_check<Simulation>(L, 1);
	std::string name = luaU_check<std::string>(L, 2);
	std::string infile = luaU_check<std::string>(L, 3);
	s->addMaterial(name, infile);
	return 1;
}

// this function wraps setMaterial(const std::string name, const double** epsilon, const std::string type)
// @how to use:
// SetMaterial(material name, new epsilon)
int MESH_SetMaterial(lua_State* L){
  Simulation* s = luaW_check<Simulation>(L, 1);
	std::string name = luaU_check<std::string>(L, 2);
  int numOfOmega = lua_rawlen(L, 3);
  bool set = false;
  double** epsilon = new double*[numOfOmega];
  std::string type = "scalar";
  int len;
  for(int i = 0; i < numOfOmega; i++){
    lua_pushinteger(L, i+1);
    lua_gettable(L, 3);
    len = lua_rawlen(L, -1);
    if(!set){
      for(int j = 0; j < numOfOmega; j++){
        epsilon[j] = new double[len];
      }
      if(len == 6) type = "diagonal";
      if(len == 10) type = "tensor";
      set = true;
    }
    for(int j = 0; j < len; j++){
      lua_pushinteger(L, j + 1);
      lua_gettable(L, 4);
      epsilon[i][j] = luaU_check<double>(L, -1);
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }

  s->setMaterial(name, epsilon, type);

  for(int i = 0; i < numOfOmega; i++){
    delete [] epsilon[i];
  }
  delete [] epsilon;
  return 1;
}

// this function wraps addLayer(const std::string name, const double thick, const std::string materialName)
// @how to use:
// AddLayer(layer name, thickness, material name)
int MESH_AddLayer(lua_State *L){
	Simulation* s = luaW_check<Simulation>(L, 1);
	std::string name = luaU_check<std::string>(L, 2);
	double thick = luaU_check<double>(L, 3);
	std::string materialName = luaU_check<std::string>(L, 4);
	s->addLayer(name, thick, materialName);
	return 1;
}

// this function wraps setLayer(const std::string name, const double thick, const std::string materialName)
// @how to use:
// SetLayer(layer name, thickness, material name)
int MESH_SetLayer(lua_State* L){
  Simulation* s = luaW_check<Simulation>(L, 1);
	std::string name = luaU_check<std::string>(L, 2);
	double thick = luaU_check<double>(L, 3);
	std::string materialName = luaU_check<std::string>(L, 4);
	s->setLayer(name, thick, materialName);
	return 1;
}

// this function wraps setLayerThickness(const std::string name, const double thick)
// @how to use
// SetLayerThickness(layer name, thickness)
int MESH_SetLayerThickness(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  double thick = luaU_check<double>(L, 3);
  s->setLayerThickness(name, thick);
  return 1;
}

// this function wraps addLayerCopy(const std::string name, const std::string originalName)
// @how to use
// AddLayerCopy(new layer name, original layer name)
int MESH_AddLayerCopy(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  std::string materialName = luaU_check<std::string>(L, 3);
  s->addLayerCopy(name, materialName);
  return 1;
}

// this function wraps deleteLayer(const std::string name)
// @how to use
// DeleteLayer(layer name)
int MESH_DeleteLayer(lua_State* L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  s->deleteLayer(name);
  return 1;
}

// this function wraps setLayerPatternGrating(const std::string layerName,
//  const std::string materialName, const double center, const double width)
// @how to use
// SetLayerPatternGrating(layer name, material name, center, width)
int MESH_SetLayerPatternGrating(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string layerName = luaU_check<std::string>(L, 2);
  std::string materialName = luaU_check<std::string>(L, 3);
  double center = luaU_check<double>(L, 4);
  double width = luaU_check<double>(L, 5);
  s->setLayerPatternGrating(layerName, materialName, center, width);
  return 1;
}

// this function wraps setLayerPatternRectangle(const std::string layerName, const std::string materialName,
//  const double centerx, const double centery, const double widthx, const double widthy)
// @how to use
// SetLayerPatternRectangle(layer name, material name, {centerx, centery}, {widthx, widthy})
int MESH_SetLayerPatternRectangle(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string layerName = luaU_check<std::string>(L, 2);
  std::string materialName = luaU_check<std::string>(L, 3);
  double vals[4];
  for(int i = 0; i < 2; i++){
    for(int j= 0; j < 2; j++){
      lua_pushinteger(L, j+1);
      lua_gettable(L, 4+i);
      vals[2*i+j] = luaU_check<double>(L, -1);
      lua_pop(L, 1);
    }
  }
  s->setLayerPatternRectangle(layerName, materialName, vals[0], vals[1], vals[2], vals[3]);
  return 1;
}

// this function wraps setLayerPatternCircle(const std::string layerName,const std::string materialName,
//  const double centerx, const double centery, const double radius)
// @how to use
// SetLayerPatternCircle(layer name, material name, {centerx, centery}, radius)
int MESH_SetLayerPatternCircle(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string layerName = luaU_check<std::string>(L, 2);
  std::string materialName = luaU_check<std::string>(L, 3);
  double radius = luaU_check<double>(L, 5);
  double vals[2];
  for(int i = 0; i < 2; i++){
    lua_pushinteger(L, i+1);
    lua_gettable(L, 4);
    vals[i] = luaU_check<double>(L, 1);
    lua_pop(L, 1);
  }
  s->setLayerPatternCircle(layerName, materialName, vals[0], vals[1], radius);
  return 1;
}

// this function wraps setGx(const int nGx)
// @how to use
// SetGx(number of Gx)
int MESH_SetGx(lua_State* L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  int n = luaU_check<int>(L, 2);
  s->setGx(n);
  return 1;
}

// this function wraps setGy(const int nGy)
// @how to use
// SetGy(number of Gy)
int MESH_SetGy(lua_State* L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  int n = luaU_check<int>(L, 2);
  s->setGy(n);
  return 1;
}

// this function wraps setSourceLayer(const std::string name)
// @how to use
// SetSourceLayer(layer name)
int MESH_SetSourceLayer(lua_State* L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  s->setSourceLayer(name);
  return 1;
}

// this function wraps setProbeLayer(const std::string name)
// @how to use
// SetProbeLayer(layer name)
int MESH_SetProbeLayer(lua_State* L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  s->setProbeLayer(name);
  return 1;
}


// this function wraps setThread(const int numThread)
// @how to use
// SetThread(number of thread)
int MESH_SetThread(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  int thread = luaU_check<int>(L, 2);
  s->setThread(thread);
  return 1;
}

// this function wraps buildRCWA()
// @how to use
// BuildRCWA()
int MESH_BuildRCWA(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  s->buildRCWA();
  return 1;
}

// this function wraps getPhi()
// @how to use
// GetPhi()
int MESH_GetPhi(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  double* phi = s->getPhi();
  int numOfOmega = s->getNumOfOmega();
  lua_createtable(L, numOfOmega, 0);
  for(int i = 0; i < numOfOmega; i++){
    lua_pushinteger(L, i+1);
    lua_pushnumber(L, phi[i]);
    lua_settable(L, -3);
  }
  return 1;
}

// this function wraps getOmega()
// @how to use
// GetOmega()
int MESH_GetOmega(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  double* omega = s->getOmega();
  int numOfOmega = s->getNumOfOmega();
  for(int i = 0; i < numOfOmega; i++){
    lua_pushinteger(L, i+1);
    lua_pushnumber(L, omega[i]);
    lua_settable(L, -3);
  }
  return 1;
}

// this function wraps getNumOfOmega()
// @how to use
// GetNumOfOmega()
int MESH_GetNumOfOmega(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  luaU_push(L, s->getNumOfOmega());
  return 1;
}

// this function wraps getPhiAtKxKy(const int omegaIndex, const double kx, const double ky = 0)
// @how to use
// GetPhiAtKxKy(omega index, kx) or
// GetPhiAtKxKy(omega index, kx, ky)
int MESH_GetPhiAtKxKy(lua_State *L){
  int n = lua_gettop(L);
	if(n != 3 && n != 4){
		return luaL_error(L, "expecting 2 or 3 arguments");
	}
  Simulation* s = luaW_check<Simulation>(L, 1);
  int omegaIdx = luaU_check<int>(L, 2);
  double kx = luaU_check<double>(L, 3);
  if(n == 3){
    luaU_push(L, s->getPhiAtKxKy(omegaIdx, kx));
  }
  else{
    double ky = luaU_check<double>(L, 4);
    luaU_push(L, s->getPhiAtKxKy(omegaIdx, kx, ky));
  }
  return 1;
}

// this function wraps getSysInfo()
// @how to use
// GetSysInfo()
int MESH_GetSysInfo(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  s->getSysInfo();
  return 1;
}

// this function wraps optUseInverseRule()
// @how to use
// OptUseInverseRule()
int MESH_OptUseInverseRule(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  s->optUseInverseRule();
  return 1;
}

// this function wraps optUseNaiveRule()
// @how to use
// OptUseNaiveRule
int MESH_OptUseNaiveRule(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  s->optUseNaiveRule();
  return 1;
}

// this function wraps optPrintIntermediate()
// @how to use
// OptPrintIntermediate()
int MESH_OptPrintIntermediate(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  s->optPrintIntermediate();
  return 1;
}

// this function wraps setKxIntegral(const int points, const double end = 0)
// @how to use
// SetKxIntegral(points) or
// SetKxIntegral(points, end)
int MESH_SetKxIntegral(lua_State *L){
  int n = lua_gettop(L);
	if(n != 2 && n != 3){
		return luaL_error(L, "expecting 1 or 2 arguments");
	}
  Simulation* s = luaW_check<Simulation>(L, 1);
  int points = luaU_check<int>(L, 2);
  if(n == 2){
    s->setKxIntegral(points);
  }
  else{
    double end = luaU_check<double>(L, 3);
    s->setKxIntegral(points, end);
  }
  return 1;
}

// this function wraps setKyIntegral(const int points, const double end = 0)
// @how to use
// SetKyIntegral(points) or
// SetKyIntegral(points, end)
int MESH_SetKyIntegral(lua_State *L){
  int n = lua_gettop(L);
  if(n != 2 && n != 3){
    return luaL_error(L, "expecting 1 or 2 arguments");
  }
  Simulation* s = luaW_check<Simulation>(L, 1);
  int points = luaU_check<int>(L, 2);
  if(n == 2){
    s->setKyIntegral(points);
  }
  else{
    double end = luaU_check<double>(L, 3);
    s->setKyIntegral(points, end);
  }
  return 1;
}

// this function wraps setKxIntegralSym(const int points, const double end = 0)
// @how to use
// SetKxIntegralSym(points) or
// SetKxIntegralSym(points, end)
int MESH_SetKxIntegralSym(lua_State *L){
  int n = lua_gettop(L);
  if(n != 2 && n != 3){
    return luaL_error(L, "expecting 1 or 2 arguments");
  }
  Simulation* s = luaW_check<Simulation>(L, 1);
  int points = luaU_check<int>(L, 2);
  if(n == 2){
    s->setKxIntegralSym(points);
  }
  else{
    double end = luaU_check<double>(L, 3);
    s->setKxIntegralSym(points, end);
  }
  return 1;
}

// this function wraps setKyIntegralSym(const int points, const double end = 0)
// @how to use
// SetKyIntegralSym(points) or
// SetKyIntegralSym(points, end)
int MESH_SetKyIntegralSym(lua_State *L){
  int n = lua_gettop(L);
  if(n != 2 && n != 3){
    return luaL_error(L, "expecting 1 or 2 arguments");
  }
  Simulation* s = luaW_check<Simulation>(L, 1);
  int points = luaU_check<int>(L, 2);
  if(n == 2){
    s->setKyIntegralSym(points);
  }
  else{
    double end = luaU_check<double>(L, 3);
    s->setKyIntegralSym(points, end);
  }
  return 1;
}

// this function wraps integrateKxKy()
// @how to use
// IntegrateKxKy()
int MESH_IntegrateKxKy(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  s->integrateKxKy();
  return 1;
}

// this function wraps integrateKxKyMPI(const int rank, const int size)
// @how to use
// IntegrateKxKyMPI(rank, size)
int MESH_IntegrateKxKyMPI(lua_State* L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  int rank = luaU_check<int>(L, 2);
  int size = luaU_check<int>(L, 3);
  s->integrateKxKyMPI(rank, size);
  return 1;
}
// this function wraps outputStructurePOVRay(const std::string outfile)
// @how to use
// OutputStructurePOVRay(file name)
int MESH_OutputStructurePOVRay(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  std::string outfile = luaU_check<std::string>(L, 2);
  s->outputStructurePOVRay(outfile);
  return 1;
}

/*======================================================*/
// constructor for the planar
/*=======================================================*/
SimulationPlanar* MESH_SimulationPlanar_New(lua_State* L){
	return new SimulationPlanar();
}

// this function wraps optUseQuadgk()
// @how to use
// OptUseQuadgk()
int MESH_OptUseQuadgk(lua_State *L){
	SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
	s->optUseQuadgk();
	return 1;
}

// this function wraps optUseQuadgl()
// @how to use
// OptUseQuadgl() or
// OptUseQuadgl(degree)
int MESH_OptUseQuadgl(lua_State *L){
  int n = lua_gettop(L);
  if(n != 1 && n != 2){
		return luaL_error(L, "expecting no or 1 argument");
	}
	SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
  if(n == 1){
  	s->optUseQuadgl();
  }
  else{
    int degree = luaU_check<int>(L, 2);
    s->optUseQuadgl(degree);
  }
	return 1;
}

// this function wraps setKParallelIntegral(const double end)
// @how to use
// SetKParallel(end)
int MESH_SetKParallel(lua_State *L){
  SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
  double end = luaU_check<double>(L, 2);
  s->setKParallelIntegral(end);
  return 1;
}

// this function wraps getPhiAtKParallel(const int omegaIndex, const double KParallel)
// @how to use
// GetPhiAtKParallel(omega index, k parallel value)
int MESH_GetPhiAtKParallel(lua_State *L){
  SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
  int omegaIdx = luaU_check<int>(L, 2);
  double k = luaU_check<double>(L, 3);
  luaU_push(L, s->getPhiAtKParallel(omegaIdx, k));
  return 1;
}

// this function wraps integrateKParallel()
// @how to use
// IntegrateKPrarallel()
int MESH_IntegrateKParallel(lua_State* L){
  SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
  s->integrateKParallel();
  return 1;
}

// this function wraps getPhi()
// @how to use
// GetPhi()
int MESH_GetPhi_Planar(lua_State* L){
  SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
  double* phi = s->getPhi();
  int numOfOmega = s->getNumOfOmega();
  lua_createtable(L, numOfOmega, 0);
  for(int i = 0; i < numOfOmega; i++){
    lua_pushinteger(L, i+1);
    lua_pushnumber(L, phi[i]);
    lua_settable(L, -3);
  }
  return 1;
}

/*======================================================*/
// constructor for the 1D grating
/*=======================================================*/
SimulationGrating* MESH_SimulationGrating_New(lua_State* L){
  return new SimulationGrating();
}

// this function wraps optUseAdaptive()
// @how to use
// OptUseAdaptive()
int MESH_OptUseAdaptive(lua_State *L){
  SimulationGrating* s = luaW_check<SimulationGrating>(L, 1);
  s->optUseAdaptive();
  return 1;
}

/*======================================================*/
// constructor for the 2D pattern
/*=======================================================*/
SimulationPattern* MESH_SimulationPattern_New(lua_State* L){
  return new SimulationPattern();
}

/*======================================================*/
// tables for the three clases
/*=======================================================*/
static luaL_Reg character_metatable_Simulation[] = {
	{ "SetPeriodicity", MESH_SetPeriodicity },
	{ "AddMaterial", MESH_AddMaterial },
  { "SetMaterial", MESH_SetMaterial, },
  { "SetLayer", MESH_SetLayer, },
	{ "AddLayer", MESH_AddLayer },
  { "SetLayerThickness", MESH_SetLayerThickness, },
  { "AddLayerCopy", MESH_AddLayerCopy, },
  { "DeleteLayer", MESH_DeleteLayer, },
  { "SetSourceLayer", MESH_SetSourceLayer },
  { "SetProbeLayer", MESH_SetProbeLayer },
  { "SetLayerPatternGrating", MESH_SetLayerPatternGrating },
  { "SetLayerPatternRectangle", MESH_SetLayerPatternRectangle },
  { "SetLayerPatternCircle", MESH_SetLayerPatternCircle },
  { "SetGx", MESH_SetGx },
  { "SetGy", MESH_SetGy },
  { "GetPhi", MESH_GetPhi },
  { "GetOmega", MESH_GetOmega },
  { "GetNumOfOmega", MESH_GetNumOfOmega },
  { "GetPhiAtKxKy", MESH_GetPhiAtKxKy },
  { "GetSysInfo", MESH_GetSysInfo },
  { "OptPrintIntermediate", MESH_OptPrintIntermediate },
  { "OptUseInverseRule", MESH_OptUseInverseRule },
  { "OptUseNaiveRule", MESH_OptUseNaiveRule },
  { "BuildRCWA", MESH_BuildRCWA },
  { "SetThread", MESH_SetThread },
  { "SetKxIntegral", MESH_SetKxIntegral },
  { "SetKxIntegralSym", MESH_SetKxIntegralSym },
  { "SetKyIntegral", MESH_SetKyIntegral },
  { "SetKyIntegralSym", MESH_SetKyIntegralSym },
  { "IntegrateKxKy", MESH_IntegrateKxKy },
  { "IntegrateKxKyMPI", MESH_IntegrateKxKyMPI },
  { "OutputStructurePOVRay", MESH_OutputStructurePOVRay },
	{NULL, NULL}
};

static luaL_Reg character_metatable_SimulationPlanar[] = {
  { "OptUseQuadgl", MESH_OptUseQuadgl },
	{ "OptUseQuadgk", MESH_OptUseQuadgk },
  { "SetKParallelIntegral", MESH_SetKParallel },
  { "GetPhiAtKParalle", MESH_GetPhiAtKParallel },
  { "GetPhi", MESH_GetPhi_Planar },
  { "IntegrateKParallel", MESH_IntegrateKParallel },
	{ NULL, NULL}
};

static luaL_Reg character_metatable_SimulationGrating[] = {
  { "OptUseAdaptive", MESH_OptUseAdaptive },
  { NULL, NULL}
};

/*======================================================*/
// special function to return the values of constants
/*=======================================================*/
static int MESH_Constants(lua_State *L){
  std::vector<const char*> properties = {"pi", "k_B", "eps_0", "m_e", "eV", "mu_0", "h", "h_bar", "c_0", "q", "sigma"};
  std::vector<double> vals = {constants.pi, constants.k_B, constants.eps_0, constants.m_e, constants.eV, constants.mu_0,
    constants.h, constants.h_bar, constants.c_0, constants.q, constants.sigma};
  int tableLen = properties.size();
  lua_createtable(L, 0, tableLen);
  for(int i = 0; i < tableLen; i++){
    lua_pushstring(L, properties[i]);
    lua_pushnumber(L, vals[i]);
    lua_settable(L, -3);
  }
  return 1;
};

/*======================================================*/
// registering the classes
/*=======================================================*/
static int luaopen_Simulation(lua_State *L){
	luaW_register<Simulation>(L, "Simulation", NULL, character_metatable_Simulation, NULL);
	luaW_register<SimulationPlanar>(L, "SimulationPlanar", NULL, character_metatable_SimulationPlanar, MESH_SimulationPlanar_New);
  luaW_register<SimulationGrating>(L, "SimulationGrating", NULL, character_metatable_SimulationGrating, MESH_SimulationGrating_New);
  luaW_register<SimulationPattern>(L, "SimulationPattern", NULL, NULL, MESH_SimulationPattern_New);
	luaW_extend<SimulationPlanar, Simulation>(L);
  luaW_extend<SimulationGrating, Simulation>(L);
  luaW_extend<SimulationPattern, Simulation>(L);

  lua_register(L, "Constants", MESH_Constants);
	return 1;
}

/*======================================================*/
// customized new state function
/*=======================================================*/
lua_State* new_MESH_lua_State(){
  lua_State *L = luaL_newstate();
  luaL_requiref(L, "MESH", &luaopen_Simulation, 1);
  lua_pop(L, 1);

  luaL_openlibs(L);
  #ifdef HAVE_MPI
    luaL_requiref(L, "MPI", &luaopen_MPI, 1); lua_pop(L, 1);
    luaL_requiref(L, "buffer", &luaopen_buffer, 1); lua_pop(L, 1);
  #endif
  return L;
}
/*======================================================*/
// information about the package
/*=======================================================*/
void usage(){
  std::cout << "========================================================" << std::endl;
  std::cout << "mesh -h for help" << std::endl;
  std::cout << "mesh -v for version information" << std::endl;
	std::cout << "mesh [input-file] to run a file" << std::endl;
    std::cout << "========================================================" << std::endl;
}
void version(){
  std::cout << "========================================================" << std::endl;
  std::cout << "Copyright (C) 2016-2018, and GNU GPL'd, by Kaifeng Chen." << std::endl;
	std::cout << "Multilayer Electromagnetic Solver for Heat transfer (MESH)" << std::endl;
	std::cout << "Version: " << PACKAGE_VERSION << std::endl;
  std::cout << "Email: " << PACKAGE_BUGREPORT << std::endl;
  #if defined(_OPENMP) && ! defined(HAVE_MPI)
	 std::cout << "With Openmp support." << std::endl;
  #endif
  #ifdef HAVE_MPI
    std::cout << "With MPI support." << std::endl;
  #endif
  std::cout << "========================================================" << std::endl;
}

/*======================================================*/
// main function initializing everything
/*=======================================================*/
int main(int argc, char *argv[]){
  if(argc <= 1){
    std::cerr << "Please type input file name!" << std::endl;
    throw UTILITY::UnknownArgException("Please type input file name!");
    return 0;
  }
  if(argc >= 3){
    std::cerr << "Please only read in one file a time!" << std::endl;
    throw UTILITY::UnknownArgException("Please only read in one file a time!");
    return 0;
  }

  int c;
  if((c = getopt(argc, argv, "vh")) != -1){
    switch(c){
      case 'v':{
        version();
        return 1;
      }
      case 'h':{
        usage();
        return 1;
      }
      case '?':{
        fprintf(stderr, "Unknown option -%c.\n", optopt);
        usage();
        return 0;
      }
      default: abort();
    }
  }

	lua_State *L = new_MESH_lua_State();

  if (luaL_dofile(L, argv[1])){
  	  cout << lua_tostring(L, -1) << endl;
  }
  lua_close(L);

  return 0;

}