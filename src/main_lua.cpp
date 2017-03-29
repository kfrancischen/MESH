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

extern "C"
{
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}
#include "luawrapper/luawrapper.hpp"
#include "luawrapper/luawrapperutil.hpp"

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
int MESH_SetPeriodicity(lua_State* L){
  int n = lua_gettop(L);
	if(n != 3 || n != 4){
		return luaL_error(L, "expecting 2 or 3 arguments");
	}
	Simulation* s = luaW_check<Simulation>(L, 1);
	double p1 = luaL_checknumber(L, 2);
	if(n == 3){
		s->setPeriodicity(p1);
	}
  else{
  	double p2 = luaL_checknumber(L, 3);
  	s->setPeriodicity(p1, p2);
  }
	return 1;
}

int MESH_AddMaterial(lua_State *L){
	Simulation* s = luaW_check<Simulation>(L, 1);
	std::string name = luaU_check<std::string>(L, 2);
	std::string infile = luaU_check<std::string>(L, 3);
	s->addMaterial(name, infile);
	return 1;
}

int MESH_SetMaterial(lua_State* L){
  // TODO
  return 1;
}

int MESH_AddLayer(lua_State *L){
	Simulation* s = luaW_check<Simulation>(L, 1);
	std::string name = luaU_check<std::string>(L, 2);
	double thick = luaL_checknumber(L, 3);
	std::string materialName = luaU_check<std::string>(L, 4);
	s->addLayer(name, thick, materialName);
	return 1;
}

int MESH_SetLayer(lua_State* L){
  Simulation* s = luaW_check<Simulation>(L, 1);
	std::string name = luaU_check<std::string>(L, 2);
	double thick = luaL_checknumber(L, 3);
	std::string materialName = luaU_check<std::string>(L, 4);
	s->setLayer(name, thick, materialName);
	return 1;
}

int MESH_SetLayerThickness(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  double thick = luaL_checknumber(L, 3);
  s->setLayerThickness(name, thick);
  return 1;
}

int MESH_AddLayerCopy(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  std::string materialName = luaU_check<std::string>(L, 3);
  s->addLayerCopy(name, materialName);
  return 1;
}

int MESH_DeleteLayer(lua_State* L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  s->deleteLayer(name);
  return 1;
}

int MESH_SetLayerPatternGrating(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string layerName = luaU_check<std::string>(L, 2);
  std::string materialName = luaU_check<std::string>(L, 3);
  double center = luaL_checknumber(L, 4);
  double width = luaL_checknumber(L, 5);
  s->setLayerPatternGrating(layerName, materialName, center, width);
  return 1;
}

int MESH_SetLayerPatternRectangle(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string layerName = luaU_check<std::string>(L, 2);
  std::string materialName = luaU_check<std::string>(L, 3);
  double vals[4];
  for(int i = 0; i < 2; i++){
    for(int j= 0; j < 2; j++){
      lua_pushinteger(L, j+1);
      lua_gettable(L, 4+i);
      vals[2*i+j] = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }
  }
  s->setLayerPatternRectangle(layerName, materialName, vals[0], vals[1], vals[2], vals[3]);
  return 1;
}

int MESH_SetLayerPatternCircle(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  std::string layerName = luaU_check<std::string>(L, 2);
  std::string materialName = luaU_check<std::string>(L, 3);
  double radius = lua_tonumber(L, 5);
  double vals[2];
  for(int i = 0; i < 2; i++){
    lua_pushinteger(L, i+1);
    lua_gettable(L, 4);
    vals[i] = lua_tonumber(L, 1);
    lua_pop(L, 1);
  }
  s->setLayerPatternCircle(layerName, materialName, vals[0], vals[1], radius);
  return 1;
}

int MESH_SetGx(lua_State* L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  int n = luaL_checkinteger(L, 2);
  s->setGx(n);
  return 1;
}

int MESH_SetGy(lua_State* L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  int n = luaL_checkinteger(L, 2);
  s->setGy(n);
  return 1;
}

int MESH_SetSourceLayer(lua_State* L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  s->setSourceLayer(name);
  return 1;
}

int MESH_SetProbeLayer(lua_State* L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  s->setProbeLayer(name);
  return 1;
}

int MESH_SaveToFile(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  std::string name = luaU_check<std::string>(L, 2);
  s->saveToFile(name);
  return 1;
}

int MESH_SetThread(lua_State *L){
  Simulation *s = luaW_check<Simulation>(L, 1);
  int thread = luaL_checkinteger(L, 2);
  s->setThread(thread);
  return 1;
}

int MESH_BuildRCWA(lua_State *L){
  Simulation* s = luaW_check<Simulation>(L, 1);
  s->buildRCWA();
  return 1;
}

int MESH_GetPhi(lua_State *L){
  // TODO
  return 1;
}

int MESH_GetOmega(lua_State *L){
  // TODO
  return 1;
}

int MESH_GetNumOfOmega(lua_State *L){
  // TODO
  return 1;
}

int MESH_GetPhiAtKxKy(lua_State *L){
  // TODO
  return 1;
}

int MESH_GetSysInfo(lua_State *L){
  // TODO
  return 1;
}

int MESH_OptUseInverseRule(lua_State *L){
  // TODO
  return 1;
}

int MESH_OptUseNaiveRule(lua_State *L){
  // TODO
  return 1;
}

int MESH_OptPrintIntermediate(lua_State *L){
  // TODO
  return 1;
}

int MESH_SetKxIntegral(lua_State *L){
  // TODO
  return 1;
}

int MESH_SetKyIntegral(lua_State *L){
  // TODO
  return 1;
}

int MESH_SetKxIntegralSym(lua_State *L){
  // TODO
  return 1;
}

int MESH_SetKyIntegralSym(lua_State *L){
  // TODO
  return 1;
}

int MESH_IntegrateKxKy(lua_State *L){
  // TODO
  return 1;
}

SimulationPlanar* MESH_SimulationPlanar_New(lua_State* L){
	return new SimulationPlanar();
}

int MESH_OptUseQuadgk(lua_State *L){
	SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
	s->optUseQuadgk();
	return 1;
}

int MESH_OptUseQuadgl(lua_State *L){
  int n = lua_gettop(L);
  if(n != 1 || n != 2){
		return luaL_error(L, "expecting 2 or 3 arguments");
	}
	SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
  if(n == 1){
  	s->optUseQuadgl();
  }
  else{
    int degree = luaL_checkinteger(L, 2);
    s->optUseQuadgl(degree);
  }
	return 1;
}

int MESH_SetKParallel(lua_State *L){
  SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
  double end = luaL_checknumber(L, 2);
  s->setKParallelIntegral(end);
  return 1;
}

int MESH_IntegrateKParallel(lua_State* L){
  SimulationPlanar* s = luaW_check<SimulationPlanar>(L, 1);
  s->integrateKParallel();
  return 1;
}

/*======================================================*/
// constructor for the 1D grating
/*=======================================================*/
SimulationGrating* MESH_SimulationGrating_New(lua_State* L){
  return new SimulationGrating();
}

// wrapper for use adaptive function
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
  { "setLayerPatternRectangle", MESH_SetLayerPatternRectangle },
  { "SetLayerPatternCircle", MESH_SetLayerPatternCircle },
  { "SetGx", MESH_SetGx },
  { "SetGy", MESH_SetGy },
  { "GetPhi", MESH_GetPhi },
  { "GetOmega", MESH_GetOmega },
  { "GetNumOfOmega", MESH_GetNumOfOmega },
  { "GetPhiAtKxKy", MESH_GetPhiAtKxKy },
  { "GetSysInfo", MESH_GetSysInfo },
  { "OptUseInverseRule", MESH_OptUseInverseRule },
  { "optUseNaiveRule", MESH_OptUseNaiveRule },
  { "SaveToFile", MESH_SaveToFile },
  { "BuildRCWA", MESH_BuildRCWA },
  { "SetThread", MESH_SetThread },
  { "SetKxIntegral", MESH_SetKxIntegral },
  { "SetKxIntegralSym", MESH_SetKxIntegralSym },
  { "SetKyIntegral", MESH_SetKyIntegral },
  { "setKyIntegralSym", MESH_SetKyIntegralSym },
  { "IntegrateKxKy", MESH_IntegrateKxKy },
	{NULL, NULL}
};

static luaL_Reg character_metatable_SimulationPlanar[] = {
  { "OptUseQuadgl", MESH_OptUseQuadgl },
	{ "OptUseQuadgk", MESH_OptUseQuadgk },
  { "SetKParallelIntegral", MESH_SetKParallel },
  { "IntegrateKParallel", MESH_IntegrateKParallel },
	{ NULL, NULL}
};

static luaL_Reg character_metatable_SimulationGrating[] = {
  { "OptUseAdaptive", MESH_OptUseAdaptive },
  { NULL, NULL}
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
	return 1;
}

/*======================================================*/
// information about the package
/*=======================================================*/
void usage(){
	std::cout << "mesh [input-file]" << std::endl;
}
void version(){
	std::cout << "Multilayer Electromagnetic Solver for Heat transfer (MESH)" << std::endl;
	std::cout << "Version " << "\t" << PACKAGE_VERSION << std::endl;
	std::cout << "With Openmp support" << std::endl;
}

/*======================================================*/
// main function initializing everything
/*=======================================================*/
int main(int argc, char *argv[]){
  if(argc <= 1){
    throw UTILITY::UnknownArgException("Please type input file name!");
    return 0;
  }
  if(argc >= 3){
    throw UTILITY::UnknownArgException("Please only read in one file a time!");
    return 0;
  }
   //lua_State *L;
	lua_State *L = luaL_newstate();
	luaL_openlibs(L);
	luaopen_Simulation(L);
	if (luaL_dofile(L, argv[1])){
	   cout << lua_tostring(L, -1) << endl;
	}
  lua_close(L);
  return 0;

}