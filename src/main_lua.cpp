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
#define LUA_COMPAT_MODULE
#include "Mesh.h"
#include "System.h"
#include "luawrapper/luawrapper.hpp"

using namespace MESH;
using namespace SYSTEM;

Ptr<Material> MESH_MaterialNew(lua_State *L){
  const char* name = luaL_checkstring(L, 1);
  return  Material::instanceNew(name);
}

static int MESH_SetOmega(lua_State* L){
  Ptr<Material> material = luaW_check<Material>(L, 1);
  int n = lua_rawlen(L, 2);
  double omegaList[n];
  for(int i = 0; i < n; i++){
    lua_pushinteger(L, i+1);
    lua_gettable(L, 2);
    lua_rawgeti(L, -1, 1);
    omegaList[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  material->setOmega(omegaList, n);
  return 0;
}

static int MESH_SetEpsilon(lua_State* L){
  Ptr<Material> material = luaW_check<Material>(L, 1);
  int n = lua_rawlen(L, 2);
  dcomplex epsilonList[n];
  for(int i = 0; i < n; i++){
    lua_pushinteger(L, i+1);
    lua_gettable(L, 2);
    lua_rawgeti(L, -1, 1);

    lua_rawgeti(L, -1, 1);
    double real = lua_tonumber(L, -1);
    lua_pop(L, 1);
    double imag = lua_tonumber(L, -1);
    lua_pop(L, 1);
    epsilonList[i] = dcomplex(real, imag);
  }

  material->setEpsilon(epsilonList, n);
  return 0;
}

int luaopen_MESH(lua_State* L){

  return 1;
}

lua_State* new_MESH_lua_state(){
  lua_State* L = luaL_newstate();

}

void usage(){
	std::cout << "MESH [input-file]" << std::endl;
}
void version(){
	std::cout << "Multilayer Electromagnetic Solver for Heat transfer (MESH)" << std::endl;
	std::cout << "Version " << "\t" << PACKAGE_VERSION << std::endl;
	std::cout << "With MPI support" << std::endl;
}

int main(int argc, char *argv[]){
  if(argc <= 1){
    throw UTILITY::UnknownArgException("Please type input file name!");
    return 0;
  }
  if(argc >= 3){
    throw UTILITY::UnknownArgException("Please only read in one file a time!");
    return 0;
  }
   lua_State *L;
   L = new_MESH_lua_state();
   return 0;

}