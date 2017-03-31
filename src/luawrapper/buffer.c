

#include <stdlib.h>
#include <string.h>

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

enum { 
  BUFFER_TYPE_INT=1,
  BUFFER_TYPE_CHAR=2,
  BUFFER_TYPE_FLOAT=3,
  BUFFER_TYPE_DOUBLE=4,
} ;

// -----------------------------------------------------------------------------
// buffer
// -----------------------------------------------------------------------------
static void buf_push_buffer(lua_State *L, const void *p, int size)
{
  void *newbuf = lua_newuserdata(L, size);
  if (newbuf == NULL) {
    luaL_error(L, "buffer: not enough memory to allocate buffer");
  }
  if (p != NULL) {
    memcpy(newbuf, p, size);
  }
  else {
    memset(newbuf, 0, size);
  }
  luaL_setmetatable(L, "buffer");
}
static int buffer_new_buffer(lua_State *L)
{
  int type = lua_type(L, 1);
  int N;
  const void *buf = NULL;
  switch (type) {
  case LUA_TNUMBER:
    N = luaL_checkinteger(L, 1);
    buf = NULL;
    break;
  case LUA_TSTRING:
    N = lua_rawlen(L, 1);
    buf = lua_tostring(L, 1);
    break;
  default:
    N = 0;
    luaL_error(L, "buffer: argument must be either number or string");
  }
  buf_push_buffer(L, buf, N);
  return 1;
}

static int buffer_sizeof(lua_State *L)
{
  int T = luaL_checkinteger(L, 1);
  int s = 0;
  switch (T) {
  case BUFFER_TYPE_CHAR: s = sizeof(char); break;
  case BUFFER_TYPE_INT: s = sizeof(int); break;
  case BUFFER_TYPE_FLOAT: s = sizeof(float); break;
  case BUFFER_TYPE_DOUBLE: s = sizeof(double); break;
  default:
    luaL_error(L, "buffer: unknown data type specification");
  }
  lua_pushnumber(L, s);
  return 1;
}

static int buffer_get_typed(lua_State *L)
{
  const char *buf = (const char *)luaL_checkudata(L, 1, "buffer");
  int T = luaL_checkinteger(L, 2); // type
  int n = luaL_checkinteger(L, 3); // index
  int N = lua_rawlen(L, 1); // buffer size
  int offset;

#define CASE(t)								\
  do {									\
    offset = n * sizeof(t);						\
    if (offset >= N) luaL_error(L, "buffer: index out of range");	\
    lua_pushnumber(L, *((t*)(buf + offset)));				\
  } while(0)								\

  switch (T) {
  case BUFFER_TYPE_CHAR: CASE(char); break;
  case BUFFER_TYPE_INT: CASE(int); break;
  case BUFFER_TYPE_FLOAT: CASE(float); break;
  case BUFFER_TYPE_DOUBLE: CASE(double); break;
  default:
    luaL_error(L, "buffer: unknown data type specification");
  }
#undef CASE

  return 1;
}
static int buffer_set_typed(lua_State *L)
{
  const char *buf = (const char *)luaL_checkudata(L, 1, "buffer");
  int T = luaL_checkinteger(L, 2); // type
  int n = luaL_checkinteger(L, 3); // index
  double v = luaL_checknumber(L, 4); // value
  int N = lua_rawlen(L, 1); // buffer size
  int offset;

#define CASE(t)								\
  do {									\
    offset = n * sizeof(t);						\
    if (offset >= N) luaL_error(L, "buffer: index out of range");	\
    *((t*)(buf + offset)) = v;						\
  } while(0)								\

  switch (T) {
  case BUFFER_TYPE_CHAR: CASE(char); break;
  case BUFFER_TYPE_INT: CASE(int); break;
  case BUFFER_TYPE_FLOAT: CASE(float); break;
  case BUFFER_TYPE_DOUBLE: CASE(double); break;
  default:
    luaL_error(L, "buffer: unknown data type specification");
  }
#undef CASE

  return 0;
}

static int buffer__index(lua_State *L)
{
  const char *buf = (const char *)luaL_checkudata(L, 1, "buffer");
  int n = luaL_checkinteger(L, 2);
  int N = lua_rawlen(L, 1);
  if (n < N) {
    lua_pushnumber(L, buf[n]);
  }
  else {
    lua_pushnil(L);
  }
  return 1;
}
static int buffer__newindex(lua_State *L)
{
  char *buf = (char *)luaL_checkudata(L, 1, "buffer"); // buffer
  int n = luaL_checkinteger(L, 2); // index
  char val = luaL_checkinteger(L, 3); // value
  int N = lua_rawlen(L, 1);  // max index
  if (n < N) {
    buf[n] = val;
  }
  else {
    luaL_error(L, "buffer: index %d out of range on buffer of length %d", n, N);
  }
  return 0;
}
static int buffer__len(lua_State *L)
{
  luaL_checkudata(L, 1, "buffer");
  lua_pushnumber(L, lua_rawlen(L, 1));
  return 1;
}
static int buffer__tostring(lua_State *L)
{
  char *buf = (char *)luaL_checkudata(L, 1, "buffer");
  lua_pushlstring(L, buf, lua_rawlen(L, 1));
  return 1;
}


int luaopen_buffer(lua_State *L)
{
  luaL_Reg buffer_types[] = {
    {"new_buffer", buffer_new_buffer},
    {"sizeof", buffer_sizeof},
    {"get_typed", buffer_get_typed},
    {"set_typed", buffer_set_typed},
    {NULL, NULL}};

  luaL_Reg buffer_meta[] = {
    {"__index", buffer__index},
    {"__newindex", buffer__newindex},
    {"__len", buffer__len},
    {"__tostring", buffer__tostring},
    {NULL, NULL}};

  luaL_newmetatable(L, "buffer");
  luaL_setfuncs(L, buffer_meta, 0);
  lua_pop(L, 1);

  lua_newtable(L);
  luaL_setfuncs(L, buffer_types, 0);

#define REG_NUMBER(s,t) lua_pushnumber(L, s); lua_setfield(L, -2, t);
  REG_NUMBER(BUFFER_TYPE_CHAR, "char");
  REG_NUMBER(BUFFER_TYPE_INT, "int");
  REG_NUMBER(BUFFER_TYPE_FLOAT, "float");
  REG_NUMBER(BUFFER_TYPE_DOUBLE, "double");
#undef REG_NUMBER

  return 1;
}
