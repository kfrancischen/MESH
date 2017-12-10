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

/**********************
 * A few word from Kaifeng Chen
 * In writing this file, I am not sure how to do the inheritance
 * If someone whoe is reading this knows how to do it, please update the whole file for me
 * Thanks!
 ********************/
#include "setup.h"
#include <Python.h>

#ifdef __cplusplus
extern "C" {
#endif

using namespace MESH;

struct module_state {
    PyObject *error;
};
/*======================================================*/
// checkers
/*=======================================================*/
static int CheckPyNumber(PyObject *obj){
	return PyFloat_Check(obj) || PyLong_Check(obj)
#if PY_MAJOR_VERSION < 3
		|| PyInt_Check(obj)
#endif
	;
}
static double AsNumberPyNumber(PyObject *obj){
	if(PyFloat_Check(obj)){
		return PyFloat_AsDouble(obj);
	}else if(PyLong_Check(obj)){
		return (double)PyLong_AsLong(obj);
	}
#if PY_MAJOR_VERSION < 3
	else if(PyInt_Check(obj)){
		return (double)PyInt_AsLong(obj);
	}
#endif
	return -1.0;
}
static int CheckPyComplex(PyObject *obj){
	return CheckPyNumber(obj) || PyComplex_Check(obj);
}
static void AsComplexPyComplex(PyObject *obj, double *re, double *im){
	if(CheckPyNumber(obj)){
		*re = AsNumberPyNumber(obj);
		*im = 0;
	}else if(PyComplex_Check(obj)){
		*re = PyComplex_RealAsDouble(obj);
		*im = PyComplex_ImagAsDouble(obj);
	}else{
		*re = -1.;
		*im = 0.;
	}
}
static PyObject *FromIntPyDefInt(int i){
#if PY_MAJOR_VERSION < 3
	return PyInt_FromLong(i);
#else
	return PyLong_FromLong(i);
#endif
}


/*======================================================*/
// converters
/*=======================================================*/
struct epsilon_converter_data{
	int eps_type; /* 0 = scalar, 1 = diagonal, 2 = tensor */
  int num_omega;
	std::vector<double> epsilon;
};
int epsilon_converter(PyObject *obj, struct epsilon_converter_data *data){
  data->num_omega = PyTuple_Size(obj);
  for(int i = 0; i < data->num_omega; i++){
    PyObject* pi = PyTuple_GetItem(obj, 0);
    // case for scalar
    if(PyTuple_Check(pi)){
      switch(PyTuple_Size(pi)){
        case 2:{
          data->eps_type = 0;
          break;
        }
        case 6:{
          data->eps_type = 1;
          break;
        }
        default:{
          data->eps_type = 2;
          break;
        }
      }
      for(int j = 0; j < PyTuple_Size(pi); j++){
        PyObject* pj = PyTuple_GetItem(pi, j);
        if(CheckPyComplex(pj)){
          double real, imag;
          AsComplexPyComplex(pj, &real, &imag);
          data->epsilon.push_back(real);
          data->epsilon.push_back(imag);
        }
        else{
          PyErr_SetString(PyExc_TypeError, "Wrong type of tensor");
          return 0;
        }
      }
    }
  }

	return 1;
}

struct interpolator_converter_data{
  int num_omega;
  std::vector<double> omega;
  std::vector<std::vector<double>> epsilon;
};

int interpolator_converter(PyObject* obj, struct interpolator_converter_data *data){
  data->num_omega = PyTuple_Size(obj);
  
  for(int i = 0; i < data->num_omega; i++){
    PyObject* pi = PyTuple_GetItem(obj, i);
    PyObject* pi_0 = PyTuple_GetItem(pi, 0);
    if(CheckPyNumber(pi_0)){
      data->omega.push_back(AsNumberPyNumber(pi_0));
    }
    else{
      PyErr_SetString(PyExc_TypeError, "Wrong type of omega value");
      return 0;
    }
    PyObject* pi_other = PyTuple_GetItem(pi, 1);
    std::vector<double> singleEpsilon;
    for(int j = 0; j < PyTuple_Size(pi_other); j++){
      PyObject* pj = PyTuple_GetItem(pi_other, j);
      if(CheckPyNumber(pj)){
        singleEpsilon.push_back(AsNumberPyNumber(pj));
      }
      else{
        PyErr_SetString(PyExc_TypeError, "Wrong type of omega value");
        return 0;
      }
    }
    data->epsilon.push_back(singleEpsilon);
  }
  return 1;
}


struct polygon_converter_data{
	int nvert;
	std::vector<double> vert;
};

int polygon_converter(PyObject *obj, struct polygon_converter_data *data){
	if(!PyTuple_Check(obj)){
		return 0;
	}
  data->nvert = PyTuple_Size(obj);
  for(int i = 0; i < data->nvert; i++){
    PyObject* pi = PyTuple_GetItem(obj, i);
    if(PyTuple_Check(pi) && (PyTuple_Size(pi) == 2)){
			for(int j = 0; j < 2; ++j){
				PyObject *pj = PyTuple_GetItem(pi, j);
				if(CheckPyNumber(pj)){
          data->vert.push_back(AsNumberPyNumber(pj));
				}
        else{
					PyErr_SetString(PyExc_TypeError, "Polygon tensor must be a list of coordinate pairs");
					return 0;
				}
			}
    }
    else{
      PyErr_SetString(PyExc_TypeError, "Polygon tensor must be a list of coordinate pairs");
		return 0;
    }
  }
	return 1;
}

struct coordinate_converter_data{
  double coor[3];
};

int coordinate_converter(PyObject *obj, struct coordinate_converter_data *data){
	if(!PyTuple_Check(obj)){
		return 0;
	}
  for(int i = 0; i < 3; i++){
    PyObject* pi = PyTuple_GetItem(obj, i);
    if(CheckPyNumber(pi)){
          data->coor[i] = AsNumberPyNumber(pi);
    }
    else{
      PyErr_SetString(PyExc_TypeError, "Wrong coordinate inputs");
      return 0;
    }
  }
  return 1;
}

/*======================================================*/
// wrapper for Interpolator
/*=======================================================*/
typedef struct {
  PyObject_HEAD;
  Interpolator* interpolator;
} MESH_Interpolator;

/* DEALLOC */
static void MESH_Interpolator_dealloc(MESH_Interpolator* self){
  delete self->interpolator;
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/* NEW */
static PyObject *MESH_Interpolator_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"vals", NULL };
	struct interpolator_converter_data interpolator_data;

  if(!PyArg_ParseTupleAndKeywords(args, kwds, "O&:Interpolator_new", kwlist, &interpolator_converter, &interpolator_data)){
    return NULL;
  }
  
  MESH_Interpolator *self;
  self = (MESH_Interpolator*)type->tp_alloc(type, 0);
  if(self != NULL){
    self->interpolator = new Interpolator(interpolator_data.omega, interpolator_data.epsilon);
  }
  return (PyObject*) self;
}

static PyObject *MESH_Interpolator_Get(MESH_Interpolator *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"x", NULL };
  double x;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "d:Get", kwlist, &x)){ 
    return NULL; 
  }
  std::vector<double> result = self->interpolator->getVal(x);
  PyObject* output = PyTuple_New(result.size());
  for(size_t i = 0; i < result.size(); i++){
    PyTuple_SetItem(output, i, PyFloat_FromDouble(result[i]));
  }
  return output;
}

/* METHOD TABLE */
static PyMethodDef MESH_Interpolator_methods[] = {
  {"Get",  (PyCFunction) MESH_Interpolator_Get, METH_VARARGS | METH_KEYWORDS, "Interpolating at a given value x"},
  {NULL}  /* Sentinel */
};

/* TYPE ... whatever */
static PyTypeObject MESH_Interpolator_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "MESH_Interpolator",          /* tp_name */
  sizeof(MESH_Interpolator),    /* tp_basicsize */
  0,                         /* tp_itemsize */
  (destructor)MESH_Interpolator_dealloc, /* tp_dealloc */
  0,                         /* tp_print */
  0,                         /* tp_getattr */
  0,                         /* tp_setattr */
  0,                         /* tp_reserved */
  0,                         /* tp_repr */
  0,                         /* tp_as_number */
  0,                         /* tp_as_sequence */
  0,                         /* tp_as_mapping */
  0,                         /* tp_hash  */
  0,                         /* tp_call */
  0,                         /* tp_str */
  0,                         /* tp_getattro */
  0,                         /* tp_setattro */
  0,                         /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT,        /* tp_flags */ 
  "MESH_Interpolator",       /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,                         /* tp_richcompare */
  0,                         /* tp_weaklistoffset */
  0,                         /* tp_iter */
  0,                         /* tp_iternext */
  MESH_Interpolator_methods, /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                          /* tp_base */ 
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  0,                         /* tp_init */
  0,                         /* tp_alloc */
  MESH_Interpolator_new, /* tp_new */
};

/*======================================================*/
// wrapper for planar
/*=======================================================*/
/* OBJECT */
typedef struct {
    PyObject_HEAD; // <----- PUTTING THIS FIRST INHERITS THE BASE PYTHON CLASS!!!
    SimulationPlanar* s;
} MESH_SimulationPlanar;


/* DEALLOC */
static void MESH_SimulationPlanar_dealloc(MESH_SimulationPlanar* self){
  // Nothing to deallocate. Use this in case the base
  // struct contains a pointer or so.
  delete self->s;
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/* NEW */
static PyObject *MESH_SimulationPlanar_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
  MESH_SimulationPlanar *self;
  self = (MESH_SimulationPlanar*)type->tp_alloc(type, 0);
   // Init vars if any ...
   if (self != NULL){
      self->s = new SimulationPlanar();
   }
   return (PyObject*) self;
}

static PyObject* MESH_SimulationPlanar_AddMaterial(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"material_name", (char*)"file_name", NULL };
	const char *materialName, *fileName;
	if(PyArg_ParseTupleAndKeywords(args, kwds, "ss:AddMaterial", kwlist, &materialName, &fileName)){ 
    std::string material_name(materialName), file_name(fileName);
    self->s->addMaterial(material_name, file_name);
  }
  else{
    static char *kwlist[] = { (char*)"material name", (char*)"data", NULL };
    struct interpolator_converter_data interpolator_data;
    const char *materialName;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:AddMaterial", kwlist, &materialName, &interpolator_converter, &interpolator_data)){
      return NULL;
    }

    std::string material_name(materialName);
    self->s->addMaterial(material_name, interpolator_data.omega, interpolator_data.epsilon);
  }
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPlanar_SetMaterial(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"material_name", (char*)"epsilon", NULL };
	const char *materialName;
	struct epsilon_converter_data epsdata;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:SetMaterial", kwlist, &materialName, &epsilon_converter, &epsdata)){ 
    return NULL; 
  }
  std::string material_name(materialName);
  double** epsilon = new double*[epsdata.num_omega];
  std::string type = "scalar";
  int eps_len = 2;
  switch(epsdata.eps_type){
    case 0:{
      // nothing
      break;
    }
    case 1:{
      eps_len = 6;
      type = "diagonal";
      break;
    }
    default:{
      eps_len = 10;
      type = "tensor";
      break;
    }
  }

  for(int i = 0; i < epsdata.num_omega;i++){
    epsilon[i] = new double[eps_len];
    for(int j = 0; j < eps_len; j++){
      epsilon[i][j] = epsdata.epsilon[eps_len*i+j];
    }
  }

  self->s->setMaterial(material_name, epsilon, type);
  for(int i = 0; i < epsdata.num_omega; i++){
    delete [] epsilon[i];
  }
  delete [] epsilon;
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPlanar_AddLayer(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"thickness", (char*)"material_name", NULL };
	const char *materialName, *layerName;
  double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayer", kwlist, &layerName, &thickness, &materialName)){ 
      return NULL; 
  }
  std::string material_name(materialName), layer_name(layerName);
  self->s->addLayer(layer_name, thickness, material_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_SetLayer(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"thickness", (char*)"material_name", NULL };
	const char *materialName, *layerName;
  double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayer", kwlist, &layerName, &thickness, &materialName)){ 
      return NULL; 
  }
  std::string material_name(materialName), layer_name(layerName);
  self->s->setLayer(layer_name, thickness, material_name);
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPlanar_SetLayerThickness(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"thickness", NULL };
	const char *layerName;
  double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sd:AddLayer", kwlist, &layerName, &thickness)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->setLayerThickness(layer_name, thickness);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_AddLayerCopy(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"copy_layer_name", NULL };
	const char *layerName, *copyLayerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss:AddLayer", kwlist, &layerName, &copyLayerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName), copy_layer_name(copyLayerName);
  self->s->addLayerCopy(layer_name, copy_layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_DeleteLayer(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", NULL };
	const char *layerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:DeleteLayer", kwlist, &layerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->deleteLayer(layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_SetSourceLayer(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", NULL };
	const char *layerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:SetSourceLayer", kwlist, &layerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->setSourceLayer(layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_SetProbeLayer(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", NULL };
	const char *layerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:SetProbeLayer", kwlist, &layerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->setProbeLayer(layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_SetProbeLayerZCoordinate(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"z", NULL };
	double z;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "d:SetProbeLayerZCoordinate", kwlist, &z)){ 
      return NULL; 
  }
  self->s->setProbeLayerZCoordinate(z);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_SetThread(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"num_thread", NULL };
	int thread;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "i:SetThread", kwlist, &thread)){ 
      return NULL; 
  }
  self->s->setThread(thread);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_InitSimulation(MESH_SimulationPlanar *self, PyObject *args){
  self->s->initSimulation();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_GetPhi(MESH_SimulationPlanar *self, PyObject *args){
  double* phi = self->s->getPhi();
  int num_omega = self->s->getNumOfOmega();
  PyObject* phi_value = PyTuple_New(num_omega);
  for(int i = 0; i < num_omega; i++){
    PyTuple_SetItem(phi_value, i, PyFloat_FromDouble(phi[i]));
  }
  return phi_value;
}

static PyObject* MESH_SimulationPlanar_GetOmega(MESH_SimulationPlanar *self, PyObject *args){
  double* omega = self->s->getOmega();
  int num_omega = self->s->getNumOfOmega();
  PyObject* omega_value = PyTuple_New(num_omega);
  for(int i = 0; i < num_omega; i++){
    PyTuple_SetItem(omega_value, i, PyFloat_FromDouble(omega[i]));
  }
  return omega_value;
}

static PyObject* MESH_SimulationPlanar_GetEpsilon(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"coordinate", NULL};
  int omega_index;
  struct coordinate_converter_data coordata;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "iO&:GetEpsilon", kwlist, &omega_index, &coordinate_converter, &coordata)){ 
    return NULL; 
  }
  double* epsilon = new double[10];
  self->s->getEpsilon(omega_index, coordata.coor, epsilon);
  
  PyObject* eps_value = PyTuple_New(10);
  for(int i = 0; i < 10; i++){
    PyTuple_SetItem(eps_value, i, PyFloat_FromDouble(epsilon[i]));
  }
  delete [] epsilon;
  return eps_value;
}

static PyObject* MESH_SimulationPlanar_OutputLayerPatternRealization(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"layer_name", (char*)"Nu", (char*)"Nv", (char*)"file_name", NULL};
  int omega_index, Nu, Nv;
  char* layerName, *fileName = (char*)"";
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "isii|s:OutputLayerPatternRealization", kwlist, &omega_index, &layerName, &Nu, &Nv, &fileName)){ 
    return NULL; 
  }
  std::string layer_name(layerName), file_name(fileName);
  if(file_name == ""){
    self->s->outputLayerPatternRealization(omega_index, layer_name, Nu, Nv);
  }
  else{
    self->s->outputLayerPatternRealization(omega_index, layer_name, Nu, Nv, file_name);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_GetNumOfOmega(MESH_SimulationPlanar *self, PyObject *args){
  int num_omega = self->s->getNumOfOmega();
  return FromIntPyDefInt(num_omega);
}

static PyObject* MESH_SimulationPlanar_GetPhiAtKxKy(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"kx", (char*)"ky", NULL};
  int omega_index;
  double kx, ky = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "id|d:GetPhiAtKxKy", kwlist, &omega_index, &kx, &ky)){ 
    return NULL; 
  }
  double value = self->s->getPhiAtKxKy(omega_index, kx, ky);
  return PyFloat_FromDouble(value);
}

static PyObject* MESH_SimulationPlanar_OutputSysInfo(MESH_SimulationPlanar *self, PyObject *args){
  self->s->outputSysInfo();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_OptPrintIntermediate(MESH_SimulationPlanar *self, PyObject *args){
  self->s->optPrintIntermediate();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_OptOnlyComputeTE(MESH_SimulationPlanar *self, PyObject *args){
  self->s->optOnlyComputeTE();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_OptOnlyComputeTM(MESH_SimulationPlanar *self, PyObject *args){
  self->s->optOnlyComputeTM();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_SetKxIntegral(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKxIntegral", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKxIntegral(num_points);
  }
  else{
    self->s->setKxIntegral(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_SetKyIntegral(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKyIntegral", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKyIntegral(num_points);
  }
  else{
    self->s->setKyIntegral(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_SetKxIntegralSym(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKxIntegralSym", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKxIntegralSym(num_points);
  }
  else{
    self->s->setKxIntegralSym(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_SetKyIntegralSym(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKyIntegralSym", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKyIntegralSym(num_points);
  }
  else{
    self->s->setKyIntegralSym(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_IntegrateKxKy(MESH_SimulationPlanar *self, PyObject *args){
  self->s->integrateKxKy();
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPlanar_IntegrateKxKyMPI(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"rank", (char*)"size", NULL};
  int rank, size;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ii:IntegrateKxKyMPI", kwlist, &rank, &size)){ 
    return NULL; 
  }
  self->s->integrateKxKyMPI(rank, size);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_OptUseQuadgk(MESH_SimulationPlanar *self, PyObject *args){
  self->s->optUseQuadgk();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_OptUseQuadgl(MESH_SimulationPlanar *self, PyObject *args){
  self->s->optUseQuadgl();
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPlanar_SetKParallel(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"integral_end", NULL};
  double integral_end;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "d:SetKParallel", kwlist, &integral_end)){ 
    return NULL; 
  }
  self->s->setKParallelIntegral(integral_end);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPlanar_GetPhiAtKParallel(MESH_SimulationPlanar *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"k_parallel", NULL};
  int omega_index;
  double k_parallel;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "id:GetPhiAtKParallel", kwlist, &omega_index, &k_parallel)){ 
    return NULL; 
  }
  return PyFloat_FromDouble(self->s->getPhiAtKParallel(omega_index, k_parallel));
}

static PyObject* MESH_SimulationPlanar_IntegrateKParallel(MESH_SimulationPlanar *self, PyObject *args){
  self->s->integrateKParallel();
  Py_RETURN_NONE;
}


/* METHOD TABLE */
static PyMethodDef MESH_SimulationPlanar_methods[] = {
  {"AddMaterial",                   (PyCFunction) MESH_SimulationPlanar_AddMaterial,                   METH_VARARGS | METH_KEYWORDS, "Adding new material to the simulation"},
  {"SetMaterial",                   (PyCFunction) MESH_SimulationPlanar_SetMaterial,                   METH_VARARGS | METH_KEYWORDS, "Setting material property"},
  {"AddLayer",                      (PyCFunction) MESH_SimulationPlanar_AddLayer,                      METH_VARARGS | METH_KEYWORDS, "Adding new layer to the simulation"},
  {"SetLayer",                      (PyCFunction) MESH_SimulationPlanar_SetLayer,                      METH_VARARGS | METH_KEYWORDS, "Setting layer property"},
  {"SetLayerThickness",             (PyCFunction) MESH_SimulationPlanar_SetLayerThickness,             METH_VARARGS | METH_KEYWORDS, "Setting layer thickness"},
  {"AddLayerCopy",                  (PyCFunction) MESH_SimulationPlanar_AddLayerCopy,                  METH_VARARGS | METH_KEYWORDS, "Making a copy of existing layer"},
  {"DeleteLayer",                   (PyCFunction) MESH_SimulationPlanar_DeleteLayer,                   METH_VARARGS | METH_KEYWORDS, "Deleting an existing layer"},
  {"SetSourceLayer",                (PyCFunction) MESH_SimulationPlanar_SetSourceLayer,                METH_VARARGS | METH_KEYWORDS, "Setting a source layer"},
  {"SetProbeLayer",                 (PyCFunction) MESH_SimulationPlanar_SetProbeLayer,                 METH_VARARGS | METH_KEYWORDS, "Setting the probe layer"},
  {"SetProbeLayerZCoordinate",      (PyCFunction) MESH_SimulationPlanar_SetProbeLayerZCoordinate,      METH_VARARGS | METH_KEYWORDS, "Setting the z-coordinate in the probe layer"},
  {"SetThread",                     (PyCFunction) MESH_SimulationPlanar_SetThread,                     METH_VARARGS | METH_KEYWORDS, "Setting the number of thread"},
  {"InitSimulation",                (PyCFunction) MESH_SimulationPlanar_InitSimulation,                METH_VARARGS | METH_KEYWORDS, "Initializing simulation"},
  {"GetPhi",                        (PyCFunction) MESH_SimulationPlanar_GetPhi,                        METH_VARARGS | METH_KEYWORDS, "Getting the value of Phi"},
  {"GetOmega",                      (PyCFunction) MESH_SimulationPlanar_GetOmega,                      METH_VARARGS | METH_KEYWORDS, "Getting all the omega values"},
  {"GetEpsilon",                    (PyCFunction) MESH_SimulationPlanar_GetEpsilon,                    METH_VARARGS | METH_KEYWORDS, "Getting epsilon at one frequency"},
  {"OutputLayerPatternRealization", (PyCFunction) MESH_SimulationPlanar_OutputLayerPatternRealization, METH_VARARGS | METH_KEYWORDS, "Outputting dielectric reconstruction"},
  {"GetNumOfOmega",                 (PyCFunction) MESH_SimulationPlanar_GetNumOfOmega,                 METH_VARARGS | METH_KEYWORDS, "Getting the number of omega"},
  {"GetPhiAtKxKy",                  (PyCFunction) MESH_SimulationPlanar_GetPhiAtKxKy,                  METH_VARARGS | METH_KEYWORDS, "Getting Phi value at a (kx,ky) pair"},
  {"OutputSysInfo",                 (PyCFunction) MESH_SimulationPlanar_OutputSysInfo,                 METH_VARARGS | METH_KEYWORDS, "Outputting system information"},
  {"OptPrintIntermediate",          (PyCFunction) MESH_SimulationPlanar_OptPrintIntermediate,          METH_VARARGS | METH_KEYWORDS, "Option to output intermediate results"},
  {"OptOnlyComputeTE",              (PyCFunction) MESH_SimulationPlanar_OptOnlyComputeTE,              METH_VARARGS | METH_KEYWORDS, "Option to only compute TE mode"},
  {"OptOnlyComputeTM",              (PyCFunction) MESH_SimulationPlanar_OptOnlyComputeTM,              METH_VARARGS | METH_KEYWORDS, "Option to only compute TM mode"},
  {"SetKxIntegral",                 (PyCFunction) MESH_SimulationPlanar_SetKxIntegral,                 METH_VARARGS | METH_KEYWORDS, "Setting kx integration range"},
  {"SetKyIntegral",                 (PyCFunction) MESH_SimulationPlanar_SetKyIntegral,                 METH_VARARGS | METH_KEYWORDS, "Setting kx integration range"},
  {"SetKxIntegralSym",              (PyCFunction) MESH_SimulationPlanar_SetKxIntegralSym,              METH_VARARGS | METH_KEYWORDS, "Setting kx integration range in symmetric case"},
  {"SetKyIntegralSym",              (PyCFunction) MESH_SimulationPlanar_SetKyIntegralSym,              METH_VARARGS | METH_KEYWORDS, "Setting ky integration range in symmetric case"},
  {"IntegrateKxKy",                 (PyCFunction) MESH_SimulationPlanar_IntegrateKxKy,                 METH_VARARGS | METH_KEYWORDS, "Action to integrate kx and ky"},
  {"IntegrateKxKyMPI",              (PyCFunction) MESH_SimulationPlanar_IntegrateKxKyMPI,              METH_VARARGS | METH_KEYWORDS, "Action to integrate kx and ky using MPI"},
  {"OptUseQuadgk",                  (PyCFunction) MESH_SimulationPlanar_OptUseQuadgk,                  METH_VARARGS | METH_KEYWORDS, "Option to use Quadgk"},
  {"OptUseQuadgl",                  (PyCFunction) MESH_SimulationPlanar_OptUseQuadgl,                  METH_VARARGS | METH_KEYWORDS, "Option to use Quadgk"},
  {"SetKParallelIntegral",          (PyCFunction) MESH_SimulationPlanar_SetKParallel,                  METH_VARARGS | METH_KEYWORDS, "Setting kParallel integration range"},
  {"GetPhiAtKParallel",             (PyCFunction) MESH_SimulationPlanar_GetPhiAtKParallel,             METH_VARARGS | METH_KEYWORDS, "Getting Phi at a particular kParallel"},
  {"IntegrateKParallel",            (PyCFunction) MESH_SimulationPlanar_IntegrateKParallel,            METH_VARARGS | METH_KEYWORDS, "Action to integration kParallel"},
  {NULL}  /* Sentinel */
};

/* TYPE ... whatever */
static PyTypeObject MESH_SimulationPlanar_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "MESH_SimulationPlanar",          /* tp_name */
  sizeof(MESH_SimulationPlanar),    /* tp_basicsize */
  0,                         /* tp_itemsize */
  (destructor)MESH_SimulationPlanar_dealloc, /* tp_dealloc */
  0,                         /* tp_print */
  0,                         /* tp_getattr */
  0,                         /* tp_setattr */
  0,                         /* tp_reserved */
  0,                         /* tp_repr */
  0,                         /* tp_as_number */
  0,                         /* tp_as_sequence */
  0,                         /* tp_as_mapping */
  0,                         /* tp_hash  */
  0,                         /* tp_call */
  0,                         /* tp_str */
  0,                         /* tp_getattro */
  0,                         /* tp_setattro */
  0,                         /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT,        /* tp_flags */ 
  "MESH_SimulationPlanar",   /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,                         /* tp_richcompare */
  0,                         /* tp_weaklistoffset */
  0,                         /* tp_iter */
  0,                         /* tp_iternext */
  MESH_SimulationPlanar_methods,    /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                          /* tp_base */ 
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  0,                         /* tp_init */
  0,                         /* tp_alloc */
  MESH_SimulationPlanar_new, /* tp_new */
};


/*======================================================*/
// wrappaer for 1D grating
/*=======================================================*/
/* OBJECT */
typedef struct {
    PyObject_HEAD; // <----- PUTTING THIS FIRST INHERITS THE BASE PYTHON CLASS!!!
    SimulationGrating* s;
} MESH_SimulationGrating;


/* DEALLOC */
static void MESH_SimulationGrating_dealloc(MESH_SimulationGrating* self){
  // Nothing to deallocate. Use this in case the base
  // struct contains a pointer or so.
  delete self->s;
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/* NEW */
static PyObject *MESH_SimulationGrating_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
  MESH_SimulationGrating* self;
  self = (MESH_SimulationGrating*)type->tp_alloc(type, 0);

   // Init vars if any ...
   if (self != NULL){
      self->s = new SimulationGrating();
   }
   return (PyObject*) self;
}

static PyObject* MESH_SimulationGrating_AddMaterial(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"material_name", (char*)"file_name", NULL };
	const char *materialName, *fileName;
	if(PyArg_ParseTupleAndKeywords(args, kwds, "ss:AddMaterial", kwlist, &materialName, &fileName)){ 
    std::string material_name(materialName), file_name(fileName);
    self->s->addMaterial(material_name, file_name);
  }
  else{
    static char *kwlist[] = { (char*)"material name", (char*)"data", NULL };
    struct interpolator_converter_data interpolator_data;
    const char *materialName;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:AddMaterial", kwlist, &materialName, &interpolator_converter, &interpolator_data)){
      return NULL;
    }

    std::string material_name(materialName);
    self->s->addMaterial(material_name, interpolator_data.omega, interpolator_data.epsilon);
  }
  Py_RETURN_NONE;

}

static PyObject* MESH_SimulationGrating_SetMaterial(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"material_name", (char*)"epsilon", NULL };
	const char *materialName;
	struct epsilon_converter_data epsdata;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:SetMaterial", kwlist, &materialName, &epsilon_converter, &epsdata)){ 
    return NULL; 
  }
  std::string material_name(materialName);
  double** epsilon = new double*[epsdata.num_omega];
  std::string type = "scalar";
  int eps_len = 2;
  switch(epsdata.eps_type){
    case 0:{
      // nothing
      break;
    }
    case 1:{
      eps_len = 6;
      type = "diagonal";
      break;
    }
    default:{
      eps_len = 10;
      type = "tensor";
      break;
    }
  }

  for(int i = 0; i < epsdata.num_omega;i++){
    epsilon[i] = new double[eps_len];
    for(int j = 0; j < eps_len; j++){
      epsilon[i][j] = epsdata.epsilon[eps_len*i+j];
    }
  }

  self->s->setMaterial(material_name, epsilon, type);
  for(int i = 0; i < epsdata.num_omega; i++){
    delete [] epsilon[i];
  }
  delete [] epsilon;
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationGrating_AddLayer(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"thickness", (char*)"material_name", NULL };
	const char *materialName, *layerName;
  double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayer", kwlist, &layerName, &thickness, &materialName)){ 
      return NULL; 
  }
  std::string material_name(materialName), layer_name(layerName);
  self->s->addLayer(layer_name, thickness, material_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetLayer(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"thickness", (char*)"material_name", NULL };
	const char *materialName, *layerName;
  double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayer", kwlist, &layerName, &thickness, &materialName)){ 
      return NULL; 
  }
  std::string material_name(materialName), layer_name(layerName);
  self->s->setLayer(layer_name, thickness, material_name);
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationGrating_SetLayerThickness(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"thickness", NULL };
	const char *layerName;
  double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sd:AddLayer", kwlist, &layerName, &thickness)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->setLayerThickness(layer_name, thickness);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_AddLayerCopy(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"copy_layer_name", NULL };
	const char *layerName, *copyLayerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss:AddLayer", kwlist, &layerName, &copyLayerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName), copy_layer_name(copyLayerName);
  self->s->addLayerCopy(layer_name, copy_layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_DeleteLayer(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", NULL };
	const char *layerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:DeleteLayer", kwlist, &layerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->deleteLayer(layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetSourceLayer(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", NULL };
	const char *layerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:SetSourceLayer", kwlist, &layerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->setSourceLayer(layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetProbeLayer(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", NULL };
	const char *layerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:SetProbeLayer", kwlist, &layerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->setProbeLayer(layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetProbeLayerZCoordinate(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"z", NULL };
	double z;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "d:SetProbeLayerZCoordinate", kwlist, &z)){ 
      return NULL; 
  }
  self->s->setProbeLayerZCoordinate(z);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetThread(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"num_thread", NULL };
	int thread;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "i:SetThread", kwlist, &thread)){ 
      return NULL; 
  }
  self->s->setThread(thread);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_InitSimulation(MESH_SimulationGrating *self, PyObject *args){
  self->s->initSimulation();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_GetPhi(MESH_SimulationGrating *self, PyObject *args){
  double* phi = self->s->getPhi();
  int num_omega = self->s->getNumOfOmega();
  PyObject* phi_value = PyTuple_New(num_omega);
  for(int i = 0; i < num_omega; i++){
    PyTuple_SetItem(phi_value, i, PyFloat_FromDouble(phi[i]));
  }
  return phi_value;
}

static PyObject* MESH_SimulationGrating_GetOmega(MESH_SimulationGrating *self, PyObject *args){
  double* omega = self->s->getOmega();
  int num_omega = self->s->getNumOfOmega();
  PyObject* omega_value = PyTuple_New(num_omega);
  for(int i = 0; i < num_omega; i++){
    PyTuple_SetItem(omega_value, i, PyFloat_FromDouble(omega[i]));
  }
  return omega_value;
}

static PyObject* MESH_SimulationGrating_GetEpsilon(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"coordinate", NULL};
  int omega_index;
  struct coordinate_converter_data coordata;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "iO&:GetEpsilon", kwlist, &omega_index, &coordinate_converter, &coordata)){ 
    return NULL; 
  }
  double* epsilon = new double[10];
  self->s->getEpsilon(omega_index, coordata.coor, epsilon);
  
  PyObject* eps_value = PyTuple_New(10);
  for(int i = 0; i < 10; i++){
    PyTuple_SetItem(eps_value, i, PyFloat_FromDouble(epsilon[i]));
  }
  delete [] epsilon;
  return eps_value;
}

static PyObject* MESH_SimulationGrating_OutputLayerPatternRealization(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"layer_name", (char*)"Nu", (char*)"Nv", (char*)"file_name", NULL};
  int omega_index, Nu, Nv;
  char* layerName, *fileName = (char*)"";
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "isii|s:OutputLayerPatternRealization", kwlist, &omega_index, &layerName, &Nu, &Nv, &fileName)){ 
    return NULL; 
  }
  std::string layer_name(layerName), file_name(fileName);
  if(file_name == ""){
    self->s->outputLayerPatternRealization(omega_index, layer_name, Nu, Nv);
  }
  else{
    self->s->outputLayerPatternRealization(omega_index, layer_name, Nu, Nv, file_name);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_GetNumOfOmega(MESH_SimulationGrating *self, PyObject *args){
  int num_omega = self->s->getNumOfOmega();
  return FromIntPyDefInt(num_omega);
}

static PyObject* MESH_SimulationGrating_GetPhiAtKxKy(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"kx", (char*)"ky", NULL};
  int omega_index;
  double kx, ky = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "id|d:GetPhiAtKxKy", kwlist, &omega_index, &kx, &ky)){ 
    return NULL; 
  }
  double value = self->s->getPhiAtKxKy(omega_index, kx, ky);
  return PyFloat_FromDouble(value);
}

static PyObject* MESH_SimulationGrating_OutputSysInfo(MESH_SimulationGrating *self, PyObject *args){
  self->s->outputSysInfo();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_OptPrintIntermediate(MESH_SimulationGrating *self, PyObject *args){
  self->s->optPrintIntermediate();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_OptOnlyComputeTE(MESH_SimulationGrating *self, PyObject *args){
  self->s->optOnlyComputeTE();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_OptOnlyComputeTM(MESH_SimulationGrating *self, PyObject *args){
  self->s->optOnlyComputeTM();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetKxIntegral(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKxIntegral", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKxIntegral(num_points);
  }
  else{
    self->s->setKxIntegral(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetKyIntegral(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKyIntegral", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKyIntegral(num_points);
  }
  else{
    self->s->setKyIntegral(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetKxIntegralSym(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKxIntegralSym", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKxIntegralSym(num_points);
  }
  else{
    self->s->setKxIntegralSym(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetKyIntegralSym(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKyIntegralSym", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKyIntegralSym(num_points);
  }
  else{
    self->s->setKyIntegralSym(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_IntegrateKxKy(MESH_SimulationGrating *self, PyObject *args){
  self->s->integrateKxKy();
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationGrating_IntegrateKxKyMPI(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"rank", (char*)"size", NULL};
  int rank, size;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ii:IntegrateKxKyMPI", kwlist, &rank, &size)){ 
    return NULL; 
  }
  self->s->integrateKxKyMPI(rank, size);
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationGrating_SetNumOfG(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_G",  NULL};
  int numG;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i:SetNumOfG", kwlist, &numG)){ 
    return NULL; 
  }
  self->s->setNumOfG(numG);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_GetNumOfG(MESH_SimulationGrating *self, PyObject *args){
  int numG = self->s->getNumOfG();
  return FromIntPyDefInt(numG);
}

static PyObject* MESH_SimulationGrating_SetLatticeGrating(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"lattice_len", NULL};
  double lattice_len;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "d:SetLattice", kwlist, &lattice_len)){ 
    return NULL; 
  }
  self->s->setLattice(lattice_len);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationGrating_SetLayerPatternGrating(MESH_SimulationGrating *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"layer_name", (char*)"material_name", (char*)"center", (char*)"width", NULL};
  char* layerName, *materialName;
  double center, width;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ssdd:SetLayerPatternGrating", kwlist, &layerName, &materialName, &center, &width)){ 
    return NULL; 
  }
  std::string layer_name(layerName), material_name(materialName);
  self->s->setLayerPatternGrating(layer_name, material_name, center, width);
  Py_RETURN_NONE;
}


/* METHOD TABLE */
static PyMethodDef MESH_SimulationGrating_methods[] = {
  {"AddMaterial",                   (PyCFunction) MESH_SimulationGrating_AddMaterial,                   METH_VARARGS | METH_KEYWORDS, "Adding new material to the simulation"},
  {"SetMaterial",                   (PyCFunction) MESH_SimulationGrating_SetMaterial,                   METH_VARARGS | METH_KEYWORDS, "Setting material property"},
  {"AddLayer",                      (PyCFunction) MESH_SimulationGrating_AddLayer,                      METH_VARARGS | METH_KEYWORDS, "Adding new layer to the simulation"},
  {"SetLayer",                      (PyCFunction) MESH_SimulationGrating_SetLayer,                      METH_VARARGS | METH_KEYWORDS, "Setting layer property"},
  {"SetLayerThickness",             (PyCFunction) MESH_SimulationGrating_SetLayerThickness,             METH_VARARGS | METH_KEYWORDS, "Setting layer thickness"},
  {"AddLayerCopy",                  (PyCFunction) MESH_SimulationGrating_AddLayerCopy,                  METH_VARARGS | METH_KEYWORDS, "Making a copy of existing layer"},
  {"DeleteLayer",                   (PyCFunction) MESH_SimulationGrating_DeleteLayer,                   METH_VARARGS | METH_KEYWORDS, "Deleting an existing layer"},
  {"SetSourceLayer",                (PyCFunction) MESH_SimulationGrating_SetSourceLayer,                METH_VARARGS | METH_KEYWORDS, "Setting a source layer"},
  {"SetProbeLayer",                 (PyCFunction) MESH_SimulationGrating_SetProbeLayer,                 METH_VARARGS | METH_KEYWORDS, "Setting the probe layer"},
  {"SetProbeLayerZCoordinate",      (PyCFunction) MESH_SimulationGrating_SetProbeLayerZCoordinate,      METH_VARARGS | METH_KEYWORDS, "Setting the z-coordinate in the probe layer"},
  {"SetThread",                     (PyCFunction) MESH_SimulationGrating_SetThread,                     METH_VARARGS | METH_KEYWORDS, "Setting the number of thread"},
  {"InitSimulation",                (PyCFunction) MESH_SimulationGrating_InitSimulation,                METH_VARARGS | METH_KEYWORDS, "Initializing simulation"},
  {"GetPhi",                        (PyCFunction) MESH_SimulationGrating_GetPhi,                        METH_VARARGS | METH_KEYWORDS, "Getting the value of Phi"},
  {"GetOmega",                      (PyCFunction) MESH_SimulationGrating_GetOmega,                      METH_VARARGS | METH_KEYWORDS, "Getting all the omega values"},
  {"GetEpsilon",                    (PyCFunction) MESH_SimulationGrating_GetEpsilon,                    METH_VARARGS | METH_KEYWORDS, "Getting epsilon at one frequency"},
  {"OutputLayerPatternRealization", (PyCFunction) MESH_SimulationGrating_OutputLayerPatternRealization, METH_VARARGS | METH_KEYWORDS, "Outputting dielectric reconstruction"},
  {"GetNumOfOmega",                 (PyCFunction) MESH_SimulationGrating_GetNumOfOmega,                 METH_VARARGS | METH_KEYWORDS, "Getting the number of omega"},
  {"GetPhiAtKxKy",                  (PyCFunction) MESH_SimulationGrating_GetPhiAtKxKy,                  METH_VARARGS | METH_KEYWORDS, "Getting Phi value at a (kx,ky) pair"},
  {"OutputSysInfo",                 (PyCFunction) MESH_SimulationGrating_OutputSysInfo,                 METH_VARARGS | METH_KEYWORDS, "Outputting system information"},
  {"OptPrintIntermediate",          (PyCFunction) MESH_SimulationGrating_OptPrintIntermediate,          METH_VARARGS | METH_KEYWORDS, "Option to output intermediate results"},
  {"OptOnlyComputeTE",              (PyCFunction) MESH_SimulationGrating_OptOnlyComputeTE,              METH_VARARGS | METH_KEYWORDS, "Option to only compute TE mode"},
  {"OptOnlyComputeTM",              (PyCFunction) MESH_SimulationGrating_OptOnlyComputeTM,              METH_VARARGS | METH_KEYWORDS, "Option to only compute TM mode"},
  {"SetKxIntegral",                 (PyCFunction) MESH_SimulationGrating_SetKxIntegral,                 METH_VARARGS | METH_KEYWORDS, "Setting kx integration range"},
  {"SetKyIntegral",                 (PyCFunction) MESH_SimulationGrating_SetKyIntegral,                 METH_VARARGS | METH_KEYWORDS, "Setting kx integration range"},
  {"SetKxIntegralSym",              (PyCFunction) MESH_SimulationGrating_SetKxIntegralSym,              METH_VARARGS | METH_KEYWORDS, "Setting kx integration range in symmetric case"},
  {"SetKyIntegralSym",              (PyCFunction) MESH_SimulationGrating_SetKyIntegralSym,              METH_VARARGS | METH_KEYWORDS, "Setting ky integration range in symmetric case"},
  {"IntegrateKxKy",                 (PyCFunction) MESH_SimulationGrating_IntegrateKxKy,                 METH_VARARGS | METH_KEYWORDS, "Action to integrate kx and ky"},
  {"IntegrateKxKyMPI",              (PyCFunction) MESH_SimulationGrating_IntegrateKxKyMPI,              METH_VARARGS | METH_KEYWORDS, "Action to integrate kx and ky using MPI"},
  {"SetLattice",                    (PyCFunction) MESH_SimulationGrating_SetLatticeGrating,             METH_VARARGS | METH_KEYWORDS, "Setting lattice constant"},
  {"SetLayerPatternGrating",        (PyCFunction) MESH_SimulationGrating_SetLayerPatternGrating,        METH_VARARGS | METH_KEYWORDS, "Setting grating pattern"},
  {"SetNumOfG",                     (PyCFunction) MESH_SimulationGrating_SetNumOfG,                     METH_VARARGS | METH_KEYWORDS, "Setting number of G"},
  {"GetNumOfG",                     (PyCFunction) MESH_SimulationGrating_GetNumOfG,                     METH_VARARGS | METH_KEYWORDS, "Getting number of G"},
  {NULL}  /* Sentinel */
};

/* TYPE ... whatever */
static PyTypeObject MESH_SimulationGrating_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "MESH_SimulationGrating",          /* tp_name */
  sizeof(MESH_SimulationGrating),    /* tp_basicsize */
  0,                         /* tp_itemsize */
  (destructor)MESH_SimulationGrating_dealloc, /* tp_dealloc */
  0,                         /* tp_print */
  0,                         /* tp_getattr */
  0,                         /* tp_setattr */
  0,                         /* tp_reserved */
  0,                         /* tp_repr */
  0,                         /* tp_as_number */
  0,                         /* tp_as_sequence */
  0,                         /* tp_as_mapping */
  0,                         /* tp_hash  */
  0,                         /* tp_call */
  0,                         /* tp_str */
  0,                         /* tp_getattro */
  0,                         /* tp_setattro */
  0,                         /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT,        /* tp_flags */ 
  "MESH_SimulationGrating",   /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,                         /* tp_richcompare */
  0,                         /* tp_weaklistoffset */
  0,                         /* tp_iter */
  0,                         /* tp_iternext */
  MESH_SimulationGrating_methods,    /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */ 
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  0,                         /* tp_init */
  0,                         /* tp_alloc */
  MESH_SimulationGrating_new, /* tp_new */
};


/*======================================================*/
// wrappaer for 2D pattern
/*=======================================================*/
/* OBJECT */
typedef struct {
  PyObject_HEAD; // <----- PUTTING THIS FIRST INHERITS THE BASE PYTHON CLASS!!!
    // Own c-variables:
    // e.g int x = 0;
    SimulationPattern *s;
} MESH_SimulationPattern;


/* DEALLOC */
static void MESH_SimulationPattern_dealloc(MESH_SimulationPattern* self){
  // Nothing to deallocate. Use this in case the base
  // struct contains a pointer or so.
  delete self->s;
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/* NEW */
static PyObject *MESH_SimulationPattern_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
  MESH_SimulationPattern* self;
  self = (MESH_SimulationPattern*)type->tp_alloc(type, 0);

   // Init vars if any ...
   if (self != NULL){
      self->s = new SimulationPattern();
   }
   return (PyObject*) self;
}

static PyObject* MESH_SimulationPattern_AddMaterial(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"material_name", (char*)"file_name", NULL };
	const char *materialName, *fileName;
	if(PyArg_ParseTupleAndKeywords(args, kwds, "ss:AddMaterial", kwlist, &materialName, &fileName)){ 
    std::string material_name(materialName), file_name(fileName);
    self->s->addMaterial(material_name, file_name);
  }
  else{
    static char *kwlist[] = { (char*)"material name", (char*)"data", NULL };
    struct interpolator_converter_data interpolator_data;
    const char *materialName;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:AddMaterial", kwlist, &materialName, &interpolator_converter, &interpolator_data)){
      return NULL;
    }

    std::string material_name(materialName);
    self->s->addMaterial(material_name, interpolator_data.omega, interpolator_data.epsilon);
  }
  Py_RETURN_NONE;

}

static PyObject* MESH_SimulationPattern_SetMaterial(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"material_name", (char*)"epsilon", NULL };
	const char *materialName;
	struct epsilon_converter_data epsdata;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO&:SetMaterial", kwlist, &materialName, &epsilon_converter, &epsdata)){ 
    return NULL; 
  }
  std::string material_name(materialName);
  double** epsilon = new double*[epsdata.num_omega];
  std::string type = "scalar";
  int eps_len = 2;
  switch(epsdata.eps_type){
    case 0:{
      // nothing
      break;
    }
    case 1:{
      eps_len = 6;
      type = "diagonal";
      break;
    }
    default:{
      eps_len = 10;
      type = "tensor";
      break;
    }
  }

  for(int i = 0; i < epsdata.num_omega;i++){
    epsilon[i] = new double[eps_len];
    for(int j = 0; j < eps_len; j++){
      epsilon[i][j] = epsdata.epsilon[eps_len*i+j];
    }
  }

  self->s->setMaterial(material_name, epsilon, type);
  for(int i = 0; i < epsdata.num_omega; i++){
    delete [] epsilon[i];
  }
  delete [] epsilon;
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPattern_AddLayer(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"thickness", (char*)"material_name", NULL };
	const char *materialName, *layerName;
  double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayer", kwlist, &layerName, &thickness, &materialName)){ 
      return NULL; 
  }
  std::string material_name(materialName), layer_name(layerName);
  self->s->addLayer(layer_name, thickness, material_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetLayer(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"thickness", (char*)"material_name", NULL };
	const char *materialName, *layerName;
  double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sds:AddLayer", kwlist, &layerName, &thickness, &materialName)){ 
      return NULL; 
  }
  std::string material_name(materialName), layer_name(layerName);
  self->s->setLayer(layer_name, thickness, material_name);
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPattern_SetLayerThickness(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"thickness", NULL };
	const char *layerName;
  double thickness;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sd:AddLayer", kwlist, &layerName, &thickness)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->setLayerThickness(layer_name, thickness);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_AddLayerCopy(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", (char*)"copy_layer_name", NULL };
	const char *layerName, *copyLayerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss:AddLayer", kwlist, &layerName, &copyLayerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName), copy_layer_name(copyLayerName);
  self->s->addLayerCopy(layer_name, copy_layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_DeleteLayer(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", NULL };
	const char *layerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:DeleteLayer", kwlist, &layerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->deleteLayer(layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetSourceLayer(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", NULL };
	const char *layerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:SetSourceLayer", kwlist, &layerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->setSourceLayer(layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetProbeLayer(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"layer_name", NULL };
	const char *layerName;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:SetProbeLayer", kwlist, &layerName)){ 
      return NULL; 
  }
  std::string layer_name(layerName);
  self->s->setProbeLayer(layer_name);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetProbeLayerZCoordinate(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"z", NULL };
	double z;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "d:SetProbeLayerZCoordinate", kwlist, &z)){ 
      return NULL; 
  }
  self->s->setProbeLayerZCoordinate(z);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetThread(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = { (char*)"num_thread", NULL };
	int thread;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "i:SetThread", kwlist, &thread)){ 
      return NULL; 
  }
  self->s->setThread(thread);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_InitSimulation(MESH_SimulationPattern *self, PyObject *args){
  self->s->initSimulation();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_GetPhi(MESH_SimulationPattern *self, PyObject *args){
  double* phi = self->s->getPhi();
  int num_omega = self->s->getNumOfOmega();
  PyObject* phi_value = PyTuple_New(num_omega);
  for(int i = 0; i < num_omega; i++){
    PyTuple_SetItem(phi_value, i, PyFloat_FromDouble(phi[i]));
  }
  return phi_value;
}

static PyObject* MESH_SimulationPattern_GetOmega(MESH_SimulationPattern *self, PyObject *args){
  double* omega = self->s->getOmega();
  int num_omega = self->s->getNumOfOmega();
  PyObject* omega_value = PyTuple_New(num_omega);
  for(int i = 0; i < num_omega; i++){
    PyTuple_SetItem(omega_value, i, PyFloat_FromDouble(omega[i]));
  }
  return omega_value;
}

static PyObject* MESH_SimulationPattern_GetEpsilon(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"coordinate", NULL};
  int omega_index;
  struct coordinate_converter_data coordata;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "iO&:GetEpsilon", kwlist, &omega_index, &coordinate_converter, &coordata)){ 
    return NULL; 
  }
  double* epsilon = new double[10];
  self->s->getEpsilon(omega_index, coordata.coor, epsilon);
  
  PyObject* eps_value = PyTuple_New(10);
  for(int i = 0; i < 10; i++){
    PyTuple_SetItem(eps_value, i, PyFloat_FromDouble(epsilon[i]));
  }
  delete [] epsilon;
  return eps_value;
}

static PyObject* MESH_SimulationPattern_OutputLayerPatternRealization(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"layer_name", (char*)"Nu", (char*)"Nv", (char*)"file_name", NULL};
  int omega_index, Nu, Nv;
  char* layerName, *fileName = (char*)"";
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "isii|s:OutputLayerPatternRealization", kwlist, &omega_index, &layerName, &Nu, &Nv, &fileName)){ 
    return NULL; 
  }
  std::string layer_name(layerName), file_name(fileName);
  if(file_name == ""){
    self->s->outputLayerPatternRealization(omega_index, layer_name, Nu, Nv);
  }
  else{
    self->s->outputLayerPatternRealization(omega_index, layer_name, Nu, Nv, file_name);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_GetNumOfOmega(MESH_SimulationPattern *self, PyObject *args){
  int num_omega = self->s->getNumOfOmega();
  return FromIntPyDefInt(num_omega);
}

static PyObject* MESH_SimulationPattern_GetPhiAtKxKy(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"omega_index", (char*)"kx", (char*)"ky", NULL};
  int omega_index;
  double kx, ky = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "id|d:GetPhiAtKxKy", kwlist, &omega_index, &kx, &ky)){ 
    return NULL; 
  }
  double value = self->s->getPhiAtKxKy(omega_index, kx, ky);
  return PyFloat_FromDouble(value);
}

static PyObject* MESH_SimulationPattern_OutputSysInfo(MESH_SimulationPattern *self, PyObject *args){
  self->s->outputSysInfo();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_OptPrintIntermediate(MESH_SimulationPattern *self, PyObject *args){
  self->s->optPrintIntermediate();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_OptOnlyComputeTE(MESH_SimulationPattern *self, PyObject *args){
  self->s->optOnlyComputeTE();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_OptOnlyComputeTM(MESH_SimulationPattern *self, PyObject *args){
  self->s->optOnlyComputeTM();
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetKxIntegral(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKxIntegral", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKxIntegral(num_points);
  }
  else{
    self->s->setKxIntegral(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetKyIntegral(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKyIntegral", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKyIntegral(num_points);
  }
  else{
    self->s->setKyIntegral(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetKxIntegralSym(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKxIntegralSym", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKxIntegralSym(num_points);
  }
  else{
    self->s->setKxIntegralSym(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetKyIntegralSym(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_points", (char*)"integral_end", NULL};
  int num_points;
  double integral_end = 0;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|d:SetKyIntegralSym", kwlist, &num_points, &integral_end)){ 
    return NULL; 
  }
  if(integral_end == 0){
    self->s->setKyIntegralSym(num_points);
  }
  else{
    self->s->setKyIntegralSym(num_points, integral_end);
  }
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_IntegrateKxKy(MESH_SimulationPattern *self, PyObject *args){
  self->s->integrateKxKy();
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPattern_IntegrateKxKyMPI(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"rank", (char*)"size", NULL};
  int rank, size;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ii:IntegrateKxKyMPI", kwlist, &rank, &size)){ 
    return NULL; 
  }
  self->s->integrateKxKyMPI(rank, size);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetNumOfG(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"num_G",  NULL};
  int numG;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i:SetNumOfG", kwlist, &numG)){ 
    return NULL; 
  }
  self->s->setNumOfG(numG);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_GetNumOfG(MESH_SimulationPattern *self, PyObject *args){
  int numG = self->s->getNumOfG();
  return FromIntPyDefInt(numG);
}

static PyObject* MESH_SimulationPattern_OptSetLatticeTruncation(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"truncation", NULL};
  char* truncation;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "s:OptSetLatticeTruncation", kwlist, &truncation)){ 
    return NULL; 
  }
  std::string trunc(truncation);
  self->s->optSetLatticeTruncation(trunc);
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPattern_SetLatticePattern(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"lattice1_len", (char*)"lattice2_len", (char*)"angle", NULL};
  double lattice1_len, lattice2_len, angle;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ddd:SetLattice", kwlist, &lattice1_len, &lattice2_len, &angle)){ 
    return NULL; 
  }
  self->s->setLattice(lattice1_len, lattice2_len, angle);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_GetReciprocalLattice(MESH_SimulationPattern *self, PyObject *args){
  double lattice[4];
  self->s->getReciprocalLattice(lattice);

  return PyTuple_Pack(2,
		PyTuple_Pack(2,
			PyFloat_FromDouble(lattice[0]), PyFloat_FromDouble(lattice[1])
		),
		PyTuple_Pack(2,
			PyFloat_FromDouble(lattice[2]), PyFloat_FromDouble(lattice[3])
		)
	);
}

static PyObject* MESH_SimulationPattern_SetLayerPatternRectangle(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"layer_name", (char*)"material_name", (char*)"center", (char*)"angle", (char*)"width", NULL};
  double center[2], width[2], angle;
  char* materialName, *layerName; 
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)d(dd):SetLayerPatternRectangle", kwlist, &layerName, &materialName, &center[0], &center[1], &angle, &width[0], &width[1])){ 
    return NULL; 
  }
  std::string layer_name(layerName), material_name(materialName);
  self->s->setLayerPatternRectangle(layer_name, material_name, center[0], center[1], angle, width[0], width[1]);
  Py_RETURN_NONE;
}


static PyObject* MESH_SimulationPattern_SetLayerPatternEllipse(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"layer_name", (char*)"material_name", (char*)"center", (char*)"angle", (char*)"half_width", NULL};
  double center[2], width[2], angle;
  char* materialName, *layerName; 
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)d(dd):SetLayerPatternEllipse", kwlist, &layerName, &materialName, &center[0], &center[1], &angle, &width[0], &width[1])){ 
    return NULL; 
  }
  std::string layer_name(layerName), material_name(materialName);
  self->s->setLayerPatternEllipse(layer_name, material_name, center[0], center[1], angle, width[0], width[1]);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetLayerPatternCircle(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"layer_name", (char*)"material_name", (char*)"center", (char*)"radius", NULL};
  double center[2], radius;
  char* materialName, *layerName; 
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)d:SetLayerPatternCircle", kwlist, &layerName, &materialName, &center[0], &center[1], &radius)){ 
    return NULL; 
  }
  std::string layer_name(layerName), material_name(materialName);
  self->s->setLayerPatternCircle(layer_name, material_name, center[0], center[1], radius);
  Py_RETURN_NONE;
}

static PyObject* MESH_SimulationPattern_SetLayerPatternPolygon(MESH_SimulationPattern *self, PyObject *args, PyObject *kwds){
  static char* kwlist[] = {(char*)"layer_name", (char*)"material_name", (char*)"center", (char*)"angle", (char*)"vertices", NULL};
  double center[2], angle;
  char* materialName, *layerName; 
  struct polygon_converter_data polydata;
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ss(dd)dO&:SetLayerPatternPolygon", kwlist, &layerName, &materialName, &center[0], &center[1], &angle, &polygon_converter, &polydata)){ 
    return NULL; 
  }
  std::string layer_name(layerName), material_name(materialName);
  int numOfPoints = polydata.nvert;
  double** edgePoints = new double*[numOfPoints];
  for(int i = 0; i < numOfPoints; i++){
    edgePoints[i] = new double[2];
    for(int j = 0; j < 2; j++){
      edgePoints[i][j] = polydata.vert[2*i+j];
    }
  }
  self->s->setLayerPatternPolygon(layer_name, material_name, center[0], center[1], angle, edgePoints, numOfPoints);
  for(int i = 0; i < numOfPoints; i++){
    delete [] edgePoints[i];
  }
  delete [] edgePoints;
  Py_RETURN_NONE;
}
/* METHOD TABLE */
static PyMethodDef MESH_SimulationPattern_methods[] = {
  {"AddMaterial",                   (PyCFunction) MESH_SimulationPattern_AddMaterial,                   METH_VARARGS | METH_KEYWORDS, "Adding new material to the simulation"},
  {"SetMaterial",                   (PyCFunction) MESH_SimulationPattern_SetMaterial,                   METH_VARARGS | METH_KEYWORDS, "Setting material property"},
  {"AddLayer",                      (PyCFunction) MESH_SimulationPattern_AddLayer,                      METH_VARARGS | METH_KEYWORDS, "Adding new layer to the simulation"},
  {"SetLayer",                      (PyCFunction) MESH_SimulationPattern_SetLayer,                      METH_VARARGS | METH_KEYWORDS, "Setting layer property"},
  {"SetLayerThickness",             (PyCFunction) MESH_SimulationPattern_SetLayerThickness,             METH_VARARGS | METH_KEYWORDS, "Setting layer thickness"},
  {"AddLayerCopy",                  (PyCFunction) MESH_SimulationPattern_AddLayerCopy,                  METH_VARARGS | METH_KEYWORDS, "Making a copy of existing layer"},
  {"DeleteLayer",                   (PyCFunction) MESH_SimulationPattern_DeleteLayer,                   METH_VARARGS | METH_KEYWORDS, "Deleting an existing layer"},
  {"SetSourceLayer",                (PyCFunction) MESH_SimulationPattern_SetSourceLayer,                METH_VARARGS | METH_KEYWORDS, "Setting a source layer"},
  {"SetProbeLayer",                 (PyCFunction) MESH_SimulationPattern_SetProbeLayer,                 METH_VARARGS | METH_KEYWORDS, "Setting the probe layer"},
  {"SetProbeLayerZCoordinate",      (PyCFunction) MESH_SimulationPattern_SetProbeLayerZCoordinate,      METH_VARARGS | METH_KEYWORDS, "Setting the z-coordinate in the probe layer"},
  {"SetThread",                     (PyCFunction) MESH_SimulationPattern_SetThread,                     METH_VARARGS | METH_KEYWORDS, "Setting the number of thread"},
  {"InitSimulation",                (PyCFunction) MESH_SimulationPattern_InitSimulation,                METH_VARARGS | METH_KEYWORDS, "Initializing simulation"},
  {"GetPhi",                        (PyCFunction) MESH_SimulationPattern_GetPhi,                        METH_VARARGS | METH_KEYWORDS, "Getting the value of Phi"},
  {"GetOmega",                      (PyCFunction) MESH_SimulationPattern_GetOmega,                      METH_VARARGS | METH_KEYWORDS, "Getting all the omega values"},
  {"GetEpsilon",                    (PyCFunction) MESH_SimulationPattern_GetEpsilon,                    METH_VARARGS | METH_KEYWORDS, "Getting epsilon at one frequency"},
  {"OutputLayerPatternRealization", (PyCFunction) MESH_SimulationPattern_OutputLayerPatternRealization, METH_VARARGS | METH_KEYWORDS, "Outputting dielectric reconstruction"},
  {"GetNumOfOmega",                 (PyCFunction) MESH_SimulationPattern_GetNumOfOmega,                 METH_VARARGS | METH_KEYWORDS, "Getting the number of omega"},
  {"GetPhiAtKxKy",                  (PyCFunction) MESH_SimulationPattern_GetPhiAtKxKy,                  METH_VARARGS | METH_KEYWORDS, "Getting Phi value at a (kx,ky) pair"},
  {"OutputSysInfo",                 (PyCFunction) MESH_SimulationPattern_OutputSysInfo,                 METH_VARARGS | METH_KEYWORDS, "Outputting system information"},
  {"OptPrintIntermediate",          (PyCFunction) MESH_SimulationPattern_OptPrintIntermediate,          METH_VARARGS | METH_KEYWORDS, "Option to output intermediate results"},
  {"OptOnlyComputeTE",              (PyCFunction) MESH_SimulationPattern_OptOnlyComputeTE,              METH_VARARGS | METH_KEYWORDS, "Option to only compute TE mode"},
  {"OptOnlyComputeTM",              (PyCFunction) MESH_SimulationPattern_OptOnlyComputeTM,              METH_VARARGS | METH_KEYWORDS, "Option to only compute TM mode"},
  {"SetKxIntegral",                 (PyCFunction) MESH_SimulationPattern_SetKxIntegral,                 METH_VARARGS | METH_KEYWORDS, "Setting kx integration range"},
  {"SetKyIntegral",                 (PyCFunction) MESH_SimulationPattern_SetKyIntegral,                 METH_VARARGS | METH_KEYWORDS, "Setting kx integration range"},
  {"SetKxIntegralSym",              (PyCFunction) MESH_SimulationPattern_SetKxIntegralSym,              METH_VARARGS | METH_KEYWORDS, "Setting kx integration range in symmetric case"},
  {"SetKyIntegralSym",              (PyCFunction) MESH_SimulationPattern_SetKyIntegralSym,              METH_VARARGS | METH_KEYWORDS, "Setting ky integration range in symmetric case"},
  {"IntegrateKxKy",                 (PyCFunction) MESH_SimulationPattern_IntegrateKxKy,                 METH_VARARGS | METH_KEYWORDS, "Action to integrate kx and ky"},
  {"IntegrateKxKyMPI",              (PyCFunction) MESH_SimulationPattern_IntegrateKxKyMPI,              METH_VARARGS | METH_KEYWORDS, "Action to integrate kx and ky using MPI"},
  {"SetLattice",                    (PyCFunction) MESH_SimulationPattern_SetLatticePattern,             METH_VARARGS | METH_KEYWORDS, "Setting lattice constants"},
  {"GetReciprocalLattice",          (PyCFunction) MESH_SimulationPattern_GetReciprocalLattice,          METH_VARARGS | METH_KEYWORDS, "Getting reciprocal lattice"},
  {"SetLayerPatternRectangle",      (PyCFunction) MESH_SimulationPattern_SetLayerPatternRectangle,      METH_VARARGS | METH_KEYWORDS, "Setting rectangle layer pattern"},
  {"SetLayerPatternEllipse",        (PyCFunction) MESH_SimulationPattern_SetLayerPatternEllipse,        METH_VARARGS | METH_KEYWORDS, "Setting ellipse layer pattern"},
  {"SetLayerPatternCircle",         (PyCFunction) MESH_SimulationPattern_SetLayerPatternCircle,         METH_VARARGS | METH_KEYWORDS, "Setting circle layer pattern"},
  {"SetLayerPatternPolygon",        (PyCFunction) MESH_SimulationPattern_SetLayerPatternPolygon,        METH_VARARGS | METH_KEYWORDS, "Setting polygon layer pattern"},
  {"SetNumOfG",                     (PyCFunction) MESH_SimulationPattern_SetNumOfG,                     METH_VARARGS | METH_KEYWORDS, "Setting number of G"},
  {"GetNumOfG",                     (PyCFunction) MESH_SimulationPattern_GetNumOfG,                     METH_VARARGS | METH_KEYWORDS, "Getting number of G"},
  {"OptSetLatticeTruncation",       (PyCFunction) MESH_SimulationPattern_OptSetLatticeTruncation,       METH_VARARGS | METH_KEYWORDS, "Setting G selection method"},
  {NULL}  /* Sentinel */
};

/* TYPE ... whatever */
static PyTypeObject MESH_SimulationPattern_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "MESH_SimulationPattern",          /* tp_name */
  sizeof(MESH_SimulationPattern),    /* tp_basicsize */
  0,                         /* tp_itemsize */
  (destructor)MESH_SimulationPattern_dealloc, /* tp_dealloc */
  0,                         /* tp_print */
  0,                         /* tp_getattr */
  0,                         /* tp_setattr */
  0,                         /* tp_reserved */
  0,                         /* tp_repr */
  0,                         /* tp_as_number */
  0,                         /* tp_as_sequence */
  0,                         /* tp_as_mapping */
  0,                         /* tp_hash  */
  0,                         /* tp_call */
  0,                         /* tp_str */
  0,                         /* tp_getattro */
  0,                         /* tp_setattro */
  0,                         /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT,         /* tp_flags */  
  "MESH_SimulationPattern",   /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,                         /* tp_richcompare */
  0,                         /* tp_weaklistoffset */
  0,                         /* tp_iter */
  0,                         /* tp_iternext */
  MESH_SimulationPattern_methods,    /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */ 
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  0,                         /* tp_init */
  0,                         /* tp_alloc */
  MESH_SimulationPattern_new, /* tp_new */
};


/*======================================================*/
// special function to return the values of constants
/*=======================================================*/
static PyObject* MESH_Constants(PyObject *self, PyObject *args){
  PyObject* consts = PyDict_New();
  std::vector<const char*> properties = {"pi", "k_B", "eps_0", "m_e", "eV", "mu_0", "h", "h_bar", "c_0", "q", "sigma"};
  std::vector<double> vals = {constants.pi, constants.k_B, constants.eps_0, constants.m_e, constants.eV, constants.mu_0,
    constants.h, constants.h_bar, constants.c_0, constants.q, constants.sigma};
  int tableLen = properties.size();
  for(int i = 0; i < tableLen; i++){
    PyDict_SetItem(consts, PyString_FromString(properties[i]), PyFloat_FromDouble(vals[i]));
  }
  return consts;
}

/* MODULE codie function table */
static PyMethodDef MESH_module_methods[] = {
  {"Constants", MESH_Constants, METH_VARARGS | METH_KEYWORDS,"Physical Constants"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

/*======================================================*/
// registering the classes
/*=======================================================*/
#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef MESH_module = {
    PyModuleDef_HEAD_INIT,
    "MESH",
    module_doc,
    sizeof(struct module_state),
    MESH_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
  };
  #define INITERROR return NULL
  PyObject * PyInit_MESH(void)
#else
  #define INITERROR return
  PyMODINIT_FUNC initMESH(void)
#endif
{   // The name after init MUST be the same name as in setup.py
  /* Initialization function for the module (*must* be called PyInit_FunctionSampler1D) */
  PyDoc_STRVAR(module_doc, "MESH: Multilayer Electromagnetic Solver for Heat transfer");

  PyObject* module;
  // Create the class types
  if (PyType_Ready(&MESH_Interpolator_Type) < 0)
    INITERROR;
  if (PyType_Ready(&MESH_SimulationPlanar_Type) < 0)
    INITERROR;
  if (PyType_Ready(&MESH_SimulationGrating_Type) < 0)
    INITERROR;
  if (PyType_Ready(&MESH_SimulationPattern_Type) < 0)
    INITERROR; 


  // Init Module
  #if PY_MAJOR_VERSION >= 3
    module = PyModule_Create(&MESH_module);
  #else
    module = Py_InitModule3("MESH", MESH_module_methods, // The name must be the same as behind the init and in the setup script
                      module_doc);
  #endif

  if (module == NULL)
    // PY3 return NULL;
    INITERROR;

  // Initialize Classes
  // Py_INCREF(&MESH_Simulation_Type);
  // PyModule_AddObject(m, "Simulation", (PyObject *)&MESH_Simulation_Type); // This name is used for the import command!!!!
  Py_INCREF(&MESH_Interpolator_Type);
  PyModule_AddObject(module, "Interpolator", (PyObject *)&MESH_Interpolator_Type);

  Py_INCREF(&MESH_SimulationPlanar_Type);
  PyModule_AddObject(module, "SimulationPlanar", (PyObject *)&MESH_SimulationPlanar_Type);

  Py_INCREF(&MESH_SimulationGrating_Type);
  PyModule_AddObject(module, "SimulationGrating", (PyObject *)&MESH_SimulationGrating_Type);

  Py_INCREF(&MESH_SimulationPattern_Type);
  PyModule_AddObject(module, "SimulationPattern", (PyObject *)&MESH_SimulationPattern_Type);

  #if PY_MAJOR_VERSION >= 3
	  return module;
  #endif
}

#ifdef __cplusplus
}
#endif
