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
#ifndef _SYSTEM_H
#define _SYSTEM_H
#include <vector>
#include <map>
#include "Common.h"
#include <string>

namespace SYSTEM{
  /*======================================================
  Implementaion of the Material class
  =======================================================*/
  class Material{
  public:
    Material(std::string name, dcomplex* epsilonList, double* omegaList, size_t size);
    Material(std::string name);
    Material(const Material& material);
    ~Material();

    std::string getName();
    dcomplex* getEpsilon();
    double* getOmega();

    void setName(const std::string name);
    void setOmega(const double* omega, size_t size);
    void setEpsilon(const dcomplex* epsilon, size_t size);

  private:
    std::string name_;
    dcomplex* epsilonList_;
    double* omegaList_;
    size_t size_;
  };

  typedef std::vector<Material*> MaterialVec;
  typedef MaterialVec::iterator MaterialIter;
  typedef MaterialVec::const_iterator const_MaterialIter;

  /*======================================================
  Implementaion of the Layer class
  =======================================================*/
  typedef std::vector< std::pair<double, double> > LayerPattern;
  typedef LayerPattern::iterator PatternIter;
  typedef LayerPattern::const_iterator const_PatternIter;

  class Layer{
  public:
    Layer(Material* material, double thickness, SOURCE source = ISNOTSOURCE_);
    Layer(Material* material);
    Layer();
    ~Layer();
    Layer(const Layer& layer);

    void setBackGround(Material* material);
    void setThickness(double thickness);
    void isSource();
    void isNotSource();

    Material* getBackGround();
    Material* getMaterialByName(std::string name);
    int getNumOfMaterial();
    double getThickness();

    const_MaterialIter getVecBegin();
    const_MaterialIter getVecEnd();

    const_PatternIter getArg1Begin();
    const_PatternIter getArg2Begin();

    const_PatternIter getArg1End();
    const_PatternIter getArg2End();


    void addPattern(Material* material, double args1[2], double args2[2], std::string pattern);
    void addPattern(Material* material, double args1[2]);

  private:

    double thickness_;
    Material* backGround_;
    MaterialVec materialVec_;
    std::string pattern_;
    LayerPattern args1_;
    LayerPattern args2_;
    SOURCE source_;
  };

  typedef std::map<int, Layer*> LayerMap;
  typedef LayerMap::iterator LayerIter;
  typedef LayerMap::const_iterator const_LayerIter;

  /*======================================================
  Implementaion of the structure class
  =======================================================*/
  class Structure{
  public:
    Structure();
    ~Structure();

    Structure(Structure& structure);

    void addLayer(Layer* layer);
    Layer* getLayerByIndex(int index);
    int getNumOfLayer();
    double* getThicknessList();

    const_LayerIter getMapBegin();
    const_LayerIter getMapEnd();

private:
    void setGx(int nGx);
    void setGy(int nGy);

    int getGx();
    int getGy();

    LayerMap layerMap_;
    int nGx_;
    int nGy_;
  };




}
#endif