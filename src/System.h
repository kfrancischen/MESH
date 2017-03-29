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
  class Material : public NamedInterface {
  public:
    static Ptr<Material> instanceNew(
      const std::string name,
      const double* omegaList,
      const EPSILON& epsilonList,
      const int numOfOmega
    );

    Material(const Material& material) = delete;
    ~Material();

    std::string getName();
    EPSTYPE getType();
    //dcomplex* getEpsilonList();
    EpsilonVal getEpsilonAtIndex(const int index);
    int getNumOfOmega();
    double* getOmegaList();

    void setOmega(const double* omegaList, const int numOfOmega);
    void setEpsilon(const EPSILON& epsilonList, const int numOfOmega);

  protected:
    Material(
      const std::string name,
      const double* omegaList,
      const EPSILON& epsilonList,
      const int numOfOmega
    );

    EPSILON epsilonList_;
    double* omegaList_;
    int numOfOmega_;
  };

  typedef std::vector< Ptr<Material> > MaterialVec;
  typedef MaterialVec::iterator MaterialIter;
  typedef MaterialVec::const_iterator const_MaterialIter;

  /*======================================================
  Implementaion of the Layer class
  =======================================================*/
  typedef std::pair<double, double> LayerPattern;
  typedef struct PATTERNWRAPPER{
    LayerPattern arg1_;
    LayerPattern arg2_;
    PATTERN type_;
  } Pattern;
  typedef std::vector< Pattern > PatternVec;
  typedef PatternVec::iterator PatternIter;
  typedef PatternVec::const_iterator const_PatternIter;

  class Layer : public NamedInterface{
  public:
    static Ptr<Layer> instanceNew(
      const string name,
      const Ptr<Material>& material,
      const double thickness
    );

    static Ptr<Layer> instanceNew(
      const string name
    );

    ~Layer();
    Layer(const Layer& layer) = delete;

    Ptr<Layer> layerCopy(const string name);

    void setBackGround(const Ptr<Material>& material);
    void setThickness(const double thickness);
    void setIsSource();
    bool checkIsSource();
    void containTensor(bool val);
    bool hasTensor();
    bool hasMaterial(const Ptr<Material>& material);

    Ptr<Material> getBackGround();
    Ptr<Material> getMaterialByName(const std::string name);
    int getNumOfMaterial();
    double getThickness();
    std::string getName();

    const_MaterialIter getMaterialsBegin();
    const_MaterialIter getMaterialsEnd();

    const_PatternIter getPatternsBegin();
    const_PatternIter getPatternsEnd();


    void addRectanlgePattern(const Ptr<Material>& material, const double args1[2], const double args2[2]);
    void addCirclePattern(const Ptr<Material>& material, const double args[2], const double radius);
    void addGratingPattern(const Ptr<Material>& material, const double center, const double width);

  private:
    enum SOURCE {ISSOURCE_, ISNOTSOURCE_};

    Layer(const string name, const Ptr<Material>& material, const double thickness);
    Layer(const string name);

    double thickness_;
    Ptr<Material> backGround_;
    MaterialVec materialVec_;
    PatternVec patternVec_;
    SOURCE source_;
    int hasTensor_ = 0;
  };

  typedef std::map<int, Ptr<Layer> > LayerMap;
  typedef LayerMap::iterator LayerIter;
  typedef LayerMap::const_iterator const_LayerIter;

  /*======================================================
  Implementaion of the structure class
  =======================================================*/
  class Structure : public PtrInterface{
  public:

    static Ptr<Structure> instanceNew();
    ~Structure();

    Structure(const Structure& structure);

    void addLayer(const Ptr<Layer>& layer);
    void deleteLayerByName(const string name);
    void deleteLayerByLayer(const Ptr<Layer>& layer);
    Ptr<Layer> getLayerByIndex(const int index);
    Ptr<Layer> getLayerByName(const std::string name);
    int getNumOfLayer();
    void getThicknessList(double* thicknessList);

    const_LayerIter getLayersBegin();
    const_LayerIter getLayersEnd();

  private:
    Structure();
    void deleteLayer(const_LayerIter it);
    void reorganizeLayers();
    LayerMap layerMap_;
  };


}
#endif