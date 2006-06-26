/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONCommonGeometryBuilder
/// \brief Class to build common materials
///
/// Geometry construction common to all stations
/// (material definition).
/// separated from AliMUONGeometryBuilder

#ifndef ALI_MUON_COMMON_GEOMETRY_BUILDER_H
#define ALI_MUON_COMMON_GEOMETRY_BUILDER_H

#include "AliMUONVGeometryBuilder.h"

class AliMUON;

class AliMUONCommonGeometryBuilder : public AliMUONVGeometryBuilder 
{
  public:
    AliMUONCommonGeometryBuilder(AliMUON* muon);
    AliMUONCommonGeometryBuilder();
    virtual  ~AliMUONCommonGeometryBuilder();

    virtual void  CreateMaterials();
    virtual void  CreateGeometry()      {}  ///< No geometry is created
    virtual void  SetSensitiveVolumes() {}  ///< No sensitive volumes are set
    virtual void  SetTransformations()  {}  ///< No transformations are set

  protected:
 
  private:
    AliMUONCommonGeometryBuilder(const AliMUONCommonGeometryBuilder& right);
    AliMUONCommonGeometryBuilder&  
                     operator = (const AliMUONCommonGeometryBuilder& right);
    // data members
    AliMUON*  fMUON; ///< the MUON detector class 

  ClassDef(AliMUONCommonGeometryBuilder,1)  // Class to build common materials 
};

#endif //ALI_MUON_COMMON_GEOMETRY_BUILDER_H







