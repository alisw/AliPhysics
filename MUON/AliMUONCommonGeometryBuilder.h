/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONCommonGeometryBuilder
// ----------------------------------
// Geometry construction common to all stations
// (material definition).
// separated from AliMUONGeometryBuilder

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
    virtual void  CreateGeometry()      {}
    virtual void  SetSensitiveVolumes() {}
    virtual void  SetTransformations()  {}

  protected:
    AliMUONCommonGeometryBuilder(const AliMUONCommonGeometryBuilder& right);
    AliMUONCommonGeometryBuilder&  
                     operator = (const AliMUONCommonGeometryBuilder& right);
 
  private:
    // data members
    AliMUON*  fMUON; // MUON detector

  ClassDef(AliMUONCommonGeometryBuilder,1)  // Common MUON geometry definitions
};

#endif //ALI_MUON_COMMON_GEOMETRY_BUILDER_H







