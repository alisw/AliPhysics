/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONSt1GeometryBuilder
// -----------------------------
// MUON Station1 coarse geometry construction class.
//
// Extracted from AliMUONv1
// by Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_ST1_GEOMETRY_BUILDER_H
#define ALI_MUON_ST1_GEOMETRY_BUILDER_H

#include "AliMUONVGeometryBuilder.h"

class AliMUON;

class AliMUONSt1GeometryBuilder : public AliMUONVGeometryBuilder
{
  public:
    AliMUONSt1GeometryBuilder(AliMUON* muon);
    AliMUONSt1GeometryBuilder(const AliMUONSt1GeometryBuilder& rhs);
    AliMUONSt1GeometryBuilder();
    virtual ~AliMUONSt1GeometryBuilder();

    // operators  
    AliMUONSt1GeometryBuilder& operator = (const AliMUONSt1GeometryBuilder& rhs);
  
    // methods
    virtual void CreateGeometry();
    virtual void SetTransformations();
    virtual void SetSensitiveVolumes();
    
  private:
     AliMUON*  fMUON; // the MUON detector class 
        
  ClassDef(AliMUONSt1GeometryBuilder,1) // MUON chamber geometry base class
};

#endif //ALI_MUON_ST1_GEOMETRY_BUILDER_H
