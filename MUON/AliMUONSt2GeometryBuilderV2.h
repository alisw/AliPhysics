/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004
//
// Class AliMUONSt2GeometryBuilderV2
// -----------------------------
// MUON Station2 geometry construction class.
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_ST2_GEOMETRY_BUILDER_V2_H
#define ALI_MUON_ST2_GEOMETRY_BUILDER_V2_H

#include "AliMUONVGeometryBuilder.h"

class AliMUON;

class AliMUONSt2GeometryBuilderV2 : public AliMUONVGeometryBuilder
{
  public:
    AliMUONSt2GeometryBuilderV2(AliMUON* muon);
    AliMUONSt2GeometryBuilderV2();
    virtual ~AliMUONSt2GeometryBuilderV2();
  
    // methods
    virtual void CreateGeometry();
    virtual void SetTransformations();
    virtual void SetSensitiveVolumes();
    
  protected:
    AliMUONSt2GeometryBuilderV2(const AliMUONSt2GeometryBuilderV2& rhs);

    // operators  
    AliMUONSt2GeometryBuilderV2& operator = (const AliMUONSt2GeometryBuilderV2& rhs);
    
  private:
     AliMUON*  fMUON; // the MUON detector class 
        
  ClassDef(AliMUONSt2GeometryBuilderV2,1) // MUON chamber geometry base class
};

#endif //ALI_MUON_ST2_GEOMETRY_BUILDER_H_V2
