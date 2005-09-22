/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONCommonGeometryBuilder
/// \brief MUON Station1 coarse geometry construction class
///
/// Extracted from AliMUONv1
/// by Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_ST1_GEOMETRY_BUILDER_H
#define ALI_MUON_ST1_GEOMETRY_BUILDER_H

#include "AliMUONVGeometryBuilder.h"

class AliMUON;

class AliMUONSt1GeometryBuilder : public AliMUONVGeometryBuilder
{
  public:
    AliMUONSt1GeometryBuilder(AliMUON* muon);
    AliMUONSt1GeometryBuilder();
    virtual ~AliMUONSt1GeometryBuilder();
  
    // methods
    virtual void CreateGeometry();
    virtual void SetTransformations();
    virtual void SetSensitiveVolumes();
    
  protected:
    AliMUONSt1GeometryBuilder(const AliMUONSt1GeometryBuilder& rhs);

    // operators  
    AliMUONSt1GeometryBuilder& operator = (const AliMUONSt1GeometryBuilder& rhs);
    
  private:
     AliMUON*  fMUON; // the MUON detector class 
        
  ClassDef(AliMUONSt1GeometryBuilder,1) //MUON Station1 coarse geometry construction class 
};

#endif //ALI_MUON_ST1_GEOMETRY_BUILDER_H
