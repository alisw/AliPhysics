/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONSt2GeometryBuilderV2
/// \brief MUON Station2 detailed geometry construction class
///
/// Authors: SANJOY PAL ,Prof. SUKALYAN CHATTOPADHAYAY  [SINP, KOLKATA]
///         &  Dr.SHAKEEL AHMAD (AMU), INDIA

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
     AliMUON*  fMUON; ///< the MUON detector class 
        
  ClassDef(AliMUONSt2GeometryBuilderV2,1) //MUON Station2 detailed geometry construction class
};

#endif //ALI_MUON_ST2_GEOMETRY_BUILDER_V2_H
