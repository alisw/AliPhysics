/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004
//
// Class AliMUONSlatGeometryBuilder
// -----------------------------
// MUON Station3 geometry construction class.
//
// Author: Eric Dumonteil (dumontei@cea.fr)


#ifndef ALI_MUON_SLAT_GEOMETRY_BUILDER_H
#define ALI_MUON_SLAT_GEOMETRY_BUILDER_H

#include "AliMUONVGeometryBuilder.h"

class AliMUON;

class AliMUONSlatGeometryBuilder : public AliMUONVGeometryBuilder
{
  public:
    AliMUONSlatGeometryBuilder(AliMUON* muon);
    AliMUONSlatGeometryBuilder();
    virtual ~AliMUONSlatGeometryBuilder();

    // methods
    virtual void CreateGeometry();
    virtual void SetTransformations();
    virtual void SetSensitiveVolumes();

  protected:
    AliMUONSlatGeometryBuilder(const AliMUONSlatGeometryBuilder& rhs);

    // operators  
    AliMUONSlatGeometryBuilder& operator = (const AliMUONSlatGeometryBuilder& rhs);
    
  private:
    Int_t  ConvertSlatNum(Int_t numslat, Int_t quadnum, Int_t fspq) const;

    AliMUON*  fMUON; // the MUON detector class 
        
  ClassDef(AliMUONSlatGeometryBuilder,1) // MUON chamber geometry base class
};

#endif //ALI_MUON_SLAT_GEOMETRY_BUILDER_H
