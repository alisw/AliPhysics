// $Id$
//
// Class AliMUONSlatGeometryBuilder
// -----------------------------
// MUON Station3 geometry construction class.
//
// Author: Eric Dumonteil, Philippe Crochet

#ifndef ALI_MUON_SLAT_GEOMETRY_BUILDER_H
#define ALI_MUON_SLAT_GEOMETRY_BUILDER_H

#include "AliMUONVGeometryBuilder.h"

class AliMUON;

class AliMUONSlatGeometryBuilder : public AliMUONVGeometryBuilder
{
  public:
    AliMUONSlatGeometryBuilder(AliMUON* muon);
    AliMUONSlatGeometryBuilder(const AliMUONSlatGeometryBuilder& rhs);
    AliMUONSlatGeometryBuilder();
    virtual ~AliMUONSlatGeometryBuilder();

    // operators  
    AliMUONSlatGeometryBuilder& operator = (const AliMUONSlatGeometryBuilder& rhs);
  
    // methods
    virtual void CreateGeometry();
    virtual void SetTransformations();
    virtual void SetSensitiveVolumes();
    
  private:
    Int_t  ConvertSlatNum(Int_t numslat, Int_t quadnum, Int_t fspq) const;

    AliMUON*  fMUON; // the MUON detector class 
        
  ClassDef(AliMUONSlatGeometryBuilder,1) // MUON chamber geometry base class
};

#endif //ALI_MUON_SLAT_GEOMETRY_BUILDER_H
