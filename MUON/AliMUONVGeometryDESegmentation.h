/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Class AliVMUONGeometryDESegmentation
// ----------------------------------
// Extension for AliSegmentation interface,
// added functions:
//  Bool_t  HasPad(Float_t x, Float_t y, Float_t z);
//  Bool_t  HasPad(Int_t ix, Int_t iy);
//
// Author:Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_V_GEOMETRY_DE_SEGMENTATION_H
#define ALI_MUON_V_GEOMETRY_DE_SEGMENTATION_H

#include "AliSegmentation.h"

class AliMUONVGeometryDESegmentation : public AliSegmentation
{
  public:
    AliMUONVGeometryDESegmentation();
    virtual ~AliMUONVGeometryDESegmentation();

    // methods
    virtual Bool_t  HasPad(Float_t x, Float_t y, Float_t z) = 0; 
                       // Returns true if a pad exists in the given position
    virtual Bool_t  HasPad(Int_t ix, Int_t iy) = 0;
                       // Returns true if a pad with given indices exists

  protected:
    AliMUONVGeometryDESegmentation(const AliMUONVGeometryDESegmentation& rhs);
  
    // operators
    AliMUONVGeometryDESegmentation& operator=(
      const AliMUONVGeometryDESegmentation & rhs);

   ClassDef(AliMUONVGeometryDESegmentation,1) // Det element segmentation interface
};

#endif //ALI_MUON_V_GEOMETRY_DE_SEGMENTATION_H








