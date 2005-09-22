/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup geometry
/// \class AliMUONVGeometryDESegmentation
/// \brief Extension for AliSegmentation interface for detection elements
///
/// Extension for AliSegmentation interface,
/// added functions:                                                         \n
/// Bool_t  HasPad(Float_t x, Float_t y, Float_t z);                         \n
/// Bool_t  HasPad(Int_t ix, Int_t iy);                                      \n
///
/// Author:Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_V_GEOMETRY_DE_SEGMENTATION_H
#define ALI_MUON_V_GEOMETRY_DE_SEGMENTATION_H

#include "AliSegmentation.h"
#include "AliMUONGeometryDirection.h"

class AliMUONSegmentManuIndex;

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

    virtual AliMUONGeometryDirection  GetDirection() = 0;
                       // Returns the direction with a constant pad size  
		       // (Direction or coordinate where the spatial resolution 
		       // is the best)
                       // Normally kDirY will correspond with cathode segmentation 
		       // for the bending plane and kDirX with cathode segmentation 
		       // for the non bending plane

  protected:
    AliMUONVGeometryDESegmentation(const AliMUONVGeometryDESegmentation& rhs);
  
    // operators
    AliMUONVGeometryDESegmentation& operator=(
      const AliMUONVGeometryDESegmentation & rhs);

   ClassDef(AliMUONVGeometryDESegmentation,1) // Det element segmentation interface
};

#endif //ALI_MUON_V_GEOMETRY_DE_SEGMENTATION_H








