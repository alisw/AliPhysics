/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup geometry
/// \class AliMUONVGeometryDESegmentation
/// \brief Extension for AliSegmentation interface for detection elements
///
/// \deprecated - To be removed with passing to maping segmentation interface
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_V_GEOMETRY_DE_SEGMENTATION_H
#define ALI_MUON_V_GEOMETRY_DE_SEGMENTATION_H

#include "AliSegmentation.h"
#include "AliMUONGeometryDirection.h"

class AliMpVSegmentation;

class AliMUONVGeometryDESegmentation : public AliSegmentation
{
  public:
    AliMUONVGeometryDESegmentation();
    virtual ~AliMUONVGeometryDESegmentation();

    // methods
                       /// Return true if a pad exists in the given position
    virtual Bool_t  HasPad(Float_t x, Float_t y, Float_t z) = 0; 
                       /// Return true if a pad with given indices exists
    virtual Bool_t  HasPad(Int_t ix, Int_t iy) = 0;

                       /// Return the direction with a constant pad size  
		       /// (Direction or coordinate where the spatial resolution 
		       /// is the best) \n
                       /// Normally kDirY will correspond with cathode segmentation 
		       /// for the bending plane and kDirX with cathode segmentation 
		       /// for the non bending plane
    virtual AliMUONGeometryDirection  GetDirection() = 0;
		       
                       /// Access to mapping
    virtual const AliMpVSegmentation* GetMpSegmentation() const = 0; 		       

  protected:
    AliMUONVGeometryDESegmentation(const AliMUONVGeometryDESegmentation& rhs);
    AliMUONVGeometryDESegmentation& operator=(
      const AliMUONVGeometryDESegmentation& rhs);

   ClassDef(AliMUONVGeometryDESegmentation,1) // Det element segmentation interface
};

#endif //ALI_MUON_V_GEOMETRY_DE_SEGMENTATION_H








