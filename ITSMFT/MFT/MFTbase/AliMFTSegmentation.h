#ifndef AliMFTSegmentation_H
#define AliMFTSegmentation_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTSegmentation
/// \brief Class for the virtual segmentation of the ALICE Muon Forward Tracker
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TNamed.h"
#include "TClonesArray.h"
class AliMFTHalfSegmentation;
class AliMFTPlane;

//====================================================================================================================================================

class AliMFTSegmentation : public TNamed {

public:
  
  enum THalfMFTType { kBottom, kTop };

  AliMFTSegmentation();
  AliMFTSegmentation(const Char_t *nameGeomFile);
  
  virtual ~AliMFTSegmentation();
  virtual void Clear(const Option_t* /*opt*/);

  /// \brief Returns pointer to the segmentation of the half-MFT
  /// \param iHalf Integer : 0 = Bottom; 1 = Top
  /// \return Pointer to a AliMFTHalfSegmentation
  AliMFTHalfSegmentation* GetHalf(Int_t iHalf) const { return ((iHalf==kTop || iHalf==kBottom) ? ( (AliMFTHalfSegmentation*) fMFTHalves->At(iHalf)) :  NULL); }

  Int_t GetDetElemLocalID(Int_t half, Int_t disk, Int_t ladder, Int_t sensor) const;
  
  Bool_t Hit2PixelID(Double_t xHit, Double_t yHit, Double_t zHit, Int_t half, Int_t disk, Int_t ladder, Int_t sensor, Int_t &xPixel, Int_t &yPixel);


private:

  TClonesArray *fMFTHalves; ///< \brief Array of pointer to AliMFTHalfSegmentation

  /// \cond CLASSIMP
  ClassDef(AliMFTSegmentation, 1);
  /// \endcond
  
};

//====================================================================================================================================================

#endif

