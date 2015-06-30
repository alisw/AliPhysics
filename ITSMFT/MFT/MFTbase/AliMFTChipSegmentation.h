#ifndef AliMFTChipSegmentation_H
#define AliMFTChipSegmentation_H 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTChipSegmentation
/// \brief Chip Segmentation description
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "AliMFTVSegmentation.h"

//====================================================================================================================================================

class AliMFTChipSegmentation : public AliMFTVSegmentation {

public:

  AliMFTChipSegmentation();
  AliMFTChipSegmentation(UInt_t uniqueID);
  
  virtual ~AliMFTChipSegmentation(){};
  virtual void Clear(const Option_t* /*opt*/) {;}
  virtual void Print(Option_t* /*option*/);

  /// \brief Transform (x,y) Hit coordinate into Pixel ID on the matrix
  Bool_t Hit2PixelID(Double_t xHit, Double_t yHit, Int_t &xPixel, Int_t &yPixel);
  
private:
  
  /// \cond CLASSIMP
  ClassDef(AliMFTChipSegmentation, 1);
  /// \endcond

};

//====================================================================================================================================================
	
#endif

