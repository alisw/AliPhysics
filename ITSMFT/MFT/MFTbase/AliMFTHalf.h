#ifndef AliMFTHalf_H
#define AliMFTHalf_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTHalf
/// \brief Class describing geometry of one half of the ALICE Muon Forward Tracker
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TNamed.h"
#include "TGeoVolume.h"

class AliMFTHalfSegmentation;
//=============================================================================================


class AliMFTHalf : public TNamed {
  
public:
  
  AliMFTHalf();
  AliMFTHalf(AliMFTHalfSegmentation *segmentation);
  
  virtual ~AliMFTHalf();
  
  /// \brief Returns the Volume holding the Half-MFT
  TGeoVolumeAssembly * GetVolume() {return fHalfVolume;};
  
protected:
  TGeoVolumeAssembly * fHalfVolume;

private:
  AliMFTHalfSegmentation * fSegmentation; ///< \brief Pointer to the half-MFT segmentation
  void CreateHalfDisks();

  
  /// \cond CLASSIMP
  ClassDef(AliMFTHalf, 1);
  /// \endcond

};

//=============================================================================================

#endif

