#ifndef AliMFTHalfDisk_H
#define AliMFTHalfDisk_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTbase
/// \class AliMFTHalfDisk
/// \brief Class Building geometry of one half of a MFT disk
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date June 9th, 2015

#include "TNamed.h"
#include "TGeoVolume.h"
#include "AliMFTConstants.h"

class AliMFTHalfDiskSegmentation;
class AliMFTSupport;
class AliMFTHeatExchanger;
//=============================================================================================


class AliMFTHalfDisk : public TNamed {
  
public:

  AliMFTHalfDisk();
  AliMFTHalfDisk(AliMFTHalfDiskSegmentation *segmentation);
  TGeoVolumeAssembly * CreateHeatExchanger();
  void CreateLadders();

  virtual ~AliMFTHalfDisk();
  
  /// \brief Returns a pointer to the Volume Assembly describing the entire half-disk
  TGeoVolumeAssembly * GetVolume() {return fHalfDiskVolume;};
  
private:

  AliMFTSupport    * fMFTSupport;             ///< \brief Disk Support
  AliMFTHeatExchanger * fMFTHeatExchanger;    ///< \brief Heat Exchanger
  TGeoVolumeAssembly * fHalfDiskVolume;       ///< \brief Half-Disk Volume
  AliMFTHalfDiskSegmentation * fSegmentation; ///< \brief Virtual Segmentation of the half-disk

  /// \cond CLASSIMP
  ClassDef(AliMFTHalfDisk, 1);
  /// \endcond
 
};

//=============================================================================================

#endif

