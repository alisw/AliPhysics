#ifndef AliMFTHalfDisk_H
#define AliMFTHalfDisk_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//=============================================================================================
//
//      Class describing geometry of one half of a MFT disk
//
//      Contact author: raphael.tieulent@cern.ch
//
//=============================================================================================

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
  
  TGeoVolumeAssembly * GetVolume() {return fHalfDiskVolume;};
  
  
protected:
  AliMFTSupport    * fMFTSupport;       // ! Disk Support
  AliMFTHeatExchanger * fMFTHeatExchanger;    // ! Heat Exchanger
  TGeoVolumeAssembly * fHalfDiskVolume;
  AliMFTHalfDiskSegmentation * fSegmentation;
  
private:

  ClassDef(AliMFTHalfDisk,1)
  
};

//=============================================================================================

#endif

