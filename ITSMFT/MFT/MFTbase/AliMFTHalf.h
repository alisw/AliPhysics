#ifndef AliMFTHalf_H
#define AliMFTHalf_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//=============================================================================================
//
//      Class describing geometry of one half of the ALICE Muon Forward Tracker
//
//      Contact author: raphael.tieulent@cern.ch
//
//=============================================================================================

#include "TNamed.h"
#include "AliLog.h"
#include "TGeoVolume.h"

class AliMFTHalfSegmentation;
//=============================================================================================


class AliMFTHalf : public TNamed {
  
public:
  
  AliMFTHalf();
  AliMFTHalf(AliMFTHalfSegmentation *segmentation);
  
  virtual ~AliMFTHalf();
  
  TGeoVolumeAssembly * GetVolume() {return fHalfVolume;};
  void Init();

  
protected:
  TGeoVolumeAssembly * fHalfVolume;

private:
  AliMFTHalfSegmentation * fSegmentation;
  void CreateHalfDisks();

  
  ClassDef(AliMFTHalf,1)
  
};

//=============================================================================================

#endif

