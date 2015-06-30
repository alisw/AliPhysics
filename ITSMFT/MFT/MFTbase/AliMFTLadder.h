#ifndef AliMFTLadder_H
#define AliMFTLadder_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//=============================================================================================
//
//      Class describing geometry of one MFT Ladder
//
//      Contact author: raphael.tieulent@cern.ch
//
//=============================================================================================

#include "TNamed.h"
#include "AliLog.h"
#include "TGeoVolume.h"

class AliMFTLadderSegmentation;
class AliMFTFlex;
//=============================================================================================


class AliMFTLadder : public TNamed {
  
public:
  
  AliMFTLadder();
  AliMFTLadder(AliMFTLadderSegmentation *segmentation);
  
  virtual ~AliMFTLadder();
  
  TGeoVolume * CreateVolume();
  void CreateSensors();

  
protected:
  const static Double_t kLadderDeltaY;
  const static Double_t kLadderDeltaZ;
  AliMFTLadderSegmentation *fSegmentation;
  AliMFTFlex      * fMFTFlex;          // ! Ladder's Flex
  TGeoVolume * fLadderVolume;

private:
  
  ClassDef(AliMFTLadder,1)
  
};

//=============================================================================================

#endif

