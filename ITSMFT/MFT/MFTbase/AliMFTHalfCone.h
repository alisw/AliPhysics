#ifndef AliMFTHalfCone_H
#define AliMFTHalfCone_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//=============================================================================================
//
//      Class describing geometry of one MFT half-cone
//
//      Contact author: sbest@pucp.pe, eric.endress@gmx.de, franck.manso@clermont.in2p3.fr
//
//=============================================================================================

#include "TNamed.h"
#include "AliLog.h"
#include "TGeoVolume.h"

//=============================================================================================


class AliMFTHalfCone : public TNamed {
  
public:
  
  AliMFTHalfCone();
  
  virtual ~AliMFTHalfCone();
  
  TGeoVolumeAssembly* CreateHalfCone(Int_t half);

  
protected:

  TGeoVolumeAssembly * fHalfCone;

private:
  
  ClassDef(AliMFTHalfCone,1)
  
};

//=============================================================================================

#endif

