#ifndef AliMFTFlex_H
#define AliMFTFlex_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//  Flex (Printed Cabled Board) class for ALICE MFT upgrade
//  This version uses TGeo
//  Authors:
//  F. Manso
//-------------------------------------------------------------------------

#include "TNamed.h"


class TGeoVolume;
class TGeoVolumeAssembly;

class AliMFTFlex : public TNamed {

public:
  AliMFTFlex();
  AliMFTFlex(AliMFTLadderSegmentation *ladder);
  virtual ~AliMFTFlex();
  TGeoVolumeAssembly*  MakeFlex(Int_t nbsensors, Double_t length);
  void Make_ElectricComponents(TGeoVolumeAssembly*  flex, Int_t nbsensors, Double_t length, Double_t zvarnish);

private:
  TGeoVolume*  Make_Lines(Int_t nbsensors, Double_t length, Double_t width, Double_t thickness);
  TGeoVolume*  Make_AGND_DGND(Double_t length, Double_t width, Double_t thickness);
  TGeoVolume*  Make_Kapton(Double_t length, Double_t width, Double_t thickness);
  TGeoVolume*  Make_Varnish(Double_t length, Double_t width, Double_t thickness, Int_t iflag);
  TGeoVolume*  Make_ElectricComponent(Double_t dx, Double_t dy, Double_t dz, Int_t iflag);

  Double_t *fFlexOrigin;
  AliMFTLadderSegmentation * fLadderSeg;

  ClassDef(AliMFTFlex,1)
};

#endif
