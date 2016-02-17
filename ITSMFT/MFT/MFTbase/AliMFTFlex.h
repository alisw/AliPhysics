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
  TGeoVolumeAssembly*   MakeFlex(Int_t nsensors);

private:
  TGeoVolumeAssembly* MakeAluTrack(Int_t nsensors);
  TGeoVolume*         MakeGlueLayer(Int_t nlayer, Double_t thickness, Double_t widthflex, Double_t lenghtflex);
  TGeoVolume*         MakeKaptonLayer(Int_t nlayer, Double_t thickness, Double_t widthflex, Double_t lenghtflex);
  TGeoVolume*         MakeVarnishLayer(Double_t thickness, Double_t widthflex, Double_t lenghtflex);
  
  Double_t  fTrackWidth;
  Double_t  fTrackThickness;
  Double_t  fTrackNb;
  Double_t  fGlueThickness;
  Double_t  fKaptonThickness;
  Double_t  fVarnishThickness;
  Int_t     fNLayer;
  Int_t     fNSensors;
  Double_t  fChipWidth;
  Double_t  fChipInterspace;
  Double_t  fChipSideOffset;
  Double_t *fFlexDimensions;
  Double_t *fFlexOrigin;
  Bool_t fIsLeftType;
  AliMFTLadderSegmentation * fLadderSeg;
  
  ClassDef(AliMFTFlex,1)
};

#endif
