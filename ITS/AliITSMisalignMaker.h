#ifndef ALIITSMISALIGNMAKER_H
#define ALIITSMISALIGNMAKER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//------------------------------------------------------------------------
//
// This class is a helper, producing ITS aligmnent objects.
// It provides also some useful functions
// See the parameters of the misalignment at the end of this script.
//
// Main author: L. Gaudichet
// Contact: andrea.dainese@lnl.infn.it
//
//------------------------------------------------------------------------

#include <TRandom3.h>
#include <TString.h>
#include <TClonesArray.h>


//-------------------------------------------------------------------------
class AliITSMisalignMaker : public TObject {
  
public:
  AliITSMisalignMaker();
  
  ~AliITSMisalignMaker() {};
  
  TClonesArray* GetArray() {return &fAlobj;}

  void  SetSeed(Int_t seed) {fRnd.SetSeed(seed); return;}
  
  Int_t AddAlignObj(char* name,Double_t dx,Double_t dy,Double_t dz,
		    Double_t dpsi,Double_t dtheta,Double_t dphi,
		    Bool_t unif);

  Int_t AddAlignObj(Int_t lay,Double_t dx,Double_t dy,Double_t dz,
		    Double_t dpsi,Double_t dtheta,Double_t dphi,
		    Bool_t unif);

  Int_t AddAlignObj(Int_t lay,Int_t ladd,Double_t dx,Double_t dy,Double_t dz,
		    Double_t dpsi,Double_t dtheta,Double_t dphi,
		    Bool_t unif);

  Int_t AddSectorAlignObj(Int_t sectMin,Int_t sectMax,
			  Double_t dx,Double_t dy,Double_t dz,
			  Double_t dpsi,Double_t dtheta,Double_t dphi,
			  Double_t xShift,Double_t yShift,Double_t zShift,
			  Double_t psiShift,Double_t thetaShift,Double_t phiShift,
			  Bool_t unif);

  TString  GetSymbName(Int_t layer) const;
  TString  GetSymbName(Int_t layer,Int_t ladd) const;
  TString  GetHalfStaveLadderSymbName(Int_t layer,Int_t ladd,Int_t halfStave) const;
  TString  GetSymbName(Int_t layer,Int_t ladd,Int_t mod) const;
  Double_t GaussCut(Double_t mean,Double_t sigma,Double_t cutn);
  
  static Int_t GetNLayers()             {return kNLayers;}
  static Int_t GetNLadders(Int_t lay)   {return fgkNLadders[lay];}
  static Int_t GetNDetectors(Int_t lay) {return fgkNDetectors[lay];}

protected:
  TRandom3     fRnd; // TRandom object
  Int_t        fInd; // index of current AliAlignObjParams in fAlobj
  TClonesArray fAlobj; // array of AliAlignObjParams
  TString      fStrSPD; // name of SPD
  TString      fStrSDD; // name of SDD
  TString      fStrSSD; // name of SSD
  TString      fStrStave; // name of SPD stave
  TString      fStrHalfStave; // name of SPD half-stave
  TString      fStrLadder; // name of SPD ladder
  TString      fStrSector; // name of SPD sector
  TString      fStrSensor; // name of sensitive volume

private:
  enum {kNLayers = 6}; // The number of layers.
  static const Int_t  fgkNLadders[kNLayers];  // Array of the number of ladders/layer(layer)
  static const Int_t  fgkNDetectors[kNLayers];// Array of the number of detector/ladder(layer)

  ClassDef(AliITSMisalignMaker,1)   //ITS Misalign Maker
};


#endif

