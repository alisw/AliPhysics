#ifndef ALIPMD_H
#define ALIPMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:PMD      //
////////////////////////////////////////////////
 
#include "AliDetector.h"

class AliLoader;
class TClonesArray;
class TFile;
class AliPMDRecPoint;
class AliPMDLoader;
class AliPMDhit;

class AliPMD : public AliDetector {

public:
  AliPMD();
  AliPMD(const char *name, const char *title);

  virtual AliLoader* MakeLoader(const char* topfoldername);

  virtual      ~AliPMD();
  virtual void  AddHit(Int_t track, Int_t* vol, Float_t* hits);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials() {}
  virtual Int_t IsVersion() const =0;
  virtual void  SetPAR(Float_t p1, Float_t p2, Float_t p3, Float_t p4);
  virtual void  SetIN(Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5);
  virtual void  SetGEO(Float_t p1, Float_t p2, Float_t p3);
  virtual void  SetPadSize(Float_t p1, Float_t p2, Float_t p3, Float_t p4);
  virtual void  StepManager();
  virtual void  MakeBranch(Option_t* option);
  virtual void  SetTreeAddress();
  
  virtual void  Hits2SDigits();
  virtual void  SDigits2Digits();
  virtual void  Hits2Digits();

  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;

  virtual void  Digits2Raw();
  virtual Bool_t Raw2SDigits(AliRawReader *rawReader);
  
 protected:
  Float_t fPar[4];           // pmdin, pmdout, thgas, thcell
  Float_t fIn[5];            // thmin, thmax, zdist, thlow, thhigh
  Float_t fGeo[3];           // wafer, edge, numqu
  Float_t fPadSize[4];       // size of the pads
  Int_t   fNumPads[4];       // number of the pads

  ClassDef(AliPMD,8)  // Base Class for Photon Multiplicity Detector
};
#endif
