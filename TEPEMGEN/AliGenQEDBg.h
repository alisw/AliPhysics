#ifndef ALIGENQEDBG_H
#define ALIGENQEDBG_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * Copyright(c) 1997, 1998, 2002, Adrian Alscher and Kai Hencken          *
 * Copyright(c) 2002 Kai Hencken, Yuri Kharlov, Serguei Sadovsky          *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Event generator of single e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 TeV/nucleon.
// Author: Yuri.Kharlov@cern.ch
// 9 October 2002

#include "AliGenMC.h"
#include "AliGenEventHeader.h"
#include "AliGenEpEmv1.h"

//-------------------------------------------------------------
class AliGenQEDBg : public AliGenEpEmv1
{
public:
  AliGenQEDBg();
  
  virtual ~AliGenQEDBg();
  virtual void Generate();
  virtual void Init();
  //
  Double_t  GetLuminosity()       const {return fLumi;}
  Double_t  GetIntegrationTime()  const {return fIntTime;}
  Double_t  GetXSection()         const {return fXSection;}
  Double_t  GetXSectionEps()      const {return fXSectionEps;}
  Double_t  GetMeanNPairs()       const {return fPairsInt;}
  //
  void      SetLumiIntTime(double lumi, double intTime);
  void      SetXSectionEps(double eps=1e-2)  {fXSectionEps = eps>0 ? eps:1e-2;}
  void      SetMinMaxXSTest(double mn,double mx);
  //
 protected:
  AliGenQEDBg(const AliGenQEDBg & gen);
  AliGenQEDBg & operator=(const AliGenQEDBg & gen);
  //
  Double_t   fLumi;         // beam luminsity
  Double_t   fXSection;     // estimated cross section in k-barns
  Double_t   fXSectionEps;  // error with wich Xsection is calculated
  Double_t   fIntTime;      // integration time in seconds
  Double_t   fPairsInt;     // estimated average number of pairs in IntTime
  Double_t   fMinXSTest;    // min number of generator calls for Xsection estimate
  Double_t   fMaxXSTest;    // max number of generator calls for Xsection estimate

  //
  ClassDef(AliGenQEDBg,1) // Generator e+e- pair background from PbPb QED interactions
};
#endif
