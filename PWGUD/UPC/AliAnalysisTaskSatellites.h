/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKSATELLITES_H
#define ALIANALYSISTASKSATELLITES_H

class TH1D;
class TH1I;
class TTree;
class TList;
class TFile;
class AliESDZDC;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSatellites : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSatellites();
  AliAnalysisTaskSatellites(const char *name);
  virtual ~AliAnalysisTaskSatellites();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  Bool_t IsSatellite(AliESDZDC *ZDCdata);

 private:
  TList *fOutputList;
   
  TTree *tOutput;

  TH1I *hDummyCounter;
  TH1I *hTriggerClassesCounter;
  TH1I *hTriggerInputsCounter;
  TH1I *hSatellitesCounter;
  Bool_t fTrgClassCINTZAC, fTrgInputV0M, fTrgInputTVX, fTrgInputV0alt, fIsSatellite;
  Bool_t fTrgClassC0V0M, fTrgClassC0TVX, fTrgClassCV0L, fTrgInputVBA, fTrgInputVBC, fTrgInputZAC;
  Int_t fRunNumber;
  UInt_t fL0inputs, fL1inputs, fTimeStamp;
  Float_t fZNATDCm[4], fZNCTDCm[4];
  
  AliAnalysisTaskSatellites(const AliAnalysisTaskSatellites&); //not implemented
  AliAnalysisTaskSatellites& operator =(const AliAnalysisTaskSatellites&); //not implemented
  
  ClassDef(AliAnalysisTaskSatellites, 1);
};

#endif
