#ifndef ALICFMUONSINGLETASK1_H
#define ALICFMUONSINGLETASK1_H

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */ 

//-----------------------------------------------------------------------
// Author : X. Lopez - LPC Clermont (fr)
//-----------------------------------------------------------------------

#include "AliAnalysisTaskSE.h"

class TH1I;
class TNtuple;
class TParticle ;
class TFile ;
class AliStack ;
class AliCFManager;
class AliESDtrack;
class AliVParticle;

class AliCFMuonSingleTask1 : public AliAnalysisTaskSE {
  public:

  AliCFMuonSingleTask1();
  AliCFMuonSingleTask1(const Char_t* name);
  AliCFMuonSingleTask1& operator= (const AliCFMuonSingleTask1& c);
  AliCFMuonSingleTask1(const AliCFMuonSingleTask1& c);
  virtual ~AliCFMuonSingleTask1();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
   
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* const io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager() const {return fCFManager;}                 // get corr manager
  void           SetQAList(TList* const list) {fQAHistList = list;}

  // Data types
  Bool_t IsReadAODData()   const {return fReadAODData;}
  void   SetReadAODData   (Bool_t flag=kTRUE) {fReadAODData=flag;}
  void SetUseMC(Bool_t isMC)       { fIsMC         = isMC;}

 protected:
  
  Bool_t          fReadAODData ;       // flag for AOD/ESD input files
  AliCFManager   *fCFManager   ;       // pointer to the CF manager
  TList          *fQAHistList  ;       // list of QA histograms
  TH1I           *fHistEventsProcessed;// simple histo for monitoring the number of events processed
  Int_t           fNevt        ;       // event countor

  Float_t Rap(Float_t e, Float_t pz);
  Float_t Phideg(Float_t phi);
  Bool_t fIsMC;                        // flag of whether the input is MC

  ClassDef(AliCFMuonSingleTask1,1);
};

#endif
