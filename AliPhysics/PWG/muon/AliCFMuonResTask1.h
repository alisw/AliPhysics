#ifndef ALICFMUONRESTASK1_H
#define ALICFMUONRESTASK1_H

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
// Author : R. Vernet, Consorzio Cometa - Catania (it)
//-----------------------------------------------------------------------
// Modification done by X. Lopez - LPC Clermont (fr)
//-----------------------------------------------------------------------

#include "AliAnalysisTaskSE.h"

class TH1I;
class TParticle ;
class TLorentzVector ;
class TFile ;
class AliStack ;
class AliCFManager;
class AliESDtrack;
class AliVParticle;

class AliCFMuonResTask1 : public AliAnalysisTaskSE {
  public:

  AliCFMuonResTask1();
  AliCFMuonResTask1(const Char_t* name);
  AliCFMuonResTask1& operator= (const AliCFMuonResTask1& c);
  AliCFMuonResTask1(const AliCFMuonResTask1& c);
  virtual ~AliCFMuonResTask1();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
   
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* const io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager() const {return fCFManager;}           // get corr manager
  void           SetQAList(TList* const list) {fQAHistList = list;}

  // Data types
  Bool_t IsReadAODData()   const {return fReadAODData;}
  void   SetReadAODData   (Bool_t flag=kTRUE) {fReadAODData=flag;}

 protected:
  
  Bool_t          fReadAODData ;   // flag for AOD/ESD input files
  AliCFManager   *fCFManager   ;   // pointer to the CF manager
  TList          *fQAHistList  ;   // list of QA histograms
  TH1I           *fHistEventsProcessed; // simple histo for monitoring the number of events processed
  Double_t        fNevt        ;   // event counter
  
  Float_t Imass(Float_t e1, Float_t px1, Float_t py1, Float_t pz1,
		Float_t e2, Float_t px2, Float_t py2, Float_t p2) const;
  Float_t Rap(Float_t e, Float_t pz) const;
  Float_t Phideg(Float_t phi) const;
  
  Double_t CostCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1,
                        Double_t px2, Double_t py2, Double_t pz2, Double_t e2, Double_t energy);
  Double_t CostHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1,
                        Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t PhiCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1,
                       Double_t px2, Double_t py2, Double_t pz2, Double_t e2, Double_t energy);
  Double_t PhiHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1,
                       Double_t px2, Double_t py2, Double_t pz2, Double_t e2, Double_t energy);
  
  ClassDef(AliCFMuonResTask1,1);
};

#endif
