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


///////////////////////////////////////////////////////////////////////////
//          ----   CORRECTION FRAMEWORK   ----
// AliCFAcceptanceCuts implementation
// Class to cut on the number of AliTrackReference's 
// for each detector
///////////////////////////////////////////////////////////////////////////
// author : R. Vernet (renaud.vernet@cern.ch)
///////////////////////////////////////////////////////////////////////////


#ifndef ALICFACCEPTANCECUTS_H
#define ALICFACCEPTANCECUTS_H

#include "AliCFCutBase.h"

class AliMCEvent;
class TH1F ;
class TH2F ;
class TBits;

class AliCFAcceptanceCuts : public AliCFCutBase
{
 public :
  AliCFAcceptanceCuts() ;
  AliCFAcceptanceCuts(const Char_t* name, const Char_t* title) ;
  AliCFAcceptanceCuts(const AliCFAcceptanceCuts& c) ;
  AliCFAcceptanceCuts& operator=(const AliCFAcceptanceCuts& c) ;
  virtual ~AliCFAcceptanceCuts() { };
  virtual Bool_t IsSelected(TObject* obj) ;
  virtual Bool_t IsSelected(TList*  /*list*/) {return kTRUE;}
  virtual void   SetEvtInfo(TObject* mcInfo) ;
  void SetMinNHitITS (Int_t nhits) {fMinNHitITS=nhits;} 
  void SetMinNHitTPC (Int_t nhits) {fMinNHitTPC=nhits;} 
  void SetMinNHitTRD (Int_t nhits) {fMinNHitTRD=nhits;} 
  void SetMinNHitTOF (Int_t nhits) {fMinNHitTOF=nhits;} 
  void SetMinNHitMUON(Int_t nhits) {fMinNHitMUON=nhits;}

  enum { 
    kCutHitsITS ,
    kCutHitsTPC ,
    kCutHitsTRD ,
    kCutHitsTOF ,
    kCutHitsMUON,
    kNCuts,           // number of single selections
    kNStepQA=2        // number of QA steps (before/after the cuts)
  };

 protected:
  AliMCEvent *fMCInfo;        // pointer to MC Information
  Int_t       fMinNHitITS ;   // min number of track references in ITS 
  Int_t       fMinNHitTPC ;   // min number of track references in TPC 
  Int_t       fMinNHitTRD ;   // min number of track references in TRD 
  Int_t       fMinNHitTOF ;   // min number of track references in TOF 
  Int_t       fMinNHitMUON ;  // min number of track references in MUON
  
  //QA histos
  TH1F*  fhCutStatistics;		// Histogram: statistics of what cuts the tracks did not survive
  TH2F*  fhCutCorrelation;		// Histogram: 2d statistics plot
  TH1F*  fhQA[kNCuts][kNStepQA];        // QA Histograms
  TBits* fBitmap ; 			// stores single selection decisions
  void SelectionBitMap(TObject* obj);
  void FillHistograms(TObject* obj, Bool_t afterCuts);
  void AddQAHistograms(TList *qaList) ;
  void DefineHistograms();

  ClassDef(AliCFAcceptanceCuts,1);
};

#endif
