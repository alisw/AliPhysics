#ifndef ALIEBYEPIDRSTIOPHY_H
#define ALIEBYEPIDRSTIOPHY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    //
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================// 

#include "THnSparse.h"
#include "TList.h"
#include "TRandom3.h"

#include "AliEbyEPidRatioBase.h"

class AliEbyEPidRatioPhy : public AliEbyEPidRatioBase {

 public:

  AliEbyEPidRatioPhy();
  virtual ~AliEbyEPidRatioPhy();
  virtual void Process();
  void SetOutList(TList* l) {fOutList = l;}
  void SetQA()     {fIsQA   = kTRUE;}  
  void SetSubRun() { fIsSub = kTRUE;}
  void SetBSRun()  {fIsBS   = kTRUE;}
  void SetIsPer()  {fIsPer  = kTRUE;}

 private:

  AliEbyEPidRatioPhy(const AliEbyEPidRatioPhy&);
  AliEbyEPidRatioPhy& operator=(const AliEbyEPidRatioPhy&); 
  virtual void Init();
  virtual void Reset();
  virtual void CreateHistograms();  
  Int_t ProcessTracks();
  Int_t ProcessParticles();  
 
  void ResetHistSet();
  void AddHistSetCent(const Char_t *name, const Char_t *title);
  void AddHistSetCentPt(const Char_t *name, const Char_t *title);
  void AddHistSetRatio(const Char_t *name, const Char_t *title);
  void AddHistInGroup(const Char_t *name, const Char_t *title, Int_t nSample);
 

  void FillHistSetCent(const Char_t *name, Int_t idx, Bool_t isMC);
  void FillHistSetRatio(const Char_t *name, Int_t idx, Bool_t isMC);
  void FillHistSetCentPt(const Char_t *name, Int_t idx, Bool_t isMC);
  void FillHistInGroup(const Char_t *name, Int_t idx, Int_t iSub, Bool_t isMC);
 
  
  TList                *fOutList;               //! Output data container
  Int_t                 fOrder;                 //  Max order of higher order distributions
  Int_t                 fNNp;                   //  N sets of arrays of particle/anti-particle counts
  Int_t              ***fNp;                    //!  Array of particle/anti-particle counts
  Int_t             ****fNpPt;                  //!  Array of particle/anti-particle counts
  Int_t                 fNMCNp;                 //  N sets of arrays of MC particle/anti-particle counts
  Int_t              ***fMCNp;                  //!  Array of MC particle/anti-particle counts
  Int_t             ****fMCNpPt;                //!  Array of MC particle/anti-particle counts
  Double_t            **fRedFactp;              //!  Array of particle/anti-particle reduced factorial
  TH1F                 *fPtBinHist;             // Hist
  Bool_t                fIsQA;                  // Check for QA
  Bool_t                fIsSub;                 //
  Bool_t                fIsBS;                  //
  Bool_t                fIsPer;                 //
  TRandom3             *fRan;                   //  
  THnSparseD           *fHnTrackUnCorrRec;      // THnSparseD : uncorrected probe particles
  THnSparseD           *fHnTrackUnCorrMc;       // THnSparseD : Original Corrected probe particles Mc


  ClassDef(AliEbyEPidRatioPhy, 1);
};

#endif
