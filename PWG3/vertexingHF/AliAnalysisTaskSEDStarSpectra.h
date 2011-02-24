#ifndef ALIANALYSISTASKSEDSTARSPECTRA_H
#define ALIANALYSISTASKSEDSTARSPECTRA_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

#include <TH2F.h>
#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliPID.h"
#include "AliRDHFCutsDStartoKpipi.h"


class AliAnalysisTaskSEDStarSpectra : public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskSEDStarSpectra();
  AliAnalysisTaskSEDStarSpectra(const Char_t* name,AliRDHFCutsDStartoKpipi* cuts);
  virtual ~AliAnalysisTaskSEDStarSpectra();

  // Implementation of interface methods  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  //Background simulation
  void     SideBandBackground(AliAODRecoCascadeHF *part, Bool_t PIDon, Int_t nSigma, AliRDHFCutsDStartoKpipi *cuts, TList *listout);
  void     WrongSignForDStar(AliAODRecoCascadeHF *part, Bool_t PIDon, Int_t nSigma, AliRDHFCutsDStartoKpipi *cuts, TList *listout);
  //cuts
  Bool_t   SingleTrackSelections(const AliAODRecoDecayHF2Prong* theD0particle, const AliAODTrack *track2);
  Bool_t   SelectPID(const AliAODTrack *track, AliPID::EParticleType pid, Double_t nsig);
  Bool_t   SelectTOFPID(const AliAODRecoDecayHF2Prong* d, const AliAODTrack *tracksoft);
  // histos
  void   FillSpectrum(AliAODRecoCascadeHF *part, Int_t isDStar, Bool_t PIDon, Int_t nSigma, AliRDHFCutsDStartoKpipi *cuts, TList *listout);
  void     DefineHistograms();
  // set analysis type
  void     SetAnalysisType(Int_t anaType) {fAnalysis = anaType;}
  void     PrintAnalysisType() {printf("Analysis type: %d\n(0: Heidelberg\t1: Utrecht)",fAnalysis);}

  // kaon PID
  void     SetPID(Bool_t usePID) {fPID = usePID;}
  Int_t    GetPID() const {return fPID;}
  // Set N sigmas for PID
  void     SetNSigmasPID(Int_t numberOfSigmasPID) {fNSigma = numberOfSigmasPID;}
  Int_t    GetNSigmasPID() const {return fNSigma;}
  // set MC usage
  void     SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t   GetMC() const {return fUseMCInfo;}
  
 private:
  
  AliAnalysisTaskSEDStarSpectra(const AliAnalysisTaskSEDStarSpectra &source);
  AliAnalysisTaskSEDStarSpectra& operator=(const AliAnalysisTaskSEDStarSpectra& source); 
  
  Int_t  fEvents;                //  n. of events
  Int_t  fAnalysis;		 //  0: HD;	1: UU;
  Double_t fD0Window;		 //  select width on D0Mass
  Double_t fPeakWindow;          //  select width on DstarMass
  Bool_t fUseMCInfo;             //  Use MC info
  TList *fOutput;                //!  User output
  TList *fOutputSpectrum;        //!  User output1
  TList *fOutputAll;             //!  User output2
  TList *fOutputPID3;            //!  User output3
  TList *fOutputPID2;            //!  User output4
  TList *fOutputPID1;            //!  User output5
  Int_t  fNSigma;                //  n sigma for kaon PID
  Bool_t fPID;                   //  PID flag
  AliRDHFCutsDStartoKpipi *fCuts; // Cuts - sent to output slot 7
  // define the histograms
  TH1F *fCEvents;             //!
  TH2F *fTrueDiff2;           //!
 
  ClassDef(AliAnalysisTaskSEDStarSpectra,7); // class for D* spectra
};

#endif

