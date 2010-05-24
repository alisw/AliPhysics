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

#include <TH2F.h>
#include "TROOT.h"
#include "TSystem.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODEvent.h"
#include "AliPID.h"


class TH2F;
class TH1I;
class TParticle;
class TFile;
class TClonesArray;
class AliCFManager;
class AliAODRecoDecay;
class AliAODRecoDecayHF2Prong;
class AliAODMCParticle;

class AliAnalysisTaskSEDStarSpectra : public AliAnalysisTaskSE {
  
 public:
  
  AliAnalysisTaskSEDStarSpectra();
  AliAnalysisTaskSEDStarSpectra(const Char_t* name);
  AliAnalysisTaskSEDStarSpectra& operator= (const AliAnalysisTaskSEDStarSpectra& c);
  AliAnalysisTaskSEDStarSpectra(const AliAnalysisTaskSEDStarSpectra& c);
  virtual ~AliAnalysisTaskSEDStarSpectra();
  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  //Background simulation
  void     SideBandBackground(Int_t ptbin, AliAODRecoCascadeHF *part, Bool_t PIDon, Int_t nSigma, AliAnalysisVertexingHF *vhf, TList *listout);
  void     WrongSignForDStar(Int_t ptbin, AliAODRecoCascadeHF *part, Bool_t PIDon, Int_t nSigma, AliAnalysisVertexingHF *vhf, TList *listout);
  //cuts
  Bool_t   SingleTrackSelections(AliAODRecoDecayHF2Prong* theD0particle, AliAODTrack *track0, AliAODTrack *track1, AliAODTrack *track2);
  Bool_t   SetSelections(Double_t pt);  
  Bool_t   SelectPID(AliAODTrack *track, AliPID::EParticleType pid, Double_t nsig);
  // histos
  void   FillSpectrum(Int_t ptbin, AliAODRecoCascadeHF *part, Int_t isDStar, Bool_t PIDon, Int_t nSigma, AliAnalysisVertexingHF *vhf, TList *listout);
  void     DefineHistograms();
  // set analysis type
  void     SetAnalysisType(Int_t anaType) {fAnalysis = anaType;}
  void     PrintAnalysisType() {printf("Analysis type: %d\n(0: Heidelberg\t1: Utrecht)",fAnalysis);}
  // set minimum ITS clusters for the analysis
  void     SetMinITSClusters(Int_t minITSClusters) {fMinITSClusters = minITSClusters;}
  Int_t    GetMinITSClusters() const {return fMinITSClusters;}
  // set minimum for soft pion pt
  void     SetMinITSClustersSoft(Int_t minITSClustersSoft) {fMinITSClustersSoft = minITSClustersSoft;}
  Int_t    GetMinITSClustersSoft() const {return fMinITSClustersSoft;}
  // kaon PID
  void     SetPID(Bool_t usePID) {fPID = usePID;}
  Int_t    GetPID() const {return fPID;}
  // Set N sigmas for PID
  void     SetNSigmasPID(Int_t numberOfSigmasPID) {fNSigma = numberOfSigmasPID;}
  Int_t    GetNSigmasPID() const {return fNSigma;}
  // set MC usage
  void     SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t   GetMC() const {return fUseMCInfo;}
  
 protected:
  
  Int_t  fEvents;                //  n. of events
  Int_t  fAnalysis;		 //  0: HD;	1: UU;
  AliAnalysisVertexingHF *fVHF;  //  Set the cuts
  AliAnalysisVertexingHF *fVHFloose;  //  Set the cuts
  Double_t fD0Window;		 //  select width on D0Mass
  Double_t fPeakWindow;          //  select width on DstarMass
  Int_t  fMinITSClusters;        //  min n. of ITS clusters for RecoDecay
  Int_t  fMinITSClustersSoft;    //  min n. of ITS clusters for RecoDecay soft pions
  Bool_t fUseMCInfo;             //  Use MC info
  TList *fOutput;                //!  User output
  TList *fOutputSpectrum;
  TList *fOutputAll;
  TList *fOutputPID3;
  TList *fOutputPID2;
  TList *fOutputPID1;
  Int_t  fNSigma;                //  n sigma for kaon PID
  Bool_t fPID;                   //  PID flag
  AliAODTrack* fAODTrack;        //!
  
  // define the histograms
  TH1F *fCEvents;             //!
  TH2F *fTrueDiff2;           //!
 
  ClassDef(AliAnalysisTaskSEDStarSpectra,6); // class for D* spectra
};

#endif

