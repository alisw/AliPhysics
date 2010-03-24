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

//-----------------------------------------------------------------------
// Author : A. Grelli, UTRECHT
//-----------------------------------------------------------------------

#include <TH2F.h>
#include "TROOT.h"
#include "TSystem.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODEvent.h"

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
  void     SideBandBackground(Double_t finvM, Double_t finvMDStar, Double_t pt, Int_t okD0, Int_t okD0bar);
  void     WrongSignForDStar(Double_t finvM, Double_t finvMDStar, Double_t pt, Int_t okD0, Int_t okD0bar);
  //cuts
  Bool_t   SetUtrechtSelections(Double_t ptD0);  
  Bool_t   SelectPID(AliAODTrack *track, Double_t nsig);
  // histos
  Bool_t   DefineHistoFroAnalysis(); 
  // set minimum ITS clusters for the analysis
  void     SetMinITSClusters(Int_t minITSClusters) {fMinITSClusters = minITSClusters;}
  Int_t    GetMinITSClusters() const {return fMinITSClusters;}
  // set minimum for soft pion pt
  void     SetMinITSClustersSoft(Int_t minITSClustersSoft) {fMinITSClustersSoft = minITSClustersSoft;}
  Int_t    GetMinITSClustersSoft() const {return fMinITSClustersSoft;}
  // kaon PID
  void     SetPID(Bool_t usePIDforKaons) {fPID = usePIDforKaons;}
  Int_t    GetPID() const {return fPID;}
  // Set N sigmas for PID
  void     SetNSigmasPID(Int_t numberOfSigmasPID) {fNSigma = numberOfSigmasPID;}
  Int_t    GetNSigmasPID() const {return fNSigma;}
  // set MC usage
  void     SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t   GetMC() const {return fUseMCInfo;}
  
 protected:
  
  Int_t  fEvents;                //  n. of events
  AliAnalysisVertexingHF *fVHF;  //  Set the cuts
  Int_t  fMinITSClusters;        //  min n. of ITS clusters for RecoDecay
  Int_t  fMinITSClustersSoft;    //  min n. of ITS clusters for RecoDecay soft pions
  Bool_t fUseMCInfo;             //  Use MC info
  TList *fOutput;                //!  User output
  Int_t  fNSigma;                //  n sigma for kaon PID
  Bool_t fPID;                   //  PID flag
  AliAODTrack* fAODTrack;        //!
  
  // define the histograms
  TH1F *fMCDStarPt;           //!    
  TH1F *fCEvents;             //!
  TH1F *fDStarMass;           //!
  TH1F *fTrueDiff;            //!
  TH2F *fTrueDiff2;           //!
  TH1F *fInvMass;             //!
  TH1F *fInvMass1;            //!
  TH1F *fInvMass2;            //!
  TH1F *fInvMass3;            //!
  TH1F *fInvMass4;            //!
  TH1F *fInvMass5;            //!
  TH1F *fPtDStar;             //!
  TH1F *fDStar;               //!
  TH1F *fDiff;                //!
  TH1F *fDiff1;               //!
  TH1F *fDiff2;               //!
  TH1F *fDiff3;               //!
  TH1F *fDiff4;               //!
  TH1F *fDiff5;               //!
  TH1F *fDiffSideBand;        //!
  TH1F *fDiffSideBand1;       //!
  TH1F *fDiffSideBand2;       //!
  TH1F *fDiffSideBand3;       //!
  TH1F *fDiffSideBand4;       //!
  TH1F *fDiffSideBand5;       //!
  TH1F *fDiffWrongSign;       //!
  TH1F *fDiffWrongSign1;      //!
  TH1F *fDiffWrongSign2;      //!
  TH1F *fDiffWrongSign3;      //!
  TH1F *fDiffWrongSign4;      //!
  TH1F *fDiffWrongSign5;      //!
 
  ClassDef(AliAnalysisTaskSEDStarSpectra,1); // class for D* spectra
};

#endif
