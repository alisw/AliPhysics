#ifndef ALIANALYSISTASKSELC2V0BACHELOR_H
#define ALIANALYSISTASKSELC2V0BACHELOR_H
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

#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliPID.h"
#include "AliAODTrack.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliNormalizationCounter.h"

/// \class AliAnalysisTaskSELc2V0bachelor

class TH1F;
class TClonesArray;
class AliAODRecoCascade;

class AliAnalysisTaskSELc2V0bachelor : public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskSELc2V0bachelor();
  AliAnalysisTaskSELc2V0bachelor(const Char_t* name, AliRDHFCutsLctoV0* cuts,
				 Bool_t useOnTheFly=kFALSE, Bool_t writeVariableTree=kTRUE, Bool_t additionalChecks=kFALSE, Bool_t trackRotation=kFALSE, Bool_t useTPCpid=kFALSE, Char_t sign=2, Bool_t checkOrigin=kFALSE);
  virtual ~AliAnalysisTaskSELc2V0bachelor();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  /// histos
  void FillLc2pK0Sspectrum(AliAODRecoCascadeHF *part, Int_t isLc,
			   Int_t &nSelectedAnal, AliRDHFCutsLctoV0 *cutsAnal,
			   TClonesArray *mcArray, Int_t originLc);

  void MakeAnalysisForLc2prK0S(AliAODEvent* aodEvent,TClonesArray *arrayLctopK0S,
			       TClonesArray *mcArray,
			       Int_t &nSelectedAnal, AliRDHFCutsLctoV0 *cutsAnal);
 
  /// set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}

  /// set flag for additional checks
  void SetAdditionalChecks(Bool_t additionalChecks) {fAdditionalChecks = additionalChecks;}
  Bool_t GetAdditionalChecks() const {return fAdditionalChecks;}

  void FillArmPodDistribution(AliAODRecoDecay *vZero,TString histoTitle, Bool_t isCandidateSelectedCuts, Bool_t isBachelorID);

  void SetUseOnTheFlyV0(Bool_t a) { fUseOnTheFlyV0=a; }
  Bool_t GetUseOnTheFlyV0() { return fUseOnTheFlyV0; }

  Int_t MatchToMClabelC(AliAODRecoCascadeHF *candidate,TClonesArray *mcArray);


  void SetNRotations(Int_t nRot=9) { fNRotations=nRot; }
  Int_t GetNRotations() { return fNRotations; }

  void SetPtMinToFillTheTree(Double_t pTmin=  0.) { fPtMinToFillTheTree=pTmin; }
  void SetPtMaxToFillTheTree(Double_t pTmax=999.) { fPtMaxToFillTheTree=pTmax; }
  Double_t GetPtMinToFillTheTree() { return fPtMinToFillTheTree; }
  Double_t GetPtMaxToFillTheTree() { return fPtMaxToFillTheTree; }

  void SetMinAngleForRot(Double_t min=5.*TMath::Pi()/6.) { fMinAngleForRot=min; }
  void SetMaxAngleForRot(Double_t max=7.*TMath::Pi()/6.) { fMaxAngleForRot=max; }
  Double_t GetMinAngleForRot() { return fMinAngleForRot; }
  Double_t GetMaxAngleForRot() { return fMaxAngleForRot; }

  void SetUseTPCPIDtoFillTree(Bool_t a) {fUseTPCPIDtoFillTree=a;}
  Bool_t GetUseTPCPIDtoFillTree() const {return fUseTPCPIDtoFillTree;}

  /// set MC usage
  void SetSign(Char_t sign) {fSign = sign;}
  Char_t GetSign() const {return fSign;}

  void SetCheckOrigin(Bool_t origin) {fCheckOrigin = origin;}
  Bool_t GetCheckOrigin() const {return fCheckOrigin;}


  void ReconstructSecVtx(Bool_t dummy=kTRUE) {fReconstructSecVtx=dummy;}
  Bool_t GetReconstructSecVtx() const {return fReconstructSecVtx;}

  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}

 private:
  
  void CheckEventSelection(AliAODEvent *aodEvent);
  void CheckEventSelectionWithCandidates(AliAODEvent *aodEvent);
  void CheckCandidatesAtDifferentLevels(AliAODRecoCascadeHF *part,AliRDHFCutsLctoV0* cutsAnal);
  void FillTheTree(AliAODRecoCascadeHF *part, AliRDHFCutsLctoV0 *cutsAnal, TClonesArray *mcArray, Int_t isLc, Int_t originLc);
  void DefineTreeVariables();

  Int_t MatchToMC(AliAODRecoCascadeHF *lc2bacV0,
		  Int_t *pdgDgLc2bacV0, Int_t *pdgDgV0,
		  TClonesArray *mcArray);

  Int_t SearchLcDaughter(TClonesArray *arrayMC, Int_t iii);

  void DefineGeneralHistograms();
  void DefineK0SHistos();
  void DefineSignalHistosSeparatedPerOrigin();
  void FillAnalysisHistograms(AliAODRecoCascadeHF *part, AliRDHFCutsLctoV0 *cutsAnal, TString appendthis);
  void TrackRotation(AliRDHFCutsLctoV0 *cutsAnal, AliAODRecoCascadeHF *part, TString appendthis);

  AliAnalysisTaskSELc2V0bachelor(const AliAnalysisTaskSELc2V0bachelor &source);
  AliAnalysisTaskSELc2V0bachelor& operator=(const AliAnalysisTaskSELc2V0bachelor& source); 
  
  Double_t Det(Double_t a00,Double_t a01,
	       Double_t a10,Double_t a11) const;
  Double_t Det(Double_t a00,Double_t a01,Double_t a02,
	       Double_t a10,Double_t a11,Double_t a12,
	       Double_t a20,Double_t a21,Double_t a22) const;
  Double_t PropagateToDCA(AliAODv0 *v, AliAODTrack *bachelor, Double_t b,
			  Double_t &xVtxLc, Double_t &yVtxLc, Double_t &zVtxLc,
			  Double_t &pxVtxLc, Double_t &pyVtxLc, Double_t &pzVtxLc);

  Double_t GetAlpha(Double_t xyz[3],Double_t pxpypz[3]);

  Int_t SearchForCommonMother(TClonesArray *mcArray,
			      Int_t dgLabels[10],Int_t ndg,
			      Int_t &ndgCk, Int_t *pdgDg, Int_t &labelMother, Int_t &nDauCand) const;

  Bool_t fUseMCInfo;          /// Use MC info
  TList *fOutput;             /// User output slot 1 // general histos
  TList *fOutputAll;          /// User output slot 4 // histos without pid and cut on V0
  TList *fOutputPIDBach;      /// User output slot 5 // histos with PID on Bachelor

  TH1F *fCEvents;                    /// Histogram to check selected events
  AliNormalizationCounter *fCounter; /// AliNormalizationCounter on output slot 2
  AliRDHFCutsLctoV0 *fAnalCuts;      /// Cuts - sent to output slot 3
  Bool_t fUseOnTheFlyV0;             /// flag to analyze also on-the-fly V0 candidates
  Int_t     fAODProtection;       /// flag to activate protection against AOD-dAOD mismatch.
  Bool_t fIsEventSelected;           /// flag for event selected

  Bool_t    fWriteVariableTree;       /// flag to decide whether to write the candidate variables on a tree variables
  TTree    *fVariablesTree;           //!<! tree of the candidate variables after track selection on output slot 4
  Float_t *fCandidateVariables;       //!<! variables to be written to the tree
  AliAODVertex *fVtx1;                /// primary vertex
  Float_t fBzkG;                      /// magnetic field value [kG]
  Bool_t fAdditionalChecks;           /// flag to fill additional histograms

  Bool_t fTrackRotation;              /// flag to check track rotation
  TList *fOutputPIDBachTR;            /// User output slot 6 // histos with PID on Bachelor and track rotation

  Double_t fMinAngleForRot; ///=5*TMath::Pi()/6;
  Double_t fMaxAngleForRot; ///=7*TMath::Pi()/6;
  Double_t fMinMass;        ///=mLcPDG-0.250;
  Double_t fMaxMass;        ///=mLcPDG+0.250;
  Int_t fNRotations;        ///=9;
  Double_t fPtMinToFillTheTree;///0.
  Double_t fPtMaxToFillTheTree;///999.
  Bool_t fUseTPCPIDtoFillTree;///kFALSE
  Char_t fSign;
  Bool_t fCheckOrigin;
  Bool_t fReconstructSecVtx;

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSELc2V0bachelor,11); /// class for Lc->p K0
  /// \endcond
};

#endif

