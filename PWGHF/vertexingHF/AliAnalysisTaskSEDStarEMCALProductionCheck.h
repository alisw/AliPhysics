#ifndef AliAnalysisTaskSEDStarEMCALProductionCheck_H
#define AliAnalysisTaskSEDStarEMCALProductionCheck_H
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

/// \class AliAnalysisTaskSEDStarEMCALProductionCheck

#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"

class AliRDHFCutsDStartoKpipi;
class AliNormalizationCounter;

class AliAnalysisTaskSEDStarEMCALProductionCheck : public AliAnalysisTaskSE
{

public:

  AliAnalysisTaskSEDStarEMCALProductionCheck();
  AliAnalysisTaskSEDStarEMCALProductionCheck(const Char_t* name, AliRDHFCutsDStartoKpipi* cuts);
  virtual ~AliAnalysisTaskSEDStarEMCALProductionCheck();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);


  void SetAODMismatchProtection(Int_t opt = 1) {fAODProtection = opt;}

  /// Background simulation
  void     SideBandBackground(AliAODRecoCascadeHF *part, AliRDHFCutsDStartoKpipi *cuts, Int_t isSel, TList *listout, TH1F** histlist);
  void     WrongSignForDStar(AliAODRecoCascadeHF *part, AliRDHFCutsDStartoKpipi *cuts, TList *listout);
  /// histos
  void   FillSpectrum(AliAODRecoCascadeHF *part, Int_t isDStar, AliRDHFCutsDStartoKpipi *cuts, Int_t isSel, TList *listout, TH1F** histlist);
  void     DefineHistograms();
  Int_t CheckOrigin(TClonesArray* arrayMC, const AliAODMCParticle *mcPartCandidate) const;
  void CreateImpactParameterHistos();

  /// set analysis type
  void     SetAnalysisType(Int_t anaType) {fAnalysis = anaType;}
  void     PrintAnalysisType() {printf("Analysis type: %d\n(0: Heidelberg\t1: Utrecht)", fAnalysis);}

  /// set MC usage
  void     SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t   GetMC() const {return fUseMCInfo;}
/// set rare mesons
  void     SetRareSearch(Bool_t theRareOn) {fDoSearch = theRareOn;}
  Bool_t   GetRareSearch() const {return fDoSearch;}
  /// impact par study
  void SetDoImpactParameterHistos(Bool_t doImp = kTRUE) {fDoImpParDstar = doImp;}
  Bool_t GetDoImpactParameterHistos() const {return fDoImpParDstar;}

  Float_t GetTrueImpactParameterD0(const AliAODMCHeader *mcHeader, TClonesArray* arrayMC, const AliAODMCParticle *partDp) const;

  void SetDoDStarVsY(Bool_t theDStarVsY) {fDoDStarVsY = theDStarVsY;}

  void     SetUseEMCalTrigger(Bool_t bUseEMCalTrigger) {fUseEMCalTrigger = bUseEMCalTrigger;}
  Bool_t   GetUseEMCalTrigger() const {return fUseEMCalTrigger;}

  void     SetTriggerSelectionString(TString nameEMCalTrigger) {fTriggerSelectionString = nameEMCalTrigger;}
  Bool_t   GetTriggerSelectionString() const {return fTriggerSelectionString;}

  void     SetCheckEMCalAcceptance(Bool_t bCheckEMCALAcceptance) {fCheckEMCALAcceptance = bCheckEMCALAcceptance;}
  Bool_t   GetCheckEMCalAcceptance() const {return fCheckEMCALAcceptance;}

  void     SetCheckEMCalAcceptanceNumber(Int_t nCheckEMCALAcceptanceNumber) {fCheckEMCALAcceptanceNumber = nCheckEMCALAcceptanceNumber;}
  Int_t    GetCheckEMCalAcceptanceNumber() const {return fCheckEMCALAcceptanceNumber;}

  void     SetApplyEMCALClusterEventCut(Bool_t bApplyEMCALClusterEventCut) {fApplyEMCALClusterEventCut = bApplyEMCALClusterEventCut;}
  Bool_t   GetApplyEMCALClusterEventCut() const {return fApplyEMCALClusterEventCut;}

private:

  AliAnalysisTaskSEDStarEMCALProductionCheck(const AliAnalysisTaskSEDStarEMCALProductionCheck &source);
  AliAnalysisTaskSEDStarEMCALProductionCheck& operator=(const AliAnalysisTaskSEDStarEMCALProductionCheck& source);

  enum {kDzMass, kDstarMass, kDeltaMass, kptMass, ketaMass, kDzSgn, kDstarSgn, kDeltaSgn, kptSgn, ketaSgn, kDzBkg, kDstarBkg, kDeltaBkg, kptBkg, ketaBkg, kSideBandMass, kWrongSignMass};

  TH1F** fAllhist;               // Histogramlist all
  TH1F** fPIDhist;               // Histogramlist with PID
  Int_t fNPtBins;                // Number of ptbins specified in the cutfile
  Int_t  fEvents;                ///  n. of events
  Int_t  fAnalysis;    ///  0: HD;  1: UU;
  Double_t fD0Window;    ///  select width on D0Mass
  Double_t fPeakWindow;          ///  select width on DstarMass
  Bool_t fUseMCInfo;             ///  Use MC info
  Bool_t fDoSearch;              ///  Rare mesons
  TList *fOutput;                //!<!  User output
  TList *fOutputAll;             //!<!  User output2
  TList *fOutputPID;             //!<!  User output3
  TList *fOutputProductionCheck; //!<!  User output
  Int_t  fNSigma;                ///  n sigma for kaon PID
  AliRDHFCutsDStartoKpipi *fCuts; /// Cuts - sent to output slot 3
  // define the histograms
  TH1F *fCEvents;             //!<!
  TH2F *fTrueDiff2;           //!<!
  TH1F *fDeltaMassD1;         //!<!
  AliNormalizationCounter *fCounter;//!<!Counter for normalization slot 4
  Int_t fAODProtection;  /// flag to activate protection against AOD-dAOD mismatch.
  /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID
  Bool_t fDoImpParDstar;  /// imppar studies
  Int_t  fNImpParBins;   /// nunber of bins in impact parameter histos
  Float_t fLowerImpPar;  /// lower limit in impact parameter (um)
  Float_t fHigherImpPar; /// higher limit in impact parameter (um)
  Bool_t  fDoDStarVsY;   /// flag to enable D* vs y
  Bool_t fUseEMCalTrigger; /// flag to use simulated EMCal trigger in MC
  TString fTriggerSelectionString; /// Level 1 name of EMCal trigger
  Bool_t fCheckEMCALAcceptance; /// flag to perform emcal acceptance check
  Int_t fCheckEMCALAcceptanceNumber; /// selection level for emcal check
  Bool_t fApplyEMCALClusterEventCut; /// flag to cut events that do not pass EMCAL cluster check

  THnSparseF *fHistMassPtImpParTCDs[5];//!<! histograms for impact paramter studies

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEDStarEMCALProductionCheck, 2); /// class for D* spectra
  /// \endcond
};

#endif

