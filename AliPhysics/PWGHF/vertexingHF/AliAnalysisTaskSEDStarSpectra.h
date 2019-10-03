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

/// \class AliAnalysisTaskSEDStarSpectra

#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"

class AliRDHFCutsDStartoKpipi;
class AliNormalizationCounter;

class AliAnalysisTaskSEDStarSpectra : public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskSEDStarSpectra();
  AliAnalysisTaskSEDStarSpectra(const Char_t* name,AliRDHFCutsDStartoKpipi* cuts);
  virtual ~AliAnalysisTaskSEDStarSpectra();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);


  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
  
  /// Background simulation
  void     SideBandBackground(AliAODRecoCascadeHF *part, AliRDHFCutsDStartoKpipi *cuts, Int_t isSel, TList *listout, TH1F** histlist);
  void     WrongSignForDStar(AliAODRecoCascadeHF *part, AliRDHFCutsDStartoKpipi *cuts, TList *listout);
    /// histos
  void   FillSpectrum(AliAODRecoCascadeHF *part, Int_t isDStar, AliRDHFCutsDStartoKpipi *cuts, Int_t isSel, TList *listout,TH1F** histlist);
  void     DefineHistograms();
  Int_t CheckOrigin(TClonesArray* arrayMC, const AliAODMCParticle *mcPartCandidate) const;
  void CreateImpactParameterHistos();

  /// set analysis type
  void     SetAnalysisType(Int_t anaType) {fAnalysis = anaType;}
  void     PrintAnalysisType() {printf("Analysis type: %d\n(0: Heidelberg\t1: Utrecht)",fAnalysis);}

  /// set MC usage
  void     SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t   GetMC() const {return fUseMCInfo;}
 /// set rare mesons
  void     SetRareSearch(Bool_t theRareOn) {fDoSearch = theRareOn;}
  Bool_t   GetRareSearch() const {return fDoSearch;}
  /// impact par study
  void SetDoImpactParameterHistos(Bool_t doImp=kTRUE){fDoImpParDstar=doImp;}
  Bool_t GetDoImpactParameterHistos() const {return fDoImpParDstar;}

  Float_t GetTrueImpactParameterD0(const AliAODMCHeader *mcHeader, TClonesArray* arrayMC, const AliAODMCParticle *partDp) const;

  void SetDoDStarVsY(Bool_t theDStarVsY) {fDoDStarVsY = theDStarVsY;}

 private:
  
  AliAnalysisTaskSEDStarSpectra(const AliAnalysisTaskSEDStarSpectra &source);
  AliAnalysisTaskSEDStarSpectra& operator=(const AliAnalysisTaskSEDStarSpectra& source); 

  enum{kDzMass, kDstarMass, kDeltaMass,kptMass, ketaMass,kDzSgn, kDstarSgn, kDeltaSgn,kptSgn, ketaSgn,kDzBkg, kDstarBkg, kDeltaBkg,kptBkg, ketaBkg, kSideBandMass, kWrongSignMass};

  TH1F** fAllhist;               // Histogramlist all
  TH1F** fPIDhist;               // Histogramlist with PID
  Int_t fNPtBins;                // Number of ptbins specified in the cutfile
  Int_t  fEvents;                ///  n. of events
  Int_t  fAnalysis;		 ///  0: HD;	1: UU;
  Double_t fD0Window;		 ///  select width on D0Mass
  Double_t fPeakWindow;          ///  select width on DstarMass
  Bool_t fUseMCInfo;             ///  Use MC info
  Bool_t fDoSearch;              ///  Rare mesons
  TList *fOutput;                //!<!  User output
  TList *fOutputAll;             //!<!  User output2
  TList *fOutputPID;             //!<!  User output3
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

  THnSparseF *fHistMassPtImpParTCDs[5];//!<! histograms for impact paramter studies

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEDStarSpectra,10); /// class for D* spectra
  /// \endcond
};

#endif

