#ifndef AliAnalysisTaskSEB0toDStarPi_H
#define AliAnalysisTaskSEB0toDStarPi_H
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

/// \class AliAnalysisTaskSEB0toDStarPi

#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <vector>
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODVertex.h"

#include "AliAnalysisTaskSE.h"

class AliRDHFCutsB0toDStarPi;
class AliNormalizationCounter;

class AliAnalysisTaskSEB0toDStarPi : public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskSEB0toDStarPi();
  AliAnalysisTaskSEB0toDStarPi(const Char_t* name,AliRDHFCutsB0toDStarPi* cuts);
  virtual ~AliAnalysisTaskSEB0toDStarPi();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  // histos
  void     DefineHistograms();

  //selection and reconstruction
  void     B0toDStarPiSignalTracksInMC(TClonesArray * mcTrackArray,AliAODEvent*  aodevent,TMatrix * B0toDStarPiLabelMatrix, TList *listout);
  void     DStarPionSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix);
  void     B0PionSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix);  
  void     D0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz,TClonesArray * mcTrackArray,TMatrix *  B0toDStarPiLabelMatrix, TClonesArray * D0TracksFromFriendFile);
  void     DStarAndB0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, TClonesArray * D0TracksFromFriendFile);

  AliAODVertex* RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField, Double_t dispersion);
  void     FillFinalTrackHistograms(AliAODRecoDecayHF2Prong * trackB0, Bool_t isDesiredCandidate,TClonesArray * mcTrackArray);

  void     FillD0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType, Int_t pdgCodeMother = -1);
  void     FillDStarAndB0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType);
  Int_t    MatchCandidateToMonteCarlo(Int_t pdgabs, AliAODRecoDecayHF2Prong * candidate, TClonesArray *mcArray, TMatrix * B0toDStarPiLabelMatrix) const;

  // set MC usage
  void     SetMC(Bool_t bUseMCInfo) {fUseMCInfo = bUseMCInfo;}
  Bool_t   GetMC() const {return fUseMCInfo;}

  Double_t DeltaInvMassDStarKpipi(AliAODRecoDecayHF2Prong * DStar) const;
  Double_t DeltaInvMassB0Kpipipi(AliAODRecoDecayHF2Prong * B0) const;

  void     SetQuickSignalAnalysis(Bool_t value){fQuickSignalAnalysis = value;}
  void     SetGetCutInfo(Bool_t value){fGetCutInfo = value;}

  void     SetShowMask(Bool_t bShowMask) {fShowMask = bShowMask;}
  Bool_t   GetShowMask() const {return fShowMask;}

  void     SetShowRejection(Bool_t bShowRejection) {fShowRejection = bShowRejection;}
  Bool_t   GetShowRejection() const {return fShowRejection;}

 private:
  
  AliAnalysisTaskSEB0toDStarPi(const AliAnalysisTaskSEB0toDStarPi &source);
  AliAnalysisTaskSEB0toDStarPi& operator=(const AliAnalysisTaskSEB0toDStarPi& source); 
  
  Int_t  fEvents;                ///  n. of events
  Bool_t fUseMCInfo;             ///  Use MC info
  Bool_t fShowMask;
  Bool_t fShowRejection;

  TList *fOutput;                //!<!  User output
  TList *fOutputD0Pion;
  TList *fOutputD0Kaon;
  TList *fOutputDStarPion;
  TList *fOutputB0Pion;
  TList *fOutputD0;
  TList *fOutputDStar;
  TList *fOutputB0;
  TList *fOutputD0_D0Pt;
  TList *fOutputD0_DStarPt;
  TList *fOutputDStar_DStarPt;
  TList *fOutputB0MC;             //!<!  User output 7


  AliRDHFCutsB0toDStarPi *fCuts; /// Cuts - sent to output slot 3

  TH1F *fCEvents;             //!<!

  AliNormalizationCounter *fCounter;//!<!Counter for normalization slot 4

  std::vector<Int_t> * fDStarPionTracks;
  std::vector<Int_t> * fB0PionTracks;
  std::vector<Int_t> * fD0Tracks;
    
  Int_t fnPtBins;
  Int_t fnPtBinLimits;
  Float_t * fPtBinLimits;                    //[fnPtBinLimits]
  Int_t fnPtBinsD0forD0ptbin;
  Int_t fnPtBinsD0forD0ptbinLimits;
  Float_t * fPtBinLimitsD0forD0ptbin;        //[fnPtBinsD0forD0ptbinLimits]
  Int_t fnPtBinsD0forDStarptbin;
  Int_t fnPtBinsD0forDStarptbinLimits;
  Float_t * fPtBinLimitsD0forDStarptbin;     //[fnPtBinsD0forDStarptbinLimits]
  Int_t fnPtBinsDStarforDStarptbin;
  Int_t fnPtBinsDStarforDStarptbinLimits;
  Float_t * fPtBinLimitsDStarforDStarptbin;  //[fnPtBinsDStarforDStarptbinLimits]

  Int_t* fDaughterHistogramArray[4][5][15];
  Int_t* fDaughterHistogramArray2D[4][6];
  Int_t* fMotherHistogramArray[6][50][46];
  Int_t* fMotherHistogramArray2D[6][16][60];
  Int_t* fMotherHistogramArrayExtra[7][10];


  Bool_t fQuickSignalAnalysis;
  Bool_t fGetCutInfo;

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEB0toDStarPi,2); /// class for B0 spectra
  /// \endcond
};

#endif

