#ifndef AliAnalysisTaskSEBPlustoD0Pi_H
#define AliAnalysisTaskSEBPlustoD0Pi_H
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
//
//                 Author Lennart van Doremalen
//           Utrecht University - l.v.r.vandoremalen@uu.nl
//
//     Several AliPhysics classes have been used as a basis for this code
//
///***********************************************************


/* $Id$ */

/// \class AliAnalysisTaskSEBPlustoD0Pi

#include <vector>
#include <TH3F.h>
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"

#include "AliAnalysisTaskSE.h"

class AliRDHFCutsBPlustoD0Pi;

class AliAnalysisTaskSEBPlustoD0Pi : public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskSEBPlustoD0Pi();
  AliAnalysisTaskSEBPlustoD0Pi(const Char_t* name,AliRDHFCutsBPlustoD0Pi* cuts);
  virtual ~AliAnalysisTaskSEBPlustoD0Pi();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  /// histos
  void     DefineHistograms();

  //selection and reconstruction
  void     BPlustoD0PiSignalTracksInMC(TClonesArray * mcTrackArray,AliAODEvent*  aodevent,TMatrix * BPlustoD0PiLabelMatrix, TList *listout);
  Bool_t   D0FirstDaughterSelection(AliAODTrack* aodTrack, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header, AliAODEvent* aodEvent);
  Bool_t   D0SecondDaughterSelection(AliAODTrack* aodTrack, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header, AliAODEvent* aodEvent);
  void     BPlusPionSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * BPlustoD0PiLabelMatrix, AliAODMCHeader * header);

  void     D0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz,TClonesArray * mcTrackArray, TMatrix * BPlustoD0PiLabelMatrix, TClonesArray * D0TracksFromFriendFile, AliAODMCHeader * header);
  void     BPlusSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * BPlustoD0PiLabelMatrix, TClonesArray * D0TracksFromFriendFile, AliAODMCHeader * header);
  Int_t    IsTrackInjected(AliAODTrack *part,AliAODMCHeader *header,TClonesArray *arrayMC);
  Bool_t   IsCandidateInjected(AliAODRecoDecayHF2Prong *part, AliAODMCHeader *header,TClonesArray *arrayMC);
  void     CutOptimizationLoop(Int_t variable, Int_t nVariables, Int_t nCuts, Int_t ptBin, Int_t fillNumber, Bool_t isDesiredCandidate, Int_t nSigmaBin);
  void     CutOptimizationVariableValues(AliAODRecoDecayHF2Prong * candidateBPlus, AliAODEvent*  aod);

  AliAODVertex* RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField, Double_t dispersion, Bool_t optUseFitter, Bool_t optPropagate, Bool_t optUseDiamondConstraint);
  void     FillFinalTrackHistograms(AliAODRecoDecayHF2Prong * selectedBPlus, AliAODVertex *primaryVertex, Double_t bz, Bool_t isDesiredCandidate,TClonesArray * mcTrackArray);

  void     FillD0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType, Int_t pdgCodeMother = -1);
  void     FillBPlusHistograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType);
  Int_t    MatchCandidateToMonteCarlo(Int_t pdgabs, AliAODRecoDecayHF2Prong * candidate, TClonesArray *mcArray, TMatrix * B0toDStarPiLabelMatrix, Bool_t bCheckLabel = kFALSE) const;

  void     TwoTrackCombinationInfo(AliExternalTrackParam * firstTrack, AliExternalTrackParam * secondTrack, AliAODVertex * primaryVertex, Double_t bz, Bool_t isDesiredCandidate, Int_t histogramNumber, UInt_t prongs[2]);
  void     ThreeTrackCombinationInfo(AliExternalTrackParam * firstTrack, AliExternalTrackParam * secondTrack, AliExternalTrackParam * thirdTrack, AliAODVertex * primaryVertex, Double_t bz, Bool_t isDesiredCandidate, Int_t histogramNumber, UInt_t prongs[3]);

  /// set MC usage
  void     SetMC(Bool_t bUseMCInfo) {fUseMCInfo = bUseMCInfo;}
  Bool_t   GetMC() const {return fUseMCInfo;}

  Double_t DeltaInvMassBPlusKpipi(AliAODRecoDecayHF2Prong * BPlus) const;

  void     SetQuickSignalAnalysis(Int_t value){fQuickSignalAnalysis = value;}

  void     SetShowMask(Bool_t bShowMask) {fShowMask = bShowMask;}
  Bool_t   GetShowMask() const {return fShowMask;}

  void     SetShowRejection(Bool_t bShowRejection) {fShowRejection = bShowRejection;}
  Bool_t   GetShowRejection() const {return fShowRejection;}

  void     SetHistMassWindow(Double_t value) {fHistMassWindow = value;}
  Double_t GetHistMassWindow() const {return fHistMassWindow;}

  void     SetDegreePerRotation(Int_t value) {fDegreePerRotation = value;}
  Int_t    GetDegreePerRotation() const {return fDegreePerRotation;}

  void     SetNumberOfRotations(Int_t value) {fNumberOfRotations = value;}
  Int_t    GetNumberOfRotations() const {return fNumberOfRotations;}

  void     SetPerformCutOptimization(Bool_t bPerformCutOptimization) {fPerformCutOptimization = bPerformCutOptimization;}
  Bool_t   GetPerformCutOptimization() const {return fPerformCutOptimization;}

  void     SetRemoveInjected(Bool_t bRemoveInjected) {fRemoveInjected = bRemoveInjected;}
  Bool_t   GetRemoveInjected() const {return fRemoveInjected;}

 private:
  
  AliAnalysisTaskSEBPlustoD0Pi(const AliAnalysisTaskSEBPlustoD0Pi &source);
  AliAnalysisTaskSEBPlustoD0Pi& operator=(const AliAnalysisTaskSEBPlustoD0Pi& source); 
  
  Int_t  fEvents;                                // 
  Bool_t fUseMCInfo;                             //  Use MC info
  Bool_t fShowMask;                              //
  Bool_t fShowRejection;                         //
  Int_t  fQuickSignalAnalysis;                   //
  Double_t fHistMassWindow;                      //  
  Int_t  fDegreePerRotation;                     //
  Int_t  fNumberOfRotations;                     //
  Bool_t fPerformCutOptimization;                //
  Bool_t fRemoveInjected;                        //

  TList *fOutput;                                //!<!  User output
  TList *fListCuts;                              //!<!  User output  
  TList *fOutputBPlusResults;                    //!<!  User output
  TList *fOutputD0FirstDaughter;                 //!<!  User output
  TList *fOutputD0SecondDaughter;                //!<!  User output
  TList *fOutputBPlusPion;                       //!<!  User output
  TList *fOutputD0;                              //!<!  User output
  TList *fOutputBPlus;                           //!<!  User output
  TList *fOutputD0_D0Pt;                         //!<!  User output

  AliRDHFCutsBPlustoD0Pi *fCuts;                 // Cuts - sent to output
  
  TH1F *fCEvents;                                //!<!
  TH1F *fCTrigger;                               //!<!
  TH1F *fCRejected;                              //!<!

  std::vector<Int_t> * fBPlusPionTracks;         //!
  std::vector<Int_t> * fD0Tracks;                //!

  Int_t fnPtBins;                                //!
  Int_t fnPtBinLimits;                           //!
  Float_t * fPtBinLimits;                        //! [fnPtBinLimits]
  Int_t fnPtBinsD0forD0ptbin;                    //!
  Int_t fnPtBinsD0forD0ptbinLimits;              //!
  Float_t * fPtBinLimitsD0forD0ptbin;            //! [fnPtBinsD0forD0ptbinLimits]
  Float_t fCutVariableValueArray[99];            //!


  TH1F* fDaughterHistogramArray[4][6][15];       //!
  TH2F* fDaughterHistogramArray2D[4][6];         //!
  TH1F* fDaughterHistogramArrayExtra[4][6];      //!
  TH1F* fMotherHistogramArray[6][99][60];        //!
  TH2F* fMotherHistogramArray2D[6][99][60];      //!
  TH1F* fMotherHistogramArrayExtra[7][10];       //!
  TH1F* fResultsHistogramArray[20][99];          //!
  TH2F* fResultsHistogramArray2D[20][99];        //!

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEBPlustoD0Pi,3); /// class for BPlus spectra
  /// \endcond
};

#endif

