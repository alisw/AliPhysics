#ifndef AliAnalysisTaskSEB0toDPi_H
#define AliAnalysisTaskSEB0toDPi_H
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

/// \class AliAnalysisTaskSEB0toDPi

#include <vector>
#include <TH3F.h>
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"

#include "AliAnalysisTaskSE.h"

class AliRDHFCutsB0toDPi;

class AliAnalysisTaskSEB0toDPi : public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskSEB0toDPi();
  AliAnalysisTaskSEB0toDPi(const Char_t* name,AliRDHFCutsB0toDPi* cuts);
  virtual ~AliAnalysisTaskSEB0toDPi();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  /// histos
  void     DefineHistograms();

  //selection and reconstruction
  void     B0toDPiSignalTracksInMC(TClonesArray * mcTrackArray,AliAODEvent*  aodevent,TMatrix * B0toDPlusPiLabelMatrix, TList *listout);
  void     DaughterSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDPiLabelMatrix, AliAODMCHeader * header, Int_t daughterType, std::vector<Int_t> * trackVector);
  Bool_t   DPlusDaughterSelection(AliAODTrack*  aodTrack, AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDPiLabelMatrix, AliAODMCHeader * header, Int_t daughterType);

  void     DPlusSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz,TClonesArray * mcTrackArray, TMatrix * B0toDPlusPiLabelMatrix, TClonesArray * DPlusTrackArray, AliAODMCHeader * header);
  void     B0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDPlusPiLabelMatrix, TClonesArray * DPlusTrackArray, AliAODMCHeader * header);
  Int_t    IsTrackInjected(AliAODTrack *part,AliAODMCHeader *header,TClonesArray *arrayMC);
  Bool_t   IsCandidateInjected(AliAODRecoDecayHF2Prong *part, AliAODMCHeader *header,TClonesArray *arrayMC);
  void     CutOptimizationLoop(Int_t variable, Int_t nVariables, Int_t nCuts, Int_t ptBin, Int_t fillNumber, Bool_t isDesiredCandidate, Int_t nSigmaBin);
  void     CutOptimizationVariableValues(AliAODRecoDecayHF2Prong * candidateB0, AliAODEvent*  aod);

  AliAODVertex* RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField, Double_t dispersion, Bool_t optUseFitter, Bool_t optPropagate, Bool_t optUseDiamondConstraint, Int_t nprongs);
  void     FillDaughterHistograms(AliAODTrack* daughterTrack, Int_t daughterType, Int_t histType, Double_t ptB0, AliAODVertex *primaryVertex, Double_t bz, Bool_t isDesiredCandidate, TClonesArray * mcTrackArray);
  void     FillFinalTrackHistograms(AliAODRecoDecayHF2Prong * selectedB0, AliAODVertex *primaryVertex, Double_t bz, Bool_t isDesiredCandidate,TClonesArray * mcTrackArray);

  void     FillDPlusHistograms(AliAODRecoDecayHF3Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType, Int_t pdgCodeMother = -1);
  void     FillB0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType);
  Int_t    MatchCandidateToMonteCarlo(Int_t pdgabs, AliAODRecoDecayHF * candidate, TClonesArray *mcArray, TMatrix * B0toDStarPiLabelMatrix, Bool_t bCheckLabel = kFALSE) const;

  void     TwoTrackCombinationInfo(AliExternalTrackParam * firstTrack, AliExternalTrackParam * secondTrack, AliAODVertex * primaryVertex, Double_t bz, Bool_t isDesiredCandidate, Int_t histogramNumber, UInt_t prongs[2]);
  void     ThreeTrackCombinationInfo(AliExternalTrackParam * firstTrack, AliExternalTrackParam * secondTrack, AliExternalTrackParam * thirdTrack, AliAODVertex * primaryVertex, Double_t bz, Bool_t isDesiredCandidate, Int_t histogramNumber, UInt_t prongs[3]);

  /// set MC usage
  void     SetMC(Bool_t bUseMCInfo) {fUseMCInfo = bUseMCInfo;}
  Bool_t   GetMC() const {return fUseMCInfo;}

  Double_t DeltaInvMassB0Kpipipi(AliAODRecoDecayHF2Prong *Bzero) const;

  void     SetQuickSignalAnalysis(Int_t value){fQuickSignalAnalysis = value;}

  void     SetShowMask(Bool_t value) {fShowMask = value;}
  Bool_t   GetShowMask() const {return fShowMask;}

  void     SetShowRejection(Bool_t value) {fShowRejection = value;}
  Bool_t   GetShowRejection() const {return fShowRejection;}

  void     SetHistMassWindow(Double_t value) {fHistMassWindow = value;}
  Double_t GetHistMassWindow() const {return fHistMassWindow;}

  void     SetDegreePerRotation(Int_t value) {fDegreePerRotation = value;}
  Int_t    GetDegreePerRotation() const {return fDegreePerRotation;}

  void     SetNumberOfRotations(Int_t value) {fNumberOfRotations = value;}
  Int_t    GetNumberOfRotations() const {return fNumberOfRotations;}

  void     SetPerformCutOptimization(Bool_t value) {fPerformCutOptimization = value;}
  Bool_t   GetPerformCutOptimization() const {return fPerformCutOptimization;}

  void     SetRemoveInjected(Bool_t value) {fRemoveInjected = value;}
  Bool_t   GetRemoveInjected() const {return fRemoveInjected;}

  void     SetUseSideBands(Bool_t value) {fUseSideBands = value;}
  Bool_t   GetUseSideBands() const {return fUseSideBands;}

  void     SetSideBandLow(Double_t value) {fSideBandLow = value;}
  Double_t GetSideBandLow() const {return fSideBandLow;}

  void     SetSideBandHigh(Double_t value) {fSideBandHigh = value;}
  Double_t GetSideBandHigh() const {return fSideBandHigh;}

  void     SetSaveTRHists(Bool_t value) {fSaveTRHists = value;}
  Bool_t   GetSaveTRHists() const {return fSaveTRHists;}


 private:
  
  AliAnalysisTaskSEB0toDPi(const AliAnalysisTaskSEB0toDPi &source);
  AliAnalysisTaskSEB0toDPi& operator=(const AliAnalysisTaskSEB0toDPi& source); 
  
  Int_t     fEvents;                             // 
  Bool_t    fUseMCInfo;                          //  Use MC info
  Bool_t    fShowMask;                           //
  Bool_t    fShowRejection;                      //
  Int_t     fQuickSignalAnalysis;                //
  Double_t  fHistMassWindow;                     //  
  Int_t     fDegreePerRotation;                  //
  Int_t     fNumberOfRotations;                  //
  Bool_t    fPerformCutOptimization;             //
  Bool_t    fRemoveInjected;                     //
  Bool_t    fUseSideBands;                       //
  Double_t  fSideBandLow;                        //
  Double_t  fSideBandHigh;                       //
  Bool_t    fSaveTRHists;                        //

  TList *fOutput;                                //!<!  User output
  TList *fListCuts;                              //!<!  User output  
  TList *fOutputB0Results;                       //!<!  User output
  TList *fOutputDPlusKaon;                       //!<!  User output
  TList *fOutputDPlusPions;                      //!<!  User output
  TList *fOutputB0Pion;                          //!<!  User output
  TList *fOutputDPlus;                           //!<!  User output
  TList *fOutputB0;                              //!<!  User output
  TList *fOutputDPlus_DPlusPt;                   //!<!  User output

  AliRDHFCutsB0toDPi *fCuts;                     // Cuts - sent to output
  
  TH1F *fCEvents;                                //!<!
  TH1F *fCTrigger;                               //!<!
  TH1F *fCRejected;                              //!<!
  TH1F *fTriggerClassesCorrelated;               //!<!
  TH1F *fTriggerSubClasses;                      //!<!

  std::vector<Int_t> * fB0PionTracks;            //!
  std::vector<Int_t> * fDPlusTracks;             //!
  // std::vector<Int_t> * fDPlusPionTracks;         //!

  Int_t fnPtBins;                                //!
  Int_t fnPtBinLimits;                           //!
  Float_t * fPtBinLimits;                        //! [fnPtBinLimits]
  Int_t fnPtBinsDPlusforDPlusptbin;              //!
  Int_t fnPtBinsDPlusforDPlusptbinLimits;        //!
  Float_t * fPtBinLimitsDPlusforDPlusptbin;      //! [fnPtBinsDPlusforDPlusptbinLimits]
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
  ClassDef(AliAnalysisTaskSEB0toDPi,1); /// class for B0 spectra
  /// \endcond
};

#endif

