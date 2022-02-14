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
//
//                 Author Lennart van Doremalen
//           Utrecht University - l.v.r.vandoremalen@uu.nl
//
//     Several AliPhysics classes have been used as a basis for this code
//
///***********************************************************


/* $Id$ */

/// \class AliAnalysisTaskSEB0toDStarPi


#include <vector>
#include <TH3F.h>
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"

#include "AliAnalysisTaskSE.h"

class AliRDHFCutsB0toDStarPi;
// class AliNormalizationCounter;

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
 
  /// histos
  void     DefineHistograms();

 /// selection and reconstruction
  void     B0toDStarPiSignalTracksInMC(TClonesArray * mcTrackArray,AliAODEvent*  aodevent,TMatrix * B0toDStarPiLabelMatrix, TList *listout);
  Bool_t   D0FirstDaughterSelection(AliAODTrack* aodTrack, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header);
  Bool_t   D0SecondDaughterSelection(AliAODTrack* aodTrack, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header);
  void     DStarPionSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header);
  void     B0PionSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header);  
  void     D0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz,TClonesArray * mcTrackArray,TMatrix *  B0toDStarPiLabelMatrix, TClonesArray * D0TracksFromFriendFile, AliAODMCHeader * header);
  void     DStarAndB0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, TClonesArray * D0TracksFromFriendFile, AliAODMCHeader * header);
  Int_t    IsTrackInjected(AliAODTrack *part,AliAODMCHeader *header,TClonesArray *arrayMC);
  Bool_t   IsCandidateInjected(AliAODRecoDecayHF2Prong *part, AliAODMCHeader *header,TClonesArray *arrayMC);

  AliAODVertex* RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField, Double_t dispersion);
  void     FillFinalTrackHistograms(AliAODRecoDecayHF2Prong * trackB0, Bool_t isDesiredCandidate,TClonesArray * mcTrackArray);

  void     FillD0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType, Int_t pdgCodeMother = -1);
  void     FillDStarAndB0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType);
  Int_t    MatchCandidateToMonteCarlo(Int_t pdgabs, AliAODRecoDecayHF2Prong * candidate, TClonesArray *mcArray, TMatrix * B0toDStarPiLabelMatrix) const;

  /// set MC usage
  void     SetMC(Bool_t bUseMCInfo) {fUseMCInfo = bUseMCInfo;}
  Bool_t   GetMC() const {return fUseMCInfo;}

  Double_t DeltaInvMassDStarKpipi(AliAODRecoDecayHF2Prong *DStar) const;
  Double_t DeltaInvMassB0Kpipipi(AliAODRecoDecayHF2Prong *Bzero) const;

  void     SetQuickSignalAnalysis(Int_t value){fQuickSignalAnalysis = value;}
  void     SetGetCutInfo(Bool_t value){fGetCutInfo = value;}

  void     SetShowMask(Bool_t bShowMask) {fShowMask = bShowMask;}
  Bool_t   GetShowMask() const {return fShowMask;}

  void     SetShowRejection(Bool_t bShowRejection) {fShowRejection = bShowRejection;}
  Bool_t   GetShowRejection() const {return fShowRejection;}

  void     SetUse3DHistograms(Bool_t bUse3DHistograms) {fUse3DHistograms = bUse3DHistograms;}
  Bool_t   GetUse3DHistograms() const {return fUse3DHistograms;}

  void     SetUpgradeSetting(Int_t nUpgradeSetting) {fUpgradeSetting = nUpgradeSetting;}
  Int_t    GetUpgradeSetting() const {return fUpgradeSetting;}

  void     SetHistMassWindow(Double_t value) {fHistMassWindow = value;}
  Double_t GetHistMassWindow() const {return fHistMassWindow;}

  void     SetDegreePerRotation(Int_t value) {fDegreePerRotation = value;}
  Int_t    GetDegreePerRotation() const {return fDegreePerRotation;}

  void     SetNumberOfRotations(Int_t value) {fNumberOfRotations = value;}
  Int_t    GetNumberOfRotations() const {return fNumberOfRotations;}

  void     SetCheckBackground(Bool_t value) {fCheckBackground = value;}
  Bool_t   GetCheckBackground() const {return fCheckBackground;}

  void     SetCheckInjected(Bool_t value) {fCheckInjected = value;}
  Bool_t   GetCheckInjected() const {return fCheckInjected;}

  void     SetRemoveInjected(Int_t value) {fRemoveInjected = value;}
  Int_t    GetRemoveInjected() const {return fRemoveInjected;}

 private:
  
  AliAnalysisTaskSEB0toDStarPi(const AliAnalysisTaskSEB0toDStarPi &source);
  AliAnalysisTaskSEB0toDStarPi& operator=(const AliAnalysisTaskSEB0toDStarPi& source); 
  
  Int_t  fEvents;                            // 
  Bool_t fUseMCInfo;                         //  Use MC info
  Bool_t fShowMask;                          //
  Bool_t fShowRejection;                     //
  Int_t  fQuickSignalAnalysis;               //
  Bool_t fGetCutInfo;                        //
  Bool_t fUse3DHistograms;                   //
  Int_t  fUpgradeSetting;                    //
  Double_t fHistMassWindow;                  //  
  Int_t  fDegreePerRotation;                 //
  Int_t  fNumberOfRotations;                 //
  Bool_t fCheckBackground;                   //
  Bool_t fCheckInjected;                     //
  Int_t  fRemoveInjected;                    //

  TList *fOutput;                            //!<!  User output
  TList *fListCuts;                          //!<!  User output
  TList *fOutputD0FirstDaughter;             //!<!  User output
  TList *fOutputD0SecondDaughter;            //!<!  User output
  TList *fOutputDStarPion;                   //!<!  User output
  TList *fOutputB0Pion;                      //!<!  User output
  TList *fOutputD0;                          //!<!  User output
  TList *fOutputDStar;                       //!<!  User output
  TList *fOutputB0;                          //!<!  User output
  TList *fOutputD0_D0Pt;                     //!<!  User output
  TList *fOutputD0_DStarPt;                  //!<!  User output
  TList *fOutputDStar_DStarPt;               //!<!  User output
  TList *fOutputB0MC;                        //!<!  User output 

  AliRDHFCutsB0toDStarPi *fCuts;             // Cuts - sent to output 

  TH1F *fCEvents;                            //!<!

  std::vector<Int_t> * fDStarPionTracks;     //!
  std::vector<Int_t> * fB0PionTracks;        //!
  std::vector<Int_t> * fD0Tracks;            //!
    
  Int_t fnPtBins;                            //!
  Int_t fnPtBinLimits;                       //!
  Float_t * fPtBinLimits;                    //! [fnPtBinLimits]
  Int_t fnPtBinsD0forD0ptbin;                //!
  Int_t fnPtBinsD0forD0ptbinLimits;          //!
  Float_t * fPtBinLimitsD0forD0ptbin;        //! [fnPtBinsD0forD0ptbinLimits]
  Int_t fnPtBinsD0forDStarptbin;             //!
  Int_t fnPtBinsD0forDStarptbinLimits;       //!
  Float_t * fPtBinLimitsD0forDStarptbin;     //! [fnPtBinsD0forDStarptbinLimits]
  Int_t fnPtBinsDStarforDStarptbin;          //!
  Int_t fnPtBinsDStarforDStarptbinLimits;    //!
  Float_t * fPtBinLimitsDStarforDStarptbin;  //! [fnPtBinsDStarforDStarptbinLimits]


  TH1F* fDaughterHistogramArray[4][6][15];   //!
  TH2F* fDaughterHistogramArray2D[4][6];     //!
  TH1F* fDaughterHistogramArrayExtra[4][6];  //!
  TH1F* fMotherHistogramArray[6][500][60];    //!
  TH2F* fMotherHistogramArray2D[6][500][60];  //!
  TH1F* fMotherHistogramArrayExtra[7][10];   //!
  TH3F* fMotherHistogramArray3D[6][500][60];  //!

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEB0toDStarPi,6);  // class for B0 spectra
  /// \endcond
};

#endif

