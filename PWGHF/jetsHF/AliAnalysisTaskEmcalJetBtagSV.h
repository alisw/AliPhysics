#ifndef ALIANALYSISTASKEMCALJETBTAGSV_H
#define ALIANALYSISTASKEMCALJETBTAGSV_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* Class AliAnalysisTaskEmcalJetBtagSV:                                    *
 * AliAnalysisTaskSE for the extraction of the B-jet spectrum        *
 * using the properties of secondary verices as tagging observables. */

/*
 Mailto: andrea.rossi@cern.ch, elena.bruna@to.infn.it,
 svallero@to.infn.it, s.lapointe@cern.ch
 ycorrale@cern.ch
 */

//--Root--
class TH1F;
class TList;
class TProfile;

//--AliRoot--
#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisUtils;
class AliParticleContainer;

//--AliHFJetsClass--
#include "AliHFJetsUtils.h"
#include "AliRDHFJetsCutsVertex.h"
#include "AliHFJetsTaggingVertex.h"

class AliHFJetsContainerVertex;

//-------------------------------------------------------------------------------------

class AliAnalysisTaskEmcalJetBtagSV : public AliAnalysisTaskEmcalJet {
  
public:
  
  AliAnalysisTaskEmcalJetBtagSV();
  
  AliAnalysisTaskEmcalJetBtagSV(const char *name);
  
  virtual ~AliAnalysisTaskEmcalJetBtagSV();
  
  // Implementation of interface methods
  virtual void    Init();
  
  virtual void    LocalInit() {Init();}
  
  virtual void    UserCreateOutputObjects();
  
  virtual void    UserExec(Option_t *option);
  
  virtual Bool_t  UserNotify();
  
  virtual void    Terminate(Option_t *option);
  
  // Setters/Getters
  void SetCorrectionMode(Bool_t mode ) {fCorrMode = mode;}
  
  void SetDoBkgRejection(Bool_t dorej) {fDoBkgRej = dorej;}
  
  void SetDoFillSecVtxQA(Bool_t doqa ) {fDoQAVtx  = doqa;}
  
  void SetDoFillV0Trks(Bool_t doV0) {fDoFillV0Trks = doV0;}
  
  void SetRecJetsBranch(const char *branch) {fRecJetsBranch = branch;}
  
  void SetGenJetsBranch(const char *branch) {fGenJetsBranch = branch;}
  
  void SetGenNamePattern(const char *patt)  {fGenNamePattern = patt; }
  
  void SetCheckMCCrossSection(Bool_t check) {fCheckMCCrossSection = check;}
  
  void SetUseWeightOn() {fUseWeight = kTRUE;}
  
  void SetJetTaggingRadius(Double_t tagradius) {fTaggingRadius = tagradius;};
  
  void SetTagger(AliHFJetsTaggingVertex *tagger)
  {
    if (fTagger) delete fTagger;
    fTagger = (AliHFJetsTaggingVertex *)tagger->Clone("fTagger");
  }
  
  void SetCuts(AliRDHFJetsCuts *cuts)
  {
    if (fCutsHFjets) delete fCutsHFjets;
    fCutsHFjets = (AliRDHFJetsCuts *)cuts->Clone("fCutsHFjets");
  }
  
  void SetJetContName(char *name)        { fJetContName = name; }
  
  void SetTrkContName(char *name)        { fTrkContName = name; }
  
  void SetRhoTaskName(char *name)        { fRhoTaskName = name; }
  
  void SetMCJetContName(char *name)      { fMCJetContName = name; }
  
  void SetMCTrkContName(char *name)      { fMCTrkContName = name; }
  
  void SetMCRhoTaskName(char *name)      { fMCRhoTaskName = name; }
  
  void SetAnalysisUtils(AliAnalysisUtils *utils) { fAnalysisUtils = utils; }
  
  void  SetDebugLevel(Int_t level) { fDebug = (AliLog::EType_t)level; }
  
  Int_t GetDebugLevel()            { return fDebug;  }
  
  virtual void     GetPythiaCrossSection();
  
  static  Bool_t   GetPythiaInfoFromFile(TString currFile,
                                         Float_t &xsec,
                                         Float_t &trials);
  
  AliAnalysisUtils *GetAnalysisUtils()           { return fAnalysisUtils; }
  
protected:
  
  // copy constructo not implemented yet
  AliAnalysisTaskEmcalJetBtagSV           (const AliAnalysisTaskEmcalJetBtagSV &);
  
  // assignment operator not implemented yet
  AliAnalysisTaskEmcalJetBtagSV &operator=(const AliAnalysisTaskEmcalJetBtagSV &);
  
  void            AnalyseCorrectionsMode();
  
  void            AnalyseDataMode();
  
  void GetFlavour2Methods(AliEmcalJet *jet,
                          Double_t (&partonnat)[2],
                          Double_t (&ptpart)[2],
                          Double_t radius);
  
  // temporary
  Bool_t    GetArrays();
  Bool_t    FillMapOfV0gTrkIDs();
  Bool_t    FillVecOfV0gTrkIDs(std::vector<Int_t> &vctrTrkIDs);
  Bool_t    IsAODtrkBelongToV0(std::vector<Int_t> &vctrTrkIDs, Int_t trkID);
  
  Double_t  GetExternalRho(Bool_t isMC = kFALSE);
  
private:
  
  Bool_t      fCorrMode;             // enable correction or data modes
  Bool_t      fDoBkgRej;
  Bool_t      fDoQAVtx;              // enable output of qa on secondary vertex
  Bool_t      fDoFillV0Trks;         // enable V0 checks
  Bool_t      fUseTriggerData;       // use emacal trigger
  
  const char *fRecJetsBranch;        // name of the AOD REC-jets branch
  const char *fGenJetsBranch;        // name of the AOD GEN-jets branch
  
  TString     fGenNamePattern;
  
  TString     fJetContName;          //  Name of the found jet array
  TString     fTrkContName;          //  Name of the found track array
  TString     fRhoTaskName;          //  Name of the rho task
  TString     fMCJetContName;        //  Name of the found jet array
  TString     fMCTrkContName;        //  Name of the found track array
  TString     fMCRhoTaskName;        //  Name of the rho task
  
  Double_t    fTaggingRadius;        // radius used in tagging the jet flavour
  
  //
  // MC weights
  //
  Double_t    fMCWeight;             ///<  pT-hard bin MC weight. It is used only internally.
  Double_t    fMCXsec;
  Double_t    fMCAvgTrials;
  
  TString     fCurrFileName;         ///<  Current file path name.
  
  Bool_t      fCheckMCCrossSection;  ///<  Retrieve from the pyxsec.root file the cross section, only if requested.
  Bool_t      fSkipWeightInfo;
  Bool_t      fUseWeight;
  
  TList                      *fOutputList;       // list of output objects
  
  // AliHFJetsContainerVertex
  AliHFJetsContainerVertex   *fhJetVtxSim;       //!<! properties of vertices within the jet MC
  AliHFJetsContainerVertex   *fhJetVtxData;      //!<! properties of vertices within the jet Data
  AliHFJetsContainerVertex   *fhQaVtx;           //!<! vertices properties
  
  TProfile                   *fhXsec;            //!<! Cross section in PYTHIA.
  TH1F                       *fhEntries;         //!<!
  TH1F                       *fhTrials;          //!<! Number of event trials in PYTHIA.
  
  AliVEvent                  *fEvent;            //! Input event
  AliAODMCHeader             *fMCHeader;         //! Input MC header
  
  AliHFJetsTaggingVertex     *fTagger;           // Jet Tagging object
  
  AliRDHFJetsCuts            *fCutsHFjets;       //  specific algo jet cut object
  AliAnalysisUtils           *fAnalysisUtils;    //  points to class with common analysis utilities
  AliParticleContainer       *fMCTracksCont;     //! MC tracks
 
  TClonesArray               *fRecJetArray;      //! Array of the found jets
  TClonesArray               *fRecTrkArray;      //! Array of PicoTracks
  TClonesArray               *fMCJetArray;       //! Array of the found mc jets
  TClonesArray               *fMCPartArray;      //! Array of MC particles for given event
  
  TClonesArray               *fHFvertexing;      //! Array of reconstructed secondary vertex (b-tagged jets)
  
  map_int_bool               *fV0gTrkMap;
  
  AliLog::EType_t             fDebug;
  
  ClassDef(AliAnalysisTaskEmcalJetBtagSV, 2);  // analysis task for MC study
};

//-------------------------------------------------------------------------------------
inline Bool_t AliAnalysisTaskEmcalJetBtagSV::UserNotify() {
  
  if (fCheckMCCrossSection) {
    
    GetPythiaCrossSection();
    
    AliDebugF(1, MSGDEBUG("MC pT-hard weight: %e"), fMCWeight);
  }
  
  fSkipWeightInfo = kFALSE;
  return kTRUE;
}

#endif
