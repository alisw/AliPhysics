#ifndef ALIANALYSISTASKEMCALJETBTAGSV_H
#define ALIANALYSISTASKEMCALJETBTAGSV_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* Class AliAnalysisTaskEmcalJetBtagSV:                                    *
 * AliAnalysisTaskSE for the extraction of the B-jet spectrum        *
 * using the properties of secondary verices as tagging observables. */

/* Mailto: andrea.rossi@cern.ch, elena.bruna@to.infn.it, svallero@to.infn.it, s.lapointe@cern.ch */

class TH1F;
class TH2F;
//class TList;
class AliAODDEvent;
class AliAODMCHeader;
class AliAODRecoDecayHF2Prong;
class AliAODRecoDecayHF;
class AliAODMCParticle;
class AliAnalysisUtils;
class AliAnalysisVertexingHF;
class AliNormalizationCounter;
class AliAODMCParticle;


class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliHFJetsContainerVertex.h"
#include "AliHFJetsTagging.h"
#include "AliHFJetsTaggingVertex.h"
#include "AliRDHFJetsCuts.h"
#include <TArrayI.h>
#include <TArrayD.h>

#include "AliHFjetsUtils.h"

class AliAnalysisTaskEmcalJetBtagSV : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetBtagSV();
  AliAnalysisTaskEmcalJetBtagSV(const char *name);
  virtual ~AliAnalysisTaskEmcalJetBtagSV();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual Bool_t  UserNotify();
  virtual void Terminate(Option_t *option);  

  void AnalyseCorrectionsMode();
  void AnalyseDataMode();

  // Setters/Getters
  void SetCorrectionsMode(Bool_t mode) {fCorrMode = mode;}
  void SetDoBkgRejection(Bool_t rej) {fDoBkgRej = rej;}
  void SetRecoJetsBranch(const char* branch) {fRecoJetsBranch = branch;}
  void SetMcJetsBranch(const char* branch) {fMcJetsBranch = branch;}
  void SetPtHardSelection(Bool_t ispthard, TString pthardmin, TString flavor) { fSelectPtHard = ispthard; fPtHardMin = pthardmin; fFlavor = flavor; }
  void SetJetTaggingRadius(Double_t tagradius) {fTaggingRadius = tagradius;};

  void DoSecondaryVertexQA(Bool_t doqa) {fDoQAVtx = doqa;}

  void SetTagger(AliHFJetsTaggingVertex *tagger)
  {
    if (fTagger) delete fTagger;
    fTagger = (AliHFJetsTaggingVertex *) tagger->Clone("fTagger");
  }
  
  void SetCuts(AliRDHFJetsCuts *cuts)
  {
    if (fCutsHFjets) delete fCutsHFjets;
    fCutsHFjets = (AliRDHFJetsCuts *)cuts->Clone("fCutsHFjets");
  }

  void SetJetContName(char *name)           { fJetContName     = name; }
  void SetTrackContName(char *name)         { fTrackContName   = name; }
  void SetExternalRhoTaskName(char *name)   { fRhoTaskName     = name; }
  void SetMcJetContName(char *name)         { fMcJetContName   = name; }
  void SetMcTrackContName(char *name)       { fMcTrackContName = name; }
  void SetMcExternalRhoTaskName(char *name) { fMcRhoTaskName   = name; } 
  
  AliAnalysisUtils *GetAnalysisUtils() { return fAnalysisUtils; }
  void   SetAnalysisUtils(AliAnalysisUtils* utils) { fAnalysisUtils = utils; }

 protected:
  Bool_t                      GetArrays();

  TH1                        *fHistTrials;                 //!trials from pyxsec.root
  TH1                        *fHistEvents;                 //!total number of events per pt hard bin
  TProfile                   *fHistXsection;               //!x section from pyxsec.root
  AliParticleContainer       *fMCTracksCont;                                   //!MC tracks
 private:
  AliAnalysisTaskEmcalJetBtagSV(const AliAnalysisTaskEmcalJetBtagSV&); // copy constructo not implemented yet
  AliAnalysisTaskEmcalJetBtagSV& operator=(const AliAnalysisTaskEmcalJetBtagSV&); // assignment operator not implemented yet

  // temporary
  void GetFlavour2Methods(AliEmcalJet *jet, Double_t (&partonnat)[2], Double_t (&ptpart)[2], Double_t radius);
  Bool_t PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials);

  Double_t GetExternalRho();   // rho
  Double_t GetMcExternalRho(); // rho  

  void StoreTrackReference(AliAODTrack *track);      // Store pointer to global track
  void ResetTrackReference();                        // Reset all pointers to global tracks
  
  /* AliHFJetsContainerVertex *fhJets;    // reco jet properties */
  AliHFJetsContainerVertex *fhQaVtx;   // vertices properties
  AliHFJetsContainerVertex *fhJetVtx;  // properties of vertices within the jet 
  AliHFJetsContainerVertex *fhJetVtxData;  // properties of vertices within the jet 
 
  AliVEvent                         *fEvent;                     //! Input event 

  Bool_t fCorrMode;                    // enable correction or data modes
  Bool_t fDoBkgRej;
  const char* fRecoJetsBranch;         // name of the AOD RECO-jets branch  
  const char* fMcJetsBranch;           // name of the AOD MC-jets branch  
  Bool_t fSelectPtHard;                // select pythia events based on their pt hard min
  Bool_t fUseTriggerData;              // use emacal trigger
  TString fPtHardMin;                  // minimum pt hard
  TString fFlavor;                     // flavour
  Double_t fTaggingRadius;             // radius used in tagging the jet flavour
  Bool_t fDoQAVtx;                     // enable output of qa on secondary vertex
  TString fFiredClass;                 // emcal trigger class
  
  TH1F *fNentries;                     //! histo for event counting and checks
  TH1F *fNtriggers;                    //! histo for event counting and checks
  AliHFJetsTaggingVertex *fTagger;     // Jet Tagging object
  AliRDHFJetsCuts *fCutsHFjets;        // specific algo jet cut object 

  AliAnalysisUtils*     fAnalysisUtils;      // points to class with common analysis utilities
   
  TClonesArray *fbJetArray;            //! b-tagged jets
  TClonesArray *fArrayMC;              // array of MC particles for given event


  // Store pointers to global tracks for pid and dca
  AliAODTrack **fGTIp;                //! Array of pointers, just nicely sorted according to the id
  const UShort_t  fTrackBuffSize;          //! Size of the above array, ~12000 for PbPb
  AliAODTrack **fGTIn;                //! Array of pointers, just nicely sorted according to the id


  TClonesArray                     *fJetArray;         //! Array of the found jets
  TClonesArray                     *fMcJetArray;       //! Array of the found mc jets
  TString                           fJetContName;      //  Name of the found jet array
  TString                           fTrackContName;    //  Name of the found track array
  TString                           fRhoTaskName;      //  Name of the rho task
  TString                           fMcJetContName;    //  Name of the found jet array
  TString                           fMcTrackContName;  //  Name of the found track array
  TString                           fMcRhoTaskName;    //  Name of the rho task
  
  TList                            *fOutputList;                  // list of output objects
  ClassDef(AliAnalysisTaskEmcalJetBtagSV, 1); // analysis task for MC study
};

#endif
