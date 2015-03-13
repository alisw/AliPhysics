#ifndef ALIANALYSISTASKSEHFJETS_H
#define ALIANALYSISTASKSEHFJETS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* Class AliAnalysisTaskSEHFJets:                                    *
 * AliAnalysisTaskSE for the extraction of the B-jet spectrum        *
 * using the properties of secondary verices as tagging observables. */

/* Mailto: andrea.rossi@cern.ch, elena.bruna@to.infn.it, svallero@to.infn.it, s.lapointe@cern.ch */


/* #include "AliAnalysisTaskSE.h"  */
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliHFJetsContainerVertex.h"
#include "AliHFJetsTagging.h"
#include "AliHFJetsTaggingVertex.h"
#include "AliRDHFJetsCuts.h"
#include <TArrayI.h>
#include <TArrayD.h>

class TH1F;
class TH2F;
//class TList;
class AliAODDEvent;
class AliAODMCHeader;
class AliAODRecoDecayHF2Prong;
class AliAODRecoDecayHF;
class AliAODMCParticle;
class AliAnalysisVertexingHF;
class AliNormalizationCounter;
class AliAODMCParticle;

class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#define RED  "\033[22;31;1m"
#define MAG  "\033[22;35;1m"
#define mage  "\033[22;35m"
#define cy  "\033[22;36m"
#define Bee  "\033[22;30m"

class AliAnalysisTaskSEHFJets : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskSEHFJets();
  AliAnalysisTaskSEHFJets(const char* name);
  virtual ~AliAnalysisTaskSEHFJets();

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
  void SetJetIdMethod(Int_t method) {fJetIdMeth = method;}
  void SetPtHardSelection(Bool_t ispthard, TString pthardmin, TString flavor) {fSelectPtHard = ispthard; fPtHardMin = pthardmin; fFlavor = flavor;}
  void SetJetTaggingRadius(Double_t tagradius) {fTaggingRadius = tagradius;};

  void SetTagger(AliHFJetsTaggingVertex *tagger){fTagger=(AliHFJetsTaggingVertex*)tagger->Clone("fTagger");}
  void SetCuts(AliRDHFJetsCuts *cuts){delete fCutsHFjets; fCutsHFjets=(AliRDHFJetsCuts*)cuts->Clone("fCutsHFjets");}

  void SetJetContName(char *j) {fJetContName=j;}
  void SetTrackContName(char *t) {fTrackContName=t;}
  void SetMcJetContName(char *mcj) {fMcJetContName=mcj;}
  void SetMcTrackContName(char *mct) {fMcTrackContName=mct;}

  AliAnalysisUtils* GetAnalysisUtils() { return fAnalysisUtils; }
  void   SetAnalysisUtils(AliAnalysisUtils* utils){ fAnalysisUtils = utils; }

 protected:
  Bool_t                         GetArrays();
  /* void                        ExecOnce(); */

  /* AliJetContainer            *fJetsCont;                   //!Jets */
  /* AliParticleContainer       *fTracksCont;                 //!Tracks */
  /* AliClusterContainer        *fCaloClustersCont;           //!Clusters   */


  TH1                        *fHistTrials;                 //!trials from pyxsec.root
  TH1                        *fHistEvents;                 //!total number of events per pt hard bin
  TProfile                   *fHistXsection;               //!x section from pyxsec.root


    AliParticleContainer       *fMCTracksCont;                                   //!MC tracks
 private:
  AliAnalysisTaskSEHFJets(const AliAnalysisTaskSEHFJets&); // copy constructo not implemented yet
  AliAnalysisTaskSEHFJets& operator=(const AliAnalysisTaskSEHFJets&); // assignment operator not implemented yet

  // temporary
  void GetFlavour2Methods(AliEmcalJet *jet, Double_t (&partonnat)[2], Double_t (&ptpart)[2], Double_t &contribution, Double_t radius);
  Bool_t PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials);

  void StoreTrackReference(AliAODTrack *track);      // Store pointer to global track
  void ResetTrackReference();                        // Reset all pointers to global tracks

  void GetJetMatching(const TList *genJetsList, const Int_t &kGenJets, TClonesArray *fTrackArraymc, const TList *recJetsList,  const Int_t &kRecJets, TClonesArray *fTrackArrayreco,  TArrayI &iMatchIndex, TArrayF &fPtFraction, Int_t iDebug, Float_t maxDist, Int_t mode);

  Double_t GetFractionOfJet(const AliEmcalJet *recJet, const AliEmcalJet *genJet, TClonesArray *fTrackArrayrec, TClonesArray *fTrackArraymc,Int_t mode);
  
  AliHFJetsContainerVertex *fhJets;    // reco jet properties
  AliHFJetsContainerVertex *fhQaVtx;   // vertices properties
  AliHFJetsContainerVertex *fhBJets;   // B-jet properties
  AliHFJetsContainerVertex *fhJetVtx;  // properties of vertices within the jet 
  AliHFJetsContainerVertex *fhJetVtxData;  // properties of vertices within the jet 
 
 AliVEvent                         *fEvent;                     //! Input event 


  Bool_t fCorrMode;                    // enable correction or data modes
  Int_t fJetIdMeth;                    // method for jet flavour type selection
  Bool_t fSelectPtHard;                // select pythia events based on their pt hard min
  Bool_t fUseTriggerData;              // use emacal trigger
  TString fPtHardMin;                  // minimum pt hard
  TString fFlavor;                     // flavour
  Double_t fTaggingRadius;             // radius used in tagging the jet flavour
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


  TClonesArray                      *fJetArray;                  //! Array of the found jets
    TClonesArray                      *fMcJetArray;                  //! Array of the found mc jets
    TString                           fJetContName;                 //  Name of the found jet array
  TString                           fTrackContName;                 //  Name of the found track array
      TString                           fMcJetContName;                 //  Name of the found jet array
  TString                           fMcTrackContName;                 //  Name of the found track array


  TList *fOutputList;                  // list of output objects
  ClassDef(AliAnalysisTaskSEHFJets,1); // analysis task for MC study
};

#endif
