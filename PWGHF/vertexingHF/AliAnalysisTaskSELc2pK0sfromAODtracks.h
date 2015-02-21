#ifndef ALIANALYSISTASKSELC2PK0SFROMAODTRACKS_H
#define ALIANALYSISTASKSELC2PK0SFROMAODTRACKS_H

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
#include "AliRDHFCutsLctopK0sfromAODtracks.h"

class THnSparse;
class TH1F;
class TClonesArray;
class AliAODRecoCascadeHF;
class AliESDVertex;
class AliAODMCParticle;

class AliAnalysisTaskSELc2pK0sfromAODtracks : public AliAnalysisTaskSE 
{
 public:
  AliAnalysisTaskSELc2pK0sfromAODtracks();
  AliAnalysisTaskSELc2pK0sfromAODtracks(const Char_t* name, AliRDHFCutsLctopK0sfromAODtracks* cuts, Bool_t writeVariableTree=kTRUE);
  virtual ~AliAnalysisTaskSELc2pK0sfromAODtracks();

  // Implementation of interface methods  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void FillROOTObjects(AliAODRecoCascadeHF *lcobj, AliAODMCParticle *mcpart, AliAODMCParticle *mcdau1, AliAODMCParticle *mcdau2, Int_t mcnused);
  void MakeAnalysis(AliAODEvent *aod, TClonesArray *mcArray);


  // set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}

  void SetReconstructPrimVert(Bool_t a) { fReconstructPrimVert=a; }

  AliAODRecoCascadeHF* MakeCascadeHF(AliAODv0 *casc, AliAODTrack *trk, AliAODEvent *aod, AliAODVertex *vert);
  AliAODVertex* ReconstructSecondaryVertex(AliAODv0 *casc, AliAODTrack *trk, AliAODEvent *aod);


 private:

  AliAnalysisTaskSELc2pK0sfromAODtracks(const AliAnalysisTaskSELc2pK0sfromAODtracks &source);
  AliAnalysisTaskSELc2pK0sfromAODtracks& operator=(const AliAnalysisTaskSELc2pK0sfromAODtracks& source); 

  void DefineTreeVariables();
  void DefineGeneralHistograms();
  void DefineAnalysisHistograms();

  AliAODVertex *CallPrimaryVertex(AliAODv0 *v0, AliAODTrack *trk, AliAODEvent *evt);
  AliAODVertex* PrimaryVertex(const TObjArray *trkArray,AliVEvent *event);

  Bool_t fUseMCInfo;                 // Use MC info
  TList *fOutput;                    //! User output slot 1 // general histos
  TList *fOutputAll;                 //! User Output slot 3  //analysis histograms 
  TList *fListCuts;                  //! User output slot 2 // Cuts 
  TH1F *fCEvents;                    //! Histogram to check selected events
  TH1F *fHTrigger;                   //! Histogram to check Trigger
  TH1F *fHCentrality;                //! Histogram to check Centrality
  AliRDHFCutsLctopK0sfromAODtracks *fAnalCuts;// Cuts - sent to output slot 2
  Bool_t fIsEventSelected;          // flag for event selected
  Bool_t    fWriteVariableTree;     // flag to decide whether to write the candidate variables on a tree variables
  TTree    *fVariablesTree;         //! tree of the candidate variables after track selection on output slot 4
  Bool_t fReconstructPrimVert;       //Reconstruct primary vertex excluding candidate tracks
  Bool_t fIsMB;       //MB trigger event
  Bool_t fIsSemi;     //SemiCentral trigger event
  Bool_t fIsCent;     //Central trigger event
  Bool_t fIsINT7;     //INT7 trigger event
  Bool_t fIsEMC7;     //EMC7 trigger event
  Float_t *fCandidateVariables;   //! variables to be written to the tree
  AliAODVertex *fVtx1;            // primary vertex
  AliESDVertex *fV1;              // primary vertex
  Double_t fBzkG;                 // magnetic field value [kG]
  Float_t  fCentrality;           //Centrality
  Float_t  fTriggerCheck;         //Stores trigger information

  //--------------------- My histograms ------------------
  THnSparse* fHistoLcK0SpMass;         //Lc mass spectra

  TH1F* fHistoBachPt;      //! Bachelor pT histogram
  TH1F* fHistod0Bach;      //! Bachelor d0 histogram
  TH1F* fHistod0V0;        //! V0 d0 histogram
  TH1F* fHistod0d0;        //! Bachelor d0 * V0 d0 histogram
  TH1F* fHistoV0CosPA;     //! V0 cosine pointing angle to primary vertex
  TH1F* fHistoProbProton;  //! Probability to be proton histogram
  TH1F* fHistoDecayLength; //! Decay length histogram
  TH1F* fHistoK0SMass;     //! K0s mass histogram


  ClassDef(AliAnalysisTaskSELc2pK0sfromAODtracks,2); // class for Lc->p K0
};
#endif

