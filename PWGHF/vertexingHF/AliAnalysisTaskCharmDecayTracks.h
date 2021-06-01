#ifndef ALIANALYSISTASKCHARMDECAYTRACKS_H
#define ALIANALYSISTASKCHARMDECAYTRACKS_H

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///*************************************************************************
/// \class Class AliAnalysisTaskCharmDecayTracks
//////////////////////////////////////////////////////////////

class TTree;
class TH1F;
class AliAODTrack;

#include "AliNeutralTrackParam.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCharmDecayTracks : public AliAnalysisTaskSE
{
public:
  
  AliAnalysisTaskCharmDecayTracks();
  virtual ~AliAnalysisTaskCharmDecayTracks();
  
  virtual void UserCreateOutputObjects();
  virtual void Init(){};
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetSelectedHadron(Int_t pdg){fSelSpecies=pdg;}
  void SelectLcpK0s(){
    fSelSpecies=4122;
    fDecayMode=1;
  }
  void SetUseCandidatesFromDeltaAOD(){fMethod=1;}
  void SetUseCharmedHadronsFromKine(){fMethod=0;}
  
  void SetReadMC(Bool_t read){fReadMC=read;}
  void SetKeepNegIDtracks(Bool_t nid){fKeepNegID=nid;}
  void SetTrackCuts(AliESDtrackCuts* cuts){
    if(fTrCuts) delete fTrCuts;
    fTrCuts=new AliESDtrackCuts(*cuts);
  }
  void SetFilterMask(UInt_t mask=16){fFilterMask=mask;}
  
  void SetUsePhysicsSelection(Bool_t opt=kTRUE){
    fUsePhysSel=opt;
  }
  void SetTriggerMask(Int_t mask){
    fTriggerMask=mask;
  }
  void SetUsePileupCut(Bool_t opt=kTRUE){
    fUsePileupCut=kTRUE;
  }

  Bool_t IsTrackSelected(AliAODTrack* track);
  Bool_t IsV0Selected(AliAODv0* v0);
  enum EMesonSpecies {kDzero, kDplus, kDstar, kDs, kLc};
  
private:
  
  enum EVarsTree {kNumOfIntVar=5, kNumOfFloatVar=12};
  enum {kMaxLabel=1000000};

  AliAnalysisTaskCharmDecayTracks(const AliAnalysisTaskCharmDecayTracks &source);
  AliAnalysisTaskCharmDecayTracks& operator=(const AliAnalysisTaskCharmDecayTracks& source);
  void MapTrackLabels(AliAODEvent* aod);
  void MapV0Labels(AliAODEvent* aod, TClonesArray* arrayMC);
  Bool_t PrepareTreeVars(AliAODMCParticle* partD, TClonesArray* arrayMC, AliAODMCHeader* mcHeader);

  TList* fOutput;                  //!<! list send on output slot 0
  TH1F*  fHistNEvents;             //!<! histo with number of events
  TH1F*  fHistNCand;               //!<! histo with number of candidates
  TH1F*  fHistTrLab;               //!<! histo for tracks with same label
  TH2F*  fHistCluTPCDupLab;        //!<! histo for tracks with same label
  TH2F*  fHistCluTPCDupLabCorrel;  //!<! histo for tracks with same label
  TH2F*  fHistCluITSDupLabCorrel;  //!<! histo for tracks with same label
  TH2F*  fHistMomDupLab;           //!<! histo for tracks with same label
  TTree* fTrackTree;               //!<! output tree
  Int_t*   fTreeVarInt;            //!<! variables to be written to the tree
  Float_t* fTreeVarFloat;          //!<! variables to be written to the tree
  AliExternalTrackParam fTrPar1;   //!<! first track
  AliExternalTrackParam fTrPar2;   //!<! second track
  AliExternalTrackParam fTrPar3;   //!<! third track
  AliNeutralTrackParam fTrParV0;   //!<! V0
  AliAODVertex *fPVertexTrk;       //!<! primary vertex
  Int_t  fSelSpecies;              /// Charmed hadron species to analyse
  Int_t  fDecayMode;               /// 0=charged hadrons; 1=V0+bachelor
  UInt_t fFilterMask;              /// FilterMask
  AliESDtrackCuts* fTrCuts;        /// track selection
  Int_t fMapTrLabel[kMaxLabel];    /// map of track labels
  Int_t fMapV0Label[kMaxLabel];    /// map of V0 labels
  Bool_t fReadMC;                  ///  flag for access to MC
  Bool_t fUsePhysSel;              /// flag use/not use phys sel
  Bool_t fUsePileupCut;            /// flag use/not use phys pileup cut
  Int_t  fTriggerMask;             /// mask used in physics selection
  Bool_t fGoUpToQuark;             /// flag for definition of c,b origin
  Bool_t fKeepNegID;               /// flag to keep also track with negative ID
  Int_t fMethod;                   /// analysis from kine or from deltaAOD
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCharmDecayTracks,5);
  /// \endcond
};

#endif
