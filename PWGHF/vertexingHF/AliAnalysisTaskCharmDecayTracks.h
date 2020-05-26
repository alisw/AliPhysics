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

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCuts.h"

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
  enum EMesonSpecies {kDzero, kDplus, kDstar, kDs, kLc};
  
private:
  
  enum EVarsTree {kNumOfIntVar=1, kNumOfFloatVar=9};
  enum {kMaxLabel=1000000};

  AliAnalysisTaskCharmDecayTracks(const AliAnalysisTaskCharmDecayTracks &source);
  AliAnalysisTaskCharmDecayTracks& operator=(const AliAnalysisTaskCharmDecayTracks& source);
  void MapTrackLabels(AliAODEvent* aod);

  TList* fOutput;                  //!<! list send on output slot 0
  TH1F*  fHistNEvents;             //!<! histo with number of events
  TTree* fTrackTree;               //!<! output tree
  Int_t*   fTreeVarInt;            //!<! variables to be written to the tree
  Float_t* fTreeVarFloat;          //!<! variables to be written to the tree
  AliExternalTrackParam fTrPar1;   //!<! first track
  AliExternalTrackParam fTrPar2;   //!<! second track
  AliExternalTrackParam fTrPar3;   //!<! third track
  Int_t  fSelSpecies;              /// Charmed hadron species to analyse
  UInt_t fFilterMask;              /// FilterMask
  AliESDtrackCuts* fTrCuts;        /// track selection
  Int_t fMapTrLabel[kMaxLabel];              /// map of track labels
  Bool_t fReadMC;                  ///  flag for access to MC
  Bool_t fUsePhysSel;              /// flag use/not use phys sel
  Bool_t fUsePileupCut;            /// flag use/not use phys pileup cut
  Int_t  fTriggerMask;             /// mask used in physics selection
  Bool_t fGoUpToQuark;             /// flag for definition of c,b origin
  Bool_t fKeepNegID;               /// flag to keep also track with negative ID

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCharmDecayTracks,1); 
  /// \endcond
};

#endif
