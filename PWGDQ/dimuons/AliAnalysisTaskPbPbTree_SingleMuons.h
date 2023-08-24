#ifndef AliAnalysisTaskPbPbTree_SingleMuons_H
#define AliAnalysisTaskPbPbTree_SingleMuons_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "TTreeStream.h"

class TObjArray;
class TLorentzVector;
class AliVParticle;
class AliAODEvent;
class AliMuonTrackCuts;
class AliTriggerAnalysis;
class AliVMultiplicity;
class AliAODTrack;


class AliAnalysisTaskPbPbTree_SingleMuons: public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskPbPbTree_SingleMuons();
  AliAnalysisTaskPbPbTree_SingleMuons(const char *name);
  virtual ~AliAnalysisTaskPbPbTree_SingleMuons();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  virtual void NotifyRun();

 private:
  AliAnalysisTaskPbPbTree_SingleMuons(const AliAnalysisTaskPbPbTree_SingleMuons&);
  AliAnalysisTaskPbPbTree_SingleMuons& operator=(const AliAnalysisTaskPbPbTree_SingleMuons&);

 //protected:

  TTree     *fOutputTree;      //! tree output
  TH1D      *fhNEv;            //! histo

  Int_t fCountCINT7;               // counter
  Int_t fCountCMUL7;               //counter
  Int_t fCountCMLL7;               //counter
  Int_t fCountCMSL7;               //counter
  Int_t fCountCMSH7;               //counter

  Int_t		fNMuons;		              // muon tracks in the event
  AliAODEvent*  fAODEvent;          //! AOD event
  char   	fTrigClass[800];	        // fired trigger classes

  AliMuonTrackCuts* fMuonTrackCuts;
  Float_t   fPercentV0M;             //! percentile V0

  TObjArray *fMuonTracks;      //! array of muon tracks

  ClassDef(AliAnalysisTaskPbPbTree_SingleMuons,1);
};

#endif
