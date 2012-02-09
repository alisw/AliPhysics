#ifndef ALIESDJETTASK_H
#define ALIESDJETTASK_H

// $Id$

class TClonesArray;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliEsdJetTask : public AliAnalysisTaskSE {
 public:
  AliEsdJetTask(const char *name=0);
  virtual ~AliEsdJetTask();

  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *option);

  void         SetAlgo(Int_t a)                      { fAlgo           = a; }
  void         SetCaloName(const char *n)            { fCaloName       = n; }
  void         SetHadCorr(Double_t c)                { fHadCorr        = c; }
  void         SetJetsName(const char *n)            { fJetsName       = n; }
  void         SetPrimTrackCuts(AliESDtrackCuts *c)  { fPrimTrCuts     = c; }
  void         SetPrimTracksName(const char *n)      { fPrimTracksName = n; }
  void         SetRadius(Double_t r)                 { fRadius         = r; }
  void         SetType(Int_t t)                      { fType           = t; }

 protected:
  void FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius);

  AliESDtrackCuts       *fPrimTrCuts;             // track cuts
  TString                fPrimTracksName;         // name of track collection (if "" use branch)
  TString                fJetsName;               // name of jet collection
  TString                fCaloName;               // name of calo collection
  Int_t                  fAlgo;                   // algo (0==kt, 1==antikt)
  Double_t               fRadius;                 // jet radius
  Int_t                  fType;                   // jet type (0=all, 1=ch, 2=neutral)
  Double_t               fHadCorr;                // hadronic correction
  TClonesArray          *fJets;                   //!jet collection
  TClonesArray          *fNeutrals;               //!neutrals collection

 private:
  AliEsdJetTask(const AliEsdJetTask&);            // not implemented
  AliEsdJetTask &operator=(const AliEsdJetTask&); // not implemented

  ClassDef(AliEsdJetTask, 1) // Jet producing task
};
#endif
