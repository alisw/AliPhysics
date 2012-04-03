#ifndef ALIESDJETTASK_H
#define ALIESDJETTASK_H

// $Id$

class TClonesArray;
class AliESDtrackCuts;
class TList;
class TH1;
class TH2;

#include "AliAnalysisTaskSE.h"

class AliEsdJetTask : public AliAnalysisTaskSE {
 public:
  AliEsdJetTask(const char *name); //was set to 0
  AliEsdJetTask();//=0);
  virtual ~AliEsdJetTask();

  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *option);

  void         SetAlgo(Int_t a)                 { fAlgo          = a;  }
  void         SetClusName(const char *n)       { fCaloName      = n;  }
  void         SetJetsName(const char *n)       { fJetsName      = n;  }
  void         SetMinJetTrackPt(Double_t min)   { fMinJetTrackPt = min;}
  void         SetRadius(Double_t r)            { fRadius        = r;  }
  void         SetTracksName(const char *n)     { fTracksName    = n;  }
  void         SetType(Int_t t)                 { fType          = t;  }

 protected:
  void FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius, Float_t fCent);

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fJetsName;               // name of jet collection
  Int_t                  fAlgo;                   // algo (0==kt, 1==antikt)
  Double_t               fRadius;                 // jet radius
  Int_t                  fType;                   // jet type (0=all, 1=ch, 2=neutral)
  Double_t               fMinJetTrackPt;          // min jet track momentum
  TClonesArray          *fJets;                   //!jet collection

 private:
  AliEsdJetTask(const AliEsdJetTask&);            // not implemented
  AliEsdJetTask &operator=(const AliEsdJetTask&); // not implemented

  ClassDef(AliEsdJetTask, 2) // Jet producing task
};
#endif
