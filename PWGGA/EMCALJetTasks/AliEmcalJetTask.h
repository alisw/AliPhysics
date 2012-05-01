#ifndef ALIEMCALJETTASK_H
#define ALIEMCALJETTASK_H

// $Id$

class TClonesArray;
class TList;
class TH1;
class TH2;

#include "AliAnalysisTaskSE.h"

class AliEmcalJetTask : public AliAnalysisTaskSE {
 public:
  AliEmcalJetTask();
  AliEmcalJetTask(const char *name);
  virtual ~AliEmcalJetTask();

  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *option);

  void         SetAlgo(Int_t a)                 { fAlgo          = a;  }
  void         SetClusName(const char *n)       { fCaloName      = n;  }
  void         SetJetsName(const char *n)       { fJetsName      = n;  }
  void         SetMinJetClusPt(Double_t min)    { fMinJetClusPt  = min;}
  void         SetMinJetTrackPt(Double_t min)   { fMinJetTrackPt = min;}
  void         SetRadius(Double_t r)            { fRadius        = r;  }
  void         SetTracksName(const char *n)     { fTracksName    = n;  }
  void         SetType(Int_t t)                 { fType          = t;  }

 protected:

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fJetsName;               // name of jet collection
  Int_t                  fAlgo;                   // algo (0==kt, 1==antikt)
  Double_t               fRadius;                 // jet radius
  Int_t                  fType;                   // jet type (0=all, 1=ch, 2=neutral)
  Double_t               fMinJetTrackPt;          // min jet track momentum
  Double_t               fMinJetClusPt;           // min jet cluster momentum
  TClonesArray          *fJets;                   //!jet collection

 private:
  AliEmcalJetTask(const AliEmcalJetTask&);            // not implemented
  AliEmcalJetTask &operator=(const AliEmcalJetTask&); // not implemented

  ClassDef(AliEmcalJetTask, 1) // Jet producing task
};
#endif
