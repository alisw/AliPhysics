#ifndef ALIEMCALJETTASK_H
#define ALIEMCALJETTASK_H

// $Id$

class TClonesArray;
class AliVEvent;

#include "AliAnalysisTaskSE.h"

class AliEmcalJetTask : public AliAnalysisTaskSE {
 public:
  AliEmcalJetTask();
  AliEmcalJetTask(const char *name);
  virtual ~AliEmcalJetTask();

  void                   UserCreateOutputObjects();
  void                   UserExec(Option_t *option);
  void                   Terminate(Option_t *option);

  void                   SetTracksName(const char *n)     { fTracksName    = n     ; }
  void                   SetClusName(const char *n)       { fCaloName      = n     ; }
  void                   SetJetsName(const char *n)       { fJetsName      = n     ; }
  void                   SetMC(Bool_t mc = kTRUE)         { fMC            = mc    ; }
  void                   SetAlgo(Int_t a)                 { fAlgo          = a     ; }
  void                   SetRadius(Double_t r)            { fRadius        = r     ; }
  void                   SetType(Int_t t)                 { fType          = t     ; }
  void                   SetMinJetClusPt(Double_t min)    { fMinJetClusPt  = min   ; }
  void                   SetMinJetTrackPt(Double_t min)   { fMinJetTrackPt = min   ; }
  void                   SetMinJetArea(Double_t a)        { fMinJetArea    = a     ; }
  void                   SetMinJetPt(Double_t j)          { fMinJetPt      = j     ; }


 protected:
  void                   FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius, Float_t /*cent*/);

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fJetsName;               // name of jet collection
  Bool_t                 fMC;                     // true = MC analysis
  Int_t                  fAlgo;                   // algo (0==kt, 1==antikt)
  Double_t               fRadius;                 // jet radius
  Int_t                  fType;                   // jet type (0=all, 1=ch, 2=neutral)
  Double_t               fMinJetTrackPt;          // min jet track momentum   (applied before clustering)
  Double_t               fMinJetClusPt;           // min jet cluster momentum (applied before clustering)
  Double_t               fMinJetArea;             // min area to keep jet in output
  Double_t               fMinJetPt;               // min jet pt to keep jet in output
  TClonesArray          *fJets;                   //!jet collection
  AliVEvent             *fEvent;                  //!current event

 private:
  AliEmcalJetTask(const AliEmcalJetTask&);            // not implemented
  AliEmcalJetTask &operator=(const AliEmcalJetTask&); // not implemented

  ClassDef(AliEmcalJetTask, 3) // Jet producing task
};
#endif
