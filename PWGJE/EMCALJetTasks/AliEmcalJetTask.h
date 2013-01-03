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

  void                   SetAlgo(Int_t a)                 { fAlgo          = a     ; }
  void                   SetClusName(const char *n)       { fCaloName      = n     ; }
  void                   SetJetsName(const char *n)       { fJetsName      = n     ; }
  void                   SetMinJetArea(Double_t a)        { fMinJetArea    = a     ; }
  void                   SetMinJetClusPt(Double_t min)    { fMinJetClusPt  = min   ; }
  void                   SetMinJetPt(Double_t j)          { fMinJetPt      = j     ; }
  void                   SetMinJetTrackPt(Double_t min)   { fMinJetTrackPt = min   ; }
  void                   SetRadius(Double_t r)            { fRadius        = r     ; }
  void                   SetTracksName(const char *n)     { fTracksName    = n     ; }
  void                   SetType(Int_t t)                 { fType          = t     ; }
  void                   SetEtaRange(Double_t emi, Double_t ema) {fEtaMin = emi; fEtaMax = ema; }
  void                   SetPhiRange(Double_t pmi, Double_t pma) {fPhiMin = pmi; fPhiMax = pma; }
  void                   SetGhostArea(Double_t gharea)    { fGhostArea      = gharea;  }

 protected:
  void                   FindJets();
  Bool_t                 DoInit();

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fJetsName;               // name of jet collection
  Int_t                  fAlgo;                   // algo (0==kt, 1==antikt)
  Double_t               fRadius;                 // jet radius
  Int_t                  fType;                   // jet type (0=all, 1=ch, 2=neutral)
  Double_t               fMinJetTrackPt;          // min jet track momentum   (applied before clustering)
  Double_t               fMinJetClusPt;           // min jet cluster momentum (applied before clustering)
  Double_t               fPhiMin;                 // minimum phi for constituents (applied before clustering)
  Double_t               fPhiMax;                 // maximum phi for constituents (applied before clustering)
  Double_t               fEtaMin;                 // minimum eta for constituents (applied before clustering)
  Double_t               fEtaMax;                 // maximum eta for constituents (applied before clustering)
  Double_t               fMinJetArea;             // min area to keep jet in output
  Double_t               fMinJetPt;               // min jet pt to keep jet in output
  Double_t               fGhostArea;              // ghost area
  Bool_t                 fIsInit;                 //!=true if already initialized
  Bool_t                 fIsMcPart;               //!=true if MC particles are given as input
  Bool_t                 fIsEmcPart;              //!=true if emcal particles are given as input (for clusters)
  TClonesArray          *fJets;                   //!jet collection
  AliVEvent             *fEvent;                  //!current event
  TClonesArray          *fTracks;                 //!tracks collection
  TClonesArray          *fClus;                   //!cluster collection

 private:
  AliEmcalJetTask(const AliEmcalJetTask&);            // not implemented
  AliEmcalJetTask &operator=(const AliEmcalJetTask&); // not implemented

  ClassDef(AliEmcalJetTask, 5) // Jet producing task
};
#endif
