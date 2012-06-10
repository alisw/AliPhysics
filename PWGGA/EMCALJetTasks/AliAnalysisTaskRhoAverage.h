#ifndef ALIANALYSISTASKRHOAVERAGE_H
#define ALIANALYSISTASKRHOAVERAGE_H

// $Id$

class TClonesArray;
class TList;
class AliEmcalJet;

#include "AliAnalysisTaskRhoBase.h"

class AliAnalysisTaskRhoAverage : public AliAnalysisTaskRhoBase {

 public:
  AliAnalysisTaskRhoAverage();
  AliAnalysisTaskRhoAverage(const char *name);
  virtual ~AliAnalysisTaskRhoAverage() {}
  
  void                   UserExec(Option_t*);

  void                   SetClustersName(const char *n)                        { fClustersName  = n    ; }
  void                   SetEtaLimits(Double_t emin, Double_t emax)            { fEtaMin        = emin ; fEtaMax = emax  ; }
  void                   SetJetsName(const char *n)                            { fJetsName      = n    ; }  
  void                   SetPhiLimits(Double_t pmin, Double_t pmax)            { fPhiMin        = pmin ; fPhiMax = pmax  ; }
  void                   SetPtMin(Double_t pt)                                 { fPtMin         = pt   ; }
  void                   SetTracksName(const char *n)                          { fTracksName    = n    ; }
  
 protected:
  void                   ExecOnce();
  Bool_t                 IsJetCluster(AliEmcalJet* jet, Int_t iclus) const;
  Bool_t                 IsJetTrack(AliEmcalJet* jet, Int_t itrack)  const;

  TString                fTracksName;                    // name of track collection
  TString                fClustersName;                  // name of clusters collection
  TString                fJetsName;                      // name of jet collection
  Double_t               fEtaMin;                        // minimum eta
  Double_t               fEtaMax;                        // maximum eta
  Double_t               fPhiMin;                        // minimum phi
  Double_t               fPhiMax;                        // maximum phi
  Double_t               fPtMin;                         // minimum pt
  TClonesArray          *fClusters;                      //!input clusters
  TClonesArray          *fJets;                          //!input jets
  TClonesArray          *fTracks;                        //!input tracks

  AliAnalysisTaskRhoAverage(const AliAnalysisTaskRhoAverage&);             // not implemented
  AliAnalysisTaskRhoAverage& operator=(const AliAnalysisTaskRhoAverage&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoAverage, 2); // Rho task, method: sum of all particle pt / full acceptance area
};
#endif
