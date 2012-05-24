#ifndef ALIANALYSISTASKRHOAVERAGE_cxx
#define ALIANALYSISTASKRHOAVERAGE_cxx

// $Id$

class TList;
class TH1F;
class TH2F;
class TClonesArray;
class TString;
class TF1;

#include <TParameter.h>

#include "AliAnalysisTaskRhoBase.h"

class AliAnalysisTaskRhoAverage : public AliAnalysisTaskRhoBase {

 public:
  AliAnalysisTaskRhoAverage();
  AliAnalysisTaskRhoAverage(const char *name);
  AliAnalysisTaskRhoAverage(const char *name, Bool_t histo);
  virtual ~AliAnalysisTaskRhoAverage() {}
  
  virtual void           UserExec(Option_t*);
  virtual void           Terminate(Option_t*);

  void                   SetTracksName(const char *n)                          { fTracksName    = n    ; }
  void                   SetClustersName(const char *n)                        { fClustersName  = n    ; }
  void                   SetJetsName(const char *n)                            { fJetsName      = n    ; }  
  void                   SetEtaLimits(Double_t emin, Double_t emax)            { fEtaMin        = emin ; fEtaMax = emax  ; }
  void                   SetPhiLimits(Double_t pmin, Double_t pmax)            { fPhiMin        = pmin ; fPhiMax = pmax  ; }
  void                   SetPtMin(Double_t pt)                                 { fPtMin         = pt   ; }
  
 protected:
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

  AliAnalysisTaskRhoAverage(const AliAnalysisTaskRhoAverage&);             // not implemented
  AliAnalysisTaskRhoAverage& operator=(const AliAnalysisTaskRhoAverage&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoAverage, 0); // Rho task, method: sum of all particle pt / full acceptance area
};
#endif
