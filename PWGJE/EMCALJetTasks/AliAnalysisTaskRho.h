#ifndef ALIANALYSISTASKRHO_H
#define ALIANALYSISTASKRHO_H

// $Id$

#include "AliAnalysisTaskRhoBase.h"

class AliAnalysisTaskRho : public AliAnalysisTaskRhoBase {

 public:
  AliAnalysisTaskRho();
  AliAnalysisTaskRho(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRho() {}

  void             SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }

  static AliAnalysisTaskRho* AddTaskRhoNew (
    const char    *nTracks                        = "usedefault",
    const char    *nClusters                      = "usedefault",
    const char    *nRho                           = "Rho",
    Double_t       jetradius                      = 0.2,
    UInt_t         acceptance                     = AliEmcalJet::kTPCfid,
    AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
    const Bool_t   histo                          = kFALSE,
    AliJetContainer::ERecoScheme_t rscheme        = AliJetContainer::pt_scheme,
    const char    *suffix                         = ""
);

 protected:
  Bool_t           Run();

  UInt_t           fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation

  AliAnalysisTaskRho(const AliAnalysisTaskRho&);             // not implemented
  AliAnalysisTaskRho& operator=(const AliAnalysisTaskRho&);  // not implemented
  
  ClassDef(AliAnalysisTaskRho, 10); // Rho task
};
#endif
