#ifndef ALIEMCALCLUSTERMAKER_H
#define ALIEMCALCLUSTERMAKER_H

// $Id$

class TClonesArray;
class AliEMCALRecoUtils;

#include "AliAnalysisTaskEmcal.h"

class AliEmcalClusterMaker : public AliAnalysisTaskEmcal {

 public:
  AliEmcalClusterMaker();
  AliEmcalClusterMaker(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliEmcalClusterMaker();

  void                   UserCreateOutputObjects();
  void                   SetOutClusName(const char *n)           { fOutCaloName    = n    ; }
  void                   SetRecoUtils(AliEMCALRecoUtils *ru) { fRecoUtils = ru; }

 protected:
  Bool_t                 Run()                                          ;
  void                   ExecOnce()                                     ;

  TString                fOutCaloName;               // name of output clusters
  AliEMCALRecoUtils     *fRecoUtils;                 // pointer to reco utils
  Bool_t                 fEsdMode;                   //!ESD/AOD mode
  TClonesArray          *fOutClusters;               //!output cluster collection

 private:
  AliEmcalClusterMaker(const AliEmcalClusterMaker&);            // not implemented
  AliEmcalClusterMaker &operator=(const AliEmcalClusterMaker&); // not implemented

  ClassDef(AliEmcalClusterMaker, 1) // Emcal cluster maker
};
#endif
