#ifndef ALIEMCALCLUSTERMAKER_H
#define ALIEMCALCLUSTERMAKER_H

class TClonesArray;
class TH1F;
class TH2F;
class AliEMCALRecoUtils;

#include "AliAnalysisTaskEmcal.h"

class AliEmcalClusterMaker : public AliAnalysisTaskEmcal {

 public:
  AliEmcalClusterMaker();
  AliEmcalClusterMaker(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliEmcalClusterMaker();

  void                   UserCreateOutputObjects();
  void                   SetOutClusName(const char *n)       { fOutCaloName    = n;  }
  void                   SetRecoUtils(AliEMCALRecoUtils *ru) { fRecoUtils      = ru; }

 protected:
  Bool_t                 Run()                                          ;
  void                   ExecOnce()                                     ;

  TString                fOutCaloName;               // name of output clusters; if empty updates old clusters instead of creating a new collection
  AliEMCALRecoUtils     *fRecoUtils;                 // pointer to reco utils
  Bool_t                 fEsdMode;                   //!ESD/AOD mode
  TClonesArray          *fOutClusters;               //!output cluster collection
  TH1F                  *fEnergyDistBefore;          //!energy distribution before
  TH2F                  *fEtaPhiDistBefore;          //!eta/phi distribution before
  TH2F                  *fEnergyTimeHistBefore;      //!energy/time distribution before
  TH1F                  *fEnergyDistAfter;           //!energy distribution after
  TH2F                  *fEtaPhiDistAfter;           //!eta/phi distribution after
  TH2F                  *fEnergyTimeHistAfter;       //!energy/time distribution after
  TH1F                  *fEnergyExoticClusters;      //!energy of exotic clusters

 private:
  AliEmcalClusterMaker(const AliEmcalClusterMaker&);            // not implemented
  AliEmcalClusterMaker &operator=(const AliEmcalClusterMaker&); // not implemented

  ClassDef(AliEmcalClusterMaker, 2) // Emcal cluster maker
};
#endif
