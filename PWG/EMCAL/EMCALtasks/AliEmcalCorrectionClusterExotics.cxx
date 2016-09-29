// AliEmcalCorrectionClusterExotics
//

#include "AliEmcalCorrectionClusterExotics.h"

#include <TH2F.h>
#include <TList.h>

#include "AliClusterContainer.h"
#include "AliEMCALRecoUtils.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterExotics);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterExotics> AliEmcalCorrectionClusterExotics::reg("AliEmcalCorrectionClusterExotics");

//________________________________________________________________________
AliEmcalCorrectionClusterExotics::AliEmcalCorrectionClusterExotics() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterExotics"),
  fEtaPhiDistBefore(0),
  fEtaPhiDistAfter(0),
  fEnergyExoticClusters(0)
{
  // Default constructor
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
}

//________________________________________________________________________
AliEmcalCorrectionClusterExotics::~AliEmcalCorrectionClusterExotics()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionClusterExotics::Initialize()
{
  // Initialization
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Initialize();
  // Do base class initializations and if it fails -> bail out
  //AliAnalysisTaskEmcal::ExecOnce();
  //if (!fInitialized) return;
  
  GetProperty("createHistos", fCreateHisto);

  /*AddContainer(kCluster);
  Double_t minE = 0.;
  GetProperty("clusterEMin", minE);
  Double_t minPt = 0.;
  GetProperty("clusterPtMin", minPt);
  
  // Settings from sample run macro
  fClusCont->SetClusECut(minE);
  fClusCont->SetClusPtCut(minPt);*/
  
  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;
  fRecoUtils->SwitchOnRejectExoticCluster();
  if (fRecoUtils)
    fRecoUtils->Print("");
  
  // Create my user objects.
  if (fCreateHisto){
    fEtaPhiDistBefore = new TH2F("hEtaPhiDistBefore","hEtaPhiDistBefore;#eta;#phi",280,-0.7,0.7,800,1.3,3.3);
    fOutput->Add(fEtaPhiDistBefore);
    fEtaPhiDistAfter = new TH2F("hEtaPhiDistAfter","hEtaPhiDistAfter;#eta;#phi",280,-0.7,0.7,800,1.3,3.3);
    fOutput->Add(fEtaPhiDistAfter);
    fEnergyExoticClusters = new TH1F("fEnergyExoticClusters","fEnergyExoticClusters;E_{ex clus} (GeV)",1500,0,150);
    fOutput->Add(fEnergyExoticClusters);
    
    // Take ownership of output list
    fOutput->SetOwner(kTRUE);
  }
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionClusterExotics::Run()
{
  // Run
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Run();
  
  if (!fClusCont) return kFALSE;
  
  // loop over clusters
  fClusCont->ResetCurrentID();
  AliVCluster *clus = 0;
  while ((clus = fClusCont->GetNextCluster())) {
    if (!clus->IsEMCAL()) continue;
    
    if (fCreateHisto) {
      Float_t pos[3] = {0.};
      clus->GetPosition(pos);
      TVector3 vec(pos);
      fEtaPhiDistBefore->Fill(vec.Eta(),vec.Phi());
    }
    
    Bool_t exResult = kFALSE;
    
    if (fRecoUtils) {
      if (fRecoUtils->IsRejectExoticCluster()) {
        Bool_t exRemoval = fRecoUtils->IsRejectExoticCell();
        fRecoUtils->SwitchOnRejectExoticCell();                  //switch on temporarily
        exResult = fRecoUtils->IsExoticCluster(clus, fCaloCells);
        if (!exRemoval) fRecoUtils->SwitchOffRejectExoticCell(); //switch back off
        
        clus->SetIsExotic(exResult);
      }
    }
    
    if (fCreateHisto) {
      if (exResult) {
        fEnergyExoticClusters->Fill(clus->E());
      }
      else {
        Float_t pos[3] = {0.};
        clus->GetPosition(pos);
        TVector3 vec(pos);
        fEtaPhiDistAfter->Fill(vec.Eta(), vec.Phi());
      }
    }
  }
  
  return kTRUE;
}
