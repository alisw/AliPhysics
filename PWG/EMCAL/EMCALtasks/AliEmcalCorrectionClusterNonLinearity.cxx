// AliEmcalCorrectionClusterNonLinearity
//

#include "AliEmcalCorrectionClusterNonLinearity.h"

#include <TList.h>

#include "AliClusterContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterNonLinearity);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterNonLinearity> AliEmcalCorrectionClusterNonLinearity::reg("AliEmcalCorrectionClusterNonLinearity");

//________________________________________________________________________
AliEmcalCorrectionClusterNonLinearity::AliEmcalCorrectionClusterNonLinearity() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterNonLinearity"),
  fEnergyDistBefore(0),
  fEnergyTimeHistBefore(0),
  fEnergyDistAfter(0),
  fEnergyTimeHistAfter(0)

{
  // Default constructor
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  
}

//________________________________________________________________________
AliEmcalCorrectionClusterNonLinearity::~AliEmcalCorrectionClusterNonLinearity()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionClusterNonLinearity::Initialize()
{
  // Initialization
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Initialize();
  // Do base class initializations and if it fails -> bail out
  //AliAnalysisTaskEmcal::ExecOnce();
  //if (!fInitialized) return;
  
  GetProperty("createHistos", fCreateHisto);

  std::string nonLinFunctStr = "";
  GetProperty("nonLinFunct", nonLinFunctStr);
  UInt_t nonLinFunct = nonlinearityFunctionMap[nonLinFunctStr];

  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;
  fRecoUtils->SetNonLinearityFunction(nonLinFunct);
  
  if (fRecoUtils) {
    fRecoUtils->InitNonLinearityParam();
    fRecoUtils->Print("");
  }

  return kTRUE;
}

//________________________________________________________________________
void AliEmcalCorrectionClusterNonLinearity::UserCreateOutputObjects()
{   
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
  
  // Create my user objects.
  if (fCreateHisto){
    fEnergyDistBefore = new TH1F("hEnergyDistBefore","hEnergyDistBefore;E_{clus} (GeV)",1500,0,150);
    fOutput->Add(fEnergyDistBefore);
    fEnergyTimeHistBefore = new TH2F("hEnergyTimeDistBefore","hEnergyTimeDistBefore;E_{clus} (GeV);time",1500,0,150,500,0,1e-6);
    fOutput->Add(fEnergyTimeHistBefore);
    fEnergyDistAfter = new TH1F("hEnergyDistAfter","hEnergyDistAfter;E_{clus} (GeV)",1500,0,150);
    fOutput->Add(fEnergyDistAfter);
    fEnergyTimeHistAfter = new TH2F("hEnergyTimeDistAfter","hEnergyTimeDistAfter;E_{clus} (GeV);time",1500,0,150,500,0,1e-6);
    fOutput->Add(fEnergyTimeHistAfter);
    
    // Take ownership of output list
    fOutput->SetOwner(kTRUE);
  }
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionClusterNonLinearity::Run()
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
      fEnergyDistBefore->Fill(clus->E());
      fEnergyTimeHistBefore->Fill(clus->E(), clus->GetTOF());
    }
    
    if (fRecoUtils) {
      if (fRecoUtils->GetNonLinearityFunction() != AliEMCALRecoUtils::kNoCorrection) {
        Double_t energy = fRecoUtils->CorrectClusterEnergyLinearity(clus);
        clus->SetNonLinCorrEnergy(energy);
      }
    }
    
    // Fill histograms only if cluster is not exotic, as in ClusterMaker (the clusters are flagged, not removed)
    if (fCreateHisto && !clus->GetIsExotic()) {
      fEnergyDistAfter->Fill(clus->GetNonLinCorrEnergy());
      fEnergyTimeHistAfter->Fill(clus->GetNonLinCorrEnergy(), clus->GetTOF());
    }
  }
  
  return kTRUE;
}
