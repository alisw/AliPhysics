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

/**
 * Default constructor
 */
AliEmcalCorrectionClusterExotics::AliEmcalCorrectionClusterExotics() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterExotics"),
  fEtaPhiDistBefore(0),
  fEtaPhiDistAfter(0),
  fEnergyExoticClusters(0),
  fExoticMinCellAmplitude(0),
  fMaxFcross(0),
  fCellCrossMaxTimeDiff(0)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionClusterExotics::~AliEmcalCorrectionClusterExotics()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionClusterExotics::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  GetProperty("createHistos", fCreateHisto);
  
  Float_t fExoticMinCellAmplitude = 4.;
  GetProperty("fExoticMinCellAmplitude", fExoticMinCellAmplitude);
  Float_t fMaxFcross = 0.97;
  GetProperty("fMaxFcross", fMaxFcross);
  Float_t fCellCrossMaxTimeDiff = 1e6;
  GetProperty("fCellCrossMaxTimeDiff", fCellCrossMaxTimeDiff);
  
  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;
  fRecoUtils->SwitchOnRejectExoticCluster();
  fRecoUtils->SetExoticCellMinAmplitudeCut(fExoticMinCellAmplitude);
  fRecoUtils->SetExoticCellFractionCut(fMaxFcross);
  fRecoUtils->SetExoticCellDiffTimeCut(fCellCrossMaxTimeDiff);
  if (fRecoUtils)
    fRecoUtils->Print("");

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionClusterExotics::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
  
  // Create my user objects.
  if (fCreateHisto){
    fEtaPhiDistBefore = new TH2F("hEtaPhiDistBefore","hEtaPhiDistBefore;#eta;#phi",280,-0.7,0.7,200*3.14,0,2*3.14);
    fOutput->Add(fEtaPhiDistBefore);
    fEtaPhiDistAfter = new TH2F("hEtaPhiDistAfter","hEtaPhiDistAfter;#eta;#phi",280,-0.7,0.7,200*3.14,0,2*3.14);
    fOutput->Add(fEtaPhiDistAfter);
    fEnergyExoticClusters = new TH1F("fEnergyExoticClusters","fEnergyExoticClusters;E_{ex clus} (GeV)",1500,0,150);
    fOutput->Add(fEnergyExoticClusters);
    
    // Take ownership of output list
    fOutput->SetOwner(kTRUE);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionClusterExotics::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  // loop over clusters
  AliVCluster *clus = 0;
  AliClusterContainer * clusCont = 0;
  TIter nextClusCont(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {

    if (!clusCont) continue;
    auto clusItCont = clusCont->all_momentum();

    for (AliClusterIterableMomentumContainer::iterator clusIterator = clusItCont.begin(); clusIterator != clusItCont.end(); ++clusIterator) {
      clus = static_cast<AliVCluster *>(clusIterator->second);

      if (!clus->IsEMCAL()) continue;

      if (fCreateHisto) {
        Float_t pos[3] = {0.};
        clus->GetPosition(pos);
        TVector3 vec(pos);
        // Phi needs to be in 0 to 2 Pi
        fEtaPhiDistBefore->Fill(vec.Eta(), TVector2::Phi_0_2pi(vec.Phi()));
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
          // Phi needs to be in 0 to 2 Pi
          fEtaPhiDistAfter->Fill(vec.Eta(), TVector2::Phi_0_2pi(vec.Phi()));
        }
      }
    }
  }
  
  return kTRUE;
}
