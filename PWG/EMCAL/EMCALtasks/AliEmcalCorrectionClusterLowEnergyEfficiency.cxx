// AliEmcalCorrectionClusterLowEnergyEfficiency
//

#include "AliEmcalCorrectionClusterLowEnergyEfficiency.h"

// #include <TList.h>
#include "AliInputEventHandler.h"

#include "AliClusterContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterLowEnergyEfficiency);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterLowEnergyEfficiency> AliEmcalCorrectionClusterLowEnergyEfficiency::reg("AliEmcalCorrectionClusterLowEnergyEfficiency");

const std::map <std::string, AliEMCALRecoUtils::NCellEfficiencyFunctions> AliEmcalCorrectionClusterLowEnergyEfficiency::fgkNCellEfficiencyFunctionMap = {
    { "kNoCorrection", AliEMCALRecoUtils::kNCeNoCorrection },
    { "kAllClusters", AliEMCALRecoUtils::kNCeAllClusters },
    { "kTestBeam", AliEMCALRecoUtils::kNCeTestBeam },
    { "kGammaAndElec", AliEMCALRecoUtils::kNCeGammaAndElec }
};

/**
 * Default constructor
 */
AliEmcalCorrectionClusterLowEnergyEfficiency::AliEmcalCorrectionClusterLowEnergyEfficiency() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterLowEnergyEfficiency"),
  fNCellDistBefore(0),
  fNCellDistAfter(0),
  fRejectNextToClus(0)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionClusterLowEnergyEfficiency::~AliEmcalCorrectionClusterLowEnergyEfficiency()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionClusterLowEnergyEfficiency::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  std::string nCellFunctStr = "";
  GetProperty("nCellFunct", nCellFunctStr);
  UInt_t nCellEffiFunct = fgkNCellEfficiencyFunctionMap.at(nCellFunctStr);

  fRejectNextToClus = false;
  GetProperty("setRejectNextToClus", fRejectNextToClus);

  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;

  if (fRecoUtils) {
    fRecoUtils->SetNCellEfficiencyFunction(nCellEffiFunct);
    fRecoUtils->InitNCellEfficiencyParam();
    fRecoUtils->Print("");
  }

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionClusterLowEnergyEfficiency::UserCreateOutputObjects()
{
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  // Create my user objects.
  if (fCreateHisto){
    fNCellDistBefore = new TH1F("hNCellDistBefore","hNCellDistBefore;N_{cells} in cluster",100,-0.5, 99.5);
    fOutput->Add(fNCellDistBefore);
    fNCellDistAfter = new TH1F("hNCellDistAfter","hNCellDistAfter;N_{cells} in cluster",100,-0.5, 99.5);
    fOutput->Add(fNCellDistAfter);

    // Take ownership of output list
    fOutput->SetOwner(kTRUE);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionClusterLowEnergyEfficiency::Run()
{
  AliEmcalCorrectionComponent::Run();

  // only apply on MC clusters
  if(!fMCEvent) return kTRUE;

  // loop over clusters
  AliVCluster *clus = 0;
  AliClusterContainer * clusCont = 0;
  TIter nextClusCont(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {

    if (!clusCont) return kFALSE;
    auto clusItCont = clusCont->all_momentum();

    for (AliClusterIterableMomentumContainer::iterator clusIterator = clusItCont.begin(); clusIterator != clusItCont.end(); ++clusIterator) {
      clus = static_cast<AliVCluster *>(clusIterator->second);

      if (!clus->IsEMCAL()) continue;

      if (fCreateHisto) {
        fNCellDistBefore->Fill(clus->GetNCells());
      }
      if (fRecoUtils) {
        if(clus->GetNCells() == 1){
          if (fRecoUtils->GetNCellEfficiencyFunction() != AliEMCALRecoUtils::kNCeNoCorrection) {
            Bool_t isAccepted = fRecoUtils->GetIsNCellCorrected(clus, (AliVCaloCells*) fEventManager.InputEvent()->GetEMCALCells());
            if ( isAccepted ) clus->SetChi2(1);
          }
        }
      }

      if (fCreateHisto) {
        // case if a cluster has only one cell but its rtificielly widened
        if(clus->Chi2() == 1. && clus->GetNCells() == 1){
          fNCellDistAfter->Fill(2);
        } else {
          fNCellDistAfter->Fill(clus->GetNCells());
        }
      }
    }
  }

  return kTRUE;
}
