// AliEmcalCorrectionClusterPositionCorrection
//

#include "AliEmcalCorrectionClusterPositionCorrection.h"

#include "AliInputEventHandler.h"
#include "AliClusterContainer.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterPositionCorrection);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterPositionCorrection> AliEmcalCorrectionClusterPositionCorrection::reg("AliEmcalCorrectionClusterPositionCorrection");

/**
 * Default constructor
 */
AliEmcalCorrectionClusterPositionCorrection::AliEmcalCorrectionClusterPositionCorrection() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterPositionCorrection"),
  fClusterPositionDiffEta(0),
  fClusterPositionDiffPhi(0),
  fSMEtaValues(0),
  fSMPhiValues(0),
  fApplyToMC(0)
{
  for(unsigned int i = 0; i < 20; ++i){
    fSMEtaValues.push_back(0); // set default values to 0 --> no shift is applied
    fSMPhiValues.push_back(0); // set default values to 0 --> no shift is applied
  }
}

/**
 * Destructor
 */
AliEmcalCorrectionClusterPositionCorrection::~AliEmcalCorrectionClusterPositionCorrection()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionClusterPositionCorrection::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();

  fApplyToMC = false;
  GetProperty("applyToMC", fApplyToMC);

  GetProperty("setSMEtaValues", fSMEtaValues);

  GetProperty("setSMPhiValues", fSMPhiValues);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionClusterPositionCorrection::UserCreateOutputObjects()
{
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  // Create my user objects.
  if (fCreateHisto){
    fClusterPositionDiffEta = new TH1F("fClusterPositionDiffEta","fClusterPositionDiffEta;#Delta#eta",400,-0.02, 0.02);
    fOutput->Add(fClusterPositionDiffEta);
    fClusterPositionDiffPhi = new TH1F("fClusterPositionDiffPhi","fClusterPositionDiffPhi;#Delta#eta",400,-0.02, 0.02);
    fOutput->Add(fClusterPositionDiffPhi);

    // Take ownership of output list
    fOutput->SetOwner(kTRUE);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionClusterPositionCorrection::Run()
{
  AliEmcalCorrectionComponent::Run();

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

      Float_t pos[3] = {0};
      clus->GetPosition(pos);
      // printf("cluster pos: x=%.04f y=%.04f z=%.04f\n",pos[0],pos[1],pos[2]);

      TVector3 vecPos(pos);

      Double_t clusPt = vecPos.Pt();
      Double_t clusEta = vecPos.Eta();
      Double_t clusPhi = vecPos.Phi();
      Int_t iSM = fGeom->GetSuperModuleNumber(clus->GetCellAbsId(0));
      if(fApplyToMC){
        clusEta-= fSMEtaValues[iSM];
        clusPhi-= fSMPhiValues[iSM];
      } else {
        clusEta+= fSMEtaValues[iSM];
        clusPhi+= fSMPhiValues[iSM];
      }
      // shift the position of the cluster
      vecPos.SetPtEtaPhi(clusPt, clusEta, clusPhi);

      // Set position after correction
      vecPos.GetXYZ(pos);
      clus->SetPosition(pos);
      // printf("TVector pos: x=%.04f y=%.04f z=%.04f\n",pos[0],pos[1],pos[2]);

      Double_t clusEtaAfter = vecPos.Eta();
      Double_t clusPhiAfter = vecPos.Phi();

      if (fCreateHisto) {
        fClusterPositionDiffEta->Fill(clusEtaAfter - clusEta);
        fClusterPositionDiffPhi->Fill(clusPhiAfter - clusPhi);
      }
    }
  }

  return kTRUE;
}
