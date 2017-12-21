// AliEmcalCorrectionPHOSCorrections
//

#include "AliEmcalCorrectionPHOSCorrections.h"

#include "AliClusterContainer.h"
#include "AliPHOSTenderSupply.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionPHOSCorrections);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionPHOSCorrections> AliEmcalCorrectionPHOSCorrections::reg("AliEmcalCorrectionPHOSCorrections");

/**
 * Default constructor
 */
AliEmcalCorrectionPHOSCorrections::AliEmcalCorrectionPHOSCorrections() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionPHOSCorrections"),
  fIsMC(kFALSE),
  fOptions()
{
}

/**
 * Destructor
 */
AliEmcalCorrectionPHOSCorrections::~AliEmcalCorrectionPHOSCorrections()
{
  if (fPHOSTender) {
    fPHOSTender->Delete();
  }
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionPHOSCorrections::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  GetProperty("isMC", fIsMC);
  
  std::string options = "";
  GetProperty("options", options);
  fOptions = options.c_str();
  
  Bool_t isEmbedded = kFALSE;
  GetProperty("isEmbedded", isEmbedded);
  
  // Create PHOS Tender Supply and configure it
  fPHOSTender = new AliPHOSTenderSupply("PHOStenderEmb");
  if(fIsMC) { //handle MC data
    fPHOSTender->SetMCProduction(fOptions.Data());
  }
  
  // Check if main event or embedded event clusters should be corrected (only applicable for embedding analyses)
  if (isEmbedded) {
    SetUsingInputEvent(kFALSE);
    fPHOSTender->SetTask(&fEventManager);
  }

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionPHOSCorrections::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
  
  // Nothing to be done.

}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionPHOSCorrections::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  CheckIfRunChanged();
  
  fPHOSTender->ProcessEvent();
  
  return kTRUE;
}

/**
 * This function is called if the run changes (it inherits from the base component).
 */
Bool_t AliEmcalCorrectionPHOSCorrections::CheckIfRunChanged()
{
  Bool_t runChanged = AliEmcalCorrectionComponent::CheckIfRunChanged();
  
  if (runChanged) {
    
    TString passNum(fFilepass(4,5));
    Int_t pass = passNum.Atoi();
    fPHOSTender->SetReconstructionPass(pass);
    
    fPHOSTender->InitTender();
  }
  return runChanged;
}
