//====================================================================
// 
// Base class for classes that calculate the multiplicity in the
// forward regions event-by-event
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODForwardMult 
// 
// Histograms 
//   
// Corrections used 
#include "AliForwardMultiplicityBase.h"
#include "AliForwardCorrectionManager.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliFMDEventInspector.h"
#include "AliFMDEnergyFitter.h"
#include "AliFMDSharingFilter.h"
#include "AliFMDDensityCalculator.h"
#include "AliFMDCorrector.h"
#include "AliFMDHistCollector.h"
#include "AliESDEvent.h"
#include <TROOT.h>
#include <iostream>
#include <iomanip>

//====================================================================
AliForwardMultiplicityBase::AliForwardMultiplicityBase(const char* name) 
  : AliAnalysisTaskSE(name), 
    fEnableLowFlux(false), 
    fFirstEvent(true),
    fCorrManager(0)
{
  // Set our persistent pointer 
  fCorrManager = &AliForwardCorrectionManager::Instance();
}

//____________________________________________________________________
Bool_t 
AliForwardMultiplicityBase::CheckCorrections(UInt_t what) const
{
  // 
  // Check if all needed corrections are there and accounted for.  If not,
  // do a Fatal exit 
  // 
  // Parameters:
  //    what Which corrections is needed
  // 
  // Return:
  //    true if all present, false otherwise
  //  

  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  // Check that we have the energy loss fits, needed by 
  //   AliFMDSharingFilter 
  //   AliFMDDensityCalculator 
  if (what & AliForwardCorrectionManager::kELossFits && !fcm.GetELossFit()) { 
    AliFatal(Form("No energy loss fits"));
    return false;
  }
  // Check that we have the double hit correction - (optionally) used by 
  //  AliFMDDensityCalculator 
  if (what & AliForwardCorrectionManager::kDoubleHit && !fcm.GetDoubleHit()) {
    AliFatal("No double hit corrections"); 
    return false;
  }
  // Check that we have the secondary maps, needed by 
  //   AliFMDCorrector 
  //   AliFMDHistCollector
  if (what & AliForwardCorrectionManager::kSecondaryMap && 
      !fcm.GetSecondaryMap()) {
    AliFatal("No secondary corrections");
    return false;
  }
  // Check that we have the vertex bias correction, needed by 
  //   AliFMDCorrector 
  if (what & AliForwardCorrectionManager::kVertexBias && 
      !fcm.GetVertexBias()) { 
    AliFatal("No event vertex bias corrections");
    return false;
  }
  // Check that we have the merging efficiencies, optionally used by 
  //   AliFMDCorrector 
  if (what & AliForwardCorrectionManager::kMergingEfficiency && 
      !fcm.GetMergingEfficiency()) {
    AliFatal("No merging efficiencies");
    return false;
  }
  // Check that we have the acceptance correction, needed by 
  //   AliFMDCorrector 
  if (what & AliForwardCorrectionManager::kAcceptance && 
      !fcm.GetAcceptance()) { 
    AliFatal("No acceptance corrections");
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t
AliForwardMultiplicityBase::ReadCorrections(const TAxis*& pe, 
					    const TAxis*& pv, 
					    Bool_t        mc)
{
  UInt_t what = AliForwardCorrectionManager::kAll;
  if (!fEnableLowFlux)
    what ^= AliForwardCorrectionManager::kDoubleHit;
  if (!GetCorrections().IsUseVertexBias())
    what ^= AliForwardCorrectionManager::kVertexBias;
  if (!GetCorrections().IsUseAcceptance())
    what ^= AliForwardCorrectionManager::kAcceptance;
  if (!GetCorrections().IsUseMergingEfficiency())
    what ^= AliForwardCorrectionManager::kMergingEfficiency;

  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  if (!fcm.Init(GetEventInspector().GetCollisionSystem(),
		GetEventInspector().GetEnergy(),
		GetEventInspector().GetField(),
		mc,
		what)) return false;
  if (!CheckCorrections(what)) return false;

  // Sett our persistency pointer 
  // fCorrManager = &fcm;

  // Get the eta axis from the secondary maps - if read in
  if (!pe) {
    pe = fcm.GetEtaAxis();
    if (!pe) AliFatal("No eta axis defined");
  }
  // Get the vertex axis from the secondary maps - if read in
  if (!pv) {
    pv = fcm.GetVertexAxis();
    if (!pv) AliFatal("No vertex axis defined");
  }

  return true;
}
//____________________________________________________________________
AliESDEvent*
AliForwardMultiplicityBase::GetESDEvent()
{
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliWarning("No ESD event found for input event");
    return 0;
  }

  // On the first event, initialize the parameters
  if (fFirstEvent && esd->GetESDRun()) {
    GetEventInspector().ReadRunDetails(esd);

    AliInfo(Form("Initializing with parameters from the ESD:\n"
                 "         AliESDEvent::GetBeamEnergy()   ->%f\n"
                 "         AliESDEvent::GetBeamType()     ->%s\n"
                 "         AliESDEvent::GetCurrentL3()    ->%f\n"
                 "         AliESDEvent::GetMagneticField()->%f\n"
                 "         AliESDEvent::GetRunNumber()    ->%d\n",
                 esd->GetBeamEnergy(),
                 esd->GetBeamType(),
                 esd->GetCurrentL3(),
                 esd->GetMagneticField(),
                 esd->GetRunNumber()));

    fFirstEvent = false;

    InitializeSubs();
  }
  return esd;
}
//____________________________________________________________________
void
AliForwardMultiplicityBase::MarkEventForStore() const
{
  // Make sure the AOD tree is filled 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah)  
    AliFatal("No AOD output handler set in analysis manager");

  ah->SetFillAOD(kTRUE);
}

//____________________________________________________________________
void
AliForwardMultiplicityBase::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  
  std::cout << "AliForwardMultiplicityBase: " << GetName() << "\n" 
	    << "  Enable low flux code:   " << (fEnableLowFlux ? "yes" : "no") 
	    << "\n"
	    << "  Off-line trigger mask:  0x" 
	    << std::hex     << std::setfill('0') 
	    << std::setw (8) << fOfflineTriggerMask 
	    << std::dec     << std::setfill (' ') << std::endl;
  gROOT->IncreaseDirLevel();
  if (fCorrManager) fCorrManager->Print();
  else  
    std::cout << "  Correction manager not set yet" << std::endl;
  GetEventInspector()   .Print(option);
  GetEnergyFitter()     .Print(option);    
  GetSharingFilter()    .Print(option);
  GetDensityCalculator().Print(option);
  GetCorrections()      .Print(option);
  GetHistCollector()    .Print(option);
  gROOT->DecreaseDirLevel();
}


//
// EOF
//
