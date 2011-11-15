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
#include "AliForwardUtil.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliFMDEventInspector.h"
#include "AliFMDSharingFilter.h"
#include "AliFMDDensityCalculator.h"
#include "AliFMDCorrector.h"
#include "AliFMDHistCollector.h"
#include "AliESDEvent.h"
#include <TROOT.h>
#include <TAxis.h>
#include <THStack.h>
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
  fBranchNames = 
    "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,"
    "AliESDFMD.,SPDVertex.,PrimaryVertex.";
}

//____________________________________________________________________
AliForwardMultiplicityBase& 
AliForwardMultiplicityBase::operator=(const AliForwardMultiplicityBase& o)
{
  if (&o == this) return *this;
  fEnableLowFlux = o.fEnableLowFlux;
  fFirstEvent    = o.fFirstEvent;
  fCorrManager   = o.fCorrManager;
  return *this;
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
  //
  // Read corrections
  //
  //
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
  //
  // Get the ESD event. IF this is the first event, initialise
  //
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
AliForwardMultiplicityBase::MakeRingdNdeta(const TList* input, 
					   const char*  inName,
					   TList*       output,
					   const char*  outName,
					   Int_t        style) const
{
  // Make dN/deta for each ring found in the input list.  
  // 
  // A stack of all the dN/deta is also made for easy drawing. 
  // 
  // Note, that the distributions are normalised to the number of
  // observed events only - they should be corrected for 
  if (!input) return;
  TList* list = static_cast<TList*>(input->FindObject(inName));
  if (!list) { 
    AliWarning(Form("No list %s found in %s", inName, input->GetName()));
    return;
  }
  
  TList* out = new TList;
  out->SetName(outName);
  out->SetOwner();
  output->Add(out);

  THStack*     dndetaRings = new THStack("all", "dN/d#eta per ring");
  const char*  names[]     = { "FMD1I", "FMD2I", "FMD2O", "FMD3I", "FMD3O", 0 };
  const char** ptr         = names;
  
  while (*ptr) { 
    TList* thisList = new TList;
    thisList->SetOwner();
    thisList->SetName(*ptr);
    out->Add(thisList);

    TH2D* h = static_cast<TH2D*>(list->FindObject(Form("%s_cache", *ptr)));
    if (!h) { 
      AliWarning(Form("Didn't find %s_cache in %s", *ptr, list->GetName()));
      ptr++;
      continue;
    }
    TH2D* copy = static_cast<TH2D*>(h->Clone("sum"));
    copy->SetDirectory(0);
    thisList->Add(copy);
    
    TH1D* norm =static_cast<TH1D*>(h->ProjectionX("norm", 0, 0, ""));
    for (Int_t i = 1; i <= copy->GetNbinsX(); i++) { 
      for (Int_t j = 1; j <= copy->GetNbinsY(); j++) { 
	Double_t c = copy->GetBinContent(i, j);
	Double_t e = copy->GetBinError(i, j);
	Double_t a = norm->GetBinContent(i);
	copy->SetBinContent(i, j, a <= 0 ? 0 : c / a);
	copy->SetBinError(i, j, a <= 0 ? 0 : e / a);
      }
    }

    TH1D* res  =static_cast<TH1D*>(copy->ProjectionX("dndeta",1,
						     h->GetNbinsY(),"e"));
    TH1D* proj =static_cast<TH1D*>(h->ProjectionX("proj",1,h->GetNbinsY(),"e"));
    res->SetTitle(*ptr);
    res->Scale(1., "width");
    copy->Scale(1., "width");
    proj->Scale(1. / norm->GetMaximum(), "width");
    norm->Scale(1. / norm->GetMaximum());

    res->SetMarkerStyle(style);
    norm->SetDirectory(0);
    res->SetDirectory(0);
    proj->SetDirectory(0);
    thisList->Add(norm);
    thisList->Add(res);
    thisList->Add(proj);
    dndetaRings->Add(res);
    ptr++;
  }
  out->Add(dndetaRings);
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
  
  std::cout << ClassName() << ": " << GetName() << "\n" 
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
  GetSharingFilter()    .Print(option);
  GetDensityCalculator().Print(option);
  GetCorrections()      .Print(option);
  GetHistCollector()    .Print(option);
  gROOT->DecreaseDirLevel();
}


//
// EOF
//
