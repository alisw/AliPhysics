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
#include "AliFMDEventPlaneFinder.h"
#include "AliESDEvent.h"
#include <TROOT.h>
#include <TSystem.h>
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
  DGUARD(fDebug, 3,"Named CTOR of AliForwardMultiplicityBase %s",name);
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
  DGUARD(fDebug,2,"Assignment to AliForwardMultiplicityBase");
  if (&o == this) return *this;
  fEnableLowFlux = o.fEnableLowFlux;
  fFirstEvent    = o.fFirstEvent;
  fCorrManager   = o.fCorrManager;
  return *this;
}
//____________________________________________________________________
Bool_t 
AliForwardMultiplicityBase::Configure(const char* macro)
{
  // --- Configure the task ------------------------------------------
  TString macroPath(gROOT->GetMacroPath());
  if (!macroPath.Contains("$(ALICE_ROOT)/PWGLF/FORWARD/analysis2")) { 
    macroPath.Append(":$(ALICE_ROOT)/PWGLF/FORWARD/analysis2");
    gROOT->SetMacroPath(macroPath);
  }
  const char* config = gSystem->Which(gROOT->GetMacroPath(), macro);
  if (!config) {
    AliWarningF("%s not found in %s", macro, gROOT->GetMacroPath());
    return false;
  }

  AliInfoF("Loading configuration of '%s' from %s",  ClassName(), config);
  gROOT->Macro(Form("%s((AliForwardMultiplicityBase*)%p)", config, this));
  delete config;
 
 return true;
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
  DGUARD(fDebug,1,"Checking corrections 0x%x", what);

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
  DGUARD(fDebug,1,"Read corrections 0x%x", what);

  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  if (!fcm.Init(GetEventInspector().GetCollisionSystem(),
		GetEventInspector().GetEnergy(),
		GetEventInspector().GetField(),
		mc,
		what)) { 
    AliWarning("Failed to read in some corrections, making task zombie");
    return false;
  }
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
  DGUARD(fDebug,1,"Get the ESD event");
  if (IsZombie()) return 0;
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliWarning("No ESD event found for input event");
    return 0;
  }

  // --- Load the data -----------------------------------------------
  LoadBranches();

  // On the first event, initialize the parameters
  if (fFirstEvent && esd->GetESDRun()) {
    GetEventInspector().ReadRunDetails(esd);
    
    AliInfo(Form("Initializing with parameters from the ESD:\n"
		 "         AliESDEvent::GetBeamEnergy()   ->%f\n"
		 "         AliESDEvent::GetBeamType()     ->%s\n"
		 "         AliESDEvent::GetCurrentL3()    ->%f\n"
		 "         AliESDEvent::GetMagneticField()->%f\n"
		 "         AliESDEvent::GetRunNumber()    ->%d",
		 esd->GetBeamEnergy(),
		 esd->GetBeamType(),
		 esd->GetCurrentL3(),
		 esd->GetMagneticField(),
		 esd->GetRunNumber()));
    
    fFirstEvent = false;
    
    GetEventPlaneFinder().SetRunNumber(esd->GetRunNumber());
    if (!SetupForData()) { 
      AliError("Failed to initialize sub-algorithms, making this a zombie");
      esd = 0; // Make sure we do nothing on this event
      Info("GetESDEvent", "ESD event pointer %p", esd);
      SetZombie(true);
      // return 0;
    }
  }
  return esd;
}
//____________________________________________________________________
void
AliForwardMultiplicityBase::MarkEventForStore() const
{
  // Make sure the AOD tree is filled 
  DGUARD(fDebug,3,"Mark AOD event for storage");
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (ah)    	
	ah->SetFillAOD(kTRUE);

}

//____________________________________________________________________
Bool_t
AliForwardMultiplicityBase::MakeSimpledNdeta(const TList* input, 
					     TList*       output,
					     Double_t&    nTr, 
					     Double_t&    nTrVtx, 
					     Double_t&    nAcc)
{
  // Get our histograms from the container 
  TH1I* hEventsTr    = 0;
  TH1I* hEventsTrVtx = 0;
  TH1I* hEventsAcc   = 0;
  TH1I* hTriggers    = 0;
  if (!GetEventInspector().FetchHistograms(input, 
					   hEventsTr, 
					   hEventsTrVtx, 
					   hEventsAcc,
					   hTriggers)) { 
    AliError(Form("Didn't get histograms from event selector "
		  "(hEventsTr=%p,hEventsTrVtx=%p,hEventsAcc=%p,hTriggers=%p)", 
		  hEventsTr, hEventsTrVtx, hEventsAcc, hTriggers));
    input->ls();
    return false;
  }
  nTr             = hEventsTr->Integral();
  nTrVtx          = hEventsTrVtx->Integral();
  nAcc            = hEventsAcc->Integral();
  Double_t vtxEff = nTrVtx / nTr;
  TH2D*   hData   = static_cast<TH2D*>(input->FindObject("d2Ndetadphi"));
  if (!hData) { 
    AliError(Form("Couldn't get our summed histogram from output "
		  "list %s (d2Ndetadphi=%p)", input->GetName(), hData));
    input->ls();
    return false;
  }

  Int_t nY      = hData->GetNbinsY();
  TH1D* dNdeta  = hData->ProjectionX("dNdeta",  1,     nY, "e");
  TH1D* dNdeta_ = hData->ProjectionX("dNdeta_", 1,     nY, "e");
  TH1D* norm    = hData->ProjectionX("norm",    0,     0,  "");
  TH1D* phi     = hData->ProjectionX("phi",     nY+1,  nY+1,  "");
  dNdeta->SetTitle("dN_{ch}/d#eta in the forward regions");
  dNdeta->SetYTitle("#frac{1}{N}#frac{dN_{ch}}{d#eta}");
  dNdeta->SetMarkerColor(kRed+1);
  dNdeta->SetMarkerStyle(20);
  dNdeta->SetDirectory(0);

  dNdeta_->SetTitle("dN_{ch}/d#eta in the forward regions");
  dNdeta_->SetYTitle("#frac{1}{N}#frac{dN_{ch}}{d#eta}");
  dNdeta_->SetMarkerColor(kMagenta+1);
  dNdeta_->SetMarkerStyle(21);
  dNdeta_->SetDirectory(0);

  norm->SetTitle("Normalization to #eta coverage");
  norm->SetYTitle("#eta coverage");
  norm->SetMarkerColor(kBlue+1);
  norm->SetMarkerStyle(21);
  norm->SetFillColor(kBlue+1);
  norm->SetFillStyle(3005);
  norm->SetDirectory(0);

  phi->SetTitle("Normalization to #phi acceptance");
  phi->SetYTitle("#phi acceptance");
  phi->SetMarkerColor(kGreen+1);
  phi->SetMarkerStyle(20);
  phi->SetFillColor(kGreen+1);
  phi->SetFillStyle(3004);
  // phi->Scale(1. / nAcc);
  phi->SetDirectory(0);

  // dNdeta->Divide(norm);
  dNdeta->Divide(phi);
  dNdeta->SetStats(0);
  dNdeta->Scale(vtxEff,	"width");

  dNdeta_->Divide(norm);
  dNdeta_->SetStats(0);
  dNdeta_->Scale(vtxEff, "width");

  output->Add(dNdeta);
  output->Add(dNdeta_);
  output->Add(norm);
  output->Add(phi);

  return true;
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
  DGUARD(fDebug,3,"Make first-shot ring dN/deta");

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
    TH2D* sumPhi = static_cast<TH2D*>(h->Clone("sum_phi"));
    sumPhi->SetDirectory(0);
    thisList->Add(sumPhi);

    TH2D* sumEta = static_cast<TH2D*>(h->Clone("sum_eta"));
    sumEta->SetDirectory(0);
    thisList->Add(sumEta);
    
    Int_t nY   = sumEta->GetNbinsY();
    TH1D* etaCov =static_cast<TH1D*>(h->ProjectionX("etaCov", 0,    0,    ""));
    TH1D* phiAcc =static_cast<TH1D*>(h->ProjectionX("phiAcc", nY+1, nY+1, ""));

    etaCov->SetTitle("Normalization to #eta coverage");
    etaCov->SetYTitle("#eta coverage");
    etaCov->SetMarkerColor(kBlue+1);
    etaCov->SetFillColor(kBlue+1);
    etaCov->SetFillStyle(3005);
    etaCov->SetDirectory(0);
    
    phiAcc->SetTitle("Normalization to #phi acceptance");
    phiAcc->SetYTitle("#phi acceptance");
    phiAcc->SetMarkerColor(kGreen+1);
    phiAcc->SetFillColor(kGreen+1);
    phiAcc->SetFillStyle(3004);
    // phiAcc->Scale(1. / nAcc);
    phiAcc->SetDirectory(0);

    // Double_t s = (etaCov->GetMaximum() > 0 ? 1. / etaCov->GetMaximum() : 1);
    for (Int_t i = 1; i <= sumEta->GetNbinsX(); i++) { 
      for (Int_t j = 1; j <= nY; j++) { 
	Double_t c = sumEta->GetBinContent(i, j);
	Double_t e = sumEta->GetBinError(i, j);
	Double_t a = etaCov->GetBinContent(i);
	Double_t p = phiAcc->GetBinContent(i);
	// Double_t t = p; // * a
	sumEta->SetBinContent(i, j, a <= 0 ? 0 : c / a);
	sumEta->SetBinError(  i, j, a <= 0 ? 0 : e / a);
	sumPhi->SetBinContent(i, j, p <= 0 ? 0 : c / p);
	sumPhi->SetBinError(  i, j, p <= 0 ? 0 : e / p);
      }
    }
    // etaCov->Scale(s);
    // phiAcc->Scale(s);

    TH1D* resPhi  =static_cast<TH1D*>(sumPhi->ProjectionX("dndeta_phi",
							  1,nY,"e"));
    resPhi->SetMarkerStyle(style);
    resPhi->SetDirectory(0);
    resPhi->Scale(1, "width");

    TH1D* resEta  =static_cast<TH1D*>(sumEta->ProjectionX("dndeta_eta",
							  1,nY,"e"));
    resEta->SetMarkerStyle(style);
    resEta->SetDirectory(0);
    resEta->Scale(1, "width");

    thisList->Add(resEta);
    thisList->Add(etaCov);
    thisList->Add(resPhi);
    thisList->Add(phiAcc);
    dndetaRings->Add(resPhi);
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
  GetEventPlaneFinder() .Print(option);
  gROOT->DecreaseDirLevel();
}


//
// EOF
//
