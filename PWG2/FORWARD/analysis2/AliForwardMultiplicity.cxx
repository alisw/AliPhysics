#include "AliForwardMultiplicity.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "AliFMDAnaParameters.h"
#include "AliESDEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include <TH1.h>
#include <TDirectory.h>
#include <TTree.h>

//====================================================================
AliForwardMultiplicity::AliForwardMultiplicity()
  : AliAnalysisTaskSE(),
    fHData(0),
    fFirstEvent(true),
    fESDFMD(),
    fHistos(),
    fAODFMD(),
    fEventInspector(),
    fEnergyFitter(),
    fSharingFilter(),
    fDensityCalculator(),
    fCorrections(),
    fHistCollector(),
    fList(0)
{
}

//____________________________________________________________________
AliForwardMultiplicity::AliForwardMultiplicity(const char* name)
  : AliAnalysisTaskSE(name), 
    fHData(0),
    fFirstEvent(true),
    fESDFMD(),
    fHistos(),
    fAODFMD(kTRUE),
    fEventInspector("event"),
    fEnergyFitter("energy"),
    fSharingFilter("sharing"), 
    fDensityCalculator("density"),
    fCorrections("corrections"),
    fHistCollector("collector"),
    fList(0)
{
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________
AliForwardMultiplicity::AliForwardMultiplicity(const AliForwardMultiplicity& o)
  : AliAnalysisTaskSE(o),
    fHData(o.fHData),
    fFirstEvent(true),
    fESDFMD(o.fESDFMD),
    fHistos(o.fHistos),
    fAODFMD(o.fAODFMD),
    fEventInspector(o.fEventInspector),
    fEnergyFitter(o.fEnergyFitter),
    fSharingFilter(o.fSharingFilter),
    fDensityCalculator(o.fDensityCalculator),
    fCorrections(o.fCorrections),
    fHistCollector(o.fHistCollector),
    fList(o.fList) 
{
}

//____________________________________________________________________
AliForwardMultiplicity&
AliForwardMultiplicity::operator=(const AliForwardMultiplicity& o)
{
  AliAnalysisTaskSE::operator=(o);

  fHData             = o.fHData;
  fFirstEvent        = o.fFirstEvent;
  fEventInspector    = o.fEventInspector;
  fEnergyFitter      = o.fEnergyFitter;
  fSharingFilter     = o.fSharingFilter;
  fDensityCalculator = o.fDensityCalculator;
  fCorrections       = o.fCorrections;
  fHistCollector     = o.fHistCollector;
  fHistos            = o.fHistos;
  fAODFMD            = o.fAODFMD;
  fList              = o.fList;

  return *this;
}

//____________________________________________________________________
void
AliForwardMultiplicity::SetDebug(Int_t dbg)
{
  fEventInspector.SetDebug(dbg);
  fEnergyFitter.SetDebug(dbg);
  fSharingFilter.SetDebug(dbg);
  fDensityCalculator.SetDebug(dbg);
  fCorrections.SetDebug(dbg);
  fHistCollector.SetDebug(dbg);
}

//____________________________________________________________________
void
AliForwardMultiplicity::Init()
{
  fFirstEvent = true;
}

//____________________________________________________________________
void
AliForwardMultiplicity::InitializeSubs()
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  pars->Init(kTRUE);


  TAxis e(pars->GetNetaBins(),  pars->GetEtaMin(),  pars->GetEtaMax());
  TAxis v(pars->GetNvtxBins(), -pars->GetVtxCutZ(), pars->GetVtxCutZ());
			
  fHistos.Init(e);
  fAODFMD.Init(e);

  fHData = static_cast<TH2D*>(fAODFMD.GetHistogram().Clone("d2Ndetadphi"));
  fHData->SetStats(0);
  fHData->SetDirectory(0);
  fList->Add(fHData);

  fEnergyFitter.Init(e);
  fEventInspector.Init(v);
  fHistCollector.Init(v);
}

//____________________________________________________________________
void
AliForwardMultiplicity::UserCreateOutputObjects()
{
  fList = new TList;

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) AliFatal("No AOD output handler set in analysis manager");
    
    
  TObject* obj = &fAODFMD;
  ah->AddBranch("AliAODForwardMult", &obj);

  fEventInspector.DefineOutput(fList);
  fEnergyFitter.DefineOutput(fList);
  fSharingFilter.DefineOutput(fList);
  fDensityCalculator.DefineOutput(fList);
  fCorrections.DefineOutput(fList);
}
//____________________________________________________________________
void
AliForwardMultiplicity::UserExec(Option_t*)
{
  // Get the input data 
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) { 
    AliWarning("No ESD event found for input event");
    return;
  }

  // On the first event, initialize the parameters 
  if (fFirstEvent && esd->GetESDRun()) { 
    AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
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
    pars->SetParametersFromESD(esd);
    pars->PrintStatus();
    fFirstEvent = false;

    InitializeSubs();
  }
  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  fAODFMD.Clear();

  Bool_t   lowFlux  = kFALSE;
  UInt_t   triggers = 0;
  Int_t    ivz      = -1;
  Double_t vz       = 0;
  UInt_t   found    = fEventInspector.Process(esd, triggers, lowFlux, ivz, vz);
  if (found & AliFMDEventInspector::kNoEvent)    return;
  if (found & AliFMDEventInspector::kNoTriggers) return;
 
  // Set trigger bits, and mark this event for storage 
  fAODFMD.SetTriggerBits(triggers);
  MarkEventForStore();

  if (found & AliFMDEventInspector::kNoSPD)     return;
  if (found & AliFMDEventInspector::kNoFMD)     return;
  if (found & AliFMDEventInspector::kNoVertex)  return;
  fAODFMD.SetIpZ(vz);

  if (found & AliFMDEventInspector::kBadVertex) return;

  // Get FMD data 
  AliESDFMD* esdFMD = esd->GetFMDData();
  // Apply the sharing filter (or hit merging or clustering if you like)
  if (!fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD)) { 
    AliWarning("Sharing filter failed!");
    return;
  }

  // Do the energy stuff 
  if (!fEnergyFitter.Accumulate(*esdFMD, triggers & AliAODForwardMult::kEmpty)){
    AliWarning("Energy fitter failed");
    return;
  }

  // Calculate the inclusive charged particle density 
  if (!fDensityCalculator.Calculate(fESDFMD, fHistos, ivz, lowFlux)) { 
    AliWarning("Density calculator failed!");
    return;
  }
  
  // Do the secondary and other corrections. 
  if (!fCorrections.Correct(fHistos, ivz)) { 
    AliWarning("Corrections failed");
    return;
  }

  if (!fHistCollector.Collect(fHistos, ivz, fAODFMD.GetHistogram())) {
    AliWarning("Histogram collector failed");
    return;
  }

  if (fAODFMD.IsTriggerBits(AliAODForwardMult::kInel))
    fHData->Add(&(fAODFMD.GetHistogram()));

  PostData(1, fList);
}

//____________________________________________________________________
void
AliForwardMultiplicity::Terminate(Option_t*)
{
  TList* list = dynamic_cast<TList*>(GetOutputData(1));
  if (!list) {
    AliError(Form("No output list defined (%p)", GetOutputData(1)));
    if (GetOutputData(1)) GetOutputData(1)->Print();
    return;
  }
  
  // Get our histograms from the container 
  TH1I* hEventsTr    = 0;//static_cast<TH1I*>(list->FindObject("nEventsTr"));
  TH1I* hEventsTrVtx = 0;//static_cast<TH1I*>(list->FindObject("nEventsTrVtx"));
  TH1I* hTriggers    = 0;
  if (!fEventInspector.FetchHistograms(list, hEventsTr, 
				       hEventsTrVtx, hTriggers)) { 
    AliError(Form("Didn't get histograms from event selector "
		  "(hEventsTr=%p,hEventsTrVtx=%p)", 
		  hEventsTr, hEventsTrVtx));
    list->ls();
    return;
  }

  TH2D* hData        = static_cast<TH2D*>(list->FindObject("d2Ndetadphi"));
  if (!hData) { 
    AliError(Form("Couldn't get our summed histogram from output "
		  "list %s (d2Ndetadphi=%p)", list->GetName(), hData));
    list->ls();
    return;
  }
  
  // TH1D* dNdeta = fHData->ProjectionX("dNdeta", 0, -1, "e");
  TH1D* dNdeta = hData->ProjectionX("dNdeta", 1, -1, "e");
  TH1D* norm   = hData->ProjectionX("dNdeta", 0, 1,  "");
  dNdeta->SetTitle("dN_{ch}/d#eta in the forward regions");
  dNdeta->SetYTitle("#frac{1}{N}#frac{dN_{ch}}{d#eta}");
  dNdeta->Divide(norm);
  dNdeta->SetStats(0);
  dNdeta->Scale(Double_t(hEventsTrVtx->GetEntries())/hEventsTr->GetEntries(),
		"width");
  list->Add(dNdeta);

  fEnergyFitter.Fit(list);
  fSharingFilter.ScaleHistograms(list,hEventsTr->Integral());
  fDensityCalculator.ScaleHistograms(list,hEventsTrVtx->Integral());
  fCorrections.ScaleHistograms(list,hEventsTrVtx->Integral());
}

//____________________________________________________________________
void
AliForwardMultiplicity::MarkEventForStore() const
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
AliForwardMultiplicity::Print(Option_t*) const
{
}

//
// EOF
//
