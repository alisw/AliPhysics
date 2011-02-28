// 
// Calculate the multiplicity in the forward regions event-by-event 
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
//
#include "AliForwardMultiplicityTask.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliForwardCorrectionManager.h"
#include "AliAnalysisManager.h"
#include <TH1.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TROOT.h>


//====================================================================
AliForwardMultiplicityTask::AliForwardMultiplicityTask()
  : AliForwardMultiplicityBase(),
    fHData(0),
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
  // 
  // Constructor
  //
}

//____________________________________________________________________
AliForwardMultiplicityTask::AliForwardMultiplicityTask(const char* name)
  : AliForwardMultiplicityBase(name),
    fHData(0),
    fESDFMD(),
    fHistos(),
    fAODFMD(false),
    fEventInspector("event"),
    fEnergyFitter("energy"),
    fSharingFilter("sharing"), 
    fDensityCalculator("density"),
    fCorrections("corrections"),
    fHistCollector("collector"),
    fList(0)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________
AliForwardMultiplicityTask::AliForwardMultiplicityTask(const AliForwardMultiplicityTask& o)
  : AliForwardMultiplicityBase(o),
    fHData(o.fHData),
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
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________
AliForwardMultiplicityTask&
AliForwardMultiplicityTask::operator=(const AliForwardMultiplicityTask& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object 
  //
  AliForwardMultiplicityBase::operator=(o);

  fHData             = o.fHData;
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
AliForwardMultiplicityTask::SetDebug(Int_t dbg)
{
  // 
  // Set debug level 
  // 
  // Parameters:
  //    dbg Debug level
  //
  fEventInspector.SetDebug(dbg);
  fEnergyFitter.SetDebug(dbg);
  fSharingFilter.SetDebug(dbg);
  fDensityCalculator.SetDebug(dbg);
  fCorrections.SetDebug(dbg);
  fHistCollector.SetDebug(dbg);
}

//____________________________________________________________________
void
AliForwardMultiplicityTask::InitializeSubs()
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
  const TAxis* pe = 0;
  const TAxis* pv = 0;

  if (!ReadCorrections(pe,pv)) return;

  fHistos.Init(*pe);
  fAODFMD.Init(*pe);

  fHData = static_cast<TH2D*>(fAODFMD.GetHistogram().Clone("d2Ndetadphi"));
  fHData->SetStats(0);
  fHData->SetDirectory(0);
  fList->Add(fHData);

  fEnergyFitter.Init(*pe);
  fEventInspector.Init(*pv);
  fDensityCalculator.Init(*pe);
  fCorrections.Init(*pe);
  fHistCollector.Init(*pv);

  this->Print();
}

//____________________________________________________________________
void
AliForwardMultiplicityTask::UserCreateOutputObjects()
{
  // 
  // Create output objects 
  // 
  //
  fList = new TList;

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) AliFatal("No AOD output handler set in analysis manager");
    
    
  TObject* obj = &fAODFMD;
  ah->AddBranch("AliAODForwardMult", &obj);

  fEventInspector.DefineOutput(fList);
  fEnergyFitter.DefineOutput(fList);
  fSharingFilter.DefineOutput(fList);
  fDensityCalculator.DefineOutput(fList);
  fCorrections.DefineOutput(fList);

  PostData(1, fList);
}
//____________________________________________________________________
void
AliForwardMultiplicityTask::UserExec(Option_t*)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  

  // static Int_t cnt = 0;
  // cnt++;
  // Get the input data 
  AliESDEvent* esd = GetESDEvent();

  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  fAODFMD.Clear();
  
  Bool_t   lowFlux  = kFALSE;
  UInt_t   triggers = 0;
  UShort_t ivz      = 0;
  Double_t vz       = 0;
  Double_t cent     = -1;
  UInt_t   found    = fEventInspector.Process(esd, triggers, lowFlux, 
					      ivz, vz, cent);
  
  if (found & AliFMDEventInspector::kNoEvent)    return;
  if (found & AliFMDEventInspector::kNoTriggers) return;

  // Set trigger bits, and mark this event for storage 
  fAODFMD.SetTriggerBits(triggers);
  fAODFMD.SetSNN(fEventInspector.GetEnergy());
  fAODFMD.SetSystem(fEventInspector.GetCollisionSystem());
  fAODFMD.SetCentrality(cent);
  MarkEventForStore();
  
  if (found & AliFMDEventInspector::kNoSPD)      return;
  if (found & AliFMDEventInspector::kNoFMD)      return;
  if (found & AliFMDEventInspector::kNoVertex)   return;
  
  if (triggers & AliAODForwardMult::kPileUp) return;
  
  fAODFMD.SetIpZ(vz);

  if (found & AliFMDEventInspector::kBadVertex) return;

  // We we do not want to use low flux specific code, we disable it here. 
  if (!fEnableLowFlux) lowFlux = false;

  // Get FMD data 
  AliESDFMD* esdFMD = esd->GetFMDData();
  // Apply the sharing filter (or hit merging or clustering if you like)
  if (!fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD)) { 
    AliWarning("Sharing filter failed!");
    return;
  }

  // Do the energy stuff 
  if (!fEnergyFitter.Accumulate(*esdFMD, cent, 
				triggers & AliAODForwardMult::kEmpty)){
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
AliForwardMultiplicityTask::Terminate(Option_t*)
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //

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
  TH1D* norm   = hData->ProjectionX("norm",   0,  1,  "");
  dNdeta->SetTitle("dN_{ch}/d#eta in the forward regions");
  dNdeta->SetYTitle("#frac{1}{N}#frac{dN_{ch}}{d#eta}");
  dNdeta->Divide(norm);
  dNdeta->SetStats(0);
  dNdeta->Scale(Double_t(hEventsTrVtx->GetEntries())/hEventsTr->GetEntries(),
		"width");
  list->Add(dNdeta);
  list->Add(norm);

  fEnergyFitter.Fit(list);
  fSharingFilter.ScaleHistograms(list,Int_t(hEventsTr->Integral()));
  fDensityCalculator.ScaleHistograms(list,Int_t(hEventsTrVtx->Integral()));
  fCorrections.ScaleHistograms(list,Int_t(hEventsTrVtx->Integral()));
}

//
// EOF
//
