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
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fHTriggers(0),
    fHData(0),
    fFirstEvent(true),
    fLowFluxCut(1000),
    fESDFMD(),
    fHistos(),
    fAODFMD(),
    fSharingFilter(),
    fDensityCalculator(),
    fCorrections(),
    fHistCollector(),
    fList(0), 
    fTree(0)
{
}

//____________________________________________________________________
AliForwardMultiplicity::AliForwardMultiplicity(const char* name)
  : AliAnalysisTaskSE(name), 
    fHEventsTr(0), 
    fHEventsTrVtx(0), 
    fHTriggers(0),
    fHData(0),
    fFirstEvent(true),
    fLowFluxCut(1000),
    fESDFMD(),
    fHistos(),
    fAODFMD(kTRUE),
    fSharingFilter("sharing"), 
    fDensityCalculator("density"),
    fCorrections("corrections"),
    fHistCollector("collector"),
    fList(0), 
    fTree(0)
{
  DefineOutput(1, TList::Class());
  // DefineOutput(2, TTree::Class());
}

//____________________________________________________________________
AliForwardMultiplicity::AliForwardMultiplicity(const AliForwardMultiplicity& o)
  : AliAnalysisTaskSE(o),
    fHEventsTr(o.fHEventsTr), 
    fHEventsTrVtx(o.fHEventsTrVtx), 
    fHTriggers(o.fHTriggers),
    fHData(o.fHData),
    fFirstEvent(true),
    fLowFluxCut(1000),
    fESDFMD(o.fESDFMD),
    fHistos(o.fHistos),
    fAODFMD(o.fAODFMD),
    fSharingFilter(o.fSharingFilter),
    fDensityCalculator(o.fDensityCalculator),
    fCorrections(o.fCorrections),
    fHistCollector(o.fHistCollector),
    fList(o.fList), 
    fTree(o.fTree)
{
}

//____________________________________________________________________
AliForwardMultiplicity&
AliForwardMultiplicity::operator=(const AliForwardMultiplicity& o)
{
  fHEventsTr         = o.fHEventsTr;
  fHEventsTrVtx      = o.fHEventsTrVtx;
  fHTriggers         = o.fHTriggers;
  fHData             = o.fHData;
  fFirstEvent        = o.fFirstEvent;
  fSharingFilter     = o.fSharingFilter;
  fDensityCalculator = o.fDensityCalculator;
  fCorrections       = o.fCorrections;
  fHistCollector     = o.fHistCollector;
  fHistos            = o.fHistos;
  fAODFMD            = o.fAODFMD;
  fList              = o.fList;
  fTree              = o.fTree;

  return *this;
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

  fHEventsTr = new TH1I("nEventsTr", "Number of events w/trigger", 
		      pars->GetNvtxBins(), 
		      -pars->GetVtxCutZ(), 
		      pars->GetVtxCutZ());
  fHEventsTr->SetXTitle("v_{z} [cm]");
  fHEventsTr->SetYTitle("# of events");
  fHEventsTr->SetFillColor(kRed+1);
  fHEventsTr->SetFillStyle(3001);
  fHEventsTr->SetDirectory(0);
  fList->Add(fHEventsTr);
  // fHEventsTr->Sumw2();

  fHEventsTrVtx = new TH1I("nEventsTrVtx", 
			   "Number of events w/trigger and vertex", 
			   pars->GetNvtxBins(), 
			   -pars->GetVtxCutZ(), 
			   pars->GetVtxCutZ());
  fHEventsTrVtx->SetXTitle("v_{z} [cm]");
  fHEventsTrVtx->SetYTitle("# of events");
  fHEventsTrVtx->SetFillColor(kBlue+1);
  fHEventsTrVtx->SetFillStyle(3001);
  fHEventsTrVtx->SetDirectory(0);
  fList->Add(fHEventsTrVtx);
  // fHEventsTrVtx->Sumw2();

      
  fHTriggers = new TH1I("triggers", "Triggers", 10, 0, 10);
  fHTriggers->SetFillColor(kRed+1);
  fHTriggers->SetFillStyle(3001);
  fHTriggers->SetStats(0);
  fHTriggers->SetDirectory(0);
  fHTriggers->GetXaxis()->SetBinLabel(1,"INEL");
  fHTriggers->GetXaxis()->SetBinLabel(2,"INEL>0");
  fHTriggers->GetXaxis()->SetBinLabel(3,"NSD");
  fHTriggers->GetXaxis()->SetBinLabel(4,"Empty");
  fHTriggers->GetXaxis()->SetBinLabel(5,"A");
  fHTriggers->GetXaxis()->SetBinLabel(6,"B");
  fHTriggers->GetXaxis()->SetBinLabel(7,"C");
  fHTriggers->GetXaxis()->SetBinLabel(8,"E");
  fList->Add(fHTriggers);

  TAxis e(pars->GetNetaBins(), pars->GetEtaMin(), pars->GetEtaMax());
  fHistos.Init(e);
  fAODFMD.Init(e);

  fHData = static_cast<TH2D*>(fAODFMD.GetHistogram().Clone("d2Ndetadphi"));
  fHData->SetStats(0);
  fHData->SetDirectory(0);
  fList->Add(fHData);

  fHistCollector.Init(*(fHEventsTr->GetXaxis()));
}

//____________________________________________________________________
void
AliForwardMultiplicity::UserCreateOutputObjects()
{
  fList = new TList;

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah)  
    AliFatal("No AOD output handler set in analysis manager");
    
    
  TObject* obj = &fAODFMD;
  ah->AddBranch("AliAODForwardMult", &obj);

  fSharingFilter.DefineOutput(fList);
  fDensityCalculator.DefineOutput(fList);
  fCorrections.DefineOutput(fList);
}
//____________________________________________________________________
void
AliForwardMultiplicity::UserExec(Option_t*)
{
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();

  // Get the input data 
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) { 
    AliWarning("No ESD event found for input event");
    return;
  }

#if 0
  static Int_t nEvents = 0;
  nEvents++;
  if (nEvents % 100 == 0) AliInfo(Form("Event # %6d", nEvents));
#endif

  // On the first event, initialize the parameters 
  if (fFirstEvent) { 
    AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
    pars->SetParametersFromESD(esd);
    pars->PrintStatus();
    fFirstEvent = false;

    InitializeSubs();
  }
  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  fAODFMD.Clear();

  // Read trigger information from the ESD and store in AOD object
  UInt_t triggers = 0;
  if (!AliForwardUtil::ReadTriggers(esd, triggers, fHTriggers)) { 
    if (am->GetDebugLevel() > 1) 
      AliWarning("Failed to read triggers from ESD");
    return;
  }
  fAODFMD.SetTriggerBits(triggers);

  // Mark this event for storage 
  MarkEventForStore();

  // Check if this is a high-flux event 
  const AliMultiplicity* testmult = esd->GetMultiplicity();
  if (!testmult) {
    if (am->GetDebugLevel() > 1) 
      AliWarning("No central multiplicity object found");
    return;
  }
  Bool_t lowFlux = testmult->GetNumberOfTracklets() < fLowFluxCut;
  if (am->GetDebugLevel() > 1) 
    AliInfo(Form("Event has %d SPD tracklets, cut is %d, this is a %s event",
		 testmult, fLowFluxCut, (lowFlux ? "low" : "high")));

  // Get the FMD ESD data 
  AliESDFMD* esdFMD = esd->GetFMDData();
  if (!esdFMD) { 
    if (am->GetDebugLevel() > 1) 
      AliWarning("No FMD data found in ESD");
    return;
  }

  // Get the vertex information 
  Double_t vz   = 0;
  Bool_t   vzOk = AliForwardUtil::ReadVertex(esd, vz);

  fHEventsTr->Fill(vz);
  if (!vzOk) { 
    if (am->GetDebugLevel() > 1) 
      AliWarning("Failed to read vertex from ESD");
    return;
  }
  fHEventsTrVtx->Fill(vz);

  // Get the vertex bin 
  Int_t ivz = fHEventsTr->GetXaxis()->FindBin(vz)-1;
  fAODFMD.SetIpZ(vz);
  if (ivz < 0 || ivz >= fHEventsTr->GetXaxis()->GetNbins()) { 
    if (am->GetDebugLevel() > 1) 
      AliWarning(Form("Vertex @ %f outside of range [%f,%f]", 
		      vz, fHEventsTr->GetXaxis()->GetXmin(), 
		      fHEventsTr->GetXaxis()->GetXmax()));
    return;
  }
  if (am->GetDebugLevel() > 2) 
    AliInfo(Form("Events vertex @ %f (bin %d), in range", vz, ivz));
  

  // Apply the sharing filter (or hit merging or clustering if you like)
  if (!fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD, vz)) { 
    AliWarning("Sharing filter failed!");
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
  TH1I* hEventsTr    = static_cast<TH1I*>(list->FindObject("nEventsTr"));
  TH1I* hEventsTrVtx = static_cast<TH1I*>(list->FindObject("nEventsTrVtx"));
  TH2D* hData        = static_cast<TH2D*>(list->FindObject("d2Ndetadphi"));
  if (!hData || !hEventsTr || !hEventsTrVtx) { 
    AliError(Form("one or more histograms could not be found in output "
		  "list %s (hEventsTr=%p,hEventsTrVtx=%p,d2Ndetadphi=%p)", 
		  list->GetName(), hEventsTr, hEventsTrVtx, hData));
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
