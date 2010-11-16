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

  fHEventsTr = new TH1I("nEvents", "Number of events w/trigger", 
		      pars->GetNvtxBins(), 
		      -pars->GetVtxCutZ(), 
		      pars->GetVtxCutZ());
  fHEventsTr->SetXTitle("v_{z} [cm]");
  fHEventsTr->SetYTitle("# of events");
  fHEventsTr->SetFillColor(kRed+1);
  fHEventsTr->SetFillStyle(3001);
  fHEventsTr->SetDirectory(0);
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

  TAxis e(pars->GetNetaBins(), pars->GetEtaMin(), pars->GetEtaMax());
  fHistos.Init(e);
  fAODFMD.Init(e);

  fHData = static_cast<TH2D*>(fAODFMD.GetHistogram().Clone("d2Ndetadphi"));
  fHData->SetStats(0);
  fHData->SetDirectory(0);
  fSharingFilter.Init();
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

  // fTree = new TTree("T", "T");
  // fTree->Branch("forward", &fAODFMD);

  PostData(1, fList);
  // PostData(2, fTree);
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
  if (!ReadTriggers(esd)) { 
#ifdef VERBOSE
    AliWarning("Failed to read triggers from ESD");
#endif
    return;
  }

  // Mark this event for storage 
  MarkEventForStore();

  // Check if this is a high-flux event 
  const AliMultiplicity* testmult = esd->GetMultiplicity();
  if (!testmult) { 
#ifdef VERBOSE
    AliWarning("No central multiplicity object found");
#endif
    return;
  }
  Bool_t lowFlux = testmult->GetNumberOfTracklets() < fLowFluxCut;

  // Get the FMD ESD data 
  AliESDFMD* esdFMD = esd->GetFMDData();
  if (!esdFMD) { 
#ifdef VERBOSE
    AliWarning("No FMD data found in ESD");
#endif
    return;
  }

  // Get the vertex information 
  Double_t vz   = 0;
  Bool_t   vzOk = ReadVertex(esd, vz);

  fHEventsTr->Fill(vz);
  if (!vzOk) { 
#ifdef VERBOSE
    AliWarning("Failed to read vertex from ESD");
#endif
    return;
  }
  fHEventsTrVtx->Fill(vz);

  // Get the vertex bin 
  Int_t ivz = fHEventsTr->GetXaxis()->FindBin(vz)-1;
  fAODFMD.SetIpZ(vz);
  if (ivz < 0 || ivz >= fHEventsTr->GetXaxis()->GetNbins()) { 
#if 0
    AliWarning(Form("Vertex @ %f outside of range [%f,%f]", 
		    vz, fHEventsTr->GetXaxis()->GetXmin(), 
		    fHEventsTr->GetXaxis()->GetXmax()));
#endif
    return;
  }

  // Apply the sharing filter (or hit merging or clustering if you like)
  if (!fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD, vz)) { 
#ifdef VERBOSE
    AliWarning("Sharing filter failed!");
#endif
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
}

//____________________________________________________________________
void
AliForwardMultiplicity::Terminate(Option_t*)
{
  TList* list = dynamic_cast<TList*>(GetOutputData(1));
  if (!list) {
    AliError("No output list defined");
    return;
  }
  // TH1D* dNdeta = fHData->ProjectionX("dNdeta", 0, -1, "e");
  TH1D* dNdeta = fHData->ProjectionX("dNdeta", 1, -1, "e");
  TH1D* norm   = fHData->ProjectionX("dNdeta", 0, 1,  "");
  dNdeta->SetTitle("dN_{ch}/d#eta in the forward regions");
  dNdeta->SetYTitle("#frac{1}{N}#frac{dN_{ch}}{d#eta}");
  dNdeta->Divide(norm);
  dNdeta->SetStats(0);
  dNdeta->Scale(Double_t(fHEventsTrVtx->GetEntries())/fHEventsTr->GetEntries(),
		"width");

  list->Add(fHEventsTr);
  list->Add(fHEventsTrVtx);
  list->Add(fHTriggers);
  list->Add(fHData);
  list->Add(dNdeta);
  
  TList* last = new TList;
  last->SetName("LastEvent");
  list->Add(last);
  last->Add(&fAODFMD.GetHistogram());
  last->Add(fHistos.fFMD1i);
  last->Add(fHistos.fFMD2i);
  last->Add(fHistos.fFMD2o);
  last->Add(fHistos.fFMD3i);
  last->Add(fHistos.fFMD3o);


  fSharingFilter.ScaleHistograms(fHEventsTr->Integral());
  fSharingFilter.Output(list);

  fDensityCalculator.ScaleHistograms(fHEventsTrVtx->Integral());
  fDensityCalculator.Output(list);

  fCorrections.ScaleHistograms(fHEventsTrVtx->Integral());
  fCorrections.Output(list);
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
Bool_t
AliForwardMultiplicity::ReadTriggers(AliESDEvent* esd)
{
  // Get the analysis manager - should always be there 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  if (!am) { 
    AliWarning("No analysis manager defined!");
    return kFALSE;
  }

  // Get the input handler - should always be there 
  AliInputEventHandler* ih = 
    static_cast<AliInputEventHandler*>(am->GetInputEventHandler());
  if (!ih) { 
    AliWarning("No input handler");
    return kFALSE;
  }
  
  // Get the physics selection - add that by using the macro 
  // AddTaskPhysicsSelection.C 
  AliPhysicsSelection* ps = 
    static_cast<AliPhysicsSelection*>(ih->GetEventSelection());
  if (!ps) { 
    AliWarning("No physics selection");
    return kFALSE;
  }
  
  // Check if this is a collision candidate (INEL)
  Bool_t inel = ps->IsCollisionCandidate(esd);
  if (inel) { 
    fAODFMD.SetTriggerBits(AliAODForwardMult::kInel);
    fHTriggers->Fill(.5);
  }
  

  // IF this is inel, see if we have a tracklet 
  if (inel) { 
    const AliMultiplicity* spdmult = esd->GetMultiplicity();
    if (!spdmult) {
      AliWarning("No SPD multiplicity");
    }
    else { 
      Int_t n = spdmult->GetNumberOfTracklets();
      for (Int_t j = 0; j < n; j++) { 
	if(TMath::Abs(spdmult->GetEta(j)) < 1) { 
	  fAODFMD.SetTriggerBits(AliAODForwardMult::kInelGt0);
	  fHTriggers->Fill(1.5);
	  break;
	}
      }
    }
  }

  // Analyse some trigger stuff 
  AliTriggerAnalysis ta;
  if (ta.IsOfflineTriggerFired(esd, AliTriggerAnalysis::kNSD1)) {
    fAODFMD.SetTriggerBits(AliAODForwardMult::kNSD);
    fHTriggers->Fill(2.5);
  }

  // Get trigger stuff 
  TString triggers = esd->GetFiredTriggerClasses();
  if (triggers.Contains("CBEAMB-ABCE-NOPF-ALL")) {
    fAODFMD.SetTriggerBits(AliAODForwardMult::kEmpty);
    fHTriggers->Fill(3.5);
  }

  if (triggers.Contains("CINT1A-ABCE-NOPF-ALL")) {
    fAODFMD.SetTriggerBits(AliAODForwardMult::kA);
    fHTriggers->Fill(4.5);
  }

  if (triggers.Contains("CINT1B-ABCE-NOPF-ALL")) {
    fAODFMD.SetTriggerBits(AliAODForwardMult::kB);
    fHTriggers->Fill(5.5);
  }


  if (triggers.Contains("CINT1C-ABCE-NOPF-ALL")) {
    fAODFMD.SetTriggerBits(AliAODForwardMult::kC);
    fHTriggers->Fill(6.5);
  }

  if (triggers.Contains("CINT1-E-NOPF-ALL")) {
    fAODFMD.SetTriggerBits(AliAODForwardMult::kE);
    fHTriggers->Fill(7.5);
  }

  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliForwardMultiplicity::ReadVertex(AliESDEvent* esd, Double_t& vz)
{
  // Get the vertex 
  const AliESDVertex* vertex = esd->GetPrimaryVertexSPD();
  if (!vertex) { 
#ifdef VERBOSE
    AliWarning("No SPD vertex found in ESD");
#endif
    return kFALSE;
  }

  // Check that enough tracklets contributed 
  if(vertex->GetNContributors() <= 0) {
#ifdef VERBOSE
    AliWarning(Form("Number of contributors to vertex is %d<=0",
		    vertex->GetNContributors()));
#endif
    return kFALSE;
  }

  // Check that the uncertainty isn't too large 
  if (vertex->GetZRes() > 0.1) { 
#ifdef VERBOSE
    AliWarning(Form("Uncertaintity in Z of vertex is too large %f > 0.1", 
		    vertex->GetZRes()));
#endif
    return kFALSE;
  }

  // Get the z coordiante 
  vz = vertex->GetZ();
	       
  return kTRUE;
}

//____________________________________________________________________
void
AliForwardMultiplicity::Print(Option_t*) const
{
}

//
// EOF
//
