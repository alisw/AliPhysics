#include "AliFMDEventInspector.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliAODForwardMult.h"
#include <TH1.h>
#include <TList.h>
#include <TDirectory.h>

//====================================================================
AliFMDEventInspector::AliFMDEventInspector()
  : TNamed(),
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fHTriggers(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.1),
    fList(0),
    fDebug(0)
{
}

//____________________________________________________________________
AliFMDEventInspector::AliFMDEventInspector(const char* name)
  : TNamed("fmdEventInspector", name),
    fHEventsTr(0), 
    fHEventsTrVtx(0), 
    fHTriggers(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.1),
    fList(0),
    fDebug(0)
{
}

//____________________________________________________________________
AliFMDEventInspector::AliFMDEventInspector(const AliFMDEventInspector& o)
  : TNamed(o), 
    fHEventsTr(o.fHEventsTr), 
    fHEventsTrVtx(o.fHEventsTrVtx), 
    fHTriggers(o.fHTriggers),
    fLowFluxCut(o.fMaxVzErr),
    fMaxVzErr(o.fMaxVzErr),
    fList(o.fList),
    fDebug(0)
{
}

//____________________________________________________________________
AliFMDEventInspector::~AliFMDEventInspector()
{
  if (fHEventsTr)    delete fHEventsTr;
  if (fHEventsTrVtx) delete fHEventsTrVtx;
  if (fHTriggers)    delete fHTriggers;  
  if (fList)         delete fList;
}
//____________________________________________________________________
AliFMDEventInspector&
AliFMDEventInspector::operator=(const AliFMDEventInspector& o)
{
  TNamed::operator=(o);
  fHEventsTr         = o.fHEventsTr;
  fHEventsTrVtx      = o.fHEventsTrVtx;
  fHTriggers         = o.fHTriggers;
  fLowFluxCut        = o.fLowFluxCut;
  fMaxVzErr          = o.fMaxVzErr;
  fDebug             = o.fDebug;
  fList              = (o.fList ? new TList : 0);
  if (fList) { 
    fList->SetName(GetName());
    if (fHEventsTr)    fList->Add(fHEventsTr);
    if (fHEventsTrVtx) fList->Add(fHEventsTrVtx);
    if (fHTriggers)    fList->Add(fHTriggers);
  }
  return *this;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::FetchHistograms(TList* d, 
				      TH1I*& hEventsTr, 
				      TH1I*& hEventsTrVtx, 
				      TH1I*& hTriggers) const
{
  hEventsTr    = 0;
  hEventsTrVtx = 0;
  hTriggers    = 0;
  TList* dd    = dynamic_cast<TList*>(d->FindObject(GetName()));
  if (!dd) return kFALSE;
  
  hEventsTr    = dynamic_cast<TH1I*>(dd->FindObject("nEventsTr"));
  hEventsTrVtx = dynamic_cast<TH1I*>(dd->FindObject("nEventsTrVtx"));
  hTriggers    = dynamic_cast<TH1I*>(dd->FindObject("triggers"));

  if (!hEventsTr || !hEventsTrVtx || !hTriggers) return kFALSE;
  return kTRUE;
}
//____________________________________________________________________
void
AliFMDEventInspector::Init(const TAxis& vtxAxis)
{
  fHEventsTr = new TH1I("nEventsTr", "Number of events w/trigger", 
			vtxAxis.GetNbins(), 
			vtxAxis.GetXmin(), 
			vtxAxis.GetXmax());
  fHEventsTr->SetXTitle("v_{z} [cm]");
  fHEventsTr->SetYTitle("# of events");
  fHEventsTr->SetFillColor(kRed+1);
  fHEventsTr->SetFillStyle(3001);
  fHEventsTr->SetDirectory(0);
  // fHEventsTr->Sumw2();
  fList->Add(fHEventsTr);

  fHEventsTrVtx = new TH1I("nEventsTrVtx", 
			   "Number of events w/trigger and vertex", 
			   vtxAxis.GetNbins(), 
			   vtxAxis.GetXmin(), 
			   vtxAxis.GetXmax());
  fHEventsTrVtx->SetXTitle("v_{z} [cm]");
  fHEventsTrVtx->SetYTitle("# of events");
  fHEventsTrVtx->SetFillColor(kBlue+1);
  fHEventsTrVtx->SetFillStyle(3001);
  fHEventsTrVtx->SetDirectory(0);
  // fHEventsTrVtx->Sumw2();
  fList->Add(fHEventsTrVtx);

      
  fHTriggers = new TH1I("triggers", "Triggers", 10, 0, 10);
  fHTriggers->SetFillColor(kRed+1);
  fHTriggers->SetFillStyle(3001);
  fHTriggers->SetStats(0);
  fHTriggers->SetDirectory(0);
  fHTriggers->GetXaxis()->SetBinLabel(kInel   +1,"INEL");
  fHTriggers->GetXaxis()->SetBinLabel(kInelGt0+1,"INEL>0");
  fHTriggers->GetXaxis()->SetBinLabel(kNSD    +1,"NSD");
  fHTriggers->GetXaxis()->SetBinLabel(kEmpty  +1,"Empty");
  fHTriggers->GetXaxis()->SetBinLabel(kA      +1,"A");
  fHTriggers->GetXaxis()->SetBinLabel(kB      +1,"B");
  fHTriggers->GetXaxis()->SetBinLabel(kC      +1,"C");
  fHTriggers->GetXaxis()->SetBinLabel(kE      +1,"E");
  fHTriggers->GetXaxis()->SetBinLabel(9,         "spare1");
  fHTriggers->GetXaxis()->SetBinLabel(10,        "spare2");
  fList->Add(fHTriggers);
}

//____________________________________________________________________
void
AliFMDEventInspector::DefineOutput(TList* dir)
{
  fList = new TList;
  fList->SetName(GetName());
  dir->Add(fList);
}

//____________________________________________________________________
UInt_t
AliFMDEventInspector::Process(const AliESDEvent* event, 
			      UInt_t&            triggers,
			      Bool_t&            lowFlux,
			      Int_t&             ivz, 
			      Double_t&          vz)
{
  // Check that we have an event 
  if (!event) { 
    AliWarning("No ESD event found for input event");
    return kNoEvent;
  }

  // Read trigger information from the ESD and store in AOD object
  if (!ReadTriggers(event, triggers)) { 
    if (fDebug > 2) 
      AliWarning("Failed to read triggers from ESD");
    return kNoTriggers;
  }

  // Check if this is a high-flux event 
  const AliMultiplicity* testmult = event->GetMultiplicity();
  if (!testmult) {
    if (fDebug > 3) 
      AliWarning("No central multiplicity object found");
    return kNoSPD;
  }
  lowFlux = testmult->GetNumberOfTracklets() < fLowFluxCut;

  // Check the FMD ESD data 
  if (!event->GetFMDData()) { 
    if (fDebug > 3) 
      AliWarning("No FMD data found in ESD");
    return kNoFMD;
  }

  // Get the vertex information 
  vz          = 0;
  Bool_t vzOk = ReadVertex(event, vz);

  fHEventsTr->Fill(vz);
  if (!vzOk) { 
    if (fDebug > 3) 
      AliWarning("Failed to read vertex from ESD");
    return kNoVertex;
  }
  fHEventsTrVtx->Fill(vz);

  // Get the vertex bin 
  ivz = fHEventsTr->GetXaxis()->FindBin(vz)-1;
  if (ivz < 0 || ivz >= fHEventsTr->GetXaxis()->GetNbins()) { 
    if (fDebug > 3) 
      AliWarning(Form("Vertex @ %f outside of range [%f,%f]", 
		      vz, fHEventsTr->GetXaxis()->GetXmin(), 
		      fHEventsTr->GetXaxis()->GetXmax()));
    ivz = -1;
    return kBadVertex;
  }
  return kOk;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadTriggers(const AliESDEvent* esd, UInt_t& triggers)
{
  triggers = 0;

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
    triggers |= AliAODForwardMult::kInel;
    fHTriggers->Fill(kInel+0.5);
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
	  triggers |= AliAODForwardMult::kInelGt0;
	  fHTriggers->Fill(kInelGt0+.5);
	  break;
	}
      }
    }
  }

  // Analyse some trigger stuff 
  AliTriggerAnalysis ta;
  if (ta.IsOfflineTriggerFired(esd, AliTriggerAnalysis::kNSD1)) {
    triggers |= AliAODForwardMult::kNSD;
    fHTriggers->Fill(kNSD+.5);
  }

  // Get trigger stuff 
  TString trigStr = esd->GetFiredTriggerClasses();
  if (trigStr.Contains("CBEAMB-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kEmpty;
    fHTriggers->Fill(kEmpty+.5);
  }

  if (trigStr.Contains("CINT1A-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kA;
    fHTriggers->Fill(kA+.5);
  }

  if (trigStr.Contains("CINT1B-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kB;
    fHTriggers->Fill(kB+.5);
  }


  if (trigStr.Contains("CINT1C-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kC;
    fHTriggers->Fill(kC+.5);
  }

  if (trigStr.Contains("CINT1-E-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kE;
    fHTriggers->Fill(kE+.5);
  }

  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadVertex(const AliESDEvent* esd, Double_t& vz)
{
  vz = 0;
  // Get the vertex 
  const AliESDVertex* vertex = esd->GetPrimaryVertexSPD();
  if (!vertex) { 
    if (fDebug > 2) 
      AliWarning("No SPD vertex found in ESD");
    return kFALSE;
  }
  
  // Check that enough tracklets contributed 
  if(vertex->GetNContributors() <= 0) {
    if (fDebug > 2)
      AliWarning(Form("Number of contributors to vertex is %d<=0",
		      vertex->GetNContributors()));
    vz = 0;
    return kFALSE;
  }

  // Check that the uncertainty isn't too large 
  if (vertex->GetZRes() > fMaxVzErr) { 
    if (fDebug > 2)
      AliWarning(Form("Uncertaintity in Z of vertex is too large %f > %f", 
		      vertex->GetZRes(), fMaxVzErr));
    return kFALSE;
  }

  // Get the z coordiante 
  vz = vertex->GetZ();
  return kTRUE;
}
//
// EOF
//

