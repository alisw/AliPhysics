// 
// Calculate the corrections in the central regions
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODCentralMult 
// 
// Histograms 
//   
// Corrections used 
// 
#include "AliCentralMCCorrectionsTask.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliESDEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODForwardMult.h"
#include "AliCentralCorrSecondaryMap.h"
#include <TH1.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TList.h>
#include <TROOT.h>
#include <iostream>

//====================================================================
namespace {
  const char* GetEventName(Bool_t tr, Bool_t vtx) 
  {
    return Form("nEvents%s%s", (tr ? "Tr" : ""), (vtx ? "Vtx" : ""));
  }
  const char* GetHitsName(UShort_t d, Char_t r) 
  {
    return Form("hitsSPD%d%c", d, r);
  }
  const char* GetStripsName(UShort_t d, Char_t r)
  {
    return Form("stripsSPD%d%c", d, r);
  }
  const char* GetPrimaryName(Char_t r, Bool_t trVtx)
  {
    return Form("primaries%s%s", 
		((r == 'I' || r == 'i') ? "Inner" : "Outer"), 
		(trVtx ? "TrVtx" : "All"));
  }
}

//====================================================================
AliCentralMCCorrectionsTask::AliCentralMCCorrectionsTask()
  : AliAnalysisTaskSE(),
    fInspector(),
    fTrackDensity(),
    fVtxBins(0),
    fFirstEvent(true),
    fHEvents(0), 
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fVtxAxis(),
    fEtaAxis(),
    fList(),
    fNPhiBins(20)    
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
}

//____________________________________________________________________
AliCentralMCCorrectionsTask::AliCentralMCCorrectionsTask(const char* name)
  : AliAnalysisTaskSE(name),
    fInspector("eventInspector"), 
    fTrackDensity("trackDensity"),
    fVtxBins(0),
    fFirstEvent(true),
    fHEvents(0), 
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fVtxAxis(10,-10,10), 
    fEtaAxis(200,-4,6),
    fList(),
    fNPhiBins(20)    
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________________
AliCentralMCCorrectionsTask::AliCentralMCCorrectionsTask(const AliCentralMCCorrectionsTask& o)
  : AliAnalysisTaskSE(o),
    fInspector(o.fInspector),
    fTrackDensity(),
    fVtxBins(0),
    fFirstEvent(o.fFirstEvent),
    fHEvents(o.fHEvents), 
    fHEventsTr(o.fHEventsTr), 
    fHEventsTrVtx(o.fHEventsTrVtx),
    fVtxAxis(10,-10,10), 
    fEtaAxis(200,-4,6),
    fList(o.fList),
    fNPhiBins(o.fNPhiBins)    
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
}

//____________________________________________________________________
AliCentralMCCorrectionsTask&
AliCentralMCCorrectionsTask::operator=(const AliCentralMCCorrectionsTask& o)
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
  if (&o == this) return *this; 
  fInspector         = o.fInspector;
  fTrackDensity      = o.fTrackDensity;
  fVtxBins           = o.fVtxBins;
  fFirstEvent        = o.fFirstEvent;
  fHEvents           = o.fHEvents;
  fHEventsTr         = o.fHEventsTr;
  fHEventsTrVtx      = o.fHEventsTrVtx;
  SetVertexAxis(o.fVtxAxis);
  SetEtaAxis(o.fEtaAxis);
  fNPhiBins          = o.fNPhiBins;

  return *this;
}

//____________________________________________________________________
void
AliCentralMCCorrectionsTask::Init()
{
  // 
  // Initialize the task 
  // 
  //
}

//____________________________________________________________________
void
AliCentralMCCorrectionsTask::SetVertexAxis(Int_t nBin, Double_t min, 
					   Double_t max)
{
  // 
  // Set the vertex axis to use
  // 
  // Parameters:
  //    nBins Number of bins
  //    vzMin Least @f$z@f$ coordinate of interation point
  //    vzMax Largest @f$z@f$ coordinate of interation point
  //
  if (max < min) max = -min;
  if (min < max) { 
    Double_t tmp = min;
    min          = max;
    max          = tmp;
  }
  if (min < -10) 
    AliWarning(Form("Minimum vertex %f < -10, make sure you want this",min));
  if (max > +10) 
    AliWarning(Form("Minimum vertex %f > +10, make sure you want this",max));
  fVtxAxis.Set(nBin, min, max);
}
//____________________________________________________________________
void
AliCentralMCCorrectionsTask::SetVertexAxis(const TAxis& axis)
{
  // 
  // Set the vertex axis to use
  // 
  // Parameters:
  //    axis Axis
  //
  SetVertexAxis(axis.GetNbins(),axis.GetXmin(),axis.GetXmax());
}

//____________________________________________________________________
void
AliCentralMCCorrectionsTask::SetEtaAxis(Int_t nBin, Double_t min, Double_t max)
{
  // 
  // Set the eta axis to use
  // 
  // Parameters:
  //    nBins Number of bins
  //    vzMin Least @f$\eta@f$ 
  //    vzMax Largest @f$\eta@f$ 
  //
  if (max < min) max = -min;
  if (min < max) { 
    Double_t tmp = min;
    min          = max;
    max          = tmp;
  }
  if (min < -4) 
    AliWarning(Form("Minimum eta %f < -4, make sure you want this",min));
  if (max > +6) 
    AliWarning(Form("Minimum eta %f > +6, make sure you want this",max));
  fVtxAxis.Set(nBin, min, max);
}
//____________________________________________________________________
void
AliCentralMCCorrectionsTask::SetEtaAxis(const TAxis& axis)
{
  // 
  // Set the eta axis to use
  // 
  // Parameters:
  //    axis Axis
  //
  SetEtaAxis(axis.GetNbins(),axis.GetXmin(),axis.GetXmax());
}

//____________________________________________________________________
void
AliCentralMCCorrectionsTask::DefineBins(TList* l)
{
  if (!fVtxBins) fVtxBins = new TObjArray(fVtxAxis.GetNbins(), 1);
  if (fVtxBins->GetEntries() > 0) return;

  fVtxBins->SetOwner();
  for (Int_t i = 1; i <= fVtxAxis.GetNbins(); i++) { 
    Double_t low  = fVtxAxis.GetBinLowEdge(i);
    Double_t high = fVtxAxis.GetBinUpEdge(i);
    VtxBin*  bin  = new VtxBin(low, high, fEtaAxis, fNPhiBins);
    fVtxBins->AddAt(bin, i);
    bin->DefineOutput(l);
  }
}

//____________________________________________________________________
void
AliCentralMCCorrectionsTask::UserCreateOutputObjects()
{
  // 
  // Create output objects 
  // 
  //
  fList = new TList;
  fList->SetOwner();
  fList->SetName(Form("%sSums", GetName()));

  DefineBins(fList);

  fHEvents = new TH1I(GetEventName(false,false),
		      "Number of all events", 
		      fVtxAxis.GetNbins(), 
		      fVtxAxis.GetXmin(), 
		      fVtxAxis.GetXmax());
  fHEvents->SetXTitle("v_{z} [cm]");
  fHEvents->SetYTitle("# of events");
  fHEvents->SetFillColor(kBlue+1);
  fHEvents->SetFillStyle(3001);
  fHEvents->SetDirectory(0);
  fList->Add(fHEvents);

  fHEventsTr = new TH1I(GetEventName(true, false), 
			"Number of triggered events",
			fVtxAxis.GetNbins(), 
			fVtxAxis.GetXmin(), 
			fVtxAxis.GetXmax());
  fHEventsTr->SetXTitle("v_{z} [cm]");
  fHEventsTr->SetYTitle("# of events");
  fHEventsTr->SetFillColor(kRed+1);
  fHEventsTr->SetFillStyle(3001);
  fHEventsTr->SetDirectory(0);
  fList->Add(fHEventsTr);

  fHEventsTrVtx = new TH1I(GetEventName(true,true),
			   "Number of events w/trigger and vertex", 
			   fVtxAxis.GetNbins(), 
			   fVtxAxis.GetXmin(), 
			   fVtxAxis.GetXmax());
  fHEventsTrVtx->SetXTitle("v_{z} [cm]");
  fHEventsTrVtx->SetYTitle("# of events");
  fHEventsTrVtx->SetFillColor(kBlue+1);
  fHEventsTrVtx->SetFillStyle(3001);
  fHEventsTrVtx->SetDirectory(0);
  fList->Add(fHEventsTrVtx);

  // Copy axis objects to output 
  TH1* vtxAxis = new TH1D("vtxAxis", "Vertex axis", 
			  fVtxAxis.GetNbins(), 
			  fVtxAxis.GetXmin(), 
			  fVtxAxis.GetXmax());
  TH1* etaAxis = new TH1D("etaAxis", "Eta axis", 
			  fEtaAxis.GetNbins(), 
			  fEtaAxis.GetXmin(), 
			  fEtaAxis.GetXmax());
  fList->Add(vtxAxis);
  fList->Add(etaAxis);

  AliInfo(Form("Initialising sub-routines: %p, %p", 
	       &fInspector, &fTrackDensity));
  fInspector.DefineOutput(fList);
  fInspector.Init(fVtxAxis);
  fTrackDensity.DefineOutput(fList);

  PostData(1, fList);
}
//____________________________________________________________________
void
AliCentralMCCorrectionsTask::UserExec(Option_t*)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  

  // Get the input data - MC event
  AliMCEvent*  mcEvent = MCEvent();
  if (!mcEvent) { 
    AliWarning("No MC event found");
    return;
  }

  // Get the input data - ESD event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) { 
    AliWarning("No ESD event found for input event");
    return;
  }

  //--- Read run information -----------------------------------------
  if (fFirstEvent && esd->GetESDRun()) {
    fInspector.ReadRunDetails(esd);
    
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

    Print();
    fFirstEvent = false;
  }

  // Some variables 
  UInt_t   triggers; // Trigger bits
  Bool_t   lowFlux;  // Low flux flag
  UShort_t iVz;      // Vertex bin from ESD
  Double_t vZ;       // Z coordinate from ESD
  Double_t cent;     // Centrality 
  UShort_t iVzMc;    // Vertex bin from MC
  Double_t vZMc;     // Z coordinate of IP vertex from MC
  Double_t b;        // Impact parameter
  Int_t    nPart;    // Number of participants 
  Int_t    nBin;     // Number of binary collisions 
  Double_t phiR;     // Reaction plane from MC
  UShort_t nClusters;// Number of SPD clusters 
  // Process the data 
  UInt_t retESD = fInspector.Process(esd, triggers, lowFlux, iVz, vZ, 
				     cent, nClusters);
  fInspector.ProcessMC(mcEvent, triggers, iVzMc, vZMc, 
		       b, nPart, nBin, phiR);

  Bool_t isInel   = triggers & AliAODForwardMult::kInel;
  Bool_t hasVtx   = retESD == AliFMDMCEventInspector::kOk;

  // Fill the event count histograms 
  if (isInel)           fHEventsTr->Fill(vZMc);
  if (isInel && hasVtx) fHEventsTrVtx->Fill(vZMc);
  fHEvents->Fill(vZMc);

  // Now find our vertex bin object 
  VtxBin* bin = 0;
  if (iVzMc > 0 && iVzMc <= fVtxAxis.GetNbins()) 
    bin = static_cast<VtxBin*>(fVtxBins->At(iVzMc));
  if (!bin) { 
    // AliError(Form("No vertex bin object @ %d (%f)", iVzMc, vZMc));
    return;
  }

  // Now process our input data and store in own ESD object 
  fTrackDensity.Calculate(*mcEvent, vZMc, *bin->fHits, bin->fPrimary);
  bin->fCounts->Fill(0.5);
}

//____________________________________________________________________
void
AliCentralMCCorrectionsTask::Terminate(Option_t*)
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //

  fList = dynamic_cast<TList*>(GetOutputData(1));
  if (!fList) {
    AliError("No output list defined");
    return;
  }

  // Output list 
  TList* output = new TList;
  output->SetOwner();
  output->SetName(Form("%sResults", GetName()));

  // --- Fill correction object --------------------------------------
  AliCentralCorrSecondaryMap* corr = new AliCentralCorrSecondaryMap;
  corr->SetVertexAxis(fVtxAxis);
  // corr->SetEtaAxis(fEtaAxis);

  TIter     next(fVtxBins);
  VtxBin*   bin = 0;
  UShort_t  iVz = 1;
  while ((bin = static_cast<VtxBin*>(next()))) 
    bin->Finish(fList, output, iVz++, corr);

  output->Add(corr);

  PostData(2, output);
}

//____________________________________________________________________
void
AliCentralMCCorrectionsTask::Print(Option_t* option) const
{
  std::cout << ClassName() << "\n"
	    << "  Vertex bins:      " << fVtxAxis.GetNbins() << '\n'
	    << "  Vertex range:     [" << fVtxAxis.GetXmin() 
	    << "," << fVtxAxis.GetXmax() << "]\n"
	    << "  Eta bins:         " << fEtaAxis.GetNbins() << '\n'
	    << "  Eta range:        [" << fEtaAxis.GetXmin() 
	    << "," << fEtaAxis.GetXmax() << "]\n"
	    << "  # of phi bins:    " << fNPhiBins 
	    << std::endl;
  gROOT->IncreaseDirLevel();
  fInspector.Print(option);
  fTrackDensity.Print(option);
  gROOT->DecreaseDirLevel();
}

//====================================================================
const char*
AliCentralMCCorrectionsTask::VtxBin::BinName(Double_t low, 
					     Double_t high) 
{
  static TString buf;
  buf = Form("vtx%+05.1f_%+05.1f", low, high);
  buf.ReplaceAll("+", "p");
  buf.ReplaceAll("-", "m");
  buf.ReplaceAll(".", "d");
  return buf.Data();
}


//____________________________________________________________________
AliCentralMCCorrectionsTask::VtxBin::VtxBin()
  : fHits(0), 
    fPrimary(0),
    fCounts(0)
{
}
//____________________________________________________________________
AliCentralMCCorrectionsTask::VtxBin::VtxBin(Double_t     low, 
					    Double_t     high, 
					    const TAxis& axis,
					    UShort_t     nPhi)
  : TNamed(BinName(low, high), 
	   Form("%+5.1fcm<v_{z}<%+5.1fcm", low, high)),
    fHits(0), 
    fPrimary(0),
    fCounts(0)
{
  fPrimary = new TH2D("primary", "Primaries", 
		      axis.GetNbins(), axis.GetXmin(), axis.GetXmax(), 
		      nPhi, 0, 2*TMath::Pi());
  fPrimary->SetXTitle("#eta");
  fPrimary->SetYTitle("#varphi [radians]");
  fPrimary->Sumw2();
  fPrimary->SetDirectory(0);

  fHits = static_cast<TH2D*>(fPrimary->Clone("hits"));
  fHits->SetTitle("Hits");
  fHits->SetDirectory(0);

  fCounts = new TH1D("counts", "Counts", 1, 0, 1);
  fCounts->SetXTitle("Events");
  fCounts->SetYTitle("# of Events");
  fCounts->SetDirectory(0);
}

//____________________________________________________________________
AliCentralMCCorrectionsTask::VtxBin::VtxBin(const VtxBin& o)
  : TNamed(o),
    fHits(0),
    fPrimary(0), 
    fCounts(0)
{
  if (o.fHits) {
    fHits = static_cast<TH2D*>(o.fHits->Clone());
    fHits->SetDirectory(0);
  }
  if (o.fPrimary) {
    fPrimary = static_cast<TH2D*>(o.fPrimary->Clone());
    fPrimary->SetDirectory(0);
  }
  if (o.fCounts) {
    fCounts = static_cast<TH1D*>(o.fCounts->Clone());
    fCounts->SetDirectory(0);
  }
}

//____________________________________________________________________
AliCentralMCCorrectionsTask::VtxBin&
AliCentralMCCorrectionsTask::VtxBin::operator=(const VtxBin& o)
{
  if (&o == this) return *this; 
  TNamed::operator=(o);
  fHits    = 0;
  fPrimary = 0;
  fCounts  = 0;
  if (o.fHits) {
    fHits = static_cast<TH2D*>(o.fHits->Clone());
    fHits->SetDirectory(0);
  }
  if (o.fPrimary) {
    fPrimary = static_cast<TH2D*>(o.fPrimary->Clone());
    fPrimary->SetDirectory(0);
  }
  if (o.fCounts) {
    fCounts = static_cast<TH1D*>(o.fCounts->Clone());
    fCounts->SetDirectory(0);
  }
  return *this;
}

//____________________________________________________________________
void
AliCentralMCCorrectionsTask::VtxBin::DefineOutput(TList* l)
{
  TList* d = new TList;
  d->SetName(GetName());
  d->SetOwner();
  l->Add(d);

  d->Add(fHits);
  d->Add(fPrimary);
  d->Add(fCounts);
}
//____________________________________________________________________
void
AliCentralMCCorrectionsTask::VtxBin::Finish(const TList* input, 
					    TList* output, 
					    UShort_t iVz, 
					    AliCentralCorrSecondaryMap* map)
{
  TList* out = new TList;
  out->SetName(GetName());
  out->SetOwner();
  output->Add(out);
  
  TList* l = static_cast<TList*>(input->FindObject(GetName()));
  if (!l) { 
    AliError(Form("List %s not found in %s", GetName(), input->GetName()));
    return;
  }


  TH2D*   hits  = static_cast<TH2D*>(l->FindObject("hits"));
  TH2D*   prim  = static_cast<TH2D*>(l->FindObject("primary"));
  if (!hits || !prim) {
    AliError(Form("Missing histograms: %p, %p", hits, prim));
    return;
  }

  TH2D* h = static_cast<TH2D*>(hits->Clone("bgCorr"));
  h->SetDirectory(0);
  h->Divide(prim);

  map->SetCorrection(iVz, h);

  out->Add(h);
}

//
// EOF
//
