// 
// Calculate the corrections in the forward regions
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
#include "AliForwardMCCorrectionsTask.h"
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
#include "AliFMDStripIndex.h"
#include "AliFMDCorrSecondaryMap.h"
#include <TH1.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TList.h>
#include <TROOT.h>
#include <iostream>

//====================================================================
namespace {
  const char* GetEventName(Bool_t tr, Bool_t vtx) 
  {
    return Form("nEvents%s%s", (tr ? "Tr" : ""), (vtx ? "Vtx" : ""));
  }
#if 0
  const char* GetHitsName(UShort_t d, Char_t r) 
  {
    return Form("hitsFMD%d%c", d, r);
  }
  const char* GetStripsName(UShort_t d, Char_t r)
  {
    return Form("stripsFMD%d%c", d, r);
  }
  const char* GetPrimaryName(Char_t r, Bool_t trVtx)
  {
    return Form("primaries%s%s", 
		((r == 'I' || r == 'i') ? "Inner" : "Outer"), 
		(trVtx ? "TrVtx" : "All"));
  }
#endif
}

//====================================================================
AliForwardMCCorrectionsTask::AliForwardMCCorrectionsTask()
  : AliAnalysisTaskSE(),
    fInspector(),
    fTrackDensity(),
    fESDFMD(),
    fVtxBins(0),
    fFirstEvent(true),
    fHEvents(0), 
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fVtxAxis(),
    fEtaAxis(),
    fList()
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
}

//____________________________________________________________________
AliForwardMCCorrectionsTask::AliForwardMCCorrectionsTask(const char* name)
  : AliAnalysisTaskSE(name),
    fInspector("eventInspector"), 
    fTrackDensity("trackDensity"),
    fESDFMD(),
    fVtxBins(0),
    fFirstEvent(true),
    fHEvents(0), 
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fVtxAxis(10,-10,10), 
    fEtaAxis(200,-4,6),
    fList()
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
AliForwardMCCorrectionsTask::AliForwardMCCorrectionsTask(const AliForwardMCCorrectionsTask& o)
  : AliAnalysisTaskSE(o),
    fInspector(o.fInspector),
    fTrackDensity(),
    fESDFMD(o.fESDFMD),
    fVtxBins(0),
    fFirstEvent(o.fFirstEvent),
    fHEvents(o.fHEvents), 
    fHEventsTr(o.fHEventsTr), 
    fHEventsTrVtx(o.fHEventsTrVtx),
    fVtxAxis(10,-10,10), 
    fEtaAxis(200,-4,6),
    fList(o.fList)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
}

//____________________________________________________________________
AliForwardMCCorrectionsTask&
AliForwardMCCorrectionsTask::operator=(const AliForwardMCCorrectionsTask& o)
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
  fESDFMD            = o.fESDFMD;
  fVtxBins           = o.fVtxBins;
  fFirstEvent        = o.fFirstEvent;
  fHEvents           = o.fHEvents;
  fHEventsTr         = o.fHEventsTr;
  fHEventsTrVtx      = o.fHEventsTrVtx;
  SetVertexAxis(o.fVtxAxis);
  SetEtaAxis(o.fEtaAxis);

  return *this;
}

//____________________________________________________________________
void
AliForwardMCCorrectionsTask::Init()
{
  // 
  // Initialize the task 
  // 
  //
}

//____________________________________________________________________
void
AliForwardMCCorrectionsTask::SetVertexAxis(Int_t nBin, Double_t min, 
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
AliForwardMCCorrectionsTask::SetVertexAxis(const TAxis& axis)
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
AliForwardMCCorrectionsTask::SetEtaAxis(Int_t nBin, Double_t min, Double_t max)
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
AliForwardMCCorrectionsTask::SetEtaAxis(const TAxis& axis)
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
AliForwardMCCorrectionsTask::DefineBins(TList* l)
{
  if (!fVtxBins) fVtxBins = new TObjArray(fVtxAxis.GetNbins(), 1);
  if (fVtxBins->GetEntries() > 0) return;

  fVtxBins->SetOwner();
  for (Int_t  i    = 1; i <= fVtxAxis.GetNbins(); i++) { 
    Double_t  low  = fVtxAxis.GetBinLowEdge(i);
    Double_t  high = fVtxAxis.GetBinUpEdge(i);
    VtxBin*   bin  = new VtxBin(low, high, fEtaAxis);
    fVtxBins->AddAt(bin, i);
    bin->DefineOutput(l);
  }
}

//____________________________________________________________________
void
AliForwardMCCorrectionsTask::UserCreateOutputObjects()
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
AliForwardMCCorrectionsTask::UserExec(Option_t*)
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

  // Clear our ESD object 
  fESDFMD.Clear();

  // Get FMD data 
  AliESDFMD* esdFMD = esd->GetFMDData();
  
  // Now process our input data and store in own ESD object 
  fTrackDensity.Calculate(*esdFMD, *mcEvent, vZMc, fESDFMD, bin->fPrimary);
  bin->fCounts->Fill(0.5);

  // And then bin the data in our vtxbin 
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      UShort_t    ns= (q == 0 ?  20 :  40);
      UShort_t    nt= (q == 0 ? 512 : 256);
      TH2D*       h = bin->fHists.Get(d,r);

      for (UShort_t s=0; s<ns; s++) { 
	for (UShort_t t=0; t<nt; t++) {
	  Float_t mult = fESDFMD.Multiplicity(d,r,s,t);
	  
	  if (mult == 0 || mult > 20) continue;

	  Float_t phi = fESDFMD.Phi(d,r,s,t) / 180 * TMath::Pi();
	  Float_t eta = fESDFMD.Eta(d,r,s,t);
	  h->Fill(eta,phi,mult);
	} // for t
      } // for s 
    } // for q 
  } // for d
}

//____________________________________________________________________
void
AliForwardMCCorrectionsTask::Terminate(Option_t*)
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
  AliFMDCorrSecondaryMap* corr = new AliFMDCorrSecondaryMap;
  corr->SetVertexAxis(fVtxAxis);
  corr->SetEtaAxis(fEtaAxis);
  
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
AliForwardMCCorrectionsTask::Print(Option_t* option) const
{
  std::cout << ClassName() << "\n"
	    << "  Vertex bins:      " << fVtxAxis.GetNbins() << '\n'
	    << "  Vertex range:     [" << fVtxAxis.GetXmin() 
	    << "," << fVtxAxis.GetXmax() << "]\n"
	    << "  Eta bins:         " << fEtaAxis.GetNbins() << '\n'
	    << "  Eta range:        [" << fEtaAxis.GetXmin() 
	    << "," << fEtaAxis.GetXmax() << "]"
	    << std::endl;
  gROOT->IncreaseDirLevel();
  fInspector.Print(option);
  fTrackDensity.Print(option);
  gROOT->DecreaseDirLevel();
}

//====================================================================
const char*
AliForwardMCCorrectionsTask::VtxBin::BinName(Double_t low, 
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
AliForwardMCCorrectionsTask::VtxBin::VtxBin()
  : fHists(), 
    fPrimary(0),
    fCounts(0)
{
}
//____________________________________________________________________
AliForwardMCCorrectionsTask::VtxBin::VtxBin(Double_t low, 
					    Double_t high, 
					    const TAxis& axis)
  : TNamed(BinName(low, high), 
	   Form("%+5.1fcm<v_{z}<%+5.1fcm", low, high)),
    fHists(), 
    fPrimary(0),
    fCounts(0)
{
  fHists.Init(axis);

  fPrimary = new TH2D("primary", "Primaries", 
		      axis.GetNbins(), axis.GetXmin(), axis.GetXmax(), 
		      40, 0, 2*TMath::Pi());
  fPrimary->SetXTitle("#eta");
  fPrimary->SetYTitle("#varphi [radians]");
  fPrimary->Sumw2();
  fPrimary->SetDirectory(0);

  fCounts = new TH1D("counts", "Counts", 1, 0, 1);
  fCounts->SetXTitle("Events");
  fCounts->SetYTitle("# of Events");
  fCounts->SetDirectory(0);
}

//____________________________________________________________________
AliForwardMCCorrectionsTask::VtxBin::VtxBin(const VtxBin& o)
  : TNamed(o),
    fHists(o.fHists),
    fPrimary(0), 
    fCounts(0)
{
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
AliForwardMCCorrectionsTask::VtxBin&
AliForwardMCCorrectionsTask::VtxBin::operator=(const VtxBin& o)
{
  if (&o == this) return *this;
  TNamed::operator=(o);
  fHists   = o.fHists;
  fPrimary = 0;
  fCounts  = 0;
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
AliForwardMCCorrectionsTask::VtxBin::DefineOutput(TList* l)
{
  TList* d = new TList;
  d->SetName(GetName());
  d->SetOwner();
  l->Add(d);

  d->Add(fHists.fFMD1i);
  d->Add(fHists.fFMD2i);
  d->Add(fHists.fFMD2o);
  d->Add(fHists.fFMD3i);
  d->Add(fHists.fFMD3o);

  d->Add(fPrimary);
  d->Add(fCounts);
}

//____________________________________________________________________
TH2D*
AliForwardMCCorrectionsTask::VtxBin::MakeBg(const TH2D* hits, 
					    const TH2D* primary) const
{
  TH2D* h = static_cast<TH2D*>(hits->Clone());
  h->SetDirectory(0);
  TString n(h->GetName());
  n.ReplaceAll("_cache", "");
  h->SetName(n);
  h->Divide(primary);
  
  return h;
}
  
//____________________________________________________________________
void
AliForwardMCCorrectionsTask::VtxBin::Finish(const TList* input, 
					    TList* output, 
					    UShort_t iVz,
					    AliFMDCorrSecondaryMap* map)
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

  TH2D*   fmd1i = static_cast<TH2D*>(l->FindObject("FMD1I_cache"));
  TH2D*   fmd2i = static_cast<TH2D*>(l->FindObject("FMD2I_cache"));
  TH2D*   fmd2o = static_cast<TH2D*>(l->FindObject("FMD2O_cache"));
  TH2D*   fmd3i = static_cast<TH2D*>(l->FindObject("FMD3I_cache"));
  TH2D*   fmd3o = static_cast<TH2D*>(l->FindObject("FMD3O_cache"));
  TH2D*   primO = static_cast<TH2D*>(l->FindObject("primary"));
  if (!fmd1i || !fmd2i || !fmd2o || !fmd3i || !fmd3o || !primO) {
    AliError(Form("Missing histogram(s): %p,%p,%p,%p,%p,%p",
		  fmd1i, fmd2i, fmd2o, fmd3i, fmd3o, primO));
    return;
  }

  // Half coverage in phi for inners
  TH2D*   primI = static_cast<TH2D*>(primO->Clone());
  primI->SetDirectory(0);
  primI->RebinY(2); 

  TH2D* bg1i = MakeBg(fmd1i, primI);
  TH2D* bg2i = MakeBg(fmd2i, primI);
  TH2D* bg2o = MakeBg(fmd2o, primO);
  TH2D* bg3i = MakeBg(fmd3i, primI);
  TH2D* bg3o = MakeBg(fmd3o, primO);
  map->SetCorrection(1, 'I', iVz, bg1i);
  map->SetCorrection(2, 'I', iVz, bg2i);
  map->SetCorrection(2, 'O', iVz, bg2o);
  map->SetCorrection(3, 'I', iVz, bg3i);
  map->SetCorrection(3, 'O', iVz, bg3o);
  out->Add(bg1i);
  out->Add(bg2i);
  out->Add(bg2o);
  out->Add(bg3i);
  out->Add(bg3o);
 
}

//
// EOF
//
