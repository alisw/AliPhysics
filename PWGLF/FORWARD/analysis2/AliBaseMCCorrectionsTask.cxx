// 
// Calculate the corrections in the base regions
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODBaseMult 
// 
// Histograms 
//   
// Corrections used 
// 
#include "AliBaseMCCorrectionsTask.h"
#include "AliBaseMCTrackDensity.h"
#include "AliCorrectionManagerBase.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODForwardMult.h"
#include "AliForwardUtil.h"
#include <TH1.h>
#include <TH2D.h>
#include <TH2F.h>
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
}

//====================================================================
AliBaseMCCorrectionsTask::AliBaseMCCorrectionsTask()
  : AliBaseESDTask(),
    fInspector(),
    fVtxBins(0),
    fHEvents(0), 
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fVtxAxis(),
    fEtaAxis(),
    fUseESDVertex(false),
    fAfterEventSel(false)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
}

//____________________________________________________________________
AliBaseMCCorrectionsTask::AliBaseMCCorrectionsTask(const char* name,
						   AliCorrectionManagerBase* m)
  : AliBaseESDTask(name, "", m),
    fInspector("eventInspector"), 
    fVtxBins(0),
    fHEvents(0), 
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fVtxAxis(10,-10,10), 
    fEtaAxis(200,-4,6),
    fUseESDVertex(false),
    fAfterEventSel(false)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  fBranchNames = 
    "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,"
    "AliESDFMD.,SPDVertex.,PrimaryVertex.";
}

//____________________________________________________________________
void
AliBaseMCCorrectionsTask::SetSatellite(Bool_t sat)
{
  if (!sat) return;
  SetVertexAxis(*AliForwardUtil::MakeFullIpZAxis(20));
  SetEtaAxis(240, -6, 6);
  SetIPzMethod("displaced");
}
//____________________________________________________________________
void
AliBaseMCCorrectionsTask::SetVertexAxis(Int_t nBin, Double_t min, 
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
  if (max < min) { 
    Double_t tmp = min;
    min          = max;
    max          = tmp;
  }
  /*
    if (min < -10) 
    AliWarning(Form("Minimum vertex %f < -10, make sure you want this",min));
    if (max > +10) 
    AliWarning(Form("Minimum vertex %f > +10, make sure you want this",max));
  */
  fVtxAxis.Set(nBin, min, max);
}
//____________________________________________________________________
void
AliBaseMCCorrectionsTask::SetVertexAxis(const TAxis& axis)
{
  // 
  // Set the vertex axis to use
  // 
  // Parameters:
  //    axis Axis
  //
  if (axis.GetXbins() && axis.GetXbins()->GetArray())
    fVtxAxis.Set(axis.GetNbins(),axis.GetXbins()->GetArray());
  else 
    SetVertexAxis(axis.GetNbins(),axis.GetXmin(),axis.GetXmax());
}

//____________________________________________________________________
void
AliBaseMCCorrectionsTask::SetEtaAxis(Int_t nBin, Double_t min, Double_t max)
{
  // 
  // Set the eta axis to use
  // 
  // Parameters:
  //    nBins Number of bins
  //    vzMin Least @f$\eta@f$ 
  //    vzMax Largest @f$\eta@f$ 
  //
  if (max < min) { 
    Double_t tmp = min;
    min          = max;
    max          = tmp;
  }
  if (min < -4) 
    AliWarning(Form("Minimum eta %f < -4, make sure you want this",min));
  if (max > +6) 
    AliWarning(Form("Minimum eta %f > +6, make sure you want this",max));
  fEtaAxis.Set(nBin, min, max);
}
//____________________________________________________________________
void
AliBaseMCCorrectionsTask::SetEtaAxis(const TAxis& axis)
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
AliBaseMCCorrectionsTask::DefineBins(TList* l)
{
  if (!fVtxBins) fVtxBins = new TObjArray(fVtxAxis.GetNbins(), 1);
  if (fVtxBins->GetEntries() > 0) return;

  fVtxBins->SetOwner();
  for (Int_t  i    = 1; i <= fVtxAxis.GetNbins(); i++) { 
    Double_t  low  = fVtxAxis.GetBinLowEdge(i);
    Double_t  high = fVtxAxis.GetBinUpEdge(i);
    VtxBin*   bin  = CreateVtxBin(low, high);
    fVtxBins->AddAt(bin, i);
    bin->CreateOutputObjects(l);
  }
}

//____________________________________________________________________
Bool_t
AliBaseMCCorrectionsTask::Book()
{
  // 
  // Create output objects 
  // 
  //
  DefineBins(fList);
  fNeededCorrections = 0;
  fExtraCorrections  = 0;

  fHEvents = 0;
  TH1* vtxAxis = 0;
  if (fVtxAxis.GetXbins() && fVtxAxis.GetXbins()->GetArray()) {
    vtxAxis = new TH1D("vtxAxis", "Vertex axis", 
		       fVtxAxis.GetNbins(), fVtxAxis.GetXbins()->GetArray()); 
    fHEvents = new TH1I(GetEventName(false,false),
			"Number of all events", 
			fVtxAxis.GetNbins(), fVtxAxis.GetXbins()->GetArray()); 
  }
  else {
    vtxAxis = new TH1D("vtxAxis", "Vertex axis", 
		       fVtxAxis.GetNbins(), 
		       fVtxAxis.GetXmin(), 
		       fVtxAxis.GetXmax());
    fHEvents = new TH1I(GetEventName(false,false),
			"Number of all events", 
			fVtxAxis.GetNbins(), 
			fVtxAxis.GetXmin(), 
			fVtxAxis.GetXmax());
  }
  fHEvents->SetXTitle("v_{z} [cm]");
  fHEvents->SetYTitle("# of events");
  fHEvents->SetFillColor(kBlue+1);
  fHEvents->SetFillStyle(3001);
  fHEvents->SetDirectory(0);
  fList->Add(fHEvents);

  fHEventsTr = static_cast<TH1I*>(fHEvents->Clone(GetEventName(true, false)));
  fHEventsTr->SetTitle("Number of triggered events");
  fHEventsTr->SetFillColor(kRed+1);
  fHEventsTr->SetDirectory(0);
  fList->Add(fHEventsTr);
  
  fHEventsTrVtx = static_cast<TH1I*>(fHEvents->Clone(GetEventName(true,true)));
  fHEventsTrVtx->SetTitle("Number of events w/trigger and vertex");
  fHEventsTrVtx->SetFillColor(kMagenta+1);
  fHEventsTrVtx->SetDirectory(0);
  fList->Add(fHEventsTrVtx);

  // Copy axis objects to output 
  TH1* etaAxis = new TH1D("etaAxis", "Eta axis", 
			  fEtaAxis.GetNbins(), 
			  fEtaAxis.GetXmin(), 
			  fEtaAxis.GetXmax());
  vtxAxis->SetXTitle("IP_{#it{z}}");
  vtxAxis->SetYTitle("dummy");
  etaAxis->SetXTitle("#eta");
  etaAxis->SetYTitle("dummy");
  
  fList->Add(vtxAxis);
  fList->Add(etaAxis);

  AliInfo(Form("Initialising sub-routines: %p, %p", 
	       &fInspector, &GetTrackDensity()));
  GetTrackDensity().CreateOutputObjects(fList);
  
  return true;
}

//____________________________________________________________________
Bool_t
AliBaseMCCorrectionsTask::Event(AliESDEvent& esd)
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
    return false;
  }

  // Some variables 
  UInt_t   triggers  = 0;    // Trigger bits
  Bool_t   lowFlux   = true; // Low flux flag
  UShort_t iVz       = 0;    // Vertex bin from ESD
  TVector3 ip(1024,1024,1000);
  Double_t cent      = -1;   // Centrality 
  UShort_t iVzMc     = 0;    // Vertex bin from MC
  TVector3 ipMc(1024,1024,1000);
  Double_t b         = -1;   // Impact parameter
  Double_t cMC       = -1;   // Centrality estimate from b
  Int_t    nPart     = -1;   // Number of participants 
  Int_t    nBin      = -1;   // Number of binary collisions 
  Double_t phiR      = 100;  // Reaction plane from MC
  UShort_t nClusters = 0;    // Number of SPD clusters 
  // Process the data 
  UInt_t retESD = fInspector.Process(&esd, triggers, lowFlux, iVz, ip, 
				     cent, nClusters);

  Bool_t isAccepted = true;
  if(fAfterEventSel) {
    if (retESD & AliFMDEventInspector::kNoEvent)    isAccepted = false; 
    if (retESD & AliFMDEventInspector::kNoTriggers) isAccepted = false; 
    if (retESD & AliFMDEventInspector::kNoVertex)   isAccepted = false;
    if (retESD & AliFMDEventInspector::kNoFMD)      isAccepted = false;  
    if (!isAccepted) return false;  
    // returns if there is not event , does not pass phys selection ,
    // has no veretx and lack of FMD data.
    // with good veretx outside z range it will continue
  }		


  fInspector.ProcessMC(mcEvent, triggers, iVzMc, ipMc, 
		       b, cMC, nPart, nBin, phiR);
  
  fInspector.CompareResults(ip.Z(),    ipMc.Z(), 
			    cent,      cMC, 
			    b, nPart,  nBin);
  // Only allow NSD events to contribute 
  // if (!(triggers & AliAODForwardMult::kMCNSD)) return false;
  Bool_t isInel   = triggers & AliAODForwardMult::kInel;
  Bool_t hasVtx   = retESD == AliFMDMCEventInspector::kOk;

  
  // Fill the event count histograms 
  if (isInel)           fHEventsTr->Fill(ipMc.Z());
  if (isInel && hasVtx) fHEventsTrVtx->Fill(ipMc.Z());
  fHEvents->Fill(ipMc.Z());

  // Now find our vertex bin object 
  VtxBin*  bin      = 0;
  UShort_t usedZbin = iVzMc; 
  if (fUseESDVertex) usedZbin = iVz;
  // if (!hasVtx)       usedZbin = 0;


  if (usedZbin > 0 && usedZbin <= fVtxAxis.GetNbins()) 
    bin = static_cast<VtxBin*>(fVtxBins->At(usedZbin));
  if (!bin) { 
    AliErrorF("No vertex bin object @ %d (%f)", iVzMc, ipMc.Z());
    return false;
  }

  return ProcessESD(esd, *mcEvent, *bin, ipMc);
}

//____________________________________________________________________
Bool_t
AliBaseMCCorrectionsTask::Finalize()
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  DefineBins(fList);
  CreateCorrections(fResults);

  TIter     next(fVtxBins);
  VtxBin*   bin = 0;
  UShort_t  iVz = 1;
  while ((bin = static_cast<VtxBin*>(next()))) 
    FinalizeVtxBin(bin, iVz++);

  return true;
}

//____________________________________________________________________
void
AliBaseMCCorrectionsTask::Print(Option_t* option) const
{
  AliBaseESDTask::Print(option);
  std::cout << "  Vertex bins:      " << fVtxAxis.GetNbins() << '\n'
	    << "  Vertex range:     [" << fVtxAxis.GetXmin() 
	    << "," << fVtxAxis.GetXmax() << "]\n"
	    << "  Eta bins:         " << fEtaAxis.GetNbins() << '\n'
	    << "  Eta range:        [" << fEtaAxis.GetXmin() 
	    << "," << fEtaAxis.GetXmax() << "]"
	    << std::endl;
}

//====================================================================
const char*
AliBaseMCCorrectionsTask::VtxBin::BinName(Double_t low, 
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
AliBaseMCCorrectionsTask::VtxBin::VtxBin()
  : TNamed(),
    fPrimary(0),
    fCounts(0)
{
}
//____________________________________________________________________
AliBaseMCCorrectionsTask::VtxBin::VtxBin(Double_t     low, 
					 Double_t     high, 
					 const TAxis& axis,
					 UShort_t     nPhi)
  : TNamed(BinName(low, high), 
	   Form("%+5.1fcm<v_{z}<%+5.1fcm", low, high)),
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

  fCounts = new TH1D("counts", "Counts", 1, 0, 1);
  fCounts->SetXTitle("Events");
  fCounts->SetYTitle("# of Events");
  fCounts->SetDirectory(0);
}


//____________________________________________________________________
TList*
AliBaseMCCorrectionsTask::VtxBin::CreateOutputObjects(TList* l)
{
  TList* d = new TList;
  d->SetName(GetName());
  d->SetOwner();
  l->Add(d);

  d->Add(fPrimary);
  d->Add(fCounts);

  return d;
}


//
// EOF
//
