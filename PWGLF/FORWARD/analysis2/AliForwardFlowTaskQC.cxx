//
// Calculate flow in the forward and central regions using the Q cumulants method.
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root
//
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include <TMath.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TParameter.h>
#include <TGraph.h>
#include "AliLog.h"
#include "AliForwardFlowTaskQC.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"
#include "AliForwardUtil.h"

ClassImp(AliForwardFlowTaskQC)
#if 0
; // For emacs 
#endif

AliForwardFlowTaskQC::AliForwardFlowTaskQC()
  : AliAnalysisTaskSE(),
    fVtxAxis(),         // Axis to contorl vertex binning
    fFMDCut(-1),        // FMD sigma cut
    fSPDCut(-1),        // SPD sigma cut
    fFlowFlags(kSymEta),// Flow flags
    fEtaGap(2.),        // Eta gap value
    fBinsFMD(),         // List with FMD flow histos
    fBinsSPD(),         // List with SPD flow histos
    fSumList(0),	// Event sum list
    fOutputList(0),	// Result output list
    fAOD(0),		// AOD input event
    fV(),               // Flow moments
    fVtx(1111),	        // Z vertex coordinate
    fCent(-1),		// Centrality
    fHistCent(),        // Histo for centrality
    fHistVertexSel()    // Histo for selected vertices
{
  // 
  // Default constructor
  //
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const char* name) 
  : AliAnalysisTaskSE(name),
    fVtxAxis(),         // Axis to contorl vertex binning
    fFMDCut(-1),        // FMD sigma cut
    fSPDCut(-1),        // SPD sigma cut
    fFlowFlags(kSymEta),// Flow flags
    fEtaGap(2.),        // Eta gap value
    fBinsFMD(),         // List with FMD flow histos
    fBinsSPD(),         // List with SPD flow histos
    fSumList(0),        // Event sum list           
    fOutputList(0),     // Result output list       
    fAOD(0),	        // AOD input event          
    fV(),               // Flow moments
    fVtx(1111),         // Z vertex coordinate      
    fCent(-1),          // Centrality               
    fHistCent(),        // Histo for centrality
    fHistVertexSel()    // Histo for selected vertices
{
  // 
  // Constructor
  //
  // Parameters:
  //  name: Name of task
  //
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const AliForwardFlowTaskQC& o)
  : AliAnalysisTaskSE(o),
    fVtxAxis(o.fVtxAxis),              // Axis to contorl vertex binning
    fFMDCut(o.fFMDCut),                // FMD sigma cut
    fSPDCut(o.fSPDCut),                // SPD sigma cut
    fFlowFlags(o.fFlowFlags),          // Flow flags
    fEtaGap(o.fEtaGap),                // Eta gap value
    fBinsFMD(),                        // List with FMD flow histos
    fBinsSPD(),                        // List with SPD flow histos
    fSumList(o.fSumList),              // Event sum list           
    fOutputList(o.fOutputList),        // Result output list       
    fAOD(o.fAOD),	               // AOD input event          
    fV(o.fV),                          // Flow moments
    fVtx(o.fVtx),                      // Z vertex coordinate      
    fCent(o.fCent),	               // Centrality
    fHistCent(o.fHistCent),            // Histo for centrality
    fHistVertexSel(o.fHistVertexSel)   // Histo for selected vertices
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
}
//_____________________________________________________________________
AliForwardFlowTaskQC&
AliForwardFlowTaskQC::operator=(const AliForwardFlowTaskQC& o)
{
  // 
  // Assignment operator 
  //
  if (&o == this) return *this;
  fVtxAxis       = o.fVtxAxis;
  fFMDCut        = o.fFMDCut;
  fSPDCut        = o.fSPDCut;
  fFlowFlags     = o.fFlowFlags;
  fEtaGap        = o.fEtaGap;
  fSumList       = o.fSumList;
  fOutputList    = o.fOutputList;
  fAOD           = o.fAOD;
  fV             = o.fV;
  fVtx           = o.fVtx;
  fCent          = o.fCent;
  fHistCent      = o.fHistCent;
  fHistVertexSel = o.fHistVertexSel;

  return *this;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::UserCreateOutputObjects()
{
  //
  // Create output objects
  //
  InitVertexBins();
  InitHists();
  PrintFlowSetup();

  PostData(1, fSumList);
  PostData(2, fOutputList);

}
//_____________________________________________________________________
void AliForwardFlowTaskQC::InitVertexBins()
{
  // 
  // Init vertexbin objects for FMD and SPD, and add them to the lists
  //
  Int_t moment = 0;
  for(UShort_t n = 0; n < fV.GetSize(); n++) {
    moment = fV.At(n);
    for (Int_t v = 1; v <= fVtxAxis->GetNbins(); v++) {
      Int_t vL = Int_t(fVtxAxis->GetBinLowEdge(v));
      Int_t vH = Int_t(fVtxAxis->GetBinUpEdge(v));
      fBinsFMD.Add(new VertexBin(vL, vH, moment, "FMD", fFlowFlags, fFMDCut, fEtaGap));
      fBinsSPD.Add(new VertexBin(vL, vH, moment, "SPD", fFlowFlags, fSPDCut, fEtaGap));
    }
  }
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::InitHists()
{
  //
  // Init histograms and add vertex bin histograms to the sum list
  //
  
  if (!fSumList)
    fSumList = new TList();
  fSumList->SetName("Sums");
  fSumList->SetOwner();

  if (!fVtxAxis) fVtxAxis = new TAxis(20, -10, 10);
  fVtxAxis->SetName("VtxAxis");
  fHistCent         = new TH1D("hCent", "Centralities", 100, 0, 100);
  fHistVertexSel    = new TH1D("hVertexSel", "Selected vertices", fVtxAxis->GetNbins(), fVtxAxis->GetXmin(), fVtxAxis->GetXmax());

  TList* dList = new TList();
  dList->SetName("Diagnostics");
//  dList->Add(fVtxAxis);
  dList->Add(fHistCent);
  dList->Add(fHistVertexSel);
  fSumList->Add(dList);

  TIter nextFMD(&fBinsFMD);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(nextFMD()))) {
    bin->AddOutput(fSumList);
  }
  TIter nextSPD(&fBinsSPD);
  while ((bin = static_cast<VertexBin*>(nextSPD()))) {
    bin->AddOutput(fSumList);
  }
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::UserExec(Option_t */*option*/)
{
  //
  // Calls the analyze function - called every event
  //
  // Parameters:
  //  option: Not used
  //
  
  Analyze();
  
  PostData(1, fSumList);

  return;
}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::Analyze()
{
  // 
  // Load FMD and SPD objects from aod tree and call the cumulants 
  // calculation for the correct vertexbin
  //

  // Reset data members
  fCent = -1;
  fVtx = 1111;

  // Get input event
//  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fAOD = dynamic_cast<AliAODEvent*>(AliForwardUtil::GetAODEvent(this));
  if (!fAOD) return kFALSE;

  const AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
  const AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClusters"));
  if (!aodfmult) return kFALSE;
  
  // Check event for triggers, get centrality, vtx etc.
  if (!CheckEvent(aodfmult)) return kFALSE;
  Int_t vtx = fVtxAxis->FindBin(fVtx)-1;

  // If everything is OK: get histos and run analysis
  const TH2D& fmddNdetadphi = aodfmult->GetHistogram();
  if ((fFlowFlags & kEtaGap)) {
    FillVtxBinListEtaGap(fBinsFMD, fmddNdetadphi, fmddNdetadphi, vtx);
  } else {
    FillVtxBinList(fBinsFMD, fmddNdetadphi, vtx);
  }

  if (aodcmult) {
    const TH2D& spddNdetadphi = aodcmult->GetHistogram();
    if ((fFlowFlags & kEtaGap)) {
      FillVtxBinListEtaGap(fBinsSPD, fmddNdetadphi, spddNdetadphi, vtx);
    } else {
      FillVtxBinList(fBinsSPD, spddNdetadphi, vtx);
    }
  }

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::FillVtxBinList(const TList& list, const TH2D& h, Int_t vtx) const
{
  //
  // Loops over list of VtxBins, fills hists of bins for current vertex
  // and runs analysis on those bins
  //
  // Parameters:
  //  list: list of VtxBins
  //  h:    dN/detadphi histogram
  //  vBin: current vertex bin
  //
  // return true on success
  //
  VertexBin* bin = 0;
  Int_t i = 0;
  Int_t nVtxBins = fVtxAxis->GetNbins();

  while ((bin = static_cast<VertexBin*>(list.At(vtx+(nVtxBins*i))))) {
    Bool_t skipFourP = !bin->FillHists(h, fCent, kFillBoth);
    bin->CumulantsAccumulate(fCent, skipFourP);
    i++;
  }

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::FillVtxBinListEtaGap(const TList& list, const TH2D& href, 
                                                  const TH2D& hdiff, Int_t vtx) const
{
  //
  // Loops over list of VtxBins, fills hists of bins for current vertex
  // and runs analysis on those bins
  //
  // Parameters:
  //  list: list of VtxBins
  //  h:    dN/detadphi histogram
  //  vBin: current vertex bin
  //
  // return true on success
  //
  VertexBin* bin = 0;
  Int_t i = 0;
  Int_t nVtxBins = fVtxAxis->GetNbins();

  while ((bin = static_cast<VertexBin*>(list.At(vtx+(nVtxBins*i))))) {
    bin->FillHists(href, fCent, kFillRef);
    bin->FillHists(hdiff, fCent, kFillDiff);
    bin->CumulantsAccumulate(fCent);
    i++;
  }

  return kTRUE;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::Terminate(Option_t */*option*/)
{
  //
  // Calls the finalize function, done at the end of the analysis
  //
  // Parameters:
  //  option: Not used
  //

  // Make sure pointers are set to the correct lists
  fSumList = dynamic_cast<TList*> (GetOutputData(1));
  if(!fSumList) {
    AliError("Could not retrieve TList fSumList"); 
    return; 
  }
  if (!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("Results");
  fOutputList->SetOwner();

  if ((fFlowFlags & kEtaGap)) {
    TParameter<Double_t>* etaGap = new TParameter<Double_t>("EtaGap", fEtaGap);
    fOutputList->Add(etaGap);
  }

  // We make summary histograms of accepted events.
  TList* list = 0;
  TH3D* hist = 0;
  TH1D* cent = 0;
  for (Int_t i = 1; i < fSumList->GetEntries(); i++) {
    list = dynamic_cast<TList*>(fSumList->At(i));
    if (!list) continue;
    hist = dynamic_cast<TH3D*>(list->At(0));
    if (!hist) continue;
    const char* histname = hist->GetName();
    TString name = "";
    for (Int_t j = 0; ; j++) {
      if (histname[j] == 'v') break;
      name += histname[j];
    }
    if ((fFlowFlags & kEtaGap)) name += "_etaGap";
    cent = (TH1D*)fOutputList->FindObject(name.Data());
    if (!cent) {
      cent = new TH1D(name.Data(), name.Data(), hist->GetNbinsY(), hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax());
      cent->GetXaxis()->Set(hist->GetNbinsY(), hist->GetYaxis()->GetXbins()->GetArray());
      fOutputList->Add(cent);
    }
    for (Int_t k = 1; k <= hist->GetNbinsY(); k++) {
      Double_t centrality = hist->GetYaxis()->GetBinCenter(k);
      Double_t events = hist->GetBinContent(0, k, 0);
      cent->Fill(centrality, events);
    }
  }

  // Run finalize on VertexBins
  Finalize();

  // Collect centralities
  MakeCentralityHists(fOutputList);
  TList* vtxList = (TList*)fOutputList->FindObject("vtxList");
  if (vtxList) MakeCentralityHists(vtxList);

  PostData(2, fOutputList);

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::AddFlowMoment(Short_t n)
{
    //
    // Add a flow moment to be calculated
    //
    Int_t size = fV.GetSize();
    fV.Set(size+1);
    fV.AddAt(n, size);
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::Finalize()
{
  //
  // Finalize command, called by Terminate()
  //

  // Reinitiate vertex bins if Terminate is called separately!
  if (fBinsFMD.GetEntries() == 0) InitVertexBins();

  // Iterate over all vertex bins objects and finalize cumulants
  // calculations
  EndVtxBinList(fBinsFMD);
  EndVtxBinList(fBinsSPD);

  return;
} 
//_____________________________________________________________________
void AliForwardFlowTaskQC::EndVtxBinList(const TList& list) const
{
  //
  // Loop over VertexBin list and call terminate on each 
  //
  // Parameters:
  //  list VertexBin list
  //
  TIter next(&list);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(next()))) {
    bin->CumulantsTerminate(fSumList, fOutputList);
  }
  return;
}
// _____________________________________________________________________
void AliForwardFlowTaskQC::MakeCentralityHists(TList* list)
{
  //
  // Loop over a list containing a TProfile2D with flow results and project
  // and project to TH1's in specific centrality bins
  //
  // Parameters:
  //  list Flow results list
  //
  TProfile2D* hist2D = 0;
  TList* centList = 0;
  TH1D* hist1D = 0;
  TObject* helper = 0;
  TIter nextProfile(list);
  while ((helper = dynamic_cast<TObject*>(nextProfile()))) {
    if (!(hist2D = dynamic_cast<TProfile2D*>(helper))) continue;
    for (Int_t cBin = 1; cBin <= hist2D->GetNbinsY(); cBin++) {
      Int_t cMin = Int_t(hist2D->GetYaxis()->GetBinLowEdge(cBin));
      Int_t cMax = Int_t(hist2D->GetYaxis()->GetBinUpEdge(cBin));
      TString name = Form("cent_%d-%d", cMin, cMax);
      centList = (TList*)list->FindObject(name.Data());
      if (!centList) { 
	centList = new TList();
	centList->SetName(name.Data());
	list->Add(centList);
      }
      hist1D = hist2D->ProjectionX(Form("%s_%s", hist2D->GetName(), name.Data()), 
                                         cBin, cBin, "E");
      hist1D->SetTitle(hist1D->GetName());
      centList->Add(hist1D);
    }
  }
}
// _____________________________________________________________________
Bool_t AliForwardFlowTaskQC::CheckEvent(const AliAODForwardMult* aodfm) 
{
  // 
  // Function to check that and AOD event meets the cuts
  //
  // Parameters: 
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  // Returns false if there is no trigger or if the centrality or vertex
  // is out of range. Otherwise true.
  //

  // First check for trigger
  if (!CheckTrigger(aodfm)) return kFALSE;

  // Then check for centrality
  if (!GetCentrality(aodfm)) return kFALSE;

  // And finally check for vertex
  if (!GetVertex(aodfm)) return kFALSE;

  // Ev. accepted - filling diag. hists
  fHistCent->Fill(fCent);
  fHistVertexSel->Fill(fVtx);
  
  return kTRUE;
}
// _____________________________________________________________________
Bool_t AliForwardFlowTaskQC::CheckTrigger(const AliAODForwardMult* aodfm) const 
{
  //
  // Function to look for a trigger string in the event.
  //
  // Parameters: 
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  // Returns true if offline trigger is present
  //
  return aodfm->IsTriggerBits(AliAODForwardMult::kOffline);
}
// _____________________________________________________________________
Bool_t AliForwardFlowTaskQC::GetCentrality(const AliAODForwardMult* aodfm) 
{
  //
  // Function to look get centrality of the event.
  //
  // Parameters: 
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  // Returns true if centrality determination is present
  //
  if (aodfm->HasCentrality()) {
    fCent = (Double_t)aodfm->GetCentrality();
    if (0. >= fCent || fCent >= 100.) return kFALSE;
  }
  else fCent = 97.5;

  return kTRUE;
}
// _____________________________________________________________________
Bool_t AliForwardFlowTaskQC::GetVertex(const AliAODForwardMult* aodfm) 
{
  //
  // Function to look for vertex determination in the event.
  //
  // Parameters: 
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  // Returns true if vertex is determined
  //
  fVtx = aodfm->GetIpZ();
  if (fVtx < fVtxAxis->GetXmin() || fVtx > fVtxAxis->GetXmax()) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________
AliForwardFlowTaskQC::VertexBin::VertexBin()
  : TNamed(),
    fMoment(0),      // Flow moment for this vertexbin
    fVzMin(0),       // Vertex z-coordinate min
    fVzMax(0),       // Vertex z-coordinate max
    fType(),         // Data type name e.g., FMD/SPD/FMDTR/SPDTR/MC
    fFlags(0),       // Use forward-backward symmetry, if detector allows it
    fSigmaCut(-1),   // Sigma cut to remove outlier events
    fEtaGap(2.),     // Eta gap value
    fCumuRef(),      // Histogram for reference flow
    fCumuDiff(),     // Histogram for differential flow
    fCumuHist(),     // Sum histogram for cumulants
    fdNdedpAcc(),    // Diagnostics histogram to make acc. maps
    fOutliers(),     // Histogram for sigma distribution
    fDebug()         // Debug level
{
  //
  // Default constructor
  //
}
//_____________________________________________________________________
AliForwardFlowTaskQC::VertexBin::VertexBin(Int_t vLow, Int_t vHigh, 
                                           UShort_t moment, TString name,
                                           UShort_t flags, Double_t cut,
                                           Double_t etaGap)
  : TNamed("", ""),
    fMoment(moment),  // Flow moment for this vertexbin
    fVzMin(vLow),     // Vertex z-coordinate min
    fVzMax(vHigh),    // Vertex z-coordinate max
    fType(name),      // Data type name e.g., FMD/SPD/FMDTR/SPDTR/MC
    fFlags(flags),    // Use forward-backward symmetry, if detector allows it
    fSigmaCut(cut),   // Sigma cut to remove outlier events
    fEtaGap(etaGap),  // Eta gap value
    fCumuRef(),       // Histogram for reference flow
    fCumuDiff(),      // Histogram for differential flow
    fCumuHist(),      // Sum histogram for cumulants
    fdNdedpAcc(),     // Diagnostics histogram to make acc. maps
    fOutliers(),      // Histogram for sigma distribution
    fDebug(0)         // Debug level
{
  //
  // Constructor
  //
  // Parameters
  //  vLow: min z-coordinate
  //  vHigh: max z-coordinate
  //  moment: flow moment
  //  name: data type name (FMD/SPD/FMDTR/SPDTR/MC)
  //  sym: data is symmetric in eta
  //
  fType.ToUpper();

  SetName(Form("%svertexBin%d_%d_%d%s", fType.Data(), moment, vLow, vHigh, ((fFlags & kEtaGap) ? "_etaGap" : "")));
  SetTitle(Form("%svertexBin%d_%d_%d%s", fType.Data(), moment, vLow, vHigh, ((fFlags & kEtaGap) ? "_etaGap" : "")));

  fDebug = AliAnalysisManager::GetAnalysisManager()->GetDebugLevel();

}
//_____________________________________________________________________
AliForwardFlowTaskQC::VertexBin&
AliForwardFlowTaskQC::VertexBin::operator=(const AliForwardFlowTaskQC::VertexBin& o)
{
  //
  // Assignment operator
  //
  // Parameters
  //  o: AliForwardFlowTaskQC::VertexBin
  //
  if (&o == this) return *this;
  fType         = o.fType;
  fCumuRef      = o.fCumuRef;
  fCumuDiff     = o.fCumuDiff;
  fCumuHist     = o.fCumuHist;
  fdNdedpAcc    = o.fdNdedpAcc;
  fOutliers     = o.fOutliers;
  fDebug        = o.fDebug;

  return *this;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::AddOutput(TList* outputlist)
{
  //
  // Add histograms to outputlist
  //
  // Parameters
  //  outputlist: list of histograms
  //

  // First we try to find an outputlist for this vertexbin
  TList* list = (TList*)outputlist->FindObject(Form("%svertex_%d_%d%s", fType.Data(), fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")));

  // If it doesn't exist we make one
  if (!list) {
    list = new TList();
    list->SetName(Form("%svertex_%d_%d%s", fType.Data(), fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")));
    outputlist->Add(list);
  }

  // We initiate the reference histogram
  fCumuRef = new TH2D(Form("%s_v%d_%d_%d%s_ref", fType.Data(), fMoment, fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")),
                        Form("%s_v%d_%d_%d%s_ref", fType.Data(), fMoment, fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")),
                        2, -6., 6., 5, 0.5, 5.5);
  fCumuRef->Sumw2();
  //list->Add(fCumuRef);

  // We initiate the differential histogram
  fCumuDiff = new TH2D(Form("%s_v%d_%d_%d%s_diff", fType.Data(), fMoment, fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")),
                       Form("%s_v%d_%d_%d%s_diff", fType.Data(), fMoment, fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")),
                       48, -6., 6., 5, 0.5, 5.5);
  fCumuDiff->Sumw2();
  //list->Add(fCumuDiff);

  // Initiate the cumulant sum histogram
  fCumuHist = new TH3D(Form("%sv%d_vertex_%d_%d%s_cumu", fType.Data(), fMoment, fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")),
                       Form("%sv%d_vertex_%d_%d%s_cumu", fType.Data(), fMoment, fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")),
                       48, -6., 6., 20, 0., 100., 26, 0.5, 26.5);
  fCumuHist->Sumw2();
  SetupCentAxis(fCumuHist->GetYaxis());

  list->Add(fCumuHist);

  // We check for diagnostics histograms (only done per type and moment, not vertexbin)
  // If they are not found we create them.
  TList* dList = (TList*)outputlist->FindObject("Diagnostics");
  if (!dList) AliFatal("No diagnostics list found, what kind of game are you running here?!?!");

  // Acceptance hists are shared over all moments
  fdNdedpAcc = (TH2F*)dList->FindObject(Form("h%sdNdedpAcc_%d_%d%s", fType.Data(), fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")));
  if (!fdNdedpAcc) {
    fdNdedpAcc = new TH2F(Form("h%sdNdedpAcc_%d_%d%s", fType.Data(), fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")), 
                          Form("%s acceptance map for %d cm < v_{z} < %d cm", fType.Data(), fVzMin, fVzMax),
                          48, -6, 6, 20, 0, TMath::TwoPi());
    fdNdedpAcc->Sumw2();
    dList->Add(fdNdedpAcc);
  }

  if (!(fFlags & kEtaGap)) {
    fOutliers = new TH2F(Form("hOutliers_%s_v%d_%d_%d", fType.Data(), fMoment, fVzMin, fVzMax),
			Form("Maximum #sigma from mean N_{ch} pr. bin - %s v_{%d}, %d < v_{z} < %d",
			fType.Data(), fMoment, fVzMin, fVzMax), 
			20, 0., 100., 500, 0., (fType.Contains("MC") ? 15. : 5.));
    dList->Add(fOutliers);
  }

}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::VertexBin::FillHists(const TH2D& dNdetadphi, Double_t cent, EFillFlow mode) 
{
  // 
  // Fill reference and differential eta-histograms
  //
  // Parameters:
  //  dNdetadphi: 2D histogram with input data
  //

  if (!fCumuRef) AliFatal("You have not called AddOutput() - Terminating!");

  Bool_t useEvent = kTRUE;

  // Fist we reset histograms
  if ((mode & kFillRef))  fCumuRef->Reset();
  if ((mode & kFillDiff)) fCumuDiff->Reset();

  // Numbers to cut away bad events and acceptance.
  Double_t runAvg = 0;
  Double_t max = 0;
  Int_t nInAvg = 0;
  Double_t avgSqr = 0;
  Int_t nBins = (dNdetadphi.GetNbinsX() * 6) / (fCumuDiff->GetNbinsX() * 5);
  Int_t nInBin = 0;
  Int_t nCurBin = 0, nPrevBin = 0;
  Int_t nCurRefBin = 0, nPrevRefBin = 0;
  Int_t nBadBins = 0;

  // Then we loop over the input and calculate sum cos(k*n*phi)
  // and fill it in the reference and differential histograms
  Double_t eta, phi, weight;
  Double_t dQnRe = 0, dQ2nRe = 0, dQnIm = 0, dQ2nIm = 0;
  for (Int_t etaBin = 1; etaBin <= dNdetadphi.GetNbinsX(); etaBin++) {
    eta = dNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    nCurBin = fCumuDiff->GetXaxis()->FindBin(eta);
    nCurRefBin = fCumuRef->GetXaxis()->FindBin(eta);
    // If we have moved to a new bin in the flow hist, and less than half the eta
    // region has been covered by it we cut it away.
    if (nPrevBin == 0) nPrevBin = nCurBin;
    if (nPrevRefBin == 0) nPrevRefBin = nCurRefBin;
    if (nCurBin != nPrevBin) {
      if (nInBin <= nBins/2) {
	for (Int_t qBin = 1; qBin <= fCumuDiff->GetNbinsY(); qBin++) {
	  Double_t removeContent = fCumuDiff->GetBinContent(nPrevBin, qBin);
	  Double_t removeEta = fCumuDiff->GetXaxis()->GetBinCenter(nPrevBin);
    	  if (nCurRefBin != nPrevRefBin) {
	    if (!(fFlags & kEtaGap)) fCumuRef->Fill(removeEta, qBin, -removeContent);
	    if ((fFlags & kSymEta)) fCumuRef->Fill(-1.*removeEta, qBin, -removeContent);
	  }
	  fCumuDiff->SetBinContent(nPrevBin, qBin, 0);
	  fCumuDiff->SetBinError(nPrevBin, qBin, 0);
	}
      }
      nInBin = 0;
      nPrevBin = nCurBin;
      if (nCurRefBin != nPrevRefBin)  nPrevRefBin = nCurRefBin;
    }
    
    Bool_t data = kFALSE;
    for (Int_t phiBin = 0; phiBin <= dNdetadphi.GetNbinsY(); phiBin++) {
      if (phiBin == 0) {
	if (dNdetadphi.GetBinContent(etaBin, phiBin) == 0) break;
	if ((fFlags & kEtaGap) && (mode & kFillRef) && TMath::Abs(eta) < fEtaGap) break;
	else data = kTRUE;
	continue;
      }
      if (!(fFlags & kSatVtx) && (nCurBin == 34 || nCurBin == 35) && (phiBin == 17 || phiBin == 18)) continue;
      if ( (fFlags & kSatVtx) && eta < 0 && (phiBin == 17 || phiBin == 18)) continue;
      phi = dNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      weight = dNdetadphi.GetBinContent(etaBin, phiBin);

      // We calculate the average Nch per. bin
      avgSqr += weight*weight;
      runAvg += weight;
      nInAvg++;
      if (weight == 0) continue;
      if (weight > max) max = weight;

      dQnRe = weight*TMath::Cos(fMoment*phi);
      dQnIm = weight*TMath::Sin(fMoment*phi);
      dQ2nRe = weight*TMath::Cos(2.*fMoment*phi);
      dQ2nIm = weight*TMath::Sin(2.*fMoment*phi);

      // if we do not have an eta gap we fill the ref sums in eta
      if (!(fFlags & kEtaGap)) {
	if ((mode & kFillRef)){
	  fCumuRef->Fill(eta, kHmult, weight);
	  fCumuRef->Fill(eta, kHQnRe, dQnRe);
	  fCumuRef->Fill(eta, kHQnIm, dQnIm);
	  
	  fCumuRef->Fill(eta, kHQ2nRe, dQ2nRe);
	  fCumuRef->Fill(eta, kHQ2nIm, dQ2nIm);
	}
      }
      // if we have an eta gap or we symmetrize around eta = 0
      // we fill in -eta
      if ((fFlags & kEtaGap) || (fFlags & kSymEta)){
	if ((mode & kFillRef)){
	  fCumuRef->Fill(-eta, kHmult, weight);
	  fCumuRef->Fill(-eta, kHQnRe, dQnRe);
	  fCumuRef->Fill(-eta, kHQnIm, dQnIm);
	  
	  fCumuRef->Fill(-eta, kHQ2nRe, dQ2nRe);
	  fCumuRef->Fill(-eta, kHQ2nIm, dQ2nIm);
	}
      }

      // If we fill diff flow, we always fill it in eta
      if ((mode & kFillDiff)) {
	fCumuDiff->Fill(eta, kHmult, weight);
	fCumuDiff->Fill(eta, kHQnRe, dQnRe);
	fCumuDiff->Fill(eta, kHQnIm, dQnIm);
	fCumuDiff->Fill(eta, kHQ2nRe, dQ2nRe);
	fCumuDiff->Fill(eta, kHQ2nIm, dQ2nIm);
      }

      // Fill acc. map
      fdNdedpAcc->Fill(eta, phi, weight);
    }
    if (data) {
      nInBin++;
    }
    // Outlier cut calculations
    if (nInAvg > 3 && !(fFlags & kEtaGap)) {
      runAvg /= nInAvg;
      avgSqr /= nInAvg;
      Double_t stdev = TMath::Sqrt(nInAvg/(nInAvg-1))*TMath::Sqrt(avgSqr - runAvg*runAvg);
      Double_t nSigma = (stdev == 0 ? 0 : (max-runAvg)/stdev);
      if (fSigmaCut > 0. && nSigma >= fSigmaCut && cent <= 60) nBadBins++;
      else nBadBins = 0;
      fOutliers->Fill(cent, nSigma);
      // We still finish the loop, for fOutliers to make sense, 
      // but we do no keep the event for analysis
      if (nBadBins > 3) useEvent = kFALSE;
    }
    runAvg = 0;
    avgSqr = 0;
    nInAvg = 0;
    max = 0;
  }

  return useEvent;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::CumulantsAccumulate(Double_t cent, Bool_t skipFourP) 
{
  // 
  // Calculate the Q cumulant of order fMoment
  //
  // Parameters:
  //  cent: Centrality of event
  //

  if (!fCumuRef) AliFatal("You have not called AddOutput() - Terminating!");

  // We create the objects needed for the analysis
  Double_t dQnRe = 0, dQ2nRe = 0, dQnIm = 0, dQ2nIm = 0, mult = 0;
  Double_t pnRe = 0, p2nRe = 0, qnRe = 0, q2nRe = 0, pnIm = 0, p2nIm = 0, qnIm = 0, q2nIm = 0;
  Double_t two = 0, four = 0, twoPrime = 0, fourPrime = 0;
  Double_t cosPhi1Phi2 = 0, cosPhi1Phi2Phi3m = 0;
  Double_t sinPhi1Phi2 = 0, sinPhi1Phi2Phi3m = 0;
  Double_t cosPsi1Phi2 = 0, cosPsi1Phi2Phi3m = 0, cosPsi1Phi2Phi3p = 0;
  Double_t sinPsi1Phi2 = 0, sinPsi1Phi2Phi3m = 0, sinPsi1Phi2Phi3p = 0;
  Double_t eta = 0;
  Double_t mp = 0, mq = 0;
  Double_t w2 = 0, w4 = 0, w2p = 0, w4p = 0;
  Int_t refEtaBin = 0;

  // We loop over the data 1 time!
  for (Int_t etaBin = 1; etaBin <= fCumuDiff->GetNbinsX(); etaBin++) {
    eta = fCumuDiff->GetXaxis()->GetBinCenter(etaBin);
    refEtaBin = fCumuRef->GetXaxis()->FindBin(eta);
    // The values for each individual etaBin bins are reset
    mult = 0;
    dQnRe = 0;
    dQnIm = 0;
    dQ2nRe = 0;
    dQ2nIm = 0;
    
    mp = 0;
    pnRe = 0;
    p2nRe = 0;
    pnIm = 0;
    p2nIm = 0;

    mq = 0;
    qnRe = 0;
    qnIm = 0;
    q2nRe = 0;
    q2nIm = 0;

    // Reference flow
    mult   = fCumuRef->GetBinContent(refEtaBin, kHmult);
    dQnRe  = fCumuRef->GetBinContent(refEtaBin, kHQnRe);
    dQnIm  = fCumuRef->GetBinContent(refEtaBin, kHQnIm);
    dQ2nRe = fCumuRef->GetBinContent(refEtaBin, kHQ2nRe);
    dQ2nIm = fCumuRef->GetBinContent(refEtaBin, kHQ2nIm);
    
    // For each etaBin bin the necessary values for differential flow
    // is calculated. .
    mp = fCumuDiff->GetBinContent(etaBin, kHmult);
    pnRe = fCumuDiff->GetBinContent(etaBin, kHQnRe);
    pnIm = fCumuDiff->GetBinContent(etaBin, kHQnIm);
    p2nRe = fCumuDiff->GetBinContent(etaBin, kHQ2nRe);
    p2nIm = fCumuDiff->GetBinContent(etaBin, kHQ2nIm);
    
    if (mult <= 3) continue; 

    if (mp == 0) continue; 
    // The reference flow is calculated 
    // 2-particle
    w2 = mult * (mult - 1.);
    two = dQnRe*dQnRe + dQnIm*dQnIm - mult;
    
    fCumuHist->Fill(eta, cent, kW2Two, two);
    fCumuHist->Fill(eta, cent, kW2, w2);

    fCumuHist->Fill(eta, cent, kQnRe, dQnRe);
    fCumuHist->Fill(eta, cent, kQnIm, dQnIm);
    fCumuHist->Fill(eta, cent, kM, mult);

    // Differential flow calculations for each eta bin bin is done:
    // 2-particle differential flow
    if (!(fFlags & kEtaGap)) {
      mq = mp;
      qnRe = pnRe;
      qnIm = pnIm;
      q2nRe = p2nRe;
      q2nIm = p2nIm;
    }

    w2p = mp * mult - mq;
    twoPrime = pnRe*dQnRe + pnIm*dQnIm - mq;
    
    fCumuHist->Fill(eta, cent, kw2two, twoPrime);
    fCumuHist->Fill(eta, cent, kw2, w2p);

    fCumuHist->Fill(eta, cent, kpnRe, pnRe);
    fCumuHist->Fill(eta, cent, kpnIm, pnIm);
    fCumuHist->Fill(eta, cent, kmp, mp);

    if ((fFlags & kEtaGap) || skipFourP) continue;
    // The reference flow is calculated 
    // 4-particle
    w4 = mult * (mult - 1.) * (mult - 2.) * (mult - 3.);
  
    four = 2.*mult*(mult-3.) + TMath::Power((TMath::Power(dQnRe,2.)+TMath::Power(dQnIm,2.)),2.)
             -4.*(mult-2.)*(TMath::Power(dQnRe,2.) + TMath::Power(dQnIm,2.))
             -2.*(TMath::Power(dQnRe,2.)*dQ2nRe+2.*dQnRe*dQnIm*dQ2nIm-TMath::Power(dQnIm,2.)*dQ2nRe)
             +(TMath::Power(dQ2nRe,2.)+TMath::Power(dQ2nIm,2.));

    fCumuHist->Fill(eta, cent, kW4Four, four);
    fCumuHist->Fill(eta, cent, kW4, w4);

    cosPhi1Phi2 = dQnRe*dQnRe - dQnIm*dQnIm - dQ2nRe;
    sinPhi1Phi2 = 2.*dQnRe*dQnIm - dQ2nIm;
      
    cosPhi1Phi2Phi3m = dQnRe*(TMath::Power(dQnRe,2)+TMath::Power(dQnIm,2))-dQnRe*dQ2nRe-dQnIm*dQ2nIm-2.*(mult-1)*dQnRe;

    sinPhi1Phi2Phi3m = -dQnIm*(TMath::Power(dQnRe,2)+TMath::Power(dQnIm,2))+dQnRe*dQ2nIm-dQnIm*dQ2nRe+2.*(mult-1)*dQnIm; 

    fCumuHist->Fill(eta, cent, kCosphi1phi2, cosPhi1Phi2);
    fCumuHist->Fill(eta, cent, kSinphi1phi2, sinPhi1Phi2);
    fCumuHist->Fill(eta, cent, kCosphi1phi2phi3m, cosPhi1Phi2Phi3m);
    fCumuHist->Fill(eta, cent, kSinphi1phi2phi3m, sinPhi1Phi2Phi3m);
    fCumuHist->Fill(eta, cent, kMm1m2, mult*(mult-1.)*(mult-2.));

    // Differential flow calculations for each eta bin bin is done:
    // 4-particle differential flow
    w4p = (mp * mult - 3.*mq)*(mult - 1.)*(mult - 2.);
 
    fourPrime = (TMath::Power(dQnRe,2.)+TMath::Power(dQnIm,2.))*(pnRe*dQnRe+pnIm*dQnIm)
                      - q2nRe*(TMath::Power(dQnRe,2.)-TMath::Power(dQnIm,2.))
                      - 2.*q2nIm*dQnRe*dQnIm
                      - pnRe*(dQnRe*dQ2nRe+dQnIm*dQ2nIm)
                      + pnIm*(dQnIm*dQ2nRe-dQnRe*dQ2nIm)
                      - 2.*mult*(pnRe*dQnRe+pnIm*dQnIm)
                      - 2.*(TMath::Power(dQnRe,2.)+TMath::Power(dQnIm,2.))*mq                      
                      + 6.*(qnRe*dQnRe+qnIm*dQnIm)                                            
                      + 1.*(q2nRe*dQ2nRe+q2nIm*dQ2nIm)                      
                      + 2.*(pnRe*dQnRe+pnIm*dQnIm)                       
                      + 2.*mq*mult                      
                      - 6.*mq; 

    fCumuHist->Fill(eta, cent, kw4four, fourPrime);
    fCumuHist->Fill(eta, cent, kw4, w4p);

    cosPsi1Phi2 = pnRe*dQnRe - pnIm*dQnIm - q2nRe;
    sinPsi1Phi2 = pnRe*dQnIm + pnIm*dQnRe - q2nIm;

    cosPsi1Phi2Phi3p = pnRe*(TMath::Power(dQnIm,2.)+TMath::Power(dQnRe,2.)-mult)
                          - 1.*(q2nRe*dQnRe+q2nIm*dQnIm)  
                          - mq*dQnRe+2.*qnRe;
 
    sinPsi1Phi2Phi3p = pnIm*(TMath::Power(dQnIm,2.)+TMath::Power(dQnRe,2.)-mult)
                          - 1.*(q2nIm*dQnRe-q2nRe*dQnIm)  
                          - mq*dQnIm+2.*qnIm; 

    cosPsi1Phi2Phi3m = pnRe*(TMath::Power(dQnRe,2.)-TMath::Power(dQnIm,2.))+2.*pnIm*dQnRe*dQnIm
                          - 1.*(pnRe*dQ2nRe+pnIm*dQ2nIm)  
                          - 2.*mq*dQnRe+2.*qnRe;
 
    sinPsi1Phi2Phi3m = pnIm*(TMath::Power(dQnRe,2.)-TMath::Power(dQnIm,2.))-2.*pnRe*dQnRe*dQnIm
                          - 1.*(pnIm*dQ2nRe-pnRe*dQ2nIm)
                          + 2.*mq*dQnIm-2.*qnIm;

    fCumuHist->Fill(eta, cent, kCospsi1phi2, cosPsi1Phi2);
    fCumuHist->Fill(eta, cent, kSinpsi1phi2, sinPsi1Phi2);
    fCumuHist->Fill(eta, cent, kCospsi1phi2phi3m, cosPsi1Phi2Phi3m);
    fCumuHist->Fill(eta, cent, kSinpsi1phi2phi3m, sinPsi1Phi2Phi3m);
    fCumuHist->Fill(eta, cent, kmpmq, (mp*mult-2.*mq)*(mult-1.));
    fCumuHist->Fill(eta, cent, kCospsi1phi2phi3p, cosPsi1Phi2Phi3p);
    fCumuHist->Fill(eta, cent, kSinpsi1phi2phi3p, sinPsi1Phi2Phi3p); 
  }
  // Event count
  fCumuHist->Fill(-7., cent, -0.5, 1.);

}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::CumulantsTerminate(TList* inlist, TList* outlist) 
{
  // 
  //  Finalizes the Q cumulant calculations
  // 
  //  Parameters:
  //   inlist: input sumlist
  //   outlist: output result list 
  //

  // Re-find cumulants hist if Terminate is called separately
  if (!fCumuHist) {
    TList* list = (TList*)inlist->FindObject(Form("%svertex_%d_%d%s", fType.Data(), fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")));
    fCumuHist = (TH3D*)list->FindObject(Form("%sv%d_vertex_%d_%d%s_cumu", fType.Data(), fMoment, 
                                        fVzMin, fVzMax, ((fFlags & kEtaGap) ? "_etaGap" : "")));
  }

  // Create result profiles
  TList* vtxList = (TList*)outlist->FindObject("vtxList");
  if (!vtxList) {
    vtxList = new TList();
    vtxList->SetName("vtxList");
    outlist->Add(vtxList);
  }
  TH1I* quality = (TH1I*)outlist->FindObject(Form("hQCquality%s%s", fType.Data(), ((fFlags & kEtaGap) ? "_etaGap" : "")));
  if (!quality) {
    quality = new TH1I(Form("hQCquality%s%s", fType.Data(), ((fFlags & kEtaGap) ? "_etaGap" : "")), 
                       Form("hQCquality%s%s", fType.Data(), ((fFlags & kEtaGap) ? "_etaGap" : "")),
                       40, 1, 41);
    quality->GetXaxis()->SetBinLabel( 1, "QC_{2}{2} > 0");
    quality->GetXaxis()->SetBinLabel( 2, "QC_{2}{2} <= 0");
    quality->GetXaxis()->SetBinLabel( 3, "QC'_{2}{2} > 0");
    quality->GetXaxis()->SetBinLabel( 4, "QC'_{2}{2} <= 0");
    quality->GetXaxis()->SetBinLabel( 5, "QC_{2}{4} < 0");
    quality->GetXaxis()->SetBinLabel( 6, "QC_{2}{4} >= 0");
    quality->GetXaxis()->SetBinLabel( 7, "QC'_{2}{4} < 0");
    quality->GetXaxis()->SetBinLabel( 8, "QC'_{2}{4} >= 0");
    quality->GetXaxis()->SetBinLabel( 9, "QC_{3}{2} > 0");
    quality->GetXaxis()->SetBinLabel(10, "QC_{3}{2} <= 0");
    quality->GetXaxis()->SetBinLabel(11, "QC'_{3}{2} > 0");
    quality->GetXaxis()->SetBinLabel(12, "QC'_{3}{2} <= 0");
    quality->GetXaxis()->SetBinLabel(13, "QC_{3}{4} < 0");
    quality->GetXaxis()->SetBinLabel(14, "QC_{3}{4} >= 0");
    quality->GetXaxis()->SetBinLabel(15, "QC'_{3}{4} < 0");
    quality->GetXaxis()->SetBinLabel(16, "QC'_{3}{4} >= 0");
    quality->GetXaxis()->SetBinLabel(17, "QC_{4}{2} > 0");
    quality->GetXaxis()->SetBinLabel(18, "QC_{4}{2} <= 0");
    quality->GetXaxis()->SetBinLabel(19, "QC'_{4}{2} > 0");
    quality->GetXaxis()->SetBinLabel(20, "QC'_{4}{2} <= 0");
    quality->GetXaxis()->SetBinLabel(21, "QC_{4}{4} < 0");
    quality->GetXaxis()->SetBinLabel(22, "QC_{4}{4} >= 0");
    quality->GetXaxis()->SetBinLabel(23, "QC'_{4}{4} < 0");
    quality->GetXaxis()->SetBinLabel(24, "QC'_{4}{4} >= 0");
    quality->GetXaxis()->SetBinLabel(25, "QC_{5}{2} > 0");
    quality->GetXaxis()->SetBinLabel(26, "QC_{5}{2} <= 0");
    quality->GetXaxis()->SetBinLabel(27, "QC'_{5}{2} > 0");
    quality->GetXaxis()->SetBinLabel(28, "QC'_{5}{2} <= 0");
    quality->GetXaxis()->SetBinLabel(29, "QC_{5}{4} < 0");
    quality->GetXaxis()->SetBinLabel(30, "QC_{5}{4} >= 0");
    quality->GetXaxis()->SetBinLabel(31, "QC'_{5}{4} < 0");
    quality->GetXaxis()->SetBinLabel(32, "QC'_{5}{4} >= 0");
    quality->GetXaxis()->SetBinLabel(33, "QC_{6}{2} > 0");
    quality->GetXaxis()->SetBinLabel(34, "QC_{6}{2} <= 0");
    quality->GetXaxis()->SetBinLabel(35, "QC'_{6}{2} > 0");
    quality->GetXaxis()->SetBinLabel(36, "QC'_{6}{2} <= 0");
    quality->GetXaxis()->SetBinLabel(37, "QC_{6}{4} < 0");
    quality->GetXaxis()->SetBinLabel(38, "QC_{6}{4} >= 0");
    quality->GetXaxis()->SetBinLabel(39, "QC'_{6}{4} < 0");
    quality->GetXaxis()->SetBinLabel(40, "QC'_{6}{4} >= 0");
    outlist->Add(quality);
  }

  TProfile2D* cumu2Sum = (TProfile2D*)outlist->FindObject(Form("%sQC2_v%d_unCorr%s", fType.Data(), fMoment, ((fFlags & kEtaGap) ? "_etaGap" : ""))); 
  TProfile2D* cumu4Sum = (TProfile2D*)outlist->FindObject(Form("%sQC4_v%d_unCorr%s", fType.Data(), fMoment, ((fFlags & kEtaGap) ? "_etaGap" : ""))); 
  if (!cumu2Sum) {
    cumu2Sum = new TProfile2D(Form("%sQC2_v%d_unCorr%s", fType.Data(), fMoment, ((fFlags & kEtaGap) ? "_etaGap" : "")),
                           Form("%sQC2_v%d_unCorr%s", fType.Data(), fMoment, ((fFlags & kEtaGap) ? "_etaGap" : "")),
	      fCumuHist->GetNbinsX(), fCumuHist->GetXaxis()->GetXmin(), fCumuHist->GetXaxis()->GetXmax(), 
	      fCumuHist->GetNbinsY(), fCumuHist->GetYaxis()->GetXmin(), fCumuHist->GetYaxis()->GetXmax());
    SetupCentAxis(cumu2Sum->GetYaxis());
    outlist->Add(cumu2Sum);
  }
  if (!cumu4Sum && !(fFlags & kEtaGap)) {
    cumu4Sum = new TProfile2D(Form("%sQC4_v%d_unCorr%s", fType.Data(), fMoment, ((fFlags & kEtaGap) ? "_etaGap" : "")),
                           Form("%sQC4_v%d_unCorr%s", fType.Data(), fMoment, ((fFlags & kEtaGap) ? "_etaGap" : "")),
	      fCumuHist->GetNbinsX(), fCumuHist->GetXaxis()->GetXmin(), fCumuHist->GetXaxis()->GetXmax(), 
	      fCumuHist->GetNbinsY(), fCumuHist->GetYaxis()->GetXmin(), fCumuHist->GetYaxis()->GetXmax());
    SetupCentAxis(cumu4Sum->GetYaxis());
    outlist->Add(cumu4Sum);
  }
 
  TProfile2D* cumu2 = new TProfile2D(Form("%sQC2_v%d_unCorr%s_vtx%3.1f", fType.Data(), fMoment, 
                         ((fFlags & kEtaGap) ? "_etaGap" : ""), (fVzMin+fVzMax)/2.),
			 Form("%sQC2_v%d_unCorr%s_vtx%3.1f", fType.Data(), fMoment,
			 ((fFlags & kEtaGap) ? "_etaGap" : ""), (fVzMin+fVzMax)/2.),
	    fCumuHist->GetNbinsX(), fCumuHist->GetXaxis()->GetXmin(), fCumuHist->GetXaxis()->GetXmax(), 
	    fCumuHist->GetNbinsY(), fCumuHist->GetYaxis()->GetXmin(), fCumuHist->GetYaxis()->GetXmax());
  SetupCentAxis(cumu2->GetYaxis());
  vtxList->Add(cumu2);
  TProfile2D* cumu4 = 0;
  if (!(fFlags & kEtaGap)) { 
    cumu4 = new TProfile2D(Form("%sQC4_v%d_unCorr%s_vtx%3.1f", fType.Data(), fMoment, 
      			    ((fFlags & kEtaGap) ? "_etaGap" : ""), (fVzMin+fVzMax)/2.),
     			    Form("%sQC4_v%d_unCorr%s_vtx%3.1f", fType.Data(), fMoment, 
   			    ((fFlags & kEtaGap) ? "_etaGap" : ""), (fVzMin+fVzMax)/2.),
  	      fCumuHist->GetNbinsX(), fCumuHist->GetXaxis()->GetXmin(), fCumuHist->GetXaxis()->GetXmax(), 
  	      fCumuHist->GetNbinsY(), fCumuHist->GetYaxis()->GetXmin(), fCumuHist->GetYaxis()->GetXmax());
    SetupCentAxis(cumu4->GetYaxis());
    vtxList->Add(cumu4);
  }

  // For flow calculations
  Double_t two = 0, qc2 = 0, vnTwo = 0, four = 0, qc4 = 0, vnFour = 0; 
  Double_t twoPrime = 0, qc2Prime = 0, vnTwoDiff = 0, fourPrime = 0, qc4Prime = 0, vnFourDiff = 0;
  Double_t w2 = 0, w4 = 0, w2p = 0, w4p = 0;
  Double_t w2Two = 0, w2pTwoPrime = 0, w4Four = 0, w4pFourPrime = 0;
  Double_t cosP1nPhi = 0, sinP1nPhi = 0, mult = 0, cosP1nPhi1P1nPhi2 = 0, sinP1nPhi1P1nPhi2 = 0;
  Double_t cosP1nPhi1M1nPhi2M1nPhi3 = 0, sinP1nPhi1M1nPhi2M1nPhi3 = 0, multm1m2 = 0;
  Double_t cosP1nPsi = 0, sinP1nPsi = 0, mp = 0, cosP1nPsi1P1nPhi2 = 0, sinP1nPsi1P1nPhi2 = 0;
  Double_t cosP1nPsi1M1nPhi2M1nPhi3 = 0, sinP1nPsi1M1nPhi2M1nPhi3 = 0, mpqMult = 0;
  Double_t cosP1nPsi1P1nPhi2M1nPhi3 = 0, sinP1nPsi1P1nPhi2M1nPhi3 = 0;

  // Loop over cumulant histogram for final calculations   
  // Centrality loop
  for (Int_t cBin = 1; cBin <= fCumuHist->GetNbinsY(); cBin++) {
    // for weighted avg.
    Double_t cent = fCumuHist->GetYaxis()->GetBinCenter(cBin);
    if (fDebug > 0) AliInfo(Form("%s - v_%d: centrality %3.1f:..", fType.Data(), fMoment, cent));
    // Eta loop
    for (Int_t etaBin = 1; etaBin <= fCumuHist->GetNbinsX(); etaBin++) {
      Double_t eta = fCumuHist->GetXaxis()->GetBinCenter(etaBin);
      // 2-particle reference flow
      w2Two = fCumuHist->GetBinContent(etaBin, cBin, kW2Two);
      w2 = fCumuHist->GetBinContent(etaBin, cBin, kW2);
      mult = fCumuHist->GetBinContent(etaBin, cBin, kM);
      if (!w2 || !mult) continue;
      cosP1nPhi = fCumuHist->GetBinContent(etaBin, cBin, kQnRe);
      sinP1nPhi = fCumuHist->GetBinContent(etaBin, cBin, kQnIm);
        
      cosP1nPhi /= mult;
      sinP1nPhi /= mult;
      two = w2Two / w2;
      qc2 = two - TMath::Power(cosP1nPhi, 2) - TMath::Power(sinP1nPhi, 2);
      if (qc2 <= 0) { 
	if (fDebug > 0) AliInfo(Form("%s: QC_%d{2} = %1.3f for eta = %1.2f and centrality %3.1f - skipping", fType.Data(), fMoment, qc2, eta, cent));
	quality->Fill((fMoment-2)*8+2);
	continue;
      }
       vnTwo = TMath::Sqrt(qc2);
      if (!TMath::IsNaN(vnTwo*mult)) 
	quality->Fill((fMoment-2)*8+1);
 //       cumu2->Fill(eta, cent, vnTwo, fCumuHist->GetBinContent(0,cBin,0)); 

      // 2-particle differential flow
      w2pTwoPrime = fCumuHist->GetBinContent(etaBin, cBin, kw2two);
      w2p = fCumuHist->GetBinContent(etaBin, cBin, kw2);
      mp = fCumuHist->GetBinContent(etaBin, cBin, kmp);
      if (!w2p || !mp) continue;
      cosP1nPsi = fCumuHist->GetBinContent(etaBin, cBin, kpnRe);
      sinP1nPsi = fCumuHist->GetBinContent(etaBin, cBin, kpnIm);

      cosP1nPsi /= mp;
      sinP1nPsi /= mp;
      twoPrime = w2pTwoPrime / w2p;
      qc2Prime = twoPrime - sinP1nPsi*sinP1nPhi - cosP1nPsi*cosP1nPhi;

      vnTwoDiff = qc2Prime / TMath::Sqrt(qc2);
      if (!TMath::IsNaN(vnTwoDiff*mp)) {
	cumu2->Fill(eta, cent, vnTwoDiff);
	quality->Fill((fMoment-2)*8+3);
      }
      else 
	quality->Fill((fMoment-2)*8+4);

      if (fDebug > 1) AliInfo(Form("%s: v_%d{2} = %1.3f for eta = %1.2f and centrality %3.1f", fType.Data(), fMoment, vnTwoDiff, eta, cent));

      if ((fFlags & kEtaGap)) continue;
      // 4-particle reference flow
      w4Four = fCumuHist->GetBinContent(etaBin, cBin, kW4Four);
      w4 = fCumuHist->GetBinContent(etaBin, cBin, kW4);
      multm1m2 = fCumuHist->GetBinContent(etaBin, cBin, kMm1m2);
      if (!w4 || !multm1m2) continue;
      cosP1nPhi1P1nPhi2 = fCumuHist->GetBinContent(etaBin, cBin, kCosphi1phi2);
      sinP1nPhi1P1nPhi2 = fCumuHist->GetBinContent(etaBin, cBin, kSinphi1phi2);
      cosP1nPhi1M1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, cBin, kCosphi1phi2phi3m);
      sinP1nPhi1M1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, cBin, kSinphi1phi2phi3m);

      cosP1nPhi1P1nPhi2 /= w2;
      sinP1nPhi1P1nPhi2 /= w2;
      cosP1nPhi1M1nPhi2M1nPhi3 /= multm1m2;
      sinP1nPhi1M1nPhi2M1nPhi3 /= multm1m2;
      four = w4Four / w4;
      qc4 = four-2.*TMath::Power(two,2.)
         - 4.*cosP1nPhi*cosP1nPhi1M1nPhi2M1nPhi3
         + 4.*sinP1nPhi*sinP1nPhi1M1nPhi2M1nPhi3-TMath::Power(cosP1nPhi1P1nPhi2,2.)-TMath::Power(sinP1nPhi1P1nPhi2,2.)
         + 4.*cosP1nPhi1P1nPhi2*(TMath::Power(cosP1nPhi,2.)-TMath::Power(sinP1nPhi,2.))
         + 8.*sinP1nPhi1P1nPhi2*sinP1nPhi*cosP1nPhi
         + 8.*two*(TMath::Power(cosP1nPhi,2.)+TMath::Power(sinP1nPhi,2.))
         - 6.*TMath::Power((TMath::Power(cosP1nPhi,2.)+TMath::Power(sinP1nPhi,2.)),2.);
      
      if (qc4 >= 0) {
	if (fDebug > 0) AliInfo(Form("%s: QC_%d{4} = %1.3f for eta = %1.2f and centrality %3.1f - skipping", fType.Data(), fMoment, qc2, eta, cent));
	quality->Fill((fMoment-2)*8+6);
   	continue;
      }
      vnFour = TMath::Power(-qc4, 0.25);
      if (!TMath::IsNaN(vnFour*mult)) 
	quality->Fill((fMoment-2)*8+5);
 //         cumu4->Fill(eta, cent, vnFour, fCumuHist->GetBinContent(0,cBin,0));

      // 4-particle differential flow
      w4pFourPrime = fCumuHist->GetBinContent(etaBin, cBin, kw4four);
      w4p = fCumuHist->GetBinContent(etaBin, cBin, kw4);
      mpqMult = fCumuHist->GetBinContent(etaBin, cBin, kmpmq);
      if (!w4p || !mpqMult) continue;
      cosP1nPsi1P1nPhi2 = fCumuHist->GetBinContent(etaBin, cBin, kCospsi1phi2);
      sinP1nPsi1P1nPhi2 = fCumuHist->GetBinContent(etaBin, cBin, kSinpsi1phi2);
      cosP1nPsi1M1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, cBin, kCospsi1phi2phi3m);
      sinP1nPsi1M1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, cBin, kSinpsi1phi2phi3m);
      cosP1nPsi1P1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, cBin, kCospsi1phi2phi3p);
      sinP1nPsi1P1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, cBin, kSinpsi1phi2phi3p); 
      
      cosP1nPsi1P1nPhi2 /= w2p;
      sinP1nPsi1P1nPhi2 /= w2p;
      cosP1nPsi1M1nPhi2M1nPhi3 /= mpqMult;
      sinP1nPsi1M1nPhi2M1nPhi3 /= mpqMult;
      cosP1nPsi1P1nPhi2M1nPhi3 /= mpqMult;
      sinP1nPsi1P1nPhi2M1nPhi3 /= mpqMult;

      fourPrime = w4pFourPrime / w4p;

      qc4Prime = fourPrime-2.*twoPrime*two
                - cosP1nPsi*cosP1nPhi1M1nPhi2M1nPhi3
                + sinP1nPsi*sinP1nPhi1M1nPhi2M1nPhi3
                - cosP1nPhi*cosP1nPsi1M1nPhi2M1nPhi3
                + sinP1nPhi*sinP1nPsi1M1nPhi2M1nPhi3
                - 2.*cosP1nPhi*cosP1nPsi1P1nPhi2M1nPhi3
                - 2.*sinP1nPhi*sinP1nPsi1P1nPhi2M1nPhi3
                - cosP1nPsi1P1nPhi2*cosP1nPhi1P1nPhi2
                - sinP1nPsi1P1nPhi2*sinP1nPhi1P1nPhi2
                + 2.*cosP1nPhi1P1nPhi2*(cosP1nPsi*cosP1nPhi-sinP1nPsi*sinP1nPhi)
                + 2.*sinP1nPhi1P1nPhi2*(cosP1nPsi*sinP1nPhi+sinP1nPsi*cosP1nPhi)
                + 4.*two*(cosP1nPsi*cosP1nPhi+sinP1nPsi*sinP1nPhi)
                + 2.*cosP1nPsi1P1nPhi2*(TMath::Power(cosP1nPhi,2.)-TMath::Power(sinP1nPhi,2.))
                + 4.*sinP1nPsi1P1nPhi2*cosP1nPhi*sinP1nPhi
                + 4.*twoPrime*(TMath::Power(cosP1nPhi,2.)+TMath::Power(sinP1nPhi,2.))
                - 6.*(TMath::Power(cosP1nPhi,2.)-TMath::Power(sinP1nPhi,2.)) 
                * (cosP1nPsi*cosP1nPhi-sinP1nPsi*sinP1nPhi)
                - 12.*cosP1nPhi*sinP1nPhi
                * (sinP1nPsi*cosP1nPhi+cosP1nPsi*sinP1nPhi);

      vnFourDiff = - qc4Prime / TMath::Power(-qc4, 0.75);
      if (!TMath::IsNaN(vnFourDiff*mp) && vnFourDiff > 0) {
	cumu4->Fill(eta, cent, vnFourDiff);
	quality->Fill((fMoment-2)*8+7);
      }
      else 
	quality->Fill((fMoment-2)*8+8);

      if (fDebug > 1) AliInfo(Form("%s: v_%d{4} = %1.3f for eta = %1.2f and centrality %3.1f", fType.Data(), fMoment, vnFourDiff, eta, cent));
    } // End of eta loop
  } // End of centrality loop
 
  cumu2Sum->Add(cumu2);
  if (!(fFlags & kEtaGap)) cumu4Sum->Add(cumu4);

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::SetupCentAxis(TAxis* axis)  
{
  // 
  // Setup centrality axis for histogram
  //
  // Parameters:
  //  axis: centrality axis
  //
  if (!axis) {
    AliError("Null pointer passed for axis");
    return;
  }

  if ((fFlags & kSatVtx)) {
    Double_t cent[3] = {0, 40, 100};
    axis->Set(2, cent);
  } else {
    Double_t cent[13] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100 };
    axis->Set(12, cent);
  }

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::PrintFlowSetup() const 
{
  //
  // Print the setup of the flow task
  //
  Printf("AliForwardFlowTaskQC::Print");
  Printf("Number of bins in vertex axis   :\t%d", fVtxAxis->GetNbins());
  Printf("Range of vertex axis            :\t[%3.1f,%3.1f]", 
			  fVtxAxis->GetXmin(), fVtxAxis->GetXmax());
  printf("Doing flow analysis for         :\t");
  for (Int_t n  = 0; n < fV.GetSize(); n++) printf("v%d ", fV.At(n));
  printf("\n");
  Printf("Satellite vertex flag           :\t%s", ((fFlowFlags & kSatVtx) ? "true" : "false"));
  Printf("Symmetrize ref. flow wrt eta = 0:\t%s", ((fFlowFlags & kSymEta) ? "true" : "false"));
  Printf("Use an eta-gap for ref. flow    :\t%s", ((fFlowFlags & kEtaGap) ? "true" : "false"));
  Printf("FMD sigma cut:                  :\t%f", fFMDCut);
  Printf("SPD sigma cut:                  :\t%f", fSPDCut);
  if ((fFlowFlags & kEtaGap)) 
    Printf("Eta gap:                        :\t%f", fEtaGap);

}
//_____________________________________________________________________
//
//
// EOF
