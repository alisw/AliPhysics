
//
// Calculate flow in the forward and central regions using the Q cumulants method.
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root or forward_flow.root
//
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TParameter.h>
#include <TMatrixD.h>
#include <TVectorD.h>
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
#include "AliVVZERO.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAnalysisFilter.h"

ClassImp(AliForwardFlowTaskQC)
#if 0
; // For emacs 
#endif

AliForwardFlowTaskQC::AliForwardFlowTaskQC()
  : AliAnalysisTaskSE(),
    fVtxAxis(),          // Axis to control vertex binning
    fCentAxis(),         // Axis to control centrality/multiplicity binning
    fFMDCut(-1),         // FMD sigma cut
    fSPDCut(-1),         // SPD sigma cut
    fFlowFlags(0),       // Flow flags
    fEtaGap(0.),         // Eta gap value
    fBinsForward(),      // List with forward flow hists
    fBinsCentral(),      // List with central flow hists
    fSumList(0),	 // Event sum list
    fOutputList(0),	 // Result output list
    fAOD(0),		 // AOD input event
    fTrackCuts(0),    // ESD track cuts
    fMaxMoment(0),       // Max flow moment
    fVtx(1111),	         // Z vertex coordinate
    fCent(-1),		 // Centrality
    fHistdNdedpV0(),     // Hist for v0
    fHistdNdedp3Cor(),   // Hist for combining detectors
    fHistFMDSPDCorr(),   // FMD SPD correlation
    fHistCent(),         // Hist for centrality
    fHistVertexSel(),    // Hist for selected vertices
    fHistEventSel()      // Hist for event selection
{
  // 
  //  Default constructor
  //
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const char* name) 
  : AliAnalysisTaskSE(name),
    fVtxAxis(),         // Axis to control vertex binning
    fCentAxis(),        // Axis to control centrality/multiplicity binning
    fFMDCut(-1),        // FMD sigma cut
    fSPDCut(-1),        // SPD sigma cut
    fFlowFlags(0x0),    // Flow flags
    fEtaGap(0.),        // Eta gap value
    fBinsForward(),     // List with forward flow hists
    fBinsCentral(),     // List with central flow hists
    fSumList(0),        // Event sum list           
    fOutputList(0),     // Result output list       
    fAOD(0),	        // AOD input event          
    fTrackCuts(0),      // ESD track cuts
    fMaxMoment(4),      // Max flow moment
    fVtx(1111),         // Z vertex coordinate      
    fCent(-1),          // Centrality               
    fHistdNdedpV0(),    // Histo for v0
    fHistdNdedp3Cor(),  // Histo for combining detectors
    fHistFMDSPDCorr(),  // FMD SPD correlation
    fHistCent(),        // Hist for centrality
    fHistVertexSel(),   // Hist for selected vertices
    fHistEventSel()     // Hist for event selection
{
  // 
  //  Constructor
  //
  //  Parameters:
  //   name: Name of task
  //
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const AliForwardFlowTaskQC& o)
  : AliAnalysisTaskSE(o),
    fVtxAxis(o.fVtxAxis),              // Axis to control vertex binning
    fCentAxis(o.fCentAxis),            // Array to control centrality/multiplicity binning
    fFMDCut(o.fFMDCut),                // FMD sigma cut
    fSPDCut(o.fSPDCut),                // SPD sigma cut
    fFlowFlags(o.fFlowFlags),          // Flow flags
    fEtaGap(o.fEtaGap),                // Eta gap value
    fBinsForward(),                    // List with forward flow hists
    fBinsCentral(),                    // List with central flow hists
    fSumList(o.fSumList),              // Event sum list           
    fOutputList(o.fOutputList),        // Result output list       
    fAOD(o.fAOD),	               // AOD input event          
    fTrackCuts(o.fTrackCuts),          // ESD track cuts
    fMaxMoment(o.fMaxMoment),          // Flow moments
    fVtx(o.fVtx),                      // Z vertex coordinate      
    fCent(o.fCent),	               // Centrality
    fHistdNdedpV0(o.fHistdNdedpV0),    // Histo for v0
    fHistdNdedp3Cor(o.fHistdNdedp3Cor),// Histo for combining detectors
    fHistFMDSPDCorr(o.fHistFMDSPDCorr),// FMD SPD correlation
    fHistCent(o.fHistCent),            // Hist for centrality
    fHistVertexSel(o.fHistVertexSel),  // Hist for selected vertices
    fHistEventSel(o.fHistEventSel)     // Hist for event selection
{
  // 
  //  Copy constructor 
  // 
  //  Parameters:
  //   o: Object to copy from 
  //
}
//_____________________________________________________________________
AliForwardFlowTaskQC&
AliForwardFlowTaskQC::operator=(const AliForwardFlowTaskQC& o)
{
  // 
  //  Assignment operator 
  //
  if (&o == this) return *this;
  fVtxAxis        = o.fVtxAxis;
  fCentAxis       = o.fCentAxis;
  fFMDCut         = o.fFMDCut;
  fSPDCut         = o.fSPDCut;
  fFlowFlags      = o.fFlowFlags;
  fEtaGap         = o.fEtaGap;
  fSumList        = o.fSumList;
  fOutputList     = o.fOutputList;
  fAOD            = o.fAOD;
  fTrackCuts      = o.fTrackCuts;
  fMaxMoment      = o.fMaxMoment;
  fVtx            = o.fVtx;
  fCent           = o.fCent;
  fHistdNdedpV0   = o.fHistdNdedpV0;
  fHistdNdedp3Cor = o.fHistdNdedp3Cor;
  fHistFMDSPDCorr = o.fHistFMDSPDCorr;
  fHistCent       = o.fHistCent;
  fHistVertexSel  = o.fHistVertexSel;
  fHistEventSel   = o.fHistEventSel;

  return *this;
}
//____________________________________________________________________
namespace {
  Double_t GetEdge(const TString& str, Int_t start, Int_t end)
  {
    TString sub(str(start, end));
    return sub.Atof();
  }
  Bool_t ExtractBins(const TString& spec, TArrayD& edges)
  {
    TArrayD tmp(200);
    Int_t   start = 0;
    Int_t   cnt   = 0;
    for (Int_t i=1; i<spec.Length(); i++) {
      if (spec[i] == '-' || spec[i] == ':') {
  Double_t c = GetEdge(spec, start, i);
  if (cnt > 0 && c < tmp[cnt-1]) {
    Warning("ExtractBins",
      "Invalid edge @ %d: %f (< %f)", cnt, c, tmp[cnt-1]);
  tmp.Set(0);
  return false;
  }
  tmp[cnt] = c;
  i++;
  start = i;
  cnt++;
      }
    }
    if (start+1 != spec.Length()) {
      Double_t c = GetEdge(spec, start, spec.Length());
      tmp[cnt] = c;
      cnt++;
    }
    edges.Set(cnt, tmp.GetArray());
    return true;
  }
}

//________________________________________________________________________
void
AliForwardFlowTaskQC::SetCentralityAxis(const char* bins)
{
  if (!bins || bins[0] == '\0') return;

  TString     spec(bins);
  if (spec.EqualTo("none", TString::kIgnoreCase))
    return;

  TArrayD edges;
  if (spec.EqualTo("default", TString::kIgnoreCase) ||
      spec.EqualTo("pbpb", TString::kIgnoreCase)) {
    //                 1  2  3   4   5   6   7   8   9   10  11
    Double_t tmp[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100 };
    edges.Set(11, tmp);
  }
  else if (spec.EqualTo("ppb", TString::kIgnoreCase) ||
     spec.EqualTo("pbp", TString::kIgnoreCase)) {
    //                 1  2  3   4   5   6   7   8
    Double_t tmp[] = { 0, 5, 10, 20, 40, 60, 80, 100 };
    edges.Set(8, tmp);
  }
  else {
    ExtractBins(spec, edges);
  }
  TAxis* centAxis = new TAxis(edges.GetSize()-1, edges.GetArray());
  SetCentralityAxis(centAxis);
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::SetFlowFlags(UShort_t flags)
{
  //
  //  Set flow flags, making sure the detector setup is right
  //
  //  Parameters:
  //   flags: Flow flags
  //
  if ((flags & kFMD) && (flags & kVZERO)) 
    AliFatal("Cannot do analysis on more than one forward detector!");
  else if (!(flags & kFMD) && !(flags & kVZERO)) 
    AliFatal("You need to add a forward detector!");
  else fFlowFlags = flags;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::UserCreateOutputObjects()
{
  //
  //  Create output objects
  //
  InitVertexBins();
  InitHists();
  if ((fFlowFlags & kTracks) && !fTrackCuts) AliFatal("No track cuts set!");
  PrintFlowSetup();

  PostData(1, fSumList);
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::InitVertexBins()
{
  // 
  //  Init vertexbin objects for forward and central detectors, and add them to the lists
  //
  for (Int_t v = 1; v <= fVtxAxis->GetNbins(); v++) {
    Int_t vL = Int_t(fVtxAxis->GetBinLowEdge(v));
    Int_t vH = Int_t(fVtxAxis->GetBinUpEdge(v));
    if ((fFlowFlags & kFMD)) {
      fBinsForward.Add(new VertexBin(vL, vH, fMaxMoment, "FMD", fFlowFlags, fFMDCut, fEtaGap));
      if (!(fFlowFlags & k3Cor)) 
	fBinsCentral.Add(new VertexBin(vL, vH, fMaxMoment, "SPD-FMD", fFlowFlags|kNUAcorr|kSPD, fSPDCut, fEtaGap));
    }
    else if ((fFlowFlags & kVZERO)) {
      fBinsForward.Add(new VertexBin(vL, vH, fMaxMoment, "VZERO", fFlowFlags, 0, fEtaGap));
      if ((fFlowFlags & kEtaGap) && !(fFlowFlags & kTracks)) 
	fBinsCentral.Add(new VertexBin(vL, vH, fMaxMoment, "SPD-VZERO", fFlowFlags|kNUAcorr|kSPD, fSPDCut, fEtaGap));
    }
  }
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::InitHists()
{
  //
  //  Init histograms and add vertex bin histograms to the sum list
  //
  if (!fSumList)
    fSumList = new TList();
  fSumList->SetName("Sums");
  fSumList->SetOwner();

  if (!fVtxAxis) fVtxAxis = new TAxis(20, -10, 10);
  fVtxAxis->SetName("VtxAxis");
  if (!fCentAxis) fCentAxis = new TAxis(20, 0, 100);
  fVtxAxis->SetName("CentAxis");
  
  fHistCent         = new TH1D("hCent", "Centralities", 100, 0, 100);
  fHistVertexSel    = new TH1D("hVertexSel", "Selected vertices", fVtxAxis->GetNbins(), fVtxAxis->GetXmin(), fVtxAxis->GetXmax());
  fHistEventSel     = new TH1I("hEventSel", "Event Selection", kOK, 0.5, kOK+0.5);
  fHistEventSel->GetXaxis()->SetBinLabel(kNoEvent, "No AOD event");
  fHistEventSel->GetXaxis()->SetBinLabel(kNoForward, "No forward det");
  fHistEventSel->GetXaxis()->SetBinLabel(kNoCentral, "No central det");
  fHistEventSel->GetXaxis()->SetBinLabel(kNoTrigger, "Not triggered");
  fHistEventSel->GetXaxis()->SetBinLabel(kNoCent, "No centrality");
  fHistEventSel->GetXaxis()->SetBinLabel(kInvCent, "Centrality outside range");
  fHistEventSel->GetXaxis()->SetBinLabel(kNoVtx, "No vertex");
  fHistEventSel->GetXaxis()->SetBinLabel(kInvVtx, "Vtx outside range");
  fHistEventSel->GetXaxis()->SetBinLabel(kOK, "OK!");

  fHistFMDSPDCorr = new TH2D("hFMDSPDCorr", "hFMDSPCCorr", 200, 0., 20000., 200, 0, 7500);

  TList* dList = new TList();
  dList->SetName("Diagnostics");
  dList->Add(fHistCent);
  dList->Add(fHistVertexSel);
  dList->Add(fHistEventSel);
  dList->Add(fHistFMDSPDCorr);
  fSumList->Add(dList);

  fHistdNdedp3Cor = TH2D(Form("hdNdedpCombined_%s", GetQCType(fFlowFlags)), Form("hdNdedpCombined_%s", GetQCType(fFlowFlags)), 
                           200, -4., 6., 20, 0., TMath::TwoPi());
  if ((fFlowFlags & kVZERO)) {
    Double_t bins[12] = { -6, -3.7, -3.2, -2.7, -2.2, -1.7, 
			  2.8, 3.4,  3.9,  4.5,  5.1, 6 };
    fHistdNdedpV0 = TH2D(Form("hdNdedpv0%s", GetQCType(fFlowFlags)), Form("hdNdedpv0%s", GetQCType(fFlowFlags)),
		    11, -6, 6, 8, 0, TMath::TwoPi());
    fHistdNdedpV0.GetXaxis()->Set(11, bins);
    if ((fFlowFlags & k3Cor)) {
      Double_t bins2[20] = { -6, -3.7, -3.2, -2.7, -2.2, // VZERO
                             -2.0, -1.5, -1.0, -0.5 , 0., 0.5, 1.0, 1.5, 2.0, // SPD
	  		    2.8, 3.4,  3.9,  4.5,  5.1, 6 }; // VZERO
      fHistdNdedp3Cor.GetXaxis()->Set(19, bins2);
      fHistdNdedp3Cor.GetYaxis()->Set(8, 0., TMath::TwoPi());
    }
  }

  TIter nextForward(&fBinsForward);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(nextForward()))) {
    bin->AddOutput(fSumList, fCentAxis);
  }
  TIter nextCentral(&fBinsCentral);
  while ((bin = static_cast<VertexBin*>(nextCentral()))) {
    bin->AddOutput(fSumList, fCentAxis);
  }
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::UserExec(Option_t */*option*/)
{
  //
  //  Calls the analyze function - called every event
  //
  //  Parameters:
  //   option: Not used
  //
  
  // Reset data members
  fCent = -1;
  fVtx = 1111;

  Analyze();
  
  PostData(1, fSumList);

  return;
}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::Analyze()
{
  // 
  //  Load forward and central detector objects from aod tree and call the 
  //  cumulants calculation for the correct vertex bin
  //
  //  Return: true on success
  //

  // Get input event
  fAOD = dynamic_cast<AliAODEvent*>(AliForwardUtil::GetAODEvent(this));
  if (!fAOD) {
    fHistEventSel->Fill(kNoEvent);
    return kFALSE;
  }

  // Get detector objects
  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
  AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClusters"));
  AliVVZERO* vzero = GetVZERO();
  if ((fFlowFlags & kVZERO)) {
    if (vzero) {
      fHistdNdedpV0.Reset();
      FillVZEROHist(vzero);
    }
  }

  // We make sure that the necessary forward object is there
  if ((fFlowFlags & kFMD) && !aodfmult) {
    fHistEventSel->Fill(kNoForward);
    return kFALSE; 
  }
  else if ((fFlowFlags & kVZERO) && !vzero) {
    fHistEventSel->Fill(kNoForward);
    return kFALSE; 
  }
  if (!aodcmult) fHistEventSel->Fill(kNoCentral);

   // Check event for triggers, get centrality, vtx etc.
  if (!CheckEvent(aodfmult)) return kFALSE;
  Int_t vtx = fVtxAxis->FindBin(fVtx)-1;
  
  // Then we assign a reference to the forward histogram of interest
  TH2D& forwarddNdedp = ((fFlowFlags & kFMD) ? aodfmult->GetHistogram() : fHistdNdedpV0);
  TH2D& spddNdedp = aodcmult->GetHistogram();
  if ((fFlowFlags & kStdQC)) {
    FillVtxBinList(fBinsForward, forwarddNdedp, vtx);
  } else if ((fFlowFlags & kEtaGap)) {
    FillVtxBinListEtaGap(fBinsForward, forwarddNdedp, forwarddNdedp, vtx);
  }
  // At the moment only clusters are supported for the central region (some day add tracks?)
  // So no extra checks necessary
  if (aodcmult) {
    if ((fFlowFlags & kStdQC)) {
      FillVtxBinList(fBinsCentral, spddNdedp, vtx);
    } else if ((fFlowFlags & kEtaGap)) {
      FillVtxBinListEtaGap(fBinsCentral, forwarddNdedp, spddNdedp, vtx);
    } else if ((fFlowFlags & k3Cor)) {
      FillVtxBinList3Cor(fBinsForward, spddNdedp, forwarddNdedp, vtx);
    }
    // Diagnostics
    if (aodfmult) {
      Double_t totForward = forwarddNdedp.Integral(1, forwarddNdedp.GetNbinsX(), 1, forwarddNdedp.GetNbinsY());
      Double_t totSPD = spddNdedp.Integral(1, spddNdedp.GetNbinsX(), 1, spddNdedp.GetNbinsY());
      fHistFMDSPDCorr->Fill(totForward, totSPD);
    }
  }

  return kTRUE;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::FillVtxBinList(const TList& list, TH2D& h, Int_t vtx, UShort_t flags) const
{
  //
  //  Loops over list of VtxBins, fills hists of bins for current vertex
  //  and runs analysis on those bins
  //
  //  Parameters:
  //   list:  list of VtxBins
  //   h:     dN/detadphi histogram
  //   vtx:   current vertex bin
  //   flags: extra flags to handle calculations.
  // 
  //   Note: The while loop is used in this function and the next 2 for historical reasons,
  //         as originally each moment had it's own VertexBin object.
  VertexBin* bin = 0;
  Int_t i = 0;
  Int_t nVtxBins = fVtxAxis->GetNbins();
  
  while ((bin = static_cast<VertexBin*>(list.At(vtx+(nVtxBins*i))))) {
    i++;
    // If no tracks do things normally
    if (!(fFlowFlags & kTracks) || (flags & kMC)) {
      if (!bin->FillHists(h, fCent, kFillBoth|flags|kReset)) continue;
    }
    // if tracks things are more complicated
    else if ((fFlowFlags & kTracks)) {
      if (!FillTracks(bin, kFillRef|kReset|flags)) continue;
      if (!bin->FillHists(h, fCent, kFillDiff|kReset|flags)) continue;
    }
    bin->CumulantsAccumulate(fCent);
  }

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::FillVtxBinListEtaGap(const TList& list, TH2D& href, 
                                                  TH2D& hdiff, Int_t vtx, UShort_t flags) const
{
  //
  //  Loops over list of VtxBins, fills hists of bins for current vertex
  //  and runs analysis on those bins
  //
  //  Parameters:
  //   list:  list of VtxBins
  //   href:  dN/detadphi histogram for ref. flow
  //   hdiff: dN/detadphi histogram for diff. flow
  //   vBin:  current vertex bin
  //   flags: extra flags to handle calculations.
  //
  VertexBin* bin = 0;
  Int_t i = 0;
  Int_t nVtxBins = fVtxAxis->GetNbins();

  while ((bin = static_cast<VertexBin*>(list.At(vtx+(nVtxBins*i))))) {
    i++;
    if (!(fFlowFlags & kTracks) || (flags & kMC)) {
      if(!bin->FillHists(href, fCent, kFillRef|flags|kReset)) continue;
    }
    else if ((fFlowFlags & kTracks)) {
      if (!FillTracks(bin, kFillRef|kReset|flags)) continue;
    }
    if (!bin->FillHists(hdiff, fCent, kFillDiff|kReset|flags)) continue;
    bin->CumulantsAccumulate(fCent);
  }

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::FillVtxBinList3Cor(const TList& list, TH2D& hcent, 
                                              TH2D& hfwd, Int_t vtx, UShort_t flags)
{
  //
  //  Loops over list of VtxBins, fills hists of bins for current vertex
  //  and runs analysis on those bins
  //
  //  Parameters:
  //   list:  list of VtxBins
  //   hcent: dN/detadphi histogram for central coverage
  //   hfwd:  dN/detadphi histogram for forward coverage
  //   vBin:  current vertex bin
  //   flags: extra flags to handle calculations.
  //
  VertexBin* bin = 0;
  Int_t i = 0;
  Int_t nVtxBins = fVtxAxis->GetNbins();

  TH2D& h = CombineHists(hcent, hfwd);

  while ((bin = static_cast<VertexBin*>(list.At(vtx+(nVtxBins*i))))) {
    i++;
    if (!bin->FillHists(h, fCent, kFillBoth|flags|kReset)) continue;
    bin->CumulantsAccumulate3Cor(fCent);
  }

  return;
}
//_____________________________________________________________________
TH2D& AliForwardFlowTaskQC::CombineHists(TH2D& hcent, TH2D& hfwd)
{
  // 
  //  Combines a forward and central d^2N/detadphi histogram.
  //  At some point it might need a flag to choose which histogram gets
  //  priority when there is an overlap, at the moment the average is chosen
  //
  //  Parameters:
  //    hcent: Central barrel detector
  //    hfwd: Forward detector
  //
  //  Return: reference to combined hist
  //

  // If hists are the same (MC input) don't do anything
  if (&hcent == &hfwd) return hcent;

  fHistdNdedp3Cor.Reset();
  // FMD, SPD input
  if ((fFlowFlags & kFMD)) {
    for (Int_t e = 1; e <= fHistdNdedp3Cor.GetNbinsX(); e++)  {
      Double_t eta = fHistdNdedp3Cor.GetXaxis()->GetBinCenter(e);
      Bool_t fwdCov = (hfwd.GetBinContent(e, 0) != 0);
      Bool_t centCov = (hcent.GetBinContent(e, 0) != 0);
      if (!fwdCov && !centCov) continue;
      else fHistdNdedp3Cor.SetBinContent(e, 0, 1);
      for (Int_t p = 1; p <= fHistdNdedp3Cor.GetNbinsY(); p++) {
	Double_t phi = fHistdNdedp3Cor.GetYaxis()->GetBinCenter(p);
	Int_t n = 0;
	Double_t cont = 0.;
	if (fwdCov) {
	  cont += hfwd.GetBinContent(e, p);
	  n++;
	}
	if (centCov) {
	  cont += hcent.GetBinContent(e, p);
	  n++;
	}
	if (cont == 0 || n == 0) continue;
	cont /= n;
	fHistdNdedp3Cor.Fill(eta, phi, cont);
      }
    }
  // VZERO, SPD input, here we do not average but cut to avoid
  // (too much) overlap.
  } else if ((fFlowFlags & kVZERO)) {
    // VZERO loop
    for (Int_t eV = 1; eV <= hfwd.GetNbinsX(); eV++) {
      Double_t eta = hfwd.GetXaxis()->GetBinLowEdge(eV)+0.1;
      if (hfwd.GetBinContent(eV, 0) == 0) continue;
      else { 
        Int_t he = fHistdNdedp3Cor.GetXaxis()->FindBin(eta);
	fHistdNdedp3Cor.SetBinContent(he, 0, 1);
      }
      for (Int_t p = 1; p <= hfwd.GetNbinsY(); p++) { 
        Double_t phi = hfwd.GetYaxis()->GetBinCenter(p);
        Double_t cont = hfwd.GetBinContent(eV, p);
        fHistdNdedp3Cor.Fill(eta, phi, cont);
      }
    }
    // SPD loop
    Int_t eSs = hcent.GetXaxis()->FindBin(-1.99);
    Int_t eSe = hcent.GetXaxis()->FindBin(1.99);
    for (Int_t eS = eSs; eS <= eSe; eS++) {
      Double_t eta = hcent.GetXaxis()->GetBinCenter(eS);
      if (hcent.GetBinContent(eS, 0) == 0) continue;
      else {
        Int_t he = fHistdNdedp3Cor.GetXaxis()->FindBin(eta);
	fHistdNdedp3Cor.SetBinContent(he, 0, 1);
      }
      for (Int_t p = 1; p <= hcent.GetNbinsY(); p++) {
        Double_t phi = hcent.GetYaxis()->GetBinCenter(p);
        Double_t cont = hcent.GetBinContent(eS, p);
        fHistdNdedp3Cor.Fill(eta, phi, cont);
      }
    }
  }
  return fHistdNdedp3Cor;
}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::FillTracks(VertexBin* bin, UShort_t mode) const
{
  // 
  //  Get TPC tracks to use for reference flow.
  //
  //  Return: TObjArray with tracks
  //
  TObjArray* trList = 0;
  AliESDEvent* esdEv = 0;
  if (AliForwardUtil::CheckForAOD() == 1) // AOD tracks
    trList = static_cast<TObjArray*>(fAOD->GetTracks());
  else 
    esdEv = dynamic_cast<AliESDEvent*>(InputEvent());
  
  Bool_t useEvent = bin->FillTracks(trList, esdEv, fTrackCuts, mode);
  return useEvent;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::Terminate(Option_t */*option*/)
{
  //
  //  Calls the finalize function, done at the end of the analysis
  //
  //  Parameters:
  //   option: Not used
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

  if ((fFlowFlags & kEtaGap) || (fFlowFlags & kTracks)) {
    TParameter<Double_t>* etaGap = new TParameter<Double_t>("EtaGap", fEtaGap);
    fOutputList->Add(etaGap);
  }
  // We only add axes in terminate, as TAxis object do not merge well,
  // and so we get a mess when running on the grid if we put them in the sum list...
  fVtxAxis->SetName("VtxAxis");
  fOutputList->Add(fVtxAxis);
  fCentAxis->SetName("CentAxis");
  fOutputList->Add(fCentAxis);

  // Run finalize on VertexBins
  Finalize();

  // Loop over output to get dN/deta hists - used for diagnostics
  TIter next(fOutputList);
  TObject* o = 0;
  TString name;
  TH2D* dNdeta = 0;
  TH1D* cent = 0;
  while ((o = next())) {
    name = o->GetName();
    if (name.Contains("dNdeta")) {
      dNdeta = dynamic_cast<TH2D*>(o);
      name.ReplaceAll("dNdeta", "cent");
      name.ReplaceAll("Ref", "");
      name.ReplaceAll("Diff", "");
      cent = dynamic_cast<TH1D*>(fOutputList->FindObject(name.Data()));
      if (!dNdeta || !cent) continue;
      for (Int_t cBin = 1; cBin <= dNdeta->GetNbinsY(); cBin++) {
        Double_t nEvents = cent->GetBinContent(cBin);
        if (nEvents == 0) continue;
	for (Int_t eBin = 1; eBin <= dNdeta->GetNbinsX(); eBin++) {
	  dNdeta->SetBinContent(eBin, cBin, dNdeta->GetBinContent(eBin, cBin)/nEvents);
	  dNdeta->SetBinError(eBin, cBin, dNdeta->GetBinError(eBin, cBin)/nEvents);
	}
      }
    }
  }   

  // Loop over output and make 1D projections for fast look at results
  MakeCentralityHists(fOutputList);
  TList* vtxList = (TList*)fOutputList->FindObject("vtxList");
  if (vtxList) MakeCentralityHists(vtxList);
  TList* nuaList = (TList*)fOutputList->FindObject("NUATerms");
  TIter nextNua(nuaList);
  o = 0;
  TH2D* h = 0;
  while ((o = nextNua())) {
    if (!(h = dynamic_cast<TH2D*>(o))) continue;
    Double_t evts = h->GetBinContent(0, 0);
    if (evts != 0) h->Scale(1./evts);
  }
  if (nuaList) MakeCentralityHists(nuaList);

  PostData(2, fOutputList);

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::Finalize()
{
  //
  //  Finalize command, called by Terminate()
  //

  // Reinitiate vertex bins if Terminate is called separately!
  if (fBinsForward.GetEntries() == 0) InitVertexBins();

  // Iterate over all vertex bins objects and finalize cumulants
  // calculations
  EndVtxBinList(fBinsForward);
  EndVtxBinList(fBinsCentral);

  return;
} 
//_____________________________________________________________________
void AliForwardFlowTaskQC::EndVtxBinList(const TList& list) const
{
  //
  //  Loop over VertexBin list and call terminate on each 
  //
  //  Parameters:
  //   list: VertexBin list
  //
  TIter next(&list);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(next()))) {
    bin->CumulantsTerminate(fSumList, fOutputList);
  }
  return;
}
// _____________________________________________________________________
void AliForwardFlowTaskQC::MakeCentralityHists(TList* list) const
{
  //
  // Loop over a list containing a TH2D with flow results
  // and project to TH1's in specific centrality bins
  //
  // Parameters:
  //  list: Flow results list
  //
  TH2D* hist2D = 0;
  TList* centList = 0;
  TH1D* hist1D = 0;
  TObject* helper = 0;
  TIter nextHist(list);
  while ((helper = dynamic_cast<TObject*>(nextHist()))) {
    if (!(hist2D = dynamic_cast<TH2D*>(helper))) continue;
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
  //  Function to check that an AOD event meets the cuts
  //
  //  Parameters: 
  //   AliAODForwardMult: forward mult object with trigger and vertex info
  //
  //  Return: false if there is no trigger or if the centrality or vertex
  //  is out of range. Otherwise true.
  //

  // First check for trigger
  if (!CheckTrigger(aodfm)) {
    fHistEventSel->Fill(kNoTrigger);
    return kFALSE;
  }

  // Then check for centrality
  if (!GetCentrality(aodfm)) {
    return kFALSE;
  }

  // And finally check for vertex
  if (!GetVertex(aodfm)) {
    return kFALSE;
  }

  // Ev. accepted - filling diag. hists
  fHistCent->Fill(fCent);
  fHistVertexSel->Fill(fVtx);
  fHistEventSel->Fill(kOK);
  
  return kTRUE;
}
// _____________________________________________________________________
Bool_t AliForwardFlowTaskQC::CheckTrigger(const AliAODForwardMult* aodfm) const 
{
  //
  //  Function to look for a trigger string in the event.
  //  First check for info in forward mult object, if not there, use the AOD header
  //
  //  Parameters: 
  //   AliAODForwardMult: forward mult object with trigger and vertex info
  //
  //  Return: true if offline trigger is present
  //
  if (aodfm) return aodfm->IsTriggerBits(AliAODForwardMult::kOffline);
  // this may need to be changed for 2011 data to handle kCentral and so on...
  else return (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
                 ->IsEventSelected() & AliVEvent::kMB);
}
// _____________________________________________________________________
Bool_t AliForwardFlowTaskQC::GetCentrality(const AliAODForwardMult* aodfm) 
{
  //
  //  Function to look get centrality of the event.
  //  First check for info in forward mult object, if not there, use the AOD header
  //
  //  Parameters: 
  //   AliAODForwardMult: forward mult object with trigger and vertex info
  //
  //  Return: true if centrality determination is present
  //
  if (aodfm) {
    if (aodfm->HasCentrality()) {
      fCent = (Double_t)aodfm->GetCentrality();
      if (fCentAxis->GetXmin() > fCent || fCent >= fCentAxis->GetXmax()) {
	fHistEventSel->Fill(kInvCent);
	return kFALSE;
      }
    }
    else {
      fCent = 97.5;
      fHistEventSel->Fill(kNoCent);
    }
    return kTRUE;
  } else {
    AliCentrality* aodCent = fAOD->GetCentrality();
    if (aodCent) {
      fCent = (Double_t)aodCent->GetCentralityPercentile("V0M");
      if (fCentAxis->GetXmin() > fCent || fCent >= fCentAxis->GetXmax()) {
	fHistEventSel->Fill(kInvCent);
	return kFALSE;
      }
    }
    else {
      fCent = 97.5;
      fHistEventSel->Fill(kNoCent);
    }
    return kTRUE;
  } 
}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::GetVertex(const AliAODForwardMult* aodfm) 
{
  //
  //  Function to look for vertex determination in the event.
  //  First check for info in forward mult object, if not there, use the AOD header
  //
  //  Parameters: 
  //   AliAODForwardMult: forward mult object with trigger and vertex info
  //
  //  Return: true if vertex is determined
  //
  if (aodfm) {
    if (aodfm->HasIpZ()) {
      fVtx = aodfm->GetIpZ();
      if (fVtx < fVtxAxis->GetXmin() || fVtx >= fVtxAxis->GetXmax()) {
	fHistEventSel->Fill(kInvVtx);
	return kFALSE;
      }
    } else {
      fVtx = 9999;
      fHistEventSel->Fill(kNoVtx);
      return kFALSE;
    }
    return kTRUE;
  } else {
    AliAODVertex* aodVtx = fAOD->GetPrimaryVertex();
    if (aodVtx) {
      fVtx = aodVtx->GetZ();
      if (fVtx < fVtxAxis->GetXmin() || fVtx >= fVtxAxis->GetXmax()) {
	fHistEventSel->Fill(kInvVtx);
	return kFALSE;
      }
    } else {
      fVtx = 9999;
      fHistEventSel->Fill(kNoVtx);
      return kFALSE;
    }
    return kTRUE;
  }
}
// _____________________________________________________________________
AliVVZERO* AliForwardFlowTaskQC::GetVZERO() const
{
  //
  //  Get VZERO object from ESD or AOD
  //
  //  Return: VZERO data object
  //
  AliVVZERO* vzero = 0;
  // Get input type
  UShort_t input = AliForwardUtil::CheckForAOD();
  switch (input) {
    // If AOD input, simply get the track array from the event
    case 1: vzero = (AliVVZERO*)fAOD->GetVZEROData();
	    break;
    case 2: {
    // If ESD input get event, apply track cuts
	      AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
	      if (!esd) return 0;
	      vzero = (AliVVZERO*)esd->GetVZEROData();
	      break;
	    }
    default: AliFatal("Neither ESD or AOD input. This should never happen");
    	    break;
  }
  return vzero;
}
// _____________________________________________________________________
void AliForwardFlowTaskQC::FillVZEROHist(AliVVZERO* vzero)
{
  //
  //  Loops over VZERO data object and fill up d^2N/detadphi histogram for flow analysis
  //
  //  Parameters:
  //   vzero: VZERO AOD data object
  //
  Int_t ring = 0;
  Int_t bin = 0;
  Double_t eta = 0;
  // Sort of right for 2010 data, do not use for 2011!
  Double_t eq[64] = { 1.43536, 1.45727, 1.44993, 1.30051, 1.17425, 1.2335, 1.22247, 1.14362, 
  		      1.14647, 1.25208, 1.17681, 1.21642, 1.16604, 1.05532, 1.03212, 1.1032, 
		      1.22941, 1.36986, 1.14652, 1.20056, 0.927086, 1.10809, 1.03343, 1.29472, 
		      1.21204, 1.29217, 1.2003, 2.10382, 1.28513, 1.40558, 1.25784, 1.21848, 
		      0.475162, 0.50421, 0.503617, 0.512471, 0.515276, 0.39831, 0.415199, 0.444664, 
		      0.521922, 0.785915, 0.703658, 0.832479, 0.77461, 0.73129, 0.778697, 0.710265, 
		      0.89686, 0.967688, 0.974225, 0.873445, 0.811096, 0.828493, 0.889609, 0.586056, 
		      1.15877, 0.954656, 0.914557, 0.979028, 1.04907, 0.748518, 0.928043, 0.98175 };

  for (Int_t i = 0; i < 64; i++) {
    if (i % 8 == 0) {
      ring++;
      bin = (ring < 5 ? 11-ring : ring-3);
      eta = fHistdNdedpV0.GetXaxis()->GetBinCenter(bin);
      fHistdNdedpV0.SetBinContent(bin, 0, 1);
    }
    Float_t amp = vzero->GetMultiplicity(i);
    amp /= eq[i];
    Double_t phi = TMath::Pi()/8.+TMath::TwoPi()*i/8.;
    while (phi > TMath::TwoPi()) phi -= TMath::TwoPi();
    fHistdNdedpV0.Fill(eta, phi, amp);
  }
}
//_____________________________________________________________________
AliForwardFlowTaskQC::VertexBin::VertexBin()
  : TNamed(),
    fMaxMoment(0),   // Max flow moment for this vertexbin
    fVzMin(0),       // Vertex z-coordinate min [cm]
    fVzMax(0),       // Vertex z-coordinate max [cm]
    fType(),         // Data type name e.g., FMD/SPD/FMDTR/SPDTR/MC
    fFlags(0),       // Flow flags
    fSigmaCut(-1),   // Sigma cut to remove outlier events
    fEtaGap(-1),     // Eta gap value
    fEtaLims(),      // Limits for binning in 3Cor method
    fCumuRef(),      // Histogram for reference flow
    fCumuDiff(),     // Histogram for differential flow
    fCumuHists(0,0), // CumuHists object for keeping track of results
    fCumuNUARef(),   // Histogram for ref NUA terms
    fCumuNUADiff(),  // Histogram for diff NUA terms
    fdNdedpRefAcc(), // Diagnostics histogram for acc. maps
    fdNdedpDiffAcc(),// Diagnostics histogram for acc. maps
    fOutliers(),     // Histogram for sigma distribution
    fDebug()         // Debug level
{
  //
  //  Default constructor
  //
}
//_____________________________________________________________________
AliForwardFlowTaskQC::VertexBin::VertexBin(Int_t vLow, Int_t vHigh, 
                                           UShort_t moment, TString name,
                                           UShort_t flags, Double_t cut,
                                           Double_t etaGap)
  : TNamed("", ""),
    fMaxMoment(moment),  // Max flow moment for this vertexbin
    fVzMin(vLow),        // Vertex z-coordinate min [cm]
    fVzMax(vHigh),       // Vertex z-coordinate max [cm]
    fType(name),         // Data type name e.g., FMD/SPD/FMDTR/SPDTR/MC
    fFlags(flags),       // Flow flags 
    fSigmaCut(cut),      // Sigma cut to remove outlier events
    fEtaGap(etaGap),     // Eta gap value
    fEtaLims(),          // Limits for binning in 3Cor method
    fCumuRef(),          // Histogram for reference flow
    fCumuDiff(),         // Histogram for differential flow
    fCumuHists(moment,0),// CumuHists object for keeping track of results
    fCumuNUARef(),       // Histogram for ref NUA terms
    fCumuNUADiff(),      // Histogram for diff NUA terms
    fdNdedpRefAcc(),     // Diagnostics histogram for acc. maps
    fdNdedpDiffAcc(),    // Diagnostics histogram for acc. maps
    fOutliers(),         // Histogram for sigma distribution
    fDebug(0)            // Debug level
{
  //
  //  Constructor
  //
  //  Parameters
  //   vLow: min z-coordinate
  //   vHigh: max z-coordinate
  //   moment: max flow moment
  //   name: data type name (FMD/SPD/FMDTR/SPDTR/MC)
  //   flags: flow flags
  //   cut: sigma cut
  //   etaGap: eta-gap value
  //
  fType.ToUpper();

  SetName(Form("%svertexBin%d_%d_%d%s", fType.Data(), moment, vLow, vHigh, GetQCType(fFlags)));
  SetTitle(Form("%svertexBin%d_%d_%d%s", fType.Data(), moment, vLow, vHigh, GetQCType(fFlags)));
  
  fDebug = AliAnalysisManager::GetAnalysisManager()->GetDebugLevel();
  if (fDebug > 0) Printf("AliForwardFlowTaskQC::VertexBin()\tDebugMode: %d", fDebug);

  // Set limits for 3 correlator method
  if ((fFlags & kFMD)) {
    fEtaLims[0] = -6.;
    fEtaLims[1] = -1.5;
    fEtaLims[2] = -0.5;
    fEtaLims[3] = 2.;
    fEtaLims[4] = 3.;
    fEtaLims[5] = 6.;
  } else if ((fFlags & kVZERO)) {
    fEtaLims[0] = -6;
    fEtaLims[1] = -2.7;
    fEtaLims[2] = -2.0;
    fEtaLims[3] = 2.0;
    fEtaLims[4] = 3.9;
    fEtaLims[5] = 6;
  }
}
//_____________________________________________________________________
AliForwardFlowTaskQC::VertexBin&
AliForwardFlowTaskQC::VertexBin::operator=(const AliForwardFlowTaskQC::VertexBin& o)
{
  //
  //  Assignment operator
  //
  //  Parameters
  //   o: AliForwardFlowTaskQC::VertexBin
  //
  if (&o == this) return *this;
  fMaxMoment     = o.fMaxMoment;
  fVzMin         = o.fVzMin;
  fVzMax         = o.fVzMax;
  fType          = o.fType;
  fFlags         = o.fFlags;
  fSigmaCut      = o.fSigmaCut;
  fEtaGap        = o.fEtaGap;
  fCumuRef       = o.fCumuRef;
  fCumuDiff      = o.fCumuDiff;
  fCumuHists     = o.fCumuHists;
  fCumuNUARef    = o.fCumuNUARef;
  fCumuNUADiff   = o.fCumuNUADiff;
  fdNdedpRefAcc  = o.fdNdedpRefAcc;
  fdNdedpDiffAcc = o.fdNdedpDiffAcc;
  fOutliers      = o.fOutliers;
  fDebug         = o.fDebug;
  for (UInt_t i = 0; i < sizeof(fEtaLims)/sizeof(Double_t); i++) fEtaLims[i] = o.fEtaLims[i];

  return *this;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::AddOutput(TList* outputlist, TAxis* centAxis)
{
  //
  //  Add histograms to outputlist
  //
  //  Parameters
  //   outputlist: list of histograms
  //   centAxis: centrality axis
  //

  // First we try to find an outputlist for this vertexbin
  TList* list = (TList*)outputlist->FindObject(Form("%svertex_%d_%d%s", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)));
  // If it doesn't exist we make one
  if (!list) {
    list = new TList();
    list->SetName(Form("%svertex_%d_%d%s", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)));
    outputlist->Add(list);
  }

  // Get bin numbers and binning defined
  Int_t nHBins = GetBinNumberSin();
  Int_t nEtaBins = 24; 
  if ((fFlags & k3Cor)) {
    if ((fFlags & kFMD)) nEtaBins = 24;
    else if ((fFlags & kVZERO)) nEtaBins = 19;
  }
  else if ((fFlags & kVZERO) && (fFlags & kEtaGap)) nEtaBins = 19;
  else if ((fFlags & kVZERO)) nEtaBins = 11;

  Double_t vzeroBins[12] = { -6, -3.7, -3.2, -2.7, -2.2, -1.7, 
                             2.8, 3.4,  3.9,  4.5,  5.1, 6 };
  Double_t vzeroBins2[20] = { -6, -3.7, -3.2, -2.7, -2.2, // VZERO
                             -2.0, -1.5, -1.0, -0.5 , 0., 0.5, 1.0, 1.5, 2.0, // SPD
	  		    2.8, 3.4,  3.9,  4.5,  5.1, 6 }; // VZERO

  Int_t nRefBins = nEtaBins; // needs to be something as default
  if ((fFlags & kStdQC)) {
    if ((fFlags & kSymEta) && !((fFlags & kTracks) && (fFlags & kSPD))) nRefBins = 1;
    else nRefBins = 2;
  } else if ((fFlags & kEtaGap )) {
    nRefBins = 2;
  } else if ((fFlags & k3Cor)) {
    nRefBins = 24;
  }

  // We initiate the reference histogram
  fCumuRef = new TH2D(Form("%s_%d_%d%s_ref", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)), 
                      Form("%s_%d_%d%s_ref", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)), 
			nRefBins, -6., 6., 
                        nHBins, 0.5, nHBins+0.5);
  if ((fFlags & kVZERO) && (fFlags & k3Cor)) fCumuRef->GetXaxis()->Set(nEtaBins, vzeroBins2);
  SetupNUALabels(fCumuRef->GetYaxis());
  list->Add(fCumuRef);
  // We initiate the differential histogram
  fCumuDiff = new TH2D(Form("%s_%d_%d%s_diff", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)),
                       Form("%s_%d_%d%s_diff", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)),
                       nEtaBins, -6., 6., nHBins, 0.5, nHBins+0.5);
  if ((fFlags & kVZERO)) {
    if ((fFlags & kSPD) || (fFlags & k3Cor) || (fFlags & kEtaGap)) 
      fCumuDiff->GetXaxis()->Set(nEtaBins, vzeroBins2);
    else fCumuDiff->GetXaxis()->Set(nEtaBins, vzeroBins);
  }
  SetupNUALabels(fCumuDiff->GetYaxis());
  list->Add(fCumuDiff);

  // Cumulants sum hists 
  Int_t cBins = centAxis->GetNbins();
  fCumuHists.ConnectList(Form("%sCumu_%d_%d%s", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)), list);
  TH3D* cumuHist = 0;
  Int_t nC2Bins = ((fFlags & kEtaGap) || (fFlags & k3Cor) ? kW2 : k3pWeight);
  Int_t nC4Bins = ((fFlags & kEtaGap) ? kW2 : ((fFlags & k3Cor) ? kW4 : kSinphi1phi2phi3p));
  for (Int_t n = 2; n <= fMaxMoment; n++) {
    // Initiate the ref cumulant sum histogram
    cumuHist = new TH3D(Form("%sv%d_vertex_%d_%d%s_cumuRef", fType.Data(), n, fVzMin, fVzMax, GetQCType(fFlags)),
			Form("%sv%d_vertex_%d_%d%s_cumuRef", fType.Data(), n, fVzMin, fVzMax, GetQCType(fFlags)), 
                        nRefBins, -6., 6., 
			cBins, 0., 100., nC2Bins, 0.5, nC2Bins+0.5);
    if ((fFlags & kVZERO) && (fFlags & k3Cor)) cumuHist->GetXaxis()->Set(nEtaBins, vzeroBins2);
    cumuHist->GetYaxis()->Set(cBins, centAxis->GetXbins()->GetArray());
    fCumuHists.Add(cumuHist);
    // Initiate the diff cumulant sum histogram
    cumuHist = new TH3D(Form("%sv%d_vertex_%d_%d%s_cumuDiff", fType.Data(), n, fVzMin, fVzMax, GetQCType(fFlags)),
			Form("%sv%d_vertex_%d_%d%s_cumuDiff", fType.Data(), n, fVzMin, fVzMax, GetQCType(fFlags)),
			nEtaBins, -6., 6., cBins, 0., 100., nC4Bins, 0.5, nC4Bins+0.5);
    if ((fFlags & kVZERO)) {
      if ((fFlags & kSPD) || (fFlags & k3Cor) || (fFlags & kEtaGap)) 
        cumuHist->GetXaxis()->Set(nEtaBins, vzeroBins2);
      else cumuHist->GetXaxis()->Set(nEtaBins, vzeroBins);
    }
    cumuHist->GetYaxis()->Set(cBins, centAxis->GetXbins()->GetArray());
    fCumuHists.Add(cumuHist);
  }

  // Common NUA histograms
  fCumuNUARef = new TH3D(Form("%s_vertex_%d_%d%s_cumuNUARef", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)),
			 Form("%s_vertex_%d_%d%s_cumuNUARef", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)),
			 nRefBins, -6., 6., 
                          cBins, 0., 100., nHBins, 0.5, nHBins+0.5);
  if ((fFlags & kVZERO) && (fFlags & k3Cor)) fCumuNUARef->GetXaxis()->Set(nEtaBins, vzeroBins2);
  fCumuNUARef->GetYaxis()->Set(cBins, centAxis->GetXbins()->GetArray());
  SetupNUALabels(fCumuNUARef->GetZaxis());
  fCumuNUARef->Sumw2();
  list->Add(fCumuNUARef);
  
  fCumuNUADiff = new TH3D(Form("%s_vertex_%d_%d%s_cumuNUADiff", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)),
                          Form("%s_vertex_%d_%d%s_cumuNUADiff", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)),
                          nEtaBins, -6., 6., cBins, 0., 100., nHBins, 0.5, nHBins+0.5);
  if ((fFlags & kVZERO)) {
    if ((fFlags & kSPD) || (fFlags & k3Cor) || (fFlags & kEtaGap)) 
      fCumuNUADiff->GetXaxis()->Set(nEtaBins, vzeroBins2);
    else fCumuNUADiff->GetXaxis()->Set(nEtaBins, vzeroBins);
  }
  fCumuNUADiff->GetYaxis()->Set(cBins, centAxis->GetXbins()->GetArray());
  SetupNUALabels(fCumuNUADiff->GetZaxis());
  fCumuNUADiff->Sumw2();
  list->Add(fCumuNUADiff);

  // We create diagnostic histograms.
  TList* dList = (TList*)outputlist->FindObject("Diagnostics");
  if (!dList) AliFatal("No diagnostics list found");

  // Acceptance hist
  Int_t nPhiBins = ((fFlags & kFMD) ? 20 : 8);
  fdNdedpRefAcc = new TH2F(Form("h%sdNdedpRefAcc_%d_%d%s", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)),
    Form("%s reference flow acceptance map for %d cm < v_{z} < %d cm", fType.Data(), fVzMin, fVzMax),
    nEtaBins, -6., 6., nPhiBins, 0, TMath::TwoPi());
  if ((fFlags & kVZERO)) {
    if ((fFlags & kSPD) || (fFlags & k3Cor) || (fFlags & kEtaGap)) 
      fdNdedpRefAcc->GetXaxis()->Set(nEtaBins, vzeroBins2);
    else fdNdedpRefAcc->GetXaxis()->Set(nEtaBins, vzeroBins);
  }
  dList->Add(fdNdedpRefAcc);

  fdNdedpDiffAcc = new TH2F(Form("h%sdNdedpDiffAcc_%d_%d%s", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)),
    Form("%s differential flow acceptance map for %d cm < v_{z} < %d cm", fType.Data(), fVzMin, fVzMax),
    nEtaBins, -6., 6., nPhiBins, 0, TMath::TwoPi());
  if ((fFlags & kVZERO)) {
    if ((fFlags & kSPD) || (fFlags & k3Cor) || (fFlags & kEtaGap)) 
      fdNdedpDiffAcc->GetXaxis()->Set(nEtaBins, vzeroBins2);
    else fdNdedpDiffAcc->GetXaxis()->Set(nEtaBins, vzeroBins);
  }
  dList->Add(fdNdedpDiffAcc);
  
  fOutliers = new TH2F(Form("hOutliers_%s_%d_%d%s", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)),
		       Form("Maximum #sigma from mean N_{ch} pr. bin - %s, %d < v_{z} < %d",
		       fType.Data(), fVzMin, fVzMax), 
		       20, 0., 100., 500, 0., ((fFlags & kMC) ? 15. : 5.));
  dList->Add(fOutliers);
  
  return;
}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::VertexBin::FillHists(TH2D& dNdetadphi, Double_t cent, UShort_t mode) 
{
  // 
  //  Fill reference and differential eta-histograms
  //
  //  Parameters:
  //   dNdetadphi: 2D histogram with input data
  //   cent: centrality
  //   mode: filling mode: kFillRef/kFillDiff/kFillBoth
  //
  if (!fCumuRef) AliFatal("You have not called AddOutput() - Terminating!");
  Bool_t useEvent = kTRUE;

  // Fist we reset histograms
  if ((mode & kReset)) {
    if ((mode & kFillRef))  fCumuRef->Reset();
    if ((mode & kFillDiff)) fCumuDiff->Reset();
  }
  // Then we loop over the input and calculate sum cos(k*n*phi)
  // and fill it in the reference and differential histograms
  Int_t nBadBins = 0;
  Double_t limit = 9999.;
  for (Int_t etaBin = 1; etaBin <= dNdetadphi.GetNbinsX(); etaBin++) {
    Double_t eta = dNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    // Numbers to cut away bad events
    Double_t runAvg = 0;
    Double_t max = 0;
    Int_t nInAvg = 0;
    Double_t avgSqr = 0;
    for (Int_t phiBin = 0; phiBin <= dNdetadphi.GetNbinsY(); phiBin++) {
      // Check for acceptance
      if (phiBin == 0) {
        if (dNdetadphi.GetBinContent(etaBin, 0) == 0) break;
        // Central limit for eta gap break for reference flow
	if ((fFlags & kEtaGap) && (mode & kFillRef) && 
	     TMath::Abs(eta) < fEtaGap) break;
	// Backward and forward eta gap break for reference flow
	if ((fFlags & kEtaGap) && (mode & kFillRef) && TMath::Abs(eta) > TMath::Abs(limit)) break;
	if ((fFlags & kStdQC) && (fFlags & kMC) && !(fFlags & kTracks)) {
	  if (!(fFlags & kSPD) && TMath::Abs(eta) < 1.75) break; 
	  if ((fFlags & kSPD) && TMath::Abs(eta) > 2.00) break;
	}
	if (limit > 1e3) limit = dNdetadphi.GetXaxis()->GetBinLowEdge(etaBin);
	continue;
      } // End of phiBin == 0
      Double_t phi = dNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      Double_t weight = dNdetadphi.GetBinContent(etaBin, phiBin);
        
      // We calculate the average Nch per. bin
      avgSqr += weight*weight;
      runAvg += weight;
      nInAvg++;
      if (weight == 0) continue;
      if (weight > max) max = weight;
      // Fill into Cos() and Sin() hists
      if ((mode & kFillRef) && !((fFlags & kTracks) && (fFlags & kMC) && TMath::Abs(eta) > 0.75)) {
	fCumuRef->Fill(eta, 0., weight);// mult goes in underflowbin - no visual, but not needed?
        fdNdedpRefAcc->Fill(eta, phi, weight);
      }
      if ((mode & kFillDiff)) {
	fCumuDiff->Fill(eta, 0., weight);
        fdNdedpDiffAcc->Fill(eta, phi, weight);
      }
      for (Int_t n = 1; n <= 2*fMaxMoment; n++) {
	Double_t cosBin = fCumuDiff->GetYaxis()->GetBinCenter(GetBinNumberCos(n));
	Double_t sinBin = fCumuDiff->GetYaxis()->GetBinCenter(GetBinNumberSin(n));
	Double_t cosnPhi = weight*TMath::Cos(n*phi);
	Double_t sinnPhi = weight*TMath::Sin(n*phi);
        // fill ref
	if ((mode & kFillRef) && !((fFlags & kTracks) && (fFlags & kMC) && TMath::Abs(eta) > 0.75)) {
	  fCumuRef->Fill(eta, cosBin, cosnPhi);
	  fCumuRef->Fill(eta, sinBin, sinnPhi);
	}
	// fill diff
	if ((mode & kFillDiff)) {
	  fCumuDiff->Fill(eta, cosBin, cosnPhi);
	  fCumuDiff->Fill(eta, sinBin, sinnPhi);
	}
      } // End of NUA loop
    } // End of phi loop
    // Outlier cut calculations
    if (nInAvg > 0) {
      runAvg /= nInAvg;
      avgSqr /= nInAvg;
      Double_t stdev = (nInAvg > 1 ? TMath::Sqrt(nInAvg/(nInAvg-1))*TMath::Sqrt(avgSqr - runAvg*runAvg) : 0);
      Double_t nSigma = (stdev == 0 ? 0 : (max-runAvg)/stdev);
      if (fSigmaCut > 0. && nSigma >= fSigmaCut && cent < 60) nBadBins++;
      else nBadBins = 0;
      fOutliers->Fill(cent, nSigma);
      // We still finish the loop, for fOutliers to make sense, 
      // but we do no keep the event for analysis 
      if (nBadBins > 3) useEvent = kFALSE;
    }
  } // End of eta bin

  return useEvent;
}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::VertexBin::FillTracks(TObjArray* trList, AliESDEvent* esd,
                                                   AliAnalysisFilter* trFilter, UShort_t mode) 
{
  // 
  //  Fill reference and differential eta-histograms
  //
  //  Parameters:
  //   trList: list of tracks
  //   mode: filling mode: kFillRef/kFillDiff/kFillBoth
  //
  if (!fCumuRef) AliFatal("You have not called AddOutput() - Terminating!");
  if (!trList && !esd) {
    AliError("FillTracks: No AOD track list or ESD event - something might be wrong!");
    return kFALSE;
  }

  // Fist we reset histograms
  if ((mode & kReset)) {
    if ((mode & kFillRef))  fCumuRef->Reset();
    if ((mode & kFillDiff)) fCumuDiff->Reset();
  }

  // Then we loop over the input and calculate sum cos(k*n*phi)
  // and fill it in the reference and differential histograms
  Int_t nTr = 0;
  if (trList) nTr = trList->GetEntries();
  if (esd) nTr = esd->GetNumberOfTracks();
  if (nTr == 0) return kFALSE;
  AliVTrack* tr = 0;
  // Cuts for AOD tracks (have already been applied to ESD tracks) - except dEdx
//  const tpcdEdx = 10;
  for (Int_t i = 0; i < nTr; i++) { // track loop
    tr = (trList ? (AliVTrack*)trList->At(i) : (AliVTrack*)esd->GetTrack(i));
    if (!tr) continue;
    if (esd) {
      AliESDtrack* esdTr = (AliESDtrack*)tr;
      if (!trFilter->IsSelected(esdTr)) continue;
    }
    else if (trList) { // If AOD input
      Double_t pTMin = 0, pTMax = 0, etaMin = 0, etaMax = 0, minNCl = 0;
      UInt_t bit = 0;
      if ((fFlags & kTPC) == kTPC)    pTMin = 0.2, pTMax = 5., etaMin = -0.8, etaMax = 0.8, minNCl = 70, bit = 128;
      if ((fFlags & kHybrid) == kHybrid) pTMin = 0.2, pTMax = 5., etaMin = -0.8, etaMax = 0.8, minNCl = 70, bit = 272;

      AliAODTrack* aodTr = (AliAODTrack*)tr;
      if (aodTr->GetID() > -1) continue;
      if (!aodTr->TestFilterBit(bit) || aodTr->Pt() > pTMax || aodTr->Pt() < pTMin || 
	aodTr->Eta() > etaMax || aodTr->Eta() < etaMin || aodTr->GetTPCNcls() < minNCl) continue;
    }

//    if (tr->GetTPCsignal() < tpcdEdx) continue;
    // Track accepted
    Double_t eta = tr->Eta();
    if (((fFlags & kSPD) || (fFlags & kEtaGap)) && TMath::Abs(eta) < fEtaGap) continue;
    Double_t phi = tr->Phi();
    Double_t weight = 1.;

    if ((mode & kFillRef)) {
      fCumuRef->Fill(eta, 0., weight);// mult goes in underflowbin - no visual, but not needed?
      fdNdedpRefAcc->Fill(eta, phi, weight);
    }
    if ((mode & kFillDiff)) {
      fCumuDiff->Fill(eta, 0., weight);
      fdNdedpDiffAcc->Fill(eta, phi, weight);
    }
    for (Int_t n = 1; n <= 2*fMaxMoment; n++) {
      Double_t cosBin = fCumuDiff->GetYaxis()->GetBinCenter(GetBinNumberCos(n));
      Double_t sinBin = fCumuDiff->GetYaxis()->GetBinCenter(GetBinNumberSin(n));
      Double_t cosnPhi = weight*TMath::Cos(n*phi);
      Double_t sinnPhi = weight*TMath::Sin(n*phi);
      // fill ref
      if ((mode & kFillRef)) {
	fCumuRef->Fill(eta, cosBin, cosnPhi);
	fCumuRef->Fill(eta, sinBin, sinnPhi);
      }
      // fill diff
      if ((mode & kFillDiff)) {
	fCumuDiff->Fill(eta, cosBin, cosnPhi);
	fCumuDiff->Fill(eta, sinBin, sinnPhi);
      }
    } // End of NUA loop
  } // End of track loop
  return kTRUE;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::CumulantsAccumulate(Double_t cent) 
{
  // 
  //  Calculate the Q cumulant up to order fMaxMoment
  //
  //  Parameters:
  //   cent: centrality of event
  //
  if (!fCumuRef) AliFatal("You have not called AddOutput() - Terminating!");

  // Fill out NUA hists
  for (Int_t etaBin = 1; etaBin <= fCumuRef->GetNbinsX(); etaBin++) {
    Double_t eta = fCumuRef->GetXaxis()->GetBinCenter(etaBin);
    if (fCumuRef->GetBinContent(etaBin, 0) <= 3) continue;
    if ((fFlags & kTracks) && (fFlags && kSPD) && !(fFlags & kEtaGap)) eta = -eta;
    for (Int_t qBin = 0; qBin <= fCumuRef->GetNbinsY(); qBin++) {
      fCumuNUARef->Fill(eta, cent, Double_t(qBin), fCumuRef->GetBinContent(etaBin, qBin));
    }
  }
  for (Int_t etaBin = 1; etaBin <= fCumuDiff->GetNbinsX(); etaBin++) {
    Double_t eta = fCumuDiff->GetXaxis()->GetBinCenter(etaBin);
    Int_t    refetaBin = fCumuRef->GetXaxis()->FindBin(eta);
    if (fCumuRef->GetBinContent(refetaBin, 0) <= 3) continue;
    if (fCumuDiff->GetBinContent(etaBin, 0) == 0) continue;
    for (Int_t qBin = 0; qBin <= fCumuDiff->GetNbinsY(); qBin++) {
      fCumuNUADiff->Fill(eta, cent, Double_t(qBin), fCumuDiff->GetBinContent(etaBin, qBin));
    }
  }

  // We create the objects needed for the analysis
  TH3D* cumuRef = 0; 
  TH3D* cumuDiff = 0; 
  // For each n we loop over the hists
  for (Int_t n = 2; n <= fMaxMoment; n++) {
    cumuRef  = (TH3D*)fCumuHists.Get('r',n);
    cumuDiff = (TH3D*)fCumuHists.Get('d',n);
    Int_t prevRefEtaBin = 0;

    // Per mom. quantities
    Double_t dQnReA = 0, dQnImA = 0, multA = 0; 
    Double_t dQnReB = 0, dQnImB = 0, multB = 0;
    Double_t dQ2nReA = 0, dQ2nImA = 0;
    Double_t two = 0, w2 = 0, four = 0, w4 = 0;
    Double_t cosPhi1Phi2 = 0, cosPhi1Phi2Phi3m = 0;
    Double_t sinPhi1Phi2 = 0, sinPhi1Phi2Phi3m = 0;
    for (Int_t etaBin = 1; etaBin <= fCumuDiff->GetNbinsX(); etaBin++) {
      Double_t eta = fCumuDiff->GetXaxis()->GetBinCenter(etaBin);
      Double_t refEta = eta;
      if ((fFlags & kTracks) && (fFlags && kSPD) && !(fFlags & kEtaGap)) refEta = -eta;
      Int_t refEtaBinA = fCumuRef->GetXaxis()->FindBin(refEta);
      if ((fFlags & kEtaGap)) refEta = -eta;
      Int_t refEtaBinB = fCumuRef->GetXaxis()->FindBin(refEta);
      if (refEtaBinA != prevRefEtaBin) {
	// Reference flow
	multA  = fCumuRef->GetBinContent(refEtaBinA, 0);
	dQnReA = fCumuRef->GetBinContent(refEtaBinA, GetBinNumberCos(n));
	dQnImA = fCumuRef->GetBinContent(refEtaBinA, GetBinNumberSin(n));
	dQ2nReA = fCumuRef->GetBinContent(refEtaBinA, GetBinNumberCos(n*2));
	dQ2nImA = fCumuRef->GetBinContent(refEtaBinA, GetBinNumberSin(n*2));
	
	multB  = fCumuRef->GetBinContent(refEtaBinB, 0);
	dQnReB = fCumuRef->GetBinContent(refEtaBinB, GetBinNumberCos(n));
	dQnImB = fCumuRef->GetBinContent(refEtaBinB, GetBinNumberSin(n));

	if (multA <= 3 || multB <= 3) continue; 
	// The reference flow is calculated 
	// 2-particle
	if ((fFlags & kStdQC)) {
	  w2 = multA * (multA - 1.);
	  two = dQnReA*dQnReA + dQnImA*dQnImA - multA;
	} else {
	  w2 = multA * multB;
	  two = dQnReA*dQnReB + dQnImA*dQnImB;
	}
      	cumuRef->Fill(eta, cent, kW2Two, two);
  	cumuRef->Fill(eta, cent, kW2, w2);

	// The reference flow is calculated 
	// 4-particle
	if ((fFlags & kStdQC)) {
	  w4 = multA * (multA - 1.) * (multA - 2.) * (multA - 3.);
	
	  four = 2.*multA*(multA-3.) + TMath::Power((TMath::Power(dQnReA,2.)+TMath::Power(dQnImA,2.)),2.)
		   -4.*(multA-2.)*(TMath::Power(dQnReA,2.) + TMath::Power(dQnImA,2.))
		   -2.*(TMath::Power(dQnReA,2.)*dQ2nReA+2.*dQnReA*dQnImA*dQ2nImA-TMath::Power(dQnImA,2.)*dQ2nReA)
		   +(TMath::Power(dQ2nReA,2.)+TMath::Power(dQ2nImA,2.));

	  cumuRef->Fill(eta, cent, kW4Four, four);
	  cumuRef->Fill(eta, cent, kW4, w4);

	  // NUA
	  cosPhi1Phi2 = dQnReA*dQnReA - dQnImA*dQnImA - dQ2nReA;
	  sinPhi1Phi2 = 2.*dQnReA*dQnImA - dQ2nImA;

	  cosPhi1Phi2Phi3m = dQnReA*(TMath::Power(dQnReA,2)+TMath::Power(dQnImA,2))
			      -dQnReA*dQ2nReA-dQnImA*dQ2nImA-2.*(multA-1)*dQnReA;

	  sinPhi1Phi2Phi3m = -dQnImA*(TMath::Power(dQnReA,2)+TMath::Power(dQnImA,2))
			      +dQnReA*dQ2nImA-dQnImA*dQ2nReA+2.*(multA-1)*dQnImA; 

	  cumuRef->Fill(eta, cent, kCosphi1phi2, cosPhi1Phi2);
	  cumuRef->Fill(eta, cent, kSinphi1phi2, sinPhi1Phi2);
	  cumuRef->Fill(eta, cent, kCosphi1phi2phi3m, cosPhi1Phi2Phi3m);
	  cumuRef->Fill(eta, cent, kSinphi1phi2phi3m, sinPhi1Phi2Phi3m);
	  cumuRef->Fill(eta, cent, k3pWeight, multA*(multA-1.)*(multA-2.));
	} // End of QC{4}
	prevRefEtaBin = refEtaBinA;
      } // End of reference flow
      // For each etaBin bin the necessary values for differential flow is calculated
      Double_t mp = fCumuDiff->GetBinContent(etaBin, 0);
      Double_t pnRe = fCumuDiff->GetBinContent(etaBin, GetBinNumberCos(n));
      Double_t pnIm = fCumuDiff->GetBinContent(etaBin, GetBinNumberSin(n));
      Double_t p2nRe = fCumuDiff->GetBinContent(etaBin, GetBinNumberCos(n*2));
      Double_t p2nIm = fCumuDiff->GetBinContent(etaBin, GetBinNumberSin(n*2));
      if (mp == 0) continue;
      Double_t mq = 0;
      Double_t qnRe = 0;
      Double_t qnIm = 0;
      Double_t q2nRe = 0;
      Double_t q2nIm = 0;

      // Differential flow calculations for each eta bin is done:
      // 2-particle differential flow
      if ((fFlags & kStdQC) && (!(fFlags & kTracks) || ((fFlags & kTracks) && (fFlags & kMC) && !(fFlags & kSPD) && TMath::Abs(eta) < 0.75))) {
	mq = mp;
	qnRe = pnRe;
	qnIm = pnIm;
	q2nRe = p2nRe;
	q2nIm = p2nIm;
      }

      Double_t w2p = mp * multB - mq;
      Double_t twoPrime = pnRe*dQnReB + pnIm*dQnImB - mq;
      
      cumuDiff->Fill(eta, cent, kW2Two, twoPrime);
      cumuDiff->Fill(eta, cent, kW2, w2p);

      if ((fFlags & kEtaGap)) continue;
      // Differential flow calculations for each eta bin bin is done:
      // 4-particle differential flow
      Double_t w4p = (mp * multA - 3.*mq)*(multA - 1.)*(multA - 2.);
   
      Double_t fourPrime = (TMath::Power(dQnReA,2.)+TMath::Power(dQnImA,2.))*(pnRe*dQnReA+pnIm*dQnImA)
	  		  - q2nRe*(TMath::Power(dQnReA,2.)-TMath::Power(dQnImA,2.))
	  		  - 2.*q2nIm*dQnReA*dQnImA
	  		  - pnRe*(dQnReA*dQ2nReA+dQnImA*dQ2nImA)
	  		  + pnIm*(dQnImA*dQ2nReA-dQnReA*dQ2nImA)
	  		  - 2.*multA*(pnRe*dQnReA+pnIm*dQnImA)
	  		  - 2.*(TMath::Power(dQnReA,2.)+TMath::Power(dQnImA,2.))*mq 
	  		  + 6.*(qnRe*dQnReA+qnIm*dQnImA)
	  		  + 1.*(q2nRe*dQ2nReA+q2nIm*dQ2nImA)
	  		  + 2.*(pnRe*dQnReA+pnIm*dQnImA) 
	  		  + 2.*mq*multA 
	  		  - 6.*mq; 

      cumuDiff->Fill(eta, cent, kW4Four, fourPrime);
      cumuDiff->Fill(eta, cent, kW4, w4p);

      // NUA
      Double_t cosPsi1Phi2 = pnRe*dQnReA - pnIm*dQnImA - q2nRe;
      Double_t sinPsi1Phi2 = pnRe*dQnImA + pnIm*dQnReA - q2nIm;

      Double_t cosPsi1Phi2Phi3p = pnRe*(TMath::Power(dQnImA,2.)+TMath::Power(dQnReA,2.)-multA)
			    - 1.*(q2nRe*dQnReA+q2nIm*dQnImA)  
			    - mq*dQnReA+2.*qnRe;
   
      Double_t sinPsi1Phi2Phi3p = pnIm*(TMath::Power(dQnImA,2.)+TMath::Power(dQnReA,2.)-multA)
			    - 1.*(q2nIm*dQnReA-q2nRe*dQnImA)  
			    - mq*dQnImA+2.*qnIm; 

      Double_t cosPsi1Phi2Phi3m = pnRe*(TMath::Power(dQnReA,2.)-TMath::Power(dQnImA,2.))+2.*pnIm*dQnReA*dQnImA
			    - 1.*(pnRe*dQ2nReA+pnIm*dQ2nImA)  
			    - 2.*mq*dQnReA+2.*qnRe;
   
      Double_t sinPsi1Phi2Phi3m = pnIm*(TMath::Power(dQnReA,2.)-TMath::Power(dQnImA,2.))-2.*pnRe*dQnReA*dQnImA
			    - 1.*(pnIm*dQ2nReA-pnRe*dQ2nImA)
			    + 2.*mq*dQnImA-2.*qnIm;

      cumuDiff->Fill(eta, cent, kCosphi1phi2, cosPsi1Phi2);
      cumuDiff->Fill(eta, cent, kSinphi1phi2, sinPsi1Phi2);
      cumuDiff->Fill(eta, cent, kCosphi1phi2phi3m, cosPsi1Phi2Phi3m);
      cumuDiff->Fill(eta, cent, kSinphi1phi2phi3m, sinPsi1Phi2Phi3m);
      cumuDiff->Fill(eta, cent, k3pWeight, (mp*multA-2.*mq)*(multA-1.));
      cumuDiff->Fill(eta, cent, kCosphi1phi2phi3p, cosPsi1Phi2Phi3p);
      cumuDiff->Fill(eta, cent, kSinphi1phi2phi3p, sinPsi1Phi2Phi3p); 
    } // End of eta loop
    // Event count
    cumuRef->Fill(-7., cent, -0.5, 1.);
  } // End of moment loop
  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::GetLimits(Int_t bin, Int_t& aLow, Int_t& aHigh,
                                                Int_t& bLow, Int_t& bHigh) const
{
  //
  //  Get the limits for the 3 correlator method
  //
  //  Parameters: 
  //   bin  : reference bin #
  //   aLow : Lowest bin to be included in v_A calculations
  //   aHigh: Highest bin to be included in v_A calculations
  //   bLow : Lowest bin to be included in v_B calculations
  //   bHigh: Highest bin to be included in v_B calculations
  //
  if ((fFlags & kFMD)) {
    switch(bin) {
      case 0:
	aLow = 14; aHigh = 15;
	bLow = 20; bHigh = 22;
	break;
      case 1:
	aLow = 16; aHigh = 16;
	bLow = 21; bHigh = 22;
	break;
      case 2:
	aLow =  6; aHigh =  7;
	bLow = 21; bHigh = 22;
	break;
      case 3:
	aLow =  6; aHigh =  7;
	bLow = 12; bHigh = 12; 
	break;
      case 4:
	aLow =  6; aHigh =  8;
	bLow = 13; bHigh = 14;
	break;
      default:
	AliFatal(Form("No limits for this eta region! (%d)", bin));
    }
  } 
  else if ((fFlags & kVZERO)) {
    switch(bin) {
      case 0:
        aLow =  6; aHigh = 13;
        bLow = 17; bHigh = 18;
	break;
      case 1:
        aLow =  6; aHigh =  9;
        bLow = 17; bHigh = 18;
	break;
      case 2:
        aLow =  2; aHigh =  3;
        bLow = 17; bHigh = 18;
	break;
      case 3:
        aLow =  2; aHigh =  3;
        bLow =  6; bHigh =  9;
	break;
      case 4:
        aLow =  2; aHigh =  3;
        bLow =  6; bHigh = 13;
	break;
      default:
	AliFatal(Form("No limits for this eta region! (%d)", bin));
    }
  }
  // Try to catch cases where fEtaLimits and these values do not correspond to each other
  if (aHigh > fCumuNUARef->GetNbinsX() || bHigh > fCumuNUARef->GetNbinsX()) 
    AliFatal(Form("Limits outside vtx range! (%d) - aHigh = %d, bHigh = %d, Nbins = %d", 
      bin, aHigh, bHigh, fCumuNUARef->GetNbinsX()));
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::CumulantsAccumulate3Cor(Double_t cent) 
{
  // 
  //  Calculate the Q cumulant up to order fMaxMoment
  //
  //  Parameters:
  //   cent: centrality of event
  //
  if (!fCumuRef) AliFatal("You have not called AddOutput() - Terminating!");

  // Fill out NUA hists
  for (Int_t etaBin = 1; etaBin <= fCumuRef->GetNbinsX(); etaBin++) {
    Double_t eta = fCumuRef->GetXaxis()->GetBinCenter(etaBin);
    if (fCumuRef->GetBinContent(etaBin, 0) == 0) continue;
    for (Int_t qBin = 0; qBin <= fCumuRef->GetNbinsY(); qBin++) {
      fCumuNUARef->Fill(eta, cent, Double_t(qBin), fCumuRef->GetBinContent(etaBin, qBin));
    }
  }
  for (Int_t etaBin = 1; etaBin <= fCumuDiff->GetNbinsX(); etaBin++) {
    Double_t eta = fCumuDiff->GetXaxis()->GetBinCenter(etaBin);
    if (fCumuDiff->GetBinContent(etaBin, 0) == 0) continue;
    for (Int_t qBin = 0; qBin <= fCumuDiff->GetNbinsY(); qBin++) {
      fCumuNUADiff->Fill(eta, cent, Double_t(qBin), fCumuDiff->GetBinContent(etaBin, qBin));
    }
  }

  // We create the objects needed for the analysis
  TH3D* cumuRef = 0; 
  TH3D* cumuDiff = 0; 
  // For each n we loop over the hists
  for (Int_t n = 2; n <= fMaxMoment; n++) {
    cumuRef  = (TH3D*)fCumuHists.Get('r',n);
    cumuDiff = (TH3D*)fCumuHists.Get('d',n);

    // Per mom. quantities
    Int_t prevLim = 0;
    Int_t aLow = 0, aHigh = 0, bLow = 0, bHigh = 0;
    Double_t dQnReA = 0, dQnImA = 0, multA = 0; 
    Double_t dQnReB = 0, dQnImB = 0, multB = 0;
    Double_t two = 0, w2 = 0;
    for (Int_t etaBin = 1; etaBin <= fCumuDiff->GetNbinsX(); etaBin++) {
      Double_t eta = fCumuDiff->GetXaxis()->GetBinCenter(etaBin);
      if (fEtaLims[prevLim] < eta) {
        GetLimits(prevLim, aLow, aHigh, bLow, bHigh);
	prevLim++;
	multA = 0; dQnReA = 0; dQnImA = 0;
	multB = 0; dQnReB = 0; dQnImB = 0;
	// Reference flow
        for (Int_t a = aLow; a <= aHigh; a++) {
	  multA  += fCumuRef->GetBinContent(a, 0);
	  dQnReA += fCumuRef->GetBinContent(a, GetBinNumberCos(n));
	  dQnImA += fCumuRef->GetBinContent(a, GetBinNumberSin(n));
	}
	for (Int_t b = bLow; b <= bHigh; b++) {
	  multB  += fCumuRef->GetBinContent(b, 0);
	  dQnReB += fCumuRef->GetBinContent(b, GetBinNumberCos(n));
	  dQnImB += fCumuRef->GetBinContent(b, GetBinNumberSin(n));
	}
	// The reference flow is calculated 
	// 2-particle
	w2 = multA * multB;
      	two = dQnReA*dQnReB + dQnImA*dQnImB;
      } // End of reference flow
      cumuRef->Fill(eta, cent, kW2Two, two);
      cumuRef->Fill(eta, cent, kW2, w2);

      // For each etaBin bin the necessary values for differential flow is calculated
      Double_t mp = fCumuDiff->GetBinContent(etaBin, 0);
      Double_t pnRe = fCumuDiff->GetBinContent(etaBin, GetBinNumberCos(n));
      Double_t pnIm = fCumuDiff->GetBinContent(etaBin, GetBinNumberSin(n));
      if (mp == 0) continue;

      // Differential flow calculations for each eta bin is done:
      // 2-particle differential flow
      Double_t w2pA = mp * multA;
      Double_t twoPrimeA = pnRe*dQnReA + pnIm*dQnImA;
      cumuDiff->Fill(eta, cent, kW2Two, twoPrimeA);
      cumuDiff->Fill(eta, cent, kW2, w2pA);

      Double_t w2pB = mp * multB;
      Double_t twoPrimeB = pnRe*dQnReB + pnIm*dQnImB;
      cumuDiff->Fill(eta, cent, kW4Four, twoPrimeB);
      cumuDiff->Fill(eta, cent, kW4, w2pB);
     } // End of eta loop
    // Event count
    cumuRef->Fill(-7., cent, -0.5, 1.);
  } // End of moment loop
  return;

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
  if (!fCumuHists.IsConnected()) {
    TList* list = (TList*)inlist->FindObject(Form("%svertex_%d_%d%s", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)));
    fCumuHists.ConnectList(Form("%sCumu_%d_%d%s", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)), list);

    if (!fCumuNUARef) 
      fCumuNUARef = (TH3D*)list->FindObject(Form("%s_vertex_%d_%d%s_cumuNUARef", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)));
    if (!fCumuNUADiff) 
      fCumuNUADiff = (TH3D*)list->FindObject(Form("%s_vertex_%d_%d%s_cumuNUADiff", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)));
  }
  // Clone to avoid normalization problems when redoing terminate locally
  fCumuNUARef = (TH3D*)fCumuNUARef->Clone(Form("%s_vertex_%d_%d%s_cumuNUARefNorm", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)));
  fCumuNUADiff = (TH3D*)fCumuNUADiff->Clone(Form("%s_vertex_%d_%d%s_cumuNUADiffNorm", fType.Data(), fVzMin, fVzMax, GetQCType(fFlags)));

  // Diagnostics histograms
  TH2I* quality = (TH2I*)outlist->FindObject(Form("hQCQuality%s%s", fType.Data(), GetQCType(fFlags)));
  if (!quality) { 
    quality = MakeQualityHist(Form("hQCQuality%s%s", fType.Data(), GetQCType(fFlags)));
    outlist->Add(quality);
  }
  TH1D* cent = (TH1D*)outlist->FindObject(Form("%s%s_cent", fType.Data(), GetQCType(fFlags)));
  if (!cent) {
    cent = new TH1D(Form("%s%s_cent", fType.Data(), GetQCType(fFlags)), 
                    Form("%s%s_cent", fType.Data(), GetQCType(fFlags)), 
                fCumuNUARef->GetNbinsY(), fCumuNUARef->GetYaxis()->GetXmin(), fCumuNUARef->GetYaxis()->GetXmax());
    cent->GetXaxis()->Set(fCumuNUARef->GetNbinsY(), fCumuNUARef->GetYaxis()->GetXbins()->GetArray());
    outlist->Add(cent);
  }
  TH2D* dNdetaRef = (TH2D*)outlist->FindObject(Form("%s%s_dNdetaRef", fType.Data(), GetQCType(fFlags)));
  if (!dNdetaRef) {
    dNdetaRef = new TH2D(Form("%s%s_dNdetaRef", fType.Data(), GetQCType(fFlags)), 
                               Form("%s%s_dNdetaRef", fType.Data(), GetQCType(fFlags)), 
                fCumuNUARef->GetNbinsX(), fCumuNUARef->GetXaxis()->GetXmin(), fCumuNUARef->GetXaxis()->GetXmax(),
                fCumuNUARef->GetNbinsY(), fCumuNUARef->GetYaxis()->GetXmin(), fCumuNUARef->GetYaxis()->GetXmax());
    dNdetaRef->GetYaxis()->Set(fCumuNUARef->GetNbinsY(), fCumuNUARef->GetYaxis()->GetXbins()->GetArray());
    dNdetaRef->Sumw2();
    outlist->Add(dNdetaRef);
  }
  TH2D* dNdetaDiff = (TH2D*)outlist->FindObject(Form("%s%s_dNdetaDiff", fType.Data(), GetQCType(fFlags)));
  if (!dNdetaDiff) {
    dNdetaDiff = new TH2D(Form("%s%s_dNdetaDiff", fType.Data(), GetQCType(fFlags)), 
                                Form("%s%s_dNdetaDiff", fType.Data(), GetQCType(fFlags)), 
                fCumuNUADiff->GetNbinsX(), fCumuNUADiff->GetXaxis()->GetXmin(), fCumuNUADiff->GetXaxis()->GetXmax(),
                fCumuNUADiff->GetNbinsY(), fCumuNUADiff->GetYaxis()->GetXmin(), fCumuNUADiff->GetYaxis()->GetXmax());
    dNdetaDiff->GetYaxis()->Set(fCumuNUADiff->GetNbinsY(), fCumuNUADiff->GetYaxis()->GetXbins()->GetArray());
    dNdetaDiff->Sumw2();
    outlist->Add(dNdetaDiff);
  }
  
  // Setting up outputs
  // Create output lists and diagnostics
  TList* vtxList = (TList*)outlist->FindObject("vtxList");
  if (!vtxList) {
    vtxList = new TList();
    vtxList->SetName("vtxList");
    outlist->Add(vtxList);
  }
  vtxList->Add(fCumuNUARef);
  vtxList->Add(fCumuNUADiff);
 
  // Setup output profiles
  CumuHistos cumu2(fMaxMoment, ((fFlags & kNUAcorr) ? 2 : 0));
  CumuHistos cumu4(fMaxMoment, ((fFlags & kNUAcorr) ? 1 : 0));

  cumu2.ConnectList(Form("%sQC2_Cumu%s_vtx_%d_%d", fType.Data(), GetQCType(fFlags), fVzMin, fVzMax), vtxList);
  if ((fFlags & kStdQC)) 
    cumu4.ConnectList(Form("%sQC4_Cumu%s_vtx_%d_%d", fType.Data(), GetQCType(fFlags), fVzMin, fVzMax), vtxList);

  for (Int_t n = 2; n <= fMaxMoment; n++) {
    // 2-particle 
    cumu2.Add(MakeOutputHist(2, n, "Ref", CumuHistos::kNoNUA));
    if ((fFlags & k3Cor)){
      cumu2.Add(MakeOutputHist(2, n, "DiffA", CumuHistos::kNoNUA));
      cumu2.Add(MakeOutputHist(2, n, "DiffB", CumuHistos::kNoNUA));
    } else {
      cumu2.Add(MakeOutputHist(2, n, "Diff", CumuHistos::kNoNUA));
    }
    // 4-particle      
    if ((fFlags & kStdQC)) {
      cumu4.Add(MakeOutputHist(4, n, "Ref", CumuHistos::kNoNUA));
      cumu4.Add(MakeOutputHist(4, n, "Diff", CumuHistos::kNoNUA));
    }
  } // End of v_n result loop
  // NUA corrected
  if ((fFlags & kNUAcorr)) {
    for (Int_t n = 2; n <= fMaxMoment; n++) {
      // 2-particle 
      cumu2.Add(MakeOutputHist(2, n, "Ref", CumuHistos::kNUAOld));
      if ((fFlags & k3Cor)) {
        cumu2.Add(MakeOutputHist(2, n, "DiffA", CumuHistos::kNUAOld));
        cumu2.Add(MakeOutputHist(2, n, "DiffB", CumuHistos::kNUAOld));
      } else {
	cumu2.Add(MakeOutputHist(2, n, "Diff", CumuHistos::kNUAOld));
      }
      // 4-particle      
      if ((fFlags & kStdQC)) {
	cumu4.Add(MakeOutputHist(4, n, "Ref", CumuHistos::kNUAOld));
	cumu4.Add(MakeOutputHist(4, n, "Diff", CumuHistos::kNUAOld));
      }   
    }
    for (Int_t n = 2; n <= fMaxMoment; n++) {
      // 2-particle 
      cumu2.Add(MakeOutputHist(2, n, "Ref", CumuHistos::kNUA));
      if ((fFlags & k3Cor)) {
	cumu2.Add(MakeOutputHist(2, n, "DiffA", CumuHistos::kNUA));
	cumu2.Add(MakeOutputHist(2, n, "DiffB", CumuHistos::kNUA));
      } else {
	cumu2.Add(MakeOutputHist(2, n, "Diff", CumuHistos::kNUA));
      }
    }
  }

  // Calculating the cumulants
  if ((fFlags & k3Cor)) {
    Calculate3CorFlow(cumu2, quality, cent, dNdetaRef, dNdetaDiff);
  } else {
    CalculateReferenceFlow(cumu2, cumu4, quality, cent, dNdetaRef);
    CalculateDifferentialFlow(cumu2, cumu4, quality, dNdetaDiff);
  }
  if ((fFlags & kNUAcorr)) {
    SolveCoupledFlowEquations(cumu2, 'r');
    if ((fFlags & k3Cor)) {
      SolveCoupledFlowEquations(cumu2, 'a');
      SolveCoupledFlowEquations(cumu2, 'b');
    } else {
      SolveCoupledFlowEquations(cumu2, 'd');
    }
  }
 
  // Add to output for immediate viewing - individual vtx bins are used for final results
  AddVertexBins(cumu2, outlist, ((fFlags & kNUAcorr) ? 2 : 0));
  if ((fFlags & kStdQC)) AddVertexBins(cumu4, outlist, ((fFlags & kNUAcorr) ? 1 : 0));

  // Setup NUA diagnoastics histograms
  TList* nualist = (TList*)outlist->FindObject("NUATerms");
  if (!nualist) {
    nualist = new TList();
    nualist->SetName("NUATerms");
    outlist->Add(nualist);
  }
  // Reference
  TH2D* nuaRef = (TH2D*)nualist->FindObject(Form("%sReferenceNUA%s", fType.Data(), GetQCType(fFlags)));
  TH2D* temp = 0;
  if (!nuaRef) {
    nuaRef = (TH2D*)fCumuNUARef->Project3D("yz");
    nuaRef->Scale(1./fCumuNUARef->GetNbinsX());
    nuaRef->SetName(Form("%sReferenceNUA%s", fType.Data(), GetQCType(fFlags)));
    nuaRef->SetTitle(Form("%sReferenceNUA%s", fType.Data(), GetQCType(fFlags)));
    nualist->Add(nuaRef);
  } else {
    temp = (TH2D*)fCumuNUARef->Project3D("yz");
    temp->Scale(1./fCumuNUARef->GetNbinsX());
    nuaRef->Add(temp);
    delete temp;
  }
  // Filling in underflow to make scaling possible in Terminate()
  nuaRef->Fill(0., -1., 1.);
  // Differential
  TH2D* nuaDiff = (TH2D*)nualist->FindObject(Form("%sDifferentialNUA%s", fType.Data(), GetQCType(fFlags)));
  if (!nuaDiff) {
    nuaDiff = (TH2D*)fCumuNUADiff->Project3D("yz");
    nuaDiff->SetName(Form("%sDifferentialNUA%s", fType.Data(), GetQCType(fFlags)));
    nuaDiff->SetTitle(Form("%sDifferentialNUA%s", fType.Data(), GetQCType(fFlags)));
    nualist->Add(nuaDiff);
  } else {
    temp = (TH2D*)fCumuNUADiff->Project3D("yz");
    nuaDiff->Add(temp);
    delete temp;
  }
  // Filling in underflow to make scaling possible in Terminate()
  nuaDiff->Fill(0., -1., 1.);

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::CalculateReferenceFlow(CumuHistos& cumu2h, CumuHistos& cumu4h, TH2I* quality, 
                                                             TH1D* chist, TH2D* dNdetaRef) const  
{
  // 
  //  Calculates the reference flow
  //
  //  Parameters:
  //   cumu2h: CumuHistos object with QC{2} cumulants
  //   cumu4h: CumuHistos object with QC{4} cumulants
  //   quality: Histogram for success rate of cumulants
  //   chist: Centrality histogram
  //   dNdetaRef: dN/deta histogram for estimating multiplicity used for ref calculations
  //

  // Normalizing common NUA hists
  for (Int_t cBin = 1; cBin <= fCumuNUARef->GetNbinsY(); cBin++) {
    Double_t cent = fCumuNUARef->GetYaxis()->GetBinCenter(cBin);
    for (Int_t eBin = 1; eBin <= fCumuNUARef->GetNbinsX(); eBin++) {
      Double_t eta = fCumuNUARef->GetXaxis()->GetBinCenter(eBin);
      Double_t mult = fCumuNUARef->GetBinContent(eBin, cBin, 0);
      if (mult == 0) continue;
      for (Int_t qBin = 1; qBin <= fCumuNUARef->GetNbinsZ(); qBin++) {
	fCumuNUARef->SetBinContent(eBin, cBin, qBin, fCumuNUARef->GetBinContent(eBin, cBin, qBin)/mult);
	fCumuNUARef->SetBinError(eBin, cBin, qBin, fCumuNUARef->GetBinError(eBin, cBin, qBin)/mult);
      }
      // Fill dN/deta diagnostics
      dNdetaRef->Fill(eta, cent, mult);
    }
  }

  // For flow calculations
  TH3D* cumuRef = 0; 
  TH2D* cumu2 = 0;
  TH2D* cumu2NUAold = 0;
  TH2D* cumu4 = 0;
  TH2D* cumu4NUA = 0;
  Int_t qualityFactor = ((fFlags & kStdQC) ? 8 : 4); // used for correctly filling in quality hists
  // Loop over cumulant histogram for final calculations 
  for (Int_t n = 2; n <= fMaxMoment; n++) { // Moment loop begins
    cumu2 = (TH2D*)cumu2h.Get('r', n, CumuHistos::kNoNUA);
    if ((fFlags & kNUAcorr)) {
      cumu2NUAold = (TH2D*)cumu2h.Get('r', n, CumuHistos::kNUAOld);
    }
    if ((fFlags & kStdQC)) {
      cumu4 = (TH2D*)cumu4h.Get('r', n, CumuHistos::kNoNUA);
      if ((fFlags & kNUAcorr)) cumu4NUA = (TH2D*)cumu4h.Get('r', n, CumuHistos::kNUAOld);
    }
    cumuRef  = (TH3D*)fCumuHists.Get('r', n);
    // Begin loops
    for (Int_t cBin = 1; cBin <= cumuRef->GetNbinsY(); cBin++) { // Centrality loop begins
      Double_t cent = cumuRef->GetYaxis()->GetBinCenter(cBin);
      if (n == 2) chist->Fill(cent, cumuRef->GetBinContent(0, cBin, 0));
      if (fDebug > 0) AliInfo(Form("%s - v_%d: centrality %3.1f:..", fType.Data(), n, cent));
      for (Int_t etaBin = 1; etaBin <= cumuRef->GetNbinsX(); etaBin++) { // Eta loop begins
	Double_t eta = cumuRef->GetXaxis()->GetBinCenter(etaBin);
	Double_t refEta = eta;
	Int_t refEtaBinA = fCumuNUARef->GetXaxis()->FindBin(refEta);
	if ((fFlags & kEtaGap)) refEta = -eta;
	Int_t refEtaBinB = fCumuNUARef->GetXaxis()->FindBin(refEta);
	// 2-particle reference flow
	Double_t w2Two = cumuRef->GetBinContent(refEtaBinA, cBin, kW2Two);
	Double_t w2 = cumuRef->GetBinContent(refEtaBinA, cBin, kW2);
	if (w2 == 0) continue;
	Double_t cosP1nPhiA = fCumuNUARef->GetBinContent(refEtaBinA, cBin, GetBinNumberCos(n));
	Double_t sinP1nPhiA = fCumuNUARef->GetBinContent(refEtaBinA, cBin, GetBinNumberSin(n));
	Double_t cosP1nPhiB = fCumuNUARef->GetBinContent(refEtaBinB, cBin, GetBinNumberCos(n));
	Double_t sinP1nPhiB = fCumuNUARef->GetBinContent(refEtaBinB, cBin, GetBinNumberSin(n));
	Double_t cos2nPhiA = fCumuNUARef->GetBinContent(refEtaBinA, cBin, GetBinNumberCos(2*n));
	Double_t sin2nPhiA = fCumuNUARef->GetBinContent(refEtaBinA, cBin, GetBinNumberSin(2*n));
	Double_t cos2nPhiB = fCumuNUARef->GetBinContent(refEtaBinB, cBin, GetBinNumberCos(2*n));
	Double_t sin2nPhiB = fCumuNUARef->GetBinContent(refEtaBinB, cBin, GetBinNumberSin(2*n));  
	Double_t two = w2Two / w2; 
	Double_t qc2 = two;
        if (qc2 >= 0) cumu2->Fill(eta, cent, TMath::Sqrt(qc2));

	if ((fFlags & kNUAcorr)) {
	  // Old NUA
	  // With no eta gap the last two terms are <<cos(phi)>>^2 and <<sin(phi)>>^2,
	  // with eta gap the different coverage is taken into account. 
	  // The next line covers both cases.
	  qc2 -= cosP1nPhiA*cosP1nPhiB + sinP1nPhiA*sinP1nPhiB;
	  // Extra NUA term from 2n cosines and sines
	  Double_t den = 1+(cos2nPhiA*cos2nPhiB + sin2nPhiA*sin2nPhiB);
	  if (den != 0) qc2 /= den;
	  else qc2 = 0;
	}
	if (qc2 <= 0) { 
	  if (fDebug > 0) 
	    AliInfo(Form("%s: QC_%d{2} = %1.3f for eta = %1.2f and centrality %3.1f - skipping", 
	    fType.Data(), n, qc2, eta, cent));
   	  quality->Fill((n-2)*qualityFactor+2, Int_t(cent));
	  continue;
	}
	Double_t vnTwo = TMath::Sqrt(qc2);
	if (!TMath::IsNaN(vnTwo)) { 
    	  quality->Fill((n-2)*qualityFactor+1, Int_t(cent));
	  if ((fFlags & kNUAcorr)) cumu2NUAold->Fill(eta, cent, vnTwo);
	}

	if (!(fFlags & kStdQC)) continue;
	// 4-particle reference flow
	Double_t w4Four = cumuRef->GetBinContent(refEtaBinA, cBin, kW4Four);
	Double_t w4 = cumuRef->GetBinContent(refEtaBinA, cBin, kW4);
	Double_t multm1m2 = cumuRef->GetBinContent(refEtaBinA, cBin, k3pWeight);
	if (w4 == 0 || multm1m2 == 0) continue;
	Double_t cosP1nPhi1P1nPhi2 = cumuRef->GetBinContent(refEtaBinA, cBin, kCosphi1phi2);
	Double_t sinP1nPhi1P1nPhi2 = cumuRef->GetBinContent(refEtaBinA, cBin, kSinphi1phi2);
	Double_t cosP1nPhi1M1nPhi2M1nPhi3 = cumuRef->GetBinContent(refEtaBinA, cBin, kCosphi1phi2phi3m);
	Double_t sinP1nPhi1M1nPhi2M1nPhi3 = cumuRef->GetBinContent(refEtaBinA, cBin, kSinphi1phi2phi3m);

	cosP1nPhi1P1nPhi2 /= w2;
	sinP1nPhi1P1nPhi2 /= w2;
	cosP1nPhi1M1nPhi2M1nPhi3 /= multm1m2;
	sinP1nPhi1M1nPhi2M1nPhi3 /= multm1m2;
	Double_t four = w4Four / w4;
	Double_t qc4 = four-2.*TMath::Power(two,2.);
        if (qc4 < 0) cumu4->Fill(eta, cent, TMath::Power(-qc4, 0.25));
	
	if ((fFlags & kNUAcorr)) {
	   qc4 += - 4.*cosP1nPhiA*cosP1nPhi1M1nPhi2M1nPhi3
		  + 4.*sinP1nPhiA*sinP1nPhi1M1nPhi2M1nPhi3-TMath::Power(cosP1nPhi1P1nPhi2,2.)-TMath::Power(sinP1nPhi1P1nPhi2,2.)
		  + 4.*cosP1nPhi1P1nPhi2*(TMath::Power(cosP1nPhiA,2.)-TMath::Power(sinP1nPhiA,2.))
		  + 8.*sinP1nPhi1P1nPhi2*sinP1nPhiA*cosP1nPhiA
		  + 8.*two*(TMath::Power(cosP1nPhiA,2.)+TMath::Power(sinP1nPhiA,2.))
		  - 6.*TMath::Power((TMath::Power(cosP1nPhiA,2.)+TMath::Power(sinP1nPhiA,2.)),2.);
	}
	if (qc4 >= 0) {
	  if (fDebug > 0) 
	    AliInfo(Form("%s: QC_%d{4} = %1.3f for eta = %1.2f and centrality %3.1f - skipping", 
	    fType.Data(), n, qc2, eta, cent));
	  quality->Fill((n-2)*qualityFactor+6, Int_t(cent));
	  continue;
	}
	Double_t vnFour = TMath::Power(-qc4, 0.25);
	if (!TMath::IsNaN(vnFour*multm1m2)) {
	  quality->Fill((n-2)*qualityFactor+5, Int_t(cent));
	  if ((fFlags & kNUAcorr)) cumu4NUA->Fill(eta, cent, vnFour);
	}
      } // End of eta
    } // End of cent
  } // End of moment

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::CalculateDifferentialFlow(CumuHistos& cumu2h, CumuHistos& cumu4h, 
                                                                TH2I* quality, TH2D* dNdetaDiff) const  
{
  // 
  //  Calculates the differential flow
  //
  //  Parameters:
  //   cumu2h: CumuHistos object with QC{2} cumulants
  //   cumu4h: CumuHistos object with QC{4} cumulants
  //   quality: Histogram for success rate of cumulants
  //   dNdetaDiff: dN/deta histogram for estimating multiplicity used for diff calculations
  //

  for (Int_t cBin = 1; cBin <= fCumuNUADiff->GetNbinsY(); cBin++) {
    Double_t cent = fCumuNUADiff->GetYaxis()->GetBinCenter(cBin);
    for (Int_t eBin = 1; eBin <= fCumuNUADiff->GetNbinsX(); eBin++) {
      Double_t eta = fCumuNUADiff->GetXaxis()->GetBinCenter(eBin);
      Double_t mult = fCumuNUADiff->GetBinContent(eBin, cBin, 0);
      if (mult == 0) continue;
      for (Int_t qBin = 1; qBin <= fCumuNUADiff->GetNbinsZ(); qBin++) {
	fCumuNUADiff->SetBinContent(eBin, cBin, qBin, fCumuNUADiff->GetBinContent(eBin, cBin, qBin)/mult);
	fCumuNUADiff->SetBinError(eBin, cBin, qBin, fCumuNUADiff->GetBinError(eBin, cBin, qBin)/mult);
      }
      dNdetaDiff->Fill(eta, cent, mult);
    }
  }

  // For flow calculations
  TH3D* cumuRef = 0; 
  TH3D* cumuDiff = 0; 
  TH2D* cumu2 = 0;
  TH2D* cumu2NUAold = 0;
  TH2D* cumu4 = 0;
  TH2D* cumu4NUA = 0;
  Int_t qualityFactor = ((fFlags & kStdQC) ? 8 : 4); // used for correctly filling in quality hists
  // Loop over cumulant histogram for final calculations 
  for (Int_t n = 2; n <= fMaxMoment; n++) { // Moment loop begins
    cumu2 = (TH2D*)cumu2h.Get('d', n, CumuHistos::kNoNUA);
    if ((fFlags & kNUAcorr)) {
      cumu2NUAold = (TH2D*)cumu2h.Get('d', n, CumuHistos::kNUAOld);
    }
    if ((fFlags & kStdQC)) {
      cumu4 = (TH2D*)cumu4h.Get('d',n);
      if ((fFlags & kNUAcorr)) cumu4NUA = (TH2D*)cumu4h.Get('d', n, CumuHistos::kNUAOld);
    }
    cumuRef  = (TH3D*)fCumuHists.Get('r',n);
    cumuDiff = (TH3D*)fCumuHists.Get('d',n);
    for (Int_t cBin = 1; cBin <= cumuDiff->GetNbinsY(); cBin++) { // Centrality loop begins
      Double_t cent = cumuDiff->GetYaxis()->GetBinCenter(cBin);
      if (fDebug > 0) AliInfo(Form("%s - v_%d: centrality %3.1f:..", fType.Data(), n, cent));
      for (Int_t etaBin = 1; etaBin <= cumuDiff->GetNbinsX(); etaBin++) { // Eta loop begins
	Double_t eta = cumuDiff->GetXaxis()->GetBinCenter(etaBin);
	Double_t refEta = eta;
	Int_t refEtaBinA = fCumuNUARef->GetXaxis()->FindBin(refEta);
	if ((fFlags & kEtaGap)) refEta = -eta;
	Int_t refEtaBinB = fCumuNUARef->GetXaxis()->FindBin(refEta);

        // Reference objects
	Double_t w2 = cumuRef->GetBinContent(refEtaBinA, cBin, kW2);
	if (w2 == 0) continue;
	Double_t two = cumuRef->GetBinContent(refEtaBinA, cBin, kW2Two);
        two /= w2;
 	Double_t cosP1nPhiA = fCumuNUARef->GetBinContent(refEtaBinA, cBin, GetBinNumberCos(n));
	Double_t sinP1nPhiA = fCumuNUARef->GetBinContent(refEtaBinA, cBin, GetBinNumberSin(n));
	Double_t cosP1nPhiB = fCumuNUARef->GetBinContent(refEtaBinB, cBin, GetBinNumberCos(n));
	Double_t sinP1nPhiB = fCumuNUARef->GetBinContent(refEtaBinB, cBin, GetBinNumberSin(n));       
	Double_t cos2nPhiB = fCumuNUARef->GetBinContent(refEtaBinB, cBin, GetBinNumberCos(2*n));
	Double_t sin2nPhiB = fCumuNUARef->GetBinContent(refEtaBinB, cBin, GetBinNumberSin(2*n));  
	
	// 2-particle differential flow
	Double_t w2pTwoPrime = cumuDiff->GetBinContent(etaBin, cBin, kW2Two);
	Double_t w2p = cumuDiff->GetBinContent(etaBin, cBin, kW2);
	if (w2p == 0) continue;
	Double_t cosP1nPsi = fCumuNUADiff->GetBinContent(etaBin, cBin, GetBinNumberCos(n));
	Double_t sinP1nPsi = fCumuNUADiff->GetBinContent(etaBin, cBin, GetBinNumberSin(n));
	Double_t cos2nPsi = fCumuNUADiff->GetBinContent(etaBin, cBin, GetBinNumberCos(2*n));
	Double_t sin2nPsi = fCumuNUADiff->GetBinContent(etaBin, cBin, GetBinNumberSin(2*n));
	Double_t twoPrime = w2pTwoPrime / w2p;

	Double_t qc2Prime = twoPrime;
	cumu2->Fill(eta, cent, qc2Prime);
	if ((fFlags & kNUAcorr)) {
	  // Old nua
	  qc2Prime -= cosP1nPsi*cosP1nPhiB + sinP1nPsi*sinP1nPhiB;
	  // Extra NUA term from 2n cosines and sines
	  qc2Prime /= (1.+(cos2nPsi*cos2nPhiB + sin2nPsi*sin2nPhiB));
	}
	if (!TMath::IsNaN(qc2Prime)) {
	  quality->Fill((n-2)*qualityFactor+3, Int_t(cent));
	  if ((fFlags & kNUAcorr)) cumu2NUAold->Fill(eta, cent, qc2Prime);
	}
	else 
	  quality->Fill((n-2)*qualityFactor+4, Int_t(cent));
	if (fDebug > 1) 
	  AliInfo(Form("%s: QC'_%d{2} = %1.3f for eta = %1.2f and centrality %3.1f", 
	  fType.Data(), n, qc2Prime, eta, cent));

        if (!(fFlags & kStdQC)) continue;
        // Reference objects
	Double_t cosP1nPhi1P1nPhi2 = cumuRef->GetBinContent(refEtaBinA, cBin, kCosphi1phi2);
	Double_t sinP1nPhi1P1nPhi2 = cumuRef->GetBinContent(refEtaBinA, cBin, kSinphi1phi2);
	Double_t cosP1nPhi1M1nPhi2M1nPhi3 = cumuRef->GetBinContent(refEtaBinA, cBin, kCosphi1phi2phi3m);
	Double_t sinP1nPhi1M1nPhi2M1nPhi3 = cumuRef->GetBinContent(refEtaBinA, cBin, kSinphi1phi2phi3m);
	Double_t multm1m2 = cumuRef->GetBinContent(refEtaBinA, cBin, k3pWeight);
	cosP1nPhi1P1nPhi2 /= w2;
	sinP1nPhi1P1nPhi2 /= w2;
	cosP1nPhi1M1nPhi2M1nPhi3 /= multm1m2;
	sinP1nPhi1M1nPhi2M1nPhi3 /= multm1m2;

	// 4-particle differential flow
	Double_t w4pFourPrime = cumuDiff->GetBinContent(etaBin, cBin, kW4Four);
	Double_t w4p = cumuDiff->GetBinContent(etaBin, cBin, kW4);
	Double_t mpqMult = cumuDiff->GetBinContent(etaBin, cBin, k3pWeight);
	if (w4p == 0 || mpqMult == 0) continue;
	Double_t cosP1nPsi1P1nPhi2 = cumuDiff->GetBinContent(etaBin, cBin, kCosphi1phi2);
	Double_t sinP1nPsi1P1nPhi2 = cumuDiff->GetBinContent(etaBin, cBin, kSinphi1phi2);
	Double_t cosP1nPsi1M1nPhi2M1nPhi3 = cumuDiff->GetBinContent(etaBin, cBin, kCosphi1phi2phi3m);
	Double_t sinP1nPsi1M1nPhi2M1nPhi3 = cumuDiff->GetBinContent(etaBin, cBin, kSinphi1phi2phi3m);
	Double_t cosP1nPsi1P1nPhi2M1nPhi3 = cumuDiff->GetBinContent(etaBin, cBin, kCosphi1phi2phi3p);
	Double_t sinP1nPsi1P1nPhi2M1nPhi3 = cumuDiff->GetBinContent(etaBin, cBin, kSinphi1phi2phi3p); 
	
	cosP1nPsi1P1nPhi2 /= w2p;
	sinP1nPsi1P1nPhi2 /= w2p;
	cosP1nPsi1M1nPhi2M1nPhi3 /= mpqMult;
	sinP1nPsi1M1nPhi2M1nPhi3 /= mpqMult;
	cosP1nPsi1P1nPhi2M1nPhi3 /= mpqMult;
	sinP1nPsi1P1nPhi2M1nPhi3 /= mpqMult;
       
	Double_t fourPrime = w4pFourPrime / w4p;
	Double_t qc4Prime = fourPrime-2.*twoPrime*two; 
	if (cumu4) cumu4->Fill(eta, cent, qc4Prime);

	if ((fFlags & kNUAcorr)) {
	  qc4Prime += - cosP1nPsi*cosP1nPhi1M1nPhi2M1nPhi3
		      + sinP1nPsi*sinP1nPhi1M1nPhi2M1nPhi3
		      - cosP1nPhiA*cosP1nPsi1M1nPhi2M1nPhi3
		      + sinP1nPhiA*sinP1nPsi1M1nPhi2M1nPhi3
		      - 2.*cosP1nPhiA*cosP1nPsi1P1nPhi2M1nPhi3
		      - 2.*sinP1nPhiA*sinP1nPsi1P1nPhi2M1nPhi3
		      - cosP1nPsi1P1nPhi2*cosP1nPhi1P1nPhi2
		      - sinP1nPsi1P1nPhi2*sinP1nPhi1P1nPhi2
		      + 2.*cosP1nPhi1P1nPhi2*(cosP1nPsi*cosP1nPhiA-sinP1nPsi*sinP1nPhiA)
		      + 2.*sinP1nPhi1P1nPhi2*(cosP1nPsi*sinP1nPhiA+sinP1nPsi*cosP1nPhiA)
		      + 4.*two*(cosP1nPsi*cosP1nPhiA+sinP1nPsi*sinP1nPhiA) 
		      + 2.*cosP1nPsi1P1nPhi2*(TMath::Power(cosP1nPhiA,2.)-TMath::Power(sinP1nPhiA,2.))
		      + 4.*sinP1nPsi1P1nPhi2*cosP1nPhiA*sinP1nPhiA
		      + 4.*twoPrime*(TMath::Power(cosP1nPhiA,2.)+TMath::Power(sinP1nPhiA,2.))
		      - 6.*(TMath::Power(cosP1nPhiA,2.)-TMath::Power(sinP1nPhiA,2.)) 
		      * (cosP1nPsi*cosP1nPhiA-sinP1nPsi*sinP1nPhiA)
		      - 12.*cosP1nPhiA*sinP1nPhiA
		      * (sinP1nPsi*cosP1nPhiA+cosP1nPsi*sinP1nPhiA);
	}
//	Double_t vnFourDiff = - qc4Prime / TMath::Power(-qc4, 0.75);
	if (!TMath::IsNaN(qc4Prime*mpqMult)) {
	  quality->Fill((n-2)*qualityFactor+7, Int_t(cent));
	  if (cumu4NUA) cumu4NUA->Fill(eta, cent, qc4Prime);
	}
	else 
	  quality->Fill((n-2)*qualityFactor+8, Int_t(cent));
	if (fDebug > 1) 
	  AliInfo(Form("%s: v_%d{4} = %1.3f for eta = %1.2f and centrality %3.1f", 
	  fType.Data(), n, qc4Prime, eta, cent));
      } // End of eta loop
    } // End of centrality loop
  } // End of moment

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::Calculate3CorFlow(CumuHistos& cumu2h, TH2I* quality, TH1D* chist,
                                                        TH2D* dNdetaRef, TH2D* dNdetaDiff) const  
{
  // 
  //  Calculates the 3 sub flow
  //
  //  Parameters:
  //   cumu2h: CumuHistos object with QC{2} cumulants
  //   quality: Histogram for success rate of cumulants
  //   chist: Centrality histogram
  //   dNdetaDiff: dN/deta histogram for estimating multiplicity used for diff calculations
  //

  // For flow calculations
  TH3D* cumuRef = 0; 
  TH3D* cumuDiff = 0; 
  TH2D* cumu2r = 0;
  TH2D* cumu2rNUAold = 0;
  TH2D* cumu2a = 0;
  TH2D* cumu2aNUAold = 0;
  TH2D* cumu2b = 0;
  TH2D* cumu2bNUAold = 0;
  // Loop over cumulant histogram for final calculations 
  for (Int_t n = 2; n <= fMaxMoment; n++) { // Moment loop begins
    cumu2r = (TH2D*)cumu2h.Get('r', n, CumuHistos::kNoNUA);
    cumu2a = (TH2D*)cumu2h.Get('a', n, CumuHistos::kNoNUA);
    cumu2b = (TH2D*)cumu2h.Get('b', n, CumuHistos::kNoNUA);
    if ((fFlags & kNUAcorr)) {
      cumu2rNUAold = (TH2D*)cumu2h.Get('r', n, CumuHistos::kNUAOld);
      cumu2aNUAold = (TH2D*)cumu2h.Get('a', n, CumuHistos::kNUAOld);
      cumu2bNUAold = (TH2D*)cumu2h.Get('b', n, CumuHistos::kNUAOld);
    }
    cumuRef  = (TH3D*)fCumuHists.Get('r',n);
    cumuDiff = (TH3D*)fCumuHists.Get('d',n);
    for (Int_t cBin = 1; cBin <= cumuRef->GetNbinsY(); cBin++) { // Centrality loop begins
      Double_t cent = cumuRef->GetYaxis()->GetBinCenter(cBin);
      if (n == 2) chist->Fill(cent, cumuRef->GetBinContent(0, cBin, 0));
      if (fDebug > 0) AliInfo(Form("%s - v_%d: centrality %3.1f:..", fType.Data(), n, cent));
      // Here it starts!
      Int_t prevLim = 0;
      Int_t aLow = 0, aHigh = 0, bLow = 0, bHigh = 0;
      Double_t cosP1nPhiA = 0;
      Double_t sinP1nPhiA = 0;
      Double_t cos2nPhiA = 0;
      Double_t sin2nPhiA = 0;
      Double_t cosP1nPhiB = 0;
      Double_t sinP1nPhiB = 0;
      Double_t cos2nPhiB = 0;
      Double_t sin2nPhiB = 0;
      Double_t multA = 0;
      Double_t multB = 0;

      for (Int_t etaBin = 1; etaBin <= cumuDiff->GetNbinsX(); etaBin++) { // Eta loop begins
	Double_t eta = cumuDiff->GetXaxis()->GetBinCenter(etaBin);
	// 2-particle reference flow
	Double_t w2Two = cumuRef->GetBinContent(etaBin, cBin, kW2Two);
	Double_t w2 = cumuRef->GetBinContent(etaBin, cBin, kW2);
	if (w2 == 0) continue;

	// Update NUA for new range!
	if (fEtaLims[prevLim] < eta) {
	  GetLimits(prevLim, aLow, aHigh, bLow, bHigh);
	  prevLim++;
	  cosP1nPhiA = 0; sinP1nPhiA = 0; cos2nPhiA = 0; sin2nPhiA = 0; multA = 0;
	  cosP1nPhiB = 0; sinP1nPhiB = 0; cos2nPhiB = 0; sin2nPhiB = 0; multB = 0;
	  for (Int_t a = aLow; a <= aHigh; a++) {
	    cosP1nPhiA += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberCos(n));
	    sinP1nPhiA += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberSin(n));
	    cos2nPhiA += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberCos(2*n));
	    sin2nPhiA += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberSin(2*n));
	    multA += fCumuNUARef->GetBinContent(a, cBin, 0);
	  }
	  for (Int_t b = bLow; b <= bHigh; b++) {
	    cosP1nPhiB += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberCos(n));
	    sinP1nPhiB += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberSin(n));
	    cos2nPhiB += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberCos(2*n));
	    sin2nPhiB += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberSin(2*n));
	    multB += fCumuNUARef->GetBinContent(b, cBin, 0);
	  }
	  if (multA == 0 || multB == 0) {
	    AliWarning(Form("Empty NUA values for 3Cor! (%s)", cumuRef->GetName()));
	    continue;
	  }
	  cosP1nPhiA /= multA;
	  sinP1nPhiA /= multA;
	  cos2nPhiA /= multA;
	  sin2nPhiA /= multA;
	  cosP1nPhiB /= multB;
	  sinP1nPhiB /= multB;
	  cos2nPhiB /= multB;
	  sin2nPhiB /= multB;

	  dNdetaRef->Fill(eta, cent, multA+multB);
	}
	Double_t two = w2Two / w2;
  
	Double_t qc2 = two;
        if (qc2 >= 0) cumu2r->Fill(eta, cent, TMath::Sqrt(qc2));

	if ((fFlags & kNUAcorr)) {
	  // Old nua
	  qc2 -= cosP1nPhiA*cosP1nPhiB + sinP1nPhiA*sinP1nPhiB;
	  // Extra NUA term from 2n cosines and sines
	  qc2 /= (1+(cos2nPhiA*cos2nPhiB + sin2nPhiA*sin2nPhiB));
	}
	if (qc2 <= 0) { 
	  if (fDebug > 0) 
	    AliInfo(Form("%s: QC_%d{2} = %1.3f for eta = %1.2f and centrality %3.1f - skipping", 
	    fType.Data(), n, qc2, eta, cent));
   	  quality->Fill((n-2)*4+2, Int_t(cent));
	  continue;
	}
	Double_t vnTwo = TMath::Sqrt(qc2);
	if (!TMath::IsNaN(vnTwo)) { 
    	  quality->Fill((n-2)*4+1, Int_t(cent));
	  if ((fFlags & kNUAcorr)) cumu2rNUAold->Fill(eta, cent, vnTwo);
	}

	// 2-particle differential flow
	Double_t w2pTwoPrimeA = cumuDiff->GetBinContent(etaBin, cBin, kW2Two);
	Double_t w2pA = cumuDiff->GetBinContent(etaBin, cBin, kW2);
	Double_t w2pTwoPrimeB = cumuDiff->GetBinContent(etaBin, cBin, kW4Four);
	Double_t w2pB = cumuDiff->GetBinContent(etaBin, cBin, kW4);
	if (w2pA == 0 || w2pB == 0) continue;
	Double_t cosP1nPsi = fCumuNUADiff->GetBinContent(etaBin, cBin, GetBinNumberCos(n));
	Double_t sinP1nPsi = fCumuNUADiff->GetBinContent(etaBin, cBin, GetBinNumberSin(n));
	Double_t cos2nPsi = fCumuNUADiff->GetBinContent(etaBin, cBin, GetBinNumberCos(2*n));
	Double_t sin2nPsi = fCumuNUADiff->GetBinContent(etaBin, cBin, GetBinNumberSin(2*n));
        Double_t mult = fCumuNUADiff->GetBinContent(etaBin, cBin, 0);
        if (mult == 0) continue;
        cosP1nPsi /= mult;
        sinP1nPsi /= mult;
        cos2nPsi /= mult;
        sin2nPsi /= mult;
	Double_t twoPrimeA = w2pTwoPrimeA / w2pA;
	Double_t twoPrimeB = w2pTwoPrimeB / w2pB;
	dNdetaDiff->Fill(eta, cent, mult);

	Double_t qc2PrimeA = twoPrimeA;
	Double_t qc2PrimeB = twoPrimeB;
	if (qc2PrimeA*qc2PrimeB >= 0) {
	  cumu2a->Fill(eta, cent, qc2PrimeA);
	  cumu2b->Fill(eta, cent, qc2PrimeB);
	}
	if ((fFlags & kNUAcorr)) {
	  // Old nua
	  qc2PrimeA -= cosP1nPsi*cosP1nPhiA + sinP1nPsi*sinP1nPhiA;
	  qc2PrimeB -= cosP1nPsi*cosP1nPhiB + sinP1nPsi*sinP1nPhiB; // Is this OK?
	  // Extra NUA term from 2n cosines and sines
	  if (cos2nPsi*cos2nPhiA + sin2nPsi*sin2nPhiA != -1.) qc2PrimeA /= (1.+(cos2nPsi*cos2nPhiA + sin2nPsi*sin2nPhiA));
	  if (cos2nPsi*cos2nPhiB + sin2nPsi*sin2nPhiB != -1.) qc2PrimeB /= (1.+(cos2nPsi*cos2nPhiB + sin2nPsi*sin2nPhiB));
	}
	if (!TMath::IsNaN(qc2PrimeA) && !TMath::IsNaN(qc2PrimeB) && qc2 != 0) {
	if (qc2PrimeA*qc2PrimeB >= 0) {
	  quality->Fill((n-2)*4+3, Int_t(cent));
	  if ((fFlags & kNUAcorr)) cumu2aNUAold->Fill(eta, cent, qc2PrimeA);
	  if ((fFlags & kNUAcorr)) cumu2bNUAold->Fill(eta, cent, qc2PrimeB);
	}
      }
      else 
	quality->Fill((n-2)*4+4, Int_t(cent));
	if (fDebug > 1) 
	  AliInfo(Form("%s: QC'a_%d{2} = %1.3f, QC'b_%d{2} = %1.3f for eta = %1.2f and centrality %3.1f", 
	  fType.Data(), n, qc2PrimeA, n, qc2PrimeB, eta, cent));
      } // End of eta loop
    } // End of centrality loop
  } // End of moment

  return;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::SolveCoupledFlowEquations(CumuHistos& cumu, const Char_t type) const  
{
  // 
  //  Function to solve the coupled flow equations
  //  We solve it by using matrix calculations:
  //  A*v_n = V => v_n = A^-1*V
  //  First we set up a TMatrixD container to make ROOT
  //  do the inversions in an efficient way, we multiply the current <<2>> estimates.
  //  Then we fill new TH2D's if the corrected <<2>>'s (cumuNUA object).
  //
  //  Parameters:
  //   cumu: CumuHistos object - uncorrected
  //   type: Reference ('r') or differential ('d') or ('a' or 'b') for 3 correlator
  //

  // We start by declaring Matrix and vector objects, as their constructors are quite heavy
  TMatrixD mNUA(fMaxMoment-1, fMaxMoment-1);
  TVectorD vQC2(fMaxMoment-1);

  for (Int_t cBin = 1; cBin <= cumu.Get(type, 2, CumuHistos::kNUAOld)->GetNbinsY(); cBin++) { // cent loop
    Double_t cent = cumu.Get(type, 2, CumuHistos::kNUAOld)->GetYaxis()->GetBinCenter(cBin);
    for (Int_t eBin = 1; eBin <= cumu.Get(type, 2, CumuHistos::kNUAOld)->GetNbinsX(); eBin++) { // eta loop
      Double_t eta = cumu.Get(type, 2, CumuHistos::kNUAOld)->GetXaxis()->GetBinCenter(eBin);
      mNUA.Zero(); // reset matrix
      vQC2.Zero(); // reset vector
      for (Int_t n = 0; n < fMaxMoment-1; n++) { // moment loop
      	vQC2(n) = static_cast<TH2D*>(cumu.Get(type, n+2, CumuHistos::kNUAOld))->GetBinContent(eBin, cBin);
      	if (type == 'r' || type == 'R') vQC2(n) *= vQC2(n); // going back to <<2>>
      	for (Int_t m = 0; m < fMaxMoment-1; m++) { // cross moment
  	  mNUA(n,m) = CalculateNUAMatrixElement(n, m, type, eBin, cBin);
    	} // End of cross moment loop
      } // End of moment loop
      // Invert matrix
      Double_t det = 0;
      mNUA.Invert(&det);
      // If determinant is non-zero we go with corrected results
      if (det != 0 ) vQC2 = mNUA*vQC2;
      else AliWarning(Form("Determinant == 0 - cent: %d-%d, eta: %f, type: '%c', data: %s, vtx: %d-%d%s", 
		      Int_t(cumu.Get(type, 2, CumuHistos::kNUAOld)->GetYaxis()->GetBinLowEdge(cBin)), 
		      Int_t(cumu.Get(type, 2, CumuHistos::kNUAOld)->GetYaxis()->GetBinUpEdge(cBin)),
		      cumu.Get(type, 2, CumuHistos::kNUAOld)->GetXaxis()->GetBinCenter(eBin), 
		      type, fType.Data(), fVzMin, fVzMax, 
		      GetQCType(fFlags, kTRUE)));
      // Go back to v_n for ref. keep <<2'>> for diff. flow).
      for (Int_t n = 0; n < fMaxMoment-1; n++) {
	Double_t vnTwo = 0;
        if (type == 'r' || type == 'R') 
          vnTwo = (vQC2(n) > 0. ? TMath::Sqrt(vQC2(n)) : 0.);
        else {
          // is really more <<2'>> in this case
	  vnTwo = vQC2(n);
	}
        // Fill in corrected v_n
	if (vnTwo != 0) static_cast<TH2D*>(cumu.Get(type, n+2, CumuHistos::kNUA))->Fill(eta, cent, vnTwo);
      } // End of moment loop
    } // End of eta loop
  } // End of centrality loop
  return;
}
//_____________________________________________________________________
Double_t AliForwardFlowTaskQC::VertexBin::CalculateNUAMatrixElement(Int_t n, Int_t m, Char_t type, Int_t binA, Int_t cBin) const  
{
  //  
  //  Calculates the (n,m)-th element in the NUA matrix: 1 if n == m, otherwise:
  //               <<cos[(n-m)phi1]>>*<<cos[(n-m)phi2]>> + <<sin[(n-m)phi1]>>*<<sin[(n-m)phi2]>>
  //    NUA(n,m) = -----------------------------------------------------------------------------
  //                   1 + <<cos[2nphi1]>>*<<cos[2nphi2]>> + <<sin[2nphi1]>>*<<sin[2nphi2]>>
  //
  //               <<cos[(n+m)phi1]>>*<<cos[(n+m)phi2]>> + <<sin[(n+m)phi1]>>*<<sin[(n+m)phi2]>>
  //             + -----------------------------------------------------------------------------
  //                   1 + <<cos[2nphi1]>>*<<cos[2nphi2]>> + <<sin[2nphi1]>>*<<sin[2nphi2]>>
  //
  //  Parameters:
  //   n: row
  //   m: coumn
  //   type: Reference ('r') or differential ('d') or ('a' or 'b')
  //   binA: eta bin of phi1
  //   cBin: centrality bin
  //
  //  Return: NUA(n,m)
  //
  if (n == m) return 1.;
  n += 2;
  m += 2;

  Double_t cosnMmPhi1 = 0, cosnMmPhi2 = 0, sinnMmPhi1 = 0, sinnMmPhi2 = 0;
  Double_t cosnPmPhi1 = 0, cosnPmPhi2 = 0, sinnPmPhi1 = 0, sinnPmPhi2 = 0;
  Double_t cos2nPhi1 = 0, cos2nPhi2 = 0, sin2nPhi1 = 0, sin2nPhi2 = 0;

  // reference flow
  if (type == 'r' || type == 'R') {
    if ((fFlags & k3Cor)) {
      Double_t eta = fCumuNUARef->GetXaxis()->GetBinCenter(binA);
      Int_t i = 0;
      while (fEtaLims[i] < eta) i++;
      Int_t aLow = 0, aHigh = 0, bLow = 0, bHigh = 0;
      GetLimits(i-1, aLow, aHigh, bLow, bHigh);
      Double_t multA = 0, multB = 0;
      for (Int_t a = aLow; a <= aHigh; a++) {
  	cosnMmPhi1 += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberCos(n-m));
    	sinnMmPhi1 += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberSin(n-m));
      	cosnPmPhi1 += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberCos(n+m));
	sinnPmPhi1 += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberSin(n+m));
	cos2nPhi1 += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberCos(2*n));
	sin2nPhi1 += fCumuNUARef->GetBinContent(a, cBin, GetBinNumberSin(2*n));
	multA += fCumuNUARef->GetBinContent(a, cBin, 0);
      }
       for (Int_t b = bLow; b <= bHigh; b++) {
  	cosnMmPhi2 += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberCos(n-m));
    	sinnMmPhi2 += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberSin(n-m));
      	cosnPmPhi2 += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberCos(n+m));
	sinnPmPhi2 += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberSin(n+m));
	cos2nPhi2 += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberCos(2*n));
	sin2nPhi2 += fCumuNUARef->GetBinContent(b, cBin, GetBinNumberSin(2*n));
	multB += fCumuNUARef->GetBinContent(b, cBin, 0);
      }
      if (multA == 0 || multB == 0) {
        if (fDebug > 0) AliWarning("multA or multB == 0 in matrix elements, aborting NUA");
        return 0.;
      }
      cosnMmPhi1 /= multA;
      sinnMmPhi1 /= multA;
      cosnPmPhi1 /= multA;
      sinnPmPhi1 /= multA;
      cos2nPhi1 /= multA; 
      sin2nPhi1 /= multA;
      cosnMmPhi2 /= multB;
      sinnMmPhi2 /= multB;
      cosnPmPhi2 /= multB;
      sinnPmPhi2 /= multB;
      cos2nPhi2 /= multB; 
      sin2nPhi2 /= multB;
    } else {
      Int_t binB = fCumuNUARef->GetXaxis()->FindBin(-fCumuNUARef->GetXaxis()->GetBinCenter(binA));
      cosnMmPhi1 = fCumuNUARef->GetBinContent(binA, cBin, GetBinNumberCos(n-m));
      sinnMmPhi1 = fCumuNUARef->GetBinContent(binA, cBin, GetBinNumberSin(n-m));
      cosnPmPhi1 = fCumuNUARef->GetBinContent(binA, cBin, GetBinNumberCos(n+m));
      sinnPmPhi1 = fCumuNUARef->GetBinContent(binA, cBin, GetBinNumberSin(n+m));
      cos2nPhi1 = fCumuNUARef->GetBinContent(binA, cBin, GetBinNumberCos(2*n));
      sin2nPhi1 = fCumuNUARef->GetBinContent(binA, cBin, GetBinNumberSin(2*n));
      cosnMmPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberCos(n-m));
      sinnMmPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberSin(n-m));
      cosnPmPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberCos(n+m));
      sinnPmPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberSin(n+m));
      cos2nPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberCos(2*n));
      sin2nPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberSin(2*n));
    }
  } // differential flow
  else if (type == 'd' || type == 'D') {
    Int_t binB = fCumuNUARef->GetXaxis()->FindBin(-fCumuNUADiff->GetXaxis()->GetBinCenter(binA));
    cosnMmPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberCos(n-m));
    sinnMmPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberSin(n-m));
    cosnPmPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberCos(n+m));
    sinnPmPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberSin(n+m));
    cos2nPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberCos(2*n));
    sin2nPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberSin(2*n));
    cosnMmPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberCos(n-m));
    sinnMmPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberSin(n-m));
    cosnPmPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberCos(n+m));
    sinnPmPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberSin(n+m));
    cos2nPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberCos(2*n));
    sin2nPhi2 = fCumuNUARef->GetBinContent(binB, cBin, GetBinNumberSin(2*n));
  } // 3 correlator part a or b
  else if (type == 'a' || type == 'A' || type == 'b' || type == 'B') {
    Double_t mult1 = 0, mult2 = 0;
    // POIs
    cosnMmPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberCos(n-m));
    sinnMmPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberSin(n-m));
    cosnPmPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberCos(n+m));
    sinnPmPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberSin(n+m));
    cos2nPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberCos(2*n));
    sin2nPhi1 = fCumuNUADiff->GetBinContent(binA, cBin, GetBinNumberSin(2*n));
    mult1 = fCumuNUADiff->GetBinContent(binA, cBin, 0);
    // RPs
    Double_t eta = fCumuNUADiff->GetXaxis()->GetBinCenter(binA);
    Int_t i = 0;
    while (fEtaLims[i] < eta) i++;
    Int_t aLow = 0, aHigh = 0, bLow = 0, bHigh = 0;
    GetLimits(i-1, aLow, aHigh, bLow, bHigh);
    Int_t lLow = ((type == 'a' || type == 'A') ? aLow : bLow);
    Int_t lHigh = ((type == 'a' || type == 'A') ? aHigh : bHigh);
    for (Int_t l = lLow; l <= lHigh; l++) {
      cosnMmPhi2 += fCumuNUARef->GetBinContent(l, cBin, GetBinNumberCos(n-m));
      sinnMmPhi2 += fCumuNUARef->GetBinContent(l, cBin, GetBinNumberSin(n-m));
      cosnPmPhi2 += fCumuNUARef->GetBinContent(l, cBin, GetBinNumberCos(n+m));
      sinnPmPhi2 += fCumuNUARef->GetBinContent(l, cBin, GetBinNumberSin(n+m));
      cos2nPhi2 += fCumuNUARef->GetBinContent(l, cBin, GetBinNumberCos(2*n));
      sin2nPhi2 += fCumuNUARef->GetBinContent(l, cBin, GetBinNumberSin(2*n));
      mult2 += fCumuNUARef->GetBinContent(l, cBin, 0);
    }
    if (mult1 == 0 || mult2 == 0) { 
      if (fDebug > 0) AliWarning("mult1 or mult2 == 0 in matrix elements, aborting NUA");
      return 0.;
    }
    cosnMmPhi1 /= mult1;
    sinnMmPhi1 /= mult1;
    cosnPmPhi1 /= mult1;
    sinnPmPhi1 /= mult1;
    cos2nPhi1 /= mult1; 
    sin2nPhi1 /= mult1;
    cosnMmPhi2 /= mult2;
    sinnMmPhi2 /= mult2;
    cosnPmPhi2 /= mult2;
    sinnPmPhi2 /= mult2;
    cos2nPhi2 /= mult2; 
    sin2nPhi2 /= mult2;
  }

  // Actual calculation
  Double_t e = cosnMmPhi1*cosnMmPhi2 + sinnMmPhi1*sinnMmPhi2 + cosnPmPhi1*cosnPmPhi2 + sinnPmPhi1*sinnPmPhi2;
  Double_t den = 1 + cos2nPhi1*cos2nPhi2 + sin2nPhi1*sin2nPhi2;
  if (den != 0) e /= den;
  else return 0.;

  return e;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::AddVertexBins(CumuHistos& cumu, TList* list, UInt_t nNUA) const  
{
  //
  //  Add up vertex bins with flow results
  //
  //  Parameters:
  //   cumu: CumuHistos object with vtxbin results
  //   list: Outout list with added results
  //   nNUA: # of NUA histograms to loop over
  //
  TH2D* vtxHist = 0;
  TProfile2D* avgProf = 0;
  TString name;
  Int_t nT = ((fFlags & k3Cor) ? 3 : 2);
  Char_t ct = '\0';
  for (UInt_t nua = 0; nua <= nNUA; nua++) { // NUA loop
    for (Int_t n = 2; n <= fMaxMoment; n++) { // moment loop
      for (Int_t t = 0; t < nT; t++) { // type loop (r/d/a/b)
        // Find type
        switch (t) {
          case 0: ct = 'r'; break;
          case 1: ct = ((fFlags & k3Cor) ? 'a' : 'd'); break;
          case 2: ct = 'b'; break;
          default: ct = '\0'; break;
	}
	vtxHist = static_cast<TH2D*>(cumu.Get(ct,n,nua));
	if (!vtxHist) {
	  AliWarning("VertexBin::AddVertexBins: vtxHist not found!");
	  continue;
	}
	name = vtxHist->GetName();
	// Strip name of vtx info
	Ssiz_t l = name.Last('x')-3;
	name.Resize(l);
	avgProf = (TProfile2D*)list->FindObject(name.Data());
	// if no output profile yet, make one
	if (!avgProf) {
	  avgProf = new TProfile2D(name.Data(), name.Data(), 
	      vtxHist->GetNbinsX(), vtxHist->GetXaxis()->GetXmin(), vtxHist->GetXaxis()->GetXmax(),
	      vtxHist->GetNbinsY(), vtxHist->GetYaxis()->GetXmin(), vtxHist->GetYaxis()->GetXmax());
	  if (vtxHist->GetXaxis()->IsVariableBinSize()) 
	    avgProf->GetXaxis()->Set(vtxHist->GetNbinsX(), vtxHist->GetXaxis()->GetXbins()->GetArray());
	  if (vtxHist->GetYaxis()->IsVariableBinSize()) 
	    avgProf->GetYaxis()->Set(vtxHist->GetNbinsY(), vtxHist->GetYaxis()->GetXbins()->GetArray());
	  list->Add(avgProf);
	}
	// Fill in, cannot be done with Add function.
	for (Int_t e = 1; e <= vtxHist->GetNbinsX(); e++) { // eta loop
	  Double_t eta = vtxHist->GetXaxis()->GetBinCenter(e);
	  for (Int_t c = 1; c <= vtxHist->GetNbinsY(); c++) { // cent loop
	    Double_t cent = vtxHist->GetYaxis()->GetBinCenter(c);
	    Double_t cont = vtxHist->GetBinContent(e, c);
	    if (cont == 0) continue;
	    avgProf->Fill(eta, cent, cont);
	  } // End of cent loop
	} //  End of eta loop
      } // End of type loop
    } // End of moment loop
  } // End of nua loop
}
//_____________________________________________________________________
Int_t AliForwardFlowTaskQC::VertexBin::GetBinNumberCos(Int_t n) const  
{
  //
  //  Get the bin number of <<cos(nphi)>>
  //
  //  Parameters:
  //   n: moment
  //
  //  Return: bin number
  //
  Int_t bin = 0;
  n = TMath::Abs(n);
  
  if (n == 0) bin = fMaxMoment*4-1;
  else        bin = n*2-1;
  
  return bin;
}
//_____________________________________________________________________
Int_t AliForwardFlowTaskQC::VertexBin::GetBinNumberSin(Int_t n) const  
{
  //
  //  Get the bin number of <<sin(nphi)>>
  //
  //  Parameters:
  //   n: moment
  //
  //  Return: bin number
  //
  Int_t bin = 0;
  n = TMath::Abs(n);
  
  if (n == 0) bin = fMaxMoment*4;
  else        bin = n*2;

  return bin;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::SetupNUALabels(TAxis* a) const
{
  // 
  //  Setup NUA labels on axis
  //
  //  Parameters:
  //   a: Axis to set up
  //
  if (a->GetNbins() != GetBinNumberSin()) AliFatal("SetupNUALabels: Wrong number of bins on axis");

  Int_t i = 1, j= 1;
  while (i <= a->GetNbins()) {
    a->SetBinLabel(i++, Form("<<cos(%d#varphi)>>", j));
    a->SetBinLabel(i++, Form("<<sin(%d#varphi)>>", j++));
  }

  return;
}
//_____________________________________________________________________
TH2I* AliForwardFlowTaskQC::VertexBin::MakeQualityHist(const Char_t* name) const 
{
  //
  //  Add a histogram for checking the analysis quality
  //
  //  Parameters:
  //   const char*: name of data type
  //
  //  Return: histogram for tracking successful calculations
  //
  Int_t nBins = ((fFlags & kStdQC) ? 8 : 4);
  TH2I* quality = new TH2I(name, name, (fMaxMoment-1)*nBins, 1, (fMaxMoment-1)*nBins+1, fCumuNUARef->GetNbinsY(), 
                           fCumuNUARef->GetYaxis()->GetXmin(), fCumuNUARef->GetYaxis()->GetXmax());
  quality->GetYaxis()->Set(fCumuNUARef->GetNbinsY(), fCumuNUARef->GetYaxis()->GetXbins()->GetArray());
  for (Int_t i = 2, j = 1; i <= fMaxMoment; i++) {
    quality->GetXaxis()->SetBinLabel(j++, Form("QC_{%d}{2} > 0", i));
    quality->GetXaxis()->SetBinLabel(j++, Form("QC_{%d}{2} <= 0", i));
    quality->GetXaxis()->SetBinLabel(j++, Form("QC'_{%d}{2} > 0", i));
    quality->GetXaxis()->SetBinLabel(j++, Form("QC'_{%d}{2} <= 0", i));
    if ((fFlags & kStdQC)) {
      quality->GetXaxis()->SetBinLabel(j++, Form("QC_{%d}{4} < 0", i));
      quality->GetXaxis()->SetBinLabel(j++, Form("QC_{%d}{4} >= 0", i));
      quality->GetXaxis()->SetBinLabel(j++, Form("QC'_{%d}{4} < 0", i));
      quality->GetXaxis()->SetBinLabel(j++, Form("QC'_{%d}{4} >= 0", i));
    }
  }

  return quality;
}
//_____________________________________________________________________
TH2D* AliForwardFlowTaskQC::VertexBin::MakeOutputHist(Int_t qc, Int_t n, const Char_t* ctype, UInt_t nua) const
{
  //
  //  Setup a TH2D for the output
  //
  //  Parameters:
  //   qc   # of particle correlations
  //   n    flow moment
  //   ref  Type: r/d/a/b
  //   nua  For nua corrected hists?
  //
  //  Return: 2D hist for results
  //
  Bool_t ref = ((ctype[0] == 'R') || (ctype[0] == 'r'));
  TAxis* xAxis = (ref ? fCumuNUARef->GetXaxis() : fCumuNUADiff->GetXaxis());
  TAxis* yAxis = (ref ? fCumuNUARef->GetYaxis() : fCumuNUADiff->GetYaxis());
  TString ntype = "";
  switch (nua) {
    case CumuHistos::kNoNUA: ntype = "Un"; break;
    case CumuHistos::kNUAOld: ntype = "NUAOld"; break;
    case CumuHistos::kNUA: ntype = "NUA"; break;
    default: break;
  }
  TH2D* h = new TH2D(Form("%sQC%d_v%d_%s_%sCorr%s_vtx_%d_%d", 
			fType.Data(), qc, n, ctype, ntype.Data(),
			GetQCType(fFlags), fVzMin, fVzMax),
		      Form("%sQC%d_v%d_%s_%sCorr%s_vtx_%d_%d", 
			fType.Data(), qc, n, ctype, ntype.Data(),
			GetQCType(fFlags), fVzMin, fVzMax),
		      xAxis->GetNbins(), xAxis->GetXmin(), xAxis->GetXmax(),
		      yAxis->GetNbins(), yAxis->GetXmin(), yAxis->GetXmax());
  if (xAxis->IsVariableBinSize()) h->GetXaxis()->Set(xAxis->GetNbins(), xAxis->GetXbins()->GetArray());
  h->GetYaxis()->Set(yAxis->GetNbins(), yAxis->GetXbins()->GetArray());

  return h;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::PrintFlowSetup() const 
{
  //
  //  Print the setup of the flow task
  //
  Printf("=======================================================");
  Printf("%s::Print", this->IsA()->GetName());
  Printf("Forward detector:                 :\t%s", ((fFlowFlags & kFMD) ? "FMD" : "VZERO"));
  Printf("Number of bins in vertex axis     :\t%d", fVtxAxis->GetNbins());
  Printf("Range of vertex axis              :\t[%3.1f,%3.1f]", 
			  fVtxAxis->GetXmin(), fVtxAxis->GetXmax());
  Printf("Number of bins in centrality axis :\t%d", fCentAxis->GetNbins());
  printf("Centrality binning                :\t");
  for (Int_t cBin = 1; cBin <= fCentAxis->GetNbins(); cBin++) {
    printf("%02d-%02d%% ", Int_t(fCentAxis->GetBinLowEdge(cBin)), Int_t(fCentAxis->GetBinUpEdge(cBin)));
    if (cBin == fCentAxis->GetNbins()) printf("\n");
    else if (cBin % 4 == 0) printf("\n\t\t\t\t\t");
  }
  printf("Doing flow analysis for           :\t");
  for (Int_t n  = 2; n <= fMaxMoment; n++) printf("v%d ", n);
  printf("\n");
  TString type = "Standard QC{2} and QC{4} calculations.";
  if ((fFlowFlags & kEtaGap)) type = "QC{2} with a rapidity gap.";
  if ((fFlowFlags & k3Cor))   type = "QC{2} with 3 correlators.";
  if ((fFlowFlags & kTPC) == kTPC)    type.ReplaceAll(".", " with TPC tracks for reference.");
  if ((fFlowFlags & kHybrid) == kHybrid) type.ReplaceAll(".", " with hybrid tracks for reference.");
  Printf("QC calculation type               :\t%s", type.Data());
  Printf("Symmetrize ref. flow wrt. eta = 0 :\t%s", ((fFlowFlags & kSymEta) ? "true" : "false"));
  Printf("Apply NUA correction terms        :\t%s", ((fFlowFlags & kNUAcorr) ? "true" : "false"));
  Printf("Satellite vertex flag             :\t%s", ((fFlowFlags & kSatVtx) ? "true" : "false"));
  Printf("FMD sigma cut:                    :\t%f", fFMDCut);
  Printf("SPD sigma cut:                    :\t%f", fSPDCut);
  if ((fFlowFlags & kEtaGap) || (fFlowFlags & kTracks)) 
    Printf("Eta gap:                          :\t%f", fEtaGap);
  Printf("=======================================================");
}
//_____________________________________________________________________
const Char_t* AliForwardFlowTaskQC::GetQCType(UShort_t flags, Bool_t prependUS) 
{
  // 
  //  Get the type of the QC calculations
  //  Used for naming of objects in the VertexBin class, 
  //  important to avoid memory problems when running multiple
  //  initializations of the class with different setups
  // 
  //  Parameters:
  //   flags: Flow flags for type determination
  //   prependUS: Prepend an underscore (_) to the name
  // 
  //  Return: QC calculation type
  //
  static TString type = "";
  if      ((flags & kStdQC))  type = "StdQC";
  else if ((flags & kEtaGap)) type = "EtaGap";
  else if ((flags & k3Cor))   type = "3Cor";
  else type = "UNKNOWN";
  if (prependUS) type.Prepend("_");
  if ((flags & kTPC) == kTPC)    type.Append("TPCTr");
  if ((flags & kHybrid) == kHybrid) type.Append("HybTr");
  
  return type.Data();
}
//_____________________________________________________________________
TH1* AliForwardFlowTaskQC::CumuHistos::Get(Char_t t, Int_t n, UInt_t nua) const
{
  //
  //  Get histogram/profile for type t and moment n
  //
  //  Parameters:
  //   t: type = 'r'/'d'
  //   n: moment
  //   nua: NUA type
  //
  n = GetPos(n, nua);
  if (n < 0) AliFatal(Form("CumuHistos out of range: (%c,%d)", t, n));

  TH1* h = 0;
  Int_t pos = -1;
  if      (t == 'r' || t == 'R') pos = n;    
  else if (t == 'd' || t == 'D') pos = n;    
  else if (t == 'a' || t == 'A') pos = 2*n;  
  else if (t == 'b' || t == 'B') pos = 2*n+1;
  else AliFatal(Form("CumuHistos wrong type: %c", t));

  if (t == 'r' || t == 'R') {
    if (pos < fRefHists->GetEntries()) {
      h = (TH1*)fRefHists->At(pos);
    }
  } else if (pos < fDiffHists->GetEntries()) {
    h = (TH1*)fDiffHists->At(pos);
  }
  if (!h) AliFatal(Form("No hist found in list %c at pos %d", t, pos));

  return h;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::CumuHistos::ConnectList(TString name, TList* l)
{
  //
  //  Connect an output list with this object
  //
  //  Parameters:
  //   name: base name
  //   l: list to keep outputs in
  //
  TString ref = name;
  ref.ReplaceAll("Cumu","CumuRef");
  fRefHists = (TList*)l->FindObject(ref.Data());
  if (!fRefHists) {
    fRefHists = new TList();
    fRefHists->SetName(ref.Data());
    l->Add(fRefHists);
  } else {
    // Check that the list correspond to fMaxMoments
    if (fRefHists->GetEntries() != (fMaxMoment-1.)*(fNUA+1)) 
      AliError("CumuHistos::ConnectList Wrong number of hists in ref list,"
               "you are doing something wrong!");
  }
  TString diff = name;
  diff.ReplaceAll("Cumu","CumuDiff");
  fDiffHists = (TList*)l->FindObject(diff.Data());
  if (!fDiffHists) {
    fDiffHists = new TList();
    fDiffHists->SetName(diff.Data());
    l->Add(fDiffHists);
  } else {
    // Check that the list correspond to fMaxMoment
    if ((fDiffHists->GetEntries() != (fMaxMoment-1.)*(fNUA+1)) &&  
        (fDiffHists->GetEntries() != 2*(fMaxMoment-1.)*(fNUA+1)))  
      AliError(Form("CumuHistos::ConnectList Wrong number of hists in diff list,"
               "you are doing something wrong! (%s)",name.Data()));
  }

}
//_____________________________________________________________________
void AliForwardFlowTaskQC::CumuHistos::Add(TH1* h) const
{
  //
  //  Add a histogram to this objects lists
  //
  //  Parameters:
  //   h: histogram/profile to add
  //
  TString name = h->GetName();
  if (name.Contains("Ref")) {
    if (fRefHists) fRefHists->Add(h);
    else AliFatal("CumuHistos::Add() - fRefHists does not exist");
    // Check that the list correspond to fMaxMoments
    if (fRefHists->GetEntries() > (fNUA+1)*(fMaxMoment-1.)) 
      AliError("CumuHistos::Add wrong number of hists in ref list, "
               "you are doing something wrong!");
  }
  else if (name.Contains("Diff")) {
    if (fDiffHists) fDiffHists->Add(h);
    else AliFatal("CumuHistos::Add() - fDiffHists does not exist");
    // Check that the list correspond to fMaxMoment
    if (fDiffHists->GetEntries() > 2*(fNUA+1)*(fMaxMoment-1.)) 
      AliError("CumuHistos::Add wrong number of hists in diff list, "
	 "you are doing something wrong!");
  }
  return;
}
//_____________________________________________________________________
Int_t AliForwardFlowTaskQC::CumuHistos::GetPos(Int_t n, UInt_t nua) const
{
  //
  //  Get position in list of moment n objects
  //  To take care of different numbering schemes
  //
  //  Parameters:
  //   n: moment
  //   nua: # of nua corrections applied
  //
  //  Return: position #
  //
  if (n > fMaxMoment) return -1;
  else return (n-2)+nua*(fMaxMoment-1);
}
//_____________________________________________________________________
//
//
// EOF
