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
#include "TH3D.h"
#include "TProfile2D.h"
#include "AliLog.h"
#include "AliForwardFlowTaskQC.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"

ClassImp(AliForwardFlowTaskQC)
#if 0
; // For emacs 
#endif

AliForwardFlowTaskQC::AliForwardFlowTaskQC()
  : AliAnalysisTaskSE(),
    fBinsFMD(),         // List with FMD flow histos
    fBinsSPD(),         // List with SPD flow histos
    fSumList(0),	// Event sum list
    fOutputList(0),	// Result output list
    fAOD(0),		// AOD input event
    fZvertex(1111),	// Z vertex coordinate
    fCent(-1),		// Centrality
    fHistCent(),        // Histo for centrality
    fHistVertexSel(),   // Histo for selected vertices
    fHistVertexAll()    // Histo for all vertices
{
  // 
  // Default constructor
  //
  for (Int_t n = 0; n <= 6; n++) fv[n] = kTRUE;
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const char* name) 
  : AliAnalysisTaskSE(name),
    fBinsFMD(),         // List with FMD flow histos
    fBinsSPD(),         // List with SPD flow histos
    fSumList(0),        // Event sum list           
    fOutputList(0),     // Result output list       
    fAOD(0),	        // AOD input event          
    fZvertex(1111),     // Z vertex coordinate      
    fCent(-1),          // Centrality               
    fHistCent(),        // Histo for centrality
    fHistVertexSel(),   // Histo for selected vertices
    fHistVertexAll()    // Histo for all vertices      
{
  // 
  // Constructor
  //
  // Parameters:
  //  name: Name of task
  //
  for (Int_t n = 0; n <= 6; n++) fv[n] = kTRUE;

  fHistCent      = new TH1D("hCent", "Centralities", 100, 0, 100);
  fHistVertexSel = new TH1D("hVertexSel", "Selected vertices", 40, -20, 20);
  fHistVertexAll = new TH1D("hVertexAll", "All vertices", 40, -20, 20);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const AliForwardFlowTaskQC& o)
  : AliAnalysisTaskSE(o),
    fBinsFMD(),                        // List with FMD flow histos
    fBinsSPD(),                        // List with SPD flow histos
    fSumList(o.fSumList),              // Event sum list           
    fOutputList(o.fOutputList),        // Result output list       
    fAOD(o.fAOD),	               // AOD input event          
    fZvertex(o.fZvertex),              // Z vertex coordinate      
    fCent(o.fCent),	               // Centrality               
    fHistCent(o.fHistCent),            // Histo for centrality
    fHistVertexSel(o.fHistVertexSel),  // Histo for selected vertices
    fHistVertexAll(o.fHistVertexAll)   // Histo for all vertices      
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  for (Int_t n = 0; n <= 6; n++) fv[n] = o.fv[n];
}
//_____________________________________________________________________
AliForwardFlowTaskQC&
AliForwardFlowTaskQC::operator=(const AliForwardFlowTaskQC& o)
{
  // 
  // Assignment operator 
  //
  fSumList       = o.fSumList;
  fOutputList    = o.fOutputList;
  fAOD           = o.fAOD;
  fZvertex       = o.fZvertex;
  fCent          = o.fCent;
  fHistCent      = o.fHistCent;
  fHistVertexSel = o.fHistVertexSel;
  fHistVertexAll = o.fHistVertexAll;
  fHistVertexAll = o.fHistVertexAll;

  for (Int_t n = 0; n <= 6; n++) fv[n] = o.fv[n];
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

  PostData(1, fSumList);
  PostData(2, fOutputList);

}
//_____________________________________________________________________
void AliForwardFlowTaskQC::InitVertexBins()
{
  // 
  // Init vertexbin objects for FMD and SPD, and add them to the lists
  //
  for(Int_t n = 1; n <= 6; n++) {
    if (!fv[n]) continue;
    for (Int_t v = -10; v < 10; v++) {
      fBinsFMD.Add(new VertexBin(v, v+1, n, "FMD"));
      fBinsSPD.Add(new VertexBin(v, v+1, n, "SPD"));
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

  TList* dList = new TList();
  dList->SetName("Diagnostics");
  dList->Add(fHistCent);
  dList->Add(fHistVertexSel);
  dList->Add(fHistVertexAll);
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
  fZvertex = 1111;

  // Get input event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) return kFALSE;

  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
  if (!aodfmult) return kFALSE;
  if (!AODCheck(aodfmult)) return kFALSE;
  TH2D fmddNdetadphi = aodfmult->GetHistogram();

  AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClusters"));
  if (!aodcmult) return kFALSE;
  TH2D spddNdetadphi = aodcmult->GetHistogram();

  // TODO: remove me!
//  fCent = 0.5;

  // Run analysis on FMD and SPD
  TIter nextFMD(&fBinsFMD);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(nextFMD()))) {
    if (bin->CheckVertex(fZvertex)) {
      if (!bin->FillHists(fmddNdetadphi)) return kFALSE; 
      bin->CumulantsAccumulate(fCent);
    }
  }

  TIter nextSPD(&fBinsSPD);
  while ((bin = static_cast<VertexBin*>(nextSPD()))) {
    if (bin->CheckVertex(fZvertex)) {
      if (!bin->FillHists(spddNdetadphi)) return kFALSE;
      bin->CumulantsAccumulate(fCent);
    }
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

  // Run finalize on VertexBins
  Finalize();

  // Collect centralities
  TProfile2D* hist2D = 0;
  TList* centList = 0;
  TH1D* hist1D = 0;
  TIter nextProfile(fOutputList);
  while ((hist2D = dynamic_cast<TProfile2D*>(nextProfile()))) {
    for (Int_t cBin = 1; cBin <= 100; ) {
      Int_t cMin = cBin - 1;
      Int_t cMax = (cMin < 80 ? (cMin < 20 ? cMin + 5 : cMin + 10) : cMin + 20);
      TString name = Form("cent_%d-%d", cMin, cMax);
      centList = (TList*)fOutputList->FindObject(name.Data());
      if (!centList) { 
	centList = new TList();
	centList->SetName(name.Data());
	fOutputList->Add(centList);
      }
      hist1D = hist2D->ProjectionX(Form("%s_%s", hist2D->GetName(), name.Data()), 
				   cMin, cMax, "E");
      hist1D->SetTitle(hist1D->GetName());
      hist1D->Scale(1./(cMax-cMin));
      centList->Add(hist1D);

      cBin = cMax+1;
    }
  }

  PostData(2, fOutputList);

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
  TIter nextFMD(&fBinsFMD);
  VertexBin* bin = 0;
  while ((bin = static_cast<VertexBin*>(nextFMD()))) {
    bin->CumulantsTerminate(fSumList, fOutputList);
  }
  TIter nextSPD(&fBinsSPD);
  while ((bin = static_cast<VertexBin*>(nextSPD()))) {
    bin->CumulantsTerminate(fSumList, fOutputList);
  }

} 
// _____________________________________________________________________
Bool_t AliForwardFlowTaskQC::AODCheck(const AliAODForwardMult* aodfm) 
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
  if (!aodfm->IsTriggerBits(AliAODForwardMult::kOffline)) return kFALSE;

  // Then check for centrality
  fCent = (Double_t)aodfm->GetCentrality();
  if (0. >= fCent || fCent >= 100.) return kFALSE;
  fHistCent->Fill(fCent);

  // And finally check for vertex
  fZvertex = aodfm->GetIpZ();
  fHistVertexAll->Fill(fZvertex);
  if (TMath::Abs(fZvertex) >= 10.) return kFALSE;
  fHistVertexSel->Fill(fZvertex);

  return kTRUE;

}
//_____________________________________________________________________
AliForwardFlowTaskQC::VertexBin::VertexBin()
  : TNamed(),
    fMoment(0),      // Flow moment for this vertexbin
    fVzMin(0),       // Vertex z-coordinate min
    fVzMax(0),       // Vertex z-coordinate max
    fType(),         // Data type (FMD/SPD/FMDTR/SPDTR/MC)
    fCumuRef(),      // Histogram for reference flow
    fCumuDiff(),     // Histogram for differential flow
    fCumuHist(),     // Sum histogram for cumulants
    fHistTwoCorr(),  // Diagnostics histogram for <<2>>
    fHistW2(),       // Diagnostics histogram for w_<<2>>
    fHistFourCorr(), // Diagnostics histogram for <<4>>
    fHistW4(),       // Diagnostics histogram for w_<<4>>      
    fdNdedpAcc()     // Diagnostics histogram to make acc. maps
{
  //
  // Default constructor
  //
}
//_____________________________________________________________________
AliForwardFlowTaskQC::VertexBin::VertexBin(const Int_t vLow, const Int_t vHigh, 
                                           const Int_t moment, const Char_t* name)
  : TNamed("", ""),
    fMoment(moment),  // Flow moment for this vertexbin
    fVzMin(vLow),     // Vertex z-coordinate min
    fVzMax(vHigh),    // Vertex z-coordinate max
    fType(name),      // Data type (FMD/SPD/FMDTR/SPDTR/MC)
    fCumuRef(),       // Histogram for reference flow
    fCumuDiff(),      // Histogram for differential flow
    fCumuHist(),      // Sum histogram for cumulants
    fHistTwoCorr(),   // Diagnostics histogram for <<2>>
    fHistW2(),        // Diagnostics histogram for w_<<2>>
    fHistFourCorr(),  // Diagnostics histogram for <<4>>
    fHistW4(),        // Diagnostics histogram for w_<<4>> 
    fdNdedpAcc()      // Diagnostics histogram to make acc. maps
{
  //
  // Constructor
  //
  // Parameters
  //  vLow: min z-coordinate
  //  vHigh: max z-coordinate
  //  moment: flow moment
  //  name: data type name (FMD/SPD/FMDTR/SPDTR/MC)
  //
  SetName(Form("%svertexBin%d_%d_%d", name, moment, vLow, vHigh));
  SetTitle(Form("%svertexBin%d_%d_%d", name, moment, vLow, vHigh));

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
  fType         = o.fType;
  fCumuRef      = o.fCumuRef;
  fCumuDiff     = o.fCumuDiff;
  fCumuHist     = o.fCumuHist;
  fHistTwoCorr  = o.fHistTwoCorr;
  fHistW2       = o.fHistW2;
  fHistFourCorr = o.fHistFourCorr;
  fHistW4       = o.fHistW4;
  fdNdedpAcc    = o.fdNdedpAcc;

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
  TList* list = (TList*)outputlist->FindObject(Form("%svertex_%d_%d", fType, fVzMin, fVzMax));

  // If it doesn't exist we make one
  if (!list) {
    list = new TList();
    list->SetName(Form("%svertex_%d_%d", fType, fVzMin, fVzMax));
    outputlist->Add(list);
  }

  // We initiate the reference histogram according to an acceptance correction map,
  // so we don't shift the SPD coverage within a reference bin
  fCumuRef = new TH2D(Form("%s_v%d_%d_%d_ref", fType, fMoment, fVzMin, fVzMax),
                        Form("%s_v%d_%d_%d_ref", fType, fMoment, fVzMin, fVzMax),
                        24, -6., 6., 5, 0.5, 5.5);
  TFile acc("$ALICE_ROOT/PWGLF/FORWARD/corrections/FlowCorrections/FlowAccMap.root", "READ");
  TH1D* accMap = (TH1D*)acc.Get(Form("%saccVertex_%d_%d", fType, fVzMin, fVzMax));
  if (accMap) {
    Int_t nBins = accMap->GetNbinsX();
    Double_t eta[48] = { 0. };
    Int_t n = 0;
    Double_t newOcc[48] = { 0. };
    Double_t prev = -1;
    for (Int_t i = 0; i < nBins; i++) {
      Double_t occ = accMap->GetBinContent(i+1);
      if (prev != occ && (occ > 0.6 || occ == 0)) {
        eta[n] = i*0.25-6.;
        newOcc[n] = occ;
        n++;
//        printf("eta: %f \t occ: %f \t Vertex: %d \n", eta[n-1], occ, fVzMin);
      }
      prev = occ;
    }
    eta[n] = 6.;
    fCumuRef->GetXaxis()->Set(n, eta);
  }
  acc.Close();

  fCumuRef->Sumw2();
  list->Add(fCumuRef);

  // We initiate the differential histogram
  fCumuDiff = new TH2D(Form("%s_v%d_%d_%d_diff", fType, fMoment, fVzMin, fVzMax),
                        Form("%s_v%d_%d_%d_diff", fType, fMoment, fVzMin, fVzMax),
                        48, -6., 6., 5, 0.5, 5.5);
  fCumuDiff->Sumw2();
  list->Add(fCumuDiff);

  // Initiate the cumulant sum histogram
  fCumuHist = new TH3D(Form("%sv%d_vertex_%d_%d", fType, fMoment, fVzMin, fVzMax),
                       Form("%sv%d_vertex_%d_%d", fType, fMoment, fVzMin, fVzMax),
                       48, -6., 6., 100, 0., 100., 26, 0.5, 26.5);
  fCumuHist->Sumw2();

  list->Add(fCumuHist);

  // We check for diagnostics histograms (only done per type and moment, not vertexbin)
  // If they are not found we create them.
  TList* dList = (TList*)outputlist->FindObject("Diagnostics");
  if (!dList) AliFatal("No diagnostics list found, what kind of game are you running here?!?!");

  // Corr. hists are shared over all vertex bins...
  fHistTwoCorr    = (TH3D*)dList->FindObject(Form("hHistTwoCorr_%s_v%d", fType, fMoment));
  if (fHistTwoCorr) {
    fHistW2       = (TH3D*)dList->FindObject(Form("hHistW2_%s_v%d", fType, fMoment));
    fHistFourCorr = (TH3D*)dList->FindObject(Form("hHistFourCorr_%s_v%d", fType, fMoment));
    fHistW4       = (TH3D*)dList->FindObject(Form("hHistW4_%s_v%d", fType, fMoment));
  } else {
    fHistTwoCorr  = new TH3D(Form("hHistTwoCorr_%s_v%d", fType, fMoment), 
                            "Two particle correlator: w_{2}<<2>>", 
                            48, -6., 6., 100, 0., 150000, 100, 0., 100.);
    fHistTwoCorr->Sumw2();
    fHistW2       = new TH3D(Form("hHistW2_%s_v%d", fType, fMoment),
                            "Two particle event weight: w_{2}",
                            48, -6., 6., 100, 0., 2e+7, 100, 0., 100.);
    fHistW2->Sumw2();
    fHistFourCorr = new TH3D(Form("hHistFourCorr_%s_v%d", fType, fMoment),  
                            "Four particle correlator: w_{4}<<4>>", 
                            48, -6., 6., 100, 0., 1e+10, 100, 0., 100.);
    fHistFourCorr->Sumw2();
    fHistW4       = new TH3D(Form("hHistW4_%s_v%d", fType, fMoment), 
                            "Four particle event weight: w_{4}",
                            48, -6., 6., 100, 0., 4e+14, 100, 0., 100.);
    fHistW4->Sumw2();

    dList->Add(fHistTwoCorr);
    dList->Add(fHistW2);
    dList->Add(fHistFourCorr);
    dList->Add(fHistW4);
  }
  
  // Acceptance hists are shared over all moments
  fdNdedpAcc = (TH2D*)dList->FindObject(Form("h%sdNdedpAcc_%d_%d", fType, fVzMin, fVzMax));
  if (!fdNdedpAcc) {
    fdNdedpAcc = new TH2D(Form("h%sdNdedpAcc_%d_%d", fType, fVzMin, fVzMax), 
                          Form("%s acceptance map for %d cm < v_{z} < %d cm", fType, fVzMin, fVzMax),
                          48, -6, 6, 20, 0, TMath::TwoPi());
    fdNdedpAcc->Sumw2();
    dList->Add(fdNdedpAcc);
  }


}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::VertexBin::CheckVertex(Double_t vz)
{
  //
  // We check if this is the correct bin for the current event's vertex
  //
  // Parameters:
  //  vZ: Current event vertex
  //
  // Returns false if out of range, true otherwise
  //
  if ((Double_t)fVzMin > vz) return kFALSE;
  if ((Double_t)fVzMax <= vz) return kFALSE;

  return kTRUE; 
}
//_____________________________________________________________________
Bool_t AliForwardFlowTaskQC::VertexBin::FillHists(TH2D dNdetadphi) 
{
  // 
  // Fill reference and differential eta-histograms
  //
  // Parameters:
  //  dNdetadphi: 2D histogram with input data
  //

  if (!fCumuRef) AliFatal("You have not called AddOutput() - Terminating!");

  // Fist we reset histograms
  fCumuRef->Reset();
  fCumuDiff->Reset();

  // Numbers to cut away bad events and acceptance.
  Double_t runAvg = 0;
  Double_t max = 0;
  Int_t nInAvg = 0;
  Int_t nBadBins = 0;
  Int_t nBins = (dNdetadphi.GetNbinsX() * 6) / (fCumuDiff->GetNbinsX() * 5);
  Int_t nInBin = 0;
  Int_t nCurBin = 0, nPrevBin = 0;

  // Then we loop over the input and calculate sum cos(k*n*phi)
  // and fill it in the reference and differential histograms
  Double_t eta, phi, weight;
  Double_t dQnRe = 0, dQ2nRe = 0, dQnIm = 0, dQ2nIm = 0;
  for (Int_t etaBin = 1; etaBin <= dNdetadphi.GetNbinsX(); etaBin++) {
    eta = dNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    nCurBin = fCumuDiff->GetXaxis()->FindBin(eta);
    // If we have moved to a new bin in the flow hist, and less than half the eta
    // region has been covered by it we cut it away.
    if (nPrevBin && nCurBin != nPrevBin) {
      if (nInBin <= nBins/2) {
	for (Int_t pBin = 1; pBin <= fCumuDiff->GetNbinsY(); pBin++) {
	  fCumuDiff->SetBinContent(nPrevBin, pBin, 0);
	  fCumuDiff->SetBinError(nPrevBin, pBin, 0);
	}
      }
      nInBin = 0;
      nPrevBin = nCurBin;
      runAvg = 0;
      nInAvg = 0;
      max = 0;
    }
    Bool_t data = kFALSE;
    for (Int_t phiBin = 1; phiBin <= dNdetadphi.GetNbinsY(); phiBin++) {
      phi = dNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      weight = dNdetadphi.GetBinContent(etaBin, phiBin);
      if (!weight) continue;
      if (!data) data = kTRUE;
      // We calculate the running average Nch per. bin
      runAvg *= nInAvg;
      runAvg += weight;
      nInAvg++;
      runAvg /= nInAvg;
      // if the current bin has more than write the avg we count a bad bin
      if (weight > max) max = weight;

      dQnRe = weight*TMath::Cos(fMoment*phi);
      dQnIm = weight*TMath::Sin(fMoment*phi);
      dQ2nRe = weight*TMath::Cos(2.*fMoment*phi);
      dQ2nIm = weight*TMath::Sin(2.*fMoment*phi);

      fCumuRef->Fill(eta, kHmult, weight);
      fCumuRef->Fill(eta, kHQnRe, dQnRe);
      fCumuRef->Fill(eta, kHQnIm, dQnIm);
      fCumuRef->Fill(eta, kHQ2nRe, dQ2nRe);
      fCumuRef->Fill(eta, kHQ2nIm, dQ2nIm);

      fCumuDiff->Fill(eta, kHmult, weight);
      fCumuDiff->Fill(eta, kHQnRe, dQnRe);
      fCumuDiff->Fill(eta, kHQnIm, dQnIm);
      fCumuDiff->Fill(eta, kHQ2nRe, dQ2nRe);
      fCumuDiff->Fill(eta, kHQ2nIm, dQ2nIm);

      // Fill acc. map
      fdNdedpAcc->Fill(eta, phi, weight);
    }
    if (data) {
      nInBin++;
      if (max > 35) nBadBins++;
    }
    // If there are too many bad bins we throw the event away!
    if (nBadBins > 2) return kFALSE;
  }
  return kTRUE;
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::VertexBin::CumulantsAccumulate(Double_t cent) 
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
  Double_t multi = 0, multp = 0, mp = 0, mq = 0;
  Double_t w2 = 0, w4 = 0, w2p = 0, w4p = 0;
  Int_t refEtaBin = 0;
  Bool_t eventCount = kFALSE;
  
  // We loop over the data 1 time!
  for (Int_t etaBin = 1; etaBin <= fCumuDiff->GetNbinsX(); etaBin++) {
    eta = fCumuDiff->GetXaxis()->GetBinCenter(etaBin);
    refEtaBin = fCumuRef->GetXaxis()->FindBin(eta);
    // The values for each individual etaBin bins are reset
    mp = 0;
    pnRe = 0;
    p2nRe = 0;
    pnIm = 0;
    p2nIm = 0;

    mult = 0;
    dQnRe = 0;
    dQnIm = 0;
    dQ2nRe = 0;
    dQ2nIm = 0;

    // Reference flow
    multi = fCumuRef->GetBinContent(refEtaBin, kHmult);
    dQnRe = fCumuRef->GetBinContent(refEtaBin, kHQnRe);
    dQnIm = fCumuRef->GetBinContent(refEtaBin, kHQnIm);
    dQ2nRe = fCumuRef->GetBinContent(refEtaBin, kHQ2nRe);
    dQ2nIm = fCumuRef->GetBinContent(refEtaBin, kHQ2nIm);
    mult += multi;
    
    // For each etaBin bin the necessary values for differential flow
    // is calculated. Here is the loop over the phi's.
    multp = fCumuDiff->GetBinContent(etaBin, kHmult);
    pnRe = fCumuDiff->GetBinContent(etaBin, kHQnRe);
    pnIm = fCumuDiff->GetBinContent(etaBin, kHQnIm);
    p2nRe = fCumuDiff->GetBinContent(etaBin, kHQ2nRe);
    p2nIm = fCumuDiff->GetBinContent(etaBin, kHQ2nIm);
    mp += multp;
    
    if (mult <= 3) continue; 

    if (!eventCount) {
     // Count number of events
      fCumuHist->Fill(-7., cent, -0.5, 1.);
      eventCount = kTRUE;
    } 
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

    // Diagnostics are filled
    fHistTwoCorr->Fill(eta, two, cent);
    fHistW2->Fill(eta, w2, cent);
    fHistFourCorr->Fill(eta, four, cent);
    fHistW4->Fill(eta, w4, cent);

    // Differential flow calculations for each eta bin bin is done:
    mq = mp;
    qnRe = pnRe;
    qnIm = pnIm;
    q2nRe = p2nRe;
    q2nIm = p2nIm;

    // 2-particle differential flow
    w2p = mp * mult - mq;
    twoPrime = pnRe*dQnRe + pnIm*dQnIm - mq;
    
    fCumuHist->Fill(eta, cent, kw2two, twoPrime);
    fCumuHist->Fill(eta, cent, kw2, w2p);

    fCumuHist->Fill(eta, cent, kpnRe, pnRe);
    fCumuHist->Fill(eta, cent, kpnIm, pnIm);
    fCumuHist->Fill(eta, cent, kmp, mp);

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
    TList* list = (TList*)inlist->FindObject(Form("%svertex_%d_%d", fType, fVzMin, fVzMax));
    fCumuHist = (TH3D*)list->FindObject(Form("%sv%d_vertex_%d_%d", fType, fMoment, fVzMin, fVzMax));
  }

  // Create result profiles
  TProfile2D* cumu2 = (TProfile2D*)outlist->FindObject(Form("%sQC2_v%d_unCorr", fType, fMoment)); 
  TProfile2D* cumu4 = (TProfile2D*)outlist->FindObject(Form("%sQC4_v%d_unCorr", fType, fMoment)); 
  if (!cumu2) {
    cumu2 = new TProfile2D(Form("%sQC2_v%d_unCorr", fType, fMoment),
                           Form("%sQC2_v%d_unCorr", fType, fMoment),
                           48, -6., 6., 100, 0., 100);
    outlist->Add(cumu2);
  }
  if (!cumu4) {
    cumu4 = new TProfile2D(Form("%sQC4_v%d_unCorr", fType, fMoment),
                           Form("%sQC4_v%d_unCorr", fType, fMoment),
                           48, -6., 6., 100, 0., 100);
    outlist->Add(cumu4);
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
  for (Int_t c = 1; c <= 80; c++) {
    Double_t nEv = 0;
    // Eta loop
    for (Int_t etaBin = 1; etaBin <= fCumuHist->GetNbinsX(); etaBin++) {
      Double_t eta = fCumuHist->GetXaxis()->GetBinCenter(etaBin);
      // 2-particle reference flow
      w2Two = fCumuHist->GetBinContent(etaBin, c, kW2Two);
      w2 = fCumuHist->GetBinContent(etaBin, c, kW2);
      mult = fCumuHist->GetBinContent(etaBin, c, kM);
      if (!w2 || !mult) continue;
      cosP1nPhi = fCumuHist->GetBinContent(etaBin, c, kQnRe);
      sinP1nPhi = fCumuHist->GetBinContent(etaBin, c, kQnIm);
        
      cosP1nPhi /= mult;
      sinP1nPhi /= mult;
      two = w2Two / w2;
      qc2 = two - TMath::Power(cosP1nPhi, 2) - TMath::Power(sinP1nPhi, 2);
      if (qc2 <= 0) continue;
      vnTwo = TMath::Sqrt(qc2);
 //     if (!TMath::IsNaN(vnTwo*mult)) 
 //       cumu2->Fill(eta, vnTwo, fCumuHist->GetBinContent(0,c,0)); 

      // 4-particle reference flow
      w4Four = fCumuHist->GetBinContent(etaBin, c, kW4Four);
      w4 = fCumuHist->GetBinContent(etaBin, c, kW4);
      multm1m2 = fCumuHist->GetBinContent(etaBin, c, kMm1m2);
      if (!w4 || !multm1m2) continue;
      cosP1nPhi1P1nPhi2 = fCumuHist->GetBinContent(etaBin, c, kCosphi1phi2);
      sinP1nPhi1P1nPhi2 = fCumuHist->GetBinContent(etaBin, c, kSinphi1phi2);
      cosP1nPhi1M1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, c, kCosphi1phi2phi3m);
      sinP1nPhi1M1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, c, kSinphi1phi2phi3m);

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
      
      if (qc4 >= 0) continue;
      vnFour = TMath::Power(-qc4, 0.25);
 //     if (!TMath::IsNaN(vnFour*mult)) 
 //         cumu4->Fill(eta, vnFour, fCumuHist->GetBinContent(0,c,0));

      // 2-particle differential flow
      w2pTwoPrime = fCumuHist->GetBinContent(etaBin, c, kw2two);
      w2p = fCumuHist->GetBinContent(etaBin, c, kw2);
      mp = fCumuHist->GetBinContent(etaBin, c, kmp);
      if (!w2p || !mp) continue;
      cosP1nPsi = fCumuHist->GetBinContent(etaBin, c, kpnRe);
      sinP1nPsi = fCumuHist->GetBinContent(etaBin, c, kpnIm);

      cosP1nPsi /= mp;
      sinP1nPsi /= mp;
      twoPrime = w2pTwoPrime / w2p;
      qc2Prime = twoPrime - sinP1nPsi*sinP1nPhi - cosP1nPsi*cosP1nPhi;

      vnTwoDiff = qc2Prime / TMath::Sqrt(qc2);
      if (!TMath::IsNaN(vnTwoDiff*mp)) cumu2->Fill(eta, (Double_t)c-1., vnTwoDiff, fCumuHist->GetBinContent(0,c,0));

      // 4-particle differential flow
      w4pFourPrime = fCumuHist->GetBinContent(etaBin, c, kw4four);
      w4p = fCumuHist->GetBinContent(etaBin, c, kw4);
      mpqMult = fCumuHist->GetBinContent(etaBin, c, kmpmq);
      if (!w4p || !mpqMult) continue;
      cosP1nPsi1P1nPhi2 = fCumuHist->GetBinContent(etaBin, c, kCospsi1phi2);
      sinP1nPsi1P1nPhi2 = fCumuHist->GetBinContent(etaBin, c, kSinpsi1phi2);
      cosP1nPsi1M1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, c, kCospsi1phi2phi3m);
      sinP1nPsi1M1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, c, kSinpsi1phi2phi3m);
      cosP1nPsi1P1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, c, kCospsi1phi2phi3p);
      sinP1nPsi1P1nPhi2M1nPhi3 = fCumuHist->GetBinContent(etaBin, c, kSinpsi1phi2phi3p); 
      
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
      if (!TMath::IsNaN(vnFourDiff*mp)) cumu4->Fill(eta, (Double_t)c-1., vnFourDiff, fCumuHist->GetBinContent(0,c,0));
    } // End of eta loop
    // Number of events:
    nEv += fCumuHist->GetBinContent(0,c,0);
    cumu2->Fill(7., (Double_t)c-1., nEv);
    cumu4->Fill(7., (Double_t)c-1., nEv);
  } // End of centrality loop
  
  return;
}
//_____________________________________________________________________
//
//
// EOF
