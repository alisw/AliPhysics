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
#include "THnSparse.h"
#include "TH3D.h"
#include "TProfile2D.h"
#include "AliLog.h"
#include "AliForwardFlowTaskQC.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODForwardMult.h"
#include "AliAODEvent.h"

//
// Enumeration for adding and retrieving stuff from the histogram
//
enum { kW2Two = 1, kW2, kW4Four, kW4, kQnRe, kQnIm, kM,
       kCosphi1phi2, kSinphi1phi2, kCosphi1phi2phi3m, kSinphi1phi2phi3m, kMm1m2, 
       kw2two, kw2, kw4four, kw4, kpnRe, kpnIm, kmp, 
       kCospsi1phi2, kSinpsi1phi2, kCospsi1phi2phi3m, kSinpsi1phi2phi3m,
       kmpmq, kCospsi1phi2phi3p, kSinpsi1phi2phi3p };

ClassImp(AliForwardFlowTaskQC)
#if 0
; // For emacs 
#endif

AliForwardFlowTaskQC::AliForwardFlowTaskQC()
  : fOutputList(0),	// Output list
    fFlowUtil(0),	// AliForwardFlowUtil
    fAOD(0),		// AOD input event
    fMC(kFALSE),	// MC flag
    fEtaBins(48),	// # of etaBin bins in histograms
    fEtaRef(12),        // # of Eta bins for reference flow
    fAddFlow(0),	// Add flow string
    fAddType(0),	// Add flow type #
    fAddOrder(0),	// Add flow order
    fZvertex(1111),	// Z vertex range
    fCent(-1)		// Centrality
{
  // 
  // Default constructor
  //
  for (Int_t n = 0; n <= 6; n++) fv[n] = kTRUE;
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const char* name) :
  AliAnalysisTaskSE(name),
  fOutputList(0),	// Output list
  fFlowUtil(0),		// AliForwardFlowUtil
  fAOD(0),		// AOD input event
  fMC(kFALSE),		// MC flag
  fEtaBins(48),		// # of Eta bins
  fEtaRef(12),          // # of Eta bins for reference flow
  fAddFlow(0),		// Add flow string
  fAddType(0),		// Add flow type #
  fAddOrder(0),		// Add flow order
  fZvertex(1111),	// Z vertex range
  fCent(-1)		// Centrality
{
  // 
  // Constructor
  //
  // Parameters:
  //  name: Name of task
  //
  for (Int_t n = 0; n <= 6; n++) fv[n] = kTRUE;
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const AliForwardFlowTaskQC& o) :
  AliAnalysisTaskSE(o),
  fOutputList(o.fOutputList),	// Output list
  fFlowUtil(o.fFlowUtil),	// AliForwardFlowUtil
  fAOD(o.fAOD),			// AOD input event
  fMC(o.fMC),			// MC flag
  fEtaBins(o.fEtaBins),		// # of Eta bins
  fEtaRef(o.fEtaRef),           // # of Eta bins for reference flow
  fAddFlow(o.fAddFlow),		// Add flow string
  fAddType(o.fAddType),		// Add flow type #
  fAddOrder(o.fAddOrder),	// Add flow order
  fZvertex(o.fZvertex),		// Z vertex range
  fCent(o.fCent)		// Centrality
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  for (Int_t n = 0; n <= 6; n++) fv[n] = o.fv[n];
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::UserCreateOutputObjects()
{
  //
  // Create output objects
  // 
  if (!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("QCumulants");
  fOutputList->SetOwner();
  fFlowUtil = new AliForwardFlowUtil(fOutputList);
  
  Double_t x[101] = { 0. };
  // Double_t etaMin = -6;
  // Double_t etaMax = 6;

  // First we have a number of options for eta binning, also if it is
  // not really eta, but centrality or pt we want to do flow as a
  // function of, then this is possible:
  if (fEtaBins == 5) {
    x[0] = 0.;
    x[1] = 1.;
    x[2] = 2.;
    x[3] = 3.;
    x[4] = 4.5;
    x[5] = 6.0;
    // etaMin = 0;
    // etaMax = 6;
  }

  else if (fEtaBins == 100) { 
    for (Int_t e = 0; e<= 100; e++) {
      x[e] = e;
    }
    // etaMin = 0;
    // etaMax = 100;
  }
  
  else {
    if (fEtaBins % 6) fEtaBins = 48;
    for (Int_t e = 0; e <=fEtaBins; e++) {
      x[e] = -6. + e*(12./(Double_t)fEtaBins);
    }
  }
 
  // Reference flow binning
  Double_t xR[101] = { 0. };
  for (Int_t r = 0; r <= fEtaRef; r++) {
    xR[r] = -6 + r*(12./(Double_t)fEtaRef);
  }
 
  // Phi binning
  Double_t phi[21] =  { 0. };
  for (Int_t p = 0; p <= 20; p++) {
    phi[p] = p*2.*TMath::Pi() / 20.;
  }
  Double_t phiMC[201] = { 0. };
  for (Int_t p = 0; p <= 200; p++) {
    phiMC[p] = p*2.*TMath::Pi() / 200.;
  }

  // Histograms for cumulants analysis

  // We loop over flow histograms here to add different orders of harmonics
  TString type = "";
  for (Int_t loop = 1; loop <= 5; loop++) {
    
    if (loop == 1) type = "FMD";
    if (loop == 2) type = "SPD";
    if (loop == 3) type = "MC";
    if (loop == 4) type = "FMDTR";
    if (loop == 5) type = "SPDTR";

    TList* tList = new TList();
    tList->SetName(type.Data());
    fOutputList->Add(tList);

    for (Int_t n = 1; n <= 6; n++) {
      if (!fv[n]) continue;
      
      TList* vList = new TList();
      vList->SetName(Form("v%d", n));
      tList->Add(vList);
      
      TString tag = TString();
      for (Int_t c = 0; c < 100; c++) {
        // Output histograms
        tag = Form("%d_%d", c, c+1);
          
        TH3D* hFlowHist = new TH3D(Form("hQ%dCumuHist%s_%s", n, type.Data(), tag.Data()), 
                   Form("hQ%dCumuHist%s_%s", n, type.Data(), tag.Data()), fEtaBins, -6, 6, 10, -5, 5, 26, 0.5, 26.5);
        TAxis* xAxis = hFlowHist->GetXaxis();
        xAxis->Set(fEtaBins, x);
        hFlowHist->Sumw2();
        vList->Add(hFlowHist);

        TProfile* hCumulant2DiffFlow = new TProfile(Form("hQ%dCumulant2DiffFlow%s_%s", n, type.Data(), tag.Data()), 
                     Form("hQ%dCumulant2DiffFlow%s_%s", n, type.Data(), tag.Data()), fEtaBins, x);
        hCumulant2DiffFlow->Sumw2();
        vList->Add(hCumulant2DiffFlow);
   
        TProfile* hCumulant4DiffFlow = new TProfile(Form("hQ%dCumulant4DiffFlow%s_%s", n, type.Data(), tag.Data()), 
                     Form("hQ%dCumulant4DiffFlow%s_%s", n, type.Data(), tag.Data()), fEtaBins, x);
        hCumulant4DiffFlow->Sumw2();
        vList->Add(hCumulant4DiffFlow);
      } // end of centrality loop

    } // end of v_{n} loop
 
    // Single Event histograms
    TH2D* hdNdphiSE = new TH2D(Form("hdNdphiSE%s", type.Data()),
                 Form("hdNdphiSE%s", type.Data()), fEtaRef, xR, (type.Contains("MC") ? 200 : 20), (type.Contains("MC") ? phiMC : phi));
    hdNdphiSE->Sumw2();
    fOutputList->Add(hdNdphiSE);

    TH2D* hdNdetadphiSE = new TH2D(Form("hdNdetadphiSE%s", type.Data()), 
                 Form("hdNdetadphiSE%s", type.Data()), fEtaBins, x, (type.Contains("MC") ? 200 : 20), (type.Contains("MC") ? phiMC : phi));
    hdNdetadphiSE->Sumw2();
    fOutputList->Add(hdNdetadphiSE);
  } // end of type loop

  TProfile2D* pMCTruth = new TProfile2D("pMCTruth", "pMCTruth", 48, -6, 6, 100, 0, 100);
  TAxis* xAxis = pMCTruth->GetXaxis();
  xAxis->Set(fEtaBins, x);
  pMCTruth->Sumw2();
  fOutputList->Add(pMCTruth);

  // Monitoring plots
  TH1D* cent = new TH1D("Centralities", "Centralities", 100, 0, 100);
  fOutputList->Add(cent);

  TH1D* vertexSel = new TH1D("VertexSelected", "VertexSelected", 50, -10, 10);
  fOutputList->Add(vertexSel);
  TH1D* vertexAll = new TH1D("VertexAll", "VertexAll", 50, -10, 10);
  fOutputList->Add(vertexAll);
  
  TH2D* vertex2D = new TH2D("CoverageVsVertex", "CoverageVsVertex", fEtaBins, -6, 6, 20, -10, 10);
  fOutputList->Add(vertex2D);
 
  PostData(1, fOutputList);
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::UserExec(Option_t */*option*/)
{
  // 
  // Process each event
  //
  // Parameters:
  //  option: Not used
  //

  // Reset data members
  fCent = -1;
  fZvertex = 1111;

  // Get input event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) return;

  // Fill histograms
  if (!fFlowUtil->LoopAODFMD(fAOD)) return;
  if (!fFlowUtil->LoopAODSPD(fAOD)) return;

  // Get centrality and vertex from flow utility and fill monitoring histograms
  fCent = fFlowUtil->GetCentrality();
  fZvertex = fFlowUtil->GetVertex();

  TH1D* cent = (TH1D*)fOutputList->FindObject("Centralities");

  cent->Fill(fCent);
  TH1D* vertex = (TH1D*)fOutputList->FindObject("VertexSelected");
  vertex->Fill(fZvertex);

  // Run analysis on FMD and SPD
  for (Int_t n = 1; n <= 6; n++) {
    if (fv[n]) {
      CumulantsMethod("FMD", n);
      CumulantsMethod("SPD", n);
    }
  }

  // And do analysis if there is.
  if (fMC) { 
    ProcessPrimary();
  }

  PostData(1, fOutputList);
}
//_____________________________________________________________________
void AliForwardFlowTaskQC::CumulantsMethod(TString type = "", 
					   Int_t harmonic = 2)
{
  // 
  // Calculate the Q cumulant of order n
  //
  // Parameters:
  //  type: Determines which histograms should be used
  //        - "FMD" or "SPD" = data histograms
  //        - "FMDTR" or "SPDTR" = track reference histograms
  //        - "MC" = MC truth histograms
  //  harmonic: Which harmonic to calculate
  //
  Double_t n = harmonic;
  
  TList* tList = (TList*)fOutputList->FindObject(type.Data());
  TList* vList = (TList*)tList->FindObject(Form("v%d", harmonic));
  // We get the histograms
  TH3D* flowHist   = (TH3D*)vList->FindObject(Form("hQ%dCumuHist%s_%d_%d", harmonic, type.Data(), (Int_t)fCent, (Int_t)fCent+1));
  TH2D* dNdetadphi = (TH2D*)fOutputList->FindObject(Form("hdNdetadphiSE%s", type.Data()));
  TH2D* dNdphi = (TH2D*)fOutputList->FindObject(Form("hdNdphiSE%s", /*"FMD"*/type.Data()));
  
  // We create the coordinate array used to fill the THnSparse, centrality and vertex is set from the beginning.
  Double_t coord[4] = { 0., fCent, fZvertex, 0. };
  // We create the objects needed for the analysis
  Double_t dQnRe = 0, dQ2nRe = 0, dQnIm = 0, dQ2nIm = 0, mult = 0;
  Double_t pnRe = 0, p2nRe = 0, qnRe = 0, q2nRe = 0, pnIm = 0, p2nIm = 0, qnIm = 0, q2nIm = 0;
  Double_t two = 0, four = 0, twoPrime = 0, fourPrime = 0;
  Double_t cosPhi1Phi2 = 0, cosPhi1Phi2Phi3m = 0;
  Double_t sinPhi1Phi2 = 0, sinPhi1Phi2Phi3m = 0;
  Double_t cosPsi1Phi2 = 0, cosPsi1Phi2Phi3m = 0, cosPsi1Phi2Phi3p = 0;
  Double_t sinPsi1Phi2 = 0, sinPsi1Phi2Phi3m = 0, sinPsi1Phi2Phi3p = 0;
  Double_t phi = 0, eta = 0;
  Double_t multi = 0, multp = 0, mp = 0, mq = 0;
  Double_t w2 = 0, w4 = 0, w2p = 0, w4p = 0;
  Int_t refEtaBin = 0;
  Bool_t kEventCount = kFALSE;

  // We loop over the data 1 time!
  for (Int_t etaBin = 1; etaBin <= dNdetadphi->GetNbinsX(); etaBin++) {
    eta = dNdetadphi->GetXaxis()->GetBinCenter(etaBin);
    refEtaBin = dNdphi->GetXaxis()->FindBin(eta);
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

    for (Int_t phiBin = 1; phiBin <= dNdphi->GetNbinsY(); phiBin++) {
      phi = dNdphi->GetYaxis()->GetBinCenter(phiBin);
      
      // Reference flow
      multi = dNdphi->GetBinContent(refEtaBin, phiBin);
      dQnRe += multi*TMath::Cos(n*phi);
      dQnIm += multi*TMath::Sin(n*phi);
      dQ2nRe += multi*TMath::Cos(2.*n*phi);
      dQ2nIm += multi*TMath::Sin(2.*n*phi);
      mult += multi;
      
      // For each etaBin bin the necessary values for differential flow
      // is calculated. Here is the loop over the phi's.
      multp = dNdetadphi->GetBinContent(etaBin, phiBin);
      mp += multp;
      pnRe += multp*TMath::Cos(n*phi);
      pnIm += multp*TMath::Sin(n*phi);
      p2nRe += multp*TMath::Cos(2.*n*phi);
      p2nIm += multp*TMath::Sin(2.*n*phi);
    }
    if (mult == 0) continue; 

    if (!kEventCount) {
     // Count number of events
      coord[0] = -7;
      coord[3] = -0.5;
      flowHist->Fill(coord[0], coord[2], coord[3], 1.);
      kEventCount = kTRUE;
    } 

    // The reference flow is calculated 
    coord[0] = eta;
    
    // 2-particle
    w2 = mult * (mult - 1.);
    two = dQnRe*dQnRe + dQnIm*dQnIm - mult;
    coord[3] = kW2Two;
    flowHist->Fill(coord[0], coord[2], coord[3], two);
    coord[3] = kW2;
    flowHist->Fill(coord[0], coord[2], coord[3], w2);

    coord[3] = kQnRe;
    flowHist->Fill(coord[0], coord[2], coord[3], dQnRe);
    coord[3] = kQnIm;
    flowHist->Fill(coord[0], coord[2], coord[3], dQnIm);
    coord[3] = kM;
    flowHist->Fill(coord[0], coord[2], coord[3], mult);
    
    // 4-particle
    w4 = mult * (mult - 1.) * (mult - 2.) * (mult - 3.);
  
    four = 2.*mult*(mult-3.) + TMath::Power((TMath::Power(dQnRe,2.)+TMath::Power(dQnIm,2.)),2.)
             -4.*(mult-2.)*(TMath::Power(dQnRe,2.) + TMath::Power(dQnIm,2.))
             -2.*(TMath::Power(dQnRe,2.)*dQ2nRe+2.*dQnRe*dQnIm*dQ2nIm-TMath::Power(dQnIm,2.)*dQ2nRe)
             +(TMath::Power(dQ2nRe,2.)+TMath::Power(dQ2nIm,2.));

    coord[3] = kW4Four;
    flowHist->Fill(coord[0], coord[2], coord[3], four);
    coord[3] = kW4;
    flowHist->Fill(coord[0], coord[2], coord[3], w4);

    cosPhi1Phi2 = dQnRe*dQnRe - dQnIm*dQnIm - dQ2nRe;
    sinPhi1Phi2 = 2.*dQnRe*dQnIm - dQ2nIm;
      
    cosPhi1Phi2Phi3m = dQnRe*(TMath::Power(dQnRe,2)+TMath::Power(dQnIm,2))-dQnRe*dQ2nRe-dQnIm*dQ2nIm-2.*(mult-1)*dQnRe;

    sinPhi1Phi2Phi3m = -dQnIm*(TMath::Power(dQnRe,2)+TMath::Power(dQnIm,2))+dQnRe*dQ2nIm-dQnIm*dQ2nRe+2.*(mult-1)*dQnIm; 

    coord[3] = kCosphi1phi2;
    flowHist->Fill(coord[0], coord[2], coord[3], cosPhi1Phi2);
    coord[3] = kSinphi1phi2;
    flowHist->Fill(coord[0], coord[2], coord[3], sinPhi1Phi2);
    coord[3] = kCosphi1phi2phi3m;
    flowHist->Fill(coord[0], coord[2], coord[3], cosPhi1Phi2Phi3m);
    coord[3] = kSinphi1phi2phi3m;
    flowHist->Fill(coord[0], coord[2], coord[3], sinPhi1Phi2Phi3m);
    coord[3] = kMm1m2;
    flowHist->Fill(coord[0], coord[2], coord[3], mult*(mult-1.)*(mult-2.));

    // Differential flow calculations for each etaBin bin is done:
    mq = mp;
    qnRe = pnRe;
    qnIm = pnIm;
    q2nRe = p2nRe;
    q2nIm = p2nIm;

    // 2-particle differential flow
    w2p = mp * mult - mq;
    twoPrime = pnRe*dQnRe + pnIm*dQnIm - mq;
    
    coord[3] = kw2two;
    flowHist->Fill(coord[0], coord[2], coord[3], twoPrime);
    coord[3] = kw2;
    flowHist->Fill(coord[0], coord[2], coord[3], w2p);

    coord[3] = kpnRe;
    flowHist->Fill(coord[0], coord[2], coord[3], pnRe);
    coord[3] = kpnIm;
    flowHist->Fill(coord[0], coord[2], coord[3], pnIm);
    coord[3] = kmp;
    flowHist->Fill(coord[0], coord[2], coord[3], mp);

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
   
    coord[3] = kw4four;
    flowHist->Fill(coord[0], coord[2], coord[3], fourPrime);
    coord[3] = kw4;
    flowHist->Fill(coord[0], coord[2], coord[3], w4p);

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

    coord[3] = kCospsi1phi2;
    flowHist->Fill(coord[0], coord[2], coord[3], cosPsi1Phi2);
    coord[3] = kSinpsi1phi2;
    flowHist->Fill(coord[0], coord[2], coord[3], sinPsi1Phi2);
    coord[3] = kCospsi1phi2phi3m;
    flowHist->Fill(coord[0], coord[2], coord[3], cosPsi1Phi2Phi3m);
    coord[3] = kSinpsi1phi2phi3m;
    flowHist->Fill(coord[0], coord[2], coord[3], sinPsi1Phi2Phi3m);
    coord[3] = kmpmq;
    flowHist->Fill(coord[0], coord[2], coord[3], (mp*mult-2.*mq)*(mult-1.));
    coord[3] = kCospsi1phi2phi3p;
    flowHist->Fill(coord[0], coord[2], coord[3], cosPsi1Phi2Phi3p);
    coord[3] = kSinpsi1phi2phi3p;
    flowHist->Fill(coord[0], coord[2], coord[3], sinPsi1Phi2Phi3p); 

  }

}
//_____________________________________________________________________
void AliForwardFlowTaskQC::Terminate(Option_t */*option*/) 
{
  // 
  //  End of job
  //  Finalizes the Q cumulant calculations
  // 
  //  Parameters:
  //    option Not used 
  //

  TH3D* cumulantsHist = 0;
  TProfile* cumulant2diffHist = 0;
  TProfile* cumulant4diffHist = 0;
  TList* tList = 0;
  TList* vList = 0;

  // For flow calculations
  Double_t two = 0, qc2 = 0, /* vnTwo = 0, */ four = 0, qc4 = 0 /*, vnFour = 0*/; 
  Double_t twoPrime = 0, qc2Prime = 0, vnTwoDiff = 0, fourPrime = 0, qc4Prime = 0, vnFourDiff = 0;
  Double_t w2 = 0, w4 = 0, w2p = 0, w4p = 0;
  Double_t w2Two = 0, w2pTwoPrime = 0, w4Four = 0, w4pFourPrime = 0;
  Double_t cosP1nPhi = 0, sinP1nPhi = 0, mult = 0, cosP1nPhi1P1nPhi2 = 0, sinP1nPhi1P1nPhi2 = 0;
  Double_t cosP1nPhi1M1nPhi2M1nPhi3 = 0, sinP1nPhi1M1nPhi2M1nPhi3 = 0, multm1m2 = 0;
  Double_t cosP1nPsi = 0, sinP1nPsi = 0, mp = 0, cosP1nPsi1P1nPhi2 = 0, sinP1nPsi1P1nPhi2 = 0;
  Double_t cosP1nPsi1M1nPhi2M1nPhi3 = 0, sinP1nPsi1M1nPhi2M1nPhi3 = 0, mpqMult = 0;
  Double_t cosP1nPsi1P1nPhi2M1nPhi3 = 0, sinP1nPsi1P1nPhi2M1nPhi3 = 0;

  Int_t nLoops = (fMC ? 5 : 2);
  TString type = "";

  // Do a loop over the difference analysis types, calculating flow
  // 2 loops for real data, 3 for MC data
  // inside each is a nested loop over each harmonic (1, 2, 3 and 4 at the moment)
  for (Int_t loop = 1; loop <= nLoops; loop++) {

    if (loop == 1) type = "FMD";
    if (loop == 2) type = "SPD";
    if (loop == 3) type = "MC";
    if (loop == 4) type = "FMDTR";
    if (loop == 5) type = "SPDTR";
    
    for (Int_t n = 1; n <= 6; n++) {
      if (!fv[n]) continue;
     
      tList = (TList*)fOutputList->FindObject(type.Data());
      vList = (TList*)tList->FindObject(Form("v%d", n));
      TString tag = "";
 
      // Centrality loop
      for (Int_t c = 0; c < 100; c++) {

        Double_t nEv = 0;
        tag = Form("%d_%d", c, c+1);
        cumulantsHist = (TH3D*)vList->FindObject(Form("hQ%dCumuHist%s_%s", n, type.Data(), tag.Data()));
        
        cumulant2diffHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant2DiffFlow%s_%s", n, type.Data(), tag.Data()));
        cumulant4diffHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant4DiffFlow%s_%s", n, type.Data(), tag.Data()));
        
        for (Int_t vertexBin = 1; vertexBin <= cumulantsHist->GetNbinsY(); vertexBin++) {
        
          for (Int_t etaBin = 1; etaBin <= cumulantsHist->GetNbinsX(); etaBin++) {
            Double_t eta = cumulantsHist->GetXaxis()->GetBinCenter(etaBin);
            // 2-particle reference flow
            w2Two = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2Two);
            if (!w2Two) continue;
            w2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2);
            cosP1nPhi = cumulantsHist->GetBinContent(etaBin, vertexBin, kQnRe);
            sinP1nPhi = cumulantsHist->GetBinContent(etaBin, vertexBin, kQnIm);
            mult = cumulantsHist->GetBinContent(etaBin, vertexBin, kM);
              
            cosP1nPhi /= mult;
            sinP1nPhi /= mult;
            two = w2Two / w2;
            qc2 = two - TMath::Power(cosP1nPhi, 2) - TMath::Power(sinP1nPhi, 2);
            if (qc2 <= 0) continue;
            // vnTwo = TMath::Sqrt(qc2);
       //     if (!TMath::IsNaN(vnTwo*mult)) 
       //       cumulant2diffHist->Fill(eta, vnTwo, cumulantsHist->GetBinContent(0,vertexBin,0)); 

            // 4-particle reference flow
            w4Four = cumulantsHist->GetBinContent(etaBin, vertexBin, kW4Four);
            w4 = cumulantsHist->GetBinContent(etaBin, vertexBin, kW4);
            cosP1nPhi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kCosphi1phi2);
            sinP1nPhi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kSinphi1phi2);
            cosP1nPhi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kCosphi1phi2phi3m);
            sinP1nPhi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kSinphi1phi2phi3m);
            multm1m2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kMm1m2);

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
            // vnFour = TMath::Power(-qc4, 0.25);
       //     if (!TMath::IsNaN(vnFour*mult)) 
       //         cumulant4diffHist->Fill(eta, vnFour, cumulantsHist->GetBinContent(0,vertexBin,0));
  
            // 2-particle differential flow
            w2pTwoPrime = cumulantsHist->GetBinContent(etaBin, vertexBin, kw2two);
            if (!w2pTwoPrime) continue;
            w2p = cumulantsHist->GetBinContent(etaBin, vertexBin, kw2);
            cosP1nPsi = cumulantsHist->GetBinContent(etaBin, vertexBin, kpnRe);
            sinP1nPsi = cumulantsHist->GetBinContent(etaBin, vertexBin, kpnIm);
            mp = cumulantsHist->GetBinContent(etaBin, vertexBin, kmp);

            cosP1nPsi /= mp;
            sinP1nPsi /= mp;
            twoPrime = w2pTwoPrime / w2p;
            qc2Prime = twoPrime - sinP1nPsi*sinP1nPhi - cosP1nPsi*cosP1nPhi;

            vnTwoDiff = qc2Prime / TMath::Sqrt(qc2);
            if (!TMath::IsNaN(vnTwoDiff*mp) && vnTwoDiff > 0) cumulant2diffHist->Fill(eta, vnTwoDiff, cumulantsHist->GetBinContent(0,vertexBin,0));

            // 4-particle differential flow
            w4pFourPrime = cumulantsHist->GetBinContent(etaBin, vertexBin, kw4four);
            w4p = cumulantsHist->GetBinContent(etaBin, vertexBin, kw4);
            cosP1nPsi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kCospsi1phi2);
            sinP1nPsi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kSinpsi1phi2);
            cosP1nPsi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kCospsi1phi2phi3m);
            sinP1nPsi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kSinpsi1phi2phi3m);
            cosP1nPsi1P1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kCospsi1phi2phi3p);
            sinP1nPsi1P1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kSinpsi1phi2phi3p); 
            mpqMult = cumulantsHist->GetBinContent(etaBin, vertexBin, kmpmq);

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
            if (!TMath::IsNaN(vnFourDiff*mp) && vnFourDiff > 0) cumulant4diffHist->Fill(eta, vnFourDiff, cumulantsHist->GetBinContent(0,vertexBin,0));
          } // End of etaBin loop
          nEv += cumulantsHist->GetBinContent(0,vertexBin,0);
        } // End of vertexBin loop
        // Number of events:
        cumulant2diffHist->Fill(7., nEv);
        cumulant4diffHist->Fill(7., nEv);
  
      } // End of centrality loop
    } // End of harmonics loop
  }  // End of type loop

  PostData(1, fOutputList);
  
  return;
}
// _____________________________________________________________________
void AliForwardFlowTaskQC::ProcessPrimary() 
{
  //
  // If fMC == kTRUE this function takes care of organizing the input 
  // Monte Carlo data and histograms so AliForwardFlowTaskQC::QCumulants 
  // can be run on it.
  //

  if (fFlowUtil->LoopAODFMDtrrefHits(fAOD)) {
   // Run analysis on TrackRefs
    for (Int_t n = 1; n <= 6; n++) {
      if (fv[n])
        CumulantsMethod("FMDTR", n);
    }
  }

  if (fFlowUtil->LoopAODSPDtrrefHits(fAOD)) {
    // Run analysis on SPD TrackRefs
    for (Int_t n = 1; n <= 6; n++) {
      if (fv[n])
        CumulantsMethod("SPDTR", n);
    }
  }

  if (fFlowUtil->LoopAODmc(fAOD, fAddFlow, fAddType, fAddOrder)) {
    // Run analysis on MC truth
    for (Int_t n = 1; n <= 6; n++) {
      if (fv[n])
        CumulantsMethod("MC", n);
    }
  }

}
//_____________________________________________________________________
//
//
// EOF
