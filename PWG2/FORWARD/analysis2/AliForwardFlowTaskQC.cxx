//
// Calculate flow in the forward regions using the Q cumulants method
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
enum { kW2Two = 1, kW2, kW2Twosq, kW2w2pTwoTwoPrime, kW2w2p, kW4Four,
       kW4, kW4Foursq, kW2w4, kW2w4TwoFour, kW2w4p, kW2w4pTwoFourPrime, 
       kW4w2p, kW4w2pFourTwoPrime, kW4w4p, kW4w4pFourFourPrime, kQnRe, kQnIm,
       kM, kCosphi1phi2, kSinphi1phi2, kCosphi1phi2phi3m, kSinphi1phi2phi3m,
       kMm1m2, kCosphi1phi2phi3p, kSinphi1phi2phi3p,
       kRW2Two, kRW2, kRW2Twosq, kRM, kRW4Four, kRW4, kRW4Foursq,
       kRW2w4, kRW2w4TwoFour, kRCosphi1phi2, kRSinphi1phi2, 
       kRCosphi1phi2phi3m, kRSinphi1phi2phi3m, kRMm1m2 };

ClassImp(AliForwardFlowTaskQC)
#if 0
; // For emacs 
#endif

AliForwardFlowTaskQC::AliForwardFlowTaskQC()
  : fOutputList(0),	// Output list
    fFlowUtil(0),	// AliForwardFlowUtil
    fAOD(0),		// AOD input event
    fMC(kFALSE),	// MC flag
    fEtaBins(20),	// # of etaBin bins in histograms
    fAddFlow(0),	// Add flow string
    fAddType(0),	// Add flow type #
    fAddOrder(0),	// Add flow order
    fZvertex(0),	// Z vertex range
    fCent(0)		// Centrality
{
  // 
  // Default constructor
  //
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const char* name) :
  AliAnalysisTaskSE(name),
  fOutputList(0),	// Output list
  fFlowUtil(0),		// AliForwardFlowUtil
  fAOD(0),		// AOD input event
  fMC(kFALSE),		// MC flag
  fEtaBins(20),		// # of Eta bins
  fAddFlow(0),		// Add flow string
  fAddType(0),		// Add flow type #
  fAddOrder(0),		// Add flow order
  fZvertex(0),		// Z vertex range
  fCent(0)		// Centrality
{
  // 
  // Constructor
  //
  // Parameters:
  //  name: Name of task
  //
  for (Int_t n = 0; n <= 4; n++) fv[n] = kTRUE;
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
  for (Int_t n = 0; n <= 4; n++) fv[n] = o.fv[n];
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

  if (fEtaBins % 6) fEtaBins = 48;

  Double_t x[49] = { 0. };
  for (Int_t e = 0; e <=fEtaBins; e++) {
    x[e] = -6. + e*(12./(Double_t)fEtaBins);
  }
//      Double_t x[6] = { 0.0, 1.0, 2.0, 3.0, 4.5, 6.0 };
  Double_t y[11] = { 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 100. };
  Int_t bins[4] = { fEtaBins, 10, 20, 40 };
  Double_t xmin[4] = { -6., 0., -10., 0.5 };
  Double_t xmax[4] = { 6., 100., 10., 40.5 };


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

    for (Int_t n = 1; n <= 4; n++) {
      if (!fv[n]) continue;
      
      TList* vList = new TList();
      vList->SetName(Form("v%d", n));
      tList->Add(vList);
      
      // Only one flow histogram is needed for each type of data;
      // z-axis bin 1:  (w_<2> * <2>).Re()
      // z-axis bin 2:  w_<2> = mqM-mp
      // z-axis bin 3:  (w_<2> * <2> * <2>).Re()
      // z-axis bin 4:  (w_<2> * w_<2'> * <2> * <2'>).Re()
      // z-axis bin 5:  w_<2> * w_<2'>
      // z-axis bin 6:  (w_<4> * <4>).Re()
      // z-axis bin 7:  w_<4>
      // z-axis bin 8:  (w_<4> * <4> * <4>).Re()
      // z-axis bin 9:  w_<2> * w_<4>
      // z-axis bin 10:  (w_<2> * w_<4> * <2> * <4>).Re()
      // z-axis bin 11:  w_<2> * w_<4'>
      // z-axis bin 12:  (w_<2> * w_<4'> * <2> * <4'>).Re()
      // z-axis bin 13:  w_<4> * w_<2'>
      // z-axis bin 14:  (w_<4> * w_<2'> * <4> * <2'>).Re()
      // z-axis bin 15:  w_<4> * w_<4'>
      // z-axis bin 16:  (w_<4>  * w_<4'> * <4> * <4'>).Re()
      // z-axis bin 17: Qn or pn.Re() = <<cos2phi or psi>>
      // z-axis bin 18: Qn or pn.Im() = <<sin2phi or psi>>
      // z-axis bin 19: M or mp
      // z-axis bin 20: (Qn*Qn-Q2n).Re() = <<cos(2(phi1 or psi1+phi2))>>
      // z-axis bin 21: (Qn*Qn-Q2n).Im() = <<sin(2(phi1 or psi1+phi2))>>
      // z-axis bin 22: <<cos(2(phi1 or psi1-phi2-phi3))>>
      // z-axis bin 23: <<sin(2(phi1 or psi1-phi2-phi3))>>
      // z-axis bin 24: M*(M-1)*(M-2) or similar for diff
      // z-axis bin 25: <<cos(2(psi1+phi2-phi3>>
      // z-axis bin 26: <<sin(2(psi1+phi2-phi3>>
      // z-axis bin 27: ref: W_<2> * <2>
      // z-axis bin 28: ref: W_<2>
      // z-axis bin 29: ref: W_<2> * <2> * <2>
      // z-axis bin 30: ref: Mult
      // z-axis bin 31: ref: W_<4> * <4>
      // z-axis bin 32: ref: W_<4>
      // z-axis bin 33: ref: W_<4> * <4> * <4>
      // z-axis bin 34: ref: W_<2> * W_<4>
      // z-axis bin 35: ref: W_<2> * W_<4> * <2> * <4>
      // z-axis bin 36: ref: <cos(n(phi1phi2))>
      // z-axis bin 37: ref: <sin(n(phi1phi2))>
      // z-axis bin 38: ref: <cos(n(phi1-phi2-phi3))>
      // z-axis bin 39: ref: <sin(n(phi1-phi2-phi3))>
      // z-axis bin 40: ref: M(M-1)(M-2)
      
      THnSparseD* hFlowHist = new THnSparseD(Form("hQ%dCumuHist%s", n, type.Data()), 
                   Form("hQ%dCumuHist%s", n, type.Data()), 4, bins, xmin, xmax);
      TAxis* yAxis = hFlowHist->GetAxis(1);
      yAxis->Set(10, y);
      hFlowHist->Sumw2();
      vList->Add(hFlowHist);
      
      TString tag = TString();
      for (Int_t c = -2; c < 10; c++) {
        // Output histograms
        if (c == -2)      tag = "mb";
        else if (c == -1) tag = "0_40";
        else              tag = Form("%d_%d", (Int_t)y[c], (Int_t)y[c+1]);
    
        TProfile* hCumulant2RefFlow = new TProfile(Form("hQ%dCumulant2RefFlow%s_%s", n, type.Data(), tag.Data()), 
                     Form("hQ%dCumulant2RefFlow%s_%s", n, type.Data(), tag.Data()), fEtaBins, x);
        hCumulant2RefFlow->Sumw2();
        vList->Add(hCumulant2RefFlow);
   
        TProfile* hCumulant4RefFlow = new TProfile(Form("hQ%dCumulant4RefFlow%s_%s", n, type.Data(), tag.Data()), 
                     Form("hQ%dCumulant4RefFlow%s_%s", n, type.Data(), tag.Data()), fEtaBins, x);
        hCumulant4RefFlow->Sumw2();
        vList->Add(hCumulant4RefFlow);
 
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
    TH1D* hdNdphiSE = new TH1D(Form("hdNdphiSE%s", type.Data()),
                 Form("hdNdphiSE%s", type.Data()), 20, 0, 2*TMath::Pi());
    hdNdphiSE->Sumw2();
    fOutputList->Add(hdNdphiSE);

    TH2D* hdNdetadphiSE = new TH2D(Form("hdNdetadphiSE%s", type.Data()), 
                 Form("hdNdetadphiSE%s", type.Data()), fEtaBins, -6, 6, 20, 0, 2*TMath::Pi());
    hdNdetadphiSE->Sumw2();
    fOutputList->Add(hdNdetadphiSE);
  } // end of type loop

  TProfile2D* pMCTruth = new TProfile2D("pMCTruth", "pMCTruth", fEtaBins, x, 10, y);
  pMCTruth->Sumw2();
  fOutputList->Add(pMCTruth);

  TH1D* cent = new TH1D("Centralities", "Centralities", 100, 0, 100);
  fOutputList->Add(cent);

  TH2D* vertex = new TH2D("CoverageVsVertex", "CoverageVsVertex", fEtaBins, -6, 6, 20, -10, 10);
  fOutputList->Add(vertex);
 
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

  // Get input event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) return;

  // Fill histograms
  if (!fFlowUtil->LoopAODFMD(fAOD)) return;

  fCent = fFlowUtil->GetCentrality();
  fZvertex = fFlowUtil->GetVertex();
  TH1D* cent = (TH1D*)fOutputList->FindObject("Centralities");
  cent->Fill(fCent);
 
  // Run analysis on FMD
  for (Int_t n = 1; n <= 4; n++) {
    if (fv[n])
      CumulantsMethod("FMD", n);
  }
  
  // Fill SPD hist
  if (!fFlowUtil->LoopAODSPD(fAOD)) return;

  // Run analysis on SPD
  for (Int_t n = 1; n <= 4; n++) {
    if (fv[n])
      CumulantsMethod("SPD", n);
  }
 
  // Find out if there's any MC data present
  if (!fMC) {
    TClonesArray* mcArray = 0;
    mcArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
    if (mcArray) fMC = kTRUE;
  }
  if (fMC) 
    ProcessPrimary();

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
  //        - "" = data histograms
  //        - "TrRef" = track reference histograms
  //        - "MC" = MC truth histograms
  //  harmonic: Which harmonic to calculate
  //
  Double_t n = harmonic;
  
  TList* tList = (TList*)fOutputList->FindObject(type.Data());
  TList* vList = (TList*)tList->FindObject(Form("v%d", harmonic));
  // We get histograms depending on if it's real data or MC truth data
  THnSparseD* flowHist   = (THnSparseD*)vList->FindObject(Form("hQ%dCumuHist%s", harmonic, type.Data()));
  TH2D* dNdetadphi = (TH2D*)fOutputList->FindObject(Form("hdNdetadphiSE%s", type.Data()));
  TString type2 = type;
  if (type.Contains("SPD")) type2 = type2.ReplaceAll("SPD", "FMD"); 
  TH1D* dNdphi = (TH1D*)fOutputList->FindObject(Form("hdNdphiSE%s", type2.Data()));
  
  // We create the objects needed for the analysis
  Double_t mult = dNdphi->GetBinContent(0);
  if (mult <= 3) return;

  Double_t coord[4] = { 0., 0., 0., 0. };
  Double_t dQnRe = 0, dQ2nRe = 0, dQnIm = 0, dQ2nIm = 0;
  Double_t pnRe = 0, p2nRe = 0, qnRe = 0, q2nRe = 0, pnIm = 0, p2nIm = 0, qnIm = 0, q2nIm = 0;
  Double_t two = 0, four = 0, twoPrime = 0, fourPrime = 0;
//  Double_t w2Twosq = 0, w2pTwoPrimesq = 0, w4Foursq = 0, w4pFourPrimesq = 0;
//  Double_t w2w2pTwoTwoPrime = 0, w2w4TwoFour = 0, w2pw4pTwoPrimeFourPrime = 0;
//  Double_t w2w4pTwoFourPrime = 0, w4w2pFourTwoPrime = 0, w4w4pFourFourPrime = 0;
  Double_t cosPhi1Phi2 = 0, cosPhi1Phi2Phi3m = 0, cosPhi1Phi2Phi3p = 0;
  Double_t sinPhi1Phi2 = 0, sinPhi1Phi2Phi3m = 0, sinPhi1Phi2Phi3p = 0;
  Double_t phi = 0;
  Double_t multi = 0, multp = 0, mp = 0, mq = 0;
  Double_t w2 = 0, w4 = 0, w2p = 0, w4p = 0;

  // We loop over the data 1 time!
  for (Int_t etaBin = 1; etaBin <= dNdetadphi->GetNbinsX(); etaBin++) {
    // The values for each individual etaBin bins are reset
    mp = 0;
    pnRe = 0;
    p2nRe = 0;
    pnIm = 0;
    p2nIm = 0;

    for (Int_t phiN = 1; phiN <= dNdphi->GetNbinsX(); phiN++) {
      phi = dNdphi->GetXaxis()->GetBinCenter(phiN);
      multi = dNdphi->GetBinContent(phiN);

      // In the phi loop on the first etaBin loop the integrated flow
      // is calculated from the dNdphi histogram
      if(etaBin == 1) {
        dQnRe += multi*TMath::Cos(n*phi);
        dQnIm += multi*TMath::Sin(n*phi);
        dQ2nRe += multi*TMath::Cos(2.*n*phi);
        dQ2nIm += multi*TMath::Sin(2.*n*phi);
      }
      
      // For each etaBin bin the necesarry values for differential flow
      // is calculated. Here is the loop over the phi's.
      multp = dNdetadphi->GetBinContent(etaBin, phiN);
      mp += multp;
      pnRe += multp*TMath::Cos(n*phi);
      pnIm += multp*TMath::Sin(n*phi);
      p2nRe += multp*TMath::Cos(2.*n*phi);
      p2nIm += multp*TMath::Sin(2.*n*phi);    
    }

    // The reference flow is calculated for the differential flow histograms
    if (etaBin == 1) {
      Double_t eta = flowHist->GetAxis(0)->GetBinCenter(0);
      coord[0] = eta;
      coord[1] = fCent;
      coord[2] = fZvertex;

      // 2-particle
      w2 = mult * (mult - 1.);
      two = dQnRe*dQnRe + dQnIm*dQnIm - mult;
      two /= w2;
//      w2Twosq = w2 * two * two; 

      coord[3] = kRW2Two;
      flowHist->Fill(coord, w2 * two);
      coord[3] = kRW2;
      flowHist->Fill(coord, w2);
//      coord[3] = kRW2Twosq;
//      flowHist->Fill(coord, w2Twosq);

      coord[3] = kQnRe;
      flowHist->Fill(coord, dQnRe);
      coord[3] = kQnIm;
      flowHist->Fill(coord, dQnIm);
      coord[3] = kRM;
      flowHist->Fill(coord, mult);

      // 4-particle
      w4 = mult * (mult - 1.) * (mult - 2.) * (mult - 3.);
      Double_t real = dQ2nRe*dQnRe*dQnRe - dQ2nRe*dQnIm*dQnIm + 2.*dQ2nIm*dQnRe*dQnIm;

      four = TMath::Power(dQnRe*dQnRe + dQnIm*dQnIm, 2); 
      four += dQ2nRe*dQ2nRe + dQ2nIm*dQ2nIm - 2.*real;
      four -= 4.*(mult - 2.)*(dQnRe*dQnRe + dQnIm*dQnIm) - 2.*mult*(mult - 3.);  
      four /= w4;

//      w4Foursq = w4 * four * four;
//      w2w4TwoFour = w2 * w4 * two * four;

      coord[3] = kRW4Four;
      flowHist->Fill(coord, w4 * four);
      coord[3] = kRW4;
      flowHist->Fill(coord, w4);
/*      coord[3] = kRW4Foursq;
      flowHist->Fill(coord, w4Foursq);
      coord[3] = kRW2w4;
      flowHist->Fill(coord, w2 * w4);
      coord[3] = kRW2w4TwoFour;
      flowHist->Fill(coord, w2w4TwoFour);
*/
      cosPhi1Phi2 = dQnRe*dQnRe - dQnIm*dQnIm - dQ2nRe;
      sinPhi1Phi2 = 2.*dQnRe*dQnIm - dQ2nIm;
      
      cosPhi1Phi2Phi3m = TMath::Power(dQnRe, 3) + dQnRe*dQnIm*dQnIm; 
      cosPhi1Phi2Phi3m -= dQnRe*dQ2nRe + dQnIm*dQ2nIm + 2.*(mult - 1.)*dQnRe;

      sinPhi1Phi2Phi3m = -TMath::Power(dQnIm, 3) - dQnRe*dQnRe*dQnIm; 
      sinPhi1Phi2Phi3m -= dQnIm*dQ2nRe - dQnRe*dQ2nIm - 2.*(mult - 1.)*dQnIm;

      coord[3] = kRCosphi1phi2;
      flowHist->Fill(coord, cosPhi1Phi2);
      coord[3] = kRSinphi1phi2;
      flowHist->Fill(coord, sinPhi1Phi2);
      coord[3] = kRCosphi1phi2phi3m;
      flowHist->Fill(coord, cosPhi1Phi2Phi3m);
      coord[3] = kRSinphi1phi2phi3m;
      flowHist->Fill(coord, sinPhi1Phi2Phi3m);
      coord[3] = kRMm1m2;
      flowHist->Fill(coord, mult*(mult-1.)*(mult-2.));

      // Count number of events
      coord[3] = 0.;
      flowHist->Fill(coord, 1.);
    } 

    // Differential flow calculations for each etaBin bin is done:
    if (mp <= 3) continue;
    Double_t eta = dNdetadphi->GetXaxis()->GetBinCenter(etaBin);
//    eta = TMath::Abs(eta);
    coord[0] = eta;

    mq = mp;
    qnRe = pnRe;
    qnIm = pnIm;
    q2nRe = p2nRe;
    q2nIm = p2nIm;
    if (type.Contains("SPD") || (type.Contains("FMD") && eta >= 4.)) {
      mq = 0;
      qnRe = 0;
      qnIm = 0;
      q2nRe = 0;
      q2nIm = 0;
    }

    // Then the reference flow is calculated for each etaBin bin also,
    // TODO: Find smart way to implement in above calculations...

    // 2-particle
    w2p = mp * (mp - 1.);
    twoPrime = pnRe*pnRe + pnIm*pnIm - mp;
    twoPrime /= w2p;
//    w2pTwoPrimesq = w2p * twoPrime * twoPrime; 

    coord[3] = kRW2Two;
    flowHist->Fill(coord, w2p * twoPrime);
    coord[3] = kRW2;
    flowHist->Fill(coord, w2p);
//    coord[3] = kRW2Twosq;
//    flowHist->Fill(coord, w2pTwoPrimesq);

    coord[3] = kRM;
    flowHist->Fill(coord, mp);

    // 4-particle
    w4p = mp * (mp - 1.) * (mp - 2.) * (mp - 3.);
    Double_t real = p2nRe*pnRe*pnRe - p2nRe*pnIm*pnIm + 2.*p2nIm*pnRe*pnIm;

    fourPrime = TMath::Power(pnRe*pnRe + pnIm*pnIm, 2); 
    fourPrime += p2nRe*p2nRe + p2nIm*p2nIm - 2.*real;
    fourPrime -= 4.*(mp - 2.)*(pnRe*pnRe + pnIm*pnIm) - 2.*mp*(mp - 3.);
  
    fourPrime /= w4p;
 //   w4pFourPrimesq = w4p * fourPrime * fourPrime;
 //   w2w4pTwoFourPrime = w2p * w4p * twoPrime * fourPrime;

    coord[3] = kRW4Four;
    flowHist->Fill(coord, w4p * fourPrime);
    coord[3] = kRW4;
    flowHist->Fill(coord, w4p);
/*    coord[3] = kRW4Foursq;
    flowHist->Fill(coord, w4pFourPrimesq);
    coord[3] = kRW2w4;
    flowHist->Fill(coord, w2p * w4p);
    coord[3] = kRW2w4TwoFour;
    flowHist->Fill(coord, w2w4pTwoFourPrime);
*/
    cosPhi1Phi2 = pnRe*pnRe - pnIm*pnIm - p2nRe;
    sinPhi1Phi2 = 2.*pnRe*pnIm - p2nIm;
      
    cosPhi1Phi2Phi3m = TMath::Power(pnRe, 3) + pnRe*pnIm*pnIm; 
    cosPhi1Phi2Phi3m -= pnRe*p2nRe + pnIm*p2nIm + 2.*(mp - 1.)*pnRe;

    sinPhi1Phi2Phi3m = -TMath::Power(pnIm, 3) - pnRe*pnRe*pnIm; 
    sinPhi1Phi2Phi3m -= pnIm*p2nRe - pnRe*p2nIm - 2.*(mp - 1.)*pnIm;

    coord[3] = kRCosphi1phi2;
    flowHist->Fill(coord, cosPhi1Phi2);
    coord[3] = kRSinphi1phi2;
    flowHist->Fill(coord, sinPhi1Phi2);
    coord[3] = kRCosphi1phi2phi3m;
    flowHist->Fill(coord, cosPhi1Phi2Phi3m);
    coord[3] = kRSinphi1phi2phi3m;
    flowHist->Fill(coord, sinPhi1Phi2Phi3m);
    coord[3] = kRMm1m2;
    flowHist->Fill(coord, mp*(mp-1.)*(mp-2.));

    // 2-particle differential flow
    w2p = mp * mult - mq;
    twoPrime = pnRe*dQnRe + pnIm*dQnIm - mq;
    twoPrime /= w2p;
//    w2pTwoPrimesq = w2p * twoPrime * twoPrime;
//    w2w2pTwoTwoPrime = w2 * w2p * two * twoPrime;
    
    coord[3] = kW2Two;
    flowHist->Fill(coord, w2p * twoPrime);
    coord[3] = kW2;
    flowHist->Fill(coord, w2p);
/*    coord[3] = kW2Twosq;
    flowHist->Fill(coord, w2pTwoPrimesq);
    coord[3] = kW2w2pTwoTwoPrime;
    flowHist->Fill(coord, w2w2pTwoTwoPrime);
    coord[3] = kW2w2p;
    flowHist->Fill(coord, w2 * w2p);
*/
    coord[3] = kQnRe;
    flowHist->Fill(coord, pnRe);
    coord[3] = kQnIm;
    flowHist->Fill(coord, pnIm);
    coord[3] = kM;
    flowHist->Fill(coord, mp);

    // 4-particle differential flow
    w4p = (mp * mult - 3.*mq)*(mult - 1.)*(mult - 2.);

    fourPrime =  pnRe*dQnRe*(dQnRe*dQnRe + dQnIm*dQnIm) + pnIm*dQnIm*(dQnRe*dQnRe + dQnIm*dQnIm);
    fourPrime -= q2nRe*dQnRe*dQnRe - q2nRe*dQnIm*dQnIm + 2.*q2nIm*dQnRe*dQnIm;
    fourPrime -= pnRe*dQnRe*dQ2nRe - pnRe*dQnIm*dQ2nIm + pnIm*dQnRe*dQ2nIm + pnIm*dQnIm*dQ2nRe;
    fourPrime -= 2.*mult*(pnRe*dQnRe + pnIm*dQnIm);

    fourPrime += - 2.*mq*(dQnRe*dQnRe + dQnIm*dQnIm) + 7.*(qnRe*dQnRe + qnIm*dQnIm); 
    fourPrime += - (dQnRe*qnRe + dQnIm*qnIm) + (q2nRe*dQ2nRe + q2nIm*dQ2nIm);
    fourPrime += 2.*(pnRe*dQnRe + pnIm*dQnIm) + 2.*mq*mult - 6.*mq;
    fourPrime /= w4p;
/*
    w4pFourPrimesq = w4p * fourPrime * fourPrime;
    w2w4pTwoFourPrime = w2 * w4p * two * fourPrime;
    w4w2pFourTwoPrime = w4 * w2p * four * twoPrime;
    w4w4pFourFourPrime = w4 * w4p * four * fourPrime;
    w2pw4pTwoPrimeFourPrime = w2p * w4p * twoPrime * fourPrime;
*/
    coord[3] = kW4Four;
    flowHist->Fill(coord, w4p * fourPrime);
    coord[3] = kW4;
    flowHist->Fill(coord, w4p);
/*    coord[3] = kW4Foursq;
    flowHist->Fill(coord, w4pFourPrimesq);
    coord[3] = kW2w4;
    flowHist->Fill(coord, w2p * w4p);
    coord[3] = kW2w4TwoFour;
    flowHist->Fill(coord, w2pw4pTwoPrimeFourPrime);
    coord[3] = kW2w4p;
    flowHist->Fill(coord, w2 * w4p);
    coord[3] = kW2w4pTwoFourPrime;
    flowHist->Fill(coord, w2w4pTwoFourPrime);
    coord[3] = kW4w2p;
    flowHist->Fill(coord, w4 * w2p);
    coord[3] = kW4w2pFourTwoPrime;
    flowHist->Fill(coord, w4w2pFourTwoPrime);
    coord[3] = kW4w4p;
    flowHist->Fill(coord, w4 * w4p);
    coord[3] = kW4w4pFourFourPrime;
    flowHist->Fill(coord, w4w4pFourFourPrime);
*/
    cosPhi1Phi2 = pnRe*dQnRe - pnIm*dQnIm - q2nRe;
    sinPhi1Phi2 = pnRe*dQnIm + pnIm*dQnRe - q2nIm;

    cosPhi1Phi2Phi3p =  pnRe*(dQnRe*dQnRe + dQnIm*dQnIm - mult);
    cosPhi1Phi2Phi3p -= q2nRe*dQnRe - q2nIm*dQnIm + mq*dQnRe - 2.*qnRe;
    sinPhi1Phi2Phi3p =  pnIm*(dQnRe*dQnRe + dQnIm*dQnIm - mult);
    sinPhi1Phi2Phi3p -= q2nIm*dQnRe - q2nRe*dQnIm + mq*dQnIm - 2.*qnIm;

    cosPhi1Phi2Phi3m =  pnRe*(dQnRe*dQnRe - dQnIm*dQnIm) + 2.*pnIm*dQnRe*dQnIm;
    cosPhi1Phi2Phi3m -= pnRe*dQ2nRe + pnIm*dQ2nIm + 2.*mq*dQnRe - 2.*qnRe;
    sinPhi1Phi2Phi3m =  pnIm*(dQnRe*dQnRe - dQnIm*dQnIm) - 2.*pnRe*dQnRe*dQnIm;
    sinPhi1Phi2Phi3m += - pnIm*dQ2nRe + pnRe*dQ2nIm + 2.*mq*dQnIm - 2.*qnIm;

    coord[3] = kCosphi1phi2;
    flowHist->Fill(coord, cosPhi1Phi2);
    coord[3] = kSinphi1phi2;
    flowHist->Fill(coord, sinPhi1Phi2);
    coord[3] = kCosphi1phi2phi3m;
    flowHist->Fill(coord, cosPhi1Phi2Phi3m);
    coord[3] = kSinphi1phi2phi3m;
    flowHist->Fill(coord, sinPhi1Phi2Phi3m);
    coord[3] = kMm1m2;
    flowHist->Fill(coord, (mp*mult-2.*mq)*(mult-1.));
    coord[3] = kCosphi1phi2phi3p;
    flowHist->Fill(coord, cosPhi1Phi2Phi3p);
    coord[3] = kSinphi1phi2phi3p;
    flowHist->Fill(coord, sinPhi1Phi2Phi3p); 

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

//  Double_t x[6] = { 0., 1., 2., 3., 4.5, 6. };

  THnSparseD* hcumulantsHist = 0;
  TProfile* cumulant2refHist = 0;
  TProfile* cumulant4refHist = 0;
  TProfile* cumulant2diffHist = 0;
  TProfile* cumulant4diffHist = 0;
  TList* tList = 0;
  TList* vList = 0;

  // For flow calculations
  Double_t two = 0, qc2 = 0, vnTwo = 0, four = 0, qc4 = 0, vnFour = 0; 
  Double_t twoPrime = 0, qc2Prime = 0, vnTwoDiff = 0, fourPrime = 0, qc4Prime = 0, vnFourDiff = 0;
  Double_t w2 = 0, w4 = 0, w2p = 0, w4p = 0;
//  Double_t sqrtW2sq = 0, sqrtW2psq = 0, w2w2p = 0;
//  Double_t w2w4 = 0, w2w4p = 0, w4w2p = 0, w4w4p = 0, w2pW4p = 0;
//  Double_t sqrtW4sq = 0, sqrtW4psq = 0;
  Double_t w2Two = 0, w2pTwoPrime = 0, w4Four = 0, w4pFourPrime = 0;
  Double_t cosP1nPhi = 0, sinP1nPhi = 0, mult = 0, cosP1nPhi1P1nPhi2 = 0, sinP1nPhi1P1nPhi2 = 0;
  Double_t cosP1nPhi1M1nPhi2M1nPhi3 = 0, sinP1nPhi1M1nPhi2M1nPhi3 = 0, multm1m2 = 0;
  Double_t cosP1nPsi = 0, sinP1nPsi = 0, mp = 0, cosP1nPsi1P1nPhi2 = 0, sinP1nPsi1P1nPhi2 = 0;
  Double_t cosP1nPsi1M1nPhi2M1nPhi3 = 0, sinP1nPsi1M1nPhi2M1nPhi3 = 0, mpqMult = 0;
  Double_t cosP1nPsi1P1nPhi2M1nPhi3 = 0, sinP1nPsi1P1nPhi2M1nPhi3 = 0;
/*
  // For error calculations
  Double_t w2Twosq = 0, w2w2pTwoTwoPrime = 0, w2pTwoPrimesq = 0;
  Double_t w4Foursq = 0, w2w4TwoFour = 0, w4w4pFourFourPrime = 0, w4pFourPrimesq = 0;
  Double_t w2w4pTwoFourPrime = 0, w4w2pFourTwoPrime = 0;
  Double_t twoErrorSq = 0, twoPrimeErrorSq = 0, fourErrorSq = 0, fourPrimeErrorSq = 0;
  Double_t vnTwoError = 0, vnTwoDiffErr = 0, vnFourError = 0, vnFourDiffErr = 0;
  Double_t cov22p = 0, cov24 = 0, cov24p = 0, cov42p = 0, cov44p = 0, cov2p2p = 0;
  Double_t w2pW4pTwoPrimeFourPrime = 0;
*/

  Int_t nLoops = (fMC ? 5 : 2);
  TString type = "";

  Int_t y[11] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100 };
  
  // Do a loop over the difference analysis types, calculating flow
  // 2 loops for real data, 3 for MC data
  // inside each is a nested loop over each harmonic (1, 2, 3 and 4 at the moment)
  for (Int_t loop = 1; loop <= nLoops; loop++) {

    if (loop == 1) type = "FMD";
    if (loop == 2) type = "SPD";
    if (loop == 3) type = "MC";
    if (loop == 4) type = "FMDTR";
    if (loop == 5) type = "SPDTR";
    
    for (Int_t n = 1; n <= 4; n++) {
      if (!fv[n]) continue;
     
      tList = (TList*)fOutputList->FindObject(type.Data());
      vList = (TList*)tList->FindObject(Form("v%d", n));

      hcumulantsHist = (THnSparseD*)vList->FindObject(Form("hQ%dCumuHist%s", n, type.Data()));
 
      // Centrality loop
      for (Int_t c = -2; c < hcumulantsHist->GetAxis(1)->GetNbins(); c++) {
        // Test
//        Double_t multInt[41] = { 0 }, vnTwoInt[41] = { 0 }, vnFourInt[41] = { 0 };
//        Double_t vnTwoErrInt[41] = { 0 }, vnFourErrInt[41] = { 0 };
//        Double_t multDiffInt[41] = { 0 }, vnTwoDiffInt[41] = { 0 }, vnFourDiffInt[41] = { 0 };
//        Double_t vnTwoDiffErrInt[41] = { 0 }, vnFourDiffErrInt[41] = { 0 };

        Double_t nEv = 0;
        TH3D* cumulantsHist = 0;
        if (c == -2) {
          hcumulantsHist->GetAxis(1)->SetRange(1, 9);
          cumulantsHist = (TH3D*)hcumulantsHist->Projection(0, 2, 3, "E");
          
          cumulant2refHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant2RefFlow%s_mb", n, type.Data()));
          cumulant4refHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant4RefFlow%s_mb", n, type.Data()));
          cumulant2diffHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant2DiffFlow%s_mb", n, type.Data()));
          cumulant4diffHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant4DiffFlow%s_mb", n, type.Data()));
        }
        else if (c == -1) {
          hcumulantsHist->GetAxis(1)->SetRange(1, 5);
          cumulantsHist = (TH3D*)hcumulantsHist->Projection(0, 2, 3, "E");
          
          cumulant2refHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant2RefFlow%s_0_40", n, type.Data()));
          cumulant4refHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant4RefFlow%s_0_40", n, type.Data()));
          cumulant2diffHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant2DiffFlow%s_0_40", n, type.Data()));
          cumulant4diffHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant4DiffFlow%s_0_40", n, type.Data()));
        }
        else {
          hcumulantsHist->GetAxis(1)->SetRange(c+1, c+1);
          cumulantsHist = (TH3D*)hcumulantsHist->Projection(0, 2, 3, "E");
          cumulant2refHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant2RefFlow%s_%d_%d", n, type.Data(), y[c], y[c+1]));
          cumulant4refHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant4RefFlow%s_%d_%d", n, type.Data(), y[c], y[c+1]));
          cumulant2diffHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant2DiffFlow%s_%d_%d", n, type.Data(), y[c], y[c+1]));
          cumulant4diffHist = (TProfile*)vList->FindObject(Form("hQ%dCumulant4DiffFlow%s_%d_%d", n, type.Data(), y[c], y[c+1]));
        }
        
        for (Int_t vertexBin = 1; vertexBin <= cumulantsHist->GetNbinsY(); vertexBin++) {
        
          Bool_t refDone = kFALSE;
          for (Int_t etaBin = 1; etaBin <= cumulantsHist->GetNbinsX(); etaBin++) {
          Double_t eta = cumulantsHist->GetXaxis()->GetBinCenter(etaBin);
            if (refDone == kFALSE || etaBin == 0) {
              // 2-particle reference flow
              w2Two = cumulantsHist->GetBinContent(etaBin, vertexBin, kRW2Two);
              if (etaBin == cumulantsHist->GetNbinsX()) {
                refDone = kTRUE;
                etaBin = -1;
                continue;
              }
              if (!w2Two) continue;
              w2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kRW2);
              cosP1nPhi = cumulantsHist->GetBinContent(etaBin, vertexBin, kQnRe);
              sinP1nPhi = cumulantsHist->GetBinContent(etaBin, vertexBin, kQnIm);
              mult = cumulantsHist->GetBinContent(etaBin, vertexBin, kRM);
/*
              if (!mult || !w2) {
                AliWarning(Form("mult = %f and w2 = %f at eta bin %d in histogram for v%d of type %s_%d_%d"
                   , mult, w2, etaBin, n, type.Data(), y[c], y[c+1]));
                continue;
              }
*/
              cosP1nPhi /= mult;
              sinP1nPhi /= mult;
              two = w2Two / w2;
              qc2 = two  - TMath::Power(cosP1nPhi, 2) - TMath::Power(sinP1nPhi, 2);
/*            
              if (qc2<=0) {
                AliWarning(Form("Negative or zero qc2(%f) at eta bin %d in histogram for v%d of type %s_%d_%d"
                   , qc2, etaBin, n, type.Data(), y[c], y[c+1]));
                continue;
              }
*/              
              vnTwo = TMath::Sqrt(qc2);
              if (!TMath::IsNaN(vnTwo*mult)) {
                if (etaBin == 0) cumulant2diffHist->Fill(eta, vnTwo, cumulantsHist->GetBinContent(0,vertexBin,0)); 
                else cumulant2refHist->Fill(eta, vnTwo, cumulantsHist->GetBinContent(0,vertexBin,0));
              }
//              if (!TMath::IsNaN(mult))  multInt[etaBin] += mult;

              // 4-particle reference flow
              w4Four = cumulantsHist->GetBinContent(etaBin, vertexBin, kRW4Four);
              w4 = cumulantsHist->GetBinContent(etaBin, vertexBin, kRW4);
              cosP1nPhi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kRCosphi1phi2);
              sinP1nPhi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kRSinphi1phi2);
              cosP1nPhi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kRCosphi1phi2phi3m);
              sinP1nPhi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kRSinphi1phi2phi3m);
              multm1m2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kRMm1m2);
/*
              if (!w4 || !multm1m2) {
                AliWarning(Form("w4 = %f and  multm1m2 = %f at eta bin %d in histogram for v%d of type %s_%d_%d"
                   , w4, multm1m2, etaBin, n, type.Data(), y[c], y[c+1]));
                continue;
              }
*/
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
/*
              if (qc4>=0) {
                AliWarning(Form("Positive or zero qc4(%f) at eta bin %d in histogram for v%d of type %s_%d_%d"
                   , qc4, etaBin, n, type.Data(), y[c], y[c+1]));
                continue;
              }
*/
              vnFour = TMath::Power(-qc4, 0.25);
              if (!TMath::IsNaN(vnFour*mult)) {
                if (etaBin == 0) cumulant4diffHist->Fill(eta, vnFour, cumulantsHist->GetBinContent(0,vertexBin,0)); 
                else cumulant4refHist->Fill(eta, vnFour, cumulantsHist->GetBinContent(0,vertexBin,0));
              }
/*
              // 2-particle reference flow error calculations
              w2Twosq = cumulantsHist->GetBinContent(etaBin, vertexBin, kRW2Twosq);
              sqrtW2sq = cumulantsHist->GetBinError(etaBin, vertexBin, kRW2);
     
              twoErrorSq = VarSQ(w2Twosq, two, w2, w2Two, sqrtW2sq);
              
              if (two<=0) {
                AliWarning(Form("two = %f at eta bin %d in histogram for v%d of type %s_%d_%d"
                   , two, etaBin, n, type.Data(), y[c], y[c+1]));
                continue;
              }
 
              vnTwoError = sqrtW2sq * TMath::Sqrt(twoErrorSq) / (2. * TMath::Sqrt(two) * w2);
              if (!TMath::IsNaN(vnTwoError)) vnTwoErrInt[etaBin] += vnTwoError*vnTwoError*mult*mult;
      
              // 4-particle reference flow error calculations
              w4Foursq = cumulantsHist->GetBinContent(etaBin, vertexBin, kRW4Foursq);
              sqrtW4sq = cumulantsHist->GetBinError(etaBin, vertexBin, kRW4);
              w2w4 = cumulantsHist->GetBinContent(etaBin, vertexBin, kRW2w4);
              w2w4TwoFour = cumulantsHist->GetBinContent(etaBin, vertexBin, kRW2w4TwoFour);
   
              fourErrorSq = VarSQ(w4Foursq, four, w4, w4Four, sqrtW4sq);
              cov24 = CovXY(w2w4TwoFour, w2w4, two*four, w2, w4);
  
              vnFourError = two*two * TMath::Power(sqrtW2sq, 2) * twoErrorSq / (w2*w2)
                            + TMath::Power(sqrtW4sq, 2) * fourErrorSq / (16. * w4*w4)
                            - two * w2w4 * cov24 / (2. * w2 * w4);
              vnFourError /= TMath::Power(2. * two*two - four, 1.5);
            
              if (vnFourError<=0) {
                AliWarning(Form("vnFourError = %f at eta bin %d in histogram for v%d of type %s_%d_%d"
                   , vnFourError, etaBin, n, type.Data(), y[c], y[c+1]));
                continue;
              }

              vnFourError = TMath::Sqrt(vnFourError);
              if (!TMath::IsNaN(vnFourError)) vnFourErrInt[etaBin] += vnFourError*vnFourError*mult*mult;
*/ 
              continue; // Skip differential flow until all reference flow is done.
            }
  
            // 2-particle differential flow
            w2pTwoPrime = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2Two);
            if (!w2pTwoPrime) continue;
            w2p = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2);
            cosP1nPsi = cumulantsHist->GetBinContent(etaBin, vertexBin, kQnRe);
            sinP1nPsi = cumulantsHist->GetBinContent(etaBin, vertexBin, kQnIm);
            mp = cumulantsHist->GetBinContent(etaBin, vertexBin, kM);
/*  
            if (!mp || !w2p) {
              AliWarning(Form("mp = %f and w2p = %f at eta bin %d in histogram for v%d of type %s_%d_%d"
                 , mp, w2p, etaBin, n, type.Data(), y[c], y[c+1]));
              continue;
            }
*/
            cosP1nPsi /= mp;
            sinP1nPsi /= mp;
            twoPrime = w2pTwoPrime / w2p;
            qc2Prime = twoPrime - sinP1nPsi*sinP1nPhi - cosP1nPsi*cosP1nPhi;

            vnTwoDiff = qc2Prime / TMath::Sqrt(qc2);
            if (!TMath::IsNaN(vnTwoDiff*mp)) cumulant2diffHist->Fill(eta, vnTwoDiff, cumulantsHist->GetBinContent(0,vertexBin,0));
//            if (!TMath::IsNaN(mp))        multDiffInt[etaBin] += mp;

            // 4-particle differential flow
            w4pFourPrime = cumulantsHist->GetBinContent(etaBin, vertexBin, kW4Four);
            w4p = cumulantsHist->GetBinContent(etaBin, vertexBin, kW4);
            cosP1nPsi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kCosphi1phi2);
            sinP1nPsi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, vertexBin, kSinphi1phi2);
            cosP1nPsi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kCosphi1phi2phi3m);
            sinP1nPsi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kSinphi1phi2phi3m);
            cosP1nPsi1P1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kCosphi1phi2phi3p);
            sinP1nPsi1P1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, vertexBin, kSinphi1phi2phi3p);
          
            mpqMult = cumulantsHist->GetBinContent(etaBin, vertexBin, kMm1m2);
            cosP1nPsi1P1nPhi2 /= w2p;
            sinP1nPsi1P1nPhi2 /= w2p;
            cosP1nPsi1M1nPhi2M1nPhi3 /= mpqMult;
            sinP1nPsi1M1nPhi2M1nPhi3 /= mpqMult;
            cosP1nPsi1P1nPhi2M1nPhi3 /= mpqMult;
            sinP1nPsi1P1nPhi2M1nPhi3 /= mpqMult;
            fourPrime = w4pFourPrime / w4p;
/*
            if (!mpqMult || !w4p) {
              AliWarning(Form("mpqMult = %f and w4p = %f at eta bin %d in histogram for v%d of type %s_%d_%d"
                 , mpqMult, w4p, etaBin, n, type.Data(), y[c], y[c+1]));
              continue;
            }
*/
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
            if (!TMath::IsNaN(vnFourDiff*mp)) cumulant4diffHist->Fill(eta, vnFourDiff, cumulantsHist->GetBinContent(0,vertexBin,0));
/*
            // 2-particle differential flow error calculations
            w2pTwoPrimesq = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2Twosq);
            sqrtW2psq = cumulantsHist->GetBinError(etaBin, vertexBin, kW2);
            w2w2pTwoTwoPrime = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2w2pTwoTwoPrime);
            w2w2p = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2w2p);
      
            cov22p = CovXY(w2w2pTwoTwoPrime, w2w2p, two*twoPrime, w2, w2p);
            twoPrimeErrorSq = VarSQ(w2pTwoPrimesq, twoPrime, w2p, w2pTwoPrime, sqrtW2psq);
     
            vnTwoDiffErr = twoPrime*twoPrime*TMath::Power(sqrtW2sq, 2)*twoErrorSq/(w2*w2)
                          + 4.*two*two*TMath::Power(sqrtW2psq, 2)*twoPrimeErrorSq/(w2p*w2p)
                          - 4.*two*twoPrime*w2w2p*cov22p/(w2*w2p);
  
            vnTwoDiffErr /= (4. * TMath::Power(two, 3));
  
            if (vnTwoDiffErr<0) {
              AliWarning(Form("Negative vnTwoDiffErr(%f) at eta bin %d in histogram for v%d of type %s_%d_%d"
                 , vnTwoDiffErr, etaBin, n, type.Data(), y[c], y[c+1]));
              continue;
            }
  
            vnTwoDiffErr = TMath::Sqrt(vnTwoDiffErr);
            if (!TMath::IsNaN(vnTwoDiffErr)) vnTwoDiffErrInt[etaBin] += vnTwoDiffErr*vnTwoDiffErr*mp*mp;
    
            // 4-particle differential flow error calculations
            sqrtW4psq = cumulantsHist->GetBinError(etaBin, vertexBin, kW4);
            w4pFourPrimesq = cumulantsHist->GetBinContent(etaBin, vertexBin, kW4Foursq);
            w2pW4p = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2w4);
            w2pW4pTwoPrimeFourPrime = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2w4TwoFour);
            w2w4p = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2w4p);
            w2w4pTwoFourPrime = cumulantsHist->GetBinContent(etaBin, vertexBin, kW2w4pTwoFourPrime);
            w4w2p = cumulantsHist->GetBinContent(etaBin, vertexBin, kW4w2p);
            w4w2pFourTwoPrime = cumulantsHist->GetBinContent(etaBin, vertexBin, kW4w2pFourTwoPrime);
            w4w4p = cumulantsHist->GetBinContent(etaBin, vertexBin, kW4w4p);
            w4w4pFourFourPrime = cumulantsHist->GetBinContent(etaBin, vertexBin, kW4w4pFourFourPrime);
         
            fourPrimeErrorSq = VarSQ(w4pFourPrimesq, fourPrime, w4p, w4pFourPrime, sqrtW4psq);
            cov24p = CovXY(w2w4pTwoFourPrime, w2w4p, two*fourPrime, w2, w4p);
            cov42p = CovXY(w4w2pFourTwoPrime, w4w2p, four*twoPrime, w4, w2p);
            cov44p = CovXY(w4w4pFourFourPrime, w4w4p, four*fourPrime, w4, w4p);
            cov2p2p = CovXY(w2pW4pTwoPrimeFourPrime, w2pW4p, twoPrime*fourPrime, w2p, w4p);
    
            // Numbers on the side reference term number in paper (cite needed) loosely
            vnFourDiffErr =  TMath::Power(2.*two*two*twoPrime - 3.*two*fourPrime + 2.*four*twoPrime, 2) 
                            * TMath::Power(sqrtW2sq, 2) * twoErrorSq / (w2*w2)
                           + 9. * TMath::Power(2.*two*twoPrime - fourPrime, 2) * TMath::Power(sqrtW4sq, 2)
                            * fourErrorSq / (16. * w4*w4)
                           + 4. * two*two * TMath::Power(2.*two*two - four, 2) * TMath::Power(sqrtW2psq, 2)
                            * twoPrimeErrorSq / (w2p*w2p)
                           + TMath::Power(2.*two*two - four, 2) * TMath::Power(sqrtW4psq, 2) 
                            * fourPrimeErrorSq
                            / (w4p*w4p)
                           - 1.5 * (2.*two*twoPrime - fourPrime) 
                            * (2.*two*two*twoPrime - 3.*two*fourPrime + 2.*four*twoPrime)
                            * w2w4 * cov24 / (w2*w4)
                           - 4. * two * (2.*two*two - four) 
                            * (2.*two*two*twoPrime - 3.*two*fourPrime + 2.*four*twoPrime)
                            * w2w2p * cov22p / (w2 * w2p)
                           + 2. * (2.*two*two - four)
                            * (2.*two*two*twoPrime - 3.*two*fourPrime + 2.*four*twoPrime)
                            * w2w4p * cov24p / (w2 * w4p)
                           + 3.*two*(2.*two*two - four)*(2.*two*twoPrime - fourPrime)
                            * w4w2p * cov42p / (w4*w2p)
                           - 1.5 * (2.*two*two - four)*(2.*two*twoPrime - fourPrime)
                            * w4w4p * cov44p / (w4 * w4p)
                           - 4.*two*TMath::Power(2.*two*two - four, 2)
                            * w2pW4p * cov2p2p / (w2p * w4p);
  
            if (2*two*two-four == 0) {
              AliWarning(Form("Dividing by zero in vnFourDiffError at eta bin %d in histogram for v%d of type %s_%d_%d"
                 , etaBin, n, type.Data(), y[c], y[c+1]));
              continue;
            }  

           vnFourDiffErr /= TMath::Power(2.*two*two - four, 3.5);
  
           if (vnFourDiffErr<0) {
              AliWarning(Form("Negative vnFourDiffErr(%f) at eta bin %d in histogram for v%d of type %s_%d_%d"
                 , vnFourDiffErr, etaBin, n, type.Data(), y[c], y[c+1]));
              continue;
            }

            vnFourDiffErr = TMath::Sqrt(vnFourDiffErr);
            if (!TMath::IsNaN(vnFourDiffErr)) vnFourDiffErrInt[etaBin] += vnFourDiffErr*vnFourDiffErr*mp*mp; 
*/
          } // End of etaBin loop
          nEv += cumulantsHist->GetBinContent(0,vertexBin,0);
        } // End of vertexBin loop
/*
        for (Int_t etaBin = 0; etaBin <= fEtaBins; etaBin++) {
          if (!multInt[etaBin]) continue;

          vnTwo = vnTwoInt[etaBin] / multInt[etaBin];
          if (etaBin == 0) cumulant2diffHist->SetBinContent(etaBin, vnTwo);
          if (etaBin > 0)  cumulant2refHist->SetBinContent(etaBin, vnTwo);
            
          vnFour = vnFourInt[etaBin] / multInt[etaBin];
          if (etaBin == 0) cumulant4diffHist->SetBinContent(etaBin, vnFour);
          if (etaBin > 0) cumulant4refHist->SetBinContent(etaBin, vnFour);
           
          vnTwoError = TMath::Sqrt(vnTwoErrInt[etaBin]) / multInt[etaBin];
          if (etaBin == 0) cumulant2diffHist->SetBinError(etaBin, vnTwoError);
          if (etaBin > 0)  cumulant2refHist->SetBinError(etaBin, vnTwoError);
          
          vnFourError = TMath::Sqrt(vnFourErrInt[etaBin]) / multInt[etaBin];
          if (etaBin == 0) cumulant4diffHist->SetBinError(etaBin, vnFourError);
          if (etaBin > 0)  cumulant4refHist->SetBinError(etaBin, vnFourError);
           
          if (!multDiffInt[etaBin] || !etaBin) continue;
          vnTwoDiff = vnTwoDiffInt[etaBin] / multDiffInt[etaBin];
          cumulant2diffHist->SetBinContent(etaBin, vnTwoDiff);

          vnFourDiff = vnFourDiffInt[etaBin] / multDiffInt[etaBin];
          cumulant4diffHist->SetBinContent(etaBin, vnFourDiff);

          vnTwoDiffErr = TMath::Sqrt(vnTwoDiffErrInt[etaBin]) / multDiffInt[etaBin];
          cumulant2diffHist->SetBinError(etaBin, vnTwoDiffErr);

          vnFourDiffErr = TMath::Sqrt(vnFourDiffErrInt[etaBin]) / multDiffInt[etaBin];
          cumulant4diffHist->SetBinError(etaBin, vnFourDiffErr);
        }
*/
        // Number of events:
        cumulant2refHist->Fill(7., nEv);
        cumulant4refHist->Fill(7., nEv);
        cumulant2diffHist->Fill(7., nEv);
        cumulant4diffHist->Fill(7., nEv);
  
        delete cumulantsHist;

      } // End of centrality loop
    } // End of harmonics loop
  fOutputList->Remove(hcumulantsHist);
  }  // End of type loop

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
    for (Int_t n = 1; n <= 4; n++) {
      if (fv[n])
        CumulantsMethod("FMDTR", n);
    }
  }

  if (fFlowUtil->LoopAODSPDtrrefHits(fAOD)) {
    // Run analysis on SPD TrackRefs
    for (Int_t n = 1; n <= 4; n++) {
      if (fv[n])
        CumulantsMethod("SPDTR", n);
    }
  }

  if (fFlowUtil->LoopAODmc(fAOD, fAddFlow, fAddType, fAddOrder)) {
    // Run analysis on MC truth
    for (Int_t n = 1; n <= 4; n++) {
      if (fv[n])
        CumulantsMethod("MC", n);
    }
  }

 }
//_____________________________________________________________________
Double_t AliForwardFlowTaskQC::VarSQ(Double_t wxx2, 
				     Double_t x, 
				     Double_t wx, 
				     Double_t wxx, 
				     Double_t sqrtwx2) const
{
  //
  // Small function to compute the variance squared - used by Terminte()
  //
  Double_t sx;
  
  sx = wxx2 + x*x*wx - 2.*x*wxx;
  sx /= wx;
  sx /= (1. - TMath::Power(sqrtwx2, 2) / (wx*wx));

  return sx;
}
//_____________________________________________________________________
Double_t AliForwardFlowTaskQC::CovXY(Double_t wxwyxy, 
				     Double_t wxwy, 
				     Double_t xy, 
				     Double_t wx, 
				     Double_t wy) const
{
  //
  // Small function to compute the covariance between two numbers
  // - used by Terminate()
  //
  Double_t cov, denominator, numerator;

  denominator = (wxwyxy / wxwy) - xy;
  numerator = 1. - (wxwy / (wx * wy));
/* 
  if (numerator == 0) {
    AliWarning("A numerator just gave 0 in CovXY");
    return 1.;
  }
*/
  cov = denominator / numerator;
  return cov;
}
//_____________________________________________________________________
//
//
// EOF
