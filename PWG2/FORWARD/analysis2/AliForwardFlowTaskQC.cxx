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
#include "TH3D.h"
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
enum { kW2avg2 = 1, kW2, kW2avg2sq, kW2w2pavg2avg2p, kW2w2p, kW4avg4,
       kW4, kW4avg4sq, kW2w4, kW2w4avg2avg4, kW2w4p, kW2w4pavg2avg4p, 
       kW4w2p, kW4w2pavg4avg2p, kW4w4p, kW4w4pavg4avg4p, kQnRe, kQnIm,
       kM, kCosphi1phi2, kSinphi1phi2, kCosphi1phi2phi3m, kSinphi1phi2phi3m,
       kMm1m2, kCosphi1phi2phi3p, kSinphi1phi2phi3p,
       kRW2avg2, kRW2, kRW2avg2sq, kRM, kRW4avg4, kRW4, kRW4avg4sq,
       kRW2w4, kRW2w4avg2avg4, kRCosphi1phi2, kRSinphi1phi2, 
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
  fCent(0)			// Centrality
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
void AliForwardFlowTaskQC::CreateOutputObjects()
{
  //
  // Create output objects
  // 
  if (!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("QCumulants");
  fOutputList->SetOwner();
  fFlowUtil = new AliForwardFlowUtil(fOutputList);
  fFlowUtil->SetVertexRange(fZvertex);

  if (fEtaBins % 20) fEtaBins = 20;

  // Histograms for cumulants analysis

  // We loop over flow histograms here to add different orders of harmonics
  TString type = "";
  for (Int_t loop = 1; loop <= 4; loop++) {
    
    if (loop == 1) type = "";
    if (loop == 2) type = "SPD";
    if (loop == 3) type = "MC";
    if (loop == 4) type = "TrRef";

    for (Int_t n = 1; n <= 4; n++) {
      if (!fv[n]) continue;
      // Only one flow histogram is needed for each type of data;
      // x-axis is etaBin-bins with differential flow, integrated is in underflowbin
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
      Double_t x[41] = { 0. };
      for (Int_t e = 0; e <=fEtaBins; e++) {
        x[e] = -4. + e*(10./(Double_t)fEtaBins);
      }
//      Double_t x[6] = { 0.0, 1.0, 2.0, 3.0, 4.5, 6.0 };
      Double_t y[11] = { 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 100. };
      Double_t z[41] = { 0. };
      for (Int_t k = 0; k <= 40; k++) {
        z[k] = 0.5 + k*1.;
      }
      TH3D* hFlowHist = new TH3D(Form("hQ%dCumuHist%s", n, type.Data()), 
                   Form("hQ%dCumuHist%s", n, type.Data()), fEtaBins, x, 10, y, 40, z);
//      hFlowHist->RebinAxis(40/fEtaBins, hFlowHist->GetXaxis());
      hFlowHist->Sumw2();
      fOutputList->Add(hFlowHist);
      TString tag = TString();
      for (Int_t c = -2; c < 10; c++) {
        // Output histograms
        if (c == -2)      tag = "mb";
        else if (c == -1) tag = "0_40";
        else              tag = Form("%d_%d", (Int_t)y[c], (Int_t)y[c+1]);
    
        TH1D* hCumulant2RefFlow = new TH1D(Form("hQ%dCumulant2RefFlow%s_%s", n, type.Data(), tag.Data()), 
                     Form("hQ%dCumulant2RefFlow%s_%s", n, type.Data(), tag.Data()), fEtaBins, x);
        hCumulant2RefFlow->Sumw2();
        fOutputList->Add(hCumulant2RefFlow);
   
        TH1D* hCumulant4RefFlow = new TH1D(Form("hQ%dCumulant4RefFlow%s_%s", n, type.Data(), tag.Data()), 
                     Form("hQ%dCumulant4RefFlow%s_%s", n, type.Data(), tag.Data()), fEtaBins, x);
        hCumulant4RefFlow->Sumw2();
        fOutputList->Add(hCumulant4RefFlow);
 
        TH1D* hCumulant2DiffFlow = new TH1D(Form("hQ%dCumulant2DiffFlow%s_%s", n, type.Data(), tag.Data()), 
                     Form("hQ%dCumulant2DiffFlow%s_%s", n, type.Data(), tag.Data()), fEtaBins, x);
        hCumulant2DiffFlow->Sumw2();
        fOutputList->Add(hCumulant2DiffFlow);
   
        TH1D* hCumulant4DiffFlow = new TH1D(Form("hQ%dCumulant4DiffFlow%s_%s", n, type.Data(), tag.Data()), 
                     Form("hQ%dCumulant4DiffFlow%s_%s", n, type.Data(), tag.Data()), fEtaBins, x);
        hCumulant4DiffFlow->Sumw2();
        fOutputList->Add(hCumulant4DiffFlow);
      } // end of centrality loop

    } // end of v_{n} loop
 
    // Single Event histograms
    TH1D* hdNdphiSE = new TH1D(Form("hdNdphiSE%s", type.Data()),
                 Form("hdNdphiSE%s", type.Data()), 20, 0, 2*TMath::Pi());
    hdNdphiSE->Sumw2();
    fOutputList->Add(hdNdphiSE);

    TH2D* hdNdetaBindphiSE = new TH2D(Form("hdNdetaBindphiSE%s", type.Data()), 
                 Form("hdNdetaBindphiSE%s", type.Data()), fEtaBins, -4, 6, 20, 0, 2*TMath::Pi());
    hdNdetaBindphiSE->Sumw2();
    fOutputList->Add(hdNdetaBindphiSE);
  } // end of type loop

  TH1D* cent = new TH1D("cent", "cent", 100, 0, 100);
  fOutputList->Add(cent);

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


  // fill histograms
  if (!fFlowUtil->LoopAODFMD(fAOD)) return;

  // Run analysis on FMD
  for (Int_t n = 1; n <= 4; n++) {
    if (fv[n])
      CumulantsMethod("", n);
  }
  
  // Find out if there's any MC data present
  if (!fMC) {
    TClonesArray* mcArray = 0;
    mcArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
    if (mcArray) fMC = kTRUE;
  }
  if (fMC) 
    ProcessPrimary();

  if (!fFlowUtil->LoopAODSPD(fAOD)) return;
  for (Int_t n = 1; n <= 4; n++) {
    if (fv[n])
      CumulantsMethod("SPD", n);
  }
   
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
  
  // We get histograms depending on if it's real data or MC truth data
  TH3D* flowHist   = (TH3D*)fOutputList->FindObject(Form("hQ%dCumuHist%s", harmonic, type.Data()));
  TH2D* dNdetaBindphi = (TH2D*)fOutputList->FindObject(Form("hdNdetaBindphiSE%s", type.Data()));
  TH1D* dNdphi = 0;
  if (!type.Contains("SPD")) dNdphi = (TH1D*)fOutputList->FindObject(Form("hdNdphiSE%s", type.Data()));
  if ( type.Contains("SPD")) dNdphi = (TH1D*)fOutputList->FindObject("hdNdphiSE");

  // We create the objects needed for the analysis
  Double_t mult = dNdphi->GetBinContent(0);
  if (type.Length() <= 1) fCent = dNdphi->GetBinContent(dNdphi->GetNbinsX()+1);

  TH1D* cent = (TH1D*)fOutputList->FindObject("cent");
  if (type.Length() <= 1) cent->Fill(fCent);

  Double_t dQnRe = 0, dQ2nRe = 0, dQnIm = 0, dQ2nIm = 0;
  Double_t pnRe = 0, p2nRe = 0, qnRe = 0, q2nRe = 0, pnIm = 0, p2nIm = 0, qnIm = 0, q2nIm = 0;
  Double_t avg2 = 0, avg4 = 0, avg2p = 0, avg4p = 0;
  Double_t w2avg2sq = 0, w2pavg2psq = 0, w4avg4sq = 0, w4pavg4psq = 0;
  Double_t w2w2pavg2avg2p = 0, w2w4avg2avg4 = 0, w2pw4pavg2pavg4p = 0;
  Double_t w2w4pavg2avg4p = 0, w4w2pavg4avg2p = 0, w4w4pavg4avg4p = 0;
  Double_t cosPhi1Phi2 = 0, cosPhi1Phi2Phi3m = 0, cosPhi1Phi2Phi3p = 0;
  Double_t sinPhi1Phi2 = 0, sinPhi1Phi2Phi3m = 0, sinPhi1Phi2Phi3p = 0;
  Double_t phi = 0;
  Double_t multi = 0, multp = 0, mp = 0, mq = 0;
  Double_t w2 = 0, w4 = 0, w2p = 0, w4p = 0;

  // We loop over the data 1 time!
  for (Int_t etaBin = 1; etaBin <= dNdetaBindphi->GetNbinsX(); etaBin++) {
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
      multp = dNdetaBindphi->GetBinContent(etaBin, phiN);
      mp += multp;
      pnRe += multp*TMath::Cos(n*phi);
      pnIm += multp*TMath::Sin(n*phi);
      p2nRe += multp*TMath::Cos(2.*n*phi);
      p2nIm += multp*TMath::Sin(2.*n*phi);    
    }

    // The reference flow is calculated for the differential flow histograms
    if (etaBin == 1) {
      Double_t eta = flowHist->GetXaxis()->GetBinCenter(0);

      // 2-particle
      w2 = mult * (mult - 1.);
      avg2 = dQnRe*dQnRe + dQnIm*dQnIm - mult;
      avg2 /= w2;
      w2avg2sq = w2 * avg2 * avg2; 

      flowHist->Fill(eta, fCent, kRW2avg2, w2 * avg2);
      flowHist->Fill(eta, fCent, kRW2, w2);
      flowHist->Fill(eta, fCent, kRW2avg2sq, w2avg2sq);

      flowHist->Fill(eta, fCent, kQnRe, dQnRe);
      flowHist->Fill(eta, fCent, kQnIm, dQnIm);
      flowHist->Fill(eta, fCent, kRM, mult);

      // 4-particle
      w4 = mult * (mult - 1.) * (mult - 2.) * (mult - 3.);
      Double_t real = dQ2nRe*dQnRe*dQnRe - dQ2nRe*dQnIm*dQnIm + 2.*dQ2nIm*dQnRe*dQnIm;

      avg4 = TMath::Power(dQnRe*dQnRe + dQnIm*dQnIm, 2); 
      avg4 += dQ2nRe*dQ2nRe + dQ2nIm*dQ2nIm - 2.*real;
      avg4 -= 4.*(mult - 2.)*(dQnRe*dQnRe + dQnIm*dQnIm) - 2.*mult*(mult - 3.);
  
      avg4 /= w4;
      w4avg4sq = w4 * avg4 * avg4;
      w2w4avg2avg4 = w2 * w4 * avg2 * avg4;

      flowHist->Fill(eta, fCent, kRW4avg4, w4 * avg4);
      flowHist->Fill(eta, fCent, kRW4, w4);
      flowHist->Fill(eta, fCent, kRW4avg4sq, w4avg4sq);
      flowHist->Fill(eta, fCent, kRW2w4, w2 * w4);
      flowHist->Fill(eta, fCent, kRW2w4avg2avg4, w2w4avg2avg4);

      cosPhi1Phi2 = dQnRe*dQnRe - dQnIm*dQnIm - dQ2nRe;
      sinPhi1Phi2 = 2.*dQnRe*dQnIm - dQ2nIm;
      
      cosPhi1Phi2Phi3m = TMath::Power(dQnRe, 3) + dQnRe*dQnIm*dQnIm; 
      cosPhi1Phi2Phi3m -= dQnRe*dQ2nRe + dQnIm*dQ2nIm + 2.*(mult - 1.)*dQnRe;

      sinPhi1Phi2Phi3m = -TMath::Power(dQnIm, 3) - dQnRe*dQnRe*dQnIm; 
      sinPhi1Phi2Phi3m -= dQnIm*dQ2nRe - dQnRe*dQ2nIm - 2.*(mult - 1.)*dQnIm;

      flowHist->Fill(eta, fCent, kRCosphi1phi2, cosPhi1Phi2);
      flowHist->Fill(eta, fCent, kRSinphi1phi2, sinPhi1Phi2);
      flowHist->Fill(eta, fCent, kRCosphi1phi2phi3m, cosPhi1Phi2Phi3m);
      flowHist->Fill(eta, fCent, kRSinphi1phi2phi3m, sinPhi1Phi2Phi3m);
      flowHist->Fill(eta, fCent, kRMm1m2, mult*(mult-1.)*(mult-2.));

      // Count number of events
      flowHist->Fill(eta, fCent, 0., 1.);
    } 

    // Differential flow calculations for each etaBin bin is done:
    if (mp == 0) continue;
    Double_t eta = dNdetaBindphi->GetXaxis()->GetBinCenter(etaBin);
//    eta = TMath::Abs(eta);

    mq = mp;
    qnRe = pnRe;
    qnIm = pnIm;
    q2nRe = p2nRe;
    q2nIm = p2nIm;
    if (type.Contains("SPD") || (type.Contains("") && eta >= 4.)) {
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
    avg2p = pnRe*pnRe + pnIm*pnIm - mp;
    avg2p /= w2p;
    w2pavg2psq = w2p * avg2p * avg2p; 

    flowHist->Fill(eta, fCent, kRW2avg2, w2p * avg2p);
    flowHist->Fill(eta, fCent, kRW2, w2p);
    flowHist->Fill(eta, fCent, kRW2avg2sq, w2pavg2psq);

    flowHist->Fill(eta, fCent, kRM, mp);

    // 4-particle
    w4p = mp * (mp - 1.) * (mp - 2.) * (mp - 3.);
    Double_t real = p2nRe*pnRe*pnRe - p2nRe*pnIm*pnIm + 2.*p2nIm*pnRe*pnIm;

    avg4p = TMath::Power(pnRe*pnRe + pnIm*pnIm, 2); 
    avg4p += p2nRe*p2nRe + p2nIm*p2nIm - 2.*real;
    avg4p -= 4.*(mp - 2.)*(pnRe*pnRe + pnIm*pnIm) - 2.*mp*(mp - 3.);
  
    avg4p /= w4p;
    w4pavg4psq = w4p * avg4p * avg4p;
    w2w4pavg2avg4p = w2p * w4p * avg2p * avg4p;

    flowHist->Fill(eta, fCent, kRW4avg4, w4p * avg4p);
    flowHist->Fill(eta, fCent, kRW4, w4p);
    flowHist->Fill(eta, fCent, kRW4avg4sq, w4pavg4psq);
    flowHist->Fill(eta, fCent, kRW2w4, w2p * w4p);
    flowHist->Fill(eta, fCent, kRW2w4avg2avg4, w2w4pavg2avg4p);

    cosPhi1Phi2 = pnRe*pnRe - pnIm*pnIm - p2nRe;
    sinPhi1Phi2 = 2.*pnRe*pnIm - p2nIm;
      
    cosPhi1Phi2Phi3m = TMath::Power(pnRe, 3) + pnRe*pnIm*pnIm; 
    cosPhi1Phi2Phi3m -= pnRe*p2nRe + pnIm*p2nIm + 2.*(mp - 1.)*pnRe;

    sinPhi1Phi2Phi3m = -TMath::Power(pnIm, 3) - pnRe*pnRe*pnIm; 
    sinPhi1Phi2Phi3m -= pnIm*p2nRe - pnRe*p2nIm - 2.*(mp - 1.)*pnIm;

    flowHist->Fill(eta, fCent, kRCosphi1phi2, cosPhi1Phi2);
    flowHist->Fill(eta, fCent, kRSinphi1phi2, sinPhi1Phi2);
    flowHist->Fill(eta, fCent, kRCosphi1phi2phi3m, cosPhi1Phi2Phi3m);
    flowHist->Fill(eta, fCent, kRSinphi1phi2phi3m, sinPhi1Phi2Phi3m);
    flowHist->Fill(eta, fCent, kRMm1m2, mp*(mp-1.)*(mp-2.));

    // 2-particle differential flow
    w2p = mp * mult - mq;
    avg2p = pnRe*dQnRe + pnIm*dQnIm - mq;
    avg2p /= w2p;
    w2pavg2psq = w2p * avg2p * avg2p;
    w2w2pavg2avg2p = w2 * w2p * avg2 * avg2p;
    
    flowHist->Fill(eta, fCent, kW2avg2, w2p * avg2p);
    flowHist->Fill(eta, fCent, kW2, w2p);
    flowHist->Fill(eta, fCent, kW2avg2sq, w2pavg2psq);
    flowHist->Fill(eta, fCent, kW2w2pavg2avg2p, w2w2pavg2avg2p);
    flowHist->Fill(eta, fCent, kW2w2p, w2 * w2p);

    flowHist->Fill(eta, fCent, kQnRe, pnRe);
    flowHist->Fill(eta, fCent, kQnIm, pnIm);
    flowHist->Fill(eta, fCent, kM, mp);

    // 4-particle differential flow
    w4p = (mp * mult - 3.*mq)*(mult - 1.)*(mult - 2.);

    avg4p =  pnRe*dQnRe*(dQnRe*dQnRe + dQnIm*dQnIm) + pnIm*dQnIm*(dQnRe*dQnRe + dQnIm*dQnIm);
    avg4p -= q2nRe*dQnRe*dQnRe - q2nRe*dQnIm*dQnIm + 2.*q2nIm*dQnRe*dQnIm;
    avg4p -= pnRe*dQnRe*dQ2nRe - pnRe*dQnIm*dQ2nIm + pnIm*dQnRe*dQ2nIm + pnIm*dQnIm*dQ2nRe;
    avg4p -= 2.*mult*(pnRe*dQnRe + pnIm*dQnIm);

    avg4p += - 2.*mq*(dQnRe*dQnRe + dQnIm*dQnIm) + 7.*(qnRe*dQnRe + qnIm*dQnIm); 
    avg4p += - (dQnRe*qnRe + dQnIm*qnIm) + (q2nRe*dQ2nRe + q2nIm*dQ2nIm);
    avg4p += 2.*(pnRe*dQnRe + pnIm*dQnIm) + 2.*mq*mult - 6.*mq;
    avg4p /= w4p;

    w4pavg4psq = w4p * avg4p * avg4p;
    w2w4pavg2avg4p = w2 * w4p * avg2 * avg4p;
    w4w2pavg4avg2p = w4 * w2p * avg4 * avg2p;
    w4w4pavg4avg4p = w4 * w4p * avg4 * avg4p;
    w2pw4pavg2pavg4p = w2p * w4p * avg2p * avg4p;

    flowHist->Fill(eta, fCent, kW4avg4, w4p * avg4p);
    flowHist->Fill(eta, fCent, kW4, w4p);
    flowHist->Fill(eta, fCent, kW4avg4sq, w4pavg4psq);
    flowHist->Fill(eta, fCent, kW2w4, w2p * w4p);
    flowHist->Fill(eta, fCent, kW2w4avg2avg4, w2pw4pavg2pavg4p);
    flowHist->Fill(eta, fCent, kW2w4p, w2 * w4p);
    flowHist->Fill(eta, fCent, kW2w4pavg2avg4p, w2w4pavg2avg4p);
    flowHist->Fill(eta, fCent, kW4w2p, w4 * w2p);
    flowHist->Fill(eta, fCent, kW4w2pavg4avg2p, w4w2pavg4avg2p);
    flowHist->Fill(eta, fCent, kW4w4p, w4 * w4p);
    flowHist->Fill(eta, fCent, kW4w4pavg4avg4p, w4w4pavg4avg4p);

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

    flowHist->Fill(eta, fCent, kCosphi1phi2, cosPhi1Phi2);
    flowHist->Fill(eta, fCent, kSinphi1phi2, sinPhi1Phi2);
    flowHist->Fill(eta, fCent, kCosphi1phi2phi3m, cosPhi1Phi2Phi3m);
    flowHist->Fill(eta, fCent, kSinphi1phi2phi3m, sinPhi1Phi2Phi3m);
    flowHist->Fill(eta, fCent, kMm1m2, (mp*mult-2.*mq)*(mult-1.));
    flowHist->Fill(eta, fCent, kCosphi1phi2phi3p, cosPhi1Phi2Phi3p);
    flowHist->Fill(eta, fCent, kSinphi1phi2phi3p, sinPhi1Phi2Phi3p); 

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

  TH3D* hcumulantsHist = 0;
  TH2D* cumulantsHist = new TH2D("tmphist", "tmphist",  fEtaBins, -4, 6, 40, 0.5, 40.5);
  TH1D* cumulant2refHist = 0;
  TH1D* cumulant4refHist = 0;
  TH1D* cumulant2diffHist = 0;
  TH1D* cumulant4diffHist = 0;

  // For flow calculations
  Double_t two = 0, qc2 = 0, vTwo2 = 0, four = 0, qc4 = 0, vTwo4 = 0; 
  Double_t twoPrime = 0, qc2Prime = 0, vTwo2diff = 0, fourPrime = 0, qc4Prime = 0, vTwo4diff = 0;
  Double_t w2 = 0, w4 = 0, w2p = 0, w4p = 0, sqrtW2sq = 0, sqrtW2psq = 0, w2W2p = 0;
  Double_t w2W4 = 0, w2W4p = 0, w4W2p = 0, w4W4p = 0, w2pW4p = 0;
  Double_t sqrtW4sq = 0, sqrtW4psq = 0;
  Double_t w2avg2 = 0, w2pavg2p = 0, w4avg4 = 0, w4pavg4p = 0;
  Double_t cosP1nPhi = 0, sinP1nPhi = 0, mult = 0, cosP1nPhi1P1nPhi2 = 0, sinP1nPhi1P1nPhi2 = 0;
  Double_t cosP1nPhi1M1nPhi2M1nPhi3 = 0, sinP1nPhi1M1nPhi2M1nPhi3 = 0, multm1m2 = 0;
  Double_t cosP1nPsi = 0, sinP1nPsi = 0, mp = 0, cosP1nPsi1P1nPhi2 = 0, sinP1nPsi1P1nPhi2 = 0;
  Double_t cosP1nPsi1M1nPhi2M1nPhi3 = 0, sinP1nPsi1M1nPhi2M1nPhi3 = 0, mpqMult = 0;
  Double_t cosP1nPsi1P1nPhi2M1nPhi3 = 0, sinP1nPsi1P1nPhi2M1nPhi3 = 0;

  // For error calculations
  Double_t w2avg2sq = 0, w2W2pavg2avg2p = 0, w2pavg2psq = 0;
  Double_t w4avg4sq = 0, w2W4avg2avg4 = 0, w4W4pavg4avg4p = 0, w4pavg4psq = 0;
  Double_t w2W4pavg2avg4p = 0, w4W2pavg4avg2p = 0;
  Double_t stwosq = 0, stwoPrimesq = 0, sfoursq = 0, sfourPrimesq = 0;
  Double_t vTwo2err = 0, vTwo2diffErr = 0, vTwo4err = 0, vTwo4diffErr = 0;
  Double_t cov22p = 0, cov24 = 0, cov24p = 0, cov42p = 0, cov44p = 0, cov2p2np = 0;
  Double_t w2pW4pavg2pavg4p = 0;
  

  Int_t nLoops = (fMC ? 4 : 2);
  TString type = "";

  Int_t y[11] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100 };
  
  // Do a loop over the difference analysis types, calculating flow
  // 2 loops for real data, 3 for MC data
  // inside each is a nested loop over each harmonic (1, 2, 3 and 4 at the moment)
  for (Int_t loop = 1; loop <= nLoops; loop++) {

    if (loop == 1) type = "";
    if (loop == 2) type = "SPD";
    if (loop == 3) type = "MC";
    if (loop == 4) type = "TrRef";
    
    for (Int_t n = 1; n <= 4; n++) {
      if (!fv[n]) continue;

      hcumulantsHist = (TH3D*)fOutputList->FindObject(Form("hQ%dCumuHist%s", n, type.Data()));
 
      // Centrality loop
      for (Int_t c = -2; c < hcumulantsHist->GetNbinsY(); c++) {
        if (c == -2) {
          hcumulantsHist->GetYaxis()->SetRange(0, 9);
          cumulantsHist = (TH2D*)hcumulantsHist->Project3D("zx oe");
          
          cumulant2refHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant2RefFlow%s_mb", n, type.Data()));
          cumulant4refHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant4RefFlow%s_mb", n, type.Data()));
          cumulant2diffHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant2DiffFlow%s_mb", n, type.Data()));
          cumulant4diffHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant4DiffFlow%s_mb", n, type.Data()));
        }
        else if (c == -1) {
          hcumulantsHist->GetYaxis()->SetRange(1, 5);
          cumulantsHist = (TH2D*)hcumulantsHist->Project3D("zx oe");
          
          cumulant2refHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant2RefFlow%s_0_40", n, type.Data()));
          cumulant4refHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant4RefFlow%s_0_40", n, type.Data()));
          cumulant2diffHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant2DiffFlow%s_0_40", n, type.Data()));
          cumulant4diffHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant4DiffFlow%s_0_40", n, type.Data()));
        }
        else {
          hcumulantsHist->GetYaxis()->SetRange(c+1, c+1);
          cumulantsHist = (TH2D*)hcumulantsHist->Project3D("zx e");
          cumulant2refHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant2RefFlow%s_%d_%d", n, type.Data(), y[c], y[c+1]));
          cumulant4refHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant4RefFlow%s_%d_%d", n, type.Data(), y[c], y[c+1]));
          cumulant2diffHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant2DiffFlow%s_%d_%d", n, type.Data(), y[c], y[c+1]));
          cumulant4diffHist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant4DiffFlow%s_%d_%d", n, type.Data(), y[c], y[c+1]));
        }

        Bool_t refDone = kFALSE;

        for (Int_t etaBin = 1; etaBin <= cumulantsHist->GetNbinsX(); etaBin++) {
          if (refDone == kFALSE || etaBin == 0) {
            // 2-particle reference flow
            w2avg2 = cumulantsHist->GetBinContent(etaBin, kRW2avg2);
            if (etaBin == cumulantsHist->GetNbinsX()) {
              refDone = kTRUE;
              etaBin = -1;
              continue;
            }
            if (!w2avg2) continue;
            w2 = cumulantsHist->GetBinContent(etaBin, kRW2);
            cosP1nPhi = cumulantsHist->GetBinContent(etaBin, kQnRe);
            sinP1nPhi = cumulantsHist->GetBinContent(etaBin, kQnIm);
            mult = cumulantsHist->GetBinContent(etaBin, kRM);
            cosP1nPhi /= mult;
            sinP1nPhi /= mult;
            two = w2avg2 / w2;
            qc2 = two  - TMath::Power(cosP1nPhi, 2) - TMath::Power(sinP1nPhi, 2); 
            vTwo2 = TMath::Sqrt(qc2);
            if (etaBin == 0) cumulant2diffHist->SetBinContent(etaBin, vTwo2);
            if (etaBin > 0)  cumulant2refHist->SetBinContent(etaBin, vTwo2);

            // 4-particle reference flow
            w4avg4 = cumulantsHist->GetBinContent(etaBin, kRW4avg4);
            w4 = cumulantsHist->GetBinContent(etaBin, kRW4);
            cosP1nPhi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, kRCosphi1phi2);
            sinP1nPhi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, kRSinphi1phi2);
            cosP1nPhi1P1nPhi2 /= w2;
            sinP1nPhi1P1nPhi2 /= w2;
            cosP1nPhi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, kRCosphi1phi2phi3m);
            sinP1nPhi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, kRSinphi1phi2phi3m);
            multm1m2 = cumulantsHist->GetBinContent(etaBin, kRMm1m2);
            cosP1nPhi1M1nPhi2M1nPhi3 /= multm1m2;
            sinP1nPhi1M1nPhi2M1nPhi3 /= multm1m2;
            four = w4avg4 / w4;
            qc4 = four-2.*pow(two,2.)
               - 4.*cosP1nPhi*cosP1nPhi1M1nPhi2M1nPhi3
               + 4.*sinP1nPhi*sinP1nPhi1M1nPhi2M1nPhi3-pow(cosP1nPhi1P1nPhi2,2.)-pow(sinP1nPhi1P1nPhi2,2.)
               + 4.*cosP1nPhi1P1nPhi2*(pow(cosP1nPhi,2.)-pow(sinP1nPhi,2.))+8.*sinP1nPhi1P1nPhi2*sinP1nPhi*cosP1nPhi
               + 8.*two*(pow(cosP1nPhi,2.)+pow(sinP1nPhi,2.))-6.*pow((pow(cosP1nPhi,2.)+pow(sinP1nPhi,2.)),2.);


            vTwo4 = TMath::Power(-qc4, 0.25);
            if (etaBin == 0) cumulant4diffHist->SetBinContent(etaBin, vTwo4);
            if (etaBin > 0) cumulant4refHist->SetBinContent(etaBin, vTwo4);
   
            // 2-particle reference flow error calculations
            w2avg2sq = cumulantsHist->GetBinContent(etaBin, kRW2avg2sq);
            sqrtW2sq = cumulantsHist->GetBinError(etaBin, kRW2);
   
            stwosq = VarSQ(w2avg2sq, two, w2, w2avg2, sqrtW2sq);
            vTwo2err = sqrtW2sq * TMath::Sqrt(stwosq) / (2. * TMath::Sqrt(two) * w2);
            if (etaBin == 0) cumulant2diffHist->SetBinError(etaBin, vTwo2err);
            if (etaBin > 0)  cumulant2refHist->SetBinError(etaBin, vTwo2err);
    
            // 4-particle reference flow error calculations
            w4avg4sq = cumulantsHist->GetBinContent(etaBin, kRW4avg4sq);
            sqrtW4sq = cumulantsHist->GetBinError(etaBin, kRW4);
            w2W4 = cumulantsHist->GetBinContent(etaBin, kRW2w4);
            w2W4avg2avg4 = cumulantsHist->GetBinContent(etaBin, kRW2w4avg2avg4);
   
            sfoursq = VarSQ(w4avg4sq, four, w4, w4avg4, sqrtW4sq);
            cov24 = CovXY(w2W4avg2avg4, w2W4, two*four, w2, w4);
  
            vTwo4err = two*two * TMath::Power(sqrtW2sq, 2) * stwosq / (w2*w2);
            vTwo4err += TMath::Power(sqrtW4sq, 2) * sfoursq / (16. * w4*w4);
            vTwo4err -= two * w2W4 * cov24 / (2. * w2 * w4);
            vTwo4err /= TMath::Power(2. * two*two - four, 1.5);
            vTwo4err = TMath::Sqrt(vTwo4err);
            if (etaBin == 0) cumulant4diffHist->SetBinError(etaBin, vTwo4err);
            if (etaBin > 0)  cumulant4refHist->SetBinError(etaBin, vTwo4err);

            continue;
          }

          // 2-particle differential flow
          w2pavg2p = cumulantsHist->GetBinContent(etaBin, kW2avg2);
          if (!w2pavg2p) continue;
          w2p = cumulantsHist->GetBinContent(etaBin, kW2);
          cosP1nPsi = cumulantsHist->GetBinContent(etaBin, kQnRe);
          sinP1nPsi = cumulantsHist->GetBinContent(etaBin, kQnIm);
          mp = cumulantsHist->GetBinContent(etaBin, kM);
          cosP1nPsi /= mp;
          sinP1nPsi /= mp;
          twoPrime = w2pavg2p / w2p;
          qc2Prime = twoPrime - sinP1nPsi*sinP1nPhi - cosP1nPsi*cosP1nPhi;
          vTwo2diff = qc2Prime / TMath::Sqrt(qc2);
          cumulant2diffHist->SetBinContent(etaBin, vTwo2diff);
 
          // 4-particle differential flow
          w4pavg4p = cumulantsHist->GetBinContent(etaBin, kW4avg4);
          w4p = cumulantsHist->GetBinContent(etaBin, kW4);
          cosP1nPsi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, kCosphi1phi2);
          sinP1nPsi1P1nPhi2 = cumulantsHist->GetBinContent(etaBin, kSinphi1phi2);
          cosP1nPsi1P1nPhi2 /= w2p;
          sinP1nPsi1P1nPhi2 /= w2p;
          cosP1nPsi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, kCosphi1phi2phi3m);
          sinP1nPsi1M1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, kSinphi1phi2phi3m);
          mpqMult = cumulantsHist->GetBinContent(etaBin, kMm1m2);
          cosP1nPsi1M1nPhi2M1nPhi3 /= mpqMult;
          sinP1nPsi1M1nPhi2M1nPhi3 /= mpqMult;
          cosP1nPsi1P1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, kCosphi1phi2phi3p);
          sinP1nPsi1P1nPhi2M1nPhi3 = cumulantsHist->GetBinContent(etaBin, kSinphi1phi2phi3p);
          cosP1nPsi1P1nPhi2M1nPhi3 /= mpqMult;
          sinP1nPsi1P1nPhi2M1nPhi3 /= mpqMult;
  
          fourPrime = w4pavg4p / w4p;
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
                    + 2.*cosP1nPsi1P1nPhi2*(pow(cosP1nPhi,2.)-pow(sinP1nPhi,2.))
                    + 4.*sinP1nPsi1P1nPhi2*cosP1nPhi*sinP1nPhi
                    + 4.*twoPrime*(pow(cosP1nPhi,2.)+pow(sinP1nPhi,2.))
                    - 6.*(pow(cosP1nPhi,2.)-pow(sinP1nPhi,2.)) 
                    * (cosP1nPsi*cosP1nPhi-sinP1nPsi*sinP1nPhi)
                    - 12.*cosP1nPhi*sinP1nPhi
                    * (sinP1nPsi*cosP1nPhi+cosP1nPsi*sinP1nPhi);
 
          vTwo4diff = - qc4Prime / TMath::Power(-qc4, 0.75);
          cumulant4diffHist->SetBinContent(etaBin, vTwo4diff);    
       
          // 2-particle differential flow error calculations
          w2pavg2psq = cumulantsHist->GetBinContent(etaBin, kW2avg2sq);
          sqrtW2psq = cumulantsHist->GetBinError(etaBin, kW2);
          w2W2pavg2avg2p = cumulantsHist->GetBinContent(etaBin, kW2w2pavg2avg2p);
          w2W2p = cumulantsHist->GetBinContent(etaBin, kW2w2p);
    
          cov22p = CovXY(w2W2pavg2avg2p, w2W2p, two*twoPrime, w2, w2p);
          stwoPrimesq = VarSQ(w2pavg2psq, twoPrime, w2p, w2pavg2p, sqrtW2psq);
   
          vTwo2diffErr = twoPrime*twoPrime*TMath::Power(sqrtW2sq, 2)*stwosq/(w2*w2);
          vTwo2diffErr += 4.*two*two*TMath::Power(sqrtW2psq, 2)*stwoPrimesq/(w2p*w2p);
          vTwo2diffErr -= 4.*two*twoPrime*w2W2p*cov22p/(w2*w2p);
          vTwo2diffErr /= (4. * TMath::Power(two, 3));
          vTwo2diffErr = TMath::Sqrt(vTwo2diffErr);
          cumulant2diffHist->SetBinError(etaBin, vTwo2diffErr);
  
          // 4-particle differential flow error calculations
          sqrtW4psq = cumulantsHist->GetBinError(etaBin, kW4);
          w4pavg4psq = cumulantsHist->GetBinContent(etaBin, kW4avg4sq);
          w2pW4p = cumulantsHist->GetBinContent(etaBin, kW2w4);
          w2pW4pavg2pavg4p = cumulantsHist->GetBinContent(etaBin, kW2w4avg2avg4);
          w2W4p = cumulantsHist->GetBinContent(etaBin, kW2w4p);
          w2W4pavg2avg4p = cumulantsHist->GetBinContent(etaBin, kW2w4pavg2avg4p);
          w4W2p = cumulantsHist->GetBinContent(etaBin, kW4w2p);
          w4W2pavg4avg2p = cumulantsHist->GetBinContent(etaBin, kW4w2pavg4avg2p);
          w4W4p = cumulantsHist->GetBinContent(etaBin, kW4w4p);
          w4W4pavg4avg4p = cumulantsHist->GetBinContent(etaBin, kW4w4pavg4avg4p);
        
          sfourPrimesq = VarSQ(w4pavg4psq, fourPrime, w4p, w4pavg4p, sqrtW4psq);
          cov24p = CovXY(w2W4pavg2avg4p, w2W4p, two*fourPrime, w2, w4p);
          cov42p = CovXY(w4W2pavg4avg2p, w4W2p, four*twoPrime, w4, w2p);
          cov44p = CovXY(w4W4pavg4avg4p, w4W4p, four*fourPrime, w4, w4p);
          cov2p2np = CovXY(w2pW4pavg2pavg4p, w2pW4p, twoPrime*fourPrime, w2p, w4p);
   
          // Numbers on the side reference term number in paper (cite needed) loosely
    /*1*/ vTwo4diffErr =  TMath::Power(2.*two*two*twoPrime - 3.*two*fourPrime + 2.*four*twoPrime, 2) 
                          * TMath::Power(sqrtW2sq, 2) * stwosq / (w2*w2);
    /*2*/ vTwo4diffErr += 9. * TMath::Power(2.*two*twoPrime - fourPrime, 2) * TMath::Power(sqrtW4sq, 2)
                          * sfoursq / (16. * w4*w4);
    /*3*/ vTwo4diffErr += 4. * two*two * TMath::Power(2.*two*two - four, 2) * TMath::Power(sqrtW2psq, 2)
                          * stwoPrimesq / (w2p*w2p);
    /*4*/ vTwo4diffErr += TMath::Power(2.*two*two - four, 2) * TMath::Power(sqrtW4psq, 2) * sfourPrimesq
                          / (w4p*w4p);
    /*5*/ vTwo4diffErr -= 1.5 * (2.*two*twoPrime - fourPrime) * (2.*two*two*twoPrime - 3.*two*fourPrime + 2.*four*twoPrime)
                          * w2W4 * cov24 / (w2*w4);
    /*6*/ vTwo4diffErr -= 4. * two * (2.*two*two - four) 
                          * (2.*two*two*twoPrime - 3.*two*fourPrime + 2.*four*twoPrime)
                          * w2W2p * cov22p / (w2 * w2p);
    /*7*/ vTwo4diffErr += 2. * (2.*two*two - four)
                          * (2.*two*two*twoPrime - 3.*two*fourPrime + 2.*four*twoPrime)
                          * w2W4p * cov24p / (w2 * w4p);
    /*8*/ vTwo4diffErr += 3.*two*(2.*two*two - four)*(2.*two*twoPrime - fourPrime)
                          * w4W2p * cov42p / (w4*w2p);
    /*9*/ vTwo4diffErr -= 1.5 * (2.*two*two - four)*(2.*two*twoPrime - fourPrime)
                          * w4W4p * cov44p / (w4 * w4p);
    /*10*/vTwo4diffErr -= 4.*two*TMath::Power(2.*two*two - four, 2)
                          * w2pW4p * cov2p2np / (w2p * w4p);
    /*11*/vTwo4diffErr /= TMath::Power(2.*two*two - four, 3.5);
          vTwo4diffErr = TMath::Sqrt(vTwo4diffErr);
  
          cumulant4diffHist->SetBinError(etaBin, vTwo4diffErr);
        } // End of etaBin loop
        // Number of events:
        Int_t nEv = (Int_t)cumulantsHist->GetBinContent(0,0);
        cumulant2diffHist->SetBinContent(cumulant2diffHist->GetNbinsX() + 1, nEv);
        cumulant4diffHist->SetBinContent(cumulant4diffHist->GetNbinsX() + 1, nEv);
      } // End of centrality loop
    } // End of harmonics loop
  }  // End of type loop

  delete cumulantsHist;

}
// _____________________________________________________________________
void AliForwardFlowTaskQC::ProcessPrimary() 
{
  //
  // If fMC == kTRUE this function takes care of organizing the input 
  // Monte Carlo data and histograms so AliForwardFlowTaskQC::QCumulants 
  // can be run on it.
  //

  if (fFlowUtil->LoopAODtrrefHits(fAOD)) {
   // Run analysis on TrackRefs
    for (Int_t n = 1; n <= 4; n++) {
      if (fv[n])
        CumulantsMethod("TrRef", n);
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
  
  sx = wxx2 + x*x*wx - (Double_t)2*x*wxx;
  sx *= (Double_t)1 / wx;
  sx *= (Double_t)1 / (1 - TMath::Power(sqrtwx2, 2) / (wx*wx));

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
  numerator = 1 - (wxwy / (wx * wy));
  
  cov = denominator / numerator;
  return cov;
}
//_____________________________________________________________________
//
//
// EOF
