//
// Calculate flow in the forward regions using the Q cumulants method
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root
//
// TODO!
// - Add centrality stuff
// END OF TODO!
//
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>
#include <TMath.h>
#include "AliLog.h"
#include "AliForwardFlowTaskQC.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliForwardFlowBase.h"
#include "AliAODForwardMult.h"

ClassImp(AliForwardFlowTaskQC)

AliForwardFlowTaskQC::AliForwardFlowTaskQC()
: fDebug(0),  		// Debug flag
  fOutputList(0),	// Output list
  fAOD(0),		// AOD input event
  fMC(kFALSE),		// MC flag
  fEtaBins(20)		// # of eta bins in histograms
{
  // 
  // Default constructor
  //
}
//_____________________________________________________________________
AliForwardFlowTaskQC::AliForwardFlowTaskQC(const char* name) :
  AliAnalysisTaskSE(name),
  fDebug(0),		// Debug flag
  fOutputList(0),	// Output list
  fAOD(0),		// AOD input event
  fMC(kFALSE),		// MC flag
  fEtaBins(20)		// # of Eta bins
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
  fDebug(o.fDebug),		// Debug flag
  fOutputList(o.fOutputList),	// Output list
  fAOD(o.fAOD),			// AOD input event
  fMC(o.fMC),			// MC flag
  fEtaBins(o.fEtaBins)		// # of Eta bins
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

  if (fEtaBins % 20) fEtaBins = 20;

  // Histograms for cumulants analysis

  // We loop over flow histograms here to add different orders of harmonics
  for (Int_t n = 1; n <= 4; n++) {
    if (!fv[n]) continue;
    // Only one flow histogram is needed for each type of data;
    // x-axis is eta-bins with differential flow, integrated is in underflowbin
    // y-axis bin 1:  (w_<2> * <2>).Re()
    // y-axis bin 2:  (w_<2> * <2>).Im()
    // y-axis bin 3:  w_<2> = M(M-1)
    // y-axis bin 4:  (w_<2> * <2> * <2>).Re()
    // y-axis bin 5:  (w_<2> * <2> * <2>).Im()
    // y-axis bin 6:  (w_<2> * w_<2'> * <2> * <2'>).Re()
    // y-axis bin 7:  (w_<2> * w_<2'> * <2> * <2'>).Im()
    // y-axis bin 8:  w_<2> * w_<2'>
    // y-axis bin 9:  (w_<4> * <4>).Re()
    // y-axis bin 10:  (w_<4> * <4>).Im()
    // y-axis bin 11:  w_<4>
    // y-axis bin 12:  (w_<4> * <4> * <4>).Re()
    // y-axis bin 13:  (w_<4> * <4> * <4>).Im()
    // y-axis bin 14:  w_<2> * w_<4>
    // y-axis bin 15:  (w_<2> * w_<4> * <2> * <4>).Re()
    // y-axis bin 16:  (w_<2> * w_<4> * <2> * <4>).Im()
    // y-axis bin 17:  w_<2> * w_<4'>
    // y-axis bin 18:  (w_<2> * w_<4'> * <2> * <4'>).Re()
    // y-axis bin 19:  (w_<2> * w_<4'> * <2> * <4'>).Im()
    // y-axis bin 20:  w_<4> * w_<2'>
    // y-axis bin 21:  (w_<4> * w_<2'> * <4> * <2'>).Re()
    // y-axis bin 22:  (w_<4> * w_<2'> * <4> * <2'>).Im()
    // y-axis bin 23:  w_<4> * w_<4'>
    // y-axis bin 24:  (w_<4>  * w_<4'> * <4> * <4'>).Re()
    // y-axis bin 25:  (w_<4>  * w_<4'> * <4> * <4'>).Im()
    // y-axis bin 26: Qn or pn.Re() = <<cos2phi or psi>>
    // y-axis bin 27: Qn or pn.Im() = <<sin2phi or psi>>
    // y-axis bin 28: M or mp
    // y-axis bin 29: (Qn*Qn-Q2n).Re() = <<cos(2(phi1 or psi1+phi2))>>
    // y-axis bin 30: (Qn*Qn-Q2n).Im() = <<sin(2(phi1 or psi1+phi2))>>
    // y-axis bin 31: <<cos(2(phi1 or psi1-phi2-phi3))>>
    // y-axis bin 32: <<sin(2(phi1 or psi1-phi2-phi3))>>
    // y-axis bin 33: M*(M-1)*(M-2) or similar for diff
    // y-axis bin 34: <<cos(2(psi1+phi2-phi3>>
    // y-axis bin 35: <<sin(2(psi1+phi2-phi3>>
    TH2D* hFlowHist = new TH2D(Form("hQ%dCumuHist", n), Form("hQ%dCumuHist", n), fEtaBins, -4, 6, 35, 0.5, 35.5);
    hFlowHist->Sumw2();
    fOutputList->Add(hFlowHist);

    TH2D* hFlowHistMC = new TH2D(Form("hQ%dCumuHistMC", n), Form("hQ%dCumuHistMC", n), fEtaBins, -4, 6, 35, 0.5, 35.5);
    hFlowHistMC->Sumw2();
    fOutputList->Add(hFlowHistMC);

    TH2D* hFlowHistTrRef = new TH2D(Form("hQ%dCumuHistTrRef", n), Form("hQ%dCumuHistTrRef", n), fEtaBins, -4, 6, 35, 0.5, 35.5);
    hFlowHistTrRef->Sumw2();
    fOutputList->Add(hFlowHistTrRef);

    // Output histograms
    TH1D* hCumulant2Flow = new TH1D(Form("hQ%dCumulant2Flow", n), Form("hQ%dCumulant2Flow", n), fEtaBins, -4, 6);
    hCumulant2Flow->Sumw2();
    fOutputList->Add(hCumulant2Flow);
 
    TH1D* hCumulant2FlowMC = new TH1D(Form("hQ%dCumulant2FlowMC", n),Form("hQ%dCumulant2FlowMC", n), fEtaBins, -4, 6);
    hCumulant2FlowMC->Sumw2();
    fOutputList->Add(hCumulant2FlowMC);
  
    TH1D* hCumulant2FlowTrRef = new TH1D(Form("hQ%dCumulant2FlowTrRef", n), Form("hQ%dCumulant2FlowTrRef", n), fEtaBins, -4, 6);
    hCumulant2FlowTrRef->Sumw2();
    fOutputList->Add(hCumulant2FlowTrRef);


    TH1D* hCumulant4Flow = new TH1D(Form("hQ%dCumulant4Flow", n), Form("hQ%dCumulant4Flow", n), fEtaBins, -4, 6);
    hCumulant4Flow->Sumw2();
    fOutputList->Add(hCumulant4Flow);
  
    TH1D* hCumulant4FlowMC = new TH1D(Form("hQ%dCumulant4FlowMC", n), Form("hQ%dCumulant4FlowMC", n), fEtaBins, -4, 6);
    hCumulant4FlowMC->Sumw2();
    fOutputList->Add(hCumulant4FlowMC);
  
    TH1D* hCumulant4FlowTrRef = new TH1D(Form("hQ%dCumulant4FlowTrRef", n), Form("hQ%dCumulant4FlowTrRef", n), fEtaBins, -4, 6);
    hCumulant4FlowTrRef->Sumw2();
    fOutputList->Add(hCumulant4FlowTrRef);
  }

  // Single Event histograms
  TH1D* hdNdphiSE = new TH1D("hdNdphiSE","hdNdphiSE", 20, 0, 2*TMath::Pi());
  hdNdphiSE->Sumw2();
  fOutputList->Add(hdNdphiSE);

  TH2D* hdNdetadphiSE = new TH2D("hdNdetadphiSE", "hdNdetadphiSE", fEtaBins, -4, 6, 20, 0, 2*TMath::Pi());
  hdNdetadphiSE->Sumw2();
  fOutputList->Add(hdNdetadphiSE);

  TH1D* hdNdphiSEMC = new TH1D("hdNdphiSEMC","hdNdphiSEMC", 20, 0, 2*TMath::Pi());
  hdNdphiSEMC->Sumw2();
  fOutputList->Add(hdNdphiSEMC);

  TH2D* hdNdetadphiSEMC = new TH2D("hdNdetadphiSEMC", "hdNdetadphiSEMC", fEtaBins, -4, 6, 20, 0, 2*TMath::Pi());
  hdNdetadphiSEMC->Sumw2();
  fOutputList->Add(hdNdetadphiSEMC);

  TH1D* hdNdphiSETrRef = new TH1D("hdNdphiSETrRef","hdNdphiSETrRef", 20, 0, 2*TMath::Pi());
  hdNdphiSETrRef->Sumw2();
  fOutputList->Add(hdNdphiSETrRef);

  TH2D* hdNdetadphiSETrRef = new TH2D("hdNdetadphiSETrRef", "hdNdetadphiSETrRef", fEtaBins, -4, 6, 20, 0, 2*TMath::Pi());
  hdNdetadphiSETrRef->Sumw2();
  fOutputList->Add(hdNdetadphiSETrRef);

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

  // Load histograms and reset from last event
  TH1D* dNdphi     = (TH1D*)fOutputList->FindObject("hdNdphiSE");
  TH2D* dNdetadphi  = (TH2D*)fOutputList->FindObject("hdNdetadphiSE");

  dNdphi->Reset();
  dNdetadphi->Reset();

  // Initiate FlowCommon and fill histograms
  AliForwardFlowBase* common = new AliForwardFlowBase(fOutputList);

  if (!common->LoopAODFMD(fAOD)) return;
//  else if (!common->LoopAODSPD(fAOD)) return;
//  if (!common->LoopAODFMDandSPD(fAOD)) return;

  // Run analysis
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

}
//_____________________________________________________________________
void AliForwardFlowTaskQC::CumulantsMethod(TString type = "", Int_t harmonic = 2)
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
  TH2D* flowHist   = (TH2D*)fOutputList->FindObject(Form("hQ%dCumuHist%s", harmonic, type.Data()));
  TH1D* dNdphi     = (TH1D*)fOutputList->FindObject(Form("hdNdphiSE%s", type.Data()));
  TH2D* dNdetadphi = (TH2D*)fOutputList->FindObject(Form("hdNdetadphiSE%s", type.Data()));

  // We create the objects needed for the analysis
  Double_t Mult = dNdphi->GetBinContent(0);

  Double_t QnRe = 0, Q2nRe = 0, QnIm = 0, Q2nIm = 0;
  Double_t pnRe = 0, p2nRe = 0, qnRe = 0, qnnRe = 0, pnIm = 0, p2nIm = 0, qnIm = 0, qnnIm = 0;
  Double_t avg2 = 0, avg4 = 0, avg2p = 0, avg4p = 0;
  Double_t w2avg2sq = 0, w2pavg2psq = 0, w4avg4sq = 0, w4pavg4psq = 0;
  Double_t w2w2pavg2avg2p = 0, w2w4avg2avg4 = 0, w2pw4pavg2pavg4p = 0;
  Double_t w2w4pavg2avg4p = 0, w4w2pavg4avg2p = 0, w4w4pavg4avg4p = 0;
  Double_t CosPhi1Phi2 = 0, CosPhi1Phi2Phi3m = 0, CosPhi1Phi2Phi3p = 0;
  Double_t SinPhi1Phi2 = 0, SinPhi1Phi2Phi3m = 0, SinPhi1Phi2Phi3p = 0;
  Double_t Phii = 0;
  Double_t multi = 0, multp = 0, mp = 0, mq = 0;
  Double_t W2 = 0, W4 = 0, W2p = 0, W4p = 0;

  // We loop over the data 1 time!
  for (Int_t eta = 1; eta <= dNdetadphi->GetNbinsX(); eta++) {
    // The values for each individual eta bins are reset
    mp = 0;
    pnRe = 0;
    p2nRe = 0;
    pnIm = 0;
    p2nIm = 0;

    for (Int_t phii = 1; phii <= dNdphi->GetNbinsX()+1; phii++) {
      Phii = dNdphi->GetXaxis()->GetBinCenter(phii);
      multi = dNdphi->GetBinContent(phii);

      // In the phi loop on the first eta loop the integrated flow
      // is calculated from the dNdphi histogram
      if(eta == 1) {
        QnRe += multi*TMath::Cos(n*Phii);
        QnIm += multi*TMath::Sin(n*Phii);
        Q2nRe += multi*TMath::Cos(2.*n*Phii);
        Q2nIm += multi*TMath::Sin(2.*n*Phii);
      }
      
      // For each eta bin the necesarry values for differential flow
      // is calculated. Here is the loop over the phi's.
      multp = dNdetadphi->GetBinContent(eta, phii);
      mp += multp;
      pnRe += multp*TMath::Cos(n*Phii);
      pnIm += multp*TMath::Sin(n*Phii);
      p2nRe += multp*TMath::Cos(2.*n*Phii);
      p2nIm += multp*TMath::Sin(2.*n*Phii);    
    }

    // The integrated flow is calculated
    if (eta == 1) {
      Double_t Eta = flowHist->GetXaxis()->GetBinCenter(0);

      // 2-particle
      W2 = Mult * (Mult - 1.);
      avg2 = QnRe*QnRe + QnIm*QnIm - Mult;
      avg2 /= W2;
      w2avg2sq = W2 * avg2 * avg2; 

      flowHist->Fill(Eta, 1, W2 * avg2);
      flowHist->Fill(Eta, 3, W2);
      flowHist->Fill(Eta, 4, w2avg2sq);

      flowHist->Fill(Eta, 26, QnRe);
      flowHist->Fill(Eta, 27, QnIm);
      flowHist->Fill(Eta, 28, Mult);

      // 4-particle
      W4 = Mult * (Mult - 1.) * (Mult - 2.) * (Mult - 3.);
      Double_t real = Q2nRe*QnRe*QnRe - Q2nRe*QnIm*QnIm + 2.*Q2nIm*QnRe*QnIm;

      avg4 = TMath::Power(QnRe*QnRe + QnIm*QnIm, 2); 
      avg4 += Q2nRe*Q2nRe + Q2nIm*Q2nIm - 2.*real;
      avg4 -= 4.*(Mult - 2.)*(QnRe*QnRe + QnIm*QnIm) - 2.*Mult*(Mult - 3.);
  
      avg4 /= W4;
      w4avg4sq = W4 * avg4 * avg4;
      w2w4avg2avg4 = W2 * W4 * avg2 * avg4;

      flowHist->Fill(Eta, 9, W4 * avg4);
      flowHist->Fill(Eta, 11, W4);
      flowHist->Fill(Eta, 12, w4avg4sq);
      flowHist->Fill(Eta, 14, W2 * W4);
      flowHist->Fill(Eta, 15, w2w4avg2avg4);

      CosPhi1Phi2 = QnRe*QnRe - QnIm*QnIm - Q2nRe;
      SinPhi1Phi2 = 2.*QnRe*QnIm - Q2nIm;
      
      CosPhi1Phi2Phi3m = TMath::Power(QnRe, 3) + QnRe*QnIm*QnIm; 
      CosPhi1Phi2Phi3m -= QnRe*Q2nRe + QnIm*Q2nIm + 2.*(Mult - 1.)*QnRe;
      SinPhi1Phi2Phi3m = -TMath::Power(QnIm, 3) - QnRe*QnRe*QnIm; 
      SinPhi1Phi2Phi3m -= QnIm*Q2nRe - QnRe*Q2nIm + 2.*(Mult - 1.)*QnIm;

      flowHist->Fill(Eta, 29, CosPhi1Phi2);
      flowHist->Fill(Eta, 30, SinPhi1Phi2);
      flowHist->Fill(Eta, 31, CosPhi1Phi2Phi3m);
      flowHist->Fill(Eta, 32, SinPhi1Phi2Phi3m);
      flowHist->Fill(Eta, 33, Mult*(Mult-1.)*(Mult-2.));

      // Count number of events
      flowHist->Fill(Eta, 0., 1.);
    } // end of harmonics loop

    // Differential flow calculations for each eta bin is done:
    if (mp == 0) continue;
    Double_t Eta = dNdetadphi->GetXaxis()->GetBinCenter(eta);

    mq = mp;
    qnRe = pnRe;
    qnIm = pnIm;
    qnnRe = p2nRe;
    qnnIm = p2nIm;

    // 2-particle differential flow
    W2p = mp * Mult - mq;
    avg2p = pnRe*QnRe + pnIm*QnIm - mq;
    avg2p /= W2p;
    w2pavg2psq = W2p * avg2p * avg2p;
    w2w2pavg2avg2p = W2 * W2p * avg2 * avg2p;
    
    flowHist->Fill(Eta, 1, W2p * avg2p);
    flowHist->Fill(Eta, 3, W2p);
    flowHist->Fill(Eta, 4, w2pavg2psq);
    flowHist->Fill(Eta, 6, w2w2pavg2avg2p);
    flowHist->Fill(Eta, 8, W2 * W2p);

    flowHist->Fill(Eta, 26, pnRe);
    flowHist->Fill(Eta, 27, pnIm);
    flowHist->Fill(Eta, 28, mp);

    // 4-particle differential flow
    W4p = (mp * Mult - 3.*mq)*(Mult - 1.)*(Mult - 2.);

    avg4p =  pnRe*QnRe*(QnRe*QnRe + QnIm*QnIm) + pnIm*QnIm*(QnRe*QnRe + QnIm*QnIm);
    avg4p -= qnnRe*QnRe*QnRe - qnnRe*QnIm*QnIm + 2.*qnnIm*QnRe*QnIm;
    avg4p -= pnRe*QnRe*Q2nRe - pnRe*QnIm*Q2nIm + pnIm*QnRe*Q2nIm + pnIm*QnIm*Q2nRe;
    avg4p -= 2.*Mult*(pnRe*QnRe + pnIm*QnIm);

    avg4p += - 2.*mq*(QnRe*QnRe + QnIm*QnIm) + 7.*(qnRe*QnRe + qnIm*QnIm); 
    avg4p += - (QnRe*qnRe + QnIm*qnIm) + (qnnRe*Q2nRe + qnnIm*Q2nIm);
    avg4p += 2.*(pnRe*QnRe + pnIm*QnIm) + 2.*mq*Mult - 6.*mq;
    avg4p /= W4p;

    w4pavg4psq = W4p * avg4p * avg4p;
    w2w4pavg2avg4p = W2 * W4p * avg2 * avg4p;
    w4w2pavg4avg2p = W4 * W2p * avg4 * avg2p;
    w4w4pavg4avg4p = W4 * W4p * avg4 * avg4p;
    w2pw4pavg2pavg4p = W2p * W4p * avg2p * avg4p;

    flowHist->Fill(Eta, 9, W4p * avg4p);
    flowHist->Fill(Eta, 11, W4p);
    flowHist->Fill(Eta, 12, w4pavg4psq);
    flowHist->Fill(Eta, 14, W2p * W4p);
    flowHist->Fill(Eta, 15, w2pw4pavg2pavg4p);
    flowHist->Fill(Eta, 17, W2 * W4p);
    flowHist->Fill(Eta, 18, w2w4pavg2avg4p);
    flowHist->Fill(Eta, 20, W4 * W2p);
    flowHist->Fill(Eta, 21, w4w2pavg4avg2p);
    flowHist->Fill(Eta, 23, W4 * W4p);
    flowHist->Fill(Eta, 24, w4w4pavg4avg4p);

    CosPhi1Phi2 = pnRe*QnRe - pnIm*QnIm - qnnRe;
    SinPhi1Phi2 = pnRe*QnIm + pnIm*QnRe - qnnIm;

    CosPhi1Phi2Phi3p =  pnRe*(QnRe*QnRe + QnIm*QnIm - Mult);
    CosPhi1Phi2Phi3p -= qnnRe*QnRe - qnnIm*QnIm + mq*QnRe - 2.*qnRe;
    SinPhi1Phi2Phi3p =  pnIm*(QnRe*QnRe + QnIm*QnIm - Mult);
    SinPhi1Phi2Phi3p -= qnnIm*QnRe - qnnRe*QnIm + mq*QnIm - 2.*qnIm;

    CosPhi1Phi2Phi3m =  pnRe*(QnRe*QnRe - QnIm*QnIm) + 2.*pnIm*QnRe*QnIm;
    CosPhi1Phi2Phi3m -= pnRe*Q2nRe + pnIm*Q2nIm + 2.*mq*QnRe - 2.*qnRe;
    SinPhi1Phi2Phi3m =  pnIm*(QnRe*QnRe - QnIm*QnIm) - 2.*pnRe*QnRe*QnIm;
    SinPhi1Phi2Phi3m += - pnIm*Q2nRe + pnRe*Q2nIm + 2.*mq*QnIm - 2.*qnIm;

    flowHist->Fill(Eta, 29, CosPhi1Phi2);
    flowHist->Fill(Eta, 30, SinPhi1Phi2);
    flowHist->Fill(Eta, 31, CosPhi1Phi2Phi3m);
    flowHist->Fill(Eta, 32, SinPhi1Phi2Phi3m);
    flowHist->Fill(Eta, 33, (mp*Mult-2.*mq)*(Mult-1.));
    flowHist->Fill(Eta, 34, CosPhi1Phi2Phi3p);
    flowHist->Fill(Eta, 35, SinPhi1Phi2Phi3p); 

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

  TH2D* cumulantsHist;
  TH1D* cumulant2Hist; 
  TH1D* cumulant4Hist; 

  Int_t nLoops = (fMC ? 3 : 1);

  // Do a loop over the difference analysis types, calculating flow
  // 1 loop for real data, 3 for MC data
  // inside each is a nested loop over each harmonic (1, 2, 3 and 4 at the moment)
  for (Int_t loop = 1; loop <= nLoops; loop++) {

    TString type;
    if (loop == 1) type = "";
    if (loop == 2) type = "MC";
    if (loop == 3) type = "TrRef";
    
    for (Int_t n = 1; n <= 4; n++) {
      if (!fv[n]) continue;

      cumulantsHist = (TH2D*)fOutputList->FindObject(Form("hQ%dCumuHist%s", n, type.Data()));
      cumulant2Hist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant2Flow%s", n, type.Data()));
      cumulant4Hist = (TH1D*)fOutputList->FindObject(Form("hQ%dCumulant4Flow%s", n, type.Data()));
  
      // For flow calculations
      Double_t Avg2 = 0, c2 = 0, vTwo2 = 0, Avg4 = 0, c4 = 0, vTwo4 = 0; 
      Double_t Avg2p = 0, d2 = 0, vTwo2diff = 0, Avg4p = 0, d4 = 0, vTwo4diff = 0;
      Double_t W2 = 0, W4 = 0, W2p = 0, W4p = 0, sqrtW2sq = 0, sqrtW2psq = 0, W2W2p = 0;
      Double_t W2W4 = 0, W2W4p = 0, W4W2p = 0, W4W4p = 0, W2pW4p = 0;
      Double_t sqrtW4sq = 0, sqrtW4psq = 0;
      Double_t W2avg2 = 0, W2pavg2p = 0, W4avg4 = 0, W4pavg4p = 0;
      Double_t AvgCos2Phi = 0, AvgSin2Phi = 0, Mult = 0, AvgCos2Phi1Phi2 = 0, AvgSin2Phi1Phi2 = 0;
      Double_t AvgCos2Phi1Phi2Phi3m = 0, AvgSin2Phi1Phi2Phi3m = 0, Multm1m2 = 0;
      Double_t AvgCos2Psi = 0, AvgSin2Psi = 0, mp = 0, AvgCos2Psi1Phi2 = 0, AvgSin2Psi1Phi2 = 0;
      Double_t AvgCos2Psi1Phi2Phi3m = 0, AvgSin2Psi1Phi2Phi3m = 0, mpqMult = 0;
      Double_t AvgCos2Psi1Phi2Phi3p = 0, AvgSin2Psi1Phi2Phi3p = 0;

      // For error calculations
      Double_t W2avg2sq = 0, W2W2pavg2avg2p = 0, W2pavg2psq = 0;
      Double_t W4avg4sq = 0, W2W4avg2avg4 = 0, W4W4pavg4avg4p = 0, W4pavg4psq = 0;
      Double_t W2W4pavg2avg4p = 0, W4W2pavg4avg2p = 0;
      Double_t sAvg2sq = 0, sAvg2psq = 0, sAvg4sq = 0, sAvg4psq = 0;
      Double_t vTwo2err = 0, vTwo2diffErr = 0, vTwo4err = 0, vTwo4diffErr = 0;
      Double_t Cov22p = 0, Cov24 = 0, Cov24p = 0, Cov42p = 0, Cov44p = 0, Cov2p2np = 0;
      Double_t W2pW4pavg2pavg4p = 0;

      for (Int_t eta = 0; eta <= cumulantsHist->GetNbinsX(); eta++) {
        if (eta == 0) {
          // 2-particle reference flow
          W2avg2 = cumulantsHist->GetBinContent(eta, 1);
          if (!W2avg2) continue;
          W2 = cumulantsHist->GetBinContent(eta, 3);
          AvgCos2Phi = cumulantsHist->GetBinContent(eta, 26);
          AvgSin2Phi = cumulantsHist->GetBinContent(eta, 27);
          Mult = cumulantsHist->GetBinContent(eta, 28);
          AvgCos2Phi /= Mult;
          AvgSin2Phi /= Mult;
          Avg2 = W2avg2 / W2;
          c2 = Avg2 - TMath::Power(AvgCos2Phi, 2) - TMath::Power(AvgSin2Phi, 2); 
          vTwo2 = TMath::Sqrt(c2);
          cumulant2Hist->SetBinContent(eta, vTwo2);

          // 4-particle reference flow
          W4avg4 = cumulantsHist->GetBinContent(eta, 9);
          W4 = cumulantsHist->GetBinContent(eta, 11);
          AvgCos2Phi1Phi2 = cumulantsHist->GetBinContent(eta, 29);
          AvgSin2Phi1Phi2 = cumulantsHist->GetBinContent(eta, 30);
          AvgCos2Phi1Phi2 /= W2;
          AvgSin2Phi1Phi2 /= W2;
          AvgCos2Phi1Phi2Phi3m = cumulantsHist->GetBinContent(eta, 31);
          AvgSin2Phi1Phi2Phi3m = cumulantsHist->GetBinContent(eta, 32);
          Multm1m2 = cumulantsHist->GetBinContent(eta, 33);
          AvgCos2Phi1Phi2Phi3m /= Multm1m2;
          AvgSin2Phi1Phi2Phi3m /= Multm1m2;
          Avg4 = W4avg4 / W4;
          c4 = Avg4 - 2. * Avg2 * Avg2;
          c4 -= 4.*AvgCos2Phi*AvgCos2Phi1Phi2Phi3m;
          c4 += 4.*AvgSin2Phi*AvgSin2Phi1Phi2Phi3m; 
          c4 -= TMath::Power(AvgCos2Phi1Phi2, 2) + TMath::Power(AvgSin2Phi1Phi2 , 2);
          c4 += 4.*AvgCos2Phi1Phi2*(TMath::Power(AvgCos2Phi, 2) - TMath::Power(AvgSin2Phi, 2));
          c4 += 8.*AvgSin2Phi1Phi2*AvgSin2Phi*AvgCos2Phi;  
          c4 += 8.*Avg2*(TMath::Power(AvgCos2Phi, 2) + TMath::Power(AvgSin2Phi, 2));
          c4 -= 6.*TMath::Power(TMath::Power(AvgCos2Phi, 2)+TMath::Power(AvgSin2Phi, 2), 2);

          vTwo4 = TMath::Power(-c4, 0.25);
          cumulant4Hist->SetBinContent(eta, vTwo4);
 
          // 2-particle reference flow error calculations
          W2avg2sq = cumulantsHist->GetBinContent(eta, 4);
          sqrtW2sq = cumulantsHist->GetBinError(eta, 3);
  
          sAvg2sq = VarSQ(W2avg2sq, Avg2, W2, W2avg2, sqrtW2sq);
          vTwo2err = sqrtW2sq * TMath::Sqrt(sAvg2sq) / (2. * TMath::Sqrt(Avg2) * W2);
          cumulant2Hist->SetBinError(eta, vTwo2err);
  
          // 4-particle reference flow error calculations
          W4avg4sq = cumulantsHist->GetBinContent(eta, 12);
          sqrtW4sq = cumulantsHist->GetBinError(eta, 11);
          W2W4 = cumulantsHist->GetBinContent(eta, 14);
          W2W4avg2avg4 = cumulantsHist->GetBinContent(eta, 15);
  
          sAvg4sq = VarSQ(W4avg4sq, Avg4, W4, W4avg4, sqrtW4sq);
          Cov24 = CovXY(W2W4avg2avg4, W2W4, Avg2*Avg4, W2, W4);
  
          vTwo4err = Avg2*Avg2 * TMath::Power(sqrtW2sq, 2) * sAvg2sq / (W2*W2);
          vTwo4err += TMath::Power(sqrtW4sq, 2) * sAvg4sq / (16. * W4*W4);
          vTwo4err -= Avg2 * W2W4 * Cov24 / (2. * W2 * W4);
          vTwo4err /= TMath::Power(2. * Avg2*Avg2 - Avg4, 1.5);
          vTwo4err = TMath::Sqrt(vTwo4err);
          cumulant4Hist->SetBinError(eta, vTwo4err);
 
          continue;
        }

        // 2-particle differential flow
        W2pavg2p = cumulantsHist->GetBinContent(eta, 1);
        if (!W2pavg2p) continue;
        W2p = cumulantsHist->GetBinContent(eta, 3);
        AvgCos2Psi = cumulantsHist->GetBinContent(eta, 26);
        AvgSin2Psi = cumulantsHist->GetBinContent(eta, 27);
        mp = cumulantsHist->GetBinContent(eta, 28);
        AvgCos2Psi /= mp;
        AvgSin2Psi /= mp;
        Avg2p = W2pavg2p / W2p;
        d2 = Avg2p - AvgCos2Psi*AvgCos2Phi - AvgSin2Psi*AvgSin2Phi; 
        vTwo2diff = d2 / TMath::Sqrt(c2);
        cumulant2Hist->SetBinContent(eta, vTwo2diff);
 
        // 4-particle differential flow
        W4pavg4p = cumulantsHist->GetBinContent(eta, 9);
        W4p = cumulantsHist->GetBinContent(eta, 11);
        AvgCos2Psi1Phi2 = cumulantsHist->GetBinContent(eta, 29);
        AvgSin2Psi1Phi2 = cumulantsHist->GetBinContent(eta, 30);
        AvgCos2Psi1Phi2 /= W2p;
        AvgSin2Psi1Phi2 /= W2p;
        AvgCos2Psi1Phi2Phi3m = cumulantsHist->GetBinContent(eta, 31);
        AvgSin2Psi1Phi2Phi3m = cumulantsHist->GetBinContent(eta, 32);
        mpqMult = cumulantsHist->GetBinContent(eta, 33);
        AvgCos2Psi1Phi2Phi3m /= mpqMult;
        AvgSin2Psi1Phi2Phi3m /= mpqMult;
        AvgCos2Psi1Phi2Phi3p = cumulantsHist->GetBinContent(eta, 34);
        AvgSin2Psi1Phi2Phi3p = cumulantsHist->GetBinContent(eta, 35);
        AvgCos2Psi1Phi2Phi3p /= mpqMult;
        AvgSin2Psi1Phi2Phi3p /= mpqMult;

        Avg4p = W4pavg4p / W4p;
        d4 = Avg4p - 2. * Avg2p * Avg2;
        d4 -= AvgCos2Psi*AvgCos2Phi1Phi2Phi3m; 
        d4 += AvgSin2Psi*AvgSin2Phi1Phi2Phi3m; 
        d4 -= AvgCos2Phi*AvgCos2Psi1Phi2Phi3m; 
        d4 += AvgSin2Phi*AvgSin2Psi1Phi2Phi3m; 
        d4 -= 2.*AvgCos2Phi*AvgCos2Psi1Phi2Phi3p;
        d4 -= 2.*AvgSin2Phi*AvgSin2Psi1Phi2Phi3p; 
        d4 -= AvgCos2Psi1Phi2*AvgCos2Phi1Phi2; 
        d4 -= AvgSin2Psi1Phi2*AvgSin2Phi1Phi2; 
        d4 += 2.*AvgCos2Phi1Phi2*(AvgCos2Psi*AvgCos2Phi - AvgSin2Psi*AvgSin2Phi);  
        d4 += 2.*AvgSin2Phi1Phi2*(AvgCos2Psi*AvgSin2Phi + AvgSin2Psi*AvgCos2Phi); 
        d4 += 4.*Avg2*(AvgCos2Psi*AvgCos2Phi + AvgSin2Psi*AvgSin2Phi);
        d4 += 2.*AvgCos2Psi1Phi2*(TMath::Power(AvgCos2Phi, 2) - TMath::Power(AvgSin2Phi, 2)); 
        d4 += 4.*AvgSin2Psi1Phi2*AvgCos2Phi*AvgSin2Phi;
        d4 += 4.*Avg2p*(TMath::Power(AvgCos2Phi, 2) + TMath::Power(AvgSin2Phi, 2)); 
        d4 -= 6.*(TMath::Power(AvgCos2Phi, 2) - TMath::Power(AvgSin2Phi, 2))
              *(AvgCos2Psi*AvgCos2Phi-AvgSin2Psi*AvgSin2Phi); 
        d4 -= 12.*AvgCos2Phi*AvgSin2Phi*(AvgSin2Psi*AvgCos2Phi+AvgCos2Psi*AvgSin2Phi); 
 
        vTwo4diff = - d4 / TMath::Power(-c4, 0.75);
        cumulant4Hist->SetBinContent(eta, vTwo4diff);    
      
        // 2-particle differential flow error calculations
        W2pavg2psq = cumulantsHist->GetBinContent(eta, 4);
        sqrtW2psq = cumulantsHist->GetBinError(eta, 3);
        W2W2pavg2avg2p = cumulantsHist->GetBinContent(eta, 6);
        W2W2p = cumulantsHist->GetBinContent(eta, 8);
  
        Cov22p = CovXY(W2W2pavg2avg2p, W2W2p, Avg2*Avg2p, W2, W2p);
        sAvg2psq = VarSQ(W2pavg2psq, Avg2p, W2p, W2pavg2p, sqrtW2psq);
 
        vTwo2diffErr = Avg2p*Avg2p*TMath::Power(sqrtW2psq, 2)*sAvg2sq/(W2*W2);
        vTwo2diffErr += 4.*Avg2*Avg2*TMath::Power(sqrtW2psq, 2)*sAvg2psq/(W2p*W2p);
        vTwo2diffErr -= 4.*Avg2*Avg2p*W2W2p*Cov22p/(W2*W2p);
        vTwo2diffErr /= (4. * TMath::Power(Avg2, 3));
        vTwo2diffErr = TMath::Sqrt(vTwo2diffErr);
        cumulant2Hist->SetBinError(eta, vTwo2diffErr);

        // 4-particle differential flow error calculations
        sqrtW4psq = cumulantsHist->GetBinError(eta, 11);
        W4pavg4psq = cumulantsHist->GetBinContent(eta, 12);
        W2pW4p = cumulantsHist->GetBinContent(eta, 14);
        W2pW4pavg2pavg4p = cumulantsHist->GetBinContent(eta, 15);
        W2W4p = cumulantsHist->GetBinContent(eta, 17);
        W2W4pavg2avg4p = cumulantsHist->GetBinContent(eta, 18);
        W4W2p = cumulantsHist->GetBinContent(eta, 20);
        W4W2pavg4avg2p = cumulantsHist->GetBinContent(eta, 21);
        W4W4p = cumulantsHist->GetBinContent(eta, 23);
        W4W4pavg4avg4p = cumulantsHist->GetBinContent(eta, 24);
      
        sAvg4psq = VarSQ(W4pavg4psq, Avg4p, W4p, W4pavg4p, sqrtW4psq);
        Cov24p = CovXY(W2W4pavg2avg4p, W2W4p, Avg2*Avg4p, W2, W4p);
        Cov42p = CovXY(W4W2pavg4avg2p, W4W2p, Avg4*Avg2p, W4, W2p);
        Cov44p = CovXY(W4W4pavg4avg4p, W4W4p, Avg4*Avg4p, W4, W4p);
        Cov2p2np = CovXY(W2pW4pavg2pavg4p, W2pW4p, Avg2p*Avg4p, W2p, W4p);
 
        // Numbers on the side reference term number in paper (cite needed) loosely
/*1*/ vTwo4diffErr =  TMath::Power(2.*Avg2*Avg2*Avg2p - 3.*Avg2*Avg4p + 2.*Avg4*Avg2p, 2) 
                      * TMath::Power(sqrtW2sq, 2) * sAvg2sq / (W2*W2);
/*2*/ vTwo4diffErr += 9. * TMath::Power(2.*Avg2*Avg2p - Avg4p, 2) * TMath::Power(sqrtW4sq, 2)
                      * sAvg4sq / (16. * W4*W4);
/*3*/ vTwo4diffErr += 4. * Avg2*Avg2 * TMath::Power(2.*Avg2*Avg2 - Avg4, 2) * TMath::Power(sqrtW2psq, 2)
                      * sAvg2psq / (W2p*W2p);
/*4*/ vTwo4diffErr += TMath::Power(2.*Avg2*Avg2 - Avg4, 2) * TMath::Power(sqrtW4psq, 2) * sAvg4psq
                      / (W4p*W4p);
/*5*/ vTwo4diffErr -= 1.5 * (2.*Avg2*Avg2p - Avg4p) * (2.*Avg2*Avg2*Avg2p - 3.*Avg2*Avg4p + 2.*Avg4*Avg2p)
                      * W2W4 * Cov24 / (W2*W4);
/*6*/ vTwo4diffErr -= 4. * Avg2 * (2.*Avg2*Avg2 - Avg4) 
                      * (2.*Avg2*Avg2*Avg2p - 3.*Avg2*Avg4p + 2.*Avg4*Avg2p)
                      * W2W2p * Cov22p / (W2 * W2p);
/*7*/ vTwo4diffErr += 2. * (2.*Avg2*Avg2 - Avg4)
                      * (2.*Avg2*Avg2*Avg2p - 3.*Avg2*Avg4p + 2.*Avg4*Avg2p)
                      * W2W4p * Cov24p / (W2 * W4p);
/*8*/ vTwo4diffErr += 3.*Avg2*(2.*Avg2*Avg2 - Avg4)*(2.*Avg2*Avg2p - Avg4p)
                      * W4W2p * Cov42p / (W4*W2p);
/*9*/ vTwo4diffErr -= 1.5 * (2.*Avg2*Avg2 - Avg4)*(2.*Avg2*Avg2p - Avg4p)
                      * W4W4p * Cov44p / (W4 * W4p);
/*10*/vTwo4diffErr -= 4.*Avg2*TMath::Power(2.*Avg2*Avg2 - Avg4, 2)
                      * W2pW4p * Cov2p2np / (W2p * W4p);
/*11*/vTwo4diffErr /= TMath::Power(2.*Avg2*Avg2 - Avg4, 3.5);
      vTwo4diffErr = TMath::Sqrt(vTwo4diffErr);

        cumulant4Hist->SetBinError(eta, vTwo4diffErr);
      } // End of eta loop

      // Number of events:
      Int_t nEv = cumulantsHist->GetBinContent(0,0);
      cumulant2Hist->SetBinContent(cumulant2Hist->GetNbinsX() + 1, nEv);
      cumulant4Hist->SetBinContent(cumulant4Hist->GetNbinsX() + 1, nEv);
    } // End of harmonics loop
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

  // Histograms are loaded and reset
  TH1D* dNdphi          = (TH1D*)fOutputList->FindObject("hdNdphiSEMC");
  TH2D* dNdetadphi      = (TH2D*)fOutputList->FindObject("hdNdetadphiSEMC");
  TH1D* dNdphiTrRef     = (TH1D*)fOutputList->FindObject("hdNdphiSETrRef");
  TH2D* dNdetadphiTrRef = (TH2D*)fOutputList->FindObject("hdNdetadphiSETrRef");

  dNdphi->Reset();
  dNdetadphi->Reset();
  dNdphiTrRef->Reset();
  dNdetadphiTrRef->Reset();

  // Loads AliFMDFlowCommon and fills histograms and runs analysis.
  // AOD events also get a TrackRef histogram
  AliForwardFlowBase* common = new AliForwardFlowBase(fOutputList);

  if (fAOD) {
    if (!common->LoopAODmc(fAOD)) return;
    if (!common->LoopAODtrrefHits(fAOD)) return;
//    if (!common->LoopMCaddptFlow(fAOD)) return;
//    if (!common->LoopMCaddpdgFlow(fAOD)) return;
//    if (!common->LoopMCaddetaFlow(fAOD)) return;
  }

  // Run analysis on MC truth
  for (Int_t n = 1; n <= 4; n++) {
    if (fv[n])
      CumulantsMethod("MC", n);
  }
  
  // Run analysis on TrackRefs
  for (Int_t n = 1; n <= 4; n++) {
    if (fv[n])
      CumulantsMethod("TrRef", n);
  }

}
//_____________________________________________________________________
Double_t AliForwardFlowTaskQC::VarSQ(Double_t wxx2, Double_t x, Double_t wx, Double_t wxx, Double_t sqrtwx2)
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
Double_t AliForwardFlowTaskQC::CovXY(Double_t wxwyxy, Double_t wxwy, Double_t xy, Double_t wx, Double_t wy)
{
  //
  // Small function to compute the covariance between two numbers
  // - used by Terminate()
  //
  Double_t Cov, denominator, numerator;

  denominator = (wxwyxy / wxwy) - xy;
  numerator = 1 - (wxwy / (wx * wy));
  
  Cov = denominator / numerator;
  return Cov;
}
//_____________________________________________________________________
//
//
// EOF
