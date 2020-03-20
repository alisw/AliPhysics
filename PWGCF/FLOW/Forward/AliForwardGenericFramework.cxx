#include "AliForwardGenericFramework.h"
#include "TMath.h"
#include <iostream>
#include "TRandom.h"
#include "AliForwardSettings.h"
#include "TH1D.h"
#include <complex>
#include <cmath>
#include "TFile.h"

using namespace std;

//_____________________________________________________________________
AliForwardGenericFramework::AliForwardGenericFramework(Int_t refbins):
  fSettings(),
  fQvector(),
  fpvector(),
  fqvector(),
  cumu_rW2(),
  cumu_rW2Two(),
  cumu_rW4(),
  cumu_rW4Four(),
  cumu_dW2B(),
  cumu_dW2TwoB(),
  cumu_dW4(),
  cumu_dW4Four(),
  cumu_dW2TwoTwoN(),
  cumu_dW2TwoTwoD(),
  cumu_dW4ThreeTwo(),
  cumu_dW4FourTwo(),
  cumu_wSC(),
  cumu_dW2A(),
  cumu_dW2TwoA(),
  cumu_dW22TwoTwoN(),
  cumu_dW22TwoTwoD()
{
  Int_t rbins[4] = {2, 6, 4, 2} ; // kind (real or imaginary), n, p, eta
  Int_t dimensions = 4;
  rbins[3] = refbins; // two bins in eta for gap, one for standard

  Double_t xmin[4] = {-1.0, -0.5, 0.5, -fSettings.fEtaUpEdge}; // kind (real or imaginary), n, p, eta
  Double_t xmax[4] = { 1,   5.5, 4.5,   fSettings.fEtaUpEdge}; // kind (real or imaginary), n, p, eta SKAL VAERE -6 - 6

  fQvector = new THnD("Qvector", "Qvector", dimensions, rbins, xmin, xmax);

  Int_t dbins[4] = {2, 6, 4, fSettings.fNDiffEtaBins} ; // kind (real or imaginary), n, p, eta

  Double_t dxmin[4] = {-1.0, -0.5, 0.5, fSettings.fEtaLowEdge}; // kind (real or imaginary), n, p, eta
  Double_t dxmax[4] = { 1,   5.5, 4.5,  fSettings.fEtaUpEdge}; // kind (real or imaginary), n, p, eta SKAL VAERE -6 - 6


  fpvector = new THnD("pvector", "pvector", dimensions, dbins, dxmin, dxmax);
  fqvector = new THnD("qvector", "qvector", dimensions, dbins, dxmin, dxmax);

  fAutoRef =  new TH1D("fAutoRef", "fAutoRef", refbins, fSettings.fEtaLowEdge, fSettings.fEtaUpEdge);
  fAutoDiff = new TH1D("fAutoDiff","fAutoDiff", fSettings.fNDiffEtaBins, fSettings.fEtaLowEdge, fSettings.fEtaUpEdge);
  fAutoRef->SetDirectory(0);
  fAutoDiff->SetDirectory(0);
}


//_____________________________________________________________________
void AliForwardGenericFramework::CumulantsAccumulate(TH2D*& dNdetadphi, double cent, double zvertex, Bool_t useFMD, Bool_t doRefFlow, Bool_t doDiffFlow)
{
  for (Int_t etaBin = 1; etaBin <= dNdetadphi->GetNbinsX(); etaBin++) {
    if ((!fSettings.use_primaries_fwd && !fSettings.esd) && useFMD){
      if (dNdetadphi->GetBinContent(etaBin, 0) == 0) continue; // No data expected for this eta 
    }
    Double_t eta = dNdetadphi->GetXaxis()->GetBinCenter(etaBin);
    if (eta > fSettings.fEtaUpEdge || eta < fSettings.fEtaLowEdge) continue;
    Double_t difEtaBin = fpvector->GetAxis(3)->FindBin(eta);
    Double_t difEta = fpvector->GetAxis(3)->GetBinCenter(difEtaBin);

    Double_t refEtaBin = fQvector->GetAxis(3)->FindBin(eta);
    Double_t refEta = fQvector->GetAxis(3)->GetBinCenter(refEtaBin);

    for (Int_t phiBin = 1; phiBin <= dNdetadphi->GetNbinsY(); phiBin++) {

      Double_t phi = dNdetadphi->GetYaxis()->GetBinCenter(phiBin);
      Double_t mult_weight = dNdetadphi->GetBinContent(etaBin, phiBin);
      Double_t weight = mult_weight;

      if (fSettings.doNUA){
        if (useFMD) weight = applyNUAforward(dNdetadphi, etaBin, phiBin, eta, phi, zvertex, weight);
        else weight = applyNUAcentral(eta, phi, zvertex, weight);
      }
      
      if (!weight || weight == 0) continue;
      for (Int_t n = 0; n <= 4; n++) {

        if ((useFMD && fSettings.sec_corr) && n >=2) weight = applySecondaryCorr(n, eta, zvertex, cent, weight);

      if (!weight || weight == 0) continue;

        for (Int_t p = 1; p <= 4; p++) {
          Double_t realPart = TMath::Power(weight, p)*TMath::Cos(n*phi);
          Double_t imPart =   TMath::Power(weight, p)*TMath::Sin(n*phi);

          Double_t re[4] = {0.5, Double_t(n), Double_t(p), difEta};
          Double_t im[4] = {-0.5, Double_t(n), Double_t(p), difEta};

          if (doDiffFlow){
            fpvector->Fill(re, realPart);
            fpvector->Fill(im, imPart);

            if ((useFMD & !(fSettings.etagap)) ||
                (!(useFMD) && (fSettings.ref_mode & fSettings.kTPCref)) ||
                (useFMD && (fSettings.ref_mode & fSettings.kFMDref))) {
              fqvector->Fill(re, realPart);
              fqvector->Fill(im, imPart);
              if ((mult_weight > 1.0) & !fSettings.etagap) fAutoDiff->Fill(refEta,weight*(weight - 1));
            }
          }

          if (doRefFlow){
            if ((fSettings.etagap) && TMath::Abs(eta)<=fSettings.gap) continue;

            if (fSettings.ref_mode & fSettings.kTPCref) {           
              if ((fSettings.TPC_maxeta > 0) & (TMath::Abs(eta) > fSettings.TPC_maxeta)) continue;
            }
            if (fSettings.ref_mode & fSettings.kFMDref) {
              if (TMath::Abs(eta) < fSettings.fmdlowcut) continue;
              if (TMath::Abs(eta) > fSettings.fmdhighcut) continue;
            }

            Double_t req[4] = {0.5, static_cast<Double_t>(n), static_cast<Double_t>(p), refEta};
            Double_t imq[4] = {-0.5, static_cast<Double_t>(n), static_cast<Double_t>(p), refEta};
            if ((mult_weight > 1.0) & !fSettings.etagap) fAutoRef->Fill(refEta,weight*(weight - 1));  
            fQvector->Fill(req, realPart);
            fQvector->Fill(imq, imPart);
          }
        } // end p loop
      } // End of n loop
    } // End of phi loop
  } // end of eta
  return;
}


void AliForwardGenericFramework::saveEvent(double cent, double zvertex,UInt_t r, Int_t ptn){
  // For each n we loop over the hists
  Double_t sample = static_cast<Double_t>(r);
  Int_t refEtaBinA;
  Int_t refEtaBinB;
  Int_t etaBinB;
  Double_t refEtaA;

  for (Int_t n = 2; n <= 4; n++) {
    Int_t prevRefEtaBin = kTRUE;
    Int_t prevbin = 0;

    for (Int_t etaBin = 1; etaBin <= fpvector->GetAxis(3)->GetNbins(); etaBin++) {
      Double_t eta = fpvector->GetAxis(3)->GetBinCenter(etaBin);

      refEtaBinA = fQvector->GetAxis(3)->FindBin(eta);
      refEtaA    = fQvector->GetAxis(3)->GetBinCenter(refEtaBinA);
      refEtaBinB = refEtaBinA;
      etaBinB    = etaBin;

      if ((fSettings.etagap)) {
        refEtaBinB = fQvector->GetAxis(3)->FindBin(-eta);
        etaBinB    = fpvector->GetAxis(3)->FindBin(-eta);
      }

      // index to get sum of weights
      Int_t index1[4] = {2, 1, 1, refEtaBinB};
      if (!(fQvector->GetBinContent(index1) > 0)) continue;
      // REFERENCE FLOW --------------------------------------------------------------------------------
      if (prevRefEtaBin & (prevbin != refEtaBinB)){ // only used once

        // two-particle cumulant
        if ((fSettings.normal_analysis || fSettings.SC_analysis) || fSettings.second_analysis){
          double two = Two(n, -n, refEtaBinA, refEtaBinB).Re();
          if (!fSettings.etagap) two += fAutoRef->GetBinContent(refEtaBinA);
          fill(cumu_rW2Two, n, ptn, sample, zvertex, refEtaA, cent, two);
          if (n==2){
            double dn2 = Two(0,0, refEtaBinA, refEtaBinB).Re();
            if (!fSettings.etagap) dn2 += fAutoRef->GetBinContent(refEtaBinA);
            fill(cumu_rW2, -n, ptn, sample, zvertex, refEtaA, cent, dn2);
          }
        }
        if (fSettings.normal_analysis){
          // four-particle cumulant
          double four = Four(n, n, -n, -n, refEtaBinA, refEtaBinB).Re();

          fill(cumu_rW4Four, n, ptn, sample, zvertex, refEtaA, cent, four);
          if (n==2){
            double dn4 = Four(0,0,0,0 , refEtaBinA, refEtaBinB).Re();
            fill(cumu_rW4, -n, ptn, sample, zvertex, refEtaA, cent, dn4);
          }
        }
        if (fSettings.etagap) prevRefEtaBin = kFALSE;
        prevbin = refEtaBinB;
      }

      // DIFFERENTIAL FLOW -----------------------------------------------------------------------------
      if ((fSettings.normal_analysis || fSettings.decorr_analysis) || fSettings.second_analysis){
        if (n==2){
          double dn2diff = TwoDiff(0,0, refEtaBinB, etaBin).Re();
          if (!fSettings.etagap) dn2diff += fAutoDiff->GetBinContent(etaBin);
          fill(cumu_dW2B, -n, ptn, sample, zvertex, eta, cent, dn2diff);          
        }
        double twodiff = TwoDiff(n, -n, refEtaBinB, etaBin).Re();
        if (!fSettings.etagap) twodiff += fAutoDiff->GetBinContent(etaBin);
        fill(cumu_dW2TwoB, n, ptn, sample, zvertex, eta, cent, twodiff);
      }
      if (fSettings.normal_analysis){
        if (n==2){
          double dn2diff = TwoDiff(0,0, refEtaBinA, etaBin).Re();
          fill(cumu_dW2A, -n, ptn, sample, zvertex, eta, cent, dn2diff);
        }
        // A side
        double twodiff = TwoDiff(n, -n, refEtaBinA, etaBin).Re();
        fill(cumu_dW2TwoA, n, ptn, sample, zvertex, eta, cent, twodiff);
      }
      if (fSettings.normal_analysis){
        if (n==2){
          double dn4diff = FourDiff(0,0,0,0, refEtaBinA, refEtaBinB, etaBin,etaBin).Re();
          fill(cumu_dW4, -n, ptn, sample, zvertex, eta, cent, dn4diff);
        }
        // four-particle cumulant
        double fourdiff = FourDiff(n, n, -n, -n, refEtaBinA, refEtaBinB, etaBin,etaBin).Re(); // A is same side
        fill(cumu_dW4Four, n, ptn, sample, zvertex, eta, cent, fourdiff);
      }

      if (eta < 0.0) {
        if (fSettings.decorr_analysis){
          // R_{n,n; 2} numerator
          double over =  (TwoDiff(-n,n,refEtaBinA, etaBinB)*TwoDiff(n,-n,refEtaBinB, etaBin)).Re();
          double under = (TwoDiff(-n,n,refEtaBinA, etaBin)*TwoDiff(n,-n,refEtaBinB, etaBinB)).Re();
          fill(cumu_dW2TwoTwoN, n, ptn, sample, zvertex, eta, cent, over);
          fill(cumu_dW2TwoTwoD, n, ptn, sample, zvertex, eta, cent, under);
        }
      }

      if (fSettings.SC_analysis & (n==2)){
        double twotwodiffN = TwoTwoDiff(2,-2, etaBin,etaBinB).Re();
        fill(cumu_dW22TwoTwoN, -n, ptn, sample, zvertex, eta, cent, twotwodiffN);
        double twotwodiffD = TwoTwoDiff(0,0, etaBin,etaBinB).Re();
        fill(cumu_dW22TwoTwoD, -n, ptn, sample, zvertex, eta, cent, twotwodiffD);

        // four-particle cumulant SC(4,2)
        double wSC = FourDiff_SC(0,0,0,0, etaBin, refEtaBinA,etaBinB,refEtaBinB).Re();
        fill(cumu_wSC, -n, ptn, sample, zvertex, eta, cent, wSC);

        // four-particle cumulant SC(4,2)
        double fourtwodiff = FourDiff_SC(2,4,-2,-4, etaBin, refEtaBinA,etaBinB,refEtaBinB).Re();
        fill(cumu_dW4FourTwo, -n, ptn, sample, zvertex, eta, cent, fourtwodiff);

        // four-particle cumulant SC(3,2)
        double threetwodiff = FourDiff_SC(2,3,-2,-3,  etaBin, refEtaBinA,etaBinB,refEtaBinB).Re();
        fill(cumu_dW4ThreeTwo, -n, ptn, sample, zvertex, eta, cent, threetwodiff);
      } 
    } //eta
  } // moment
  return;
}


TComplex AliForwardGenericFramework::Q(Int_t n, Int_t p, Int_t etabin)
{
  double sign = (n < 0) ? -1 : 1;

  Int_t imindex[4] = {1, TMath::Abs(n)+1, p, etabin};
  Int_t reindex[4] = {2, TMath::Abs(n)+1, p, etabin};

  return TComplex(fQvector->GetBinContent(reindex),sign*fQvector->GetBinContent(imindex));;
}

TComplex AliForwardGenericFramework::p(Int_t n, Int_t p, Int_t etabin)
{
  double sign = (n > 0) ? 1 : ((n < 0) ? -1 : 1);

  Int_t imindex[4] = {1, TMath::Abs(n)+1, p, etabin};
  Int_t reindex[4] = {2, TMath::Abs(n)+1, p, etabin};

  return TComplex(fpvector->GetBinContent(reindex),sign*fpvector->GetBinContent(imindex));;
}


TComplex AliForwardGenericFramework::q(Int_t n, Int_t p, Int_t etabin)
{
  double sign = (n > 0) ? 1 : ((n < 0) ? -1 : 1);
  Int_t imindex[4] = {1, TMath::Abs(n)+1, p, etabin};
  Int_t reindex[4] = {2, TMath::Abs(n)+1, p, etabin};

  return TComplex(fqvector->GetBinContent(reindex),sign*fqvector->GetBinContent(imindex));;
}


TComplex AliForwardGenericFramework::Two(Int_t n1, Int_t n2, Int_t eta1, Int_t eta2)
{
  TComplex formula = 0;
  if (!fSettings.etagap) {
     formula = Q(n1,1,eta1)*Q(n2,1,eta1) - Q(n1+n2,2,eta1);
  }
  else{
     formula = Q(n1,1,eta1)*Q(n2,1,eta2);
  }
  return formula;
}


TComplex AliForwardGenericFramework::TwoDiff(Int_t n1, Int_t n2, Int_t refetabin, Int_t diffetabin)
{
  if (fSettings.etagap){
    return p(n1,1, diffetabin)*Q(n2,1, refetabin);
  }
  else return p(n1,1, diffetabin)*Q(n2,1, refetabin) - q(n1+n2,1, diffetabin);
}


TComplex AliForwardGenericFramework::TwoTwoDiff(Int_t n1, Int_t n2, Int_t diffetabin1, Int_t diffetabin2)
{
  return p(n1,1, diffetabin1)*p(n2,1, diffetabin2);
}


TComplex AliForwardGenericFramework::Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4,Int_t eta1, Int_t eta2)
{
  TComplex formula = 0;
  if (eta1 != eta2) {
    formula = Two(n1,n2,eta1,eta1)*Two(n3,n4,eta2,eta2);
  }
  else{

   formula = Q(n1,1,eta1)*Q(n2,1,eta1)*Q(n3,1,eta1)*Q(n4,1,eta1)-Q(n1+n2,2,eta1)*Q(n3,1,eta1)*Q(n4,1,eta1)-Q(n2,1,eta1)*Q(n1+n3,2,eta1)*Q(n4,1,eta1)
                    - Q(n1,1,eta1)*Q(n2+n3,2,eta1)*Q(n4,1,eta1)+2.*Q(n1+n2+n3,3,eta1)*Q(n4,1,eta1)-Q(n2,1,eta1)*Q(n3,1,eta1)*Q(n1+n4,2,eta1)
                    + Q(n2+n3,2,eta1)*Q(n1+n4,2,eta1)-Q(n1,1,eta1)*Q(n3,1,eta1)*Q(n2+n4,2,eta1)+Q(n1+n3,2,eta1)*Q(n2+n4,2,eta1)
                    + 2.*Q(n3,1,eta1)*Q(n1+n2+n4,3,eta1)-Q(n1,1,eta1)*Q(n2,1,eta1)*Q(n3+n4,2,eta1)+Q(n1+n2,2,eta1)*Q(n3+n4,2,eta1)
                    + 2.*Q(n2,1,eta1)*Q(n1+n3+n4,3,eta1)+2.*Q(n1,1,eta1)*Q(n2+n3+n4,3,eta1)-6.*Q(n1+n2+n3+n4,4,eta1);
  }
  return formula;
}


TComplex AliForwardGenericFramework::FourDiff(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t refetabinA, Int_t refetabinB, Int_t diffetabin,Int_t qetabin)
{
  TComplex formula = 0;
  if (!fSettings.etagap){
    formula = p(n1,1,diffetabin)*Q(n2,1,refetabinA)*Q(n3,1,refetabinA)*Q(n4,1,refetabinA)-q(n1+n2,2,qetabin)*Q(n3,1,refetabinA)*Q(n4,1,refetabinA)-Q(n2,1,refetabinA)*q(n1+n3,2,qetabin)*Q(n4,1,refetabinA)
            - p(n1,1,diffetabin)*Q(n2+n3,2,refetabinA)*Q(n4,1,refetabinA)+2.*q(n1+n2+n3,3,qetabin)*Q(n4,1,refetabinA)-Q(n2,1,refetabinA)*Q(n3,1,refetabinA)*q(n1+n4,2,qetabin)
            + Q(n2+n3,2,refetabinA)*q(n1+n4,2,qetabin)-p(n1,1,diffetabin)*Q(n3,1,refetabinA)*Q(n2+n4,2,refetabinA)+q(n1+n3,2,qetabin)*Q(n2+n4,2,refetabinA)
            + 2.*Q(n3,1,refetabinA)*q(n1+n2+n4,3,qetabin)-p(n1,1,diffetabin)*Q(n2,1,refetabinA)*Q(n3+n4,2,refetabinA)+q(n1+n2,2,qetabin)*Q(n3+n4,2,refetabinA)
            + 2.*Q(n2,1,refetabinA)*q(n1+n3+n4,3,qetabin)+2.*p(n1,1,diffetabin)*Q(n2+n3+n4,3,refetabinA)-6.*q(n1+n2+n3+n4,4,qetabin);
  }
  else{
    formula = p(n1,1,diffetabin)*Q(n2,1,refetabinA)*Q(n3,1,refetabinB)*Q(n4,1,refetabinB)
            - q(n1+n2,2,qetabin)*Q(n3,1,refetabinB)*Q(n4,1,refetabinB)
            - p(n1,1,diffetabin)*Q(n2,1,refetabinA)*Q(n3+n4,2,refetabinB)
            + q(n1+n2,2,qetabin)*Q(n3+n4,2,refetabinB);
  }
  return formula;
}


TComplex AliForwardGenericFramework::FourDiff_SC(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t diffetabinA, Int_t refetabinA, Int_t diffetabinB,Int_t refetabinB)
{
  TComplex formula = 0;
  formula = p(n1,1,diffetabinA)*Q(n2,1,refetabinA)*p(n3,1,diffetabinB)*Q(n4,1,refetabinB)
          - p(n1+n2,2,diffetabinA)*p(n3,1,diffetabinB)*Q(n4,1,refetabinB)
          - p(n1,1,diffetabinA)*Q(n2,1,refetabinA)*p(n3+n4,2,diffetabinB)
          + p(n1+n2,2,diffetabinA)*p(n3+n4,2,diffetabinB);
  return formula;
}

void AliForwardGenericFramework::reset() {
  fQvector->Reset();
  fpvector->Reset();
  fqvector->Reset();
  fAutoRef->Reset();
  fAutoDiff->Reset();
}
