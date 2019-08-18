#include "AliForwardGenericFramework.h"
#include "TMath.h"
#include <iostream>
#include "TRandom.h"
#include "AliForwardSettings.h"
#include "TH1D.h"
#include <complex>
#include <cmath>
#include "TH2F.h"
#include "TFile.h"
#include "AliForwardNUATask.h"

using namespace std;

//_____________________________________________________________________
AliForwardGenericFramework::AliForwardGenericFramework()
{
  Int_t rbins[4] = {2, 6, 4, 2} ; // kind (real or imaginary), n, p, eta
  Int_t dimensions = 4;
  if (fSettings.etagap) {
     rbins[3] = 2; // two bins in eta for gap, one for standard
  }

  Double_t xmin[4] = {-1.0, -0.5, 0.5, -6}; // kind (real or imaginary), n, p, eta
  Double_t xmax[4] = { 1,   5.5, 4.5,  6}; // kind (real or imaginary), n, p, eta SKAL VAERE -6 - 6

  fQvector = new THnD("Qvector", "Qvector", dimensions, rbins, xmin, xmax);

  Int_t dbins[4] = {2, 6, 4, fSettings.fNDiffEtaBins} ; // kind (real or imaginary), n, p, eta
  xmin[3] = -4.0;

  fpvector = new THnD("pvector", "pvector", dimensions, dbins, xmin, xmax);
  fqvector = new THnD("qvector", "qvector", dimensions, dbins, xmin, xmax);
  //fpvector->SetDirectory(0);
  //fqvector->SetDirectory(0);

  // fAutoRef = TH1F("fAutoRef","fAutoRef", 2, fSettings.fEtaLowEdge, fSettings.fEtaUpEdge);
  // fAutoDiff = TH1F("fAutoDiff","fAutoDiff", fSettings.fNDiffEtaBins, fSettings.fEtaLowEdge, fSettings.fEtaUpEdge);
  // fAutoRef.SetDirectory(0);
  // fAutoDiff.SetDirectory(0);
}


//_____________________________________________________________________
void AliForwardGenericFramework::CumulantsAccumulate(TH2D*& dNdetadphi, double cent, double zvertex, Bool_t useFMD, Bool_t doRefFlow, Bool_t doDiffFlow)
{

  for (Int_t etaBin = 1; etaBin <= dNdetadphi->GetNbinsX(); etaBin++) {

    Double_t eta = dNdetadphi->GetXaxis()->GetBinCenter(etaBin);
    Double_t difEtaBin = fpvector->GetAxis(3)->FindBin(eta);
    Double_t difEta = fpvector->GetAxis(3)->GetBinCenter(difEtaBin);

    Double_t refEtaBin = fQvector->GetAxis(3)->FindBin(eta);
    Double_t refEta = fQvector->GetAxis(3)->GetBinCenter(refEtaBin);

    for (Int_t phiBin = 1; phiBin <= dNdetadphi->GetNbinsY(); phiBin++) {
      Double_t phi = dNdetadphi->GetYaxis()->GetBinCenter(phiBin);
      Double_t weight = dNdetadphi->GetBinContent(etaBin, phiBin);
      //if (weight > 0) std::cout << "weight = " << weight << std::endl;

      if (!fSettings.use_primaries_fwd && !fSettings.esd){
        if (!fSettings.mc){
          if (dNdetadphi->GetBinContent(etaBin, 0) == 0 && useFMD) break;
        }
      }

      if (fSettings.doNUA){
        // holes in the FMD
        if ((fSettings.nua_mode & fSettings.kFill) && useFMD){
          if (etaBin >= 125 && etaBin <=137){
            if (phiBin == 17 || phiBin == 18) weight = 1.;
          }
          if (etaBin >= 168 && etaBin <=185){
            if (phiBin == 14) weight = 1.;
          }
        }

        if ((fSettings.nua_mode & fSettings.kInterpolate) && useFMD) weight = AliForwardNUATask::InterpolateWeight(dNdetadphi,phiBin,etaBin,weight);


        if (!useFMD && !fSettings.use_primaries_cen) {
          Int_t nuaeta = fSettings.nuacentral->GetXaxis()->FindBin(eta);
          Int_t nuaphi = fSettings.nuacentral->GetYaxis()->FindBin(phi);
          Int_t nuavtz = fSettings.nuacentral->GetZaxis()->FindBin(zvertex);
          weight = weight*fSettings.nuacentral->GetBinContent(nuaeta,nuaphi,nuavtz);
        }

        if (useFMD && !fSettings.use_primaries_fwd) {
          Int_t nuaeta = fSettings.nuaforward->GetXaxis()->FindBin(eta);
          Int_t nuaphi = fSettings.nuaforward->GetYaxis()->FindBin(phi);
          Int_t nuavtz = fSettings.nuaforward->GetZaxis()->FindBin(zvertex);
          weight = weight*fSettings.nuaforward->GetBinContent(nuaeta,nuaphi,nuavtz);
        }
      }

    if (weight == 0) continue; 
    for (Int_t n = 0; n <= 4; n++) {

      if (useFMD && fSettings.sec_corr){

        if ((fSettings.ref_mode & fSettings.kFMDref)){ //doRefFlow && 
          if (!fSettings.use_primaries_fwd && n>=2 && n<=4) {
            Int_t seceta = fSettings.seccorr_fwd->GetXaxis()->FindBin(eta);
            Int_t secvtz = fSettings.seccorr_fwd->GetYaxis()->FindBin(zvertex);
            Int_t secn = n-1;//fSettings.seccorr_fwd->GetZaxis()->FindBin(n-2);
            weight = weight*fSettings.seccorr_fwd->GetBinContent(seceta,secvtz,secn);
          }
        }
      }

      if (weight == 0) continue;

      for (Int_t p = 1; p <= 4; p++) {
        Double_t realPart = TMath::Power(weight, p)*TMath::Cos(n*phi);
        Double_t imPart =   TMath::Power(weight, p)*TMath::Sin(n*phi);

        Double_t re[4] = {0.5, Double_t(n), Double_t(p), difEta};
        Double_t im[4] = {-0.5, Double_t(n), Double_t(p), difEta};

        if (doDiffFlow){
          fpvector->Fill(re, realPart);
          fpvector->Fill(im, imPart);

          //if (!(fSettings.etagap) && doRefFlow){
          if (!useFMD){
            fqvector->Fill(re, realPart);
            fqvector->Fill(im, imPart);
             // if (weight > 1.00001 ){
             //   fAutoDiff.Fill(eta, TMath::Gamma(weight+1)/TMath::Gamma(weight-1));
             // }
          }
        }

        if (doRefFlow){
          if ((fSettings.etagap) && TMath::Abs(eta)<=fSettings.gap) continue;
          if (fSettings.etagap && TMath::Abs(eta)>3.0) continue;
          if (fSettings.ref_mode & fSettings.kFMDref) {
            if (TMath::Abs(eta) < 2.0) continue;
          }
          Double_t req[4] = {0.5, static_cast<Double_t>(n), static_cast<Double_t>(p), refEta};
          Double_t imq[4] = {-0.5, static_cast<Double_t>(n), static_cast<Double_t>(p), refEta};
          fQvector->Fill(req, realPart);
          fQvector->Fill(imq, imPart);
           // if (weight > 1.00001){
           //   fAutoRef.Fill(eta, TMath::Gamma(weight+1)/TMath::Gamma(weight-1));
           // }
        }
      } // end p loop
    } // End of n loop
  } // End of phi loop
} // end of eta

  return;
}



void AliForwardGenericFramework::saveEvent(TList* outputList, double cent, double zvertex,UInt_t r, Int_t ptn){
  TList* analysisList = static_cast<TList*>(outputList->FindObject("cumulants"));

  THnD* cumuRef = static_cast<THnD*>(analysisList->At(0));
  THnD* cumuDiff = static_cast<THnD*>(analysisList->At(1));

  // For each n we loop over the hists
  Double_t noSamples = static_cast<Double_t>(r);

  for (Int_t n = 2; n <= 4; n++) {
    Int_t prevRefEtaBin = kTRUE;

    for (Int_t etaBin = 1; etaBin <= fpvector->GetAxis(3)->GetNbins(); etaBin++) {
      Double_t eta = fpvector->GetAxis(3)->GetBinCenter(etaBin);

      Int_t refEtaBinA = fQvector->GetAxis(3)->FindBin(eta);
      Int_t refEtaBinB = refEtaBinA;
      Int_t etaBinB = etaBin;

      if ((fSettings.etagap)) {
        refEtaBinB = fQvector->GetAxis(3)->FindBin(-eta);
        etaBinB = fpvector->GetAxis(3)->FindBin(-eta);
      }

      Double_t refEtaA = fQvector->GetAxis(3)->GetBinCenter(refEtaBinA);


      // index to get sum of weights
      Int_t index1[4] = {2, 1, 1, refEtaBinB};
      if (fQvector->GetBinContent(index1) > 0){
        // REFERENCE FLOW --------------------------------------------------------------------------------
        if (prevRefEtaBin){ // only used once

          // if (!(fSettings.etagap)){
          //   Double_t z[5] = {noSamples, zvertex, refEtaA, cent, Double_t(fSettings.kW2Two)};
          //   fQcorrfactor->Fill(z, fAutoRef.GetBinContent(etaBin));
          // }

          // two-particle cumulant
          double two = Two(n, -n, refEtaBinA, refEtaBinB).Re();
          double dn2 = Two(0,0, refEtaBinA, refEtaBinB).Re();

          Double_t x[7] = {Double_t(n-2),Double_t(ptn),noSamples, zvertex, refEtaA, cent, Double_t(fSettings.kW2TwoA)};//kW4FourA

          //x[4] = Double_t(fSettings.kW2Two);//kW2TwoA

          cumuRef->Fill(x, two);
          x[6] = Double_t(fSettings.kW2A);//kW2A
          cumuRef->Fill(x, dn2);

          // four-particle cumulant
          double four = Four(n, n, -n, -n, refEtaBinA, refEtaBinB).Re();
          double dn4 = Four(0,0,0,0 , refEtaBinA, refEtaBinB).Re();

          x[6] = Double_t(fSettings.kW4FourA);//kW4FourA
          cumuRef->Fill(x, four);
          x[6] = Double_t(fSettings.kW4A);//kW4A
          cumuRef->Fill(x, dn4);

          prevRefEtaBin = kFALSE;
        }
        // DIFFERENTIAL FLOW -----------------------------------------------------------------------------
        // if (n == 2 && (!(fSettings.etagap))){
        //   Double_t k[5] = {noSamples, zvertex, eta, cent, Double_t(fSettings.kW2Two)};
        //   fpcorrfactor->Fill(k, fAutoDiff.GetBinContent(etaBin));
        // }

        // two-particle cumulant
        double twodiff = TwoDiff(n, -n, refEtaBinB, etaBin).Re();
        double dn2diff = TwoDiff(0,0, refEtaBinB, etaBin).Re();

        Double_t y[7] = {Double_t(n-2),Double_t(ptn),noSamples, zvertex, eta, cent, Double_t(fSettings.kW2TwoB)};//kW2TwoB

        cumuDiff->Fill(y, twodiff);
        y[6] = Double_t(fSettings.kW2B);//kW2B
        cumuDiff->Fill(y, dn2diff);

        // A side
        
        twodiff = TwoDiff(n, -n, refEtaBinA, etaBin).Re();
        dn2diff = TwoDiff(0,0, refEtaBinA, etaBin).Re();

        y[6] = Double_t(fSettings.kW2TwoA);
        cumuDiff->Fill(y, twodiff);
        y[6] = Double_t(fSettings.kW2A);
        cumuDiff->Fill(y, dn2diff);
        
        if (eta > 0.0) {
          // R_{n,n; 2} numerator
          double over = (TwoDiff(-2,2,refEtaBinB, etaBin)*TwoDiff(2,-2,refEtaBinA, etaBinB)).Re();
          y[6] = Double_t(fSettings.kWTwoTwoN);
          cumuDiff->Fill(y, over);

          // R_{n,n; 2} denominator
          double under = (TwoDiff(-2,2,refEtaBinB, etaBinB)*TwoDiff(2,-2,refEtaBinA, etaBin)).Re();
          y[6] = Double_t(fSettings.kWTwoTwoD);
          cumuDiff->Fill(y, under);
        }

        // four-particle cumulant
        double fourdiff = FourDiff(n, n, -n, -n, refEtaBinA, refEtaBinB, etaBin,etaBin).Re(); // A is same side
        double dn4diff = FourDiff(0,0,0,0, refEtaBinA, refEtaBinB, etaBin,etaBin).Re();

        y[6] = Double_t(fSettings.kW4FourB);//kW4FourB
        cumuDiff->Fill(y, fourdiff);
        y[6] = Double_t(fSettings.kW4B);//kW4B
        cumuDiff->Fill(y, dn4diff);

        // four-particle cumulant SC(4,2)
        
        double fourtwodiff = FourDiff(4, 2, -4, -2, refEtaBinA,refEtaBinB, etaBin,etaBin).Re();

        y[6] = Double_t(fSettings.kW4FourTwoB);//kW4FourTwoB
        cumuDiff->Fill(y, fourtwodiff);

        // four-particle cumulant SC(3,2)
        double threetwodiff = FourDiff(3, 2, -3, -2, refEtaBinA,refEtaBinB, etaBin,etaBin).Re();

        y[6] = Double_t(fSettings.kW4ThreeTwoB);//kW4ThreeTwoB
        cumuDiff->Fill(y, threetwodiff);
        


        // A side
        // four-particle cumulant
        
        fourdiff = FourDiff(n, n, -n, -n, refEtaBinB, refEtaBinA, etaBin,etaBin).Re();
        dn4diff = FourDiff(0,0,0,0, refEtaBinB,refEtaBinA,  etaBin,etaBin).Re();

        y[6] = Double_t(fSettings.kW4FourA);
        cumuDiff->Fill(y, fourdiff);
        y[6] = Double_t(fSettings.kW4A);
        cumuDiff->Fill(y, dn4diff);

        
        // four-particle cumulant SC(4,2)
        fourtwodiff = FourDiff(4, 2, -4, -2, refEtaBinB,refEtaBinA, etaBin,etaBin).Re();

        y[6] = Double_t(fSettings.kW4FourTwoA);
        cumuDiff->Fill(y, fourtwodiff);

        // four-particle cumulant SC(3,2)
        threetwodiff = FourDiff(3, 2, -3, -2, refEtaBinB,refEtaBinA, etaBin,etaBin).Re();

        y[6] = Double_t(fSettings.kW4ThreeTwoA);
        cumuDiff->Fill(y, threetwodiff);
        

      } // if w2 > 0
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
  if (eta1 == eta2) {
     formula = Q(n1,1,eta1)*Q(n2,1,eta1) - Q(n1+n2,2,eta1);
  }
  else{
     formula = Q(n1,1,eta1)*Q(n2,1,eta2);
  }
  return formula;
}



TComplex AliForwardGenericFramework::TwoDiff(Int_t n1, Int_t n2, Int_t refetabin, Int_t diffetabin)
{
  return p(n1,1, diffetabin)*Q(n2,1, refetabin);// - q(n1+n2,1, diffetabin);
}


// TComplex AliForwardGenericFramework::TwoTwoDiff1(Int_t n1, Int_t n2, Int_t ref1, Int_t ref2, Int_t diff1,Int_t diff2)
// {
//   // 
//   return Q(n2,1, ref1)*p(n1,1, diff2)*Q(n2,1, ref2)*p(n1,1, diff1);// - q(n1+n2,1, diffetabin);
// }

// TComplex AliForwardGenericFramework::TwoTwoDiff2(Int_t n1, Int_t n2, Int_t ref1, Int_t ref2, Int_t diff1,Int_t diff2)
// {
//   // 
//   return Q(n2,1, ref1)*p(n1,1, diff2)*Q(n2,1, ref2)*p(n1,1, diff1);// - q(n1+n2,1, diffetabin);
// }

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
  // if (refetabinPos == refetabinNeg){
  //   formula = p(n1,1,diffetabin)*Q(n2,1,refetabin)*Q(n3,1,refetabin)*Q(n4,1,refetabin)-q(n1+n2,2,qetabin)*Q(n3,1,refetabin)*Q(n4,1,refetabin)-Q(n2,1,refetabin)*q(n1+n3,2,qetabin)*Q(n4,1,refetabin)
  //           - p(n1,1,diffetabin)*Q(n2+n3,2,refetabin)*Q(n4,1,refetabin)+2.*q(n1+n2+n3,3,qetabin)*Q(n4,1,refetabin)-Q(n2,1,refetabin)*Q(n3,1,refetabin)*q(n1+n4,2,qetabin)
  //           + Q(n2+n3,2,refetabin)*q(n1+n4,2,qetabin)-p(n1,1,diffetabin)*Q(n3,1,refetabin)*Q(n2+n4,2,refetabin)+q(n1+n3,2,qetabin)*Q(n2+n4,2,refetabin)
  //           + 2.*Q(n3,1,refetabin)*q(n1+n2+n4,3,qetabin)-p(n1,1,diffetabin)*Q(n2,1,refetabin)*Q(n3+n4,2,refetabin)+q(n1+n2,2,qetabin)*Q(n3+n4,2,refetabin)
  //           + 2.*Q(n2,1,refetabin)*q(n1+n3+n4,3,qetabin)+2.*p(n1,1,diffetabin)*Q(n2+n3+n4,3,refetabin)-6.*q(n1+n2+n3+n4,4,qetabin);
  // }
  // else{
    formula = p(n1,1,diffetabin)*Q(n2,1,refetabinA)*Q(n3,1,refetabinB)*Q(n4,1,refetabinB)
            - q(n1+n2,2,qetabin)*Q(n3,1,refetabinB)*Q(n4,1,refetabinB)
            - p(n1,1,diffetabin)*Q(n2,1,refetabinA)*Q(n3+n4,2,refetabinB)
            + q(n1+n2,2,qetabin)*Q(n3+n4,2,refetabinB);
  // }
  return formula;
}


void AliForwardGenericFramework::reset() {
  fQvector->Reset();
  fpvector->Reset();
  fqvector->Reset();
}
