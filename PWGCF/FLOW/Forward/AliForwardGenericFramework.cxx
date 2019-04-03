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
  Int_t rbins[4] = {2, 6, 5, 2} ; // kind (real or imaginary), n, p, eta
  Int_t dimensions = 4;
  if (fSettings.etagap) {
     rbins[3] = 2 ; // two bins in eta for gap, one for standard
  }

  Double_t xmin[4] = {-1.0, -0.5, 0.5, -6}; // kind (real or imaginary), n, p, eta
  Double_t xmax[4] = { 1,   5.5, 4.5,  6}; // kind (real or imaginary), n, p, eta SKAL VAERE -6 - 6

  fQvector = new THnD("Qvector", "Qvector", dimensions, rbins, xmin, xmax);

  Int_t dbins[4] = {2, 6, 4, fSettings.fNDiffEtaBins} ; // kind (real or imaginary), n, p, eta
  xmin[3] = -4.0;

  fpvector = new THnD("pvector", "pvector", dimensions, dbins, xmin, xmax);
  fqvector = new THnD("qvector", "qvector", dimensions, dbins, xmin, xmax);

  fAutoRef = TH1F("fAutoRef","fAutoRef", fSettings.fNRefEtaBins, fSettings.fEtaLowEdge, fSettings.fEtaUpEdge);
  fAutoDiff = TH1F("fAutoDiff","fAutoDiff", fSettings.fNDiffEtaBins, fSettings.fEtaLowEdge, fSettings.fEtaUpEdge);
  fAutoRef.SetDirectory(0);
  fAutoDiff.SetDirectory(0);
}


//_____________________________________________________________________
void AliForwardGenericFramework::CumulantsAccumulate(TH2D& dNdetadphi, TList* outputList, double cent, double zvertex, TString detType, Bool_t doRefFlow, Bool_t doDiffFlow)
{
  //TList* eventList = static_cast<TList*>(outputList->FindObject("EventInfo"));
  //TH1D* fFMDHits = static_cast<TH1D*>(eventList->FindObject("FMDHits"));

  for (Int_t etaBin = 1; etaBin <= dNdetadphi.GetNbinsX(); etaBin++) {

    Double_t eta = dNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    Double_t difEtaBin = fpvector->GetAxis(3)->FindBin(eta);
    Double_t difEta = fpvector->GetAxis(3)->GetBinCenter(difEtaBin);

    Double_t refEtaBin = fQvector->GetAxis(3)->FindBin(eta);
    Double_t refEta = fQvector->GetAxis(3)->GetBinCenter(refEtaBin);

    for (Int_t phiBin = 1; phiBin <= dNdetadphi.GetNbinsY(); phiBin++) {
      Double_t phi = dNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      Double_t weight = dNdetadphi.GetBinContent(etaBin, phiBin);

      if (!fSettings.use_primaries_fwd && !fSettings.esd){
        if (dNdetadphi.GetBinContent(etaBin, 0) == 0 && detType == "forward") break;
      }



      if (fSettings.doNUA){
        // holes in the FMD
        if ((fSettings.nua_mode & fSettings.kFill) && detType == "forward"){
          if (etaBin >= 125 && etaBin <=137){
            if (phiBin == 17 || phiBin == 18) weight = 1.;
          }
          if (etaBin >= 168 && etaBin <=185){
            if (phiBin == 14) weight = 1.;
          }
        }

        if ((fSettings.nua_mode & fSettings.kInterpolate) && detType == "forward") weight = AliForwardNUATask::InterpolateWeight(dNdetadphi,phiBin,etaBin,weight);

        if (doRefFlow){
          if (detType == "central" && !fSettings.use_primaries_cen) {
            Double_t nuaeta = fSettings.nuacentral_ref->GetXaxis()->FindBin(eta);
            Double_t nuaphi = fSettings.nuacentral_ref->GetYaxis()->FindBin(phi);
            Double_t nuavtz = fSettings.nuacentral_ref->GetZaxis()->FindBin(zvertex);
            weight = weight*fSettings.nuacentral_ref->GetBinContent(nuaeta,nuaphi,nuavtz);
          }

          if (detType == "forward" && !fSettings.use_primaries_fwd) {
            Double_t nuaeta = fSettings.nuaforward_ref->GetXaxis()->FindBin(eta);
            Double_t nuaphi = fSettings.nuaforward_ref->GetYaxis()->FindBin(phi);
            Double_t nuavtz = fSettings.nuaforward_ref->GetZaxis()->FindBin(zvertex);
            weight = weight*fSettings.nuaforward_ref->GetBinContent(nuaeta,nuaphi,nuavtz);
          }
        }
        else{
          if (detType == "central" && !fSettings.use_primaries_cen) {
            Double_t nuaeta = fSettings.nuacentral->GetXaxis()->FindBin(eta);
            Double_t nuaphi = fSettings.nuacentral->GetYaxis()->FindBin(phi);
            Double_t nuavtz = fSettings.nuacentral->GetZaxis()->FindBin(zvertex);
            weight = weight*fSettings.nuacentral->GetBinContent(nuaeta,nuaphi,nuavtz);
          }

          if (detType == "forward" && !fSettings.use_primaries_fwd) {
            Double_t nuaeta = fSettings.nuaforward->GetXaxis()->FindBin(eta);
            Double_t nuaphi = fSettings.nuaforward->GetYaxis()->FindBin(phi);
            Double_t nuavtz = fSettings.nuaforward->GetZaxis()->FindBin(zvertex);
            weight = weight*fSettings.nuaforward->GetBinContent(nuaeta,nuaphi,nuavtz);
          }
        }
      }

    if (weight == 0) continue; 
    for (Int_t n = 0; n <= 5; n++) {

      /*
      if (doRefFlow && (fSettings.ref_mode & fSettings.kFMDref)){
        if (!fSettings.use_primaries_fwd && n>=2 && n<=4) {
          Double_t seceta = fSettings.seccorr_fwd->GetXaxis()->FindBin(eta);
          Double_t secvtz = fSettings.seccorr_fwd->GetYaxis()->FindBin(zvertex);
          Double_t secn = fSettings.seccorr_fwd->GetZaxis()->FindBin(n-2);
          weight = weight*fSettings.seccorr_fwd->GetBinContent(seceta,secvtz,secn);
        }
      }
      if (doRefFlow && (fSettings.ref_mode & fSettings.kITSref)){
        if (!fSettings.use_primaries_cen && n>=2 && n<=4) {
          Double_t seceta = fSettings.seccorr_cen->GetXaxis()->FindBin(eta);
          Double_t secvtz = fSettings.seccorr_cen->GetYaxis()->FindBin(zvertex);
          Double_t secn = fSettings.seccorr_cen->GetZaxis()->FindBin(n-2);
          weight = weight*fSettings.seccorr_cen->GetBinContent(seceta,secvtz,secn);
        }
      }
      */


      for (Int_t p = 1; p <= 4; p++) {
        Double_t realPart = TMath::Power(weight, p)*TMath::Cos(n*phi);
        Double_t imPart =   TMath::Power(weight, p)*TMath::Sin(n*phi);

        Double_t re[4] = {0.5, Double_t(n), Double_t(p), difEta};
        Double_t im[4] = {-0.5, Double_t(n), Double_t(p), difEta};

        if (doDiffFlow){
          fpvector->Fill(re, realPart);
          fpvector->Fill(im, imPart);

          if (!(fSettings.etagap) && doRefFlow){
            fqvector->Fill(re, realPart);
            fqvector->Fill(im, imPart);
            if (weight > 1.00001 ){
              fAutoDiff.Fill(eta, TMath::Gamma(weight+1)/TMath::Gamma(weight-1));
            }
          }
        }

        if (doRefFlow){
          if ((fSettings.etagap) && fabs(eta)<=fSettings.gap) continue;

          Double_t req[4] = {0.5, static_cast<Double_t>(n), static_cast<Double_t>(p), refEta};
          Double_t imq[4] = {-0.5, static_cast<Double_t>(n), static_cast<Double_t>(p), refEta};
          fQvector->Fill(req, realPart);
          fQvector->Fill(imq, imPart);
          if (!(fSettings.etagap) && weight > 1.00001){
            fAutoRef.Fill(eta, TMath::Gamma(weight+1)/TMath::Gamma(weight-1));
          }
        }
      } // end p loop
    } // End of n loop
  } // End of phi loop
} // end of eta

  return;
}



void AliForwardGenericFramework::saveEvent(TList* outputList, double cent, double zvertex,UInt_t r, Int_t ptn){
  TList* analysisList = static_cast<TList*>(outputList->FindObject("Analysis"));
  TList* refList = static_cast<TList*>(analysisList->FindObject("Reference"));
  TList* autoList = static_cast<TList*>(analysisList->FindObject("AutoCorrection"));
  THnD*  fQcorrfactor = static_cast<THnD*>(autoList->FindObject("fQcorrfactor"));
  THnD*  fpcorrfactor = static_cast<THnD*>(autoList->FindObject("fpcorrfactor"));
  TList* difList = static_cast<TList*>(analysisList->FindObject("Differential"));

  THnD* cumuRef = 0;
  THnD* cumuDiff = 0;

  // For each n we loop over the hists
  Double_t noSamples = static_cast<Double_t>(r);

  for (Int_t n = 2; n <= 5; n++) {
    Int_t prevRefEtaBin = kTRUE;

    cumuRef = static_cast<THnD*>(refList->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)));
    cumuDiff = static_cast<THnD*>(difList->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)));

    for (Int_t etaBin = 1; etaBin <= fpvector->GetAxis(3)->GetNbins(); etaBin++) {
      Double_t eta = fpvector->GetAxis(3)->GetBinCenter(etaBin);

      Double_t refEtaBinA = fQvector->GetAxis(3)->FindBin(eta);
      Double_t refEtaBinB = refEtaBinA;

      if ((fSettings.etagap)) {
        refEtaBinB = fQvector->GetAxis(3)->FindBin(-eta);
      }

      Double_t refEtaA = fQvector->GetAxis(3)->GetBinCenter(refEtaBinA);

      Int_t n_0 = 0;
      Int_t p_1 = 1;
      // index to get sum of weights
      Int_t index1[4] = {fQvector->GetAxis(0)->FindBin(0.5), fQvector->GetAxis(1)->FindBin(n_0), fQvector->GetAxis(2)->FindBin(p_1), static_cast<Int_t>(refEtaBinB)};
      if (fQvector->GetBinContent(index1) > 0){
        // REFERENCE FLOW --------------------------------------------------------------------------------
        if (prevRefEtaBin){ // only used once

          if (!(fSettings.etagap)){
            Double_t z[5] = {noSamples, zvertex, refEtaA, cent, Double_t(fSettings.kW2Two)};

            fQcorrfactor->Fill(z, fAutoRef.GetBinContent(etaBin));
          }
          // two-particle cumulant
          double two = Two(n, -n, refEtaBinA, refEtaBinB).Re();
          double dn2 = Two(0,0, refEtaBinA, refEtaBinB).Re();

          Double_t x[5] = {noSamples, zvertex, refEtaA, cent, Double_t(fSettings.kW4Four)};
          x[4] = Double_t(fSettings.kW2Two);

          cumuRef->Fill(x, two);
          x[4] = Double_t(fSettings.kW2);
          cumuRef->Fill(x, dn2);

          // four-particle cumulant
          double four = Four(n, n, -n, -n, refEtaBinA, refEtaBinB).Re();
          double dn4 = Four(0,0,0,0 , refEtaBinA, refEtaBinB).Re();

          x[4] = Double_t(fSettings.kW4Four);
          cumuRef->Fill(x, four);
          x[4] = Double_t(fSettings.kW4);
          cumuRef->Fill(x, dn4);

          prevRefEtaBin = kFALSE;
        }
        // DIFFERENTIAL FLOW -----------------------------------------------------------------------------
        if (n == 2 && (!(fSettings.etagap))){
          Double_t k[5] = {noSamples, zvertex, eta, cent, Double_t(fSettings.kW2Two)};

          fpcorrfactor->Fill(k, fAutoDiff.GetBinContent(etaBin));
        }
        // two-particle cumulant
        double twodiff = TwoDiff(n, -n, refEtaBinB, etaBin).Re();
        double dn2diff = TwoDiff(0,0, refEtaBinB, etaBin).Re();

        Double_t y[5] = {noSamples, zvertex, eta, cent, Double_t(fSettings.kW2Two)};

        cumuDiff->Fill(y, twodiff);
        y[4] = Double_t(fSettings.kW2);
        cumuDiff->Fill(y, dn2diff);

        // four-particle cumulant
        double fourdiff = FourDiff(n, n, -n, -n, refEtaBinB, etaBin,etaBin).Re();
        double dn4diff = FourDiff(0,0,0,0, refEtaBinB, etaBin,etaBin).Re();

        y[4] = Double_t(fSettings.kW4Four);
        cumuDiff->Fill(y, fourdiff);
        y[4] = Double_t(fSettings.kW4);
        cumuDiff->Fill(y, dn4diff);
      } // if w2 > 0
    } //eta
  } // moment

  return;
}



TComplex AliForwardGenericFramework::Q(int n, int p, int etabin)
{
  double sign = (n < 0) ? -1 : 1;

  Int_t imindex[4] = {fQvector->GetAxis(0)->FindBin(-0.5), fQvector->GetAxis(1)->FindBin(fabs(n)), fQvector->GetAxis(2)->FindBin(p), etabin};
  Int_t reindex[4] = {fQvector->GetAxis(0)->FindBin(0.5), fQvector->GetAxis(1)->FindBin(fabs(n)), fQvector->GetAxis(2)->FindBin(p), etabin};

  TComplex Q_temp = TComplex(fQvector->GetBinContent(reindex),sign*fQvector->GetBinContent(imindex));

  return Q_temp;
}

TComplex AliForwardGenericFramework::p(int n, int p, int etabin)
{
  double sign = (n > 0) ? 1 : ((n < 0) ? -1 : 1);
  Int_t imindex[4] = {fpvector->GetAxis(0)->FindBin(-0.5), fpvector->GetAxis(1)->FindBin(fabs(n)), fpvector->GetAxis(2)->FindBin(p), etabin};
  Int_t reindex[4] = {fpvector->GetAxis(0)->FindBin(0.5), fpvector->GetAxis(1)->FindBin(fabs(n)), fpvector->GetAxis(2)->FindBin(p), etabin};

  TComplex p_temp = TComplex(fpvector->GetBinContent(reindex),sign*fpvector->GetBinContent(imindex));

  return p_temp;
}


TComplex AliForwardGenericFramework::q(int n, int p, int etabin)
{
  double sign = (n > 0) ? 1 : ((n < 0) ? -1 : 1);
  Int_t imindex[4] = {fqvector->GetAxis(0)->FindBin(-0.5), fqvector->GetAxis(1)->FindBin(fabs(n)), fqvector->GetAxis(2)->FindBin(p), etabin};
  Int_t reindex[4] = {fqvector->GetAxis(0)->FindBin(0.5), fqvector->GetAxis(1)->FindBin(fabs(n)), fqvector->GetAxis(2)->FindBin(p), etabin};

  TComplex q_temp = TComplex(fqvector->GetBinContent(reindex),sign*fqvector->GetBinContent(imindex));

  return q_temp;
}

TComplex AliForwardGenericFramework::Two(int n1, int n2, int eta1, int eta2)
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

TComplex AliForwardGenericFramework::TwoDiff(int n1, int n2, int refetabin, int diffetabin)
{
  TComplex formula =0;

  formula = p(n1,1, diffetabin)*Q(n2,1, refetabin) - q(n1+n2,1, diffetabin);
  return formula;
}

TComplex AliForwardGenericFramework::Four(int n1, int n2, int n3, int n4,int eta1, int eta2)
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

TComplex AliForwardGenericFramework::FourDiff(int n1, int n2, int n3, int n4, int refetabin, int diffetabin,int qetabin)
{

  TComplex formula = p(n1,1,diffetabin)*Q(n2,1,refetabin)*Q(n3,1,refetabin)*Q(n4,1,refetabin)-q(n1+n2,2,qetabin)*Q(n3,1,refetabin)*Q(n4,1,refetabin)-Q(n2,1,refetabin)*q(n1+n3,2,qetabin)*Q(n4,1,refetabin)
                    - p(n1,1,diffetabin)*Q(n2+n3,2,refetabin)*Q(n4,1,refetabin)+2.*q(n1+n2+n3,3,qetabin)*Q(n4,1,refetabin)-Q(n2,1,refetabin)*Q(n3,1,refetabin)*q(n1+n4,2,qetabin)
                    + Q(n2+n3,2,refetabin)*q(n1+n4,2,qetabin)-p(n1,1,diffetabin)*Q(n3,1,refetabin)*Q(n2+n4,2,refetabin)+q(n1+n3,2,qetabin)*Q(n2+n4,2,refetabin)
                    + 2.*Q(n3,1,refetabin)*q(n1+n2+n4,3,qetabin)-p(n1,1,diffetabin)*Q(n2,1,refetabin)*Q(n3+n4,2,refetabin)+q(n1+n2,2,qetabin)*Q(n3+n4,2,refetabin)
                    + 2.*Q(n2,1,refetabin)*q(n1+n3+n4,3,qetabin)+2.*p(n1,1,diffetabin)*Q(n2+n3+n4,3,refetabin)-6.*q(n1+n2+n3+n4,4,qetabin);
  return formula;
}


void AliForwardGenericFramework::reset() {
  delete fQvector;
  delete fpvector;
  delete fqvector;
}
