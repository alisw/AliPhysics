/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*   Origin: Alberica Toia, CERN, Alberica.Toia@cern.ch                   */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class to determine centrality percentiles from 1D distributions          // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TH1D.h>
#include <TString.h>
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TH2F.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <vector>
#include "AliCentralityGlauberFit.h"


ClassImp(AliCentralityGlauberFit)  
 
//______________________________________________________________________________


AliCentralityGlauberFit::AliCentralityGlauberFit() {
  // standard constructor
  // glauber file
  TFile *f = TFile::Open("hj-unquenched.root");
  glau_ntuple = (TNtuple*) f->Get("gnt");
  glau_ntuple->SetBranchAddress("npart",&npart);
  glau_ntuple->SetBranchAddress("ncoll",&ncoll);
  glau_ntuple->SetBranchAddress("b",&b);
}

//--------------------------------------------------------------------------------------------------

AliCentralityGlauberFit::~AliCentralityGlauberFit() {
  // destructor
}

//--------------------------------------------------------------------------------------------------

void AliCentralityGlauberFit::AddHisto(TString name) {
  histnames.push_back(name);
}

//--------------------------------------------------------------------------------------------------

void AliCentralityGlauberFit::SetOutputFile(TString filename) {
  outrootfilename = filename;
}

//--------------------------------------------------------------------------------------------------

void AliCentralityGlauberFit::SetGlauberParam(
					      Int_t   fNmu,
					      Float_t fmulow,
					      Float_t fmuhigh,
					      Int_t   fNk,
					      Float_t fklow,
					      Float_t fkhigh,
					      Int_t fNalpha,
					      Float_t falphalow,
					      Float_t falphahigh)

{
  Nmu=fNmu;
  mulow=fmulow;
  muhigh=fmuhigh;
  Nk=fNk;
  klow=fklow;
  khigh=fkhigh;
  Nalpha=fNalpha;
  alphalow=falphalow;
  alphahigh=falphahigh;
}

//--------------------------------------------------------------------------------------------------

Float_t AliCentralityGlauberFit::GetPercentileCrossSection() {
  return percentXsec;
}

//--------------------------------------------------------------------------------------------------

void AliCentralityGlauberFit::MakeFits(TString infilename) {
  TH1D *hDATA;
  TH1D *thistGlau;
  
  float efflow = 0.60;
  float effhigh = 1.00;

  TFile *outrootfile;
  
  // open inrootfile, outrootfile
  inrootfile  = new TFile(infilename);
  outrootfile = new TFile(outrootfilename,"RECREATE");
  
  // loop over all distribution names  
  vector<TString>::const_iterator hni;
  for(hni=histnames.begin(); hni!=histnames.end(); hni++) {
    hDATA  = (TH1D*) (inrootfile->Get(*hni)); 
    hDATA->Rebin(10);
    TH1D *hGLAU = new TH1D("hGLAU","hGLAU",hDATA->GetNbinsX(),0,hDATA->GetNbinsX()*hDATA->GetBinWidth(1));
    int chi2min = 9999;
    Float_t alpha_min=-1;
    Float_t mu_min=-1;
    Float_t k_min=-1;
    Float_t eff_min=-1;
    Float_t alpha, mu, k, eff, chi2;   

    for (int nalpha=0;nalpha<Nalpha;nalpha++) {
      alpha = alphalow + ((float) nalpha) * 0.05;
      mu=0.0;
      for (int nmu=0; nmu<Nmu; nmu++) {
	//	mu = (mulow*(1-nalpha*0.1*mulow) )  + ((float) nmu ) * (muhigh  - mulow ) / Nmu;
	mu = mulow   + ((float) nmu ) * (muhigh  - mulow ) / Nmu;
	
	for (int nk=0; nk<Nk; nk++) {
	  k = klow  + ((float) nk ) * (khigh  - klow ) / Nk;
	  
	  for (int neff=0; neff<20; neff++) {
	    eff = efflow + ((float) neff) * (effhigh - efflow) / 20.0;
	    
	    thistGlau = GlauberHisto(mu,k,eff,alpha,hDATA,kFALSE);
	    chi2 = CalculateChi2(hDATA,thistGlau,eff);
	    if ( chi2 < chi2min ) {
	      chi2min = chi2;
	      alpha_min=alpha;
	      mu_min=mu;
	      k_min=k;
	      eff_min=eff;
	    }
	  }
	}
      }
    }

    thistGlau = GlauberHisto(mu_min,k_min,eff_min,alpha_min,hDATA,kTRUE);
    hGLAU = (TH1D *) thistGlau->Clone("hGLAU");
    hGLAU->SetName( ((TString)hDATA->GetName()).Append(Form("_GLAU_%d_%d_%d_%d",(int)(100*mu_min),(int)(100*k_min),(int)(100*alpha_min),(int)(100*eff_min))));
    hGLAU->SetTitle( ((TString)hDATA->GetName()).Append(Form("_GLAU_%.2f_%.2f_%.2f_%.2f",mu_min,k_min,alpha_min,eff_min)));

    heffi = GetTriggerEfficiencyFunction(hDATA, hGLAU);
    SaveHisto(hDATA,hGLAU,heffi,outrootfile);

  }
  // close inrootfile, outrootfile
  inrootfile->Close();
  outrootfile->Close();
  
}

//--------------------------------------------------------------------------------------------------

TH1D * AliCentralityGlauberFit::GlauberHisto(Float_t mu, Float_t k, Float_t eff, Float_t alpha, TH1D *hDATA, Bool_t save) {
    
  TH1D *hSample = NBDhist(mu,k);
  TH1D *h1 = new TH1D(Form("fit_%.2f_%.2f_%.2f_%.2f",mu,k,eff,alpha),"",hDATA->GetNbinsX(),0,hDATA->GetNbinsX()*hDATA->GetBinWidth(1));
  TFile *outFile = NULL;
  TNtuple *ntuple = NULL;
 
  if (save) {
    outFile = TFile::Open("pippo.root", "RECREATE");
    ntuple = new TNtuple("gnt", "Glauber ntuple", "npart:ncoll:b:ntot");
  } 
  
  Int_t nents = glau_ntuple->GetEntries();
  for (Int_t i=0;i<nents;++i) {
    glau_ntuple->GetEntry(i);
    Int_t n = TMath::Nint(TMath::Power(npart,alpha));
    //Int_t n = TMath::Nint(TMath::Power(ncoll,alpha));
    Int_t ntot=0;
    for(Int_t j = 0; j<n; ++j) 
      ntot+=hSample->GetRandom();
    h1->Fill(ntot);
    
    if (save) 
      ntuple->Fill(npart,ncoll,b,ntot);
   
  }
  
  if (save) {
    ntuple->Write();
    outFile->Close();
  }
  
  
  return h1;
  
}

//--------------------------------------------------------------------------------------------------

Float_t AliCentralityGlauberFit::CalculateChi2(TH1D *hDATA, TH1D *thistGlau, Float_t eff) {
  // note that for different values of neff the mc distribution is identical, just re-scaled below
  // normalize the two histogram
  // scale MC up but leave data alone - that preserves error bars !!!!
  
  float mcintegral = thistGlau->Integral();
  for (int i=1; i<=thistGlau->GetNbinsX(); i++) {
    float scale = (hDATA->Integral()/mcintegral) / ((float) eff);
    float temp = scale * thistGlau->GetBinContent(i);
    thistGlau->SetBinContent(i,temp);
  }
  
  // calculate the chi2 between MC and real data over some range ????
  int lowchibin =   10;
  int highchibin =  80;
  float chi2 = 0.0;
  for (int i=lowchibin; i<=highchibin; i++) {
    if (hDATA->GetBinContent(i) < 1.0) continue;
    float diff = pow( (thistGlau->GetBinContent(i) - hDATA->GetBinContent(i)) , 2);
    diff = diff / pow(hDATA->GetBinError(i) , 2);
    chi2 += diff;
  }
  chi2 = chi2 / (highchibin - lowchibin + 1);
  return chi2;
}

//--------------------------------------------------------------------------------------------------

TH1D * AliCentralityGlauberFit::GetTriggerEfficiencyFunction(TH1D *hist1, TH1D *hist2) {
  heffi = (TH1D*)hist1->Clone("heffi");
  heffi->Divide(hist2);
  return heffi;
}

//--------------------------------------------------------------------------------------------------

Float_t AliCentralityGlauberFit::GetTriggerEfficiencyIntegral(TH1D *hist1, TH1D *hist2) {
  effi = hist1->Integral()/hist2->Integral();
  return effi;
}

//--------------------------------------------------------------------------------------------------

void AliCentralityGlauberFit::SaveHisto(TH1D *hist1, TH1D *hist2, TH1D *heffi, TFile *outrootfile) {
  outrootfile->cd();
  hist1->Write();
  hist2->Write();
  heffi->Write();
}

//--------------------------------------------------------------------------------------------------

Double_t AliCentralityGlauberFit::NBD(Int_t n, Double_t mu, Double_t k)
{
  Double_t mudk = mu/k;
  Double_t ret = TMath::Gamma(n+k) / TMath::Gamma(k) / TMath::Factorial(n) *
                 TMath::Power(mudk,n) / TMath::Power(1+mudk,n+k);
  return ret;
}

//--------------------------------------------------------------------------------------------------

TH1D *AliCentralityGlauberFit::NBDhist(Double_t mu, Double_t k)
{
  TH1D *h = new TH1D("htemp","",100,0,100);
  h->SetName(Form("nbd_%f_%f",mu,k));
  h->SetDirectory(0);
  for (Int_t i=0;i<100;++i) {
    Double_t val = NBD(i,mu,k);
    h->Fill(i,val);
  }
  return h;
}


