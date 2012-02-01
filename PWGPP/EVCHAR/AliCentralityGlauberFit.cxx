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
#include <TMinuit.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include "AliCentralityGlauberFit.h"
#include "AliLog.h"

ClassImp(AliCentralityGlauberFit)  
 
//______________________________________________________________________________
AliCentralityGlauberFit::AliCentralityGlauberFit(const char *filename) : 
  fNmu(0),    
  fMulow(0),  
  fMuhigh(0), 
  fNk(0),     
  fKlow(0),   
  fKhigh(0),  
  fNalpha(0), 
  fAlphalow(0),   
  fAlphahigh(0),  
  fRebinFactor(0),
  fScalemin(0),    
  fMultmin(0),    
  fMultmax(0),    
  fGlauntuple(0), 
  fNpart(0),       
  fNcoll(0),       
  fB(0),           
  fTaa(0),         
  fEffi(0),        
  fhEffi(0),      
  fTempHist(0),   
  fGlauberHist(0), 
  fFastFit(0),      
  fAncestor(2),      
  fNBD(0),         
  fUseChi2(kTRUE),      
  fUseAverage(kFALSE),
  fhAncestor(0),
  fNevents(100000),
  fNtrials(100),
  fInrootfilename(0),
  fInntuplename(0),
  fOutrootfilename(0),
  fOutntuplename(0),
  fAncfilename("ancestor_hists.root"),
  fHistnames()
{
  // Standard constructor.
  TFile *f = 0;
  if (filename) {  // glauber file
    f = TFile::Open(filename);
    fGlauntuple = (TNtuple*) f->Get("nt_Pb_Pb");
    fGlauntuple->SetBranchAddress("Npart",&fNpart);
    fGlauntuple->SetBranchAddress("Ncoll",&fNcoll);
    fGlauntuple->SetBranchAddress("B",&fB);
    fGlauntuple->SetBranchAddress("tAA",&fTaa);
  }

  fNBD = new TF1 ("fNBD", AliCentralityGlauberFit::NBDFunc, 0, 100,2);
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::SetRangeToFit(Double_t fmultmin, Double_t fmultmax)
{
  // Set fit range.

  fMultmin=fmultmin;
  fMultmax=fmultmax;
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::SetRangeToScale(Double_t fscalemin)
{
  // Set range where to scale simulated histo to data

  fScalemin=fscalemin;
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::SetGlauberParam(
					      Int_t   Nmu,
					      Double_t mulow,
					      Double_t muhigh,
					      Int_t   Nk,
					      Double_t klow,
					      Double_t khigh,
					      Int_t   Nalpha,
					      Double_t alphalow,
					      Double_t alphahigh)
{
  // Set Glauber parameters.

  fNmu=Nmu;
  fMulow=mulow;
  fMuhigh=muhigh;
  fNk=Nk;
  fKlow=klow;
  fKhigh=khigh;
  fNalpha=Nalpha;
  fAlphalow=alphalow;
  fAlphahigh=alphahigh;
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::MakeFits() 
{
  // Make fits.

  TH1F *hDATA;
  TH1F *thistGlau; 
  TFile *inrootfile;
  TFile *outrootfile;  
  //FILE* fTxt = fopen ("parameters.txt","w");
  
  // open inrootfile, outrootfile
  std::cout << "input file "  << fInrootfilename  << std::endl;
  std::cout << "output file " << fOutrootfilename << std::endl;
  inrootfile  = new TFile(fInrootfilename,"OPEN");
  outrootfile = new TFile(fOutrootfilename,"RECREATE");
  
  // loop over all distribution names  
  std::vector<TString>::const_iterator hni;
  for(hni=fHistnames.begin(); hni!=fHistnames.end(); hni++) {
    hDATA  = (TH1F*) (inrootfile->Get(*hni)); 
    if (!hDATA) {
      TList *list  = (TList*) (inrootfile->Get("CentralityStat")); 
      hDATA  = (TH1F*) (list->FindObject(*hni));
    } 
    hDATA->Rebin(fRebinFactor);
    //TH1F *hGLAU = new TH1F("hGLAU","hGLAU",hDATA->GetNbinsX(),0,hDATA->GetNbinsX()*hDATA->GetBinWidth(1));
    Double_t chi2min = 9999999.0;
    Double_t alpha_min=-1;
    Double_t mu_min=-1;
    Double_t k_min=-1;
    Double_t alpha, mu, k, chi2;   

    for (Int_t nalpha=0;nalpha<fNalpha;nalpha++) {
      alpha = fAlphalow   + ((Double_t) nalpha ) * (fAlphahigh  - fAlphalow ) / fNalpha;

      mu=0.0;
      for (Int_t nmu=0; nmu<fNmu; nmu++) {
	mu = fMulow   + ((Double_t) nmu ) * (fMuhigh  - fMulow ) / fNmu;

	for (Int_t nk=0; nk<fNk; nk++) {
	  k = fKlow  + ((Double_t) nk ) * (fKhigh  - fKlow ) / fNk;
	  	    
	  thistGlau = GlauberHisto(mu,k,alpha,hDATA,kFALSE);
	  chi2 = CalculateChi2(hDATA,thistGlau);
	  //fprintf(fTxt, "%3.3f %3.3f %3.3f %3.3f\n", (float) mu, (float) k, (float) alpha, chi2);
	  if ( chi2 < chi2min ) {
	    chi2min = chi2;
	    alpha_min=alpha;
	    mu_min=mu;
	    k_min=k;
	  }
	}
      }
    }
    
    thistGlau = GlauberHisto(mu_min,k_min,alpha_min,hDATA,kTRUE);
    
    TH1F * hGLAU = 0x0;
    hGLAU = (TH1F *) thistGlau->Clone("hGLAU");
    hGLAU->SetName( ((TString)hDATA->GetName()).Append(Form("_GLAU")));
    hGLAU->SetTitle( ((TString)hDATA->GetName()).Append(Form("_GLAU_%.3f_%.3f_%.3f",
                                                             mu_min,k_min,alpha_min)));

    Double_t mcintegral = hGLAU->Integral(hGLAU->FindBin(fScalemin),hGLAU->GetNbinsX());
    Double_t scale = (hDATA->Integral(hDATA->FindBin(fScalemin),hDATA->GetNbinsX())/mcintegral);
    hGLAU->Scale(scale);

    fhEffi = GetTriggerEfficiencyFunction(hDATA, hGLAU);
    SaveHisto(hDATA,hGLAU,fhEffi,outrootfile);
    //fclose (fTxt);

    std::cout << "chi2 min is " << chi2min << std::endl;
    std::cout << "fitted " << hGLAU->Integral(hGLAU->FindBin(fMultmin),
                                              hGLAU->FindBin(fMultmax))/hGLAU->Integral() 
              << " of the total cross section" << std::endl; 
    fTempHist=hDATA;
    fGlauberHist=hGLAU;
  }

  // close inrootfile, outrootfile
  inrootfile->Close();
  outrootfile->Close();
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::MakeFitsMinuitNBD(Double_t alpha, Double_t mu, Double_t k) 
{
  // Make fits using Minuit.

  if (alpha<0) 
    alpha = fAlphalow;
  if (mu<0)
    mu = fMulow;
  if (k<0)
    k = fKlow;
  printf("Calling Minuit with starting values: %f %f %f\n", alpha, mu, k);

  TH1F *hDATA;
  TH1F *thistGlau; 
  TFile *inrootfile;
  TFile *outrootfile;
  
  // open inrootfile, outrootfile
  std::cout << "input file "  << fInrootfilename  << std::endl;
  std::cout << "output file " << fOutrootfilename << std::endl;
  inrootfile  = new TFile(fInrootfilename,"OPEN");
  outrootfile = new TFile(fOutrootfilename,"RECREATE");

  // loop over all distribution names  
  std::vector<TString>::const_iterator hni;
  for(hni=fHistnames.begin(); hni!=fHistnames.end(); hni++) {
    hDATA  = (TH1F*) (inrootfile->Get(*hni)); 
    hDATA->Rebin(fRebinFactor);
    fTempHist=hDATA;
    TH1F *hGLAU = new TH1F("hGLAU","hGLAU",hDATA->GetNbinsX(),0,hDATA->GetNbinsX()*hDATA->GetBinWidth(1));
    hGLAU->Sumw2();

    // Minimize here
    if(gMinuit) delete gMinuit;
    new TMinuit(3);
    gMinuit->mncler();
    gMinuit->SetFCN(AliCentralityGlauberFit::MinuitFcnNBD);
    gMinuit->SetObjectFit(this); 

    Double_t arglist[2]={0};
    Int_t ierflg;
    if (fUseChi2) arglist[0] = 1; // should be 1 if you want to minimize chi2
    else arglist[0] = 0.5;        // should be 0.5 for ll
    gMinuit->mnexcm("SET ERR",arglist, 1, ierflg);
    gMinuit->mnparm(0,"alpha", alpha,  (fAlphahigh-fAlphalow)/fNalpha,  fAlphalow, fAlphahigh, ierflg);
    gMinuit->mnparm(1,"mu"   , mu,     (fMuhigh-fMulow)/fNmu, fMulow, fMuhigh, ierflg);
    gMinuit->mnparm(2,"k"    , k,      (fKhigh-fKlow)/fNk, fKlow, fKhigh, ierflg);

    // Call migrad
    arglist[0] = 100; // max calls
    arglist[1] = 0.1; // tolerance    
    gMinuit->mnexcm("SIMPLEX",arglist, 2, ierflg);
    arglist[0] = 1000; // max calls
    arglist[1] = 0.1;  // tolerance    
    gMinuit->mnexcm("MIGrad",arglist, 2, ierflg);
    //gMinuit->mnexcm("IMPROVE",arglist, 0, ierflg);

    if (ierflg != 0) {
      AliWarning("Abnormal termination of minimization.");
    }
    
    // ______________________ Get chi2 and Fit Status __________
    Double_t amin,edm,errdef;
    Int_t    nvpar, nparx, icstat;
    
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    Double_t chi2min = amin;
    std::cout << "Fit status " << icstat << std::endl;

    Double_t alpha_min, mu_min, k_min;
    Double_t alpha_mine, mu_mine, k_mine;
    gMinuit->GetParameter(0, alpha_min , alpha_mine );
    gMinuit->GetParameter(1, mu_min    , mu_mine    );
    gMinuit->GetParameter(2, k_min     , k_mine     );

    // print some infos
    std::cout << "chi2 min is " << chi2min << ", " << alpha_min << ", "<< mu_min<< ", "  
              << k_min << std::endl;

    thistGlau = GlauberHisto(mu_min,k_min,alpha_min,hDATA,kTRUE);

    hGLAU = (TH1F *) thistGlau->Clone("hGLAU");
    hGLAU->SetName( ((TString)hDATA->GetName()).Append(Form("_GLAU")));
    hGLAU->SetTitle( ((TString)hDATA->GetName()).Append(Form("_GLAU_%.3f_%.3f_%.3f",
                                                             mu_min,k_min,alpha_min)));

    
    std::cout << "fitted " << hGLAU->Integral(hGLAU->FindBin(fMultmin),
                                              hGLAU->FindBin(fMultmax))/hGLAU->Integral() 
              << " of the total cross section" << std::endl; 

    Double_t mcintegral = hGLAU->Integral(hGLAU->FindBin(fScalemin),hGLAU->GetNbinsX());
    Double_t scale = (hDATA->Integral(hDATA->FindBin(fScalemin),hDATA->GetNbinsX())/mcintegral);
    hGLAU->Scale(scale);

    std::cout << "Chi2 final " << CalculateChi2(hDATA,hGLAU) << std::endl;

    fhEffi = GetTriggerEfficiencyFunction(hDATA, hGLAU);
    SaveHisto(hDATA,hGLAU,fhEffi,outrootfile);
    fGlauberHist=hGLAU;
  }
  // close inrootfile, outrootfile
  inrootfile->Close();
  outrootfile->Close();
}

//--------------------------------------------------------------------------------------------------
TH1F *AliCentralityGlauberFit::GlauberHisto(Double_t mu, Double_t k, Double_t alpha, 
                                            TH1F *hDATA, Bool_t save) 
{
  // Get Glauber histogram.
  static TH1F *h1 = (TH1F*)hDATA->Clone();
  h1->Reset();
  h1->SetName(Form("fit_%.3f_%.3f_%.3f",mu,k,alpha));

  if (fUseAverage) {
    fhAncestor = MakeAncestor(alpha);
    for (Int_t np=1; np<=fhAncestor->GetNbinsX(); ++np) {  
      Double_t nanc = fhAncestor->GetBinCenter(np);
      Double_t weights = fhAncestor->GetBinContent(np);

      if (weights <= 0) continue;
      Int_t trials = (Int_t) (20 * nanc * (int) mu);
      if (trials <=0) continue;
      for (Int_t j=0; j<trials; j++) {
       	double nbdvalue = NBD(j, mu * nanc, k * nanc);
       	h1->Fill((double) j, nbdvalue * weights);
      }
    }
    return h1;
  }

  TH1F *hSample = fFastFit != 2 ? NBDhist(mu,k) : 0;

  TFile *outFile = NULL;
  TNtuple *ntuple = NULL;
 
  if (save) {
    outFile = new TFile(fOutntuplename,"RECREATE");
    ntuple = new TNtuple("gnt", "Glauber ntuple", "Npart:Ncoll:B:tAA:ntot");
  } 

  Int_t nents = 0;
  if (fGlauntuple)
    nents = fGlauntuple->GetEntries(); 

  for (Int_t i=0;i<fNevents;++i) {
    if (fGlauntuple)
      fGlauntuple->GetEntry(i % nents);
    else {
      fNpart = 2;
      fNcoll = 1;
    }
    Int_t n=0;
    //if (fAncestor == 1)      n = (Int_t)(TMath::Power(fNpart,alpha));
    if (fAncestor == 1)      n = (Int_t)(TMath::Power(fNcoll,alpha));
    else if (fAncestor == 2) n = (Int_t)(alpha * fNpart + (1-alpha) * fNcoll);
    else if (fAncestor == 3) n = (Int_t)((1-alpha) * fNpart/2 + alpha * fNcoll);

    Int_t ntot=0;
    if (fFastFit == 1) {
      ntot = (Int_t)(n*hSample->GetRandom()); // NBD
    }
    else if (fFastFit == 2) {
      Double_t sigma = k*TMath::Sqrt(n*mu);    
      ntot = (Int_t)(gRandom->Gaus(n*mu,sigma)); // Gaussian
    }
    else {
      for(Int_t j = 0; j<(Int_t)n; ++j) 
	ntot += (Int_t)hSample->GetRandom();
    }
    h1->Fill(ntot);    

    if (save) 
      ntuple->Fill(fNpart,fNcoll,fB,fTaa,ntot);
  }

  if (save) {
    ntuple->Write();
    outFile->Close();
  }
  
  if (fFastFit != 2) delete hSample;
  return h1;
}

//--------------------------------------------------------------------------------------------------
Double_t AliCentralityGlauberFit::CalculateChi2(TH1F *hDATA, TH1F *thistGlau) 
{
  // Note that for different values of neff the mc distribution is identical, just re-scaled below
  // normalize the two histogram, scale MC up but leave data alone - that preserves error bars !!!!
  
  Int_t lowchibin =   hDATA->FindBin(fMultmin);
  Int_t highchibin =  hDATA->FindBin(fMultmax);

  Double_t mcintegral = thistGlau->Integral(1,thistGlau->GetNbinsX());
  Double_t scale = (hDATA->Integral(1,hDATA->GetNbinsX())/mcintegral);
  thistGlau->Scale(scale);

  // calculate the chi2 between MC and real data over some range ????
  if (fUseChi2) {
    Double_t chi2 = 0.0;
    for (Int_t i=lowchibin; i<=highchibin; i++) {
      if (hDATA->GetBinContent(i) < 1.0) continue;
      Double_t diff = TMath::Power((thistGlau->GetBinContent(i) - hDATA->GetBinContent(i)),2);
      diff = diff / (TMath::Power(hDATA->GetBinError(i),2) + TMath::Power(thistGlau->GetBinError(i),2)); 
      chi2 += diff;
    }
    chi2 = chi2 / (highchibin - lowchibin + 1);
    return chi2;
  }
  // "-2 log likelihood ratio(mu;n) = 2[(mu - n) + n ln(n/mu)]"
  else {
    std::cout << "LL" << std::endl;

    Double_t ll = 0.0;
    for (Int_t i=lowchibin; i<=highchibin; i++) {
      Double_t data = hDATA    ->GetBinContent(i);
      Double_t mc   = thistGlau->GetBinContent(i);
      Int_t idata = TMath::Nint(data);
      if (mc < 1.e-9) mc = 1.e-9;
      Double_t fsub = - mc + idata * TMath::Log(mc);
      Double_t fobs = 0;
      if (idata > 0) {
	for(Int_t istep = 0; istep < idata; istep++){
	  if (istep > 1) 
	    fobs += TMath::Log(istep);
	}
      }
      fsub -= fobs;
      ll -=  fsub ;
    }
    return 2*ll;
  }
}

//--------------------------------------------------------------------------------------------------
TH1F *AliCentralityGlauberFit::GetTriggerEfficiencyFunction(TH1F *hist1, TH1F *hist2) 
{
  // Get efficiency.

  fhEffi = (TH1F*)hist1->Clone("heffi");
  fhEffi->Divide(hist2);
  return fhEffi;
}

//--------------------------------------------------------------------------------------------------
Double_t AliCentralityGlauberFit::GetTriggerEfficiencyIntegral(TH1F *hist1, TH1F *hist2) 
{
  // Get eff integral.

  fEffi = hist1->Integral()/hist2->Integral();
  return fEffi;
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::SaveHisto(TH1F *hist1, TH1F *hist2, TH1F *heffi, TFile *outrootfile) 
{
  // Save histograms.

  outrootfile->cd();
  hist1->Write();
  hist2->Write();
  heffi->Write();
}

//--------------------------------------------------------------------------------------------------
Double_t AliCentralityGlauberFit::NBD(Int_t n, Double_t mu, Double_t k) const
{
  // Compute NBD.
  double F;
  double f;

  if (n+k > 100.0) {
    // log method for handling large numbers
    F  = TMath::LnGamma(n + k)- TMath::LnGamma(n + 1.)- TMath::LnGamma(k);
    f  = n * TMath::Log(mu/k) - (n + k) * TMath::Log(1.0 + mu/k);
    F = F+f;
    F = TMath::Exp(F);
  } else {
    F  = TMath::Gamma(n + k) / ( TMath::Gamma(n + 1.) * TMath::Gamma(k) );
    f  = n * TMath::Log(mu/k) - (n + k) * TMath::Log(1.0 + mu/k);
    f  = TMath::Exp(f);
    F *= f;
  }

  return F;

}

//--------------------------------------------------------------------------------------------------
TH1F *AliCentralityGlauberFit::NBDhist(Double_t mu, Double_t k)
{
  // Interface for TH1F.

  TH1F *h = new TH1F("htemp","",100,-0.5,299.5);
  h->SetName(Form("nbd_%f_%f",mu,k));
  h->SetDirectory(0);
  for (Int_t i=0;i<300;++i) {
    Double_t val = NBD(i,mu,k);
    if (val>1e-20) h->Fill(i,val);
  }
  return h;
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::MinuitFcnNBD(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // FCN for minuit.

  Double_t alpha = par[0];
  Double_t mu    = par[1];
  Double_t k     = par[2];

  if (0) { //avoid warning
    gin=gin;
    npar=npar;
    iflag=iflag;
  }
  AliCentralityGlauberFit * obj = (AliCentralityGlauberFit *) gMinuit->GetObjectFit();
  TH1F * thistGlau = obj->GlauberHisto(mu,k,alpha,obj->GetTempHist(),kFALSE);
  f = obj->CalculateChi2(obj->GetTempHist(),thistGlau);
  Printf("Minuit step: chi2=%f, alpha=%f,  mu=%f, k=%f\n",f,alpha,mu,k);
}

//--------------------------------------------------------------------------------------------------
Double_t AliCentralityGlauberFit::NBDFunc(const Double_t *x, const Double_t *par) 
{
  // TF1  interface.

  Double_t mu = par[0];
  Double_t k  = par[1];
  Double_t n  = x[0];
  Double_t ret = exp( TMath::LnGamma(n+k) - TMath::LnGamma(k) - TMath::LnGamma(n+1) ) * 
    TMath::Power(mu/(mu+k),n) * TMath::Power(1-mu/(mu+k),k);
  return ret;
}

//--------------------------------------------------------------------------------------------------
TH1F *AliCentralityGlauberFit::MakeAncestor(Double_t alpha)
{
  // Make ancestor histogram.

  TString hname(Form("fhAncestor_%.3f",alpha));
  if (fhAncestor) {
    if (hname.CompareTo(fhAncestor->GetName())==0)
      return fhAncestor;
  }

  delete fhAncestor;
  fhAncestor = 0;
  TFile *ancfile = TFile::Open(fAncfilename,"read");
  if (ancfile && ancfile->IsOpen()) {
    fhAncestor = dynamic_cast<TH1F*>(ancfile->Get(hname));
    if (fhAncestor) {
      fhAncestor->SetDirectory(0);
      delete ancfile;
      return fhAncestor;
    }
  }
  delete ancfile;
  
  fhAncestor = new TH1F(hname,hname,3000,0,3000);
  fhAncestor->SetDirectory(0);
  Int_t nents = fGlauntuple->GetEntries(); 
  for (Int_t i=0;i<nents;++i) {
    fGlauntuple->GetEntry(i % nents);
    Int_t n=0;
    if (fAncestor == 1)    n = (Int_t) (TMath::Power(fNpart,alpha));
    //if (fAncestor == 1)      n = (Int_t) (TMath::Power(fNcoll,alpha));
    else if (fAncestor == 2) n = (Int_t) (alpha * fNpart + (1-alpha) * fNcoll);
    else if (fAncestor == 3) n = (Int_t) ((1-alpha) * fNpart/2 + alpha * fNcoll);
    fhAncestor->Fill(n);
  }

  ancfile = TFile::Open(fAncfilename,"update");
  if (ancfile && ancfile->IsOpen()) {
    fhAncestor->Write();
  }
  delete ancfile;
  return fhAncestor;
}
