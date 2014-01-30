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
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include "AliCentralityGlauberFit.h"
#include "AliLog.h"

ClassImp(AliCentralityGlauberFit)  
 
//______________________________________________________________________________
AliCentralityGlauberFit::AliCentralityGlauberFit(const char *filename) : 
  fNk(0),     
  fKlow(0),   
  fKhigh(0),  
  fNalpha(0), 
  fAlphalow(0),   
  fAlphahigh(0),  
  fNsigma(0), 
  fSigmalow(0),   
  fSigmahigh(0),  
  fNbog(0), 
  fBoglow(0),   
  fBoghigh(0),  
  fNCP(0), 
  fCPlow(0),   
  fCPhigh(0),  
  fRebinFactor(0),
  fScalemin(0),    
  fMultmin(0),    
  fMultmax(0),    
  fGlauntuple(0), 
  fNpart(0),       
  fNcoll(0),       
  fB(0),           
  fTaa(0),         
  fTempHist(0),   
  fGlauberHist(0), 
  fUseChi2(kTRUE),      
  fNevents(100000),
  fInrootfilename(0),
  fInntuplename(0),
  fOutrootfilename(0),
  fOutntuplename(0),
  fAncfilename("ancestor_hists.root"),
  fHistnames(),
  fIsZN(kTRUE)
{
  // Standard constructor.
  TFile *f = 0;
  if (filename) {  // glauber file
    f = TFile::Open(filename);
    fGlauntuple = (TNtuple*) f->Get("nt_p_Pb");
    fGlauntuple->SetBranchAddress("Npart",&fNpart);
    fGlauntuple->SetBranchAddress("Ncoll",&fNcoll);
    fGlauntuple->SetBranchAddress("B",&fB);
    fGlauntuple->SetBranchAddress("tAA",&fTaa);
  }

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
					      Int_t   Nk,
					      Double_t klow,
					      Double_t khigh,
					      Int_t   Nalpha,
					      Double_t alphalow,
					      Double_t alphahigh,
					      Int_t   Nsigma,
					      Double_t sigmalow,
					      Double_t sigmahigh,
					      Int_t   Nbog,
					      Double_t boglow,
					      Double_t boghigh,
					      Int_t   NCP,
					      Double_t CPlow,
					      Double_t CPhigh)
{
  // Set Glauber parameters.
  fNk=Nk;
  fKlow=klow;
  fKhigh=khigh;
  fNalpha=Nalpha;
  fAlphalow=alphalow;
  fAlphahigh=alphahigh;
  fNsigma=Nsigma;
  fSigmalow=sigmalow;
  fSigmahigh=sigmahigh;
  fNbog=Nbog;
  fBoglow=boglow;
  fBoghigh=boghigh;
  fNCP=NCP;
  fCPlow=CPlow;
  fCPhigh=CPhigh;
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
    std::cout << " ->  Getting histogram " << *hni << std::endl;
    if (!hDATA) {
      TList *list  = (TList*) (inrootfile->Get("CentralityStat")); 
      //TList *list  = (TList*) (inrootfile->Get("VZEROEquaFactorStat")); 
      hDATA  = (TH1F*) (list->FindObject(*hni));
    } 
    //hDATA->Rebin(fRebinFactor);
    //TH1F *hGLAU = new TH1F("hGLAU","hGLAU",hDATA->GetNbinsX(),0,hDATA->GetNbinsX()*hDATA->GetBinWidth(1));

    Double_t chi2min = 9999999.0;
    Double_t CP_min=-1;
    Double_t bog_min=-1;
    Double_t sigma_min=-1;
    Double_t alpha_min=-1;
    Double_t k_min=-1;
    Double_t alpha, k, sigma, bog, CP, chi2;   

    for (Int_t nCP=0;nCP<fNCP;nCP++) {
      CP = fCPlow   + ((Double_t) nCP ) * (fCPhigh  - fCPlow ) / fNCP;

      for (Int_t nbog=0;nbog<fNbog;nbog++) {
	bog = fBoglow   + ((Double_t) nbog ) * (fBoghigh  - fBoglow ) / fNbog;
	
	for (Int_t nsigma=0;nsigma<fNsigma;nsigma++) {
	  sigma = fSigmalow   + ((Double_t) nsigma ) * (fSigmahigh  - fSigmalow ) / fNsigma;
	  
	  for (Int_t nalpha=0;nalpha<fNalpha;nalpha++) {
	    alpha = fAlphalow   + ((Double_t) nalpha ) * (fAlphahigh  - fAlphalow ) / fNalpha;
	    
	    for (Int_t nk=0; nk<fNk; nk++) {
	      k = fKlow  + ((Double_t) nk ) * (fKhigh  - fKlow ) / fNk;
	      
	      thistGlau = GlauberHisto(k,alpha,sigma,bog,CP,hDATA,kFALSE);
	      chi2 = CalculateChi2(hDATA,thistGlau);
	      //cout << "chi2 " << chi2<< endl;
	      if ( chi2 < chi2min ) {
		chi2min = chi2;
		alpha_min=alpha;
		sigma_min=sigma;
		bog_min=bog;
		CP_min=CP;
		k_min=k;
	      }
	    }
	  }
	}
      }
    }
    thistGlau->Reset();
    thistGlau = GlauberHisto(k_min,alpha_min,sigma_min,bog_min,CP_min, hDATA,kTRUE);

    TH1F * hGLAU = 0x0;
    hGLAU = (TH1F *) thistGlau->Clone("hGLAU");
    hGLAU->SetName( ((TString)hDATA->GetName()).Append(Form("_GLAU")));
    hGLAU->SetTitle( ((TString)hDATA->GetName()).Append(Form("_GLAU_%.3f_%.3f_%.3f_%.3f_%.3f",
                                                             k_min,alpha_min,sigma_min,bog_min,CP_min)));

    Int_t lastBin=0;
    if(fIsZN) lastBin = 100;//ZN
    else lastBin = hDATA->GetNbinsX(); // ZP
    Double_t mcintegral = hGLAU->Integral(2, lastBin,"width");
    Double_t dataintegral = hDATA->Integral(2, lastBin,"width");
    Double_t scale = (dataintegral/mcintegral);
    //
    std::cout << "hGLAU -> Integral in bin range:" << "1-" << lastBin << " = " << 
    	hGLAU->Integral(2, lastBin) << std::endl;
    std::cout << "hDATA -> Integral in bin range:" << "1-" << lastBin << " = " <<
    	hDATA->Integral(2, lastBin) << std::endl;
    //
    //printf(" binwidth: hGLAU %f  hDATA %f\n", hGLAU->GetBinWidth(100), hDATA->GetBinWidth(100));
    printf(" scale %f \n", scale);
    
    hGLAU->Scale(scale);
    //
    std::cout << "hGLAU -> Integral (whole range): " << hGLAU->Integral(1, hGLAU->GetNbinsX()) << std::endl;
    std::cout << "hDATA -> Integral (whole range): " << hDATA->Integral(1, hDATA->GetNbinsX()) << std::endl;

    SaveHisto(hDATA, hGLAU, outrootfile);
    //fclose (fTxt);

    printf("\n \t k = %1.2f  alpha = %1.2f  sigma = %1.2f  bog = %1.2f  CP = %1.2f \n\n", k,alpha, sigma, bog, CP);
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
TH1F *AliCentralityGlauberFit::GlauberHisto(Double_t k, Double_t alpha, Double_t sigma, Double_t bog, Double_t CP, 
                                            TH1F *hDATA, Bool_t save) 
{
  // Get Glauber histogram.
  static TH1F *h1 = (TH1F*)hDATA->Clone("h1");
  h1->Reset();
  h1->SetName(Form("fit_%.3f_%.3f",k,alpha));

  TFile *outFile = NULL;
  TNtuple *ntuple = NULL;
 
  if (save) {
    outFile = new TFile(fOutntuplename,"RECREATE");
    ntuple = new TNtuple("gnt", "Glauber ntuple", "Npart:Ncoll:B:tAA:ntot:nbn:ngn:Etot");
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

    // Slow Nucleon Model from Chiara
    Double_t ntot=0., n=0.;
    Double_t nbn=0., ngn=0., nbp=0., ngp=0.;
    MakeSlowNucleons2(fNcoll,alpha,k,bog,CP, nbn,ngn,nbp,ngp);
    //MakeSlowNucleons2s(fNcoll,alpha,k,bog,CP, nbn,ngn,nbp,ngp);

    // acceptance
    //
    if(fIsZN) n = alpha*ngn+alpha*nbn; // ZNA
    else n =  alpha*ngp+alpha*nbp; // ZPA
    //----------------------------------------

    //------ experimental resolution -> Gaussian smearing ------------------------------------
    Double_t resexp=0;
    //if (n>0) resexp = sigma*gRandom->Gaus(0.0,1.0)/TMath::Sqrt(n);
    //else resexp=0;
    //ntot = n*(1+resexp);
        
    //if(n>0) resexp = sigma*TMath::Sqrt(n);    
    //else resexp=0;
    ntot = (int) (gRandom->Gaus(n, resexp));
    //----------------------------------------

    // non-lineary effect -------------------------------
    //ntot = k*ntot;
    //ntot = ntot + k*ntot*ntot;

    //cout << ngn << " " << nbn << " "  << ntot << endl;

    Double_t nFact = 1.577;
    Double_t Etot = nFact*ntot;
    //

    if (n>0)
    resexp = 1./TMath::Sqrt(n)*sigma*gRandom->Gaus(0.0,1.0);
    
    Etot = Etot*(1+resexp);
    //printf("  ntot %f  Etot %f \n", ntot, Etot);
    
    if(ntot>0) {
      h1->Fill(Etot);     
      if (save) 
	ntuple->Fill(fNpart, fNcoll, fB, fTaa, ntot, nbn, ngn, Etot);
    }
  }

  if (save) {
    ntuple->Write();
    outFile->Close();
  }
  return h1;
}

//--------------------------------------------------------------------------------------------------
Double_t AliCentralityGlauberFit::CalculateChi2(TH1F *hDATA, TH1F *thistGlau) 
{
  // Note that for different values of neff the mc distribution is identical, just re-scaled below
  // normalize the two histogram, scale MC up but leave data alone - that preserves error bars !!!!
  
  Int_t lowchibin =   hDATA->FindBin(fMultmin);
  Int_t highchibin =  hDATA->FindBin(fMultmax);

  Double_t mcintegral = thistGlau->Integral(1, thistGlau->GetNbinsX());
  Double_t scale = (hDATA->Integral(1, hDATA->GetNbinsX())/mcintegral);
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
void AliCentralityGlauberFit::SaveHisto(TH1F *hist1, TH1F *hist2, TFile *outrootfile) 
{
  // Save histograms.

  outrootfile->cd();
  hist1->Write();
  hist2->Write();
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::MakeSlowNucleons(Int_t ncoll, Double_t alpha, Double_t k, Double_t &nbn, Double_t &ngn)
{
// from AliGenSlowNucleonModelExp (Chiara Oppedisano)
// Return the number of black and gray nucleons
//
// Number of collisions

  Int_t fP = 82;
  Int_t fN = 126;
  Float_t nu = (Float_t) ncoll; 
  //
// Mean number of gray nucleons 
    Float_t nGray         = alpha * nu;
    Float_t nGrayNeutrons = nGray * fN / (fN + fP);
    Float_t nGrayProtons  = nGray - nGrayNeutrons;

// Mean number of black nucleons 
    Float_t nBlack  = 0.;
    if(!fApplySaturation || (fApplySaturation && nGray<fnGraySaturation)) nBlack = k * nu;
    else if(fApplySaturation && nGray>=fnGraySaturation) nBlack = fnBlackSaturation;
    Float_t nBlackNeutrons = nBlack * 0.84;
    Float_t nBlackProtons  = nBlack - nBlackNeutrons;

// Actual number (including fluctuations) from binomial distribution
    Double_t p, ngp, nbp;
    
//  gray neutrons
    p =  nGrayNeutrons/fN;
    ngn = gRandom->Binomial((Int_t) fN, p);
    
//  gray protons
    p =  nGrayProtons/fP;
    ngp = gRandom->Binomial((Int_t) fP, p); 

//  black neutrons
    p =  nBlackNeutrons/fN;
    nbn = gRandom->Binomial((Int_t) fN, p);
    
//  black protons
    p =  nBlackProtons/fP;
    nbp = gRandom->Binomial((Int_t) fP, p);

}


//--------------------------------------------------------------------------------------------------
 void AliCentralityGlauberFit::MakeSlowNucleons2(Int_t ncoll, Double_t alpha, Double_t k, Double_t bog, Double_t CP, Double_t &nbn, Double_t &ngn, Double_t &nbp, Double_t &ngp)
{
// from AliGenSlowNucleonModelExp (Chiara Oppedisano)
//
// Return the number of black and gray nucleons
//
// Number of collisions

   // based on E910 model ================================================================

  Int_t fP = 82;
  Int_t fN = 126;
  
  Float_t nu = (Float_t) ncoll; 
  //
  nu = gRandom->Gaus(nu,0.5); 
  if(nu<1.) nu = 1.;
  //
  //float sdp = gRandom->Rndm();
  //if(nu==1. && sdp<=0.2) return;

  //  gray protons
  Float_t  poverpd = 0.843; 
  Float_t  zAu2zPb = 82./79.;
  Float_t  nGrayp = (-0.27 + 0.63 * nu - 0.0008 *nu *nu)*poverpd*zAu2zPb;

  Double_t p;
  p =  nGrayp/fP;
  ngp = gRandom->Binomial((Int_t) fP, p);
  //ngp = gRandom->Gaus(nGrayp, TMath::Sqrt(fP*p*(1-p)));
  if(nGrayp<0.) ngp=0;

  //  black protons
  //Float_t blackovergray = 3./7.;// from spallation
  //Float_t blackovergray = 0.65; // from COSY
  Float_t blackovergray = bog;
   // Float_t blackovergray = 0.65; // from COSY
  Float_t nBlackp  = blackovergray*nGrayp; 

  p =  nBlackp/fP;
  nbp = gRandom->Binomial((Int_t) fP, p);
  //nbp = gRandom->Gaus(nBlackp, TMath::Sqrt(fP*p*(1-p)));
  if(nBlackp<0.) nbp=0;
  
  if(nu<3.){
    nGrayp = -0.836 + 0.9112 *nu - 0.05381 *nu *nu;
    nBlackp  = blackovergray*nGrayp; 
  }
  

  //  gray neutrons
  Float_t nGrayNeutrons = 0.;
  Float_t nBlackNeutrons = 0.;
  Float_t cp = (nGrayp+nBlackp)/CP;

  // if(cp>0.){
  //   //Float_t nSlow  = 51.5+469.2/(-8.762-cp);
  //   Float_t nSlow  = 60.0+469.2/(-8.762-cp);
  //   //if(cp<2.5) nSlow = 1+(9.9-1)/(2.5-0)*(cp-0);
  //   if(cp<3.) nSlow = 0.+(11.6-0.)/(3.-0.)*(cp-0.);
    
  //   nGrayNeutrons = nSlow * 0.1; 
  //   nBlackNeutrons = nSlow - nGrayNeutrons;
  // }
  if(cp>0.){
    Float_t paramnSlow[3] = {61., 470., 7.};
    Float_t nSlow  = paramnSlow[0]+paramnSlow[1]/(-paramnSlow[2]-cp);
    float paramRetta = paramnSlow[0]+paramnSlow[1]/(-paramnSlow[2]-3);
    if(cp<3.) nSlow = 0.+(paramRetta-0.)/(3.-0.)*(cp-0.);
    
    nGrayNeutrons = nSlow * 0.1;
    nBlackNeutrons = nSlow - nGrayNeutrons;
  }
  else{
    // Sikler "pasturato"
    nGrayNeutrons = 0.47 * 2.2 *  nu; // fAlphaGray=2.3
    nBlackNeutrons = 0.88 * 3.6 * nu; // fAlphaBlack=3.6  
    if(nGrayNeutrons<0.) nGrayNeutrons=0.;
    if(nBlackNeutrons<0.) nBlackNeutrons=0.;
    //printf("nslowp=0 -> ncoll = %1.0f -> ngrayn = %1.0f  nblackn = %1.0f \n", nu, nGrayNeutrons, nBlackNeutrons);
  }

  p =  nGrayNeutrons/fN;
  //ngn = gRandom->Binomial((Int_t) fN, p);
  ngn = gRandom->Gaus(nGrayNeutrons, TMath::Sqrt(fN*p*(1-p)));

  //  black neutrons
  p =  nBlackNeutrons/fN;
  //nbn = gRandom->Binomial((Int_t) fN, p);
  nbn = gRandom->Gaus(nBlackNeutrons, TMath::Sqrt(fN*p*(1-p)));
}

//--------------------------------------------------------------------------------------------------
 void AliCentralityGlauberFit::MakeSlowNucleons2s(Int_t ncoll, Double_t alpha, Double_t k, Double_t bog, Double_t CP, Double_t &nbn, Double_t &ngn, Double_t &nbp, Double_t &ngp)
{
// from AliGenSlowNucleonModelExp (Chiara Oppedisano)
//
// Return the number of black and gray nucleons
//
// Number of collisions

   // based on E910 model ================================================================

  Int_t fP = 82;
  Int_t fN = 126;
  
  Float_t nu = (Float_t) ncoll; 
  Float_t sigmap = 0.25;

  //  gray protons
  Float_t  poverpd = 0.843; 
  Float_t  zAu2zPb = 82./79.;
  Float_t  grayp = (-0.27 + 0.63 * nu - 0.0008 *nu *nu)*poverpd*zAu2zPb;
  Float_t  nGrayp = gRandom->Gaus(grayp, sigmap);
  if(nGrayp<0.) nGrayp=0.;

  Double_t p=0.;
  p =  nGrayp/fP;
  ngp = gRandom->Binomial((Int_t) fP, p);
  //ngp = gRandom->Gaus(nGrayp, TMath::Sqrt(fP*p*(1-p)));
  if(nGrayp<0.) ngp=0;

  //  black protons
  //Float_t blackovergray = 3./7.;// from spallation
  Float_t blackovergray = 0.65; // from COSY
  //Float_t blackovergray = bog;
  Float_t blackp  = blackovergray*nGrayp; 
  Float_t nBlackp = gRandom->Gaus(blackp, sigmap);
  if(nBlackp<0.) nBlackp=0.;

  p =  nBlackp/fP;
  nbp = gRandom->Binomial((Int_t) fP, p);
  //nbp = gRandom->Gaus(nBlackp, TMath::Sqrt(fP*p*(1-p)));
  if(nBlackp<0.) nbp=0;
  
  //  gray neutrons
  Float_t nGrayNeutrons = 0.;
  Float_t nBlackNeutrons = 0.;
  Float_t cp = (nGrayp+nBlackp)/CP;

  if(cp>0.){
    Float_t paramnSlow[3] = {61., 470., 7.};
    Float_t nSlow  = paramnSlow[0]+paramnSlow[1]/(-paramnSlow[2]-cp);
    //float paramRetta = paramnSlow[0]+paramnSlow[1]/(-paramnSlow[2]-3);
    //if(cp<3.) nSlow = 0.+(paramRetta-0.)/(3.-0.)*(cp-0.);
    
    nGrayNeutrons = nSlow * 0.1;
    nBlackNeutrons = nSlow - nGrayNeutrons;
  }
  else{
    // Sikler "pasturato"
    nGrayNeutrons = 0.47 * 2.2 *  nu; // fAlphaGray=2.3
    nBlackNeutrons = 0.88 * 3.6 * nu; // fAlphaBlack=3.6  
    //printf("nslowp=0 -> ncoll = %1.0f -> ngrayn = %1.0f  nblackn = %1.0f \n", nu, nGrayNeutrons, nBlackNeutrons);
  }
  //
  if(nGrayNeutrons<0.) nGrayNeutrons=0.;
  if(nBlackNeutrons<0.) nBlackNeutrons=0.;

  p =  nGrayNeutrons/fN;
  //ngn = gRandom->Binomial((Int_t) fN, p);
  ngn = gRandom->Gaus(nGrayNeutrons, TMath::Sqrt(fN*p*(1-p)));
  if(p<0.) ngn=0;

  //  black neutrons
  p =  nBlackNeutrons/fN;
  //nbn = gRandom->Binomial((Int_t) fN, p);
  nbn = gRandom->Gaus(nBlackNeutrons, TMath::Sqrt(fN*p*(1-p)));
  if(p<0.) nbn=0;
}

Double_t AliCentralityGlauberFit::ConvertToEnergy(Double_t T)
{ 
   TDatabasePDG * pdg = TDatabasePDG::Instance();
    const Double_t kMassNeutron = pdg->GetParticle(kNeutron)->Mass();
    Double_t fPmax =10;   

    Double_t m, p, pmax, ekin, f;
    m = kMassNeutron;
    pmax = TMath::Sqrt(2*T*(T+TMath::Sqrt(T*T+m*m)));
    
    do
      {
	p = gRandom->Rndm() * fPmax;
	f = Maxwell(m, p, T) / Maxwell(m , pmax, T);
      }
    while(f < gRandom->Rndm());

    return 13.5*TMath::Sqrt(p*p+m*m);
}

Double_t AliCentralityGlauberFit::Maxwell(Double_t m, Double_t p, Double_t T)
{
/* Relativistic Maxwell-distribution */
    Double_t ekin;
    ekin = TMath::Sqrt(p*p+m*m)-m;
    return (p*p * exp(-ekin/T));
}

