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

#include <TNamed.h>
#include <TH1D.h>
#include <TString.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom1.h>
#include <TROOT.h>
#include <TH2F.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <vector>
#include <TMinuit.h>
#include <TCanvas.h>
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
  fNgamma(0), 
  fgammalow(0),   
  fgammahigh(0),  
  fNsigmap(0), 
  fsigmaplow(0),   
  fsigmaphigh(0),  
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
  fIsZN(kTRUE),
  fNTrials(0),
  fChi2Min(2000.),
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
    fGlauntuple = (TNtuple*) f->Get("nt_p_Pb");
    fGlauntuple->SetBranchAddress("Npart",&fNpart);
    fGlauntuple->SetBranchAddress("Ncoll",&fNcoll);
    fGlauntuple->SetBranchAddress("B",&fB);
    fGlauntuple->SetBranchAddress("tAA",&fTaa);
  }
  for(int i=0; i<4; i++) fParamnSlow[i] = 0.;

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
					      Int_t   Ngamma,
					      Double_t gammalow,
					      Double_t gammahigh,
					      Int_t   Nsigmap,
					      Double_t sigmaplow,
					      Double_t sigmaphigh)
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
  fNgamma=Ngamma;
  fgammalow=gammalow;
  fgammahigh=gammahigh;
  fNsigmap=Nsigmap;
  fsigmaplow=sigmaplow;
  fsigmaphigh=sigmaphigh;
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::SetNParam(Double_t a, Double_t b, Double_t c, Double_t d)
{
   // Setting parameters that are adjusted
   fParamnSlow[0] = a;
   fParamnSlow[1] = b;
   fParamnSlow[2] = c;
   fParamnSlow[3] = d;
   for(int i=0; i<4; i++) printf("  INITIAL PARAMETER VALUES : paramnSlow[%d] = %f\n", i, fParamnSlow[i]);
}

//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::MakeFits() 
{
  // Make fits.

  TH1F *hDATA = 0x0;
  TH1F *thistGlau; 
  TFile *inrootfile;
  TFile *outrootfile;  
  
  // open inrootfile, outrootfile
  std::cout << "input file "  << fInrootfilename  << std::endl;
  std::cout << "output file " << fOutrootfilename << std::endl << std::endl;
  inrootfile  = new TFile(fInrootfilename,"OPEN");
  outrootfile = new TFile(fOutrootfilename,"RECREATE");
  
  // loop over all distribution names  
  std::vector<TString>::const_iterator hni;
  for(hni=fHistnames.begin(); hni!=fHistnames.end(); hni++) {
    hDATA  = (TH1F*) (inrootfile->Get(*hni)); 
    std::cout << " ->  Getting histogram " << *hni << std::endl << std::endl;
    if (!hDATA) {
      printf("   @@@@@@ No DATA histo in input -> EXITING!!!!!\n\n");
      return;
      //TList *list  = (TList*) (inrootfile->Get("CentralityStat")); 
      //hDATA  = (TH1F*) (list->FindObject(*hni));
    } 
    //hDATA->Rebin(fRebinFactor);
    //TH1F *hGLAU = new TH1F("hGLAU","hGLAU",hDATA->GetNbinsX(),0,hDATA->GetNbinsX()*hDATA->GetBinWidth(1));

    Double_t gamma_min=-1;
    Double_t sigmap_min=-1;
    Double_t bog_min=-1;
    Double_t sigma_min=-1;
    Double_t alpha_min=-1;
    Double_t k_min=-1, chi2=0.;
    Double_t alpha=0., k=0., sigma=0., bog=0., gamma=0., sigmap=0.;   
    
    for(Int_t nsigmap=0; nsigmap<fNsigmap; nsigmap++) {
      sigmap = fsigmaplow   + ((Double_t) nsigmap ) * (fsigmaphigh  - fsigmaplow ) / fNsigmap;
    
     for(Int_t ngamma=0; ngamma<fNgamma; ngamma++) {
      gamma = fgammalow   + ((Double_t) ngamma ) * (fgammahigh  - fgammalow ) / fNgamma;

      for(Int_t nbog=0; nbog<fNbog; nbog++) {
	bog = fBoglow   + ((Double_t) nbog ) * (fBoghigh  - fBoglow ) / fNbog;
	
	for(Int_t nsigma=0; nsigma<fNsigma; nsigma++) {
	  sigma = fSigmalow   + ((Double_t) nsigma ) * (fSigmahigh  - fSigmalow ) / fNsigma;
	  
	  for(Int_t nalpha=0; nalpha<fNalpha; nalpha++) {
	    alpha = fAlphalow   + ((Double_t) nalpha ) * (fAlphahigh  - fAlphalow ) / fNalpha;
	    
	    for(Int_t nk=0; nk<fNk; nk++) {
	      k = fKlow  + ((Double_t) nk ) * (fKhigh  - fKlow ) / fNk;


	       for(int i=0; i<fNTrials; i++){
	         printf("  **** calling GlauberHisto with acc = %1.3f R = %1.4f k = %1.4f bog = %1.3f  gamma = %1.3f sigmap = %1.3f\n",
		 	alpha, sigma, k, bog, gamma, sigmap);
	         if(fIsZN){
		   // Updating parameters for slow n 
		   fParamnSlow[0] = fParamnSlow[0]; // + 0.5*i;
		   fParamnSlow[1] = fParamnSlow[1]; // + 10.*i;
		   fParamnSlow[2] = fParamnSlow[2] + 0.05*i; 
		   fParamnSlow[3] = fParamnSlow[3];// + 0.01*i; 
		   //
		   printf("  \t\t trial #%d  n slow param.: %1.2f %1.2f %1.2f %1.2f\n",
		 	i,fParamnSlow[0], fParamnSlow[1], fParamnSlow[2], fParamnSlow[3]); 
                 }
		 
	         //thistGlau = GlauberHisto(k,alpha,sigma,bog,gamma,sigmap,hDATA,kFALSE);
	         thistGlau = GlauberHisto(k,alpha,sigma,bog,gamma,sigmap,hDATA,kTRUE);
		 //
	         chi2 = CalculateChi2(hDATA, thistGlau);
	         printf("   ->  Chi2 %f \n", chi2);
		 if(chi2 < fChi2Min){
		   fChi2Min = chi2;
		   alpha_min = alpha;
		   sigma_min = sigma;
		   bog_min = bog;
		   gamma_min = gamma;
		   sigmap_min = sigmap;
		   k_min = k;
		 }
	       }
	    }
	  }
	}
      }
     }
    }
    if(fNsigmap>1 || fNgamma>1 || fNbog>1 || fNsigma>1 || fNalpha>1 || fNk>1 || fNTrials>1){
      thistGlau->Reset();
      printf("  **** calling GlauberHisto with acc = %1.3f R = %1.4f bog = %1.3f  gamma = %1.3f sigmap = %1.3f\n",
	 	alpha, sigma, bog, gamma, sigmap);
      thistGlau = GlauberHisto(k_min,alpha_min,sigma_min,bog_min,gamma_min, sigmap_min, hDATA,kTRUE);
    }
    
    TH1F * hGLAU = 0x0;
    hGLAU = (TH1F *) thistGlau->Clone("hGLAU");
    
    hGLAU->SetName( ((TString)hDATA->GetName()).Append(Form("_GLAU")));
    if(fIsZN) hGLAU->SetTitle(((TString)hDATA->GetName()).Append(Form("_GLAU_%.3f_%.3f_%.3f_%.2f-%.2f-%.3f-%.2f",
                                      bog,sigmap,gamma,fParamnSlow[0], fParamnSlow[1], fParamnSlow[2], fParamnSlow[3])));
    else hGLAU->SetTitle(((TString)hDATA->GetName()).Append(Form("_GLAU_%.3f_%.3f_%.3f",
                                      alpha,bog,sigmap)));

    Int_t firstBin=1;
    if(fIsZN) firstBin = 8;//ZN
    Int_t lastBin=0;
    if(fIsZN) lastBin = 800;//ZN
    else lastBin = hDATA->GetNbinsX(); // ZP
    //
    Double_t mcintegral = hGLAU->Integral(firstBin, lastBin,"width");
    Double_t dataintegral = hDATA->Integral(firstBin, lastBin,"width");
    Double_t scale = 1.;
    if(mcintegral>0.) (scale = dataintegral/mcintegral);
    printf("\n scale %f \n", scale);
    //
    std::cout << "hGLAU -> Integral in bin range:" << firstBin << "-" << lastBin << " = " << 
    	mcintegral << std::endl;
    std::cout << "hDATA -> Integral in bin range:" << firstBin << "-" << lastBin << " = " <<
    	dataintegral << std::endl;
    
    hGLAU->Scale(scale);
    //
    std::cout << "hGLAU -> Integral (whole range): " << hGLAU->Integral(1, hGLAU->GetNbinsX(),"width") << std::endl;
    std::cout << "hDATA -> Integral (whole range): " << hDATA->Integral(1, hDATA->GetNbinsX(),"width") << std::endl;
    
    TCanvas *c2 = new TCanvas("c2","checkGLAU",100,0,700,700);
    c2->cd();
    hGLAU->SetLineColor(kPink);
    hGLAU->Draw("");
    hDATA->Draw("same");

    SaveHisto(hDATA, hGLAU, outrootfile);

    printf("\n \t k = %1.4f  alpha = %1.3f  sigma = %1.2f  bog = %1.4f  gamma = %1.4f sigmap = %.3f \n", k,alpha, sigma, bog, gamma, sigmap);
    if(fIsZN) printf(" \t a = %1.2f  b = %1.3f  c = %1.2f  d = %1.2f\n\n", fParamnSlow[0], fParamnSlow[1], fParamnSlow[2], fParamnSlow[3]);
    std::cout << "chi2 min is " << fChi2Min << std::endl << std::endl;
    std::cout << "fitted " << hGLAU->Integral(hGLAU->FindBin(fMultmin),
                                              hGLAU->FindBin(fMultmax))/hGLAU->Integral() 
              << " of the total cross section" << std::endl << std::endl; 
    fTempHist = hDATA;
    fGlauberHist = hGLAU;
  }

  // close inrootfile, outrootfile
  inrootfile->Close();
  outrootfile->Close();
}


//--------------------------------------------------------------------------------------------------
TH1F *AliCentralityGlauberFit::GlauberHisto(Double_t k, Double_t alpha, Double_t sigma, Double_t bog, 
		Double_t gamma, Double_t sigmap, TH1F *hDATA, Bool_t save) 
{
  // Get Glauber histogram.
  static TH1F *h1 = (TH1F*) hDATA->Clone("h1");
  h1->Reset();
  h1->SetName(Form("fit_%.3f_%.3f",k,alpha));
  
  TFile *outFile = NULL;
  TNtuple *ntuple = NULL;
  
  Double_t ntot=0, nblack=0, ngray=0, Etot=0, lcp=0, nslowp=0, nslown=0;
  
  if(save){
    outFile = new TFile(fOutntuplename,"RECREATE");
    ntuple = new TNtuple("gnt", "Glauber ntuple", "fNpart:fNcoll:fB:fTaa:ntot:nblack:ngray:Etot:lcp:nslowp:nslown");
  } 

  Int_t nents=0, evtot=0, evn=0, evzeron=0, evzerop=0, evenzero=0;
  if(fGlauntuple) nents = fGlauntuple->GetEntries(); 

  for (Int_t i=0;i<fNevents;++i) {
    if(fGlauntuple) fGlauntuple->GetEntry(i % nents);
    else{
      printf(" NOT GETTING ENTRY FROM TNTUPLE for ev. %d -> setting Ncoll=1 \n\n",i);
      fNpart = 2;
      fNcoll = 1;
    }
    //printf(" Getting entry %d from ntuple  Ncoll = %f\n",i%nents, fNcoll);

    // Slow Nucleon Model from Chiara
    Double_t n=0., nbn=0., ngn=0., nbp=0., ngp=0.;
    //
    //MakeSlowNucleons(fNcoll,nbn,ngn,nbp,ngp); // pure Sikler
    //MakeSlowNucleons2(fNcoll,bog,gamma, nbn,ngn,nbp,ngp,lcp); // smear in Ncoll
    MakeSlowNucleons2s(fNcoll,bog,gamma,sigmap, nbn,ngn,nbp,ngp,lcp);  // smear in Nslow
    
    evtot++; 
    if(fNcoll==1) evn++;
    //
    if(fIsZN){
      n = alpha*(ngn+nbn); 
      nblack = nbn;
      ngray = ngn;
    }
    else{
      n = alpha*(ngp+nbp); // ZPA
      nblack = nbp;
      ngray = ngp;
    }    
    nslown = nbn+ngn;
    nslowp = nbp+ngp;

    if((alpha*(ngn+nbn)<0.5)) evzeron++;
    if((alpha*(ngp+nbp)<0.5)) evzerop++;

    //------ experimental resolution -> Gaussian smearing ------------------------------------
    Double_t resexp = 0.;
    //if (n>0) resexp = sigma*gRandom->Gaus(0.0,1.0)/TMath::Sqrt(n);
    //else resexp=0;
    //ntot = n*(1+resexp);
        
    if(n>=0){
      resexp = sigma*TMath::Sqrt(n);
      ntot = (gRandom->Gaus(n, resexp));
    }
    //else n=0.;
    
    // non-lineary effect -------------------------------
    //ntot = ntot + k*ntot*ntot;

    Double_t nFact = 4.*82./208.;
    Etot = nFact*n;
        
    // USE ONLY IF ZDC TDC IS REQUESTED TO CONDITION THE SPECTRUM!!!!!!!!! 
    float resexpe=0.;
    if(n>0.){
      resexpe = 1./TMath::Sqrt(n)*sigma*gRandom->Gaus(0.0, 1.0);
      //resexpe = 1./TMath::Sqrt(n)*sigma;
      //
      Etot = Etot*(1+resexpe);
    }
    // non-lineary effect -------------------------------
    //Etot = Etot + k*Etot*Etot;

    if(Etot<0.5) evenzero++;
     
    h1->Fill(Etot);	
    if(save) ntuple->Fill(fNpart, fNcoll, fB, fTaa, ntot, nblack, ngray, Etot, lcp, nslowp, nslown);
  }

  if(save) {
    //printf("\n ******** Fraction of events with Ncoll=1 in the ntuple %f\n",(float) (1.*evn/evtot));
    if(fIsZN){
      printf(" ******** Fraction of events with no neutrons %f\n",(float) (1.*evzeron/evtot));
      printf(" ******** Fraction of events with no neutron energy %f\n",(float) (1.*evenzero/evtot));
      printf(" DATA: fraction of events with no ZN signal = %f\n\n", 1-alpha);
    }
    else printf("\n ******** Fraction of events with no protons %f \n",(float) (1.*evzerop/evtot));
    
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
  
  Int_t lowchibin  =  8;
  Int_t highchibin =  0;
  if(fIsZN) highchibin = 800;
  else highchibin = 25;

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
void AliCentralityGlauberFit::MakeSlowNucleons(Int_t ncoll, Double_t &nbn, Double_t &ngn,Double_t &nbp, Double_t &ngp)
{
// from AliGenSlowNucleonModelExp (Chiara Oppedisano)
// Return the number of black and gray nucleons
//
// Number of collisions

  Float_t fP = 82;
  Float_t fN = 126;
  Float_t nu = (Float_t) ncoll; 
  //
// Mean number of gray nucleons 
    Float_t nGray         = 2. * nu;
    Float_t nGrayNeutrons = nGray * fN / (fN + fP);
    Float_t nGrayProtons  = nGray - nGrayNeutrons;

// Mean number of black nucleons 
    Float_t nBlack  = 0.;
    if(nGray>=15.) nBlack = 28.;
    Float_t nBlackNeutrons = nBlack * 0.84;
    Float_t nBlackProtons  = nBlack - nBlackNeutrons;

// Actual number (including fluctuations) from binomial distribution    
//  gray neutrons
    ngn = (double) gRandom->Binomial((Int_t) fN, nGrayNeutrons/fN);
    
//  gray protons
    ngp = (double) gRandom->Binomial((Int_t) fP, nGrayProtons/fP); 

//  black neutrons
    nbn = (double) gRandom->Binomial((Int_t) fN, nBlackNeutrons/fN);
    
//  black protons
    nbp = (double) gRandom->Binomial((Int_t) fP, nBlackProtons/fP);
    
   // printf(" Ncoll %d  ngp %f nbp %f  ngn %f nbn %f \n",ncoll,ngp, nbp, ngn, nbn);

}


//--------------------------------------------------------------------------------------------------
void AliCentralityGlauberFit::MakeSlowNucleons2(Int_t ncoll, Double_t bog, Double_t gamma, 
 	Double_t &nbn, Double_t &ngn, Double_t &nbp, Double_t &ngp, Double_t &lcp)
{
// from AliGenSlowNucleonModelExp (Chiara Oppedisano)
//
// Return the number of black and gray nucleons
//

 // PROTONS ----------------------------------------------------------------

  Int_t fP = 82;
  Int_t fN = 126;
  
  Float_t nu = (Float_t) ncoll; 
  //
  nu = gRandom->Gaus(nu,0.5); 
  if(nu<1.) nu = 1.;

  // PROTONS ----------------------------------------------------------------
  //  gray protons
  Float_t  poverpd = 0.843; 
  Float_t  zAu2zPb = 82./79.;
  Float_t  nGrayp = (-0.27 + 0.63 * nu - 0.0008 *nu *nu)*poverpd*zAu2zPb;
  if(nGrayp<0.) nGrayp=0;

  Double_t p;
  p =  nGrayp/fP;
  ngp = (double) gRandom->Binomial(fP, p);

  //  black protons
  //Float_t blackovergray = 3./7.;// from spallation
  //Float_t blackovergray = 0.65; // from COSY
  Float_t blackovergray = bog;
  Float_t nBlackp  = blackovergray*nGrayp; 
  if(nBlackp<0.) nBlackp=0;

  p =  nBlackp/(fP-ngp);
  nbp = (double) gRandom->Binomial((Int_t) (fP-ngp), p);
  // SATURATION in no. of black protons
  //if(ngp>=9.) nbp = 15;
  
 // NEUTRONS ----------------------------------------------------------------

  /*if(nu<3.){
    nGrayp = -0.836 + 0.9112 *nu - 0.05381 *nu *nu;
    nBlackp  = blackovergray*nGrayp; 
  } */ 

  //  gray neutrons
  Float_t nGrayNeutrons = 0.;
  Float_t nBlackNeutrons = 0.;
  lcp = (nGrayp+nBlackp)/gamma;

  Float_t nSlow  = 0.;
  if(lcp>0.){
    nSlow  = fParamnSlow[0]-fParamnSlow[1]/(fParamnSlow[2]+lcp)+fParamnSlow[3]*lcp;
    //
//changed
//    float xconj = fParamnSlow[1]/fParamnSlow[0]-fParamnSlow[2];
    float xconj = 1.;
    float yconj = fParamnSlow[0]-fParamnSlow[1]/(fParamnSlow[2]+xconj)+fParamnSlow[3]*xconj;
    if(lcp<xconj) nSlow = 0.+(yconj-0.)*(lcp-0.)/(xconj-0.);
    //
    // Factor to scale COSY to obtain ALICE n mltiplicities
    //    nSlow = 1.68*nSlow;  // now integrated in fParamnSlow[0], fParamnSlow[1]
    //
    nBlackNeutrons = nSlow * 0.9;
    //nBlackNeutrons = nSlow * 0.5;
    nGrayNeutrons  = nSlow - nBlackNeutrons;
  }
  else{
    // slightly adjusted Sikler corrected for neutron yield...
    nBlackNeutrons = 126./208. * 4. * nu; 
    //nGrayNeutrons  = 126./208. * 2. * nu; // a la Sikler (<Ngn>~0.50*<Nbn>)
    //nGrayNeutrons  = 0.47 * 2. * nu; // a la DPMJET (<Ngn>~0.39*<Nbn>)
    nGrayNeutrons = nBlackNeutrons/9.;   // a la spallation (<Ngn>~0.11*<Nbn>)
  }
  //
  if(nGrayNeutrons<0.) nGrayNeutrons=0.;
  if(nBlackNeutrons<0.) nBlackNeutrons=0.;

  p =  nGrayNeutrons/fN;
  ngn = gRandom->Binomial(fN, p);
//changed
  if(nGrayNeutrons>=10) ngn = gRandom->Gaus(nGrayNeutrons+0.5, TMath::Sqrt(fN*p*(1-p)));
  //else ngn = gRandom->PoissonD(nGrayNeutrons);
  //else ngn = gRandom->Binomial((Int_t) fN, p);
  if(ngn<0.) ngn=0.;

  //  black neutrons
  p =  nBlackNeutrons/(fN-ngn);
  nbn = gRandom->Binomial((Int_t) (fN-ngn), p);
//changed
  if(nBlackNeutrons>=10) nbn = gRandom->Gaus(nBlackNeutrons+0.5, TMath::Sqrt(fN*p*(1-p)));
  else nbn = gRandom->PoissonD(nBlackNeutrons);
  //else nbn = gRandom->Binomial((Int_t) fN, p);
  if(nbn<0.) nbn=0.;
  
}

//--------------------------------------------------------------------------------------------------
 void AliCentralityGlauberFit::MakeSlowNucleons2s(Float_t ncoll, Double_t bog, Double_t gamma, Double_t sigmap, 
 	Double_t &nbn, Double_t &ngn, Double_t &nbp, Double_t &ngp, Double_t &lcp)
{
// from AliGenSlowNucleonModelExp (Chiara Oppedisano)
//
// Return the number of black and gray nucleons
//

  float fP = 82;
  float fN = 126;
  
  Float_t nu = ncoll;
  //Float_t nu = ncoll*(gRandom->Rndm()+0.5);

 // PROTONS ----------------------------------------------------------------
  //  gray protons
  Float_t  poverpd = 0.843; 
  Float_t  zAu2zPb = 82./79.;
  Float_t  grayp = (-0.27 + 0.63 * nu - 0.0008 *nu *nu)*poverpd*zAu2zPb;
  //if(nu<3.) grayp = -0.836 + 0.9112 *nu - 0.05381 *nu *nu;
  //
//changed
  //Float_t  nGrayp = grayp;
  Float_t  nGrayp = gRandom->Gaus(grayp, sigmap);
//changed
//  Float_t  nGrayp = gRandom->Gaus(grayp, (0.8*sigmap-sigmap)*grayp/90+sigmap);
  if(nGrayp<0.) nGrayp=0;

  Double_t p=0.;
  p =  nGrayp/fP;
  ngp = (double) gRandom->Binomial((int) fP, p);

  //  black protons
  //Float_t blackovergray = 3./7.;// =0.43 from spallation
  //Float_t blackovergray = 0.65; // from COSY
  Float_t blackovergray = bog;
  //
  Float_t nBlackp = blackovergray*nGrayp;
  //if(nBlackp<0.) nBlackp=0;

  p =  nBlackp/(fP-ngp);
  nbp = (double) gRandom->Binomial((int) (fP-ngp), p);
  // SATURATION in no. of black protons
  //if(ngp>=9.) nbp = 15;

 // NEUTRONS ----------------------------------------------------------------
  // Relevant ONLY for n
  // NB-> DOESN'T WORK FOR SMEARING IN Nslow!!!!!!!!!!
  /*if(nu<=3.){
    grayp = -0.836 + 0.9112 *nu - 0.05381 *nu *nu;
    // smearing
    //nGrayp = gRandom->Gaus(grayp, sigmap);
    nGrayp = grayp;
    nBlackp  = blackovergray*nGrayp; 
  }*/
    
  //  gray neutrons
  lcp = (nGrayp+nBlackp)/gamma;
  
  Float_t nGrayNeutrons = 0.;
  Float_t nBlackNeutrons = 0.;
  Float_t nSlow  = 0.;
  if(lcp>0.){
    nSlow  = fParamnSlow[0]-fParamnSlow[1]/(fParamnSlow[2]+lcp)+fParamnSlow[3]*lcp;
    //
//changed
//    float xconj = fParamnSlow[1]/fParamnSlow[0]-fParamnSlow[2];
    /*float xconj = 1.;
    float yconj = fParamnSlow[0]-fParamnSlow[1]/(fParamnSlow[2]+xconj)+fParamnSlow[3]*xconj;
    if(lcp<xconj) nSlow = 0.+(yconj-0.)*(lcp-0.)/(xconj-0.);*/
    //
    // Factor to scale COSY to obtain ALICE n mltiplicities
    //    nSlow = 1.68*nSlow;  // now integrated in fParamnSlow[0], fParamnSlow[1]
    //
    nBlackNeutrons = nSlow * 0.9;
    //nBlackNeutrons = nSlow * 0.5;
    nGrayNeutrons  = nSlow - nBlackNeutrons;
  }
  else{
   // taking into account the prob4bility p to have 0 neutrons when we have 0 protons  
   // p = 1.-(0.82*0.96) ~ 0.22 (assuming that 100% of ev. without n are also without p)
   //Float_t r = gRandom->Rndm();
   //if(r>0.1){
     // Sikler corrected for neutron yield...
     nBlackNeutrons = 126./208. * 4. * nu; 
     //nGrayNeutrons  = 126./208. * 2. * nu; // a la Sikler (<Ngn>~0.50*<Nbn>)
     //nGrayNeutrons  = 0.47 * 2. * nu; // a la DPMJET (<Ngn>~0.39*<Nbn>)
     //nGrayNeutrons = nBlackNeutrons/9.;   // a la spallation (<Ngn>~0.11*<Nbn>)
     // smearing!
//changed
     nBlackNeutrons = gRandom->Gaus(nBlackNeutrons, sigmap);
     //nGrayNeutrons  = nBlackNeutrons/2.; // a la Sikler (<Ngn>~0.50*<Nbn>)
     nGrayNeutrons = nBlackNeutrons/9.;   // a la spallation (<Ngn>~0.11*<Nbn>)
     //printf("  Sikler -> Ncoll %f <Ngrayn> %f  <Nblackn> %f \n",nu,nGrayNeutrons, nBlackNeutrons);
   //}
  }
  //
  if(nGrayNeutrons<0.) nGrayNeutrons=0.;
  if(nBlackNeutrons<0.) nBlackNeutrons=0.;

  //  black neutrons
  p =  nBlackNeutrons/fN;
  nbn = gRandom->Binomial((Int_t) fN, p);
  TRandom1 r;
  //nbn = (double) r.Binomial((int) fN, p);
//changed
//  if(nBlackNeutrons>=10) nbn = r.Gaus(nBlackNeutrons+0.5, TMath::Sqrt(fN*p*(1-p)));
  if(nbn<0.) nbn=0.;

  //  gray neutrons
  p =  nGrayNeutrons/(fN-nbn);
  ngn = gRandom->Binomial((Int_t) (fN-nbn), p);
  //ngn = (double) (r.Binomial((Int_t) (fN-nbn), p));
//changed
  /*if(nGrayNeutrons>=10) ngn = gRandom->Gaus(nGrayNeutrons+0.5, TMath::Sqrt(fN*p*(1-p)));
  else ngn = gRandom->Binomial((Int_t) fN, p);*/
  if(ngn<0.) ngn=0.;

  //printf(" Ncoll %1.2f   Ngp %1.2f Nbp %1.2f   Nbn %1.2f Ngn %1.2f\n",nu, ngp,nbp,nbn,ngn);

}


Double_t AliCentralityGlauberFit::ConvertToEnergy(Double_t T)
{ 
   TDatabasePDG * pdg = TDatabasePDG::Instance();
    const Double_t kMassNeutron = pdg->GetParticle(kNeutron)->Mass();
    Double_t fPmax =10;   

    Double_t m, p, pmax, f;
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

