#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TH1D.h"
#include "TH1.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#endif

/******************************************************
 *
 *  Macro to compute fprompt for pPb and PbPb data
 *    using the outputs of the HFPtSpectrumRaa maco
 *  Questions and complains to Z.Conesa del Valle
 *
 *******************************************************/

Bool_t elossFDQuadSum = true;

TGraphAsymmErrors * DrawFpromptVsRaaElossHypo(const char *infile="", //const char *outfile="",
					      Double_t MinHypo=1./3., Double_t MaxHypo=3.0, Double_t CentralHypo = 1.0, Bool_t isYieldUnc=true);


//_______________________________________________________________________
void DrawFpromptVsRaaElossHypoCombined(const char *infileNb="", const char *infilefc="", const char *outfile="", Double_t MinHypo=1./3., Double_t MaxHypo=3.0, Double_t CentralHypo = 1.0, Bool_t isYieldUnc=true)
{

  TH2F *hdraw = new TH2F("hdraw","fprompt vs pt; p_{T} (GeV/c); f_{prompt}",36,0.,36.,10.,0.,1);
  hdraw->Draw();

  cout<<endl<<endl<<" ********* Checking Nb method ********* "<<endl<<endl;
  TGraphAsymmErrors *gfPromptNb = DrawFpromptVsRaaElossHypo(infileNb,MinHypo,MaxHypo,CentralHypo,isYieldUnc);
  gfPromptNb->SetNameTitle("gfPromptNb","fprompt Nb method, Eloss range considered");
  gfPromptNb->SetLineColor(kBlack);
  gfPromptNb->SetMarkerColor(kGreen+2);
  gfPromptNb->SetMarkerStyle(25);
  gfPromptNb->Draw("P");
					     
  cout<<endl<<endl<<" ********* Checking fc method ********* "<<endl<<endl;
  TGraphAsymmErrors *gfPromptfc = DrawFpromptVsRaaElossHypo(infilefc,MinHypo,MaxHypo,CentralHypo,isYieldUnc);
  gfPromptfc->SetNameTitle("gfPromptfc","fprompt fc method, Eloss range considered");
  gfPromptfc->SetLineColor(kBlack);
  gfPromptfc->SetMarkerColor(kRed);
  gfPromptfc->SetMarkerStyle(21);
  gfPromptfc->Draw("P");

  cout<<endl<<endl<<" ********* Computing the envelope ******** "<<endl<<endl;
  TGraphAsymmErrors *gFPromptCombined = new TGraphAsymmErrors(0);
  gFPromptCombined->SetNameTitle("gFPromptCombined","fprompt Nb + fc, Eloss range considered");

  Int_t nbins=gfPromptNb->GetN();
  for(Int_t i=0; i<nbins; i++) {
    Double_t pt=0.,value=0.,ptfc=0.,valuefc=0.,exlow=0.,eylow=0.,eyNblow=0.,eyfclow=0.,eyhigh=0.,eyNbhigh=0.,eyfchigh=0.;
    gfPromptNb->GetPoint(i,pt,value);
    gfPromptfc->GetPoint(i,ptfc,valuefc);
    exlow = gfPromptNb->GetErrorXlow(i);
    eyNblow = gfPromptNb->GetErrorYlow(i);
    eyNbhigh = gfPromptNb->GetErrorYhigh(i);
    eyfclow = gfPromptfc->GetErrorYlow(i);
    eyfchigh = gfPromptfc->GetErrorYhigh(i);
    Double_t uncertainties[5]={ value, value-eyNblow, value+eyNbhigh, valuefc-eyfclow, valuefc+eyfchigh };
    eylow = value - TMath::MinElement(5,uncertainties);
    eyhigh = TMath::MaxElement(5,uncertainties) - value;
    gFPromptCombined->SetPoint(i,pt,value);
    gFPromptCombined->SetPointError(i,exlow,exlow,eylow,eyhigh);
  }
  
  gFPromptCombined->SetFillStyle(0);
  gFPromptCombined->SetLineColor(kBlack);
  gFPromptCombined->SetMarkerColor(kBlack);
  gFPromptCombined->SetMarkerStyle(20);
  gFPromptCombined->Draw("2P");

  cout<<endl<<endl;

  TFile* fout= new TFile(outfile,"recreate");
  gfPromptNb->Write();
  gfPromptfc->Write();
  gFPromptCombined->Write();
  fout->Write();
}

//_______________________________________________________________________
TGraphAsymmErrors * DrawFpromptVsRaaElossHypo(const char *infile, 
					      //const char *outfile="",
					      Double_t MinHypo, Double_t MaxHypo, Double_t CentralHypo, Bool_t isYieldUnc)
{

  Float_t pt=0.,TAB=0.,sigmaPP=0.,invyieldAB=0.,RaaCharm=0.,RaaBeauty=0., fc=0.;
  Float_t yieldFdHigh=0., yieldFdLow=0.,RaaCharmFdHigh=0.,RaaCharmFdLow=0.;

  TFile *fRaa = new TFile(infile,"read");
  TH1D *hRABvsPt = (TH1D*)fRaa->Get("hRABvsPt");
  TNtuple *ntRaa = (TNtuple*)fRaa->Get("ntupleRAB");
  ntRaa->SetBranchAddress("pt",&pt);
  ntRaa->SetBranchAddress("TAB",&TAB);
  ntRaa->SetBranchAddress("sigmaPP",&sigmaPP);
  ntRaa->SetBranchAddress("invyieldAB",&invyieldAB);
  ntRaa->SetBranchAddress("RABCharm",&RaaCharm);
  ntRaa->SetBranchAddress("RABBeauty",&RaaBeauty);
  ntRaa->SetBranchAddress("RABCharmFDHigh",&RaaCharmFdHigh);
  ntRaa->SetBranchAddress("RABCharmFDLow",&RaaCharmFdLow);
  ntRaa->SetBranchAddress("invyieldABFDHigh",&yieldFdHigh);
  ntRaa->SetBranchAddress("invyieldABFDLow",&yieldFdLow);
  ntRaa->SetBranchAddress("fc",&fc);

  // define the binning
  Int_t nbins = hRABvsPt->GetNbinsX();

  //
  //
  // Search the central value of the energy loss hypothesis Rb = Rc (bin)
  //
  Double_t ElossCentral[nbins+1];
  Double_t fcV[nbins+1], fcVmax[nbins+1], fcVmin[nbins+1];
  Double_t fdMax[nbins+1], fdMin[nbins+1];
  for(Int_t i=0; i<=nbins; i++) { 
    ElossCentral[i]=0.; 
    fcV[i]=0.; fcVmax[i]=0.; fcVmin[i]=6.;
    fdMax[i]=0.; fdMin[i]=6.;
  }
  //
  Int_t netries = ntRaa->GetEntries();
  //
  for(Int_t ientry=0; ientry<=netries; ientry++){
    ntRaa->GetEntry(ientry);
    Double_t ElossHypo =  1. / (RaaCharm / RaaBeauty) ;
    //
    // Find the bin for the central Eloss hypo
    //
    if( TMath::Abs( ElossHypo - CentralHypo ) < 0.075 ){
      Int_t hABbin = hRABvsPt->FindBin( pt );
      Double_t DeltaIni = TMath::Abs( ElossCentral[ hABbin ] - CentralHypo );
      Double_t DeltaV = TMath::Abs( ElossHypo - CentralHypo );
      //      cout << " pt " << ptTot << " ECentral Tot " << ElossCentralTot[ hABbin ] << " Ehypo "<< ElossHypo ;
      if ( DeltaV < DeltaIni ) ElossCentral[ hABbin ] = ElossHypo;
      //      cout << " final ECentral " << ElossCentralTot[ hABbin ] << endl;
    }
  }

  //
  for(Int_t ientry=0; ientry<=netries; ientry++){

    ntRaa->GetEntry(ientry);
    Double_t ElossHypo =  1. / (RaaCharm / RaaBeauty) ;
    Int_t hABbin = hRABvsPt->FindBin( pt );

    // Get the central value of Raa
    if ( ElossHypo == ElossCentral[ hABbin ] ) {

      fcV[hABbin] = fc;

      // Now look at FONLL band for the central value
      Double_t elems[3]= {0., 0., 0.};
      // 	fdMax[hABbin]= (RaaCharmFdHigh-RaaCharm)/RaaCharm;
      // 	fdMin[hABbin]= (RaaCharm-RaaCharmFdLow)/RaaCharm;
      if(isYieldUnc) {
	elems[0]=invyieldAB;
	elems[1]=yieldFdHigh;
	elems[2]=yieldFdLow;
      } else {
	elems[0]=RaaCharm;
	elems[1]=RaaCharmFdHigh;
	elems[2]=RaaCharmFdLow;
      }
      fdMax[hABbin]= (TMath::MaxElement(3,elems)-elems[0])/elems[0];
      fdMin[hABbin]= (elems[0]-TMath::MinElement(3,elems))/elems[0];
      //      cout<<" Check FONLL band giving for pt="<< pt <<" +"<< fdMax[hABbin]<< " -"<< fdMin[hABbin]<<" (relative variation)"<<endl;
    }

    // Here check the Eloss band
    if( RaaCharm>0 && ElossHypo >= MinHypo && ElossHypo <=MaxHypo ) {
      if( fc < fcVmin[hABbin] ) fcVmin[hABbin] = fc;
      if( fc > fcVmax[hABbin] ) fcVmax[hABbin] = fc;
    }
  }

  TGraphAsymmErrors *gFPrompt = new TGraphAsymmErrors(0);
  Int_t j=0;
  for(Int_t i=1; i<=nbins; i++) {
    Double_t xpt = hRABvsPt->GetBinCenter(i);
    Double_t width  = hRABvsPt->GetBinWidth(i)/2;
    Double_t elossLow = (fcV[i]-fcVmin[i])/fcV[i];
    Double_t elossHigh = (fcVmax[i]-fcV[i])/fcV[i];
    Double_t fcLow = fcV[i]-fcVmin[i];
    Double_t fcHigh = fcVmax[i]-fcV[i];
    //    cout<<" pt="<<xpt<<" , fprompt="<<fcV[i]<<" ["<< fcVmin[i] <<","<<fcVmax[i] <<"] ( Eloss) "<<endl;

    if(elossFDQuadSum){ // add quadratically
      fcLow = TMath::Sqrt( fdMin[i]*fdMin[i] + elossLow*elossLow )*fcV[i];
      fcHigh = TMath::Sqrt( fdMax[i]*fdMax[i] + elossHigh*elossHigh )*fcV[i];
    } else { // add linearly
      fcLow = ( fdMin[i] + elossLow )*fcV[i];
      fcHigh = ( fdMax[i] + elossHigh )*fcV[i];
    }

    gFPrompt->SetPoint(j,xpt,fcV[i]);
    //    gFPrompt->SetPointError(j,width,width,fcV[i]-fcVmin[i],fcVmax[i]-fcV[i]);
    gFPrompt->SetPointError(j,width,width,fcLow,fcHigh);
    j++;

    cout<<" pt="<<xpt<<" , fprompt="<<fcV[i]<<" +"<<fcHigh<<" -"<<fcLow<<" ( splitted it is +"<<elossLow*fcV[i]<<" -"<<elossHigh*fcV[i]<<" Eloss, +"<<fdMax[i]*fcV[i]<<" -"<<fdMin[i]*fcV[i]<<" FONLL-band )"<<endl;
  }
  //  gFPrompt->Draw("AP");

  // if(strcmp(outfile,"")!=0) {
  //   TFile *out = new TFile(outfile,"recreate");
  //   gFPrompt->Write();
  // }

  fRaa->Close();

  return gFPrompt;
  
}
