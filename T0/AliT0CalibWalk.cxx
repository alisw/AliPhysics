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
 
/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class T0 walk correction      Alla Maevskaya alla@inr.ru 20.11.2007        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliT0CalibWalk.h"
#include "AliLog.h"

#include <TObjArray.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH2F.h> 
#include <TMath.h>
#include <TSystem.h>
#include <Riostream.h>
#include <TSpectrum.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TF1.h>

ClassImp(AliT0CalibWalk)

//________________________________________________________________
  AliT0CalibWalk::AliT0CalibWalk():   TNamed(),
				      fWalk(0),
				      fAmpLEDRec(0),
				      fQTC(0),
				      fAmpLED(0)
{
  //
}

//________________________________________________________________
AliT0CalibWalk::AliT0CalibWalk(const char* name):TNamed(),
				      fWalk(0),
				      fAmpLEDRec(0),				  				      fQTC(0),
				      fAmpLED(0)
    
{
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());

}

//________________________________________________________________
AliT0CalibWalk::AliT0CalibWalk(const AliT0CalibWalk& calibda) :
  TNamed(calibda),		
  fWalk(0),
  fAmpLEDRec(0),
  fQTC(0),
  fAmpLED(0)
 
{
// copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());


}

//________________________________________________________________
AliT0CalibWalk &AliT0CalibWalk::operator =(const AliT0CalibWalk& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
 
  return *this;
}

//________________________________________________________________
AliT0CalibWalk::~AliT0CalibWalk()
{
  //
}
//________________________________________________________________  

Bool_t AliT0CalibWalk::MakeWalkCorrGraph(const char *laserFile)
{
  //make walk corerction for preprocessor
  Int_t npeaks = 20;
  Int_t sigma=3;
  Bool_t down=false;
  Int_t index[20];
  Bool_t ok=true;
  Float_t   mips[100];
  Int_t nfound=0;
  /*  
  //061110 
 Float_t timefitth[24] = { 5594, 5618, 5630, 5662,  5623, 5645,
			    5645, 5602, 5691, 5669,  5636, 5602,
			    5651, 5654, 5612, 5625,  5669, 5678,
			    5652, 5673, 5565, 5717,  5692, 5738} ;
  */
 
  //051210
  Float_t timefitth[24] = { 5587, 5609, 5621, 5653, 5615, 5635, 
			    5637, 5595, 5683, 5659, 5629, 5601,
			    5643, 5654, 5601, 5615, 5661, 5667, 
			    5645, 5665, 5663, 5709, 5685, 5725 };
 
  Float_t xdata[4000], ydata[4000];
  Int_t npoints =0;
  TH1F* www[24];
  TFile *ff = new TFile("calibwalk0.root");

  Int_t ipmt, jpoint; 
  
  gFile = TFile::Open(laserFile);
  if(!gFile) {
    AliError("No input laser data found ");
  }
  else
    {
      //      gFile->ls();
      TH1F* hAmp = (TH1F*) gFile->Get("hAmpLaser");
      Int_t nmips=0;
      for (Int_t ibin=0; ibin<1200; ibin++) { 
	Float_t bincont = hAmp->GetBinContent(ibin);
	if(bincont>0){ 
	  mips[nmips] = hAmp->GetXaxis()->GetBinCenter(ibin);
	  cout<<ibin<<" bincont "<<bincont<<" amp "<< mips[nmips]<<endl;
	  nmips++;
	}	
      }    
      
      /*      if (nmips<17) {
	ok=false;
	return ok;
	} */
      //  TCanvas *c1 = new TCanvas("c1", "LED-CFD C side",0,48,1280,951);
      // c1->Divide(4,3);
   //   TCanvas *c2 = new TCanvas("c1", "LED-CFD C side",0,48,1280,951);
   //   c2->Divide(4,3);
        
      Float_t x1[50], y1[50]; 
      Float_t x2[50], xx2[50],y2[50];
      Float_t xx1[50],yy1[50], xx[50];
      
      Float_t cfd0[24];
      Int_t startim = 0;
      
      
      for (Int_t i=0; i<24; i++)
	{
	  //	 if(i==13) continue;
	  //	  else c2->cd(i+1-12);
	  for (Int_t im=startim; im<nmips; im++)
	    {	      

	      x1[im] = y1[im] = x2[im] = y2[im] = 0; 
	      if(i==23 && im<2) continue;
	      xx[im]=mips[im];
	      TString cfd = Form("hCFD%i_%i",i+1,im+1);
	      TString qtc = Form("hQTC%i_%i",i+1,im+1);
	      TString led = Form("hLED%i_%i",i+1,im+1);
	      
	      TH1F *hCFD = (TH1F*) gFile->Get(cfd.Data()) ;
	      TH1F *hLED = (TH1F*) gFile->Get(led.Data());
	      TH1F *hQTC = (TH1F*) gFile->Get(qtc.Data()) ;
	      if(!hCFD )
 	      	AliWarning(Form(" no CFD data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
	      if(!hQTC )
 	      	AliWarning(Form(" no QTC correction data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
	      if(!hLED)	      
	      	AliWarning(Form(" no LED correction data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
      	      
	      if(hCFD )	{
		hCFD->SetDirectory(0);
		hCFD->GetXaxis()->SetRangeUser(timefitth[i]-100,timefitth[i]+100);

		TF1 *fg = new TF1("fg","gaus",timefitth[i]-50,timefitth[i]+50);
		hCFD->Fit("fg","RQ", " ",timefitth[i]-50,timefitth[i]+50);
		Double_t par[3];
		fg->GetParameters(&par[0]);

		TSpectrum *s = new TSpectrum(2*npeaks,1);
		nfound = s->Search(hCFD,sigma," ",0.1);
		if(nfound!=0){
		  Float_t *xpeak = s->GetPositionX();
		  TMath::Sort(nfound, xpeak, index,down);
		  Float_t xp = xpeak[index[0]];
		  Double_t hmax = xp+3*sigma;
		  Double_t hmin = xp-3*sigma;
		  hCFD->GetXaxis()->SetRangeUser(hmin-10,hmax+10);
		  //	if (i==23) 
		  //		  cout<<"PMT "<<i<<" MIPS "<<mips[im]<<" mean "<<hCFD->GetMean()<<
		  //     " xp "<<xp<<" fit "<<par[1]<<endl;
		  //		  cout<<"PMT "<<i<<" MIPS "<<im<<" xp "<<xp<<endl;
		}
		else
		  {
		    hCFD->Rebin(2);
		    TSpectrum *s1 = new TSpectrum(2*npeaks,1);
		    nfound = s1->Search(hCFD,sigma," ",0.1);
		    if(nfound!=0){
		      Float_t *xpeak = s1->GetPositionX();
		      TMath::Sort(nfound, xpeak, index,down);
		      Float_t xp = xpeak[index[0]];
		      Double_t hmax = xp+3*sigma;
		      Double_t hmin = xp-3*sigma;
		      hCFD->GetXaxis()->SetRangeUser(hmin-10,hmax+10);
		    }
		    else 
		      {
			ok=true;
			//			ok=false;
			printf("no peak in CFD spectrum for PMT %i amplitude %i\n",i,im);
			return ok;
		      }
		
		  }
		// } //if no C12 80MIPs
		//	  if (i<12) cout<<"!!!! PMT "<<i<<" MIPS "<<im<<" fit  "<<par[1]<<
		//" mean "<<hCFD->GetMean()<<endl;
		 if (im == 0) cfd0[i] = par[1];
		if(i==23 && im==2) cfd0[i] = par[1]; 
		   y1[im] =  par[1] - cfd0[i];
		   if ((i == 5 && im == 15) ||(i == 8 && im == 19) || (i==11 && im==0) || i==23)
		       y1[im] =  hCFD->GetMean() - cfd0[i];
		  
		//		if (im == 0) cfd0[i] = hCFD->GetMean();
		//		y1[im] =  hCFD->GetMean() - cfd0[i];
		  
		   //	if (i==23 && im<2) y1[im]=0;
		/*		if( i == 11 && im==10 ){	y1[im] = y1[im-1];
		cout<<im<<" "<<mips[im]<<" mip "<<" pmt "<<i<<" mean "<<hCFD->GetMean()<<" walk "<<	y1[im]<<endl;
		}*/ 

	      }
	      if( hQTC) {
		x1[im] = hQTC->GetMean();
		if( x1[im] == 0) {
		  ok=false;
		  printf("no peak in QTC signal for PMT %i amplitude %i\n",i,im);
		  return ok;
		  }
	      }
	      
	      if( hLED){
		hLED->SetDirectory(0);
		TSpectrum *s = new TSpectrum(2*npeaks,1);
		nfound = s->Search(hLED,sigma," ",0.1);
		if(nfound!=0){
		  Float_t *xpeak = s->GetPositionX();
		  TMath::Sort(nfound, xpeak, index,down);
		  Float_t xp = xpeak[index[0]];
		  Double_t hmax = xp+10*sigma;
		  Double_t hmin = xp-10*sigma;
		  hLED->GetXaxis()->SetRangeUser(hmin-10,hmax+10);
		}
		else
		  { 
		    ok=false;
		    printf("no peak in LED spectrum for PMT %i amplitude %i\n",i,im);
		    return ok;
		  }
		x2[im] = hLED->GetMean();
		xx2[im] = x2[nmips-im-1];
	      }
	      
	      if (hQTC) delete  hQTC;
	      if (hCFD) delete  hCFD;
	      if (hLED) delete  hLED;
	    }
	  
	  for (Int_t imi=0; imi<nmips; imi++)
	    {
	      yy1[imi] = Float_t (mips[nmips-imi-1]);
	      xx1[imi] = x2[nmips-imi-1]; 
	      //	      cout<<" LED rev "<<i<<" "<<imi<<" "<<xx1[imi]<<" "<< yy1[imi]<<"nmips-imi-1 = "<< nmips-imi-1<<endl;
	    }
	  
	  if(i==0) cout<<"Making graphs..."<<endl;
	  // TGraph *grwalkqtc = new TGraph (nmips-4,x1,y1);
	  // grwalkqtc->SetTitle(Form("PMT%i",i));
	  TGraph *grwalkled = new TGraph (nmips-4,x2,y1);
	  grwalkled->SetTitle(Form("PMT%i",i));
	  // fWalk.AddAtAndExpand(grwalkqtc,i);
	  fAmpLEDRec.AddAtAndExpand(grwalkled,i);
	  //	  cout<<" add walk "<<i<<endl;
	  
	  //fit amplitude graphs to make comparison wth new one	  
	  TGraph *grampled = new TGraph (nmips,xx1,yy1);
	  TGraph *grqtc = new TGraph (nmips,x1,xx);
	  fQTC.AddAtAndExpand(grqtc,i);	 
	  fAmpLED.AddAtAndExpand(grampled,i);
	  //	  cout<<" add amp "<<i<<endl;
	  
	  if(i==23)
	    cout<<"Graphs created..."<<endl;
	  //	   if(i==13) continue;
	  
	  www[i] = (TH1F*)ff->Get(Form("hy%i",i) );
	  www[i]->Print();
	  npoints =0;
	  for(Int_t ibin=1; ibin< (www[i]->GetNbinsX() -2); ibin++) {
	     if( www[i]->GetBinContent(ibin)>0 && 
		www[i]->GetBinContent(ibin) !=1000 ) {
	      ydata[ibin-1] = www[i]->GetBinContent(ibin) - 1000;
	      xdata[ibin-1] = ibin+799;
	      cout<<ibin<<" "<<xdata[ibin]<<" "<<ydata[ibin]<<endl;
	      npoints++;
	     }
	  } 


	  TGraph *grwalkqtc = new TGraph (npoints-3,xdata,ydata);
	  grwalkqtc->SetTitle(Form("PMT%i",i));
	  fWalk.AddAtAndExpand(grwalkqtc,i);
	}
	  
      
    } //if gFile exits

  return ok;
}
void AliT0CalibWalk::GetMeanAndSigma(TH1F* hist, Float_t &mean, Float_t &sigma) 
{

  const double window = 2.;  //fit window 
 
  double meanEstimate, sigmaEstimate; 
  int maxBin;
  maxBin        =  hist->GetMaximumBin(); //position of maximum
  meanEstimate  =  hist->GetBinCenter( maxBin); // mean of gaussian sitting in maximum
  sigmaEstimate = hist->GetRMS();
  TF1* fit= new TF1("fit","gaus", meanEstimate - window*sigmaEstimate, meanEstimate + window*sigmaEstimate);
  fit->SetParameters(hist->GetBinContent(maxBin), meanEstimate, sigmaEstimate);
  hist->Fit("fit","RQ","Q");

  mean  = (Float_t) fit->GetParameter(1);
  sigma = (Float_t) fit->GetParameter(2);
  printf(" mean %f max %f sigma %f \n",mean, meanEstimate , sigma);  

  delete fit;
}



