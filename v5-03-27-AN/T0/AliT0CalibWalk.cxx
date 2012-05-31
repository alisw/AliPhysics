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
  Float_t sigma,cfdmean, qtmean, ledmean;
  Bool_t ok=true;
  Float_t   mips[50];
  
  gFile = TFile::Open(laserFile);
  if(!gFile) {
    AliError("No input laser data found ");
  }
  else
    {
      //      gFile->ls();
      TH1F* hAmp = (TH1F*) gFile->Get("hAmpLaser");
      Int_t nmips=0;
      for (Int_t ibin=0; ibin<2000; ibin++) {
	Float_t bincont = hAmp->GetBinContent(ibin);
	if(bincont>0){ 
	  mips[nmips] = hAmp->GetXaxis()->GetBinCenter(ibin);
	  cout<<ibin<<" bincont "<<bincont<<" amp "<< mips[nmips]<<endl;
	  nmips++;
	}	
      }    
      
      if (nmips<17) {
	ok=false;
	return ok;
      } 
       
      Float_t x1[50], y1[50]; 
      Float_t x2[50], xx2[50],y2[50];
      Float_t xx1[50],yy1[50], xx[50];
      
      Float_t cfd0 = 0;
      
      for (Int_t ii=0; ii<nmips; ii++)
	x1[ii] = y1[ii] = x2[ii] = y2[ii] = 0; 
      
      for (Int_t i=0; i<24; i++)
	{
	  cfd0 = 0;
	  for (Int_t im=0; im<nmips; im++)
	    {	      
	      TString cfd = Form("hCFD%i_%i",i+1,im+1);
	      TString qtc = Form("hQTC%i_%i",i+1,im+1);
	      TString led = Form("hLED%i_%i",i+1,im+1);
	      
	      TH1F *hCFD = (TH1F*) gFile->Get(cfd.Data()) ;
	      TH1F *hLED = (TH1F*) gFile->Get(led.Data());
	      TH1F *hQTC = (TH1F*) gFile->Get(qtc.Data()) ;
	      //	      hCFD->SetDirectory(0);
	      //	      hQTC->SetDirectory(0);
	      //	      hLED->SetDirectory(0);
	      if(!hCFD )
 	      	AliWarning(Form(" no CFD data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
	      if(!hQTC )
 	      	AliWarning(Form(" no QTC correction data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
	      if(!hLED)	      
	      	AliWarning(Form(" no LED correction data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
	      if( hCFD && hCFD->GetEntries()<500 ) {
		ok=false;
		printf("no peak in CFD spectrum for PMT %i amplitude %i\n",i,im);
		return ok;
	      }
	      if(hCFD && hCFD->GetEntries()>500 ) {
		if( hCFD->GetRMS() >= 1.5) 
		  GetMeanAndSigma(hCFD, cfdmean, sigma);
		else
		  cfdmean = hCFD->GetMean();
		
		Int_t   maxBin = hCFD->GetMaximumBin(); 
		Double_t  meanEstimate = hCFD->GetBinCenter( maxBin); 
		if(TMath::Abs(meanEstimate - cfdmean) > 20 ) cfdmean = meanEstimate; 
		if (im == 0) cfd0 = cfdmean;
		y1[im] =  cfdmean - cfd0;
	      }	
	      if(hQTC && hQTC->GetEntries()>500) {
		GetMeanAndSigma(hQTC, qtmean, sigma);
		
		x1[im] = qtmean;
		if( x1[im] == 0) {
		  ok=false;
		  printf("no peak in QTC signal for PMT %i amplitude %i\n",i,im);
		  return ok;
		}
	      }
		
	      if( hLED && hLED->GetEntries()>500) {
		GetMeanAndSigma(hLED, ledmean, sigma);
	      }				   
	      else
		{ 
		  ok=false;
		  printf("no peak in LED spectrum for PMT %i amplitude %i\n",i,im);
		  return ok;
		}
	      x2[im] = ledmean;
	      xx2[im] = x2[nmips-im-1];
		
		xx[im]=mips[im];
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
	  
      TGraph *grwalkqtc = new TGraph (nmips,x1,y1);
      grwalkqtc->SetTitle(Form("PMT%i",i));
      TGraph *grwalkled = new TGraph (nmips,x2,y1);
      grwalkled->SetTitle(Form("PMT%i",i));
      fWalk.AddAtAndExpand(grwalkqtc,i);
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



