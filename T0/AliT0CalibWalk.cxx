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
  Int_t npeaks = 20;
  Int_t sigma=3.;
  Bool_t down=false;
  Int_t index[20];
  Bool_t ok=true;
  Float_t   mips[20];
  Int_t nfound=0;
  
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
      
      if (nmips<15) ok=false; 
       
      Float_t x1[50], y1[50]; 
      Float_t x2[50], xx2[50],y2[50];
      Float_t xx1[50],yy1[50], xx[50];
      
      Float_t cfd0[24];
      Int_t startim = 0;
      
      for (Int_t ii=0; ii<nmips; ii++)
	x1[ii] = y1[ii] = x2[ii] = y2[ii] = 0; 
      
      for (Int_t i=0; i<24; i++)
	{
	  for (Int_t im=startim; im<nmips; im++)
	    {	      
	      TString cfd = Form("hCFD%i_%i",i+1,im+1);
	      TString qtc = Form("hQTC%i_%i",i+1,im+1);
	      TString led = Form("hLED%i_%i",i+1,im+1);
	      
	      TH1F *hCFD = (TH1F*) gFile->Get(cfd.Data()) ;
	      TH1F *hLED = (TH1F*) gFile->Get(led.Data());
	      TH1F *hQTC = (TH1F*) gFile->Get(qtc.Data()) ;
	      hCFD->SetDirectory(0);
	      hQTC->SetDirectory(0);
	      hLED->SetDirectory(0);
	      if(!hCFD )
 	      	AliWarning(Form(" no CFD data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
	      if(!hQTC )
 	      	AliWarning(Form(" no QTC correction data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
	      if(!hLED)	      
	      	AliWarning(Form(" no LED correction data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
      	      
	      if(hCFD )	{
		TSpectrum *s = new TSpectrum(2*npeaks,1);
		nfound = s->Search(hCFD,sigma," ",0.1);
		if(nfound!=0){
		  Float_t *xpeak = s->GetPositionX();
		  TMath::Sort(nfound, xpeak, index,down);
		  Float_t xp = xpeak[index[0]];
		  Double_t hmax = xp+3*sigma;
		  Double_t hmin = xp-3*sigma;
		  hCFD->GetXaxis()->SetRangeUser(hmin-10,hmax+10);
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
		      ok=false;
		  }

		if (im == 0) cfd0[i] = hCFD->GetMean();
		y1[im] =  hCFD->GetMean() - cfd0[i];
	      }
	      if( hQTC) x1[im] = hQTC->GetMean();
	      
	      if( hLED){
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
		  ok=false;
		x2[im] = hLED->GetMean();
		xx2[im] = x2[nmips-im-1];
	      }
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
	  
	  
	  /*
	  
	  
	  Float_t x1[50], y1[50]; 
	  Float_t x2[50], xx2[50],y2[50];
	  Float_t xx1[50],yy1[50], xx[50];
	  
	  Int_t nmips=20;
	  for (Int_t i=0; i<24; i++)
	  {
 
	  for (Int_t im=0; im<nmips; im++)
	  {
	  x1[im]=xx[im]=500+im*200;
	  y1[im]=0;
	  x2[im]=xx1[im]=yy1[im]=260+20*im;
	  } 
	  */	  
      Bool_t ok=true;
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
  ok=true;
  return ok;
}




