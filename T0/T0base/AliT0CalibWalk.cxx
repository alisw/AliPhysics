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

using std::cout;
using std::endl;
using std::ifstream;
ClassImp(AliT0CalibWalk)

//________________________________________________________________
  AliT0CalibWalk::AliT0CalibWalk():   TNamed(),
				      fWalk(0),
				      fAmpLEDRec(0),
				      fQTC(0),
                                      fAmpLED(0), 
                                      fCalibByData(kFALSE)
{
  //
}

//________________________________________________________________
AliT0CalibWalk::AliT0CalibWalk(const char* name):TNamed(),
				      fWalk(0),
				      fAmpLEDRec(0),				  				      fQTC(0),
				      fAmpLED(0), 
                                      fCalibByData(kFALSE)
    
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
  fAmpLED(0) , 
  fCalibByData(kFALSE)
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
  cout<<" @@@@ fCalibByData "<<fCalibByData<<endl;

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
      /*     
      if (nmips<17) {
	ok=false;
	return ok;
      } 
      */     
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
		  //	  ok=false;
		  printf("no peak in LED spectrum for PMT %i amplitude %i\n",i,im);
		  ledmean=0;
		  //  return ok;
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
      if(!fCalibByData)
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
    }
  Float_t xpoint, ypoint, xdata[250], ydata[250];
  Int_t ipmt;
  if(fCalibByData) {
    cout<<" read ingraph "<<endl;
    ifstream ingraph ("calibfit.txt");
    for (Int_t i=0; i<24; i++) 
      {
	for (Int_t ip=0; ip<200; ip++) 
	  {
	    ingraph>>ipmt>>xpoint>>ypoint; 
	    // cout<<i<<" "<<ipmt<<" "<<ip<<" "<< xpoint<<" "<<ypoint<<endl;  
	    xdata[ip]=xpoint;
	    ydata[ip]=ypoint;
	  }
	for (Int_t ip=200; ip<250; ip++) {
	    xdata[ip] =xdata[ip-1]+10;
	    ydata[ip]=ydata[199];
	}
       
	TGraph *grwalkqtc = new TGraph (250,xdata,ydata);
	grwalkqtc->Print();
	grwalkqtc->SetTitle(Form("PMT%i",i) );
	fWalk.AddAtAndExpand(grwalkqtc,i);
	for (Int_t ip=0; ip<250; ip++)
	  {
	    xdata[ip]=0;
	    ydata[ip]=0;
	  }
      }
  }

  
return ok;
}
void AliT0CalibWalk::SetWalk2015(TString filename)
{
  printf("AliT0CalibWalk::SetWalk2015\n");
  //set zero LED correction
  Float_t ampled[200], walkled[200], ampqtc[200];;
  for (int  i=0; i<200; i++) {
    ampled[i] = Float_t(i*10);  
    walkled[i]=0;
    ampqtc[i] = Float_t (i*100);
  }
  Int_t nbins=0;
  Double_t xgr[350], ygr[350];
  Float_t ampcut10[24] = {1563, 1544, 1555, 1645, 1745, 1555, 
			  1526, 1455, 1739, 1413, 1614, 1194,
			  1693, 1454, 1566, 1670, 1516, 1810, 
			  1607, 1516, 1675, 1566, 1443, 1545 }; 
  
  Float_t meanmax[24]= {3077.13, 3069.43, 3070.27, 3123.97, 3150.55, 3127.22,
			3110.7,  3088.24, 3164.09, 3154.94, 3142.93, 3103.78, 
			3091.14, 3055.07, 3043.82, 3051.3,  3062.58, 3087.11,
			3055.03, 3072.61, 3064.25, 3120.34, 3115.8, 3095.93};
  TProfile * prQTC_CFD[24];  
  TFile *f = new TFile(filename.Data());
  if (!f) {
    printf (" no file \n");
    return;
  }
  for (int i=0; i<24; i++) {
    prQTC_CFD[i]= (TProfile*) f->Get(Form("Slew%i",i+1) );
  }
  // collect graph
  for (Int_t i=0; i<24; i++)
    {   
      nbins=0;
      cout<<" nachalo "<<i<<endl;
      for (Int_t j=10; j<350; j++) {
	Int_t nentr = prQTC_CFD[i]->GetBinEntries(j);
	Float_t prqt = prQTC_CFD[i]->GetBinContent(j);
	xgr[nbins]=prQTC_CFD[i]->GetBinCenter(j);
	cout<<" bin "<<j<<" cont "<<prqt<<endl;
	if (prQTC_CFD[i]->GetBinCenter(j)>ampcut10[i]) {
	  if(nentr>1000) {
	    ygr[nbins]=prqt-meanmax[i];
	  }
	  else 
	    ygr[nbins]= ygr[nbins-1];
	  cout<<nbins<<" "<< xgr[nbins]<<" "<< ygr[nbins]<<endl;
	  nbins++;
	}
      }
      TGraph *grwalkqtc = new TGraph (nbins,xgr,ygr);
      grwalkqtc->SetTitle(Form("PMT%i",i+1));
      fWalk.AddAtAndExpand(grwalkqtc,i);
      TGraph *grwalkled = new TGraph (100,ampled,walkled);
      fAmpLEDRec.AddAtAndExpand(grwalkled,i);
      grwalkled->SetTitle(Form("PMT%i",i+1));
      
      //fit amplitude graphs to make comparison wth new one	  
      TGraph *grampled = new TGraph (200,ampled,ampled);
      TGraph *grqtc = new TGraph (200,ampqtc,ampqtc);
      fQTC.AddAtAndExpand(grqtc,i);	 
      fAmpLED.AddAtAndExpand(grampled,i);
      cout<<" add amp "<<i<<endl;
      
      
    }
  printf(" AliT0CalibWalk return \n");
  
}
//__________________________________________________________________________________________
void AliT0CalibWalk::SetWalkDima(TString filename)
{
  Float_t ampled[200], walkled[200], ampqtc[200];;
  for (int  i=0; i<200; i++) {
    ampled[i] = Float_t(i*10);  
    walkled[i]=0;
    ampqtc[i] = Float_t (i*100);
  }
 TFile *file = new TFile(filename.Data());
 file->ls();
  if (!file) {
    printf (" no file \n");
    return;
  }
  TDirectoryFile *dr = (TDirectoryFile*)file->Get("resultGraphs");
  TGraph *currGraph ;
  TString aPMTname;
  for(Int_t iPMT = 0; iPMT < 24; iPMT++)
    {
      if(iPMT<12)  aPMTname = Form("C_%02d_QTCCFDgraph", iPMT+1);
      else aPMTname = Form("A_%02d_QTCCFDgraph", iPMT-12+1);
      printf("  %s  \n",aPMTname.Data());
      currGraph = (TGraph*)dr->Get(aPMTname.Data());
      //     currGraph = (TGraph*)file->FindObjectAny(aPMTname.Data());
      currGraph->SetTitle(Form("PMT%i",iPMT+1));
      fWalk.AddAtAndExpand(currGraph,iPMT);
      //      fWalk.At(iPMT)->Print();
      TGraph *grwalkled = new TGraph (100,ampled,walkled);
      fAmpLEDRec.AddAtAndExpand(grwalkled,iPMT);
      
      TGraph *grampled = new TGraph (200,ampled,ampled);
      TGraph *grqtc = new TGraph (200,ampqtc,ampqtc);
      fQTC.AddAtAndExpand(grqtc,iPMT);	 
       fAmpLED.AddAtAndExpand(grampled,iPMT);
      cout<<" add amp "<<iPMT<<endl;
    }
      

}
//__________________________________________________________________________________________
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
  //  printf(" mean %f max %f sigma %f \n",mean, meanEstimate , sigma);  

  delete fit;
}

//______________________________________________________________________________
void AliT0CalibWalk::SetWalkZero() 
{
  Float_t amp[100], walk[100];
  for (int  i=0; i<100; i++) {
    amp[i] = Float_t(i*100);  
    walk[i]=0;
  }
  TGraph *gramp, *grwalk;
  for (int igr=0; igr<24; igr++) {
    gramp = new TGraph(100, amp,amp);
    grwalk = new TGraph(100, amp,walk);
    grwalk->SetTitle(Form("PMT%i",igr+1));
    gramp->SetTitle(Form("PMT%i",igr+1));
    fWalk.AddAtAndExpand(grwalk,igr);
    fAmpLEDRec.AddAtAndExpand(grwalk,igr);
    fQTC.AddAtAndExpand(gramp,igr);	 
    fAmpLED.AddAtAndExpand(gramp,igr);
     
  }
}
