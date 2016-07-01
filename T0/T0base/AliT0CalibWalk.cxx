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
  // 2016 1MIP field 0.5T
  /*Float_t cfd1mip[24] = {9415.1, 9410.65, 9410.9, 9459.59, 9458.37, 9464.44,
			 9449.28, 9433.68, 9498.78, 9494.35, 9470.42, 9443.67, 
			 9451.75, 9404.85, 9368.4, 9389.45, 9394.14, 9476.02, 
			 9409.17, 9411.73, 9428.91, 9434.13, 9455.88, 9495.24};*/
   Float_t cfd1mip[24] =  {9448.56, 9440.15, 9443.97, 9492.87, //field 0.5T
			  9489.09, 9494.43, 9482.78, 9470.24, 
			  9528.77, 9529.98, 9498.3, 9478.93,
			  9487.18, 9437.12, 9410.01, 9416.12, 
			  9425.38, 9506.05, 9443.7, 9444.27, 
			  9459.62, 9471, 9482.63, 9528.97 };

 
  Float_t sigma,meancfd, meanqtcold, meanqtc0new, meanqtc1new;
  Bool_t ok=true;
  gFile = TFile::Open(laserFile);
  if(!gFile) {
    AliError("No input laser data found ");
  }
  else
    {
      //      gFile->ls();
    }    
  Float_t x1[50], y1[50],  x2[50], x3[50];
  Float_t cfd0 = 0;
  
  for (Int_t ii=0; ii<50; ii++)
    x1[ii] = y1[ii] = x2[ii] = x3[ii] = 0; 
  
  for (Int_t ipmt=0; ipmt<24; ipmt++)
    {
      cfd0 = 0;
      int nmips=30;
      for (Int_t im=0; im<nmips; im++)
	{	      
	  TString cfd = Form("hCFD%i_%i",ipmt+1,im+1);
	  TString qtc = Form("hQTC%i_%i",ipmt+1,im+1);
	  TString qtc0new = Form("hQTC0new%i_%i",ipmt+1,im+1);
	  TString qtc1new = Form("hQTC1new%i_%i",ipmt+1,im+1);
	  
	  TH1F *hCFD = (TH1F*) gFile->Get(cfd.Data()) ;
	  TH1F *hQTC = (TH1F*) gFile->Get(qtc.Data());
	  TH1F *hQTC0new = (TH1F*) gFile->Get(qtc0new.Data()) ;
	  TH1F *hQTC1new = (TH1F*) gFile->Get(qtc1new.Data()) ;
	  //	  if( hCFD && hCFD->GetEntries()<2 ) {
	  //	    ok=false;
	    //	    printf("no peak in CFD spectrum for PMT %i amplitude %i\n",ipmt,im);
	  //	    return ok;
	  //	  }
	  if (hCFD->GetEntries()>200)
	    GetMeanAndSigma(hCFD, meancfd,  sigma);
	  else
	    {
	    meancfd = cfd1mip[ipmt];
	    cout<<ipmt<<" "<<im<<" !!! no entried "<<meancfd<<endl;
	    }
	  if (hQTC->GetEntries()>200)
	    GetMeanAndSigma(hQTC, meanqtcold,  sigma);
	  else
	    meanqtcold=hQTC->GetMean();
	  if (hQTC0new->GetEntries()>200)
	    GetMeanAndSigma(hQTC0new, meanqtc0new,  sigma);
	  else 
	    meanqtc0new=hQTC0new->GetMean();
	  if (hQTC1new->GetEntries()>200)
	    GetMeanAndSigma(hQTC1new, meanqtc1new,  sigma);
	  else 
	    meanqtc1new=hQTC1new->GetMean();
	  	  
	  //	  cout<<"@@@ "<<im<<" "<<ipmt<<"  "<<meancfd<<" ("<<hCFD->GetEntries()<<") "
	  //	      <<" "<< meanqtcold<<" ("<< hQTC->GetEntries()<<") "  <<endl;
	  
	  Int_t   maxBin = hCFD->GetMaximumBin(); 
	  Double_t  meanEstimate = hCFD->GetBinCenter( maxBin); 
	  if(TMath::Abs(meanEstimate - meancfd) > 20 ) {
	    cout<<" razoshlos' "<<ipmt<<" "<<im<<" meanEstimate "<<meanEstimate<<
	      " "<<meancfd<<endl;
	    meancfd = meanEstimate; 
	  }
	  y1[im+2] =  meancfd - cfd1mip[ipmt];
	  x1[im+2] = meanqtcold;
	  //	  if (ipmt==2 && im==5) x1[im+2]=x1[im+1];
	  x2[im+2] = meanqtc0new;
	  x3[im+2] = meanqtc1new;
	  if(ipmt==22)	  cout<<im<<" "<<ipmt<<"  "<<meancfd<<" ("<<hCFD->GetEntries()<<") "<<y1[im+2]  <<" "<< meanqtcold<<" ("<< hQTC->GetEntries()<<") "<<x1[im+2]  <<endl;
	  
	  
	  if (hQTC) delete  hQTC;
	  if (hCFD) delete  hCFD;
	  if (hQTC0new) delete  hQTC0new;
	  if (hQTC1new) delete  hQTC1new;
	  
	}
      
      if ( ipmt==0) {
	y1[0]= y1[1]=y1[2]=3;
	for (int ia=3; ia<6; ia++) y1[ia]=5;
	for (int ia=5; ia<10; ia++) y1[ia]=3;
      }
      if(ipmt==1 ) for (int ia=3; ia<4; ia++) y1[ia]=0;
      if(ipmt==2 ) y1[0]= y1[1]=y1[2]=-4;
      if(ipmt==3 ) {
	y1[0]= 	y1[1]=10;
	y1[2]=-5;
	for (int ia=3; ia<10; ia++) y1[ia]=0;
	y1[22]=0;
      }
     if( ipmt==4) {
        x1[0]=x1[3]-40;
	x1[1]=x1[3]-20;
	x1[2]=x1[3]-10;
	y1[0]=y1[1]=y1[2]=y1[3]=0;
       }
      x1[0]=x1[2]-40;
      x1[1]=x1[2]-20;
      //  cout<<" vpered "<<x1[2]<<" "<<x1[1]<<" "<<x1[0]<<endl;
      x2[0]=x2[1]=x2[2];
      x3[0]=x3[1]=x3[2];
      y1[0]=y1[1]=y1[2];
      
      if(ipmt==5 ) {  y1[0]=y1[1]=y1[2]=y1[3]=5; y1[4]=3;}
      if(ipmt==6 ) {  y1[0]=y1[1]=5;  y1[2]=0;}
      if(ipmt==7 )  	for (int iii=0; iii<9; iii++)  y1[iii]=0; 
      if(ipmt==8 ) {  y1[0]=y1[1]=y1[2]=y1[3]=0; }
      if(ipmt==9 ) { 
	for (int ia=0; ia<5; ia++) y1[ia]=-5;
	for (int ia=20; ia<23; ia++) y1[ia]=-1;
	for (int ia=23; ia<27; ia++) y1[ia]=-5;
     }
      
     if(ipmt==10 ) {  y1[0]=y1[1]=y1[2]=5; y1[3]=y1[4]=2; }
     if(ipmt==11 ) {  y1[0]=y1[1]=y1[2]=-1; y1[3]=1;y1[4]=2; }
     if(ipmt==15 ) {  y1[0]=y1[1]=y1[2]=6;  }
     if(ipmt==16 ) {  y1[0]=y1[1]=y1[2]=0;  }
     if(ipmt==17 ) {  y1[0]=y1[1]=y1[2]=0;  }
     if(ipmt==18 ) {  y1[0]=y1[1]=y1[2]=y1[4]=0;  }
     if(ipmt==19 ) {  y1[0]=y1[1]=y1[2]=2;  }
     if(ipmt==20 ) {  y1[0]=y1[1]=y1[2]=4;  }
     if(ipmt==21 ) {  y1[0]=y1[1]=y1[2]=0;  }
     if(ipmt==22 ) {  y1[0]=y1[1]=2;y1[2]=1;  }
     if(ipmt==23 ) {  y1[0]=y1[1]=y1[2]=0;  }

 
      if(ipmt==18 ) y1[15]=9443 -  cfd1mip[ipmt];
      
       if(ipmt==0) cout<<"Making graphs..."<<endl;
       cout<<"________________________________"<<ipmt<<"_____________"<<endl;   
      TGraph *grwalkqtc = new TGraph (nmips+2,x1,y1);
      grwalkqtc->SetTitle(Form("PMT%i",ipmt));
      if(ipmt<12) grwalkqtc->Print();
      fWalk.AddAtAndExpand(grwalkqtc,ipmt);
       cout<<"________________________________"<<ipmt<<"_____________"<<endl;   
      TGraph *grwalkled = new TGraph (nmips+2,x2,y1);
      grwalkled->SetTitle(Form("PMT%i",ipmt));
      fAmpLEDRec.AddAtAndExpand(grwalkled,ipmt);
      //   grwalkled->Print();
      //	  cout<<" add walk "<<i<<endl;
      
      //fit amplitude graphs to make comparison wth new one	  
      TGraph *grampled = new TGraph (nmips+2,x3,y1);
      TGraph *grqtc = new TGraph (nmips+2,x1,x1);
      fQTC.AddAtAndExpand(grqtc,ipmt);	 
      fAmpLED.AddAtAndExpand(grampled,ipmt);
      //	  cout<<" add amp "<<i<<endl;
      
      if(ipmt==23) cout<<"Graphs created..."<<endl;   
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
