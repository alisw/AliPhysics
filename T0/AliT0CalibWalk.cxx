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
void AliT0CalibWalk::SetWalk(Int_t ipmt)
{
  //read QTC walk graph from external file

  Int_t mv, ps; 
  Int_t x[70000], y[70000], index[70000];
  Float_t time[10000],amplitude[10000];
  string buffer;
  Bool_t down=false;
  
  const char * filename = gSystem->ExpandPathName("$ALICE_ROOT/T0/data/CFD-Amp.txt");
  ifstream inFile(filename);
  if(!inFile) AliError(Form("Cannot open file %s !",filename));
  
  Int_t i=0;
  while(getline(inFile,buffer)){
    inFile >> ps >> mv;

    x[i]=ps; y[i]=mv;
    i++;
  }
  inFile.close();
  //  cout<<" number of data "<<i<<endl;
 
  TMath::Sort(i, y, index,down);
  Int_t amp=0, iin=0, isum=0, sum=0;
  Int_t ind=0;
  for (Int_t ii=0; ii<i; ii++)
    {
      ind=index[ii];
      if(y[ind] == amp)
	{
	  sum +=x[ind];
	  iin++;
	}
      else
	{
	  if(iin>0)
	    time[isum] = Float_t (sum/(iin));
	  else
	    time[isum] =Float_t (x[ind]);
	  amplitude[isum] = Float_t (amp);
	  amp=y[ind];
	  iin=0;
	  isum++;
	  sum=0;
	}
    }

  inFile.close();

  TGraph* gr = new TGraph(isum, amplitude, time);
  fWalk.AddAtAndExpand(gr,ipmt);

  //should be change to real
  Double_t xq[10] = { 1220, 1370, 1542, 1697, 1860, 2023,2171,2331,2495,2684};
  Double_t yq[10] = {1,2,3,4,5,6,7,8,9,10};
  TGraph* gr1 = new TGraph(10, xq, yq);
  fQTC.AddAtAndExpand(gr1,ipmt);
}

//________________________________________________________________

void AliT0CalibWalk::SetAmpLEDRec(Int_t ipmt)
{
  // read LED walk from external file
  Float_t mv, ps; 
  Float_t x[100], y[100];
  string buffer;
  
  const char * filename = gSystem->ExpandPathName("$ALICE_ROOT/T0/data/CFD-LED.txt");
   ifstream inFile(filename);
  if(!inFile) {AliError(Form("Cannot open file %s !",filename));}
  
  inFile >> mv>>ps;
  Int_t i=0;
  
  while(getline(inFile,buffer)){
    x[i]=mv; y[i]=ps;	
    inFile >> mv >> ps;
    i++;
  }
  inFile.close();
  Float_t y1[100], x1[100];
  for (Int_t ir=0; ir<i; ir++){
    y1[ir]=y[i-ir]; x1[ir]=x[i-ir];}
  TGraph* gr = new TGraph(i,y1,x1);
  fAmpLEDRec.AddAtAndExpand(gr,ipmt);
  //should be change to real
  Double_t xq[10] = { 411, 412,413,415,417,419,422,428,437,452};
  Double_t yq[10] = {1,2,3,4,5,6,7,8,9,10};
  TGraph* gr1 = new TGraph(10, xq, yq);
  fAmpLED.AddAtAndExpand(gr1,ipmt);
  
}

//________________________________________________________________

void AliT0CalibWalk::MakeWalkCorrGraph(const char *laserFile)
{
  //make walk corerction for preprocessor
  Int_t nmips=14;
  gFile = TFile::Open(laserFile);
  if(!gFile) {
    AliError("No input laser data found ");
  }
  else
    {
      //      gFile->ls();

      Float_t x1[14], y1[14]; 
      Float_t x2[14], y2[14];
  
      Float_t mips[14] = {1.,2,3,4,5,6,7,8,9,10,20,30,40,50};
     
      Float_t xx1[14],yy1[14], xx[14];
      for (Int_t ii=0; ii<14; ii++)
	{      x1[ii]=0; y1[ii]=0; x2[ii]=0; y2[ii]=0; }
      
      
      TH2F*  hCFDvsQTC[24][14]; TH2F*  hCFDvsLED[24][14];
      
      for (Int_t i=0; i<24; i++)
	{
	  for (Int_t im=0; im<nmips; im++)
	    {	      
	      TString qtc = Form("CFDvsQTC%i_%i",i+1,im+1);
	      TString led = Form("CFDvsLED%i_%i",i+1,im+1);
	      hCFDvsQTC[i][im] = (TH2F*) gFile->Get(qtc.Data()) ;
	      hCFDvsLED[i][im] = (TH2F*) gFile->Get(led.Data());
	      
	      if(!hCFDvsQTC[i][im] || !hCFDvsLED[i][im]) 
	      	AliWarning(Form(" no walk correction data in LASER DA for channel %i for amplitude %f MIPs",i,mips[im]));
	      
	      Float_t qtc0 = hCFDvsQTC[i][0]->GetMean(2);	      
	      if(hCFDvsQTC[i][im])
		{	  
		  x1[im] = hCFDvsQTC[i][im]->GetMean(1);
		  y1[im] =  hCFDvsQTC[i][im]->GetMean(2) - qtc0;
		}
	      Float_t led0 = hCFDvsLED[i][0]->GetMean(2);
	      if( hCFDvsLED[i][im]){
		  x2[im] = hCFDvsLED[i][im]->GetMean(1);
		  y2[im] =  hCFDvsLED[i][im]->GetMean(2) - led0;
		}
	      xx[im]=mips[im];
	    }
	  for (Int_t imi=0; imi<nmips; imi++)
	    {
	      yy1[imi] = Float_t (mips[nmips-imi-1]);
	      xx1[imi] = x2[nmips-imi-1]; 
	    }
	  if(i==0) cout<<"Making graphs..."<<endl;
	  
	  TGraph *grwalkqtc = new TGraph (nmips,x1,y1);
	  TGraph *grwalkled = new TGraph (nmips,x2,y2);
	  fWalk.AddAtAndExpand(grwalkqtc,i);
	  fAmpLEDRec.AddAtAndExpand(grwalkled,i);
	 
	  //fit amplitude graphs to make comparison wth new one	  
	  TGraph *grampled = new TGraph (nmips,xx1,yy1);
	  TGraph *grqtc = new TGraph (nmips,x1,xx);
	  fQTC.AddAtAndExpand(grqtc,i);	 
	  fAmpLED.AddAtAndExpand(grampled,i);

	  //fit amplitude graphs
	  Double_t parled [2], parqtc[2];

	  TF1 *fuled = new TF1("fuled","expo",x1[0],x1[nmips]);
	  TF1 *fuqtc = new TF1 ("fuqtc","pol1",x1[0],x1[nmips]);
	  grqtc->Fit("fuqtc","Q"," ",x1[0],x1[nmips]);
	  grampled->Fit("fuled","Q"," ",xx1[0],xx1[nmips]);
	  fuled->GetParameters(&parled[0]);
	  fuqtc->GetParameters(&parqtc[0]);

	  for(Int_t ipar=0; ipar<2; ipar++) {
	    SetQTCpar(i,ipar,parqtc[ipar]);
	    SetAmpLEDpar(i,ipar,parled[ipar]);
	    //	    cout<<" pars :::: "<<i<<" "<<ipar<<" qtc "<<parqtc[ipar]<<" led "<<parled[ipar]<<endl;
	  }

	  if(i==23){
	    cout<<"Graphs created..."<<endl;
	  }
	}
    } //if gFile exits
}




