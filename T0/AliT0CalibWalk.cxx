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


ClassImp(AliT0CalibWalk)

//________________________________________________________________
  AliT0CalibWalk::AliT0CalibWalk():   TNamed(),
				      fWalk(0),
				      fAmpLEDRec(0)
{
  //
}

//________________________________________________________________
AliT0CalibWalk::AliT0CalibWalk(const char* name):TNamed(),
				      fWalk(0),
				      fAmpLEDRec(0)				      
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
  fAmpLEDRec(0)

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
  cout<<" number of data "<<i<<endl;
 
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
  
}

//________________________________________________________________

void AliT0CalibWalk::MakeWalkCorrGraph(const char *laserFile)
{
  //make walk corerction for preprocessor

  TFile *gFile = TFile::Open(laserFile);
  //gSystem->Load("libSpectrum");
  TGraph *gr[24];
  TGraph *grLED[24];
  Int_t npeaks = 20;
  Int_t sigma=3;
  Bool_t down=false;

  Int_t index[20];
  Char_t buf1[10], buf2[10], buf3[10], title[10], title2[10], titleLED[10], title2LED[10];
  TSpectrum *s = new TSpectrum(2*npeaks,1.);

  for (Int_t d=0; d<2; d++)
  {
    for (Int_t i=0; i<12; i++)
    {
      sprintf(buf1,"T0_C_%i_CFD",i+1);
      if (d==0)
      {		
        sprintf(buf2,"CFD_QTC%i",i+1);
      }
      else
      {
	sprintf(buf3,"CFD_LED%i",i+1);
      }
      TH2F *qtccfd = (TH2F*) gFile->Get(buf2);
      TH2F *ledcfd = (TH2F*) gFile->Get(buf3);
      // cout<<buf1<<" "<<buf2<<endl;
      TH1F *cfd = (TH1F*) gFile->Get(buf1);
      Int_t nfound = s->Search(cfd,sigma,"goff",0.05);
      // cout<<"Found "<<nfound<<" peaks sigma "<<sigma<<endl;
      if(nfound!=0)
      {
        Float_t *xpeak = s->GetPositionX();
        TMath::Sort(nfound, xpeak, index,down);
        Float_t xp = xpeak[index[0]];
        Float_t hmax = xp+10*sigma;
        Float_t hmin = xp-10*sigma;
	if (d==0)
	{
          Int_t nbins= qtccfd->GetXaxis()->GetNbins();
          TProfile *prY = qtccfd->ProfileX();
          prY->SetMaximum(hmax);
          prY->SetMinimum(hmin);
          Int_t np=nbins/20;
          Double_t *xx = new Double_t[np];
          Double_t *yy = new Double_t[np];
          Int_t ng=0;
          Double_t yg=0;
          for (Int_t ip=1; ip<nbins; ip++)
          {
            if(ip%20 != 0 )
            {
              if (prY->GetBinContent(ip) !=0)
              {
                yg +=prY->GetBinContent(ip);
              }
              ng++;
            }
            else
            {
              xx[ip/20] = Float_t (prY->GetBinCenter(ip));
              yy[ip/20] = yg/ng;
              yg=0;
              ng=0;
            }
          }
          sprintf(title,"Walk %i",i+1);
          gr[i] = new TGraph(np,xx,yy);
          gr[i]->SetTitle(title);
          gr[i]->SetMinimum(hmin);
          gr[i]->SetMaximum(hmax);
          gr[i]->SetMarkerStyle(7);
          sprintf(title2,"Walk %i",i+13);
          gr[i+12] = new TGraph(np,xx,yy);
          gr[i+12]->SetTitle(title2);
          gr[i+12]->SetMinimum(hmin);
          gr[i+12]->SetMaximum(hmax);
          gr[i+12]->SetMarkerStyle(7);

          fWalk.AddAtAndExpand(gr[i],i);	  
          fWalk.AddAtAndExpand(gr[i+12],i+12);
          delete [] xx;
          delete [] yy;
	  delete prY;
	}
	else
	{
	  Int_t nbinsLED= ledcfd->GetXaxis()->GetNbins();
          TProfile *prYLED = ledcfd->ProfileX();
          prYLED->SetMaximum(hmax);
          prYLED->SetMinimum(hmin);
          Int_t npLED=nbinsLED/20;
          Double_t *xxLED = new Double_t[npLED];
          Double_t *yyLED = new Double_t[npLED];
          Int_t ngLED=0;
          Double_t ygLED=0;
          for (Int_t ip=1; ip<nbinsLED; ip++)
          {
            if(ip%20 != 0 )
            {
              if (prYLED->GetBinContent(ip) !=0)
              {
                ygLED +=prYLED->GetBinContent(ip);
              }
              ngLED++;
            }
            else
            {
              xxLED[ip/20] = Float_t (prYLED->GetBinCenter(ip));
              yyLED[ip/20] = ygLED/ngLED;
              ygLED=0;
              ngLED=0;
            }
          }
          sprintf(titleLED,"Walk LED %i",i+1);
          grLED[i] = new TGraph(npLED,xxLED,yyLED);
          grLED[i]->SetTitle(titleLED);
          grLED[i]->SetMinimum(hmin);
          grLED[i]->SetMaximum(hmax);
          grLED[i]->SetMarkerStyle(7);
          sprintf(title2LED,"Walk LED%i",i+13);
          grLED[i+12] = new TGraph(npLED,xxLED,yyLED);
          grLED[i+12]->SetTitle(title2LED);
          grLED[i+12]->SetMinimum(hmin);
          grLED[i+12]->SetMaximum(hmax);
          grLED[i+12]->SetMarkerStyle(7);

          fAmpLEDRec.AddAtAndExpand(grLED[i],i);
          fAmpLEDRec.AddAtAndExpand(grLED[i+12],i+12);
          delete [] xxLED;
          delete [] yyLED;
	  delete prYLED;
	}
      }
      delete cfd;
      delete qtccfd;
      delete ledcfd;			
    }
  }

}

