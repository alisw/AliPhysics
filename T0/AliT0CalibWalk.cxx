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
#include "AliT0LookUpValue.h"
#include "AliLog.h"
#include "AliRun.h"

#include <TObjArray.h>
#include <TGraph.h>
#include <TFile.h>
#include <TAxis.h>
#include <TH2F.h> 
#include <TMath.h>
#include <TSystem.h>
#include <Riostream.h>
#include <TSpectrum.h>
#include <TVirtualFitter.h>
#include <TProfile.h>

#include <string>

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
  TFile *gFile = TFile::Open(laserFile);
  gSystem->Load("libSpectrum");
  TGraph *gr[24];
  Int_t npeaks = 20;
  Int_t sigma=3;
  Bool_t down=false;

  Int_t index[20];
  Char_t buf1[10], buf2[10],title[10];

  for (Int_t i=0; i<24; i++)
  {
    sprintf(buf1,"T0_C_%i_CFD",i+1);
    sprintf(buf2,"CFD_QTC%i",i+1);
    // cout<<buf1<<" "<<buf2<<endl;
    TH2F *qtc_cfd = (TH2F*) gFile->Get(buf2);
    TH1F *cfd = (TH1F*) gFile->Get(buf1);
    TSpectrum *s = new TSpectrum(2*npeaks,1.);
    Int_t nfound = s->Search(cfd,sigma," ",0.05);
    // cout<<"Found "<<nfound<<" peaks sigma "<<sigma<<endl;
    if(nfound!=0)
    {
      Float_t *xpeak = s->GetPositionX();
      TMath::Sort(nfound, xpeak, index,down);
      Float_t xp = xpeak[index[0]];
      Float_t hmax = xp+10*sigma;
      Float_t hmin = xp-10*sigma;
      Int_t nbins= qtc_cfd->GetXaxis()->GetNbins();
      TProfile *pr_y = qtc_cfd->ProfileX();
      pr_y->SetMaximum(hmax);
      pr_y->SetMinimum(hmin);
      Int_t np=nbins/20;
      Double_t *xx = new Double_t[np];
      Double_t *yy = new Double_t[np];
      Int_t ng=0;
      Double_t yg=0;
      for (Int_t ip=1; ip<nbins; ip++)
      {
        if(ip%20 != 0 )
        {
          if (pr_y->GetBinContent(ip) !=0)
          {
            yg +=pr_y->GetBinContent(ip);
          }
          ng++;
        }
        else
        {
          xx[ip/20] = Float_t (pr_y->GetBinCenter(ip));
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
      fWalk.AddAtAndExpand(gr[i],i);	  
       delete [] xx;
      delete [] yy;
      delete gr[i];
    }
  }

  gFile->Close();
  delete gFile;
}

