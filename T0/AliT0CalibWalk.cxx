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
 
  Int_t npeaks = 20;
  Int_t sigma=3.;
  Bool_t down = false;

  Int_t index[20];
  Char_t buf1[10], buf2[10],buf3[20];
  
  
  for (Int_t i=0; i<24; i++)
    {
      if(i>11)	{  sprintf(buf1,"T0_A_%i_CFD",i+1-12); }

      if(i<12)	{ sprintf(buf1,"T0_C_%i_CFD",i+1); }

      sprintf(buf2,"CFD_QTC%i",i+1);
      sprintf(buf3,"CFD_LED%i",i+1);
      cout<<buf1<<" "<<buf2<<" "<<buf3<<endl;
      TH2F *qtc_cfd = (TH2F*) gFile->Get(buf2);
      TH2F *led_cfd = (TH2F*) gFile->Get(buf3);
      TH1F *cfd = (TH1F*) gFile->Get(buf1);
      //       cfd->Draw();
      TSpectrum *s = new TSpectrum(2*npeaks,1);
      Int_t nfound = s->Search(cfd,sigma," ",0.05);
      cout<<"Found "<<nfound<<" peaks sigma "<<sigma<<endl;;
      if(nfound!=0){
	Float_t *xpeak = s->GetPositionX();
  	TMath::Sort(nfound, xpeak, index,down);
  	Float_t xp = xpeak[index[0]];
	Int_t xbin = cfd->GetXaxis()->FindBin(xp);
	Float_t yp = cfd->GetBinContent(xbin);
	cout<<"xbin = "<<xbin<<"\txpeak = "<<xpeak[1]<<"\typeak = "<<yp<<endl;
	Float_t hmax = xp+10*sigma;
	Float_t hmin = xp-10*sigma;
	cout<<hmin<< " "<<hmax<<endl;
	//	    cfd->GetXaxis()->SetRange(hmin,hmax);
	//	    TF1 *g1 = new TF1("g1", "gaus", hmin, hmax);
	// cfd->Fit("g1","R");
	//	Int_t hminbin=  qtc_cfd->GetXaxis()->GetFirst();
	//	Int_t hmaxbin=  qtc_cfd->GetXaxis()->GetLast();
	Int_t nbins= qtc_cfd->GetXaxis()->GetNbins();
	cout<<" nbins "<<nbins<<endl;
	//	cout<<"  qtc_cfd "<<hminbin<<" "<<hmaxbin<<" "<<nbins<<endl;
	//  qtc_cfd->Draw();
	TProfile *pr_y = qtc_cfd->ProfileX();
	//	Float_t maxHr=pr_y->GetBinCenter(hmaxbin);
	
	pr_y->SetMaximum(hmax);
	pr_y->SetMinimum(hmin);
	Int_t np=nbins/20;
	Double_t *xx = new Double_t[np];
	Double_t *yy = new Double_t[np];
	Int_t ng=0;
	Double_t yg=0;
	for (Int_t ip=1; ip<nbins; ip++)
	  {
	    if(ip%20 != 0 ) {
	      if (pr_y->GetBinContent(ip) !=0)
		yg +=pr_y->GetBinContent(ip);
	      //	cout<<ng<<" "<<pr_y->GetBinContent(ip)<<" "<<yg<<endl;
	      ng++;}
	    else {
	      xx[ip/20] = Float_t (pr_y->GetBinCenter(ip));
	      yy[ip/20] = yg/ng;
	      yg=0;
	      ng=0;
	      //	     	      cout<<ip<<" "<<ip/20<<" "<< xx[ip/20]<<" "<< yy[ip/20]<<endl;
	    }
	  }
	TGraph *gr = new TGraph(np,xx,yy);
	gr->SetMinimum(hmin);
	gr->SetMaximum(hmax);
	fWalk.AddAtAndExpand(gr,i);	  
	//	fWalk.AddAtAndExpand(gr,i+12);

	Int_t nbinsled= led_cfd->GetXaxis()->GetNbins();
	cout<<" nbins led "<<nbinsled<<endl;
	//	cout<<"  qtc_cfd "<<hminbin<<" "<<hmaxbin<<" "<<nbins<<endl;
	//  qtc_cfd->Draw();
	TProfile *prled_y = led_cfd->ProfileX();
	//	Float_t maxHr=pr_y->GetBinCenter(hmaxbin);
	
	prled_y->SetMaximum(hmax);
	prled_y->SetMinimum(hmin);
	Int_t npled=nbinsled/20;
	Double_t *xxled = new Double_t[np];
	Double_t *yyled = new Double_t[np];
	Int_t ngled=0;
	Double_t ygled=0;
	for (Int_t ip=1; ip<nbinsled; ip++)
	  {
	    if(ip%20 != 0 ) {
	      if (prled_y->GetBinContent(ip) !=0)
		ygled +=prled_y->GetBinContent(ip);
	      //	cout<<ng<<" "<<pr_y->GetBinContent(ip)<<" "<<yg<<endl;
	      ngled++;}
	    else {
	      xxled[ip/20] = Float_t (prled_y->GetBinCenter(ip));
	      yyled[ip/20] = ygled/ngled;
	      ygled=0;
	      ngled=0;
	    }
	  }
	TGraph *grled = new TGraph(npled,xxled,yyled);
	grled->SetMinimum(hmin);
	grled->SetMaximum(hmax);
	fAmpLEDRec.AddAtAndExpand(grled,i);	  
 	//	delete [] xx;
	//	delete [] yy;
	
      }
      
    }
  


}

