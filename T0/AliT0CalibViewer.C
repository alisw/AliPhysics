/* $Id$ */


#if !defined(__CINT__)
#include "TControlBar.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "AliRun.h"
#include "AliT0CalibData.h"
#include "AliT0AlignData.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif
TFile *gFile=0;

void AliT0CalibViewer()
{
  TControlBar *menu = new TControlBar("vertical","T0 CalibViewer");
  menu->AddButton("Open file","OpenFile()",
		  "Open file");
  menu->AddButton("Draw CFD","DrawCFD()",
		  "Draw all CFD");
  menu->AddButton("Draw LED","DrawLED()",
		  "Draw all LED");
  menu->AddButton("Draw LED-CFD","DrawLEDminCFD()",
		  "Draw LED-CFD");
  menu->AddButton("Draw QTC","DrawQTC()",
                  "Draw all QTC");
  menu->AddButton("Draw CFDvsQTC","DrawCFDvsQTC()",
		  "Draw CFD vs QTC");
  menu->AddButton("Draw CFDvsLED","DrawCFDvsLED()",
		  "Draw CFD vs LED-CFD");
  menu->Show();
}
void OpenFile()
{

const char *ft[]={"T0 raw files","t0tree*.root","All files","*",0,0};
  TString dir(".");
  TGFileInfo fi; fi.fFileTypes=ft; fi.fIniDir=StrDup(dir);
  new TGFileDialog(gClient->GetRoot(), 0x0, kFDOpen, &fi);

  if(!fi.fFilename) return;
  if(gFile){ gFile->Close(); gFile=0;}
  
  gFile=TFile::Open(fi.fFilename); 

}

//------------------------------------------------------------------------
void DrawCFD()
{
  Int_t npeaks = 20;
  Int_t sigma=3.;
  Bool_t down=false;
  Int_t index[20];
   
  TCanvas *c1 = new TCanvas("c1", "c1",0,48,1280,951);
  c1->Divide(4,3);
  gStyle->SetOptFit(1111);
  //c1->Divide(2,2);
  Char_t buf1[10];
  for (Int_t i=0; i<12; i++)
    {
      c1->cd(i+1);
      sprintf(buf1,"T0_C_%i_CFD",i+1);
      TH1F *cfd = (TH1F*) gFile->Get(buf1);
      //     cfd->Draw();
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
	Float_t hmax = xp+3*sigma;
	Float_t hmin = xp-3*sigma;
	cout<<hmin<< " "<<hmax<<endl;
	cfd->GetXaxis()->SetRange(hmin-20,hmax+20);
	//	cout<<" cfd range "<<mean<<" rms "<<rms<<" "<<hmin<<" "<<hmax<<endl;
	TF1 *g1 = new TF1("g1", "gaus", hmin, hmax);
	//	cfd->Fit("g1","R");
	cfd->Draw();
      }
    }
 
}
//------------------------------------------------------------------------
void DrawLED()
{
  Int_t npeaks = 10;
  Int_t sigma=10.;
  Bool_t down=false;
  Int_t index[20];
  
  TCanvas *c1 = new TCanvas("c1", "c1",0,48,1280,951);
  c1->Divide(4,3);
  
  Char_t buf1[20];
  for (Int_t i=0; i<12; i++)
    {
      c1->cd(i+1);
      sprintf(buf1,"T0_C_%i_LED",i+1);
      TH1F *cfd = (TH1F*) gFile->Get(buf1);
      cfd->Draw();
      TSpectrum *s = new TSpectrum(2*npeaks,1);
      Int_t nfound = s->Search(cfd,sigma," ",0.2);
      cout<<"Found "<<nfound<<" peaks sigma "<<sigma<<endl;;
      if(nfound!=0) {
	Float_t *xpeak = s->GetPositionX();
	TMath::Sort(nfound, xpeak, index,down);
	Float_t xp = xpeak[index[0]];
        Int_t xbin = cfd->GetXaxis()->FindBin(xp);
        Float_t yp = cfd->GetBinContent(xbin);
        cout<<"xbin = "<<xbin<<"\txpeak = "<<xpeak[index[0]]<<"\typeak = "<<yp<<endl;
        Float_t hmin=xp-3*sigma;
        Float_t hmax =xp+3*sigma;
        cfd->GetXaxis()->SetRange(hmin,hmax);
 	TF1 *g1 = new TF1("g1", "gaus", hmin, hmax);
	cfd->Fit("g1","R");
	//      cfd->Draw();
	
      }
      
    }
}

//------------------------------------------------------------------------
void DrawQTC()
{
  
  TCanvas *c1 = new TCanvas("c1", "QTC",0,48,1280,951);
  // c1->Divide(4,3);
  c1->Divide(2,3);
  
  Char_t buf1[10];
  for (Int_t i=0; i<12; i++)
    {
      c1->cd(i+1);
      sprintf(buf1,"QTC%i",i+1);
      TH1F *qtc = (TH1F*) gFile->Get(buf1);
      Float_t mean = qtc->GetMean();
      Float_t rms = qtc->GetRMS();
      Float_t hminR=mean - 0.1*mean;
      Float_t hmaxR =mean + 0.1*mean;
      qtc->GetXaxis()->SetRange(hminR,hmaxR);
      Float_t hmin=mean - 3*rms;
      Float_t hmax =mean + 3*rms;
      qtc->GetXaxis()->SetRange(hmin,hmax);
      // TF1 *g2 = new TF1("g2", "gaus", hmin, hmax);
       //             qtc->Fit("g2","RQ");
       qtc->Draw();
    }
  
  
}

//------------------------------------------------------------------------
void DrawLEDminCFD()
{
  Int_t npeaks = 10;
  Int_t sigma=10.;
  Bool_t down=false;
  Int_t index[20];
  
  TCanvas *c1 = new TCanvas("c1", "c1",0,48,1280,951);
  c1->Divide(4,3);
  
  Char_t buf1[20];
  for (Int_t i=0; i<12; i++)
    {
      c1->cd(i+1);
      sprintf(buf1,"LED-CFD%i",i+1);
      TH1F *cfd = (TH1F*) gFile->Get(buf1);
      cfd->Draw();
      TSpectrum *s = new TSpectrum(2*npeaks,1);
      Int_t nfound = s->Search(cfd,sigma," ",0.2);
      cout<<"Found "<<nfound<<" peaks sigma "<<sigma<<endl;;
      if(nfound!=0) {
	Float_t *xpeak = s->GetPositionX();
	TMath::Sort(nfound, xpeak, index,down);
	Float_t xp = xpeak[index[0]];
        Int_t xbin = cfd->GetXaxis()->FindBin(xp);
        Float_t yp = cfd->GetBinContent(xbin);
        cout<<"xbin = "<<xbin<<"\txpeak = "<<xpeak[index[0]]<<"\typeak = "<<yp<<endl;
        Float_t hmin=xp-3*sigma;
        Float_t hmax =xp+3*sigma;
        cfd->GetXaxis()->SetRange(hmin,hmax);
	TF1 *g1 = new TF1("g1", "gaus", hmin, hmax);
        cfd->Fit("g1","RQ");

	
      }
      
    }
  /*
  TCanvas *c1 = new TCanvas("c1", "c1",0,48,1280,951);
  c1->Divide(2,2);
  Char_t buf1[10];
  for (Int_t i=0; i<4; i++)
    {
      c1->cd(i+1);
      sprintf(buf1,"LED-CFD%i",i+1);
      TH1F *cfd = (TH1F*) file->Get(buf1);
      //  cout<<buf1<<" "<<cfd<<endl;
      //     cfd->Draw();
      //   cfd->GetXaxis()->SetRange(0,100);
      Float_t mean = cfd->GetMean();
      Float_t rms = cfd->GetRMS();
      Float_t hmin=mean - 3*rms;
      Float_t hmax =mean + 3*rms;
      cfd->GetXaxis()->SetRange(hmin-10,hmax+10);
      cout<<" cfd range "<<mean<<" rms "<<rms<<" "<<hmin<<" "<<hmax<<endl;
      //     TF1 *g1 = new TF1("g1", "gaus", hmin, hmax);
      //  cfd->Fit("g1","RQ");
      cfd->Draw();
    }
  */

}
//------------------------------------------------------------------------
void DrawCFDvsQTC()
{
  Int_t npeaks = 20;
  Int_t sigma=3.;
  Bool_t down=false;

  Int_t index[20];
  Char_t buf1[10], buf2[10];
  
  TCanvas *c1 = new TCanvas("c1", "c1",0,48,1280,951);
  gStyle->SetOptStat(0);
  c1->Divide(3,2);
  
  for (Int_t i=0; i<5; i++)
    {
      c1->cd(i+1);
      sprintf(buf1,"T0_C_%i_CFD",i+1);
      sprintf(buf2,"CFDvsQTC%i",i+1);
      cout<<buf1<<" "<<buf2<<endl;
      TH2F *qtc_cfd = (TH2F*) gFile->Get(buf2);
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
	Int_t hminbin=  qtc_cfd->GetXaxis()->GetFirst();
	Int_t hmaxbin=  qtc_cfd->GetXaxis()->GetLast();
	Int_t nbins= qtc_cfd->GetXaxis()->GetNbins();
	cout<<"  qtc_cfd "<<hminbin<<" "<<hmaxbin<<" "<<nbins<<endl;
	//  qtc_cfd->Draw();
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
	    if(ip%20 != 0 ) {
	      if (pr_y->GetBinContent(ip) !=0)
		yg +=pr_y->GetBinContent(ip);
	      //	cout<<ng<<" "<<pr_y->GetBinContent(ip)<<" "<<yg<<endl;
	      ng++;}
	    else {
	      xx[ip/20] = Float_t (ip);
	      yy[ip/20] = yg/ng;
	      yg=0;
	      ng=0;
	      //	      cout<<ip<<" "<<ip/20<<" "<< xx[ip/20]<<" "<< yy[ip/20]<<endl;
	    }
	  }
	TH2F *hr = new TH2F("hr"," ",np,0,nbins, np, hmin, hmax);
	hr->Draw();
	TGraph *gr = new TGraph(np,xx,yy);
	gr->SetMinimum(hmin);
	gr->SetMaximum(hmax);
	gr->SetMarkerStyle(20);
	gr->Draw("P");
	//	delete [] xx;
	//	delete [] yy;
	
	// pr_y->Rebin(10);
		   //  pr_y->Draw();
	    
      }
      
    }
  

}
//------------------------------------------------------------------------
void DrawCFDvsLED()
{
  Int_t runNumber=1098;
  
  Int_t npeaks = 20;
  Int_t sigma=2.;
  Char_t buf1[10], buf2[10];
  
  for (Int_t i=0; i<1; i++)
    {
      //   c1->cd(i+1);
      sprintf(buf1,"T0_C_%i_CFD",i+1);
      sprintf(buf2,"CFD_LED%i",i+1);
      
      TH2F *qtc_cfd = (TH2F*) gFile->Get(buf2);
      TH1F *cfd = (TH1F*) gFile->Get(buf1);
      //       cfd->Draw();
      TSpectrum *s = new TSpectrum(2*npeaks,1);
      Int_t nfound = s->Search(cfd,sigma," ",0.05);
      cout<<"Found "<<nfound<<" peaks sigma "<<sigma<<endl;;
      if(nfound!=0){
	Double_t max=0.0; 
	Double_t tabmax[2] = {0.0, 0.0};
	Float_t *xpeak = s->GetPositionX();
	for(Int_t k=0; k<1 ;k++)
	  {
	    Float_t xp = xpeak[k];
	    Int_t xbin = cfd->GetXaxis()->FindBin(xp);
	    Float_t yp = cfd->GetBinContent(xbin);
	    cout<<"xbin = "<<xbin<<"\txpeak = "<<xpeak[k]<<"\typeak = "<<yp<<endl;
	  }
	Float_t hmin=xp-10*sigma;
	Float_t hmax =xp+10*sigma;
	cout<<hmin<< " "<<hmax<<endl;
	//	    cfd->GetXaxis()->SetRange(hmin,hmax);
	//	    TF1 *g1 = new TF1("g1", "gaus", hmin, hmax);
	// cfd->Fit("g1","R");
	Int_t hminbin=  qtc_cfd->GetXaxis()->GetFirst();
	Int_t hmaxbin=  qtc_cfd->GetXaxis()->GetLast();
	Int_t nbins= qtc_cfd->GetXaxis()->GetNbins();
	cout<<"  qtc_cfd "<<hminbin<<" "<<hmaxbin<<" "<<nbins<<endl;
	//  qtc_cfd->Draw();
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
	    if(ip%20 != 0 ) {
	      if (pr_y->GetBinContent(ip) !=0)
		yg +=pr_y->GetBinContent(ip);
	      //	cout<<ng<<" "<<pr_y->GetBinContent(ip)<<" "<<yg<<endl;
	      ng++;}
	    else {
	      xx[ip/20] = Float_t (ip);
	      yy[ip/20] = yg/ng;
	      yg=0;
	      ng=0;
	      cout<<ip<<" "<<ip/20<<" "<< xx[ip/20]<<" "<< yy[ip/20]<<endl;
	    }
	  }
	TH2F *hr = new TH2F("hr"," ",np,0,nbins, np, hmin, hmax);
	hr->Draw();
	TGraph *gr = new TGraph(np,xx,yy);
	gr->SetMinimum(hmin);
	gr->SetMaximum(hmax);
	gr->Draw("P");
	delete [] xx;
	delete [] yy;
	
	// pr_y->Rebin(10);
	//  pr_y->Draw();
	
      }
      
    }
  

}
