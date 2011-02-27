#include "TMath.h"
#include "TList.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TCanvas.h"
#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include "TStyle.h"
#include "TGrid.h"
 
 
class TList;
class TH2F;
class TH1F;
class TH1D;
class TProfile;
class TFile;
class TCanvas;
class TStyle;
class TF1;
class TGrid;





Double_t convolution(Double_t* x,Double_t *par)
{
	  //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 200.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) 
      {
        	 xx = xlow + (i-.5) * step;
         	fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         	sum += fland * TMath::Gaus(x[0],xx,par[3]);

         	xx = xupp - (i-.5) * step;
         	fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         	sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}




void GetGainModuleLevel(TString filename,Bool_t normal=1,Int_t ntofit=500, Bool_t grid=0)
{
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1,0);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1);	
	if(grid)
		TGrid::Connect("alien://");
	
	TFile* file_data=TFile::Open(filename);
	if(!file_data)
		return;
	TList *listin=0x0;
	listin=(TList*)file_data->Get("output");
	if(!listin)
		listin=(TList*)file_data->Get("PWG1dEdxSSDQA/output");
	if(!listin)	
		listin=(TList*)file_data->Get("PWG1dEdxSSDQA/SSDdEdxQA");
	if(!listin)	
		listin=(TList*)file_data->Get("SSDdEdxQA");
	if(!listin)	
		return;
	TH2F* fHistQ=0x0;
	if (normal)		
		fHistQ=(TH2F*)listin ->FindObject("QACharge");
	else
		fHistQ=(TH2F*)listin ->FindObject("QAChargeCorrected");
	if(!fHistQ)
		return;
	TH2F* fHistCR=(TH2F*)listin ->FindObject("QAChargeRatio");
	if(!fHistCR)
		return;
	TList *listout1=new TList();
	
	TList *listout2=new TList();
	
	TH1F* fHistMPVs=new TH1F("HistMPVS","HistMPVs;MPV;N",75,70,95);
	fHistMPVs->SetDirectory(0);
	listout2->Add(fHistMPVs);
	
	TH1F* fHistSL=new TH1F("HistSL","#sigma_{Landau};#sigma_{Landau};N",40,0,16);
	fHistSL->SetDirectory(0);
	listout2->Add(fHistSL);
	
	TH1F* fHistSG=new TH1F("HistSG","#sigma_{Gaus};#sigma_{Gaus};N",40,0,16);
	fHistSG->SetDirectory(0);
	listout2->Add(fHistSG);
	
	TH1F* fHistCRmean=new TH1F("HistCRmean","HistCRmean;CRmean;N",200,-1,1);
	fHistCRmean->SetDirectory(0);
	listout2->Add(fHistCRmean);
	
	TH1F* fHistCRRMS=new TH1F("HistCRRMS","HistCRRMS;CRRMS;N",100,0,1);
	fHistCRRMS->SetDirectory(0);
	listout2->Add(fHistCRRMS);
	
	TH1F* fHistGainP=new TH1F("HistGainP","HistGainP;CRGainPcorr;N",120,0.5,2.0);
	fHistGainP->SetDirectory(0);
	listout2->Add(fHistGainP);
	
	TH1F* fHistGainN=new TH1F("HistGainN","HistGainN;CRGainNcorr;N",120,0.5,2.0);
	fHistGainN->SetDirectory(0);
	listout2->Add(fHistGainN);
	
	TH1F *fMPVGraph = new TH1F("MPVgraph","MPVgraph;Module number;MPV",1698,-0.5,1697.5);
	fMPVGraph->SetMarkerColor(kRed);
	fMPVGraph->SetMarkerStyle(22);
	listout2->Add(fMPVGraph);
	
	TH1F *fCRmeanGraph = new TH1F("CRmeangraph","CRmeangraph;Module number;MPV",1698,-0.5,1697.5);
	fCRmeanGraph->SetMarkerColor(kBlue);
	fCRmeanGraph->SetMarkerStyle(23);
	listout2->Add(fCRmeanGraph);
	
	Float_t gainP[1698];
	Float_t gainN[1698];
	Float_t mpv[1698];
	Int_t flag[1698];	
	
	ofstream outfiletxt;
	outfiletxt.open("gain.txt");
	 outfiletxt.width(10) ;
	outfiletxt.setf(outfiletxt.left);
	outfiletxt<<"MODULE"<<"\t";
	outfiletxt.width(10);
	outfiletxt.setf(outfiletxt.left);
	outfiletxt<<"FLAG"<<"\t";
	outfiletxt.width(10) ;
	outfiletxt.setf(outfiletxt.left);
	outfiletxt<<"GainPcorr"<<"\t";
	outfiletxt.width(10) ;
	outfiletxt.setf(outfiletxt.left);
	outfiletxt<<"GainNcorr"<<"\t";
	outfiletxt.width(10) ;
	outfiletxt.setf(outfiletxt.left);
	outfiletxt<<"MPV"<<endl;
	
	
	ofstream outfiletxtbad;
	outfiletxtbad.open("badModules.txt");
	
	
	
	
	for (int i =0;i<1698;i++)
	{
		cout<<i<<endl;
		TString tmpQ("Q");
		tmpQ+=i;
		TString tmpCR("CR");
		tmpCR+=i;
		TH1D* fHist1DCR= fHistCR->ProjectionY(tmpCR,i+1,i+1);
		Double_t mean=fHist1DCR->GetMean();
		if(!(TMath::Abs(mean)<1.0)||fHist1DCR->GetEntries()<10)
		{		
			flag[i]=-2;
			gainN[i]=1.0;
			gainP[i]=1.0;
			mpv[i]=1.0;
			continue;
		}
		fHistCRmean->Fill(mean);
		fHistCRRMS->Fill(fHist1DCR->GetRMS());
		gainN[i]=1.0/(1.0+mean);
		gainP[i]=1.0/(1.0-mean);
		fHistGainP->Fill(gainP[i]);
		fHistGainN->Fill(gainN[i]);
		fCRmeanGraph->SetBinContent(i+1,mean);
		fCRmeanGraph->SetBinError(i+1,fHist1DCR->GetRMS());
		
		
		TH1D* fHist1DQ=fHistQ->ProjectionY(tmpQ,i+1,i+1);
		 fHist1DQ->SetDirectory(0);
		 listout1->Add(fHist1DQ);
		if(fHist1DQ->GetEntries()<ntofit)
		{
			flag[i]=-1;
			mpv[i]=1.0;
			 outfiletxtbad<<"Low statistic \t module= "<<i<<" netries="<<fHist1DQ->GetEntries()<<endl;
			continue;
		}
		else
		{
			tmpQ+="fit";
			Float_t range=fHist1DQ->GetBinCenter(fHist1DQ->GetMaximumBin());
			TF1 *f1 = new TF1(tmpQ,convolution,range*0.45,range*3.0,4);
			f1->SetParameters(7.0,range,1.0,5.5);
			Float_t normalization=fHist1DQ->GetEntries()*fHist1DQ->GetXaxis()->GetBinWidth(2)/f1->Integral(range*0.45,range*3.0);
			f1->SetParameters(7.0,range,normalization,5.5);
			//f1->SetParameters(7.0,range,fHist1DQ->GetMaximum(),5.5);
			f1->SetParNames("sigma Landau","MPV","N","sigma Gaus");
			f1->SetParLimits(0,2.0,100.0);
			f1->SetParLimits(3,0.0,100.0);
			if(fHist1DQ->Fit(tmpQ,"BRQ")==0)
			{
				mpv[i]=f1->GetParameter(1);
				fHistMPVs->Fill(mpv[i]);
				fHistSL->Fill(f1->GetParameter(0));
				fHistSG->Fill(f1->GetParameter(3));
				flag[i]=1;		
				fMPVGraph->SetBinContent(i+1,f1->GetParameter(1));
				fMPVGraph->SetBinError(i+1,f1->GetParError(1));
				if(mpv[i]<75.0)
				{
					outfiletxtbad<<"MPV lower than 75 \t module="<<i<<endl;
					flag[i]=0;
				}	
				if(mpv[i]>100.0)
				{
					outfiletxtbad<<"MPV higher than 100 \t module="<<i<<endl;
					flag[i]=0;
					
				}
				if(f1->GetParError(1)>1.0)
				{
					outfiletxtbad<<"MPV high error on MPV  \t module="<<i<<endl;				
					flag[i]=0;
				}
			}
			else
			{
				mpv[i]=1;
				flag[i]=0;
				 outfiletxtbad<<"BAD FIT \t module="<<i<<endl;
				//continue;
			}	
		}	
	}	
	
	for (int i=0;i<1698;i++)
	{	
		outfiletxt.setf(outfiletxt.scientific);
		outfiletxt.precision(2);	
		 outfiletxt.width(10) ;
		outfiletxt.setf(outfiletxt.left);
		outfiletxt<<i<<"\t";
		 outfiletxt.width(10) ;
		outfiletxt.setf(outfiletxt.left);
		outfiletxt<<flag[i]<<"\t";
		 outfiletxt.width(10) ;
		outfiletxt.setf(outfiletxt.left);
		outfiletxt<<gainP[i]<<"\t";
		 outfiletxt.width(10) ;
		outfiletxt.setf(outfiletxt.left);
		outfiletxt<<gainN[i]<<"\t";
		 outfiletxt.width(10) ;
		outfiletxt.setf(outfiletxt.left);
		outfiletxt<<mpv[i]<<endl;
	}
		 
	TCanvas *c1 = new TCanvas("1","1",1200,800);
	c1->Divide(2,1);
	c1->cd(1);
	fHistQ->Draw("colz");
	c1->cd(2);
	fHistCR->Draw("colz");
	
	
	
	TFile* fout1=TFile::Open("gain_all_fits.root","recreate");
	listout1->Write("output",TObject::kSingleKey);	
	fout1->Close();
	
	TFile* fout2=TFile::Open("gain_control_plots.root","recreate");
	listout2->Write("output",TObject::kSingleKey);	
	fout2->Close();
	
	
	
}
