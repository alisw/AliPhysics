#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TError.h>
#include <TROOT.h>
#include <TKey.h>
#include <TH2.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TEnv.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TTree.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliTriggerInput.h"
#include "AliTriggerConfiguration.h"
#endif

Int_t DrawTrendingV0QA(TString mergedTrendFile ="trending.root")
{
	if(!mergedTrendFile)
	{
		printf("Cannot open merged trend file with V0 QA");
		return 1;
	}
	char outfilename[200]="ProductionQA.hist.root";
	TString plotDir(".");
	
	Int_t runNumber=0, NumberVoieOff=0, numberBadOffset=0;
	TFile*fin=TFile::Open(mergedTrendFile.Data());
	if(!fin)
	{
		Printf("ERROR: trending file not found. Exiting ...\n");
		return -1;
	}
	TTree*ttree=(TTree*)fin->Get("trending");
	if(!ttree)
	{
		printf("Invalid trending tree");
		return 2;
	}
	Float_t TimesA=0.,TimesC=0., BB_BG=0.,BB_EE=0.,AdcA=0.;
    	Float_t AdcC=0.,MultA=0.,MultC=0.;
    	
   	
   	ttree->SetBranchAddress("run",&runNumber);
   	ttree->SetBranchAddress("TimesA",&TimesA);
	ttree->SetBranchAddress("TimesC",&TimesC);
	ttree->SetBranchAddress("BB_BG",&BB_BG);
  	ttree->SetBranchAddress("BB_EE",&BB_EE);
  	ttree->SetBranchAddress("AdcA" ,&AdcA );
  	ttree->SetBranchAddress("AdcC" ,&AdcC );
  	ttree->SetBranchAddress("MultA",&MultA);
  	ttree->SetBranchAddress("MultC",&MultC);
  	ttree->SetBranchAddress("NumberVoieOff",&NumberVoieOff);
  	ttree->SetBranchAddress("numberBadOffset",&numberBadOffset);
  	
  	Int_t nRuns=ttree->GetEntries();
  	TList list;
  	
  	TH1F * hTimeA = new TH1F("hTimeA","BB Leading time;;Time (ns)",nRuns,-0.5,nRuns-0.5);
  	TH1F * hTimeC = new TH1F("hTimeC","BB Leading time;;Time (ns)",nRuns,-0.5,nRuns-0.5);
  	TH1F * hBB_BG = new TH1F("hBB_BG","Trigger ratio",nRuns,-0.5,nRuns-0.5);
  	TH1F * hBB_EE = new TH1F("hBB_EE","Trigger ratio",nRuns,-0.5,nRuns-0.5);
  	TH1F * hAdcA = new TH1F("hAdcA","Average Charge",nRuns,-0.5,nRuns-0.5);
  	TH1F * hAdcC = new TH1F("hAdcC","Average Charge",nRuns,-0.5,nRuns-0.5);
  	TH1F * hMultA = new TH1F("hMultA","Average number of Fired cell",nRuns,-0.5,nRuns-0.5);
  	TH1F * hMultC = new TH1F("hMultC","Average number of Fired cell",nRuns,-0.5,nRuns-0.5);
  	TH1F * hNumberVoieOff=new TH1F("hNumberVoieOff","Number of chanel off",nRuns,-0.5,nRuns-0.5);
  	hNumberVoieOff->SetMaximum(70);
  	TH1F * hNumberBadOffset=new TH1F("hNumberBadOffset","Number of pdestal",nRuns,-0.5,nRuns-0.5);
  	hNumberBadOffset->SetMaximum(70);
  	
  	list.Add(hTimeA);
  	list.Add(hTimeC);
  	list.Add(hBB_BG);
  	list.Add(hBB_EE);
  	list.Add(hAdcA );
  	list.Add(hAdcC );
  	list.Add(hMultA);
  	list.Add(hMultC);
  	list.Add(hNumberVoieOff);
  	list.Add(hNumberBadOffset);
  	char runlabel[6];
  	
  	
  	for(Int_t irun=0;irun<nRuns;irun++)
  	{
  		ttree->GetEntry(irun);
  		sprintf(runlabel,"%i",runNumber);
  		
  		hTimeA->SetBinContent(irun+1,TimesA);
  		hTimeA->GetXaxis()->SetBinLabel(irun+1,runlabel);
  		
    		hTimeC->SetBinContent(irun+1,TimesC);
    		hTimeC->GetXaxis()->SetBinLabel(irun+1,runlabel);
    		
    		hBB_BG->SetBinContent(irun+1,BB_BG);
    		hBB_BG->GetXaxis()->SetBinLabel(irun+1,runlabel);
    		
    		hBB_EE->SetBinContent(irun+1,BB_EE);
    		hBB_EE->GetXaxis()->SetBinLabel(irun+1,runlabel);
    		
    		hAdcA->SetBinContent(irun+1,AdcA);
    		hAdcA ->GetXaxis()->SetBinLabel(irun+1,runlabel);
    		
    		hAdcC->SetBinContent(irun+1,AdcC);
    		hAdcC ->GetXaxis()->SetBinLabel(irun+1,runlabel);
    		
    		hMultA->SetBinContent(irun+1,MultA);
    		hMultA->GetXaxis()->SetBinLabel(irun+1,runlabel);
    		
    		hMultC->SetBinContent(irun+1,MultC);
    		hMultC->GetXaxis()->SetBinLabel(irun+1,runlabel);
    		
    		hNumberVoieOff->SetBinContent(irun+1,NumberVoieOff);
    		hNumberVoieOff->GetXaxis()->SetBinLabel(irun+1,runlabel);
    		
    		hNumberBadOffset->SetBinContent(irun+1,numberBadOffset);
    		hNumberBadOffset->GetXaxis()->SetBinLabel(irun+1,runlabel);
  	}
  	
  	TFile*fout=new TFile(outfilename,"recreate");
  	fout->cd();
  	list.Write();
  	fout->Close();
  	int maxRun =runNumber;
  	ttree->GetEntry(0);
  	int minRun = runNumber;
  	
  	gStyle->SetOptStat(0);
  	hTimeA->SetMarkerStyle(20);
  	hTimeA->SetMarkerColor(2);
  
  	hTimeC->SetMarkerStyle(20);
  	hTimeC->SetMarkerColor(4);

  	
  	
  	
  	TCanvas * c = new TCanvas("c","Leading time versus run number");
  	hTimeA->GetYaxis()->SetRange(0,10);
  	hTimeA->Draw("P");
  	hTimeA->SetMinimum(TMath::Min(hTimeA->GetMinimum(),hTimeC->GetMinimum())-1.);
  	hTimeA->SetMaximum(TMath::Max(hTimeA->GetMaximum(),hTimeC->GetMaximum())+1.);
  
  	hTimeC->GetYaxis()->SetRange(0,10);
  	hTimeC->Draw("Psame");
  	TLegend * lg = new TLegend(0.8,0.9,1,1);
  	lg->AddEntry(hTimeA,"V0A - 8 ns","p");
  	lg->AddEntry(hTimeC,"V0C","p");
  	lg->Draw("same");
  	float shiftA=8.0;
  	TPave * pavA = new TPave(-0.5,TMath::Max(hTimeA->GetMinimum(),1.5-shiftA),nRuns-0.5,TMath::Min(hTimeA->GetMaximum(),33.5-shiftA),0);
  	pavA->SetFillStyle(3004);
  	pavA->SetFillColor(2);
  	TPave * pavC = new TPave(-0.5,TMath::Max(hTimeC->GetMinimum(),0.5),nRuns-0.5,TMath::Min(hTimeC->GetMaximum(),25.5),0);
  	pavC->SetFillStyle(3005);
  	pavC->SetFillColor(4);
  
  	pavA->Draw("same");
  	pavC->Draw("same");

  	
  	c->Print(Form("%s/QA_Resume_%d_%d.pdf(",plotDir.Data(),minRun,maxRun));
  	c->Write();
  
  	TCanvas * c2 = new TCanvas("c2","Trigger ratios");
  	c2->SetGridy();
  
  	hBB_BG->SetMarkerStyle(20);
  	hBB_BG->SetMarkerColor(2);
  	hBB_EE->SetMarkerStyle(20);
  	hBB_EE->SetMarkerColor(4);
  	
  	hBB_BG->Draw("P");
  	hBB_EE->Draw("Psame");
  	TLegend * lg2 = new TLegend(0.8,0.9,1,1);
  	lg2->AddEntry(hBB_BG,"BG / BB","p");
  	lg2->AddEntry(hBB_EE,"EE / BB","p");
  	lg2->Draw("same");
  	
  	c2->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  	c2->Write();
  	
  	
  	TCanvas * c3 = new TCanvas("c3","Average Charge");
  	c3->SetGridy();
  	
  	hAdcA->SetMarkerStyle(20);
  	hAdcA->SetMarkerColor(2);
  	hAdcC->SetMarkerStyle(20);
  	hAdcC->SetMarkerColor(4);
  	hAdcA->SetMinimum(0);
  	hAdcA->SetMaximum(100);
  	
  	hAdcA->Draw("P");
  	hAdcC->Draw("Psame");
  	TLegend * lg3 = new TLegend(0.8,0.9,1,1);
  	lg3->AddEntry(hAdcA,"V0A","p");
  	lg3->AddEntry(hAdcC,"V0C","p");
  	lg3->Draw("same");
  	
  	c3->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  	c3->Write();
  	
  
  	TCanvas * c4 = new TCanvas("c4","Average number of cell");
  	c4->SetGridy();
  	
  	hMultA->SetMarkerStyle(20);
  	hMultA->SetMarkerColor(2);
  	hMultC->SetMarkerStyle(20);
  	hMultC->SetMarkerColor(4);
  	hMultA->SetMinimum(0);
  	hMultA->SetMaximum(32);
  	
  	hMultA->Draw("P");
  	hMultC->Draw("Psame");
  	TLegend * lg4 = new TLegend(0.8,0.9,1,1);
  	lg4->AddEntry(hMultA,"V0A","p");
  	lg4->AddEntry(hMultC,"V0C","p");
  	lg4->Draw("same");
  	
  	c4->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  	c4->Write();

  	TCanvas * c5 = new TCanvas("c5","");
  	c5->cd();
  	hNumberVoieOff->Draw();
  	TCanvas * c6 = new TCanvas("c6","");
  	c6->cd();
  	hNumberBadOffset->Draw();
	c5->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	c6->Print(Form("%s/QA_Resume_%d_%d.pdf)",plotDir.Data(),minRun,maxRun));
  	c5->Write();
	return 0;
}
