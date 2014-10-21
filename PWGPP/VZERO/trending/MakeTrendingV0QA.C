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


Int_t MakeTrendingV0QA(TString qafilename,Int_t runNumber,TString ocdbStorage = "raw://",Bool_t IsOnGrid = kFALSE,Bool_t canvasE = kFALSE)
{
	if (!qafilename) 
	{
     		Printf("Error - Invalid input file");
  		return 1;
	}
	gStyle->SetPalette(1);
	
	TString treePostFileName=Form("trending_%i.root",runNumber);
	if(IsOnGrid)
		TGrid::Connect("alien://");
	TFile*fin=TFile::Open(qafilename,"r");

	if(!fin)
	{
		Printf("ERROR: QA output not found. Exiting ...\n");
		return -1;
	}
	else
	{
		Printf("INFO: QA output file %s open. \n",fin->GetName());
	}
	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage(ocdbStorage);
	man->SetRun(runNumber);
    	AliCDBEntry *entry2=0;
    	entry2 = man->Get("GRP/GRP/Data");
    	AliGRPObject* fGRPData=0;
    	if (entry2) 
    	{
      		printf("Found an AliGRPObject in GRP/GRP/Data, reading it\n");
      		fGRPData = dynamic_cast<AliGRPObject*>(entry2->GetObject()); 
     		entry2->SetOwner(0);
   	}
    	TString activeDetList(AliDAQ::ListOfTriggeredDetectors(fGRPData->GetDetectorMask()));
    	TString runType(fGRPData->GetRunType());
    	TString beamType(fGRPData->GetBeamType());
    	TString machineMode(fGRPData->GetMachineMode());
    	TString lhcState(fGRPData->GetLHCState());
    	printf("activeDetList %s\nrunType %s\nbeamType %s\nmachineMode %s\nlhcState %s\n",
    		activeDetList.Data(),runType.Data(),beamType.Data(),
    		machineMode.Data(),lhcState.Data());
    
    	time_t duration = fGRPData->GetTimeEnd() - fGRPData->GetTimeStart();
    
   	if(!activeDetList.Contains("VZERO"))
    	{ 
       		printf("RUN WITH VZERO NOT ACTIVE\n");
      		return 0;
    	}
    	if(!runType.Contains("PHYSICS"))
    	{ 
      		printf("RUN NO PHYSICS\n");
      		return 0;
    	}
    	if(duration<120)
    	{ 
      		printf("RUNS SHORTER THAN 2 MIN\n");
     		return 0;
    	}
    	Float_t TimesA=-9999.,TimesC=-9999., BB_BG=-9999.,BB_EE=-9999.,AdcA=-9999.;
    	Float_t AdcC=-9999.,MultA=-9999.,MultC=-9999.;
    	Int_t NumberVoieOff=0, numberBadOffset=0;
   	TTree *ttree=new TTree("trending","tree of trending variables");

   	ttree->Branch("run",&runNumber,"run/I");
   	ttree->Branch("TimesA",&TimesA,"BB Leading time;;Time (ns)/F");
	ttree->Branch("TimesC",&TimesC,"BB Leading time;;Time (ns)/F");
	ttree->Branch("BB_BG",&BB_BG,"Trigger ratio/F");
  	ttree->Branch("BB_EE",&BB_EE,"Trigger ratio/F");
  	ttree->Branch("AdcA" ,&AdcA ,"Average Charge/F");
  	ttree->Branch("AdcC" ,&AdcC ,"Average Charge/F");
  	ttree->Branch("MultA",&MultA,"Average number of Fired cell/F");
  	ttree->Branch("MultC",&MultC,"Average number of Fired cell/F");
  	ttree->Branch("NumberVoieOff",&NumberVoieOff,"Number of path off/I");
  	ttree->Branch("numberBadOffset",&numberBadOffset,"Number of bad offset /I");
  	
  	char v0QAdirName[20]="VZERO_Performance";
  	TDirectoryFile * v0QAdir=(TDirectoryFile*)fin->Get(v0QAdirName);
  	if(!v0QAdir)
  	{
  		printf("ERROR: VZERO QA directory not present in input file.\n");
  		return -1;
  	}
  	TList *list = (TList*)v0QAdir->Get("QAVZEROHists");
    	if(!list) 
    	{
      		cout << "ERROR: No list found" << endl;
		return -1;
    	}
    	TH2F *hEvents = (TH2F*)list->FindObject("hEvents");
    	TH1F *hAdcNoTimeA = (TH1F*)list->FindObject("hAdcNoTimeA");
    	TH1F *hAdcWithTimeA = (TH1F*)list->FindObject("hAdcWithTimeA");
   	TH1F *hAdcNoTimeC = (TH1F*)list->FindObject("hAdcNoTimeC");
    	TH1F *hAdcWithTimeC = (TH1F*)list->FindObject("hAdcWithTimeC");
    	TH2F *hadcpmtwithtime = (TH2F*)list->FindObject("hadcpmtwithtime");	
    	TH1F *htimepmtA = (TH1F*)list->FindObject("htimepmtA");
    	TH1F *htimepmtC = (TH1F*)list->FindObject("htimepmtC");
    	TH1F *hwidthA = (TH1F*)list->FindObject("hwidthA");
    	TH1F *hwidthC = (TH1F*)list->FindObject("hwidthC");
    	TH1F *hV0ampl = (TH1F*)list->FindObject("hV0ampl");
    	TH2F *htimepmt = (TH2F*)list->FindObject("htimepmt");	
   	TH2F *hwidthpmt = (TH2F*)list->FindObject("hwidthpmt");	
    	TH2F *hadcwidthA = (TH2F*)list->FindObject("hadcwidthA");	
    	TH2F *hadcwidthC = (TH2F*)list->FindObject("hadcwidthC");	
    	TH2F *hAdcTimeA = (TH2F*)list->FindObject("hAdcTimeA");	
    	TH2F *hAdcTimeC = (TH2F*)list->FindObject("hAdcTimeC");	
    	TH2F *htimecorr = (TH2F*)list->FindObject("htimecorr");	
    	TH1F *hV0A = (TH1F*)list->FindObject("hV0a");
    	TH1F *hV0C = (TH1F*)list->FindObject("hV0c");
    	TH1F *hV0multA = (TH1F*)list->FindObject("hV0multA");
    	TH1F *hV0multC = (TH1F*)list->FindObject("hV0multC");
    	TH2F* hVtxXYBB  =(TH2F*) list->FindObject("fhVtxXYBB");
    	TH1F* hVtxZBB   =(TH1F*) list->FindObject("fhVtxZBB");
    	TH2F* hVtxXYBGA =(TH2F*) list->FindObject("fhVtxXYBGA");
    	TH1F* hVtxZBGA  =(TH1F*) list->FindObject("fhVtxZBGA");
    	TH2F* hVtxXYBGC =(TH2F*) list->FindObject("fhVtxXYBGC");
    	TH1F* hVtxZBGC  =(TH1F*) list->FindObject("fhVtxZBGC");
    	
    	float BB = hEvents->GetBinContent(2,2);
    	float EE = hEvents->GetBinContent(1,1);
    	float BGA = hEvents->GetBinContent(3,2);
    	float BGC = hEvents->GetBinContent(2,3);
    	
    	if(hAdcWithTimeA->GetEntries()==0)
    		return 0;
    	
    	{
    		TSpectrum s;
        	float shiftA = 8.;
    		Int_t nPeaksFound = s.Search(htimepmtA);
    		Float_t *peaks = s.GetPositionY();
		Float_t *posiX = s.GetPositionX();
    		Float_t maxY = 0.;
    		Int_t index = -1;
	
    		for(int i=0;i<nPeaksFound;i++) 
    		{
    			if(peaks[i]>maxY && posiX[i]>0.) 
    			{
				maxY = peaks[i];
				index = i;
			}	
    		}
    		Float_t maxX = (index >= 0) ? s.GetPositionX()[index] : -11111;	
		
		TF1 *fgaus = new TF1("gausbbbb","gaus",maxX-1.,maxX+1.);
    		htimepmtA->Fit(fgaus,"","",maxX-1.,maxX+1.);
    		TimesA=fgaus->GetParameter(1)-shiftA;
    		delete fgaus;
    	}
    	{
   		TSpectrum s;
    		Int_t nPeaksFound = s.Search(htimepmtC);
		Float_t *peaks = s.GetPositionY();
		Float_t maxY = 0.;
		Int_t index = -1;
    		for(int i=0;i<nPeaksFound;i++) 
    		{
			if(peaks[i]>maxY) 
			{
			maxY = peaks[i];
			index = i;
			}
		}	
    		Float_t maxX = (index >= 0) ? s.GetPositionX()[index] : -11111;	
		TF1 *fgaus = new TF1("gausffff","gaus",maxX-1.,maxX+1.);
    		htimepmtC->Fit(fgaus,"","",maxX-1.,maxX+1.);
    		TimesC=fgaus->GetParameter(1);
    		delete fgaus;
    	}
    	if(BB) 
    	{
      		BB_BG=(BGA+BGC)/BB;
      		BB_EE=EE/BB;
    	}else
    	{
      		BB_BG=0;
      		BB_EE=0;
    	}
    	
    	MultA=hV0A->GetMean();
    	MultC=hV0C->GetMean();
    
    	AdcA=hAdcWithTimeA->GetMean();
    	AdcC=hAdcWithTimeC->GetMean();
    	
	double valBin=0;
	TH1D*hadcXFull=hadcpmtwithtime->ProjectionX("hadcXFull",1,hadcpmtwithtime->GetYaxis()->GetLast());
	TH1D*hadcX=hadcpmtwithtime->ProjectionX("hadcX",10,20);
	for(Int_t i=0;i<64;i++)
	{
		valBin=hadcXFull->GetBinContent(i+1);
		if(valBin==0)
			NumberVoieOff++;
		valBin=hadcX->GetBinContent(i+1);
		if(valBin==0)
			numberBadOffset++;
	}
	    	
    	TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
    	ttree->Fill();
	trendFile->cd();
	ttree->Write();
	trendFile->Close();
	
    	if(canvasE)
    	{
   	
    		TCanvas * cOut = new TCanvas("cOut",Form("Run %d",runNumber));
    		cOut->Divide(2,2);
    		cOut->cd(1); cOut->GetPad(1)->SetLogy();
    		hAdcNoTimeA->Draw("l");
    		hAdcWithTimeA->Draw("same"); hAdcWithTimeA->SetLineColor(2);
    
     		cOut->cd(2); cOut->GetPad(2)->SetLogy();
    		hAdcNoTimeC->Draw("l");
    		hAdcWithTimeC->Draw("same"); hAdcWithTimeC->SetLineColor(2);
    	
       		cOut->cd(3); cOut->GetPad(3)->SetLogz();
    		hadcpmtwithtime->Draw("colz");
    	
    		cOut->cd(4); cOut->GetPad(4)->SetLogz();
    		hEvents->Draw("colz text");
    	
    		cOut->Print(Form("QA_Run_%d.pdf(",runNumber));
	
	
    		cOut->cd(1); cOut->GetPad(1)->SetLogy();
    		htimepmtA->GetXaxis()->SetRangeUser(-25.,25.); htimepmtA->Draw();
    	
    		cOut->cd(2); cOut->GetPad(2)->SetLogy();
    		htimepmtC->GetXaxis()->SetRangeUser(-25.,25.); htimepmtC->Draw();
    	
       		cOut->cd(3); cOut->GetPad(3)->SetLogy();cOut->GetPad(3)->SetLogz(0);
    		hwidthA->GetXaxis()->SetRangeUser(0.,50.); hwidthA->Draw();
    	
    		cOut->cd(4); cOut->GetPad(4)->SetLogy();cOut->GetPad(4)->SetLogz(0);
    		hwidthC->GetXaxis()->SetRangeUser(0.,50.); hwidthC->Draw();
    	
    		cOut->Print(Form("QA_Run_%d.pdf",runNumber));
    		
    	
    		cOut->cd(1); cOut->GetPad(1)->SetLogy(0);cOut->GetPad(1)->SetLogz();
    		htimepmt->Draw("colz");
    	
    		cOut->cd(2); cOut->GetPad(2)->SetLogy(0);cOut->GetPad(2)->SetLogz();
    		hwidthpmt->GetYaxis()->SetRangeUser(0.,50.); hwidthpmt->Draw("colz");
    		
    		cOut->cd(3); cOut->GetPad(3)->SetLogy(0);cOut->GetPad(3)->SetLogz();
    		hadcwidthA->GetYaxis()->SetRangeUser(0.,50.); hadcwidthA->Draw("colz");
    	
    		cOut->cd(4); cOut->GetPad(4)->SetLogy(0);cOut->GetPad(4)->SetLogz();
    		hadcwidthC->GetYaxis()->SetRangeUser(0.,50.); hadcwidthC->Draw("colz");
    	
    		cOut->Print(Form("QA_Run_%d.pdf",runNumber));
    		
	
    		cOut->cd(1); cOut->GetPad(1)->SetLogy(0);cOut->GetPad(1)->SetLogz();
    		hAdcTimeA->Draw("colz");
	
    		cOut->cd(2); cOut->GetPad(2)->SetLogy(0);cOut->GetPad(2)->SetLogz();
    		hAdcTimeC->Draw("colz");
 	
    		cOut->cd(3); cOut->GetPad(3)->SetLogy(); cOut->GetPad(3)->SetLogz(0);
    		hV0ampl->Draw();
	
    		cOut->cd(4); cOut->GetPad(4)->SetLogy(0); cOut->GetPad(4)->SetLogz(0);
    		htimecorr->Draw("colz");
    	
    		cOut->Print(Form("QA_Run_%d.pdf",runNumber));
	
		
    	
    		cOut->cd(1);  cOut->GetPad(1)->SetLogy(1);cOut->GetPad(1)->SetLogz(0);
    		hV0A->GetXaxis()->SetRangeUser(0.,33.);hV0A->Draw();
    	
    		cOut->cd(2); cOut->GetPad(2)->SetLogy(1);cOut->GetPad(2)->SetLogz(0);
    		hV0C->GetXaxis()->SetRangeUser(0.,33.);hV0C->Draw();
    	
    		cOut->cd(3); cOut->GetPad(3)->SetLogy(); cOut->GetPad(3)->SetLogz(0);
    		hV0multA->Draw();
		
    		cOut->cd(4); cOut->GetPad(4)->SetLogy(); cOut->GetPad(3)->SetLogz(0);
    		hV0multC->Draw();
    	
    		cOut->Print(Form("QA_Run_%d.pdf",runNumber));
	
		
    		cOut->Clear();
    		cOut->Divide(2,3);
    	
    		cOut->cd(1);  cOut->GetPad(1)->SetLogy(0);cOut->GetPad(1)->SetLogz(1);
    		hVtxXYBB->Draw("colz");
    	
    		cOut->cd(2); cOut->GetPad(2)->SetLogy(1);cOut->GetPad(2)->SetLogz(0);
    		hVtxZBB->Draw();
    	
      		cOut->cd(3); cOut->GetPad(3)->SetLogy(0); cOut->GetPad(3)->SetLogz(1);
    		hVtxXYBGA->Draw("colz");
    	
    		cOut->cd(4); cOut->GetPad(4)->SetLogy(); cOut->GetPad(3)->SetLogz(0);
    		hVtxZBGA->Draw();
    	
    		cOut->cd(5); cOut->GetPad(5)->SetLogy(0); cOut->GetPad(5)->SetLogz(1);
    		hVtxXYBGC->Draw("colz");
    	
    		cOut->cd(6); cOut->GetPad(6)->SetLogy(); cOut->GetPad(6)->SetLogz(0);
    		hVtxZBGC->Draw();
    		
    		cOut->Print(Form("QA_Run_%d.pdf)",runNumber));
    		delete cOut;
    					
 	}
	delete v0QAdir;
 	return 0;
}
