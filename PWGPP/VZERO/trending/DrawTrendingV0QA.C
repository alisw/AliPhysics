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
  Int_t badrun=0;
  Int_t goodrun=0;
  Double_t maxVoieOff=0;
  Double_t maxBadOffset=0;
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
  Float_t AdcC=0.,MultA=0.,MultC=0.,AdcAError=0.,AdcCError=0.,ChargeCh46=0.,ChargeAllNo46=0.;;
  Float_t TriggerEff_CVLN=0.,TriggerEff_CVHN=0.,TriggerEff_CVHN2=0.;
  Float_t TriggerEff_CVLN_Error=0.,TriggerEff_CVHN_Error=0.,TriggerEff_CVHN2_Error=0.;
  Float_t PMTEdges[64]={0.};
  Float_t PMTEdgesError[64]={0.};
  Float_t PMTEdges[64]={0.};
  Float_t PMTEdgesError[64]={0.};
  Float_t PMmeanAdc[64]={0.};
  Float_t RingmeanAdc[8]={0.};
  Float_t VOXmeanAdc[2]={0.};
  Float_t VOmeanAdc=0.;
  Bool_t isPP;
  Int_t higtVoltage=0,invalidInput=0,qaNotFound=0,v0active=0,v0qaNotfound=0,noEntries=0;
  
  TH1I*HigtVoltage= new TH1I("HigtVoltage","Run with higt voltage is off",0,0,1);
  TH1I*InvalidInput= new TH1I("InvalidInput","run with invalid input",0,0,1);
  TH1I*QaNotFound= new TH1I("QaNotFound","QA not found",0,0,1);
  TH1I*V0active= new TH1I("V0active","v0 not active",0,0,1);
  TH1I*V0qaNotfound= new TH1I("V0qaNotfound","V0 QA not found",0,0,1);
  TH1I*NoEntrie= new TH1I("NoEntrie","no entries in qa file",0,0,1);
  
  ttree->SetBranchAddress("HigtVoltage",&higtVoltage);
  ttree->SetBranchAddress("invalidInput",&invalidInput);
  ttree->SetBranchAddress("qaNotFound",&qaNotFound);
  ttree->SetBranchAddress("v0active",&v0active);
  ttree->SetBranchAddress("v0qaNotfound",&v0qaNotfound);
  ttree->SetBranchAddress("noEntries",&noEntries);
  ttree->SetBranchAddress("run",&runNumber);
  ttree->SetBranchAddress("isPP",&isPP);
  ttree->GetEntry(0);
  if(isPP)
    {
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
      ttree->SetBranchAddress("ChargeCh46",&ChargeCh46);
      ttree->SetBranchAddress("ChargeAllNo46",&ChargeAllNo46);
      for(int i = 0; i < 64; i++)
	  ttree->SetBranchAddress(Form("PMmeanAdc[%d]",i),&PMmeanAdc[i]);
      for(int i=0;i<8;i++)
	ttree->SetBranchAddress(Form("RingmeanAdc[%d]",i),&RingmeanAdc[i]);
      ttree->SetBranchAddress("VOXmeanAdc[0]",&VOXmeanAdc[0]);
      ttree->SetBranchAddress("VOXmeanAdc[1]",&VOXmeanAdc[1]);
      ttree->SetBranchAddress("VOmeanAdc"    ,&VOmeanAdc);
      
      Int_t nRuns=ttree->GetEntries();
      TList list;
      
      TH1F * hTimeA = new TH1F("hTimeA","BB Leading time;;Time (ns)",nRuns,-0.5,nRuns-0.5);
      TH1F * hTimeC = new TH1F("hTimeC","BB Leading time;;Time (ns)",nRuns,-0.5,nRuns-0.5);
      TH1F * hBB_BG = new TH1F("hBB_BG","Trigger ratio",nRuns,-0.5,nRuns-0.5);
      TH1F * hBB_EE = new TH1F("hBB_EE","Trigger ratio",nRuns,-0.5,nRuns-0.5);
      TH1F * hAdcA  = new TH1F("hAdcA","Average Charge",nRuns,-0.5,nRuns-0.5);
      TH1F * hAdcC  = new TH1F("hAdcC","Average Charge",nRuns,-0.5,nRuns-0.5);
      TH1F * hMultA = new TH1F("hMultA","Average number of Fired cell",nRuns,-0.5,nRuns-0.5);
      TH1F * hMultC = new TH1F("hMultC","Average number of Fired cell",nRuns,-0.5,nRuns-0.5);
      TH1F * hNumberVoieOff   =new TH1F("hNumberVoieOff","Number of chanel off",nRuns,-0.5,nRuns-0.5);
      TH1F * hNumberBadOffset =new TH1F("hNumberBadOffset","Number of pdestal",nRuns,-0.5,nRuns-0.5);
      TH1F *hChargeCh46       = new TH1F("hChargeCh46","Integreted charge of chanel 46",nRuns,-0.5,nRuns-0.5);
      TH1F *hChargeChAllNo46  = new TH1F("hChargeChAllNo46","Integreted charge of all chanel without 46",nRuns,-0.5,nRuns-0.5);
      TH1F *hChargeRatio      = new TH1F("hChargeRatio","Ratio of integreted charge of chanel 46 over all chanel",nRuns,-0.5,nRuns-0.5);
      TH1F **hPMmeanAdc       = new TH1F*[64];
      TH1F **hRingmeanAdc     = new TH1F*[32];
      TH1F **hVOXmeanAdc      = new TH1F*[2];
      for (int i=0;i<64;i++)
	hPMmeanAdc[i]   = new TH1F(Form("hPMmeanAdc[%d]",i)  ,Form("Mean ADC for chanel %d",i),nRuns,-0.5,nRuns-0.5);
      for (int i=0;i<32;i++)
	hRingmeanAdc[i] = new TH1F(Form("hRingmeanAdc[%d]",i),Form("Mean ADC for Ring  %d",i) ,nRuns,-0.5,nRuns-0.5);
      hVOXmeanAdc[0]    = new TH1F("hVOXmeanAdc[0]","Mean ADC for V0C",nRuns,-0.5,nRuns-0.5);
      hVOXmeanAdc[1]    = new TH1F("hVOXmeanAdc[1]","Mean ADC for V0A",nRuns,-0.5,nRuns-0.5);
      TH1F *hVOmeanAdc  = new TH1F("hVOmeanAdc"    ,"mean ADC for V0" ,nRuns,-0.5,nRuns-0.5);
      //PMmeanAdc[64] RingmeanAdc[8] VOXmeanAdc[2] VOmeanAdc
      hNumberBadOffset->SetMaximum(70);
      hNumberVoieOff->SetMaximum(70);
      
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
      list.Add(hChargeCh46);
      list.Add(hChargeChAllNo46);
      list.Add(hChargeRatio);
      char runlabel[6];
      
      
      for(Int_t irun=0;irun<nRuns;irun++)
	{
	  ttree->GetEntry(irun);
	  sprintf(runlabel,"%i",runNumber);
	  //	  if(runNumber==130601)
	  //	    printf("%d\t%d\t%d\n",qaNotFound,invalidInput,runNumber);
	  if(invalidInput==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);
	      
	      InvalidInput->SetBinContent(badrun+1,invalidInput);
	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);
	      
	    }
	  else if(qaNotFound==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);
	      
	      QaNotFound->SetBinContent(badrun+1,qaNotFound);
	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);
	      
	    }
	  else if(v0active==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);
	      
	      V0active->SetBinContent(badrun+1,v0active);
	      
	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);
	      
	    }
	  else if(v0qaNotfound==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);
	      
	      V0qaNotfound->SetBinContent(badrun+1,v0qaNotfound);
	      
	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);
	      
	    }
	  else if(noEntries==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);
	      
	      NoEntrie->SetBinContent(badrun+1,noEntries);
	      
	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);
	      
	    }
	  else 
	    {
	      hTimeA->SetBins(goodrun+1,0,goodrun+1);	
	      hTimeA->SetBinContent(goodrun+1,TimesA);
	      hTimeA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      
	      hTimeC->SetBins(goodrun+1,0,goodrun+1);
	      hTimeC->SetBinContent(goodrun+1,TimesC);
	      hTimeC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hBB_BG->SetBins(goodrun+1,0,goodrun+1);
	      hBB_BG->SetBinContent(goodrun+1,BB_BG);
	      hBB_BG->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      
	      hBB_EE->SetBins(goodrun+1,0,goodrun+1);
	      hBB_EE->SetBinContent(goodrun+1,BB_EE);
	      hBB_EE->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      
	      hAdcA->SetBins(goodrun+1,0,goodrun+1);
	      hAdcA->SetBinContent(goodrun+1,AdcA);
	      hAdcA ->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      
	      hAdcC->SetBins(goodrun+1,0,goodrun+1);
	      hAdcC->SetBinContent(goodrun+1,AdcC);
	      hAdcC ->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      
	      hMultA->SetBins(goodrun+1,0,goodrun+1);
	      hMultA->SetBinContent(goodrun+1,MultA);
	      hMultA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      
	      hMultC->SetBins(goodrun+1,0,goodrun+1);
	      hMultC->SetBinContent(goodrun+1,MultC);
	      hMultC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      
	      hNumberVoieOff->SetBins(goodrun+1,0,goodrun+1);
	      hNumberVoieOff->SetBinContent(goodrun+1,NumberVoieOff);
	      hNumberVoieOff->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      
	      hNumberBadOffset->SetBins(goodrun+1,0,goodrun+1);
	      hNumberBadOffset->SetBinContent(goodrun+1,numberBadOffset);
	      hNumberBadOffset->GetXaxis()->SetBinLabel(goodrun+1,runlabel);

	      maxVoieOff=NumberVoieOff;
	      maxBadOffset=numberBadOffset;

	      hChargeCh46->SetBins(goodrun+1,0,goodrun+1);
	      hChargeCh46->SetBinContent(goodrun+1,ChargeCh46);
	      hChargeCh46->GetXaxis()->SetBinLabel(goodrun+1,runlabel);

	      hChargeChAllNo46->SetBins(goodrun+1,0,goodrun+1);
	      hChargeChAllNo46->SetBinContent(goodrun+1,ChargeAllNo46);
	      hChargeChAllNo46->GetXaxis()->SetBinLabel(goodrun+1,runlabel);

	      Double_t ratio=-5;
	      if (ChargeAllNo46!=0)
		ratio=ChargeCh46/ChargeAllNo46;

	      hChargeRatio->SetBins(goodrun+1,0,goodrun+1);
	      hChargeRatio->SetBinContent(goodrun+1,ratio);
	      hChargeRatio->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      /*
      list.Add(hChargeCh46);
      list.Add(hChargeChAllNo46);
      list.Add(hChargeRatio);
*/
	      for(int i=0;i<64;i++)
		{
		  hPMmeanAdc[i]->SetBins(goodrun+1,0,goodrun+1);
		  hPMmeanAdc[i]->SetBinContent(goodrun+1,PMmeanAdc[i]);
		  hPMmeanAdc[i]->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		}
	      for(int i=0;i<8;i++)
		{
		  hRingmeanAdc[i]->SetBins(goodrun+1,0,goodrun+1);
		  hRingmeanAdc[i]->SetBinContent(goodrun+1,RingmeanAdc[i]);
		  hRingmeanAdc[i]->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		}
	      for(int i=0;i<2;i++)
		{
		  hVOXmeanAdc[i]->SetBins(goodrun+1,0,goodrun+1);
		  hVOXmeanAdc[i]->SetBinContent(goodrun+1,VOXmeanAdc[i]);
		  hVOXmeanAdc[i]->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		}
	      hVOmeanAdc->SetBins(goodrun+1,0,goodrun+1);
	      hVOmeanAdc->SetBinContent(goodrun+1,VOmeanAdc);
	      hVOmeanAdc->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      goodrun++;
	    }
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
      if(hTimeA->GetEntries()!=0)
	{
	  hTimeA->SetMinimum(TMath::Min(hTimeA->GetMinimum(),hTimeC->GetMinimum())-1.);
	  hTimeA->SetMaximum(TMath::Max(hTimeA->GetMaximum(),hTimeC->GetMaximum())+1.);
	}
      hTimeC->GetYaxis()->SetRange(0,10);
      hTimeC->Draw("Psame");
      TLegend * lg = new TLegend(0.8,0.9,1,1);
      lg->AddEntry(hTimeA,"V0A - 8 ns","p");
      lg->AddEntry(hTimeC,"V0C","p");
      lg->Draw("same");
      float shiftA=8.0;
      TPave * pavA = new TPave(0,TMath::Max(hTimeA->GetMinimum(),1.5-shiftA),goodrun,TMath::Min(hTimeA->GetMaximum(),33.5-shiftA),0);
      pavA->SetFillStyle(3004);
      pavA->SetFillColor(2);
      TPave * pavC = new TPave(0,TMath::Max(hTimeC->GetMinimum(),0.5),goodrun,TMath::Min(hTimeC->GetMaximum(),25.5),0);
      pavC->SetFillStyle(3005);
      pavC->SetFillColor(4);
      
      pavA->Draw("same");
      pavC->Draw("same");
      
	    
      c->Print(Form("%s/QA_Resume_%d_%d.pdf(",plotDir.Data(),minRun,maxRun));
      c->SaveAs(Form("%s/V0QA__Leading_time.png",plotDir.Data()));
      c->Write();
      
      TCanvas * c2 = new TCanvas("c2","Trigger ratios");
      c2->SetGridy();
      Double_t max=TMath::Max(hBB_BG->GetMaximum(),hBB_EE->GetMaximum());
      Double_t min=TMath::Min(hBB_BG->GetMinimum(0.000001),hBB_EE->GetMinimum(0.0001));
      printf("max: %.5f\tmin: %.5f\n",max,min);
      if(max==0) max=1.0;
      Double_t relatMean=(max-min)/max;
      hBB_BG->SetMarkerStyle(20);
      hBB_BG->SetMarkerColor(2);
      hBB_EE->SetMarkerStyle(20);
      hBB_EE->SetMarkerColor(4);
      hBB_BG->GetYaxis()->SetRangeUser(min*0.8,max*1.5);
      hBB_BG->Draw("P");
      hBB_EE->Draw("Psame");
      TLegend * lg2 = new TLegend(0.8,0.9,1,1);
      lg2->AddEntry(hBB_BG,"BG / BB","p");
      lg2->AddEntry(hBB_EE,"EE / BB","p");
      lg2->Draw("same");
      printf("diff max - min: %.3f\n",relatMean);
      if(relatMean>0.1)
	c2->SetLogy();
      c2->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      c2->SaveAs(Form("%s/V0QA__Trigger_ratios.png",plotDir.Data()));
      c2->Write();
      
      
      //PMmeanAdc[64] RingmeanAdc[8] VOXmeanAdc[2] VOmeanAdc
      Int_t Color[8]={kRed,kMagenta,kBlue,kCyan,kGreen,kYellow,kOrange,kPink};
      /*     Int_t color=kRed;
      TCanvas * cPMmeanAdc = new TCanvas("cPMmeanAdc","mean ADC PM");
      cPMmeanAdc->SetGridy();
      for(int i=0;i<64;i++)
	{
	  if(!(i%10))
	    color=Color[i/10]+i;
	  hPMmeanAdc[i]->SetMarkerStyle(2);
	  hPMmeanAdc[i]->SetMarkerColor(color-i);
	  hPMmeanAdc[i]->Draw("Psame");
	}
      cPMmeanAdc->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cPMmeanAdc->SaveAs(Form("%s/V0QA__PM_mean_ADC.png",plotDir.Data()));
      cPMmeanAdc->Write();
      */

      printf("toto\n\n");
      TCanvas *cRing  = new TCanvas("cRing","mean ADC Ring");
      cRing->SetGridy();
      TCanvas ** cPMmeanAdc = new TCanvas*[8];//("cPMmeanAdc","mean ADC PM");
      for(int i=0;i<8;i++)
	{
	  cPMmeanAdc[i] = new TCanvas(Form("cPMmeanAdc[%d]",i),Form("mean ADC PM for ring %d",i));
	  cPMmeanAdc[i] ->SetGridy();
	  for(int j=0;j<8;j++)
	    {
	      hPMmeanAdc[j+8*i]->SetMarkerStyle(2);
	      hPMmeanAdc[j+8*i]->SetMarkerColor(Color[j]);
	      hPMmeanAdc[j+8*i]->Draw("Psame");
	    }
	  cPMmeanAdc[i] ->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  cPMmeanAdc[i] ->SaveAs(Form("%s/V0QA__PM_mean_%i_ADC.png",plotDir.Data(),i));
	  cPMmeanAdc[i] ->Write();
	  cRing->cd();
	  hRingmeanAdc[i]->SetMarkerStyle(2);
	  hRingmeanAdc[i]->SetMarkerColor(Color[i]);
	  hRingmeanAdc[i]->Draw("Psame");
	}
      cRing->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cRing->SaveAs(Form("%s/V0QA__Ring_mean_ADC.png",plotDir.Data()));
      cRing->Write();

      TFile*eageV0=new TFile("eageV0.root","recreate");
      for(int i=0;i<64;i++)hPMmeanAdc[i]->Write();
      eageV0->Close();

      TCanvas *cV0X  = new TCanvas("cV0X","mean ADC V0X");
      cV0X->SetGridy();
      hVOXmeanAdc[0]->SetMarkerStyle(2);
      hVOXmeanAdc[0]->SetMarkerColor(kRed);
      hVOXmeanAdc[0]->Draw("P");
      hVOXmeanAdc[1]->SetMarkerStyle(2);
      hVOXmeanAdc[1]->SetMarkerColor(kBlue);
      hVOXmeanAdc[1]->Draw("Psame");
      cV0X->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cV0X->SaveAs(Form("%s/V0QA__V0X_mean_ADC.png",plotDir.Data()));
      cV0X->Write();
      //VOmeanAdc

      TCanvas *cV0  = new TCanvas("cV0","mean ADC V0");
      cV0->SetGridy();
      hVOmeanAdc[0]->SetMarkerStyle(20);
      hVOmeanAdc[0]->SetMarkerColor(kRed);
      hVOmeanAdc[0]->Draw("P");
      cV0->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cV0->SaveAs(Form("%s/V0QA__V0_mean_ADC.png",plotDir.Data()));
      cV0->Write();


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
      c3->SaveAs(Form("%s/V0QA__Average_charge.png",plotDir.Data()));
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
      c4->SaveAs(Form("%s/V0QA__Average_number_of_cell.png",plotDir.Data()));
      c4->Write();
	    
      TCanvas * c5 = new TCanvas("c5","");
      c5->cd();
      hNumberVoieOff->GetYaxis()->SetRangeUser(0,maxVoieOff+1.);
      hNumberVoieOff->Draw();
      c5->SaveAs(Form("%s/V0QA__NumberVoieOff.png",plotDir.Data()));
      TCanvas * c6 = new TCanvas("c6","");
      c6->cd();
      hNumberBadOffset->GetYaxis()->SetRangeUser(0,maxBadOffset+1.);
      hNumberBadOffset->Draw();
      c6->SaveAs(Form("%s/V0QA__Number_Bad_Offset.png",plotDir.Data()));



      TCanvas *CchargeInteger = new TCanvas("CchargeInteger","");
      CchargeInteger->cd();
      hChargeChAllNo46->Draw("P");
      hChargeChAllNo46->SetLineColor(2);
      hChargeChAllNo46->SetMarkerStyle(20);
      hChargeChAllNo46->SetMarkerColor(2);
      hChargeCh46->Draw("Psame");
      hChargeCh46->SetLineColor(1);
      hChargeCh46->SetMarkerStyle(20);
      hChargeCh46->SetMarkerColor(1);
      CchargeInteger->SaveAs(Form("%s/V0QA__Charge_Chanel_46.png",plotDir.Data()));

      TCanvas *CratioCharge  = new TCanvas("CratioCharge","");
      CratioCharge->cd();
      hChargeRatio->Draw("P");
      hChargeRatio->SetLineColor(3);
      hChargeRatio->SetMarkerStyle(20);
      hChargeRatio->SetMarkerColor(3);
      CratioCharge->SaveAs(Form("%s/V0QA__Charge_Chanel_Ratio.png",plotDir.Data()));

      if(HigtVoltage->GetEntries())
	{
	  TCanvas * c7 = new TCanvas("c7","");
	  c7->cd();
	  HigtVoltage->Draw();
	  c7->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c7->SaveAs(Form("%s/V0QA__Hight_Voltage.png",plotDir.Data()));
	}
      if(InvalidInput->GetEntries())
	{
	  TCanvas * c8 = new TCanvas("c8","");
	  c8->cd();
	  InvalidInput->Draw();
	  c8->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c8->SaveAs(Form("%s/V0QA__InvalidInput_%d_%d.png",plotDir.Data(),minRun,maxRun));
	}
      if(QaNotFound->GetEntries())
	{
	  TCanvas * c9 = new TCanvas("c9","");
	  c9->cd();
	  QaNotFound->Draw();
	  c9->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c9->SaveAs(Form("%s/V0QA__QA_Not_Found_%d_%d.png",plotDir.Data()));
	}
      if(V0active->GetEntries())
	{
	  TCanvas * c10 = new TCanvas("c10","");
	  c10->cd();
	  V0active->Draw();
	  c10->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c10->SaveAs(Form("%s/V0QA__V0_not_Active.png",plotDir.Data()));
	}
      if(V0qaNotfound->GetEntries())
	{
	  TCanvas * c11 = new TCanvas("c11","");
	  c11->cd();
	  V0qaNotfound->Draw();
	  c11->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c11->SaveAs(Form("%s/V0QA__v0qa_Not_Found.png",plotDir.Data()));
	}
      if(NoEntrie->GetEntries())
	{
	  TCanvas * c12 = new TCanvas("c12","");
	  c12->cd();
	  NoEntrie->Draw();
	  c12->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c12->SaveAs(Form("%s/V0QA__NoEntrie.png",plotDir.Data()));
	}

      c5->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      c6->Print(Form("%s/QA_Resume_%d_%d.pdf)",plotDir.Data(),minRun,maxRun));
      c5->Write();
    }
  else
    {
      ttree->SetBranchAddress("TimesA",&TimesA);
      ttree->SetBranchAddress("TimesC",&TimesC);
      ttree->SetBranchAddress("BB_BG",&BB_BG);
      ttree->SetBranchAddress("BB_EE",&BB_EE);
      ttree->SetBranchAddress("AdcA" ,&AdcA );
      ttree->SetBranchAddress("AdcAError",&AdcAError);
      ttree->SetBranchAddress("AdcC" ,&AdcC );
      ttree->SetBranchAddress("AdcCError",&AdcCError);
      ttree->SetBranchAddress("MultA",&MultA);
      ttree->SetBranchAddress("MultC",&MultC);
      ttree->SetBranchAddress("TriggerEff_CVLN",&TriggerEff_CVLN);
      ttree->SetBranchAddress("TriggerEff_CVHN",&TriggerEff_CVHN);
      ttree->SetBranchAddress("TriggerEff_CVHN2",&TriggerEff_CVHN2);
      ttree->SetBranchAddress("TriggerEff_CVLN_Error",&TriggerEff_CVLN_Error);
      ttree->SetBranchAddress("TriggerEff_CVHN_Error",&TriggerEff_CVHN_Error);
      ttree->SetBranchAddress("TriggerEff_CVHN2_Error",&TriggerEff_CVHN2_Error);
      for(int i = 0; i < 64; ++i)
	{
	  ttree->SetBranchAddress(Form("PMTEdges[%d]",i),&PMTEdges[i]);
	  ttree->SetBranchAddress(Form("PMTEdgesError[%d]",i),&PMTEdgesError[i]);
	}
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
      TH1F * hTriggerEff_CVLN = new TH1F("hTriggerEff_CVLN","CVLN / CVBN",nRuns,-0.5,nRuns-0.5);
      TH1F * hTriggerEff_CVHN = new TH1F("hTriggerEff_CVHN","CVHN / CVBN",nRuns,-0.5,nRuns-0.5);
      TH1F * hTriggerEff_CVHN2 = new TH1F("hTriggerEff_CVHN2","CVHN / CVLN",nRuns,-0.5,nRuns-0.5);
      TH1F * hPMTEdges[64];
      for(int i = 0; i < 64; ++i)
	hPMTEdges[i] = new  TH1F(Form("hPMTEdges%d",i),Form("Multiplicity edge Cell %d",i),nRuns,-0.5,nRuns-0.5);
	    
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
      list.Add(hTriggerEff_CVLN);
      list.Add(hTriggerEff_CVHN);
      list.Add(hTriggerEff_CVHN2);
      for(int i=0;i<64;i++)
	list.Add(hPMTEdges[i]);
      list.Add(hNumberVoieOff);
      list.Add(hNumberBadOffset);
      char runlabel[6];
      Int_t badrun=0,goodrun=0;
	    
      for(Int_t irun=0;irun<nRuns;irun++)
	{
	  ttree->GetEntry(irun);
	  sprintf(runlabel,"%i",runNumber);

	  if(invalidInput==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);

	      InvalidInput->SetBinContent(badrun+1,invalidInput);
	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);

	    }
	  else if(qaNotFound==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);
		    
	      QaNotFound->SetBinContent(badrun+1,qaNotFound);

	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);

	    }
	  else if(v0active==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);

	      V0active->SetBinContent(badrun+1,v0active);

	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);

	    }
	  else if(v0qaNotfound==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);

	      V0qaNotfound->SetBinContent(badrun+1,v0qaNotfound);

	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);

	    }
	  else if(noEntries==1)
	    {
	      InvalidInput->SetBins(badrun+1,0,badrun+1);	
	      InvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      QaNotFound->SetBins(badrun+1,0,badrun+1);
	      QaNotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0active->SetBins(badrun+1,0,badrun+1);
	      V0active->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      V0qaNotfound->SetBins(badrun+1,0,badrun+1);
	      V0qaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      NoEntrie->SetBins(badrun+1,0,badrun+1);
	      NoEntrie->GetXaxis()->SetBinLabel(badrun+1,runlabel);
	      HigtVoltage->SetBins(badrun,0,badrun);
	      HigtVoltage->GetXaxis()->SetBinLabel(badrun,runlabel);

	      NoEntrie->SetBinContent(badrun+1,noEntries);

	      badrun++;
	      if(higtVoltage==1)
		HigtVoltage->SetBinContent(badrun,higtVoltage);
	    }
	  else 
	    {
	      hTimeA->SetBins(goodrun+1,0,goodrun+1);
	      hTimeA->SetBinContent(goodrun+1,TimesA);
	      hTimeA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hTimeC->SetBins(goodrun+1,0,goodrun+1);
	      hTimeC->SetBinContent(goodrun+1,TimesC);
	      hTimeC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hBB_BG->SetBins(goodrun+1,0,goodrun+1);
	      hBB_BG->SetBinContent(goodrun+1,BB_BG);
	      hBB_BG->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hBB_EE->SetBins(goodrun+1,0,goodrun+1);
	      hBB_EE->SetBinContent(goodrun+1,BB_EE);
	      hBB_EE->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hAdcA->SetBins(goodrun+1,0,goodrun+1);
	      hAdcA->SetBinContent(goodrun+1,AdcA);
	      hAdcA->SetBinError(goodrun+1,AdcAError);
	      hAdcA ->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hAdcC->SetBins(goodrun+1,0,goodrun+1);
	      hAdcC->SetBinContent(goodrun+1,AdcC);
	      hAdcC->SetBinError(goodrun+1,AdcCError);
	      hAdcC ->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hMultA->SetBins(goodrun+1,0,goodrun+1);
	      hMultA->SetBinContent(goodrun+1,MultA);
	      hMultA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hMultC->SetBins(goodrun+1,0,goodrun+1);
	      hMultC->SetBinContent(goodrun+1,MultC);
	      hMultC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hTriggerEff_CVLN->SetBins(goodrun+1,0,goodrun+1);
	      hTriggerEff_CVLN ->SetBinContent(goodrun+1,TriggerEff_CVLN );
	      hTriggerEff_CVLN ->SetBinError(goodrun+1,TriggerEff_CVLN_Error);
	      hTriggerEff_CVLN ->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hTriggerEff_CVHN->SetBins(goodrun+1,0,goodrun+1);
	      hTriggerEff_CVHN ->SetBinContent(goodrun+1,TriggerEff_CVHN );
	      hTriggerEff_CVHN ->SetBinError(goodrun+1,TriggerEff_CVHN_Error);
	      hTriggerEff_CVHN ->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      hTriggerEff_CVHN2->SetBins(goodrun+1,0,goodrun+1);
	      hTriggerEff_CVHN2->SetBinContent(goodrun+1,TriggerEff_CVHN2);
	      hTriggerEff_CVHN2 ->SetBinError(goodrun+1,TriggerEff_CVHN2_Error);
	      hTriggerEff_CVHN2->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		    
	      for(int i=0;i<64;i++)
		{
		  hPMTEdges[i]->SetBins(goodrun+1,0,goodrun+1);
		  hPMTEdges[i]->SetBinContent(goodrun+1,PMTEdges[i]);
		  hPMTEdges[i]->SetBinError(goodrun+1,PMTEdgesError[i]);
		  hPMTEdges[i]->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
		}
		    
		    
	      hNumberVoieOff->SetBins(goodrun+1,0,goodrun+1);
	      hNumberVoieOff->SetBinContent(goodrun+1,NumberVoieOff);
	      hNumberVoieOff->GetXaxis()->SetBinLabel(goodrun+1,runlabel);

	      hNumberBadOffset->SetBins(goodrun+1,0,goodrun+1);
	      hNumberBadOffset->SetBinContent(goodrun+1,numberBadOffset);
	      hNumberBadOffset->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	      maxVoieOff=NumberVoieOff;
	      maxBadOffset=numberBadOffset;

	      goodrun++;
	    }
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
      TPave * pavA = new TPave(0,TMath::Max(hTimeA->GetMinimum(),1.5-shiftA),goodrun,TMath::Min(hTimeA->GetMaximum(),33.5-shiftA),0);
      pavA->SetFillStyle(3004);
      pavA->SetFillColor(2);
      TPave * pavC = new TPave(0,TMath::Max(hTimeC->GetMinimum(),0.5),goodrun,TMath::Min(hTimeC->GetMaximum(),25.5),0);
      pavC->SetFillStyle(3005);
      pavC->SetFillColor(4);
	    
      pavA->Draw("same");
      pavC->Draw("same");
  
	    
      c->Print(Form("%s/QA_Resume_%d_%d.pdf(",plotDir.Data(),minRun,maxRun));
      c->SaveAs(Form("%s/V0QA__Leading_time.png",plotDir.Data()));
      c->Write();
	    
      TCanvas * c2 = new TCanvas("c2","Trigger ratios");
      c2->SetGridy();


      Double_t max=TMath::Max(hBB_BG->GetMaximum(),hBB_EE->GetMaximum());
      Double_t min=TMath::Min(hBB_BG->GetMinimum(0.000001),hBB_EE->GetMinimum(0.0001));
      printf("max: %.5f\tmin: %.5f\n",max,min);
      if(max==0) max=1.0;
      Double_t relatMean=(max-min)/max;
      hBB_BG->SetMarkerStyle(20);
      hBB_BG->SetMarkerColor(2);
      hBB_EE->SetMarkerStyle(20);
      hBB_EE->SetMarkerColor(4);
      hBB_BG->GetYaxis()->SetRangeUser(min*0.8,max*1.5);
      hBB_BG->Draw("P");
      hBB_EE->Draw("Psame");
      TLegend * lg2 = new TLegend(0.8,0.9,1,1);
      lg2->AddEntry(hBB_BG,"BG / BB","p");
      lg2->AddEntry(hBB_EE,"EE / BB","p");
      lg2->Draw("same");

      c2->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      c2->SaveAs(Form("%s/V0QA__Trigger_ratios.png",plotDir.Data()));
      c2->Write();
	    
	    
      TCanvas * c3 = new TCanvas("c3","Average Charge");
      c3->SetGridy();
	    
      hAdcA->SetMarkerStyle(20);
      hAdcA->SetMarkerColor(2);
      hAdcA->SetLineColor(2);
      hAdcC->SetMarkerStyle(20);
      hAdcC->SetMarkerColor(4);
      hAdcC->SetLineColor(4);
      hAdcA->SetMinimum(0);
      max=TMath::Max(hAdcA->GetMaximum(),hAdcC->GetMaximum());
      hAdcA->SetMaximum(max*1.1);

      hAdcA->Draw("P");
      hAdcC->Draw("Psame");
      TLegend * lg3 = new TLegend(0.8,0.9,1,1);
      lg3->AddEntry(hAdcA,"V0A","p");
      lg3->AddEntry(hAdcC,"V0C","p");
      lg3->Draw("same");
	    
      c3->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      c3->SaveAs(Form("%s/V0QA__Average_charge.png",plotDir.Data()));
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
      c4->SaveAs(Form("%s/V0QA__Average_number_of_Cell.png",plotDir.Data()));
      c4->Write();
	    

      TCanvas * c5 = new TCanvas("c5","Trigger Efficiency");
      c5->SetGridy();
	    
      hTriggerEff_CVLN->SetMarkerStyle(20);
      hTriggerEff_CVLN->SetMarkerColor(2);
      hTriggerEff_CVHN->SetMarkerStyle(20);
      hTriggerEff_CVHN->SetMarkerColor(4);
      hTriggerEff_CVHN2->SetMarkerStyle(4);
      hTriggerEff_CVHN2->SetMarkerColor(1);
      hTriggerEff_CVLN->SetMinimum(0);
      hTriggerEff_CVLN->SetMaximum(15.);
      hTriggerEff_CVLN->SetTitle("Centrality triggers fraction");
	    
      hTriggerEff_CVLN->Draw("P");
      hTriggerEff_CVHN->Draw("Psame");
      hTriggerEff_CVHN2->Draw("Psame");
      TLegend * lg5 = new TLegend(0.7,0.8,1,1);
      lg5->AddEntry(hTriggerEff_CVLN,"(CPBI2/100) / CVLN","p");
      lg5->AddEntry(hTriggerEff_CVHN,"(CPBI2/100) / CVHN","p");
      lg5->AddEntry(hTriggerEff_CVHN2,"CVLN / CVHN","p");
      lg5->Draw("same");
	    
      c5->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      c5->SaveAs(Form("%s/V0QA__PbPb_Trigger_Efficiency.png",plotDir.Data()));
      c5->Write();
	    
      TCanvas * cedge[8];
      for(int i = 0; i < 8; ++i)
	{
	  cedge[i] = new TCanvas(Form("cedge%d",i),Form("Edge Ring %d",i));
	  cedge[i]->SetGridy();
	  cedge[i]->Divide(3,3);
	  for(int iCh = 0; iCh < 8; ++iCh)
	    {
	      cedge[i]->cd(iCh+1);
	      hPMTEdges[iCh+i*8]->SetMarkerStyle(20);
	      hPMTEdges[iCh+i*8]->Draw("P");
	    }
	  cedge[i]->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  cedge[i]->SaveAs(Form("%s/V0QA__Edge_Ring_%d.png",plotDir.Data(),i));
	  cedge[i]->Write();
	}
	    
      TCanvas * cT = new TCanvas("cT","");
      cT->cd();
      hNumberVoieOff->GetYaxis()->SetRangeUser(0,maxVoieOff+1);
      hNumberVoieOff->Draw();
      cT->SaveAs(Form("%s/V0QA__Number_of_Voie_Off.png",plotDir.Data()));
      TCanvas * cT2 = new TCanvas("cT2","");
      cT2->cd();
      hNumberBadOffset->GetYaxis()->SetRangeUser(0,maxBadOffset+1.);
      hNumberBadOffset->Draw();
      cT2->SaveAs(Form("%s/V0QA__Number_Bad_Offset.png",plotDir.Data()));
      if(HigtVoltage->GetEntries())
	{
	  TCanvas * c7 = new TCanvas("c7","");
	  c7->cd();
	  HigtVoltage->Draw();
	  c7->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c7->SaveAs(Form("%s/V0QA__Hight_Voltage.png",plotDir.Data()));
	}
      if(InvalidInput->GetEntries())
	{
	  TCanvas * c8 = new TCanvas("c8","");
	  c8->cd();
	  InvalidInput->Draw();
	  c8->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c8->SaveAs(Form("%s/V0QA__InvalidInput_%d_%d.png",plotDir.Data(),minRun,maxRun));
	}
      if(QaNotFound->GetEntries())
	{
	  TCanvas * c9 = new TCanvas("c9","");
	  c9->cd();
	  QaNotFound->Draw();
	  c9->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c9->SaveAs(Form("%s/V0QA__QA_Not_Found_%d_%d.png",plotDir.Data()));
	}
      if(V0active->GetEntries())
	{
	  TCanvas * c10 = new TCanvas("c10","");
	  c10->cd();
	  V0active->Draw();
	  c10->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c10->SaveAs(Form("%s/V0QA__V0_not_Active.png",plotDir.Data()));
	}
      if(V0qaNotfound->GetEntries())
	{
	  TCanvas * c11 = new TCanvas("c11","");
	  c11->cd();
	  V0qaNotfound->Draw();
	  c11->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c11->SaveAs(Form("%s/V0QA__v0qa_Not_Found.png",plotDir.Data()));
	}
      if(NoEntrie->GetEntries())
	{
	  TCanvas * c12 = new TCanvas("c12","");
	  c12->cd();
	  NoEntrie->Draw();
	  c12->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
	  c12->SaveAs(Form("%s/V0QA__NoEntrie.png",plotDir.Data()));
	}
      cT->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cT2->Print(Form("%s/QA_Resume_%d_%d.pdf)",plotDir.Data(),minRun,maxRun));
      cT->Write();
      cT2->Write();
    }
  return 0;
}
