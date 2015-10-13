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

Int_t DrawTrendingADQA(TString mergedTrendFile ="trending.root")
{
  Int_t badrun=0;
  Int_t goodrun=0;
  if(!mergedTrendFile)
    {
      printf("Cannot open merged trend file with AD QA");
      return 1;
    }
  char outfilename[200]="ProductionQA.hist.root";
  TString plotDir(".");
  
  TFile *treandFile = TFile::Open(mergedTrendFile.Data());
  if(!treandFile)
    {
      Printf("ERROR: trending file not found. Exiting ...\n");
      return -1;
    }
  TTree*ttree=(TTree*)treandFile->Get("trending");
  if(!ttree)
    {
      printf("Invalid trending tree");
      return 2;
    }
  //----------------------Take trending tree------------------------------------ 
  Int_t runNumber=0;
  Int_t adReady=0,invalidInput=0,qaNotFound=0,adActive=0,adQANotfound=0,noEntries=0;
  Int_t LHCstate = -1; Float_t runTime = 0;
  Float_t meanTotalChargeADA = -1024, meanTotalChargeADC = -1024;
  Float_t meanTimeADA = -1024, meanTimeADC = -1024;
  Float_t meanTimeSigmaADA = -1024, meanTimeSigmaADC = -1024;
  Float_t meanTimeErrADA = -1024, meanTimeErrADC = -1024;
  Float_t meanTimeSigmaErrADA = -1024, meanTimeSigmaErrADC = -1024;
  Float_t rateUBA = -1024, rateUBC = -1024, rateUGA = -1024, rateUGC = -1024;
  Float_t rateADAND = -1024, rateADOR = -1024, rateErr = -1024;
  Float_t rateRatioADV0AND = -1024, rateRatioADV0OR = -1024;
  Float_t saturationADA = -1024, saturationADC = -1024;
  Float_t MPV[16], MPVErr[16];
  Float_t meanPedestal[32],widthPedestal[32];
  Float_t slewingChi2ADA = -1024, slewingChi2ADC = -1024; 
  Float_t ratePhysADAND = -1024, ratePhysADOR = -1024;
  Float_t ratePhysBBA = -1024, ratePhysBBC = -1024;
  Float_t ratePhysBGA = -1024, ratePhysBGC = -1024;
  Float_t channelTimeMean[16], channelTimeSigma[16];
  Float_t flagNoTimeFraction[16];
  Float_t thresholdData[16], thresholdOCDB[16];
  
  ttree->SetBranchAddress("adReady",&adReady);
  ttree->SetBranchAddress("invalidInput",&invalidInput);
  ttree->SetBranchAddress("qaNotFound",&qaNotFound);
  ttree->SetBranchAddress("adActive",&adActive);
  ttree->SetBranchAddress("adQANotfound",&adQANotfound);
  ttree->SetBranchAddress("noEntries",&noEntries);
  ttree->SetBranchAddress("run",&runNumber);
  
  ttree->SetBranchAddress("LHCstate",&LHCstate);
  ttree->SetBranchAddress("runTime",&runTime);
  ttree->SetBranchAddress("meanTotalChargeADA",&meanTotalChargeADA);
  ttree->SetBranchAddress("meanTotalChargeADC",&meanTotalChargeADC);
  ttree->SetBranchAddress("meanTimeADA",&meanTimeADA);
  ttree->SetBranchAddress("meanTimeADC",&meanTimeADC);
  ttree->SetBranchAddress("meanTimeErrADA",&meanTimeErrADA);
  ttree->SetBranchAddress("meanTimeErrADC",&meanTimeErrADC);
  ttree->SetBranchAddress("meanTimeSigmaADA",&meanTimeSigmaADA);
  ttree->SetBranchAddress("meanTimeSigmaADC",&meanTimeSigmaADC);
  ttree->SetBranchAddress("meanTimeSigmaErrADA",&meanTimeSigmaErrADA);
  ttree->SetBranchAddress("meanTimeSigmaErrADC",&meanTimeSigmaErrADC);
  ttree->SetBranchAddress("rateUBA",&rateUBA);
  ttree->SetBranchAddress("rateUBC",&rateUBC);
  ttree->SetBranchAddress("rateUGA",&rateUGA);
  ttree->SetBranchAddress("rateUGC",&rateUGC);
  ttree->SetBranchAddress("rateADAND",&rateADAND);
  ttree->SetBranchAddress("rateADOR",&rateADOR);
  ttree->SetBranchAddress("rateRatioADV0AND",&rateRatioADV0AND);
  ttree->SetBranchAddress("rateRatioADV0OR",&rateRatioADV0OR);
  ttree->SetBranchAddress("rateErr",&rateErr);
  ttree->SetBranchAddress("MPV",&MPV);
  ttree->SetBranchAddress("MPVErr",&MPVErr);
  ttree->SetBranchAddress("meanPedestal",&meanPedestal);
  ttree->SetBranchAddress("widthPedestal",&widthPedestal);
  ttree->SetBranchAddress("slewingChi2ADA",&slewingChi2ADA);
  ttree->SetBranchAddress("slewingChi2ADC",&slewingChi2ADC);
  ttree->SetBranchAddress("saturationADA",&saturationADA);
  ttree->SetBranchAddress("saturationADC",&saturationADC);
  ttree->SetBranchAddress("ratePhysADAND",&ratePhysADAND);
  ttree->SetBranchAddress("ratePhysADOR",&ratePhysADOR);
  ttree->SetBranchAddress("ratePhysBBA",&ratePhysBBA);
  ttree->SetBranchAddress("ratePhysBBC",&ratePhysBBC);
  ttree->SetBranchAddress("ratePhysBGA",&ratePhysBGA);
  ttree->SetBranchAddress("ratePhysBGC",&ratePhysBGC);
  ttree->SetBranchAddress("channelTimeMean", &channelTimeMean[0]);
  ttree->SetBranchAddress("channelTimeSigma", &channelTimeSigma[0]);
  ttree->SetBranchAddress("flagNoTimeFraction", &flagNoTimeFraction[0]);
  ttree->SetBranchAddress("thresholdData", &thresholdData[0]);
  ttree->SetBranchAddress("thresholdOCDB", &thresholdOCDB[0]);
  
  //----------------------Make trending histos------------------------------------
  ttree->GetEntry(0);
  Int_t nRuns=ttree->GetEntries();
  TList fListHist;
  
  TH1I *hADready= new TH1I("hADready","Run with higt voltage off",0,0,1);
  TH1I *hInvalidInput= new TH1I("hInvalidInput","run with invalid input",0,0,1);
  TH1I *hQANotFound= new TH1I("hQANotFound","QA not found",0,0,1);
  TH1I *hADactive= new TH1I("hADactive","AD not active",0,0,1);
  TH1I *hADqaNotfound= new TH1I("hADqaNotfound","AD QA not found",0,0,1);
  TH1I *hNoEntries= new TH1I("hNoEntries","no entries in qa file",0,0,1);
  
  TH1F *hRunTime = new TH1F("hRunTime","Run duration;;Time (min)",nRuns,-0.5,nRuns-0.5);
  TH1I *hLHCstate = new TH1I("hLHCstate","LHC state;;State",nRuns,-0.5,nRuns-0.5);
  TH1I *hNEvents = new TH1I("hNEvents","Number of events;;Number of events",nRuns,-0.5,nRuns-0.5);
      
  TH1F *hMeanTotalChargeADA = new TH1F("hMeanTotalChargeADA","Mean total charge;;Charge (ADC counts)",nRuns,-0.5,nRuns-0.5);
  TH1F *hMeanTotalChargeADC = new TH1F("hMeanTotalChargeADC","Mean total charge;;Charge (ADC counts)",nRuns,-0.5,nRuns-0.5);
  
  TH1F *hMeanTimeADA = new TH1F("hMeanTimeADA","Mean time;;Time (ns)",nRuns,-0.5,nRuns-0.5);
  TH1F *hMeanTimeADC = new TH1F("hMeanTimeADC","Mean time;;Time (ns)",nRuns,-0.5,nRuns-0.5);
  TH1F *hMeanTimeSigmaADA = new TH1F("hMeanTimeSigmaADA","Mean time Sigma;;#sigma (ns)",nRuns,-0.5,nRuns-0.5);
  TH1F *hMeanTimeSigmaADC = new TH1F("hMeanTimeSigmaADC","Mean time Sigma;;#sigma (ns)",nRuns,-0.5,nRuns-0.5);
  
  TH1F *hRateUBA = new TH1F("hRateUBA","Mean trigger rate UB;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  TH1F *hRateUBC = new TH1F("hRateUBC","Mean trigger rate UB;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  TH1F *hRateUGA = new TH1F("hRateUGA","Mean trigger rate UG;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  TH1F *hRateUGC = new TH1F("hRateUGC","Mean trigger rate UG;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  
  TH1F *hRateADAND = new TH1F("hRateADAND","Mean trigger rate AD;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  TH1F *hRateADOR = new TH1F("hRateADOR","Mean trigger rate AD;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  
  TH1F *hRatioVZEROADAND = new TH1F("hRatioVZEROADAND","Trigger rate ratio VZERO/AD",nRuns,-0.5,nRuns-0.5);
  TH1F *hRatioVZEROADOR = new TH1F("hRatioVZEROADOR","Trigger rate ratio VZERO/AD",nRuns,-0.5,nRuns-0.5);
  
  TH2F *hMPV = new TH2F("hMPV","MIP position per channel;Channel;MPV (ADC counts)",16,-0.5,15.5,nRuns,-0.5,nRuns-0.5);
  TH2F *hMeanPedestal = new TH2F("hMeanPedestal","Pedestal mean per channel;;Pedestal (ADC counts)",32,-0.5,31.5,nRuns,-0.5,nRuns-0.5);
  TH2F *hWidthPedestal = new TH2F("hWidthPedestal","Pedestal width per channel;;Width (ADC counts)",32,-0.5,31.5,nRuns,-0.5,nRuns-0.5);
  
  TH1F *hSlewingChi2ADA = new TH1F("hSlewingChi2ADA","Time slewing Chi2;;Chi2",nRuns,-0.5,nRuns-0.5);
  TH1F *hSlewingChi2ADC = new TH1F("hSlewingChi2ADC","Time slewing Chi2;;Chi2",nRuns,-0.5,nRuns-0.5);
  
  TH1F *hSaturationADA = new TH1F("hSaturationADA","Saturation;;Saturation",nRuns,-0.5,nRuns-0.5);
  TH1F *hSaturationADC = new TH1F("hSaturationADC","Saturation;;Saturation",nRuns,-0.5,nRuns-0.5);
  
  TH2F *hChannelTimeMean = new TH2F("hChannelTimeMean","Mean time per channel;Channel;Time (ns)",16,-0.5,15.5,nRuns,-0.5,nRuns-0.5);
  TH2F *hChannelTimeSigma = new TH2F("hChannelTimeSigma","Time resolution per channel;Channel;Time (ns)",16,-0.5,15.5,nRuns,-0.5,nRuns-0.5);
  
  TH1F *hRatePhysADAND = new TH1F("hRatePhysADAND","Physics selection rate AD;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  TH1F *hRatePhysADOR = new TH1F("hRatePhysADOR","Physics selection rate AD;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  TH1F *hRatePhysBBA = new TH1F("hRatePhysBBA","Physics selection rate BB;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  TH1F *hRatePhysBBC = new TH1F("hRatePhysBBC","Physics selection rate BB;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  TH1F *hRatePhysBGA = new TH1F("hRatePhysBGA","Physics selection rate BG;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  TH1F *hRatePhysBGC = new TH1F("hRatePhysBGC","Physics selection rate BG;;Trigger rate",nRuns,-0.5,nRuns-0.5);
  
  TH2F *hFlagNoTime = new TH2F("hFlagNoTime","Fraction of events with BB/BG flag but not time;Channel;Fraction",16,-0.5,15.5,nRuns,-0.5,nRuns-0.5);
  TH2F *hThresholdData = new TH2F("hThresholdData","Threshold from data fit;Channel;Threshold",16,-0.5,15.5,nRuns,-0.5,nRuns-0.5);
  TH2F *hThresholdOCDB = new TH2F("hThresholdOCDB","Threshold from OCDB;Channel;Threshold",16,-0.5,15.5,nRuns,-0.5,nRuns-0.5);
  
  char ChannelInt[10];
  for(Int_t iPM=0; iPM<16; iPM++){
  	for(Int_t iInt=0; iInt<2; iInt++){
		sprintf(ChannelInt,"Ch %d Int%d",iPM,iInt);
  		hMeanPedestal->GetXaxis()->SetBinLabel(1+iPM+16*iInt,ChannelInt);
  		hWidthPedestal->GetXaxis()->SetBinLabel(1+iPM+16*iInt,ChannelInt);
		}
	}
  
  fListHist.Add(hRunTime);
  fListHist.Add(hLHCstate); 
  fListHist.Add(hNEvents);   
  fListHist.Add(hMeanTotalChargeADA);
  fListHist.Add(hMeanTotalChargeADC);
  fListHist.Add(hMeanTimeADA);
  fListHist.Add(hMeanTimeADC);
  fListHist.Add(hMeanTimeSigmaADA);
  fListHist.Add(hMeanTimeSigmaADC);
  fListHist.Add(hRateUBA);
  fListHist.Add(hRateUBC);
  fListHist.Add(hRateUGA);
  fListHist.Add(hRateUGC);
  fListHist.Add(hRateADAND);
  fListHist.Add(hRateADOR);
  fListHist.Add(hRatioVZEROADAND);
  fListHist.Add(hRatioVZEROADOR);
  fListHist.Add(hMPV);
  fListHist.Add(hMeanPedestal);
  fListHist.Add(hWidthPedestal);
  fListHist.Add(hSlewingChi2ADA);
  fListHist.Add(hSlewingChi2ADC);
  fListHist.Add(hSaturationADA);
  fListHist.Add(hSaturationADC);
  fListHist.Add(hChannelTimeMean);
  fListHist.Add(hChannelTimeSigma);
  fListHist.Add(hRatePhysADAND);
  fListHist.Add(hRatePhysADOR);
  fListHist.Add(hRatePhysBBA);
  fListHist.Add(hRatePhysBBC);
  fListHist.Add(hRatePhysBGA);
  fListHist.Add(hRatePhysBGC);
  fListHist.Add(hFlagNoTime);
  fListHist.Add(hThresholdData);
  fListHist.Add(hThresholdOCDB);
  
  
  //----------------------Loop over runs in tree------------------------------------
  char runlabel[6];      
  for(Int_t irun=0;irun<nRuns;irun++){
  
      ttree->GetEntry(irun);
      sprintf(runlabel,"%i",runNumber);
      
      //----------------------Bad runs------------------------------------
      if(invalidInput==1 || qaNotFound==1 || adActive==1 || adQANotfound==1 || noEntries==1){
      	hInvalidInput->SetBins(badrun+1,0,badrun+1);	  
      	hInvalidInput->GetXaxis()->SetBinLabel(badrun+1,runlabel);
      	hQANotFound->SetBins(badrun+1,0,badrun+1);
      	hQANotFound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
      	hADactive->SetBins(badrun+1,0,badrun+1);
      	hADactive->GetXaxis()->SetBinLabel(badrun+1,runlabel);
      	hADqaNotfound->SetBins(badrun+1,0,badrun+1);
      	hADqaNotfound->GetXaxis()->SetBinLabel(badrun+1,runlabel);
      	hNoEntries->SetBins(badrun+1,0,badrun+1);
      	hNoEntries->GetXaxis()->SetBinLabel(badrun+1,runlabel);
      	hADready->SetBins(badrun,0,badrun);
      	hADready->GetXaxis()->SetBinLabel(badrun,runlabel);
      	if(adReady==1) hADready->SetBinContent(badrun,adReady);
	}
      if(invalidInput==1){
    	hInvalidInput->SetBinContent(badrun+1,invalidInput);
    	badrun++;
    	}
      else if(qaNotFound==1){
	hQANotFound->SetBinContent(badrun+1,qaNotFound);
	badrun++;
	}
      else if(adActive==1){ 
	hADactive->SetBinContent(badrun+1,adActive);
	badrun++;
	}
      else if(adQANotfound==1){
	hADqaNotfound->SetBinContent(badrun+1,adQANotfound);
	badrun++;
	}
      else if(noEntries==1){ 
	hNoEntries->SetBinContent(badrun+1,noEntries);
	badrun++;
	}
      //----------------------Good runs------------------------------------
      else 
	{
	hRunTime->SetBins(goodrun+1,0,goodrun+1);
	hRunTime->SetBinContent(goodrun+1,runTime);
	hRunTime->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hLHCstate->SetBins(goodrun+1,0,goodrun+1);
	hLHCstate->SetBinContent(goodrun+1,LHCstate);
	hLHCstate->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hNEvents->SetBins(goodrun+1,0,goodrun+1);
	hNEvents->SetBinContent(goodrun+1,(Int_t)2*TMath::Power(rateErr,-2));
	hNEvents->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hMeanTotalChargeADA->SetBins(goodrun+1,0,goodrun+1);
	hMeanTotalChargeADA->SetBinContent(goodrun+1,meanTotalChargeADA);
	hMeanTotalChargeADA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hMeanTotalChargeADC->SetBins(goodrun+1,0,goodrun+1);
	hMeanTotalChargeADC->SetBinContent(goodrun+1,meanTotalChargeADC);
	hMeanTotalChargeADC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hMeanTimeADA->SetBins(goodrun+1,0,goodrun+1);
	hMeanTimeADA->SetBinContent(goodrun+1,meanTimeADA);
	hMeanTimeADA->SetBinError(goodrun+1,meanTimeErrADA);
	hMeanTimeADA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hMeanTimeADC->SetBins(goodrun+1,0,goodrun+1);
	hMeanTimeADC->SetBinContent(goodrun+1,meanTimeADC);
	hMeanTimeADC->SetBinError(goodrun+1,meanTimeErrADC);
	hMeanTimeADC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);

	hMeanTimeSigmaADA->SetBins(goodrun+1,0,goodrun+1);
	hMeanTimeSigmaADA->SetBinContent(goodrun+1,meanTimeSigmaADA);
	hMeanTimeSigmaADA->SetBinError(goodrun+1,meanTimeSigmaErrADA);
	hMeanTimeSigmaADA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hMeanTimeSigmaADC->SetBins(goodrun+1,0,goodrun+1);
	hMeanTimeSigmaADC->SetBinContent(goodrun+1,meanTimeSigmaADC);
	hMeanTimeSigmaADC->SetBinError(goodrun+1,meanTimeSigmaErrADC);
	hMeanTimeSigmaADC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);

	hRateUBA->SetBins(goodrun+1,0,goodrun+1);
	hRateUBA->SetBinContent(goodrun+1,rateUBA);
	if(rateUBA!=0)hRateUBA->SetBinError(goodrun+1,rateErr/rateUBA);
	hRateUBA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRateUBC->SetBins(goodrun+1,0,goodrun+1);
	hRateUBC->SetBinContent(goodrun+1,rateUBC);
	if(rateUBC!=0)hRateUBC->SetBinError(goodrun+1,rateErr/rateUBC);
	hRateUBC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRateUGA->SetBins(goodrun+1,0,goodrun+1);
	hRateUGA->SetBinContent(goodrun+1,rateUGA);
	if(rateUGA!=0)hRateUGA->SetBinError(goodrun+1,rateErr/rateUGA);
	hRateUGA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRateUGC->SetBins(goodrun+1,0,goodrun+1);
	hRateUGC->SetBinContent(goodrun+1,rateUGC);
	if(rateUGC!=0)hRateUGC->SetBinError(goodrun+1,rateErr/rateUGC);
	hRateUGC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);

	hRateADAND->SetBins(goodrun+1,0,goodrun+1);
	hRateADAND->SetBinContent(goodrun+1,rateADAND);
	if(rateADAND!=0)hRateADAND->SetBinError(goodrun+1,rateErr/rateADAND);
	hRateADAND->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRateADOR->SetBins(goodrun+1,0,goodrun+1);
	hRateADOR->SetBinContent(goodrun+1,rateADOR);
	if(rateADOR!=0)hRateADOR->SetBinError(goodrun+1,rateErr/rateADOR);
	hRateADOR->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRatioVZEROADAND->SetBins(goodrun+1,0,goodrun+1);
	hRatioVZEROADAND->SetBinContent(goodrun+1,rateRatioADV0AND);
	if(hRatioVZEROADAND!=0)hRatioVZEROADAND->SetBinError(goodrun+1,rateErr/rateRatioADV0AND);
	hRatioVZEROADAND->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRatioVZEROADOR->SetBins(goodrun+1,0,goodrun+1);
	hRatioVZEROADOR->SetBinContent(goodrun+1,rateRatioADV0OR);
	if(hRatioVZEROADOR!=0)hRatioVZEROADOR->SetBinError(goodrun+1,rateErr/rateRatioADV0OR);
	hRatioVZEROADOR->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRatePhysADAND->SetBins(goodrun+1,0,goodrun+1);
	hRatePhysADAND->SetBinContent(goodrun+1,rateADAND);
	if(rateADAND!=0)hRatePhysADAND->SetBinError(goodrun+1,rateErr/ratePhysADAND);
	hRatePhysADAND->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRatePhysADOR->SetBins(goodrun+1,0,goodrun+1);
	hRatePhysADOR->SetBinContent(goodrun+1,ratePhysADOR);
	if(ratePhysADOR!=0)hRatePhysADOR->SetBinError(goodrun+1,rateErr/ratePhysADOR);
	hRatePhysADOR->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRatePhysBBA->SetBins(goodrun+1,0,goodrun+1);
	hRatePhysBBA->SetBinContent(goodrun+1,ratePhysBBA);
	if(ratePhysBBA!=0)hRatePhysBBA->SetBinError(goodrun+1,rateErr/ratePhysBBA);
	hRatePhysBBA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRatePhysBBC->SetBins(goodrun+1,0,goodrun+1);
	hRatePhysBBC->SetBinContent(goodrun+1,ratePhysBBC);
	if(ratePhysBBC!=0)hRatePhysBBC->SetBinError(goodrun+1,rateErr/ratePhysBBC);
	hRatePhysBBC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRatePhysBGA->SetBins(goodrun+1,0,goodrun+1);
	hRatePhysBGA->SetBinContent(goodrun+1,ratePhysBGA);
	if(ratePhysBGA!=0)hRatePhysBGA->SetBinError(goodrun+1,rateErr/ratePhysBGA);
	hRatePhysBGA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hRatePhysBGC->SetBins(goodrun+1,0,goodrun+1);
	hRatePhysBGC->SetBinContent(goodrun+1,ratePhysBGC);
	if(ratePhysBGC!=0)hRatePhysBGC->SetBinError(goodrun+1,rateErr/ratePhysBGC);
	hRatePhysBGC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hMPV->SetBins(16,-0.5,15.5,goodrun+1,0,goodrun+1);
	for(Int_t i=0; i<16; i++){
		hMPV->SetBinContent(i+1,goodrun+1,MPV[i]);
		hMPV->SetBinError(i+1,goodrun+1,MPVErr[i]);
		}
	hMPV->GetYaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hMeanPedestal->SetBins(32,-0.5,31.5,goodrun+1,0,goodrun+1);
	for(Int_t i=0; i<32; i++)hMeanPedestal->SetBinContent(i+1,goodrun+1,meanPedestal[i]);
	hMeanPedestal->GetYaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hWidthPedestal->SetBins(32,-0.5,31.5,goodrun+1,0,goodrun+1);
	for(Int_t i=0; i<32; i++)hWidthPedestal->SetBinContent(i+1,goodrun+1,widthPedestal[i]);
	hWidthPedestal->GetYaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hSlewingChi2ADA->SetBins(goodrun+1,0,goodrun+1);
	hSlewingChi2ADA->SetBinContent(goodrun+1,slewingChi2ADA);
	hSlewingChi2ADA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hSlewingChi2ADC->SetBins(goodrun+1,0,goodrun+1);
	hSlewingChi2ADC->SetBinContent(goodrun+1,slewingChi2ADC);
	hSlewingChi2ADC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hSaturationADA->SetBins(goodrun+1,0,goodrun+1);
	hSaturationADA->SetBinContent(goodrun+1,saturationADA);
	hSaturationADA->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hSaturationADC->SetBins(goodrun+1,0,goodrun+1);
	hSaturationADC->SetBinContent(goodrun+1,saturationADC);
	hSaturationADC->GetXaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hChannelTimeMean->SetBins(16,-0.5,15.5,goodrun+1,0,goodrun+1);
	for(Int_t i=0; i<16; i++) hChannelTimeMean->SetBinContent(i+1,goodrun+1,channelTimeMean[i]);
	hChannelTimeMean->GetYaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hChannelTimeSigma->SetBins(16,-0.5,15.5,goodrun+1,0,goodrun+1);
	for(Int_t i=0; i<16; i++) hChannelTimeSigma->SetBinContent(i+1,goodrun+1,channelTimeSigma[i]);
	hChannelTimeSigma->GetYaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hFlagNoTime->SetBins(16,-0.5,15.5,goodrun+1,0,goodrun+1);
	for(Int_t i=0; i<16; i++) hFlagNoTime->SetBinContent(i+1,goodrun+1,flagNoTimeFraction[i]);
	hFlagNoTime->GetYaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hThresholdData->SetBins(16,-0.5,15.5,goodrun+1,0,goodrun+1);
	for(Int_t i=0; i<16; i++) hThresholdData->SetBinContent(i+1,goodrun+1,thresholdData[i]);
	hThresholdData->GetYaxis()->SetBinLabel(goodrun+1,runlabel);
	
	hThresholdOCDB->SetBins(16,-0.5,15.5,goodrun+1,0,goodrun+1);
	for(Int_t i=0; i<16; i++) hThresholdOCDB->SetBinContent(i+1,goodrun+1,thresholdOCDB[i]);
	hThresholdOCDB->GetYaxis()->SetBinLabel(goodrun+1,runlabel);
	
	goodrun++;
        }
    }
    
  TFile*fout=new TFile(outfilename,"recreate");
  fout->cd();
  fListHist.Write();
  fout->Close();
  
  //----------------------Print trending plots------------------------------------   
  int maxRun =runNumber;
  ttree->GetEntry(0);
  int minRun = runNumber;
    	
  myOptions();
  gROOT->ForceStyle();
  
  TDatime now;
  int iDate = now.GetDate();
  int iYear=iDate/10000;
  int iMonth=(iDate%10000)/100;
  int iDay=iDate%100;
  char* cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
  		    "Jul","Aug","Sep","Oct","Nov","Dec"};
  char cStamp1[25],cStamp2[25];
  sprintf(cStamp1,"%i %s %i",iDay, cMonth[iMonth-1], iYear);
  sprintf(cStamp2,"%i/%.2d/%i",iDay, iMonth, iYear);

  

  TCanvas *c1 = new TCanvas("MeanTotalCharge"," ",800,400); 
  c1->Draw();								 
  c1->cd();
  TPad *myPad1 = new TPad("myPad1", "The pad",0,0,1,1);
  myPadSetUp(myPad1,0.15,0.1,0.04,0.15);
  myPad1->SetGridy();
  myPad1->Draw();
  myPad1->cd();

  myHistSetUp(hMeanTotalChargeADA);
  myHistSetUp(hMeanTotalChargeADC);
  hMeanTotalChargeADA->SetMarkerColor(kBlue);
  hMeanTotalChargeADC->SetMarkerColor(kRed);
  hMeanTotalChargeADA->SetMarkerStyle(kFullCircle);
  hMeanTotalChargeADC->SetMarkerStyle(kFullCircle); 
  myScaleSetUp(hMeanTotalChargeADA,hMeanTotalChargeADC); 
  hMeanTotalChargeADA->Draw("P");
  hMeanTotalChargeADC->Draw("Psame");

  TLegend *myLegend1 = new TLegend(0.70,0.67,0.97,0.82);
  myLegendSetUp(myLegend1,0.04,1);
  myLegend1->AddEntry(hMeanTotalChargeADA,"ADA","p");
  myLegend1->AddEntry(hMeanTotalChargeADC,"ADC","p");
  myLegend1->Draw();
    	
  c1->Print(Form("%s/QA_Resume_%d_%d.pdf(",plotDir.Data(),minRun,maxRun));
  c1->SaveAs(Form("%s/ADQA__Mean_charge.png",plotDir.Data()));
  
  TCanvas *c2 = new TCanvas("MeanTime"," ",800,400); 
  c2->Draw();						
  c2->cd();
  TPad *myPad2 = new TPad("myPad2", "The pad",0,0,1,1);
  myPadSetUp(myPad2,0.15,0.1,0.04,0.15);
  myPad2->SetGridy();
  myPad2->Draw();
  myPad2->cd();

  myHistSetUp(hMeanTimeADA);
  myHistSetUp(hMeanTimeADC);
  hMeanTimeADA->SetMarkerColor(kBlue);
  hMeanTimeADC->SetMarkerColor(kRed);
  hMeanTimeADA->SetMarkerStyle(kFullCircle);
  hMeanTimeADC->SetMarkerStyle(kFullCircle);  
  //myScaleSetUp(hMeanTimeADA,hMeanTimeADC); 
  hMeanTimeADA->GetYaxis()->SetRangeUser(50,70); 
  hMeanTimeADA->Draw("E");
  hMeanTimeADC->Draw("Esame");
  myLegend1->Draw();
  
  c2->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c2->SaveAs(Form("%s/ADQA__Mean_time.png",plotDir.Data()));
  
  TCanvas *c3 = new TCanvas("MeanTimeSigma"," ",800,400); 
  c3->Draw();						
  c3->cd();
  TPad *myPad3 = new TPad("myPad3", "The pad",0,0,1,1);
  myPadSetUp(myPad3,0.15,0.1,0.04,0.15);
  myPad3->SetGridy();
  myPad3->Draw();
  myPad3->cd();

  myHistSetUp(hMeanTimeSigmaADA);
  myHistSetUp(hMeanTimeSigmaADC);
  hMeanTimeSigmaADA->SetMarkerColor(kBlue);
  hMeanTimeSigmaADC->SetMarkerColor(kRed);
  hMeanTimeSigmaADA->SetMarkerStyle(kFullCircle);
  hMeanTimeSigmaADC->SetMarkerStyle(kFullCircle);
  //myScaleSetUp(hMeanTimeSigmaADA,hMeanTimeSigmaADC);
  hMeanTimeSigmaADA->GetYaxis()->SetRangeUser(0.1,1.0);
  hMeanTimeSigmaADA->Draw("E");
  hMeanTimeSigmaADC->Draw("Esame");
  myLegend1->Draw();
  
  c3->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c3->SaveAs(Form("%s/ADQA__Mean_time_Sigma.png",plotDir.Data()));
  
  TCanvas *c21 = new TCanvas("ChannelTime"," ",800,400); 
  c21->Draw();						
  c21->cd();
  TPad *myPad21 = new TPad("myPad21", "The pad",0,0,1,1);
  myPadSetUp(myPad21,0.15,0.1,0.04,0.15);
  myPad21->SetGridy();
  myPad21->Draw();
  myPad21->cd();

  myHistSetUp(hChannelTimeMean); 
  hChannelTimeMean->Draw("COLZTEXT");
  
  c21->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  TCanvas *c22 = new TCanvas("Channel time 1D"," ",1200,800); 
  c22->Draw();						
  c22->cd();
  TPad *myPad22 = new TPad("myPad22", "The pad",0,0,1,1);
  myPad22->Divide(4,4);
  myPad22->Draw();
  
  TH1D *hChannelSlice;
  for(Int_t i = 0; i<16; i++){
  	myPadSetUp(myPad22->cd(i+1),0.15,0.00,0.05,0.15);
  	myPad22->cd(i+1);
  	gPad->SetGridy();
	hChannelSlice = hChannelTimeMean->ProjectionY("hChannelSlice",i+1,i+1);
	myHistSetUp(hChannelSlice);
	hChannelSlice->SetMarkerStyle(kFullCircle);
	hChannelSlice->SetMarkerSize(0.5);
	hChannelSlice->GetXaxis()->SetTitle("");
	hChannelSlice->GetYaxis()->SetTitle("Time mean (ns)");
	hChannelSlice->SetTitle(Form("Mean time channel %d",i));
	if(i<8)hChannelSlice->GetYaxis()->SetRangeUser(62,68);
	else hChannelSlice->GetYaxis()->SetRangeUser(53,59);
	hChannelSlice->DrawCopy("P");
	}
 
  c22->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  TCanvas *c31 = new TCanvas("ChannelSigma"," ",800,400); 
  c31->Draw();						
  c31->cd();
  TPad *myPad31 = new TPad("myPad31", "The pad",0,0,1,1);
  myPadSetUp(myPad31,0.15,0.1,0.04,0.15);
  myPad31->SetGridy();
  myPad31->Draw();
  myPad31->cd();

  myHistSetUp(hChannelTimeSigma); 
  hChannelTimeSigma->Draw("COLZTEXT");
  
  c31->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  TCanvas *c32 = new TCanvas("Channel sigma 1D"," ",1200,800); 
  c32->Draw();						
  c32->cd();
  TPad *myPad32 = new TPad("myPad32", "The pad",0,0,1,1);
  myPad32->Divide(4,4);
  myPad32->Draw();
  
  for(Int_t i = 0; i<16; i++){
  	myPadSetUp(myPad32->cd(i+1),0.15,0.00,0.05,0.15);
  	myPad32->cd(i+1);
  	gPad->SetGridy();
	hChannelSlice = hChannelTimeSigma->ProjectionY("hChannelSlice",i+1,i+1);
	myHistSetUp(hChannelSlice);
	hChannelSlice->SetMarkerStyle(kFullCircle);
	hChannelSlice->SetMarkerSize(0.5);
	hChannelSlice->GetXaxis()->SetTitle("");
	hChannelSlice->GetYaxis()->SetTitle("Time resolution (ns)");
	hChannelSlice->SetTitle(Form("Time resolution channel %d",i));
	hChannelSlice->GetYaxis()->SetRangeUser(0.1,1.1);
	hChannelSlice->DrawCopy("P");
	}
  c32->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  TCanvas *c4 = new TCanvas("Rate UB"," ",800,400); 
  c4->Draw();						
  c4->cd();
  TPad *myPad4 = new TPad("myPad4", "The pad",0,0,1,1);
  myPadSetUp(myPad4,0.15,0.1,0.04,0.15);
  myPad4->SetGridy();
  myPad4->Draw();
  myPad4->cd();

  myHistSetUp(hRateUBA);
  myHistSetUp(hRateUBC);
  hRateUBA->SetMarkerColor(kBlue);
  hRateUBC->SetMarkerColor(kRed);
  hRateUBA->SetMarkerStyle(kFullCircle);
  hRateUBC->SetMarkerStyle(kFullCircle);
  hRateUBA->GetYaxis()->SetRangeUser(0,1.4);
  hRateUBA->Draw("E");
  hRateUBC->Draw("Esame");
  myLegend1->Draw();
  
  c4->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c4->SaveAs(Form("%s/ADQA__RateUB.png",plotDir.Data()));
  
  TCanvas *c41 = new TCanvas("Physics selection rate BB"," ",800,400); 
  c41->Draw();						
  c41->cd();
  TPad *myPad41 = new TPad("myPad41", "The pad",0,0,1,1);
  myPadSetUp(myPad41,0.15,0.1,0.04,0.15);
  myPad41->SetGridy();
  myPad41->Draw();
  myPad41->cd();

  myHistSetUp(hRatePhysBBA);
  myHistSetUp(hRatePhysBBC);
  hRatePhysBBA->SetMarkerColor(kBlue);
  hRatePhysBBC->SetMarkerColor(kRed);
  hRatePhysBBA->SetMarkerStyle(kFullCircle);
  hRatePhysBBC->SetMarkerStyle(kFullCircle);
  hRatePhysBBA->GetYaxis()->SetRangeUser(0,1.4);
  hRatePhysBBA->Draw("E");
  hRatePhysBBC->Draw("Esame");
  myLegend1->Draw();
  
  c41->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c41->SaveAs(Form("%s/ADQA__RatePhysBB.png",plotDir.Data()));
  
  TCanvas *c5 = new TCanvas("Rate UG"," ",800,400); 
  c5->Draw();						
  c5->cd();
  TPad *myPad5 = new TPad("myPad5", "The pad",0,0,1,1);
  myPadSetUp(myPad5,0.15,0.1,0.04,0.15);
  myPad5->SetGridy();
  myPad5->Draw();
  myPad5->cd();

  myHistSetUp(hRateUGA);
  myHistSetUp(hRateUGC);
  hRateUGA->SetMarkerColor(kBlue);
  hRateUGC->SetMarkerColor(kRed);
  hRateUGA->SetMarkerStyle(kFullCircle);
  hRateUGC->SetMarkerStyle(kFullCircle);
  myScaleSetUp(hRateUGA,hRateUGC);
  hRateUGA->Draw("E");
  hRateUGC->Draw("Esame");
  myLegend1->Draw();
  
  c5->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c5->SaveAs(Form("%s/ADQA__RateUG.png",plotDir.Data()));
  
  TCanvas *c51 = new TCanvas("Physics selection rate BG"," ",800,400); 
  c51->Draw();						
  c51->cd();
  TPad *myPad51 = new TPad("myPad51", "The pad",0,0,1,1);
  myPadSetUp(myPad51,0.15,0.1,0.04,0.15);
  myPad51->SetGridy();
  myPad51->Draw();
  myPad51->cd();

  myHistSetUp(hRatePhysBGA);
  myHistSetUp(hRatePhysBGC);
  hRatePhysBGA->SetMarkerColor(kBlue);
  hRatePhysBGC->SetMarkerColor(kRed);
  hRatePhysBGA->SetMarkerStyle(kFullCircle);
  hRatePhysBGC->SetMarkerStyle(kFullCircle);
  myScaleSetUp(hRatePhysBGA,hRatePhysBGC);
  hRatePhysBGA->Draw("E");
  hRatePhysBGC->Draw("Esame");
  myLegend1->Draw();
  
  c51->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c51->SaveAs(Form("%s/ADQA__RatePhysBG.png",plotDir.Data()));
  
  TCanvas *c6 = new TCanvas("Rate AD"," ",800,400); 
  c6->Draw();						
  c6->cd();
  TPad *myPad6 = new TPad("myPad6", "The pad",0,0,1,1);
  myPadSetUp(myPad6,0.15,0.1,0.04,0.15);
  myPad6->SetGridy();
  myPad6->Draw();
  myPad6->cd();

  myHistSetUp(hRateADAND);
  myHistSetUp(hRateADOR);
  hRateADAND->SetMarkerColor(kBlue);
  hRateADOR->SetMarkerColor(kRed);
  hRateADAND->SetMarkerStyle(kFullCircle);
  hRateADOR->SetMarkerStyle(kFullCircle);
  hRateADAND->GetYaxis()->SetRangeUser(0,1.4);
  hRateADAND->Draw("E");
  hRateADOR->Draw("Esame");
  
  TLegend *myLegend2 = new TLegend(0.70,0.67,0.97,0.82);
  myLegendSetUp(myLegend2,0.04,1);
  myLegend2->AddEntry(hRateADAND,"ADAND","p");
  myLegend2->AddEntry(hRateADOR,"ADOR","p");
  myLegend2->Draw();
  
  c6->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c6->SaveAs(Form("%s/ADQA__RateAD.png",plotDir.Data()));
  
  TCanvas *c61 = new TCanvas("Physics selection rate AD"," ",800,400); 
  c61->Draw();						
  c61->cd();
  TPad *myPad61 = new TPad("myPad61", "The pad",0,0,1,1);
  myPadSetUp(myPad61,0.15,0.1,0.04,0.15);
  myPad61->SetGridy();
  myPad61->Draw();
  myPad61->cd();

  myHistSetUp(hRatePhysADAND);
  myHistSetUp(hRatePhysADOR);
  hRatePhysADAND->SetMarkerColor(kBlue);
  hRatePhysADOR->SetMarkerColor(kRed);
  hRatePhysADAND->SetMarkerStyle(kFullCircle);
  hRatePhysADOR->SetMarkerStyle(kFullCircle);
  hRatePhysADAND->GetYaxis()->SetRangeUser(0,1.4);
  hRatePhysADAND->Draw("E");
  hRatePhysADOR->Draw("Esame");

  myLegend2->Draw();
  
  c61->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c61->SaveAs(Form("%s/ADQA__RatePhysAD.png",plotDir.Data()));
  
  TCanvas *c11 = new TCanvas("Rate VZERO_AD"," ",800,400); 
  c11->Draw();						
  c11->cd();
  TPad *myPad11 = new TPad("myPad11", "The pad",0,0,1,1);
  myPadSetUp(myPad11,0.15,0.1,0.04,0.15);
  myPad11->SetGridy();
  myPad11->Draw();
  myPad11->cd();

  myHistSetUp(hRatioVZEROADAND);
  myHistSetUp(hRatioVZEROADOR);
  hRatioVZEROADAND->SetMarkerColor(kBlue);
  hRatioVZEROADOR->SetMarkerColor(kRed);
  hRatioVZEROADAND->SetMarkerStyle(kFullCircle);
  hRatioVZEROADOR->SetMarkerStyle(kFullCircle);
  myScaleSetUp(hRatioVZEROADAND,hRatioVZEROADOR);
  hRatioVZEROADAND->Draw("E");
  hRatioVZEROADOR->Draw("Esame");
  
  TLegend *myLegend3 = new TLegend(0.70,0.67,0.97,0.82);
  myLegendSetUp(myLegend3,0.04,1);
  myLegend3->AddEntry(hRatioVZEROADAND,"VZERO_AND/AD_AND","p");
  myLegend3->AddEntry(hRatioVZEROADOR,"VZERO_OR/AD_OR","p");
  myLegend3->Draw();
  
  c11->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c11->SaveAs(Form("%s/ADQA__RateRatio.png",plotDir.Data()));
  
  TCanvas *c15 = new TCanvas("FlagNoTime"," ",800,400); 
  c15->Draw();						
  c15->cd();
  TPad *myPad15 = new TPad("myPad15", "The pad",0,0,1,1);
  myPadSetUp(myPad15,0.15,0.1,0.04,0.15);
  myPad15->SetGridy();
  myPad15->Draw();
  myPad15->cd();

  myHistSetUp(hFlagNoTime);
  hFlagNoTime->Draw("COLZTEXT");
 
  c15->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c15->SaveAs(Form("%s/ADQA__FlagNoTime.png",plotDir.Data()));
  
  TCanvas *c151 = new TCanvas("FlagNoTime 1D"," ",1200,800); 
  c151->Draw();						
  c151->cd();
  TPad *myPad151 = new TPad("myPad151", "The pad",0,0,1,1);
  myPad151->Divide(4,4);
  myPad151->Draw();
  
  for(Int_t i = 0; i<16; i++){
  	myPadSetUp(myPad151->cd(i+1),0.15,0.00,0.05,0.15);
  	myPad151->cd(i+1);
  	gPad->SetGridy();
	hChannelSlice = hFlagNoTime->ProjectionY("hChannelSlice",i+1,i+1);
	myHistSetUp(hChannelSlice);
	hChannelSlice->SetLineWidth(1);
	myScaleSetUp(hChannelSlice);
	hChannelSlice->GetXaxis()->SetTitle("");
	hChannelSlice->GetYaxis()->SetTitle("Fraction with flag but no time");
	hChannelSlice->SetTitle(Form("Flag but no time fraction channel %d",i));
	//hChannelSlice->GetYaxis()->SetRangeUser(2,30);
	hChannelSlice->DrawCopy("HIST");
	}
 
  c151->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  TCanvas *c7 = new TCanvas("MPV"," ",800,400); 
  c7->Draw();						
  c7->cd();
  TPad *myPad7 = new TPad("myPad7", "The pad",0,0,1,1);
  myPadSetUp(myPad7,0.15,0.1,0.04,0.15);
  myPad7->SetGridy();
  myPad7->Draw();
  myPad7->cd();

  myHistSetUp(hMPV);
  hMPV->Draw("COLZTEXT");
 
  c7->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c7->SaveAs(Form("%s/ADQA__MPV.png",plotDir.Data()));
  
  TCanvas *c71 = new TCanvas("MPV 1D"," ",1200,800); 
  c71->Draw();						
  c71->cd();
  TPad *myPad71 = new TPad("myPad71", "The pad",0,0,1,1);
  myPad71->Divide(4,4);
  myPad71->Draw();
  
  for(Int_t i = 0; i<16; i++){
  	myPadSetUp(myPad71->cd(i+1),0.15,0.00,0.05,0.15);
  	myPad71->cd(i+1);
  	gPad->SetGridy();
	hChannelSlice = hMPV->ProjectionY("hChannelSlice",i+1,i+1);
	myHistSetUp(hChannelSlice);
	hChannelSlice->SetLineWidth(1);
	myScaleSetUp(hChannelSlice);
	hChannelSlice->GetXaxis()->SetTitle("");
	hChannelSlice->GetYaxis()->SetTitle("MPV (ADC counts)");
	hChannelSlice->SetTitle(Form("MPV channel %d",i));
	hChannelSlice->GetYaxis()->SetRangeUser(2,30);
	hChannelSlice->DrawCopy("HIST");
	}
 
  c71->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  
  
  TCanvas *c13 = new TCanvas("SlewingChi2"," ",800,400); 
  c13->Draw();						
  c13->cd();
  TPad *myPad13 = new TPad("myPad13", "The pad",0,0,1,1);
  myPadSetUp(myPad13,0.15,0.1,0.04,0.15);
  myPad13->SetGridy();
  myPad13->Draw();
  myPad13->cd();

  myHistSetUp(hSlewingChi2ADA);
  myHistSetUp(hSlewingChi2ADC);
  hSlewingChi2ADA->SetMarkerColor(kBlue);
  hSlewingChi2ADC->SetMarkerColor(kRed);
  hSlewingChi2ADA->SetMarkerStyle(kFullCircle);
  hSlewingChi2ADC->SetMarkerStyle(kFullCircle);  
  //myScaleSetUp(hSlewingChi2ADA,hSlewingChi2ADC); 
  hSlewingChi2ADA->GetYaxis()->SetRangeUser(-0.1,1); 
  hSlewingChi2ADA->Draw("P");
  hSlewingChi2ADC->Draw("Psame");
  myLegend1->Draw();
  
  c13->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c13->SaveAs(Form("%s/ADQA__SlewingChi2.png",plotDir.Data()));
  
  TCanvas *c14 = new TCanvas("Saturation"," ",800,400); 
  c14->Draw();						
  c14->cd();
  TPad *myPad14 = new TPad("myPad14", "The pad",0,0,1,1);
  myPadSetUp(myPad14,0.15,0.1,0.04,0.15);
  myPad14->SetGridy();
  myPad14->Draw();
  myPad14->cd();

  myHistSetUp(hSaturationADA);
  myHistSetUp(hSaturationADC);
  hSaturationADA->SetMarkerColor(kBlue);
  hSaturationADC->SetMarkerColor(kRed);
  hSaturationADA->SetMarkerStyle(kFullCircle);
  hSaturationADC->SetMarkerStyle(kFullCircle);  
  myScaleSetUp(hSaturationADA,hSaturationADC); 
  hSaturationADA->GetYaxis()->SetRange(0,hSaturationADA->GetYaxis()->GetLast()); 
  hSaturationADA->Draw("P");
  hSaturationADC->Draw("Psame");
  myLegend1->Draw();
  
  c14->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c14->SaveAs(Form("%s/ADQA__Saturation.png",plotDir.Data()));
  
  if(hADready->GetEntries())
    {
      TCanvas * cBad1 = new TCanvas("cBad1","");
      cBad1->cd();
      hADready->Draw();
      cBad1->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cBad1->SaveAs(Form("%s/ADQA__adReady.png",plotDir.Data()));
    }
  if(hInvalidInput->GetEntries())
    {
      TCanvas * cBad2 = new TCanvas("cBad2","");
      cBad2->cd();
      hInvalidInput->Draw();
      cBad2->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cBad2->SaveAs(Form("%s/ADQA__InvalidInput_%d_%d.png",plotDir.Data(),minRun,maxRun));
    }
  if(hQANotFound->GetEntries())
    {
      TCanvas * cBad3 = new TCanvas("cBad3","");
      cBad3->cd();
      hQANotFound->Draw();
      cBad3->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cBad3->SaveAs(Form("%s/ADQA__QA_Not_Found_%d_%d.png",plotDir.Data()));
    }
  if(hADactive->GetEntries())
    {
      TCanvas * cBad4 = new TCanvas("cBad4","");
      cBad4->cd();
      hADactive->Draw();
      cBad4->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cBad4->SaveAs(Form("%s/ADQA__AD_not_Active.png",plotDir.Data()));
    }
  if(hADqaNotfound->GetEntries())
    {
      TCanvas * cBad5 = new TCanvas("cBad5","");
      cBad5->cd();
      hADqaNotfound->Draw();
      cBad5->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cBad5->SaveAs(Form("%s/ADQA__adQA_Not_Found.png",plotDir.Data()));
    }
  if(hNoEntries->GetEntries())
    {
      TCanvas * cBad6 = new TCanvas("cBad6","");
      cBad6->cd();
      hNoEntries->Draw();
      cBad6->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
      cBad6->SaveAs(Form("%s/ADQA__NoEntries.png",plotDir.Data()));
    }
    
  TCanvas *c16 = new TCanvas("Thresholds"," ",1200,800); 
  c16->Draw();						
  c16->cd();
  TPad *myPad16 = new TPad("myPad16", "The pad",0,0,1,1);
  myPad16->Divide(4,4);
  myPad16->Draw();
  
  TH1D *hChannelSliceBlue;
  TH1D *hChannelSliceRed;
  for(Int_t i = 0; i<16; i++){
  	myPadSetUp(myPad16->cd(i+1),0.15,0.00,0.05,0.15);
  	myPad16->cd(i+1);
  	gPad->SetGridy();
	hChannelSliceRed = hThresholdData->ProjectionY("hChannelSliceRed",i+1,i+1);
	myHistSetUp(hChannelSliceRed);
	hChannelSliceRed->SetLineWidth(1);
	hChannelSliceRed->SetLineColor(kRed);
	hChannelSliceRed->GetXaxis()->SetTitle("");
	hChannelSliceRed->GetYaxis()->SetTitle("Threshold [ADC counts]");
	hChannelSliceRed->SetTitle(Form("Threshold, channel %d",i));
	
	hChannelSliceBlue = hThresholdOCDB->ProjectionY("hChannelSliceBlue",i+1,i+1);
	myHistSetUp(hChannelSliceBlue);
	hChannelSliceBlue->SetLineWidth(1);
	hChannelSliceBlue->SetLineColor(kBlue);
	myScaleSetUp(hChannelSliceRed,hChannelSliceBlue);
	hChannelSliceRed->DrawCopy("HIST");
	hChannelSliceBlue->DrawCopy("HISTSAME");
	
	TLegend *myLegend4 = new TLegend(0.36,0.62,0.62,0.84);
        myLegendSetUp(myLegend4,0.08,1);
        myLegend4->AddEntry(hChannelSliceRed,"Data","l");
        myLegend4->AddEntry(hChannelSliceBlue,"OCDB","l");
	myLegend4->Draw();
	}
 
  c16->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));

  TCanvas *c8 = new TCanvas("Pedestal mean"," ",1200,400); 
  c8->Draw();						
  c8->cd();
  TPad *myPad8 = new TPad("myPad8", "The pad",0,0,1,1);
  myPadSetUp(myPad8,0.15,0.1,0.04,0.15);
  myPad8->SetGridy();
  myPad8->Draw();
  myPad8->cd();

  myHistSetUp(hMeanPedestal);
  hMeanPedestal->Draw("COLZTEXT");
 
  c8->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c8->SaveAs(Form("%s/ADQA__PedMean.png",plotDir.Data())); 
  
  TCanvas *c81 = new TCanvas("Pedestal mean 1D Int0"," ",1200,800); 
  c81->Draw();						
  c81->cd();
  TPad *myPad81 = new TPad("myPad81", "The pad",0,0,1,1);
  myPad81->Divide(4,4);
  myPad81->Draw();
  
  for(Int_t i = 0; i<16; i++){
  	myPadSetUp(myPad81->cd(i+1),0.15,0.00,0.05,0.15);
  	myPad81->cd(i+1);
  	gPad->SetGridy();
	hChannelSlice = hMeanPedestal->ProjectionY("hChannelSlice",i+1,i+1);
	myHistSetUp(hChannelSlice);
	hChannelSlice->SetLineWidth(1);
	hChannelSlice->GetXaxis()->SetTitle("");
	hChannelSlice->GetYaxis()->SetTitle("Pedestal mean");
	myScaleSetUp(hChannelSlice);
	hChannelSlice->SetTitle(Form("Pedestal mean channel %d, Int0",i));
	hChannelSlice->DrawCopy("HIST");
	}
 
  c81->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  TCanvas *c82 = new TCanvas("Pedestal mean 1D Int1"," ",1200,800); 
  c82->Draw();						
  c82->cd();
  TPad *myPad82 = new TPad("myPad82", "The pad",0,0,1,1);
  myPad82->Divide(4,4);
  myPad82->Draw();

  for(Int_t i = 0; i<16; i++){
  	myPadSetUp(myPad82->cd(i+1),0.15,0.00,0.05,0.15);
  	myPad82->cd(i+1);
  	gPad->SetGridy();
	hChannelSlice = hMeanPedestal->ProjectionY("hChannelSlice",i+17,i+17);
	myHistSetUp(hChannelSlice);
	hChannelSlice->SetLineWidth(1);
	hChannelSlice->GetXaxis()->SetTitle("");
	hChannelSlice->GetYaxis()->SetTitle("Pedestal mean");
	myScaleSetUp(hChannelSlice);
	hChannelSlice->SetTitle(Form("Pedestal mean channel %d, Int1",i));
	hChannelSlice->DrawCopy("HIST");
	}
 
  c82->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  TCanvas *c9 = new TCanvas("Pedestal width"," ",1200,400); 
  c9->Draw();						
  c9->cd();
  TPad *myPad9 = new TPad("myPad9", "The pad",0,0,1,1);
  myPadSetUp(myPad9,0.15,0.1,0.04,0.15);
  myPad9->SetGridy();
  myPad9->Draw();
  myPad9->cd();

  myHistSetUp(hWidthPedestal);
  hWidthPedestal->Draw("COLZTEXT");
 
  c9->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c9->SaveAs(Form("%s/ADQA__PedWidth.png",plotDir.Data())); 
  
  TCanvas *c91 = new TCanvas("Pedestal width 1D Int0"," ",1200,800); 
  c91->Draw();						
  c91->cd();
  TPad *myPad91 = new TPad("myPad91", "The pad",0,0,1,1);
  myPad91->Divide(4,4);
  myPad91->Draw();
  
  for(Int_t i = 0; i<16; i++){
  	myPadSetUp(myPad91->cd(i+1),0.15,0.00,0.05,0.15);
  	myPad91->cd(i+1);
  	gPad->SetGridy();
	hChannelSlice = hWidthPedestal->ProjectionY("hChannelSlice",i+1,i+1);
	myHistSetUp(hChannelSlice);
	hChannelSlice->SetLineWidth(1);
	hChannelSlice->GetXaxis()->SetTitle("");
	hChannelSlice->GetYaxis()->SetTitle("Pedestal width");
	myScaleSetUp(hChannelSlice);
	hChannelSlice->SetTitle(Form("Pedestal width channel %d, Int0",i));
	hChannelSlice->DrawCopy("HIST");
	}
 
  c91->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  TCanvas *c92 = new TCanvas("Pedestal width 1D Int1"," ",1200,800); 
  c92->Draw();						
  c92->cd();
  TPad *myPad92 = new TPad("myPad92", "The pad",0,0,1,1);
  myPad92->Divide(4,4);
  myPad92->Draw();

  for(Int_t i = 0; i<16; i++){
  	myPadSetUp(myPad92->cd(i+1),0.15,0.00,0.05,0.15);
  	myPad92->cd(i+1);
  	gPad->SetGridy();
	hChannelSlice = hWidthPedestal->ProjectionY("hChannelSlice",i+17,i+17);
	myHistSetUp(hChannelSlice);
	hChannelSlice->SetLineWidth(1);
	hChannelSlice->GetXaxis()->SetTitle("");
	hChannelSlice->GetYaxis()->SetTitle("Pedestal width");
	myScaleSetUp(hChannelSlice);
	hChannelSlice->SetTitle(Form("Pedestal width channel %d, Int1",i));
	hChannelSlice->DrawCopy("HIST");
	}
 
  c92->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  
  TCanvas *c10 = new TCanvas("Run time"," ",800,400); 
  c10->Draw();						
  c10->cd();
  TPad *myPad10 = new TPad("myPad10", "The pad",0,0,1,1);
  myPadSetUp(myPad10,0.15,0.1,0.04,0.15);
  myPad10->SetGridy();
  myPad10->Draw();
  myPad10->cd();

  myHistSetUp(hRunTime);
  hRunTime->SetMarkerStyle(kFullCircle);
  hRunTime->Draw("P");
 
  c10->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c10->SaveAs(Form("%s/ADQA__RunTime.png",plotDir.Data()));
  
  TCanvas *c11 = new TCanvas("LHC state"," ",800,400); 
  c11->Draw();						
  c11->cd();
  TPad *myPad11 = new TPad("myPad11", "The pad",0,0,1,1);
  myPadSetUp(myPad11,0.15,0.1,0.04,0.15);
  myPad11->SetGridy();
  myPad11->Draw();
  myPad11->cd();

  myHistSetUp(hLHCstate);
  hLHCstate->SetMarkerStyle(kFullCircle);
  hLHCstate->GetYaxis()->Set(5,-1.5,3.5);
  hLHCstate->GetYaxis()->SetRangeUser(-0.5,3.5);
  
  hLHCstate->GetYaxis()->SetBinLabel(2,"No Beam");
  hLHCstate->GetYaxis()->SetBinLabel(3,"Squeeze");
  hLHCstate->GetYaxis()->SetBinLabel(4,"Adjust");
  hLHCstate->GetYaxis()->SetBinLabel(5,"Stable beams");
  
  hLHCstate->Draw("P");
 
  c11->Print(Form("%s/QA_Resume_%d_%d.pdf",plotDir.Data(),minRun,maxRun));
  c11->SaveAs(Form("%s/ADQA__LHCstate.png",plotDir.Data())); 
  
  TCanvas *c12 = new TCanvas("Number of events"," ",800,400); 
  c12->Draw();						
  c12->cd();
  TPad *myPad12 = new TPad("myPad12", "The pad",0,0,1,1);
  myPadSetUp(myPad12,0.15,0.1,0.04,0.15);
  myPad12->SetGridy();
  myPad12->Draw();
  myPad12->cd();
  gPad->SetLogy();

  myHistSetUp(hNEvents);
  hNEvents->SetMarkerStyle(kFullCircle);
  hNEvents->Draw("P");
 
  c12->Print(Form("%s/QA_Resume_%d_%d.pdf)",plotDir.Data(),minRun,maxRun));
  c12->SaveAs(Form("%s/ADQA__NEvents.png",plotDir.Data()));
  
  
  
  
  return 0;
}

void myScaleSetUp(TH1* histoBlue, TH1* histoRed){

  Float_t min[2], max[2];

  max[0] = histoBlue->GetBinContent(histoBlue->GetMaximumBin());
  max[1] = histoRed->GetBinContent(histoRed->GetMaximumBin());

  min[0] = histoBlue->GetBinContent(histoBlue->GetMinimumBin());
  min[1] = histoRed->GetBinContent(histoRed->GetMinimumBin());
	
  histoBlue->GetYaxis()->SetRangeUser(0.5*TMath::MinElement(2,min),1.5*TMath::MaxElement(2,max));  

}

void myScaleSetUp(TH1* histo){

  Float_t min, max;

  max = histo->GetBinContent(histo->GetMaximumBin());

  min = histo->GetBinContent(histo->GetMinimumBin());
	
  histo->GetYaxis()->SetRangeUser(0.5*min,1.5*max);  

}


void myPadSetUp(TVirtualPad *currentPad, float currentLeft=0.31, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.11, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07,int columns=2){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  currentLegend->SetNColumns(columns);
  return;
}


void myHistSetUp(TH1 *currentGraph=0){
 
  currentGraph->SetLabelSize(0.05,"xyz");
  currentGraph->SetLabelFont(42,"xyz"); 
  currentGraph->SetLabelOffset(0.01,"xyz");
  currentGraph->SetTitleFont(42,"xyz");   
  currentGraph->GetXaxis()->SetTitleOffset(1.1);
  currentGraph->GetXaxis()->LabelsOption("v");
  currentGraph->GetYaxis()->SetTitleOffset(1.1);  
  currentGraph->SetTitleSize(0.06,"xyz");
  currentGraph->SetStats(kFALSE); 
  return;
}
void myHistSetUp(TH2 *currentGraph=0){
 
  currentGraph->SetLabelSize(0.03,"xyz");
  currentGraph->SetLabelFont(42,"xyz"); 
  currentGraph->SetLabelOffset(0.01,"xyz");
  currentGraph->SetTitleFont(42,"xyz"); 
  currentGraph->GetXaxis()->SetTitleOffset(1.3);
  currentGraph->GetYaxis()->SetTitleOffset(1.3);
  currentGraph->GetYaxis()->LabelsOption("h");  
  currentGraph->SetTitleSize(0.05,"xyz");
  currentGraph->SetStats(kFALSE);
  return;
}

void myOptions(Int_t lStat=0){
  // Set gStyle
  int font = 42;
  // From plain
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.0,"xyz");  
  gStyle->SetTitleSize(0.06,"xyz");  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  gStyle->SetErrorX(0);
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111111);
    gStyle->SetOptFit(1111111);
    }
  else {
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}
