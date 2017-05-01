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
#include "AliADCalibData.h"
#endif

bool IsADReady(Int_t run);

Int_t MakeTrendingADQA(TString QAfilename ="QAresults.root",Int_t runNumber = 244619,TString ocdbStorage = "raw://",Bool_t IsOnGrid = kFALSE,Bool_t printResults = kTRUE, Bool_t doRebin = kTRUE)
{
  
  TTree *ttree=new TTree("trending","tree of trending variables");
  Int_t adReady=0,invalidInput=0,qaNotFound=0,adActive=0,adQANotfound=0,noEntries=0;
  Int_t LHCstate = -1; Float_t runTime = 0; Int_t fillNumber = 0;
  Float_t meanTotalChargeADA = -1024, meanTotalChargeADC = -1024;
  Float_t meanChargeChannelTime[16];
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
  
  Float_t integratedChargeChannel[16];
  Float_t triggerChargeChannel[16];
  Float_t tailChargeChannel[16];
  Float_t integratedChargeChannel_Weighted[16];
  Float_t triggerChargeChannel_Weighted[16];
  Float_t tailChargeChannel_Weighted[16];

  TString treePostFileName="trending.root";

  ttree->Branch("adReady",&adReady,"adReady/I");
  ttree->Branch("invalidInput",&invalidInput,"invalidInput/I");
  ttree->Branch("qaNotFound",&qaNotFound,"qaNotFound/I");
  ttree->Branch("adActive",&adActive,"adActive/I");
  ttree->Branch("adQANotfound",&adQANotfound,"adQANotfound/I");
  ttree->Branch("noEntries",&noEntries,"noEntries/I");
  ttree->Branch("run",&runNumber,"run/I");
  ttree->Branch("fill",&fillNumber,"fill/I");
  if (!QAfilename) 
    {
    if(!IsADReady(runNumber)) adReady=1;
    invalidInput=1;
    Printf("Error - Invalid input file");
    TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
    ttree->Fill();
    trendFile->cd();
    ttree->Write();
    trendFile->Close();
    return 1;
    }
  gStyle->SetPalette(1);
  
  if(IsOnGrid) TGrid::Connect("alien://");
  TFile *QAfile = TFile::Open(QAfilename,"r");

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbStorage);
  man->SetRun(runNumber);
  AliCDBEntry *entry2=0;
  entry2 = man->Get("GRP/GRP/Data");
  AliGRPObject* fGRPData=0;
  if (!entry2) 
    return 0;
  printf("Found an AliGRPObject in GRP/GRP/Data, reading it\n");
  fGRPData = dynamic_cast<AliGRPObject*>(entry2->GetObject()); 
  entry2->SetOwner(0);
    
  TString activeDetList(AliDAQ::ListOfTriggeredDetectors(fGRPData->GetDetectorMask()));
  TString runType(fGRPData->GetRunType());
  TString beamType(fGRPData->GetBeamType());
  TString machineMode(fGRPData->GetMachineMode());
  TString lhcState(fGRPData->GetLHCState());
  
  AliCDBEntry *entry3=0;
  entry3 = man->Get("GRP/GRP/LHCData");
  AliLHCData* fLHCData=0;
  if (!entry3) 
    return 0;
  printf("Found an AliLHCData in GRP/GRP/LHCData, reading it\n");
  fLHCData = dynamic_cast<AliLHCData*>(entry3->GetObject()); 
  entry3->SetOwner(0);
  fillNumber = fLHCData->GetFillNumber();
  
  //-----------------------Trending variables------------------------------------------------
  ttree->Branch("LHCstate",&LHCstate,"LHC state;;State/I");
  ttree->Branch("runTime",&runTime,"runTime;;Duration(min)/F");
  ttree->Branch("meanTotalChargeADA",&meanTotalChargeADA,"Mean total charge ADA;;Charge (ADC counts)/F");
  ttree->Branch("meanTotalChargeADC",&meanTotalChargeADC,"Mean total charge ADC;;Charge (ADC counts)/F");
  ttree->Branch("meanTimeADA",&meanTimeADA,"Mean time ADA;;Time (ns)/F");
  ttree->Branch("meanTimeADC",&meanTimeADC,"Mean time ADC;;Time (ns)/F");
  ttree->Branch("meanChargeChannelTime", &meanChargeChannelTime[0], "meanChargeChannelTime[16]/F");
  ttree->Branch("meanTimeErrADA",&meanTimeErrADA,"Mean time err ADA;;Time (ns)/F");
  ttree->Branch("meanTimeErrADC",&meanTimeErrADC,"Mean time err ADC;;Time (ns)/F");
  ttree->Branch("meanTimeSigmaADA",&meanTimeSigmaADA,"Mean time Sigma ADA;;Time Sigma (ns)/F");
  ttree->Branch("meanTimeSigmaADC",&meanTimeSigmaADC,"Mean time Sigma ADC;;Time Sigma (ns)/F");
  ttree->Branch("meanTimeSigmaErrADA",&meanTimeSigmaErrADA,"Mean time Sigma err ADA;;Time Sigma (ns)/F");
  ttree->Branch("meanTimeSigmaErrADC",&meanTimeSigmaErrADC,"Mean time Sigma err ADC;;Time Sigma (ns)/F");
  ttree->Branch("rateUBA",&rateUBA,"Trigger rate UBA;;Rate/F");
  ttree->Branch("rateUBC",&rateUBC,"Trigger rate UBC;;Rate/F");
  ttree->Branch("rateUGA",&rateUGA,"Trigger rate UGA;;Rate/F");
  ttree->Branch("rateUGC",&rateUGC,"Trigger rate UGC;;Rate/F");
  ttree->Branch("rateADAND",&rateADAND,"Trigger rate ADAND;;Rate/F");
  ttree->Branch("rateADOR",&rateADOR,"Trigger rate ADOR;;Rate/F");
  ttree->Branch("rateRatioADV0AND",&rateRatioADV0AND,"Trigger rate ratio VZEROAND ADAND;;Ratio/F");
  ttree->Branch("rateRatioADV0OR",&rateRatioADV0OR,"Trigger rate ratio VZEROOR ADOR;;Ratio/F");
  ttree->Branch("rateErr",&rateErr,"Trigger rate err;;Rate/F");
  ttree->Branch("MPV", &MPV[0], "MPV[16]/F");
  ttree->Branch("MPVErr", &MPVErr[0], "MPVErr[16]/F");
  ttree->Branch("meanPedestal", &meanPedestal[0], "meanPedestal[32]/F");
  ttree->Branch("widthPedestal", &widthPedestal[0], "widthPedestal[32]/F");
  ttree->Branch("slewingChi2ADA",&slewingChi2ADA,"Time slewing Chi2 ADA;;Chi2perNDF/F");
  ttree->Branch("slewingChi2ADC",&slewingChi2ADC,"Time slewing Chi2 ADC;;Chi2perNDF/F");
  ttree->Branch("saturationADA",&saturationADA,"Saturation ADA;;Saturated fraction/F");
  ttree->Branch("saturationADC",&saturationADC,"Saturation ADC;;Saturated fraction/F");
  ttree->Branch("ratePhysADAND",&ratePhysADAND,"Offline rate ADAND;;Rate/F");
  ttree->Branch("ratePhysADOR",&ratePhysADOR,"Offline rate ADOR;;Rate/F");
  ttree->Branch("ratePhysBBA",&ratePhysBBA,"Offline rate BBA;;Rate/F");
  ttree->Branch("ratePhysBBC",&ratePhysBBC,"Offline rate BBC;;Rate/F");
  ttree->Branch("ratePhysBGA",&ratePhysBGA,"Offline rate BGA;;Rate/F");
  ttree->Branch("ratePhysBGC",&ratePhysBGC,"Offline rate BGC;;Rate/F");
  ttree->Branch("channelTimeMean", &channelTimeMean[0], "channelTimeMean[16]/F");
  ttree->Branch("channelTimeSigma", &channelTimeSigma[0], "channelTimeSigma[16]/F");
  ttree->Branch("channelTimeSigma", &channelTimeSigma[0], "channelTimeSigma[16]/F");
  ttree->Branch("flagNoTimeFraction", &flagNoTimeFraction[0], "flagNoTimeFraction[16]/F");
  ttree->Branch("thresholdData", &thresholdData[0], "thresholdData[16]/F");
  ttree->Branch("thresholdOCDB", &thresholdOCDB[0], "thresholdOCDB[16]/F");
  
  ttree->Branch("integratedChargeChannel", &integratedChargeChannel[0], "integratedChargeChannel[16]/F");
  ttree->Branch("triggerChargeChannel", &triggerChargeChannel[0], "triggerChargeChannel[16]/F");
  ttree->Branch("tailChargeChannel", &tailChargeChannel[0], "tailChargeChannel[16]/F");
  
  ttree->Branch("integratedChargeChannel_Weighted", &integratedChargeChannel_Weighted[0], "integratedChargeChannel_Weighted[16]/F");
  ttree->Branch("triggerChargeChannel_Weighted", &triggerChargeChannel_Weighted[0], "triggerChargeChannel_Weighted[16]/F");
  ttree->Branch("tailChargeChannel_Weighted", &tailChargeChannel_Weighted[0], "tailChargeChannel_Weighted[16]/F");
  

  if(!QAfile)
    {
    if(!IsADReady(runNumber)) adReady=1;
    qaNotFound=1;
    Printf("ERROR: QA output not found. Exiting ...\n");
    TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
    ttree->Fill();
    trendFile->cd();
    ttree->Write();
    trendFile->Close();
    return -1;
    }
  else
    {
    Printf("INFO: QA output file %s open. \n",QAfile->GetName());
    }

  printf("activeDetList %s\nrunType %s\nbeamType %s\nmachineMode %s\nlhcState %s\n",
	 activeDetList.Data(),runType.Data(),beamType.Data(),
	 machineMode.Data(),lhcState.Data());
  
  time_t duration = fGRPData->GetTimeEnd() - fGRPData->GetTimeStart();
  runTime = (Float_t)(duration/60);

  if(!activeDetList.Contains("AD"))
    { 
    if(!IsADReady(runNumber)) adReady=1;
    adActive=1;
    printf("RUN WITH AD NOT ACTIVE\n");
    TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
    ttree->Fill();
    trendFile->cd();
    ttree->Write();
    trendFile->Close();
    
    return 0;
    }
  if(!runType.Contains("PHYSICS"))
    { 
    printf("RUN NO PHYSICS\n");
    return 0;
    }
  if(lhcState.Contains("NO BEAM"))LHCstate = 0;
  if(lhcState.Contains("SQUEEZE"))LHCstate = 1;
  if(lhcState.Contains("ADJUST"))LHCstate = 2;
  if(lhcState.Contains("STABLE BEAMS"))LHCstate = 3;
      
  char adQAdirName[5]="ADQA";
  QAfile->Print();
  TDirectoryFile* adQAdir=(TDirectoryFile*)QAfile->Get(adQAdirName);
  if(!adQAdir)
    {
    if(!IsADReady(runNumber)) adReady=1;
    else adQANotfound=1;
    printf("ERROR: AD QA directory not present in input file.\n");
    TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
    ttree->Fill();
    trendFile->cd();
    ttree->Write();
    trendFile->Close();
  
    return -1;
    }
  TList *flistQA = (TList*)adQAdir->Get("ADQAListHist");
  if(!flistQA) 
    {
    cout << "ERROR: No list found" << endl;
    return -1;
    }
  //-----------------------QA histograms------------------------------------------------
  TH1F *fHistTotalChargePerEventADA = dynamic_cast<TH1F*> (flistQA->FindObject("fHistTotalChargePerEventADA"));
  TH1F *fHistTotalChargePerEventADC = dynamic_cast<TH1F*> (flistQA->FindObject("fHistTotalChargePerEventADC"));

  TH2F *fHistChargePerPM_All = dynamic_cast<TH2F*> (flistQA->FindObject("fHistChargePerPM_All"));
  TH2F *fHistChargePerPM_BB = dynamic_cast<TH2F*> (flistQA->FindObject("fHistChargePerPM_BB"));
  TH2F *fHistChargePerPM_BG = dynamic_cast<TH2F*> (flistQA->FindObject("fHistChargePerPM_BG"));
  TH2F *fHistChargePerPM_Time = dynamic_cast<TH2F*> (flistQA->FindObject("fHistChargePerPM_Time"));

  TH2F *fHistTimePerPM_Corr = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTimePerPM_Corr"));
  TH2F *fHistTimePerPM_UnCorr = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTimePerPM_UnCorr"));

  TH2F *fHistWidthPerPM = dynamic_cast<TH2F*> (flistQA->FindObject("fHistWidthPerPM"));

  TH2F *fHistTimeVsChargeADA_Cut = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTimeVsChargeADA_Cut"));
  TH2F *fHistTimeVsChargeADA_Corr = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTimeVsChargeADA_Corr"));
  TH2F *fHistTimeVsChargeADA_UnCorr = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTimeVsChargeADA_UnCorr"));
  TH2F *fHistTimeVsChargeADC_Cut = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTimeVsChargeADC_Cut"));
  TH2F *fHistTimeVsChargeADC_Corr = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTimeVsChargeADC_Corr"));
  TH2F *fHistTimeVsChargeADC_UnCorr = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTimeVsChargeADC_UnCorr"));

  TH2F *fHistWidthVsCharge = dynamic_cast<TH2F*> (flistQA->FindObject("fHistWidthVsCharge"));

  TH1F *fHistNBBflagsADA = dynamic_cast<TH1F*> (flistQA->FindObject("fHistNBBflagsADA"));
  TH1F *fHistNBBflagsADC = dynamic_cast<TH1F*> (flistQA->FindObject("fHistNBBflagsADC"));
  TH2F *fHistNBBflagsADAVsADC = dynamic_cast<TH2F*> (flistQA->FindObject("fHistNBBflagsADAVsADC"));

  TH1F *fHistNBGflagsADA = dynamic_cast<TH1F*> (flistQA->FindObject("fHistNBGflagsADA"));
  TH1F *fHistNBGflagsADC = dynamic_cast<TH1F*> (flistQA->FindObject("fHistNBGflagsADC"));
  TH2F *fHistNBGflagsADAVsADC = dynamic_cast<TH2F*> (flistQA->FindObject("fHistNBGflagsADAVsADC"));

  TH1F *fHistNBBCoincidencesADA = dynamic_cast<TH1F*> (flistQA->FindObject("fHistNBBCoincidencesADA"));
  TH1F *fHistNBBCoincidencesADC = dynamic_cast<TH1F*> (flistQA->FindObject("fHistNBBCoincidencesADC"));
  TH2F *fHistNBBCoincidencesADAVsADC = dynamic_cast<TH2F*> (flistQA->FindObject("fHistNBBCoincidencesADAVsADC"));

  TH1F *fHistNBGCoincidencesADA = dynamic_cast<TH1F*> (flistQA->FindObject("fHistNBGCoincidencesADA"));
  TH1F *fHistNBGCoincidencesADC = dynamic_cast<TH1F*> (flistQA->FindObject("fHistNBGCoincidencesADC"));
  TH2F *fHistNBGCoincidencesADAVsADC = dynamic_cast<TH2F*> (flistQA->FindObject("fHistNBGCoincidencesADAVsADC"));

  TH1F *fHistChargeNoFlag = dynamic_cast<TH1F*> (flistQA->FindObject("fHistChargeNoFlag"));
  TH2F *fHistTimeNoFlag = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTimeNoFlag"));
  TH2F *fHistChargeNoTime = dynamic_cast<TH2F*> (flistQA->FindObject("fHistChargeNoTime"));
  TH1F *fHistFlagNoTime = dynamic_cast<TH1F*> (flistQA->FindObject("fHistFlagNoTime"));

  TH2F *fHistChargePerCoincidence = dynamic_cast<TH2F*> (flistQA->FindObject("fHistChargePerCoincidence"));

  TH1F *fHistMeanTimeADA = dynamic_cast<TH1F*> (flistQA->FindObject("fHistMeanTimeADA"));
  TH1F *fHistMeanTimeADC = dynamic_cast<TH1F*> (flistQA->FindObject("fHistMeanTimeADC"));
  TH1F *fHistMeanTimeDifference = dynamic_cast<TH1F*> (flistQA->FindObject("fHistMeanTimeDifference"));

  TH2F *fHistMeanTimeCorrelation = dynamic_cast<TH2F*> (flistQA->FindObject("fHistMeanTimeCorrelation"));
  TH2F *fHistMeanTimeSumDiff = dynamic_cast<TH2F*> (flistQA->FindObject("fHistMeanTimeSumDiff"));

  TH2F *fHistDecision = dynamic_cast<TH2F*> (flistQA->FindObject("fHistDecision"));

  TH1F *fHistTriggerMasked = dynamic_cast<TH1F*> (flistQA->FindObject("fHistTriggerMasked"));
  TH1F *fHistTriggerOthers = dynamic_cast<TH1F*> (flistQA->FindObject("fHistTriggerOthers"));

  TH2F *fHistChargeVsClockInt0 = dynamic_cast<TH2F*> (flistQA->FindObject("fHistChargeVsClockInt0"));
  TH2F *fHistChargeVsClockInt1 = dynamic_cast<TH2F*> (flistQA->FindObject("fHistChargeVsClockInt1"));
  TH2F *fHistBBFlagVsClock = dynamic_cast<TH2F*> (flistQA->FindObject("fHistBBFlagVsClock"));
  TH2F *fHistBGFlagVsClock = dynamic_cast<TH2F*> (flistQA->FindObject("fHistBGFlagVsClock"));
  TH2F *fHistBBFlagPerChannel = dynamic_cast<TH2F*> (flistQA->FindObject("fHistBBFlagPerChannel"));
  TH2F *fHistBGFlagPerChannel = dynamic_cast<TH2F*> (flistQA->FindObject("fHistBGFlagPerChannel"));
  TH2F *fHistMaxChargeClock = dynamic_cast<TH2F*> (flistQA->FindObject("fHistMaxChargeClock"));
  TH2F *fHistMaxChargeValueInt0 = dynamic_cast<TH2F*> (flistQA->FindObject("fHistMaxChargeValueInt0"));
  TH2F *fHistMaxChargeValueInt1 = dynamic_cast<TH2F*> (flistQA->FindObject("fHistMaxChargeValueInt1"));
  
  TH3F *fHistTimeVsChargePerPM_UnCorr = dynamic_cast<TH3F*> (flistQA->FindObject("fHistTimeVsChargePerPM_UnCorr"));
  
  //Gain mon cosmics
  TH2F *fHistChargeTriggerPerChannel_PF = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTriggerChargePerChannel_PF"));
  //Gain mon collisions
  TH2F *fHistChargeTriggerPerChannel_PF_TVX = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTriggerChargePerChannel_TVX"));
  TH2F *fHistChargeTailPerChannel_PF_TVX = dynamic_cast<TH2F*> (flistQA->FindObject("fHistTailChargePerChannel_TVX"));
  TH2F *fHistChargePerPM_BB_PF_TVX = dynamic_cast<TH2F*> (flistQA->FindObject("fHistIntegratedChargePerChannel_TVX"));
  TH3F *fHistChargeTriggerPerPMPerV0Flag = dynamic_cast<TH3F*> (flistQA->FindObject("fHistTriggerChargePerPMPerV0Flag"));
  TH3F *fHistChargeTailPerPMPerV0Flag = dynamic_cast<TH3F*> (flistQA->FindObject("fHistTailChargePerPMPerV0Flag"));
  TH3F *fHistChargeBBPerPMPerV0Flag = dynamic_cast<TH3F*> (flistQA->FindObject("fHistIntegratedChargePerPMPerV0Flag"));

  if(fHistTriggerMasked->GetEntries()==0)
    {
    if(!IsADReady(runNumber)) adReady=1;
    noEntries=1;
    TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
    ttree->Fill();
    trendFile->cd();
    ttree->Write();
    trendFile->Close();
  
    printf("fHistTriggerMasked no entries\n");
    return 0;
    }  
  //-----------------------Fill trending variables------------------------------------------------
  meanTotalChargeADA = fHistTotalChargePerEventADA->GetMean();
  meanTotalChargeADC = fHistTotalChargePerEventADC->GetMean();
    
  Int_t minFitRange = 50;
  Int_t maxFitRange = 65;
  Int_t fitStatus; 
  
  fitStatus = fHistMeanTimeADA->Fit("gaus","R+","",minFitRange,maxFitRange);
  if(fitStatus ==0){
  	TF1 *fitTimeADA = (TF1*) fHistMeanTimeADA->GetFunction("gaus");
  	meanTimeADA = fitTimeADA->GetParameter(1);
  	meanTimeSigmaADA = fitTimeADA->GetParameter(2);
 	meanTimeErrADA = fitTimeADA->GetParError(1);
  	meanTimeSigmaErrADA = fitTimeADA->GetParError(2);
	}
	
  TH1D *hTimeChannelSlice[16];
  for(Int_t i = 8; i<16; i++){
  	TString sliceName = "hTimeChannelSlice";
	sliceName += i;
	TString sliceTitle = "Corrected time, ch";
	sliceTitle += i;
  	hTimeChannelSlice[i] = fHistTimePerPM_Corr->ProjectionY(sliceName.Data(),i+1,i+1);
	hTimeChannelSlice[i]->SetTitle(sliceTitle.Data());
	fitStatus = hTimeChannelSlice[i]->Fit("gaus","R+","",minFitRange,maxFitRange);
	if(fitStatus ==0){
  		TF1 *fitTimeChannel = (TF1*) hTimeChannelSlice[i]->GetFunction("gaus");
  		channelTimeMean[i] = fitTimeChannel->GetParameter(1);
  		channelTimeSigma[i] = fitTimeChannel->GetParameter(2);
	}
  }
  
  minFitRange = 62;
  maxFitRange = 68;
  
  fitStatus = fHistMeanTimeADC->Fit("gaus","R+","",minFitRange,maxFitRange);
  if(fitStatus ==0){
  	TF1 *fitTimeADC = (TF1*) fHistMeanTimeADC->GetFunction("gaus");
  	meanTimeADC = fitTimeADC->GetParameter(1);
  	meanTimeSigmaADC = fitTimeADC->GetParameter(2);
  	meanTimeErrADC = fitTimeADC->GetParError(1);
  	meanTimeSigmaErrADC = fitTimeADC->GetParError(2);
	}
	
  for(Int_t i = 0; i<8; i++){
  	TString sliceName = "hTimeChannelSlice";
	sliceName += i;
	TString sliceTitle = "Corrected time, ch";
	sliceTitle += i;
  	hTimeChannelSlice[i] = fHistTimePerPM_Corr->ProjectionY(sliceName.Data(),i+1,i+1);
	hTimeChannelSlice[i]->SetTitle(sliceTitle.Data());
	fitStatus = hTimeChannelSlice[i]->Fit("gaus","R+","",minFitRange,maxFitRange);
	if(fitStatus ==0){
  		TF1 *fitTimeChannel = (TF1*) hTimeChannelSlice[i]->GetFunction("gaus");
  		channelTimeMean[i] = fitTimeChannel->GetParameter(1);
  		channelTimeSigma[i] = fitTimeChannel->GetParameter(2);
	}
  }
  
  const Double_t fTOFADA = 56.63;
  TH1D *hSlewingSlice;
  slewingChi2ADA = 0;
  for(Int_t i = 0; i<550; i++){
  	hSlewingSlice = fHistTimeVsChargeADA_Cut->ProjectionX("hSlewingSlice",i+1,i+1);
	if(hSlewingSlice->Integral() < 100) continue;
	slewingChi2ADA += TMath::Power(hSlewingSlice->GetMean() - fTOFADA, 2);
  	}
  slewingChi2ADA = slewingChi2ADA/550;
  
  const Double_t fTOFADC = 65.21;
  TH1D *hSlewingSlice;
  slewingChi2ADC = 0;
  for(Int_t i = 0; i<550; i++){
  	hSlewingSlice = fHistTimeVsChargeADC_Cut->ProjectionX("hSlewingSlice",i+1,i+1);
	if(hSlewingSlice->Integral() < 100) continue;
	slewingChi2ADC += TMath::Power(hSlewingSlice->GetMean() - fTOFADC, 2);
  	}
  slewingChi2ADC = slewingChi2ADC/550;
  
  Float_t nEvents = fHistTotalChargePerEventADA->GetEntries();
  rateUBA = fHistTriggerMasked->GetBinContent(7)/nEvents;
  rateUBC = fHistTriggerMasked->GetBinContent(8)/nEvents;
  rateUGA = fHistTriggerMasked->GetBinContent(4)/nEvents;
  rateUGC = fHistTriggerMasked->GetBinContent(6)/nEvents;
  rateADAND = fHistTriggerMasked->GetBinContent(1)/nEvents;
  rateADOR = fHistTriggerMasked->GetBinContent(2)/nEvents;
  if(fHistTriggerMasked->GetBinContent(1) != 0)rateRatioADV0AND = fHistTriggerOthers->GetBinContent(1)/fHistTriggerMasked->GetBinContent(1);
  if(fHistTriggerMasked->GetBinContent(2) != 0)rateRatioADV0OR = fHistTriggerOthers->GetBinContent(2)/fHistTriggerMasked->GetBinContent(2);
  rateErr = 2/TMath::Sqrt(nEvents);
  
  //cout<<"Rate ADand = "<<rateADAND<<" Rate ADor = "<<rateADOR<<" Ratio AND/OR = "<<rateADAND/rateADOR<<endl;
  
  ratePhysADAND = fHistDecision->GetBinContent(2,2)/nEvents;
  ratePhysADOR = 0;
  for(Int_t i=1;i<=4;i++){ 
  	ratePhysADOR += fHistDecision->GetBinContent(i,2)/nEvents;
	ratePhysADOR += fHistDecision->GetBinContent(2,i)/nEvents;
	}
  ratePhysADOR -= ratePhysADAND;
  
  ratePhysBBA = 0; ratePhysBBC = 0; ratePhysBGA = 0; ratePhysBGC = 0;
  for(Int_t i=1;i<=4;i++){ 
  	ratePhysBBA += fHistDecision->GetBinContent(2,i)/nEvents;
	ratePhysBGA += fHistDecision->GetBinContent(3,i)/nEvents;
	ratePhysBBC += fHistDecision->GetBinContent(i,2)/nEvents;
	ratePhysBGC += fHistDecision->GetBinContent(i,3)/nEvents;
	}
  
  //cout<<"Phys sel ADand = "<<ratePhysADAND<<" Phys sel ADor = "<<ratePhysADOR<<" Ratio AND/OR = "<<ratePhysADAND/ratePhysADOR<<endl;
  
  for(Int_t i = 0; i<16; i++) flagNoTimeFraction[i] = fHistFlagNoTime->GetBinContent(i+1)/nEvents;
  
  TH1D *hChargeSliceAll[16];
  TH1D *hChargeSliceTime[16];
  TH1D *hChargeSliceBB[16];
  TH1D *hThresholdShape[16];
  
  for(Int_t i = 0; i<16; i++){
  	TString channelNameAll = "hChargeSliceAll";
	channelNameAll += i;
  	hChargeSliceAll[i] = fHistChargePerPM_All->ProjectionY(channelNameAll.Data(),i+1,i+1);
	
	TString channelNameTime = "hChargeSliceTime";
	channelNameTime += i;
  	hChargeSliceTime[i] = fHistChargePerPM_Time->ProjectionY(channelNameTime.Data(),i+1,i+1);
	meanChargeChannelTime[i] = hChargeSliceTime[i]->GetMean();
	
	TString channelNameBB = "hChargeSliceBB";
	channelNameBB += i;
  	hChargeSliceBB[i] = fHistChargePerPM_BB->ProjectionY(channelNameBB.Data(),i+1,i+1);
	
	TString channelTitle = "Integrated charge, PM ";
	channelTitle += i;
	hChargeSliceAll[i]->SetTitle(channelTitle.Data());
	
	TString thresholdShapeName = "hThresholdShape";
	thresholdShapeName += i;
	hThresholdShape[i] = (TH1D*)hChargeSliceTime[i]->Clone(thresholdShapeName.Data());
	hThresholdShape[i]->Divide(hChargeSliceAll[i]);
	
	TString thresholdTitle = "Threshold shape, PM ";
	thresholdTitle += i;
	hThresholdShape[i]->SetTitle(thresholdTitle.Data());
	}

  AliCDBEntry *entCD = man->Get("AD/Calib/Data");
  AliADCalibData *fCalibData = (AliADCalibData*)entCD->GetObject();
	
  TF1 *thrFit = new TF1("thrFit","1/(1+TMath::Exp(-x + [0]))",1,50);
  for(Int_t i = 0; i<16; i++){
  	thrFit->SetParameter(0,fCalibData->GetCalibDiscriThr(i));
  	fitStatus = hThresholdShape[i]->Fit("thrFit","+R");
	if(fitStatus == 0) thresholdData[i] = thrFit->GetParameter(0);
	thresholdOCDB[i] = fCalibData->GetCalibDiscriThr(i);
  	}
  
  for(Int_t i = 0; i<32; i++){
  	meanPedestal[i] = fCalibData->GetPedestal(i);
	widthPedestal[i] = fCalibData->GetSigma(i);
	}

  TH1D *hMaxChargeSlice;
  
  Float_t sat=0, nonsat=0;
  for(Int_t i = 0; i<8; i++){
  	hMaxChargeSlice = fHistMaxChargeValueInt0->ProjectionY("hMaxChargeSlice",i+1,i+1);
	sat += hMaxChargeSlice->Integral(1000,1025);
	nonsat += hMaxChargeSlice->Integral(meanPedestal[i]+3*widthPedestal[i],1025);
	hMaxChargeSlice = fHistMaxChargeValueInt1->ProjectionY("hMaxChargeSlice",i+1,i+1);
	sat += hMaxChargeSlice->Integral(1000,1025);
	nonsat += hMaxChargeSlice->Integral(meanPedestal[i]+3*widthPedestal[i],1025);
	}
  if(nonsat != 0)saturationADC = sat/nonsat;
  
  sat=0; nonsat=0;
  for(Int_t i = 8; i<16; i++){
  	hMaxChargeSlice = fHistMaxChargeValueInt0->ProjectionY("hMaxChargeSlice",i+1,i+1);
	sat += hMaxChargeSlice->Integral(1000,1025);
	nonsat += hMaxChargeSlice->Integral(meanPedestal[i]+3*widthPedestal[i],1025);
	hMaxChargeSlice = fHistMaxChargeValueInt1->ProjectionY("hMaxChargeSlice",i+1,i+1);
	sat += hMaxChargeSlice->Integral(1000,1025);
	nonsat += hMaxChargeSlice->Integral(meanPedestal[i]+3*widthPedestal[i],1025);
	}
  if(nonsat != 0)saturationADA = sat/nonsat;
  
  
  //Gain mon cosmics
  TF1 *fitMPV = new TF1("fitMPV","landau",0,5);
  Double_t xq[1] = {0.80};
  Double_t yq[1];
  
  TH1D* hCTSlice;
  for(Int_t i = 0; i<16; i++){
	hCTSlice = fHistChargeTriggerPerChannel_PF->ProjectionY("hCTSlice",i+1,i+1);	

	//hCTSlice->Fit(fitMPV);
  	//MPV[i] = fitMPV->GetParameter(1);
	//MPVErr[i] = fitMPV->GetParError(1);
	//MPV[i] = hCTSlice->GetMean();
	//hCTSlice->GetQuantiles(1,yq,xq);
	//MPV[i] = yq[0];
	MPV[i] = fCalibData->GetADCperMIP(i);
	MPVErr[i] = 0.0;
  	}
	
  //Gain mon collisions	 
  for(Int_t i = 0; i<16; i++){
	hCTSlice = fHistChargeTriggerPerChannel_PF_TVX->ProjectionY("hCTSlice",i+1,i+1);	
	hCTSlice->GetQuantiles(1,yq,xq);
	triggerChargeChannel[i] = yq[0];
	//cout<<"Quantile unweighted = "<<triggerChargeChannel[i]<<endl;
  	}
  for(Int_t i = 0; i<16; i++){
	hCTSlice = fHistChargeTailPerChannel_PF_TVX->ProjectionY("hCTSlice",i+1,i+1);	
	hCTSlice->GetQuantiles(1,yq,xq);
	tailChargeChannel[i] = yq[0];
	//cout<<"Quantile unweighted = "<<tailChargeChannel[i]<<endl;
  	}
   for(Int_t i = 0; i<16; i++){
	hCTSlice = fHistChargePerPM_BB_PF_TVX->ProjectionY("hCTSlice",i+1,i+1);	
	hCTSlice->GetQuantiles(1,yq,xq);
	integratedChargeChannel[i] = yq[0];
	//cout<<"Quantile unweighted = "<<integratedChargeChannel[i]<<endl;
  	}

   TH1D* hVZEROFlags;
   TH1D* hWeighted;	
   for(Int_t i = 0; i<16; i++){
   	hVZEROFlags = dynamic_cast<TH1D*>fHistChargeTriggerPerPMPerV0Flag->Project3D("z");
	if(hVZEROFlags->Integral() != 0)hVZEROFlags->Scale(1/hVZEROFlags->Integral());
	hWeighted = dynamic_cast<TH1D*>fHistChargeTriggerPerPMPerV0Flag->Project3D("y");
	hWeighted->SetName("hWeighted");
	hWeighted->Reset("ICES");
	for(Int_t j = 0; j<64; j++){
		hCTSlice = fHistChargeTriggerPerPMPerV0Flag->ProjectionY("hCTSlice",i+1,i+1,j+1,j+1);
		Float_t weight = 0.0;
		if(hVZEROFlags->GetBinContent(j+1) != 0) weight = 1/(64*hVZEROFlags->GetBinContent(j+1));
		//if(i==0)cout<<"Weight "<<j<<" = "<<weight<<endl;
		hWeighted->Add(hCTSlice,weight);
		}
	
	hWeighted->GetQuantiles(1,yq,xq);
	triggerChargeChannel_Weighted[i] = yq[0];
	//cout<<"Quantile weighted = "<<triggerChargeChannel_Weighted[i]<<endl;
  	}
   for(Int_t i = 0; i<16; i++){
   	hVZEROFlags = dynamic_cast<TH1D*>fHistChargeTailPerPMPerV0Flag->Project3D("z");
	if(hVZEROFlags->Integral() != 0)hVZEROFlags->Scale(1/hVZEROFlags->Integral());
	hWeighted = dynamic_cast<TH1D*>fHistChargeTailPerPMPerV0Flag->Project3D("y");
	hWeighted->SetName("hWeighted");
	hWeighted->Reset("ICES");
	for(Int_t j = 0; j<64; j++){
		hCTSlice = fHistChargeTailPerPMPerV0Flag->ProjectionY("hCTSlice",i+1,i+1,j+1,j+1);
		Float_t weight = 0.0;
		if(hVZEROFlags->GetBinContent(j+1) != 0) weight = 1/(64*hVZEROFlags->GetBinContent(j+1));
		//if(i==0)cout<<"Weight "<<j<<" = "<<weight<<endl;
		hWeighted->Add(hCTSlice,weight);
		}
	
	hWeighted->GetQuantiles(1,yq,xq);
	tailChargeChannel_Weighted[i] = yq[0];
	//cout<<"Quantile weighted = "<<tailChargeChannel_Weighted[i]<<endl;
  	}
   for(Int_t i = 0; i<16; i++){
   	hVZEROFlags = dynamic_cast<TH1D*>fHistChargeBBPerPMPerV0Flag->Project3D("z");
	if(hVZEROFlags->Integral() != 0)hVZEROFlags->Scale(1/hVZEROFlags->Integral());
	hWeighted = dynamic_cast<TH1D*>fHistChargeBBPerPMPerV0Flag->Project3D("y");
	hWeighted->SetName("hWeighted");
	hWeighted->Reset("ICES");
	for(Int_t j = 0; j<64; j++){
		hCTSlice = fHistChargeBBPerPMPerV0Flag->ProjectionY("hCTSlice",i+1,i+1,j+1,j+1);
		Float_t weight = 0.0;
		if(hVZEROFlags->GetBinContent(j+1) != 0) weight = 1/(64*hVZEROFlags->GetBinContent(j+1));
		//if(i==0)cout<<"Weight "<<j<<" = "<<weight<<endl;
		hWeighted->Add(hCTSlice,weight);
		}
	
	hWeighted->GetQuantiles(1,yq,xq);
	integratedChargeChannel_Weighted[i] = yq[0];
	//cout<<"Quantile weighted = "<<integratedChargeChannel_Weighted[i]<<endl;
  	}		
  
  
  TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
  ttree->Fill();
  trendFile->cd();
  ttree->Write();
  trendFile->Close();
  
  //-----------------------Print QA histos---------------------------------------------------------
  if(printResults)
    { 
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

    TCanvas *c1 = new TCanvas("TotalChargePerEvent"," ",800,400); 
    c1->Draw(); 							   
    c1->cd();
    TPad *myPad1 = new TPad("myPad1", "The pad",0,0,1,1);
    myPadSetUp(myPad1,0.15,0.1,0.04,0.15);
    myPad1->SetLogy();
    myPad1->Draw();
    myPad1->cd();
    myHistSetUp(fHistTotalChargePerEventADA);
    myHistSetUp(fHistTotalChargePerEventADC);
    if(doRebin)fHistTotalChargePerEventADA->Rebin(10);
    if(doRebin)fHistTotalChargePerEventADC->Rebin(10);
    fHistTotalChargePerEventADA->SetLineColor(kBlue);
    fHistTotalChargePerEventADC->SetLineColor(kRed);
    fHistTotalChargePerEventADA->SetLineWidth(2);
    fHistTotalChargePerEventADC->SetLineWidth(2);
    fHistTotalChargePerEventADA->Draw();
    fHistTotalChargePerEventADC->Draw("same");

    TLegend *myLegend1 = new TLegend(0.70,0.67,0.97,0.82);
    myLegendSetUp(myLegend1,0.04,1);
    myLegend1->AddEntry(fHistTotalChargePerEventADA,"ADA","l");
    myLegend1->AddEntry(fHistTotalChargePerEventADC,"ADC","l");
    myLegend1->Draw();

    c1->Print(Form("ADQA_Run_%d.pdf(",runNumber));
    
    myHistSetUp(fHistChargePerPM_All);
    myHistSetUp(fHistChargePerPM_Time);
    myHistSetUp(fHistChargePerPM_BB);
    myHistSetUp(fHistChargePerPM_BG);

    TCanvas *c21 = new TCanvas("ChargePerPM_All"," ",1500,500);
    c21->Draw();
    c21->cd();
    TPad *myPad21 = new TPad("myPad21", "The pad",0,0,1,1);
    myPad21->Divide(2,1);
    myPad21->Draw();
    myPadSetUp(myPad21->cd(1),0.15,0.15,0.15,0.15);
    myPadSetUp(myPad21->cd(2),0.15,0.15,0.15,0.15);
    myPad21->cd(2);
    gPad->SetLogz();
    fHistChargePerPM_All->GetYaxis()->SetRangeUser(1,50);
    fHistChargePerPM_All->DrawCopy("COLZ");
    myPad21->cd(1);
    gPad->SetLogz();
    fHistChargePerPM_All->GetYaxis()->SetRangeUser(-1,-1);
    if(doRebin)fHistChargePerPM_All->RebinY(10);
    fHistChargePerPM_All->DrawCopy("COLZ");
    c21->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    TCanvas *c25 = new TCanvas("ChargePerPM_Time"," ",1500,500);
    c25->Draw();
    c25->cd();
    TPad *myPad25 = new TPad("myPad25", "The pad",0,0,1,1);
    myPad25->Divide(2,1);
    myPad25->Draw();
    myPadSetUp(myPad25->cd(1),0.15,0.15,0.15,0.15);
    myPadSetUp(myPad25->cd(2),0.15,0.15,0.15,0.15);
    myPad25->cd(2);
    gPad->SetLogz();
    fHistChargePerPM_Time->GetYaxis()->SetRangeUser(1,50);
    fHistChargePerPM_Time->DrawCopy("COLZ");
    myPad25->cd(1);
    gPad->SetLogz();
    fHistChargePerPM_Time->GetYaxis()->SetRangeUser(-1,-1);
    if(doRebin)fHistChargePerPM_Time->RebinY(10);
    fHistChargePerPM_Time->DrawCopy("COLZ"); 
    c25->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    TCanvas *c22 = new TCanvas("ChargePerPM_BB"," ",1500,500);
    c22->Draw();
    c22->cd();
    TPad *myPad22 = new TPad("myPad22", "The pad",0,0,1,1);
    myPad22->Divide(2,1);
    myPad22->Draw();
    myPadSetUp(myPad22->cd(1),0.15,0.15,0.15,0.15);
    myPadSetUp(myPad22->cd(2),0.15,0.15,0.15,0.15);
    myPad22->cd(2);
    gPad->SetLogz();
    fHistChargePerPM_BB->GetYaxis()->SetRangeUser(1,50);
    fHistChargePerPM_BB->DrawCopy("COLZ");
    myPad22->cd(1);
    gPad->SetLogz();
    fHistChargePerPM_BB->GetYaxis()->SetRangeUser(-1,-1);
    if(doRebin)fHistChargePerPM_BB->RebinY(10);
    fHistChargePerPM_BB->DrawCopy("COLZ");
    c22->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    TCanvas *c23 = new TCanvas("ChargePerPM_BG"," ",1500,500);
    c23->Draw();
    c23->cd();
    TPad *myPad23 = new TPad("myPad23", "The pad",0,0,1,1);
    myPad23->Divide(2,1);
    myPad23->Draw();
    myPadSetUp(myPad23->cd(1),0.15,0.15,0.15,0.15);
    myPadSetUp(myPad23->cd(2),0.15,0.15,0.15,0.15);
    myPad23->cd(2);
    gPad->SetLogz();
    fHistChargePerPM_BG->GetYaxis()->SetRangeUser(1,50);
    fHistChargePerPM_BG->DrawCopy("COLZ");
    myPad23->cd(1);
    gPad->SetLogz();
    fHistChargePerPM_BG->GetYaxis()->SetRangeUser(-1,-1);
    if(doRebin)fHistChargePerPM_BG->RebinY(10);
    fHistChargePerPM_BG->DrawCopy("COLZ");

    c23->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    TCanvas *c24 = new TCanvas("ChargeLandau"," ",1500,800);
    c24->Draw();
    c24->cd();
    TPad *myPad24 = new TPad("myPad24", "The pad",0,0,1,1);
    myPad24->Divide(4,4);
    myPad24->Draw();
    
    TLegend *myLegend3 = new TLegend(0.36,0.62,0.62,0.84);
    myLegendSetUp(myLegend3,0.08,1);
    myLegend3->AddEntry(hChargeSliceAll[0],"All events","f");
    myLegend3->AddEntry(hChargeSliceTime[0],"Events with time","f");
    myLegend3->AddEntry(hChargeSliceBB[0],"Events with BB flag","f");
    
    for(Int_t i=0; i<16; i++){
    	myPadSetUp(myPad24->cd(i+1),0.00,0.00,0.00,0.15);
    	myPad24->cd(i+1);
    	gPad->SetLogy();
	myHistSetUp(hChargeSliceAll[i]);
	myHistSetUp(hChargeSliceTime[i]);
	myHistSetUp(hChargeSliceBB[i]);
	myScaleSetUp(hChargeSliceAll[i],hChargeSliceTime[i]);
    	hChargeSliceAll[i]->GetXaxis()->SetRangeUser(0,50);
	hChargeSliceAll[i]->SetFillStyle(3001);
	hChargeSliceAll[i]->SetFillColor(11);
	hChargeSliceTime[i]->SetFillStyle(3001);
	hChargeSliceTime[i]->SetFillColor(kYellow);
	hChargeSliceBB[i]->SetFillStyle(3001);
	hChargeSliceBB[i]->SetFillColor(kRed);
	hChargeSliceAll[i]->Draw();
	hChargeSliceTime[i]->Draw("same");
	hChargeSliceBB[i]->Draw("same");
	myLegend3->Draw();
	}
 
    c24->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    /*/
    TCanvas *c25 = new TCanvas("ThresholdShape"," ",1500,800);
    c25->Draw();
    c25->cd();
    TPad *myPad25 = new TPad("myPad25", "The pad",0,0,1,1);
    myPad25->Divide(4,4);
    myPad25->Draw();
    
    TF1 *thrModel = new TF1("thrModel","1/(1+TMath::Exp(-x + [0]))",1,50);
    thrModel->SetLineColor(kGreen);
    
    TLegend *myLegend4 = new TLegend(0.48,0.26,0.73,0.47);
    myLegendSetUp(myLegend4,0.08,1);
    myLegend4->AddEntry(hThresholdShape[0],"With_Time/All","l");
    myLegend4->AddEntry(thrModel,"Model","l");
    
    for(Int_t i=0; i<16; i++){
    	myPadSetUp(myPad25->cd(i+1),0.00,0.00,0.00,0.15);
    	myPad25->cd(i+1);
    	gPad->SetLogy();
	myHistSetUp(hThresholdShape[i]);
    	hThresholdShape[i]->GetXaxis()->SetRangeUser(0,50);
	hThresholdShape[i]->SetFillStyle(3001);
	hThresholdShape[i]->SetFillColor(kYellow);
	hThresholdShape[i]->Draw();
	thrModel->SetParameter(0,fCalibData->GetCalibDiscriThr(i));
	thrModel->DrawCopy("same");
	myLegend4->Draw();
	}
 
    c25->Print(Form("ADQA_Run_%d.pdf",runNumber));
    /*/
    
    myHistSetUp(fHistTimePerPM_Corr);
    myHistSetUp(fHistTimePerPM_UnCorr);

    TCanvas *c3 = new TCanvas("LeadingTimePerPM"," ",1000,500);
    c3->Draw();
    c3->cd();
    TPad *myPad3 = new TPad("myPad3", "The pad",0,0,1,1);
    myPadSetUp(myPad3,0.15,0.1,0.1,0.15);
    myPad3->Divide(2,1);
    myPad3->Draw();
    myPadSetUp(myPad3->cd(1),0.15,0.15,0.15,0.15);
    myPadSetUp(myPad3->cd(2),0.15,0.15,0.15,0.15);
    myPad3->cd(1);
    gPad->SetLogz();
    if(doRebin)fHistTimePerPM_Corr->RebinY(2);
    fHistTimePerPM_Corr->Draw("COLZ");
    myPad3->cd(2);
    gPad->SetLogz();
    if(doRebin)fHistTimePerPM_UnCorr->RebinY(2);
    fHistTimePerPM_UnCorr->Draw("COLZ");

    c3->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    TCanvas *c31 = new TCanvas("TimeFitChannels"," ",1500,800);
    c31->Draw();
    c31->cd();
    TPad *myPad31 = new TPad("myPad31", "The pad",0,0,1,1);
    myPad31->Divide(4,4);
    myPad31->Draw();
    TLatex latex;
    
    for(Int_t i=0; i<16; i++){
    	myPadSetUp(myPad31->cd(i+1),0.00,0.00,0.00,0.15);
    	myPad31->cd(i+1);
	myHistSetUp(hTimeChannelSlice[i]);

    	if(i<8)hTimeChannelSlice[i]->GetXaxis()->SetRangeUser(62,68);
	if(i>7)hTimeChannelSlice[i]->GetXaxis()->SetRangeUser(52,62);
	hTimeChannelSlice[i]->Draw();
	
	latex.SetNDC();
	latex.SetTextFont(42);
	latex.SetTextSize(0.1);
	latex.DrawLatex(0.8,0.9,TString::Format("#mu = %2.1f",channelTimeMean[i]));
	latex.DrawLatex(0.8,0.8,TString::Format("#sigma = %0.2f",channelTimeSigma[i]));
	
	}
 
    c31->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    
    //--------------------------Time slewing---------------------------------    
    TH2D *hChannelSlice[16] = new TH2D(); 
    for(Int_t i=0; i<16; i++){
	TString channelName = "hChannelSlice";
	channelName += i;
	TString channelTitle = "Raw time Vs Charge, PM ";
	channelTitle += i;
	fHistTimeVsChargePerPM_UnCorr->GetZaxis()->SetRange(i+1,i+1);
	hChannelSlice[i] = dynamic_cast<TH2D*>fHistTimeVsChargePerPM_UnCorr->Project3D("xy");
	hChannelSlice[i]->SetName(channelName.Data());
	hChannelSlice[i]->SetTitle(channelTitle.Data());
	hChannelSlice[i]->GetXaxis()->SetTitle("Log10(1/charge) [ADC counts]");
	hChannelSlice[i]->GetYaxis()->SetTitle("Leading time [TDC counts]");
	}
	
    TSpline3 *fTimeSlewingSpline[16];
    AliCDBEntry *entTS = man->Get("AD/Calib/TimeSlewing");
    TList *fListSplines = (TList*)entTS->GetObject();
    for(Int_t i=0; i<16; i++)fTimeSlewingSpline[i] = dynamic_cast<TSpline3*> (fListSplines->At(i));

    TCanvas *c41 = new TCanvas("TimeSlewingADA"," ",1800,600);
    c41->Draw();
    c41->cd();
    TPad *myPad41 = new TPad("myPad41", "The pad",0,0,1,1);
    myPad41->Divide(4,2);
    myPad41->Draw();

    for(Int_t i=8; i<16; i++){
    	    myPadSetUp(myPad41->cd(i-7),0.15,0.1,0.1,0.15);
    	    myPad41->cd(i-7);
    	    gPad->SetLogz();
	    if(doRebin){hChannelSlice[i]->RebinX(2); hChannelSlice[i]->RebinY(2);}
	    hChannelSlice[i]->GetYaxis()->SetRangeUser(1600,2400);
    	    hChannelSlice[i]->DrawCopy("COLZ");
    	    fTimeSlewingSpline[i]->Draw("Psame");
    	    }
    c41->Print(Form("ADQA_Run_%d.pdf",runNumber)); 

    TCanvas *c42 = new TCanvas("TimeSlewingADC"," ",1800,600);
    c42->Draw();
    c42->cd();
    TPad *myPad42 = new TPad("myPad42", "The pad",0,0,1,1);
    myPad42->Divide(4,2);
    myPad42->Draw();

    for(Int_t i=0; i<8; i++){
    	    myPadSetUp(myPad42->cd(i+1),0.15,0.1,0.1,0.15);
    	    myPad42->cd(i+1);
    	    gPad->SetLogz();
	    if(doRebin){hChannelSlice[i]->RebinX(2); hChannelSlice[i]->RebinY(2);}
	    hChannelSlice[i]->GetYaxis()->SetRangeUser(1600,2400);
    	    hChannelSlice[i]->DrawCopy("COLZ");
    	    fTimeSlewingSpline[i]->Draw("Psame");
    	    }
    c42->Print(Form("ADQA_Run_%d.pdf",runNumber));
    //-------------------------------------------------------------------
    
    myHistSetUp(fHistTimeVsChargeADA_Cut);
    myHistSetUp(fHistTimeVsChargeADA_Corr);
    myHistSetUp(fHistTimeVsChargeADA_UnCorr);
    myHistSetUp(fHistTimeVsChargeADC_Cut);
    myHistSetUp(fHistTimeVsChargeADC_Corr);
    myHistSetUp(fHistTimeVsChargeADC_UnCorr);

    TCanvas *c5 = new TCanvas("TimeVsCharge"," ",1000,1000);
    c5->Draw();
    c5->cd();
    TPad *myPad5 = new TPad("myPad5", "The pad",0,0,1,1);
    myPad5->Divide(3,2);
    myPad5->Draw();
    myPadSetUp(myPad5->cd(1),0.15,0.1,0.1,0.15);
    myPadSetUp(myPad5->cd(2),0.15,0.1,0.1,0.15);
    myPadSetUp(myPad5->cd(3),0.15,0.1,0.1,0.15);
    myPadSetUp(myPad5->cd(4),0.15,0.1,0.1,0.15);
    myPadSetUp(myPad5->cd(5),0.15,0.1,0.1,0.15);
    myPadSetUp(myPad5->cd(6),0.15,0.1,0.1,0.15);
    myPad5->cd(1);
    gPad->SetLogz();
    if(doRebin){fHistTimeVsChargeADA_Cut->RebinX(2); fHistTimeVsChargeADA_Cut->RebinY(5);}
    fHistTimeVsChargeADA_Cut->GetXaxis()->SetRangeUser(40,80);
    fHistTimeVsChargeADA_Cut->Draw("COLZ");
    myPad5->cd(2);
    gPad->SetLogz();
    if(doRebin){fHistTimeVsChargeADA_Corr->RebinX(2); fHistTimeVsChargeADA_Corr->RebinY(5);}
    fHistTimeVsChargeADA_Corr->GetXaxis()->SetRangeUser(40,80);
    fHistTimeVsChargeADA_Corr->Draw("COLZ");
    myPad5->cd(3);
    gPad->SetLogz();
    if(doRebin){fHistTimeVsChargeADA_UnCorr->RebinX(2); fHistTimeVsChargeADA_UnCorr->RebinY(5);}
    fHistTimeVsChargeADA_UnCorr->GetXaxis()->SetRangeUser(170,220);
    fHistTimeVsChargeADA_UnCorr->Draw("COLZ");
    myPad5->cd(4);
    gPad->SetLogz();
    if(doRebin){fHistTimeVsChargeADC_Cut->RebinX(2); fHistTimeVsChargeADC_Cut->RebinY(5);}
    fHistTimeVsChargeADC_Cut->GetXaxis()->SetRangeUser(40,80);
    fHistTimeVsChargeADC_Cut->Draw("COLZ");
    myPad5->cd(5);
    gPad->SetLogz();
    if(doRebin){fHistTimeVsChargeADC_Corr->RebinX(2); fHistTimeVsChargeADC_Corr->RebinY(5);}
    fHistTimeVsChargeADC_Corr->GetXaxis()->SetRangeUser(40,80);
    fHistTimeVsChargeADC_Corr->Draw("COLZ");
    myPad5->cd(6);
    gPad->SetLogz();
    if(doRebin){fHistTimeVsChargeADC_UnCorr->RebinX(2); fHistTimeVsChargeADC_UnCorr->RebinY(5);}
    fHistTimeVsChargeADC_UnCorr->GetXaxis()->SetRangeUser(170,220);
    fHistTimeVsChargeADC_UnCorr->Draw("COLZ");

    c5->Print(Form("ADQA_Run_%d.pdf",runNumber));

    myHistSetUp(fHistNBBCoincidencesADAVsADC);
    fHistNBBCoincidencesADAVsADC->GetYaxis()->SetNdivisions(505);
    fHistNBBCoincidencesADAVsADC->GetXaxis()->SetNdivisions(505);
    myHistSetUp(fHistNBBCoincidencesADA);
    myHistSetUp(fHistNBBCoincidencesADC);
    fHistNBBCoincidencesADA->GetXaxis()->SetNdivisions(505);
    fHistNBBCoincidencesADA->GetXaxis()->SetTitle("Number of BB coincidences");
    fHistNBBCoincidencesADA->SetTitle("BB coincidences per side");
    fHistNBBCoincidencesADC->GetXaxis()->SetNdivisions(505);
    fHistNBBCoincidencesADA->SetLineColor(kBlue);
    fHistNBBCoincidencesADC->SetLineColor(kRed);
    fHistNBBCoincidencesADA->SetLineWidth(2);
    fHistNBBCoincidencesADC->SetLineWidth(2);

    TCanvas *c6 = new TCanvas("Coincidences"," ",1000,500);
    c6->Draw();
    c6->cd();
    TPad *myPad6 = new TPad("myPad6", "The pad",0,0,1,1);
    myPad6->Divide(2,1);
    myPad6->Draw();
    myPadSetUp(myPad6->cd(1),0.15,0.1,0.1,0.15);
    myPadSetUp(myPad6->cd(2),0.15,0.1,0.1,0.15);
    myPad6->cd(1);
    gPad->SetLogy();
    fHistNBBCoincidencesADA->Draw();
    fHistNBBCoincidencesADC->Draw("same");
    myLegend1->Draw();
    myPad6->cd(2);
    gPad->SetLogz();
    fHistNBBCoincidencesADAVsADC->Draw("COLZ");
    
    c6->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    TCanvas *c8 = new TCanvas("Noflag"," ",1000,500);
    c8->Draw();
    c8->cd();
    TPad *myPad8 = new TPad("myPad8", "The pad",0,0,1,1);
    myPad8->Divide(2,1);
    myPad8->Draw();
    myPadSetUp(myPad8->cd(1),0.15,0.1,0.1,0.15);
    myPadSetUp(myPad8->cd(2),0.15,0.1,0.1,0.15);
    myPad8->cd(1);
    gPad->SetLogz();
    if(doRebin)fHistChargeNoTime->RebinY(10);
    fHistChargeNoTime->Draw("COLZ");
    myPad8->cd(2);
    gPad->SetLogz();
    if(doRebin)fHistTimeNoFlag->RebinY(2);
    fHistTimeNoFlag->Draw("COLZ");

    c8->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    myHistSetUp(fHistMeanTimeADA);
    myHistSetUp(fHistMeanTimeADC);
    myHistSetUp(fHistMeanTimeDifference);
    
    TCanvas *c9 = new TCanvas("MeanTime"," ",1500,500);
    c9->Draw();
    c9->cd();
    TPad *myPad9 = new TPad("myPad9", "The pad",0,0,1,1);
    myPadSetUp(myPad9,0.15,0.1,0.1,0.15);
    myPad9->Divide(3,1);
    myPad9->Draw();
    myPadSetUp(myPad9->cd(1),0.15,0.15,0.15,0.15);
    myPadSetUp(myPad9->cd(2),0.15,0.15,0.15,0.15);
    myPadSetUp(myPad9->cd(3),0.15,0.15,0.15,0.15);
    myPad9->cd(1);
    //gPad->SetLogy();
    fHistMeanTimeADA->GetXaxis()->SetRangeUser(52,62);
    fHistMeanTimeADA->Draw();
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.65,0.75,TString::Format("#mu = %2.1f",meanTimeADA));
    latex.DrawLatex(0.65,0.65,TString::Format("#sigma = %0.2f",meanTimeSigmaADA));
    myPad9->cd(2);
    //gPad->SetLogy();
    fHistMeanTimeADC->GetXaxis()->SetRangeUser(62,68);
    fHistMeanTimeADC->Draw();
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.65,0.75,TString::Format("#mu = %2.1f",meanTimeADC));
    latex.DrawLatex(0.65,0.65,TString::Format("#sigma = %0.2f",meanTimeSigmaADC));
    myPad9->cd(3);
    //gPad->SetLogy();
    fHistMeanTimeDifference->GetXaxis()->SetRangeUser(-16,0);
    fHistMeanTimeDifference->Draw();

    c9->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    myHistSetUp(fHistMeanTimeCorrelation);
    myHistSetUp(fHistMeanTimeSumDiff);
    
    TCanvas *c10 = new TCanvas("MeanTimeCorr"," ",1000,500);
    c10->Draw();
    c10->cd();
    TPad *myPad10 = new TPad("myPad10", "The pad",0,0,1,1);
    myPadSetUp(myPad10,0.15,0.1,0.1,0.15);
    myPad10->Divide(2,1);
    myPad10->Draw();
    myPadSetUp(myPad10->cd(1),0.15,0.15,0.15,0.15);
    myPadSetUp(myPad10->cd(2),0.15,0.15,0.15,0.15);
    myPad10->cd(1);
    gPad->SetLogz();
    if(doRebin){fHistMeanTimeCorrelation->RebinX(2);fHistMeanTimeCorrelation->RebinY(2);}
    fHistMeanTimeCorrelation->Draw("COLZ");
    myPad10->cd(2);
    gPad->SetLogz();
    if(doRebin){fHistMeanTimeSumDiff->RebinX(2);fHistMeanTimeSumDiff->RebinY(2);} 
    fHistMeanTimeSumDiff->Draw("COLZ");

    c10->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    TCanvas *c11 = new TCanvas("TriggerInputs"," ",1000,500);
    c11->Draw();
    c11->cd();
    TPad *myPad11 = new TPad("myPad11", "The pad",0,0,1,1);
    myPadSetUp(myPad11,0.15,0.15,0.15,0.15);
    myPad11->Draw();
    gPad->SetLogy();
    fHistTriggerMasked->Draw("HIST");
    fHistTriggerMasked->Draw("TEXT0SAME");
    
    c11->Print(Form("ADQA_Run_%d.pdf",runNumber));
    
    TCanvas *c12 = new TCanvas("Decisions"," ",1000,500);
    c12->Draw();
    c12->cd();
    TPad *myPad12 = new TPad("myPad12", "The pad",0,0,1,1);
    myPadSetUp(myPad12,0.15,0.15,0.15,0.15);
    myPad12->Draw();
    gPad->SetLogz();
    fHistDecision->Draw("COLZTEXT");

    c12->Print(Form("ADQA_Run_%d.pdf)",runNumber));
    
    if(runNumber > 251861){
        AliCDBEntry *entTrends = 0x0;
        entTrends = man->Get("AD/Calib/PMTrends");
	
    	TClonesArray *fGraphs = (TClonesArray*)entTrends->GetObject();
    	TGraph *fGraphPM[32];
  
    	TCanvas *cHV = new TCanvas("AD0_HV","AD0_HV",1000,500);
    	cHV->Divide(4,4);
   
    	TCanvas *cImon = new TCanvas("AD0_Imon","AD0_Imon",1000,500);
    	cImon->Divide(4,4);
  
    	for(int iPM = 0; iPM<16 ; iPM++){
  		myPadSetUp(cHV->cd(iPM+1),0.2,0.10,0.02,0.2);
  		cHV->cd(iPM+1);

		fGraphPM[iPM] = (TGraph*)fGraphs->At(iPM);
	
		fGraphPM[iPM]->SetMarkerStyle(20);
		fGraphPM[iPM]->SetMarkerSize(0.5);
		fGraphPM[iPM]->SetTitle(Form("PM%d",iPM));
		fGraphPM[iPM]->GetYaxis()->SetTitle("Voltage [V]");
		fGraphPM[iPM]->GetXaxis()->SetTitle("DayTime");
		fGraphPM[iPM]->GetYaxis()->SetTitleSize(0.08);
		fGraphPM[iPM]->GetXaxis()->SetTitleSize(0.08);
		fGraphPM[iPM]->GetYaxis()->SetLabelSize(0.08);
		fGraphPM[iPM]->GetXaxis()->SetLabelSize(0.08);
		fGraphPM[iPM]->GetXaxis()->SetNdivisions(502);
		fGraphPM[iPM]->GetXaxis()->SetTimeFormat("%H:%M");
		fGraphPM[iPM]->Draw("APL");
  		}
		cHV->Print(Form("HVtrends_Run_%d.pdf(",runNumber));
	
    	for(int iPM = 16; iPM<32 ; iPM++){
  		cImon->cd(iPM-15);
		myPadSetUp(cImon->cd(iPM-15),0.2,0.10,0.02,0.2);
		fGraphPM[iPM] = (TGraph*)fGraphs->At(iPM);
		fGraphPM[iPM]->SetMarkerStyle(20);
		fGraphPM[iPM]->SetMarkerSize(0.5);
		fGraphPM[iPM]->SetTitle(Form("PM%d",iPM-16));
		fGraphPM[iPM]->GetYaxis()->SetTitle("Current [mA]");
		fGraphPM[iPM]->GetXaxis()->SetTitle("DayTime");
		fGraphPM[iPM]->GetYaxis()->SetTitleSize(0.08);
		fGraphPM[iPM]->GetXaxis()->SetTitleSize(0.08);
		fGraphPM[iPM]->GetYaxis()->SetLabelSize(0.08);
		fGraphPM[iPM]->GetXaxis()->SetLabelSize(0.08);
		fGraphPM[iPM]->GetXaxis()->SetNdivisions(502);
		fGraphPM[iPM]->GetXaxis()->SetTimeFormat("%H:%M");
		fGraphPM[iPM]->Draw("APL");		
  		}
	cImon->Print(Form("HVtrends_Run_%d.pdf)",runNumber));
	}

        
    
    }
  delete adQAdir;
  
 
  return 0;
}

bool IsADReady(Int_t run)
{

  if(run<200000) return false;
  bool result=true;
  
  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://");
  man->SetRun(run);

  AliCDBEntry *ent = man->Get("AD/Calib/Data");
  AliADCalibData *calData = (AliADCalibData*)ent->GetObject();

  int HV=calData->GetMeanHV(0);
  if(HV<100)
    result=false;
    
  return result;
}

void myScaleSetUp(TH1* histoBlue, TH1* histoRed){

  Float_t min[2], max[2];

  max[0] = histoBlue->GetBinContent(histoBlue->GetMaximumBin());
  max[1] = histoRed->GetBinContent(histoRed->GetMaximumBin());

  min[0] = histoBlue->GetBinContent(histoBlue->GetMinimumBin());
  min[1] = histoRed->GetBinContent(histoRed->GetMinimumBin());
	
  histoBlue->GetYaxis()->SetRangeUser(TMath::MinElement(2,min)+0.9,1000*TMath::MaxElement(2,max));  

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
  currentGraph->GetYaxis()->SetTitleOffset(1.1);  
  currentGraph->SetTitleSize(0.06,"xyz");
  currentGraph->SetStats(kFALSE); 
  return;
}
void myHistSetUp(TH2 *currentGraph=0){
 
  currentGraph->SetLabelSize(0.05,"xyz");
  currentGraph->SetLabelFont(42,"xyz"); 
  currentGraph->SetLabelOffset(0.01,"xyz");
  currentGraph->SetTitleFont(42,"xyz"); 
  currentGraph->GetXaxis()->SetTitleOffset(1.3);
  currentGraph->GetYaxis()->SetTitleOffset(1.3);  
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
