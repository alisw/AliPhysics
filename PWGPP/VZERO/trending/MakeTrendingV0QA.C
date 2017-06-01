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
#include "AliVZEROCalibData.h"
#endif

bool HigtVoltage(Int_t run);

Int_t MakeTrendingV0QA(TString qafilename,Int_t runNumber,TString ocdbStorage = "raw://",Bool_t IsOnGrid = kFALSE,Bool_t canvasE = kFALSE)
{
  TTree *ttree=new TTree("trending","tree of trending variables");
  Int_t higtVoltage=0,invalidInput=0,qaNotFound=0,v0active=0,v0qaNotfound=0,noEntries=0;
  Int_t NumberVoieOff=0, numberBadOffset=0;
  Float_t TimesA=-9999.,TimesC=-9999., BB_BG=-9999.,BB_EE=-9999.,AdcA=-9999.;
  Float_t AdcC=-9999.,MultA=-9999.,MultC=-9999.,ChargeCh46=-9999.,ChargeAllNo46=-9999.;
  Float_t AdcAError=-9999.,AdcCError=-9999.;
  Float_t TriggerEff_CVLN=-9999.,TriggerEff_CVHN=-9999.,TriggerEff_CVHN2=-9999.;
  Float_t TriggerEff_CVLN_Error=-9999.,TriggerEff_CVHN_Error=-9999.,TriggerEff_CVHN2_Error=-9999.;
  Float_t PMTEdges[64]={-9999.};
  Float_t PMTEdgesError[64]={-9999.};
  Float_t PMmeanAdc[64]={-9999.};
  Float_t RingmeanAdc[8]={-9999.};
  Float_t VOXmeanAdc[2]={-9999.};
  Float_t VOmeanAdc=-9999.;
  Bool_t isPP;
  TString treePostFileName="trending.root";

  ttree->Branch("HigtVoltage",&higtVoltage,"higtVoltage/I");
  ttree->Branch("invalidInput",&invalidInput,"invalidInput/I");
  ttree->Branch("qaNotFound",&qaNotFound,"qaNotFound/I");
  ttree->Branch("v0active",&v0active,"v0active/I");
  ttree->Branch("v0qaNotfound",&v0qaNotfound,"v0qaNotfound/I");
  ttree->Branch("noEntries",&noEntries,"noEntries/I");
  ttree->Branch("NumberVoieOff",&NumberVoieOff,"Number of path off/I");
  ttree->Branch("numberBadOffset",&numberBadOffset,"Number of bad offset /I");
  ttree->Branch("run",&runNumber,"run/I");
  if (!qafilename) 
    {
      if(!HigtVoltage(runNumber))
	higtVoltage=1;
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
  
  if(IsOnGrid)
    TGrid::Connect("alien://");
  TFile*fin=TFile::Open(qafilename,"r");

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
  if(beamType!="A-A")
    isPP=kTRUE;
  else
    isPP=kFALSE;
  if(beamType!="A-A")
    {
      ttree->Branch("isPP",&isPP,"isPP/O");
      ttree->Branch("TimesA",&TimesA,"BB Leading time;;Time (ns)/F");
      ttree->Branch("TimesC",&TimesC,"BB Leading time;;Time (ns)/F");
      ttree->Branch("BB_BG",&BB_BG,"Trigger ratio/F");
      ttree->Branch("BB_EE",&BB_EE,"Trigger ratio/F");
      ttree->Branch("AdcA" ,&AdcA ,"Average Charge/F");
      ttree->Branch("AdcC" ,&AdcC ,"Average Charge/F");
      ttree->Branch("MultA",&MultA,"Average number of Fired cell/F");
      ttree->Branch("MultC",&MultC,"Average number of Fired cell/F");
      ttree->Branch("ChargeCh46",&ChargeCh46,"Integreted charge of chanel 46/F");
      ttree->Branch("ChargeAllNo46",&ChargeAllNo46,"Integreted charge of all chanel without 46/F");
      for(int i = 0; i < 64; i++)
	  ttree->Branch(Form("PMmeanAdc[%d]",i),&PMmeanAdc[i],Form("Mean ADC Cell %d/F",i));
      for(int i=0;i<8;i++)
	ttree->Branch(Form("RingmeanAdc[%d]",i),&RingmeanAdc[i],Form("Mean ADC Ring %d/F",i));
      ttree->Branch("VOXmeanAdc[0]",&VOXmeanAdc[0],"Mean ADC V0C/F");
      ttree->Branch("VOXmeanAdc[1]",&VOXmeanAdc[1],"Mean ADC V0A/F");
      ttree->Branch("VOmeanAdc",&VOmeanAdc,"Mean ADC V0/F");

      /*
  //PMmeanAdc[64] RingmeanAdc[8] VOXmeanAdc[2] VOmeanAdc
      */

      // ChargeCh46 ChargeAllNo46
    }
  else
    {
      ttree->Branch("isPP",&isPP,"isPP/O");
      ttree->Branch("TimesA",&TimesA,"BB Leading time;;Time (ns)/F");
      ttree->Branch("TimesC",&TimesC,"BB Leading time;;Time (ns)/F");
      ttree->Branch("BB_BG",&BB_BG,"Trigger ratio/F");
      ttree->Branch("BB_EE",&BB_EE,"Trigger ratio/F");
      ttree->Branch("AdcA" ,&AdcA ,"Average Charge/F");
      ttree->Branch("AdcAError",&AdcAError,"AdcAError/F");
      ttree->Branch("AdcC" ,&AdcC ,"Average Charge/F");
      ttree->Branch("AdcCError",&AdcCError,"AdcCError/F");
      ttree->Branch("MultA",&MultA,"Average number of Fired cell/F");
      ttree->Branch("MultC",&MultC,"Average number of Fired cell/F");
      ttree->Branch("TriggerEff_CVLN",&TriggerEff_CVLN,"CVLN CVBN/F");
      ttree->Branch("TriggerEff_CVHN",&TriggerEff_CVHN,"CVHN CVBN/F");
      ttree->Branch("TriggerEff_CVHN2",&TriggerEff_CVHN2,"CVHN CVLN/F");
      ttree->Branch("TriggerEff_CVLN_Error",&TriggerEff_CVLN_Error,"CVLN CVBN/F");
      ttree->Branch("TriggerEff_CVHN_Error",&TriggerEff_CVHN_Error,"CVHN CVBN/F");
      ttree->Branch("TriggerEff_CVHN2_Error",&TriggerEff_CVHN2_Error,"CVHN CVLN/F");
      for(int i = 0; i < 64; ++i)
	{
	  ttree->Branch(Form("PMTEdges[%d]",i),&PMTEdges[i],Form("Multiplicity edge Cell %d/F",i));
	  ttree->Branch(Form("PMTEdgesError[%d]",i),&PMTEdgesError[i],Form("Multiplicity edge Cell %d/F",i));
	}
    }

  if(!fin)
    {
      if(!HigtVoltage(runNumber))
	higtVoltage=1;
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
      Printf("INFO: QA output file %s open. \n",fin->GetName());
    }

  printf("activeDetList %s\nrunType %s\nbeamType %s\nmachineMode %s\nlhcState %s\n",
	 activeDetList.Data(),runType.Data(),beamType.Data(),
	 machineMode.Data(),lhcState.Data());
  
  time_t duration = fGRPData->GetTimeEnd() - fGRPData->GetTimeStart();



  if(!activeDetList.Contains("VZERO"))
    { 
      if(!HigtVoltage(runNumber))
	higtVoltage=1;
      v0active=1;
      printf("RUN WITH VZERO NOT ACTIVE\n");
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
  if(beamType!="A-A")
    {
      if(duration<120)
	{ 
	  printf("RUNS SHORTER THAN 2 MIN\n");
	  return 0;
	}
      

      char v0QAdirName[30]="VZERO_Performance";
      fin->Print();
      printf("\n\n");
      fin->ls();
      TDirectoryFile * v0QAdir=(TDirectoryFile*)fin->Get(v0QAdirName);
      if(!v0QAdir)
	{
	  if(!HigtVoltage(runNumber))
	    higtVoltage=1;
	  else
	    v0qaNotfound=1;
	  printf("ERROR: VZERO QA directory not present in input file.\n");
	  TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
	  ttree->Fill();
	  trendFile->cd();
	  ttree->Write();
	  trendFile->Close();
      
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
      TH2F *hadcpmtnotime = (TH2F*)list->FindObject("hadcpmtnotime");
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
	{
	  if(!HigtVoltage(runNumber))
	    higtVoltage=1;
	  noEntries=1;
	  TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
	  ttree->Fill();
	  trendFile->cd();
	  ttree->Write();
	  trendFile->Close();
      
	  printf("hAdcWithTimeA no entries\n");
	  return 0;
	}
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
      
      char v0QAdirNameTrig[30]="VZERO_Performance_Trig";
      TDirectoryFile * v0QAdirTrig=(TDirectoryFile*)fin->Get(v0QAdirNameTrig);
      TH2F *hadcpmtwithtimeTrig ;
      if(v0QAdirTrig)
	{
	  TList *listTrig = (TList*)v0QAdirTrig->Get("QAVZEROHistsTrig");
	  hadcpmtwithtimeTrig = (TH2F*)listTrig->FindObject("hadcpmtwithtime");
	}
      else
	{
	  hadcpmtwithtimeTrig=new TH2F("hadcpmtwithtimeTrig","hadcpmtwithtimeTrig",64,-0.5,63.5,1,-10000,-9998);
	  for(Int_t i=0;i<64;i++)
	    {
	      hadcpmtwithtimeTrig->SetBinContent(i+1,1,-9999.);
	      hadcpmtwithtimeTrig->SetBinError(i+1,1,0);
	    }
	}
      double valBin=0;
      double Max=hadcpmtwithtime->GetYaxis()->GetLast(), Min=1;
      double mean=Max-Min,meanRing=0,meanV0X=0,meanV0=0;
      TH1D*hadcXFull=hadcpmtwithtime->ProjectionX("hadcXFull",1,Max);
      TH1D*hadcX=hadcpmtwithtime->ProjectionX("hadcX",5,15);
      ChargeAllNo46=0;
      for(Int_t i=0;i<64;i++)
	{
	  TH1D*MeanValue=(TH1D*)hadcpmtwithtimeTrig->ProjectionY("MeanValue",i+1,i+1);
	  valBin=hadcX->GetBinContent(i+1);
	  if(valBin==0)
	    numberBadOffset++;
	  valBin=hadcXFull->GetBinContent(i+1);
	  if(valBin==0)
	    NumberVoieOff++;
	  if(i==46)
	    ChargeCh46=valBin;
	  else
	    ChargeAllNo46+=valBin;

	  PMmeanAdc[i]=MeanValue->GetMean();
	  printf("PMmeanAdc[%d]: %.3f\n",i,PMmeanAdc[i]);
	  meanRing+=PMmeanAdc[i];
	  meanV0X +=PMmeanAdc[i];
	  meanV0  +=PMmeanAdc[i];
	  if(!((i+1)%8))
	    {
	      RingmeanAdc[i/8]=meanRing;
	      meanRing=0;
	    }
	  if(!((i+1)%32))
	    {
	      VOXmeanAdc[i/32]=meanV0X;
	      meanV0X=0;
	    }
	}
      VOmeanAdc=meanV0;
      if(NumberVoieOff>=63)
	ChargeAllNo46=-9999.;
      else
	ChargeAllNo46=ChargeAllNo46/(63.-NumberVoieOff);
      TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
      ttree->Fill();
      trendFile->cd();
      ttree->Write();
      trendFile->Close();
      TFile * chargeFile = new TFile("chargeADC.root","recreate");
      chargeFile->cd();
      Double_t val=0;
      for(Int_t channel=0;channel<=65;channel++)
	{
	  for(Int_t adc=0;adc<=201;adc++)
	    {
	      val=0;
	      val=hadcpmtnotime->GetBinContent(channel,adc);
	      val+=hadcpmtwithtime->GetBinContent(channel,adc);
	      hadcpmtnotime->SetBinContent(channel,adc,val);
	    }
	}
      hadcpmtnotime->Write();
      chargeFile->Close();

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
	  

	  cOut->Clear("D");
	  cOut->cd(1);
	  cOut->GetPad(1)->SetLogy(0);//cOut->GetPad(1)->SetLogz();
	  hAdcTimeA->Draw("colz");

	  
	  cOut->cd(2);
	  cOut->GetPad(2)->SetLogy(0);//cOut->GetPad(2)->SetLogz();
	  hAdcTimeC->Draw("colz");

	  cOut->cd(3);
	  cOut->GetPad(3)->SetLogy();// cOut->GetPad(3)->SetLogz(0);
	  hV0ampl->Draw();

	  cOut->cd(4);
	  cOut->GetPad(4)->SetLogy(0);// cOut->GetPad(4)->SetLogz(0);
	  htimecorr->Draw("colz");


	  printf("%p\n",cOut);
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
    }
  else
    {
      printf("Pb-Pb run \n\n\n");
      if(duration<600)
	printf("RUNS SHORTER THAN 10 MIN\n");
      
      char v0QAdirName[30]="VZERO_PbPb_Performance";
      TDirectoryFile * v0QAdir=(TDirectoryFile*)fin->Get(v0QAdirName);
      if(!v0QAdir)
	{
	  printf("ERROR: VZERO QA directory not present in input file.\n");
	  return -1;
	}
  
      TList *list = (TList*)v0QAdir->Get("PbPbVZEROHists");
      if(!list) 
	{
	  cout << "ERROR: No list found" << endl;
	  return -1;
	}
      
 
      TH1F *hL2Triggers=(TH1F*)list->FindObject("hL2Triggers");
      TH2F *hTriggerDecision = (TH2F*)list->FindObject("hTriggerDecision");
      TH1F *hAdcNoTimeA = (TH1F*)list->FindObject("hAdcNoTimeV0A");
      TH1F *hAdcWithTimeA = (TH1F*)list->FindObject("hAdcWithTimeV0A");
      TH1F *hAdcNoTimeC = (TH1F*)list->FindObject("hAdcNoTimeV0C");
      TH1F *hAdcWithTimeC = (TH1F*)list->FindObject("hAdcWithTimeV0C");
      TH2F *hadcpmtwithtime = (TH2F*)list->FindObject("hadcpmtwithtime");	
      TH1F *htimepmtA = (TH1F*)list->FindObject("htimepmtV0A");
      TH1F *htimepmtC = (TH1F*)list->FindObject("htimepmtV0C");
      TH1F *hwidthA = (TH1F*)list->FindObject("hwidthV0A");
      TH1F *hwidthC = (TH1F*)list->FindObject("hwidthV0C");
      //	TH1F *hV0ampl = (TH1F*)list->FindObject("hV0ampl");
      TH2F *htimepmt = (TH2F*)list->FindObject("htimepmt");	
      TH2F *hwidthpmt = (TH2F*)list->FindObject("hwidthpmt");	
      TH2F *hadcwidthA = (TH2F*)list->FindObject("hadcwidthV0A");	
      TH2F *hadcwidthC = (TH2F*)list->FindObject("hadcwidthV0C");	
      TH2F *hAdcTimeA = (TH2F*)list->FindObject("hAdcTimeV0A");	
      TH2F *hAdcTimeC = (TH2F*)list->FindObject("hAdcTimeV0C");	
      TH2F *htimecorr = (TH2F*)list->FindObject("htimecorr");	
      TH2F *hNFlags   = (TH2F*)list->FindObject("hNFlags");
      TH1F *hV0A = (TH1F*)hNFlags->ProjectionX("hV0A",1,hNFlags->GetNbinsY());
      TH1F *hV0C = (TH1F*)hNFlags->ProjectionY("hV0C",1,hNFlags->GetNbinsX());
      TH2F* hVtxXYBB  =(TH2F*) list->FindObject("fhVtxXYBB");
      TH1F* hVtxZBB   =(TH1F*) list->FindObject("fhVtxZBB");
      TH2F* hVtxXYBGA =(TH2F*) list->FindObject("fhVtxXYBGA");
      TH1F* hVtxZBGA  =(TH1F*) list->FindObject("fhVtxZBGA");
      TH2F* hVtxXYBGC =(TH2F*) list->FindObject("fhVtxXYBGC");
      TH1F* hVtxZBGC  =(TH1F*) list->FindObject("fhVtxZBGC");
  
      Int_t last=hL2Triggers->GetXaxis()->GetLast();
      if (last==0)
	printf("Arg, hL2Triggers n'ai pas cool");
      Bool_t CPBI2_B1=kFALSE, CSEMI_R1=kFALSE, CCENT_R2=kFALSE, CVLN_R1=kFALSE, CVHN_R2=kFALSE;
      Bool_t CVLN_B2 =kFALSE, CPBI1   =kFALSE, CVLN    =kFALSE, CVHN   =kFALSE, CINT7=kFALSE;
      for(Int_t i=0;i<last;i++)
	{
	  TString trigerName=hL2Triggers->GetXaxis()->GetBinLabel(i);
	  if(trigerName.Contains("CINT7"))
	    CINT7=kTRUE;
	  if(trigerName.Contains("CPBI2_B1"))
	    CPBI2_B1=kTRUE;
	  if(trigerName.Contains("CSEMI_R1"))
	    CSEMI_R1=kTRUE;
	  if(trigerName.Contains("CCENT_R2"))
	    CCENT_R2=kTRUE;
	  if(trigerName.Contains("CVLN_R1"))
	    CVLN_R1=kTRUE;
	  if(trigerName.Contains("CVHN_R2"))
	    CVHN_R2=kTRUE;
	  if(trigerName.Contains("CVLN_B2"))
	    CVLN_B2=kTRUE;
	  if(trigerName.Contains("CPBI1"))
	    CPBI1=kTRUE;
	  if(trigerName.Contains("CVLN"))
	    CVLN=kTRUE;
	  if(trigerName.Contains("CVHN"))
	    CVHN=kTRUE;
	}
      TString trigMB,trigCVLN,trigCVHN;
      if("CINT7")
	trigMB = "CINT7";
      else
	{
	  if(CPBI2_B1)
	    trigMB = "CPBI2_B1";
	  else
	    trigMB   = "CPBI1";
	  if(CSEMI_R1)
	    trigCVLN = "CSEMI_R1";
	  else if(CVLN_R1)
	    trigCVLN = "CVLN_R1";
	  else if(CVLN_B2)
	    trigCVLN = "CVLN_B2";
	  else
	    trigCVLN = "CVLN";
	  if(CCENT_R2)
	    trigCVHN = "CCENT_R2";
	  else if(CVHN_R2)
	    trigCVHN = "CVHN_R2";
	  else
	    trigCVHN = "CVHN";
	}
      printf("Relou: %s\t%s\t%s\t%s\n",trigMB.Data(), trigCVLN.Data(),trigCVHN.Data());
      TH2F* hRecoMult = (TH2F*) list->FindObject(Form("hRecoMult_%s-",trigMB.Data()));
      TH2F* hRecoMultPMT = (TH2F*) list->FindObject(Form("hRecoMultPMT_%s-",trigMB.Data()));
      TH1F* hTotRecoMult = (TH1F*) list->FindObject(Form("hTotRecoMult_%s-",trigMB.Data()));
      TH1F* hTotRecoMult_CVLN = (TH1F*) list->FindObject(Form("hTotRecoMult_%s-",trigCVLN.Data()));
      TH1F* hTotRecoMult_CVHN = (TH1F*) list->FindObject(Form("hTotRecoMult_%s-",trigCVHN.Data()));
      TH2F* hEqualizedMult = (TH2F*) list->FindObject(Form("hEqualizedMult_%s-",trigMB.Data()));

      Double_t BB  = hTriggerDecision->GetBinContent(2,2);
      Double_t EE  = hTriggerDecision->GetBinContent(1,1);
      Double_t BGA = hTriggerDecision->GetBinContent(3,2);
      Double_t BGC = hTriggerDecision->GetBinContent(2,3);
 
      if(hAdcWithTimeA->GetEntries()==0)
	{
	  TFile * trendFile = new TFile(treePostFileName.Data(),"recreate");
	  ttree->Fill();
	  trendFile->cd();
	  ttree->Write();
	  trendFile->Close();
      
	  printf("hAdcWithTimeA no entries\n");
	  return 0;
	}


      Double_t cVLN =0, cVHN = 0,cVBN=0;
      if (hTotRecoMult_CVLN)
	cVLN = hTotRecoMult_CVLN->GetEntries();
      if (hTotRecoMult_CVHN)
        cVHN = hTotRecoMult_CVHN->GetEntries();
      if(hTotRecoMult)
        cVBN = hTotRecoMult->GetEntries();

      man->SetRun(runNumber);
      AliCDBEntry *entryCTP = man->Get("GRP/CTP/Config");
      AliTriggerConfiguration *configCTP = (AliTriggerConfiguration*)entryCTP->GetObject();
      TObjArray  inputsArray = configCTP->GetInputs();

      Double_t rnd1=1., rnd2=1., bc1=1., bc2=1.;
      AliTriggerInput * input;
      input = (AliTriggerInput*)(inputsArray.FindObject("RND1"));
      if(input) rnd1 =  (input->GetSignature())/(double)(0x7fffffff );
  
      input = (AliTriggerInput*)(inputsArray.FindObject("RND2"));
      if(input) rnd2 =  (input->GetSignature())/(double)(0x7fffffff );
  
      input = (AliTriggerInput*)(inputsArray.FindObject("BC1"));
      if(input) bc1 =  1./(input->GetSignature()+1.);
  
      input = (AliTriggerInput*)(inputsArray.FindObject("BC2"));
      if(input) bc2 =  1./(input->GetSignature()+1.);

      Double_t scaleVBN = 1., scaleVLN = 1., scaleVHN = 1.;
      if(trigMB.Contains("B1")) 
	{
	  if(bc1>0.) 
	    scaleVBN=1./bc1;
	} else if(trigMB.Contains("B2"))
	{
	  if(bc2>0.)
	    scaleVBN=1./bc2;
	} else if(trigMB.Contains("R1"))
	{
	  if(rnd1>0.)
	    scaleVBN=1./rnd1;
 	} else if(trigMB.Contains("R2"))
	{
	  if(rnd2>0.) 
	    scaleVBN=1./rnd2;
	}
  
      if(trigCVLN.Contains("B1")) 
	{
	  if(bc1>0.) 
	    scaleVLN=1./bc1;
	} else if(trigCVLN.Contains("B2"))
	{
	  if(bc2>0.)
	    scaleVLN=1./bc2;
	} else if(trigCVLN.Contains("R1")) 
	{
	  if(rnd1>0.) 
	    scaleVLN=1./rnd1;
	} else if(trigCVLN.Contains("R2")) 
	{
	  if(rnd2>0.) 
	    scaleVLN=1./rnd2;
	}
    
      if(trigCVHN.Contains("B1")) 
	{
	  if(bc1>0.)
	    scaleVHN=1./bc1;
	} else if(trigCVHN.Contains("B2"))
	{
	  if(bc2>0.)
	    scaleVHN=1./bc2;
	} else if(trigCVHN.Contains("R1"))
	{
	  if(rnd1>0.)
	    scaleVHN=1./rnd1;
	} else if(trigCVHN.Contains("R2")) 
	{
	  if(rnd2>0.)
	    scaleVHN=1./rnd2;
	}
  

      if (hTotRecoMult_CVLN)
	hTotRecoMult_CVLN->Scale(scaleVLN);
      if (hTotRecoMult_CVHN)
	hTotRecoMult_CVHN->Scale(scaleVHN);
      if(hTotRecoMult)
	hTotRecoMult->Scale(scaleVBN);
      cVBN *= scaleVBN/100.;
      cVLN *= scaleVLN;
      cVHN *= scaleVHN;
      if(cVBN >0.)
	{   
	  if(cVLN>0.)
	    TriggerEff_CVLN=cVBN/cVLN;
	  if((cVLN>0.)&&(cVHN>0.)) 
	    TriggerEff_CVLN_Error=cVBN/cVLN*(TMath::Sqrt(1./cVLN+1./cVBN));
	  if(cVHN>0.)
	    TriggerEff_CVHN=cVBN/cVHN;
	  if(cVHN>0.)
	    TriggerEff_CVHN_Error=cVBN/cVHN*(TMath::Sqrt(1./(cVHN/scaleVHN)+1./(cVBN/scaleVBN*100.)));
	  if((cVLN>0.)&&(cVHN>0.))
	    TriggerEff_CVHN2=cVLN/cVHN;
	  if(cVHN>0. && cVLN>0.) 
	    TriggerEff_CVHN2_Error=cVLN/cVHN*(TMath::Sqrt(1./(cVHN/scaleVHN)+1./(cVLN/scaleVLN)));
	}

      Double_t beta1[64], beta2[64];
      Double_t q = 1. - 1.e-4;
      Double_t q2 = 1. - 2.e-4;
      if(hRecoMultPMT)
	{
	  if(hRecoMultPMT->GetEntries()!=0)
	    {
	      for(int i = 0; i < 64; ++i)
		{
		  ((TH1D*)hRecoMultPMT->ProjectionY(Form("hRecoMultPMT%d",i),i+1,i+1))->GetQuantiles(1,&beta1[i],&q);
		  ((TH1D*)hRecoMultPMT->ProjectionY(Form("hRecoMultPMT%d",i),i+1,i+1))->GetQuantiles(1,&beta2[i],&q2);
		  PMTEdges[i]=(beta1[i]+beta2[i])/2.;
		  PMTEdgesError[i]=beta1[i] - beta2[i];
		  printf("%.3f\t%.3f\n", PMTEdges[i],PMTEdgesError[i]);
		}
	    }
	}
      Double_t betaSide1[2];
      Double_t betaSide2[2];
      if(hRecoMultPMT)
	{
	  if(hRecoMult->GetEntries()!=0)
	    {
	      for(int i = 0; i < 2; ++i)
		{
		  if(i==0)
		    {
		      ((TH1D*)hRecoMult->ProjectionY(Form("hRecoMult1%d",i)))->GetQuantiles(1,&betaSide1[i],&q);
		      ((TH1D*)hRecoMult->ProjectionY(Form("hRecoMult2%d",i)))->GetQuantiles(1,&betaSide2[i],&q2);
		    }
		  else
		    {
		      ((TH1D*)hRecoMult->ProjectionX(Form("hRecoMult1%d",i)))->GetQuantiles(1,&betaSide1[i],&q);
		      ((TH1D*)hRecoMult->ProjectionX(Form("hRecoMult2%d",i)))->GetQuantiles(1,&betaSide2[i],&q2);
		    }
		}
	    }
	  printf("hRecoMult empty\n");
	}
      AdcA=(betaSide1[1] + betaSide2[1])/2.;
      AdcC=(betaSide1[0] + betaSide2[0])/2.;
      AdcAError=betaSide1[1] - betaSide2[1];
      AdcCError=betaSide1[0] - betaSide2[0];
      printf("AdcA: %.3f\tAdcC: %.3f\n",AdcA,AdcC);
      TSpectrum s;
      int nPeaksFound;
      float * peaks;
      float max;
      float shiftA = 8.;
  
      nPeaksFound = s.Search(htimepmtA); peaks = s.GetPositionX(); max = -25.;
      for(int i=0;i<nPeaksFound;i++)
	{
	  if(peaks[i]>max)
	    max = peaks[i];
	}
      htimepmtA->Fit("gaus","","",max-1.,max+1.);
      TimesA=htimepmtA->GetFunction("gaus")->GetParameter(1)-shiftA;
  
      nPeaksFound = s.Search(htimepmtC); peaks = s.GetPositionX(); max = -25.;
      for(int i=0;i<nPeaksFound;i++) 
	{
	  if(peaks[i]>max) 
	    max = peaks[i];
	}	
      htimepmtC->Fit("gaus","","",max-1.,max+1.);
      TimesC=htimepmtC->GetFunction("gaus")->GetParameter(1);
  
      if(BB) 
	{
	  BB_BG=((BGA+BGC)/BB);
	  BB_EE=(EE/BB);
	} 
      else
	{
	  BB_BG=(0);
	  BB_EE=(0);
	}
      //printf("%f\t%f\n\n",TimesA,TimesC);
      MultA=(hV0A->GetMean());
      MultC=(hV0C->GetMean());

      double valBin=-9999;
      TH1D*hadcXFull=hadcpmtwithtime->ProjectionX("hadcXFull",1,hadcpmtwithtime->GetYaxis()->GetLast());
      TH1D*hadcX=hadcpmtwithtime->ProjectionX("hadcX",5,15);
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
	


      //-------------
      cout << "DEBUG: Processing finished. Now painting..." << endl;
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
      
	  if(hEqualizedMult)
	    {
	      cOut->cd(4); cOut->GetPad(4)->SetLogy(0);cOut->GetPad(4)->SetLogz(1);
	      hEqualizedMult->Draw("colz");
	    }
	  cOut->Print(Form("QA_Run_%d.pdf(",runNumber));
	  //-------------
      
	  cOut->Clear();
	  cOut->Divide(2,2);
	  cOut->cd(1); cOut->GetPad(1)->SetLogy();
	  htimepmtA->GetXaxis()->SetRangeUser(-25.,25.); htimepmtA->Draw();
      
	  cOut->cd(2); cOut->GetPad(2)->SetLogy();
	  htimepmtC->GetXaxis()->SetRangeUser(-25.,25.); htimepmtC->Draw();
      
	  cOut->cd(3); cOut->GetPad(3)->SetLogy();cOut->GetPad(3)->SetLogz(0);
	  //hwidthA->GetXaxis()->SetRangeUser(0.,50.); 
	  hwidthA->Draw();
      
	  cOut->cd(4); cOut->GetPad(4)->SetLogy();cOut->GetPad(4)->SetLogz(0);
	  //hwidthC->GetXaxis()->SetRangeUser(0.,50.); 
	  hwidthC->Draw();
      
	  cOut->Print(Form("QA_Run_%d.pdf",runNumber));
	  //-------------
      
	  cOut->Clear();
	  cOut->Divide(2,2);
	  cOut->cd(1); cOut->GetPad(1)->SetLogy(0);cOut->GetPad(1)->SetLogz();
	  htimepmt->Draw("colz");
      
	  cOut->cd(2); cOut->GetPad(2)->SetLogy(0);cOut->GetPad(2)->SetLogz();
	  //hwidthpmt->GetYaxis()->SetRangeUser(0.,50.); 
	  hwidthpmt->Draw("colz");
      
	  cOut->cd(3); cOut->GetPad(3)->SetLogy(0);cOut->GetPad(3)->SetLogz();
	  //hadcwidthA->GetYaxis()->SetRangeUser(0.,50.); 
	  hadcwidthA->Draw("colz");
      
	  cOut->cd(4); cOut->GetPad(4)->SetLogy(0);cOut->GetPad(4)->SetLogz();
	  //hadcwidthC->GetYaxis()->SetRangeUser(0.,50.); 
	  hadcwidthC->Draw("colz");
      
	  cOut->Print(Form("QA_Run_%d.pdf",runNumber));
	  //-------------
      
	  cOut->Clear();
	  cOut->Divide(2,2);
	  cOut->cd(1); cOut->GetPad(1)->SetLogy(0);cOut->GetPad(1)->SetLogz();
	  hAdcTimeA->Draw("colz");
      
	  cOut->cd(2); cOut->GetPad(2)->SetLogy(0);cOut->GetPad(2)->SetLogz();
	  hAdcTimeC->Draw("colz");
      
	  cOut->cd(3); cOut->GetPad(3)->SetLogz();
	  hTriggerDecision->Draw("text");
      
	  cOut->cd(4); cOut->GetPad(4)->SetLogz(); cOut->GetPad(4)->SetLogz(0);
	  htimecorr->Draw("colz");
      
	  cOut->Print(Form("QA_Run_%d.pdf",runNumber));
	  //-------------
  
	  cOut->Clear();
	  cOut->Divide(2,2);
	  cOut->cd(1);  cOut->GetPad(1)->SetLogy(1);cOut->GetPad(1)->SetLogz(0);
	  hV0A->GetXaxis()->SetRangeUser(0.,33.);hV0A->Draw();
  
	  cOut->cd(2); cOut->GetPad(2)->SetLogy(1);cOut->GetPad(2)->SetLogz(0);
	  hV0C->GetXaxis()->SetRangeUser(0.,33.);hV0C->Draw();
      
	  if(hTotRecoMult)
	    {
	      cOut->cd(3); cOut->GetPad(3)->SetLogy(0);cOut->GetPad(3)->SetLogz(0); cOut->GetPad(3)->SetGridy(1);
	      hTotRecoMult->Sumw2();
	      if(hTotRecoMult_CVLN)
		{
		  hTotRecoMult_CVLN->Sumw2();
		  hTotRecoMult_CVHN->Sumw2();
		  hTotRecoMult_CVLN->Divide(hTotRecoMult);
		  hTotRecoMult_CVHN->Divide(hTotRecoMult);
		  hTotRecoMult_CVLN->Draw("e"); hTotRecoMult_CVLN->SetLineColor(4);hTotRecoMult_CVLN->SetMarkerStyle(0);hTotRecoMult_CVLN->SetMaximum(1.2);
		  hTotRecoMult_CVLN->SetTitle("Multiplicity efficiency CVLN (blue) and CHVN (red)");
		  hTotRecoMult_CVHN->Draw("esame"); hTotRecoMult_CVHN->SetLineColor(2);hTotRecoMult_CVHN->SetMarkerStyle(0);
		}
	      cOut->cd(4); cOut->GetPad(4)->SetLogy(0);cOut->GetPad(4)->SetLogz(1); cOut->GetPad(4)->SetGridx(1); cOut->GetPad(4)->SetGridy(1);
	      hRecoMult->Draw("colz");
	    }
      
	  cOut->Print(Form("QA_Run_%d.pdf",runNumber));
	  //-------------
      
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
	  //-------------
  
	  delete cOut;
	}
    }
  return 0;
}

bool HigtVoltage(Int_t run)
{
  bool result=true;
  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://");
  man->SetRun(run);

  AliCDBEntry *ent = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calData = (AliVZEROCalibData*)ent->GetObject();

  int dead=calData->IsChannelDead(0);
  if(dead==1)
    result=false;
  return result;
}
