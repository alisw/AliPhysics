#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include "AliITSgeomTGeo.h"
#endif

Double_t LangausFun(Double_t *x, Double_t *par);
void MakePlots(TString ntupleFileName);


void TrendQAtrainSDD(TString period,
		     TString recoPass,
		     TString qaTrain="QA",
		     Int_t firstRun=0,
		     Int_t lastRun=999999999,
		     TString fileName="QAresults.root"){

  gStyle->SetOptStat(0);
  

  TGrid::Connect("alien:");
  Int_t year=0;
  if(period.Contains("LHC09")) year=2009;
  else if(period.Contains("LHC10")) year=2010;
  else if(period.Contains("LHC11")) year=2011;

  TString outFilNam=Form("TrendingSDD_%s_%s_%s.root",period.Data(),recoPass.Data(),qaTrain.Data());


  const Int_t nVariables=35;
  TNtuple* ntsdd=new TNtuple("ntsdd","SDD trending","nrun:fracTrackWithClu1:errfracTrackWithClu1:fracTrackWithClu2:errfracTrackWithClu2:fracTrackWithClu3:errfracTrackWithClu3:fracTrackWithClu4:errfracTrackWithClu4:fracTrackWithClu5:errfracTrackWithClu5:fracTrackWithClu6:errfracTrackWithClu6:meanTrPts3:errmeanTrPts3:meanTrPts4:errmeanTrPts4:minDrTime:errminDrTime:meanDrTime:errmeanDrTime:fracExtra:errfracExtra:meandEdxLay3:errmeandEdxLay3:meandEdxLay4:errmeandEdxLay4:meandEdxTB0:errmeandEdxTB0:meandEdxTB5:errmeandEdxTB5:nMod95:nMod80:nMod60:nModEmpty");
  Float_t xnt[nVariables];

  TBits* readRun=new TBits(999999);
  readRun->ResetAllBits();
  if(!gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",outFilNam.Data()))){
    TFile* oldfil=new TFile(outFilNam.Data());
    TNtuple* ntmp=(TNtuple*)oldfil->Get("ntsdd");
    Bool_t isOK=kFALSE;
    if(ntmp){
      if(ntmp->GetNvar()==ntsdd->GetNvar()){
	isOK=kTRUE;
	TObjArray* arr1=(TObjArray*)ntsdd->GetListOfBranches();
	TObjArray* arr2=(TObjArray*)ntmp->GetListOfBranches();
	for(Int_t iV=0; iV<ntmp->GetNvar(); iV++){
	  TString vnam1=arr1->At(iV)->GetName();
	  TString vnam2=arr2->At(iV)->GetName();
	  if(vnam1!=vnam2) isOK=kFALSE;
	  ntmp->SetBranchAddress(vnam2.Data(),&xnt[iV]);
	}
	if(isOK){
	  for(Int_t nE=0; nE<ntmp->GetEntries(); nE++){
	    ntmp->GetEvent(nE);
	    Int_t theRun=(Int_t)(xnt[0]+0.0001);
	    readRun->SetBitNumber(theRun);
	    ntsdd->Fill(xnt);
	  }
	}
      }
    }
    if(!isOK){
      printf("Ntuple in local file not OK -> will be recreated\n");
    }
    oldfil->Close();
    delete oldfil;
  }

  if(!gGrid||!gGrid->IsConnected()) {
    printf("gGrid not found! exit macro\n");
    return;
  }

  TString  path=Form("/alice/data/%d/%s/",year,period.Data());
  TGridResult *gr = gGrid->Query(path,fileName);
  Int_t nFiles = gr->GetEntries();
  printf("================>%d files found\n", nFiles);
  if (nFiles < 1) return;


  if (nFiles > 1){
    for (Int_t iFil = 0; iFil <nFiles ; iFil++) { 
      TString fileNameLong=Form("%s",gr->GetKey(iFil,"turl"));
      if(!fileNameLong.Contains(recoPass.Data())) continue;
      if(!fileNameLong.Contains(qaTrain.Data())) continue;
      TString runNumber=fileNameLong;
      runNumber.ReplaceAll(Form("alien:///alice/data/%d/%s/",year,period.Data()),"");
      runNumber.Remove(9,runNumber.Sizeof());
   
      Int_t iRun=atoi(runNumber.Data());
      if(readRun->TestBitNumber(iRun)){ 
	printf("Run %d aleady in local ntuple -> skipping it\n",iRun);
	continue;
      }
      if(iRun<firstRun) continue;
      if(iRun>lastRun) continue;

      TString isMerged=fileNameLong;
      isMerged.Remove(isMerged.Sizeof()-16); 
      isMerged.Remove(0,isMerged.Sizeof()-5);
      if(!isMerged.Contains("QA")) continue;
      printf("Open File %s  Run %d\n",fileNameLong.Data(),iRun);



      TFile* f=TFile::Open(fileNameLong.Data());  

      TDirectoryFile* df=(TDirectoryFile*)f->Get("SDD_Performance");
      if(!df){
	printf("Run %d SDD_Performance MISSING -> Exit\n",iRun);
	continue;
      }
      TList* l=(TList*)df->Get("coutputRP");
      if(!df){
	printf("Run %d coutputRP TList MISSING -> Exit\n",iRun);
	continue;
      }  
      TH1F* hcllay=(TH1F*)l->FindObject("hCluInLay");
      Float_t fracT[6]={0.,0.,0.,0.,0.,0.};
      Float_t efracT[6]={0.,0.,0.,0.,0.,0.};
      if(hcllay->GetBinContent(1)>0){
	for(Int_t iLay=0; iLay<6; iLay++){
	  fracT[iLay]=hcllay->GetBinContent(iLay+2)/hcllay->GetBinContent(1);
	  efracT[iLay]=TMath::Sqrt(fracT[iLay]*(1-fracT[iLay])/hcllay->GetBinContent(1));
	}
      }
      TH1F* hmodT=(TH1F*)l->FindObject("hTPMod");
      TH1F* hgamod=(TH1F*)l->FindObject("hGAMod");
      Int_t bestMod=0;
      for(Int_t iMod=0; iMod<260;iMod++){
	Int_t gda=(Int_t)hgamod->GetBinContent(iMod+1);
	if(gda>bestMod) bestMod=gda;
      }
      Int_t nChunks=1;
      if(bestMod>512){
	nChunks=(Int_t)(bestMod/512.+0.5);
      }
      hgamod->Scale(1./nChunks);

      TH1F* hev=(TH1F*)l->FindObject("hNEvents");
      Int_t nTotEvents=hev->GetBinContent(2);
      Int_t nTrigEvents=hev->GetBinContent(3);
      Int_t nEvents=nTotEvents;
      printf("Run %d Number of Events = %d Triggered=%d\n",iRun,nTotEvents,nTrigEvents);
      if(nTrigEvents>0){ 
	nEvents=nTrigEvents;
      }
      if(nTotEvents==0) continue;
      Int_t nModGood3=0;
      Int_t nModGood4=0;
      Int_t nModBadAn=0;
      Float_t sumtp3=0;
      Float_t sumtp4=0;
      Float_t sumEtp3=0;
      Float_t sumEtp4=0;
      for(Int_t iMod=0; iMod<260; iMod++){
	Float_t tps=hmodT->GetBinContent(iMod+1);
	Float_t ga=hgamod->GetBinContent(iMod+1);
	if(ga<500) nModBadAn++;
	Float_t tpsN=0.;
	Float_t etpsN=0.;
	if(ga>0){
	  tpsN=tps/ga/(Float_t)nEvents;
	  etpsN=TMath::Sqrt(tps)/ga/(Float_t)nEvents;
	  if(iMod<84){
	    sumtp3+=tpsN;
	    sumEtp3+=(etpsN*etpsN);
	    nModGood3++;
	  }
	  else{
	    sumtp4+=tpsN;
	    sumEtp4+=(etpsN*etpsN);
	    nModGood4++;
	  }
	}
      }

      TH1F* hapmod=(TH1F*)l->FindObject("hAllPmod");
      TH1F* hgpmod=(TH1F*)l->FindObject("hGoodPmod");
      //     TH1F* hmpmod=(TH1F*)l->FindObject("hMissPmod");
      TH1F* hbrmod=(TH1F*)l->FindObject("hBadRegmod");
      TH1F* hskmod=(TH1F*)l->FindObject("hSkippedmod");
      TH1F* hoamod=(TH1F*)l->FindObject("hOutAccmod");
      TH1F* hnrmod=(TH1F*)l->FindObject("hNoRefitmod");
      Int_t nBelow95=0;
      Int_t nBelow80=0;
      Int_t nBelow60=0;
      Int_t nZeroP=0;
      for(Int_t imod=0; imod<260;imod++){
	Float_t numer=hgpmod->GetBinContent(imod+1)+hbrmod->GetBinContent(imod+1)+hoamod->GetBinContent(imod+1)+hnrmod->GetBinContent(imod+1)+hskmod->GetBinContent(imod+1);
	Float_t denom=hapmod->GetBinContent(imod+1);
	if(denom>0){
	  Float_t eff=numer/denom;
	  if(eff<0.95) nBelow95++;
	  if(eff<0.80) nBelow80++;
	  if(eff<0.60) nBelow60++;
	}
	if(hmodT->GetBinContent(imod+1)<1.){
	  nZeroP++;
	}	
      }

      TH1F* htimT=(TH1F*)l->FindObject("hDrTimTPAll");
      TH1F* htimTe=(TH1F*)l->FindObject("hDrTimTPExtra");
      
      Double_t fracExtra=0.;
      Double_t errFracExtra=0.;
      if(htimT->GetEntries()>0){
	fracExtra=htimTe->GetEntries()/htimT->GetEntries();
	errFracExtra=TMath::Sqrt(htimTe->GetEntries())/htimT->GetEntries();
      }
      Double_t averPoints=0.;
      Double_t cntBins=0.;
      for(Int_t iBin=1; iBin<=htimT->GetNbinsX(); iBin++){
	Float_t tim=htimT->GetBinCenter(iBin);
	if(tim>2000. && tim<4000.){
	  averPoints+=htimT->GetBinContent(iBin);
	  cntBins+=1;
	}
      }
      Double_t minTime=-999.;
      Double_t errMinTime=0.;
      if(cntBins>0){ 
	averPoints/=cntBins;
	for(Int_t iBin=1; iBin<=htimT->GetNbinsX(); iBin++){
	  if(htimT->GetBinContent(iBin)>0.5*averPoints){
	    minTime=htimT->GetBinCenter(iBin);
	    errMinTime=0.5*htimT->GetBinWidth(iBin);
	    break;
	  }
	}
      }
      TH2F* hdedxmod=(TH2F*)l->FindObject("hdEdxVsMod");
      TH1D* hdedxLay3=hdedxmod->ProjectionY("hdedxLay3",1,84);
      TH1D* hdedxLay4=hdedxmod->ProjectionY("hdedxLay4",85,260);
      
      TH1F* hSigTim0=(TH1F*)l->FindObject("hSigTimeInt0");
      TH1F* hSigTim5=(TH1F*)l->FindObject("hSigTimeInt5");

      Int_t index=0;
      xnt[index++]=(Float_t)iRun;
      xnt[index++]=fracT[0];
      xnt[index++]=efracT[0];
      xnt[index++]=fracT[1];
      xnt[index++]=efracT[1];
      xnt[index++]=fracT[2];
      xnt[index++]=efracT[2];
      xnt[index++]=fracT[3];
      xnt[index++]=efracT[3];
      xnt[index++]=fracT[4];
      xnt[index++]=efracT[4];
      xnt[index++]=fracT[5];
      xnt[index++]=efracT[5];
      xnt[index++]=sumtp3/nModGood3;
      xnt[index++]=TMath::Sqrt(sumEtp3)/nModGood3;
      xnt[index++]=sumtp4/nModGood4;
      xnt[index++]=TMath::Sqrt(sumEtp4)/nModGood4;
      xnt[index++]=minTime;
      xnt[index++]=errMinTime;
      xnt[index++]=htimT->GetMean();
      xnt[index++]=htimT->GetMeanError();
      xnt[index++]=fracExtra;
      xnt[index++]=errFracExtra;
      xnt[index++]=hdedxLay3->GetMean();
      xnt[index++]=hdedxLay3->GetMeanError();
      xnt[index++]=hdedxLay4->GetMean();
      xnt[index++]=hdedxLay4->GetMeanError();
      xnt[index++]=hSigTim0->GetMean();
      xnt[index++]=hSigTim0->GetMeanError();
      xnt[index++]=hSigTim5->GetMean();
      xnt[index++]=hSigTim5->GetMeanError();
      xnt[index++]=(Float_t)nBelow95;
      xnt[index++]=(Float_t)nBelow80;
      xnt[index++]=(Float_t)nBelow60;
      xnt[index++]=(Float_t)nZeroP;
      ntsdd->Fill(xnt);
    }
  }

  TFile* outfil=new TFile(outFilNam.Data(),"recreate");
  outfil->cd();
  ntsdd->Write();
  outfil->Close();

  MakePlots(outFilNam);

}

void MakePlots(TString ntupleFileName){
  TFile* fil=new TFile(ntupleFileName.Data(),"read");
  if(!fil){
    printf("File with ntuple does not exist\n");
    return;
  }
  TNtuple* ntsdd=(TNtuple*)fil->Get("ntsdd");

  Float_t nrun;
  Float_t meanTrPts3,errmeanTrPts3,meanTrPts4,errmeanTrPts4;
  Float_t minDrTime,errminDrTime,meanDrTime,errmeanDrTime;
  Float_t fracTrackWithClu1,fracTrackWithClu2,errfracTrackWithClu1,errfracTrackWithClu2;
  Float_t fracTrackWithClu3,fracTrackWithClu4,errfracTrackWithClu3,errfracTrackWithClu4;
  Float_t fracTrackWithClu5,fracTrackWithClu6,errfracTrackWithClu5,errfracTrackWithClu6;
  Float_t fracExtra,errfracExtra;
  Float_t meandEdxLay3,errmeandEdxLay3,meandEdxLay4,errmeandEdxLay4;
  Float_t meandEdxTB0,errmeandEdxTB0,meandEdxTB5,errmeandEdxTB5;
  Float_t nMod95,nMod80,nMod60,nModEmpty;
  
  ntsdd->SetBranchAddress("nrun",&nrun);
  ntsdd->SetBranchAddress("fracTrackWithClu1",&fracTrackWithClu1);
  ntsdd->SetBranchAddress("errfracTrackWithClu1",&errfracTrackWithClu1);
  ntsdd->SetBranchAddress("fracTrackWithClu2",&fracTrackWithClu2);
  ntsdd->SetBranchAddress("errfracTrackWithClu2",&errfracTrackWithClu2);
  ntsdd->SetBranchAddress("fracTrackWithClu3",&fracTrackWithClu3);
  ntsdd->SetBranchAddress("errfracTrackWithClu3",&errfracTrackWithClu3);
  ntsdd->SetBranchAddress("fracTrackWithClu4",&fracTrackWithClu4);
  ntsdd->SetBranchAddress("errfracTrackWithClu4",&errfracTrackWithClu4);
  ntsdd->SetBranchAddress("fracTrackWithClu5",&fracTrackWithClu5);
  ntsdd->SetBranchAddress("errfracTrackWithClu5",&errfracTrackWithClu5);
  ntsdd->SetBranchAddress("fracTrackWithClu6",&fracTrackWithClu6);
  ntsdd->SetBranchAddress("errfracTrackWithClu6",&errfracTrackWithClu6);
  ntsdd->SetBranchAddress("nMod95",&nMod95);
  ntsdd->SetBranchAddress("nMod80",&nMod80);
  ntsdd->SetBranchAddress("nMod60",&nMod60);
  ntsdd->SetBranchAddress("nModEmpty",&nModEmpty);

  ntsdd->SetBranchAddress("meanTrPts3",&meanTrPts3);
  ntsdd->SetBranchAddress("errmeanTrPts3",&errmeanTrPts3);
  ntsdd->SetBranchAddress("meanTrPts4",&meanTrPts4);
  ntsdd->SetBranchAddress("errmeanTrPts4",&errmeanTrPts4);
  ntsdd->SetBranchAddress("minDrTime",&minDrTime);
  ntsdd->SetBranchAddress("errminDrTime",&errminDrTime);
  ntsdd->SetBranchAddress("meanDrTime",&meanDrTime);
  ntsdd->SetBranchAddress("errmeanDrTime",&errmeanDrTime);
  ntsdd->SetBranchAddress("fracExtra",&fracExtra);
  ntsdd->SetBranchAddress("errfracExtra",&errfracExtra);
  ntsdd->SetBranchAddress("meandEdxTB0",&meandEdxTB0);
  ntsdd->SetBranchAddress("errmeandEdxTB0",&errmeandEdxTB0);
  ntsdd->SetBranchAddress("meandEdxTB5",&meandEdxTB5);
  ntsdd->SetBranchAddress("errmeandEdxTB5",&errmeandEdxTB5);
  ntsdd->SetBranchAddress("meandEdxLay3",&meandEdxLay3);
  ntsdd->SetBranchAddress("errmeandEdxLay3",&errmeandEdxLay3);
  ntsdd->SetBranchAddress("meandEdxLay4",&meandEdxLay4);
  ntsdd->SetBranchAddress("errmeandEdxLay4",&errmeandEdxLay4);

  TH1F* histotrp3=new TH1F("histotrp3","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histotrp4=new TH1F("histotrp4","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histominTime=new TH1F("histominTime","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histomeanTime=new TH1F("histomeanTime","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histofracExtra=new TH1F("histofracExtra","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histodEdxTB0=new TH1F("histodEdxTB0","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histodEdxTB5=new TH1F("histodEdxTB5","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histodEdxLay3=new TH1F("histodEdxLay3","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histodEdxLay4=new TH1F("histodEdxLay4","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu1=new TH1F("histoTrackClu1","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu2=new TH1F("histoTrackClu2","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu3=new TH1F("histoTrackClu3","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu4=new TH1F("histoTrackClu4","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu5=new TH1F("histoTrackClu5","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoTrackClu6=new TH1F("histoTrackClu6","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());

  TH1F* histoNmodEffBelow95=new TH1F("histoNmodEffBelow95","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoNmodEffBelow80=new TH1F("histoNmodEffBelow80","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoNmodEffBelow60=new TH1F("histoNmodEffBelow60","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());
  TH1F* histoNmodEmpty=new TH1F("histoNmodEmpty","",(Int_t)ntsdd->GetEntries(),0.,ntsdd->GetEntries());

  for(Int_t i=0; i<ntsdd->GetEntries();i++){
    ntsdd->GetEvent(i);
    histoTrackClu1->SetBinContent(i+1,fracTrackWithClu1);
    histoTrackClu1->SetBinError(i+1,errfracTrackWithClu1);
    histoTrackClu1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoTrackClu2->SetBinContent(i+1,fracTrackWithClu2);
    histoTrackClu2->SetBinError(i+1,errfracTrackWithClu2);
    histoTrackClu2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoTrackClu3->SetBinContent(i+1,fracTrackWithClu3);
    histoTrackClu3->SetBinError(i+1,errfracTrackWithClu3);
    histoTrackClu3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoTrackClu4->SetBinContent(i+1,fracTrackWithClu4);
    histoTrackClu4->SetBinError(i+1,errfracTrackWithClu4);
    histoTrackClu4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoTrackClu5->SetBinContent(i+1,fracTrackWithClu5);
    histoTrackClu5->SetBinError(i+1,errfracTrackWithClu5);
    histoTrackClu5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoTrackClu6->SetBinContent(i+1,fracTrackWithClu6);
    histoTrackClu6->SetBinError(i+1,errfracTrackWithClu6);
    histoTrackClu6->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histominTime->SetBinContent(i+1,minDrTime);
    histominTime->SetBinError(i+1,errminDrTime);
    histominTime->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histomeanTime->SetBinContent(i+1,meanDrTime);
    histomeanTime->SetBinError(i+1,errmeanDrTime);
    histomeanTime->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histotrp3->SetBinContent(i+1,meanTrPts3);
    histotrp3->SetBinError(i+1,errmeanTrPts3);
    histotrp3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histotrp4->SetBinContent(i+1,meanTrPts4);
    histotrp4->SetBinError(i+1,errmeanTrPts3);
    histotrp4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histofracExtra->SetBinContent(i+1,fracExtra);
    histofracExtra->SetBinError(i+1,errfracExtra);
    histofracExtra->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histodEdxTB0->SetBinContent(i+1,meandEdxTB0);
    histodEdxTB0->SetBinError(i+1,errmeandEdxTB0);
    histodEdxTB0->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histodEdxTB5->SetBinContent(i+1,meandEdxTB5);
    histodEdxTB5->SetBinError(i+1,errmeandEdxTB5);
    histodEdxTB5->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histodEdxLay3->SetBinContent(i+1,meandEdxLay3);
    histodEdxLay3->SetBinError(i+1,errmeandEdxLay3);
    histodEdxLay3->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histodEdxLay4->SetBinContent(i+1,meandEdxLay4);
    histodEdxLay4->SetBinError(i+1,errmeandEdxLay4);
    histodEdxLay4->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));

    histoNmodEffBelow95->SetBinContent(i+1,nMod95);
    histoNmodEffBelow95->SetBinError(i+1,0.0001);
    histoNmodEffBelow95->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoNmodEffBelow80->SetBinContent(i+1,nMod80);
    histoNmodEffBelow80->SetBinError(i+1,0.0001);
    histoNmodEffBelow80->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoNmodEffBelow60->SetBinContent(i+1,nMod60);
    histoNmodEffBelow60->SetBinError(i+1,0.0001);
    histoNmodEffBelow60->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    histoNmodEmpty->SetBinContent(i+1,nModEmpty);
    histoNmodEmpty->SetBinError(i+1,0.0001);
    histoNmodEmpty->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
  }

  gStyle->SetOptStat(0);


  TCanvas* c1=new TCanvas("c1","Occupancy");
  histotrp3->SetLineColor(1);
  histotrp3->SetMarkerStyle(20);
  histotrp3->Draw();
  histotrp3->SetMinimum(0.);
  histotrp4->SetLineColor(2);
  histotrp4->SetMarkerColor(2);
  histotrp4->SetMarkerStyle(22);
  histotrp4->Draw("same");
  histotrp3->GetYaxis()->SetTitle("Track Point Occupancy");
  TLegend* leg=new TLegend(0.7,0.15,0.88,0.35);
  TLegendEntry* ent=leg->AddEntry(histotrp3,"Layer3","PL");
  ent=leg->AddEntry(histotrp4,"Layer4","PL");
  ent->SetTextColor(histotrp4->GetMarkerColor());
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Update();

  TCanvas* c2=new TCanvas("c2","DriftTime",600,900);
  c2->Divide(1,2);
  c2->cd(1);
  histominTime->Draw();
  histominTime->SetMinimum(250);
  histominTime->SetMaximum(1000);
  histominTime->GetYaxis()->SetTitle("Minimum Drift Time (ns)");
  c2->cd(2);
  histomeanTime->Draw();
  histomeanTime->GetYaxis()->SetTitle("Average Drift Time (ns)");
  c2->Update();

  TCanvas* c3=new TCanvas("c3","ExtraClusters");
  histofracExtra->Draw();
  histofracExtra->GetYaxis()->SetTitle("Fraction of Extra Clusters");
  c3->Update();


  TCanvas* c4=new TCanvas("c4","Charge");
  histodEdxTB0->SetLineColor(1);
  histodEdxTB0->SetMarkerStyle(20);
  histodEdxTB0->Draw();
  histodEdxTB0->SetMinimum(0.);
  histodEdxTB5->SetLineColor(4);
  histodEdxTB5->SetMarkerColor(4);
  histodEdxTB5->SetMarkerStyle(23);
  histodEdxTB5->Draw("same");
  histodEdxTB0->GetYaxis()->SetTitle("<dE/dx> (keV/300 #mum)");
  TLegend* leg2=new TLegend(0.6,0.15,0.88,0.35);
  ent=leg2->AddEntry(histodEdxTB0,"Small drift time","PL");
  ent=leg2->AddEntry(histodEdxTB5,"Large drift time","PL");
  ent->SetTextColor(histodEdxTB5->GetMarkerColor());
  leg2->SetFillStyle(0);
  leg2->Draw();
  c4->Update();

  TCanvas* c4b=new TCanvas("c4b","Charge per Layer");
  histodEdxLay3->SetLineColor(1);
  histodEdxLay3->SetMarkerStyle(20);
  histodEdxLay3->Draw();
  histodEdxLay3->SetMinimum(0.);
  histodEdxLay4->SetLineColor(4);
  histodEdxLay4->SetMarkerColor(4);
  histodEdxLay4->SetMarkerStyle(23);
  histodEdxLay4->Draw("same");
  histodEdxLay3->GetYaxis()->SetTitle("<dE/dx> (keV/300 #mum)");
  TLegend* leg2b=new TLegend(0.6,0.15,0.88,0.35);
  ent=leg2b->AddEntry(histodEdxLay3,"Layer 3","PL");
  ent=leg2b->AddEntry(histodEdxLay4,"Layer 4","PL");
  ent->SetTextColor(histodEdxLay4->GetMarkerColor());
  leg2b->SetFillStyle(0);
  leg2b->Draw();
  c4b->Update();

  TCanvas* c5=new TCanvas("c5","TrackWithSDD");
  histoTrackClu3->Draw();
  histoTrackClu3->SetLineColor(1);
  histoTrackClu3->SetMarkerStyle(20);
  histoTrackClu3->Draw();
  histoTrackClu3->SetMinimum(0.);
  histoTrackClu4->SetLineColor(2);
  histoTrackClu4->SetMarkerColor(2);
  histoTrackClu4->SetMarkerStyle(22);
  histoTrackClu4->Draw("same");
  histoTrackClu1->SetLineColor(kGray+1);
  histoTrackClu1->SetMarkerColor(kGray+1);
  histoTrackClu1->SetMarkerStyle(24);
  histoTrackClu1->Draw("same");
  histoTrackClu2->SetLineColor(kGray+2);
  histoTrackClu2->SetMarkerColor(kGray+2);
  histoTrackClu2->SetMarkerStyle(26);
  histoTrackClu2->Draw("same");
  histoTrackClu5->SetLineColor(4);
  histoTrackClu5->SetMarkerColor(4);
  histoTrackClu5->SetMarkerStyle(29);
  histoTrackClu5->Draw("same");
  histoTrackClu6->SetLineColor(kBlue+1);
  histoTrackClu6->SetMarkerColor(kBlue+1);
  histoTrackClu6->SetMarkerStyle(30);
  histoTrackClu6->Draw("same");
  histoTrackClu3->GetYaxis()->SetTitle("Fraction of Tracks with Point in SDD Layers");
  TLegend* leg3=new TLegend(0.7,0.15,0.88,0.35);
  ent=leg3->AddEntry(histoTrackClu1,"Layer1","PL");
  ent->SetTextColor(histoTrackClu1->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu2,"Layer2","PL");
  ent->SetTextColor(histoTrackClu2->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu3,"Layer3","PL");
  ent->SetTextColor(histoTrackClu3->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu4,"Layer4","PL");
  ent->SetTextColor(histoTrackClu4->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu5,"Layer5","PL");
  ent->SetTextColor(histoTrackClu5->GetMarkerColor());
  ent=leg3->AddEntry(histoTrackClu6,"Layer6","PL");
  ent->SetTextColor(histoTrackClu6->GetMarkerColor());

  leg3->SetFillStyle(0);
  leg3->Draw();
  c5->Update();

  TCanvas* c6=new TCanvas("c6","Modules with low eff",800,1000);
  c6->Divide(1,4);
  c6->cd(1);
  histoNmodEffBelow95->Draw("E");
  c6->cd(2);
  histoNmodEffBelow80->Draw("E");
  c6->cd(3);
  histoNmodEffBelow60->Draw("E");
  c6->cd(4);
  histoNmodEmpty->Draw("E");

}

Double_t LangausFun(Double_t *x, Double_t *par) {

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
  Double_t np = 100.0;      // number of convolution steps
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
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}
