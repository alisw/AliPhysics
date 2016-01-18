/*
  
  find    /hera/alice/miranov/alice-tpc-notes/reconstruction/dataProductionPreparation/ATO-240/data/benchmarkFull/testATO240_07_31_15_Gitv5-06-35-228-gbf86468/testATO240_07_31_15_Gitv5-06-35-228-gbf86468/2015/LHC15f/000225106/cpass1/ -name AliESDs_Barrel.root > esd.list
  
   Int_t run=225106
   .L $ALICE_PHYSICS/../src/PWGPP/CalibMacros/CPass0/ConfigCalibTrain.C
   ConfigCalibTrain(run,"local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB/");

   .x $NOTES/aux/rootlogon.C

   .L $ALICE_PHYSICS/../src/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolationMacro.C+
   // cp   $ALICE_PHYSICS/../src/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolationMacro.C AliTPCcalibAlignInterpolationMacro.C 

*/


#include  "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include "TSystem.h"
#include "THn.h"
#include "TStatToolkit.h"
#include "TTreeStream.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"
//
#include "AliSysInfo.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliExternalTrackParam.h"
#include "AliTPCcalibAlignInterpolation.h"
#include "AliNDLocalRegression.h"



TTree * treeDist =0;   // global tree used to check the content of the fit

void AliTPCcalibAlignInterpolationLoop(Int_t nchunks, Int_t neventsMax);
void  CreateDistortionMapsFromFile(Int_t type, const char * inputFile, const char *outputFile, Int_t delta=1);
void MakeNDFit(const char * inputFile, const char * inputTree, Float_t sector0,  Float_t sector1,  Float_t theta0, Float_t theta1);


void AliTPCcalibAlignInterpolationMacro(Int_t action,Int_t param0=0, Int_t param1=0){
  if (action==0)   AliTPCcalibAlignInterpolationLoop(100000,1000000000);
  if (action==1)    {
    Int_t startTime=TString(gSystem->Getenv("mapStartTime")).Atoi();
    Int_t stopTime =TString(gSystem->Getenv("mapStopTime" )).Atoi();
    AliTPCcalibAlignInterpolation::FillHistogramsFromChain("residual.list", 6,4,startTime, stopTime, param1,param0);
  }
  if (action==4)    AliTPCcalibAlignInterpolation::FillHistogramsFromStreamers("residual.list",1,1,1);

  if (action==2) {
    Int_t startTime=TString(gSystem->Getenv("mapStartTime")).Atoi();
    TString outFile = TString::Format("ResidualHistograms_His%d.root",param0);
    if (startTime>0)  outFile = TString::Format("ResidualHistograms_His%d_Time%d.root",param0,startTime);
    CreateDistortionMapsFromFile(param0, outFile.Data(),TString::Format("ResidualMapFull_%d.root",param0).Data(), param1);
  }
  if (action==3){ //make distortion fit
    TString inputFile=gSystem->Getenv("inputFile");
    TString inputTree=gSystem->Getenv("inputTree");
    Int_t runNumber=TString(gSystem->Getenv("runNumber")).Atoi();
    Int_t sec0=TString(gSystem->Getenv("varSec0")).Atoi();
    Int_t sec1=TString(gSystem->Getenv("varSec1")).Atoi();
    Float_t theta0=TString(gSystem->Getenv("varTheta0")).Atof();
    Float_t theta1=TString(gSystem->Getenv("varTheta1")).Atof();
    ::Info("AliTPCcalibAlignInterpolationM","MakeFit");
    printf("Env varaible:\n%s\t%s\t%s\t%s\t%s\t%s\t\n",inputFile.Data(), inputTree.Data(), gSystem->Getenv("varSec0") , gSystem->Getenv("varSec1"), 	gSystem->Getenv("varTheta0"), gSystem->Getenv("varTheta1"));	          
    printf("%s\t%s\t%d\t%d\t%f\t%f\t\n",inputFile.Data(), inputTree.Data(), sec0,sec1,theta0,theta1);		      

    MakeNDFit(inputFile, inputTree, sec0,sec1,theta0,theta1);		      
  } 

  if (action==5) {
    Int_t runNumber=TString(gSystem->Getenv("runNumber")).Atoi();
     AliTPCcalibAlignInterpolation::MakeEventStatInfo("cat residual.list",300,param0,1);
  }
 
}

void AliTPCcalibAlignInterpolationLoop(Int_t nchunks, Int_t neventsMax){
  //
  //
  //
  AliTPCcalibAlignInterpolation * calibInterpolation = new  AliTPCcalibAlignInterpolation("calibInterpolation","calibInterpolation",kFALSE);
  calibInterpolation->SetStreamLevel( AliTPCcalibAlignInterpolation::kStremInterpolation);
  calibInterpolation->SetSyswatchStep(200);
  //
  TString  esdList0 = gSystem->GetFromPipe("cat esd.list");
  TObjArray *esdArray= esdList0.Tokenize("\n");  
  Int_t nesd = esdArray->GetEntriesFast();
  
  AliESDEvent *esdEvent = new AliESDEvent();
  AliESDfriend *esdFriend = 0;
  
  //
  for (Int_t iesd=0; iesd<TMath::Min(nesd,nchunks); iesd++){
    printf("%d\n",iesd);
    TFile *esdFile = TFile::Open(esdArray->At(iesd)->GetName(),"read");
    if (!esdFile) return;
    TTree *esdTree = (TTree*)esdFile->Get("esdTree");
    esdEvent->ReadFromTree(esdTree);    
    esdTree->AddFriend("esdFriendTree", TString(gSystem->DirName(esdArray->At(iesd)->GetName()))+"/AliESDfriends_Barrel.root");
    
    esdTree->SetBranchStatus("ESDfriend.", 1);
    esdFriend = (AliESDfriend*)(esdEvent->FindListObject("AliESDfriend"));
    if (!esdFriend) return;
    if (esdFriend) esdTree->SetBranchAddress("ESDfriend.", &esdFriend);
    Int_t nEvents=TMath::Min(Int_t(esdTree->GetEntries()), neventsMax);
    for (Int_t iEvent=0; iEvent<nEvents; iEvent++){
      esdTree->GetEvent(iEvent);
      printf("%d\t%d\n",iesd,iEvent);
      calibInterpolation->Process(esdEvent);
    }    
 }
  delete  calibInterpolation;  
}

void TestSpeedMapsFromFile(){
  //
  // Aim of test - check the speed of the projection makeing depending on the order of projections
  //
  const char * inputFile="ResidualHistograms.root";
  const char * outputFile="ResidualMapTest.root";
  TFile *fHistos  = TFile::Open(inputFile);
  THnF *histoITSTRD= (THnF*) fHistos->Get("deltaRPhiTPCITSTRD");
  //
  TTreeSRedirector * pcstream = new TTreeSRedirector(outputFile,"recreate");
  TMatrixD projectionInfo(5,5);
  projectionInfo(0,0)=0;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;
  projectionInfo(1,0)=1;  projectionInfo(1,1)=1;  projectionInfo(1,2)=8; 
  projectionInfo(2,0)=2;  projectionInfo(2,1)=1;  projectionInfo(2,2)=4;
  projectionInfo(3,0)=3;  projectionInfo(3,1)=1;  projectionInfo(3,2)=2;
  projectionInfo(4,0)=4;  projectionInfo(4,1)=1;  projectionInfo(4,2)=2;

  AliSysInfo::AddStamp("Test1234_Begin",-1,1);
  TStatToolkit::MakeDistortionMap(4, histoITSTRD, pcstream, projectionInfo); 
  AliSysInfo::AddStamp("Test1234_End",-1,1);  
  printf("Test1234_End\n");  
  projectionInfo(0,0)=0;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;
  projectionInfo(1,0)=2;  projectionInfo(1,1)=1;  projectionInfo(1,2)=8; 
  projectionInfo(2,0)=3;  projectionInfo(2,1)=1;  projectionInfo(2,2)=4;
  projectionInfo(3,0)=4;  projectionInfo(3,1)=1;  projectionInfo(3,2)=2;
  projectionInfo(4,0)=1;  projectionInfo(4,1)=1;  projectionInfo(4,2)=2;
  AliSysInfo::AddStamp("Test2341_Begin",-1,2);
  TStatToolkit::MakeDistortionMap(4, histoITSTRD, pcstream, projectionInfo); 
  AliSysInfo::AddStamp("Test2341_End",-1,2);
  printf("Test2341_End\n");
  //
  projectionInfo(0,0)=0;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;
  projectionInfo(1,0)=3;  projectionInfo(1,1)=1;  projectionInfo(1,2)=8; 
  projectionInfo(2,0)=4;  projectionInfo(2,1)=1;  projectionInfo(2,2)=4;
  projectionInfo(3,0)=1;  projectionInfo(3,1)=1;  projectionInfo(3,2)=2;
  projectionInfo(4,0)=2;  projectionInfo(4,1)=1;  projectionInfo(4,2)=2;
  AliSysInfo::AddStamp("Test3412_Begin",-1,3);
  TStatToolkit::MakeDistortionMap(4, histoITSTRD, pcstream, projectionInfo); 
  AliSysInfo::AddStamp("Test3412_End",-1,3);
  printf("Test3412_End\n");
  //
  projectionInfo(0,0)=0;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;
  projectionInfo(1,0)=4;  projectionInfo(1,1)=1;  projectionInfo(1,2)=8; 
  projectionInfo(2,0)=1;  projectionInfo(2,1)=1;  projectionInfo(2,2)=4;
  projectionInfo(3,0)=2;  projectionInfo(3,1)=1;  projectionInfo(3,2)=2;
  projectionInfo(4,0)=3;  projectionInfo(4,1)=1;  projectionInfo(4,2)=2;
  AliSysInfo::AddStamp("Test4123_Begin",-1,3);
  TStatToolkit::MakeDistortionMap(4, histoITSTRD, pcstream, projectionInfo); 
  AliSysInfo::AddStamp("Test4123_End",-1,3);
  printf("Test4123_End\n");
  delete pcstream;
}

void  CreateDistortionMapsFromFile(Int_t type, const char * inputFile, const char *outputFile, Int_t delta){
  //
  // Create distortion maps from residual histograms
  // TPC cluster to ITS, ITS-TRD and ITS-TOF track fits
  //
  /*
    inputFile="ResidualHistograms.root"
    outputFile="ResidualMap.root"

   */
  TFile *fHistos  = TFile::Open(inputFile);
  
  //  THnF *histoITS = (THnF*) fHistos->Get("deltaRPhiTPCITS");
  //THnF *histoITSTRD= (THnF*) fHistos->Get("deltaRPhiTPCITSTRD");
  //THnF *histoITSTOF = (THnF*) fHistos->Get("deltaRPhiTPCITSTOF");
  //THnF *histoITSZ = (THnF*) fHistos->Get("deltaZTPCITS");
  //THnF *histoITSTRDZ= (THnF*) fHistos->Get("deltaZTPCITSTRD");
  //THnF *histoITSTOFZ = (THnF*) fHistos->Get("deltaZTPCITSTOF");
  
  TTreeSRedirector * pcstream = new TTreeSRedirector(outputFile,"recreate");  
//   TMatrixD projectionInfo(5,5);
//   projectionInfo(0,0)=4;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;
//   projectionInfo(1,0)=0;  projectionInfo(1,1)=1;  projectionInfo(1,2)=delta; 
//   projectionInfo(2,0)=3;  projectionInfo(2,1)=1;  projectionInfo(2,2)=delta;
//   projectionInfo(3,0)=2;  projectionInfo(3,1)=1;  projectionInfo(3,2)=delta;
//   projectionInfo(4,0)=1;  projectionInfo(4,1)=1;  projectionInfo(4,2)=delta;
  //
  TMatrixD projectionInfo(5,5);
  projectionInfo(0,0)=4;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;
  projectionInfo(1,0)=0;  projectionInfo(1,1)=0;  projectionInfo(1,2)=delta; 
  projectionInfo(2,0)=3;  projectionInfo(2,1)=0;  projectionInfo(2,2)=delta;
  projectionInfo(3,0)=2;  projectionInfo(3,1)=0;  projectionInfo(3,2)=delta;
  projectionInfo(4,0)=1;  projectionInfo(4,1)=0;  projectionInfo(4,2)=delta;

  const char *hisNames[6]={"deltaRPhiTPCITS","deltaRPhiTPCITSTRD", "deltaRPhiTPCITSTOF", "deltaZTPCITS","deltaZTPCITSTRD", "deltaZTPCITSTOF"};
  if (type<6){
    printf("Processing-  TStatToolkit::MakeDistortionMap(4, %s, pcstream, projectionInfo)\n", hisNames[type]);
    THnF *histo= (THnF*) fHistos->Get(hisNames[type]);
    TTree * treeMeta = (TTree*)fHistos->Get("metaData");
    if (histo) {
      //TStatToolkit::MakeDistortionMap(4, histo, pcstream, projectionInfo); 
      TStatToolkit::MakeDistortionMapFast(histo, pcstream, projectionInfo,kFALSE); 
      if (treeMeta){
	pcstream->GetFile()->cd();
	TTree * treeMeta2 = treeMeta->CopyTree("1");
	treeMeta2->Write("metaData");
	TTree * treeMap = ((*pcstream)<<hisNames[type]).GetTree();
	const char *query[20]={"startTime","stopTime","meanTime", "minTime", "maxTime","ntracksUsed"};
	const char *aliases[20]={"startTime","stopTime","meanTime", "minTime", "maxTime","ntracksUsed"};
	for (Int_t i=0; i<6; i++){
	  Int_t entries=treeMeta2->Draw(query[i],"","goff");
	  Double_t mean=TMath::Mean(entries,treeMeta2->GetV1());
	  treeMap->SetAlias(aliases[i], TString::Format("(0+%f)",mean));
	}
	for (Int_t icurrent=0; icurrent<8; icurrent++){
	  Int_t entries=treeMeta2->Draw(TString::Format("meanNcl.fElements[%d]",72+icurrent).Data(),"","goff");
	  Double_t mean=TMath::Mean(entries,treeMeta2->GetV1());
	  treeMap->SetAlias(TString::Format("meanNcl%d",icurrent).Data(), TString::Format("(0+%f)",mean).Data());
	  entries=treeMeta2->Draw(TString::Format("meanNclUsed.fElements[%d]",72+icurrent).Data(),"","goff");
	  mean=TMath::Mean(entries,treeMeta2->GetV1());
	  treeMap->SetAlias(TString::Format("meanNclUsed%d",icurrent).Data(), TString::Format("(0+%f)",mean).Data());	  
	}
	pcstream->GetFile()->cd();
	treeMap->Write();    
      }
    }else{
      fHistos->ls();
      printf("%s does not exist",hisNames[type]);
    }

    delete histo;

  }
  // TPC delta histograms  
  if (type>5) {
    TMatrixD projectionInfo(5,5);
    projectionInfo(0,0)=0;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;
    projectionInfo(1,0)=2;  projectionInfo(1,1)=2;  projectionInfo(1,2)=2; 
    projectionInfo(2,0)=3;  projectionInfo(2,1)=2;  projectionInfo(2,2)=2;
    projectionInfo(3,0)=4;  projectionInfo(3,1)=1;  projectionInfo(3,2)=2;
    projectionInfo(4,0)=1;  projectionInfo(4,1)=1;  projectionInfo(4,2)=2;

    TString axis=(type<9)? "Y":"Z";
    Int_t iter=(type-6)%3;
    TString hName= TString::Format("delta%sIter%d",axis.Data(),iter);
    THnF *histo= (THnF*) fHistos->Get(hName);
    printf("Processing-  TStatToolkit::MakeDistortionMap(4, %s, pcstream, projectionInfo)\n",hName.Data());    
    histo->SetName(hName.Data());
    TStatToolkit::MakeDistortionMap(4, histo, pcstream, projectionInfo); 
    delete histo;
  }



  delete pcstream;
  
  //TStatToolkit::MakeDistortionMap(4, histoITS,    pcstream, projectionInfo); 
  //TStatToolkit::MakeDistortionMap(4, histoITSZ,    pcstream, projectionInfo); 
  //TStatToolkit::MakeDistortionMap(4, histoITSTRDZ, pcstream, projectionInfo); 
  //delete pcstream;
  //
}


void MakeNDFit(const char * inputFile, const char * inputTree, Float_t sector0,  Float_t sector1,  Float_t theta0, Float_t theta1){
  //
  /*
    const char * inputFile="ResidualMapFull_1.root"
    const char * inputTree="deltaRPhiTPCITSTRDDist"
    Float_t sector0=2, sector1 =3;
    Float_t theta0=0, theta1=1;
  */
  TTreeSRedirector * pcstream = new TTreeSRedirector(TString::Format("%sFit_sec%d_%d_theta%d_%d.root",inputTree,Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data(),"recreate");
  TTreeSRedirector * pcstreamFit = new TTreeSRedirector(TString::Format("fitTree_%sFit_sec%d_%d_theta%d_%d.root",inputTree,Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data(),"recreate");

  Int_t runNumber=TString(gSystem->Getenv("runNumber")).Atoi();
  TFile * fdist = TFile::Open(inputFile);
  if (!fdist){
    ::Error("MakeNDFit","Intput file %s not accessible\n",inputFile);
    return;
  }
  treeDist = (TTree*)fdist->Get(inputTree);
  if (!treeDist){
    ::Error("MakeNDFit","Intput tree %s not accessible\n",inputTree);
    return;    
  }
  const Double_t pxmin=8.48499984741210938e+01; //param.GetPadRowRadii(0,0)-param.GetPadPitchLength(0,0)/2
  const Double_t pxmax=2.46600006103515625e+02; //2.46600006103515625e+02param.GetPadRowRadii(36,param.GetNRow(36)-1)+param.GetPadPitchLength(36,param.GetNRow(36)-1)/2.
  Int_t     ndim=4;
  Int_t     nbins[4]= {10,  (sector1-sector0-0.1)*15,        abs(theta1-theta0)*10,        3};  // {radius, phi bin, }
  Double_t  xmin[4] = {pxmin,  sector0+0.05,   theta0,                            -2.0};
  Double_t  xmax[4] = {pxmax, sector1-0.05,   theta1,               2.0};
  //

  THnF* hN= new THnF("exampleFit","exampleFit", ndim, nbins, xmin,xmax);
  treeDist->SetAlias("isOK","rms>0&&entries>50&&abs(rmsG/rms)<2&&abs(mean-meanG)<6");
  TCut cutFit="isOK&&abs(mean-meanG)<6.0";
  TCut cutAcceptFit=TString::Format("sectorCenter>%f&&sectorCenter<%f&&kZCenter>%f&&kZCenter<%f", sector0-0.5,sector1+0.5,theta0,theta1).Data();
  TCut cutAcceptDraw=TString::Format("sectorCenter>%f&&sectorCenter<%f&&kZCenter>%f&&kZCenter<%f", sector0,sector1,theta0,theta1).Data();
  
  AliNDLocalRegression *fitCorrs[6]={0};
  for (Int_t icorr=0; icorr<2; icorr++){
    fitCorrs[icorr]= new  AliNDLocalRegression();
    fitCorrs[icorr]->SetName(TString::Format("%sFit%d_sec%d_%d_theta%d_%d",inputTree,icorr, Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data());  
    Int_t hashIndex=fitCorrs[icorr]->GetVisualCorrectionIndex();
    fitCorrs[icorr]->SetHistogram((THn*)(hN->Clone()));  
    TStopwatch timer;
    fitCorrs[0]->SetStreamer(pcstream);
    if (icorr==0) fitCorrs[icorr]->MakeFit(treeDist,"mean:1", "RCenter:sectorCenter:kZCenter:qptCenter",cutFit+cutAcceptFit,"5:0.05:0.1:3","2:2:2:2",0.0001);
    if (icorr==1) fitCorrs[icorr]->MakeFit(treeDist,"mean:1", "RCenter:sectorCenter:kZCenter:qptCenter",cutFit+cutAcceptFit,"7.:0.075:0.15:3","2:2:2:2",0.0001);
    timer.Print();
    AliNDLocalRegression::AddVisualCorrection(fitCorrs[icorr]);
    treeDist->SetAlias(TString::Format("meanG_Fit%d",icorr).Data(),TString::Format("AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter+0)",hashIndex).Data());
  }
  //
  // Make smoothing at boundaries
  //    
  fitCorrs[2]=(AliNDLocalRegression *)fitCorrs[0]->Clone();
  fitCorrs[3]=(AliNDLocalRegression *)fitCorrs[1]->Clone();
  Int_t nDims=4;
  Int_t indexes[4]={0,1,2,3};
  Double_t relWeight0[12]={1,1,1,   1,1,1, 1,1,1, 1,1,1};
  Double_t relWeightC[12]={0.25,1,1,   0.25,1,1, 0.5,1,1, 0.25,1,1};
  for (Int_t iter=0; iter<3; iter++){
    fitCorrs[2]->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeight0, 0);
    fitCorrs[3]->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeight0, 0);
  }
  fitCorrs[4]=(AliNDLocalRegression *)fitCorrs[2]->Clone();
  fitCorrs[5]=(AliNDLocalRegression *)fitCorrs[3]->Clone();
  fitCorrs[4]->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeightC, 0, kTRUE);
  fitCorrs[5]->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeightC, 0, kTRUE);
  fitCorrs[2]->SetName(TString::Format("%s_Smooth3",fitCorrs[0]->GetName()).Data());
  fitCorrs[3]->SetName(TString::Format("%s_Smooth3",fitCorrs[1]->GetName()).Data());
  fitCorrs[4]->SetName(TString::Format("%s_SmoothConst3",fitCorrs[0]->GetName()).Data());
  fitCorrs[5]->SetName(TString::Format("%s_SmoothConst3",fitCorrs[1]->GetName()).Data());  
  //
  // Make QA and Store fit
  //
  TCanvas *canvasQA = new TCanvas("canvasQA","canvasQA",1200,1000);
  canvasQA->Divide(1,4);
  
  TH1* his=0;
  TFile * fout = pcstream->GetFile();
  pcstream->GetFile()->cd();
  for (Int_t iter=0; iter<6; iter++){
    fitCorrs[iter]->Write();
    fitCorrs[iter]->DumpToTree(4, (*pcstreamFit)<<TString::Format("tree%s", fitCorrs[iter]->GetName()).Data());
  }
  canvasQA->cd(1)->SetLogz();
  treeDist->Draw("meanG_Fit0:meanG>>hisQA2D(200,-1,1,200,-1,1)",cutFit+cutAcceptDraw,"colz");
  treeDist->GetHistogram()->Write("hisQA2D");
  canvasQA->cd(2)->SetLogy();
  //
  treeDist->SetLineColor(1); treeDist->SetMarkerColor(1);
  treeDist->Draw("meanG_Fit0-meanG>>hisQA1D(300,-0.5,0.5)",cutFit+cutAcceptDraw,"");
  his=treeDist->GetHistogram();
  Double_t mean=his->GetMean();
  Double_t rms=his->GetRMS();
  treeDist->GetHistogram()->Write("hisQA1D");
  treeDist->SetLineColor(2); treeDist->SetMarkerColor(2);
  treeDist->Draw("meanG_Fit0-meanG_Fit1>>hisQA1DifFit(300,-0.5,0.5)",cutFit+cutAcceptDraw,"same");
  his=treeDist->GetHistogram();
  Double_t meanDiffFit=his->GetMean();
  Double_t rmsDiffFit=his->GetRMS();
  treeDist->GetHistogram()->Write("hisQA1DifFit");

  treeDist->SetLineColor(1); treeDist->SetMarkerColor(1);
  canvasQA->cd(3)->SetLogy();
  treeDist->Draw("(meanG_Fit0-meanG)/(rms/sqrt(entries))>>hisQAPull(100,-10,10)",cutFit+cutAcceptDraw,"");
  his=treeDist->GetHistogram();
  Double_t meanPull=his->GetMean();
  Double_t rmsPull=his->GetRMS();
  treeDist->GetHistogram()->Write("hisQAPull");
  treeDist->SetLineColor(2); treeDist->SetMarkerColor(2);
  treeDist->Draw("(meanG_Fit0-meanG_Fit1)/(rms/sqrt(entries))>>hisQAPullDifFit(100,-10,10)",cutFit+cutAcceptDraw,"same");
  his=treeDist->GetHistogram();
  Double_t meanPullDiffFit=his->GetMean();
  Double_t rmsPullDiffFit=his->GetRMS();
  treeDist->GetHistogram()->Write("hisQAPullDiffFit");

  canvasQA->SaveAs((TString::Format("%sFit_sec%d_%d_theta%d_%dQA.png",inputTree,Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data()));
  TObjString input=inputTree;
 
  (*pcstream)<<"qa"<<
    "input.="<<&input<<
    "runNumber="<<runNumber<<
    "sec0="<<sector0<<
    "sec1="<<sector1<<
    "theta0="<<theta0<<
    "theta1="<<theta1<<
    "mean="<<mean<<
    "rms="<<rms<<
    "meanPull="<<meanPull<<
    "rmsPull="<<rmsPull<<
    //
    "meanDiffFit="<<meanDiffFit<<
    "rmsDiffFit="<<rmsDiffFit<<
    "meanPullDiffFit="<<meanPullDiffFit<<
    "rmsPullDiffFit="<<rmsPullDiffFit<<
    "\n";
  
  delete pcstream;
  delete pcstreamFit;

}

void loadDistortionTree(){
  //
  //
  //
  TTree * tree =  0;
  TObjArray* array = TString(gSystem->GetFromPipe("cat map.list")).Tokenize("\n");

  for (Int_t i=0; i<array->GetEntries(); i++){
    printf("%s\n",array->At(i)->GetName());
    TString runName=gSystem->BaseName(gSystem->DirName(array->At(i)->GetName()));
    if (TString(array->At(i)->GetName()).Contains("_0.root")){
      tree = AliTPCcalibAlignInterpolation::AddFriendDistortionTree(tree,array->At(i)->GetName(),"deltaRPhiTPCITSTRDDist",runName+"TRDY");
    }
    if (TString(array->At(i)->GetName()).Contains("_1.root")){
      tree = AliTPCcalibAlignInterpolation::AddFriendDistortionTree(tree,array->At(i)->GetName(),"deltaRPhiTPCITSTOFDist",runName+"TOFY");
    }
    if (TString(array->At(i)->GetName()).Contains("_2.root")){
      tree = AliTPCcalibAlignInterpolation::AddFriendDistortionTree(tree,array->At(i)->GetName(),"deltaZTPCITSTOFDist",runName+"TOFZ");
    }
  }
}

void makeCurrentTrend(){
  //
  //
  //
  TCut cutFit = "refCurrent!=0&&id<=240220";
  TChain * chain=  AliXRDPROOFtoolkit::MakeChainRandom("timeInfo.list","sumaryTime",0,1000);
  TStopwatch timer;
  TTree *tree = chain->CopyTree("1");
  timer.Print();
  tree->SetAlias("refCurrent","0.5*(grNcl79.fY+grNcl80.fY)");
  tree->SetAlias("time","grNcl80.fX");
  //
  // make current normalization aliases
  //
  for (Int_t iSec=0; iSec<=80; iSec++){
    Int_t entries = tree->Draw(TString::Format("grNcl%d.fY/grNcl%d.fY", iSec,72+iSec/18).Data(),cutFit,"goff");
    Double_t median=TMath::Median(entries,tree->GetV1());
    printf("isec=%d\tnorm=%f\n",iSec,median);
    tree->SetAlias(TString::Format("normNcl%d",iSec).Data(),  TString::Format("(grNcl%d.fY)/(%f+0)",iSec,median).Data());
    tree->SetAlias(TString::Format("normNclRatio%d",iSec).Data(),  TString::Format("(grNcl%d.fY/grNcl%d.fY)/(%f+0)",iSec,iSec,median).Data());   
  }

  TCanvas *pcanvasCurrent[4]={0};
  for (Int_t iType=0; iType<4; iType++){
    pcanvasCurrent[iType] = new TCanvas(TString::Format("canvasCurrentROC%d",iType).Data(),TString::Format("canvasCurrentROC%d",iType).Data(),1100,1100);
    pcanvasCurrent[iType] = new TCanvas(TString::Format("canvasCurrentROC%d",iType).Data(),TString::Format("canvasCurrentROC%d",iType).Data(),1100,1100);
    pcanvasCurrent[iType]->Divide(3,6,0);
    for (Int_t iSec=0; iSec<=18; iSec++){
      pcanvasCurrent[iType]->cd(iSec+1);
      tree->Draw(TString::Format("normNclRatio%d:time", iSec+18*iType).Data(),cutFit,"");      
    }
    
  }
}


Bool_t  FitDrift(TTree * treeIn, Int_t timeStampMin, Int_t timeStampMax, Int_t npointMax, const char delta){
  //
  // Make drift velocity+radial distortion rough linear fit 
  // see  fitDriftString string - radial components idnependent on A side and C side
  //
  //   1.) Make fit 
  //   2.) Make QA plots
  //   3.) Store fit value in the tree user info
  //
  /*
    const char * chInput="/hera/alice/alien/calib/alice/data/2015/LHC15o/000245231/cpass0_pass1/ResidualMerge/001/ResidualTrees.root";    
    TFile * finput = TFile::Open(chInput);
    TTree * treeIn = (TTree*)finput->Get("delta");
    npointMax=10000;
    delta="tof1.fElements";

  */
  const Int_t kMinPoints=100;
  treeIn->SetAlias("drift","(1-abs(vecZ.fElements)/250.)");
  treeIn->SetAlias("normR","(vecR.fElements/250.)");
  treeIn->SetAlias("normGY","(sin(vecPhi.fElements)*vecR.fElements/250.)");
  treeIn->SetAlias("sideA","(int(vecSec.fElements)%36)<18");
  treeIn->SetAlias("sideC","(int(vecSec.fElements)%36)>=18");
  treeIn->SetAlias("delta","tof1.fElements");
  //
  // 1.) Make fit
  //
  Int_t  npointsMax=2000;
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;
  TString fitDriftString="";
  fitDriftString+="drift*(2*sideA-1)++";
  fitDriftString+="drift*normGY*(2*sideA-1)++";
  fitDriftString+="sideA*normR++";         // elements to correct for the radial distortions
  fitDriftString+="sideA*normR*drift++";   // 
  fitDriftString+="sideC*normR++";
  fitDriftString+="sideC*normR*drift++";
  // 
  TCut cutFit="Iteration$<npValid&&vecZ.fElements!=0&&abs(delta)<20";  
  TString *strDrift = TStatToolkit::FitPlane(treeIn,"delta:1", fitDriftString.Data(),cutFit, chi2,npoints,param,covar,-1,0, npointsMax/5);
  if (npoints==kMinPoints) return kFALSE;
  
  //
  treeIn->SetAlias("fitDrift",strDrift->Data());
  strDrift = TStatToolkit::FitPlane(treeIn,"delta:1", fitDriftString.Data(),cutFit+"abs(fitDrift-delta)<5", chi2,npoints,param,covar,-1,0, npointsMax/2);
  treeIn->SetAlias("fitDrift",strDrift->Data());
  strDrift = TStatToolkit::FitPlane(treeIn,"delta:1", fitDriftString.Data(),cutFit+"abs(fitDrift-delta)<5", chi2,npoints,param,covar,-1,0, npointsMax);
  treeIn->SetAlias("fitDrift",strDrift->Data());
  TObjArray* tokArr = strDrift->Tokenize("++");
  tokArr->Print();
  //
  // 2.) drift QA plots
  //
  TH1 * hisQA[4]={0};
  TObjArray fitArray(3);
  treeIn->Draw("delta-fitDrift:vecR.fElements:vecPhi.fElements>>hisSecA(90,-3.14,3.14,20,85,245)",cutFit+"sideA&&abs(fitDrift-delta)<2.5","profcolzgoff",4*npointsMax);
  hisQA[0]=(TH1*)treeIn->GetHistogram()->Clone();
  treeIn->Draw("delta-fitDrift:vecR.fElements:vecPhi.fElements>>hisSecC(90,-3.14,3.14,20,85,245)",cutFit+"sideC&&abs(fitDrift-delta)<2.5","profcolzgoff",4*npointsMax);
  hisQA[1]=(TH1*)treeIn->GetHistogram()->Clone();
  //
  treeIn->Draw("delta-fitDrift:vecZ.fElements>>hisZ(44,-220,220,100,-3,3)",cutFit+"","colzgoff",4*npointsMax);
  hisQA[2]=(TH1*)treeIn->GetHistogram()->Clone();
  ((TH2*)hisQA[2])->FitSlicesY(0,0,-1,0,"QNR",&fitArray);
  //
  //
  //
  TCanvas *canvasZFit = new TCanvas("canvasZFit","canvasZFit",1000,800);
  canvasZFit->Divide(1,3);
  canvasZFit->cd(1);
  hisQA[0]->Draw("colz");
  canvasZFit->cd(2);
  hisQA[1]->Draw("colz");
  canvasZFit->cd(3);
  hisQA[2]->Draw("colz");
  fitArray.At(1)->Draw("same");
  canvasZFit->SaveAs("ATO-108_canvasZFit.png");
  //
  //gStyle->SetOptStat();
  TCanvas *canvasZFitModel = new TCanvas("canvasZFitModel","canvasZFitModel",1000,800);
  canvasZFitModel->Divide(3,1);
  canvasZFitModel->cd(1);
  treeIn->Draw("delta:vecZ.fElements",cutFit+"","",npointsMax);
  canvasZFitModel->cd(2);
  treeIn->Draw("delta:fitDrift",cutFit+"","",npointsMax);
  canvasZFitModel->cd(3);
  treeIn->Draw("fitDrift*(2*sideA-1):vecZ.fElements:vecR.fElements",cutFit+"","profcolz",4*npointsMax);
  canvasZFitModel->SaveAs("ATO-108_canvasZFitModel.png");
  //
  // 
  //

}
