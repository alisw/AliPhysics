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
void MakeEventStatInfo(const char * inputList="cat residual.list", Int_t timeInterval=300, Int_t id=0, Int_t skip=1);


void AliTPCcalibAlignInterpolationMacro(Int_t action,Int_t param0=0, Int_t param1=0){
  if (action==0)   AliTPCcalibAlignInterpolationLoop(100000,1000000000);
  if (action==1)    {
    Int_t startTime=TString(gSystem->Getenv("mapStartTime")).Atoi();
    Int_t stopTime =TString(gSystem->Getenv("mapStopTime" )).Atoi();
    AliTPCcalibAlignInterpolation::FillHistogramsFromChain("residual.list", 6,10,startTime, stopTime, param1,param0);
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
    MakeEventStatInfo("cat residual.list",300,param0,1);
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
    if (histo) {
      //TStatToolkit::MakeDistortionMap(4, histo, pcstream, projectionInfo); 
      TStatToolkit::MakeDistortionMapFast(histo, pcstream, projectionInfo,kFALSE); 
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
    const char * inputFile="ResidualMapFull_0.root"
    const char * inputTree="deltaRPhiTPCITSTRDDist"
    Float_t sector0=0,Int_t(sector1)=1;
    Float_t theta0=0, theta1=1;
  */
  Int_t runNumber=TString(gSystem->Getenv("runNumber")).Atoi();
  TFile * fdist = TFile::Open(inputFile);
  treeDist = (TTree*)fdist->Get(inputTree);
  Int_t     ndim=4;
  Int_t     nbins[4]= {10,  (sector1-sector0)*10,        abs(theta1-theta0)*10,        3};  // {radius, phi bin, }
  Double_t  xmin[4] = {84,  sector0,   theta0,                            -2.0};
  Double_t  xmax[4] = {245, sector1,   theta1,               2.0};
  //

  THnF* hN= new THnF("exampleFit","exampleFit", ndim, nbins, xmin,xmax);
  treeDist->SetAlias("isOK","rms>0&&entries>50&&abs(rmsG/rms)<2&&abs(mean-meanG)<6");
  TCut cutFit="isOK&&abs(mean-meanG)<6.0";
  TCut cutAcceptFit=TString::Format("sectorCenter>%f&&sectorCenter<%f&&kZCenter>%f&&kZCenter<%f", sector0-0.5,sector1+0.5,theta0,theta1).Data();
  TCut cutAcceptDraw=TString::Format("sectorCenter>%f&&sectorCenter<%f&&kZCenter>%f&&kZCenter<%f", sector0,sector1,theta0,theta1).Data();
  
  AliNDLocalRegression *fitCorrs[2]={0};
  for (Int_t icorr=0; icorr<2; icorr++){
    fitCorrs[icorr]= new  AliNDLocalRegression();
    fitCorrs[icorr]->SetName(TString::Format("%sFit%d_sec%d_%d_theta%d_%d",inputTree,icorr, Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data());  
    Int_t hashIndex=fitCorrs[icorr]->GetVisualCorrectionIndex();
    fitCorrs[icorr]->SetHistogram((THn*)(hN->Clone()));  
    TStopwatch timer;
    if (icorr==0) fitCorrs[icorr]->MakeFit(treeDist,"mean:1", "RCenter:sectorCenter:kZCenter:qptCenter",cutFit+cutAcceptFit,"5:0.1:0.1:3","2:2:2:2",0.00001);
    if (icorr==1) fitCorrs[icorr]->MakeFit(treeDist,"mean:1", "RCenter:sectorCenter:kZCenter:qptCenter",cutFit+cutAcceptFit,"7.:0.15:0.15:3","2:2:2:2",0.00001);
    timer.Print();
    AliNDLocalRegression::AddVisualCorrection(fitCorrs[icorr]);
    treeDist->SetAlias(TString::Format("meanG_Fit%d",icorr).Data(),TString::Format("AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter+0)",hashIndex).Data());
  }
  
  //
  // Make QA and Store fit
  //
  TCanvas *canvasQA = new TCanvas("canvasQA","canvasQA",1200,1000);
  canvasQA->Divide(1,4);
  TTreeSRedirector * pcstream = new TTreeSRedirector(TString::Format("%sFit_sec%d_%d_theta%d_%d.root",inputTree,Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data(),"recreate");
  
  TH1* his=0;
  TFile * fout = pcstream->GetFile();
  fitCorrs[0]->Write();
  fitCorrs[1]->Write();
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

void MakeEventStatInfo(const char * inputList, Int_t timeInterval, Int_t id, Int_t skip){
  //
  /// Code to query statistical event information from the ResidualTrees.root file 
  /// output written to file residualInfo.root
  ///   \param const char * inputList - ascii file with input list
  ///   \param Int_t timeInterval     - length of time interval (beginning of time intervals rounded)
  ///   \param id                     - additional ID added to the tree
  ///   \param skip                   - parameter skip file
  /// Algorithm:
  ///   1.) Cache information per files - beginTime and endTime for file
  ///   2.) Cache information per time interval

  /*
    run=240204;
    GetResidualStatInfo("cat residual.list",300,run,1);
  */
  TObjArray *array = TString(gSystem->GetFromPipe(TString::Format("%s",inputList).Data())).Tokenize("\n");
  Int_t nFiles=array->GetEntries();
  if (nFiles<=0) {
    ::Error("GetResidualStatInfo. Wrong input list",inputList);
    return;
  }
  TStopwatch timer;
  //
  // 1.) Cache information per files - beginTime and endTime for file
  //
  TStopwatch timer1;
  TTreeSRedirector * pcstream = new TTreeSRedirector("residualInfo.root", "recreate");
  for (Int_t iFile=0; iFile<nFiles; iFile+=skip){
    timer.Start();
    printf("%d\t%s\n",iFile,array->At(iFile)->GetName());
    TFile * f = TFile::Open(array->At(iFile)->GetName());
    if (f==NULL) continue;
    TTree * treeInfo = (TTree*)f->Get("eventInfo");
    if (treeInfo==NULL) continue;
    Int_t entriesInfo=treeInfo->GetEntries();
    Int_t entries=treeInfo->Draw("B1","1","goff");
    Double_t maxTime=TMath::MaxElement(entries,treeInfo->GetV1());
    Double_t minTime=TMath::MinElement(entries,treeInfo->GetV1());
    Double_t meanTime=TMath::Mean(entries,treeInfo->GetV1());
    TObjString fname(array->At(iFile)->GetName());
    (*pcstream)<<"summary1"<<
      "iFile="<<iFile<<
      "fname.="<<&fname<<
      "events="<<entriesInfo<<
      "minTime="<<minTime<<
      "maxTime="<<maxTime<<
      "meanTime="<<meanTime<<
      "\n";
    timer.Print();
  }
  delete pcstream;
  ::Info("GetResidualStatInfo","Total time");
  timer1.Print();
  //
  // 2.) Cache information per time interval
  //
  TStopwatch timer2;
  pcstream = new TTreeSRedirector("residualInfo.root", "update");
  TTree * treeSummary1=(TTree*)(pcstream->GetFile()->Get("summary1"));
  Int_t entries = treeSummary1->Draw("minTime","1","goff");
  Long64_t minTime = TMath::MinElement(entries, treeSummary1->GetV1());
  entries = treeSummary1->Draw("maxTime","1","goff");
  Long64_t maxTime = TMath::MaxElement(entries, treeSummary1->GetV1());
  minTime=timeInterval*(minTime/timeInterval);
  maxTime=timeInterval*(1+(maxTime/timeInterval));
  Int_t nIntervals=(maxTime-minTime)/timeInterval;
  Int_t nIntervalsQA=(maxTime-minTime)/15;
  //
  TH1F  * hisEvent= new TH1F("hisEvent","hisEvent",nIntervalsQA,minTime,maxTime);
  const Int_t nSec=81; // 72 sector +5 sumarry info+ 4 medians
  TProfile * profArrayNcl[nSec]={0};
  TProfile * profArrayNclUsed[nSec]={0};
  TGraphErrors * grArrayNcl[nSec]={0};
  TGraphErrors * grArrayNclUsed[nSec]={0};
  TProfile * profArrayITSNcl[3]={0};
  TGraphErrors * grArrayITSNcl[3]={0};
  
  for (Int_t isec=0; isec<nSec; isec++){
    profArrayNcl[isec]=new TProfile(TString::Format("TPCnclSec%d",isec).Data(), TString::Format("TPCnclSec%d",isec).Data(), nIntervalsQA,minTime,maxTime);
    profArrayNclUsed[isec]=new TProfile(TString::Format("TPCnclUsedSec%d",isec).Data(), TString::Format("TPCnclUsedSec%d",isec).Data(), nIntervalsQA,minTime,maxTime);
  }
   for (Int_t iits=0; iits<3; iits++){
    profArrayITSNcl[iits]=new TProfile(TString::Format("ITSnclSec%d",iits).Data(), TString::Format("ITSnclSec%d",iits).Data(), nIntervalsQA,minTime,maxTime);    
  }

  TVectorF *vecNClTPC=0;
  TVectorF *vecNClTPCused=0;
  Int_t nITS[3]={0};
  Int_t timeStamp=0;
  for (Int_t iFile=0; iFile<nFiles; iFile+=skip){
    timer.Start();
    printf("%d\t%s\n",iFile,array->At(iFile)->GetName());    
    TFile * f = TFile::Open(array->At(iFile)->GetName());
    if (f==NULL) continue;
    TTree * treeInfo = (TTree*)f->Get("eventInfo"); 
    if (treeInfo==NULL) continue;
    treeInfo->SetBranchAddress("vecNClTPC.",&vecNClTPC);
    treeInfo->SetBranchAddress("vecNClTPCused.",&vecNClTPCused);
    treeInfo->SetBranchAddress("nSPD",&nITS[0]);
    treeInfo->SetBranchAddress("nSDD",&nITS[1]);
    treeInfo->SetBranchAddress("nSSD",&nITS[2]);
    Bool_t hasTimeStamp=(treeInfo->GetBranch("timeStamp")!=NULL);
    if (hasTimeStamp) treeInfo->SetBranchAddress("timeStamp",&timeStamp);
    if (!hasTimeStamp) ((TBranch*)(treeInfo->GetListOfBranches()->At(1)))->SetAddress(&timeStamp);
    Int_t treeEntries=treeInfo->GetEntries();
    for (Int_t iEntry=0; iEntry<treeEntries; iEntry++){
      treeInfo->GetEntry(iEntry);
      hisEvent->Fill(timeStamp);
      for (Int_t isec=0; isec<72; isec++){
	profArrayNcl[isec]->Fill(timeStamp, (*vecNClTPC)[isec]);
	profArrayNclUsed[isec]->Fill(timeStamp, (*vecNClTPC)[isec]);
	if (isec<36){
	  if (isec<18) 	profArrayNcl[72]->Fill(timeStamp, (*vecNClTPC)[isec]);
	  if (isec>=18) profArrayNcl[73]->Fill(timeStamp, (*vecNClTPC)[isec]);
	  if (isec<18) 	profArrayNclUsed[72]->Fill(timeStamp, (*vecNClTPCused)[isec]);
	  if (isec>=18) profArrayNclUsed[73]->Fill(timeStamp, (*vecNClTPCused)[isec]);
	}else{
	  if ((isec%36)<18)  profArrayNcl[74]->Fill(timeStamp, (*vecNClTPC)[isec]);
	  if ((isec%36)>=18) profArrayNcl[75]->Fill(timeStamp, (*vecNClTPC)[isec]);
	  if ((isec%36)<18)  profArrayNclUsed[74]->Fill(timeStamp, (*vecNClTPCused)[isec]);
	  if ((isec%36)>=18) profArrayNclUsed[75]->Fill(timeStamp, (*vecNClTPCused)[isec]);
	}
	profArrayNcl[76]->Fill(timeStamp, (*vecNClTPC)[isec]);
	profArrayNclUsed[76]->Fill(timeStamp, (*vecNClTPCused)[isec]);
      }
      profArrayNcl[77]->Fill(timeStamp, TMath::Median(18, &((vecNClTPC->GetMatrixArray())[0])));
      profArrayNcl[78]->Fill(timeStamp, TMath::Median(18, &((vecNClTPC->GetMatrixArray())[18])));
      profArrayNcl[79]->Fill(timeStamp, TMath::Median(18, &((vecNClTPC->GetMatrixArray())[36])));
      profArrayNcl[80]->Fill(timeStamp, TMath::Median(18, &((vecNClTPC->GetMatrixArray())[54])));
      //
      profArrayNclUsed[77]->Fill(timeStamp, TMath::Median(18, &((vecNClTPCused->GetMatrixArray())[0])));
      profArrayNclUsed[78]->Fill(timeStamp, TMath::Median(18, &((vecNClTPCused->GetMatrixArray())[18])));
      profArrayNclUsed[79]->Fill(timeStamp, TMath::Median(18, &((vecNClTPCused->GetMatrixArray())[36])));
      profArrayNclUsed[80]->Fill(timeStamp, TMath::Median(18, &((vecNClTPCused->GetMatrixArray())[54])));
      for (Int_t iits=0; iits<3; iits++){
	profArrayITSNcl[iits]->Fill(timeStamp,nITS[iits]);
      }
    }
    timer.Print();
  }
  timer2.Print();
  TGraphErrors grEvent(hisEvent);
  (*pcstream)<<"sumaryTime"<<
    "id="<<id<<
    "grEvent.="<<&grEvent;
  for (Int_t isec=0; isec<nSec; isec++){
    grArrayNcl[isec] = new TGraphErrors((profArrayNcl[isec]));
    grArrayNclUsed[isec] = new TGraphErrors((profArrayNclUsed[isec]));
    (*pcstream)<<"sumaryTime"<<
      TString::Format("grNcl%d.=",isec).Data()<<grArrayNcl[isec]<<
      TString::Format("grNclUsed%d.=",isec).Data()<<grArrayNclUsed[isec];
  }
  for (Int_t iits=0; iits<3; iits++){
    grArrayITSNcl[iits] = new TGraphErrors((profArrayITSNcl[iits]));
    (*pcstream)<<"sumaryTime"<<
      TString::Format("grITSNcl%d.=",iits).Data()<<grArrayITSNcl[iits];
  }
  
  
  (*pcstream)<<"sumaryTime"<<"\n";
  for (Int_t isec=0; isec<nSec; isec++){
    delete 	profArrayNcl[isec];
    delete 	profArrayNclUsed[isec];
    delete 	grArrayNcl[isec];
    delete 	grArrayNclUsed[isec];
  }
  delete hisEvent;
  delete pcstream;

  printf("StatInfo.minTime\t%d\n",minTime);
  printf("StatInfo.maxTime\t%d\n",maxTime);
  delete array;
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
    tree->SetAlias(TString::Format("normNclRatio%d",iSec).Data(),  TString::Format("(grNcl%d.fY/grNcl%d.fY)/(%f+0)",iSec,median).Data());   
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
