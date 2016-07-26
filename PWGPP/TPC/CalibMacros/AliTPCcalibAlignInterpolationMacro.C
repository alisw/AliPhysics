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
#include "TStyle.h"
#include "TLegend.h"
//
#include "AliSysInfo.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliExternalTrackParam.h"
#include "AliTPCcalibAlignInterpolation.h"
#include "AliNDLocalRegression.h"
#include "AliMathBase.h"



TTree * treeDist =0;   // global tree used to check the content of the fit

void AliTPCcalibAlignInterpolationLoop(Int_t nchunks, Int_t neventsMax);
void  CreateDistortionMapsFromFile(Int_t type, const char * inputFile, const char *outputFile, Int_t delta=1);


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

    AliTPCcalibAlignInterpolation::MakeNDFit(inputFile, inputTree, sec0,sec1,theta0,theta1);		      
  } 

  if (action==5) {
    Int_t runNumber=TString(gSystem->Getenv("runNumber")).Atoi();
     AliTPCcalibAlignInterpolation::MakeEventStatInfo("cat residual.list",300,param0,1);
  }

  if (action==6){  //fit drift velocity
    Int_t deltaT=TString(gSystem->Getenv("driftDeltaT" )).Atoi();
    Int_t sigmaT=TString(gSystem->Getenv("driftSigmaT" )).Atoi();
    if (deltaT<=0 || sigmaT<=0){
      ::Error("AliTPCcalibAlignInterpolation::FitDrift FAILED ","Invalid parameter value for the deltaT %d and sigmaT", deltaT, sigmaT);
      return;
    }
    ::Info("AliTPCcalibAlignInterpolation::FitDrift","Begin");
    AliTPCcalibAlignInterpolation::FitDrift(deltaT, sigmaT);
    ::Info("AliTPCcalibAlignInterpolation::FitDrift","End");
  }
 
}

void AliTPCcalibAlignInterpolationLoop(Int_t nchunks, Int_t neventsMax){
  //
  // Invocation of the AliTPCcalibAlignInterpolation from the macro for test purposes (instead of train usage)
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




void DrawFailedFitsExample(const char * fname, const char * treeName, Int_t firstEntry=0, Int_t entriesDraw=9){
  //
  // Function to draw the TH1 where the Gaussian fit failed
  //   this function was important at the moment when we considered gaussian fit to be default position estimator
  // Later on decided LTM estimator used instead
  //
  //
  //
  /* 
    fname="/hera/alice/miranov/alice-tpc-notes/SpaceChargeDistortion/data/ATO-108/alice/data/2015/LHC15o0901/000245231/Time1448613300/his1/ResidualMapFull_1.root";
    treeName="deltaRPhiTPCITSTRDDistDump";
    firstEntry=0;
    entriesDraw=9;
   */
  TFile *ff = TFile::Open(fname);
  gStyle->SetOptFit(1);
  TH1* his=0;
  TTree * tree = (TTree*)ff->Get(treeName);
  tree->SetBranchAddress("hDump.",&his);
  Int_t entriesSel=tree->Draw("Entry$:binMedian:rms","entries>150&&isFitValid==0","goff");
  TVectorD vecSel(entriesSel, tree->GetV1());
  TVectorD vecMed(entriesSel, tree->GetV2());
  TCanvas * canvas = new TCanvas("canvasGausFailed","canvasGausFailed",900,900);
  entriesDraw=TMath::Nint(TMath::Sqrt(entriesDraw));
  canvas->Divide(entriesDraw, entriesDraw);
  entriesDraw*=entriesDraw;
  //
  TF1 fgaus0("fgaus0","gaus",-10,10);
  TF1 fgaus1("fgaus1","gaus",-10,10);
  fgaus0.SetLineColor(2);
  fgaus1.SetLineColor(4);
  for (Int_t i=firstEntry; i<firstEntry+entriesDraw && i<entriesSel; i++){
    canvas->cd(i+1);
    tree->GetEntry(TMath::Nint(vecSel[i]));
    TH1 *his2=(TH1*)his->Clone();
    his2->GetXaxis()->UnZoom();
    fgaus0.SetParameters(his2->GetEntries()/(2.5*his2->GetRMS()/his2->GetBinWidth(50)), vecMed[i], his2->GetRMS());
    fgaus1.SetParameters(his2->GetEntries()/(2.5*his2->GetRMS()/his2->GetBinWidth(50)), vecMed[i], his2->GetRMS());
    fgaus0.SetRange(his2->GetMean()-4*his2->GetRMS(),his2->GetMean()+4*his2->GetRMS());
    his2->Fit(&fgaus0,"rl");
    fgaus1.Draw("same");
  }
  canvas->SaveAs("ATO-108_FaledFitExample.png");

}  
  
void DumpDiffEstimators(const char *fname0,const char *fname1, const char *treeName0, const char *treeName1 ){
  /*
    fname0="/hera/alice/miranov/alice-tpc-notes/SpaceChargeDistortion/data/ATO-108/alice/data/2015/LHC15o0901/000245231/Time1448613300/his1/ResidualMapFull_1.root";
    fname1="/hera/alice/miranov/alice-tpc-notes/SpaceChargeDistortion/data/ATO-108/alice/data/2015/LHC15o0901/000245231/Time1448614500/his1/ResidualMapFull_1.root";
    treeName0="deltaRPhiTPCITSTRDDist";
    treeName1="deltaRPhiTPCITSTRDDist";
  */
  TFile *f0= TFile::Open(fname0);
  TFile *f1= TFile::Open(fname1);
  TTree *tree0 = (TTree*)f0->Get(treeName0);
  TTree *tree1 = (TTree*)f1->Get(treeName1);
  tree0->AddFriend(tree1,"T1");
  TCut cutDraw="abs(qptCenter)<0.1&&min(entries,T1.entries)>20";
  TObjArray arrayFit(3);
  
  TH1*hisRMS[4]={0};
  TH1*hisMean[4]={0};
  TH1*his2D[4]={0};
  const char * chLegend[4]={"mean all","mean restricted","gaus fit", "binMedian"};
  Int_t entries= tree0->Draw("1/sqrt(entries)",cutDraw,"goff");
  Double_t xmin=TMath::KOrdStat(entries, tree0->GetV1(), Int_t(entries*0.02));
  Double_t xmax=TMath::KOrdStat(entries, tree0->GetV1(), Int_t(entries*0.98));
  for (Int_t i=0; i<4; i++){
    TString query="";
    if (i==0) query=TString::Format("(mean0-T1.mean0):1/sqrt(entries)>>hisMean0(10,%f,%f,100,-0.5,0.5)",xmin,xmax);
    if (i==1) query=TString::Format("(mean-T1.mean):1/sqrt(entries)>>hisMean(10,%f,%f,100,-0.5,0.5)",xmin,xmax);
    if (i==2) query=TString::Format("(meanG-T1.meanG):1/sqrt(entries)>>hisMeanG(10,%f,%f,100,-0.5,0.5)",xmin,xmax);
    if (i==3) query=TString::Format("(binMedian-T1.binMedian):1/sqrt(entries)>>hisMedian(10,%f,%f,100,-0.5,0.5)",xmin,xmax);
    tree0->Draw(query.Data(), cutDraw,"goff");
    ((TH2*)tree0->GetHistogram())->FitSlicesY(0,0,-1,0,"QNR",&arrayFit);
    hisMean[i]=(TH1*)arrayFit.At(1);
    hisRMS[i]=(TH1*)arrayFit.At(2);
    his2D[i]=tree0->GetHistogram();
    his2D[i]->SetTitle(chLegend[i]);
    his2D[i]->GetXaxis()->SetTitle("1/#sqrt(N)");
    his2D[i]->GetYaxis()->SetTitle("#sigma (cm)");
  }
  //
  TCanvas *canvasDrawEst= new TCanvas("canvasDrawEst","canvasDrawEst",1000,800);
  canvasDrawEst->Divide(2,3);
  for (Int_t i=0; i<4; i++){
    canvasDrawEst->cd(i+1);
    his2D[i]->Draw("colz");
    hisMean[i]->Draw("same");    
  }
  canvasDrawEst->cd(5);
  TLegend *legend= new TLegend(0.1,0.1,0.7,0.7,"Estimator resoutution as function of the entries");
  legend->SetNColumns(2);
  for (Int_t i=0; i<4; i++){
    hisRMS[i]->SetMarkerStyle(21+i);
    hisRMS[i]->SetMarkerColor(1+i);
    hisRMS[0]->SetMinimum(0);
    if (i==0) hisRMS[i]->Draw();
    hisRMS[i]->Draw("same");
    legend->AddEntry(hisRMS[i], chLegend[i]);
  }
  canvasDrawEst->cd(6);
  legend->Draw();
  canvasDrawEst->SaveAs("ATO-108_ComparisonOfEstimators_000245231_Time1448613300_his1.png");


}

Bool_t LTMHisto(TH1 *his, TVectorD &params , Float_t fraction){
  // 
  // LTM : Trimmed mean on histogram - Modified version for binned data
  //       Using binned information     
  //
  // Robust statistic to estimate properties of the distribution
  // To handle binning error special treatment
  // for definition of unbinned data see:
  //     http://en.wikipedia.org/w/index.php?title=Trimmed_estimator&oldid=582847999
  //
  // Function parameters:
  //     his1D   - input histogram
  //     params  - vector with parameters
  //             - 0 - area
  //             - 1 - mean
  //             - 2 - rms 
  //
  Int_t nbins    = his->GetNbinsX();
  Int_t nentries = (Int_t)his->Integral(0,nbins);
  const Double_t kEpsilon=0.0000000001;

  if (nentries<=0) return 0;
  if (fraction>1) fraction=0.9999;
  if (fraction<0) return 0;
  
  TVectorD vec(nentries);
  Int_t all=0;
  for (Int_t ibin=1; ibin<=nbins; ibin++){
    Int_t ncont=his->GetBinContent(ibin);
    printf("%d\t%d\n",ibin,ncont);
    Double_t lowEdge= his->GetBinLowEdge(ibin);
    Double_t width  = his->GetBinWidth(ibin);
    for (Int_t icont=0; icont<ncont; icont++){
      Double_t x=lowEdge+icont*(width/ncont);
      vec[all++]=x;
    }
  }
  Double_t mean, rms;
  if (fraction*all<3) return 0;
  AliMathBase::EvaluateUni(all, vec.GetMatrixArray(), mean,rms, fraction*all);
  params[0]=all;
  params[1]=mean;
  params[2]=rms;
  params[3]=rms/TMath::Sqrt(fraction*all);
  return kTRUE;
}

