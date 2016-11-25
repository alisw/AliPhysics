/*
  Steering macro for the Detector/<action> trending

  //
  // code to check consitency
  //
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/install/include -I$ALICE_ROOT/STEER\
   -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW  -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/TPCBase  -I$ALICE_ROOT/TPC/TPCRec -I$ALICE_ROOT/TPC/TPCCalib -I$ALICE_PHYSICS/../src/PWGPP/TPC/");
  .L $ALICE_PHYSICS/../src/PWGPP/QA/scripts/qaTrending.C+
  period="LHC15o";
  pass="cpass1_pass1";
  detector="T0";
  referenceDet="ITS;TPC;TRD;TOF";
  
  treeQADet=InitTrees(detector,referenceDet);
  treeQADet->SetAlias("tagID","run");  
  InitSummaryTrending(treeQADet);

*/

#include "TTree.h"
#include "TString.h"
#include "TStatToolkit.h"
#include "AliExternalInfo.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TLegend.h"


TTree * InitTrees(const char * detector,   const char * referenceDet);    // load QA trees and Sttering tree <det>   QA.<Det>, Logbook, Lo
void    InitSummaryTrending(TTree * tree);                     // can be generic ?      
void    SetDrawStyle();                                        // setting Default draw style
TObjArray * GetTexDescription(TLatex *latex);                  // should go to base code

//
// Macro global variables:
//
char * year=0, *period=0, *pass=0;
TTree * treeQADet=0;
std::map<std::string,TTree*> qaMap;
TCanvas *canvasQA=0;          // central canvas for QA plots
TObjArray * descriptionQA=0;  // descriptionQA of process  (time/aliroot/root/pariod/pass/user)
TObjArray*   oaMultGr=0;      // status bar multigraph
TString      statusDescription[3];   // status bar description


//
//
// Drawing configuration can be changed
//
// colors & markers:
Int_t colPosA=kGreen+2; //kRed;
Int_t colPosC=kMagenta+1; //kAzure-4;
Int_t colNegA=kRed; //kRed;//kOrange;
Int_t colNegC=kAzure-4; //kAzure-4;//kGreen;
Int_t colSum=kBlack;
Int_t marA1=21;
Int_t marA2=24;
Int_t marC1=20;
Int_t marC2=25;
Int_t marSum=34;//full cross
Int_t marCorr=31;//snowflake
Float_t markerSize=0.9;
const char * defColors="1;2;4;3;856;616";
const char * defMarkers="21;22;25;26;27;28";
const char * bandColors="1;400;400;632;632;416;416"; // kBlack;kYellow;kYellow;kRed; kRed; kGreeen; kGreen
// shifting of graphs within one run for better visibility:
Float_t sh_gr0=-0.3;
Float_t sh_gr1=-0.1;
Float_t sh_gr2=+0.1;
Float_t sh_gr3=+0.3;
// maximal canvas width
Int_t kMaxCanvasWidth=2500;
// status bar height
Float_t statPadHeight=0.30; //fraction of canvas height (was 0.25 until Aug 2014)
Float_t statPadBotMar=0.40; //bottom margin of pad for run numbers
Int_t   padRightMargin=0.07;


void qaTrending(){
}

TTree * InitTrees(const char * detector,  const char *referenceDet){
  //
  // Init tree for given period
  //  all  trees stored in qaMap
  //  FriendTree added to the master treeQADet
  // Currentrly in the code we assume ID="run";
  // Bug in root tree ?  - in case more than one friend tree set - indeces looks corrupted 
  //    0.) QA tree per detector Raw+MC (if exist)
  //    1.) QA trees per refernce detectors specified by string refDet e.g "TPC;TRD;TOF"
  //    2.) Logbook tree per run
  //    3.) Logbook tree per run/detector
  //    3.) RCT table
  //    4.) CPass table 
  // 
  //  tree is created with addition of friend trees which can be used in queries
  //  queries as analog to the SQL statement

  /* 
    period="LHC15o"; pass="cpass1_pass1"
  */
  ::Info("qaTrending::InitTrees::Begin","Detector %s, RefDet=%s",detector, referenceDet);
  AliExternalInfo info;
  Int_t treeCounter=0;
  // Load trees
  TObjArray * detArray=TString(referenceDet).Tokenize(";");
  Int_t nrefDets=detArray->GetEntries();
  TVectorF runCounter(5+nrefDets*2);  // <QADet>, <Logbook>, <Logbook.Det>, <MonAlisa>, <CPass>, <QA.RefDet>xnrefDets, <Logbook.RefDet>xnrefDets
  //
  ::Info("qaTrending::InitTrees::End","Laoding trees");
  treeQADet = info.GetTreeDataQA(detector,period, pass);
  if (!treeQADet){
    ::Error("qaTrending.InitTrees", "Input QA tree %s not available", detector); 
    return 0;
  }
  runCounter[treeCounter++]=treeQADet->Draw("run","1","goff");
  qaMap[TString::Format("QA.%s",detector).Data()]=treeQADet;
  //
  qaMap["Logbook"]=info.GetTree("Logbook",period,"");
  qaMap["Logbook"]->AddFriend(treeQADet,"QADet");
  treeQADet->AddFriend(qaMap["Logbook"],"Logbook");
  runCounter[treeCounter++]=treeQADet->Draw("run","1","goff");
  //
  qaMap["MonALISA.RCT"]=info.GetTree("MonALISA.RCT",period,pass);
  qaMap["MonALISA.RCT"]->AddFriend(treeQADet,"QADet");
  treeQADet->AddFriend(qaMap["MonALISA.RCT"],"MonALISA.RCT");
  runCounter[treeCounter++]=treeQADet->Draw("run","1","goff");
  //
  TTree *treeLogbookDetector =info.GetTree("Logbook.detector",period,"");
  if (treeLogbookDetector){
    if (detArray->GetEntries()>0){
      for (Int_t idet=0; idet<detArray->GetEntries(); idet++){
	// Load Logbook.RefDet
	const char *detName=detArray->At(idet)->GetName();	
	TTree * treeLog =treeLogbookDetector->CopyTree(TString::Format("detector==\"%s\"",detName).Data());
	if (treeLog->GetEntries()<=0){
	  ::Error("qaTrending.InitTrees","Missing Tree Logbook. %s - check the syntax",detName);
	}else{
	  treeLog->BuildIndex("run");
	  qaMap[TString::Format("Logbook.%s",detName).Data()]=treeLog;
	  treeLog->AddFriend(treeQADet, "QADet");
	  treeQADet->AddFriend(treeLog,  TString::Format("Logbook.%s",detName).Data());
	  runCounter[treeCounter++]=treeQADet->Draw("run","1","goff");
	}
	// Load QA.RefDet
	TTree * treeQARefDet = info.GetTreeDataQA(detName,period, pass);
	if (treeQARefDet){
	  qaMap[TString::Format("QA.%s",detName).Data()]=treeQARefDet;
	  treeQARefDet->AddFriend(treeQADet, "QADet");
	  treeQADet->AddFriend(treeQARefDet,  TString::Format("QA.%s",detName).Data());
	  runCounter[treeCounter++]=treeQADet->Draw("run","1","goff");
	}else{
	  ::Error("qaTrending.InitTrees","Missing Tree QA.%s - check the syntax",detName);
	}
      }
    }
  }
  //
  // Check consistency of data
  // 
  ::Info("qaTrending::InitTrees::End","Checking trees");
  TList *arrFriends = treeQADet->GetListOfFriends();
  for (Int_t ifriend=0; ifriend<arrFriends->GetEntries(); ifriend++){
    Int_t entries = treeQADet->Draw(TString::Format("run-%s.run", arrFriends->At(ifriend)->GetName()).Data(),"1","goff");
    Double_t mean=0;
    if (entries>0) {
      mean=TMath::Mean(entries, treeQADet->GetV1());
    }
    if (mean==0){
      ::Info("qaTrending::InitTrees", "Friend:\t%s\t%d\t%f", arrFriends->At(ifriend)->GetName(), entries,mean);
    }else{
      ::Error("qaTrending::InitTrees", "Friend:\t%s\t%d\t%f", arrFriends->At(ifriend)->GetName(), entries,mean);
    }
  }
  delete detArray;
  ::Info("qaTrending::InitTrees::End","Detector %s, RefDet=%s",detector, referenceDet);

  return treeQADet;  
}


void SetDrawStyle(){
  //
  // Standard drawing style
  //
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetLabelSize(0.04,"x");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

 


TObjArray * GetTexDescription(TLatex *latex){
  //
  //
  //
  TObjArray * description= new TObjArray();
  TString sTimestamp  = TString::Format("Creation time:%s",gSystem->GetFromPipe("date").Data());
  TString sUser=TString::Format("User:%s",gSystem->GetFromPipe("echo $USER").Data());  
  TString sPeriod=TString::Format("period:%s",period);  
  TString sPass=TString::Format("pass:%s",pass);  
  TString sAlirootVer;
  TString sAliphysicsVer;
  //TString runStat=
  if (gSystem->GetFromPipe("echo $ALICE_VER") == "master" || gSystem->GetFromPipe("echo $ALICE_VER") == ""){
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("wdir=`pwd`; cd $ALICE_ROOT/../src; git describe; cd $wdir;");
  }else {
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("echo $ALIROOT_VERSION");
  }
  if (gSystem->GetFromPipe("echo $ALIPHYSICS_VER") == "master" || gSystem->GetFromPipe("echo $ALIPHYSICS_VER") == ""){
    sAliphysicsVer = "AliPhysics: " + gSystem->GetFromPipe("wdir=`pwd`; cd $ALICE_PHYSICS/../src; git describe; cd $wdir;");
  }else {
    sAliphysicsVer = "AliPhysics: " + gSystem->GetFromPipe("echo $ALIPHYSICS_VERSION");
  }
  //
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY(),sTimestamp.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-1.1*latex->GetTextSize(),sUser.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-2.2*latex->GetTextSize(),sPeriod.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-3.3*latex->GetTextSize(),sPass.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-4.4*latex->GetTextSize(),sAlirootVer.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-5.5*latex->GetTextSize(),sAliphysicsVer.Data()));
  return description;
}


void InitSummaryTrending(TTree * tree){
  //
  // Init drawing for the <detector> QA
  // Detector specific qaConfig() has to be called before invoking this function
  //    0.) Make descriptor 
  //    1.) Make default canvas - addopt canvas width to the number of entries to draw
  //    3.) compute detector  status graphs

  //
  //   0.) Make descriptor 
  //
  TLatex *latex= new TLatex;
  latex->SetX(0.11);
  latex->SetY(0.8);
  latex->SetTextSize(0.03);
  descriptionQA = GetTexDescription(latex);
  //
  //   1.) Make default canvas - addopt canvas width to the number of entries to draw
  //
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"resolution:tagID","");
  Int_t numberOfTags = gr->GetN();
  cout<<"number of graph entries: "<<numberOfTags<<endl;
  double SpaceForLegend = 0.;
  Int_t canvas_width  = SpaceForLegend + (numberOfTags+5)*30;
  Int_t canvas_height = 600; 
  if ( canvas_width>kMaxCanvasWidth)   canvas_width=kMaxCanvasWidth;
  canvasQA = new TCanvas("canvasQA","canvasQA",canvas_width,canvas_height);
  canvasQA->SetGrid(3);
  canvasQA->cd();
  gPad->SetTicks(1,2);
  // canvasQA->SetRightMargin(SpaceForLegend/canvas_width);
  double leftlegend  = 1 - 180./canvasQA->GetWw();
  double rightlegend = 1 - 10./canvasQA->GetWw();
  //
  // 2.) process config file qaConfig.C to initialize status aliases (outliers etc.), status bar criteria, status lines, ...
  //
  TString sStatusbarVars  = statusDescription[0];
  TString sStatusbarNames = statusDescription[1];
  TString sCriteria       = statusDescription[2];
  cout << "sStatusbarVars = " << sStatusbarVars.Data() << endl;
  cout << "sCriteria      = " << sCriteria.Data() << endl;
  //
  // 3.) compute detector status graphs
  //
  TObjArray* oaStatusbarVars = sStatusbarVars.Tokenize(";");
  TObjArray* oaStatusbarNames = sStatusbarNames.Tokenize(";");
  oaMultGr = new TObjArray();
  int igr=0;
  for (Int_t vari=oaStatusbarVars->GetEntriesFast()-1; vari>=0; vari--){ // invert the order of the status graphs
    TString sVar = Form("%s:tagID", oaStatusbarVars->At(vari)->GetName()); //e.g. -> dcar:run
    oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, sVar.Data(),  "", sCriteria.Data(), igr) );
    TString sYtitle = oaStatusbarNames->At(vari)->GetName(); // set better name for y axis of statuspad
    ((TMultiGraph*) oaMultGr->At(igr))->SetTitle(sYtitle.Data());
    igr++;
  }
}
