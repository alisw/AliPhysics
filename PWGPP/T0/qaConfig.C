/*
  //
  // macro -qaConfig to define standard statistical properties for summary T0 QA variables
  //
  // Define properties:
  //   1.) <varname>_Mean, LTS, median for variable of interes  
  //   2.) outlier (usually 6 sigma) and warning (usualy 3 sigma) bands
  //   **  <varname>_Warning, <varname>_WarningMin,   <varname>_OutlierMax
  //   **  <varname>_Outlier, <varname>_OutlierMin,  <varname>_OutlierMax, 
  //   **  <varname>_PhysAcc, <varname>_PhysAccMin,  <varname>_PhysAccMax
  // Define status bar configuartion
  //  
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/install/include -I$ALICE_ROOT/STEER\
   -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW  -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/TPCBase  -I$ALICE_ROOT/TPC/TPCRec -I$ALICE_ROOT/TPC/TPCCalib -I$ALICE_PHYSICS/../src/PWGPP/TPC/");
  .L $ALICE_PHYSICS/../src/PWGPP/T0/qaConfig.C+
  period="LHC15o";
  pass="cpass1_pass1";
  tree=InitTrees();
  trendingTree=tree;
  tree->SetAlias("tagID","run");

  TString returnString[3];
  qaConfig(trendingTree, returnString);
  DrawTrending(tree);

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


//
// Setup configuration:
//
char * year=0, *period=0, *pass=0;

//
// Drawing configuration:
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
Int_t padRightMargin=0.07;


TTree * trendingTree=0;
TTree * rctTree=0, *cpassTree=0;


TTree * InitTrees(){
  //
  // Example tree for testing
  //
  AliExternalInfo info;
  TTree * tree = info.GetTreeDataQA("T0",period, pass);
  return tree;  
}


void SetDrawStyle(){
  //
  //
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetLabelSize(0.04,"x");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}


Int_t ComputeRange(TTree* tree, const char* varname, Float_t &plotmean, Float_t &plotoutlier)
{
  //the function computes useful numbers for plot ranges from the outlier criteria
  plotmean    = (Float_t) TFormula("fcn", tree->GetAlias(Form("%s_RobustMean",varname))).Eval(0);
  plotoutlier = (Float_t) TFormula("fcn", tree->GetAlias(Form("%s_OutlierMax",varname))).Eval(0) - plotmean;
  return 1;
}

Int_t PlotStatusLines(TTree * tree, const char * expr, const char * cut)
{
  //the function plots status lines
  const char* alias = "varname_OutlierMin:varname_OutlierMax:varname_WarningMin:varname_WarningMax:varname_PhysAccMin:varname_PhysAccMax:varname_RobustMean";
  TMultiGraph* mgStatusLines = TStatToolkit::MakeStatusLines(tree,expr,cut,alias);

  if (mgStatusLines) mgStatusLines->Draw("l");
  else { cout << " no mgStatusLines available!" << endl; return 0; }

  return 1;
}
 

TTree * GetCPassTree(const char * period, const * char pass){
  //
  // try to find production information about pass OCDB export
  //
  TTree * treeProd=0;
  AliExternalInfo info;
  treeProdArray = info.GetTreeCPass();
  treeProdArray->Scan("ID:Description",TString::Format("strstr(Description,\"%s\")&&strstr(Description,\"%s\")",period,pass).Data(),"col=10:100");
  // check all candidata production and select one which exports OCDB
  Int_t entries= treeProdArray->Draw("ID:Description",TString::Format("strstr(Description,\"%s\")&&strstr(Description,\"%s\")",period,pass).Data(),"col=10:100");  
  for (Int_t ientry=0; ientry<entries; ientry++){
    TTree * treeProd0 = info.GetTreeProdCycleByID(TString::Format("%.0f",treeProdArray->GetV1()[ientry]).Data());
    treeProd0->Draw("1","strstr(outputdir,\"OCDB\")==1");
    if (treeProd0->Draw("1","strstr(outputdir,\"OCDB\")==1")==0) {
      delete treeProd0;
      continue;
    }
    treeProd=treeProd0;
  }

}

Bool_t MyRegExp(const char * chinput,Int_t p1){
  //
  //
  //
  printf("%s\n",chinput);
  return p1;
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

Int_t qaConfig(TTree* tree, TString* returnStrings)
{
  // Configure alarms and status bars
  //    0.) Define standard cut to define statisical properties
  //    1.) Define standard aliases Outlier/Warning/PhysAcc for standard variables 
  //    2.) Define custom alaises 
  //    4.) Define status bar layout
  //    5.) Define metadata descibing variables
  //
  // 0. Define standard cut  
  //
  tree->SetAlias("statisticOK", "abs(resolution)<100");
  Float_t entryFrac=0.8, nsigmaOutlier=6., nsigmaWarning=3., epsilon=1.0e-6;  
  //
  // 1. Define aliases Outlier/Warning/PhysAcc for combined variables 
  //
  TString sTrendVars="resolution;tzeroOrAPlusOrC;tzeroOrA;tzeroOrC";
  
  TObjArray* oaTrendVars = sTrendVars.Tokenize(",;");
  for (Int_t vari=0; vari<oaTrendVars->GetEntriesFast(); vari++)
  {
    TString sVar( oaTrendVars->At(vari)->GetName() );
    // outliers, warnings and robust mean are set for all variables identically.
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_RobustMean:(MeanEF+0):%f", entryFrac));  // robust mean
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_RMS:(RMSEF+0):%f", entryFrac));          // robust rms
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_OutlierMin:(MeanEF-%f*RMSEF-%f):%f", nsigmaOutlier, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_OutlierMax:(MeanEF+%f*RMSEF+%f):%f", nsigmaOutlier, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_WarningMin:(MeanEF-%f*RMSEF-%f):%f", nsigmaWarning, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_WarningMax:(MeanEF+%f*RMSEF+%f):%f", nsigmaWarning, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*RMSEF-%f):%f", nsigmaWarning, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*RMSEF+%f):%f", nsigmaWarning, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_Outlier:(varname>varname_OutlierMax||varname<varname_OutlierMin)"));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_Warning:(varname>varname_WarningMax||varname<varname_WarningMin)"));
    // here we set dummy physics accepatable range - this value should be physics drivven - 
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAcc:(varname>varname_WarningMax||varname<varname_WarningMin)"));
  }      
  //
  // 3.) custom variables 
  //
  TStatToolkit::SetStatusAlias(tree,"tzeroOrAPlusOrC" ,    "statisticOK", Form("varname_Outlier:(varname>varname_OutlierMax||varname<varname_OutlierMin)|| abs(varname)>40"));
  TStatToolkit::SetStatusAlias(tree,"tzeroOrAPlusOrC" ,    "statisticOK", Form("varname_Warning:(varname>varname_WarningMax||varname<varname_WarningMin)|| abs(varname)>20"));

  TStatToolkit::SetStatusAlias(tree, "tzeroOrAPlusOrC",    "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f+0)", 30.)); // 30 ps cut
  TStatToolkit::SetStatusAlias(tree, "tzeroOrAPlusOrC",    "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f+0)", 30.)); // 30 ps cut
  //
  // 4.) Define status bar layout
  //
  TString sStatusbarVars ("resolution;tzeroOrAPlusOrC");  // variable names  - used for query 
  TString sStatusbarNames("T0 resolution;#Delta_{T0} (A+C);");  // 
  TString sCriteria("(1):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)"); // or to just show vetos: (varname_PhysAcc&&varname_Warning)
  returnStrings[0] = sStatusbarVars;
  returnStrings[1] = sStatusbarNames;
  returnStrings[2] = sCriteria;
  //
  //    5.) Define metadata describing variables
  TStatToolkit::AddMetadata(tree,"resolution.Title","Time resolution");
  TStatToolkit::AddMetadata(tree,"resolution.AxisTitle","Time resolution (ps)");
  //
  TStatToolkit::AddMetadata(tree,"tzeroOrAPlusOrC.Title","Delta T");
  TStatToolkit::AddMetadata(tree,"tzeroOrAPlusOrC.AxisTitle","#DeltaT (ps)");
  TStatToolkit::AddMetadata(tree,"tzeroOrAPlusOrC.Legend","(A+C)");
  TStatToolkit::AddMetadata(tree,"tzeroOrA.Title","#Delta_{T}");
  TStatToolkit::AddMetadata(tree,"tzeroOrA.AxisTitle","#Delta_{T} (ps)");
  TStatToolkit::AddMetadata(tree,"tzeroOrA.Legend","A side");
  TStatToolkit::AddMetadata(tree,"tzeroOrA.Title","#Delta_{T}");
  TStatToolkit::AddMetadata(tree,"tzeroOrA.AxisTitle","#Delta_{T} (ps)");
  TStatToolkit::AddMetadata(tree,"tzeroOrC.Legend","C side");
  //
  //
  for (Int_t pmt=1; pmt<=24; pmt++){
    TStatToolkit::AddMetadata(tree,TString::Format("amplPMT%d.AxisTitle",pmt).Data(),"<Q>_{PMT} (a.u.)");
    TStatToolkit::AddMetadata(tree,TString::Format("amplPMT%d.Title",pmt).Data(),"<Q>_{PMT}");
    TStatToolkit::AddMetadata(tree,TString::Format("amplPMT%d.Legend",pmt).Data(),TString::Format("PMT %d",pmt).Data());
  }

}


void DrawTrending(TTree * tree){
  //
  // Draw trending
  // 
  tree->SetAlias("tagID","run");
  TMultiGraph * multiGraph=0;
  TMultiGraph * multiLine=0;
  TLatex *latex= new TLatex;
  latex->SetX(0.11);
  latex->SetY(0.8);
  latex->SetTextSize(0.03);
  TObjArray * description = GetTexDescription(latex);
  //
  // 1.) Make default canvas - addopt canvas width to the number of entries to draw
  //
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"resolution:tagID","");
  Int_t numberOfTags = gr->GetN();
  cout<<"number of graph entries: "<<numberOfTags<<endl;
  double SpaceForLegend = 0.;
  Int_t canvas_width  = SpaceForLegend + (numberOfTags+5)*30;
  Int_t canvas_height = 600; 
  if ( canvas_width>kMaxCanvasWidth)   canvas_width=kMaxCanvasWidth;
  TCanvas *canvasQA = new TCanvas("canvasQA","canvasQA",canvas_width,canvas_height);
  canvasQA->SetGrid(3);
  canvasQA->cd();
  gPad->SetTicks(1,2);
  // canvasQA->SetRightMargin(SpaceForLegend/canvas_width);
  double leftlegend  = 1 - 180./canvasQA->GetWw();
  double rightlegend = 1 - 10./canvasQA->GetWw();
  //
  // 1.) process config file qaConfig.C to initialize status aliases (outliers etc.), status bar criteria, status lines, ...
  //
  TString returnStrings[3];
  qaConfig(tree, returnStrings);
  TString sStatusbarVars  = returnStrings[0];
  TString sStatusbarNames = returnStrings[1];
  TString sCriteria       = returnStrings[2];
  cout << "sStatusbarVars = " << sStatusbarVars.Data() << endl;
  cout << "sCriteria      = " << sCriteria.Data() << endl;
  // compute TPC status graphs
  TObjArray* oaStatusbarVars = sStatusbarVars.Tokenize(";");
  TObjArray* oaStatusbarNames = sStatusbarNames.Tokenize(";");
  TObjArray* oaMultGr = new TObjArray();
  int igr=0;
  for (Int_t vari=oaStatusbarVars->GetEntriesFast()-1; vari>=0; vari--){ // invert the order of the status graphs
    TString sVar = Form("%s:tagID", oaStatusbarVars->At(vari)->GetName()); //e.g. -> dcar:run
    oaMultGr->Add( TStatToolkit::MakeStatusMultGr(tree, sVar.Data(),  "", sCriteria.Data(), igr) );
    TString sYtitle = oaStatusbarNames->At(vari)->GetName(); // set better name for y axis of statuspad
    ((TMultiGraph*) oaMultGr->At(igr))->SetTitle(sYtitle.Data());
    igr++;
  }
  //
  // 
  //
  
  Float_t plotmean, plotoutlier;
  //
  // Plot 1.)
  //  
  /****** Resolution  vs ID (run number or period, ...) ******/
  canvasQA->Clear(); 
  TLegend *legend1=new TLegend(0.11,0.11,0.3,0.3,"T0 default QA - Time resolution ");
  multiGraph=TStatToolkit::MakeMultGraph(tree, "T0 <T>","resolution:tagID:5","","21","1;",kTRUE,0.8,9, legend1);
  multiLine =TStatToolkit::MakeMultGraph(tree, "T0 <T>","resolution_RobustMean;resolution_WarningMin;resolution_WarningMax;resolution_OutlierMin;resolution_OutlierMax:tagID", "","0;0;0;0;0;0;0;0",bandColors,kTRUE,0,0,0);
  multiGraph->SetMinimum(-100);
  multiGraph->SetMaximum(200);
  multiGraph->Draw();
  legend1->Draw();
  multiLine->Draw();
  TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin); canvasQA->Draw("ap");
  description->DrawClone();
  canvasQA->SaveAs("resolution_vs_TagID.png");
  canvasQA->Clear();

  //
  // Plot 2.)
  //  
  /****** Delta  vs ID (run number or period, ...) ******/
  canvasQA->Clear();
  TLegend *legend2=new TLegend(0.11,0.11,0.3,0.3,"T0 default QA - #DeltaT ");
  legend2->SetNColumns(2);
  legend2->SetBorderSize(0);
  multiGraph=TStatToolkit::MakeMultGraph(tree, "T0 <T>","tzeroOrAPlusOrC;tzeroOrA;tzeroOrC:tagID:10;10;10","","21;22;25;26;27;28","1;2;4;3;856;616",kTRUE,0.8,9, legend2);
  multiLine =TStatToolkit::MakeMultGraph(tree, "T0 <T>","tzeroOrAPlusOrC_RobustMean;tzeroOrAPlusOrC_WarningMin;tzeroOrAPlusOrC_WarningMax;tzeroOrAPlusOrC_OutlierMin;tzeroOrAPlusOrC_OutlierMax:tagID", "","0;0;0;0;0;0;0;0",bandColors,kTRUE,0,0,0);
  multiGraph->Draw();
  multiLine->Draw();
  legend2->Draw();
  TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin);  canvasQA->Draw();
  description->DrawClone();

  canvasQA->SaveAs("T0QADeltaTTrending.png");
  canvasQA->Clear();

  


}
