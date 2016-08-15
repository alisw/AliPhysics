/*
 
 First version -  code to make a T0 QA trending 
   to run the code in shell: 
   aliroot -b -q  $ALICE_PHYSICS/../src/PWGPP/QA/scripts/qaTrending.C+ '$ALICE_PHYSICS/../src/PWGPP/QA/detectorQAscripts/qaTrendingT0.C("LHC15o","cpass1_pass1")'

   //
   To develop the code run following sequence and after add new metadata to  qaConfig() or new plots to   DrawSummaryTrending()

 gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/install/include -I$ALICE_ROOT/STEER\
  -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW  -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/TPCBase  -I$ALICE_ROOT/TPC/TPCRec -I$ALICE_ROOT/TPC/TPCCalib -I$ALICE_PHYSICS/../src/PWGPP/TPC/");
 
  .L $ALICE_PHYSICS/../src/PWGPP/QA/scripts/qaTrending.C+
  .L $ALICE_PHYSICS/../src/PWGPP/QA/detectorQAscripts/qaTrendingT0.C
  //

  period="LHC16f";
  pass="cpass1_pass1";
  detector="T0";
  referenceDet="T0;ITS;TPC;TRD;TOF";  
  tree=InitTrees(detector,referenceDet);
  
  qaConfig(tree,statusDescription);
  tree->SetAlias("tagID","run");
  InitSummaryTrending(tree);
  DrawSummaryTrending(tree);



*/



Int_t   qaConfig(TTree* tree, TString* returnStrings);         // actual detector configuration for alarms
void    DrawSummaryTrending(TTree * tree);                     // actual detector drawing code 
//      



void qaTrendingT0(const char *pperiod, const char *ppass){
  //
  //
  //
  period=pperiod;
  pass=ppass;
  detector="T0";
  referenceDet="T0;ITS;TPC;TRD;TOF;EVS";  
  tree=InitTrees(detector,referenceDet);
  qaConfig(tree,statusDescription);
  tree->SetAlias("tagID","run");
  InitSummaryTrending(tree);
  DrawSummaryTrending(tree);
}


Int_t qaConfig(TTree* tree, TString* returnStrings)
{
  // Configure alarms and status bars
  //    0.) Define standard cut to define statisical properties
  //    1.) Define standard aliases Outlier/Warning/PhysAcc for standard variables 
  //    2.) Define custom alaises 
  //    4.) Define status bar layout
  //    5.) Define metadata describing variables
  //
  // 0. Define standard cut  
  //
  logbookConfig(tree);
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
  TStatToolkit::AddMetadata(tree,"tzeroOrAPlusOrC.Title","#Delta_{T}");
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
    //
    TStatToolkit::AddMetadata(tree,TString::Format("timeDelayPMT%d.AxisTitle",pmt).Data(),"#Delta_{T} (a.u.)");
    TStatToolkit::AddMetadata(tree,TString::Format("timeDelayPMT%d.Title",pmt).Data(),"#Delta_{T}");
    TStatToolkit::AddMetadata(tree,TString::Format("timeDelayPMT%d.Legend",pmt).Data(),TString::Format("PMT %d",pmt).Data());
  }
}




void DrawSummaryTrending(TTree * tree){
  //
  // Make default trending plots
  //
  TMultiGraph * multiGraph=0;
  TMultiGraph * multiLine=0;
  SetDrawStyle();
  //
  // Plot 0. Draw, print basic time info
  //
  if (qaMap["Logbook.T0"]){
    canvasQA->Clear(); 
    TLegend *legendTime=new TLegend(0.11,0.03,0.3,0.3,"T0 default QA - Run time information ");  legendTime->SetNColumns(2);   legendTime->SetBorderSize(0);
    multiGraph=TStatToolkit::MakeMultGraph(tree, "T0 <T>","runDuration;pauseDuration:tagID","",defMarkers,defColors,kTRUE,0.8,9, legendTime);
    multiGraph->Draw();
    legendTime->Draw();
    descriptionQA->DrawClone();
    TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
    TStatToolkit::DrawStatusGraphs(oaMultGr);
    canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin); canvasQA->Draw("ap");
    canvasQA->SaveAs("Logbook.T0.Timing_Vs_TagID.png");
  }

  if (qaMap["Logbook.T0"]){
    canvasQA->Clear(); 
    TLegend *legendStat=new TLegend(0.11,0.03,0.3,0.3,"T0 default QA - Statistic information "); legendStat->SetNColumns(2);   legendStat->SetBorderSize(0);
    multiGraph=TStatToolkit::MakeMultGraph(tree, "T0 <T>","totalEvents;totalEventsPhysics;Logbook.T0.eventCountPhysics;Logbook.TPC.eventCountPhysics:tagID","",defMarkers,defColors,kTRUE,0.8,9, legendStat);
    multiGraph->Draw();
    descriptionQA->DrawClone();
    legendStat->Draw();    
    TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
    TStatToolkit::DrawStatusGraphs(oaMultGr);
    canvasQA->cd(1)->SetRightMargin(padRightMargin);     canvasQA->cd(2)->SetRightMargin(padRightMargin);     canvasQA->Draw("ap");
    canvasQA->SaveAs("Logbook.T0.Statistic_vs_TagID.png");
  }
  //
  // Plot 1.)
  //  
  /****** Resolution  vs ID (run number or period, ...) ******/
  canvasQA->Clear(); 
  TLegend *legend1=new TLegend(0.11,0.03,0.3,0.3,"T0 default QA - Time resolution ");
  multiGraph=TStatToolkit::MakeMultGraph(tree, "T0 <T>","resolution:tagID:5","","21","1;",kTRUE,0.8,9, legend1);
  multiLine =TStatToolkit::MakeMultGraph(tree, "T0 <T>","resolution_RobustMean;resolution_WarningMin;resolution_WarningMax;resolution_OutlierMin;resolution_OutlierMax:tagID", "","0;0;0;0;0;0;0;0",bandColors,kTRUE,0,0,0);
  multiGraph->SetMinimum(-100);
  multiGraph->SetMaximum(200);
  multiGraph->Draw();
  legend1->Draw();
  multiLine->Draw();
  descriptionQA->DrawClone();
  TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin); canvasQA->Draw("ap");
  canvasQA->SaveAs("QA.T0.Resolution_Vs_TagID.png");
  canvasQA->Clear();
  //
  // Plot 2.)
  //  
  /****** Delta  vs ID (run number or period, ...) ******/
  canvasQA->Clear();
  TLegend *legend2=new TLegend(0.11,0.03,0.3,0.3,"T0 default QA - #DeltaT "); legend2->SetNColumns(2);   legend2->SetBorderSize(0);
  multiGraph=TStatToolkit::MakeMultGraph(tree, "T0 <T>","tzeroOrAPlusOrC;tzeroOrA;tzeroOrC:tagID:10;10;10","","21;22;25;26;27;28","1;2;4;3;856;616",kTRUE,0.8,9, legend2);
  multiLine =TStatToolkit::MakeMultGraph(tree, "T0 <T>","tzeroOrAPlusOrC_RobustMean;tzeroOrAPlusOrC_WarningMin;tzeroOrAPlusOrC_WarningMax;tzeroOrAPlusOrC_OutlierMin;tzeroOrAPlusOrC_OutlierMax:tagID", "","0;0;0;0;0;0;0;0",bandColors,kTRUE,0,0,0);
  multiGraph->Draw();
  multiLine->Draw();
  legend2->Draw();
  descriptionQA->DrawClone();  
  TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin);  canvasQA->Draw();
  canvasQA->SaveAs("QA.T0.DeltaT_Vs_TagID.png");
  canvasQA->Clear();
  //
  // Plot 3.)
  //  
  /****** Delta PMT  vs ID (run number or period, ...) ******/
  for (Int_t ipmt=1; ipmt<24; ipmt+=6){
    canvasQA->Clear();
    TLegend *legendPMT=new TLegend(0.11,0.03,0.3,0.3,"T0 default QA - #DeltaT PMT ");legendPMT->SetNColumns(2);  legendPMT->SetBorderSize(0);
    TString drawStr="";
    for (Int_t idr=0; idr<6; idr++) drawStr+=TString::Format("timeDelayPMT%d;",ipmt+idr);
    drawStr+=":tagID";
    multiGraph=TStatToolkit::MakeMultGraph(tree, "T0 <T>",drawStr,"",defMarkers,defColors,kTRUE,0.8,3, legendPMT);
    multiGraph->Draw();
    legendPMT->Draw();
    TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
    TStatToolkit::DrawStatusGraphs(oaMultGr);
    descriptionQA->DrawClone();  
    canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin);  canvasQA->Draw();
    canvasQA->SaveAs(TString::Format("QA.T0.DeltaPMT%d_%d_Vs_TagID.png",ipmt,ipmt+5).Data());
    canvasQA->Clear();
  }

}
