/*
   First version -  code to make a TRD QA trending 
   to run the code in shell: 
   aliroot -b -q  $ALICE_PHYSICS/../src/PWGPP/QA/scripts/qaTrending.C+ '$ALICE_PHYSICS/../src/PWGPP/QA/detectorQAscripts/qaTrendingTRD.C("LHC15o","cpass1_pass1")'

   //
   To develop the code run following sequence and after add new metadata to  qaConfig() or new plots to   DrawSummaryTrending()
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/install/include -I$ALICE_ROOT/STEER\
  -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW  -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/TPCBase  -I$ALICE_ROOT/TPC/TPCRec -I$ALICE_ROOT/TPC/TPCCalib -I$ALICE_PHYSICS/../src/PWGPP/TPC/"); 
  .L $ALICE_PHYSICS/../src/PWGPP/QA/scripts/qaTrending.C+
  .L $ALICE_PHYSICS/../src/PWGPP/QA/detectorQAscripts/qaTrendingTRD.C
  //
  period="LHC15o";
  pass="cpass1_pass1";
  detector="TRD";
  referenceDet="T0;ITS;TPC;TRD;TOF;EVS";  
  tree=InitTrees(detector,referenceDet);
  qaConfig(tree,statusDescription);
  tree->SetAlias("tagID","run");
  InitSummaryTrending(tree);
  DrawSummaryTrending(tree);
  //
*/

void qaTrendingTRD(const char *pperiod, const char *ppass){
  //
  //
  //
  period=pperiod;
  pass=ppass;
  detector="TRD";
  referenceDet="T0;ITS;TPC;TRD;TOF;EVS";  
  tree=InitTrees(detector,referenceDet);
  qaConfig(tree,statusDescription);
  tree->SetAlias("tagID","run");
  InitSummaryTrending(tree);
  DrawSummaryTrending(tree);
}


Int_t   qaConfig(TTree* tree, TString* returnStrings);         // actual detector configuration for alarms
void    DrawSummaryTrending(TTree * tree);                     // actual detector drawing code 
//      



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
  tree->SetAlias("statisticOK", "TPCTRDmatchEffPosAll>0"); // decalre as runs to extract properties runs with TPC
  Float_t entryFrac=0.8, nsigmaOutlier=6., nsigmaWarning=3., epsilon=1.0e-6;  
  //
  // 1. Define aliases Outlier/Warning/PhysAcc for combined variables 
  //
  TString sTrendVars="AvNclsPerTrackAll;AvTRDtrkltsPerTrackAll;TPCTRDmatchEffPosAll;TPCTRDmatchEffNegAll;meanT0;meanVdrift;TRDTOFmatchEffPosAll;TRDTOFmatchEffNegAll";
  //
  //
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
    // here we set dummy physics accepatable range - this value should be physics driven - 
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAcc:(varname>varname_WarningMax||varname<varname_WarningMin)"));
  }      
  //
  // 3.) custom variables 
  //
  tree->SetAlias("TPCTRDmatchEffAll_Outlier","(TPCTRDmatchEffPosAll_Outlier||TPCTRDmatchEffNegAll_Outlier)");
  tree->SetAlias("TPCTRDmatchEffAll_Warning","(TPCTRDmatchEffPosAll_Warning||TPCTRDmatchEffNegAll_Warning)");
  tree->SetAlias("TPCTRDmatchEffAll_PhysAcc","(TPCTRDmatchEffPosAll_PhysAcc||TPCTRDmatchEffNegAll_PhysAcc)");
  //
  tree->SetAlias("TRDTOFmatchEffAll_Outlier","(TRDTOFmatchEffPosAll_Outlier||TRDTOFmatchEffNegAll_Outlier)");
  tree->SetAlias("TRDTOFmatchEffAll_Warning","(TRDTOFmatchEffPosAll_Warning||TRDTOFmatchEffNegAll_Warning)");
  tree->SetAlias("TRDTOFmatchEffAll_PhysAcc","(TRDTOFmatchEffPosAll_PhysAcc||TRDTOFmatchEffNegAll_PhysAcc)");


  //
  // 4.) Define status bar layout
  //
  TString sStatusbarVars ("AvNclsPerTrackAll;AvTRDtrkltsPerTrackAll;TPCTRDmatchEffAll;TRDTOFmatchEffAll;meanT0;meanVdrift");       // variable names  - used for query 
  TString sStatusbarNames("N_{cls};N_{tracklet};#epsilon_{TPC->TRD};#epsilon_{TRD->TOF} ;meanT0;meanVdrift");  // 
  TString sCriteria("(1):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)"); // or to just show vetos: (varname_PhysAcc&&varname_Warning)
  returnStrings[0] = sStatusbarVars;
  returnStrings[1] = sStatusbarNames;
  returnStrings[2] = sCriteria;
  //
  //    5.) Define metadata describing variables
  TStatToolkit::AddMetadata(tree,"TPCTRDmatchEffPosAll.Title","#epsilon_{TPC->TRD}");
  TStatToolkit::AddMetadata(tree,"TPCTRDmatchEffPosAll.AxisTitle","#epsilon_{TPC->TRD}");
  TStatToolkit::AddMetadata(tree,"TPCTRDmatchEffPosAll.Legend","(Q>0)");
  TStatToolkit::AddMetadata(tree,"TPCTRDmatchEffNegAll.Title","#epsilon_{TPC->TRD}");
  TStatToolkit::AddMetadata(tree,"TPCTRDmatchEffNegAll.AxisTitle","#epsilon_{TPC->TRD}");
  TStatToolkit::AddMetadata(tree,"TPCTRDmatchEffNegAll.Legend","(Q<0)");
  TStatToolkit::AddMetadata(tree,"TRDTOFmatchEffPosAll.Title","#epsilon_{TPC->TRD}");
  TStatToolkit::AddMetadata(tree,"TRDTOFmatchEffPosAll.AxisTitle","#epsilon_{TPC->TRD}");
  TStatToolkit::AddMetadata(tree,"TRDTOFmatchEffPosAll.Legend","(Q>0)");
  TStatToolkit::AddMetadata(tree,"TRDTOFmatchEffNegAll.Title","#epsilon_{TPC->TRD}");
  TStatToolkit::AddMetadata(tree,"TRDTOFmatchEffNegAll.AxisTitle","#epsilon_{TPC->TRD}");
  TStatToolkit::AddMetadata(tree,"TRDTOFmatchEffNegAll.Legend","(Q<0)");

}




void DrawSummaryTrending(TTree * tree){
  //
  // Make default trending plots
  // All branches and aliases from the main tree (qaConfig) and friend trees (defined in the qaTrendingTRD)
  // 
  //
  TMultiGraph * multiGraph=0;
  TMultiGraph * multiLine=0;
  SetDrawStyle();
  //
  // Plot 0. Draw, print basic time info
  //
  if (qaMap["Logbook.TRD"]){
    canvasQA->Clear(); 
    TLegend *legendTime=new TLegend(0.11,0.03,0.3,0.3,"T0 default QA - Run time information ");  legendTime->SetNColumns(2);   legendTime->SetBorderSize(0);
    multiGraph=TStatToolkit::MakeMultGraph(tree, "T0 <T>","runDuration;pauseDuration:tagID","",defMarkers,defColors,kTRUE,0.8,9, legendTime);
    multiGraph->Draw();
    legendTime->Draw();
    descriptionQA->DrawClone();
    TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
    TStatToolkit::DrawStatusGraphs(oaMultGr);
    canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin); canvasQA->Draw("ap");
    canvasQA->SaveAs("Logbook.TRD.Timing_Vs_TagID.png");
  }

  if (qaMap["Logbook.TRD"]){
    canvasQA->Clear(); 
    TLegend *legendStat=new TLegend(0.11,0.03,0.3,0.3,"T0 default QA - Statistic information "); legendStat->SetNColumns(2);   legendStat->SetBorderSize(0);
    multiGraph=TStatToolkit::MakeMultGraph(tree, "T0 <T>","totalEvents;totalEventsPhysics;Logbook.TRD.eventCountPhysics;Logbook.TPC.eventCountPhysics:tagID","",defMarkers,defColors,kTRUE,0.8,9, legendStat);
    multiGraph->Draw();
    descriptionQA->DrawClone();
    legendStat->Draw();    
    TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
    TStatToolkit::DrawStatusGraphs(oaMultGr);
    canvasQA->cd(1)->SetRightMargin(padRightMargin);     canvasQA->cd(2)->SetRightMargin(padRightMargin);     canvasQA->Draw("ap");
    canvasQA->SaveAs("Logbook.TRD.Statistic_vs_TagID.png");
  }
  //
  // Plot 1.)
  //  
  /****** TPC-TRD matching Efficiency  vs ID (run number or period, ...) ******/
  canvasQA->Clear(); 
  TLegend *legend1=new TLegend(0.11,0.03,0.3,0.3,"TPC->TRD matching efficiency");
  legend1->SetBorderSize(0);
  multiGraph=TStatToolkit::MakeMultGraph(tree, "TPC->TRD <#Epsilon>","TPCTRDmatchEffPosAll;TPCTRDmatchEffNegAll:tagID:TPCTRDmatchEffPosAllErr;TPCTRDmatchEffNegAllErr","",defMarkers,defColors,kTRUE,0.8,9, legend1);
  multiLine =TStatToolkit::MakeMultGraph(tree, "TPC->TRD <#Epsilon>","TPCTRDmatchEffPosAll_RobustMean;TPCTRDmatchEffPosAll_WarningMin;TPCTRDmatchEffPosAll_WarningMax;TPCTRDmatchEffPosAll_OutlierMin;TPCTRDmatchEffPosAll_OutlierMax:tagID", "","0;0;0;0;0;0;0;0",bandColors,kTRUE,0,0,0);
  multiGraph->SetMinimum(0.0);
  multiGraph->SetMaximum(1.5);
  multiGraph->Draw();
  legend1->Draw();
  multiLine->Draw();
  descriptionQA->DrawClone();
  TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin); canvasQA->Draw("ap");
  canvasQA->SaveAs("QA.TRD.MatchingEffTPCTRD_Vs_TagID.png");
  canvasQA->Clear();
  //
  //
  //
  /****** TRD->TOF matching Efficiency  vs ID (run number or period, ...) ******/
  canvasQA->Clear(); 
  TLegend *legend1=new TLegend(0.11,0.03,0.3,0.3,"TRD->TOF matching efficency");
  legend1->SetBorderSize(0);
  multiGraph=TStatToolkit::MakeMultGraph(tree, "TRD->TOF <#Epsilon>","TRDTOFmatchEffPosAll;TRDTOFmatchEffNegAll:tagID:TRDTOFmatchEffPosAllErr;TRDTOFmatchEffNegAllErr","",defMarkers,defColors,kTRUE,0.8,9, legend1);
  multiLine =TStatToolkit::MakeMultGraph(tree, "TRD->TOF <#Epsilon>","TRDTOFmatchEffPosAll_RobustMean;TRDTOFmatchEffPosAll_WarningMin;TRDTOFmatchEffPosAll_WarningMax;TRDTOFmatchEffPosAll_OutlierMin;TRDTOFmatchEffPosAll_OutlierMax:tagID", "","0;0;0;0;0;0;0;0",bandColors,kTRUE,0,0,0);
  multiGraph->SetMinimum(0.0);
  multiGraph->SetMaximum(1.5);
  multiGraph->Draw();
  legend1->Draw();
  multiLine->Draw();
  descriptionQA->DrawClone();
  TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin); canvasQA->Draw("ap");
  canvasQA->SaveAs("QA.TRD.MatchingEffTRDTOF_Vs_TagID.png");
  canvasQA->Clear();
  //
  //
   /****** TRD->TOF matching Efficiency  vs ID (run number or period, ...) ******/
  canvasQA->Clear(); 
  TLegend *legend1=new TLegend(0.11,0.03,0.3,0.3,"TRD->TOF matching efficency");
  legend1->SetBorderSize(0);
  multiGraph=TStatToolkit::MakeMultGraph(tree, "TRD->TOF <#Epsilon>","TPCTRDmatchEffPosAll-TPCTRDmatchEffNegAll;TRDTOFmatchEffPosAll-TRDTOFmatchEffNegAll:tagID:TRDTOFmatchEffPosAllErr;TRDTOFmatchEffNegAllErr","",defMarkers,defColors,kTRUE,0.8,9, legend1);
  multiGraph->SetMinimum(0.0);
  multiGraph->SetMaximum(1.5);
  multiGraph->Draw();
  legend1->Draw();
  descriptionQA->DrawClone();
  TStatToolkit::AddStatusPad(canvasQA, statPadHeight, statPadBotMar);
  TStatToolkit::DrawStatusGraphs(oaMultGr);
  canvasQA->cd(1)->SetRightMargin(padRightMargin); canvasQA->cd(2)->SetRightMargin(padRightMargin); canvasQA->Draw("ap");
  canvasQA->SaveAs("QA.TRD.MatchingEffPosMinusNeg_Vs_TagID.png");
  canvasQA->Clear();



}
