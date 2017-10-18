/*
  Append logbook metadata decribing tree structure,  and annotating varaibles
  
  Usage:
     1.) Metadata can be setup infoking macro:
         AliExternalInfo info;
         TTree * tree = info.GetTree("Logbook","LHC15o","cpass1_pass1","QA.TRD;QA.TOF;QA.TPC;Logbook.detector");     
         //
         metadataMacro="$ALICE_ROOT/STAT/Macros/logbookAddMetadata.C";
         gROOT->ProcessLine(TString::Format(".x %s((TTree*)%p,0);",metadataMacro,tree).Data());
     2.) Macro can be executed automatically if proper configuation file leaded AliExternalInfo.cfg
     

*/

void logbookAddMetadata(TTree*tree, Int_t verbose=0){
  //
  // Set metadata infomation 
  //
  if (tree==NULL) {
    ::Error("logbookAddMetadata","Start processing. Emtpy tree");
    return;
  }
  ::Info("logbookAddMetadata","Start processing Tree %s",tree->GetName());
  TObjArray * branches=tree->GetListOfBranches();
  for (Int_t ibr=0; ibr<branches->GetEntriesFast(); ibr++){
    TStatToolkit::AddMetadata(tree,TString::Format("%s.Description",branches->At(ibr)->GetName()).Data(),
			      TString::Format("Source alice logbook https.\n\t\t %s",tree->GetTitle()).Data());
    Bool_t isEventBranch= TString(branches->At(ibr)->GetName()).Contains("Events");
    Bool_t isDataBranch= TString(branches->At(ibr)->GetName()).Contains("Data");
    Bool_t isMagnet= TString(branches->At(ibr)->GetName()).Contains("magnetCurrent");
    TString brClass="Logbook";
    if (isEventBranch) {
      brClass+=" Stat";
      TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),
				TString::Format("counts").Data());      
    }

    if (isDataBranch)  {
      brClass+=" Data";
      TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),
				TString::Format("size (By)").Data());      
    }
    if (isMagnet) brClass+=" Setup";

    
    TStatToolkit::AddMetadata(tree,TString::Format("%s.class",branches->At(ibr)->GetName()).Data(),brClass.Data());    
  }
  // Fill Based and  custom metadata
  //
  // Index
  TStatToolkit::AddMetadata(tree,"run.class","Base Logbook Index");
  TStatToolkit::AddMetadata(tree,"LHCperiod.class","Base Logbook Index");
  TStatToolkit::AddMetadata(tree,"LHCFillNumber.class","Base Logbook Index");
  // Stat
  TStatToolkit::AddMetadata(tree,"totalEvents.class","Base Logbook Stat");
  TStatToolkit::AddMetadata(tree,"totalEventsPhysics.class","Base Logbook Stat");
  // Time
  TStatToolkit::AddMetadata(tree,"DAQ_time_start.class","Base Logbook Time");
  TStatToolkit::AddMetadata(tree,"DAQ_time_stop.class","Logbook Time");
  TStatToolkit::AddMetadata(tree,"runDuration.class","Base Logbook Time");
  TStatToolkit::AddMetadata(tree,"ctpDuration.class","Base Logbook Time");
  // Setup
  TStatToolkit::AddMetadata(tree,"detectorMask.class","Base Logbook Setup");
  TStatToolkit::AddMetadata(tree,"HLTmode.class","Base Logbook Setup");
  TStatToolkit::AddMetadata(tree,"numberOfDetectors.class","Base Logbook Setup");  
  //
   // Logbook html metadata
  TStatToolkit::AddMetadata(treeMC,"ctpDuration.thead","&Delta;T(s)");
  TStatToolkit::AddMetadata(treeMC,"ctpDuration.tooltip","Logbook CTP duration ");
  TStatToolkit::AddMetadata(treeMC,"ctpDuration.html","<a href=https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=rund&p_run=%d{run}&p_tab=gi&p_subtab=rs>%d{ctpDuration}</a>");
  TStatToolkit::AddMetadata(treeMC,"detectorMask.thead","mask");
  TStatToolkit::AddMetadata(treeMC,"detectorMask.tooltip","Detector mask");
  TStatToolkit::AddMetadata(treeMC,"DAQ_time_start.thead","T<sub>start</sub>");
  TStatToolkit::AddMetadata(treeMC,"LHCperiod.thead","period");
  TStatToolkit::AddMetadata(treeMC,"LHCFillNumber.thead","fill");
  TStatToolkit::AddMetadata(treeMC,"HLTmode.thead","H L T");
  TStatToolkit::AddMetadata(treeMC,"numberOfDetectors.thead","N<sub>det</sub>");
  TStatToolkit::AddMetadata(treeMC,"numberOfDetectors.tooltip","Number of detectors");
  TStatToolkit::AddMetadata(treeMC,"totalEventsPhysics.thead","N<sub>evphys</sub>");
  TStatToolkit::AddMetadata(treeMC,"interactionRate.thead","rate(Hz)");

  TStatToolkit::AddMetadata(treeMC,"DAQ_time_start.html","<a href=https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=rund&p_run=%d{run}&p_tab=gi&p_subtab=rs>%t{DAQ_time_start}</a>");
  TStatToolkit::AddMetadata(treeMC,"run.html","<a href=\"000%d{run}/index.html\">%d{run}</a>");
  TStatToolkit::AddMetadata(treeMC,"detectorMask.html","%x{detectorMask&0xFFFF}");
  TStatToolkit::AddMetadata(treeMC,"interactionRate.html","%2.2f{interactionRate}");

  //
  TList * mlist = (TList*)(tree->GetUserInfo()->FindObject("metaTable"));
  mlist->Sort();
  if (verbose==1){
    mlist->Print();
  }
  if (verbose==2){
    // AliTreeToolkit::printMetadata(treeLogbook, "[class]",0);  // not yet comitted
  }
  ::Info("logbookAddMetadata","End");

}
