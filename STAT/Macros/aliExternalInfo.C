/// \ingroup Macros
/// \file    aliExternalInfo.C
/// \brief   Demo usage of the information from the AliExternalInfo and visualization using the TStatToolkit, AliTreePlayer and AliDrawStyle
///
/// \author  Marian Ivanov
///
/// Demo usage of the information from the AliExternalInfo and visualization using the TStatToolkit, AliTreePlayer and AliDrawStyle
/// See documentation in following functions:
/// * logbook/QA/RCT  multi figure export multipad
///   * ::drawLogbook
/// * logbook/QA/RCT  multi figure export multigraph
/// * TODO:  logbook/QA/RCT  multi figure export multiselection
/// * TODO:  makeHtml
/// * TODO:  add period pass to the legend in standard plots

/// The code can be used as it is.
/// It will be running as part of UnitTest
/// ## Example usage
/// ###   Load library and define some example variables
/*!
    .L $AliRoot_SRC/STAT/Macros/aliExternalInfo.C
    // define some example selection and copy past the code below
    TString  period="LHC15o", pass="pass1", runSelection="QA.TPC.meanTPCncl>0", varSelection="run:runDuration:totalEventsPhysics:totalNumberOfFilesMigrated:QA.TPC.meanMIP";
    TString source="Logbook;QA.TPC;QA.TRD;QA.ITS;MonALISA.RCT;QA.EVS";
*/

AliExternalInfo info;
TLatex latex;

/// AliExternalInfo demo - Generic draw of the logbook information - per period
/// external source .e.g "QA.TPC;QA.TRD;QA.ITS;MonALISA.RCT" could be used in selection.
/// \param period
/// \param pass
/// \param runSelection
/// \param varSelection
/*!
 * ### Example usage:
\code
    .L $AliRoot_SRC/STAT/Macros/aliExternalInfo.C
    drawLogbook("LHC15o","pass1","QA.TPC.meanTPCncl>0","runDuration:totalEventsPhysics:totalNumberOfFilesMigrated:MonALISA.RCT.tpc_value");
    drawLogbook("LHC10h","pass2","totalEventsPhysics>1000&&totalNumberOfFilesMigrated>20","runDuration:totalEventsPhysics:totalNumberOfFilesMigrated:MonALISA.RCT.tpc_value");
    drawLogbook("LHC11h","pass2","totalEventsPhysics>1000&&totalNumberOfFilesMigrated>20","runDuration:totalEventsPhysics:totalNumberOfFilesMigrated:MonALISA.RCT.tpc_value");
\endcode
*/
void drawLogbook(TString period, TString pass,TString runSelection="", TString varSelection="runDuration:totalEventsPhysics:totalNumberOfFilesMigrated"){
  TTree * treeLogbook = info.GetTree("Logbook",period,pass,"QA.TPC;QA.TRD;QA.ITS;MonALISA.RCT");
  TObjArray *varList=0;
  if (varSelection[0]=="["){ //variable list using class selection
    // Use class selection to select variables
    varList=AliTreePlayer::selectMetadata(treeLogbook,varSelection,0);
    Int_t nvars=varList->GetEntries();
    for (Int_t i=0; i<nvars;i++){
      TString vname=varList->At(i)->GetName();
      vname.ReplaceAll(".class","");
    }
    // varList=AliTreePlayer::selectMetadata(treeLogbook, "[class==\"Logbook&&Time\"]",0);
  }else{
    varList=varSelection.Tokenize(":");
  }
  Int_t nvars=varList->GetEntries();
  TCanvas * canvas = new TCanvas("drawLogbook","drawLogbook",1200,1000);
  canvas->Divide(1,nvars);
  for (Int_t i=0; i<nvars; i++){
    canvas->cd(i+1);
    TString var=TString::Format("%s:run",varList->At(i)->GetName());
    TStatToolkit::MakeGraphSparse(treeLogbook,var,runSelection,25,1,1)->Draw("ap");
    Double_t mean = TMath::Mean(treeLogbook->GetSelectedRows(),treeLogbook->GetV1());
    Double_t sum  = treeLogbook->GetSelectedRows()*mean;
    latex.DrawLatexNDC(0.15,0.8,TString::Format("Mean=%0.0f",mean));
    latex.DrawLatexNDC(0.15,0.7,TString::Format("Sum=%0.0f",sum));
  }
  canvas->SaveAs(TString::Format("%s_%s.png",pass.Data(),period.Data()).Data());
}

/// AliExternalInfo - MultiGraph example
/*!
 * @param period        -  period - e.g LHC15o
 * @param pass          -  reconstruction pass  - e.g pass2 (needed in case QA selection)
 * @param runSelection  -  runSelection (standard tree cut )
 * @param varSelection  -  array of variables to show
*/
/*!
 * ## Example usage:
 \code
   .L $AliRoot_SRC/STAT/Macros/aliExternalInfo.C
  drawLogbookMultiExpr("LHC15o","pass1", "totalNumberOfFilesMigrated>2000", "runDuration;totalEventsPhysics/1000;totalNumberOfFilesMigrated:run");
  drawLogbookMultiExpr("LHC15o","pass1", "totalNumberOfFilesMigrated>2000&&QA.TPC.meanMIP>40", "runDuration;totalEventsPhysics/1000;totalNumberOfFilesMigrated:run");
\endcode
*/
void drawLogbookMultiExpr(TString period, TString pass,TString runSelection="", TString varSelection="runDuration;totalEventsPhysics/1000;totalNumberOfFilesMigrated:run") {
  TTree *treeLogbook = info.GetTree("Logbook", period, pass, "QA.TPC;QA.TRD;QA.ITS;MonALISA.RCT");
  TLegend * legend = new TLegend(0.15,0.7,0.4,0.8, "Trending");
  TMultiGraph *m0 = TStatToolkit::MakeMultGraph(treeLogbook, "", varSelection.Data(), runSelection.Data(), "25;21;22;24", "1;2;4", kTRUE, 1, 6, legend);
  m0->Draw("ap");
  legend->Draw();
  delete tree;
}

/// drawLogbookMultiCut - example of usage of multigraph
/// \param period
/// \param pass
/// \param runSelection
/// \param varSelection
/// TODO - implemnetation missing. Only hardwired selection
/*!
 */
void drawLogbookMultiCut(TString period, TString pass,TString runSelection="", TString varSelection="LHCperiod:runDuration:totalEventsPhysics:totalNumberOfFilesMigrated") {
  TTree *treeLogbook = info.GetTree("Logbook", period, pass, "QA.TPC;QA.TRD;QA.ITS;MonALISA.RCT");
  TMultiGraph *m1 = TStatToolkit::MakeMultGraph(treeLogbook,"","runDuration:run","totalNumberOfFilesMigrated>2000;(run<245000);abs(run-245500)<500;(run>246000)","figTemplateTRD","figTemplateTRD", kFALSE, 1,6,0);
  m1->Draw("ap");
}

/// makeHTMLPage - make a standard html page
/// \param period        - period               e.g.   : LHC15o
/// \param pass          - reconstruction pass  e.g.   : pass1
/// \param runSelection  - run selection cut    e.g.   : MonALISA.RCT.tpc_value>0&&QA.TPC.meanTPCncl>0
/// \param varSelection  - list of variables to export : run:runDuration:totalEventsPhysics:totalNumberOfFilesMigrated:QA.TPC.meanMIP:QA.TPC.meanMIP
/// \param source"       - source list          e.g:   : Logbook;QA.TPC;QA.TRD;QA.ITS;MonALISA.RCT
/// ##Example usage:
/*!
\code
 .L $AliRoot_SRC/STAT/Macros/aliExternalInfo.C
 makeHTMLPage("LHC15o","pass1", "QA.TPC.meanTPCncl>0", "run:runDuration:totalEventsPhysics:totalNumberOfFilesMigrated:QA.TPC.meanMIP:QA.TPC.meanMIP", "Logbook;QA.TPC;QA.TRD;QA.ITS;MonALISA.RCT")
 makeHTMLPage("LHC15o","pass1", "QA.TPC.meanTPCncl>0", "run:runDuration:totalEventsPhysics:ocdbStatusCounter:ocdbHVStatusCounter:TPC_Status:meanTPCncl_Status:PID_Status:DCAz_Status:DCAr_Status:tpcItsMatch_Status", "Logbook;QA.TPC;QA.TRD;QA.ITS;MonALISA.RCT")
 makeHTMLPage("LHC15o","pass1", "QA.TPC.meanTPCncl>0", "run:runDuration:totalEventsPhysics:ocdbStatusCounter:ocdbHVStatusCounter:MonALISA.RCT.tpc_value:rctMismatch:TPC_Status:meanTPCncl_Status:PID_Status:DCAz_Status:DCAr_Status:tpcItsMatch_Status", "Logbook;QA.TPC;QA.TRD;QA.ITS;MonALISA.RCT")
//
 makeHTMLPage("LHC15o","pass1", "QA.TPC.meanTPCncl>0", "run:runDuration:#%d{year}/%d{period.GetName()}/%d{pass.GetName()}/%d{run}:totalEventsPhysics:totalNumberOfFilesMigrated:QA.TPC.meanMIP", "Logbook;QA.TPC;")
 \endcode
 */
void makeHTMLPage(TString  period,TString  pass, TString runSelection, TString varSelection, TString source){
  TTree * tree = info.GetTree("Logbook",period,pass,source);
  AliTreePlayer::selectWhatWhereOrderBy(tree,varSelection.Data(), runSelection.Data(),"",0,100000,"html","table.html");
  delete tree;
  gSystem->GetFromPipe("$AliPhysics_SRC/PWGPP/scripts/makeHtmlv1.sh index.html table.html 0");

}
