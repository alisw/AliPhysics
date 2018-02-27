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
#include "TString.h"
#include "TObjArray.h"
#include "TTree.h"
#include "AliExternalInfo.h"
#include "TLatex.h"
#include "TLeaf.h"
#include "TSystem.h"
#include "TStatToolkit.h"
#include "AliTreePlayer.h"
#include "TLegend.h"

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
  if (varSelection[0]=='['){ //variable list using class selection
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
  delete treeLogbook;
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

/// \brief Cache  MC production trees, store summary information in formated text files -> root trees
/// \param dataType  -
/// \param fileList
void CacheTestMCProductions(TString dataType, const char *fileList=NULL){
  AliExternalInfo info;
  info.fLoadMetadata=kFALSE;
  TObjArray* periodList = NULL;
  TArrayI nRuns;
  if (fileList!=NULL) {
    periodList=(gSystem->GetFromPipe(TString::Format("cat %s", fileList).Data())).Tokenize("\n");
    nRuns.Set(periodList->GetEntries());

  }else{
    TTree * tree = info.GetTree("MonALISA.ProductionMC","","");
    Int_t nProd=tree->GetEntries();
    periodList = new TObjArray(nProd);
    nRuns.Set(nProd);
    TLeaf *leaf = tree->GetLeaf("Tag");
    TLeaf *leafRuns = tree->GetLeaf("Number_of_runs");
    for (Int_t iProd=0; iProd<nProd; iProd++){
      tree->GetEntry(iProd);
      TString prodName=((char*)leaf->GetValuePointer());
      if (prodName.Contains("LHC")==0) continue;
      periodList->AddAt(new TObjString(((char*)leaf->GetValuePointer())),iProd);
      nRuns[iProd]=leafRuns->GetValue();
    }
    delete tree;
  }
  for (Int_t iPeriod=0; iPeriod<periodList->GetEntriesFast(); iPeriod++){
    TTree* tree = info.GetTree(dataType.Data(),periodList->At(iPeriod)->GetName(),"passMC");
    if (tree){
      Int_t entries=tree->Draw("run","1","goff");
      TString sInfo=periodList->At(iPeriod)->GetName();
      sInfo+="\t";
      sInfo+=dataType;
      sInfo+="\t";
      sInfo+=TString::Format("%d\t",entries);
      sInfo+=TString::Format("%d\t",nRuns[iPeriod]);
      for (Int_t j=0; j<entries; j++) {
        sInfo+=TString::Format("%2.0f,",tree->GetV1()[j]);
        ::Info("CacheTestMCProductionsRun:","%s\t%s\t%d\t%d\t%d\t%2.0f",periodList->At(iPeriod)->GetName(),dataType.Data(),entries,nRuns[iPeriod],j, tree->GetV1()[j]);
      }
      sInfo+="0";
      ::Info("CacheTestMCProductionsPeriod:","%s\n",sInfo.Data());
      delete tree;
    }else{
      ::Error("CacheTestMCProductionsPeriod:","%s\t%s\t-1\t%d\t0",periodList->At(iPeriod)->GetName(), dataType.Data(),nRuns[iPeriod]);
    }
  }
}


/// Cache MC production information
void CacheTrendingProductions(TString dataType){
  AliExternalInfo info;
  info.fLoadMetadata=kFALSE;
  TObjArray* periodList = NULL, *idList=NULL;
  //
  TTree * tree = info.GetTree("MonALISA.ProductionCycle","","");
  Int_t nProd=tree->GetEntries();
  periodList = new TObjArray(nProd);
  idList= new TObjArray(nProd);
  TLeaf *leafTag = tree->GetLeaf("Tag");
  TLeaf *leafID   =  tree->GetLeaf("ID");
  for (Int_t iProd=0; iProd<nProd; iProd++){
    tree->GetEntry(iProd);
    TString prodName=((char*)leafTag->GetValuePointer());
    TString  idName =TString::Format("%d",TMath::Nint(leafID->GetValue()));
    if (prodName.Contains("LHC")==0) continue;
    periodList->AddAt(new TObjString(prodName),iProd);
    idList->AddAt(new TObjString(idName),iProd);
  }
  delete tree;
  //
  for (Int_t iPeriod=0; iPeriod<periodList->GetEntriesFast(); iPeriod++) {
    TTree* treeP = info.GetTreeProdCycleByID(idList->At(iPeriod)->GetName());
    if (treeP==NULL) continue;
    TLeaf *leafOutput = treeP->GetLeaf("outputdir");
    Int_t nRuns= treeP->GetEntries();
    treeP->GetEntry(0);
    TString path=((char*)leafOutput->GetValuePointer());
    TObjArray *pArray = path.Tokenize("/");
    if (pArray==NULL) continue;
    Int_t nElems=pArray->GetEntries();
    if (nElems<4) continue;
    TString aperiod=pArray->At(3)->GetName();
    TString apass =pArray->At(nElems-1)->GetName();
    delete pArray;
    ::Info("CacheTrendingProductions","%s\t%s\t%s\t%s\t%d",idList->At(iPeriod)->GetName(),path.Data(), aperiod.Data(),apass.Data(),nRuns);
    delete treeP;
    TTree* treeQA = info.GetTree(dataType.Data(),aperiod.Data(),apass.Data());
    if (treeQA){
      Int_t entries=treeQA->Draw("run","1","goff");
      TString sInfo=aperiod;
      sInfo+="\t";
      sInfo+=apass;
      sInfo+="\t";
      sInfo+=dataType;
      sInfo+="\t";
      sInfo+=TString::Format("%d\t",entries);
      sInfo+=TString::Format("%d\t",nRuns);
      for (Int_t j=0; j<entries; j++) {
        sInfo+=TString::Format("%2.0f,",treeQA->GetV1()[j]);
        ::Info("CacheTrendingProductionsRun:","%s\t%s\t%s\t%d\t%d\t%2.0f",aperiod.Data(),apass.Data(),dataType.Data(),entries,nRuns,treeQA->GetV1()[j]);
      }
      sInfo+="0";
      ::Info("CacheTrendingProductionsPeriod:","%s\n",sInfo.Data());
      delete treeQA;
    }else{
      ::Error("CacheTrendingProductionsPeriod:","%s\t%s\t%s\t-1\t%d\t0",aperiod.Data(),apass.Data(), dataType.Data(),nRuns);
    }
  }
}

void CheckProductions(){
  AliExternalInfo info;
  //to add there production yer
  TTree * treeRaw= info.GetTree("QA.Period","data","");
  treeRaw->SetAlias("isTPC","type==\"QA.TPC\"");
  treeRaw->SetAlias("isITS","type==\"QA.ITS\"");
  treeRaw->SetAlias("isTRD","type==\"QA.TRD\"");
  // black list for production
  treeRaw->SetAlias("isBlack","strstr(pass,\"clean\")!=0||strstr(pass,\"rec\")!=0||strstr(pass,\"its\")!=0||strstr(pass,\"cpass\")!=0||strstr(pass,\"vpass\")!=0||strstr(pass,\"muon\")!=0||strstr(pass,\"cosmic\")!=0||strstr(pass,\"align\")!=0||strstr(pass,\"FAST\")!=0||strstr(pass,\"scan\")!=0||strstr(pass,\"test\")!=0");

}
