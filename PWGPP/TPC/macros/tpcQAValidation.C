/// \ingroup PWGPP/TPC/macros/
/// \brief  RAW QA using cimpbined QA information
/// JIRA
///      - $NOTESData/JIRA/ATO-360 and ATO-361
/// \author marian  Ivanov marian.ivanov@cern.ch
///
/*!
\code
  gSystem->AddIncludePath("-I$ALICE_ROOT/include/"); //couldn't add include path in .rootr
  .L $AliPhysics_SRC/PWGPP/TPC/macros/tpcQAValidation.C+
  TString period="LHC15o";
  TString pass="pass3_lowIR_pidfix";
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  //
  InitTPCValidation("LHC15o", "pass3_lowIR_pidfix",0,0); //short period
  //InitTPCValidation("LHC17p","cpass1_pass1",0,0);
  //
  MakeReport();
 
  MakeStatusPlots();
  trendingDraw->fReport->Close();
  //
  makeHtml();

\endcode
*/

// TODO - show error message in case f missing information - instaead of failing
//      - probelem can be emulated e.g excluding QA.EVS input
//

#include <TError.h>
#include <TMultiGraph.h>
#include "TCanvas.h"
#include "TLatex.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TTreeFormula.h"
#include "AliExternalInfo.h"
#include "AliTreeTrending.h"
#include "TStatToolkit.h"
#include "TPRegexp.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "AliDrawStyle.h"
#include "AliTreePlayer.h"

AliExternalInfo *externalInfo = 0;
AliTreeTrending *trendingDraw = 0;
TTree *treeQA;
std::vector <Double_t> cRange;
std::vector <Double_t> cRange2;
std::vector <Double_t> cRange5;
TString sCriteria("(present):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)");
// current variables - used instead of local variables to enable copy paste of part of code
TString queryString, queryTitle; ///
TMultiGraph *graph = 0, *lines = 0;
const char *outputDir="./";

/// \param treeQA
/// \param doCheck
/// \param verbose
void makeTPCAlarms(TTree *treeQA, Bool_t doCheck, Int_t verbose);
Bool_t InitTPCValidation(TString period, TString pass, Int_t verbose,  Int_t doCheck);
void MakeReport();
void MakeStatusPlots();
//void tpcQAValidation(const char *period = "LHC15o", const char *pass= "pass3", const char *soutputDir = "./");
void makeHtml();
// specific functions
void makeHtmlDCA();
void MakeJSROOTHTML(TString prefix, TString outputName);

void LoadStyles(){
  AliDrawStyle::RegisterCssStyle("AliTreeTrending", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/AliTreeTrending.css", 0));
  AliDrawStyle::RegisterCssStyle("AliTreeTrending_QATPC", AliDrawStyle::ReadCSSFile("$AliPhysics_SRC/PWGPP/TPC/macros/TPCQAWebpage/QAtabs/AliTreeTrending_TPCQA.css", 0));

}


/// function to create a set of the comparison plots MC/AnchorRaw data
/// \param mcPeriod    -  MC production name
/// \param outputDir   -  directory where png and html files are stored
void tpcQAValidation(const char *period, const char *pass, const char *sOutputDir) {
  if (InitTPCValidation(period, pass, 0, 0)) {
    MakeReport();
  } else ::Error("tpcQAValidation", "InitTPCValidation returned with error -> skip plotting!");
}

/// makeTPCAlarms - Set of expression aliases/formulas  (Warning,Outlier,PhysAcc)
/// define variables for alarms for MC/raw mismatch
/// * 1.) Absolute aliases
/// * 2.) Relative aliases
/// * 3.) Combined status aliases
/// \param treeQA  - input tree
/// \param doCheck - force check of the variables
/// \param verbose - set verbosity for make alarms
void makeTPCAlarms(TTree * treeQA, Bool_t doCheck,Int_t verbose){
  ::Info("makeTPCalarms","Done with aliases");
  // Alarms for the MIP QA per sector OFF
  treeQA->SetAlias("gainCalib_Warning","Max$(abs(meanMIPvsSector.fElements/meanMIP-1)*(meanMIPvsSector.fElements/meanMIP>0.1))>0.015"); // Warning    1.5 % OFF
  treeQA->SetAlias("gainCalib_Outlier","Max$(abs(meanMIPvsSector.fElements/meanMIP-1)*(meanMIPvsSector.fElements/meanMIP>0.1))>0.03");  // Outlier    3.% % OFF
  treeQA->SetAlias("gainCalib_PhysAcc","Max$(abs(meanMIPvsSector.fElements/meanMIP-1)*(meanMIPvsSector.fElements/meanMIP>0.1))<0.02");  // Acceptable 2.% % ONN
  treeQA->SetAlias("MIPquality_Warning",TString::Format("%s||gainCalib_Warning", treeQA->GetAlias("MIPquality_Warning")));
  treeQA->SetAlias("MIPquality_Outlier",TString::Format("%s||gainCalib_Outlier", treeQA->GetAlias("MIPquality_Outlier")));
  treeQA->SetAlias("MIPquality_PhysAcc",TString::Format("%s||gainCalib_PhysAcc", treeQA->GetAlias("MIPquality_PhysAcc")));




}
///
/// \param period
/// \param pass
/// \param verbose
/// \param doCheck
/// \return
Bool_t InitTPCValidation(TString period, TString pass, Int_t verbose, Int_t doCheck) {

  cRange = std::vector < Double_t > {0.13, 0.01, 0.5, 0.35};
  cRange2 = std::vector < Double_t > {0.13, 0.01, 0.5, 0.3};
  cRange5 = std::vector < Double_t > {0.13, 0.01, 0.8, 0.3};
  externalInfo = new AliExternalInfo(".", "", verbose);
  trendingDraw = new AliTreeTrending("QA validation","QA validation");
  trendingDraw->SetDefaultStyle();
  gStyle->SetOptTitle(0);
  treeQA = externalInfo->GetTree("QA.TPC", period, pass, "QA.TPC;QA.TRD;QA.TOF;QA.ITS;Logbook;QA.EVS;QA.rawTPC;Logbook.detector:TPC:detector==\"TPC\"");
  //
  makeTPCAlarms(treeQA, doCheck, verbose);
  TString sStatusBarVars("MIPquality;dcaz;dcar;tpcItsMatch;meanTPCncl;itsTpcPulls");
  TString sStatusBarNames("MIPquality;dcaz;dcar;tpcItsMatch;meanTPCncl;itsTpcPulls");
  TString sCriteria("(present):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)"); // status bar markers
  TString statusString[3];
  statusString[0] = sStatusBarVars;
  statusString[1] = sStatusBarNames;
  statusString[2] = sCriteria;

  trendingDraw->SetTree(treeQA);
  treeQA->SetAlias("tagID", "run");
  treeQA->SetAlias("present","run>0");
  treeQA->SetAlias("defaultCut", "run==QA.TPC.run");
  Bool_t initStatus=trendingDraw->InitSummaryTrending(statusString, 0.015, "defaultCut");
  makeTPCAlarms(treeQA,doCheck,verbose);

  // treeQA->SetAlias("defaultCut", "run>0");   // TODO enable all runs - need to make reorting of the multigraphs
  if (initStatus) {
    trendingDraw->fWorkingCanvas->Print(TString(outputDir) + "/report.pdf[", "pdf");
    return kTRUE;
  }else {

    return kFALSE;
  }
}

/// ## MakeReport
/// \param outputDir
void MakeReport() {

  trendingDraw->fWorkingCanvas->Clear();
  trendingDraw->fWorkingCanvas->Print(TString(outputDir) + "/report.pdf]", "pdf");
  //
    // 1.) Event properties  ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabEvent.html)
  trendingDraw->MakePlot(outputDir, "interactionRate.png", "Interaction rate", cRange, "",
                         "Logbook.averageEventsPerSecond;QA.EVS.interactionRate:run", "defaultCut", "figTemplateTRD",
                         "figTemplateTRD", 1, 1, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "runDuration.png", "Run duration ", cRange, "", "Logbook.runDuration:run",
                         "defaultCut", "figTemplateTRD", "figTemplateTRD", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "bField.png", "Magnet current", cRange, "", "Logbook.L3_magnetCurrent:run",
                         "defaultCut", "figTemplateTRD", "figTemplateTRD", 1, 1.0, 3, kTRUE);
  trendingDraw->MakePlot(outputDir, "eventCounters.png", "Event counters ", cRange, "",
                         "Logbook.totalEvents;totalEventsPhysics;totalEventsCalibration;Logbook.detector_TPC.eventCountPhysics;TPC.Anchor.nEvents;QA.TPC.nEvents:run",
                         "defaultCut", "figTemplateTRD", "figTemplateTRD", 1, 1.0, 4, kFALSE);




}


/// TODO - alias for the html table names, hints
/// TODO - fix formatting in AliTreeFormulaF- in case of missing entry write undefined
/// TODO - tooltip - explaining status bits (detector mask, resp. status mask)
void makeHtml() {
  //
  //
  // EVS
  TStatToolkit::AddMetadata(treeQA, "interactionRate.thead", "rate(Hz)");
  TStatToolkit::AddMetadata(treeQA, "interactionRate.tooltip", "Interaction rate (source EVS)");
  TStatToolkit::AddMetadata(treeQA, "bz.thead", "B<sub>z</sub> (T)");
  TStatToolkit::AddMetadata(treeQA, "bz.tooltip", "Barrel field");
  TStatToolkit::AddMetadata(treeQA, "detectorMask.html","%b{detectorMask&0x7F}");
  TStatToolkit::AddMetadata(treeQA, "detectorMask.headerTooltip","Detector bit mask:\n . SPD\n . SDD\n . SSS\n . TPC\n . TRD");
  //MakeStatusBitMasks();
  makeHtmlDCA();
}

void makeHtmlDCA(){

  //TString metaDCA="<a title
  TString pathMC="<a href=\"http://aliqatpc.web.cern.ch/aliqatpc/sim/%d{year}/%s{period.GetName()}/passMC/000%d{run}/";
  TString pathData="<a href=\"http://aliqatpc.web.cern.ch/aliqatpc/data/%d{TPC.Anchor.year}/%s{TPC.Anchor.period.GetName()}/%{TPC.Anchor.pass.GetName()}/000%d{run}/";
  TStatToolkit::AddMetadata(treeQA, "dcarAP0.html", (pathMC+"dca_and_phi.png\">%2.4f{dcarAP0}</a>").Data());
  TStatToolkit::AddMetadata(treeQA, "TPC.Anchor.dcarAP0.html", (pathData+"dca_and_phi.png\">%2.4f{TPC.Anchor.dcarAP0}</a>").Data());
  TStatToolkit::AddMetadata(treeQA, "dcarCP0.html", (pathMC+"dca_and_phi.png\">%2.4f{dcarCP0}</a>").Data());
  TStatToolkit::AddMetadata(treeQA, "TPC.Anchor.dcarCP0.html", (pathData+"dca_and_phi.png\">%2.4f{TPC.Anchor.dcarCP0}</a>").Data());

  TStatToolkit::AddMetadata(treeQA, "dcarAP0.thead", "&sigma;<sup>A</sup><sub>MC DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeQA, "dcarAP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeQA, "dcarCP0.thead", "&sigma;<sup>C</sup><sub>MC DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeQA, "dcarCP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeQA, "TPC.Anchor.dcarAP0.thead", "&sigma;<sup>A</sup><sub>Anchor DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeQA, "TPC.Anchor.dcarAP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeQA, "TPC.Anchor.dcarCP0.thead", "&sigma;<sup>C</sup><sub>Anchor DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeQA, "TPC.Anchor.dcarCP0.tooltip", "DCA (r&phi;) at infinite pt");
   //
  TStatToolkit::AddMetadata(treeQA, "dcar_MaskWarning.thead", "DCA<sub>MCWarning</sub></sub>");
  TStatToolkit::AddMetadata(treeQA, "dcar_MaskWarning.tooltip", "DCAr status");
  TStatToolkit::AddMetadata(treeQA, "dcar_MaskWarning.html","<a href=\"dcarStatusMC.png\">%x{dcar_MaskWarning}</a>");
  //
  TString runSelection = "defaultCut";
  TString varSelection = "";
  TString logbookBase = "";
  logbookBase = "run:LHCFillNumber:LHCperiod:detectorMask:HLTmode:DAQ_time_start:ctpDuration:totalEventsPhysics:interactionRate:bz;2.2:";
  // tabDCA
  TString tpcDCA = "dcarAP0:TPC.Anchor.dcarAP0:dcarCP0:TPC.Anchor.dcarCP0:dcar_MaskWarning";//:dcarStatusString:dcarRawStatusString";
  // TString tpcDCA="dcarAP1;2.2f:dcarStatusString:dcarStatusOutlier:dcarStatusWarning";
  // TObjArray *tpcDCAArray = AliTreePlayer::selectMetadata(treeQA->GetFriend("QA.TPC"), "[class==\"TPC&&DCAr&&!Err&&!Chi2\"]",0);
  AliTreePlayer::selectWhatWhereOrderBy(treeQA, (logbookBase + tpcDCA).Data(), runSelection.Data(), "", 0, 100000,
                                        "html", "stableDCA.html");
  gSystem->GetFromPipe("$AliPhysics_SRC/PWGPP/scripts/makeHtmlv1.sh index.html stableDCA.html 600");
}



/// Test html link printing
/// \param testHtml
/*!
 Used to test correctness of the format string e.g.:
 \code
 PrintTestHtmlLink("<a href=\"http://aliqatpc.web.cern.ch/aliqatpc/data/%d{year}/%s{TPC.Anchor.period.GetName()}/{TPC.Anchor.pass.GetName()}/000%d{run}/dca_and_phi.png\">%2.2f{dcarCP0}</a>")
\endcode
 */
void PrintTestHtmlLink(TString testHtml){
  AliTreeFormulaF testLink("testLink",testHtml,treeQA,7);
  //testLink.fDebug=1;
  printf("\n %s \n",testLink.PrintValue());
}

///
/// \param path
/// \param figName
/// \param statusString
/// \param selection
/// \param friendName
void MakeStatusPlot(const char *path, TString figName, TString statusString, TString selection, TString friendName=""){
  TPRegexp suffix("_Warning$");
  TString currentString="",statusVar="",statusTitle="", statusMask="";
  Int_t counter=0;
  AliTreeTrending::DecomposeStatusAlias(treeQA, statusString,statusVar,statusTitle,suffix,counter,statusMask);
  statusVar+="run";
  trendingDraw->MakeStatusPlot(path, figName.Data(), statusVar, statusTitle, selection.Data(), sCriteria,friendName);
}


///
/// \param statusString
/// \param selection
/// \param friendName
void MakeStatusBitMask(TString statusString){
  TPRegexp suffix("_Warning$");
  TString currentString="",statusVar="",statusTitle="",statusMask="";
  Int_t counter=0;
  AliTreeTrending::DecomposeStatusAlias(treeQA, statusString,statusVar,statusTitle,suffix,counter,statusMask);
  statusString.ReplaceAll("_Warning","");
  TString statusAliasName, statusAliasMask;
  //
  statusAliasName=statusString+"_MaskWarning";
  statusAliasMask=statusMask;
  statusAliasMask.ReplaceAll("_#","_Warning");
  treeQA->SetAlias(statusAliasName.Data(),statusAliasMask.Data());

}


void MakeStatusPlots(){
  try {
    MakeStatusPlot("./", "dcarStatus.png", "dcar_Warning", "1");
    MakeStatusPlot("./", "dcazStatus.png", "dcaz_Warning", "1");
    MakeStatusPlot("./", "dEdxStatus.png", "MIPquality_Warning", "1");


//    MakeStatusPlot("./", "itsEffStatusMCToAnchor.png", "mcAnchor.itsEffStatus_Warning", "1");
////    MakeStatusPlot("./", "itsEffStatusMC.png", "itsEffStatus_Warning", "1");                  //missing corresponding alias in MC tree
////    MakeStatusPlot("./", "itsEffStatusAnchor.png", "itsEffStatus_Warning", "1","TRD.Anchor"); //missing corresponding alias in Anchor tree
//    //
//
//    MakeStatusPlot("./", "dEdxStatusAnchor.png", "MIPquality_Warning", "1","TPC.Anchor");
//    MakeStatusPlot("./", "dEdxStatusMCToAnchor.png", "mcAnchor.dEdx_Warning", "1");
//
//    MakeStatusPlot("./", "nclStatusMCToAnchor.png", "mcAnchor.ncl_Warning", "1");
////    MakeStatusPlot("./", "nclStatusAnchor.png", "ncl_Warning", "1","TPC.Anchor");    //missing corresponding alias
////    MakeStatusPlot("./", "nclStatusMC.png", "ncl_Warning", "1");                     //missing corresponding alias
  }
  catch (const std::invalid_argument& ia) {
	  std::cerr << "Invalid argument: " << ia.what() << '\n';
  }
}

void MakeStatusBitMasks(){
  MakeStatusBitMask("dcar_Warning");
  MakeStatusBitMask("mcAnchor.dcarResol_Warning");
}


/*!
 * MakeJSROOTHTML make hsroot html page - see alose READMEjsrootQA.md
 * * plots from the report.root file divided into categories using regular exression
 * @param prefix      - html prefix for jsroot see example testing jsroot on local node localhost
 * @param outputName  - name of the output html file
\code
    TString outputName="jsrootMCAnchor.html";
    TString prefix="http://localhost:90/data/jsroot/jsroot/index.htm?file=http://localhost:90/data/alice-tpc-notes/JIRA/ATO-83/test/report.root"
    MakeJSROOTHTML(prefix,outputName);
\endcode
*/
void MakeJSROOTHTML(TString prefix, TString outputName){
  TString description="";
  TString figList="";
  TString items;
  TFile *fReport= TFile::Open("report.root");
  FILE * pFile;
  pFile = fopen (outputName.Data(),"w");
  //
  items=AliTreeTrending::ArrayNameToString(fReport->GetListOfKeys(),".*tatus.*",",");
  AliTreeTrending::AddJSROOTHtmlLink(pFile,"Status", prefix, Form("items=[%s]",items.Data()));
  items=AliTreeTrending::ArrayNameToString(fReport->GetListOfKeys(),"(.*vent.*|.*Mult.*|.*Vert.*)",",");
  AliTreeTrending::AddJSROOTHtmlLink(pFile,"General", prefix, Form("items=[%s]",items.Data()));
  items=AliTreeTrending::ArrayNameToString(fReport->GetListOfKeys(),"ncl.*",",");
  AliTreeTrending::AddJSROOTHtmlLink(pFile,"Ncl status", prefix, Form("items=[%s]",items.Data()));
  items=AliTreeTrending::ArrayNameToString(fReport->GetListOfKeys(),"dcar.*",",");
  AliTreeTrending::AddJSROOTHtmlLink(pFile,"DCAr status", prefix, Form("items=[%s]",items.Data()));
  items=AliTreeTrending::ArrayNameToString(fReport->GetListOfKeys(),"dcaz.*",",");
  AliTreeTrending::AddJSROOTHtmlLink(pFile,"DCAz status", prefix, Form("items=[%s]",items.Data()));
  items=AliTreeTrending::ArrayNameToString(fReport->GetListOfKeys(),".*Eff.*",",");
  AliTreeTrending::AddJSROOTHtmlLink(pFile,"Efficiency", prefix, Form("items=[%s]",items.Data()));
  items=AliTreeTrending::ArrayNameToString(fReport->GetListOfKeys(),"(.*MIP.*|.*dEdx.*)",",");
  AliTreeTrending::AddJSROOTHtmlLink(pFile,"dEdx", prefix, Form("items=[%s]",items.Data()));
  fclose (pFile);
}

