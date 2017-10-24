/// \ingroup PWGPP/TPC/macros/
/// \brief  MC/Anchor raw QA comparison
///
/// \author marian  Ivanov marian.ivanov@cern.ch, sebastian.lehner@cern.ch
///
/*!
\code
   gSystem->AddIncludePath("-I$ALICE_ROOT/include/"); //couldn't add include path in .rootrc
  .L $AliPhysics_SRC/PWGPP/TPC/macros/tpcMCValidation.C+
  TString mcPeriod="LHC15k1a1";
  TString mcPass="passMC";
  TString anchorPeriod="LHC15o";
  TString anchorPass="pass3_lowIR_pidfix";
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  InitTPCMCValidation("LHC15k1a1","passMC","LHC15o", "pass3_lowIR_pidfix",0,0);
\endcode
*/
#include <TError.h>
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
TTree *treeMC;
std::vector <Double_t> cRange;
std::vector <Double_t> cRange2;
std::vector <Double_t> cRange5;
TString sCriteria("1-2*(varname_Warning):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)");
TString queryString, queryTitle; ///

void makeTPCMCAlarms(TTree *treeMC, Bool_t doCheck, Int_t verbose);

Bool_t InitTPCMCValidation(TString mcPeriod, TString mcPass, TString anchorPeriod, TString anchorPass, Int_t verbose,
                           Int_t doCheck);

void MakeReport(const char *outputDir);

void tpcMCValidation(const char *mcPeriod = "LHC15k1a1", const char *outputDir = "./");


/// function to create a set of the comparison plots MC/AnchorRaw data
/// \param mcPeriod    -  MC production name
/// \param outputDir   -  directory where png and html files are stored
void tpcMCValidation(const char *mcPeriod, const char *outputDir) {
  cout << "INITIALIZING TPC MC Validation" << endl;
  AliExternalInfo i;
  cout << mcPeriod << endl;
  TString anchorProdNamePass = i.GetMCPassGuess(TString::Format("%s", mcPeriod));
  cout << "Anchor Production Name and Pass: " << anchorProdNamePass << endl;
  TObjArray *subStrL;
  subStrL = TPRegexp("^([^ ]+)").MatchS(anchorProdNamePass);
  TString anchorProdName = ((TObjString *) subStrL->At(0))->GetString();
  subStrL = TPRegexp("([^ ])+$").MatchS(anchorProdNamePass);
  TString anchorPassName = ((TObjString *) subStrL->At(0))->GetString();
  if (InitTPCMCValidation(mcPeriod, "passMC", anchorProdName, anchorPassName, 0, 0)) {
    MakeReport(outputDir);
  } else ::Error("tpcMCValidation", "InitTPCMCValidation returned with error -> skip plotting!");
}

/// makeTPCMCAlarms - Set of expression aliases/formulas  (Warning,Outlier,PhysAcc)
/// define variables for alarms for MC/raw mismatch
/// * 1.) Absolute aliases
/// * 2.) Relative aliases
/// * 3.) Combined status aliases
/// \param treeMC  - input tree
/// \param doCheck - force check of the variables
/// \param verbose - set verbosity for make alarms
void makeTPCMCAlarms(TTree * treeMC, Bool_t doCheck,Int_t verbose){
  //  TODO - Sebastian - add status bar for the dEdx - similar like int the standard QA
  //  TODO - adding new combined status bar status plots of individual alarms to be provided
  //  TODO - continue
  //  ==============================================================
  //  1.)  Partial alarms  (variable, variableAnchor deltaWarning,deltaError, PhysAcc)
  //                 deltaWarning and deltaError can be an expression which is understood by TTreeFormula
  //  ==============================================================
  // 1.) Absolute aliases
  TString sTrendVars=";";
  {
    // Ncl
    sTrendVars+="QA.TPC.meanTPCncl,TPC.Anchor.meanTPCncl,10,20,5;";       // delta Ncl  warning 10 ,  error 20     (nominal ~ 100-140)
    sTrendVars+="QA.TPC.meanTPCnclF,TPC.Anchor.meanTPCnclF,0.05,0.10,0.05;"; // delta NclF  warning 5%,  error 10%    (nominal ~ 90%)
    // dcaR resolution
    sTrendVars+="QA.TPC.dcarAP0,TPC.Anchor.dcarAP0,0.02,0.05,0.02;";     // dcarAP0;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcarCP0,TPC.Anchor.dcarCP0,0.02,0.05,0.02;";     // dcarCP0;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcarAP1,TPC.Anchor.dcarAP1,0.02,0.05,0.02;";     // dcarAP1;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcarCP1,TPC.Anchor.dcarCP1,0.02,0.05,0.02;";     // dcarCP1;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
//    sTrendVars+="QA.TPC.dcarCP0,TPC.Anchor.dcarCP0,0.02,0.05,0.02;";     // dcarCP0;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
//    sTrendVars+="QA.TPC.dcarAP1,TPC.Anchor.dcarAP1,0.02,0.05,0.02;";     // dcarAP1;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
//    sTrendVars+="QA.TPC.dcarCP1,TPC.Anchor.dcarCP1,0.02,0.05,0.02;";     // dcarCP1;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
    // Eff ITS: TPC->ITS
    sTrendVars+="QA.ITS.EffoneSPDPt02,ITS.Anchor.EffoneSPDPt02,0.05,0.1,0.07;";
    sTrendVars+="QA.ITS.EffoneSPDPt1,ITS.Anchor.EffoneSPDPt1,0.05,0.1,0.07;";
    sTrendVars+="QA.ITS.EffoneSPDPt10,ITS.Anchor.EffoneSPDPt10,0.05,0.1,0.07;";
    sTrendVars+="QA.ITS.EffTOTPt02,ITS.Anchor.EffTOTPt02,0.05,0.1,0.07;";
    sTrendVars+="QA.ITS.EffTOTPt1,ITS.Anchor.EffTOTPt1,0.05,0.1,0.07;";
    sTrendVars+="QA.ITS.EffTOTPt10,ITS.Anchor.EffTOTPt10,0.05,0.1,0.07;";
    // dEdX
    sTrendVars+="QA.TPC.meanMIP,TPC.Anchor.meanMIP,1,2,1;";     // meanMIP;  warning 1; error 2; physics acceptable 1; (nominal ~ 50)
    sTrendVars+="QA.TPC.resolutionMIP,TPC.Anchor.resolutionMIP,0.02,0.05,0.02;";     // resolutionMIP;  warning 0.01; error 0.02; physics acceptable 0.01; (nominal ~ 0.06)
    sTrendVars+="QA.TPC.MIPattachSlopeA,TPC.Anchor.MIPattachSlopeA,0.02,0.05,0.02;";     // MIPattachSlopeA;  warning xxx; error xxx  (nominal ~ 0)
    sTrendVars+="QA.TPC.MIPattachSlopeC,TPC.Anchor.MIPattachSlopeC,0.02,0.05,0.02;";     // MIPattachSlopeC;  warning xxx; error xxx  (nominal ~ 0)
    sTrendVars+="QA.TPC.meanMIPele,TPC.Anchor.meanMIPele,0.02,0.05,0.02;";     // MIPattachSlopeA;  warning xxx; error xxx  (nominal ~ 0)
    sTrendVars+="QA.TPC.resolutionMIPele,TPC.Anchor.resolutionMIPele,0.02,0.05,0.02;";     // MIPattachSlopeC;  warning xxx; error xxx  (nominal ~ 0)  
 
  }
  TStatToolkit::MakeAnchorAlias(treeMC,sTrendVars, doCheck, verbose);

  // 2.) Make aliases for Mean of MC-anchor of variables
  treeMC->SetAlias("diff0.meanTPCncl" , "(QA.TPC.meanTPCncl-TPC.Anchor.meanTPCncl)");
  treeMC->SetAlias("diff0.meanTPCnclF" ,"(QA.TPC.meanTPCnclF-TPC.Anchor.meanTPCnclF)");
  treeMC->SetAlias("ratio.dcarAP0" ,    "(QA.TPC.dcarAP0/TPC.Anchor.dcarAP0)");
  treeMC->SetAlias("ratio.dcarAP1" ,    "(QA.TPC.dcarAP1/TPC.Anchor.dcarAP1)");
  treeMC->SetAlias("ratio.dcarCP0" ,    "(QA.TPC.dcarCP0/TPC.Anchor.dcarCP0)");
  treeMC->SetAlias("ratio.dcarCP1" ,    "(QA.TPC.dcarCP1/TPC.Anchor.dcarCP1)");
  
  treeMC->SetAlias("ratio.meanMIP" ,    "(QA.TPC.meanMIP/TPC.Anchor.meanMIP)");
  treeMC->SetAlias("ratio.resolutionMIP" ,    "(QA.TPC.resolutionMIP/TPC.Anchor.resolutionMIP)");
  treeMC->SetAlias("diff0.MIPattachSlopeA" ,    "(QA.TPC.MIPattachSlopeA-TPC.Anchor.MIPattachSlopeA)");
  treeMC->SetAlias("diff0.MIPattachSlopeC" ,    "(QA.TPC.MIPattachSlopeC-TPC.Anchor.MIPattachSlopeC)");
  treeMC->SetAlias("diff0.meanMIPele" ,    "(QA.TPC.meanMIPele-TPC.Anchor.meanMIPele)"); 
  treeMC->SetAlias("diff0.resolutionMIPele" ,    "(QA.TPC.resolutionMIPele-TPC.Anchor.resolutionMIPele)");
  
  treeMC->SetAlias("statisticOK", "(meanTPCncl>0)");
  TString sDiffVars="";
  sDiffVars+="diff0.meanTPCncl;diff0.meanTPCnclF;ratio.dcarAP0;ratio.dcarAP1;ratio.dcarCP0;ratio.dcarCP1;ratio.meanMIP;ratio.resolutionMIP;diff0.MIPattachSlopeA;diff0.MIPattachSlopeC;diff0.meanMIPele;diff0.resolutionMIPele;";
  TObjArray* oaTrendVars = sDiffVars.Tokenize(";");
  Float_t entryFraction=0.8, nSigmaOutlier=6., nSigmaWarning=3., epsilon=1.0e-6, rangeFactor=0.1;
  for (Int_t iVar=0; iVar<oaTrendVars->GetEntriesFast(); iVar++) {
    TString sVar( oaTrendVars->At(iVar)->GetName() );
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_OutlierMin:(MeanEF-%f*RMSEF-%f):%f", nSigmaOutlier, epsilon, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_OutlierMax:(MeanEF+%f*RMSEF+%f):%f", nSigmaOutlier, epsilon, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_WarningMin:(MeanEF-%f*RMSEF-%f):%f", nSigmaWarning, epsilon, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_WarningMax:(MeanEF+%f*RMSEF+%f):%f", nSigmaWarning, epsilon, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f",   rangeFactor, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f",   rangeFactor, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_RobustMean:(MeanEF+0):%f", entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_Outlier:(varname>varname_OutlierMax||varname<varname_OutlierMin)"));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_Warning:(varname>varname_WarningMax||varname<varname_WarningMin)"));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_PhysAcc:(varname>varname_PhysAccMin&&varname<varname_PhysAccMax)"));
  }

  // 3.) Configure combined status. Default using logical OR of problems
  TString sCombinedStatus=";";
  sCombinedStatus+="ncl,absDiff.QA.TPC.meanTPCncl,absDiff.QA.TPC.meanTPCnclF;diff0.QA.TPC.meanTPCncl,diff0.QA.TPC.meanTPCnclF;"; // Status number of clusters and findable clusters
  sCombinedStatus+="dcarResol,absDiff.QA.TPC.dcarAP0,absDiff.QA.TPC.dcarAP1,absDiff.QA.TPC.dcarCP0,absDiff.QA.TPC.dcarCP1,ratio.dcarAP0,ratio.dcarAP1;";  // Status: DCA resolution
  sCombinedStatus+="itsEffStatus,absDiff.QA.ITS.EffoneSPDPt02,absDiff.QA.ITS.EffoneSPDPt1,absDiff.QA.ITS.EffoneSPDPt10,absDiff.QA.ITS.EffTOTPt02,absDiff.QA.ITS.EffTOTPt1,absDiff.QA.ITS.EffTOTPt10;";
  sCombinedStatus+="dEdX,ratio.meanMIP,ratio.resolutionMIP,diff0.MIPattachSlopeA,diff0.MIPattachSlopeC,diff0.meanMIPele,diff0.resolutionMIPele,absDiff.QA.TPC.meanMIP,absDiff.QA.TPC.resolutionMIP,absDiff.QA.TPC.MIPattachSlopeA,absDiff.QA.TPC.MIPattachSlopeC,absDiff.QA.TPC.meanMIPele,absDiff.QA.TPC.resolutionMIPele;";

  // Status: ITS:TPC-ITS matching efficiency
  TStatToolkit::MakeCombinedAlias(treeMC,sCombinedStatus,doCheck, verbose);
  ::Info("InitTPCMCValidation","Done with aliases");
}
///
/// \param mcPeriod
/// \param mcPass
/// \param anchorPeriod
/// \param anchorPass
/// \param verbose
/// \param doCheck
/// \return
Bool_t InitTPCMCValidation(TString mcPeriod, TString mcPass, TString anchorPeriod, TString anchorPass, Int_t verbose,
                           Int_t doCheck) {

  cRange = std::vector < Double_t > {0.13, 0.01, 0.5, 0.35};
  cRange2 = std::vector < Double_t > {0.13, 0.01, 0.5, 0.3};
  cRange5 = std::vector < Double_t > {0.13, 0.01, 0.8, 0.3};
  externalInfo = new AliExternalInfo(".", "", verbose);
  trendingDraw = new AliTreeTrending;
  trendingDraw->SetDefaultStyle();
  gStyle->SetOptTitle(0);

  treeMC = externalInfo->GetTree("QA.TPC", mcPeriod, mcPass, "QA.TPC;QA.TRD;QA.TOF;QA.ITS");

  TTree *treeAnchorTPC = NULL;
  treeAnchorTPC = externalInfo->GetTree("QA.TPC", anchorPeriod, anchorPass,
                                        "Logbook;QA.EVS;Logbook.detector:TPC:detector==\"TPC\"");
  if (treeAnchorTPC == NULL) {
    ::Error("InitTPCMCValidation", "Failed to get QA.TPC tree");
    return kFALSE;
  }

  TTree *treeAnchorTRD0;
  if (externalInfo->GetTree("QA.TRD", anchorPeriod, anchorPass, "Logbook") != 0) {
    treeAnchorTRD0 = externalInfo->GetTree("QA.TRD", anchorPeriod, anchorPass, "Logbook");
    ::Info("InitTPCMCValidation", "QA.TRD tree entries: %d", Int_t(treeAnchorTRD0->GetEntries()));
  } else {
    ::Error("InitTPCMCValidation", "Failed to get QA.TRD tree");
    return kFALSE;
  }

  TTree *treeAnchorITS0;
  if (externalInfo->GetTree("QA.ITS", anchorPeriod, anchorPass, "Logbook") != 0) {
    treeAnchorITS0 = externalInfo->GetTree("QA.ITS", anchorPeriod, anchorPass, "Logbook");
    ::Info("InitTPCMCValidation", "QA.ITS tree entries: %d", Int_t(treeAnchorITS0->GetEntries()));
  } else {
    ::Error("InitTPCMCValidation", "Failed to get QA.ITS tree");
    return kFALSE;
  }
  treeMC->AddFriend(treeAnchorTPC, "TPC.Anchor");
  treeMC->AddFriend(treeAnchorTRD0, "TRD.Anchor");
  treeMC->AddFriend(treeAnchorITS0, "ITS.Anchor");
  treeMC->SetAlias("QA.TPC.nEvents", "QA.TPC.entriesVertX");
  treeMC->SetAlias("TPC.Anchor.nEvents", "TPC.Anchor.entriesVertX");
  // check the match between MC and MC anchor
  {
    Int_t entriesMatch = treeMC->Draw("QA.TPC.meanTPCncl-TPC.Anchor.meanTPCncl", "1");
    Int_t entriesMC = treeMC->Draw("QA.TPC.meanTPCncl", "1");
    Int_t entriesAnchor = treeMC->Draw("TPC.Anchor.meanTPCncl", "1");
    if (entriesMatch == 0) {
      ::Error("InitTPCMCValidation", "No run match between the MC and and anchor raw");
      ::Error("InitTPCMCValidation", "QA runs:\tMC=%d\tAnchor=%d\tMatch=%d", entriesMC, entriesAnchor, entriesMatch);
      return kFALSE;
    } else {
      ::Info("InitTPCMCValidation", "QA runs:\tMC=%d\tAnchor=%d\tMatch=%d", entriesMC, entriesAnchor, entriesMatch);
    }
  }
  //
  makeTPCMCAlarms(treeMC, doCheck, verbose);
  TString sStatusBarVars("ncl;dcarResol;itsEffStatus;dEdX");
  TString sStatusBarNames("#(cl);dcar;itsEffStatus;dEdX;");
  TString sCriteria("(1):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)"); // status bar markers
  TString statusString[3];
  statusString[0] = sStatusBarVars;
  statusString[1] = sStatusBarNames;
  statusString[2] = sCriteria;

  trendingDraw->AddUserDescription(new TNamed("MC period", mcPeriod.Data()));
  trendingDraw->AddUserDescription(new TNamed("MC pass", mcPass.Data()));
  trendingDraw->AddUserDescription(new TNamed("Anchor period", anchorPeriod.Data()));
  trendingDraw->AddUserDescription(new TNamed("Anchor pass", anchorPass.Data()));
  trendingDraw->SetTree(treeMC);
  treeMC->SetAlias("tagID", "run");
  treeMC->SetAlias("defaultCut", "run==TPC.Anchor.run");

  if (trendingDraw->InitSummaryTrending(statusString, 0.015, "defaultCut")) return kTRUE;
  else return kFALSE;
}

/// ## MakeReport
/// \param outputDir
void MakeReport(const char *outputDir) {
  TMultiGraph *graph = 0, *lines = 0;
  TString queryString = "";
  //outputDir=".";
  trendingDraw->fWorkingCanvas->Print(TString(outputDir) + "/report.pdf[", "pdf");
  //
  // DONE: In some cases different sources provided different run lists
  // DONE: MakeMultGraph should match "run graphs" and "rebinning graphs" if needed (Marian - in TStatToolkit)
  // DONE: Remove CLion convention warning (Marian)
  // TODO: Optionally add offset to the X value in order to distinguish overlapped  graphs TStatToolkit::
  // TODO: Partially done. Add TMultiGraph "class?" to specify additional options (see bellow) (Boris, Marian) - groupName to be used for that
  // TODO: Optionally draw y value on top of the markers (for single graphs) can be coded in the class (lower priority)
  // TODO: For each tab we should provide status figure (decomposition of status to the components)

  // TODO: For MC, Anchor and for MC/Anchor comparison (Marian, Sebastian)

  // 1.) Event properties  ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabEvent.html)
  trendingDraw->MakePlot(outputDir, "interactionRate.png", "Interaction rate", cRange, "",
                         "Logbook.averageEventsPerSecond;QA.EVS.interactionRate:run", "defaultCut", "figTemplateTRD",
                         "figTemplateTRD", 1, 1, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "runDuration.png", "Run duration ", cRange, "", "Logbook.runDuration:run",
                         "defaultCut", "figTemplateTRD", "figTemplateTRD", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "bField.png", "Magnet current", cRange, "", "Logbook.L3_magnetCurrent:run",
                         "defaultCut", "figTemplateTRD", "figTemplateTRD", 1, 0.75, 3, kTRUE);
  trendingDraw->MakePlot(outputDir, "eventCounters.png", "Event counters ", cRange, "",
                         "Logbook.totalEvents;totalEventsPhysics;totalEventsCalibration;Logbook.detector_TPC.eventCountPhysics;TPC.Anchor.nEvents;QA.TPC.nEvents:run",
                         "defaultCut", "figTemplateTRD", "figTemplateTRD", 1, 0.75, 4, kTRUE);
  //trendingDraw->MakePlot(outputDir,"runTime.png","Time ",cRange,"",":run","defaultCut","figTemplateTRDPair","figTemplateTRDPair",1,0.75,4,kTRUE); //TODO - add time format
  trendingDraw->MakePlot(outputDir, "meanMult.png", "Mean TPC multiplicity (|DCA|<3 cm)", cRange, "",
                         "QA.TPC.meanMult;TPC.Anchor.meanMult:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanMultPos.png", "Mean TPC multiplicity (q>0, |DCA|<3 cm)", cRange, "",
                         "QA.TPC.meanMultPos;TPC.Anchor.meanMultPos:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanMultNeg.png", "Mean TPC multiplicity (q<0, |DCA|<3 cm)", cRange, "",
                         "QA.TPC.meanMultNeg;TPC.Anchor.meanMultNeg:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanMult_comb2.png", "meanMult_comb2", cRange, "",
                         "meanMult_comb2;TPC.Anchor.meanMult_comb2:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanVertX.png", "meanVertX:run MC/Anchor", cRange, "",
                         "QA.TPC.meanVertX;TPC.Anchor.meanVertX:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 5, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanVertY.png", "meanVertY:run MC/Anchor", cRange, "",
                         "QA.TPC.meanVertY;TPC.Anchor.meanVertY:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 5, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanVertZ.png", "meanVertZ:run MC/Anchor", cRange, "",
                         "QA.TPC.meanVertZ;TPC.Anchor.meanVertZ:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 5, kTRUE);
  //
  // 2.) Number of clusters comparison ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabNcl.html)
  // TODO: Add all estimators fo missing chambers (Ncl, Voltage, RawQA, tracks)
  trendingDraw->MakeStatusPlot(outputDir, "nclStatus.png",
                               "absDiff.QA.TPC.meanTPCncl;absDiff.QA.TPC.meanTPCnclF;diff0.meanTPCncl;diff0.meanTPCnclF;run",
                               "#Delta^{A}_{ATPCncl};#Delta^{A}_{TPCnclF};#Delta^{R}_{TPCncl};#Delta^{R}_{TPCnclF};",
                               "defaultCut", sCriteria);
  trendingDraw->MakePlot(outputDir, "meanTPCncl.png", "Number of clusters", cRange, "",
                         "QA.TPC.meanTPCncl;TPC.Anchor.meanTPCncl:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 5, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanTPCnclFindable.png", "Cluster fraction #left(#frac{N_{cl}}{N_{find.}}#right)",
                         cRange, "", "QA.TPC.meanTPCnclF;TPC.Anchor.meanTPCnclF:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 5, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanTPCNclRatioMCtoAnchor.png", "Number of clusters MC/Anchor", cRange, "",
                         "meanTPCncl/TPC.Anchor.meanTPCncl;meanTPCnclF/TPC.Anchor.meanTPCnclF:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 5, kTRUE);
  trendingDraw->MakePlot(outputDir, "iroc.png", "IROC #chambers", cRange, "",
                         "iroc_A_side;TPC.Anchor.iroc_A_side;iroc_C_side;TPC.Anchor.iroc_C_side:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "oroc.png", "OROC #chambers", cRange, "",
                         "oroc_A_side;TPC.Anchor.oroc_A_side;oroc_C_side;TPC.Anchor.oroc_C_side:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  //
  // 3.) Matching efficiency ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabEff.html)
  // TODO: provide description of variables (Metadata attribute Description exist - where to place? )
  // TODO: Add ITS layer matching trending
  // TODO: Add TRD matching
  trendingDraw->MakePlot(outputDir, "matchingTPCITSEffNoPileUpCut.png", "Matching efficiency(no pileup cut):MC/Anchor",
                         cRange, "",
                         "QA.TPC.tpcItsMatchA;TPC.Anchor.tpcItsMatchA;QA.TPC.tpcItsMatchC;TPC.Anchor.tpcItsMatchC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPCITSEffPileUpCut.png", "Matching efficiency (pileup cut):MC/Anchor",
                         cRange, "",
                         "QA.ITS.EffTOTPt02;ITS.Anchor.EffTOTPt02;QA.ITS.EffTOTPt1;ITS.Anchor.EffTOTPt1:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPCITSEffPileUpCutHighPt.png",
                         "Matching efficiency (pileup cut):MC/Anchor", cRange, "",
                         "QA.ITS.EffTOTPt1;ITS.Anchor.EffTOTPt1;QA.ITS.EffTOTPt10;ITS.Anchor.EffTOTPt10:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPCITSEffACRatio.png", "Matching efficiency A/C ratio:MC/Anchor", cRange,
                         "",
                         "QA.TPC.tpcItsMatchA/QA.TPC.tpcItsMatchC;TPC.Anchor.tpcItsMatchA/TPC.Anchor.tpcItsMatchC;QA.TPC.tpcItsMatchHighPtA/QA.TPC.tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtA/TPC.Anchor.tpcItsMatchHighPtC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPCTRDEffPileUpCut.png", "Matching efficiency (pileup cut):MC/Anchor",
                         cRange, "",
                         "QA.TRD.TPCTRDmatchEffPosAll;TRD.Anchor.TPCTRDmatchEffPosAll;QA.TRD.TPCTRDmatchEffNegAll;TRD.Anchor.TPCTRDmatchEffNegAll:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  // TODO plots bellow not in the html (to be described or removed)
  trendingDraw->MakePlot(outputDir, "matchingTPC-ITSEff.png", "Matching efficiency:MC/Anchor", cRange, "",
                         "QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffTOTPt02;QA.ITS.EffTOTPt1;QA.ITS.EffTOTPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffTOTPt02;ITS.Anchor.EffTOTPt1;ITS.Anchor.EffTOTPt10:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcItsMatch.png", "TPC - ITS Match", cRange, "",
                         "tpcItsMatchHighPtA;TPC.Anchor.tpcItsMatchHighPtA;tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPC-ITSEff_1.png", "Matching efficiency x", cRange, "",
                         "QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffoneSPDPt02;QA.ITS.EffoneSPDPt1;QA.ITS.EffoneSPDPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffoneSPDPt02;ITS.Anchor.EffoneSPDPt1;ITS.Anchor.EffoneSPDPt10:run",
                         "defaultCut", "21;24;25;27;28;21;24;25;27;28", "2;2;2;figTemplateTRDPair;4;4;4", 1, 0.75, 6,
                         kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPC-ITSEff_2.png", "Matching efficiency y", cRange, "",
                         "QA.TPC.tpcItsMatchA/TPC.Anchor.tpcItsMatchA;QA.TPC.tpcItsMatchC/TPC.Anchor.tpcItsMatchC;QA.ITS.EffTOTPt02/ITS.Anchor.EffTOTPt02;QA.ITS.EffTOTPt1/ITS.Anchor.EffTOTPt1;QA.ITS.EffTOTPt10/ITS.Anchor.EffTOTPt10:run",
                         "defaultCut", "21;24;25;27;28;21;24;25;27;28", "1;2;4;3;6;2;4;4;4;4;4", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcItsMatch.png", "tpcItsMatch", cRange, "",
                         "tpcItsMatchA;TPC.Anchor.tpcItsMatchA;tpcItsMatchC;TPC.Anchor.tpcItsMatchC:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcItsMatchHighPt.png", "tpcItsMatchHighPt", cRange, "",
                         "tpcItsMatchHighPtA;TPC.Anchor.tpcItsMatchHighPtA;tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  //
  // 4.) DCA  ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabDCA.html)
  //
  trendingDraw->MakePlot(outputDir, "dcarP0.png",
                         "HighPt: DCA_{xy} #sigma_{0} (#sigma^{2}=#sigma_{0}^{2}+#sigma_{1}^{2}/p_{T}^{2}) ", cRange,
                         "", "QA.TPC.dcarAP0;TPC.Anchor.dcarAP0;QA.TPC.dcarCP0;TPC.Anchor.dcarCP0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcarP1.png",
                         "MS: DCA_{xy} #sigma_{1} (#sigma^{2}=#sigma_{0}^{2}+#sigma_{1}^{2}/p_{T}^{2})", cRange, "",
                         "QA.TPC.dcarAP1;TPC.Anchor.dcarAP1;QA.TPC.dcarCP1;TPC.Anchor.dcarCP1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(0, "offsetR.png", "Offset dR", cRange, "",
                         "offsetdRA;TPC.Anchor.offsetdRA;offsetdRC;TPC.Anchor.offsetdRC:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->AppendBand(outputDir, "offsetdR.png",
                           "offsetdRA_RobustMean;offsetdRA_OutlierMin;offsetdRA_OutlierMax;offsetdRA_WarningMin;offsetdRA_WarningMax:run",
                           "defaultCut", "1;1;1;1;1,1;2;2;3;3", "1;1;1;1;1,1;figTemplateTRDPair", 1, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_A_0.png", "dcar_A_0", cRange, "",
                         "dcar_posA_0;TPC.Anchor.dcar_posA_0;dcar_negA_0;TPC.Anchor.dcar_negA_0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_C_0.png", "dcar_C_0", cRange, "",
                         "dcar_posC_0;TPC.Anchor.dcar_posC_0;dcar_negC_0;TPC.Anchor.dcar_negC_0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_C_1.png", "dcar_C_1", cRange, "",
                         "dcar_posC_1;TPC.Anchor.dcar_posC_1;dcar_negC_1;TPC.Anchor.dcar_negC_1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_A_1.png", "dcar_A_1", cRange, "",
                         "dcar_posA_1;TPC.Anchor.dcar_posA_1;dcar_negA_1;TPC.Anchor.dcar_negA_1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_A_2.png", "dcar_A_2", cRange, "",
                         "dcar_posA_2;TPC.Anchor.dcar_posA_2;dcar_negA_2;TPC.Anchor.dcar_negA_2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_C_2.png", "dcar_C_2", cRange, "",
                         "dcar_posC_2;TPC.Anchor.dcar_posC_2;dcar_negC_2;TPC.Anchor.dcar_negC_2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_A_0.png", "dcaz_A_0", cRange, "",
                         "dcaz_posA_0;TPC.Anchor.dcaz_posA_0;dcaz_negA_0;TPC.Anchor.dcaz_negA_0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_C_0.png", "dcaz_C_0", cRange, "",
                         "dcaz_posC_0;TPC.Anchor.dcaz_posC_0;dcaz_negC_0;TPC.Anchor.dcaz_negC_0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_C_1.png", "dcaz_C_1", cRange, "",
                         "dcaz_posC_1;TPC.Anchor.dcaz_posC_1;dcaz_negC_1;TPC.Anchor.dcaz_negC_1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_A_1.png", "dcaz_A_1", cRange, "",
                         "dcaz_posA_1;TPC.Anchor.dcaz_posA_1;dcaz_negA_1;TPC.Anchor.dcaz_negA_1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_A_2.png", "dcaz_A_2", cRange, "",
                         "dcaz_posA_2;TPC.Anchor.dcaz_posA_2;dcaz_negA_2;TPC.Anchor.dcaz_negA_2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_C_2.png", "dcaz_C_2", cRange, "",
                         "dcaz_posC_2;TPC.Anchor.dcaz_posC_2;dcaz_negC_2;TPC.Anchor.dcaz_negC_2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcarFitPar_comb4.png", "dcarFitpar_comb4", cRange, "",
                         "dcarFitpar_comb4;TPC.Anchor.dcarFitpar_comb4:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  //trendingDraw->MakePlot(outputDir,"rmsDCAMultParMCtoAnchor.png","DCA Resolution mult due MS:MC/Anchor",cRange,"","dcarAP1;dcarCP1;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP1:run","defaultCut","figTemplateTRDPair","figTemplateTRDPair",1,0.75,6,kTRUE);
  //
  // 5.) dEdx ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabdEdx.html)
  // TODO: Add sub-status lines
  trendingDraw->MakePlot(outputDir, "meanMIP.png", "<Mean dEdx_{MIP}> (a.u) ", cRange, "",
                         "QA.TPC.meanMIP;TPC.Anchor.meanMIP:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanElectron.png", "<dEdx_{el}> (a.u) ", cRange, "",
                         "QA.TPC.meanMIPele;TPC.Anchor.meanMIPele:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "electronMIPSeparation.png", "<dEdx_{el}>-<dEdx_{MIP}>", cRange, "",
                         "QA.TPC.electroMIPSeparation;TPC.Anchor.electroMIPSeparation:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "electronMIPRatio.png", "<dEdx_{el}>/<dEdx_{MIP}>", cRange, "",
                         "QA.TPC.meanMIPele/QA.TPC.meanMIP;TPC.Anchor.meanMIPele/TPC.Anchor.meanMIP:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "MIPattachSlope.png", "MIPattachSlope", cRange, "",
                         "MIPattachSlopeA;TPC.Anchor.MIPattachSlopeA;MIPattachSlopeC;TPC.Anchor.MIPattachSlopeC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "MIPattachSlopeA.png", "MIPattachSlopeA", cRange, "",
                         "MIPattachSlopeA;TPC.Anchor.MIPattachSlopeA:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "MIPattachSlopeC.png", "MIPattachSlopeC", cRange, "",
                         "MIPattachSlopeC;TPC.Anchor.MIPattachSlopeC:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 4, kTRUE);
  //
  queryString = "QA.TPC.meanMIP_Warning*1.1;QA.TPC.resolutionMIP_Warning*1.2;QA.TPC.MIPattachSlopeA_Warning*1.3;";
  queryString += "QA.TPC.MIPattachSlopeC_Warning*1.4;QA.TPC.meanMIPele_Warning*1.5;QA.TPC.resolutionMIPele_Warning*1.6;QA.TPC.electroMIPSeparation_Warning*1.7:run";
  trendingDraw->MakePlot(outputDir, "dEdxStatus.png", "dEdx status", cRange, "", queryString.Data(), "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);

  { // MIP resolution  and separation
    trendingDraw->fWorkingCanvas->Clear();
    TLegend *legend = new TLegend(cRange[0], cRange[1], cRange[2], cRange[3], "DCA Resolution mult due MS:MC/Anchor");
    legend->SetBorderSize(0);
    graph = TStatToolkit::MakeMultGraph(treeMC, "", "meanMIPele/meanMIP;TPC.Anchor.meanMIPele/TPC.Anchor.meanMIP:run",
                                        "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, legend);
    lines = TStatToolkit::MakeMultGraph(treeMC, "",
                                        "meanMIPele_RobustMean/meanMIP_RobustMean;TPC.Anchor.meanMIPele_RobustMean/TPC.Anchor.meanMIP_RobustMean:run",
                                        "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, 0);
    if (!graph) {
      legend = new TLegend(cRange[0], cRange[1], cRange[2], cRange[3], "Plotting error!!");
      ::Error("tpcMCValidation", "No plot returned -> dummy plot!");
    } else {
      TStatToolkit::DrawMultiGraph(graph, "alp");
      trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
    }
    TStatToolkit::DrawMultiGraph(lines, "l");
    legend->SetFillStyle(0);
    legend->Draw();
    trendingDraw->fWorkingCanvas->SaveAs(TString(outputDir) + "/mipToEleSeparation.png");
    trendingDraw->fWorkingCanvas->Print(TString(outputDir) + "/report.pdf");
  }

  //
  // 6.) matching (tabMatching.html))
  //
  trendingDraw->MakePlot(outputDir, "offsetdZ.png", "offsetdZ", cRange, "",
                         "offsetdZA;TPC.Anchor.offsetdZA;offsetdZC;TPC.Anchor.offsetdZC:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "offsetd_comb4.png", "offsetd_comb4", cRange, "",
                         "offsetd_comb4;TPC.Anchor.offsetd_comb4:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "lambdaPull.png", "lambdaPull", cRange, "", "lambdaPull;TPC.Anchor.lambdaPull:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcConstrainPhiC.png", "tpcConstrainPhiC", cRange, "",
                         "tpcConstrainPhiC;TPC.Anchor.tpcConstrainPhiC:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "yPull.png", "yPull", cRange, "", "yPull;TPC.Anchor.yPull:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "zPull.png", "zPull", cRange, "", "zPull;TPC.Anchor.zPull:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "itsTpcPulls_comb4.png", "itsTpcPulls_comb4", cRange, "",
                         "itsTpcPulls_comb4;TPC.Anchor.itsTpcPulls_comb4:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcConstrainPhi.png", "tpcConstrainPhi", cRange, "",
                         "tpcConstrainPhiA;TPC.Anchor.tpcConstrainPhiA;tpcConstrainPhiC;TPC.Anchor.tpcConstrainPhiC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcConstrainPhi_comb2.png", "tpcConstrainPhi_comb2", cRange, "",
                         "tpcConstrainPhi_comb2;TPC.Anchor.tpcConstrainPhi_comb2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "deltaPt.png", "deltaPt", cRange, "", "deltaPt;TPC.Anchor.deltaPt:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "deltaPtAC.png", "deltaPtAC", cRange, "",
                         "deltaPtA;TPC.Anchor.deltaPtA;deltaPtC;TPC.Anchor.deltaPtC:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 0.75, 6, kTRUE);


  // Status plots
  TString statusExpression,statusTitle;

  // DCA status
  statusExpression="absDiff.QA.TPC.dcarAP0;absDiff.QA.TPC.dcarCP0;absDiff.QA.TPC.dcarAP1;absDiff.QA.TPC.dcarCP1;";
  statusTitle="absDiff.QA.TPC.dcarAP0;absDiff.QA.TPC.dcarCP0;absDiff.QA.TPC.dcarAP1;absDiff.QA.TPC.dcarCP1;";
  statusExpression+="dcarAP0;dcarAP1;dcarCP0;dcarCP1;dcar_posA_0;dcar_posA_1;dcar_posA_2;dcar_posC_0;dcar_posC_1;dcar_posC_2;";
  statusTitle+="dcarAP0;dcarAP1;dcarCP0;dcarCP1;dcar_posA_0;dcar_posA_1;dcar_posA_2;dcar_posC_0;dcar_posC_1;dcar_posC_2;";
  statusExpression+="TPC.Anchor.dcarAP0;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP0;TPC.Anchor.dcarCP1;TPC.Anchor.dcar_posA_0;TPC.Anchor.dcar_posA_1;TPC.Anchor.dcar_posA_2;TPC.Anchor.dcar_posC_0;TPC.Anchor.dcar_posC_1;TPC.Anchor.dcar_posC_2;run";
  statusTitle+="TPC.Anchor.dcarAP0;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP0;TPC.Anchor.dcarCP1;TPC.Anchor.dcar_posA_0;TPC.Anchor.dcar_posA_1;TPC.Anchor.dcar_posA_2;TPC.Anchor.dcar_posC_0;TPC.Anchor.dcar_posC_1;TPC.Anchor.dcar_posC_2";
  trendingDraw->MakeStatusPlot("./", "dcarStatus.png", statusExpression, statusTitle, "defaultCut",sCriteria);

  //
  trendingDraw->fWorkingCanvas->Clear();
  trendingDraw->fWorkingCanvas->Print(TString(outputDir) + "/report.pdf]", "pdf");
}

///
void MakeDCAStatusPlot(){  // Example
   // Status plots
  TString statusExpression,statusTitle;
  // DCA status
  statusExpression="absDiff.QA.TPC.dcarAP0;absDiff.QA.TPC.dcarCP0;absDiff.QA.TPC.dcarAP1;absDiff.QA.TPC.dcarCP1;";
  statusTitle="#Delta^{MC-Anchor}_{dcarAP0};#Delta^{MC-Anchor}_{dcarCP0};#Delta^{MC-Anchor}_{dcarAP1};#Delta^{MC-Anchor}_{dcarCP1};";
  statusExpression+="ratio.dcarAP0;ratio.dcarCP0;ratio.dcarAP1;ratio.dcarCP1;";
  statusTitle+="#Delta^{MC/Anchor-#mu}_{dcarAP0};#Delta^{MC/Anchor-#mu}_{dcarCP0};#Delta^{MC/Anchor-#mu}_{dcarAP1};#Delta^{MC/Anchor-#mu}_{dcarCP1};";
  statusExpression+="dcarAP0;dcarAP1;dcarCP0;dcarCP1;";
  statusTitle+="MC_{dcarAP0};MC_{dcarAP1};MC_{dcarCP0};MC_{dcarCP1};";
  statusExpression+="TPC.Anchor.dcarAP0;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP0;TPC.Anchor.dcarCP1;";
  statusTitle+="Anchor_{dcarAP0};Anchor_{dcarAP1};Anchor_{dcarCP0};Anchor_{dcarCP1};";
  statusExpression+="run;";
  trendingDraw->MakeStatusPlot("./", "dcarStatus.png", statusExpression.Data(), statusTitle.Data(), "defaultCut",sCriteria);
}

///
void MakedEdXStatusPlot() {  //TODO - Sebastian
  // To get the list of warning for standard QA treeMC->GetListOfAliases()->Print("","*MIP*arning")
  // Status plots
  TString statusExpression,statusTitle;
  // MIP status

  statusExpression="ratio.meanMIP;ratio.resolutionMIP;diff0.MIPattachSlopeA;diff0.MIPattachSlopeC;diff0.resolutionMIPele;diff0.meanMIPele;absDiff.QA.TPC.meanMIP;absDiff.QA.TPC.resolutionMIP;absDiff.QA.TPC.MIPattachSlopeA;absDiff.QA.TPC.MIPattachSlopeC;";
  statusTitle="ratio-meanMIP;ratio-resolutionMIP;diff0.MIPattachSlopeA;diff0.MIPattachSlopeC;diff0.resolutionMIPele;diff0.meanMIPele;#Delta^{MC/Anchor-#mu}_{meanMIP};#Delta^{MC/Anchor-#mu}_{resolutionMIP};#Delta^{MC/Anchor-#mu}_{MIPattachSlopeA};#Delta^{MC/Anchor-#mu}_{MIPattachSlopeC};";
  statusExpression+="run;";
  trendingDraw->MakeStatusPlot("./", "dEdXStatus.png", statusExpression.Data(), statusTitle.Data(), "defaultCut",sCriteria);

}


/// TODO - use html hyperlinks + formatting
/// TODO  - bug fix in the html
/// TODO - alias for the html table names, hints
/// TODO - table name and hints   metadata tag ".thead", ".tooltip"
/// TODO - fix formatting in AliTreeFormulaF- in case of missing entry write undefined
///
/// TODO - tooltip - explaining status bits (detector mask, resp. status mask)
void makeHtml() {
  //
  //
  // EVS
  TStatToolkit::AddMetadata(treeMC, "interactionRate.thead", "rate(Hz)");
  TStatToolkit::AddMetadata(treeMC, "interactionRate.tooltip", "Interaction rate (source EVS)");
  TStatToolkit::AddMetadata(treeMC, "bz.thead", "B<sub>z</sub> (T)");
  TStatToolkit::AddMetadata(treeMC, "bz.tooltip", "Barrel field");
  TStatToolkit::AddMetadata(treeMC, "detectorMask.html","%b{detectorMask&0x7F}");
  TStatToolkit::AddMetadata(treeMC, "detectorMask.headerTooltip","Detector bit mask:\n* SPD\n* SDD\n* SSS\n* TPC\n* TRD");
}

void makeHtmlDCA(){
  //TString metaDCA="<a title
  TString pathMC="<a href=\"http://aliqatpc.web.cern.ch/aliqatpc/sim/%d{year}/%s{period.GetName()}/passMC/000%d{run}/";
  TString pathData="<a href=\"http://aliqatpc.web.cern.ch/aliqatpc/data/%d{TPC.Anchor.year}/%s{TPC.Anchor.period.GetName()}/%{TPC.Anchor.pass.GetName()}/000%d{run}/";
  TStatToolkit::AddMetadata(treeMC, "dcarAP0.html", (pathMC+"dca_and_phi.png\">%2.4f{dcarAP0}</a>").Data());
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarAP0.html", (pathData+"dca_and_phi.png\">%2.4f{TPC.Anchor.dcarAP0}</a>").Data());
  TStatToolkit::AddMetadata(treeMC, "dcarCP0.html", (pathMC+"dca_and_phi.png\">%2.4f{dcarCP0}</a>").Data());
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarAP0.html", (pathData+"dca_and_phi.png\">%2.4f{TPC.Anchor.dcarCP0}</a>").Data());

  TStatToolkit::AddMetadata(treeMC, "dcarAP0.thead", "&sigma;<sup>A</sup><sub>MC DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeMC, "dcarAP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeMC, "dcarCP0.thead", "&sigma;<sup>C</sup><sub>MC DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeMC, "dcarCP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarAP0.thead", "&sigma;<sup>A</sup><sub>Anchor DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarAP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarCP0.thead", "&sigma;<sup>C</sup><sub>Anchor DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarCP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarCP0.html", "<a href=\"http://aliqatpc.web.cern.ch/aliqatpc/data/%d{year}/%s{TPC.Anchor.period.GetName()}/passMC/000%d{run}/dca_and_phi.png\">%2.2f{dcarCP0}</a>");


  //
  TString runSelection = "defaultCut";
  TString varSelection = "";
  TString logbookBase = "";
  logbookBase = "run:LHCFillNumber:LHCperiod:detectorMask:HLTmode:DAQ_time_start:ctpDuration:totalEventsPhysics:interactionRate:bz;2.2:";
  // tabDCA
  TString tpcDCA = "dcarAP0:TPC.Anchor.dcarAP0:dcarCP0:TPC.Anchor.dcarCP0";//:dcarStatusString:dcarRawStatusString";
  // TString tpcDCA="dcarAP1;2.2f:dcarStatusString:dcarStatusOutlier:dcarStatusWarning";
  // TObjArray *tpcDCAArray = AliTreePlayer::selectMetadata(treeMC->GetFriend("QA.TPC"), "[class==\"TPC&&DCAr&&!Err&&!Chi2\"]",0);
  AliTreePlayer::selectWhatWhereOrderBy(treeMC, (logbookBase + tpcDCA).Data(), runSelection.Data(), "", 0, 100000,
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
  AliTreeFormulaF testLink("testLink",testHtml,treeMC,7);
  //testLink.fDebug=1;
  printf("\n %s \n",testLink.PrintValue());
}

/// html useful links
///  - Monalisa - raw run details query - used on mouse over
///     * http://alimonitor.cern.ch/raw/rawrun_details.jsp?run=280897
///  - multi-line tooltip using title tag  - can be used for bitmask explanation
///      http://jsfiddle.net/rzea/vsp6840b/3/
///  usage of datalist:
///   -https://www.w3schools.com/tags/tryit.asp?filename=tryhtml5_datalist
///   https://www.w3schools.com/tags/tryit.asp?filename=tryhtml5_details


/// TODO - automatic decomposition status aliases
/// TODO - get form logbook detailed per run information  (put into html metadata string)
/// TODO - Ful name of the alias expression (maybe recursive)
/// TODO - add html preview for table header
/// tableHeader, tableHeader_Tooltip, tableHeader_Title, tableHeader_html