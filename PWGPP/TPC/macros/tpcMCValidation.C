/// \ingroup PWGPP/TPC/macros/
/// \brief  MC/Anchor raw QA comparison
///
/// \author marian  Ivanov marian.ivanov@cern.ch, sebastian.lehner@cern.ch
///
/*!
\code
  gSystem->AddIncludePath("-I$ALICE_ROOT/include/"); //couldn't add include path in .rootr
  .L $AliPhysics_SRC/PWGPP/TPC/macros/tpcMCValidation.C+
  TString mcPeriod="LHC15k1a1";
  TString mcPass="passMC";
  TString anchorPeriod="LHC15o";
  TString anchorPass="pass3_lowIR_pidfix";
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  //
  InitTPCMCValidation("LHC15k1a1","passMC","LHC15o", "pass3_lowIR_pidfix",0,0); //short period
//  InitTPCMCValidation("LHC16g1c","passMC","LHC15o", "pass1",0,0);   // long period example
  //
  MakeReport();
 
  MakeStatusPlots();
  trendingDraw->fReport->Close();
  //
  makeHtml();

\endcode
*/

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
TTree *treeMC =0;
std::vector <Double_t> cRange;
std::vector <Double_t> cRange2;
std::vector <Double_t> cRange5;
TString sCriteria("(present):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)");
// current variables - used instead of local variables to enable copy paste of part of code
TString queryString, queryTitle; ///
TMultiGraph *graph = 0, *lines = 0;
const char *outputDir="./";

/// \param treeMC
/// \param doCheck
/// \param verbose
void makeTPCMCAlarms(TTree *treeMC, Bool_t doCheck, Int_t verbose);

Bool_t InitTPCMCValidation(TString mcPeriod, TString mcPass, TString anchorPeriod, TString anchorPass, Int_t verbose,
                           Int_t doCheck);

void MakeReport();
void MakeStatusPlots();
void tpcMCValidation(const char *mcPeriod = "LHC15k1a1", const char *soutputDir = "./");
void makeHtmlDCA();


/// function to create a set of the comparison plots MC/AnchorRaw data
/// \param mcPeriod    -  MC production name
/// \param outputDir   -  directory where png and html files are stored
void tpcMCValidation(const char *mcPeriod, const char *sOutputDir) {
  outputDir=sOutputDir;
  cout << "INITIALIZING TPC MC Validation" << endl;
  AliExternalInfo i;
  cout << mcPeriod << endl;
  TString anchorProdNamePass = i.GetMCPassGuess(TString::Format("%s", mcPeriod));
  cout << "Anchor Production Name and Pass: " << anchorProdNamePass << endl;
  if(anchorProdNamePass.Contains("not found")) {
      ::Error("tpcMCValidation", "MC not found in guess -> skip plotting!");
      return;
  }
  TObjArray *subStrL;
  subStrL = TPRegexp("^([^ ]+)").MatchS(anchorProdNamePass);
  TString anchorProdName = ((TObjString *) subStrL->At(0))->GetString().ReplaceAll("?","");
  subStrL = TPRegexp("([^ ])+$").MatchS(anchorProdNamePass);
  TString anchorPassName = ((TObjString *) subStrL->At(0))->GetString();
  if (InitTPCMCValidation(mcPeriod, "passMC", anchorProdName, anchorPassName, 0, 0)) {
    MakeReport();
  }
  else ::Error("tpcMCValidation", "InitTPCMCValidation returned with error -> skip plotting!");
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
  //  ==============================================================
  //  1.)  Partial alarms  (variable, variableAnchor deltaWarning,deltaError, PhysAcc)
  //                 deltaWarning and deltaError can be an expression which is understood by TTreeFormula
  //  ==============================================================
  // 1.) Absolute aliases  alarms on |MC-Anchor|
  TString sTrendVars=";";
  {
    // Ncl
    sTrendVars+="QA.TPC.meanTPCncl,TPC.Anchor.meanTPCncl,10,20,5;";       // delta Ncl  warning 10 ,  error 20     (nominal ~ 100-140)
    sTrendVars+="QA.TPC.meanTPCnclF,TPC.Anchor.meanTPCnclF,0.05,0.10,0.05;"; // delta NclF  warning 5%,  error 10%    (nominal ~ 90%)
    // dcaR resolution
    sTrendVars+="QA.TPC.dcarAP0,TPC.Anchor.dcarAP0,0.02,0.05,0.02;";     // dcar;  warning 0.02 cm; error 0.05 cm  (nominal ~ 0.02 cm)
    sTrendVars+="QA.TPC.dcarCP0,TPC.Anchor.dcarCP0,0.02,0.05,0.02;";     // dcar;  warning 0.02 cm; error 0.05 cm  (nominal ~ 0.02 cm)
    sTrendVars+="QA.TPC.dcarAP1,TPC.Anchor.dcarAP1,0.02,0.05,0.02;";     // dcar;  warning 0.02 cm; error 0.05 cm  (nominal ~ 0.02 cm)
    sTrendVars+="QA.TPC.dcarCP1,TPC.Anchor.dcarCP1,0.02,0.05,0.02;";     // dcar;  warning 0.02 cm; error 0.05 cm  (nominal ~ 0.02 cm)
    
    sTrendVars+="QA.TPC.dcar_posA_0,TPC.Anchor.dcar_posA_0,0.2,1.0,0.2;";     // dcar;  warning 0.1 cm; error 0.15 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcar_negA_0,TPC.Anchor.dcar_negA_0,0.2,1.0,0.2;";     // dcar;  warning 0.1 cm; error 0.15 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcar_posC_0,TPC.Anchor.dcar_posC_0,0.2,1.0,0.2;";     // dcar;  warning 0.1 cm; error 0.15 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcar_negC_0,TPC.Anchor.dcar_negC_0,0.2,1.0,0.2;";     // dcar;  warning 0.1 cm; error 0.15 cm  (nominal ~ 0.2 cm)    
 
    sTrendVars+="QA.TPC.dcaz_posA_0,TPC.Anchor.dcaz_posA_0,0.2,1.0,0.2;";     // dcaz;  warning 0.1 cm; error 0.15 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcaz_negA_0,TPC.Anchor.dcaz_negA_0,0.2,1.0,0.2;";     // dcaz;  warning 0.1 cm; error 0.15 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcaz_posC_0,TPC.Anchor.dcaz_posC_0,0.2,1.0,0.2;";     // dcaz;  warning 0.1 cm; error 0.15 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcaz_negC_0,TPC.Anchor.dcaz_negC_0,0.2,1.0,0.2;";     // dcaz;  warning 0.1 cm; error 0.15 cm  (nominal ~ 0.2 cm)

    // Eff ITS: TPC->ITS //TODO add comment for each cut variable  (Sebastian)
    sTrendVars+="QA.ITS.EffoneSPDPt02,ITS.Anchor.EffoneSPDPt02,0.05,0.1,0.07;"; // ITSeff warning+-5%, error+-10%; acceptable +-7%
    sTrendVars+="QA.ITS.EffoneSPDPt1,ITS.Anchor.EffoneSPDPt1,0.05,0.1,0.07;";   // ITSeff warning+-5%, error+-10%; acceptable +-7% 
    sTrendVars+="QA.ITS.EffoneSPDPt10,ITS.Anchor.EffoneSPDPt10,0.05,0.1,0.07;"; // ITSeff warning+-5%, error+-10%; acceptable +-7% 
    sTrendVars+="QA.ITS.EffTOTPt02,ITS.Anchor.EffTOTPt02,0.05,0.1,0.07;";       // ITSeff warning+-5%, error+-10%; acceptable +-7%
    sTrendVars+="QA.ITS.EffTOTPt1,ITS.Anchor.EffTOTPt1,0.05,0.1,0.07;";         // ITSeff warning+-5%, error+-10%; acceptable +-7%
    sTrendVars+="QA.ITS.EffTOTPt10,ITS.Anchor.EffTOTPt10,0.05,0.1,0.07;";       // ITSeff warning+-5%, error+-10%; acceptable +-7%
    
    // Eff TRD: TPC->TRD
    sTrendVars+="QA.TRD.TPCTRDmatchEffPosAll,TRD.Anchor.TPCTRDmatchEffPosAll,0.05,0.1,0.07;";  //
    sTrendVars+="QA.TRD.TPCTRDmatchEffNegAll,TRD.Anchor.TPCTRDmatchEffNegAll,0.05,0.1,0.07;";  //
    
    // dEdx
    sTrendVars+="QA.TPC.meanMIP,TPC.Anchor.meanMIP,1,2,1;";     // meanMIP;  warning 1; error 2; physics acceptable 1; (nominal ~ 50)
    sTrendVars+="QA.TPC.resolutionMIP,TPC.Anchor.resolutionMIP,0.02,0.05,0.02;";     // resolutionMIP;  warning 0.01; error 0.02; physics acceptable 0.01; (nominal ~ 0.06)
    sTrendVars+="QA.TPC.MIPattachSlopeA,TPC.Anchor.MIPattachSlopeA,1,1.5,1;";     // MIPattachSlopeA;  warning 1; error 1.5, physics acceptable 1   (nominal ~ 0)
    sTrendVars+="QA.TPC.MIPattachSlopeC,TPC.Anchor.MIPattachSlopeC,1,1.5,1;";     // MIPattachSlopeC;  warning 1; error 1.5, physics acceptable 1  (nominal ~ 0)
    sTrendVars+="QA.TPC.meanMIPele,TPC.Anchor.meanMIPele,3,6,3;";                 // MIPattachSlopeA;  warning 3; error 6;  physics acceptable 3  (nominal ~ 30)
    sTrendVars+="QA.TPC.resolutionMIPele,TPC.Anchor.resolutionMIPele,0.02,0.05,0.02;";     // MIPattachSlopeC;  warning xxx; error xxx  (nominal ~ 0)  
 
  }
  TStatToolkit::MakeAnchorAlias(treeMC,sTrendVars, doCheck, verbose);

  // 2.) Make aliases for Mean of MC-anchor of variables
  //     alarms on |(MC-Anchor)-<MC-Anchor)>| (diff0)    // and |(MC/Anchor)-<MC/Anchor)>| (ratio)
  treeMC->SetAlias("diff0.meanTPCncl" , "(QA.TPC.meanTPCncl-TPC.Anchor.meanTPCncl)");
  treeMC->SetAlias("diff0.meanTPCnclF" ,"(QA.TPC.meanTPCnclF-TPC.Anchor.meanTPCnclF)");
  treeMC->SetAlias("ratio.dcarAP0" ,    "(QA.TPC.dcarAP0/TPC.Anchor.dcarAP0)");
  treeMC->SetAlias("ratio.dcarAP1" ,    "(QA.TPC.dcarAP1/TPC.Anchor.dcarAP1)");
  treeMC->SetAlias("ratio.dcarCP0" ,    "(QA.TPC.dcarCP0/TPC.Anchor.dcarCP0)");
  treeMC->SetAlias("ratio.dcarCP1" ,    "(QA.TPC.dcarCP1/TPC.Anchor.dcarCP1)");
  treeMC->SetAlias("diff0.dcar_posA_0" ,    "(QA.TPC.dcar_posA_0-TPC.Anchor.dcar_posA_0)");
  treeMC->SetAlias("diff0.dcar_posC_0" ,    "(QA.TPC.dcar_posC_0-TPC.Anchor.dcar_posC_0)");
  treeMC->SetAlias("diff0.dcar_negA_0" ,    "(QA.TPC.dcar_negA_0-TPC.Anchor.dcar_negA_0)");
  treeMC->SetAlias("diff0.dcar_negC_0" ,    "(QA.TPC.dcar_negC_0-TPC.Anchor.dcar_negC_0)");  
  treeMC->SetAlias("diff0.dcaz_posA_0" ,    "(QA.TPC.dcaz_posA_0-TPC.Anchor.dcaz_posA_0)");
  treeMC->SetAlias("diff0.dcaz_posC_0" ,    "(QA.TPC.dcaz_posC_0-TPC.Anchor.dcaz_posC_0)");
  treeMC->SetAlias("diff0.dcaz_negA_0" ,    "(QA.TPC.dcaz_negA_0-TPC.Anchor.dcaz_negA_0)");
  treeMC->SetAlias("diff0.dcaz_negC_0" ,    "(QA.TPC.dcaz_negC_0-TPC.Anchor.dcaz_negC_0)");
  treeMC->SetAlias("ratio.meanMIP" ,    "(QA.TPC.meanMIP/TPC.Anchor.meanMIP)");
  treeMC->SetAlias("ratio.resolutionMIP" ,    "(QA.TPC.resolutionMIP/TPC.Anchor.resolutionMIP)");
  treeMC->SetAlias("diff0.MIPattachSlopeA" ,    "(QA.TPC.MIPattachSlopeA-TPC.Anchor.MIPattachSlopeA)");
  treeMC->SetAlias("diff0.MIPattachSlopeC" ,    "(QA.TPC.MIPattachSlopeC-TPC.Anchor.MIPattachSlopeC)");
  treeMC->SetAlias("diff0.meanMIPele" ,    "(QA.TPC.meanMIPele-TPC.Anchor.meanMIPele)"); 
  treeMC->SetAlias("diff0.resolutionMIPele" ,    "(QA.TPC.resolutionMIPele-TPC.Anchor.resolutionMIPele)");

  TString sDiffVars="";
  sDiffVars+="diff0.meanTPCncl;diff0.meanTPCnclF;";
  sDiffVars+="ratio.dcarAP0;ratio.dcarAP1;ratio.dcarCP0;ratio.dcarCP1;diff0.dcar_posA_0;diff0.dcar_posC_0;diff0.dcar_negA_0;diff0.dcar_negC_0;";
  sDiffVars+="diff0.dcaz_posA_0;diff0.dcaz_posC_0;diff0.dcaz_negA_0;diff0.dcaz_negC_0;";
  sDiffVars+="ratio.meanMIP;ratio.resolutionMIP;diff0.MIPattachSlopeA;diff0.MIPattachSlopeC;diff0.meanMIPele;diff0.resolutionMIPele;";

  TObjArray* oaTrendVars = sDiffVars.Tokenize(";");
  Float_t entryFraction=0.8, nSigmaOutlier=6., nSigmaWarning=3., epsilon=1.0e-6, rangeFactor=0.1;
  for (Int_t iVar=0; iVar<oaTrendVars->GetEntriesFast(); iVar++) {
    TString sVar( oaTrendVars->At(iVar)->GetName() );
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_OutlierMin:(MeanEF-%f*RMSEF-%f):%f", nSigmaOutlier, epsilon, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_OutlierMax:(MeanEF+%f*RMSEF+%f):%f", nSigmaOutlier, epsilon, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_WarningMin:(MeanEF-%f*RMSEF-%f):%f", nSigmaWarning, epsilon, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_WarningMax:(MeanEF+%f*RMSEF+%f):%f", nSigmaWarning, epsilon, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", rangeFactor, entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", rangeFactor, entryFraction));
    //
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_RobustMean:(MeanEF+0):%f", entryFraction));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_Outlier:(varname>varname_OutlierMax||varname<varname_OutlierMin)"));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_Warning:(varname>varname_WarningMax||varname<varname_WarningMin)"));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_PhysAcc:(varname>varname_PhysAccMin&&varname<varname_PhysAccMax)"));
  }
  // TODO Setup physics acceptable in case default+-0.1 is not acceptable (Sebastian)
  //    All variable with units
  //    ratios for which default 10 % cut is not appropriate
  //    In some case e.g efficiency we should allow higher efficiency ... - could be asymmetric
  treeMC->SetAlias("diff0.MIPattachSlopeA_PhysAcc","abs(diff0.MIPattachSlopeA-diff0.MIPattachSlopeA_RobustMean)<0.5");      // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.MIPattachSlopeC_PhysAcc","abs(diff0.MIPattachSlopeC-diff0.MIPattachSlopeC_RobustMean)<0.5");      // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.meanMIPele_PhysAcc","abs(diff0.meanMIPele-diff0.meanMIPele_RobustMean)<0.5");                     // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.resolutionMIPele_PhysAcc","abs(diff0.resolutionMIPele-diff0.resolutionMIPele_RobustMean)<0.5");   // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.dcar_posA_0_PhysAcc","abs(diff0.dcar_posA_0-diff0.dcar_posA_0_RobustMean)<0.5");                  // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.dcar_posC_0_PhysAcc","abs(diff0.dcar_posC_0-diff0.dcar_posC_0_RobustMean)<0.5");                  // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.dcar_negA_0_PhysAcc","abs(diff0.dcar_negA_0-diff0.dcar_negA_0_RobustMean)<0.5");                  // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.dcar_negC_0_PhysAcc","abs(diff0.dcar_negC_0-diff0.dcar_negC_0_RobustMean)<0.5");                  // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.dcaz_posA_0_PhysAcc","abs(diff0.dcaz_posA_0-diff0.dcaz_posA_0_RobustMean)<0.5");                  // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.dcaz_posC_0_PhysAcc","abs(diff0.dcaz_posC_0-diff0.dcaz_posC_0_RobustMean)<0.5");                  // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.dcaz_negA_0_PhysAcc","abs(diff0.dcaz_negA_0-diff0.dcaz_negA_0_RobustMean)<0.5");                  // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.dcaz_negC_0_PhysAcc","abs(diff0.dcaz_negC_0-diff0.dcaz_negC_0_RobustMean)<0.5");                  // phys. acceptable +-0.5  
  treeMC->SetAlias("diff0.meanTPCncl_PhysAcc","abs(diff0.meanTPCncl-diff0.meanTPCncl_RobustMean)<0.5");                     // phys. acceptable +-0.5
  treeMC->SetAlias("diff0.meanTPCnclF_PhysAcc","abs(diff0.meanTPCnclF-diff0.meanTPCnclF_RobustMean)<0.5");                  // phys. acceptable +-0.5  
  treeMC->SetAlias("ratio.dcarAP0_PhysAcc","abs(ratio.dcarAP0-ratio.dcarAP0_RobustMean)<0.5");                              // phys. acceptable +-0.5
  treeMC->SetAlias("ratio.dcarAP1_PhysAcc","abs(ratio.dcarAP1-ratio.dcarAP1_RobustMean)<0.5");                              // phys. acceptable +-0.5   
  treeMC->SetAlias("ratio.dcarCP0_PhysAcc","abs(ratio.dcarCP0-ratio.dcarCP0_RobustMean)<0.5");                              // phys. acceptable +-0.5
  treeMC->SetAlias("ratio.dcarCP1_PhysAcc","abs(ratio.dcarCP1-ratio.dcarCP1_RobustMean)<0.5");                              // phys. acceptable +-0.5   
  treeMC->SetAlias("ratio.meanMIP_PhysAcc","abs(ratio.meanMIP-ratio.meanMIP_RobustMean)<0.5");                        // phys. acceptable +-0.5
  treeMC->SetAlias("ratio.resolutionMIP_PhysAcc","abs(ratio.resolutionMIP-ratio.resolutionMIP_RobustMean)<0.5");            // phys. acceptable +-0.5
  
  
  // 3.) Configure combined status. Default using logical OR of problems
  TString sCombinedStatus=";";
  sCombinedStatus+="mcAnchor.ncl,absDiff.QA.TPC.meanTPCncl,absDiff.QA.TPC.meanTPCnclF,diff0.meanTPCncl,diff0.meanTPCnclF;"; // Status number of clusters and findable clusters
  sCombinedStatus+="mcAnchor.dcarResol,absDiff.QA.TPC.dcarAP0,absDiff.QA.TPC.dcarAP1,absDiff.QA.TPC.dcarCP0,absDiff.QA.TPC.dcarCP1,ratio.dcarAP0,ratio.dcarAP1,ratio.dcarCP0,ratio.dcarCP1,diff0.dcar_posA_0,diff0.dcar_posC_0,diff0.dcar_negA_0,diff0.dcar_negC_0;";  // Status: DCAr resolution
  sCombinedStatus+="mcAnchor.dcazResol,absDiff.QA.TPC.dcaz_posA_0,absDiff.QA.TPC.dcaz_posC_0,absDiff.QA.TPC.dcaz_negA_0,absDiff.QA.TPC.dcaz_negC_0,diff0.dcaz_posA_0,diff0.dcaz_posC_0,diff0.dcaz_negA_0,diff0.dcaz_negC_0;";  // Status: DCAz resolution
  
  sCombinedStatus+="mcAnchor.itsEffStatus,absDiff.QA.ITS.EffoneSPDPt02,absDiff.QA.ITS.EffoneSPDPt1,absDiff.QA.ITS.EffoneSPDPt10,absDiff.QA.ITS.EffTOTPt02,absDiff.QA.ITS.EffTOTPt1,";
  sCombinedStatus+="absDiff.QA.TRD.TPCTRDmatchEffPosAll,absDiff.QA.TRD.TPCTRDmatchEffNegAll;";
  
  sCombinedStatus+="mcAnchor.dEdx,ratio.meanMIP,ratio.resolutionMIP,diff0.MIPattachSlopeA,diff0.MIPattachSlopeC,diff0.meanMIPele,diff0.resolutionMIPele,absDiff.QA.TPC.meanMIP,absDiff.QA.TPC.resolutionMIP,absDiff.QA.TPC.MIPattachSlopeA,absDiff.QA.TPC.MIPattachSlopeC,absDiff.QA.TPC.meanMIPele,absDiff.QA.TPC.resolutionMIPele;";

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
  if(externalInfo==0) externalInfo = new AliExternalInfo(".", "", verbose);
  if(trendingDraw==0) trendingDraw = new AliTreeTrending("mcAnchor","mcAnchor");
  trendingDraw->SetDefaultStyle();
  
  gStyle->SetOptTitle(0);

  if(treeMC==0) treeMC = externalInfo->GetTree("QA.TPC", mcPeriod, mcPass, "QA.TPC;QA.TRD;QA.TOF;QA.ITS");

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
  treeMC->SetAlias("statisticOK","QA.TPC.nEvents>100&&TPC.Anchor.nEvents>0");
  treeMC->SetAlias("present","(run==TPC.Anchor.run)");      // check presence of reference detectors for MC tree
  treeAnchorTPC->SetAlias("present","(run>0)");                                // define present alias for anchor tree
  TStatToolkit::AddMetadata(treeMC,"QA.TPC.Legend",("MC_{TPC}: "+mcPeriod).Data());
  TStatToolkit::AddMetadata(treeMC,"Logbook.Legend","Logbook");
  TStatToolkit::AddMetadata(treeMC,"QA.ITS.Legend",("MC_{ITS}: "+mcPeriod).Data());
  TStatToolkit::AddMetadata(treeMC,"QA.TRD.Legend",("MC_{TRD}: "+mcPeriod).Data());
  TStatToolkit::AddMetadata(treeMC,"TPC.Anchor.Legend",(anchorPeriod+"_{"+anchorPass+"}").Data());
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
  TString sStatusBarVars("mcAnchor.ncl;mcAnchor.dcarResol;mcAnchor.dcazResol;mcAnchor.itsEffStatus;mcAnchor.dEdx");
  TString sStatusBarNames("#(cl)_{MC Anchor};dca_{R}_{MC Anchor};dca_{z}_{MC Anchor};#epsilon_{ITS}_{MC Anchor};dEdx_{MC Anchor}");
  TString sCriteria("(present):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)"); // status bar markers
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
  // treeMC->SetAlias("defaultCut", "run>0");   // TODO enable all runs - need to make reorting of the multigraphs
  if (trendingDraw->InitSummaryTrending(statusString, 0.015, "defaultCut")) {
    trendingDraw->fWorkingCanvas->Print(TString(outputDir) + "/report.pdf[", "pdf");
    return kTRUE;
  }
  else return kFALSE;
}

/// ## MakeReport
/// \param outputDir
void MakeReport() {

  //
  // TODO: Partially done. Add TMultiGraph "class?" to specify additional options (see bellow) (Boris, Marian) - groupName to be used for that
  // OLDTODO: Optionally draw y value on top of the markers - not needed with jsroot
  // TODO - add time format to y axis - to do with class
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
                         "defaultCut", "figTemplateTRD", "figTemplateTRD", 1, 1.0, 4, kTRUE);
  //trendingDraw->MakePlot(outputDir,"runTime.png","Time ",cRange,"",":run","defaultCut","figTemplateTRDPair","figTemplateTRDPair",1,1.0,4,kTRUE);
  trendingDraw->MakePlot(outputDir, "meanMult.png", "Mean TPC multiplicity (|DCA|<3 cm)", cRange, "",
                         "QA.TPC.meanMult;TPC.Anchor.meanMult:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanMultPos.png", "Mean TPC multiplicity (q>0, |DCA|<3 cm)", cRange, "",
                         "QA.TPC.meanMultPos;TPC.Anchor.meanMultPos:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanMultNeg.png", "Mean TPC multiplicity (q<0, |DCA|<3 cm)", cRange, "",
                         "QA.TPC.meanMultNeg;TPC.Anchor.meanMultNeg:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanMult_comb2.png", "meanMult_comb2", cRange, "",
                         "meanMult_comb2;TPC.Anchor.meanMult_comb2:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanVertX.png", "meanVertX:run MC/Anchor", cRange, "",
                         "QA.TPC.meanVertX;TPC.Anchor.meanVertX:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 5, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanVertY.png", "meanVertY:run MC/Anchor", cRange, "",
                         "QA.TPC.meanVertY;TPC.Anchor.meanVertY:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 5, kTRUE);
  trendingDraw->MakePlot(outputDir, "meanVertZ.png", "meanVertZ:run MC/Anchor", cRange, "",
                         "QA.TPC.meanVertZ;TPC.Anchor.meanVertZ:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 5, kTRUE);
  //
  // 2.) Number of clusters comparison ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabNcl.html)
  // TODO: Add all estimators fo missing chambers (Ncl, Voltage, RawQA, tracks)
  trendingDraw->MakePlot(outputDir, "meanTPCncl.png", "Number of clusters", cRange, "",
                         "QA.TPC.meanTPCncl;TPC.Anchor.meanTPCncl:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 5, kTRUE);
  trendingDraw->AppendDefaultBands(outputDir,"meanTPCncl.png","QA.TPC.meanTPCncl", "absDiff.QA.TPC.meanTPCncl", "defaultCut","");
  
  trendingDraw->MakePlot(outputDir, "meanTPCnclFindable.png", "Cluster fraction #left(#frac{N_{cl}}{N_{find.}}#right)",
                         cRange, "", "QA.TPC.meanTPCnclF;TPC.Anchor.meanTPCnclF:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 5, kTRUE);
  trendingDraw->AppendDefaultBands(outputDir,"meanTPCnclFindable.png","QA.TPC.meanTPCnclF", "absDiff.QA.TPC.meanTPCnclF", "defaultCut","");
    
  trendingDraw->MakePlot(outputDir, "meanTPCNclRatioMCtoAnchor.png", "Number of clusters MC/Anchor", cRange, "",
                         "meanTPCncl/TPC.Anchor.meanTPCncl;meanTPCnclF/TPC.Anchor.meanTPCnclF:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 5, kTRUE);
  trendingDraw->MakePlot(outputDir, "iroc.png", "IROC #chambers", cRange, "",
                         "iroc_A_side;TPC.Anchor.iroc_A_side;iroc_C_side;TPC.Anchor.iroc_C_side:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "oroc.png", "OROC #chambers", cRange, "",
                         "oroc_A_side;TPC.Anchor.oroc_A_side;oroc_C_side;TPC.Anchor.oroc_C_side:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  //
  // 3.) Matching efficiency ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabEff.html)
  // TODO: Add ITS layer matching trending
  // TODO: Add TRD matching
  trendingDraw->MakePlot(outputDir, "matchingTPCITSEffNoPileUpCut.png", "Matching efficiency(no pileup cut):MC/Anchor",
                         cRange, "",
                         "QA.TPC.tpcItsMatchA;TPC.Anchor.tpcItsMatchA;QA.TPC.tpcItsMatchC;TPC.Anchor.tpcItsMatchC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPCITSEffPileUpCut.png", "Matching efficiency (pileup cut):MC/Anchor",
                         cRange, "",
                         "QA.ITS.EffTOTPt02;ITS.Anchor.EffTOTPt02;QA.ITS.EffTOTPt1;ITS.Anchor.EffTOTPt1:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPCITSEffPileUpCutHighPt.png",
                         "Matching efficiency (pileup cut):MC/Anchor", cRange, "",
                         "QA.ITS.EffTOTPt1;ITS.Anchor.EffTOTPt1;QA.ITS.EffTOTPt10;ITS.Anchor.EffTOTPt10:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPCITSEffACRatio.png", "Matching efficiency A/C ratio:MC/Anchor", cRange,
                         "",
                         "QA.TPC.tpcItsMatchA/QA.TPC.tpcItsMatchC;TPC.Anchor.tpcItsMatchA/TPC.Anchor.tpcItsMatchC;QA.TPC.tpcItsMatchHighPtA/QA.TPC.tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtA/TPC.Anchor.tpcItsMatchHighPtC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPCTRDEffPileUpCut.png", "Matching efficiency (pileup cut):MC/Anchor",
                         cRange, "",
                         "QA.TRD.TPCTRDmatchEffPosAll;TRD.Anchor.TPCTRDmatchEffPosAll;QA.TRD.TPCTRDmatchEffNegAll;TRD.Anchor.TPCTRDmatchEffNegAll:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  // TODO plots bellow not in the html (to be described or removed)
  trendingDraw->MakePlot(outputDir, "matchingTPC-ITSEff.png", "Matching efficiency:MC/Anchor", cRange, "",
                         "QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffTOTPt02;QA.ITS.EffTOTPt1;QA.ITS.EffTOTPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffTOTPt02;ITS.Anchor.EffTOTPt1;ITS.Anchor.EffTOTPt10:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcItsMatch.png", "TPC - ITS Match", cRange, "",
                         "tpcItsMatchHighPtA;TPC.Anchor.tpcItsMatchHighPtA;tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPC-ITSEff_1.png", "Matching efficiency x", cRange, "",
                         "QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffoneSPDPt02;QA.ITS.EffoneSPDPt1;QA.ITS.EffoneSPDPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffoneSPDPt02;ITS.Anchor.EffoneSPDPt1;ITS.Anchor.EffoneSPDPt10:run",
                         "defaultCut", "21;24;25;27;28;21;24;25;27;28", "2;2;2;figTemplateTRDPair;4;4;4", 1, 1.0, 6,
                         kTRUE);
  trendingDraw->MakePlot(outputDir, "matchingTPC-ITSEff_2.png", "Matching efficiency y", cRange, "",
                         "QA.TPC.tpcItsMatchA/TPC.Anchor.tpcItsMatchA;QA.TPC.tpcItsMatchC/TPC.Anchor.tpcItsMatchC;QA.ITS.EffTOTPt02/ITS.Anchor.EffTOTPt02;QA.ITS.EffTOTPt1/ITS.Anchor.EffTOTPt1;QA.ITS.EffTOTPt10/ITS.Anchor.EffTOTPt10:run",
                         "defaultCut", "21;24;25;27;28;21;24;25;27;28", "1;2;4;3;6;2;4;4;4;4;4", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcItsMatch.png", "tpcItsMatch", cRange, "",
                         "tpcItsMatchA;TPC.Anchor.tpcItsMatchA;tpcItsMatchC;TPC.Anchor.tpcItsMatchC:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcItsMatchHighPt.png", "tpcItsMatchHighPt", cRange, "",
                         "tpcItsMatchHighPtA;TPC.Anchor.tpcItsMatchHighPtA;tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  //
  // 4.) DCA  ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabDCA.html)
  //
  trendingDraw->MakePlot(outputDir, "dcarP0.png",
                         "HighPt: DCA_{xy} #sigma_{0} (#sigma^{2}=#sigma_{0}^{2}+#sigma_{1}^{2}/p_{T}^{2}) ", cRange,
                         "", "QA.TPC.dcarAP0;TPC.Anchor.dcarAP0;QA.TPC.dcarCP0;TPC.Anchor.dcarCP0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcarP1.png",
                         "MS: DCA_{xy} #sigma_{1} (#sigma^{2}=#sigma_{0}^{2}+#sigma_{1}^{2}/p_{T}^{2})", cRange, "",
                         "QA.TPC.dcarAP1;TPC.Anchor.dcarAP1;QA.TPC.dcarCP1;TPC.Anchor.dcarCP1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(0, "offsetR.png", "Offset dR", cRange, "",
                         "offsetdRA;TPC.Anchor.offsetdRA;offsetdRC;TPC.Anchor.offsetdRC:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->AppendBand(outputDir, "offsetdR.png",
                           "offsetdRA_RobustMean;offsetdRA_OutlierMin;offsetdRA_OutlierMax;offsetdRA_WarningMin;offsetdRA_WarningMax:run",
                           "defaultCut", "1;1;1;1;1,1;2;2;3;3", "1;1;1;1;1,1;figTemplateTRDPair", 1, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_A_0.png", "dcar_A_0", cRange, "",
                         "dcar_posA_0;TPC.Anchor.dcar_posA_0;dcar_negA_0;TPC.Anchor.dcar_negA_0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_C_0.png", "dcar_C_0", cRange, "",
                         "dcar_posC_0;TPC.Anchor.dcar_posC_0;dcar_negC_0;TPC.Anchor.dcar_negC_0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_C_1.png", "dcar_C_1", cRange, "",
                         "dcar_posC_1;TPC.Anchor.dcar_posC_1;dcar_negC_1;TPC.Anchor.dcar_negC_1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_A_1.png", "dcar_A_1", cRange, "",
                         "dcar_posA_1;TPC.Anchor.dcar_posA_1;dcar_negA_1;TPC.Anchor.dcar_negA_1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_A_2.png", "dcar_A_2", cRange, "",
                         "dcar_posA_2;TPC.Anchor.dcar_posA_2;dcar_negA_2;TPC.Anchor.dcar_negA_2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcar_C_2.png", "dcar_C_2", cRange, "",
                         "dcar_posC_2;TPC.Anchor.dcar_posC_2;dcar_negC_2;TPC.Anchor.dcar_negC_2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_A_0.png", "dcaz_A_0", cRange, "",
                         "dcaz_posA_0;TPC.Anchor.dcaz_posA_0;dcaz_negA_0;TPC.Anchor.dcaz_negA_0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_C_0.png", "dcaz_C_0", cRange, "",
                         "dcaz_posC_0;TPC.Anchor.dcaz_posC_0;dcaz_negC_0;TPC.Anchor.dcaz_negC_0:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_C_1.png", "dcaz_C_1", cRange, "",
                         "dcaz_posC_1;TPC.Anchor.dcaz_posC_1;dcaz_negC_1;TPC.Anchor.dcaz_negC_1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_A_1.png", "dcaz_A_1", cRange, "",
                         "dcaz_posA_1;TPC.Anchor.dcaz_posA_1;dcaz_negA_1;TPC.Anchor.dcaz_negA_1:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_A_2.png", "dcaz_A_2", cRange, "",
                         "dcaz_posA_2;TPC.Anchor.dcaz_posA_2;dcaz_negA_2;TPC.Anchor.dcaz_negA_2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcaz_C_2.png", "dcaz_C_2", cRange, "",
                         "dcaz_posC_2;TPC.Anchor.dcaz_posC_2;dcaz_negC_2;TPC.Anchor.dcaz_negC_2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "dcarFitPar_comb4.png", "dcarFitpar_comb4", cRange, "",
                         "dcarFitpar_comb4;TPC.Anchor.dcarFitpar_comb4:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  //trendingDraw->MakePlot(outputDir,"rmsDCAMultParMCtoAnchor.png","DCA Resolution mult due MS:MC/Anchor",cRange,"","dcarAP1;dcarCP1;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP1:run","defaultCut","figTemplateTRDPair","figTemplateTRDPair",1,1.0,6,kTRUE);

  // 5.) dEdx ($AliPhysic_SRC/PWGPP/TPC/macros/TPCQAWebpage/MCAnchor/tabdEdx.html)
  // TODO: Add sub-status lines
  trendingDraw->MakePlot(outputDir, "meanMIP.png", "<Mean dEdx_{MIP}> (a.u) ", cRange, "",
                         "QA.TPC.meanMIP;TPC.Anchor.meanMIP:run:fitMIP.fElements[4];TPC.Anchor.fitMIP.fElements[4]", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->AppendDefaultBands(outputDir,"meanMIP.png","QA.TPC.meanMIP", "absDiff.QA.TPC.meanMIP", "defaultCut","");
  
  trendingDraw->MakePlot(outputDir, "meanElectron.png", "<dEdx_{el}> (a.u) ", cRange, "",
                         "QA.TPC.meanMIPele;TPC.Anchor.meanMIPele:run:fitElectron.fElements[4];TPC.Anchor.fitElectron.fElements[4]", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->AppendDefaultBands(outputDir,"meanElectron.png","QA.TPC.meanMIPele", "absDiff.QA.TPC.meanMIPele", "defaultCut","");
    
  trendingDraw->MakePlot(outputDir, "electronMIPSeparation.png", "<dEdx_{el}>-<dEdx_{MIP}>", cRange, "",
                         "QA.TPC.electroMIPSeparation;TPC.Anchor.electroMIPSeparation:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "electronMIPRatio.png", "<dEdx_{el}>/<dEdx_{MIP}>", cRange, "",
                         "QA.TPC.meanMIPele/QA.TPC.meanMIP;TPC.Anchor.meanMIPele/TPC.Anchor.meanMIP:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "MIPattachSlope.png", "MIPattachSlope", cRange, "",
                         "MIPattachSlopeA;TPC.Anchor.MIPattachSlopeA;MIPattachSlopeC;TPC.Anchor.MIPattachSlopeC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "MIPattachSlopeA.png", "MIPattachSlopeA", cRange, "",
                         "MIPattachSlopeA;TPC.Anchor.MIPattachSlopeA:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "MIPattachSlopeC.png", "MIPattachSlopeC", cRange, "",
                         "MIPattachSlopeC;TPC.Anchor.MIPattachSlopeC:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "resolutionMIP.png", "resolutionMIP", cRange, "",
                         "QA.TPC.resolutionMIP;TPC.Anchor.resolutionMIP:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
  trendingDraw->MakePlot(outputDir, "resolutionMIPele.png", "resolutionMIPele", cRange, "",
                         "QA.TPC.resolutionMIPele;TPC.Anchor.resolutionMIPele:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 4, kTRUE);
//  trendingDraw->MakePlot(outputDir, "mipToEleSeparation.png", "meanMIPele/meanMIP", cRange, "",
//                         "meanMIPele/meanMIP;TPC.Anchor.meanMIPele/TPC.Anchor.meanMIP:run", "defaultCut", "figTemplateTRDPair",
//                         "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->AppendBand(outputDir,"mipToEleSeparation.png","meanMIPele_RobustMean/meanMIP_RobustMean;TPC.Anchor.meanMIPele_RobustMean/TPC.Anchor.meanMIP_RobustMean:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", kTRUE, 1.0, kTRUE);


  //
  // 6.) matching (tabMatching.html))
  //
  trendingDraw->MakePlot(outputDir, "offsetdZ.png", "offsetdZ", cRange, "",
                         "offsetdZA;TPC.Anchor.offsetdZA;offsetdZC;TPC.Anchor.offsetdZC:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "offsetd_comb4.png", "offsetd_comb4", cRange, "",
                         "offsetd_comb4;TPC.Anchor.offsetd_comb4:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "lambdaPull.png", "lambdaPull", cRange, "", "lambdaPull;TPC.Anchor.lambdaPull:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcConstrainPhiC.png", "tpcConstrainPhiC", cRange, "",
                         "tpcConstrainPhiC;TPC.Anchor.tpcConstrainPhiC:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "yPull.png", "yPull", cRange, "", "yPull;TPC.Anchor.yPull:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "zPull.png", "zPull", cRange, "", "zPull;TPC.Anchor.zPull:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "itsTpcPulls_comb4.png", "itsTpcPulls_comb4", cRange, "",
                         "itsTpcPulls_comb4;TPC.Anchor.itsTpcPulls_comb4:run", "defaultCut", "figTemplateTRDPair",
                         "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcConstrainPhi.png", "tpcConstrainPhi", cRange, "",
                         "tpcConstrainPhiA;TPC.Anchor.tpcConstrainPhiA;tpcConstrainPhiC;TPC.Anchor.tpcConstrainPhiC:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "tpcConstrainPhi_comb2.png", "tpcConstrainPhi_comb2", cRange, "",
                         "tpcConstrainPhi_comb2;TPC.Anchor.tpcConstrainPhi_comb2:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "deltaPt.png", "deltaPt", cRange, "", "deltaPt;TPC.Anchor.deltaPt:run",
                         "defaultCut", "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  trendingDraw->MakePlot(outputDir, "deltaPtAC.png", "deltaPtAC", cRange, "",
                         "deltaPtA;TPC.Anchor.deltaPtA;deltaPtC;TPC.Anchor.deltaPtC:run", "defaultCut",
                         "figTemplateTRDPair", "figTemplateTRDPair", 1, 1.0, 6, kTRUE);
  //
  trendingDraw->fWorkingCanvas->Clear();
  trendingDraw->fWorkingCanvas->Print(TString(outputDir) + "/report.pdf]", "pdf");
}


/// TODO - alias for the html table names, hints
/// TODO - fix formatting in AliTreeFormulaF- in case of missing entry write undefined
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
  TStatToolkit::AddMetadata(treeMC, "detectorMask.headerTooltip","Detector bit mask:\n . SPD\n . SDD\n . SSS\n . TPC\n . TRD");
  //MakeStatusBitMasks();
  makeHtmlDCA();
}

void makeHtmlDCA(){

  //TString metaDCA="<a title
  TString pathMC="<a href=\"http://aliqatpc.web.cern.ch/aliqatpc/sim/%d{year}/%s{period.GetName()}/passMC/000%d{run}/";
  TString pathData="<a href=\"http://aliqatpc.web.cern.ch/aliqatpc/data/%d{TPC.Anchor.year}/%s{TPC.Anchor.period.GetName()}/%{TPC.Anchor.pass.GetName()}/000%d{run}/";
  TStatToolkit::AddMetadata(treeMC, "dcarAP0.html", (pathMC+"dca_and_phi.png\">%2.4f{dcarAP0}</a>").Data());
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarAP0.html", (pathData+"dca_and_phi.png\">%2.4f{TPC.Anchor.dcarAP0}</a>").Data());
  TStatToolkit::AddMetadata(treeMC, "dcarCP0.html", (pathMC+"dca_and_phi.png\">%2.4f{dcarCP0}</a>").Data());
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarCP0.html", (pathData+"dca_and_phi.png\">%2.4f{TPC.Anchor.dcarCP0}</a>").Data());

  TStatToolkit::AddMetadata(treeMC, "dcarAP0.thead", "&sigma;<sup>A</sup><sub>MC DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeMC, "dcarAP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeMC, "dcarCP0.thead", "&sigma;<sup>C</sup><sub>MC DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeMC, "dcarCP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarAP0.thead", "&sigma;<sup>A</sup><sub>Anchor DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarAP0.tooltip", "DCA (r&phi;) at infinite pt");
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarCP0.thead", "&sigma;<sup>C</sup><sub>Anchor DCA<sub>1/pt=0</sub></sub> (cm)");
  TStatToolkit::AddMetadata(treeMC, "TPC.Anchor.dcarCP0.tooltip", "DCA (r&phi;) at infinite pt");
   //
  TStatToolkit::AddMetadata(treeMC, "dcar_MaskWarning.thead", "DCA<sub>MCWarning</sub></sub>");
  TStatToolkit::AddMetadata(treeMC, "dcar_MaskWarning.tooltip", "DCAr status");
  TStatToolkit::AddMetadata(treeMC, "dcar_MaskWarning.html","<a href=\"dcarStatusMC.png\">%x{dcar_MaskWarning}</a>");
  //
  TString runSelection = "defaultCut";
  TString varSelection = "";
  TString logbookBase = "";
  logbookBase = "run:LHCFillNumber:LHCperiod:detectorMask:HLTmode:DAQ_time_start:ctpDuration:totalEventsPhysics:interactionRate:bz;2.2:";
  // tabDCA
  TString tpcDCA = "dcarAP0:TPC.Anchor.dcarAP0:dcarCP0:TPC.Anchor.dcarCP0:dcar_MaskWarning";//:dcarStatusString:dcarRawStatusString";
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
  AliTreeTrending::DecomposeStatusAlias(treeMC, statusString,statusVar,statusTitle,suffix,counter,statusMask);
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
  AliTreeTrending::DecomposeStatusAlias(treeMC, statusString,statusVar,statusTitle,suffix,counter,statusMask);
  statusString.ReplaceAll("_Warning","");
  TString statusAliasName, statusAliasMask;
  //
  statusAliasName=statusString+"_MaskWarning";
  statusAliasMask=statusMask;
  statusAliasMask.ReplaceAll("_#","_Warning");
  treeMC->SetAlias(statusAliasName.Data(),statusAliasMask.Data());

}


void MakeStatusPlots(){
  try {
    MakeStatusPlot("./", "dcarStatusMC.png", "dcar_Warning", "1");
    MakeStatusPlot("./", "dcarStatusAnchor.png", "dcar_Warning", "1", "TPC.Anchor");
    MakeStatusPlot("./", "dcarStatusMCToAnchor.png", "mcAnchor.dcarResol_Warning", "1");

    MakeStatusPlot("./", "dcazStatusMC.png", "dcaz_Warning", "1");
    MakeStatusPlot("./", "dcazStatusAnchor.png", "dcaz_Warning", "1", "TPC.Anchor");
    MakeStatusPlot("./", "dcazStatusMCToAnchor.png", "mcAnchor.dcazResol_Warning", "1");
    
    MakeStatusPlot("./", "itsEffStatusMCToAnchor.png", "mcAnchor.itsEffStatus_Warning", "1");
//    MakeStatusPlot("./", "itsEffStatusMC.png", "itsEffStatus_Warning", "1");                  //missing corresponding alias in MC tree
//    MakeStatusPlot("./", "itsEffStatusAnchor.png", "itsEffStatus_Warning", "1","TRD.Anchor"); //missing corresponding alias in Anchor tree
    //
    MakeStatusPlot("./", "dEdxStatusMC.png", "MIPquality_Warning", "1");
    MakeStatusPlot("./", "dEdxStatusAnchor.png", "MIPquality_Warning", "1","TPC.Anchor");
    MakeStatusPlot("./", "dEdxStatusMCToAnchor.png", "mcAnchor.dEdx_Warning", "1");

    MakeStatusPlot("./", "nclStatusMCToAnchor.png", "mcAnchor.ncl_Warning", "1");
//    MakeStatusPlot("./", "nclStatusAnchor.png", "ncl_Warning", "1","TPC.Anchor");    //missing corresponding alias
//    MakeStatusPlot("./", "nclStatusMC.png", "ncl_Warning", "1");                     //missing corresponding alias    
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




/// DONE  - automatic decomposition status aliases
/// DONE  - for the status graphs with friend - we need to use run list for other status
///       - RebinMultiGraph should keep offsets
/// DONE  (TO COMMIT) - assign proper names to generated graphs
/// DONE  (TO COMMIT) - use anchor names in the legend
/// DONE - FIX  Draw MakeMultiGraph drawing graph 2 times - Fix in TStatToolkit///
/// TODO - (Ongoing) problem with the markers and y axis labels in JSROOT - trying to avoid it increasing marker size from0.75->1. - did not help
///      - way around avoid marker style 20 - dow not work for status bar
///      - jsroot responsible contacted https://alice.its.cern.ch/jira/browse/ATO-393


/// TODO - get form logbook detailed per run information  (put into html metadata string)
/// TODO - Ful name of the alias expression (maybe recursive)
/// TODO - add html preview for table header
///          tableHeader, tableHeader_Tooltip, tableHeader_Title, tableHeader_html

/// TODO - Graph rebin to be extended to X,Y,Z
/// TODO - tooltip/title- to make bullets using UTF printing
///
/// TODO  - make dEdx,dcaz, TPC-TRD eff (together with its) full QA (Sebastian)
///      - all status(MC/Anchor, MC, Anchor) plots
///      - add them to web page = modifying tabdEdx.html, tabdDCAR.html, tabDCAZ.html tabEff.html
///      - add plot for all variables contributing to alarms
///      - for TRD eff see the variables in treeMC->GetFriend("QA.TRD")->GetListOfBranches()->Print("","*ff*")
///      - plot is there can be added to the alarms
///      
///      - for plots with only one variable plot band 
///
/// TOD0  - dcar, dcaz - including fits (currently only resolution)
///       - TPC-ITS matching
///
/// TODO - add calibration trending plots and alarms (Marian)
///      - TStatToolkit::AddMetadata(treeMC,"status2.Legend","TPC ON")
///      - TStatToolkit::GetMetadata(treeMC,"status2.Legend")->GetTitle()


/// TODO - AddAppendBandDefault - adding warning outlier and phys acceptable around reference data;
