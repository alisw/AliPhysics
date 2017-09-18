/*
  .L $NOTES/JIRA/ATO-83/code/tpcMCValidationStandardQA.C+
  TString mcPeriod="LHC15k1a1";
  TString mcPass="passMC";
  TString anchorPeriod="LHC15o";
  TString anchorPass="pass3_lowIR_pidfix";
  
  InitTPCMCValidation("LHC15k1a1","passMC","LHC15o", "pass3_lowIR_pidfix",0,0);
 

*/ 
//gSystem->AddIncludePath("-I$ALICE_ROOT/include/"); //couldn't add include path in .rootrc

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

AliExternalInfo   *pinfo=0;
AliTreeTrending   *trendingDraw=0;  
TTree * treeMC;
  
void makeTPCMCAlarms(TTree * treeMC, Bool_t doCheck,Int_t verbose);
void tpcMCValidationStandard(TString mcPeriod,Int_t verbose,Int_t doCheck);
Bool_t InitTPCMCValidation(TString mcPeriod,  TString mcPass, TString anchorPeriod,  TString anchorPass, Int_t verbose,Int_t doCheck);
void MakeReport(const char* mcrddir);
void MakePlot(TTree * tree,const char* mcrddir, const char *fname, const char *LegendTitle, std::vector<Double_t>& legpos, const char *groupName, const char* expr, const char * cut, const char * markers, const char *colors, Bool_t drawSparse, Float_t msize, Float_t sigmaRange, Bool_t comp);

void tpcMCValidation(const char* MCper ="LHC15k1a1",const char* mcrddir="./"){

//gROOT->LoadMacro("tpcMCValidationStandardQA.C+");
cout<<"INITIALIZING TPC MC Validation"<<endl;
AliExternalInfo i;
cout<<MCper<<endl;
TString AnchProdNamenPass = i.GetMCPassGuess(TString::Format("%s",MCper));
//TString AnchProdName;
cout<<"Anchor Production Name and Pass: "<<AnchProdNamenPass<<endl;
TObjArray *subStrL;
subStrL = TPRegexp("^([^ ]+)").MatchS(AnchProdNamenPass);
TString AnchProdName = ((TObjString*)subStrL->At(0))->GetString();

subStrL = TPRegexp("([^ ])+$").MatchS(AnchProdNamenPass);
TString AnchPassName = ((TObjString*)subStrL->At(0))->GetString();


if(InitTPCMCValidation(MCper,"passMC",AnchProdName, AnchPassName,0,0)){

    MakeReport(mcrddir);
}
else ::Error("tpcMCValidation","InitTPCMCValidation returned with error -> skip plotting!");
}



void tpcMCValidationStandard(TString mcPeriod,  TString mcPass, TString anchorPeriod,  TString anchorPass, Int_t verbose,Int_t doCheck){
  InitTPCMCValidation( mcPeriod,mcPass, anchorPeriod, anchorPass,verbose,doCheck);
}



void makeTPCMCAlarms(TTree * treeMC, Bool_t doCheck,Int_t verbose){
  //
  //  define variables for alarms for MC/raw mismatch
  //
  //  ==============================================================
  //  1.)  Partial alarms  (variable, variableAnchor deltaWarning,deltaError, PhysAcc) 
  //                 deltaWarning and deltaErrror can be an expression which is understtod by TTreeFormula
  //  ==============================================================
  TString sTrendVars=";";   
  {
    // Ncl
    sTrendVars+="QA.TPC.meanTPCncl,TPC.Anchor.meanTPCncl,5,10,5;";       // delta Ncl  warning 5 ,  error 10     (nominal ~ 100-140)
    sTrendVars+="QA.TPC.meanTPCnclF,TPC.Anchor.meanTPCnclF,0.02,0.05,0.05;"; // delta NclF  warning 2%,  error 5%    (nominal ~ 90%)
    // dcaR resolution
    sTrendVars+="QA.TPC.dcarAP0,TPC.Anchor.dcarAP0,0.02,0.05,0.02;";     // dcarAP0;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcarCP0,TPC.Anchor.dcarCP0,0.02,0.05,0.02;";     // dcarCP0;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcarAP1,TPC.Anchor.dcarAP1,0.02,0.05,0.02;";     // dcarAP1;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
    sTrendVars+="QA.TPC.dcarCP1,TPC.Anchor.dcarCP1,0.02,0.05,0.02;";     // dcarCP1;  warning 0.02cm; error 0.05 cm  (nominal ~ 0.2 cm)
    // Eff ITS: TPC->ITS
    sTrendVars+="QA.ITS.EffoneSPDPt02,ITS.Anchor.EffoneSPDPt02,0.05,0.1,0.07;";
    sTrendVars+="QA.ITS.EffoneSPDPt1,ITS.Anchor.EffoneSPDPt1,0.05,0.1,0.07;";
    sTrendVars+="QA.ITS.EffoneSPDPt10,ITS.Anchor.EffoneSPDPt10,0.05,0.1,0.07;";        
    sTrendVars+="QA.ITS.EffTOTPt02,ITS.Anchor.EffTOTPt02,0.05,0.1,0.07;";
    sTrendVars+="QA.ITS.EffTOTPt1,ITS.Anchor.EffTOTPt1,0.05,0.1,0.07;";
    sTrendVars+="QA.ITS.EffTOTPt10,ITS.Anchor.EffTOTPt10,0.05,0.1,0.07;";    
    // Eff TRD: TPC->TRD    
  }
  TStatToolkit::MakeAnchorAlias(treeMC,sTrendVars, doCheck, verbose);
  
  
  //Make aliases for Mean of MC-anchor of variables
  treeMC->SetAlias("mcrddiff_meanTPCncl" , "QA.TPC.meanTPCncl-TPC.Anchor.meanTPCncl");
  
  
  
  treeMC->SetAlias("statisticOK", "(meanTPCncl>0)");
  
  TString sDiffVars=";"; 
  sDiffVars+=";mcrddiff_meanTPCncl;";
          
  TObjArray* oaTrendVars = sDiffVars.Tokenize(";");
  Float_t entryFrac=0.8, nsigmaOutlier=6., nsigmaWarning=3., epsilon=1.0e-6, combfac=1.;
  
  for (Int_t vari=0; vari<oaTrendVars->GetEntriesFast(); vari++)
  {
    TString sVar( oaTrendVars->At(vari)->GetName() );

    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_OutlierMin:(MeanEF-%f*RMSEF-%f):%f", nsigmaOutlier, epsilon, entryFrac));   
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_OutlierMax:(MeanEF+%f*RMSEF+%f):%f", nsigmaOutlier, epsilon, entryFrac));   
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_WarningMin:(MeanEF-%f*RMSEF-%f):%f", nsigmaWarning, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_WarningMax:(MeanEF+%f*RMSEF+%f):%f", nsigmaWarning, epsilon, entryFrac));    
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", 0.05*combfac, entryFrac));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", 0.05*combfac, entryFrac));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "statisticOK", Form("varname_RobustMean:(MeanEF+0):%f", entryFrac));
  }
  
 for (Int_t vari=0; vari<oaTrendVars->GetEntriesFast(); vari++)
  {
    TString sVar( oaTrendVars->At(vari)->GetName() );
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "", Form("varname_Outlier:(varname>varname_OutlierMax||varname<varname_OutlierMin)"));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "", Form("varname_Warning:(varname>varname_WarningMax||varname<varname_WarningMin)"));
    TStatToolkit::SetStatusAlias(treeMC, sVar.Data(),    "", Form("varname_PhysAcc:(varname>varname_PhysAccMin&&varname<varname_PhysAccMax)"));
 }
  
  //  ==============================================================
  // 2.) Configure combined status. Default using logiacal OR of problems
  //  ==============================================================
  TString sCombinedStatus=";";
  sCombinedStatus+="ncl,TPC.Anchor.meanTPCncl,TPC.Anchor.meanTPCnclF;";               // Status number of clusters and findable clusters
  sCombinedStatus+="dcarResol,TPC.Anchor.dcarAP0,TPC.Anchor.dcarAP1,TPC.Anchor.dcarCP0,TPC.Anchor.dcarCP1;";  // Status: DCA resolution
  sCombinedStatus+="itsEffStatus,ITS.Anchor.EffoneSPDPt02,ITS.Anchor.EffoneSPDPt1,ITS.Anchor.EffoneSPDPt10,ITS.Anchor.EffTOTPt02,ITS.Anchor.EffTOTPt1,ITS.Anchor.EffTOTPt10;";
  sCombinedStatus += "mcrddiff,mcrddiff_meanTPCncl,mcrddiff_meanTPCncl;"; // Status number of clusters and findable clusters
 
  // Status: ITS:TPC-ITS matching efficiency 
  TStatToolkit::MakeCombinedAlias(treeMC,sCombinedStatus,doCheck, verbose);
  ::Info("InitTPCMCValidation","Done with aliases");
}


Bool_t InitTPCMCValidation(TString mcPeriod,  TString mcPass, TString anchorPeriod,  TString anchorPass, Int_t verbose,Int_t doCheck){

  pinfo=new AliExternalInfo(".","",verbose);
  trendingDraw= new AliTreeTrending;
  trendingDraw->SetDefaultStyle();
  
  treeMC = pinfo->GetTree("QA.TPC",mcPeriod,mcPass,"QA.TPC;QA.TRD;QA.TOF;QA.ITS");
  
  TTree * treeAnchorTPC;
  if(pinfo->GetTree("QA.TPC",anchorPeriod,anchorPass,"Logbook;QA.EVS")!=0) treeAnchorTPC=pinfo->GetTree("QA.TPC",anchorPeriod,anchorPass,"Logbook;QA.EVS");
  else{ 
      ::Error("InitTPCMCValidation","Failed to get QA.TPC tree");
      return kFALSE;
  }
  
  TTree * treeAnchorTRD0; 
  if(pinfo->GetTree("QA.TRD",anchorPeriod,anchorPass,"Logbook")!=0){
      treeAnchorTRD0=pinfo->GetTree("QA.TRD",anchorPeriod,anchorPass,"Logbook");
      ::Info("InitTPCMCValidation","QA.TRD tree entries: %d",treeAnchorTRD0->GetEntries());
  }
  else{
      ::Error("InitTPCMCValidation","Failed to get QA.TRD tree");
      return kFALSE;
  }
  
  TTree * treeAnchorITS0;
  if(pinfo->GetTree("QA.ITS",anchorPeriod,anchorPass,"Logbook")!=0){ 
      treeAnchorITS0 = pinfo->GetTree("QA.ITS",anchorPeriod,anchorPass,"Logbook");
      ::Info("InitTPCMCValidation","QA.ITS tree entries: %d",treeAnchorITS0->GetEntries());
  }
  else{
      ::Error("InitTPCMCValidation","Failed to get QA.ITS tree");
      return kFALSE;
  }

  treeMC->AddFriend(treeAnchorTPC,"TPC.Anchor");
  treeMC->AddFriend(treeAnchorTRD0,"TRD.Anchor");
  treeMC->AddFriend(treeAnchorITS0,"ITS.Anchor");
  makeTPCMCAlarms(treeMC,doCheck,verbose);
  TString sStatusbarVars ("ncl;dcarResol;itsEffStatus;mcrddiff;");
  TString sStatusbarNames("#(cl);dcar;itsEffStatus;mcrddiff;");
  TString sCriteria("(1):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)"); // or to just show vetos: (varname_PhysAcc&&varname_Warning)
  TString statusString[3];
  statusString[0] = sStatusbarVars;
  statusString[1] = sStatusbarNames;
  statusString[2] = sCriteria;

  trendingDraw->AddUserDescription(new TNamed("MC period",mcPeriod.Data()));
  trendingDraw->AddUserDescription(new TNamed("MC pass",mcPass.Data()));
  trendingDraw->AddUserDescription(new TNamed("Anchor period",anchorPeriod.Data()));
  trendingDraw->AddUserDescription(new TNamed("Anchor pass",anchorPass.Data()));
  trendingDraw->SetTree(treeMC);
  treeMC->SetAlias("tagID","run");
  treeMC->SetAlias("defaultcut","run==TPC.Anchor.run");

  if(trendingDraw->InitSummaryTrending(statusString,0.015,"defaultcut")) return kTRUE;
  else return kFALSE;
}

void MakeReport(const char* mcrddir){

  std::vector<Double_t> cRange{0.13,0.01,0.5,0.35};
  std::vector<Double_t> cRange2{0.13,0.01,0.5,0.3};
  std::vector<Double_t> cRange5{0.13,0.01,0.8,0.3};
  TMultiGraph *graph=0,*lines=0;
  

  MakePlot(treeMC,mcrddir,"matchingTPC-ITSEffe.png","Number of clusters",cRange,"","meanTPCncl;TPC.Anchor.meanTPCncl:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"matchingTPC-ITSEffe.png","Matching efficiency:MC/Anchor",cRange,"","QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffTOTPt02;QA.ITS.EffTOTPt1;QA.ITS.EffTOTPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffTOTPt02;ITS.Anchor.EffTOTPt1;ITS.Anchor.EffTOTPt10:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcItsMatch.png","TPC - ITS Match",cRange,"","tpcItsMatchHighPtA;TPC.Anchor.tpcItsMatchHighPtA;tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtC:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"rmsDCAMultSpartMCtoAnchor.png","DCA Resolution mult due MS:MC/Anchor",cRange,"","dcarAP1;dcarCP1;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"matchingTPC-ITSEffe_1.png","Matching efficiencyx",cRange,"","QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffoneSPDPt02;QA.ITS.EffoneSPDPt1;QA.ITS.EffoneSPDPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffoneSPDPt02;ITS.Anchor.EffoneSPDPt1;ITS.Anchor.EffoneSPDPt10:run","defaultcut","21;24;25;27;28;21;24;25;27;28","2;2;2;2;2;4;4;4;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"meanTPCNclRatioMCtoAnchor.png","Number of clusters MC/Anchor",cRange,"","meanTPCncl/TPC.Anchor.meanTPCncl;meanTPCnclF/TPC.Anchor.meanTPCnclF:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"matchingTPC-ITSEffe_2.png","Matching efficiency",cRange,"","QA.TPC.tpcItsMatchA/TPC.Anchor.tpcItsMatchA;QA.TPC.tpcItsMatchC/TPC.Anchor.tpcItsMatchC;QA.ITS.EffTOTPt02/ITS.Anchor.EffTOTPt02;QA.ITS.EffTOTPt1/ITS.Anchor.EffTOTPt1;QA.ITS.EffTOTPt10/ITS.Anchor.EffTOTPt10:run","defaultcut","21;24;25;27;28;21;24;25;27;28","1;2;4;3;6;2;4;4;4;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"meanVertX.png","meanVertX:run MC/Anchor",cRange,"","QA.TPC.meanVertX;TPC.Anchor.meanVertX:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"meanVertY.png","meanVertY:run MC/Anchor",cRange,"","QA.TPC.meanVertY;TPC.Anchor.meanVertY:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"meanVertZ.png","meanVertZ:run MC/Anchor",cRange,"","offsetdRA;TPC.Anchor.offsetdRA:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);
  MakePlot(treeMC,mcrddir,"offsetdRA.png","offsetdRA",cRange,"","offsetdRA;TPC.Anchor.offsetdRA:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"offsetdZA.png","offsetdZA",cRange,"","offsetdZA;TPC.Anchor.offsetdZA:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"offsetdRC.png","offsetdRC",cRange,"","offsetdRC;TPC.Anchor.offsetdRC:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"offsetdZC.png","offsetdZC",cRange,"","offsetdZC;TPC.Anchor.offsetdZC:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"offsetd_comb4.png","offsetd_comb4",cRange,"","offsetd_comb4;TPC.Anchor.offsetd_comb4:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"meanMultPos.png","meanMultPos",cRange,"","meanMultPos;TPC.Anchor.meanMultPos:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"meanMultNeg.png","meanMultNeg",cRange,"","meanMultNeg;TPC.Anchor.meanMultNeg:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"meanMult_comb2.png","meanMult_comb2",cRange,"","meanMult_comb2;TPC.Anchor.meanMult_comb2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcItsMatchA.png","tpcItsMatchA",cRange,"","tpcItsMatchA;TPC.Anchor.tpcItsMatchA:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcItsMatchHighPtA.png","tpcItsMatchHighPtA",cRange,"","tpcItsMatchHighPtA;TPC.Anchor.tpcItsMatchHighPtA:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcItsMatchC.png","tpcItsMatchC",cRange,"","tpcItsMatchC;TPC.Anchor.tpcItsMatchC:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcItsMatchHighPtC.png","tpcItsMatchHighPtC",cRange,"","tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtC:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcItsMatch_comb4.png","tpcItsMatch_comb4",cRange,"","tpcItsMatch_comb4;TPC.Anchor.tpcItsMatch_comb4:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"lambdaPull.png","lambdaPull",cRange,"","lambdaPull;TPC.Anchor.lambdaPull:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcConstrainPhiC.png","tpcConstrainPhiC",cRange,"","tpcConstrainPhiC;TPC.Anchor.tpcConstrainPhiC:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"yPull.png","yPull",cRange,"","yPull;TPC.Anchor.yPull:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"zPull.png","zPull",cRange,"","zPull;TPC.Anchor.zPull:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"itsTpcPulls_comb4.png","itsTpcPulls_comb4",cRange,"","itsTpcPulls_comb4;TPC.Anchor.itsTpcPulls_comb4:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcConstrainPhiA.png","tpcConstrainPhiA",cRange,"","tpcConstrainPhiA;TPC.Anchor.tpcConstrainPhiA:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcConstrainPhiC.png","tpcConstrainPhiC",cRange,"","tpcConstrainPhiC;TPC.Anchor.tpcConstrainPhiC:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"tpcConstrainPhi_comb2.png","tpcConstrainPhi_comb2",cRange,"","tpcConstrainPhi_comb2;TPC.Anchor.tpcConstrainPhi_comb2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"deltaPt.png","deltaPt",cRange,"","deltaPt;TPC.Anchor.deltaPt:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"deltaPtA.png","deltaPtA",cRange,"","deltaPtA;TPC.Anchor.deltaPtA:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"deltaPtC.png","deltaPtC",cRange,"","deltaPtC;TPC.Anchor.deltaPtC:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcarAP0.png","dcarAP0",cRange,"","dcarAP0;TPC.Anchor.dcarAP0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcarAP1.png","dcarAP1",cRange,"","dcarAP1;TPC.Anchor.dcarAP1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcarCP0.png","dcarCP0",cRange,"","dcarCP0;TPC.Anchor.dcarCP0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcarCP1.png","dcarCP1",cRange,"","dcarCP1;TPC.Anchor.dcarCP1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcarFitpar_comb4.png","dcarFitpar_comb4",cRange,"","dcarFitpar_comb4;TPC.Anchor.dcarFitpar_comb4:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_posA_0.png","dcar_posA_0",cRange,"","dcar_posA_0;TPC.Anchor.dcar_posA_0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_negA_0.png","dcar_negA_0",cRange,"","dcar_negA_0;TPC.Anchor.dcar_negA_0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_posC_0.png","dcar_posC_0",cRange,"","dcar_posC_0;TPC.Anchor.dcar_posC_0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_negC_0.png","dcar_negC_0",cRange,"","dcar_negC_0;TPC.Anchor.dcar_negC_0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_posA_1.png","dcar_posA_1",cRange,"","dcar_posA_1;TPC.Anchor.dcar_posA_1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_negA_1.png","dcar_negA_1",cRange,"","dcar_negA_1;TPC.Anchor.dcar_negA_1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_posC_1.png","dcar_posC_1",cRange,"","dcar_posC_1;TPC.Anchor.dcar_posC_1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_negC_1.png","dcar_negC_1",cRange,"","dcar_negC_1;TPC.Anchor.dcar_negC_1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_posA_2.png","dcar_posA_2",cRange,"","dcar_posA_2;TPC.Anchor.dcar_posA_2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_negA_2.png","dcar_negA_2",cRange,"","dcar_negA_2;TPC.Anchor.dcar_negA_2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_posC_2.png","dcar_posC_2",cRange,"","dcar_posC_2;TPC.Anchor.dcar_posC_2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcar_negC_2.png","dcar_negC_2",cRange,"","dcar_negC_2;TPC.Anchor.dcar_negC_2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_posA_0.png","dcaz_posA_0",cRange,"","dcaz_posA_0;TPC.Anchor.dcaz_posA_0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_negA_0.png","dcaz_negA_0",cRange,"","dcaz_negA_0;TPC.Anchor.dcaz_negA_0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_posC_0.png","dcaz_posC_0",cRange,"","dcaz_posC_0;TPC.Anchor.dcaz_posC_0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_negC_0.png","dcaz_negC_0",cRange,"","dcaz_negC_0;TPC.Anchor.dcaz_negC_0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_posA_1.png","dcaz_posA_1",cRange,"","dcaz_posA_1;TPC.Anchor.dcaz_posA_1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_negA_1.png","dcaz_negA_1",cRange,"","dcaz_negA_1;TPC.Anchor.dcaz_negA_1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_posC_1.png","dcaz_posC_1",cRange,"","dcaz_posC_1;TPC.Anchor.dcaz_posC_1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_negC_1.png","dcaz_negC_1",cRange,"","dcaz_negC_1;TPC.Anchor.dcaz_negC_1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_posA_2.png","dcaz_posA_2",cRange,"","dcaz_posA_2;TPC.Anchor.dcaz_posA_2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_negA_2.png","dcaz_negA_2",cRange,"","dcaz_negA_2;TPC.Anchor.dcaz_negA_2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_posC_2.png","dcaz_posC_2",cRange,"","dcaz_posC_2;TPC.Anchor.dcaz_posC_2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"dcaz_negC_2.png","dcaz_negC_2",cRange,"","dcaz_negC_2;TPC.Anchor.dcaz_negC_2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"iroc_A_side.png","iroc_A_side",cRange,"","iroc_A_side;TPC.Anchor.iroc_A_side:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"oroc_A_side.png","oroc_A_side",cRange,"","oroc_A_side;TPC.Anchor.oroc_A_side:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"iroc_C_side.png","iroc_C_side",cRange,"","iroc_C_side;TPC.Anchor.iroc_C_side:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"oroc_C_side.png","oroc_C_side",cRange,"","oroc_C_side;TPC.Anchor.oroc_C_side:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);   
  MakePlot(treeMC,mcrddir,"MIPattachSlopeA.png","MIPattachSlopeA",cRange,"","MIPattachSlopeA;TPC.Anchor.MIPattachSlopeA:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"MIPattachSlopeC*(-1).png","MIPattachSlopeC*(-1)",cRange,"","MIPattachSlopeC*(-1);TPC.Anchor.MIPattachSlopeC*(-1):run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"MIPattachSlope_comb2.png","MIPattachSlope_comb2",cRange,"","MIPattachSlope_comb2;TPC.Anchor.MIPattachSlope_comb2:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE);  
  MakePlot(treeMC,mcrddir,"electroMIPSeparation.png","electroMIPSeparation",cRange,"","electroMIPSeparation;TPC.Anchor.electroMIPSeparation:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE); 
  
  MakePlot(treeMC,mcrddir,"mcrddiff_statvar.png","mcrddiff Status Variable",cRange,"","mcrddiff_Warning;mcrddiff_Outlier;mcrddiff_PhysAcc:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,kTRUE); 

 { // MIP resolution  and separation
    trendingDraw->fWorkingCanvas->Clear(); 
    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"DCA Resolution mult due MS:MC/Anchor"); legend->SetBorderSize(0);
    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanMIPele/meanMIP;TPC.Anchor.meanMIPele/TPC.Anchor.meanMIP:run","defaultcut","25;21;25;21","2;4;4;4",1,0.75,6,legend);
    lines = TStatToolkit::MakeMultGraph(treeMC,"","meanMIPele_RobustMean/meanMIP_RobustMean;TPC.Anchor.meanMIPele_RobustMean/TPC.Anchor.meanMIP_RobustMean:run","defaultcut","25;21;25;21","2;4;4;4",1,0.75,6,0);
    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!");
    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} 
    else { TStatToolkit::DrawMultiGraph(graph,"alp");         
    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
    TStatToolkit::DrawMultiGraph(lines,"l");
    legend->SetFillStyle(0);legend->Draw();
   
    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/mipToEleSeparation.png");    
  }
}

void MakePlot(TTree * tree,const char* mcrddir, const char *fname, const char *LegendTitle, std::vector<Double_t>& legpos, const char *groupName, const char* expr, const char * cut, const char * markers, const char *colors, Bool_t drawSparse, Float_t msize, Float_t sigmaRange, Bool_t comp)
{
  TMultiGraph *graph=0;
  
   trendingDraw->fWorkingCanvas->Clear();  
   TLegend legend = TLegend(legpos[0],legpos[1],legpos[2],legpos[3],LegendTitle);
   legend.SetBorderSize(0);
   graph = TStatToolkit::MakeMultGraph(tree,groupName,expr,cut,markers,colors,drawSparse,msize,sigmaRange,&legend,comp);
   if(!graph){
       ::Error("MakePlot","No plot returned -> dummy plot!");
   } 
   else {
       TStatToolkit::DrawMultiGraph(graph,"alp");
       trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
       legend.SetFillStyle(0);
       legend.Draw();
   }

   trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/"+TString(fname)); 

}