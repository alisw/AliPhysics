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

AliExternalInfo   *pinfo=0;
AliTreeTrending   *trendingDraw=0;  
TTree * treeMC;
  
void makeTPCMCAlarms(TTree * treeMC, Bool_t doCheck,Int_t verbose);
void tpcMCValidationStandard(TString mcPeriod,Int_t verbose,Int_t doCheck);
Bool_t InitTPCMCValidation(TString mcPeriod,  TString mcPass, TString anchorPeriod,  TString anchorPass, Int_t verbose,Int_t doCheck);
void MakeReport(const char* mcrddir);


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

    Double_t cRange[4]={0.13,0.01,0.5,0.35};
    Double_t cRange2[4]={0.13,0.01,0.5,0.3};
    Double_t cRange5[4]={0.13,0.01,0.8,0.3};
    TMultiGraph *graph=0,*lines=0;

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
  //  ==============================================================
  // 2.) Configure combined status. Default using logiacal OR of problems
  //  ==============================================================
  TString sCombinedStatus=";";
  sCombinedStatus+="ncl,TPC.Anchor.meanTPCncl,TPC.Anchor.meanTPCnclF;";               // Status number of clusters and findable clusters
  sCombinedStatus+="dcarResol,TPC.Anchor.dcarAP0,TPC.Anchor.dcarAP1,TPC.Anchor.dcarCP0,TPC.Anchor.dcarCP1;";  // Status: DCA resolution
  sCombinedStatus+="itsEffStatus,ITS.Anchor.EffoneSPDPt02,ITS.Anchor.EffoneSPDPt1,ITS.Anchor.EffoneSPDPt10,ITS.Anchor.EffTOTPt02,ITS.Anchor.EffTOTPt1,ITS.Anchor.EffTOTPt10;";  

  // Status: ITS:TPC-ITS matching efficiency 
  TStatToolkit::MakeCombinedAlias(treeMC,sCombinedStatus,doCheck, verbose);
}


Bool_t InitTPCMCValidation(TString mcPeriod,  TString mcPass, TString anchorPeriod,  TString anchorPass, Int_t verbose,Int_t doCheck){

  pinfo=new AliExternalInfo(".","",verbose);
  trendingDraw= new AliTreeTrending;
  trendingDraw->SetDefaultStyle();
  
  treeMC = pinfo->GetTree("QA.TPC",mcPeriod,mcPass,"QA.TPC;QA.TRD;QA.TOF;QA.ITS");
//    treeMC = pinfo->GetTree("QA.ITS", anchorPeriod, anchorPass, "Logbook");
  
  TTree * treeAnchorTPC;
  if(pinfo->GetTree("QA.TPC",anchorPeriod,anchorPass,"Logbook;QA.EVS")!=0) treeAnchorTPC=pinfo->GetTree("QA.TPC",anchorPeriod,anchorPass,"Logbook;QA.EVS");
  else{ 
      ::Error("InitTPCMCValidation","Failed to get QA.TPC tree");
      return kFALSE;
  }
  
  TTree * treeAnchorTRD0; 
  if(pinfo->GetTree("QA.TRD",anchorPeriod,anchorPass,"Logbook")!=0) treeAnchorTRD0=pinfo->GetTree("QA.TRD",anchorPeriod,anchorPass,"Logbook");
  else{
      ::Error("InitTPCMCValidation","Failed to get QA.TRD tree");  
      return kFALSE;
  }
  
  TTree * treeAnchorITS0;
  if(pinfo->GetTree("QA.ITS",anchorPeriod,anchorPass,"Logbook")!=0){ 
      treeAnchorITS0 = pinfo->GetTree("QA.ITS",anchorPeriod,anchorPass,"Logbook");
//      treeAnchorITS0 =pinfo->GetTree("QA.TPC",mcPeriod,mcPass,"QA.TPC;QA.TRD;QA.TOF;QA.ITS");
//  treeAnchorITS0->SetBranchStatus("*",0);
//  treeAnchorITS0->SetBranchStatus("Eff*",1);

  }
  else{
      ::Error("InitTPCMCValidation","Failed to get QA.ITS tree");
      return kFALSE;
  }

  treeMC->AddFriend(treeAnchorTPC,"TPC.Anchor");
  treeMC->AddFriend(treeAnchorTRD0,"TRD.Anchor");
  treeMC->AddFriend(treeAnchorITS0,"ITS.Anchor");
  ::Info("here","%d < %d",treeMC->GetEntries(),treeAnchorITS0->GetEntries());
//  treeMC->Scan("ITS.Anchor.EffoneSPDPt02");
  makeTPCMCAlarms(treeMC,doCheck,verbose);
  TString sStatusbarVars ("ncl;dcarResol;itsEffStatus;");
  TString sStatusbarNames("#(cl);dcar;itsEffStatus;");
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
  if(trendingDraw->InitSummaryTrending(statusString,0.015,"1")) return kTRUE;
  else return kFALSE;
}

void MakeReport(const char* mcrddir){

  Double_t cRange[4]={0.13,0.01,0.5,0.35};
  Double_t cRange2[4]={0.13,0.01,0.5,0.3};
  Double_t cRange5[4]={0.13,0.01,0.8,0.3};
  TMultiGraph *graph=0,*lines=0;
  

  
  { // matching efficiency
    trendingDraw->fWorkingCanvas->Clear(); 
    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Matching efficiency:MC/Anchor"); legend->SetBorderSize(0);
    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcItsMatchHighPtA;TPC.Anchor.tpcItsMatchHighPtA;tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtC:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,legend);
    
    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);
    ::Error("tpcMCValidation","No plot returned -> dummy plot!");}
    else { TStatToolkit::DrawMultiGraph(graph,"alp"); 
    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
    legend->Draw();

    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/matchingTPC-ITSEffe3.png");    
  }
  
  { //Number of clusters comparison MC/real data
    trendingDraw->fWorkingCanvas->Clear();  
    TLegend *legend = new TLegend(cRange2[0],cRange2[1],cRange2[2],cRange2[3],"Number of clusters MC/Anchor"); legend->SetBorderSize(0);
    legend->SetNColumns(2);
    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanTPCncl;TPC.Anchor.meanTPCncl:run","defaultcut","25;21","2;4",1,0.75,5,legend);
    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");}
    else { TStatToolkit::DrawMultiGraph(graph,"alp");          
    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
    legend->Draw();

    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/meanTPCNclMCtoAnchor.png");
  }
//  
//  {// matching efficiency
//   trendingDraw->fWorkingCanvas->Clear(); 
//   TLegend *legend = new TLegend(cRange5[0],cRange5[1],cRange5[2],cRange5[3],"Matching efficiency:MC/Anchor"); legend->SetBorderSize(0);
//   legend->SetNColumns(5);
//   graph = TStatToolkit::MakeMultGraph(treeMC,"","QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffTOTPt02;QA.ITS.EffTOTPt1;QA.ITS.EffTOTPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffTOTPt02;ITS.Anchor.EffTOTPt1;ITS.Anchor.EffTOTPt10:run","defaultcut","21;24;25;27;28;21;24;25;27;28","2;2;2;2;2;4;4;4;4;4",1,1.5,5,legend);
//   if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");}
//   else { TStatToolkit::DrawMultiGraph(graph,"alp");
//   trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//   legend->Draw(); 
//
//   trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/matchingTPC-ITSEffe.png");    
//  }
//  
//  { // DCA resolution at high pt
//   trendingDraw->fWorkingCanvas->Clear();  
//   TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"DCA Resolution at q/pt=0:MC/Anchor"); legend->SetBorderSize(0);
//   graph = TStatToolkit::MakeMultGraph(treeMC,"","dcarAP0;dcarCP0;TPC.Anchor.dcarAP0;TPC.Anchor.dcarCP0:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,legend);
//   if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);
//   ::Error("tpcMCValidation","No plot returned -> dummy plot!");} 
//   else { TStatToolkit::DrawMultiGraph(graph,"alp");
//   trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//   legend->Draw();
//
//   trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/rmsDCAAt1pt0MCtoAnchor.png");    
// }
//
// { // DCA resolution MS part
//   trendingDraw->fWorkingCanvas->Clear();  
//   TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"DCA Resolution mult due MS:MC/Anchor"); legend->SetBorderSize(0);
//   graph = TStatToolkit::MakeMultGraph(treeMC,"","dcarAP1;dcarCP1;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP1:run","defaultcut","25;21;25;21","2;2;4;4",1,0.75,6,legend);
//   if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);
//   ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");
//   trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//
//   trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/rmsDCAMultSpartMCtoAnchor.png");    
// }
// 
// { // matching efficiency
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange5[0],cRange5[1],cRange5[2],cRange5[3],"Matching efficiency:pass1_lowIR/pass3_lowIR_pidfix"); legend->SetBorderSize(0);
//    legend->SetNColumns(5);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffoneSPDPt02;QA.ITS.EffoneSPDPt1;QA.ITS.EffoneSPDPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffoneSPDPt02;ITS.Anchor.EffoneSPDPt1;ITS.Anchor.EffoneSPDPt10:run","defaultcut","21;24;25;27;28;21;24;25;27;28","2;2;2;2;2;4;4;4;4;4",1,1.5,5,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);
//    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp"); 
//    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}    
//    legend->Draw();
//
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/matchingTPC-ITSEffe_1.png");    
// }
// { 
//   trendingDraw->fWorkingCanvas->Clear();  
//   TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Number of clusters MC/Anchor"); legend->SetBorderSize(0);
//   graph = TStatToolkit::MakeMultGraph(treeMC,"","meanTPCncl/TPC.Anchor.meanTPCncl;meanTPCnclF/TPC.Anchor.meanTPCnclF:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//   if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//   legend->Draw();    
//  
//   trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/meanTPCNclRatioMCtoAnchor.png");    
// }
// { // matching efficiency
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange5[0],cRange5[1],cRange5[2],cRange[3],"Matching efficiency:pass1_lowIR/pass3_lowIR_pidfix"); legend->SetBorderSize(0);
//    legend->SetNColumns(2);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","QA.TPC.tpcItsMatchA/TPC.Anchor.tpcItsMatchA;QA.TPC.tpcItsMatchC/TPC.Anchor.tpcItsMatchC;QA.ITS.EffTOTPt02/ITS.Anchor.EffTOTPt02;QA.ITS.EffTOTPt1/ITS.Anchor.EffTOTPt1;QA.ITS.EffTOTPt10/ITS.Anchor.EffTOTPt10;:run","defaultcut","21;24;25;27;28;21;24;25;27;28","1;2;4;3;6;2;4;4;4;4;4",1,1.5,8,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/matchingTPC-ITSEffe_2.png");    
// }
// { // MIP resolution  and separation
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"DCA Resolution mult due MS:MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanMIPele/meanMIP;TPC.Anchor.meanMIPele/TPC.Anchor.meanMIP:run","defaultcut","25;21;25;21","2;4;4;4",1,0.75,6,legend);
//    lines = TStatToolkit::MakeMultGraph(treeMC,"","meanMIPele_RobustMean/meanMIP_RobustMean;TPC.Anchor.meanMIPele_RobustMean/TPC.Anchor.meanMIP_RobustMean:run","defaultcut","25;21;25;21","2;4;4;4",1,0.75,6,0);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    TStatToolkit::DrawMultiGraph(lines,"l");
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/mipToEleSeparation.png");    
//  }
// 
// 
// 
// 
// 
// 
// //////////////NEW PLOTS//////////////////////////////
// 
// 
// 
// 
//
// { // Separation Power
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"(meanMIPele-meanMIP)/(0.5*(resolutionMIP*meanMIP+resolutionMIPele*meanMIPele)):run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","PIDSepPow_comb2;TPC.Anchor.PIDSepPow_comb2:run","defaultcut","21;25","4;2",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/SeparationPower_vs_run.png");    
//  } 
//  
//  
//   { // VertX
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"meanVertX:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","QA.TPC.meanVertX;TPC.Anchor.meanVertX:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/meanVertX.png");    
//  } 
//  
//     { // VertY
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"meanVertY:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanVertY;TPC.Anchor.meanVertY:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/meanVertY.png");    
//  } 
//  
//  
//       { // VertZ
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"meanVertZ:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanVertZ;TPC.Anchor.meanVertZ:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/meanVertZ.png");    
//  } 
//  
//  
//     { // offsetdRA                             //Errors not working!
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"offsetdRA:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","offsetdRA;TPC.Anchor.offsetdRA:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/offsetdRA.png");    
//  } 
//  
//  
//       { // offsetdZA
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"offsetdZA:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","offsetdZA;TPC.Anchor.offsetdZA:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/offsetdZA.png");    
//  } 
//  
//  
//    { // offsetdRC
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"offsetdRC:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","offsetdRC;TPC.Anchor.offsetdRC:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/offsetdRC.png");    
//  } 
//  
//  
//    { // offsetdZC
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"offsetdZC:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","offsetdZC;TPC.Anchor.offsetdZC:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/offsetdZC.png");    
//  } 
//  
//  
//    { // offsetd_comb4
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"offsetd_comb4:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","offsetd_comb4;TPC.Anchor.offsetd_comb4:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/offsetd_comb4.png");    
//  } 
//  
//  
//  
//  
//  
//  
//      { // meanMultPos
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"meanMultPos:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanMultPos;TPC.Anchor.meanMultPos:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/meanMultPos.png");    
//  } 
//  
//      { // meanMultNeg
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"meanMultNeg:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanMultNeg;TPC.Anchor.meanMultNeg:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/meanMultNeg.png");  
//  } 
//
//  
//    { // meanMult_comb2
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"meanMult_comb2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanMult_comb2;TPC.Anchor.meanMultNeg:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/meanMult_comb2.png");  
//  }  
//  
//  
//  
//  
//  
//  
//  
//  
//  
//  
//  
//        { // tpcItsMatchA
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"tpcItsMatchA:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcItsMatchA;TPC.Anchor.tpcItsMatchA:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/tpcItsMatchA.png");    
//  } 
//  
//      { // tpcItsMatchHighPtA
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"tpcItsMatchHighPtA:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcItsMatchHighPtA;TPC.Anchor.tpcItsMatchHighPtA:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/tpcItsMatchHighPtA.png");  
//  } 
//  
//  
//          { // tpcItsMatchC
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"tpcItsMatchC:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcItsMatchC;TPC.Anchor.tpcItsMatchC:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/tpcItsMatchC.png");    
//  } 
//  
//      { // tpcItsMatchHighPtC
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"tpcItsMatchHighPtC:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtC:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/tpcItsMatchHighPtC.png");  
//  } 
//
//  
//    { // tpcItsMatch_comb4
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"tpcItsMatch_comb4:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcItsMatch_comb4;TPC.Anchor.tpcItsMatch_comb4:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/tpcItsMatch_comb4.png");  
//  }  
//  
//  
//  
//  
//  
//  
//    /****** ITS-TPC matching quality  ******/
//        { // lambdaPull
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"lambdaPull:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","lambdaPull;TPC.Anchor.lambdaPull:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/lambdaPull.png");    
//  } 
//  
//      { // tpcConstrainPhiC
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"tpcConstrainPhiC:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcConstrainPhiC;TPC.Anchor.tpcConstrainPhiC:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/tpcConstrainPhiC.png");  
//  } 
//  
//  
//          { // yPull
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"yPull:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","yPull;TPC.Anchor.yPull:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/yPull.png");    
//  } 
//  
//      { // zPull
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"zPull:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","zPull;TPC.Anchor.zPull:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/zPull.png");  
//  } 
//
//  
//    { // itsTpcPulls_comb4
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"itsTpcPulls_comb4:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","itsTpcPulls_comb4;TPC.Anchor.itsTpcPulls_comb4:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/itsTpcPulls_comb4.png");  
//  }  
//  
//  
//  
//  /****** pullPhi for TPC Constrain  ******/  
//        { // tpcConstrainPhiA
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"tpcConstrainPhiA:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcConstrainPhiA;TPC.Anchor.tpcConstrainPhiA:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/tpcConstrainPhiA.png");    
//  } 
//  
//      { // tpcConstrainPhiC
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"tpcConstrainPhiC:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcConstrainPhiC;TPC.Anchor.tpcConstrainPhiC:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/tpcConstrainPhiC.png");  
//  } 
//  
//  
//          { // tpcConstrainPhi_comb2
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"tpcConstrainPhi_comb2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcConstrainPhi_comb2;TPC.Anchor.tpcConstrainPhi_comb2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/tpcConstrainPhi_comb2.png");    
//  } 
//  
//  
//  
//    /****** 1/Pt  ******/
//   { // deltaPt
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"deltaPt:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","deltaPt;TPC.Anchor.deltaPt:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/deltaPt.png");    
//  } 
//  
//      { // deltaPtA
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"deltaPtA:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","deltaPtA;TPC.Anchor.deltaPtA:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/deltaPtA.png");  
//  } 
//  
//  
//      { // deltaPtC
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"deltaPtC:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","deltaPtC;TPC.Anchor.deltaPtC:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/deltaPtC.png");  
//  } 
//  
//  
//  
//    /****** DCAr fitting parameters  ******/
//        { // dcarAP0
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcarAP0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcarAP0;TPC.Anchor.dcarAP0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcarAP0.png");    
//  } 
//  
//      { // dcarAP1
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcarAP1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcarAP1;TPC.Anchor.dcarAP1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcarAP1.png");  
//  } 
//  
//  
//    { // dcarCP0
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcarCP0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcarCP0;TPC.Anchor.dcarCP0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcarCP0.png");    
//  } 
//  
//    { // dcarCP1
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcarCP1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcarCP1;TPC.Anchor.dcarCP1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcarCP1.png");  
//  } 
//
//  
//    { // dcarFitpar_comb4
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcarFitpar_comb4:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcarFitpar_comb4;TPC.Anchor.dcarFitpar_comb4:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcarFitpar_comb4.png");  
//  }   
//  
//  
//  
//  
//    ////////////////////////////////////////////////////////////////////////
//  //test DCAR plots
//  //DCAr first parameter _0
//  
//    { // dcar_posA_0:run
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_posA_0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_posA_0;TPC.Anchor.dcar_posA_0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_posA_0.png");    
//    } 
//  
//    { // dcar_negA_0
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_negA_0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_negA_0;TPC.Anchor.dcar_negA_0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_negA_0.png");  
//    } 
//
//  
//    { // dcar_posC_0
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_posC_0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_posC_0;TPC.Anchor.dcar_posC_0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_posC_0.png");  
//    }
//    
//    { // dcar_negC_0
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_negC_0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_negC_0;TPC.Anchor.dcar_negC_0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_negC_0.png");  
//    }    
//    
//    
//    
//  ////////////////////////////////////////////////////////////////////////
//  //DCAr second parameter _1
//  
//    
//    { // dcar_posA_1:run
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_posA_1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_posA_1;TPC.Anchor.dcar_posA_1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_posA_1.png");    
//    } 
//  
//    { // dcar_negA_1
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_negA_1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_negA_1;TPC.Anchor.dcar_negA_1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_negA_1.png");  
//    } 
//
//  
//    { // dcar_posC_1
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_posC_1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_posC_1;TPC.Anchor.dcar_posC_1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_posC_1.png");  
//    }
//    
//    { // dcar_negC_1
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_negC_1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_negC_1;TPC.Anchor.dcar_negC_1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_negC_1.png");  
//    }    
//    
//    
//////////////////////////////////////////////////////////////////////////
//  //DCAr second parameter _2
//  
//    
//    { // dcar_posA_2:run
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_posA_2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_posA_2;TPC.Anchor.dcar_posA_2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_posA_2.png");    
//    } 
//  
//    { // dcar_negA_2
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_negA_2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_negA_2;TPC.Anchor.dcar_negA_2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_negA_2.png");  
//    } 
//
//  
//    { // dcar_posC_2
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_posC_2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_posC_2;TPC.Anchor.dcar_posC_2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_posC_2.png");  
//    }
//    
//    { // dcar_negC_2
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcar_negC_2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcar_negC_2;TPC.Anchor.dcar_negC_2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcar_negC_2.png");  
//    }   
//    
//    
//  ////////////////////////////////////////////////////////////////////////
//  //test DCAR plots
//  //DCAr first parameter _0
//  
//    { // dcaz_posA_0:run
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_posA_0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_posA_0;TPC.Anchor.dcaz_posA_0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_posA_0.png");    
//    } 
//  
//    { // dcaz_negA_0
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_negA_0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_negA_0;TPC.Anchor.dcaz_negA_0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_negA_0.png");  
//    } 
//
//  
//    { // dcaz_posC_0
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_posC_0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_posC_0;TPC.Anchor.dcaz_posC_0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_posC_0.png");  
//    }
//    
//    { // dcaz_negC_0
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_negC_0:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_negC_0;TPC.Anchor.dcaz_negC_0:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_negC_0.png");  
//    }    
//    
//    
//    
//  ////////////////////////////////////////////////////////////////////////
//  //DCAr second parameter _1
//  
//    
//    { // dcaz_posA_1:run
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_posA_1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_posA_1;TPC.Anchor.dcaz_posA_1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_posA_1.png");    
//    } 
//  
//    { // dcaz_negA_1
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_negA_1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_negA_1;TPC.Anchor.dcaz_negA_1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_negA_1.png");  
//    } 
//
//  
//    { // dcaz_posC_1
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_posC_1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_posC_1;TPC.Anchor.dcaz_posC_1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_posC_1.png");  
//    }
//    
//    { // dcaz_negC_1
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_negC_1:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_negC_1;TPC.Anchor.dcaz_negC_1:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_negC_1.png");  
//    }    
//    
//    
//////////////////////////////////////////////////////////////////////////
//  //DCAr second parameter _2
//  
//    
//    { // dcaz_posA_2:run
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_posA_2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_posA_2;TPC.Anchor.dcaz_posA_2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_posA_2.png");    
//    } 
//  
//    { // dcaz_negA_2
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_negA_2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_negA_2;TPC.Anchor.dcaz_negA_2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_negA_2.png");  
//    } 
//
//  
//    { // dcaz_posC_2
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_posC_2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_posC_2;TPC.Anchor.dcaz_posC_2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_posC_2.png");  
//    }
//    
//    { // dcaz_negC_2
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"dcaz_negC_2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcaz_negC_2;TPC.Anchor.dcaz_negC_2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/dcaz_negC_2.png");  
//    }   
//    
//    
//    
//  ///////////////////////////////////////////////////////////////////////////////////////////////
//  // Plot Occupancy IROC, OROC, A,C side  
//    
//    { // iroc_A_side:run
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"iroc_A_side:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","iroc_A_side;TPC.Anchor.iroc_A_side:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/iroc_A_side.png");    
//    } 
//  
//    { // oroc_A_side
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"oroc_A_side:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","oroc_A_side;TPC.Anchor.oroc_A_side:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/oroc_A_side.png");  
//    } 
//
//  
//    { // iroc_C_side:run
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"iroc_C_side:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","iroc_C_side;TPC.Anchor.iroc_C_side:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/iroc_C_side.png");    
//    } 
//  
//    { // oroc_C_side
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"oroc_C_side:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","oroc_C_side;TPC.Anchor.oroc_C_side:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/oroc_C_side.png");  
//    }    
//    
//    
//    
//    /****** attachment parameters for A and C side ******/
//    
//    { // MIPattachSlopeA
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"MIPattachSlopeA:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","MIPattachSlopeA;TPC.Anchor.MIPattachSlopeA:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/MIPattachSlopeA.png");  
//    }      
//
//    { // MIPattachSlopeC*(-1)
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"MIPattachSlopeC*(-1):run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","MIPattachSlopeC*(-1);TPC.Anchor.MIPattachSlopeC*(-1):run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/MIPattachSlopeC.png");  
//    }     
//    
//    { // MIPattachSlope_comb2
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"MIPattachSlope_comb2:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","MIPattachSlope_comb2;TPC.Anchor.MIPattachSlope_comb2:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/MIPattachSlope_comb2.png");  
//    } 
//    
//  /****** electron and MIPs separation ******/
//    
//    
//    { // electroMIPSeparation
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"electroMIPSeparation:run MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","electroMIPSeparation;TPC.Anchor.electroMIPSeparation:run","defaultcut","25;21","2;4",1,0.75,6,legend);
//    if(!graph){ legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Plotting error!!"); legend->SetBorderSize(0);    ::Error("tpcMCValidation","No plot returned -> dummy plot!");} else { TStatToolkit::DrawMultiGraph(graph,"alp");           trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);}
//    legend->Draw();
//   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/electroMIPSeparation.png");  
//    } 
    
    
    
    
    
    
    
    
    
    
    
    
      /****** electron and MIPs separation ******/
    
    
//    { // Runs
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Run Numbers MC/Anchor"); legend->SetBorderSize(0);
//    treeMC->Draw("run");
//    treeMC->SetFillStyle(3001);
//    treeMC->SetFillColor(kBlue);
//    treeMC->Draw("TPC.Anchor.run","","same");
//    legend->Draw();
////   
//    trendingDraw->fWorkingCanvas->SaveAs(TString(mcrddir)+"/run_numbers.png");  
//    }  

}