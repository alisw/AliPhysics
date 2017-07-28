/*
  .L $NOTES/JIRA/ATO-83/code/tpcMCValidationStandardQA.C+
  TString mcPeriod="LHC15k1a1";
  TString mcPass="passMC";
  TString anchorPeriod="LHC15o";
  TString anchorPass="pass3_lowIR_pidfix";
  
  InitTPCMCValidation("LHC15k1a1","passMC","LHC15o", "pass3_lowIR_pidfix",0,0);
 

*/ 
//gSystem->AddIncludePath("-I$ALICE_ROOT/include/"); //couldn't add include path in .rootrc

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
void InitTPCMCValidation(TString mcPeriod,  TString mcPass, TString anchorPeriod,  TString anchorPass, Int_t verbose,Int_t doCheck);
void MakeReport();


void tpcMCValidation(const char* MCper ="LHC15k1a1"){

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


InitTPCMCValidation(MCper,"passMC",AnchProdName, AnchPassName,0,0);


Double_t cRange[4]={0.13,0.01,0.5,0.35};
Double_t cRange2[4]={0.13,0.01,0.5,0.3};
Double_t cRange5[4]={0.13,0.01,0.8,0.3};
TMultiGraph *graph=0,*lines=0;

MakeReport();
//trendingDraw->fWorkingCanvas->Clear();  
//TLegend *legend2 = new TLegend(cRange2[0],cRange2[1],cRange2[2],cRange2[3],"Number of clusters: MC/Anchor"); legend2->SetBorderSize(0);
//legend2->SetNColumns(2);
//graph = TStatToolkit::MakeMultGraph(treeMC,"","meanTPCncl;TPC.Anchor.meanTPCncl:run","1","25;21","2;4",1,0.75,5,legend2);
//TStatToolkit::DrawMultiGraph(graph,"alp");
//legend2->Draw();
//trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//trendingDraw->fWorkingCanvas->SaveAs("meanTPCNclMCtoAnchor.png"); 
////trendingDraw->fWorkingCanvas->Draw();
//
//{trendingDraw->fWorkingCanvas->Clear(); 
//TLegend *legend = new TLegend(cRange5[0],cRange5[1],cRange5[2],cRange5[3],"Matching efficiency: MC/Anchor"); legend->SetBorderSize(0);
//legend->SetNColumns(5);
//graph = TStatToolkit::MakeMultGraph(treeMC,"","QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffTOTPt02;QA.ITS.EffTOTPt1;QA.ITS.EffTOTPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffTOTPt02;ITS.Anchor.EffTOTPt1;ITS.Anchor.EffTOTPt10:run","1","21;24;25;27;28;21;24;25;27;28","2;2;2;2;2;4;4;4;4;4",1,1.5,5,legend);
//TStatToolkit::DrawMultiGraph(graph,"alp");
//legend->Draw(); 
//trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//trendingDraw->fWorkingCanvas->SaveAs("matchingTPC-ITSEffe.png");    
////trendingDraw->fWorkingCanvas->Draw(); 
//}
//
//{
//trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange5[0],cRange5[1],cRange5[2],cRange5[3],"Matching efficiency:pass1_lowIR/pass3_lowIR_pidfix"); legend->SetBorderSize(0);
//    legend->SetNColumns(5);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffoneSPDPt02;QA.ITS.EffoneSPDPt1;QA.ITS.EffoneSPDPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffoneSPDPt02;ITS.Anchor.EffoneSPDPt1;ITS.Anchor.EffoneSPDPt10:run","1","21;24;25;27;28;21;24;25;27;28","2;2;2;2;2;4;4;4;4;4",1,1.5,5,legend);
//    TStatToolkit::DrawMultiGraph(graph,"alp");    
//    legend->Draw();
//    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//    trendingDraw->fWorkingCanvas->SaveAs("matchingTPC-ITSEffe1.png"); 
//}
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


void InitTPCMCValidation(TString mcPeriod,  TString mcPass, TString anchorPeriod,  TString anchorPass, Int_t verbose,Int_t doCheck){
  /*
   // Example setting for the code tests:
   TString mcPeriod="LHC15k1a1";
   TString mcPass="passMC";
   TString anchorPeriod="LHC15o";
   TString anchorPass="pass3_lowIR_pidfix";
   
   tpcMCValidationStandard("LHC15k1a1","passMC","LHC15o", "pass3_lowIR_pidfix"  )
   
  */
  // 0.) Init helper classes 
  pinfo=new AliExternalInfo(".","",verbose);
  trendingDraw= new AliTreeTrending;
  trendingDraw->SetDefaultStyle();
  //
  // 1.) Configure input trees
  //
  treeMC = pinfo->GetTree("QA.TPC",mcPeriod,mcPass,"QA.TPC;QA.TRD;QA.TOF;QA.ITS");
  TTree * treeAnchorTPC = pinfo->GetTree("QA.TPC",anchorPeriod,anchorPass,"Logbook;QA.EVS");
  TTree * treeAnchorTRD0 = pinfo->GetTree("QA.TRD",anchorPeriod,anchorPass,"Logbook");
  TTree * treeAnchorITS0 = pinfo->GetTree("QA.ITS",anchorPeriod,anchorPass,"Logbook");

  treeMC->AddFriend(treeAnchorTPC,"TPC.Anchor");
  treeMC->AddFriend(treeAnchorTRD0,"TRD.Anchor");
  treeMC->AddFriend(treeAnchorITS0,"ITS.Anchor");
  //
  // 2.) Configure alarms and status bars
  //
  makeTPCMCAlarms(treeMC,doCheck,verbose);
  TString sStatusbarVars ("ncl;dcarResol;itsEffStatus;");
  TString sStatusbarNames("#(cl);dcar;itsEffStatus;");
  TString sCriteria("(1):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)"); // or to just show vetos: (varname_PhysAcc&&varname_Warning)
  TString statusString[3];
  statusString[0] = sStatusbarVars;
  statusString[1] = sStatusbarNames;
  statusString[2] = sCriteria;
  //
  trendingDraw->AddUserDescription(new TNamed("MC period",mcPeriod.Data()));
  trendingDraw->AddUserDescription(new TNamed("MC pass",mcPass.Data()));
  trendingDraw->AddUserDescription(new TNamed("Anchor period",anchorPeriod.Data()));
  trendingDraw->AddUserDescription(new TNamed("Anchor pass",anchorPass.Data()));
  trendingDraw->SetTree(treeMC);
  treeMC->SetAlias("tagID","run");
  trendingDraw->InitSummaryTrending(statusString); 
}

void MakeReport(){
  //
  // 3.) make a report
  //
  Double_t cRange[4]={0.13,0.01,0.5,0.35};
  Double_t cRange2[4]={0.13,0.01,0.5,0.3};
  Double_t cRange5[4]={0.13,0.01,0.8,0.3};
  TMultiGraph *graph=0,*lines=0;
  { //Number of clusters comparison MC/real data
    trendingDraw->fWorkingCanvas->Clear();  
    TLegend *legend2 = new TLegend(cRange2[0],cRange2[1],cRange2[2],cRange2[3],"Number of clusters MC/Anchor"); legend2->SetBorderSize(0);
    legend2->SetNColumns(2);
    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanTPCncl;TPC.Anchor.meanTPCncl:run","1","25;21","2;4",1,0.75,5,legend2);
    TStatToolkit::DrawMultiGraph(graph,"alp");
    legend2->Draw();
    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
    trendingDraw->fWorkingCanvas->SaveAs("mcrd_com/meanTPCNclMCtoAnchor.png");
  }
//  
//  { //
//    trendingDraw->fWorkingCanvas->Clear();  
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Number of clusters MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanTPCncl/TPC.Anchor.meanTPCncl;meanTPCnclF/TPC.Anchor.meanTPCnclF :run","1","25;21","2;4",1,0.75,6,legend);
//    TStatToolkit::DrawMultiGraph(graph,"alp");
//    legend->Draw();    
//    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//    trendingDraw->fWorkingCanvas->SaveAs("meanTPCNclRatioMCtoAnchor.png");    
//  }
//  
//  { // DCA resolution at high pt
//    trendingDraw->fWorkingCanvas->Clear();  
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"DCA Resolution at q/pt=0:MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcarAP0;dcarCP0;TPC.Anchor.dcarAP0;TPC.Anchor.dcarCP0:run","1","25;21;25;21","2;2;4;4",1,0.75,6,legend);
//    TStatToolkit::DrawMultiGraph(graph,"alp");
//    legend->Draw();
//    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//    trendingDraw->fWorkingCanvas->SaveAs("rmsDCAAt1pt0MCtoAnchor.png");    
//  }
//
//  { // DCA resolution MS part
//    trendingDraw->fWorkingCanvas->Clear();  
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"DCA Resolution mult due MS:MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","dcarAP1;dcarCP1;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP1:run","1","25;21;25;21","2;2;4;4",1,0.75,6,legend);
//    TStatToolkit::DrawMultiGraph(graph,"alp");
//    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//    trendingDraw->fWorkingCanvas->SaveAs("rmsDCAMultSpartMCtoAnchor.png");    
//  }
//
//  { // MIP resolution  and separation
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"DCA Resolution mult due MS:MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","meanMIPele/meanMIP;TPC.Anchor.meanMIPele/TPC.Anchor.meanMIP:run","1","25;21;25;21","2;4;4;4",1,0.75,6,legend);
//    lines = TStatToolkit::MakeMultGraph(treeMC,"","meanMIPele_RobustMean/meanMIP_RobustMean;TPC.Anchor.meanMIPele_RobustMean/TPC.Anchor.meanMIP_RobustMean:run","1","25;21;25;21","2;4;4;4",1,0.75,6,0);
//    TStatToolkit::DrawMultiGraph(graph,"alp");
//    TStatToolkit::DrawMultiGraph(lines,"l");
//    legend->Draw();
//    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//    trendingDraw->fWorkingCanvas->SaveAs("mipToEleSeparation.png");    
//  }
//
//  { // matching efficiency
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange[0],cRange[1],cRange[2],cRange[3],"Matching efficiency:MC/Anchor"); legend->SetBorderSize(0);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","tpcItsMatchHighPtA;TPC.Anchor.tpcItsMatchHighPtA;tpcItsMatchHighPtC;TPC.Anchor.tpcItsMatchHighPtC:run","1","25;21;25;21","2;2;4;4",1,0.75,6,legend);
//    TStatToolkit::DrawMultiGraph(graph,"alp");
//    legend->Draw();
//    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//    trendingDraw->fWorkingCanvas->SaveAs("matchingTPC-ITSEffe3.png");    
//  }

  { // matching efficiency
    trendingDraw->fWorkingCanvas->Clear(); 
    TLegend *legend = new TLegend(cRange5[0],cRange5[1],cRange5[2],cRange5[3],"Matching efficiency:pass1_lowIR/pass3_lowIR_pidfix"); legend->SetBorderSize(0);
    legend->SetNColumns(5);
    graph = TStatToolkit::MakeMultGraph(treeMC,"","QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffTOTPt02;QA.ITS.EffTOTPt1;QA.ITS.EffTOTPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffTOTPt02;ITS.Anchor.EffTOTPt1;ITS.Anchor.EffTOTPt10:run","1","21;24;25;27;28;21;24;25;27;28","2;2;2;2;2;4;4;4;4;4",1,1.5,5,legend);
    TStatToolkit::DrawMultiGraph(graph,"alp");
    legend->Draw(); 
    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
    trendingDraw->fWorkingCanvas->SaveAs("mcrd_com/matchingTPC-ITSEffe.png");    
  }

//  { // matching efficiency
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange5[0],cRange5[1],cRange5[2],cRange5[3],"Matching efficiency:pass1_lowIR/pass3_lowIR_pidfix"); legend->SetBorderSize(0);
//    legend->SetNColumns(5);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","QA.TPC.tpcItsMatchA;QA.TPC.tpcItsMatchC;QA.ITS.EffoneSPDPt02;QA.ITS.EffoneSPDPt1;QA.ITS.EffoneSPDPt10;TPC.Anchor.tpcItsMatchA;TPC.Anchor.tpcItsMatchC;ITS.Anchor.EffoneSPDPt02;ITS.Anchor.EffoneSPDPt1;ITS.Anchor.EffoneSPDPt10:run","1","21;24;25;27;28;21;24;25;27;28","2;2;2;2;2;4;4;4;4;4",1,1.5,5,legend);
//    TStatToolkit::DrawMultiGraph(graph,"alp");    
//    legend->Draw();
//    trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//    trendingDraw->fWorkingCanvas->SaveAs("matchingTPC-ITSEffe_1.png");    
//  }
  
  
//  { // matching efficiency
//    trendingDraw->fWorkingCanvas->Clear(); 
//    TLegend *legend = new TLegend(cRange5[0],cRange5[1],cRange5[2],cRange[3],"Matching efficiency:pass1_lowIR/pass3_lowIR_pidfix"); legend->SetBorderSize(0);
//    legend->SetNColumns(2);
//    graph = TStatToolkit::MakeMultGraph(treeMC,"","QA.TPC.tpcItsMatchA/TPC.Anchor.tpcItsMatchA;QA.TPC.tpcItsMatchC/TPC.Anchor.tpcItsMatchC;QA.ITS.EffTOTPt02/ITS.Anchor.EffTOTPt02;QA.ITS.EffTOTPt1/ITS.Anchor.EffTOTPt1;QA.ITS.EffTOTPt10/ITS.Anchor.EffTOTPt10;:run","1","21;24;25;27;28;21;24;25;27;28","1;2;4;3;6;2;4;4;4;4;4",1,1.5,8,legend);
//    TStatToolkit::DrawMultiGraph(graph,"alp");
//    legend->Draw();
//     trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);
//     trendingDraw->fWorkingCanvas->SaveAs("matchingTPC-ITSEffe_2.png");    
//  }


}



/*
void debugFormulas(){  
  //
  // Bug found in parsing of the TTreeFormulas
  //    Shortes script reproducing problem to be created and put into the UnitTest of the AliTreePplayerTest.C
  //
  //
  // In some cases TTreeFromaul does not recognize properly branches
  //      Unit test to check Tree  match
  //  Problem analogue observed in the past and Unit test Writen in  AliTreePlayerTest.C:reproduceIndexProblem()
  // 

  treeMC->GetFriend("QA.TRD").Scan("run:TPCTRDmatchEffNegAll");
  treeMC->GetFriend("Anchor.QA.TRD").Scan("run:TPCTRDmatchEffNegAll")
  // This is not fine - 
  treeMC->Scan("run:QA.TRD.run:Anchor.QA.TRD.run:QA.TRD.TPCTRDmatchEffNegAll:Anchor.QA.TRD.TPCTRDmatchEffNegAll","","");
  // This looks fine - but exact test / comparison neeede
  treeMCTRD->Scan("run:QA.TRD.run:Anchor.QA.TRD.run:QA.TRD.TPCTRDmatchEffNegAll:Anchor.QA.TRD.TPCTRDmatchEffNegAll","","");
  //
  dump
    treeMC->GetFriend("QA.TRD")->Scan("run:TPCTRDmatchEffNegAll","","colsize=30");        >treeMC.qa.trd.list0
    treeMC->GetFriend("Anchor.QA.TRD")->Scan("run:TPCTRDmatchEffNegAll","","colsize=30"); >treeMC.anchor.qa.trd.list0
    //
    treeMC->Scan("run:QA.TRD.TPCTRDmatchEffNegAll","","colsize=30"); >treeMC.qa.trd.list1
    treeMC->Scan("run:Anchor.QA.TRD.TPCTRDmatchEffNegAll","","colsize=30"); >treeMC.anchor.qa.trd.list1     
    //
    treeMCTRD->Scan("run:QA.TRD.TPCTRDmatchEffNegAll","","colsize=30"); >treeMCTRD.qa.trd.list1
    treeMCTRD->Scan("run:Anchor.QA.TRD.TPCTRDmatchEffNegAll","","colsize=30"); >treeMCTRD.anchor.qa.trd.list1 
    //
    for a in `ls *.list*`; do cat $a | gawk '{print $4"\t"$6}'> $a.Table; mv $a.Table $a; done 

    meld  treeMC.qa.trd.list0  treeMCTRD.qa.trd.list1 treeMC.qa.trd.list1    # first 2 collumns are the same 3 columns is wrong
    meld  treeMC.anchor.qa.trd.list0 treeMC.anchor.qa.trd.list1 treeMCTRD.anchor.qa.trd.list1  #  all are the same but in anchor we have more runs ...
    meld  treeMC.anchor.qa.trd.list1 treeMC.qa.trd.list1   #anchor tree used instead of the MC

  
}

*/
