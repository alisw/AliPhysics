/*
  Test suit for the AliTreePlayer:

  aliroot -b -q   $ALICE_ROOT/../src/STAT/test/AliTreePlayerTest.C+ > AliTreePlayerTest.log
  cat AliTreePlayerTest.log | grep "AliTreePlayerTest\."

  Test should be integrated to AliRoot test suit (not cmake test as it need access to the external files)

  Output should be all test OK:
  I-AliTreePlayerTest.testSelectMetadata AND invariant: Test OK: N(Logbook&&(Stat))=N(Logbook&&(Stat&&Base))+N(Logbook&&(Stat&&(!Base)))   7=2+5
  I-AliTreePlayerTest.testSelectMetadata Negatioation test: Test OK:   N(!(A||B))=N((!A)&&(!B)) 73==73
  I-AliTreePlayerTest.testselectTreeInfo: Test OK:   N(A))==N(A&&B)&&N(A&&!B)) 46=6+40
  I-AliTreePlayerTest.testConvertTree: Test OK


  To run test interacivally, or test subsets:
     .L  $ALICE_ROOT/../src/STAT/test/AliTreePlayerTest.C+
     AliTreePlayerTest();

  Tests: see desciption in idnividual test
  void testSelectMetadata();
  void testselectTreeInfo();
  //void testselectWhatWhereOrderBy();
  void testConvertTree();
  in case




  AliTreePlayer::selectMetadata(treeLogbook, "[class==\"Logbook&&Time\"]",0)->Print();
  AliTreePlayer::selectWhatWhereOrderBy(treeTPC,"run:Logbook.run:QA.TPC.run:meanMIP:meanMIPele:meanMIPvsSector.fElements:fitMIP.fElements","meanMIP>0", "", 0,10,"html","qatpc.html");
  AliTreePlayer::selectWhatWhereOrderBy(treeTRD,"run:Logbook.run:QA.TRD.run:meanMIP:meanMIPele:meanMIPvsSector.fElements:fitMIP.fElements","meanMIP>0", "", 0,10,"json","qatpc.json");  
  AliTreePlayer::selectWhatWhereOrderBy(treeTRD,"run:Logbook.run:QA.TRD.run:meanMIP:meanMIPele:meanMIPvsSector.fElements:fitMIP.fElements","meanMIP>0", "", 0,10,"csv","qatpc.csv");  
*/

#include "TStatToolkit.h"
#include "Riostream.h"
#include <iostream>
#include "TSystem.h"
#include "TNamed.h"
#include "TFile.h"
#include "TTree.h"
#include "TPRegexp.h"
#include "TFriendElement.h"
#include "AliExternalInfo.h"
#include "TTreeFormula.h"
#include "TTreeFormulaManager.h"
#include "AliTreePlayer.h"


TTree * testTree=0;
void testSelectMetadata();
void testselectTreeInfo();
//void testselectWhatWhereOrderBy();
void testConvertTree();


void AliTreePlayerTest(){
  // test all  
  // Input data are provided by AliExternalInfo 
  // 
  AliExternalInfo info;
  TTree * treeLogbook = info.GetTree("Logbook","LHC15o","cpass1_pass1","QA.TPC;QA.TRD;QA.TOF;");
  TTree * treeTPC = info.GetTree("QA.TPC","LHC15o","cpass1_pass1","QA.TRD;QA.TPC;QA.TOC;QA.TOF;Logbook");
  TTree * treeTRD = info.GetTree("QA.TRD","LHC15o","cpass1_pass1","QA.TPC;QA.TRD;QA.TOF;Logbook;");
  //
  testTree=treeLogbook;
  testSelectMetadata();  // run test
  //
  testTree=treeTPC;
  testselectTreeInfo();
  //
  testConvertTree();
}

void testSelectMetadata(){
  //
  // test selectMetadata function
  //
  // 1.) Init trees
  TObjArray *array =0;
  if (testTree==NULL){
    AliExternalInfo info;
    testTree= info.GetTree("Logbook","LHC15o","cpass1_pass1","QA.TPC;QA.TRD;QA.TOF;");
    if (testTree==NULL){
      ::Error("AliTreePlayerTest.testSelectMetadata","Input not available");
    }
  }
  // 2.) Make test
  // 2.a) AND invariant test    N(A)=N(A&&B)+N(A&&(!B))
  Int_t statAll, stat0, stat1;
  array = AliTreePlayer::selectMetadata(testTree, "[class==\"Logbook&&(Stat)\"]",0);
  statAll=array->GetEntries();
  array = AliTreePlayer::selectMetadata(testTree, "[class==\"Logbook&&(Stat&&Base)\"]",0);
  stat0=array->GetEntries();
  array = AliTreePlayer::selectMetadata(testTree, "[class==\"Logbook&&(Stat&&(!Base))\"]",0);
  stat1=array->GetEntries();
  if (statAll!=stat0+stat1 || statAll==0){
    ::Error("AliTreePlayerTest.testSelectMetadata AND invariant","Test ERROR: N(Logbook&&(Stat))!=N(Logbook&&(Stat&&Base))+N(Logbook&&(Stat&&(!Base)))  ");
  }else{
    ::Info("AliTreePlayerTest.testSelectMetadata AND invariant","Test OK: N(Logbook&&(Stat))=N(Logbook&&(Stat&&Base))+N(Logbook&&(Stat&&(!Base)))\t %d=%d+%d",statAll,stat0,stat1);
  }
  // 2.b) Negation test    N(!(A||B))=N((!A)&&(!B))
  Int_t stat2b_0=AliTreePlayer::selectMetadata(testTree, "[class==\"!(Stat||Time)\"]",0)->GetEntries();
  Int_t stat2b_1=AliTreePlayer::selectMetadata(testTree, "[class==\"(!Stat&&!Time)\"]",0)->GetEntries();
  if (stat2b_0!=stat2b_1){
    ::Error("AliTreePlayerTest.testSelectMetadata Negatioation test","Test ERROR:   N(!(A||B))!=N((!A)&&(!B)) ");
  }else{
    ::Info("AliTreePlayerTest.testSelectMetadata Negatioation test","Test OK:   N(!(A||B))=N((!A)&&(!B)) %d==%d",stat2b_0,stat2b_1);
  }
}

void testselectTreeInfo(){
  //
  //  test using as and input the QA.TPC tree
  //
  Int_t stat2a_0 = AliTreePlayer::selectTreeInfo(testTree, "([.name:meanMIP] )",0)->GetEntries();
  Int_t stat2a_1 = AliTreePlayer::selectTreeInfo(testTree, "([.name:meanMIP] && ![.name:_])",0)->GetEntries();
  Int_t stat2a_2 = AliTreePlayer::selectTreeInfo(testTree, "([.name:meanMIP] && [.name:_])",0)->GetEntries();
  if (stat2a_0!=stat2a_1+stat2a_2){
    ::Error("AliTreePlayerTest.testselectTreeInfo","Test ERROR:   N(A))!=N(A&&B)&&N(A&&!B)) ");
  }else{
    ::Info("AliTreePlayerTest.testselectTreeInfo","Test OK:   N(A))==N(A&&B)&&N(A&&!B)) %d=%d+%d",stat2a_0,stat2a_1,stat2a_2);
  }
}


void testselectWhatWhereOrderByForTRD(){
  // aim of test  -
  // 1.)  code not crashing
  // 2.)  test unique IDs  (currently ther is a bug in )
  // 3.)  define properties in case of missing information
  //      Bug in the TTree::Scan
  //
  // Export algorithm:
  //    1.) Load trees
  //    2.) Set metadata (appropriate config file)
  //    3.) Export default table
  //
  AliExternalInfo info;
  TTree * treeTRD0 = info.GetTree("QA.TRD","LHC15o","cpass1_pass1","QA.TRD;Logbook;");
  TTree * treeTRD = info.GetTree("QA.TRD","LHC15o","cpass1_pass1","QA.TPC;QA.TRD;QA.TOF;Logbook;MonALISA.RCT");
  TTree * treeTPC = info.GetTree("QA.TPC","LHC15o","cpass1_pass1","QA.TPC;QA.TRD;QA.TOF;Logbook;MonALISA.RCT");
  //   check/filter  available data :  AliTreePlayer::selectTreeInfo(treeTRD, "([.name:bz])",0)->Print();

  TString what0="";
  what0+="Entry$:QA.TRD.run:QA.TPC.run:Logbook.run:Logbook.LHCperiod:QA.TPC.pass.GetName():QA.TPC.bz:Logbook.totalEventsPhysics:";  // add run properties
  AliTreePlayer::selectWhatWhereOrderBy(treeTRD,what0.Data(),"1", "", 0,1000,"json","qatrdtest0.json");
  AliTreePlayer::selectWhatWhereOrderBy(treeTPC,what0.Data(),"1", "", 0,1000,"json","qatpctest0.json");
  //
 //what+="MonALISA.RCT.tpc_value:MonALISA.RCT.trd_value:MonALISA.RCT.hlt_mode_value:";           // add MONALISA
 // what+="TPCTRDmatchEffPosAll:TPCTRDmatchEffNegAll:TRDTOFmatchEffPosAll:TRDTOFmatchEffPosAll";  // add eff.
 // what+="";

}

void testConvertTree() {
  AliExternalInfo info;
  TTree *treeTPC = info.GetTree("QA.TPC", "LHC15o", "cpass1_pass1");
  TTree *treeTRD = info.GetTree("QA.TRD", "LHC15o", "cpass1_pass1");
  TTree *treeTRD2 = info.GetTree("QA.TRD", "LHC15o", "cpass1_pass1");

  treeTPC->BuildIndex("run");
  treeTRD->BuildIndex("run");
  treeTRD->AddFriend(treeTRD2,"QA.TRD");
  treeTRD->AddFriend(treeTPC,"QA.TPC");

  AliTreePlayer::selectWhatWhereOrderBy(treeTRD,"run:QA.TRD.run:QA.TPC.run","1", "", 0,1000,"csv","sparsetest.csv");
  // test that the 1 and 2 columns are the same - count different collumns 
  //      cat sparsetest.csv | grep -v "run" | gawk  '{ sum+=($1!=$2) } END {print sum} '
  Int_t mismatch0 = (gSystem->GetFromPipe(" cat sparsetest.csv | grep -v \"run\" | gawk  '{ sum+=($1!=$2) } END {print sum} '")).Atoi(); // should be 0 count run mismatch
  //  cat sparsetest.csv | grep -v "run" | gawk  '{sum+=($1==prev); prev=$1;} END {print sum} '  
  Int_t mismatch1=(gSystem->GetFromPipe("cat sparsetest.csv | grep -v \"run\" | gawk  '{sum+=($1==prev); prev=$1;} END {print sum} '")).Atoi();
  if (mismatch0>0||mismatch1){
    ::Error("AliTreePlayerTest.testConvertTree","Test ERROR Run mismatch for sparse trees");
  }else{
    ::Info("AliTreePlayerTest.testConvertTree","Test OK");
  }
}





void reproduceIndexProblem(){
  //
  // In case name and title correspond each other  - Freind trees corraltion works properly
  // In old implementation of the TTreeSredirector branch name and title were different
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("testSparseTree.root","recreate");
  Int_t nruns=100;
  for (Int_t irun=0; irun<nruns;irun++) {
    (*pcstream)<<"runAll"<<"run="<<irun<<"\n";
    if (irun%2==0)  (*pcstream)<<"run2"<<"run="<<irun<<"\n";
    if (irun%3==0)  (*pcstream)<<"run3"<<"run="<<irun<<"\n";
    if (irun%4==0)  (*pcstream)<<"run4"<<"run="<<irun<<"\n";
    if (gRandom->Rndm()<0.5)  (*pcstream)<<"runR05"<<"run="<<irun<<"\n";
  }
  delete pcstream;
  pcstream = new TTreeSRedirector("testSparseTree.root");
  TList * tlist = pcstream->GetFile()->GetListOfKeys();
  TTree * tree  = (TTree*)pcstream->GetFile()->Get(tlist->At(0)->GetName());
  tree->BuildIndex("run");
  TTree * treeR  = (TTree*)pcstream->GetFile()->Get("runR05");
  treeR->BuildIndex("run");
  for (Int_t ikey=0; ikey<=tlist->GetLast(); ikey++) {
    TTree *treeF = (TTree*)pcstream->GetFile()->Get(tlist->At(ikey)->GetName());
    treeF->BuildIndex("run");
    tree->AddFriend(treeF, tlist->At(ikey)->GetName());
    treeR->AddFriend(treeF, tlist->At(ikey)->GetName());
  }
  printf("printed entries should be unioue per row and the same within the row  or not defined (in case of sparse)");
  treeR->Scan("run:run2.run:run3.run:run4.run:runR05.run","","",20);

}
