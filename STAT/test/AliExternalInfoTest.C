//
// Unit test for the AliExternalInfo
// 


/*
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/install/include -I$ALICE_ROOT/STEER\
  -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW  -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/TPCBase  -I$ALICE_ROOT/TPC/TPCRec -I$ALICE_ROOT/TPC/TPCCalib -I$ALICE_PHYSICS/../src/PWGPP/TPC/  -I$ALICE_ROOT/../src/STAT/");

  
  .L $ALICE_ROOT/../src/STAT/test/AliExternalInfoTest.C+
  

  TestMCProduction(); 
  TestProductionAccess();
*/


#include "AliExternalInfo.h"
#include "TTree.h"
#include "TMath.h"

void TestMCProduction();
void TestProductionAccess();

void AliExternalInfoTest(){
  //
  //
  //
  TestMCProduction();
  TestProductionAccess();
}

void TestMCProduction(){
  //
  // Check availaibility fo the MC production infomation
  //     - MC production tree
  //     - MC tree 
  // Compatibility of the trees   
  //     
  ::Info("AliExternalInfo.TestMCProduction","Begin");
  // 1.) Load MC production information
  AliExternalInfo info;
  TTree * treeMCProd= info.GetTree("MonALISA.ProductionMC","","","MonALISA.MC");
  TTree * treeMC= info.GetTree("MonALISA.MC","","","MonALISA.ProductionMC");
  //
  // 2.) Check availaibility
  //
  Int_t nentries=treeMC->Draw("prodName!=MonALISA.ProductionMC.prodName","1");
  Double_t mean=(nentries>0) ? TMath::Mean(nentries, treeMC->GetV1()):1;
  if (nentries>0) {
    ::Info("AliExternalInfo.TestMCProduction.","TestStat(prodName==MonALISA.ProductionMC.prodName)%d-OK",nentries);
  }else{
    ::Error("AliExternalInfo.TestMCProduction.","TestStat(prodName==MonALISA.ProductionMC.prodName)%d-FAILED",nentries);
  }
  if (TMath::Abs(mean)<0.5/(1+nentries)) {
    ::Info("AliExternalInfo.TestMCProduction","TestEqual(prodName==MonALISA.ProductionMC.prodName)%6.6f-OK",mean);
  }else{
    ::Error("AliExternalInfo.TestMCProduction","TestEqual(prodName==MonALISA.ProductionMC.prodName)%6.6f-FAILED",mean);
  }
  ::Info("AliExternalInfo.TestMCProduction","End");

}


void TestProductionAccess(){
  //
  // Check availability of the AliExternalInfo
  //
  AliExternalInfo info;
  TTree * tree = info.GetTree("MonALISA.ProductionCycle", "", "");
  //
  if (tree->GetEntries()>0){
    ::Info("AliExternalInfo.TestProductionAccess","Nproductions=%d >0 - OK", tree->GetEntries());
  }else{
    ::Error("AliExternalInfo.TestProductionAccess","Nproductions=%d - FAILED", tree->GetEntries());
    return;
  };
  // check and dump some counters
  tree->Draw("strstr(Description,\"for\")");
}
