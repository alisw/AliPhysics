#if !defined( __CINT__) || defined(__MAKECINT__)

#include <vector>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TH1F.h>
#include <AliCDBManager.h>
#include <AliTRDCalibra.h>

#endif
using namespace std;

TTree *AliTRDMergeTwoTrees(const char* variablecali){
  //
  // After having simulated and reconstructed events in subrepertories 000%d
  // this macro create a tree, "sum" of the tree in 0000/TRD.calibration.root and 0001/TRD.calibration.root
  // variablecali is the type: treeCH2d, treePH2d or treePRF2d
  //



  //Variables
  AliCDBManager *man = AliCDBManager::Instance();
  man->GetStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  AliTRDCalibra *calibra = AliTRDCalibra::Instance();

   //Add
  TTree *tree = calibra->Sum2Trees("/d/alice06/bailhache/aliroottrd/productiondefaulttout/trackletefficiency/sumtest/0000/TRD.calibration.root","/d/alice06/bailhache/aliroottrd/productiondefaulttout/trackletefficiency/sumtest/0001/TRD.calibration.root",variablecali);

  return tree;

}
