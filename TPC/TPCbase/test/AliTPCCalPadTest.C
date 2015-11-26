#include "TTree.h"
#include "TFile.h"
#include "TError.h"
#include "TObjString.h"
#include "TMap.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGeoGlobalMagField.h"

#include "TStatToolkit.h"
#include "AliTPCCalPad.h"
#include "AliTPCcalibDB.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliGRPManager.h"
#include "AliMagF.h"

/** Functions to dump and compare AliTPCCalPad unit tests

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  .L AliTPCCalPad.C+g

  //dump first object before modifications
  DumpUnitTestTrees("outputFile=AliTPCCalPad_Test.old.root");

  //dump second object after modifications
  DumpUnitTestTrees("outputFile=AliTPCCalPad_Test.new.root");
  
  //run comparison
  CompareTrees("fileName1=AliTPCCalPad_Test.old.root;fileName2=AliTPCCalPad_Test.new.root")

  // === from command line:
  aliroot -b -q AliTPCCalPadTest.C+'(1,"outputFile=AliTPCCalPad_Test.old.root")'

  aliroot -b -q AliTPCCalPadTest.C+'(2,"fileName1=AliTPCCalPad_Test.old.root;fileName2=AliTPCCalPad_Test.new.root")'

*/

void DumpUnitTestTrees(AliTPCCalPad *calPad, TString fileName="");
void DumpUnitTestTrees(TString arguments);
void CompareSingle(const TTree &tree, const TString expr, const TString infoType);
void CompareTrees(TString arguments);
void ParseArguments(const TString arguments, TMap &argumentMap);
void InitCDB(Int_t run, TString storage);
TString CollectBranchNames(TTree &tree);
TString BranchListMod(const TString branchList, const TString friendName="ref", const TString modifier="-");
void ConfigOCDB(Int_t run, const char *ocdb="raw://");

void AliTPCCalPadTest(Int_t function, TString arguments)
{
  switch (function) {
    case 1:
      DumpUnitTestTrees(arguments);
      break;
    case 2:
      CompareTrees(arguments);
      break;
    default:
     Error("AliTPCCalPad","Function %d is unknown", function);
  }
}

//______________________________________________________________________________
void ParseArguments(const TString arguments, TMap &argumentMap)
{
  /// Parse the string and make a map out of it
  /// The assumed format is key1: value1; key2: value2; ...

  // clear map
  argumentMap.SetOwnerKeyValue();
  argumentMap.DeleteAll();

  if (arguments.CountChar('=')==0) return;

  TObjArray *arr=arguments.Tokenize(";");
  for (Int_t ipair=0; ipair<arr->GetEntriesFast(); ++ipair) {
    const TString &pair = ((TObjString*)arr->At(ipair))->String();
    TObjArray *arrKeyValue = pair.Tokenize("=");
    if (arrKeyValue->GetEntriesFast()!=2) {
      Error("ParseArguments","Not a real pair: '%s'", pair.Data());
      delete arrKeyValue;
      continue;
    }

    TString &key  =((TObjString*)arrKeyValue->At(0))->String();
    TString &value=((TObjString*)arrKeyValue->At(1))->String();
    key  .Remove(TString::kBoth, ' ');
    value.Remove(TString::kBoth, ' ');

    argumentMap.Add(new TObjString(key), new TObjString(value));
    delete arrKeyValue;
  }

  delete arr;
}

//______________________________________________________________________________
void DumpUnitTestTrees(AliTPCCalPad *calPad, TString fileName)
{
  /// Dump the unit test tree for 'calPad' to 'fileName'
  calPad->DumpUnitTestTrees(fileName);
}

//______________________________________________________________________________
void InitCDB(Int_t run, TString storage)
{
  ConfigOCDB(run, storage.Data());
}

//______________________________________________________________________________
void DumpUnitTestTrees(TString arguments)
{
  /// Dump the unit test trees
  /// By default use as input GetCETmean from AliTPCcalibDB
  /// Interpreted arguments:
  /// outputFile : output file name (default: "" -> given in the DumpUnitTestTrees function)
  /// runNumber  : run number to use for CDB manager (default: 0)
  /// storage    : CDB storage to use (default: 'local://$ALICE_ROOT/OCDB'
  /// 
  /// TODO: implement generic method to access any data member of AliTPCcalibDB if possible

  //============================================================================
  // ===| parse arguments |=====================================================
  //
  TMap argumentMap;
  ParseArguments(arguments, argumentMap);

  // --- parse output file name
  TObjString *outputFile = (TObjString*)argumentMap.GetValue("outputFile");
  TString fileName;
  if (outputFile) {
    fileName=outputFile->String();
  }

  // --- parse run number, 0 by default
  TObjString *runNumber  =  (TObjString*)argumentMap.GetValue("runNumber");
  Int_t run=0;
  if (runNumber) {
    run=runNumber->String().Atoi();
  }

  // --- set up CDB manager
  TObjString *ostorage    = (TObjString*)argumentMap.GetValue("storage");
  TString storage("local://$ALICE_ROOT/OCDB");
  if (ostorage) {
    storage=ostorage->GetName();
  }
  InitCDB(run, storage);

  // get calib db
  AliTPCcalibDB *db = AliTPCcalibDB::Instance();

  // get CalPad to be dumped
  AliTPCCalPad *calPad = db->GetCETmean();

  // dump CalPad
  DumpUnitTestTrees(calPad, fileName);
}

//______________________________________________________________________________
Double_t CompareSingle(TTree &tree, const TString expr, const TString infoType, 
                   const TStatToolkit::ENormType normType=TStatToolkit::kMax, const Double_t limit=1e-30)
{
  /// Compare single variables defined in 'expr' (split at ':')
  /// \param[in] tree          input tree
  /// \param[in] expr          variable expression, single variables separated by a ':'
  /// \param[in] infoType      type of information to show<br>
  ///                          can contain any of<br>
  ///                          draw: show summary canvases<br>
  ///                          info: show info massage (default)<br>
  /// \param[in] normType      Normalisation type to use
  /// \param[in] limit         distance (norm) limit below which the result is interpreted as identical
  /// \return                  fraction of good variables

  //============================================================================
  // ===| parse arguments |=====================================================
  //
  //Bool_t drawCanvas = infoType.Contains("draw");
  Bool_t dumpInfo   = infoType.Contains("info");


  // split variables
  TObjArray *arrVars = expr.Tokenize(":");

  const Int_t nvars=arrVars->GetEntriesFast();
  Int_t nGood=0;
  // loop over all variables and compare them
  for (Int_t ivar=0; ivar<nvars; ++ivar) {
    TString variable=arrVars->At(ivar)->GetName();
    Double_t norm=TMath::Abs(TStatToolkit::GetDistance(&tree, variable.Data(), "", normType));
    if (norm<limit) {
      if (dumpInfo) Info("CompareSingle","Norm of '%s' below limit (%.3g < %.3g) -> GOOD", variable.Data(), norm, limit);
      ++nGood;
    } else {
      if (dumpInfo) Error("CompareSingle","Norm of '%s' above limit (%.3g >= %.3g) -> BAD", variable.Data(), norm, limit);
    }
  }

  delete arrVars;
  return Double_t(nGood)/Double_t(nvars);
}

//______________________________________________________________________________
void CompareTrees(TString arguments)
{
  /// Compare two output CalPad trees
  /// \param[in] arguments can contain the following information:<br>
  ///                      fileName1: name of the first file<br>
  ///                      fileName2: name of the second file<br>
  ///                      infoType:  can contain any of<br>
  ///                                 draw: show summary canvases<br>
  ///                                 info: show info massage (default)<br>

  //============================================================================
  // ===| parse arguments |=====================================================
  //
  TMap argumentMap;
  ParseArguments(arguments, argumentMap);

  // --- fileName1
  TObjString *fileName1 = (TObjString*)argumentMap.GetValue("fileName1");
  if (!fileName1) {
    Error("CompareTrees","Required argument '%s' not given", "fileName1");
    return;
  }

  // --- fileName2
  TObjString *fileName2 = (TObjString*)argumentMap.GetValue("fileName2");
  if (!fileName2) {
    Error("CompareTrees","Required argument '%s' not given", "fileName2");
    return;
  }

  // --- infoType
  TObjString *ostrInfoType = (TObjString*)argumentMap.GetValue("infoType");
  TString infoType("info");
  if (ostrInfoType) infoType=ostrInfoType->String();

  //============================================================================
  // ===| get trees |===========================================================
  //

  // --- first file
  TFile f1(fileName1->GetName());
  if (!f1.IsOpen() || f1.IsZombie()) return;
  TTree *tValues1 = (TTree*)f1.Get("vars");
  TTree *tCalPad1 = (TTree*)f1.Get("calPads");

  // --- second file
  TFile f2(fileName2->GetName());
  if (!f2.IsOpen() || f2.IsZombie()) return;
  TTree *tValues2 = (TTree*)f2.Get("vars");
  TTree *tCalPad2 = (TTree*)f2.Get("calPads");

  // --- add as friend
  tValues1->AddFriend(tValues2, "ref");
  tCalPad1->AddFriend(tCalPad2, "ref");

  gROOT->cd();

  //============================================================================
  // ===| Do comparison |=======================================================
  //

  // ---| vars tree, single branches |------------------------------------------
  const TString branchesVars   = CollectBranchNames(*tValues1);
  const TString branchDiffVars = BranchListMod(branchesVars);

  const Double_t fracVarsGood  = CompareSingle(*tValues1, branchDiffVars, infoType);

  Info("CompareTrees","Single variable test of tree 'vars' has a good fraction of %.2f%%", fracVarsGood*100);

  // ---| calPad tree, selected branches |--------------------------------------
  const TString branchesCalPad  ("localFit.fElements:globalFit.fElements:globalFitRegions.fElements");
  const TString branchDiffCalPad=BranchListMod(branchesCalPad);

  const Double_t fracCalPadGood  = CompareSingle(*tCalPad1, branchDiffCalPad, infoType);
  
  Info("CompareTrees","CalPad variable test of tree 'vars' has a good fraction of %.2f%%", fracCalPadGood*100);
}

//______________________________________________________________________________
TString CollectBranchNames(TTree &tree)
{
  /// Collect all Branch Names in a draw string
  TString branchList;
  const TObjArray * branches = tree.GetListOfBranches();
  for (Int_t ibranch=0; ibranch<branches->GetEntries(); ++ibranch) {
    if (ibranch>0) branchList+=":";
    branchList+=branches->At(ibranch)->GetName();
  }
  return branchList;
}

//______________________________________________________________________________
TString BranchListMod(const TString branchList, const TString friendName/*="ref"*/, const TString modifier/*="-"*/)
{
  /// Add to each branch (separated by ':') 'modifier friendName.branch
  /// E.g.  branchList = a:b
  ///       friendName = ref
  ///       modifier   = -
  /// return a-ref.a:b-ref.b
  
  TString branchListMod;
  TObjArray *arr = branchList.Tokenize(":");
  for (Int_t ibranch=0; ibranch<arr->GetEntriesFast(); ++ibranch) {
    if (ibranch>0) branchListMod+=":";
    branchListMod+=arr->At(ibranch)->GetName();
    branchListMod+=modifier;
    branchListMod+=friendName;
    if (!friendName.EndsWith(".")) branchListMod+=".";
    branchListMod+=arr->At(ibranch)->GetName();
  }

  return branchListMod;
}

//______________________________________________________________________________
void ConfigOCDB(Int_t run, const char *ocdb/*="raw://"*/) 
{

  // OCDB
  printf("setting run to %d\n",run);
  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(run); 

  // magnetic field
  if ( !TGeoGlobalMagField::Instance()->GetField() ) {
    printf("Loading field map...\n");
    AliGRPManager grpMan;
    if( !grpMan.ReadGRPEntry() ) { 
      printf("Cannot get GRP entry\n"); 
    }
    if( !grpMan.SetMagField() ) { 
      printf("Problem with magnetic field setup\n"); 
    }
  }
  if ( !TGeoGlobalMagField::Instance()->GetField()){
    AliMagF::BMap_t smag = AliMagF::k5kG;
    Double_t bzfac = 1;
    AliMagF* magF= new AliMagF("Maps","Maps", bzfac, 1., smag);
    TGeoGlobalMagField::Instance()->SetField(magF);
  }

  // geometry
  printf("Loading geometry...\n");
  AliGeomManager::LoadGeometry();
  if( !AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC") ) {
    Error("ConfigOCDB","Problem with align objects\n"); 
  }
}

