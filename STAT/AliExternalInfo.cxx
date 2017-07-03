#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>

#include "AliLog.h"
#include "TROOT.h"

#include "AliExternalInfo.h"

#include "TObjArray.h"
#include "TString.h"
#include "TWebFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TStatToolkit.h"
#include "TLeaf.h"
#include "TVirtualIndex.h"
#include "TPRegexp.h" 
ClassImp(AliExternalInfo)

using namespace std;

const TString AliExternalInfo::fgkDefaultConfig="$ALICE_ROOT/STAT/Macros/AliExternalInfo.cfg";

AliExternalInfo::AliExternalInfo(TString localStorageDirectory, TString configLocation, Int_t verbose/*, Bool_t copyToLocalStorage*/) :
                                /*fCopyDataToLocalStorage(copyToLocalStorage),*/
                                TObject(),
				fVerbose(verbose),
                                fConfigLocation(configLocation),
                                fLocalStorageDirectory(localStorageDirectory),
                                fConfigMap(),
                                fTree(0x0),
                                fChain(new TChain()),
                                fChainMap(),
                                fMaxCacheSize(-1)
{
  // use default cache path from Env variable if specified 
  if (gSystem->Getenv("AliExternalInfoCache")!=NULL){
    if (localStorageDirectory.Length()<2)   fLocalStorageDirectory=gSystem->Getenv("AliExternalInfoCache");
  }
  ReadConfig(configLocation, fVerbose);
}

AliExternalInfo::~AliExternalInfo() {}


/// Reads the configuration files. Lines beginning with an '#' are ignored.
/// 
/// Use the format which is in the config.cfg by default. Adding ressources like the ones already
/// there should work without problems.
/// \param configLocation - semicolon separated configuration list
///                       - applied as in CSS 
/// \param verbose        - verbosity
void AliExternalInfo::ReadConfig( TString configLocation, Int_t verbose){

  // Prepend config files on first load
  if (fConfigLocation.Length()==0){
    // Check if local config exists
    if(!gSystem->AccessPathName("AliExternalInfo.cfg")){ // false if file exists
      configLocation=AliExternalInfo::fgkDefaultConfig+";AliExternalInfo.cfg;"+configLocation;
    } else{
      configLocation=AliExternalInfo::fgkDefaultConfig+configLocation;
    }
  }

  TObjArray *configArray = configLocation.Tokenize(";");
  fConfigLocation+=configLocation;
  fConfigLocation+=";";

  Int_t nConfig=configArray->GetEntries();
  if (nConfig==0){
    ::Error("AliExternalInfo::ReadConfig","Invalid configuration description %s",configLocation.Data());
    return;
  }

  for (Int_t iConfig=0; iConfig<nConfig; iConfig++){
    TString cName=configArray->At(iConfig)->GetName();
    if (cName=="default") cName="$ALICE_ROOT/STAT/Macros/AliExternalInfo.cfg";    
    if (cName.Length()==0) continue;
    TString configFileName=gSystem->ExpandPathName(cName.Data());
    if (verbose>0){
      ::Info("AliExternalInfo::ReadConfig","Path: %s\t%s",cName.Data(), configFileName.Data());
    }
    if (gSystem->AccessPathName(configFileName)!=0) {   // be aware of strange convention for the gSystem->AccessPathName - 0 mean it is OK
      AliError(TString::Format("Could not find config file '%s'", configFileName.Data()));
      const TString defaultConfigFileName=gSystem->ExpandPathName(fgkDefaultConfig);
      if (defaultConfigFileName!=configFileName) {
	AliError("Using default config file instead");
	configFileName=defaultConfigFileName;
      }
    }
    
    std::ifstream configFile(configFileName);
    if (!configFile.is_open()) {
      AliError(TString::Format("Could not open config file '%s'", configFileName.Data()));
      return;
    }
    
    //
    std::string line;
    while (std::getline(configFile, line)){
      TString temp_line(line.c_str()); // Use TString for easier handling
      
      if (temp_line.EqualTo("")) continue; // Ignore empty lines
      // temp_line = temp_line.Strip(TString::EStripType::kBoth, ' '); // Strip trailing and leading spaces
      temp_line = temp_line.Strip(TString::kBoth, ' '); // Strip trailing and leading spaces
      if (temp_line.First('#') == 0) continue; // If line starts with a '#' treat is as a comment
      
      TObjArray arrTokens = (*(temp_line.Tokenize(' ')));
      const TString key(arrTokens.At(0)->GetName());
      const TString value(arrTokens.At(1)->GetName());
      
      fConfigMap[key] = value;
    }
  }
  delete configArray;
  return;
}

/// Prints out the config which was read in previously. Useful to check if anything went wrong
void AliExternalInfo::PrintConfig(){
  // Loop through the map (Would be much easier in c++11)
  AliInfo("User defined resources are\n");
  // looping over map with const_iterator
  typedef std::map<TString,TString>::const_iterator it_type;
  for(it_type iterator = fConfigMap.begin(); iterator != fConfigMap.end(); ++iterator) {
    std::cout << iterator->first << " " << iterator->second << "\n";
  }
  return;
}


/// Sets up all variables according to period, pass and type. Extracts information from the config file
void AliExternalInfo::SetupVariables(TString& internalFilename, TString& internalLocation, Bool_t& resourceIsTree, TString& pathStructure, \
                                     TString& detector, TString& rootFileName, TString& treeName, const TString& type, const TString& period, const TString& pass, TString & indexName){
  // Check if resource is a tree in a root file or not
  TString lpass=pass;
  if (fConfigMap[type + ".nopass"].Contains("true")) lpass="";
  pathStructure = CreatePath(type, period, lpass);

  // if (fConfigMap.count(type + ".treename") > 0) resourceIsTree = kTRUE;
  // else resourceIsTree = kFALSE;
  if (type.Contains("MonALISA") == kTRUE) resourceIsTree = kFALSE;
  else resourceIsTree = kTRUE;

  // To distinguish different detector QA you have to add a <det>_ to the trending.root. Here we check the detector!
  if (type.Contains("QA")) {
   Int_t firstDotOfType(type.First('.') + 1);
   Int_t lastCharOfType(type.Length() - 1);
   detector = type(firstDotOfType, lastCharOfType) + "_";
//    std::cout << "DETECTOR: " << detector << std::endl;
  }

  rootFileName = fConfigMap[type + ".filename"];
  treeName     = fConfigMap[type + ".treename"];
  indexName= fConfigMap[type + ".indexname"];
  if (indexName.Length()<=0) indexName="run";
  // Create the local path where to store the information of the resource
  internalLocation += pathStructure;
  AliInfo(TString::Format("Information will be stored/retrieved in/from %s", internalLocation.Data()));

  if (!(period.Last('*') == period.Length() - 1) || !(pass.Last('*') == pass.Length() - 1) || period.Length() == 0){
//     std::cout << "mkdir " << internalLocation.Data() << std::endl;
    gSystem->mkdir(internalLocation.Data(), kTRUE);
  }

  // Resulting internal path to the file
  if (resourceIsTree) internalFilename = internalLocation + detector + rootFileName; // e.g data/<year>/<period>/<pass>/<det>_trending.root
  else internalFilename = internalLocation + rootFileName; // e.g data/<year>/<period>/<pass>/MC.root

  return;
};

/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// \param period Period, e.g. 'LHC15f'
/// \param pass E.g. 'pass2' or 'passMC'
/// Implements the functionality to download the ressources in the specified storage. If it is not
/// in the form of a tree it creates one. You can use this function or the functions available in
/// the class definition as an abbrevation
/// \return If downloading and creation of tree was successful true, else false
Bool_t AliExternalInfo::Cache(TString type, TString period, TString pass){
  AliInfo(TString::Format("Caching of %s %s from %s in start path %s", period.Data(), pass.Data(), type.Data(), fLocalStorageDirectory.Data()));

  // initialization of local variables
  TString internalFilename = ""; // Resulting path to the file
  TString internalLocation = fLocalStorageDirectory; // Gets expanded in this function to the right directory
  TString externalLocation = "";
  Bool_t resourceIsTree = kFALSE;
  TString detector = "";
  TString rootFileName = "";
  TString treeName = "";
  TString pathStructure = "";
  TString indexName=""; 
  TString oldIndexName= fConfigMap[type + ".indexname"];  // rename index branch to avoid incositencies (bug in ROOT - the same index branch name requeired) 

  // initialization of external variables
  externalLocation = fConfigMap[type + ".location"];

  // Setting up all the local variables
  SetupVariables(internalFilename, internalLocation, resourceIsTree, pathStructure, detector, rootFileName, treeName, type, period, pass,indexName);

  // Checking if resource needs to be downloaded
  const Bool_t downloadNeeded = IsDownloadNeeded(internalFilename, type);

  if (downloadNeeded == kTRUE){
    // Download resources in the form of .root files in a tree
    if (resourceIsTree == kTRUE){
      externalLocation += pathStructure + rootFileName;
      AliInfo(TString::Format("Information retrieved from: %s", externalLocation.Data()));

      // Check if external location is a http address or locally accessible
//       std::cout << externalLocation(0, 4) << std::endl;
      TFile *file = TFile::Open(externalLocation);
      if (file && !file->IsZombie()){ // Checks if webresource is available
        if (file->Cp(internalFilename)) {
          AliInfo("Caching successful");
          return kTRUE;
        }
        else {
          AliError("Copying to internal location failed");
          return kFALSE;
        }
      }
      else {
        AliError("Ressource not available");
        return kFALSE;
      }
      delete file;
    }
    else {
      //Set up external path with period and pass if necessary
      if (period != "" && pass != ""){
        externalLocation = TString::Format(externalLocation.Data(), period.Data(), pass.Data());
      }
      else if (period != "" && pass == ""){
        externalLocation = TString::Format(externalLocation.Data(), period.Data());
      }

      TString mifFilePath = ""; // Gets changed in Curl command
      TString command = Curl(mifFilePath, internalLocation, rootFileName, externalLocation);

      std::cout << command << std::endl;
      gSystem->Exec(command.Data());
      if (oldIndexName.Length()==0){
	gSystem->Exec(TString::Format("cat %s | sed -l 1 s/raw_run/run/ |  sed -l 1 s/RunNo/run/ > %s",mifFilePath.Data(),  (mifFilePath+"RunFix").Data())); // use standrd run number IDS
      }else{
	gSystem->Exec(TString::Format("cat %s | sed -l 1 s/%s/%s/  > %s",mifFilePath.Data(), oldIndexName.Data(), indexName.Data(),  (mifFilePath+"RunFix").Data())); // use standrd run number IDS
      } 

      gSystem->GetFromPipe(TString::Format("cat %s  | sed s_\\\"\\\"_\\\"\\ \\\"_g | sed s_\\\"\\\"_\\\"\\ \\\"_g > %s",  (mifFilePath+"RunFix").Data(),  (mifFilePath+"RunFix").Data()).Data());
      // Store it in a tree inside a root file
      TFile tempfile(internalFilename, "RECREATE");
      tempfile.cd();
      TTree tree(treeName, treeName);

      if ( (tree.ReadFile(mifFilePath, "", '\"')) > 0) {
        if (fVerbose>1) AliInfo("Successfully read in tree");
      }
      else {
        AliError("-- Error while reading tree");
        return kFALSE;
      }

      tree.Write();
      tempfile.Close();
      if (fVerbose>0) AliInfo(TString::Format("Write tree to file: %s", internalFilename.Data()));
      return kTRUE;
    }
  }
  else {// downloadIsNeeded == kFALSE
    return kTRUE;
  }
}

/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// \param period Period, e.g. 'LHC15f'
/// \param pass E.g. 'pass2' or 'passMC'
/// Returns the tree with the information from the corresponding resource
/// \return TTree* with corresponding resource
TTree* AliExternalInfo::GetTree(TString type, TString period, TString pass, Int_t buildIndex){
  TString internalFilename = ""; // Resulting path to the file
  TString internalLocation = fLocalStorageDirectory; // Gets expanded in this function to the right directory
  TString externalLocation = "";
  Bool_t resourceIsTree = kFALSE;
  TString detector = "";
  TString rootFileName = "";
  TString treeName = "";
  TString pathStructure = "";
  TString indexName=""; 
  TString metadataMacro=fConfigMap[type + ".metadataMacro"];

  TTree* tree = 0x0;

  // Setting up all the local variables
  SetupVariables(internalFilename, internalLocation, resourceIsTree, pathStructure, detector, rootFileName, treeName, type, period, pass,indexName);

  //std::cout << "internalFilename: " << internalFilename << " rootFileName: " << rootFileName << std::endl;
  if (fVerbose>1) AliInfo(TString::Format("Caching start internalFileName\t%s\trootFileName\t%s", internalFilename.Data(), rootFileName.Data()).Data());
  if (gSystem->AccessPathName(internalFilename.Data()) == kTRUE) {
    if (Cache(type, period, pass) == kFALSE) {
      AliError(TString::Format("Caching failed internalFileName\t%s\trootFileName\t%s", internalFilename.Data(), rootFileName.Data()).Data());
      return tree;
    }
  }else{
    Bool_t downloadNeeded = IsDownloadNeeded(internalFilename, type);
    if (downloadNeeded){
      AliInfo(TString::Format("Caching %s", internalFilename.Data()).Data());
      if (Cache(type, period, pass) == kFALSE) {
	AliError(TString::Format("Caching failed internalFileName\t%s\trootFileName\t%s", internalFilename.Data(), rootFileName.Data()).Data());
	return tree;
      }
    }    
  }

  // Creating and returning the tree from the file
  TFile* treefile = new TFile(internalFilename.Data());

  // ---| loop over possible tree names |---
  TObjArray *arr = treeName.Tokenize(",");
  for (Int_t iname=0; iname<arr->GetEntriesFast(); ++iname) {
    tree = dynamic_cast<TTree*>( treefile->Get(arr->At(iname)->GetName()) );
    if (tree) break;
  }
  delete arr;
  TTreeSRedirector::FixLeafNameBug(tree);
  if (tree != 0x0) {
    if (fVerbose>1) AliInfo(TString::Format("Successfully read %s/%s",internalFilename.Data(), tree->GetName()));
    if (buildIndex==1) BuildIndex(tree, type);
  } else {
    AliError("Error while reading tree");
  }

  const TString cacheSize=fConfigMap[type + ".CacheSize"];
  Long64_t cache=cacheSize.Atoll();
  if (fMaxCacheSize>0) {
    if (cache>0) {
      cache=TMath::Min(fMaxCacheSize, cache);
    } else {
      cache=fMaxCacheSize;
    }
  }
  if (cache>0) tree->SetCacheSize(cache);
  //

  if (metadataMacro.Length()>0){  // rename branch  with index if specified in configuration file
    if (fVerbose>1) printf("Processing metadata macro:\n gROOT->ProcessLine(.x %s((TTree*)%p,%d);",     metadataMacro.Data(),tree, fVerbose);
    gROOT->ProcessLine(TString::Format(".x %s((TTree*)%p,%d);",metadataMacro.Data(),tree,fVerbose).Data());
  }

  return tree;
}

/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// \param period Period, e.g. 'LHC15f'
/// \param pass E.g. 'pass2' or 'passMC'
/// \param friendList - semicolomn separated array of "friend" trees attached to the tree
/// Returns the tree with the information from the corresponding resource - trees from the friendList are added as friend trees
/// see example usage
/// \return TTree* with corresponding resource

TTree*  AliExternalInfo::GetTree(TString type, TString period, TString pass, TString friendList){
  //  
  /*
    Example usage: - loading QA.TPC together with logbook and QA.EVS - visualizing data volume as fuction of interaction rate
    
    TTree * treeTPC = info.GetTree("QA.TPC","LHC15o","pass1","QA.rawTPC;QA.ITS;QA.TPC;QA.TRD;QA.TOF;QA.EVS;Logbook;Logbook.detector:TPC:detector==\"TPC\""); 
    treeTPC->Draw("Logbook.detector_TPC.bytesInjectedPhysics/Logbook.detector_TPC.eventCountPhysics/1000000:QA.EVS.interactionRate","","*")

   */
  TTree * tree = GetTree(type, period,pass);  
  if (tree==NULL) tree=  GetTree(type, period,"");
  if (tree==NULL) tree=  GetTree(type, "","");
  if (tree==NULL){
    ::Error("AliExternalInfo::GetTree","Friend tree %s not valid or empty",type.Data()); 
    return 0;
  }
  TString indexName= fConfigMap[type + ".indexname"];
  TString oldIndexName= fConfigMap[type + ".oldindexname"];
  if (oldIndexName.Length()>0 && tree->FindBranch(oldIndexName.Data())){
    tree->FindBranch(oldIndexName.Data())->SetName(indexName.Data());
  }
  if (indexName.Length()<=0) indexName="run";
  Int_t entries = tree->Draw(indexName.Data(),"","goff");
  if (entries<=0){
    ::Error("AliExternalInfo::GetTree","Friend tree %s not valid or empty",type.Data()); 
    return 0;
  }

  TObjArray * arrFriendList= friendList.Tokenize(";");
  for (Int_t ilist=0; ilist<arrFriendList->GetEntriesFast(); ilist++){
    TString fname=arrFriendList->At(ilist)->GetName();
    TString conditionName="";   
    TString condition="";
    Int_t nDots = fname.CountChar(':');
    // in case there are more than one entry for primary index - secondary key has to be specified
    // following syntax is used in this case <treeID>:conditionName:condition
    //     e.g Logbook.detector:TPC:detector==\"TPC\"
    if (nDots!=0 && nDots!=2) continue;
    if (nDots==2){
      TObjArray * tokenArray = fname.Tokenize(":");
      fname=tokenArray->At(0)->GetName();
      conditionName=tokenArray->At(1)->GetName();
      condition=tokenArray->At(2)->GetName();
      delete tokenArray;
    }
    //
    TTree *ftree= GetTree(fname.Data(), period,pass,nDots==0);  // Standard build index if not custom selection
    if (ftree==NULL || ftree->GetEntries()<=0){
      ::Error("AliExternalInfo::GetTree","Friend tree %s not valid or empty",fname.Data()); 
      continue;
    }    
    if (nDots==2){
      tree->SetAlias(conditionName.Data(),"(1+0)");    
      ftree->SetAlias(conditionName.Data(),condition.Data());
      ftree->BuildIndex(tree->GetTreeIndex()->GetMajorName(), conditionName.Data());    
      tree->AddFriend(ftree, (fname+"_"+conditionName).Data());
    }else{
      tree->AddFriend(ftree, fname.Data());
    }
    //
    ftree->SetTitle(TString::Format("%s.%s",fname.Data(),ftree->GetTitle()).Data());
    ftree->SetName(TString::Format("%s.%s",fname.Data(),ftree->GetName()).Data());
    //ftree->AddFriend(tree, type.Data());
    Int_t fentries = tree->Draw(indexName.Data(),"","goff");
    ::Info("AliExternalInfo::GetTree","AddFriend %s+%s - entries=%d", type.Data(),  fname.Data(),fentries);
  }
  return tree;
}



/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// \param period Period, e.g. 'LHC15f'. Here you can use wildcards like in 'ls', e.g. 'LHC15*'
/// \param pass E.g. 'pass2' or 'passMC'. Here you can use wildcards like in 'ls', e.g. 'pass*'
/// Returns a chain with the information from the corresponding resources.
/// \return TChain*
TChain* AliExternalInfo::GetChain(TString type, TString period, TString pass){
  // FIXME  - here we should also fix Leave name bug
  TChain* chain = 0x0;
  TString internalFilename = ""; // Resulting path to the file
  TString internalLocation = fLocalStorageDirectory; // Gets expanded in this function to the right directory
  TString externalLocation = "";
  Bool_t resourceIsTree = kFALSE;
  TString detector = "";
  TString rootFileName = "";
  TString treeName = "";
  TString pathStructure = "";
  TString indexName=""; 
 

  // Setting up all the local variables
  SetupVariables(internalFilename, internalLocation, resourceIsTree, pathStructure, detector, rootFileName, treeName, type, period, pass,indexName);

  TString cmd = TString::Format("/bin/ls %s", internalFilename.Data());
  // std::cout << "==== cmd: " << cmd << std::endl;

  TString files=gSystem->GetFromPipe(cmd.Data());
  TObjArray *arrFiles=files.Tokenize("\n");
  AliInfo(TString::Format("Files to add to chain: %s", files.Data()));

  //function to get tree namee based on type
  chain=new TChain(treeName.Data());
   // ---| loop over possible tree names |---
  TObjArray *arrTreeName = treeName.Tokenize(",");

  for (Int_t ifile=0; ifile<arrFiles->GetEntriesFast(); ++ifile) {
    for (Int_t itree=0; itree<arrTreeName->GetEntriesFast(); itree++){
      TFile *ftemp=TFile::Open(arrFiles->At(ifile)->GetName());
      if (ftemp==NULL) continue;
      if (ftemp->GetListOfKeys()->FindObject(arrTreeName->At(itree)->GetName())){
	chain->AddFile(arrFiles->At(ifile)->GetName(),0, arrTreeName->At(itree)->GetName());
	
      }
      delete ftemp;
    }
  }
  
  const TString cacheSize=fConfigMap[type + ".CacheSize"];
  Long64_t cache=cacheSize.Atoll();
  if (fMaxCacheSize>0) {
    if (cache>0) {
      cache=TMath::Min(fMaxCacheSize, cache);
    } else {
      cache=fMaxCacheSize;
    }
  }
  if (cache>0) chain->SetCacheSize(cache);

  AddChain(type, period, pass);
  delete arrFiles;
  delete arrTreeName;
  return chain;
};

/// Every tree you create is added to a big tree acting as a friend.
/// You can have access to this tree with the GetFriendsTree() function.
/// @TODO Add 'return false' when adding to the friends tree was not successful
/// \return kTRUE
/// Comment:  function should do only Build index - adding freind tree is not sufficently general
///           
Bool_t AliExternalInfo::BuildIndex(TTree* tree, TString type){
  //
  //
  //
  TString indexName= fConfigMap[type + ".indexname"];
  TString oldIndexName= fConfigMap[type + ".oldindexname"];
  TString metadataMacro=fConfigMap[type + "metadataMacro"];
  //
  if (oldIndexName.Length()>0){  // rename branch  with index if specified in configuration file
    if (tree->GetBranch(oldIndexName.Data())) {
      tree->GetBranch(oldIndexName.Data())->SetName(indexName.Data());
    }
  }
  if (indexName.Length()<=0) { // set default index name
    if (tree->GetListOfBranches()->FindObject("run"))  indexName="run";    
  }
  if (indexName.Length()<=0) {
    ::Error("AliExternalInfo::BuildIndex","Index %s not avaible for type %s", indexName.Data(), type.Data());
  }  
  if (tree->GetBranch(indexName.Data()) && TString(tree->GetBranch(indexName.Data())->GetTitle()).Contains("/C")){
    BuildHashIndex(tree,indexName.Data(),"hashIndex");
  }else{
    tree->BuildIndex(indexName.Data());
  }
  TStatToolkit::AddMetadata(tree,"TTree.indexName",indexName.Data());

 //  TString name = "";

//   if(type.Contains("QA")){ // use TPC instead of QA.TPC
//     name = type(3, type.Length()-1);
//   }
//   else if(type.Contains("MonALISA")){
//     name = type(9, type.Length()-1);
//   }
//   else {
//     name = type;
//   }
//   tree->SetName(name);
//   if (fTree == 0x0) fTree = dynamic_cast<TTree*>(tree->Clone());

 
//   fTree->AddFriend(tree, name);  
//  AliInfo(TString::Format("Added as friend with the name: %s",name.Data()));
  ::Info("AliExternalInfo::BuildIndex", "TreeName:%s;IndexName:%s",tree->GetName(), indexName.Data());
  return kTRUE;
}
/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// \param period Period, e.g. 'LHC15f'. Here you can use wildcards like in 'ls', e.g. 'LHC15*'
/// \param pass E.g. 'pass2' or 'passMC'. Here you can use wildcards like in 'ls', e.g. 'pass*'
/// Not fully working. Should automatically add every downloaded file to a huge chain.
Bool_t AliExternalInfo::AddChain(TString type, TString period, TString pass){
  // Adds chain of trees and buildsindexs and addfriends it
  TString internalFilename = ""; // Resulting path to the file
  TString internalLocation = fLocalStorageDirectory; // Gets expanded in this function to the right directory
  TString externalLocation = "";
  Bool_t resourceIsTree = kFALSE;
  TString detector = "";
  TString rootFileName = "";
  TString treeName = "";
  TString pathStructure = "";
  TString indexName=""; 

  // Setting up all the local variables
  SetupVariables(internalFilename, internalLocation, resourceIsTree, pathStructure, detector, rootFileName, treeName, type, period, pass,indexName);
  AliInfo(TString::Format("Add to internal Chain: %s", internalFilename.Data()));
  AliInfo(TString::Format("with tree name: %s",        treeName.Data()));
  fChain->AddFile(internalFilename.Data(), TChain::kBigNumber, treeName);

  return kTRUE;
}

/// \param period Period in the form eg LHC15f
/// Generate from eg LHC15f the year 2015
/// \return Year in the form "2015"
const TString AliExternalInfo::GetYearFromPeriod(const TString &period){
  TString year(period(3,2));
  year = TString::Format("20%s", year.Data());
  return year;
}

/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// \param period Period, e.g. 'LHC15f'. Here you can use wildcards like in 'ls', e.g. 'LHC15*'
/// \param pass E.g. 'pass2' or 'passMC'. Here you can use wildcards like in 'ls', e.g. 'pass*'
/// Returns a TString containing the complete directory structure where the root
/// file should be stored/loaded from. eg './data/2015/LHC15f/pass2/'
const TString AliExternalInfo::CreatePath(TString type, TString period, TString pass){
  // Create the local path from the type, period and pass of the resource
  TString internalLocation;
  //Check if period is MC and adjust storage hierarchy
  if (period.Length() == 6 || (period == "" && type != "MonALISA.MC") || type == "MonALISA.ProductionCycleID"  || type == "MonALISA.RCT"|| type == "TriggerClasses") { // put everything in the form LHCYYx or with empty period and pass in "data"
    internalLocation.Append("/data/");
  }
  else { // everything which is not in the form "LHCYYx" put in /sim/
    internalLocation.Append("/sim/");
  }
  // std::cout << "Internal save path 1: " << internalLocation << std::endl;
  if (period != "" && type != "MonALISA.ProductionCycleID" && type != "QA.rawTPC"){
    const TString year = GetYearFromPeriod(period);
    internalLocation += year + "/";
    // std::cout << "Internal save path 2: " << internalLocation << std::endl;
    internalLocation += period + "/";
    // std::cout << "Internal save path 3: " << internalLocation << std::endl;
    if (pass != ""){
      internalLocation += pass + "/";
      // std::cout << "Internal save path 4: " << internalLocation << std::endl;
    }
  }
  else if (type == "QA.rawTPC"){ // only needs the year
    const TString year = GetYearFromPeriod(period);
    internalLocation += year + "/";
  }
  else if (type == "MonALISA.ProductionCycleID") { // Needed for ProductionCycleIDs
    internalLocation += "ID_" + period + "/";
  }
  return internalLocation;
}

/// \param file Exact location of the file which should be checked
/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// Checks if the file is older than the timeout which is specified in the config file of this
/// specific resource type
/// \returns true if download of the resource is needed and false if not
Bool_t AliExternalInfo::IsDownloadNeeded(TString file, TString type){
  Int_t timeOut = 0;
  timeOut = atoi(fConfigMap[type + ".timeout"]);

  // std::cout << "-- Check, if " << file << " is already there" << std::endl;
  if (gSystem->AccessPathName(file.Data()) == kTRUE) {
    AliInfo("-- File not found locally --> Caching from remote");
    return kTRUE;
  }
  else {
    // std::cout << "---- File already downloaded --> Check if older than timelimit" << std::endl;
    struct stat st;
    stat(file.Data(), &st);
    std::time_t timeNow = std::time(0);
    long int timeFileModified = st.st_mtime;
    // std::cout << "------ File is " << timeNow - timeFileModified << " seconds old" << std::endl;
    if (timeNow - timeFileModified < timeOut) {
      AliInfo(TString::Format("-- File is %li s old; NOT older than the set timelimit %d s",timeNow - timeFileModified, timeOut));
      return kFALSE; // if file is younger than the set time limit, it will not be downloaded again
    }
    else {
      AliInfo(TString::Format("-- File is %li s old; Older than the set timelimit %d s",timeNow - timeFileModified, timeOut));
      return kTRUE;
    }
  }
}
/// \param mifFilePath Location of the mif-file which is downloaded from Monalisa; Is changed by this function
/// \param internalLocation Directory where the root file is stored
/// \param rootFileName Location of the newly created root file
/// \param externalLocation Location specified in the config file
/// Composes the curl-command in a TString which afterwards then can be executed
/// \return curl-command in a TString
const TString AliExternalInfo::Curl(TString& mifFilePath, const TString& internalLocation, TString rootFileName, const TString& externalLocation){
  TString command = "";
  TString certificate("$HOME/.globus/usercert.pem");
  TString privateKey("$HOME/.globus/userkey.pem");

  // TString internalLocation = internalFilename(0, internalFilename.Last('/') + 1); // path to internal location
  // Create path to logfile with same name as the root file
  TString logFileName = rootFileName.ReplaceAll(".root", ".log");
  logFileName.Prepend(internalLocation);

  mifFilePath = rootFileName.ReplaceAll(".log", ".mif");
  mifFilePath.Prepend(internalLocation);

  command = TString::Format("curl -z %s -k --tlsv1 --cert %s --key %s -o %s 2>%s \"%s\"",
                                     mifFilePath.Data(), certificate.Data(), privateKey.Data(),
                                     mifFilePath.Data(), logFileName.Data(), externalLocation.Data());

  return command;
}


/// \param period  LHC period
/// \param pass    calibratio pass

TTree * AliExternalInfo::GetCPassTree(const char * period, const  char *pass){
  //
  // Try to find production information about pass OCDB export
  // To find the production description field of the overlal production table is queried
  //
  // Warnig:
  //    In some cases mif format internaly used not stable 
  //    Unit consistency test should be part of procedure
  //
  TTree * treeProdArray=0, *treeProd=0;
  AliExternalInfo info;
  treeProdArray = GetTreeCPass();
  treeProdArray->Scan("ID:Description:Tag",TString::Format("strstr(Tag,\"%s\")&&strstr(Tag,\"%s\")&& strstr(Description,\"merging\")",period,pass).Data(),"col=10:100:100");
  // check all candidata production and select one which exports OCDB
  Int_t entries= treeProdArray->Draw("ID:Description",TString::Format("strstr(Description,\"%s\")&&strstr(Description,\"%s\")&& strstr(Description,\"merging\")",period,pass).Data(),"goff");  
  for (Int_t ientry=0; ientry<entries; ientry++){
    TTree * treeProd0 = info.GetTreeProdCycleByID(TString::Format("%.0f",treeProdArray->GetV1()[ientry]).Data());
    Int_t status = treeProd0->Draw("1","strstr(outputdir,\"OCDB\")==1","goff");       // check presence of the OCDB 
    if (status==0) {
      delete treeProd0;
      continue;
    }
    treeProd=treeProd0;
  }
  return treeProd;
}

/// Add branch index branch with name chindexName calculated from input string branch chbranchName
/// and make thi branch as a primary index 
/// used in order to qury FreindTrees with string keyname (impossible with standard ROOT)
/// 
/// \param tree            - input tree
/// \param chbranchName    - branch name with index
/// \param chindexName     - name if the index branch
void AliExternalInfo::BuildHashIndex(TTree* tree, const char *chbranchName,  const char *chindexName){
  //
  //
  Int_t indexName=0;
  char  pbranchName[100];
  TBranch *brIndexMC = tree->Branch(chindexName,&indexName,TString::Format("%s/I",chindexName).Data()); // branch to fill
  TBranch *branch=tree->GetBranch(chbranchName); // branhc to get string
  if (branch!=NULL){
    branch->SetAddress(&pbranchName);
  }else{
    //;
  }
  Int_t entries= tree->GetEntries();
  for (Int_t ientry=0; ientry<entries; ientry++){
    branch->GetEntry(ientry);
    indexName=TMath::Hash(pbranchName);
    brIndexMC->Fill();
  }
  tree->BuildIndex(chindexName);
  brIndexMC->SetAddress(0);  // reset pointers to the branches to 0
  branch->SetAddress(0);     // reset pointers to the branches to 0
}

void   AliExternalInfo::PrintConfigSelected(const char *expName, const char *expValue){
  //
  //
  PrintMapSelected(fConfigMap, expName,expValue);
}


void AliExternalInfo::PrintMapSelected(std::map<TString, TString> infoMap, const char *expName, const char *expValue){
  //
  // This function to move to some util class
  // Simple printing of the configuration filtering using Reg. expression
  //     algorithm consuted with the http://stackoverflow.com/questions/17253690/finding-in-a-std-map-using-regex
  //     No support in std::map - linear loop also sugested (faster version have to make some assumptions)
  //
  // More complex filter with logical operation will be written later  
  //     expName   - regular expression to selec in tag name
  //     expValue   - regular expression to selec in tag value
  //
  TPRegexp regExpName(expName);
  TPRegexp regExpValue(expValue);
  std::map<TString, TString>::iterator it;    
  for (it=infoMap.begin(); it!=infoMap.end(); ++it){
    Bool_t isSelected=kTRUE;
    if (regExpName.GetPattern().Length()>0) isSelected&= regExpName.Match(it->first);
    if (regExpValue.GetPattern().Length()>0) isSelected&= regExpValue.Match(it->second);
    if (isSelected) printf("%s\t%s\n", it->first.Data(), it->second.Data());
  }
  
}


TTree*  AliExternalInfo::GetTreeAliVersRD(){            
    
//    returns and stores tree ("dumptree") containing the relevant information of real data productions
//    for guessing the anchor pass of MC productions
   TTree * treeProd = GetTreeProdCycle();              // getting tree with information on real data productions (list, id, tag) - id will be used to get info for each production via GetTreeProdCycleByID(TString::Format("%d",id))
   TFile* outfile;
   

   Bool_t downloadNeeded = IsDownloadNeeded(fLocalStorageDirectory+TString::Format("/dumptree_RD.root"),TString::Format("QA.TPC"));     //check if download is needed
   if(!downloadNeeded){
        outfile = TFile::Open(fLocalStorageDirectory+TString::Format("%s","/dumptree_RD.root"),"UPDATE");
        if(outfile->GetListOfKeys()->Contains("dumptree_RD")){
           AliInfo("-- dumptree_RD.root found locally and validated--> Not caching");
           return((TTree*)outfile->Get("dumptree_RD"));
        } 
    } 
        
    else  AliInfo("-- dumptree_RD.root not validated--> Caching from remote");
   
    outfile= new TFile(fLocalStorageDirectory+"/dumptree_RD.root","RECREATE");
    TTree *dumptree = treeProd->CloneTree();       //tree that will hold information for guessing
    
    Int_t id=0;
    char tag[1000];
    
    dumptree->SetBranchAddress("ID",&id);
    dumptree->SetBranchAddress("Tag",&tag);
    
    char paliroot[1000];            // variables that will hold information read of the tree: GetTreeProdCycleByID(TString::Format("%d",id))
    char paliphysics[1000];
    char poutputdir[1000];

    TObjString sprodname;           // variables that will be written into dumptree 
    TObjString spassname;
    TObjString soutputpath;
    TObjString saliroot;
    TObjString saliphysics;
    Bool_t consist;
    
    TString soutputdir;
    TObjArray *subStrL;

    TBranch* braliroot= dumptree->Branch("aliroot",&saliroot);
    TBranch* braliphys= dumptree->Branch("aliphysics",&saliphysics);
    TBranch* brprodname= dumptree->Branch("prodName",&sprodname);
    TBranch* brpassname= dumptree->Branch("passName",&spassname);
    TBranch* broutputpath= dumptree->Branch("outputPath",&soutputpath);
    TBranch* brconsist= dumptree->Branch("nameconsistency",&consist);

    Int_t entries=dumptree->GetEntries();
    for (Int_t i=0; i<entries; i++){            //loop over all IDs
      dumptree->GetEntry(i);
      AliInfo(TString::Format("Getting ProdCyle ID: %d",id));
      TTree * tree= GetTreeProdCycleByID(TString::Format("%d",id));         //get tree with production info for each ID

      if (tree==NULL) cout<<"err0"<<endl;
      if (tree->GetBranch("app_aliphysics")==NULL) continue;
      tree->SetBranchAddress("app_aliphysics",&paliphysics);      // set prod info branch addresses
      tree->SetBranchAddress("app_aliroot",&paliroot);
      tree->SetBranchAddress("outputdir",&poutputdir);

      tree->GetEntry(0);                                            // read prod info tree entry
      
      soutputdir= TString::Format("%s",poutputdir);                 // extract production name from outputdir

      subStrL = TPRegexp("(?=LHC)(.*?)(?=/)").MatchS(soutputdir);
      sprodname = *((TObjString *)subStrL->At(0)); 
      delete subStrL;
      subStrL = TPRegexp("[A-Za-z0-9]*").MatchS(sprodname.String());
      sprodname = *((TObjString *)subStrL->At(0));
      delete subStrL;

      subStrL = TPRegexp("[^/]+$").MatchS(soutputdir);              //extract pass name from utputdir
      spassname = *((TObjString *)subStrL->At(0));
      delete subStrL;      

      saliroot = TObjString(paliroot);
      saliphysics = TObjString(paliphysics);
      soutputpath= TObjString(poutputdir);
      
      if(TString::Format("%s",tag).Contains(sprodname.GetString())) consist=kTRUE;
      else consist = kFALSE;
      
      braliroot->Fill();
      braliphys->Fill();
      brprodname->Fill(); 
      brpassname->Fill(); 
      broutputpath->Fill();
      brconsist->Fill();
      
      delete tree;
    }
    outfile->cd();
    dumptree->Write("dumptree_RD");

    delete outfile;
    return dumptree; 
  }
  
  
TTree*  AliExternalInfo::GetTreeAliVersMC(){
//    returns and stores tree ("dumptree") containing the relevant information of MC productions
//    for guessing the anchor pass of MC productions
   TTree* treeMC = GetTreeMC(); 
   TTree* treeProdMC = GetTree("MonALISA.ProductionMC","","");          //tree with "Description" branch - needed to determine if the MC production is general purpose production
   treeMC->AddFriend(treeProdMC);
   TFile* outfile;

   Bool_t downloadNeeded = IsDownloadNeeded(fLocalStorageDirectory+TString::Format("/dumptree_MC.root"),TString::Format("QA.TPC"));     //check if download is needed
   if(!downloadNeeded){
        outfile = TFile::Open(fLocalStorageDirectory+TString::Format("%s","/dumptree_MC.root"),"UPDATE");
        if(outfile->GetListOfKeys()->Contains("dumptree_MC")){
           AliInfo("-- dumptree_MC.root found locally and validated--> Not caching");
           return((TTree*)outfile->Get("dumptree_MC"));
        } 
    } 
    
   AliInfo("-- dumptree_MC.root not validated--> Caching from remote");
   
   outfile= new TFile(fLocalStorageDirectory+"/dumptree_MC.root","RECREATE");
   TTree* dumptree=treeMC->CloneTree();
   
   //variable to read tree from GetTreeMC()
   char panchprodname[1000];
   char prunlist[50000];
   char pdescr[50000];
   
   dumptree->SetBranchAddress("anchorProdTag",&panchprodname);
   dumptree->SetBranchAddress("runList",&prunlist);
   dumptree->SetBranchAddress("Description",&pdescr);

   TObjString sMCanchprodname; 
   TObjString sMCdescr;
   Int_t first=-1;
   Int_t last=-1;

   TString sfirst;
   TString slast;
   TString sanprod;
   TObjArray *subStrL;
   
   TBranch* brMCanchprodname= dumptree->Branch("anchorProdTag_ForGuess",&sMCanchprodname);
   TBranch* brfirst= dumptree->Branch("First_Run",&first);
   TBranch* brlast= dumptree->Branch("Last_Run",&last);
   TBranch* brMCdescr= dumptree->Branch("Description",&sMCdescr);
   
   Int_t entries=dumptree->GetEntries();  
   for (Int_t i=0; i<entries; i++){             //loop overall MC production
     dumptree->GetEntry(i);                       //read info

     sanprod= TString::Format("%s",panchprodname);          //extract anchor production name from anchorProdTag
     subStrL = TPRegexp("[A-Za-z0-9]*").MatchS(sanprod);
     sanprod = ((TObjString *)subStrL->At(0))->GetString(); 
     delete subStrL; 

     if(TString::Format("%s",prunlist).Length()!=0){        //extract first run number
        subStrL = TPRegexp("^[^ ,]+").MatchS(TString::Format("%s",prunlist));
        sfirst = ((TObjString *)subStrL->At(0))->GetString();
        delete subStrL;
     }
     else sfirst=TString("-1");

     if(TString::Format("%s",prunlist).Length()!=0){        //extract last run number
        subStrL = TPRegexp("[^ ,]+$").MatchS(TString::Format("%s",prunlist));
        slast = ((TObjString *)subStrL->At(0))->GetString();
        delete subStrL;
     }
     else slast=TString("-1");
     
     sMCanchprodname = TObjString(sanprod);
     sMCdescr = TObjString(pdescr);
     first=sfirst.Atoi();
     last=slast.Atoi();
          
     brMCanchprodname->Fill();
     brfirst->Fill();
     brlast->Fill();
     brMCdescr->Fill();
     
   }
   outfile->cd();
   dumptree->Write("dumptree_MC");

   delete outfile;
   delete treeMC;
   delete treeProdMC;
   return dumptree;
}



class anchprod{             //class that holds information about aliphysics and aliroot version, the pass name and the production name - will be used to sort the real data and MC productions
public:
TString aliphys;
TString aliroot;
TString anchpass;
TString anchprodname;

anchprod(){
aliphys="-1";
aliroot="-1";
anchpass="-1";
anchprodname="-1";
}
bool operator< (const anchprod & otheranchprod) const       //define comparison operator for lexicographic ordering
	{
		return (aliphys+aliroot < otheranchprod.aliphys+otheranchprod.aliroot);         //compare first aliphysics, if same (i.e. empty) then aliroot, duplicates in ali versions would be stored only once also if passes diferent
	}
};


TTree*  AliExternalInfo::GetTreeMCPassGuess(){
//    returns and stores (dumptree_MC.root) the tree containing the pass guesses for each MC production
    
    TFile *MCFile;
    TTree *MCTree;
        
    TFile *RDFile;
    TTree *RDTree;
    
    Bool_t downloadNeeded = IsDownloadNeeded(fLocalStorageDirectory+TString::Format("/dumptree_RD.root"),TString::Format("QA.TPC"));
    if(!downloadNeeded){
        AliInfo("found dumptree_RD.root -> getting dumptree_RD and make guesses");
        RDFile= TFile::Open(fLocalStorageDirectory+"/dumptree_RD.root","UPDATE");
        RDTree = dynamic_cast<TTree*>(RDFile->Get("dumptree_RD"));
    }
    else{
        AliInfo("dumptree_RD.root not available -> run GetTreeAliVersRD()...");
        GetTreeAliVersRD();
        AliInfo("Got tree from GetTreeAliVersRD()");
        RDFile= TFile::Open(fLocalStorageDirectory+"/dumptree_RD.root","UPDATE");
        RDTree = dynamic_cast<TTree*>(RDFile->Get("dumptree_RD"));
    }

    TObjString* osrdprod=0;             //variables for reading RD tree
    TObjString* osrdpass=0;
    TObjString* osrdaliphys=0;
    TObjString* osrdaliroot=0;

    RDTree->GetBranch("aliphysics")->SetAddress(&osrdaliphys);
    RDTree->GetBranch("aliroot")->SetAddress(&osrdaliroot);
    RDTree->GetBranch("prodName")->SetAddress(&osrdprod);
    RDTree->GetBranch("passName")->SetAddress(&osrdpass); 
   
    downloadNeeded = IsDownloadNeeded(fLocalStorageDirectory+TString::Format("/dumptree_MC.root"),TString::Format("QA.TPC"));
    if(!downloadNeeded){
        MCFile= TFile::Open(fLocalStorageDirectory+"/dumptree_MC.root","UPDATE");
        if(MCFile->GetListOfKeys()->Contains("dumptree_MC_guess")){
        AliInfo("Guesses already available - done.");
        return(dynamic_cast<TTree*>(MCFile->Get("dumptree_MC_guess")));
        } 
        else if(MCFile->GetListOfKeys()->Contains("dumptree_MC")) {
        AliInfo("found dumptree_MC.root -> getting dumptree_MC and make guesses");
        MCTree = dynamic_cast<TTree*>(MCFile->Get("dumptree_MC"));
        }
    }
    else{
        AliInfo("dumptree_MC.root not available -> run GetTreeAliVersMC()...");
        GetTreeAliVersMC();
        AliInfo("Got tree from GetTreeAliVersMC()");
        MCFile= TFile::Open(fLocalStorageDirectory+"/dumptree_MC.root","UPDATE");
        MCTree = dynamic_cast<TTree*>(MCFile->Get("dumptree_MC"));
    }

    MCFile->cd();

    TObjString osMCaliroot;          //char arrays for reading from MCTree
    TObjString osMCaliphysics;
    TObjString osMCprodname;
    TObjString* osMCanchprodname=0; //variables for reading TObjString from MCTree
    TObjString* osMCdescr=0;
    Bool_t isgp = kFALSE;         //is MC production a general purpose production?
        
    char pMCaliroot[1000];     //variable to read tree from GetTreeMC()
    char pMCaliphysics[1000];
    char pMCprodname[1000];
     
    MCTree->GetBranch("aliphysics")->SetAddress(&pMCaliphysics);       //set branch addresses
    MCTree->GetBranch("aliroot")->SetAddress(&pMCaliroot);
    MCTree->GetBranch("prodName")->SetAddress(&pMCprodname);
    MCTree->GetBranch("anchorProdTag_ForGuess")->SetAddress(&osMCanchprodname);
    MCTree->GetBranch("Description")->SetAddress(&osMCdescr);
    
    TBranch *branchpass = MCTree->Branch("anchorPassName_guess",&osrdpass);         //adding new branches holding information about pass and ali-versions of guessed RD production
    TBranch *branchroot = MCTree->Branch("anchoraliroot_guess",&osrdaliroot);
    TBranch *branchphys = MCTree->Branch("anchoraliphys_guess",&osrdaliphys);

    int n = RDTree->GetEntries();
    int m = MCTree->GetEntries();
    Bool_t prfound= kFALSE; 

    anchprod tempprod;

    multiset<anchprod> list;                                                        //set that will hold "anchprod" instances of RD productions, that have have matching production names
    multiset<anchprod>::iterator it = list.begin();

    for (Int_t i=0; i<m; i++) {     //loop over MC productions

        MCTree->GetEntry(i);
        
        isgp = TPRegexp("General").MatchB(osMCdescr->String(),"i") && TPRegexp("Purpose").MatchB(osMCdescr->String(),"i");
        
        osMCaliroot = TObjString(pMCaliroot);          //get TObjStrings for guessing
        osMCaliphysics= TObjString(pMCaliphysics);
        osMCprodname= TObjString(pMCprodname);
    
        cout<<endl;
        cout<<i<<" of "<<m<<endl;
        cout<<"MC Production name: "<<osMCprodname.String()<<" Anchor Production name: "<<osMCanchprodname->String()<<" MC aliphys: "<<osMCaliphysics.String()<<" MC aliroot: "<<osMCaliroot.String()<<" MC description: "<<osMCdescr->String()<<" isgp: "<<isgp<<endl;

        prfound=kFALSE;         // flag to know if any RD production with matching production name was found
        list.clear();                  //reset set of RDinfos

        for (Int_t j=0; j<n; j++) {        //search for matching RD production and get pass info

            RDTree->GetEntry(j);

            if(osMCanchprodname->String()==""){
                cout<<"Complete match: no. No anchorprod found -> skip"<<endl;
                break;}

            if(osMCanchprodname->String()==osrdprod->String()){           //if RDprodname is right then save the prod infos to be later eventually able to look up what was closest in terms of aliphys/aliroot
                prfound=kTRUE;                                            //found a prodction with correct production name

                tempprod.aliphys=osrdaliphys->GetString();                //store info in list of "anchprod"
                tempprod.aliroot=osrdaliroot->GetString();
                tempprod.anchpass=osrdpass->GetString();
                tempprod.anchprodname=osrdprod->GetString();
                list.insert(tempprod);

            }            
            if(osMCanchprodname->String()==osrdprod->String() && (( osMCaliphysics.String()!="" && osMCaliphysics.String()==osrdaliphys->String())  || (osMCaliphysics.String()=="" && osrdaliroot->String()==osMCaliroot.String())  ) && !osrdpass->GetString().Contains("cpass")  && !osrdpass->GetString().Contains("cosmic") && (!isgp || TPRegexp("^pass").MatchB(osrdpass->GetString(),"i"))){      //check for perfect match, i.e. prodname and aliphys/aliroot match - exclude RDpasses that contain "comsic" or "cpass" and only allow pass names starting with "pass" when MC prod is general purpose
         
                 cout<<"Complete match: yes"<<endl;
                 cout<<"Used for guess: RDphys:"<<osrdaliphys->GetString()<<" RDroot: "<<osrdaliroot->GetString()<<" RDpass guess: "<<osrdpass->GetString()<<endl;
                 prfound=kTRUE;
                 break;             //break if perfect match found      ?? wanted ??
            }

            if(j==n-1 && prfound){            //if match not found check what was closest if matching prodname was found
                             tempprod.aliphys=osMCaliphysics.String();         //make anchprd instance with MC info and insert into its lexicographical position
                             tempprod.aliroot=osMCaliroot.String();
                             tempprod.anchprodname=osMCprodname.String();
                             tempprod.anchpass=TString("MCPass");               //dummy info to make visible in list what was MC entry
                             list.insert(tempprod);     //insert MC prod 
                             it=list.find(tempprod);         //get iterator to pointer before MC prod
                             if (it == list.begin()){
                                 prfound =kFALSE;       //no interesting ones found
                                 break;
                             }

            int l =0;
            cout<<"Complete match: no -> take closest aliversion as guess from RD productions: "<<endl;
            for (multiset<anchprod>::iterator iter=list.begin(); iter!=list.end(); ++iter){         //cout the ordered list of productions
              l++;
              cout<<l<<" RDprodname: "<<(*iter).anchprodname<<" AliPhys: "<<(*iter).aliphys<<" AliRoot: "<<(*iter).aliroot<<" RDpass: "<<(*iter).anchpass<<endl;
            }
            
            --it;                                                       //let iterator point to entry right before MC entry that is neither a "cpass" nor a "cosmics" pass
            for (multiset<anchprod>::iterator iter=it; ; --iter){       //go backwards through list and take first pass guess that is not a cpass
                cout<<"Looking for MC aliphys: "<<tempprod.aliphys<<" MC aliroot: "<<tempprod.aliroot<<endl;
                if(!((*iter).anchpass).Contains("cpass") && !((*iter).anchpass).Contains("cosmic") && (!isgp || TPRegexp("^pass").MatchB((*iter).anchpass,"i"))){
                cout<<"Used for guess: RDphys:"<<(*iter).aliphys<<" RDroot: "<<(*iter).aliroot<<" RDpass guess: "<<(*iter).anchpass<<endl;
                break;
                    }
                if(iter==list.begin()) break;
                }

            }
        }
        branchpass->Fill();
        branchphys->Fill();
        branchroot->Fill();		
        }

    MCTree->Write("dumptree_MC_guess");
    delete MCFile;    
    return(MCTree);
    
}
    
     
TString  AliExternalInfo::GetMCPassGuess(TString sMCprodname){
    
//returns string with Pass guess for a given MC production name    
 TFile* MCFile;
 TTree* guesstree; 
    
 Bool_t downloadNeeded = IsDownloadNeeded(fLocalStorageDirectory+TString::Format("/dumptree_MC.root"),TString::Format("QA.TPC"));   //check if file containing guesses is present
 if(!downloadNeeded) {
     AliInfo("dumptree_MC.root available");
     MCFile= TFile::Open(fLocalStorageDirectory+"/dumptree_MC.root");
     if( !(MCFile->GetListOfKeys()->Contains("dumptree_MC_guess"))){
        AliInfo("dumptree_MC_guess not available - get it vai GetTreeMCPassGuess()");
        guesstree = GetTreeMCPassGuess();
     }
     else{
         AliInfo("dumptree_MC_guess available");
         guesstree = dynamic_cast<TTree*>(MCFile->Get("dumptree_MC_guess"));
     }
 }
 
 else{
     AliInfo("dumptree_MC_guess not available ");
     GetTreeMCPassGuess();
     MCFile= TFile::Open(fLocalStorageDirectory+"/dumptree_MC.root");
     guesstree = dynamic_cast<TTree*>(MCFile->Get("dumptree_MC_guess"));
     AliInfo("got the tree with guesses");
}
 
 char pMCprodname[1000];
 TObjString osMCprodname=0;
 TObjString* osMCpassguess=0;
 
 guesstree->GetBranch("prodName")->SetAddress(&pMCprodname);
 guesstree->GetBranch("anchorPassName_guess")->SetAddress(&osMCpassguess);
 
 for(int i=0;i<guesstree->GetEntries();i++){        //loop over tree with guesses
     guesstree->GetEntry(i);
     osMCprodname = TObjString(pMCprodname);
     if(osMCprodname.String()==sMCprodname){       //if match found return corresponding guess
         cout<<"Anchor Pass guess for "<<osMCprodname.String()<<": "<<osMCpassguess->String()<<endl;
         return(osMCpassguess->String());
     }
 }
 cout<<osMCprodname.String()<<" was not found in list of MC productions"<<endl;
 return(TString::Format("MC production not found"));
}
