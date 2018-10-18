#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <string>
#include <map>
#include <limits>
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
#include "TError.h"
ClassImp(AliExternalInfo)

using namespace std;

const TString AliExternalInfo::fgkDefaultConfig="$ALICE_ROOT/STAT/Macros/AliExternalInfo.cfg";

AliExternalInfo::AliExternalInfo(TString localStorageDirectory, TString configLocation, Int_t verbose/*, Bool_t copyToLocalStorage*/) :
                                /*fCopyDataToLocalStorage(copyToLocalStorage),*/
                                TObject(),
				                        fVerbose(verbose),
                                fLoadMetadata(kTRUE),
                                fConfigLocation(configLocation),
                                fLocalStorageDirectory(localStorageDirectory),
                                fConfigMap(),
                                fTree(0x0),
                                fChain(new TChain()),
                                fChainMap(),
                                fMaxCacheSize(-1),
                                fLogCache(0x0)
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
      ::Error("AliExternalInfo::ReadConfig", "Could not find config file '%s'", configFileName.Data());
      const TString defaultConfigFileName=gSystem->ExpandPathName(fgkDefaultConfig);
      if (defaultConfigFileName!=configFileName) {
	      ::Error("AliExternalInfo::ReadConfig", "Using default config file instead");
	      configFileName=defaultConfigFileName;
      }
    }
    
    std::ifstream configFile(configFileName);
    if (!configFile.is_open()) {
      ::Error("AliExternalInfo::ReadConfig", "Could not open config file '%s'", configFileName.Data());
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
  ::Info("AliExternalInfo::PrintConfig", "User defined resources are\n");
  // looping over map with const_iterator
  typedef std::map<TString,TString>::const_iterator it_type;
  for(it_type iterator = fConfigMap.begin(); iterator != fConfigMap.end(); ++iterator) {
    std::cout << iterator->first << " " << iterator->second << "\n";
  }
  return;
}


///  AliExternalInfo::GetProductionTree(TString period, TString pass)
/// \param period    - period ID
/// \param pass      - pass ID
/// \return production tree ()
/// * Input data source MonALISA
///   * Production information from MonALISA web interface - querying tags
///     * PPass - https://alimonitor.cern.ch/prod/?t=1&res_path=mif
///     * CPass - https://alimonitor.cern.ch/prod/?t=2&res_path=mif
///     * MC - https://alimonitor.cern.ch/prod/?t=3&res_path=mif
///   * Tags query not fully reliable as syntax was chenging several times
///   * Only reliable information - path column in reulting tree (path to the output data)
/// * Algorithm:
///   * loop over all possible production (Prod, ProdCPassm, ProdMC)  for given period
///   * loop over all passes for given production
///   * check presence of output path
///   * save tree in cache file
/*!
  Example usage:
  \code
  AliExternalInfo info;
  TTree *  prodTree = info.GetProductionTree("LHC17f","pass1");
  prodTree->Scan("RunNo:outputdir:jobs_error:jobs_total","","col=10:50:10:10");
  TTree *  mcProdTree = info.GetProductionTree("LHC17k2","");
  mcProdTree->Scan("RunNo:outputdir:jobs_error:jobs_total","","col=10:50:10:10");
  \endcode
*/
TTree *  AliExternalInfo::GetProductionTree(TString period, TString pass){
  TTree *productionTree=NULL;
  TPRegexp regexpDir(TString::Format("/%s$",pass.Data()));
  for (Int_t productionSource=0; productionSource<3; productionSource++){
    TTree *  treeProductionAll = NULL;
    if (productionSource==0) treeProductionAll = GetTree("MonALISA.ProductionCycle","","");
    if (productionSource==1) treeProductionAll = GetTree("MonALISA.ProductionCPass","","");
    if (productionSource==2) treeProductionAll = GetTree("MonALISA.ProductionMC","","");
    treeProductionAll->SetBranchStatus("*",kFALSE);
    treeProductionAll->SetBranchStatus("Tag",kTRUE);
    treeProductionAll->SetBranchStatus("ID",kTRUE);
    Int_t entriesAll=treeProductionAll->GetEntries();
    for (Int_t i=0; i<entriesAll; i++){
      treeProductionAll->GetEntry(i);
      TString tag=((char*)treeProductionAll->GetLeaf("Tag")->GetValuePointer());
      if (tag.Contains(period.Data())==0) continue;
      Int_t id = TMath::Nint(treeProductionAll->GetLeaf("ID")->GetValue());
      if (fVerbose&0x2) ::Info("GetProductionTree","Check id %d\t%s",id, tag.Data());
      TTree * cTree = GetTreeProdCycleByID(TString::Format("%d", id));
      if (cTree==NULL) continue;
      if (productionSource==2) {  // do not check dir for the MC production  - different naming convention
        productionTree=cTree;
        break;
      }
      cTree->GetEntry(0);
      TString dir=((char*)cTree->GetLeaf("outputdir")->GetValuePointer());
      if (regexpDir.Match(dir)){
         productionTree=cTree;
         break;
      }
    }
    delete treeProductionAll;
    if (productionTree) { // save file in predefined folder
      TString outputPath = CreatePath("MonALISA.Production",period,pass);
      if (gSystem->Getenv("AliExternalInfoCache")){
        outputPath.Prepend(gSystem->Getenv("AliExternalInfoCache"));
      }else{
        outputPath.Prepend(".");
      }
      gSystem->mkdir(outputPath,kTRUE);
      outputPath+=fConfigMap["MonALISA.Production.filename"];
      TFile *f = TFile::Open(outputPath.Data(),"recreate");
      TTree *copyTree = productionTree->CopyTree("");
      copyTree->Write(fConfigMap["MonALISA.Production.treename"]);
      f->Close();
      if (fVerbose&0x1) ::Info("GetProductionTree","Make cache path %s",outputPath.Data());
      return productionTree;
    }
  }
  return NULL;
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
  ::Info("AliExternalInfo::SetupVariables", "Information will be stored/retrieved in/from %s", internalLocation.Data());

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
  ::Info("AliExternalInfo::Cache", "Caching of %s %s from %s in start path %s", period.Data(), pass.Data(), type.Data(), fLocalStorageDirectory.Data());

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

  TString mifFilePath = ""; // Gets changed in Curl command

  if (downloadNeeded == kTRUE){
    if (resourceIsTree == kTRUE && externalLocation.Contains("http")) {
      externalLocation += pathStructure + rootFileName;
      Int_t fstatus=0;
      TString command = CurlTree(internalFilename, externalLocation);
      std::cout << command << std::endl;
      gSystem->Exec(command.Data());
      TFile * fcache = TFile::Open(internalFilename);
      if (fcache!=NULL && !fcache->IsZombie()) {
        fstatus|=1;
        delete fcache;
      }
      if (fstatus==1) {
        gSystem->GetFromPipe(Form("touch %s",internalFilename.Data()));  // Update the access and modification times of each FILE to the current time
        return kTRUE;
      }else{
        ::Error("AliExternalInfo::Cache", "Curl caching failed");
        gSystem->GetFromPipe(Form("rm %s",internalFilename.Data()));
        return kFALSE;
      }
    }
    // Download resources in the form of .root files in a tree
    if (resourceIsTree == kTRUE ) {
      externalLocation += pathStructure + rootFileName;
      ::Info("AliExternalInfo::Cache", "Information retrieved from: %s", externalLocation.Data());
      // Check if external location is a http address or locally accessible
      //    std::cout << externalLocation(0, 4) << std::endl;
      TFile *file = TFile::Open(externalLocation);
      if (file && !file->IsZombie()) { // Checks if webresource is available
        ::Info("AliExternalInfo::Cache", "Resource available");
        if (file->Cp(internalFilename)) {
          ::Info("AliExternalInfo::Cache", "Caching with TFile::Cp() successful");
          return kTRUE;
        } else {
          ::Error("AliExternalInfo::Cache", "Copying to internal location failed");
          return kFALSE;
        }
      } else {
        ::Error("AliExternalInfo::Cache", "Resource not available");
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


      TString command = CurlMif(mifFilePath, internalLocation, rootFileName, externalLocation);

      std::cout << command << std::endl;
      gSystem->Exec(command.Data());
      if (oldIndexName.Length()==0){
        gSystem->Exec(TString::Format("cat %s | sed -l 1 s/raw_run/run/ |  sed -l 1 s/RunNo/run/ > %s",mifFilePath.Data(),  (mifFilePath+"RunFix0").Data())); // use standard run number IDS
      }else{
        gSystem->Exec(TString::Format("cat %s | sed -l 1 s/%s/%s/  > %s",mifFilePath.Data(), oldIndexName.Data(), indexName.Data(),  (mifFilePath+"RunFix0").Data())); // use standrad run number IDS
      }

      gSystem->GetFromPipe(TString::Format("cat %s  | sed s_\\\"\\\$_\\\"0_g | sed s_\\\"\\\"_\\\"\\ \\\"_g | sed s_\\\"\\\"_\\\"\\ \\\"_g > %s",  (mifFilePath+"RunFix0").Data(),  (mifFilePath+"RunFix").Data()).Data());
      // Store it in a tree inside a root file
      TFile tempfile(internalFilename, "RECREATE");
      tempfile.cd();
      TTree tree(treeName, treeName);

      if ( (tree.ReadFile(mifFilePath+"RunFix", "", '\"')) > 0) {
        if (fVerbose>1) ::Info("AliExternalInfo::Cache", "Successfully read in tree");
      }
      else {
        ::Error("AliExternalInfo::Cache", "Error while reading tree");
        return kFALSE;
      }

      tree.Write();
      tempfile.Close();
      if (fVerbose>0) ::Info("AliExternalInfo::Cache", "Write tree to file: %s", internalFilename.Data());
      return kTRUE;
    }
  }
  else {// downloadIsNeeded == kFALSE
    return kTRUE;
  }
}

/// Cache selected production trees.  Input production list obtained from MonALISA web interface
/// \param select      - selection mask
/// \param reject      - rejection mask
/// \param sourceList  - list of detectors to cache
/*!
   Example usage:
   \code
        AliExternalInfo::CacheProduction(TPRegexp("LHC17.*"),TPRegexp("cpass0"),"QA.TPC;QA.EVS;QA.TRD;QA.rawTPC;QA.ITS;Logbook;QA.TOF;Logbook.detector");
   \endcode
*/
void AliExternalInfo::CacheProduction(TPRegexp select, TPRegexp reject, TString sourceList){
  AliExternalInfo info;
  TTree* treeProd = info.GetTreeProdCycle();
  Int_t entries=treeProd->GetEntries();
  TObjArray * detectorArray=sourceList.Tokenize(";");
  TString rejpat=reject.GetModifiers();
  for (Int_t i=0; i<entries; i++){
    treeProd->GetEntry(i);
    char * productionTag= (char*)treeProd->GetLeaf("Tag")->GetValuePointer();
    if (select.Match(productionTag)==0) continue;
    if (!rejpat.EqualTo("") && reject.Match(productionTag)==1) continue;
    printf("Caching\t%s\n",productionTag);
    TString production(productionTag);
    Int_t pos=production.First('_');
    if (pos<0) continue;
    if (pos>production.Length()-4) continue;
    printf("Caching\t%s\n",productionTag);
    TString period( production(0,pos));
    TString pass(production(pos+1, production.Length()-pos-1));
    printf("Caching\t%s\t%s\t%s\n",productionTag,period.Data(),pass.Data());
    for (Int_t iDet=0;iDet<detectorArray->GetEntries(); iDet++) {
      info.Cache(detectorArray->At(iDet)->GetName(), period.Data(), pass.Data());
    }
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
  if (fVerbose>1) ::Info("AliExternalInfo::GetTree", "Caching start internalFileName\t%s\trootFileName\t%s", internalFilename.Data(), rootFileName.Data());
  if (gSystem->AccessPathName(internalFilename.Data()) == kTRUE) {
    if (Cache(type, period, pass) == kFALSE) {
      ::Error("AliExternalInfo::GetTree", "Caching failed internalFileName\t%s\trootFileName\t%s", internalFilename.Data(), rootFileName.Data());
      return tree;
    }
  }else{
    Bool_t downloadNeeded = IsDownloadNeeded(internalFilename, type);
    if (downloadNeeded){
      ::Info("AliExternalInfo::GetTree", "Caching %s", internalFilename.Data());
      if (Cache(type, period, pass) == kFALSE) {
	::Error("AliExternalInfo::GetTree", "Caching failed internalFileName\t%s\trootFileName\t%s", internalFilename.Data(), rootFileName.Data());
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
    if (fVerbose>1) ::Info("AliExternalInfo::GetTree", "Successfully read %s/%s",internalFilename.Data(), tree->GetName());
    if (buildIndex==1) BuildIndex(tree, type);
  } else {
    //::Error("AliExternalInfo::GetTree", "Error while reading tree: ");
    ::Error("AliExternalInfo::GetTree", "ERROR READING: %s", treeName.Data());
    delete treefile;
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

  if (fLoadMetadata && metadataMacro.Length()>0){  // rename branch  with index if specified in configuration file
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
    if (fVerbose&0x2) ::Error("AliExternalInfo::GetTree","Tree %s not valid or empty",type.Data()); 
    return 0;
  }
  if (fVerbose&0x2) ::Info("AliExternalInfo::GetTree","Tree %s has %d entries",type.Data(),(int)tree->GetEntries()); 
  TString indexName= fConfigMap[type + ".indexname"];
  TString oldIndexName= fConfigMap[type + ".oldindexname"];
  if (oldIndexName.Length()>0 && tree->FindBranch(oldIndexName.Data())){
    tree->FindBranch(oldIndexName.Data())->SetName(indexName.Data());
  }
  if (indexName.Length()<=0) indexName="run";
  Int_t entries = tree->Draw(indexName.Data(),"","goff");
  if (entries<=0){
    if (fVerbose&0x2) ::Error("AliExternalInfo::GetTree","Friend tree %s not valid or empty",type.Data()); 
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
      if (fVerbose&0x2) ::Error("AliExternalInfo::GetTree","Friend tree %s not valid or empty",fname.Data()); 
      continue;
    }
    else if (fVerbose&0x2) ::Info("AliExternalInfo::GetTree","Friend tree %s has %d entries",fname.Data(),(int)ftree->GetEntries());
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
    if (fVerbose&0x2) ::Info("AliExternalInfo::GetTree","AddFriend %s+%s - entries=%d", type.Data(),  fname.Data(),fentries);
  }
  return tree;
}



/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// \param period Period, e.g. 'LHC15f'. Here you can use wildcards like in 'ls', e.g. 'LHC15*'
/// \param pass E.g. 'pass2' or 'passMC'. Here you can use wildcards like in 'ls', e.g. 'pass*'
/// Returns a chain with the information from the corresponding resources.
/// \return TChain*
TChain* AliExternalInfo::GetChain(TString type, TString period, TString pass, Int_t buildIndex){
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
  ::Info("AliExternalInfo::GetChain", "Files to add to chain: %s", files.Data());

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
  BuildIndex(chain,type);
  TString metadataMacro=fConfigMap[type + ".metadataMacro"];
  chain->Draw("Entry$","1","goff",1);
  if (metadataMacro.Length()>0 && chain->GetTree()) {  // rename branch  with index if specified in configuration file
    if (fVerbose>1) printf("Processing metadata macro:\n gROOT->ProcessLine(.x %s((TTree*)%p,%d);",     metadataMacro.Data(),chain->GetTree(), fVerbose);
    gROOT->ProcessLine(TString::Format(".x %s((TTree*)%p,%d);",metadataMacro.Data(),chain->GetTree(),fVerbose).Data());
  }

  delete arrFiles;
  delete arrTreeName;
  return chain;
};

/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// \param period Period, e.g. 'LHC15f'. Here you can use wildcards like in 'ls', e.g. 'LHC15*'
/// \param pass E.g. 'pass2' or 'passMC'. Here you can use wildcards like in 'ls', e.g. 'pass*'
/// Returns a chain with the information from the corresponding resources.
/// \return TChain*
TChain* AliExternalInfo::GetChain(TString type, TString period, TString pass, TString friendList){
  TChain *chain = GetChain(type.Data(),period.Data(),pass.Data(),kFALSE);
  if (chain==0){
    ::Error("AliExternalInfo::GetChain", "Invalid tree description %s\t%s\t%s",type.Data(), period.Data(),pass.Data());
  }
  TObjArray * arrFriendList= friendList.Tokenize(";");
  for (Int_t ilist=0; ilist<arrFriendList->GetEntriesFast(); ilist++) {

    TString fname=arrFriendList->At(ilist)->GetName();
    TString conditionName="";
    TString condition="";
    Int_t nDots = fname.CountChar(':');
    TChain *chainF =NULL;

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

    chainF=GetChain(fname.Data(), period.Data(), pass.Data(),kTRUE);

    if (chainF){
      if (nDots!=2) {
        chain->AddFriend(chainF, arrFriendList->At(ilist)->GetName());
      }else{
        chain->SetAlias(conditionName.Data(),"(1+0)");
        chainF->SetAlias(conditionName.Data(),condition.Data());
        chainF->BuildIndex(chainF->GetTreeIndex()->GetMajorName(), conditionName.Data());
        chain->AddFriend(chainF, (fname+"_"+conditionName).Data());
      }
    }else{
      ::Error("AliExternalInfo::GetChain", "Invalid friend tree\t%s\t%s",arrFriendList->At(ilist)->GetName(), friendList.Data());
      continue;
    }
  }
  return chain;
}

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
      tree->SetAlias(indexName.Data(),oldIndexName.Data());
    }
  }
  if (indexName.Length()<=0) { // set default index name
     indexName="run";
    if (tree->GetListOfBranches()!=NULL) if (tree->GetListOfBranches()->FindObject("run"))  indexName="run";
  }
  if (indexName.Length()<=0) {
    ::Error("AliExternalInfo::BuildIndex","Index %s not available for type %s", indexName.Data(), type.Data());
  }  
  if (tree->GetBranch(indexName.Data()) && TString(tree->GetBranch(indexName.Data())->GetTitle()).Contains("/C")){
    BuildHashIndex(tree,indexName.Data(),"hashIndex");
  }else{
    tree->BuildIndex(indexName.Data());
  }
  TStatToolkit::AddMetadata(tree,"TTree.indexName",indexName.Data());

  if (fVerbose&0x2) ::Info("AliExternalInfo::BuildIndex", "TreeName:%s;IndexName:%s",tree->GetName(), indexName.Data());
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
  ::Info("AliExternalInfo::AddChain", "Add to internal Chain: %s", internalFilename.Data());
  ::Info("AliExternalInfo::AddChain", "with tree name: %s",        treeName.Data());
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
  if (period=="sim")  return "/sim/";
  if (period=="data") return "/data/";
  //Check if period is MC and adjust storage hif::erarchy
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
    ::Info("AliExternalInfo::IsDownloadNeeded", "-- File not found locally --> Caching from remote");
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
      ::Info("AliExternalInfo::IsDownloadNeeded", "-- File is %li s old; NOT older than the set timelimit %d s",timeNow - timeFileModified, timeOut);
      return kFALSE; // if file is younger than the set time limit, it will not be downloaded again
    }
    else {
      ::Info("AliExternalInfo::IsDownloadNeeded", "-- File is %li s old; Older than the set timelimit %d s",timeNow - timeFileModified, timeOut);
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
const TString AliExternalInfo::CurlMif(TString& mifFilePath, const TString& internalLocation, TString rootFileName, const TString& externalLocation){
  TString command = "";
  TString certificate("$HOME/.globus/usercert.pem");
  TString privateKey("$HOME/.globus/userkey.pem");

  if(!gSystem->AccessPathName(certificate.Data())){
    ::Error("AliExternalInfo::CurlMif","Grid certificate can not be found in: %s",certificate.Data());
    return TString("No Certificate!");
  }
  if(!gSystem->AccessPathName(privateKey.Data())) {
    ::Error("AliExternalInfo::CurlMif","Grid private key can not be found in: %s",privateKey.Data());
    return TString("No private Key!");
  }
  // TString internalLocation = internalFilename(0, internalFilename.Last('/') + 1); // path to internal location
  // Create path to logfile with same name as the root file
  TString logFileName = rootFileName.ReplaceAll(".root", ".log");
  logFileName.Prepend(internalLocation);

  mifFilePath = rootFileName.ReplaceAll(".log", ".mif");
  mifFilePath.Prepend(internalLocation);

  command = TString::Format("curl -z %s -k --tlsv1 --cert %s --key %s -o %s 2>%s \"%s\"",
                                     mifFilePath.Data(), certificate.Data(), privateKey.Data(),
                                     mifFilePath.Data(), logFileName.Data(), externalLocation.Data());
  if ((fVerbose&0x4)>0) {
    ::Info("AliExternalInfo::Curl","%s",command.Data());
  }
  return command;
}

const TString AliExternalInfo::CurlTree(const TString internalFilename, const TString& externalLocation){
  TString command = "";
  TString certificate("$HOME/.globus/usercert.pem");
  TString privateKey("$HOME/.globus/userkey.pem");
  
  if(!gSystem->AccessPathName(certificate.Data())){
    ::Error("AliExternalInfo::CurlTree","Grid certificate can not be found in: %s",certificate.Data());
    return TString("No Certificate!");
  }
  if(!gSystem->AccessPathName(privateKey.Data())) {
    ::Error("AliExternalInfo::CurlTree","Grid private key can not be found in: %s",privateKey.Data());
    return TString("No private Key!");
  }

  command = TString::Format("curl -Lk -z %s --tlsv1 --cert %s --key %s -o %s \"%s\"",     //-L option required to get files from redirected URL
                                     internalFilename.Data(),certificate.Data(), privateKey.Data(),
                                     internalFilename.Data(), externalLocation.Data());

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
  char  pbranchName[10000];
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
//    and chaches the logbooks (and QA trees) which are needed for GetMCAnchPerGuess()  
   TTree * treeProd = GetTreeProdCycle();              // getting tree with information on real data productions (list, id, tag) - id will be used to get info for each production via GetTreeProdCycleByID(TString::Format("%d",id))
   TFile* outfile;
   TTree *logTree = NULL;

   Bool_t downloadNeeded = IsDownloadNeeded(fLocalStorageDirectory+TString::Format("/dumptree_RD.root"),TString::Format("QA.TPC"));     //check if download is needed
   if(!downloadNeeded){
        outfile = TFile::Open(fLocalStorageDirectory+TString::Format("%s","/dumptree_RD.root"),"UPDATE");
        if (outfile!=NULL) if(outfile->GetListOfKeys()->Contains("dumptree_RD")){
	  ::Info("AliExternalInfo::GetTreeAliVersRD", "-- dumptree_RD.root found locally and validated--> Not caching");
           return((TTree*)outfile->Get("dumptree_RD"));
        } 
    } 
        
   else  ::Info("AliExternalInfo::GetTreeAliVersRD", "-- dumptree_RD.root not validated--> Caching from remote");
   
    outfile= new TFile(fLocalStorageDirectory+"/dumptree_RD.root","RECREATE");
    TTree *dumptree = treeProd->CloneTree();       //tree that will hold information for guessing
    
    Int_t id=0;
    char tag[1000];
    
    dumptree->SetBranchAddress("ID",&id);
    dumptree->SetBranchAddress("Tag",&tag);
    
    char paliroot[1000];            // variables that will hold information read of the tree: GetTreeProdCycleByID(TString::Format("%d",id))
    char paliphysics[1000];
    char poutputdir[1000];
    Int_t runn;
    TString sRunList;

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
    TBranch* brrunlist= dumptree->Branch("runList",&sRunList);

    Int_t entries=dumptree->GetEntries();
    for (Int_t i=0; i<entries; i++){            //loop over all IDs
      dumptree->GetEntry(i);
      ::Info("AliExternalInfo::GetTreeAliVersRD", "Getting ProdCyle ID: %d",id);
      TTree * tree= GetTreeProdCycleByID(TString::Format("%d",id));         //get tree with production info for each ID

      if (tree==NULL) {
        ::Error("AliExternalInfo::GetTreeAliVersRD","GetTreeProdCycleByID(%d) returned bad tree",id);
        continue;
      }
      if (tree->GetBranch("app_aliphysics")==NULL) {
        TFile *ftree = ( tree->GetBranch("RunNo")!=NULL) ? tree->GetBranch("RunNo")->GetFile():NULL;
        delete tree;
        delete ftree;
        continue;
      }
      tree->SetBranchAddress("app_aliphysics",&paliphysics);      // set prod info branch addresses
      tree->SetBranchAddress("app_aliroot",&paliroot);
      tree->SetBranchAddress("outputdir",&poutputdir);
      tree->SetBranchAddress("RunNo",&runn);
      sRunList = "";
      
      //Loop over ProdCycleTree to get the string of run numbers (may be used when guessing the anchor pass)
      for(Int_t k=0; k<tree->GetEntries(); k++){
        tree->GetEntry(k);
        if(k!=0) sRunList.Append(TString(", "));
        sRunList.Append(TString::Format("%d",runn));
      }
      tree->GetEntry(0);
      soutputdir= TString::Format("%s",poutputdir);                 // extract production name from outputdir
      if(soutputdir.EndsWith("/")) soutputdir=soutputdir.Remove(soutputdir.Length()-1);  //if directory ends with "/" remove this
      
      subStrL = TPRegexp("(?=LHC)(.*?)(?=/)").MatchS(soutputdir);
      sprodname = *((TObjString *)subStrL->At(0)); 
      delete subStrL;
      subStrL = TPRegexp("[A-Za-z0-9]*").MatchS(sprodname.String());
      sprodname = *((TObjString *)subStrL->At(0));
      delete subStrL;

      subStrL = TPRegexp("[^/]+$").MatchS(soutputdir);              //extract pass name from utputdir
      spassname = *((TObjString *)subStrL->At(0));
      delete subStrL;      

      //GetLoogbook for anchor period and pass

      logTree = GetTree("Logbook", sprodname.GetString(), spassname.GetString(),"");
        if (logTree == NULL) {
          ::Error("AliExternalInfo::GetTreeAliVersRD", "Failed to get LogBook for %s %s",sprodname.GetString().Data(),spassname.GetString().Data());
      }      
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
      brrunlist->Fill();
      TFile *ftree = ( tree->GetBranch("RunNo")!=NULL) ? tree->GetBranch("RunNo")->GetFile():NULL;
      delete tree;
      if (ftree!=NULL) {
        delete ftree;  /// delete tree together with file
      }else{
        ::Error("AliExternalInfo::GetTreeAliVersRD","Invalif file");
      }

    }
    outfile->cd();
    dumptree->Write("dumptree_RD");
    delete dumptree;
    delete outfile;
    return GetTreeAliVersRD();
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
        if (!outfile){
          ::Error("AliExternalInfo::GetTreeAliVersMC","File dumptree_MC.root not available or corrupted. Recreating");
        }else {
	        ::Info(" AliExternalInfo::GetTreeAliVersMC", "-- dumptree_MC.root found locally and validated--> Not caching");
	        return((TTree*)outfile->Get("dumptree_MC"));
        } 
    } 
    
   if (fVerbose&0x2) ::Info("AliExternalInfo::GetTreeAliVersMC","-- dumptree_MC.root not validated--> Caching from remote");
   
   outfile= new TFile(fLocalStorageDirectory+"/dumptree_MC.root","RECREATE");
   TTree* dumptree=treeMC->CloneTree();
   
   //variable to read tree from GetTreeMC()
   char panchprodname[1000];
   char prunlist[50000];
   char pdescr[50000];
   char pprodname[50000];
   
   dumptree->SetBranchAddress("anchorProdTag",&panchprodname);
   dumptree->SetBranchAddress("runList",&prunlist);
   dumptree->SetBranchAddress("Description",&pdescr);
   dumptree->SetBranchAddress("prodName",&pprodname);
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
   for (Int_t i=0; i<entries; i++){           //loop overall MC production
     dumptree->GetEntry(i);                       //read info
     if (fVerbose&0x2) ::Info("AliExternalInfo::GetTreeAliVersMC","%d\tout of %d\n",i,entries);
     sanprod= TString::Format("%s",panchprodname);          //extract anchor production name from anchorProdTag
     subStrL = TPRegexp("[A-Za-z0-9]+").MatchS(sanprod);
     if(subStrL->GetLast()==-1) sanprod =TString("");
     else sanprod = ((TObjString *)subStrL->At(0))->GetString(); 
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
     
     if(sanprod==""){
       if (fVerbose&0x2) ::Info("AliExternalInfo::GetTreeAliVersMC","For MC prod %s no anchor period was specified in dumptree_MC - guessing based on run list",     pprodname);
       sanprod=GetMCAnchPerGuess(prunlist)+TString("?");
     }
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
   //dumptree->SetDirectory(gROOT);
   dumptree->Write("dumptree_MC");
  // TODO - fix memory leak in code above - either keep cache or delete intermediate trees )_
  delete treeMC;
  delete treeProdMC;
  delete dumptree;
  //delete outfile;
  outfile->Close();
  return GetTreeAliVersMC();
}



class anchprod{             //class that holds information about aliphysics and aliroot version, the pass name and the production name - will be used to sort the real data and MC productions
public:
TString aliphys;
TString aliroot;
TString anchpass;
TString anchprodname;
int runNumMatch;
int runN;

anchprod(){
aliphys="-1";
aliroot="-1";
anchpass="-1";
anchprodname="-1";
runNumMatch=-1;
runN=-1;
}

bool operator< (const anchprod & otheranchprod) const       //define comparison operator for lexicographic ordering, older AliPhys/AliRoot first and lower runnumber match first
	{
            if(aliphys+aliroot != otheranchprod.aliphys+otheranchprod.aliroot){
		return (aliphys+aliroot < otheranchprod.aliphys+otheranchprod.aliroot);         //compare first aliphysics, if same (i.e. empty) then aliroot, duplicates in ali versions would be stored only once also if passes different
            }
            else return(runNumMatch<otheranchprod.runNumMatch);
        }
};
 
TTree*  AliExternalInfo::GetTreeMCPassGuess(){
    Bool_t loadMetadataBackup= fLoadMetadata;
    fLoadMetadata=kFALSE;
//    returns and stores (dumptree_MC.root) the tree containing the pass guesses for each MC production
    TFile *MCFileG=0;
    TTree *MCTreeG=0;
    TTree *MCTree=0;
        
    //TFile *RDFile;
    TTree *RDTree;
    
    //Check if fMCGuessTree has already been cached   
    Bool_t downloadNeeded = IsDownloadNeeded(fLocalStorageDirectory+TString::Format("/dumptree_MC_guess.root"),TString::Format("QA.TPC"));
    if(!downloadNeeded){
      MCFileG= TFile::Open(fLocalStorageDirectory+"/dumptree_MC_guess.root","UPDATE");
      if (MCFileG==NULL){
        ::Info("AliExternalInfo::GetTreeMCPassGuess","dumptree_MC_guess.root corrupted -> Overwrite it"); /// TODO - is it save ?
      }else {
        if (MCFileG->GetListOfKeys()->Contains("dumptree_MC_guess")) {
          if (fVerbose & 0x2) ::Info("AliExternalInfo::GetTreeMCPassGuess", "Guesses already available - done.");
          return (dynamic_cast<TTree *>(MCFileG->Get("dumptree_MC_guess")));
        }
      }
    }
    if (fVerbose&0x2) ::Info("AliExternalInfo::GetTreeMCPassGuess","dumptree_MC_guess.root not available -> get it");
   
    RDTree=GetTreeAliVersRD();   

    TObjString* osrdprod=0;             //variables for reading RD tree
    TObjString* osrdpass=0;
    TObjString* osrdaliphys=0;
    TObjString* osrdaliroot=0;

    TObjArray* arrRDrunlist=0;
    TString* sRDrunlist=0;
    
    RDTree->GetBranch("aliphysics")->SetAddress(&osrdaliphys);
    RDTree->GetBranch("aliroot")->SetAddress(&osrdaliroot);
    RDTree->GetBranch("prodName")->SetAddress(&osrdprod);
    RDTree->GetBranch("passName")->SetAddress(&osrdpass); 
    RDTree->GetBranch("runList")->SetAddress(&sRDrunlist);

    MCTree = GetTreeAliVersMC();    
    MCFileG=TFile::Open(fLocalStorageDirectory+"/dumptree_MC_guess.root","RECREATE");
    MCTreeG=new TTree("dumptree_MC_guess","dumptree_MC_guess"); 
    
    TObjString* osMCaliroot=0;          //char arrays for reading from MCTree
    TObjString* osMCaliphysics=0;
    TObjString* osMCprodname=0;
    TObjString* osMCanchprodname=0; //variables for reading TObjString from MCTree
    TObjString*osMCanchprodnameOrig=0;
    TObjString* osMCdescr=0;
    TObjString* osMCRunList=0;
    TString sMCanchprodname;      //holds anchor period name after "?" was removed ("?" indicates that anchor period was not specified, but guessed using the run list and the logbook)

    Bool_t isgp = kFALSE;         //is MC production a general purpose production?
        
    char pMCaliroot[1000];     //variable to read tree from GetTreeMC()
    char pMCaliphysics[1000];
    char pMCprodname[1000];
    char pMCrunlist[10000];
    TString sMCrunlist;
    TObjArray* arrMCrunlist=0;
    Int_t matches=0;
    Int_t runNMC=0;
    Int_t runNAnchor=0;    
    Int_t rank=0;
    
    MCTree->GetBranch("aliphysics")->SetAddress(&pMCaliphysics);       //set branch addresses
    MCTree->GetBranch("aliroot")->SetAddress(&pMCaliroot);
    MCTree->GetBranch("prodName")->SetAddress(&pMCprodname);
    MCTree->GetBranch("anchorProdTag_ForGuess")->SetAddress(&osMCanchprodname);
    MCTree->GetBranch("Description")->SetAddress(&osMCdescr);
    MCTree->GetBranch("runList")->SetAddress(&pMCrunlist); 
    MCTreeG->Branch("MCProdName",&osMCprodname);
    MCTreeG->Branch("MCAliphysics",&osMCaliphysics);       //set branch addresses
    MCTreeG->Branch("MCAliroot",&osMCaliroot);
    MCTreeG->Branch("MCDescription",&osMCdescr);
    MCTreeG->Branch("MCRunList",&osMCRunList);
    MCTreeG->Branch("AnchorProdTag",&osMCanchprodnameOrig);
    MCTreeG->Branch("AnchorPassName",&osrdpass);         //adding new branches holding information about pass and ali-versions of guessed RD production
    MCTreeG->Branch("Anchoraliroot",&osrdaliroot);
    MCTreeG->Branch("Anchoraliphys",&osrdaliphys);
    MCTreeG->Branch("rankGuess",&rank);
    MCTreeG->Branch("runNMatches",&matches);
    MCTreeG->Branch("runNMC",&runNMC);
    MCTreeG->Branch("runNAnchor",&runNAnchor);
    
    int n = RDTree->GetEntries();
    int m = MCTree->GetEntries();
    Bool_t prfound= kFALSE; 
    Bool_t onPassBlackList=kFALSE;
    
    anchprod tempprod;

    multiset<anchprod> list;                                                        //set that will hold "anchprod" instances of RD productions, that have have matching production names
    multiset<anchprod>::iterator it = list.begin();

    for (Int_t i=0; i<m; i++) {     //loop over MC productions

      MCTree->GetEntry(i);
      isgp = TPRegexp("General").MatchB(osMCdescr->String(),"i") && TPRegexp("Purpose").MatchB(osMCdescr->String(),"i");

      osMCaliroot = new TObjString(pMCaliroot);          //get TObjStrings for guessing
      osMCaliphysics= new TObjString(pMCaliphysics);
      osMCprodname= new TObjString(pMCprodname);

      osMCanchprodnameOrig=new TObjString(osMCanchprodname->String());   //Copy Anchor Production before removing "?"
      
      cout<<endl;
      cout<<i<<" of "<<m<<endl;
      cout<<"MC Production name: "<<osMCprodname->String()<<" Anchor Production name: "<<osMCanchprodname->String()<<" MC aliphys: "<<osMCaliphysics->String()<<" MC aliroot: "<<osMCaliroot->String()<<" MC description: "<<osMCdescr->String()<<" isgp: "<<isgp<<endl;
      sMCanchprodname=osMCanchprodname->String().ReplaceAll("?","");
      if(osMCanchprodname->String().Contains("?")) cout<<"MC anchor production was guessed based on run list comparison with logbook! Use "<<sMCanchprodname<<" as anchor period for pass guessing"<<endl;

      prfound=kFALSE;         // flag to know if any RD production with matching production name was found
      list.clear();                  //reset set of RDinfos

     if(sMCanchprodname==""){
     if (fVerbose&0x2) ::Warning("AliExternalInfo::GetTreeMCPassGuess", "No anchor pass name provided for MC prod %s -> Guesses set to NONE",osMCprodname->String().Data());
        osMCanchprodname->SetString("NONE");
        osrdpass->SetString("NONE");
        cout<<"Used for guess: RDpass guess: NONE"<<endl;
        rank=-1;
        matches=-1;
        runNMC=-1;
        runNAnchor=-1;
        MCTreeG->Fill();
      }        
     else{     //if AnchorProd specified

      for (Int_t j=0; j<n; j++) {        //search for matching RD production and get pass info

       RDTree->GetEntry(j);
       
        if(sMCanchprodname==osrdprod->String()){           //if RDprodname is right then save the prod infos to be later eventually able to look up what was closest in terms of aliphys/aliroot
                                                      //found a prodction with correct production name
          tempprod.aliphys=osrdaliphys->GetString();                //store info in list of "anchprod"
          tempprod.aliroot=osrdaliroot->GetString();
          tempprod.anchpass=osrdpass->GetString();
          tempprod.anchprodname=osrdprod->GetString();

          //Get number of runlist matches between current MC and RD run list
          sMCrunlist = TString(pMCrunlist);
          arrMCrunlist = sMCrunlist.Tokenize(" ,;\t");

          arrRDrunlist = sRDrunlist->Tokenize(" ,;\t");
          matches=0;
          runNMC=arrMCrunlist->GetEntries();
          tempprod.runN=arrRDrunlist->GetEntries();
          runNAnchor=arrRDrunlist->GetEntries();           
          for(int l=0; l<runNMC; l++){
            if(l==0 && j==0)    cout<<"Number of entries MC run list: "<<arrMCrunlist->GetEntries()<<endl;
            for(int k=0; k<runNAnchor; k++){
               if(((TObjString *)arrMCrunlist->At(l))->GetString().Atoi()==((TObjString *)arrRDrunlist->At(k))->GetString().Atoi()) matches++;
               }
          }
           
          tempprod.runNumMatch=matches;
          onPassBlackList=kFALSE;
          TString passBlackList = TString("cpass,vpass,vdm,cosmic,muon,fast,pass0,its,recpoint,cleanesd");
          TObjArray* arrBlack= passBlackList.Tokenize(",");
          for(int m=0; m<arrBlack->GetEntries(); m++){
              if( (tempprod.anchpass).Contains(((TObjString *)arrBlack->At(m))->GetString(),TString::kIgnoreCase)) onPassBlackList=kTRUE;
            }
           
          if(!onPassBlackList && (!isgp || TPRegexp("^pass").MatchB(tempprod.anchpass,"i"))) {
            list.insert(tempprod);
            prfound=kTRUE;
          }
         }            

         if(j==n-1 && prfound){            //if match not found check what was closest if matching prodname was found
           tempprod.aliphys=osMCaliphysics->String();         //make anchprd instance with MC info and insert into its lexicographical position
           tempprod.aliroot=osMCaliroot->String();
           tempprod.anchprodname=osMCprodname->String();
           tempprod.anchpass=TString("MCPass");               //dummy info to make visible in list what was MC entry
           tempprod.runNumMatch=9999; 
           tempprod.runN=runNMC;
           list.insert(tempprod);     //insert MC prod 
           it=list.find(tempprod);         //get iterator to pointer before MC prod
           
           if (it == list.begin()){
             prfound =kFALSE;       //no interesting RD productions listed before MC - guess NONE
             cout<<"No RD period was listed before MC production: Used for guess: RDphys:NONE, RDroot: NONE, RDpass guess: NONE"<<endl;
             *osrdpass=TObjString("NONE");
             *osrdaliphys=TObjString("NONE");
             *osrdaliroot=TObjString("NONE");
             rank=-1;
             matches=-1;
             runNMC=-1;
             runNAnchor=-1;
             MCTreeG->Fill();
             
             break;
            }     //if passing here there must be at least one guess
           
          int l =0;
          cout<<"Take closest aliversion and then most run number matches as guess from RD productions: "<<endl;
          for (multiset<anchprod>::iterator iter=list.begin(); iter!=list.end(); ++iter){         //cout the ordered list of productions
             l++;
            cout<<l<<" Prodname: "<<(*iter).anchprodname<<" AliPhys: "<<(*iter).aliphys<<" AliRoot: "<<(*iter).aliroot<<" run number matches: "<<(*iter).runNumMatch<<" (length RunList:"<<(*iter).runN<<") pass: "<<(*iter).anchpass<<endl;
          }

          --it;                                                       //let iterator point to entry right before MC entry that is neither a "cpass" nor a "cosmics" pass
          rank=0;
          for (multiset<anchprod>::iterator iter=it; ; --iter){       //go backwards through list and take first pass guess that is not a cpass
            rank++;
            cout<<"Guess #"<<rank<<" : RDphys:"<<(*iter).aliphys<<" RDroot: "<<(*iter).aliroot<<" RDpass guess: "<<(*iter).anchpass<<endl;
            *osrdpass=TObjString(iter->anchpass);
            *osrdaliphys=TObjString(iter->aliphys);
            *osrdaliroot=TObjString(iter->aliroot);

            MCTreeG->Fill();     

            if(iter==list.begin()){    //if list contains no viable guess set Anchor and AnchorPass guess to "NONE"

              break;
            }
         }
        }
        else if(j==n-1 && !prfound){      //if RD contained no matching period
            rank=-1;
            matches=-1;
            runNMC=-1;
            runNAnchor=-1;
            cout<<"No RD period found that matches the Anchor Period: Used for guess: RDphys:NONE, RDroot: NONE, RDpass guess: NONE"<<endl;
            *osrdpass=TObjString("NONE");
            *osrdaliphys=TObjString("NONE");
            *osrdaliroot=TObjString("NONE");

            MCTreeG->Fill();         
            break;   
        }
       }
      }
     }
    
    MCFileG->cd();
    MCTreeG->Write();
//    delete MCFileG;
    fLoadMetadata=loadMetadataBackup;
    return MCTreeG;
}

TTree* AliExternalInfo::GetLogbookCache(){
//  cout<<"here"<<endl;
  Bool_t downloadNeeded = IsDownloadNeeded(fLocalStorageDirectory+TString::Format("/logbook.root"),TString::Format("Logbook"));   //check if file containing guesses is present
  if(!downloadNeeded) {
     if (fVerbose&0x2)::Info("AliExternalInfo::GetLogbookCache","logbook.root available");
     TFile* logFile= TFile::Open(fLocalStorageDirectory+"/logbook.root");
     if( !(logFile->GetListOfKeys()->Contains("logbook"))){
       if (fVerbose&0x2) ::Info("AliExternalInfo::GetLogbookCache","logbook tree not available");
       logFile=new TFile(fLocalStorageDirectory+"/logbook.root","RECREATE");
       TChain* tempLog=GetChain("Logbook","LHC1**","");
       TTree* tempTree=tempLog->CloneTree(-1,"fast");
       logFile->cd();
       tempTree->Write();

       return tempTree;
     }
     else{
         if (fVerbose&0x2) ::Info("AliExternalInfo::GetLogbookCache","logbook available");
         TTree* tempTree = dynamic_cast<TTree*>(logFile->Get("logbook"));
         return tempTree;
     }
  }
  else{
      if (fVerbose&0x2)::Info("AliExternalInfo::GetLogbookCache","logbook file not available ");
      TChain* tempLog=GetChain("Logbook","LHC1**","");
      TFile* logFile=new TFile(fLocalStorageDirectory+"/logbook.root","RECREATE");
       
       TTree* tempTree=tempLog->CloneTree(-1,"fast");
       logFile->cd();
       tempTree->Write();
       return tempTree;
  }  
}



TString  AliExternalInfo::GetMCAnchPerGuess(const char* prunlist){

/*
AliExternalInfo* info = new AliExternalInfo(".", "", 2);
info->GetMCPassGuess("LHC17l1");
*/ 
  if(fLogCache==0){
    ::Info("AliExternalInfo::GetMCAnchPerGuess","Logbook not already opened - getting it!");
    fLogCache = GetLogbookCache();
  }
  else {
    if (fVerbose>2) ::Info("AliExternalInfo::GetMCAnchPerGuess","Logbook already opened!");
  }
  
  std::string * pAnchProdName = new std::string();
  Int_t run;  
  fLogCache->GetBranch("LHCperiod")->SetAddress(&pAnchProdName); 
  fLogCache->GetBranch("run")->SetAddress(&run);
  
  std::map<std::string,int> map;
  std::string sProdGuess;

  if (fVerbose&0x4) cout<<"RunList: "<<prunlist<<endl;

  TString s = TString(prunlist);
  TObjArray* arr= s.Tokenize(" ,;\t");
  Int_t nRuns=arr->GetEntries();

  for(int k=0; k<fLogCache->GetEntries(); k++){
    fLogCache->GetEntry(k);
    for(int l=0; l<nRuns; l++){
      if( std::atoi(((TObjString *)arr->UncheckedAt(l))->GetString().Data())==run){
        if (fVerbose&0x4) cout<<"Anchor Prod Guess : "<<pAnchProdName->c_str()<<" for run number: "<<((TObjString *)arr->UncheckedAt(l))->GetString()<<endl;
        ++map[*pAnchProdName];
      }
    }
  }
  if(map.size()==0)  {
    if (fVerbose&0x4) ::Warning("AliExternalInfo::GetMCAnchPerGuess","No Anchor Production found via run number matching");

    return "NONE";
  }
  else {
    int max=0;

  std::map<std::string, int>::iterator it;
  for ( it = map.begin(); it != map.end(); it++ )
  {
      if (fVerbose&0x4)std::cout << it->first<< ':'<< it->second << std::endl;
      if(it->second>max) {
        sProdGuess=it->first;
        max=it->second;
      }
    }
    if (fVerbose&0x4) ::Info("AliExternalInfo::GetMCAnchPerGuess","AnchorProd guess is %s based on %d run number hits",sProdGuess.c_str(),max);  
    return(TString(sProdGuess.c_str()));
  }
return "NONE";
}



TString  AliExternalInfo::GetMCPassGuess(TString sMCprodname){
/*
AliExternalInfo* info = new AliExternalInfo(".", "", 2);
info->GetMCPassGuess("LHC17l1");
*/     
//returns string with Pass guess for a given MC production name    

 TObjString* osMCprodname=0;
 TObjString* osAnchprodname=0;
 TObjString* osMCpassguess=0;
 Int_t rank;

 fMCGuessTree=GetTreeMCPassGuess();

 cout<<"MC Pass Guess Tree Entries: "<<fMCGuessTree->GetEntries()<<endl;
 fMCGuessTree->GetBranch("AnchorProdTag")->SetAddress(&osAnchprodname); 
 fMCGuessTree->GetBranch("MCProdName")->SetAddress(&osMCprodname);
 fMCGuessTree->GetBranch("AnchorPassName")->SetAddress(&osMCpassguess);
 fMCGuessTree->GetBranch("rankGuess")->SetAddress(&rank);
 
 for(int i=0;i<fMCGuessTree->GetEntries();i++){        //loop over tree with guesses
    fMCGuessTree->GetEntry(i);
    if(osMCprodname->String()==sMCprodname && rank==1){       //if best guess found
      if (fVerbose&0x2) ::Info("AliExternalInfo::GetMCPassGuess","Anchor Production:%s, Pass guess for %s: %s ",osAnchprodname->String().Data(),osMCprodname->String().Data(),osMCpassguess->String().Data());
         return(osAnchprodname->String()+" "+osMCpassguess->String());
   }
 }
 if (fVerbose&0x2) ::Error("AliExternalInfo::GetMCPassGuess", "%s was not found in list of MC productions", osMCprodname->String().Data());
   return(TString::Format("MC production not found"));
}
