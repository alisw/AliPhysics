#include <fstream>
#include <iostream>
#include <map>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>

#include "AliLog.h"

#include "AliExternalInfo.h"

#include "TObjArray.h"
#include "TString.h"
#include "TWebFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"

ClassImp(AliExternalInfo)

const TString AliExternalInfo::fgkDefaultConfig="$ALICE_ROOT/STAT/Macros/AliExternalInfo.cfg";

AliExternalInfo::AliExternalInfo(TString localStorageDirectory, TString configLocation/*, Bool_t copyToLocalStorage*/) :
                                /*fCopyDataToLocalStorage(copyToLocalStorage),*/
                                TObject(),
                                fConfigLocation(configLocation),
                                fLocalStorageDirectory(localStorageDirectory),
                                fLocationTimeOutMap(),
                                fTree(0x0),
                                fChain(new TChain()),
                                fChainMap(),
                                fMaxCacheSize(-1)
{
  ReadConfig();
}

AliExternalInfo::~AliExternalInfo() {}


/// Reads the configuration file. Lines beginning with an '#' are ignored.
/// Use the format which is in the config.cfg by default. Adding ressources like the ones already
/// there should work without problems.
void AliExternalInfo::ReadConfig(){
  TString configFileName=gSystem->ExpandPathName(fConfigLocation.Data());
  if (gSystem->AccessPathName(configFileName)) {
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

    fLocationTimeOutMap[key] = value;
  }
  return;
}

/// Prints out the config which was read in previously. Useful to check if anything went wrong
void AliExternalInfo::PrintConfig(){
  // Loop through the map (Would be much easier in c++11)
  std::cout << "User defined resources are:\n";
  // looping over map with const_iterator
  typedef std::map<TString,TString>::const_iterator it_type;
  for(it_type iterator = fLocationTimeOutMap.begin(); iterator != fLocationTimeOutMap.end(); ++iterator) {
    std::cout << iterator->first << " " << iterator->second << "\n";
  }
  return;
}


/// Sets up all variables according to period, pass and type. Extracts information from the config file
void AliExternalInfo::SetupVariables(TString& internalFilename, TString& internalLocation, Bool_t& resourceIsTree, TString& pathStructure, \
                                     TString& detector, TString& rootFileName, TString& treeName, const TString& type, const TString& period, const TString& pass){
  // Check if resource is a tree in a root file or not
  pathStructure = CreatePath(type, period, pass);

  // if (fLocationTimeOutMap.count(type + ".treename") > 0) resourceIsTree = kTRUE;
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

  rootFileName = fLocationTimeOutMap[type + ".filename"];
  treeName     = fLocationTimeOutMap[type + ".treename"];

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

  // initialization of external variables
  externalLocation = fLocationTimeOutMap[type + ".location"];

  // Setting up all the local variables
  SetupVariables(internalFilename, internalLocation, resourceIsTree, pathStructure, detector, rootFileName, treeName, type, period, pass);

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

      TString mifFilePath = ""; // Gets changed in Wget command
      TString command = Wget(mifFilePath, internalLocation, rootFileName, externalLocation);

      std::cout << command << std::endl;
      gSystem->Exec(command.Data());

      // Store it in a tree inside a root file
      TTree tree(treeName, treeName);

      if ( (tree.ReadFile(mifFilePath, "", '\"')) > 0) {
        AliInfo("-- Successfully read in tree");
      }
      else {
        AliError("-- Error while reading tree");
        return kFALSE;
      }

      TFile tempfile(internalFilename, "RECREATE");
      tempfile.cd();
      tree.Write();
      tempfile.Close();
      AliInfo(TString::Format("Write tree to file: %s", internalFilename.Data()));
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
TTree* AliExternalInfo::GetTree(TString type, TString period, TString pass){
  TString internalFilename = ""; // Resulting path to the file
  TString internalLocation = fLocalStorageDirectory; // Gets expanded in this function to the right directory
  TString externalLocation = "";
  Bool_t resourceIsTree = kFALSE;
  TString detector = "";
  TString rootFileName = "";
  TString treeName = "";
  TString pathStructure = "";

  TTree* tree = 0x0;

  // Setting up all the local variables
  SetupVariables(internalFilename, internalLocation, resourceIsTree, pathStructure, detector, rootFileName, treeName, type, period, pass);

  std::cout << "internalFilename: " << internalFilename << " rootFileName: " << rootFileName << std::endl;

  if (gSystem->AccessPathName(internalFilename.Data()) == kTRUE) {
    if (Cache(type, period, pass) == kFALSE) {
      std::cout << "Caching of ressource was not successful; Nullpointer is returned!\n" << std::endl;
      return tree;
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

  if (tree != 0x0) {
    AliInfo("-- Successfully read in tree");
    AddTree(tree, type);
  } else {
    AliError("Error while reading tree");
  }

  const TString cacheSize=fLocationTimeOutMap[type + ".CacheSize"];
  Long64_t cache=cacheSize.Atoll();
  if (fMaxCacheSize>0) {
    if (cache>0) {
      cache=TMath::Min(fMaxCacheSize, cache);
    } else {
      cache=fMaxCacheSize;
    }
  }
  if (cache>0) tree->SetCacheSize(cache);

  return tree;
}

/// \param type Type of the resource as described in the config file, e.g. QA.TPC, MonALISA.RCT
/// \param period Period, e.g. 'LHC15f'. Here you can use wildcards like in 'ls', e.g. 'LHC15*'
/// \param pass E.g. 'pass2' or 'passMC'. Here you can use wildcards like in 'ls', e.g. 'pass*'
/// Returns a chain with the information from the corresponding resources.
/// \return TChain*
TChain* AliExternalInfo::GetChain(TString type, TString period, TString pass){
  TChain* chain = 0x0;
  TString internalFilename = ""; // Resulting path to the file
  TString internalLocation = fLocalStorageDirectory; // Gets expanded in this function to the right directory
  TString externalLocation = "";
  Bool_t resourceIsTree = kFALSE;
  TString detector = "";
  TString rootFileName = "";
  TString treeName = "";
  TString pathStructure = "";

  // Setting up all the local variables
  SetupVariables(internalFilename, internalLocation, resourceIsTree, pathStructure, detector, rootFileName, treeName, type, period, pass);

  TString cmd = TString::Format("/bin/ls %s", internalFilename.Data());
  // std::cout << "==== cmd: " << cmd << std::endl;

  TString files=gSystem->GetFromPipe(cmd.Data());
  TObjArray *arrFiles=files.Tokenize("\n");
  AliInfo(TString::Format("Files to add to chain: %s", files.Data()));

  //function to get tree namee based on type
  chain=new TChain(treeName.Data());
  for (Int_t ifile=0; ifile<arrFiles->GetEntriesFast(); ++ifile) {
    chain->AddFile(arrFiles->At(ifile)->GetName());
  }

  const TString cacheSize=fLocationTimeOutMap[type + ".CacheSize"];
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
  return chain;
};

/// Every tree you create is added to a big tree acting as a friend.
/// You can have access to this tree with the GetFriendsTree() function.
/// @TODO Add 'return false' when adding to the friends tree was not successful
/// \return kTRUE
Bool_t AliExternalInfo::AddTree(TTree* tree, TString type){

  if (tree->BuildIndex("run") < 0) tree->BuildIndex("raw_run");
  TString name = "";

  if(type.Contains("QA")){ // use TPC instead of QA.TPC
    name = type(3, type.Length()-1);
  }
  else if(type.Contains("MonALISA")){
    name = type(9, type.Length()-1);
  }
  else {
    name = type;
  }

  tree->SetName(name);
  if (fTree == 0x0) fTree = dynamic_cast<TTree*>(tree->Clone());
  fTree->AddFriend(tree, name);

  AliInfo(TString::Format("Added as friend with the name: %s",name.Data()));
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

  // Setting up all the local variables
  SetupVariables(internalFilename, internalLocation, resourceIsTree, pathStructure, detector, rootFileName, treeName, type, period, pass);
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
  if (period.Length() == 6 || (period == "" && type != "MonALISA.MC") || type == "MonALISA.ProductionCycleID" || type == "TriggerClasses") { // put everything in the form LHCYYx or with empty period and pass in "data"
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
  timeOut = atoi(fLocationTimeOutMap[type + ".timeout"]);

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
/// Composes the wget-command in a TString which afterwards then can be executed
/// \return wget-command in a TString
const TString AliExternalInfo::Wget(TString& mifFilePath, const TString& internalLocation, TString rootFileName, const TString& externalLocation){
  TString command = "";
  TString certificate("$HOME/.globus/usercert.pem");
  TString privateKey("$HOME/.globus/userkey.pem");

  // TString internalLocation = internalFilename(0, internalFilename.Last('/') + 1); // path to internal location
  // Create path to logfile with same name as the root file
  TString logFileName = rootFileName.ReplaceAll(".root", ".log");
  logFileName.Prepend(internalLocation);

  mifFilePath = rootFileName.ReplaceAll(".log", ".mif");
  mifFilePath.Prepend(internalLocation);

  command = TString::Format("wget --no-check-certificate --secure-protocol=TLSv1 --certificate=%s --private-key=%s -o %s -O %s \"%s\"",
                                     certificate.Data(), privateKey.Data(), logFileName.Data(),
                                     mifFilePath.Data(), externalLocation.Data());
  return command;
}
