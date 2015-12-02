#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>

#include "AliExternalInfo.h"
#include "TSystem.h"
#include "TWebFile.h"
#include "TString.h"
#include "TChain.h"

ClassImp(AliExternalInfo)

AliExternalInfo::AliExternalInfo() :  fCertificate("$HOME/.globus/usercert.pem"),
                                      fPrivateKey("$HOME/.globus/userkey.pem"),
                                      fGlobalDir("."),
                                      fFormMCFilename("MC"),
                                      fFormMCDir("/sim/"),
                                      fFormMCFiletype(".mif"),
                                      fFormRCTFilename("RCT"),
                                      fFormRCTDir("/data/%s/%s/%s/"),
                                      fFormRCTFiletype(".mif"),
                                      fFormLogbookFilename("logbook"),
                                      fFormLogbookDir("/data/%s/%s/%s/"),
                                      fFormLogbookFiletype(".root"),
                                      fFormTriggerClassesFilename("trigger_classes"),
                                      fFormTriggerClassesDir("/data/%s/"),
                                      fFormTriggerClassesFiletype(".root"),
                                      fFormTrendingFilename("trending"),
                                      fFormTrendingDir("/data/%s/%s/%s/"),
                                      fFormTrendingFiletype(".root"),
                                      fFormProdCycleFilename("ProdCycle"),
                                      fFormProdCycleDir("/data/"),
                                      fFormProdCycleFiletype(".mif"),
                                      fFormProdPassesFilename("ProdPasses"),
                                      fFormProdPassesDir("/data/"),
                                      fFormProdPassesFiletype(".mif"),
                                      fTimeLimit(60 * 60 * 24)
{}
AliExternalInfo::AliExternalInfo(TString GlobalPath) : 
                                      fCertificate("$HOME/.globus/usercert.pem"),
                                      fPrivateKey("$HOME/.globus/userkey.pem"),
                                      fGlobalDir(GlobalPath),
                                      fFormMCFilename("MC"),
                                      fFormMCDir("/sim/"),
                                      fFormMCFiletype(".mif"),
                                      fFormRCTFilename("RCT"),
                                      fFormRCTDir("/data/%s/%s/%s/"),
                                      fFormRCTFiletype(".mif"),
                                      fFormLogbookFilename("logbook"),
                                      fFormLogbookDir("/data/%s/%s/%s/"),
                                      fFormLogbookFiletype(".root"),
                                      fFormTriggerClassesFilename("trigger_classes"),
                                      fFormTriggerClassesDir("/data/%s/"),
                                      fFormTriggerClassesFiletype(".root"),
                                      fFormTrendingFilename("trending"),
                                      fFormTrendingDir("/data/%s/%s/%s/"),
                                      fFormTrendingFiletype(".root"),
                                      fFormProdCycleFilename("ProdCycle"),
                                      fFormProdCycleDir("/data/"),
                                      fFormProdCycleFiletype(".mif"),
                                      fFormProdPassesFilename("ProdPasses"),
                                      fFormProdPassesDir("/data/"),
                                      fFormProdPassesFiletype(".mif"),
                                      fTimeLimit(60 * 60 * 24)
{}
AliExternalInfo::~AliExternalInfo(){ /*delete fTree;*/}


// Genereric function to download the information into a specific folder structure

Bool_t AliExternalInfo::Cache(TString period, TString pass, Int_t type){
  std::cout << "Caching started: period=" << period << " pass= " << pass << std::endl;
  // Needs to be adjusted by the later switch case for correct download and storage
  TString remotepath("");
  TString localpath("");
  TString file("");
  TString command("");
  TString year = GetYearFromPeriod(period);
  switch (type){
    case kRCT:
    {
      // sanity checks
      if ( (CheckPeriod(period) == kFALSE) || (CheckPass(pass) == kFALSE) ) return kFALSE;
      // end of sanity checks
      TString numberOfPass("");
      numberOfPass = pass[4];
      remotepath = (Form("https://alimonitor.cern.ch/configuration/index.jsp?partition=%s&pass=%s&res_path=mif", period.Data(), numberOfPass.Data()));
      SetUpForPeriodPass(period, pass, localpath, file, fFormRCTFilename, fFormMCFiletype, type);
      // std::cout << period << " " << pass << " " << file << " " << title << " " << localpath << " " << type << std::endl;
      command = wget(localpath, fFormRCTFilename, fFormMCFiletype, remotepath);
      break;
    }
    case kMC:
    {
      remotepath = "https://alimonitor.cern.ch/MC/?res_path=mif";
      SetUpForPeriodPass(period, pass, localpath, file, fFormMCFilename, fFormMCFiletype, type);
      command = wget(localpath, fFormMCFilename, fFormMCFiletype, remotepath);
      break;
    }
    case kProdCycle:
    {
      if (period == "") // for downloading whole master database
      {
        remotepath = TString::Format("https://alimonitor.cern.ch/prod/?t=1&res_path=mif");
      }
      else if ( pass == "") // for downloading a specific id
      {
        remotepath = TString::Format("https://alimonitor.cern.ch/prod/jobs.jsp?t=%s&res_path=mif", period.Data());
      }
      else {  // for searching a specific pass and period, NOT YET IMPLEMENTED
        std::cout << "!!!!!!!!!!!!!!!!      NOT IMPLEMENTED      !!!!!!!!!!!!!!!!!!!" << std::endl;
        TTree* tempTree;
        tempTree = GetTreeProdCycle();

        char* description = new char[200]();
        char* tag = new char[200]();
        Int_t ID;
        tempTree->SetBranchAddress("Description", description);
        tempTree->SetBranchAddress("Tag", tag);
        tempTree->SetBranchAddress("ID", &ID);
        Int_t nEntries = tempTree->GetEntriesFast();

        for (Int_t i = 0; i < nEntries; ++i){
          tempTree->GetEntry(i);
          std::cout << "ID: " << ID << "  Description: " << description << "\n" << "Tag: " << tag << "\n";
          // NOT YET IMPLEMENTED
          //
        }
        std::cout << std::endl;
        delete[] description;
      }
      TString filename_temp = fFormProdCycleFilename+period;
      SetUpForPeriodPass(period, pass, localpath, file, filename_temp, fFormProdCycleFiletype, type);
      command = wget(localpath, filename_temp, fFormProdCycleFiletype, remotepath);
      break;
    }
    case kProdCycleCpasses:{
      remotepath = TString::Format("https://alimonitor.cern.ch/prod/?t=3&res_path=mif");
      SetUpForPeriodPass(period, pass, localpath, file, fFormProdPassesFilename, fFormProdPassesFiletype, type);
      command = wget(localpath, fFormProdPassesFilename, fFormProdPassesFiletype, remotepath);
      break;
    }
    case kLogbook:
    {
      // sanity checks
      if (CheckPeriod(period) == kFALSE) return kFALSE;
      // end of sanity checks
      remotepath = TString::Format("http://aliqamod.web.cern.ch/aliqamod/data/%s/%s/%s%s", year.Data(), period.Data(), fFormLogbookFilename.Data(), fFormLogbookFiletype.Data());
      SetUpForPeriodPass(period, pass, localpath, file, fFormLogbookFilename, fFormLogbookFiletype, type);

      break;
    }
    case kTriggerClasses:
    {
      remotepath = TString::Format("http://aliqamod.web.cern.ch/aliqamod/data/%s/%s%s", period.Data(), fFormTriggerClassesFilename.Data(), fFormLogbookFiletype.Data()); // ATTENTION 'period' is year!
      SetUpForPeriodPass(period, pass, localpath, file, fFormTriggerClassesFilename, fFormTriggerClassesFiletype, type);
      break;
    }
    case kTPC:
    {
      if (CheckPeriod(period) == kFALSE) return kFALSE;
      // remotepath = TString::Format("http://aliqa<det>.web.cern.ch/aliqa<det>/data/%s/%s/%s/%s%s", year.Data(), period.Data(), pass.Data(), fFormTrendingFilename.Data(), fFormTrendingFiletype.Data()); // ATTENTION 'period' is year!
      remotepath = TString::Format("http://aliqatpc.web.cern.ch/aliqatpc/data/%s/%s/%s/%s%s", year.Data(), period.Data(), pass.Data(), fFormTrendingFilename.Data(), fFormTrendingFiletype.Data()); // ATTENTION 'period' is year!
      SetUpForPeriodPass(period, pass, localpath, file, fFormTrendingFilename, fFormTrendingFiletype, type);
      break;
    }
    case kRawQATPC:
    {
      if (CheckPeriod(period) == kFALSE) return kFALSE;
      // remotepath = TString::Format("http://aliqa<det>.web.cern.ch/aliqa<det>/data/%s/%s/%s/%s%s", year.Data(), period.Data(), pass.Data(), fFormTrendingFilename.Data(), fFormTrendingFiletype.Data()); // ATTENTION 'period' is year!
      remotepath = TString::Format("http://aliqatpc.web.cern.ch/aliqatpc/data/%s/%s%s", year.Data(), "OCDBscan", fFormTrendingFiletype.Data()); // ATTENTION 'period' is year!
      TString filename("OCDBscan");
      SetUpForPeriodPass(period, pass, localpath, file, filename, fFormTrendingFiletype, type);
      break;
    }
    case kTOF:
    {
      if (CheckPeriod(period) == kFALSE) return kFALSE;
      // remotepath = TString::Format("http://aliqa<det>.web.cern.ch/aliqa<det>/data/%s/%s/%s/%s%s", year.Data(), period.Data(), pass.Data(), fFormTrendingFilename.Data(), fFormTrendingFiletype.Data()); // ATTENTION 'period' is year!
      remotepath = TString::Format("http://aliqatof.web.cern.ch/aliqatof/data/%s/%s/%s/%s%s", year.Data(), period.Data(), pass.Data(), fFormTrendingFilename.Data(), fFormTrendingFiletype.Data()); // ATTENTION 'period' is year!
      SetUpForPeriodPass(period, pass, localpath, file, fFormTrendingFilename, fFormTrendingFiletype, type);
      break;
    }
    case kEVS:
    {
      if (CheckPeriod(period) == kFALSE) return kFALSE;
      remotepath = TString::Format("http://aliqaevs.web.cern.ch/aliqaevs/data/%s/%s/%s/%s%s", year.Data(), period.Data(), pass.Data(), fFormTrendingFilename.Data(), fFormTrendingFiletype.Data()); // ATTENTION 'period' is year!
      SetUpForPeriodPass(period, pass, localpath, file, fFormTrendingFilename, fFormTrendingFiletype, type);
      break;
    }
  }

  // Check if there is a file already
  std::cout << "-- Download from: " << remotepath << std::endl;
  if (IsDownloadNeeded(file) == kFALSE) return kFALSE;
  std::cout << "-- Save in: " << file << std::endl;

  // Create folder structure
  std::cout << "-- mkdir " << localpath << std::endl;
  gSystem->mkdir(localpath.Data(), kTRUE);


  // download information
  if ( type == kMC || type == kRCT || type == kProdCycle || type == kProdCycleCpasses ){
    std::cout << command << std::endl;
    gSystem->Exec(command.Data());
  }
  else if ((type == kLogbook) || (type == kTriggerClasses || (type == kTPC) || type == kTOF || type == kEVS || kRawQATPC)){
    std::cout << "-- TWebFile webfile(" << remotepath << ")" << std::endl;
    TWebFile webfile(remotepath);
    if (webfile.Cp(file)) std::cout << "Caching successful" << std::endl;
    else std::cout << "Caching FAILED" << std::endl;
  }
  return kTRUE;
}

Bool_t AliExternalInfo::CacheRCT(TString period, TString pass){
  return Cache(period, pass, kRCT);
}
Bool_t AliExternalInfo::CacheMC(){
  return Cache(static_cast<TString>(""), static_cast<TString>("passMC"), kMC);
}
Bool_t AliExternalInfo::CacheProdCycle(){
  return Cache("", "", kProdCycle); // pay attention, first parameter in "AliExternalInfo::Cache()" is normally period
}
Bool_t AliExternalInfo::CacheProdPasses(){
  return Cache("", "", kProdCycleCpasses); // pay attention, first parameter in "AliExternalInfo::Cache()" is normally period
}
Bool_t AliExternalInfo::CacheProdCycle(TString id){
  return Cache(id, "", kProdCycle); // pay attention, first parameter in "AliExternalInfo::Cache()" is normally period
}
// Bool_t AliExternalInfo::CacheProdCycle(TString period, TString pass){
//   return Cache(period, pass, kProdCycle); // pay attention, first parameter in "AliExternalInfo::Cache()" is normally period
// }
Bool_t AliExternalInfo::CacheLogbook(TString period){
  return Cache(period, "", kLogbook);
}
Bool_t AliExternalInfo::CacheTriggerClasses(TString year){
  return Cache(year, "", kTriggerClasses); // pay attention, first parameter in "AliExternalInfo::Cache()" is normally period
}
Bool_t AliExternalInfo::CacheDataQA(TString detector, TString period, TString pass){
  if      ( detector == "TPC" ) return Cache(period, pass, kTPC);
  else if ( detector == "RawQATPC" ) return Cache(period, pass, kRawQATPC);
  else if ( detector == "TOF" ) return Cache(period, pass, kTOF);
  else if ( detector == "EVS" ) return Cache(period, pass, kEVS);
  else return 1;
}
// Bool_t AliExternalInfo::CacheSimQA(TString detector, TString period){
//   // return Cache(period, pass, kTPC); // pay attention, first parameter in "AliExternalInfo::Cache()" is normally period
// }

// Bool_t AliExternalInfo::CacheAll(){
//   // CacheMC();
//   // @TODO Think about a way to see which periods are existing!!!
//   // for "all periods and passes": CacheRCT(...);
//   return kTRUE;
// }


// Generic function to get the tree according to parameters

TTree* AliExternalInfo::GetTree(TString period, TString pass, TString anchorYear, TString productionTag, Int_t type){
  std::cout << "GetTree()" << std::endl;
  TTree* tree = 0x0;
  char delimiter = 'n';
  TString localpath("");  // complete path to the directory of the file
  TString file("");       // complete path to file
  TString filename("");   // filename without type
  TString filetype("");   // filetype of the stored file
  TString name("");       // name of the tree
  TString title("");      // title of the tree
  switch (type){
    case kRCT:
    {
      // sanity checks
      if ( (CheckPeriod(period) == kFALSE) || (CheckPass(pass) == kFALSE) ) return tree;
      // end of sanity checks
      delimiter = '\"';
      name = "RCT";
      title = TString::Format("%s_%s", period.Data(), pass.Data()); // "period_pass"
      SetUpForPeriodPass(period, pass, localpath, file, fFormRCTFilename, fFormMCFiletype, type);
      // std::cout << period << " " << pass << " " << file << " " << title << " " << localpath << " " << type << std::endl;
      break;
    }
    case kMC:
    {
      anchorYear = ""; // Set to not create a warning, implemented later
      productionTag = ""; // Set to not create a warning, implemented later
      delimiter = '\"';
      name = "MC";
      title = "passMC";
      localpath = fFormMCDir;
      localpath.Prepend(fGlobalDir);
      file=localpath+fFormMCFilename+fFormMCFiletype;
      break;
    }
    case kProdCycle:
    {
      delimiter = '\"';
      name = "ProdCycle";
      title = "ProdCycle";
      localpath = fFormProdCycleDir;
      localpath.Prepend(fGlobalDir);
      file=localpath+fFormProdCycleFilename+period+fFormProdCycleFiletype;
      break;
    }
    case kProdCycleCpasses:
    {
      delimiter = '\"';
      name = "ProdPasses";
      title = "ProdPasses";
      localpath = fFormProdPassesDir;
      localpath.Prepend(fGlobalDir);
      file=localpath+fFormProdPassesFilename+fFormProdPassesFiletype;
      break;
    }
    case kLogbook:
    {
      // sanity checks
      if ( CheckPeriod(period) == kFALSE) return tree;
      // end of sanity checks
      name = "logbook";
      title = TString::Format("%s", period.Data()); // "period_pass"
      SetUpForPeriodPass(period, pass, localpath, file, fFormLogbookFilename, fFormLogbookFiletype, type);
      break;
    }
    case kTriggerClasses:
    {
      name = "trigger_classes";
      title = TString::Format("%s", period.Data()); // "year"
      SetUpForPeriodPass(period, pass, localpath, file, fFormTriggerClassesFilename, fFormTriggerClassesFiletype, type);
      break;
    }
    case kTPC:
    {
      name = "tpcQA";
      title = TString::Format("%s/%s", period.Data(), pass.Data()); // "year"
      SetUpForPeriodPass(period, pass, localpath, file, fFormTrendingFilename, fFormTrendingFiletype, type);
      break;
    }
    case kRawQATPC:
    {
      name = "dcs";
      title = TString::Format("%s/%s", period.Data(), pass.Data()); // "year"
      TString filename("OCDBscan");
      SetUpForPeriodPass(period, pass, localpath, file, filename, fFormTrendingFiletype, type);
      break;
    }
    case kTOF:
    {
      name = "trending";
      title = TString::Format("%s/%s", period.Data(), pass.Data()); // "year"
      SetUpForPeriodPass(period, pass, localpath, file, fFormTrendingFilename, fFormTrendingFiletype, type);
      break;
    }
    case kEVS:
    {
      name = "trending";
      title = TString::Format("%s/%s", period.Data(), pass.Data()); // "year"
      SetUpForPeriodPass(period, pass, localpath, file, fFormTrendingFilename, fFormTrendingFiletype, type);
      break;
    }
  }
  // Tree is actually set
  tree = new TTree(name, title);

  std::cout << "-- Read tree in " << file << std::endl;

  // Checking if tree is already cached, if not cache it automatically according to type
  if (gSystem->AccessPathName(file.Data()) == kTRUE) {
    std::cout << "---- ERROR: File not found, maybe not cached yet\n";
    std::cout << "------ Automatically Caching" << std::endl;
    if ( type == kRCT){
      Cache(period, pass, kRCT);
    }
    else if (type == kMC){
      Cache(static_cast<TString>(""), static_cast<TString>("passMC"), kMC);
    }
    else if (type == kProdCycle){
      Cache(period, "", kProdCycle);
    }
    else if (type == kProdCycleCpasses){
      Cache(period, "", kProdCycleCpasses);
    }
    else if (type == kLogbook){
      Cache(period, "", kLogbook);
    }
    else if (type == kTriggerClasses){
      Cache(period, "", kTriggerClasses);
    }
    else if (type == kTPC) {
      Cache(period, pass, kTPC);
    }
    else if (type == kRawQATPC) {
      Cache(period, pass, kRawQATPC);
    }
    else if (type == kTOF) {
      Cache(period, pass, kTOF);
    }
    else if (type == kEVS) {
      Cache(period, pass, kEVS);
    }
  }

  if ( type == kMC || type == kRCT || type == kProdCycle || type == kProdCycleCpasses){
    // tree->ReadFile(file, "", delimiter);
    if ( (tree->ReadFile(file, "", delimiter)) > 0) std::cout << "-- Successfully read in tree" << std::endl;
    else std::cout << "-- Error while reading tree" << std::endl;
    if (file.Contains(".mif")) {
      file.Chop();file.Chop();file.Chop();
      file.Append("root");
    }
    TFile tempfile(file, "RECREATE");
    tempfile.cd();
    tree->Write();
    tempfile.Close();
    std::cout << "Write tree to file: " << file << std::endl;
  }
  else if (type == kLogbook || type == kTriggerClasses || type == kTPC || type == kRawQATPC || type == kTOF || type == kEVS){
    TFile* treefile = new TFile(file.Data());
    tree = dynamic_cast<TTree*>( treefile->Get(name));
    if (tree==NULL) tree= dynamic_cast<TTree*>( treefile->Get("trending")); // MI temporary FIX  - code has to be reviewed and rewritten
    if (tree != 0x0) {std::cout << "-- Successfully read in tree" << std::endl;}
    else std::cout << "-- Error while reading tree" << std::endl;
  }
  return tree;
}

// ##################
// Returns a pointer to the MC tree*
TTree* AliExternalInfo::GetTreeMC(TString period, TString anchorYear, TString productionTag){
  return GetTree(period, "passMC", anchorYear, productionTag, kMC);
}

// ##################
// Returns a pointer to the RCT tree
TTree* AliExternalInfo::GetTreeRCT(TString period, TString pass){
  return GetTree(period, pass, "", "", kRCT);
}

// ##################
// Returns a pointer to the Production Cycle tree
TTree* AliExternalInfo::GetTreeProdCycle(){
  return GetTree("", "", "", "", kProdCycle);
}

// ##################
// Returns a pointer to the Production Cycle tree
TTree* AliExternalInfo::GetTreeProdPasses(){
  return GetTree("", "", "", "", kProdCycleCpasses);
}

// ##################
// Returns a pointer to the Production Cycle tree
TTree* AliExternalInfo::GetTreeProdCycle(TString id){
  return GetTree(id, "", "", "", kProdCycle);
}
// ##################
// Returns a pointer to the Production Cycle tree
// TTree* AliExternalInfo::GetTreeProdCycle(TString period, TString pass){
//   return GetTree(period, pass, "", "", kProdCycle);
// }

// ##################
// Returns a pointer to the Logbook tree
TTree* AliExternalInfo::GetTreeLogbook(TString period){
  return GetTree(period, "", "", "", kLogbook);
}

// ##################
// Returns a pointer to the Trigger Classes tree
TTree* AliExternalInfo::GetTreeTriggerClasses(TString year){
  return GetTree(year, "", "", "", kTriggerClasses);
}

TTree* AliExternalInfo::GetTreeDataQA(TString detector, TString period, TString pass){
  if      (detector == "TPC" ) return GetTree(period, pass, "", "", kTPC);
  else if (detector == "RawQATPC" ) return GetTree(period, pass, "", "", kRawQATPC);
  else if (detector == "TOF" ) return GetTree(period, pass, "", "", kTOF);
  else if (detector == "EVS" ) return GetTree(period, pass, "", "", kEVS);
  else{
    return new TTree;
  }
}

TChain* AliExternalInfo::GetChain(TString period, TString pass, TString anchorYear, TString productionTag, Int_t type)
{
  //blabla file name
  TString year = GetYearFromPeriod(period);
  TString name("");
  TString localpath("");
  TString file("");
  switch (type){
    case kRCT:
    {
      name = "RCT";
      SetUpForPeriodPass(period, pass, localpath, file, fFormRCTFilename, ".root", type);
      // std::cout << period << " " << pass << " " << file << " " << title << " " << localpath << " " << type << std::endl;
      break;
    }
    case kMC:
    {
      anchorYear = ""; // Set to not create a warning, implemented later
      productionTag = ""; // Set to not create a warning, implemented later
      name = "MC";
      localpath = fFormMCDir;
      localpath.Prepend(fGlobalDir);
      file=localpath+fFormMCFilename+fFormMCFiletype;
      break;
    }
    case kProdCycle:
    {
      name = "ProdCycle";
      localpath = fFormProdCycleDir;
      localpath.Prepend(fGlobalDir);
      file=localpath+fFormProdCycleFilename+period+".root";
      break;
    }
    case kProdCycleCpasses:
    {
      name = "ProdPasses";
      localpath = fFormProdPassesDir;
      localpath.Prepend(fGlobalDir);
      file=localpath+fFormProdPassesFilename+fFormProdPassesFiletype;
      break;
    }
    case kLogbook:
    {
      name = "logbook";
      SetUpForPeriodPass(period, pass, localpath, file, fFormLogbookFilename, fFormLogbookFiletype, type);
      break;
    }
    case kTriggerClasses:
    {
      name = "trigger_classes";
      SetUpForPeriodPass(period, pass, localpath, file, fFormTriggerClassesFilename, fFormTriggerClassesFiletype, type);
      break;
    }
    case kTPC:
    {
      name = "tpcQA";
      SetUpForPeriodPass(period, pass, localpath, file, fFormTrendingFilename, fFormTrendingFiletype, type);
      break;
    }
    case kRawQATPC:
    {
      name = "dcs";
      TString filename("OCDBscan");
      SetUpForPeriodPass(period, pass, localpath, file, filename, fFormTrendingFiletype, type);
      break;
    }
    case kTOF:
    {
      name = "trending";
      SetUpForPeriodPass(period, pass, localpath, file, fFormTrendingFilename, fFormTrendingFiletype, type);
      break;
    }
    case kEVS:
    {
      name = "trending";
      // TString global_properties = "global_properties"; // comment in these two lines if you want to have the global properties file with the interaction rate
      // SetUpForPeriodPass(period, pass, localpath, file, global_properties, fFormTrendingFiletype, type);
      SetUpForPeriodPass(period, pass, localpath, file, fFormTrendingFilename, fFormTrendingFiletype, type);
      break;
    }
  }

  TString cmd = TString::Format("ls %s", file.Data());
  std::cout << "==== cmd: " << cmd << std::endl;


  // TString fileName=fFormTrendingFilename;
  // TString cmd=TString::Format("ls %s/data/*/%s/%s/%s_TPC.root", fGlobalDir.Data(), period.Data(), pass.Data(), fileName.Data());
  TString files=gSystem->GetFromPipe(cmd.Data());
  TObjArray *arrFiles=files.Tokenize("\n");
  std::cout << "Files to add to chain: " << files << std::endl;

  //function to get tree namee based on type
  TChain *chain=new TChain(name.Data());
  for (Int_t ifile=0; ifile<arrFiles->GetEntriesFast(); ++ifile) {
    chain->AddFile(arrFiles->At(ifile)->GetName());
  }

  delete arrFiles;
  return chain;
}

// ##################
// Returns a pointer to the MC tree*
TChain* AliExternalInfo::GetChainMC(TString period, TString anchorYear, TString productionTag){
  return GetChain(period, "passMC", anchorYear, productionTag, kMC);
}

// ##################
// Returns a pointer to the RCT tree
TChain* AliExternalInfo::GetChainRCT(TString period, TString pass){
  return GetChain(period, pass, "", "", kRCT);
}

// // ##################
// // Returns a pointer to the Production Cycle tree
// TChain* AliExternalInfo::GetChainProdCycle(){
//   return GetChain("", "", "", "", kProdCycle);
// }

// // ##################
// // Returns a pointer to the Production Cycle tree
// TChain* AliExternalInfo::GetChainProdPasses(){
//   return GetChain("", "", "", "", kProdCycleCpasses);
// }

// ##################
// Returns a pointer to the Production Cycle tree
TChain* AliExternalInfo::GetChainProdCycle(TString id){
  return GetChain(id, "", "", "", kProdCycle);
}
// ##################
// Returns a pointer to the Production Cycle tree
// TChain* AliExternalInfo::GetChainProdCycle(TString period, TString pass){
//   return GetChain(period, pass, "", "", kProdCycle);
// }

// ##################
// Returns a pointer to the Logbook tree
TChain* AliExternalInfo::GetChainLogbook(TString period){
  return GetChain(period, "", "", "", kLogbook);
}

// ##################
// Returns a pointer to the Trigger Classes tree
TChain* AliExternalInfo::GetChainTriggerClasses(TString year){
  return GetChain(year, "", "", "", kTriggerClasses);
}

TChain* AliExternalInfo::GetChainDataQA(TString detector, TString period, TString pass){
  if      (detector == "TPC" ) return GetChain(period, pass, "", "", kTPC);
  else if (detector == "RawQATPC" ) return GetChain(period, pass, "", "", kRawQATPC);
  else if (detector == "TOF" ) return GetChain(period, pass, "", "", kTOF);
  else if (detector == "EVS" ) return GetChain(period, pass, "", "", kEVS);
  else{
    return new TChain;
  }
}

// ##################
// Gives true if Download is needed, because the file is not found or it is older than the set timelimit
Bool_t AliExternalInfo::IsDownloadNeeded(TString file){
  std::cout << "-- Check, if " << file << " is already there" << std::endl;
  if (gSystem->AccessPathName(file.Data()) == kTRUE) {
    std::cout << "---- File not found locally --> Caching from remote" << std::endl;
    return kTRUE;
  }
  else {
    std::cout << "---- File already downloaded --> Check if older than timelimit" << std::endl;
    struct stat st;
    stat(file.Data(), &st);
    std::time_t timeNow = std::time(0);
    long int timeFileModified = st.st_mtime;
    std::cout << "------ File is " << timeNow - timeFileModified << " seconds old" << std::endl;
    if (timeNow - timeFileModified < fTimeLimit) {
      std::cout << "-------- File is not older than the set timelimit " << fTimeLimit << " seconds --> no caching needed" << std::endl;
      return kFALSE; // if file is younger than the set time limit, it will not be downloaded again
    }
    else {
      std::cout << "-------- File is older than the set timelimit " << fTimeLimit << " seconds --> Caching from remote" << std::endl;
      return kTRUE;
    }
  }
}


// Setting up the correct paths and filenames
void AliExternalInfo::SetUpForPeriodPass(const TString& period, const TString& pass, TString& localpath, TString& file, const TString& filename, const TString& filetype, Int_t type){
  TString year = GetYearFromPeriod(period);
  if (type == kRCT){
    localpath = TString::Format(fFormRCTDir.Data(), year.Data(), period.Data(), pass.Data());
  }
  else if (type == kMC){
    localpath = fFormMCDir;
  }
  else if (type == kProdCycle){
    // filename.Append(period);
    localpath = fFormProdCycleDir;
  }
  else if (type == kProdCycleCpasses){
    // filename.Append(period);
    localpath = fFormProdPassesDir;
  }
  else if (type == kLogbook){
    localpath = TString::Format(fFormLogbookDir.Data(), year.Data(), period.Data(), pass.Data());
  }
  else if (type == kTriggerClasses){
    localpath = TString::Format(fFormTriggerClassesDir.Data(), period.Data());
  }
  else if (type == kTPC){
    localpath = TString::Format(fFormTrendingDir.Data(), year.Data(), period.Data(), pass.Data());
    localpath.Prepend(fGlobalDir);
    file = localpath+filename+"_TPC"+filetype;
    return;
  }
  else if (type == kRawQATPC){
    localpath = TString::Format("/data/%s/", year.Data());
    localpath.Prepend(fGlobalDir);
    file = localpath+filename+""+filetype;
    return;
  }
    else if (type == kTOF){
    localpath = TString::Format(fFormTrendingDir.Data(), year.Data(), period.Data(), pass.Data());
    localpath.Prepend(fGlobalDir);
    file = localpath+filename+"_TOF"+filetype;
    return;
  }
  else if (type == kEVS){
    localpath = TString::Format(fFormTrendingDir.Data(), year.Data(), period.Data(), pass.Data());
    localpath.Prepend(fGlobalDir);
    file = localpath+filename+"_EVS"+filetype;
    return;
  }
  localpath.Prepend(fGlobalDir);
  file = localpath+filename+filetype;
}


// ##################
// Checks if path is correct
Bool_t AliExternalInfo::CheckPeriod(TString period){
  if (period.Contains("LHC") == kFALSE) {
    std::cout << "ERROR: Period is not given correctly; Should be something like \"LHC10b\"; Input was: " << period << std::endl;
    return kFALSE;
  }
  if (period.Length() != 6) {
    std::cout << "ERROR: Period name is to long/short, should be 6 characters like \"LHC10b\"; Input was: " << period.Length() << std::endl;
    return kFALSE;
  }
  return kTRUE;
}


Bool_t AliExternalInfo::CheckPass(TString pass) {
  if (pass.Contains("pass") == kFALSE) {
    std::cout << "ERROR: Pass is not given correctly; Should be something like \"pass2\"; Input was: " << pass << std::endl;
    return kFALSE;
  }
  if (pass.Length() != 5) {
    std::cout << "ERROR: Pass name is to long/short, should be 5 characters like \"pass2\"; Input was: " << pass.Length() << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

TString AliExternalInfo::wget(const TString& localpath, const TString& filename, const TString& filetype, const TString& remotepath){
  return TString::Format("wget --no-check-certificate --secure-protocol=TLSv1 --certificate=%s --private-key=%s -o %s%s.log -O %s%s%s \"%s\"",
                                     fCertificate.Data(), fPrivateKey.Data(),
                                     localpath.Data(), filename.Data(),
                                     localpath.Data(), filename.Data(), filetype.Data(),
                                     remotepath.Data());;
}

void AliExternalInfo::CreateMapWithPassesAndPeriods(){
  // Creates a vector containing objects of the class "Period" which itself contains the different passes in a vector of TString
  TTree* tree = GetTreeProdCycle();
  std::vector<Period> vec;
  // Reading of the period names in the tree
  char* temp_tag = new char[200]();
  tree->SetBranchAddress("Tag", temp_tag);
  Int_t nEntries = tree->GetEntriesFast();

  // Looping over the tree and chopping the period_pass name into pieces and storing it in a vector
  for (Int_t i = 0; i < nEntries; ++i){
    tree->GetEntry(i);

    // std::cout << "Tag: " << temp_tag;
    TString tag(temp_tag);
    // std::cout << "    Contains \"_\" = " << tag.First('_') << "\n";

    TObjArray objArray;
    int first_occurence_ = tag.First('_');
    if (first_occurence_ != -1){
      TString period_name(tag(0,first_occurence_));
      TString pass_name(tag(first_occurence_+1, tag.Length()-first_occurence_));

      if (period_name.Length() != 6) {continue;} // If period does not look like "LHC10b" pattern
      // std::cout << "Token0 = " << period_name << "   Token1 = " << pass_name << "    Token0.Lenght = " << period_name.Length() << "\n";

      // std::cout << "-- SIZE 6" << std::endl;
      int found = 0;
      unsigned int j = 0;
      for (; j < vec.size(); ++j){
        if (vec.at(j).GetPeriod() == period_name){
          found = 1;
          break;
        }
      }
      if (found == 0){ // if period is already found
        Period p(period_name);
        p.AddPass(pass_name);

        vec.push_back(p);
      }
      else{
        vec.at(j).AddPass(pass_name);
      } // if '_' in tag
    } // for loop
  }
  std::cout << std::endl; // To flush the buffer
  delete[] temp_tag;
  for (unsigned int i = 0; i < vec.size(); ++i){
    vec.at(i).PrintPeriodPass();
  }
  return;
}

TString AliExternalInfo::GetYearFromPeriod(TString period){
  TString year(period(3,2));
  year = TString::Format("20%s", year.Data());
  return year;
}
