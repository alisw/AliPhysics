#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TFile.h"
#include "TGrid.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TKey.h"
#include "TTree.h"
#include "TParameter.h"
#include "TFileMerger.h"
#include "THashList.h"

// Aliroot includes
#include "AliAnalysisManager.h"
#include "AliAnalysisAlien.h"
#include "AliESDInputHandler.h"
#include "AliCounterCollection.h"

#define COMPILEMACRO

#endif


//_____________________________________________________________________________
void LoadLibs()
{
//  gSystem->Load("libTree");
//  gSystem->Load("libGeom");
//  gSystem->Load("libVMC");
//  gSystem->Load("libPhysics");
//  gSystem->Load("libProof");
//
//  gSystem->Load("libANALYSIS");
//  gSystem->Load("libOADB");
//  gSystem->Load("libANALYSISalice");
//  gSystem->Load("libCORRFW");
//  gSystem->Load("libPWGmuon");
  TString libName = "libPWGPPMUONlite";
  TString getLib = gSystem->GetLibraries(libName.Data(),"",kFALSE);
  if ( getLib.IsNull() ) gSystem->Load(libName.Data());
}


//_____________________________________________________________________________
AliAnalysisAlien* CreateAlienHandler()
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Set the run mode
  plugin->SetRunMode("terminate");

  // Declare all libraries
  plugin->SetAdditionalLibs("libCORRFW.so libPWGHFbase.so libPWGmuon.so libPWGPPMUONlite.so");

  plugin->SetAdditionalRootLibs("libXMLParser.so libGui.so libProofPlayer.so");

  plugin->AddIncludePath("-I.");
  plugin->AddIncludePath("-I$ALICE_PHYSICS/PWGPP/MUON/lite");

  return plugin;
}

enum {
  trackQA = 1 << 0,
  trigQA  = 1 << 1
};

//_____________________________________________________________________________
void terminateQA ( TString outfilename = "QAresults.root", Bool_t isMC = kFALSE, Bool_t usePhysicsSelection = kTRUE, UInt_t mask = (trackQA|trigQA), UInt_t force = (trackQA|trigQA) )
{
  //
  // Run terminate on QA output
  // Terminate is skipped if it was already run during the production
  // Unless the "force" option is specified
  //

  LoadLibs();

  AliAnalysisAlien* alienHandler = CreateAlienHandler();

  AliAnalysisManager* mgr = new AliAnalysisManager("testAnalysis");
  mgr->SetCommonFileName(outfilename.Data());
  mgr->SetGridHandler(alienHandler);

  // Needed to the manager (but not used in terminate mode)
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);

  TString trigOutName = "trigChEff_ANY_Apt_allTrig.root";
  if ( ( force & trigQA ) == 0 ) {
    if ( gSystem->AccessPathName(trigOutName) == 0 ) {
      printf("Terminate already done for trigger. Skip\n");
      mask &= ~trigQA;
    }
  }
  if ( ( force & trackQA ) == 0 ) {
    TFile* file = TFile::Open(outfilename.Data());
    TKey* key = file->FindKeyAny("general2");
    if ( key ) {
      printf("Terminate already done for tracker. Skip\n");
      mask &= ~trackQA;
    }
    delete file;
  }

  if ( mask == 0 ) return;

#ifndef COMPILEMACRO

  if ( mask & trigQA ) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskMTRchamberEfficiency.C");
    AliAnalysisTaskTrigChEff* trigChEffTask = AddTaskMTRchamberEfficiency(isMC);
    TString physSelName = "PhysSelPass";
    if ( ! usePhysicsSelection ) physSelName += ",PhysSelReject";
    trigChEffTask->SetTerminateOptions(physSelName,"ANY","-5_105",Form("FORCEBATCH NoSelMatchApt FromTrg %s?%s?ANY?-5_105?NoSelMatchAptFromTrg",trigOutName.Data(),physSelName.Data()));
  }
  if ( mask & trackQA ) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskMuonQA.C");
    AliAnalysisTaskMuonQA* muonQATask = AddTaskMuonQA(usePhysicsSelection);
  }

#endif

//  // Check if terminate was already performed
//  if  ( ! force ) {
//    TObject* paramContainer = mgr->GetParamOutputs()->At(0);
//    if ( paramContainer ) {
//      TFile* file = TFile::Open(outfilename);
//      if ( file->FindObjectAny(paramContainer->GetName() ) ) {
//        printf("\nTerminate was already executed!\n");
//        printf("Nothing to be done\n");
//        file->Close();
//        return;
//      }
//      file->Close();
//    }
//  }


  if ( ! mgr->InitAnalysis()) {
    printf("Fatal: Cannot initialize analysis\n");
    return;
  }
  mgr->PrintStatus();
  mgr->StartAnalysis("grid terminate");

  if ( ! gSystem->AccessPathName("outputs_valid") ) gSystem->Exec("rm outputs_valid");
}


//_____________________________________________________________________________
TString GetFullPath ( TString filename )
{
  if ( filename.BeginsWith("alien://") ) return filename;
  TString dirName = gSystem->DirName(filename);
  TString baseName = gSystem->BaseName(filename);
  TString currDir = gSystem->pwd();
  gSystem->cd(dirName);
  TString fullDir = gSystem->pwd();
  gSystem->cd(currDir);
  TString fullPath = fullDir.Data();
  if ( ! fullDir.EndsWith("/") ) fullPath.Append("/");
  fullPath += baseName;
  return fullPath;
}

//_____________________________________________________________________________
TString GetBaseName ( TString filename )
{
  TString baseName = gSystem->BaseName(filename);
  Int_t idx = baseName.Index("#");
  if ( idx > 0 ) baseName.Remove(0,idx+1);
  return baseName;
}


//_____________________________________________________________________________
void CopyDir(TDirectory *source) {
  //copy all objects and subdirs of directory source as a subdir of the current directory
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir->mkdir(source->GetName());
  adir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write(obj->GetName(),TObject::kSingleKey);
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}


//_____________________________________________________________________________
UInt_t GetQAInfo ( const char* qaFileName, TString dirNames = "MUON_QA MTR_ChamberEffMap MUON.TrigEfficiencyMap MUON.TriggerEfficiencyMap" )
{
  LoadLibs();

  UInt_t info = 0;

  TString outFilename = GetBaseName(qaFileName);
  TString inFullPath = GetFullPath(qaFileName);
  TString outFullPath = GetFullPath(outFilename);
  if ( inFullPath == outFullPath ) {
    printf("Warning: input and output are same file!\n");
    return info;
  }

  if ( inFullPath.BeginsWith("alien") && ! gGrid ) TGrid::Connect("alien://");

  TObjArray* dirList = dirNames.Tokenize(" ");
  TFile* inFile = TFile::Open(qaFileName);
  if ( ! inFile ) {
    // This might happen when checking for QAresults_outer.root
    // when an input QAresults_barrel.root is provided
    printf("Warning: file %s cannot be opened\n",qaFileName);
    return info;
  }
  TFile* outFile = TFile::Open(outFilename,"RECREATE");
  TIter next(dirList);
  TObjString* objStr = 0x0;
  while ( (objStr=static_cast<TObjString*>(next())) ) {
    inFile->cd();
    TString currDir = objStr->String();
    TObject* obj = inFile->Get(currDir.Data());
    if ( ! obj ) continue;
    if ( currDir == "MUON_QA" ) info |= trackQA;
    else info |= trigQA;
    outFile->cd();
    CopyDir(static_cast<TDirectory*>(obj));
  }
  delete outFile;
  delete inFile;
  delete dirList;

  return info;
}


//_____________________________________________________________________________
Bool_t CheckMergedOverlap ( TString fileList )
{
  LoadLibs();
  TObjArray* arr = fileList.Tokenize(" ");
  THashList triggerList;
  triggerList.SetOwner();
  Bool_t hasOverlap = kFALSE;
  for ( Int_t iarr=0; iarr<arr->GetEntries(); iarr++ ) {
    TFile* file = TFile::Open(arr->At(iarr)->GetName());
    AliCounterCollection* eventCounters = (AliCounterCollection*)file->FindObjectAny("eventCounters");
    if ( eventCounters ) {
      TString listFromContainer = eventCounters->GetKeyWords("trigger");
      TObjArray* trigArr = listFromContainer.Tokenize(",");
      for ( Int_t itrig=0; itrig<trigArr->GetEntries(); itrig++ ) {
        TString currTrig = trigArr->At(itrig)->GetName();
        if ( triggerList.FindObject(currTrig.Data()) ) {
          if ( currTrig != "ANY" ) {
            printf("Warning: duplicated trigger %s\n", currTrig.Data());
            hasOverlap = kTRUE;
          }
        }
        else triggerList.Add(new TObjString(currTrig));
      }
      delete trigArr;
    }
    delete file;
  }
  delete arr;

  return hasOverlap;
}

//_____________________________________________________________________________
Bool_t GetMergedQAInfo ( TString fileList, TString outFilename = "QAresults.root" )
{
  LoadLibs();
  TObjArray* arr = fileList.Tokenize(" ");
  TFileMerger fm;
  fm.OutputFile(outFilename.Data());
  for ( Int_t iarr=0; iarr<arr->GetEntries(); iarr++ ) {
    fm.AddFile(arr->At(iarr)->GetName());
  }
  delete arr;
  return fm.Merge();
}

//_____________________________________________________________________________
Bool_t AddTreeVariable ( TList& parList, const char* varName, char varType, Float_t val )
{
  if ( varType == 'D' ) varType = 'F';
  TString parName = Form("%s/%c",varName,varType);
  if ( varType == 'F' ) {
    parList.Add(new TParameter<float>(parName,val));
  }
  else if ( varType == 'I' ) {
    parList.Add(new TParameter<int>(parName,(Int_t)val));
  }
  else {
    printf("Error: variable type %c not accepted", varType);
    return kFALSE;
  }
  return kTRUE;
}


//_____________________________________________________________________________
void FillTree ( TTree* tree, TList &parList )
{
  Int_t nVars = parList.GetEntries();
  TArrayI varInt(nVars);
  TArrayF varFloat(nVars);
  for ( Int_t ivar=0; ivar<nVars; ivar++ ) {
    TObject* obj = parList.At(ivar);
    TString varName = obj->GetName();
    TString branchName = varName;
    branchName.Remove(varName.Length()-2);
    if ( varName.EndsWith("F") ) {
      varFloat[ivar] = ((TParameter<float>*)obj)->GetVal();
      tree->Branch(branchName.Data(),&varFloat[ivar],varName.Data());
    }
    else if ( varName.EndsWith("I") ) {
      varInt[ivar] = (Int_t)((TParameter<int>*)obj)->GetVal();
      tree->Branch(branchName.Data(),&varInt[ivar],varName.Data());
    }
  }
  tree->Fill();
}


//_____________________________________________________________________________
void AddTrigVars ( TString filename, TList &parList )
{
  TString trigOutName = "trigChEff_ANY_Apt_allTrig.root";
  if ( gSystem->AccessPathName(trigOutName.Data()) ) trigOutName = filename;
  TFile* file = TFile::Open(filename.Data());
  TList* inList = (TList*)file->FindObjectAny("triggerChamberEff");
  TString hChNames[] = {"bendPlaneCountChamber","nonBendPlaneCountChamber","allTracksCountChamber"};
  Int_t nHistos = sizeof(hChNames)/sizeof(hChNames[0]);
  for ( Int_t ihisto=0; ihisto<nHistos; ihisto++ ) {
    TH1* histo = (TH1*)inList->FindObject(hChNames[ihisto].Data());
    for ( Int_t ibin=1; ibin<=4; ibin++ ) {
      Double_t currVal = ( histo ) ? histo->GetBinContent(ibin) : 0.;
      AddTreeVariable(parList, Form("%s%i",hChNames[ihisto].Data(),ibin),'F',currVal);
    }
  }
  delete file;
}

//_____________________________________________________________________________
void MakeTrend ( const char* qaFile, Int_t runNumber, Bool_t isMC = kFALSE, Bool_t usePhysicsSelection = kTRUE, UInt_t mask = (trackQA|trigQA) )
{
  UInt_t info = GetQAInfo(qaFile);
  if ( info == 0 ) return;

  TString baseFilename = "QAresults.root";
  TString inFilename = GetBaseName(qaFile);

  UInt_t forceTerminate = 0;
  if ( inFilename.Contains("barrel") ) {
    TString outerInFilename(qaFile);
    outerInFilename.ReplaceAll("barrel","outer");
    UInt_t outerInfo = GetQAInfo(outerInFilename);
    if ( outerInfo ) {
      // Merge outer and barrel
      TString fileList = GetBaseName(outerInFilename);
      fileList += " " + inFilename;
      Bool_t isMergedOk = GetMergedQAInfo(fileList,baseFilename);
      if ( isMergedOk ) {
        printf("Merged files: %s => %s\n",fileList.Data(),baseFilename.Data());
        CheckMergedOverlap(fileList);
        gSystem->Exec(Form("rm %s",fileList.Data())); // Remove QAresults_barrel and outer
        forceTerminate = info; // Re-do terminate when merging barrel and outer
      }
    }
  }

  if ( inFilename != baseFilename && gSystem->AccessPathName(baseFilename) ) {
    gSystem->Exec(Form("mv -v %s %s",inFilename.Data(),baseFilename.Data()));
  }

  UInt_t checkedMask = mask&info;

  terminateQA(baseFilename,isMC,usePhysicsSelection,checkedMask,forceTerminate);

  TList parList;
  parList.SetOwner();
  AddTreeVariable(parList, "run", 'I', runNumber);

  // function for trigger
  if ( checkedMask & trigQA ) AddTrigVars(baseFilename.Data(),parList);

  TFile* outFile = TFile::Open("trending.root","RECREATE");
  TTree* tree = new TTree("trending","trending");

  FillTree(tree, parList);
  tree->Write();
  delete outFile;
}
