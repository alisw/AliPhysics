#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFileMerger.h"
#include "TGrid.h"
#include "TList.h"
#include "TSystem.h"
#include "TGridResult.h"
#include "TMap.h"
#include "TFile.h"
#include "TROOT.h"
#endif

Bool_t AddFileList(TString, TString&, TString);
void ReadListFromFile(TString, TString&, TString);
Int_t GetLastStage(TGridResult*);

void mergeGridFiles(TString outFilename, TString inFileList, TString addPrefix = "alien://", Int_t nFilesPerStep = 100, Bool_t copyLocal = kTRUE, TString dirsToMerge = "")
{
  TString fileList = "";
  ReadListFromFile(inFileList, fileList, addPrefix);

  if ( fileList.IsNull() ) {
    printf("List of merging files is null: Nothing done!\n");
    return;
  }

  TString logFilename = "toBeMerged.txt";

  if ( ! gGrid ) copyLocal = kFALSE;

  TString currList = fileList;
  TString currOutput = "", mergedFiles = "";
  TFileMerger* fileMerger = 0x0;
  TObjArray* dirsToMergeArray = 0x0;
  TList* sourceList = new TList();
  if ( ! dirsToMerge.IsNull() ) {
    dirsToMergeArray = dirsToMerge.Tokenize(" ");
    dirsToMergeArray->SetOwner();
  }

  Bool_t showProgressBar = ! gROOT->IsBatch();

  for ( Int_t istep = 0; istep < 100; istep++ ) {
    TObjArray* array = currList.Tokenize(" ");
    array->SetOwner();
    currList = "";
    Int_t nFiles = array->GetEntries();
    Int_t subStep = -1;
    for (Int_t ifile = 0; ifile < nFiles; ifile++ ) {
      if ( ! fileMerger ) fileMerger = new TFileMerger(copyLocal);
      TString currFilename = array->At(ifile)->GetName();
      Bool_t isFileAdded = fileMerger->AddFile(currFilename.Data(), showProgressBar);
      if ( ! isFileAdded ) continue;
      //printf("New file %s\n", gROOT->GetListOfFiles()->Last()->GetName()); // REMEMBER TO CUT
      if ( dirsToMergeArray ) {
        if ( ! sourceList ) sourceList = new TList();
        sourceList->Add(gROOT->GetListOfFiles()->Last());
      }
      
      Int_t nFilesToMerge = fileMerger->GetMergeList()->GetEntries();
      if ( nFilesToMerge % nFilesPerStep != 0 && ifile < nFiles - 1 ) 
	continue;
      // The following part is executed only at the end of each step
      currOutput = outFilename;
      if ( nFiles > nFilesPerStep ) {
	subStep++;
	currOutput.ReplaceAll(".root",Form("_%i_%i.root", istep, subStep));
	AddFileList(currOutput, currList, "");
      }
      if ( dirsToMergeArray ) {
        TFile* outFile = new TFile(currOutput.Data(), "recreate");
        for ( Int_t idir=0; idir<dirsToMergeArray->GetEntries(); idir++ ) {
          outFile->cd();
          TDirectory* currDir = outFile->mkdir(dirsToMergeArray->At(idir)->GetName());
          fileMerger->MergeRecursive(currDir, sourceList);
        }
        outFile->Close();
      }
      else {
        fileMerger->OutputFile(currOutput.Data());
        fileMerger->Merge();
      }
      printf("\nMerged in %s:\n", currOutput.Data());
      mergedFiles = "";
      for ( Int_t ientry=0; ientry<nFilesToMerge; ientry++ )
	mergedFiles += Form("%s ", fileMerger->GetMergeList()->At(ientry)->GetName());
      printf("%s\n\n", mergedFiles.Data());

      // remove merged files
      if ( istep > 0 )
	gSystem->Exec(Form("rm %s", mergedFiles.Data()));

      delete fileMerger;
      fileMerger = 0x0;
      if ( dirsToMergeArray ) {
        delete sourceList;
        sourceList = 0x0;
      }

      // Write log file to keep trace of files to be merged
      ofstream logFile(logFilename.Data());
      TString logString = "";
      for ( Int_t jfile = ifile + 1; jfile < nFiles; jfile++ ) {
	logString += Form("%s ", array->At(jfile)->GetName());
      }
      logString.Append(currList.Data());
      logString.ReplaceAll(" ", "\n");
      logFile << logString.Data() << endl;;
      logFile.close();
    } // loop on files


    delete array;
    printf("Step %i completed!\n", istep);

    if ( nFiles <= nFilesPerStep ) break;
  } // loop on steps

  gSystem->Exec(Form("rm %s", logFilename.Data()));
}


//___________________________________________________
Bool_t AddFileList(TString filename, TString& fileList, TString addPrefix)
{
  if ( filename.IsNull() || ! filename.Contains(".root") ) return kFALSE;

  if ( ! addPrefix.IsNull() && ! filename.Contains(addPrefix.Data()) )
    filename.Prepend(addPrefix.Data());

  if ( filename.Contains("alien://") && ! gGrid )
    TGrid::Connect("alien://");
  
  if ( ! fileList.IsNull() )
    fileList.Append(" ");
  fileList.Append(filename.Data());

  return kTRUE;
}

//___________________________________________________
void ReadListFromFile(TString filename, TString& fileList, TString addPrefix)
{
  ifstream inFile(filename.Data());
  TString currLine = "";
  if ( inFile.is_open() ) {
    while ( ! inFile.eof() ) {
      currLine.ReadLine(inFile,kTRUE); // Read line
      AddFileList(currLine, fileList, addPrefix);
    }
    inFile.close();
  }
}


//___________________________________________________
void completeProd(TString runListName="runList.txt", TString prodDir = "", TString baseDir="/alice/data/2010/LHC10h", TString outTaskFilename="QAresults.root", Int_t nFilesPerStep = 50, TString dirsToMerge = "MUON_QA MUON.TriggerEfficiencyMap", Bool_t mergeFast = kFALSE, Bool_t overwriteExisting = kFALSE)
{
  TString outFilename = "completeFileList.txt";

  // Get run list from file
  ifstream inFile(runListName.Data());
  TObjArray runList;
  runList.SetOwner();
  TString currRun;
  Int_t minRun = 99999999;
  Int_t maxRun = -1;
  if (inFile.is_open()) {
    while (! inFile.eof() ) {
      currRun.ReadLine(inFile,kTRUE); // Read line
      if ( currRun.IsNull() || ! currRun.IsDigit() ) continue;
      Int_t currRunInt = currRun.Atoi();
      minRun = TMath::Min(currRunInt, minRun);
      maxRun = TMath::Max(currRunInt, maxRun);
      runList.Add(new TObjString(Form("%d", currRunInt)));
    }
    inFile.close();
  }

  outFilename.ReplaceAll(".txt", Form("_%i_%i.txt", minRun, maxRun));

  //TString filePattern[3] = {"", "Stage*/","*/"};
  TString filePattern[2] = {"","*/"};

  ofstream outFile(outFilename.Data());

  // if ( outTaskFilename.Contains("QAresults.root") ) {
  const Int_t kNlibs = 5; // 1
  //TString loadLibs[kNlibs] = {"libPWG3base.so"};
  //TString loadLibs[kNlibs] = {"libANALYSIS.so", "libANALYSISalice.so", "libTENDER.so", "libPWG1.so", "libPWG3base.so"};
  TString loadLibs[kNlibs] = {"libANALYSIS.so", "libOADB.so", "libANALYSISalice.so", "libCORRFW.so", "libPWG3base.so"};
  for ( Int_t ilib=0; ilib<kNlibs; ilib++ ) {
    Int_t exitVal = gSystem->Load(loadLibs[ilib].Data());
    if ( exitVal < 0 ) {
      printf("Please run with aliroot if you're merging QA objects!\n");
      return;
    }
  }
  //}

  if ( ! gGrid )
    TGrid::Connect("alien://");

  baseDir.ReplaceAll("alien://","");

  TMap* map = 0x0;
  TString stageName = "";
  TString runsWithoutOut = "";
  for ( Int_t irun=0; irun<runList.GetEntries(); irun++ ) {
    TString currRunString = ((TObjString*)runList.At(irun))->GetString();

    TString localOut = outTaskFilename;
    localOut.ReplaceAll(".root", Form("_%s.root", currRunString.Data()));
    if ( ! gSystem->AccessPathName(localOut.Data()) ) {
      if ( overwriteExisting )
        printf("Overwriting existing file %s\n", localOut.Data());
      else {
        printf("Warning: merged file %s already exist: do not overwrite\n", localOut.Data());
        outFile << gSystem->pwd() << "/" << localOut.Data() << endl;
        continue;
      }
    }

    TString tmpFilename = Form("/tmp/mergeListRun%s.txt", currRunString.Data());
    ofstream tmpFile(tmpFilename.Data());
    TString mergeFilename = "";

    Int_t nPatterns = ( mergeFast ) ? 1 : 2;
    
    for ( Int_t ipattern=0; ipattern<nPatterns; ipattern++ ) {
      TString command = ( prodDir.Contains("private") ) ? Form("find %s/ *%s/%s%s", baseDir.Data(), currRunString.Data(), filePattern[ipattern].Data(), outTaskFilename.Data()) : Form("find %s/*%s %s/%s%s", baseDir.Data(), currRunString.Data(), prodDir.Data(), filePattern[ipattern].Data(), outTaskFilename.Data());

      printf("%s\n", command.Data());

      TGridResult* res = gGrid->Command(command);

      if ( ! res || res->GetEntries() == 0 ) continue;

      Int_t mergeStage = ( ipattern == 1 ) ? GetLastStage(res)  : -1;
      stageName = Form("Stage_%i", mergeStage);

      TIter nextmap(res);
      while ( ( map = (TMap*)nextmap() ) ) {
        // Loop 'find' results and get next LFN
        TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
        if (!objs || !objs->GetString().Length()) 
          continue;

        mergeFilename = objs->GetString();

        if ( mergeStage > 0 && ! mergeFilename.Contains(stageName.Data()) ) continue;

        tmpFile << mergeFilename.Data() << endl;
      } // loop on grid lfns

      delete res;

      tmpFile.close();

      // If the merging of specific directories is required
      // run merging also if there is 1 file
      // so that we get rid of other sub-directories
      if ( ipattern == 1 || ! dirsToMerge.IsNull() ) {
        mergeFilename = outTaskFilename;
        mergeFilename.ReplaceAll(".root", Form("_%s.root", currRunString.Data()));
        mergeGridFiles(mergeFilename, tmpFilename, "alien://", nFilesPerStep, kTRUE, dirsToMerge);
      }

      if ( ! mergeFilename.Contains("alien://") )
        outFile << gSystem->pwd() << "/";
      outFile << mergeFilename.Data() << endl;
      gSystem->Exec(Form("rm %s", tmpFilename.Data()));
      break;
    } // loop on pattern
    if ( mergeFilename.IsNull() ) runsWithoutOut += currRunString + " ";
  } // loop on runs
  if ( ! runsWithoutOut.IsNull() )
    printf("\nNo output found in runs\n%s\n",runsWithoutOut.Data());
  printf("\nOutput written in:\n%s\n", outFilename.Data());
}


//___________________________________________________
Int_t GetLastStage(TGridResult* res)
{
  Int_t lastStage = 0;

  TMap* map = 0x0;
  TIter nextmap(res);
  TString filename = "", currToken = "";
  while ( ( map = (TMap*)nextmap() ) ) {
    // Loop 'find' results and get next LFN
    TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
    if (!objs || !objs->GetString().Length()) 
      continue;

    filename = objs->GetString();
    
    if ( ! filename.Contains("Stage_") ) continue;

    TObjArray* array = filename.Tokenize("/");
    array->SetOwner();
    for ( Int_t ientry=0; ientry<array->GetEntries(); ientry++ ) {
      currToken = array->At(ientry)->GetName();
      if ( currToken.Contains("Stage_") ) {
        currToken.Remove(0,6);
        Int_t currStage = currToken.Atoi();
        lastStage = TMath::Max(currStage, lastStage);
        break;
      }
    }
    delete array;
  } // loop on grid lfns

  return lastStage;
}
