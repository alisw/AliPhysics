#include "TFileMerger.h"
#include "TSystem.h"
#include <TString.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>

void MergeSet(const TString mergeFileName,
	      const TString runSetName = "runSet.txt",
	      const TString runFileName = "runFile.txt")
{
  // Author: Henrik.Qvigstad@cern.ch, 2012
  //
  // Use to merge sets of runs defined by file by name of 'runSetName'
  //
  // For some reason, it does not work without compiling:
  //
  // .L MergeSet.C+g
  // MergeSet()

  // Make runToFill map
  std::ifstream runFile;
  runFile.open(runFileName.Data());
  
  std::map<int, std::string> runToFile; // maps: run -> fill
  while (1) {
    int run=0;
    std::string file;
    runFile >> run >> file;
    if (!runFile.good()) break;
    runToFile[run] = file;
  }


  // Make FillMergers and add run files to them.
  std::ifstream runSet;
  TFileMerger* fileMerger = new TFileMerger(); // create merger
  fileMerger->OutputFile(mergeFileName.Data());
  runSet.open(runSetName.Data());
  Int_t runNumber=0;
  while (1) {
    runSet >> runNumber;
    if (!runSet.good()) break;

//     if(169553 == runNumber ) //169094 
//       continue;
    
    if( ! runToFile.count(runNumber) ) {
      Printf("can't map run %d to file", runNumber);
      continue;
    }
   
    std::string file = runToFile[runNumber];
    printf("adding %s: ", file.data());
    fileMerger->AddFile(file.data());
  }
  
  Printf("Merging %s...", mergeFileName.Data());
  fileMerger->Merge();
  Printf("Merging %s... done! \n", mergeFileName.Data());
  delete fileMerger;
}


void MergeTpcSets()
{
  MergeSet("AnalysisResults_GoodTpc.root", "../runsGoodTpcLHC11h.txt", "runFile.txt");
  MergeSet("AnalysisResults_SemiGoodTpc.root", "../runsSemiGoodTpcLHC11h.txt", "runFile.txt");
  MergeSet("AnalysisResults_BadTpc.root", "../runsBadTpcLHC11h.txt", "runFile.txt");
}
