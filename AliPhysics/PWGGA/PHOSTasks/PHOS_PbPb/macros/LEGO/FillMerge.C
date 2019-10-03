#include "TFileMerger.h"
#include "TSystem.h"

#include <map>
#include <iostream>
#include <fstream>
#include <string>

std::map<int, int> runToFill; // maps: run -> fill
std::map<int, TFileMerger*> fillToMerger;

void makeFillMap()
{
}

void FillMerge(const TString runFile="runFile.txt",
	       const TString runFillFile = "run_fill.txt")
{
  // Use to merge run files to fill files
  // For some reason, it does not work without compiling:
  // .L FillMerge.C+g
  // FillMerge()
  //
  // run_fill.txt must be a two column txt file, 
  // column 1 being run and column 2 being fill.

  //TGrid::Connect("alien://");
  

  // Make runToFill map
  std::ifstream rfin;
  rfin.open(runFillFile.Data());
  
  while (1) {
    int run=0, fill=0;
    rfin >> run >> fill;
    if (!rfin.good()) break;
    runToFill[run] = fill;
  }


  // Make FillMergers and add run files to them.
  std::ifstream in;
  in.open(runFile.Data());
  char rootFileName[256];
  Int_t runNumber=0;
  while (1) {
    in >> runNumber >> rootFileName;
    if (!in.good()) break;

    if(169553 == runNumber ) //169094 
      continue;
    
    if( ! runToFill.count(runNumber) ) {
      Printf("can't map run %d to fill, not in file %s", runNumber, runFillFile.Data());
      continue;
    }
    
    int fill = runToFill[runNumber];
    
    if( ! fillToMerger.count( fill ) ) { // if no merger for fill
      fillToMerger[fill] = new TFileMerger(); // create merger
      gSystem->mkdir( Form("%s/fill/%d", gSystem->pwd(), fill), true );
      fillToMerger[fill]->OutputFile(Form("%s/fill/%d/AnalysisResults.root", gSystem->pwd(), fill));
    }
    
    fillToMerger[fill]->AddFile(rootFileName);
  }

  std::ofstream out("fillFile.txt");

  
  // Merge:
  std::map<int, TFileMerger*>::iterator it;
  for( it = fillToMerger.begin(); it != fillToMerger.end(); ++it ) {
    TFileMerger* fm = (*it).second;
    Printf("merging %s", fm->GetOutputFileName());
    fm->Merge();
    out << (*it).first << " " << fm->GetOutputFileName() << std::endl;
    delete fm;
  }

}
