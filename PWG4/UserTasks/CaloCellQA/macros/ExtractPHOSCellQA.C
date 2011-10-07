//------------------------------------------------------------------------
// PWG1 QA train produces QA histograms QAresults.root, where PHOS
// histograms for two events types are stored in
// TObjArray *PHOSCellsQA_AnyInt and TObjArray *PHOSCellsQA_PHI7
// As each a root file for eah run contains, by design, unique histograms
// per run, the root files from different runs cannot be merged 
// via TFileMerger.
// This drawback of the QA design is solved by extracting the PHOS
// histograms from TObjArray's of QAresults.root to separate files per run
// and per event type.
//
// Usage:
// 1) Create a list of files QAresults.root produced by the PWG1 QA train,
//    to a text file, let say QAresults.txt
// 2) Compile this macro:
//    .L ExtractPHOSCellQA.C++
// 3) Run conversion of QAresults to new root files with PHOS histograms:
//    ExtractPHOSCellQA("QAresult.txt")
// 4) On the output, the new root files will be created per run
//    with the names AnyInt_<run>.root and PHI7_<run>.root
//
// Author: Yuri Kharlov. 3-Oct-2011
//------------------------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TObjArray.h>
#include <TString.h>
#include <TFile.h>
#include <TGrid.h>
#include <Riostream.h>
#include <stdio.h>
using namespace std;
#endif

void ExtractPHOSCellQA(const TString QAfilelist="QAresult.list")
{

  ifstream in;
  in.open(QAfilelist.Data());
  char rootFileName[256];
  TString oldRootFileName, newRootFileName;
  TFile *oldRootFile, *newRootFile;
  TObjArray *histAnyInt, *histPHI7;
  Bool_t firstFile = kTRUE;

  while (1) {
    in >> rootFileName;
    if (!in.good()) break;
    printf("root file is %s",rootFileName);
    oldRootFileName = rootFileName;
    if (oldRootFileName.BeginsWith("alien://")) {
      TGrid::Connect("alien://");
      firstFile = kFALSE;
    }
    oldRootFile = TFile::Open(rootFileName,"read");
    histAnyInt = (TObjArray*)oldRootFile->Get("PHOSCellsQA_AnyInt");
    histPHI7   = (TObjArray*)oldRootFile->Get("PHOSCellsQA_PHI7");
    if (histAnyInt == 0 || histPHI7 == 0) {
      printf(" does not contain PHOSCellQA histograms\n");
      continue;
    }
    else {
      printf(" contains PHOSCellQA histograms\n");
    }

    char *runNum = strtok(rootFileName+35,"/");

    newRootFileName = Form("AnyInt_%s.root",runNum);
    newRootFile = TFile::Open(newRootFileName,"recreate");
    histAnyInt->Write(0,0);
    newRootFile->Close();

    newRootFileName = Form("PHI7_%s.root",runNum);
    newRootFile = TFile::Open(newRootFileName,"recreate");
    histPHI7  ->Write(0,0);
    newRootFile->Close();

    histAnyInt->Clear();
    histPHI7  ->Clear();
    oldRootFile  ->Close();
   }
}
