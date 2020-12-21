#include <iostream>
#include "TSystem.h"
#include "TFileMerger.h"
#include "TGrid.h"
#include "TGridCollection.h"
#include "TFile.h"
#include "TH1F.h"
#include "TList.h"
#include "TH1I.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "AliAnalysisTaskdNdetaMC.h"
#include "TROOT.h"

using namespace std;



void MergeCollectionFromGrid(const char * incollection = "test.xml", const char * outputfile= "dndeta_merged.root")
{
  // for running with root only
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");

  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  // load dndeta task
  gROOT->LoadMacro("AliAnalysisTaskdNdetaMC.cxx+");


  TFileMerger * fileMerger  = new TFileMerger(0); // dont merge local files

  TGrid::Connect("alien://");
  TGridCollection * coll = gGrid->OpenCollection(incollection);
  Int_t  ifile=0;
  while(coll->Next()){
    fileMerger->AddFile(TString("alien://")+coll->GetLFN());
    ifile++;
    //    if(ifile>2) break;
  }
  fileMerger->OutputFile("tmp.root");
  fileMerger->Merge();

  // Reopen the merged file, normalize histos and save them back.
  // SOME DUPLICATED CODE... SOME CLEAN UP WOULD BE GOOD
  AliAnalysisTaskdNdetaMC * localTask  = new AliAnalysisTaskdNdetaMC("merger", "tmp.root");
  localTask->Finalize();
  
  TFile * fout = new TFile (outputfile, "recreate");
  localTask->GetList()->Write();
  fout->Close();
}
