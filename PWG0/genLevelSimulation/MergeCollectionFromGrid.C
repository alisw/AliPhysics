#include <iostream>
#include "TSystem.h"
#include "TFileMerger.h"
#include "TGrid.h"
#include "TAlienCollection.h"
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
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 

  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  // load dndeta task
  gROOT->LoadMacro("AliAnalysisTaskdNdetaMC.cxx+");


  TFileMerger * fileMerger  = new TFileMerger(0); // dont merge local files

  TGrid::Connect("alien://");
  TGridCollection * coll = TAlienCollection::Open (incollection);
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
