#ifndef __CINT__
#include "TFileMerger.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"
#endif

void recursiveMergingRCT(TString pathNameRCT, TString pathNameWithWildCards, TString outputPath, TString outputFile,
                         Int_t nmerge = 5, Bool_t cleanUpIntermediateFiles = kTRUE)
{
  //e.g. pathNameRCT = "/hera/alice/bhess/macros/myRCTlist.txt"
  //e.g. pathNameWithWildCards = "/hera/alice/bhess/train/V010.aod.pp/myTrain/mergedRuns/aod.pp.new/LHC10*/*/bhess_PID_Jets_Inclusive_PureGauss.root";
  
  TString runListFilePathName = Form("%s/%s2.log", outputPath.Data(), outputFile.Data());
  gSystem->Exec(Form("ls %s > %s", pathNameWithWildCards.Data(), runListFilePathName.Data()));
  
  TString runListFilePathName2 = runListFilePathName;
  runListFilePathName2 = runListFilePathName2.ReplaceAll("2.log", ".log");
  gSystem->Exec(Form("for ii in $(cat %s); do grep $ii %s; done > %s", pathNameRCT.Data(), runListFilePathName.Data(), 
                runListFilePathName2.Data()));
  
  
  TString s=gSystem->GetFromPipe(Form("cat %s", runListFilePathName2.Data()));
  
  gSystem->Exec(Form("rm %s", runListFilePathName.Data()));

  TObjArray *arr=s.Tokenize("\n");
  
  if (arr->GetEntries() == 0)
    return;
  
  TString outputPathName = Form("%s/%s", outputPath.Data(), outputFile.Data());

  Int_t depth = 0;
  while (arr->GetEntries() > 1){
    printf("depth: %d\n",depth);
    for (Int_t iIter=0; iIter<TMath::Ceil((Double_t)arr->GetEntries()/((Double_t)nmerge)); ++iIter){
      printf("Iter: %d\n",iIter);
      TFileMerger m(0);
      m.OutputFile(Form("%s.%d.%d.root",outputPathName.Data(),depth,iIter));
      printf("writing output file: %s\n", Form("%s.%d.%d.root",outputPathName.Data(),depth,iIter));
      for (Int_t ifile=iIter*nmerge; ifile<(iIter+1)*nmerge; ++ifile){
        if (!arr->At(ifile)) continue;
        printf("Adding file: %s\n",arr->At(ifile)->GetName());
        m.AddFile(arr->At(ifile)->GetName());
      }
      m.Merge();
    }
    delete arr;
    arr=0x0;
    s=gSystem->GetFromPipe(Form("ls %s.%d.[0-9]*.root",outputPathName.Data(),depth));
    arr=s.Tokenize("\n");
    ++depth;
    printf("\n-----------\n");
  }
  
  gSystem->Exec(Form("mv %s.%d.0.root %s.root",outputPathName.Data(),depth-1,outputPathName.Data()));
  
  if (cleanUpIntermediateFiles) {
    delete arr;
    arr=0x0;
    s=gSystem->GetFromPipe(Form("ls %s.[0-9]*.[0-9]*.root",outputPathName.Data()));
    arr=s.Tokenize("\n");
    
    for (Int_t i=0; i<arr->GetEntries(); ++i){
      printf("Removing file: %s\n", arr->At(i)->GetName());
      gSystem->Exec(Form("rm %s", arr->At(i)->GetName()));
    }
  }
}

