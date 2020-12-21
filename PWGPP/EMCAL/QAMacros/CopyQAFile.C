#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "Riostream.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TSystem.h"
#include "TFile.h"
#include "TError.h"
#endif

///
/// \file copyQAFile.C
/// \ingroup EMCALOfflineMacros
/// \brief copy a QAresult.root file from alien to a local directory and check it.
///
/// copy a QAresult.root file from alien to a local directory and check it. 
///
/// \author Marie Germain, <Marie.Germain@subatech.in2p3.fr>, SUBATECH
///

//_____________________________________________________________________________
Int_t GetRunNumber(TString filePath)
{
  TObjArray* array = filePath.Tokenize("/");
  array->SetOwner();
  TString auxString = "";
  Int_t runNum = -1;
  for ( Int_t ientry=0; ientry<array->GetEntries(); ientry++ ) {
    auxString = array->At(ientry)->GetName();
    if ( auxString.Length() == 9 && auxString.IsDigit() ) {
      runNum = auxString.Atoi();
      break;
    }
  }
  delete array;

  if ( runNum < 0 ) {
    array = auxString.Tokenize("_");
    array->SetOwner();
    auxString = array->Last()->GetName();
    auxString.ReplaceAll(".root","");
    if ( auxString.IsDigit() )
      runNum = auxString.Atoi();
    delete array;
  }

  return runNum;
}

//_____________________________________________________________________________
Int_t CopyQAFile(TString inFilename, TString baseOutDir=".", Bool_t makeRunDir=kTRUE, TString changeFilename="", Int_t timeOut = 10)
{

  gSystem->Setenv("XRDCLIENTMAXWAIT",Form("%d",timeOut));
  gEnv->SetValue("XNet.RequestTimeout", timeOut);
  gEnv->SetValue("XNet.ConnectTimeout", timeOut);
  gEnv->SetValue("XNet.TransactionTimeout", timeOut);
  TFile::SetOpenTimeout(timeOut);


  if ( inFilename.Contains("alien://") && ! gGrid )

    if (! TGrid::Connect("alien://")) {
      Error(__FUNCTION__,"Error connecting to alien");
      return -1;
    }

  TObjArray* array = inFilename.Tokenize("/");
  array->SetOwner();
  TString outFilename = changeFilename.IsNull() ? array->Last()->GetName() : changeFilename.Data();
  delete array;

  if ( makeRunDir ) {
    Int_t runNumber = GetRunNumber(inFilename);
    if ( runNumber >= 0 ) {
      baseOutDir = Form("%s/%i", baseOutDir.Data(), runNumber);
      if ( gSystem->AccessPathName(baseOutDir.Data()) )
        gSystem->mkdir(baseOutDir.Data());
    }
    else Warning(__FUNCTION__,"run number not found!");
  }
  outFilename.Prepend(Form("%s/", baseOutDir.Data()));
  Bool_t showProgressBar = ! gROOT->IsBatch();
 
  if ( gSystem->AccessPathName(outFilename.Data())) {
    if (! TFile::Cp(inFilename.Data(), outFilename.Data(), showProgressBar)) {
      Error(__FUNCTION__,"Error copying the file from alien");
      return -2;
    }
  }

  printf("file: %s\n", inFilename.Data());
  printf("outDir: %s\n", baseOutDir.Data());
  printf("outFile: %s\n", outFilename.Data());

  gErrorIgnoreLevel = kWarning +1; 
  TFile f(outFilename.Data());
  gErrorIgnoreLevel = -1;
   
  if (f.IsZombie()) {
    Error(__FUNCTION__,"Error opening outFile");
    return -3;
  }

  if (f.TestBit(TFile::kRecovered)) {
    Info(__FUNCTION__,"The file is likely to be corrupted");
    return -4;
  }

  return 0; 
}

