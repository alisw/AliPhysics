#include "AliMultEstimator.h"
#include "AliMultSelectionCuts.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCalibrator.h"
#include "AliMultSelectionCalibratorMC.h"
#include <TString.h>
#include <TSystem.h>
#include <TF1.h>
#include <TFile.h>

void CalibratePeriodMC( TString lPeriodName         = "",
                        TString inputFileNameMC     = "",
                        TString inputFileNameData   = "",
                        TString dataOADBFile        = "",
                        Int_t defaultRunNumber      = 0,
                        Bool_t enableBufferFiles    = kFALSE,
                        Bool_t enableSuperCalib     = kFALSE
                      ){
  cout<<"Run!"<<endl;

  //Load ALICE stuff
  TString gLibs[] = {"STEER", "ANALYSIS", "ANALYSISalice", "ANALYSIScalib"};
  TString thislib = "lib";
  for(Int_t ilib = 0; ilib<4; ilib++){
    thislib="lib";
    thislib.Append(gLibs[ilib].Data());
    cout<<"Will load "<<thislib.Data()<<endl;
    gSystem->Load(thislib.Data());
  }
  gSystem->SetIncludePath("-I$ROOTSYS/include  -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  cout<<"Alive! "<<endl;

  //All fine, let's try the new MC calibrator
  AliMultSelectionCalibratorMC *lCalib = new AliMultSelectionCalibratorMC("lCalib");

  lCalib->SetupStandardInput();

  //Actual Input files
  lCalib -> SetInputFileData  ( inputFileNameData.Data() ) ;
  lCalib -> SetInputFileOADB  ( dataOADBFile.Data()) ;
  lCalib -> SetInputFileMC    ( inputFileNameMC.Data() ) ;

  // set default runnumber
  lCalib->SetRunToUseAsDefault(defaultRunNumber);
  TString addNameOutput  = "";
  if(enableSuperCalib){
    lCalib->SetUseQuadraticMapping(enableSuperCalib);
    addNameOutput   = "-SuperCalib";
  }

  //Buffer files
  if (enableBufferFiles){
    lCalib -> SetBufferFileData ( Form("buffer-data-anc%s.root", lPeriodName.Data()) );
    lCalib -> SetBufferFileMC   ( Form("buffer-MC-%s.root", lPeriodName.Data()) );
  }

  //Output OADB
  lCalib -> SetDebugFile     ( Form("debug-%s%s.root",lPeriodName.Data(),addNameOutput.Data()) );
  lCalib -> SetOutputFile     ( Form("OADB-%s%s.root",lPeriodName.Data(),addNameOutput.Data()) );
  lCalib -> Calibrate     ();

}
