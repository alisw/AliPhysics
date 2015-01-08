/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// PHOS Preprocessor class. It runs by Shuttle at the end of the run,
// calculates calibration coefficients and dead/bad channels
// to be posted in OCDB
//
// Author: Hisayuki Torii, 11 August 2009
///////////////////////////////////////////////////////////////////////////////
#include "AliLog.h"
#include "TFile.h"
#include "TKey.h"
#include "TList.h"
#include "TString.h"
#include "TObjString.h"

#include "AliPHOSPreprocessorPHYS.h"

ClassImp(AliPHOSPreprocessorPHYS)

//_______________________________________________________________________________________
AliPHOSPreprocessorPHYS::AliPHOSPreprocessorPHYS() : AliPreprocessor("PHS",0) {
  //default constructor
}
//_______________________________________________________________________________________
AliPHOSPreprocessorPHYS::AliPHOSPreprocessorPHYS(AliShuttleInterface* shuttle): AliPreprocessor("PHS",shuttle) {
  // Constructor

  AddRunType("PHYSICS");
  AddRunType("COSMICS"); // Does this exist??
}
//_______________________________________________________________________________________
UInt_t AliPHOSPreprocessorPHYS::Process(TMap* /*valueSet*/){
  // process data retrieved by the Shuttle
  //
  
  TString runType = GetRunType();
  Log(Form("Run type: %s",runType.Data()));

  if( runType=="PHYSICS" || runType=="COSMICS") {
    
    Bool_t calibEmc_OK = CalibratePhys();
    if(calibEmc_OK) return 0;
    else return 1;
  }
  Log(Form("Unknown run type %s. Do nothing and return OK.",runType.Data()));
  return 0;
}
//_______________________________________________________________________________________
Bool_t AliPHOSPreprocessorPHYS::CalibratePhys(){
  //process PHYSICS event retrieved by the Shuttle
  //
  
  TList* list = 0;
  list = GetFileSources(kDAQ,"PHOSDApi0mip");
  if(!list) {
    Log(Form("DAQ sources list not found!"));
    return kFALSE;
  }
  if(!list->GetEntries()) {
    Log(Form("Got empty sources list. It seems PHYS DA did not produce any files!"));
    return kFALSE;
  }

  TIter iter(list);
  TObjString *source;
  
  while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
    AliInfo(Form("found source %s", source->String().Data()));

    TString fileName = GetFile(kDAQ, "PHOSDApi0mip", source->GetName());
    AliInfo(Form("Got filename: %s",fileName.Data()));

    TFile* file = TFile::Open(fileName);
    if(!file) {
      Log(Form("File %s is not opened, something goes wrong!",fileName.Data()));
      return kFALSE;
    }

    file->ls();
    //TList * keylist = file->GetListOfKeys();
    Int_t nkeys   = file->GetNkeys();
    if(nkeys==0){
      Log(Form("Not enough  (%d) for calibration.",nkeys));
      return 1; // it's not fatal! May be short run..
    }
    
    file->Close();
    delete file;


  }
  return kTRUE;
}
//_______________________________________________________________________________________
     

