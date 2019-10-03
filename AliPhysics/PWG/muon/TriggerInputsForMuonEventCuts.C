#if !defined(__CINT__) || defined(__MAKECINT__)

// ROOT includes
#include "Riostream.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "THashList.h"
#include "TArrayI.h"
#include "TMath.h"

// STEER includes
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerInput.h"
#endif

/// This macro reads the CTP configuration in the OCDB
/// and returns the trigger inputs for the requested run list.
/// The output is in a format that suits the AliMuonEventCuts class
///
/// The file with the input run list must contain one run number per line
/// If the user is interested only in a subset of inputs,
/// he/she can specify them in selectedInputs, separated by commas
/// e.g "0MSL,0MSH,0MUL,0MLL"
///
/// author: Diego Stocco, Subatech, Nantes

void TriggerInputsForMuonEventCuts ( TString runListFilename, TString selectedInputs="", TString defaultStorage = "raw://" )
{
  AliCDBManager::Instance()->SetDefaultStorage(defaultStorage.Data());
  TObjArray inputsList;
  inputsList.SetOwner();
  
  TObjArray* selectedInputsList = selectedInputs.Tokenize(",");

  // Read input run list
  ifstream inFile(runListFilename.Data());
  TString srun = "";
  if ( inFile.is_open() ) {
    while ( ! inFile.eof() ) {
      srun.ReadLine(inFile,kFALSE);
      if ( ! srun.IsDigit() ) continue;
      
      // For each run, read trigger inputs from OCDB
      Int_t runNumber = srun.Atoi();
      AliCDBManager::Instance()->SetRun(runNumber);
      
      // Get trigger class configuration
      AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
      if ( ! entry ) continue;
      
      THashList* runInputs = new THashList();
      runInputs->SetOwner();
      runInputs->SetUniqueID((UInt_t)runNumber);
      AliTriggerConfiguration* trigConf = (AliTriggerConfiguration*)entry->GetObject();
      const TObjArray& trigInputsArray = trigConf->GetInputs();
      AliTriggerInput* trigInput = 0x0;
      TIter next(&trigInputsArray);
      while ( ( trigInput = static_cast<AliTriggerInput*>(next()) ) ) {
        if ( selectedInputsList->GetEntriesFast() > 0 && ! selectedInputsList->FindObject(trigInput->GetName()) ) continue;
        Int_t inputId = (Int_t)TMath::Log2(trigInput->GetMask());
        TObjString* currInput = new TObjString(trigInput->GetName());
        currInput->SetUniqueID(inputId);
        runInputs->Add(currInput);
      }
      inputsList.Add(runInputs);
    }
    inFile.close();
  }
  delete selectedInputsList;
  
  // Loop on the trigger inputs
  // and group runs with an equal list of inputs
  Int_t nentries = inputsList.GetEntries();
  TArrayI checkMask(nentries);
  checkMask.Reset(1);
  for ( Int_t irun=0; irun<nentries; irun++ ) {
    if ( checkMask[irun] == 0 ) continue;
    THashList* currList = static_cast<THashList*>(inputsList.At(irun));
    TString runRange = Form("Run range: %u", currList->GetUniqueID());
    for ( Int_t jrun=irun+1; jrun<nentries; jrun++ ) {
      if ( checkMask[jrun] == 0 ) continue;
      THashList* checkList = static_cast<THashList*>(inputsList.At(jrun));
      Bool_t isDifferent = kFALSE;
      for ( Int_t itrig=0; itrig<currList->GetEntries(); itrig++ ) {
        TObjString* currInput = static_cast<TObjString*>(currList->At(itrig));
        TObject* checkInput = checkList->FindObject(currInput->GetName());
        if ( ! checkInput || checkInput->GetUniqueID() != currInput->GetUniqueID() ) {
          isDifferent = kTRUE;
          break;
        }
      } // loop on trigger inputs
      if ( isDifferent ) continue;
      checkMask[jrun] = 0;
      runRange += Form(",%u", checkList->GetUniqueID());
    } // loop on runs
    
    TString outString = "\nSetTrigInputsMap(\"";
    for ( Int_t itrig=0; itrig<currList->GetEntries(); itrig++ ) {
      TObjString* currInput = static_cast<TObjString*>(currList->At(itrig));
      outString += Form("%s:%u,",currInput->GetString().Data(), currInput->GetUniqueID());
    }
    outString.Append("\");\n");
    outString.ReplaceAll(",\"","\"");
    outString += runRange;
    printf("%s\n", outString.Data());
  } // loop on runs
}
