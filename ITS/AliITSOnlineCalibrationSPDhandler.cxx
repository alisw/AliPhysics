//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// For easier handling of dead and noisy pixels they are kept in    //
// container maps (AliITSIntMap).                                   //
// Handling of inactive equipments,HSs,chips have been added.       //
// A pixel that is either dead or inactive is called SILENT here.   //
// The lists of single dead and noisy pixels are separated from the //
// information about which eq/hs/chip are inactive.                 //
// The TArrayS objects that are put in the AliITSCalibrationSPD     //
// objects can be obtained from the methods GetDeadArray and        //
// GetNoisyArray.                                                   //
//////////////////////////////////////////////////////////////////////   

#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliITSOnlineCalibrationSPD.h"
#include "AliITSTriggerConditions.h"
#include "AliITSIntMap.h"
#include <TObjArray.h>
#include <TArrayI.h>
#include <TArrayS.h>
#include <TFile.h>
#include <TError.h>
#include <fstream>

#ifndef SPD_DA_OFF // you may want to exclude cdb classes, if using this class outside aliroot
#include "AliITSCalibrationSPD.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#endif

/* $Id$ */

//____________________________________________________________________________________________
AliITSOnlineCalibrationSPDhandler::AliITSOnlineCalibrationSPDhandler():
  fFileLocation("."),
  fTriggerConditions(0)
{
  // constructor
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    fNrDead[gloChip]=0;
    fNrSparseDead[gloChip]=0;
    fNrNoisy[gloChip]=0;
    fDeadPixelMap[gloChip] = new AliITSIntMap();
    fSparseDeadPixelMap[gloChip] = new AliITSIntMap();
    fNoisyPixelMap[gloChip] = new AliITSIntMap();    
  }
  ActivateALL();
  UnSetDeadALL();
}
//____________________________________________________________________________________________
AliITSOnlineCalibrationSPDhandler::AliITSOnlineCalibrationSPDhandler(const AliITSOnlineCalibrationSPDhandler& handle): 
  fFileLocation("."),
  fTriggerConditions(handle.fTriggerConditions)
{
  // copy constructor
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    fNrDead[gloChip] = handle.fNrDead[gloChip];
    fNrSparseDead[gloChip] = handle.fNrSparseDead[gloChip];
    fNrNoisy[gloChip] = handle.fNrNoisy[gloChip];
    fDeadPixelMap[gloChip] = handle.fDeadPixelMap[gloChip]->Clone();
    fSparseDeadPixelMap[gloChip] = handle.fSparseDeadPixelMap[gloChip]->Clone();
    fNoisyPixelMap[gloChip] = handle.fNoisyPixelMap[gloChip]->Clone();
  }
  for (UInt_t eq=0; eq<20; eq++) {
    fActiveEq[eq] = handle.fActiveEq[eq];
    fDeadEq[eq]=handle.fDeadEq[eq];
    for (UInt_t hs=0; hs<6; hs++) {
      fActiveHS[eq][hs] = handle.fActiveHS[eq][hs];
      for (UInt_t chip=0; chip<10; chip++) {
	fActiveChip[eq][hs][chip] = handle.fActiveChip[eq][hs][chip];
      }
    }
  }
  fFileLocation = handle.fFileLocation;
}
//____________________________________________________________________________________________
AliITSOnlineCalibrationSPDhandler::~AliITSOnlineCalibrationSPDhandler() {
  //  ClearMaps();
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    delete fDeadPixelMap[gloChip];
    delete fSparseDeadPixelMap[gloChip];
    delete fNoisyPixelMap[gloChip];
  }
}
//____________________________________________________________________________________________
AliITSOnlineCalibrationSPDhandler& AliITSOnlineCalibrationSPDhandler::operator=(const AliITSOnlineCalibrationSPDhandler& handle) {
  // assignment operator
  if (this!=&handle) {
    this->ClearMaps();
    for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
      fNrDead[gloChip] = handle.fNrDead[gloChip];
      fNrSparseDead[gloChip] = handle.fNrSparseDead[gloChip];
      fNrNoisy[gloChip] = handle.fNrNoisy[gloChip];
      fDeadPixelMap[gloChip] = handle.fDeadPixelMap[gloChip]->Clone();
      fSparseDeadPixelMap[gloChip] = handle.fSparseDeadPixelMap[gloChip]->Clone();
      fNoisyPixelMap[gloChip] = handle.fNoisyPixelMap[gloChip]->Clone();
    }
    for (UInt_t eq=0; eq<20; eq++) {
      fActiveEq[eq] = handle.fActiveEq[eq];
      fDeadEq[eq] = handle.fDeadEq[eq];
      for (UInt_t hs=0; hs<6; hs++) {
	fActiveHS[eq][hs] = handle.fActiveHS[eq][hs];
	for (UInt_t chip=0; chip<10; chip++) {
	  fActiveChip[eq][hs][chip] = handle.fActiveChip[eq][hs][chip];
	}
      }
    }
    fFileLocation = handle.fFileLocation;
  }
  return *this;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ClearMaps() {
  // clear the lists of dead and noisy
  ResetDead();
  ResetNoisy();
  ActivateALL();
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ResetDead() {
  // reset the dead pixel map and inactive eq,hs,chip
  UnSetDeadALL();
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    fNrDead[gloChip]=0;
    fNrSparseDead[gloChip]=0;
    fDeadPixelMap[gloChip]->Clear();
    fSparseDeadPixelMap[gloChip]->Clear();
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ResetNoisy() {
  // clear the list of noisy pixels
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    fNrNoisy[gloChip]=0;
    fNoisyPixelMap[gloChip]->Clear();
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ResetDeadForChip(UInt_t eq, UInt_t hs, UInt_t chip) {
  // clear the dead pixels for this chip
  SetDeadChip(eq,hs,chip,kFALSE);
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      Int_t key = GetKey(eq,hs,chip,col,row);
      if (fDeadPixelMap[gloChip]->Remove(key)) {
	fNrDead[gloChip]--;
      }
      if (fSparseDeadPixelMap[gloChip]->Remove(key)) {
	fNrSparseDead[gloChip]--;
      }
    }
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ResetNoisyForChip(UInt_t eq, UInt_t hs, UInt_t chip) {
  // clear the noisy pixels for this chip
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::ResetNoisyForChip","global chip nr (%d) out of bounds\n",gloChip);
    return;
  }
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      Int_t key = GetKey(eq,hs,chip,col,row);
      if (fNoisyPixelMap[gloChip]->Remove(key)) {
	fNrNoisy[gloChip]--;
      }
    }
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ResetDeadForEq(UInt_t eq) {
  // clear the dead pixels for this eq
  if (eq>=20) {
    Error("AliITSOnlineCalibrationSPDhandler::ResetDeadForEq", "eq (%d) out of bounds.",eq);
    return;
  }
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      ResetDeadForChip(eq, hs, chip);
    }
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ResetNoisyForEq(UInt_t eq) {
  // clear the noisy pixels for this eq
  if (eq>=20) {
    Error("AliITSOnlineCalibrationSPDhandler::ResetNoisyForEq", "eq (%d) out of bounds.",eq);
    return;
  }
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      ResetNoisyForChip(eq, hs, chip);
    }
  }
}


//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadFromFiles() {
  // read files from file location (active,dead,noisy info). returns true if at least one file found
  Bool_t b1 = ReadNoisyFromFiles();
  Bool_t b2 = ReadSilentFromFiles();
  return (b1 || b2);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadSilentFromFiles() {
  // read dead,active files from file location. returns true if at least one file found
  Bool_t returnval=kFALSE;
  for (UInt_t eq=0; eq<20; eq++) {
    if (ReadSilentFromFile(eq)) {
      returnval=kTRUE;
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromFiles() {
  // read dead,active files from file location. returns true if at least one file found
  Bool_t returnval=kFALSE;
  for (UInt_t eq=0; eq<20; eq++) {
    if (ReadDeadFromFile(eq)) {
      returnval=kTRUE;
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadSilentFromFile(UInt_t eq) {
  // read dead file for eq from file location. 
  TString fileName = Form("%s/SPD_Dead_%d.root",fFileLocation.Data(),eq);
  return ReadSilentFromFileName(fileName.Data());
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromFile(UInt_t eq) {
  // read dead file for eq from file location. 
  TString fileName = Form("%s/SPD_Dead_%d.root",fFileLocation.Data(),eq);
  return ReadDeadFromFileName(fileName.Data());
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadSilentFromFileName(const char *fileName) {
  // read dead from file fileName (including inactive)
  return ReadDeadFromFileName(fileName, kTRUE);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromFileName(const char *fileName, Bool_t inactive) {
  // read dead from file fileName
  AliITSOnlineCalibrationSPD* calib;
  FILE* fp0 = fopen(fileName, "r");
  if (fp0 == NULL) {return kFALSE;}
  else {
    fclose(fp0);
    TFile file(fileName, "READ");
    if (file.IsOpen()) {
      file.GetObject("AliITSOnlineCalibrationSPD", calib);
      file.Close();
      if (calib!=NULL) {
	UInt_t nrDead=calib->GetNrBad();
	for (UInt_t index=0; index<nrDead; index++) {
	  UInt_t key = calib->GetKeyAt(index);
	  UInt_t eq = GetEqIdFromKey(key);
	  UInt_t hs = GetHSFromKey(key);
	  UInt_t chip = GetChipFromKey(key);
	  UInt_t col = GetColFromKey(key);
	  UInt_t row = GetRowFromKey(key);
	  SetDeadPixel(eq,hs,chip,col,row);
	}
	UInt_t eq1 = calib->GetEqNr();
	if (calib->IsDeadEq()) SetDeadEq(eq1);
	else                   SetDeadEq(eq1,kFALSE);
	for (UInt_t hs=0; hs<6; hs++) {
	  if (calib->IsDeadHS(hs)) SetDeadHS(eq1,hs);
	  else                     SetDeadHS(eq1,hs,kFALSE);
	  for (UInt_t chip=0; chip<10; chip++) {
	    if (calib->IsDeadChip(hs,chip)) SetDeadChip(eq1,hs,chip);
	    else                            SetDeadChip(eq1,hs,chip,kFALSE);
	  }
	}
	if (inactive) {
	  UInt_t eq = calib->GetEqNr();
	  if (calib->IsActiveEq()) ActivateEq(eq);
	  else                     ActivateEq(eq,kFALSE);
	  for (UInt_t hs=0; hs<6; hs++) {
	    if (calib->IsActiveHS(hs)) ActivateHS(eq,hs);
	    else                       ActivateHS(eq,hs,kFALSE);
	    for (UInt_t chip=0; chip<10; chip++) {
	      if (calib->IsActiveChip(hs,chip)) ActivateChip(eq,hs,chip);
	      else                              ActivateChip(eq,hs,chip,kFALSE);
	    }
	  }
	}
      }
    }
  }
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromFiles() {
  // read noisy files from file location. returns true if at least one file found
  Bool_t returnval=kFALSE;
  for (UInt_t eq=0; eq<20; eq++) {
    if (ReadNoisyFromFile(eq)) {
      returnval=kTRUE;
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromFile(UInt_t eq) {
  // read noisy file for eq from file location. 
  TString fileName = Form("%s/SPD_Noisy_%d.root",fFileLocation.Data(),eq);
  return ReadNoisyFromFileName(fileName.Data());
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromFileName(const char *fileName) {
  // read noisy from file fileName
  AliITSOnlineCalibrationSPD* calib;
  FILE* fp0 = fopen(fileName, "r");
  if (fp0 == NULL) {return kFALSE;}
  else {
    fclose(fp0);
    TFile file(fileName, "READ");
    if (file.IsOpen()) {
      file.GetObject("AliITSOnlineCalibrationSPD", calib);
      file.Close();
      if (calib!=NULL) {
	UInt_t nrNoisy=calib->GetNrBad();
	for (UInt_t index=0; index<nrNoisy; index++) {
	  UInt_t key = calib->GetKeyAt(index);
	  UInt_t eq = GetEqIdFromKey(key);
	  UInt_t hs = GetHSFromKey(key);
	  UInt_t chip = GetChipFromKey(key);
	  UInt_t col = GetColFromKey(key);
	  UInt_t row = GetRowFromKey(key);
	  SetNoisyPixel(eq,hs,chip,col,row);
	}
      }
    }
  }
  return kTRUE;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromText(const char *fileName, UInt_t module) {
  // read dead from a text file (lines with eq,hs,chip,col,row). returns nr of pixels added (not already here)
  // insert only those pixels that belong to module (or all if module=240). 
  UInt_t newNrDead=0;
  ifstream textFile;
  textFile.open(fileName, ifstream::in);
  if (textFile.fail()) {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadDeadFromText","No dead text file (%s) present.",fileName);
  }
  else {
    while(1) {
      UInt_t eq,hs,chip,col,row;
      textFile >> eq; if (textFile.eof()) break;
      textFile >> hs; if (textFile.eof()) break;
      textFile >> chip; if (textFile.eof()) break;
      textFile >> col; if (textFile.eof()) break;
      textFile >> row; 
      if (module==240 || (Int_t)module==AliITSRawStreamSPD::GetModuleNumber(eq,hs,chip)){
	if (SetDeadPixel(eq,hs,chip,col,row)) {
	  newNrDead++;
	}
      }
      if (textFile.eof()) break;
    }
    textFile.close();
  }
  return newNrDead;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromText(const char *fileName, UInt_t module) {
  // read noisy from a text file (lines with eq,hs,chip,col,row). returns nr of pixels added (not already here)
  // insert only those pixels that belong to module (or all if module=240). 
  UInt_t newNrNoisy=0;
  ifstream textFile;
  textFile.open(fileName, ifstream::in);
  if (textFile.fail()) {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadNoisyFromText","No noisy text file (%s) present.",fileName);
  }
  else {
    while(1) {
      UInt_t eq,hs,chip,col,row;
      textFile >> eq; if (textFile.eof()) break;
      textFile >> hs; if (textFile.eof()) break;
      textFile >> chip; if (textFile.eof()) break;
      textFile >> col; if (textFile.eof()) break;
      textFile >> row; 
      if (module==240 || (Int_t)module==AliITSRawStreamSPD::GetModuleNumber(eq,hs,chip)){
	if (SetNoisyPixel(eq,hs,chip,col,row)) {
	  newNrNoisy++;
	}
      }
      if (textFile.eof()) break;
    }
    textFile.close();
  }
  return newNrNoisy;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ReadPITConditionsFromText(const char *fileName) {
  // read PIT conditions file from text as printed out at P2
  //  !!! please note that the chip numbering goes from 9 to 0 in the text. In PVSS panels is the opposite.
  if(fTriggerConditions) fTriggerConditions->ResetAll();
  else fTriggerConditions = new AliITSTriggerConditions();
  fTriggerConditions->ReadFromTextFile(fileName);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadPITConditionsFromDB(Int_t runNr, const Char_t *storage){
// read PIT conditions from the OCDB

  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBEntry *cdbEntry = man->Get("TRIGGER/SPD/PITConditions", runNr);
 if(cdbEntry) {
  fTriggerConditions = (AliITSTriggerConditions*)cdbEntry->GetObject();
  return kTRUE;
  } else return kFALSE;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::WriteToFilesAlways() {
  // write the lists of active,dead,noisy to files
  for (UInt_t eq=0; eq<20; eq++) {
    WriteSilentToFile(eq);
    WriteNoisyToFile(eq);
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::WriteToFiles() {
  // write the lists of dead and noisy to files (only if there are >0 dead or noisy pixels) , returns nr of files produced
  return (WriteNoisyToFiles() + WriteSilentToFiles());
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::WriteSilentToFilesAlways() {
  // write the lists of silent to files
  for (UInt_t eq=0; eq<20; eq++) {
      WriteSilentToFile(eq);
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::WriteDeadToFilesAlways() {
  // write the lists of dead to files
  for (UInt_t eq=0; eq<20; eq++) {
      WriteDeadToFile(eq);
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::WriteNoisyToFilesAlways() {
  // write the lists of noisy to files
  for (UInt_t eq=0; eq<20; eq++) {
    WriteNoisyToFile(eq);
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::WriteSilentToFiles() {
  // write the list of silent to file (only if there are >0 silent pixels) , returns nr of files produced
  UInt_t nrFiles=0;
  for (UInt_t eq=0; eq<20; eq++) {
    if (GetNrSilentEq(eq) > 0) {
      WriteSilentToFile(eq);
      nrFiles++;
    }
  }
  return nrFiles;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::WriteDeadToFiles() {
  // write the list of dead to file (only if there are >0 dead pixels) , returns nr of files produced
  UInt_t nrFiles=0;
  for (UInt_t eq=0; eq<20; eq++) {
    if (GetNrDeadEq(eq) > 0) {
      WriteDeadToFile(eq);
      nrFiles++;
    }
  }
  return nrFiles;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::WriteNoisyToFiles() {
  // write the list of noisy to file (only if there are >0 noisy pixels) , returns nr of files produced
  UInt_t nrFiles=0;
  for (UInt_t eq=0; eq<20; eq++) {
    if (GetNrNoisyEq(eq) > 0) {
      WriteNoisyToFile(eq);
      nrFiles++;
    }
  }
  return nrFiles;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::WriteSilentToFile(UInt_t eq) {
  WriteDeadToFile(eq,kTRUE);
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::WriteDeadToFile(UInt_t eq, Bool_t inactive) {
  // write the lists of dead (and inactive if input boolean is true) for eq to file
  AliITSOnlineCalibrationSPD* calib = new AliITSOnlineCalibrationSPD();
  calib->SetEqNr(eq);
  calib->SetBadList(GetDeadArrayOnline(eq));
  calib->SetNrBad(GetNrDeadEq(eq));
  if (IsDeadEq(eq)) calib->SetDeadEq();
  else              calib->SetDeadEq(kFALSE);
  for (UInt_t hs=0; hs<6; hs++) {
    if (IsDeadHS(eq,hs)) calib->SetDeadHS(hs);
    else                 calib->SetDeadHS(hs,kFALSE);
    for (UInt_t chip=0; chip<10; chip++) {
      if (IsDeadChip(eq,hs,chip)) calib->SetDeadChip(hs,chip);
      else                        calib->SetDeadChip(hs,chip,kFALSE);
    }
  }
  if (inactive) {
    if (IsActiveEq(eq)) calib->ActivateEq();
    else                calib->ActivateEq(kFALSE);
    for (UInt_t hs=0; hs<6; hs++) {
      if (IsActiveHS(eq,hs)) calib->ActivateHS(hs);
      else                   calib->ActivateHS(hs,kFALSE);
      for (UInt_t chip=0; chip<10; chip++) {
	if (IsActiveChip(eq,hs,chip)) calib->ActivateChip(hs,chip);
	else                          calib->ActivateChip(hs,chip,kFALSE);
      }
    }
  }
  TString fileName = Form("%s/SPD_Dead_%d.root",fFileLocation.Data(),eq);
  TFile file(fileName.Data(), "RECREATE");
  file.WriteTObject(calib, "AliITSOnlineCalibrationSPD");
  file.Close();
  delete calib;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::WriteNoisyToFile(UInt_t eq) {
  // write the lists of noisy for eq to file
  AliITSOnlineCalibrationSPD* calib = new AliITSOnlineCalibrationSPD();
  calib->SetEqNr(eq);
  calib->SetBadList(GetNoisyArrayOnline(eq));
  calib->SetNrBad(GetNrNoisyEq(eq));
  TString fileName = Form("%s/SPD_Noisy_%d.root",fFileLocation.Data(),eq);
  TFile file(fileName.Data(), "RECREATE");
  file.WriteTObject(calib, "AliITSOnlineCalibrationSPD");
  file.Close();
  delete calib;
}
//____________________________________________________________________________________________
#ifndef SPD_DA_OFF
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadModuleFromDB(UInt_t module, Int_t runNr, const Char_t *storage, Bool_t treeSerial) {
  // reads dead pixels from DB for given module and runNr
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBEntry *cdbEntry = man->Get("ITS/Calib/SPDDead", runNr);
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return kFALSE;
  }
  else {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadDeadModuleFromDB","Calibration for run %d not found in database.",runNr);
    return kFALSE;
  }
  AliITSCalibrationSPD* calibSPD;
  calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);

  UInt_t nrDead = calibSPD->GetNrBadSingle();
  if (nrDead>0) {
    if (!treeSerial) RecursiveInsertDead(calibSPD,module,0,nrDead-1);
    else {
      for (UInt_t index=0; index<nrDead; index++) {
	UInt_t colM = calibSPD->GetBadColAt(index);
	UInt_t rowM = calibSPD->GetBadRowAt(index);
	SetDeadPixelM(module,colM,rowM);
      }
    }
  }
  for (UInt_t chipIndex=0; chipIndex<5; chipIndex++) {
    UInt_t eq,hs,chip,col,row;
    AliITSRawStreamSPD::OfflineToOnline(module, chipIndex*32, 0, eq, hs, chip, col, row);
    if (calibSPD->IsChipBad(chipIndex)) {
      SetDeadChip(eq,hs,chip);
    }
    else {
      SetDeadChip(eq,hs,chip,kFALSE);
    }
  }

  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyModuleFromDB(UInt_t module, Int_t runNr, const Char_t *storage, Bool_t treeSerial) {
  // reads noisy pixels from DB for given module and runNr
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBEntry *cdbEntry = man->Get("ITS/Calib/SPDNoisy", runNr);
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return kFALSE;
  }
  else {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadNoisyModuleFromDB","Calibration for run %d not found in database.",runNr);
    return kFALSE;
  }
  AliITSCalibrationSPD* calibSPD;
  calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
  UInt_t nrNoisy = calibSPD->GetNrBadSingle();
  if (nrNoisy>0) {
    if (!treeSerial) RecursiveInsertNoisy(calibSPD,module,0,nrNoisy-1);
    else {
      for (UInt_t index=0; index<nrNoisy; index++) {
	UInt_t colM = calibSPD->GetBadColAt(index);
	UInt_t rowM = calibSPD->GetBadRowAt(index);
	SetNoisyPixelM(module,colM,rowM);
      }
    }
  }
  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadFromDB(Int_t runNr, const Char_t *storage, Bool_t treeSerial) {
  // reads dead and noisy pixels from DB for given runNr
  // note that you may want to clear the lists (if they are not empty) before reading
  return (ReadNoisyFromDB(runNr,storage,treeSerial) && ReadDeadFromDB(runNr,storage,treeSerial));
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromDB(Int_t runNr, const Char_t *storage, Bool_t treeSerial) {
  // reads dead pixels from DB for given runNr
  // note that you may want to clear the list (if it is not empty) before reading
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBEntry *cdbEntry = man->Get("ITS/Calib/SPDDead", runNr);
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return kFALSE;
  }
  else {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadDeadFromDB","Calibration for run %d not found in database.",runNr);
    return kFALSE;
  }
  AliITSCalibrationSPD* calibSPD;
  for (UInt_t module=0; module<240; module++) {
    calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    UInt_t nrDead = calibSPD->GetNrBadSingle();
    if (nrDead>0) {
      if (!treeSerial) {
	RecursiveInsertDead(calibSPD,module,0,nrDead-1);
      }

      else {
	for (UInt_t index=0; index<nrDead; index++) {
	  UInt_t colM = calibSPD->GetBadColAt(index);
	  UInt_t rowM = calibSPD->GetBadRowAt(index);
	  SetDeadPixelM(module,colM,rowM);
	}
      }
    }
    for (UInt_t chipIndex=0; chipIndex<5; chipIndex++) {
      UInt_t eq,hs,chip,col,row;
      AliITSRawStreamSPD::OfflineToOnline(module, chipIndex*32, 0, eq, hs, chip, col, row);
      if (calibSPD->IsChipBad(chipIndex)) {
	SetDeadChip(eq,hs,chip);
      }
      else {
	SetDeadChip(eq,hs,chip,kFALSE);
      }
    }
  }
  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadSparseDeadFromDB(Int_t runNr, const Char_t *storage, Bool_t treeSerial) {
  // reads dead pixels from DB for given runNr
  // note that you may want to clear the list (if it is not empty) before reading
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBEntry *cdbEntry = man->Get("ITS/Calib/SPDSparseDead", runNr);
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return kFALSE;
  }
  else {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadSparseDeadFromDB","Calibration for run %d not found in database.",runNr);
    return kFALSE;
  }
  AliITSCalibrationSPD* calibSPD;
  for (UInt_t module=0; module<240; module++) {
    calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    UInt_t nrDead = calibSPD->GetNrBadSingle();
    if (nrDead>0) {
      if (!treeSerial) {
        RecursiveInsertSparseDead(calibSPD,module,0,nrDead-1);
      }

      else {
        for (UInt_t index=0; index<nrDead; index++) {
          UInt_t colM = calibSPD->GetBadColAt(index);
          UInt_t rowM = calibSPD->GetBadRowAt(index);
          SetSparseDeadPixelM(module,colM,rowM);
        }
      }
    }
    for (UInt_t chipIndex=0; chipIndex<5; chipIndex++) {
      UInt_t eq,hs,chip,col,row;
      AliITSRawStreamSPD::OfflineToOnline(module, chipIndex*32, 0, eq, hs, chip, col, row);
      if (calibSPD->IsChipBad(chipIndex)) {
        SetDeadChip(eq,hs,chip);
      }
      else {
        SetDeadChip(eq,hs,chip,kFALSE);
      }
    }
  }
  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromDB(Int_t runNr, const Char_t *storage, Bool_t treeSerial) {
  // reads noisy pixels from DB for given runNr
  // note that you may want to clear the list (if it is not empty) before reading
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBEntry *cdbEntry = man->Get("ITS/Calib/SPDNoisy", runNr);
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return kFALSE;
  }
  else {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadNoisyFromDB","Calibration for run %d not found in database.",runNr);
    return kFALSE;
  }
  AliITSCalibrationSPD* calibSPD;
  for (UInt_t module=0; module<240; module++) {
    calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    UInt_t nrNoisy = calibSPD->GetNrBadSingle();
    if (nrNoisy>0) {    
      if (!treeSerial) {
	RecursiveInsertNoisy(calibSPD,module,0,nrNoisy-1);
      }
      else {
	for (UInt_t index=0; index<nrNoisy; index++) {
	  UInt_t colM = calibSPD->GetBadColAt(index);
	  UInt_t rowM = calibSPD->GetBadRowAt(index);
	  SetNoisyPixelM(module,colM,rowM);
	}
      }
    }
  }
  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromDBasNoisy(Int_t runNr, const Char_t *storage, Bool_t treeSerial) {
  // reads dead pixels (put as noisy) from DB for given runNr
  // note that you may want to clear the list (if it is not empty) before reading
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBEntry *cdbEntry = man->Get("ITS/Calib/SPDNoisy", runNr);
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return kFALSE;
  }
  else {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadDeadFromDBasNoisy","Calibration for run %d not found in database.",runNr);
    return kFALSE;
  }
  AliITSCalibrationSPD* calibSPD;
  for (UInt_t module=0; module<240; module++) {
    calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    UInt_t nrDead = calibSPD->GetNrBadSingle();
    if (nrDead>0) {
      if (!treeSerial) {
	RecursiveInsertDead(calibSPD,module,0,nrDead-1);
      }

      else {
	for (UInt_t index=0; index<nrDead; index++) {
	  UInt_t colM = calibSPD->GetBadColAt(index);
	  UInt_t rowM = calibSPD->GetBadRowAt(index);
	  SetDeadPixelM(module,colM,rowM);
	}
      }
    }
    for (UInt_t chipIndex=0; chipIndex<5; chipIndex++) {
      UInt_t eq,hs,chip,col,row;
      AliITSRawStreamSPD::OfflineToOnline(module, chipIndex*32, 0, eq, hs, chip, col, row);
      if (calibSPD->IsChipBad(chipIndex)) {
	SetDeadChip(eq,hs,chip);
      }
      else {
	SetDeadChip(eq,hs,chip,kFALSE);
      }
    }
  }
  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromCalibObj(TObjArray* calObj) {
  // reads dead pixels from calib object
  for (UInt_t module=0; module<240; module++) {
    AliITSCalibrationSPD* calibSPD = (AliITSCalibrationSPD*)calObj->At(module);
    for (Int_t i=0; i<calibSPD->GetNrBadSingle(); i++) {
      SetDeadPixelM(module,calibSPD->GetBadColAt(i),calibSPD->GetBadRowAt(i));
    }
    for (UInt_t chipIndex=0; chipIndex<5; chipIndex++) {
      UInt_t eq,hs,chip,col,row;
      AliITSRawStreamSPD::OfflineToOnline(module, chipIndex*32, 0, eq, hs, chip, col, row);
      if (calibSPD->IsChipBad(chipIndex)) {
	SetDeadChip(eq,hs,chip);
      }
      else {
	SetDeadChip(eq,hs,chip,kFALSE);
      }
    }
  }
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromCalibObj(TObjArray* calObj) {
  // reads noisy pixels from calib object
  for (UInt_t module=0; module<240; module++) {
    AliITSCalibrationSPD* calibSPD = (AliITSCalibrationSPD*)calObj->At(module);
    for (Int_t i=0; i<calibSPD->GetNrBadSingle(); i++) {
      SetNoisyPixelM(module,calibSPD->GetBadColAt(i),calibSPD->GetBadRowAt(i));
    }
  }
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::WriteToDB(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage) {
  // writes dead and noisy pixels to DB for given runNrs
  // overwrites any previous entries
  return (WriteNoisyToDB(runNrStart,runNrEnd,storage) && WriteDeadToDB(runNrStart,runNrEnd,storage));
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::WriteDeadToDB(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage) {
  // writes dead pixels to DB for given runNrs
  // overwrites any previous entries
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Henrik Tydesjo");
  metaData->SetComment("Created by AliITSOnlineCalibrationSPDhandler.");
  AliCDBId idCalSPD("ITS/Calib/SPDDead",runNrStart,runNrEnd);
  TObjArray* spdEntry = new TObjArray(240);
  spdEntry->SetOwner(kTRUE);
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calibSPD = new AliITSCalibrationSPD();
    spdEntry->Add(calibSPD);
  }
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    calibSPD->SetNrBadSingle( GetNrDeadSingle(module) );
    calibSPD->SetBadList( GetDeadArray(module) );
    for (UInt_t chipIndex=0; chipIndex<5; chipIndex++) {
      UInt_t eq,hs,chip,col,row;
      AliITSRawStreamSPD::OfflineToOnline(module, chipIndex*32, 0, eq, hs, chip, col, row);
      if (IsSilentChip(eq,hs,chip)) {
	calibSPD->SetChipBad(chipIndex);
      }
      else {
	calibSPD->UnSetChipBad(chipIndex);
      }
    }
  }
  AliCDBEntry* cdbEntry = new AliCDBEntry((TObject*)spdEntry,idCalSPD,metaData);
  man->Put(cdbEntry);
  delete spdEntry;
  delete cdbEntry;
  delete metaData;
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::WriteSparseDeadToDB(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage) {
  // writes dead pixels to DB for given runNrs
  // overwrites any previous entries
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Annalisa Mastroserio");
  metaData->SetComment("Created by AliITSOnlineCalibrationSPDhandler.");
  AliCDBId idCalSPD("ITS/Calib/SPDSparseDead",runNrStart,runNrEnd);
  TObjArray* spdEntry = new TObjArray(240);
  spdEntry->SetOwner(kTRUE);
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calibSPD = new AliITSCalibrationSPD();
    spdEntry->Add(calibSPD);
  }
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    //printf(" AliITSOnlineCalibrationSPDhandler::WriteSparseDeadToDB :  nr Sparse dead in module %i - %i \n",module,GetNrSparseDead(module));
    calibSPD->SetNrBadSingle( GetNrSparseDead(module) );
    calibSPD->SetBadList( GetSparseDeadArray(module) );
    for (UInt_t chipIndex=0; chipIndex<5; chipIndex++) {
      UInt_t eq,hs,chip,col,row;
      AliITSRawStreamSPD::OfflineToOnline(module, chipIndex*32, 0, eq, hs, chip, col, row);
      if (IsSilentChip(eq,hs,chip)) {
        calibSPD->SetChipBad(chipIndex);
      }
      else {
        calibSPD->UnSetChipBad(chipIndex);
      }
    }
  }
  AliCDBEntry* cdbEntry = new AliCDBEntry((TObject*)spdEntry,idCalSPD,metaData);
  man->Put(cdbEntry);
  delete spdEntry;
  delete cdbEntry;
  delete metaData;
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::WriteDeadToDBasNoisy(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage) {
  // writes dead pixels to DB for given runNrs
  // overwrites any previous entries
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Henrik Tydesjo");
  metaData->SetComment("Created by AliITSOnlineCalibrationSPDhandler.");
  AliCDBId idCalSPD("ITS/Calib/SPDNoisy",runNrStart,runNrEnd);
  TObjArray* spdEntry = new TObjArray(240);
  spdEntry->SetOwner(kTRUE);
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calibSPD = new AliITSCalibrationSPD();
    spdEntry->Add(calibSPD);
  }
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    calibSPD->SetNrBadSingle( GetNrDeadSingle(module) );
    calibSPD->SetBadList( GetDeadArray(module) );
    for (UInt_t chipIndex=0; chipIndex<5; chipIndex++) {
      UInt_t eq,hs,chip,col,row;
      AliITSRawStreamSPD::OfflineToOnline(module, chipIndex*32, 0, eq, hs, chip, col, row);
      if (IsSilentChip(eq,hs,chip)) {
	calibSPD->SetChipBad(chipIndex);
      }
      else {
	calibSPD->UnSetChipBad(chipIndex);
      }
    }
  }
  AliCDBEntry* cdbEntry = new AliCDBEntry((TObject*)spdEntry,idCalSPD,metaData);
  man->Put(cdbEntry);
  delete spdEntry;
  delete cdbEntry;
  delete metaData;
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::WriteNoisyToDB(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage) {
  // writes noisy pixels to DB for given runNrs
  // overwrites any previous entries
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Henrik Tydesjo");
  metaData->SetComment("Created by AliITSOnlineCalibrationSPDhandler.");
  AliCDBId idCalSPD("ITS/Calib/SPDNoisy",runNrStart,runNrEnd);
  TObjArray* spdEntry = new TObjArray(240);
  spdEntry->SetOwner(kTRUE);
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calibSPD = new AliITSCalibrationSPD();
    spdEntry->Add(calibSPD);
  }
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    calibSPD->SetNrBadSingle( GetNrNoisySingle(module) );
    calibSPD->SetBadList( GetNoisyArray(module) );
  }
  AliCDBEntry* cdbEntry = new AliCDBEntry((TObject*)spdEntry,idCalSPD,metaData);
  man->Put(cdbEntry);
  delete spdEntry;
  delete cdbEntry;
  delete metaData;
  return kTRUE;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::RecursiveInsertDead(AliITSCalibrationSPD* calibSPD, UInt_t module, Int_t lowInd, Int_t highInd) {
  // inserts dead pixels recursively, used when reading from db
  if (lowInd>highInd) return;
  Int_t thisInd = lowInd+(highInd-lowInd)/2;
  SetDeadPixelM(module,calibSPD->GetBadColAt(thisInd),calibSPD->GetBadRowAt(thisInd));
  RecursiveInsertDead(calibSPD,module,lowInd,thisInd-1);
  RecursiveInsertDead(calibSPD,module,thisInd+1,highInd);
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::RecursiveInsertSparseDead(AliITSCalibrationSPD* calibSPD, UInt_t module, Int_t lowInd, Int_t highInd) {
  // inserts sparse dead pixels recursively, used when reading from db
  if (lowInd>highInd) return;
  Int_t thisInd = lowInd+(highInd-lowInd)/2;
  SetSparseDeadPixelM(module,calibSPD->GetBadColAt(thisInd),calibSPD->GetBadRowAt(thisInd));
  RecursiveInsertSparseDead(calibSPD,module,lowInd,thisInd-1);
  RecursiveInsertSparseDead(calibSPD,module,thisInd+1,highInd);
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::RecursiveInsertNoisy(AliITSCalibrationSPD* calibSPD, UInt_t module, Int_t lowInd, Int_t highInd) {
  // inserts noisy pixels recursively, used when reading from db
  if (lowInd>highInd) return;
  Int_t thisInd = lowInd+(highInd-lowInd)/2;
  SetNoisyPixelM(module,calibSPD->GetBadColAt(thisInd),calibSPD->GetBadRowAt(thisInd));
  RecursiveInsertNoisy(calibSPD,module,lowInd,thisInd-1);
  RecursiveInsertNoisy(calibSPD,module,thisInd+1,highInd);
}

#endif
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::GenerateDCSConfigFile(const Char_t* fileName) {
  // generates an ascii file in the format as the one produced by online da (but with dummy runNr=0)
  ofstream dcsfile;
  dcsfile.open(fileName);
  dcsfile << "[SPD SCAN]\n";
  dcsfile << "RunNumber=" << "0" << "\n"; // dummy value
  dcsfile << "Type=" << "4" << "\n";
  dcsfile << "Router=" << "0" << "\n"; // dummy value
  dcsfile << "ActualDetConfiguration=" << "0,-1,-1\n"; // dummy values
  dcsfile << "[NOISY]\n";
  for (UInt_t module=0; module<240; module++) {
    UInt_t headkey=20*10*6;
    for (UInt_t ind=0; ind<GetNrNoisy(module); ind++) {
      UInt_t newkey = GetNoisyEqIdAt(module,ind)*10*6 +
	GetNoisyHSAt(module,ind)*10 +
	GetNoisyChipAt(module,ind);
      if (newkey!=headkey) { // print eq,hs,chip header
	headkey = newkey;
	dcsfile << "-" << newkey/(6*10) << "," << (newkey%(6*10))/10 << "," << (newkey%(6*10))%10 << "\n";
      }
      dcsfile << GetNoisyColAt(module,ind) << "," << GetNoisyRowAt(module,ind) << "\n";
    }
  }
  dcsfile.close();
}
//____________________________________________________________________________________________
TArrayS AliITSOnlineCalibrationSPDhandler::GetSilentArray(UInt_t module, Bool_t treeSerial) {
  // get a TArrayS of the silent=dead+inactive pixels (format for the AliITSCalibrationSPD object)
  // NB! with new implementation of AliITSCalibrationSPD this is not needed anymore
  TArrayS returnArray;

  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t size=0;
  if ( !( IsActiveEq(eq) && IsActiveHS(eq,hs) ) ) {
    size = 8192*5;
  }
  else {
    for (UInt_t ch=0; ch<5; ch++) {
      UInt_t chip = GetChipFromOffline(module,ch*32);
      if (!(IsActiveChip(eq,hs,chip))) {
	size += 8192;
      }
      else {
	UInt_t gloChip = GetGloChip(eq,hs,chip);
	size += fNrDead[gloChip];
      }
    }
  }
  returnArray.Set(size*2);

  UInt_t gloIndex=0;
  if ( !( IsActiveEq(eq) && IsActiveHS(eq,hs) ) ) {
    for (UInt_t colM=0; colM<160; colM++) {
      for (UInt_t rowM=0; rowM<256; rowM++) {
	returnArray.AddAt(colM,gloIndex*2);
	returnArray.AddAt(rowM,gloIndex*2+1);
	gloIndex++;
      }
    }
  }
  else {
    for (UInt_t ch=0; ch<5; ch++) {
      UInt_t chip = GetChipFromOffline(module,ch*32);
      if (!(IsActiveChip(eq,hs,chip))) {
	for (UInt_t colM=ch*32; colM<ch*32+32; colM++) {
	  for (UInt_t rowM=0; rowM<256; rowM++) {
	    returnArray.AddAt(colM,gloIndex*2);
	    returnArray.AddAt(rowM,gloIndex*2+1);
	    gloIndex++;
	  }
	}
      }
      else {
	UInt_t gloChip = GetGloChip(eq,hs,chip);
	if (treeSerial) fDeadPixelMap[gloChip]->PrepareSerialize(); // for tree ordered values
	else fDeadPixelMap[gloChip]->PrepareSerializeOrdered(); // for key ordered values
	for (UInt_t index=0; index<fNrDead[gloChip]; index++) {
	  Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
	  Int_t colM = GetColMFromKey(key);
	  Int_t rowM = GetRowMFromKey(key);
	  returnArray.AddAt(colM,gloIndex*2);
	  returnArray.AddAt(rowM,gloIndex*2+1);
	  gloIndex++;
	}
      }
    }
  }
  return returnArray;
}
//____________________________________________________________________________________________
TArrayS AliITSOnlineCalibrationSPDhandler::GetDeadArray(UInt_t module, Bool_t treeSerial) {
  // get a TArrayS of the single dead pixels (format for the AliITSCalibrationSPD object)
  TArrayS returnArray;

  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t size=GetNrDeadSingle(module);
  returnArray.Set(size*2);
  UInt_t gloIndex=0;
  for (UInt_t ch=0; ch<5; ch++) {
    UInt_t chip = GetChipFromOffline(module,ch*32);
    UInt_t gloChip = GetGloChip(eq,hs,chip);
    if (treeSerial) fDeadPixelMap[gloChip]->PrepareSerialize(); // for tree ordered values
    else fDeadPixelMap[gloChip]->PrepareSerializeOrdered(); // for key ordered values
    if (!IsSilentChip(eq,hs,chip)) {
      for (UInt_t index=0; index<fNrDead[gloChip]; index++) {
	Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
	Int_t colM = GetColMFromKey(key);
	Int_t rowM = GetRowMFromKey(key);
	returnArray.AddAt(colM,gloIndex*2);
	returnArray.AddAt(rowM,gloIndex*2+1);
	gloIndex++;
      }
    }
  }
  return returnArray;
}
//____________________________________________________________________________________________
TArrayS AliITSOnlineCalibrationSPDhandler::GetSparseDeadArray(UInt_t module, Bool_t treeSerial) {
  // get a TArrayS of the single dead pixels (format for the AliITSCalibrationSPD object)
  TArrayS returnArray;

  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t size=GetNrSparseDead(module);
  returnArray.Set(size*2);
  UInt_t gloIndex=0;
  for (UInt_t ch=0; ch<5; ch++) {
    UInt_t chip = GetChipFromOffline(module,ch*32);
    UInt_t gloChip = GetGloChip(eq,hs,chip);
    if (treeSerial) fSparseDeadPixelMap[gloChip]->PrepareSerialize(); // for tree ordered values
    else fSparseDeadPixelMap[gloChip]->PrepareSerializeOrdered(); // for key ordered values
    if (!IsSilentChip(eq,hs,chip)) {
      for (UInt_t index=0; index<fNrSparseDead[gloChip]; index++) {
        Int_t key = fSparseDeadPixelMap[gloChip]->GetKeyIndex(index);
        Int_t colM = GetColMFromKey(key);
        Int_t rowM = GetRowMFromKey(key);
        returnArray.AddAt(colM,gloIndex*2);
        returnArray.AddAt(rowM,gloIndex*2+1);
        gloIndex++;
      }
    }
  }
  return returnArray;
}
//____________________________________________________________________________________________
TArrayS AliITSOnlineCalibrationSPDhandler::GetNoisyArray(UInt_t module, Bool_t treeSerial) {
  // get a TArrayS of the single noisy pixels (format for the AliITSCalibrationSPD object)
  TArrayS returnArray;

  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t size=GetNrNoisySingle(module);
  returnArray.Set(size*2);
  UInt_t gloIndex=0;
  for (UInt_t ch=0; ch<5; ch++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,ch*32));
    if (treeSerial) fNoisyPixelMap[gloChip]->PrepareSerialize(); // for tree ordered values
    else fNoisyPixelMap[gloChip]->PrepareSerializeOrdered(); // for key ordered values
    for (UInt_t index=0; index<fNrNoisy[gloChip]; index++) {
      Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(index);
      Int_t colM = GetColMFromKey(key);
      Int_t rowM = GetRowMFromKey(key);
      returnArray.AddAt(colM,gloIndex*2);
      returnArray.AddAt(rowM,gloIndex*2+1);
      gloIndex++;
    }
  }
  return returnArray;
}
//____________________________________________________________________________________________
TArrayI AliITSOnlineCalibrationSPDhandler::GetDeadArrayOnline(UInt_t eq) {
  // get a TArrayI of the single dead pixels (format for the AliITSOnlineCalibrationSPD object)
  TArrayI returnArray;
  // fix size of array
  UInt_t size=0;
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      UInt_t gloChip = GetGloChip(eq,hs,chip);
      size+=fNrDead[gloChip];
    }
  }
  returnArray.Set(size);
  // put keys in array
  UInt_t gloIndex=0;
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      UInt_t gloChip = GetGloChip(eq,hs,chip);
      fDeadPixelMap[gloChip]->PrepareSerialize(); // for tree ordered values
      for (UInt_t index=0; index<fNrDead[gloChip]; index++) {
	Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
	returnArray.AddAt(key,gloIndex);
	gloIndex++;
      }
    }
  }
  return returnArray;
}
//____________________________________________________________________________________________
TArrayI AliITSOnlineCalibrationSPDhandler::GetNoisyArrayOnline(UInt_t eq) {
  // get a TArrayI of the single noisy pixels (format for the AliITSOnlineCalibrationSPD object)
  TArrayI returnArray;
  // fix size of array
  UInt_t size=0;
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      UInt_t gloChip = GetGloChip(eq,hs,chip);
      size+=fNrNoisy[gloChip];
    }
  }
  returnArray.Set(size);
  // put keys in array
  UInt_t gloIndex=0;
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      UInt_t gloChip = GetGloChip(eq,hs,chip);
      fNoisyPixelMap[gloChip]->PrepareSerialize(); // for tree ordered values
      for (UInt_t index=0; index<fNrNoisy[gloChip]; index++) {
	Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(index);
	returnArray.AddAt(key,gloIndex);
	gloIndex++;
      }
    }
  }
  return returnArray;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::PrintEqSummary() {
// print summary (nr of dead and noisy) for each equipment
  printf("-----------\n");
  printf("Eq summary:\n");
  printf("-----------\n");
  for (UInt_t eq=0; eq<20; eq++) {
    printf("Eq %*d: %*d silent(dead+inactive) , %*d dead , %*d sparse-dead %*d noisy\n",2,eq,6,GetNrSilentEq(eq),6,GetNrDeadEq(eq),6,GetNrSparseDeadEq(eq),6,GetNrNoisyEq(eq));
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::PrintSilent() const {
  // print the inactive and dead pixels to screen
  printf("-----------------------------------------------------------\n");
  printf("Inactive or dead Equipments: (eq  |  module1 .. module12)\n");
  printf("-----------------------------------------------------------\n");
  for (UInt_t eq=0; eq<20; eq++) {
    if (IsSilentEq(eq)) {
      printf("%*d  |  ",2,eq);
      for (UInt_t hs=0; hs<6; hs++) {
	for (UInt_t chip=0; chip<10; chip+=5) {
	  UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
	  if (hs>0 || chip==5) printf(",");
	  printf("%*d",3,module);
	}
      }
      printf("\n");
    }
  }

  printf("-----------------------------------------------------------\n");
  printf("Inactive or dead Half-staves: (eq,hs  |  module1,module2)\n");
  printf("-----------------------------------------------------------\n");
  for (UInt_t eq=0; eq<20; eq++) {
    if (!IsSilentEq(eq)) {
      for (UInt_t hs=0; hs<6; hs++) {
	if (IsSilentHS(eq,hs)) {
	  printf("%*d,%*d  |  ",2,eq,1,hs);
	  for (UInt_t chip=0; chip<10; chip+=5) {
	    UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
	    if (chip==5) printf(",");
	    printf("%*d",3,module);
	  }
	  printf("\n");
	}
      }
    }
  }

  printf("-----------------------------------------------------------\n");
  printf("Inactive or dead Chips: (eq,hs,chip  |  module,colM1-colM2)\n");
  printf("-----------------------------------------------------------\n");
  for (UInt_t eq=0; eq<20; eq++) {
    if (!IsSilentEq(eq)) {
      for (UInt_t hs=0; hs<6; hs++) {
	if (!IsSilentHS(eq,hs)) {
	  for (UInt_t chip=0; chip<10; chip++) {
	    if (IsSilentChip(eq,hs,chip)) {
	      printf("%*d,%*d,%*d  |  ",2,eq,1,hs,1,chip);
	      UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
	      UInt_t colM1 = AliITSRawStreamSPD::GetOfflineColFromOnline(eq,hs,chip,0);
	      UInt_t colM2 = AliITSRawStreamSPD::GetOfflineColFromOnline(eq,hs,chip,31);
	      printf("%*d,%*d-%*d\n",3,module,3,colM1,3,colM2);
	    }
	  }
	}
      }
    }
  }

  PrintDead();

}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::PrintDead() const {
  // print the single dead pixels to screen (disregards inactive eq,hs,chip)
  printf("------------------------------------------------------\n");
  printf("Dead Pixels: (eq,hs,chip,col,row  |  module,colM,rowM)\n");
  printf("------------------------------------------------------\n");
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t index=0; index<fNrDead[gloChip]; index++) {
      Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);

      UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
      UInt_t colM = AliITSRawStreamSPD::GetOfflineColFromOnline(eq,hs,chip,col);
      UInt_t rowM = AliITSRawStreamSPD::GetOfflineRowFromOnline(eq,hs,chip,row);

      printf("%*d,%*d,%*d,%*d,%*d  |  %*d,%*d,%*d\n",2,eq,1,hs,1,chip,2,col,3,row,3,module,3,colM,3,rowM);
    }
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::PrintSparseDead() const {
  // print the single dead pixels to screen (disregards inactive eq,hs,chip)
  printf("------------------------------------------------------\n");
  printf("Sparse Dead Pixels: (eq,hs,chip,col,row  |  module,colM,rowM)\n");
  printf("------------------------------------------------------\n");
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t index=0; index<fNrSparseDead[gloChip]; index++) {
      Int_t key = fSparseDeadPixelMap[gloChip]->GetKeyIndex(index);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);

      UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
      UInt_t colM = AliITSRawStreamSPD::GetOfflineColFromOnline(eq,hs,chip,col);
      UInt_t rowM = AliITSRawStreamSPD::GetOfflineRowFromOnline(eq,hs,chip,row);

      printf("%*d,%*d,%*d,%*d,%*d  |  %*d,%*d,%*d\n",2,eq,1,hs,1,chip,2,col,3,row,3,module,3,colM,3,rowM);
    }
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::PrintNoisy() const {
  // print the dead pixels to screen
  printf("-------------------------------------------------------\n");
  printf("Noisy Pixels: (eq,hs,chip,col,row  |  module,colM,rowM)\n");
  printf("-------------------------------------------------------\n");
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t index=0; index<fNrNoisy[gloChip]; index++) {
      Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(index);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);

      UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
      UInt_t colM = AliITSRawStreamSPD::GetOfflineColFromOnline(eq,hs,chip,col);
      UInt_t rowM = AliITSRawStreamSPD::GetOfflineRowFromOnline(eq,hs,chip,row);

      printf("%*d,%*d,%*d,%*d,%*d  |  %*d,%*d,%*d\n",2,eq,1,hs,1,chip,2,col,3,row,3,module,3,colM,3,rowM);
    }
  }
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::SetDeadPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // set a dead pixel, returns false if pixel is already dead
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::SetDeadPixel", "eq,hs,chip nrs (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  if (col>=32 && row>=256) {
    Error("AliITSOnlineCalibrationSPDhandler::SetDeadPixel", "col,row nrs (%d,%d) out of bounds.",col,row);
    return kFALSE;
  }
  Int_t key = GetKey(eq,hs,chip,col,row);
  // if noisy we dont want to add it...
  if (fNoisyPixelMap[gloChip]->Find(key) != NULL) return kFALSE;
  if (fDeadPixelMap[gloChip]->Insert(key,gloChip)) {
    fNrDead[gloChip]++;
    return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::SetSparseDeadPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // set a dead pixel, returns false if pixel is already dead
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::SetSparseDeadPixel", "eq,hs,chip nrs (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  if (col>=32 && row>=256) {
    Error("AliITSOnlineCalibrationSPDhandler::SetSparseDeadPixel", "col,row nrs (%d,%d) out of bounds.",col,row);
    return kFALSE;
  }
  Int_t key = GetKey(eq,hs,chip,col,row);
  // if noisy we dont want to add it...
  if (fSparseDeadPixelMap[gloChip]->Find(key) != NULL) return kFALSE;
  if (fSparseDeadPixelMap[gloChip]->Insert(key,gloChip)) {
    fNrSparseDead[gloChip]++;
    //printf(" AliITSOnlineCalibrationSPDhandler::SetSparseDeadPixel nSparse Dead : %i \n",fNrSparseDead[gloChip]);	
    return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::SetNoisyPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // set a noisy pixel, returns false if pixel is already noisy
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::SetNoisyPixel", "eq,hs,chip nrs (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  if (col>=32 && row>=256) {
    Error("AliITSOnlineCalibrationSPDhandler::SetNoisyPixel", "col,row nrs (%d,%d) out of bounds.",col,row);
    return kFALSE;
  }
  Int_t key = GetKey(eq,hs,chip,col,row);
  // if dead before - remove from the dead list 
  if (fDeadPixelMap[gloChip]->Remove(key)) {
    fNrDead[gloChip]--;
  }
  if (fNoisyPixelMap[gloChip]->Insert(key,gloChip)) {
    fNrNoisy[gloChip]++;
    return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::SetDeadPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // set a dead pixel, returns false if pixel is already dead
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return SetDeadPixel(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::SetSparseDeadPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // set a dead pixel, returns false if pixel is already dead
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return SetSparseDeadPixel(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::SetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // set a noisy pixel, returns false if pixel is already noisy
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return SetNoisyPixel(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetDeadPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // unset a dead pixel, returns false if pixel is not dead
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::UnSetDeadPixel", "eq,hs,chip nrs (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  Int_t key = GetKey(eq,hs,chip,col,row);
  if (fDeadPixelMap[gloChip]->Remove(key)) {
    fNrDead[gloChip]--;
    return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetSparseDeadPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // unset a dead pixel, returns false if pixel is not dead
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::UnSetSparseDeadPixel", "eq,hs,chip nrs (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  Int_t key = GetKey(eq,hs,chip,col,row);
  if (fSparseDeadPixelMap[gloChip]->Remove(key)) {
    fNrSparseDead[gloChip]--;
    return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetNoisyPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // unset a noisy pixel, returns false if pixel is not noisy
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::UnSetNoisyPixel", "eq,hs,chip nrs (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  Int_t key = GetKey(eq,hs,chip,col,row);
  if (fNoisyPixelMap[gloChip]->Remove(key)) {
    fNrNoisy[gloChip]--;
    return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetDeadPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // unset a dead pixel, returns false if pixel is not dead
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return UnSetDeadPixel(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetSparseDeadPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // unset a dead pixel, returns false if pixel is not dead
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return UnSetSparseDeadPixel(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // unset a noisy pixel, returns false if pixel is not noisy
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return UnSetNoisyPixel(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelBad(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const {
  // is the pixel bad (silent or noisy)
  return (IsPixelSilent(eq,hs,chip,col,row) || IsPixelNoisy(eq,hs,chip,col,row));
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelSilent(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const {
  // is the pixel silent (dead or inactive)?
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200 || col>=32 || row>=256) {
    Error("AliITSOnlineCalibrationSPDhandler::IsPixelSilent", "eq,hs,chip,col,row nrs (%d,%d,%d,%d,%d) out of bounds.",eq,hs,chip,col,row);
    return kFALSE;
  }
  if (IsSilentChip(eq,hs,chip)) return kTRUE;
  else return IsPixelDead(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDead(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const {
  // is the pixel dead?
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200 || col>=32 || row>=256) {
    Error("AliITSOnlineCalibrationSPDhandler::IsPixelDead", "eq,hs,chip,col,row nrs (%d,%d,%d,%d,%d) out of bounds.",eq,hs,chip,col,row);
    return kFALSE;
  }
  UInt_t key = GetKey(eq,hs,chip,col,row);
  if (IsDeadEq(eq) || IsDeadHS(eq,hs) || IsDeadChip(eq,hs,chip)) return kTRUE;
  else {
    if ( fDeadPixelMap[gloChip]->Find(key) != NULL ) return kTRUE;
    else return kFALSE;
  }
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisy(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const {
  // is the pixel noisy?
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200 || col>=32 || row>=256) {
    Error("AliITSOnlineCalibrationSPDhandler::IsPixelNoisy","eq,hs,chip,col,row nrs ( %d, %d, %d, %d, %d ) out of bounds.",eq,hs,chip,col,row);
    return kFALSE;
  }
  UInt_t key = GetKey(eq,hs,chip,col,row);
  if ( fNoisyPixelMap[gloChip]->Find(key) != NULL ) return kTRUE;
  else return kFALSE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelBadM(UInt_t module, UInt_t colM, UInt_t rowM) const {
  // is the pixel bad (silent or noisy)?
  return (IsPixelSilentM(module,colM,rowM) || IsPixelNoisyM(module,colM,rowM));
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelSilentM(UInt_t module, UInt_t colM, UInt_t rowM) const {
  // is the pixel silent (dead or inactive)?
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return IsPixelSilent(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDeadM(UInt_t module, UInt_t colM, UInt_t rowM) const {
  // is the pixel dead?
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return IsPixelDead(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisyM(UInt_t module, UInt_t colM, UInt_t rowM) const  {
  // is the pixel noisy?
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return IsPixelNoisy(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelBadKey(Int_t key) const {
  // is this pixel silent (dead or inactive)?
  UInt_t eq = GetEqIdFromKey(key);
  UInt_t hs = GetHSFromKey(key);
  UInt_t chip = GetChipFromKey(key);
  UInt_t col = GetColFromKey(key);
  UInt_t row = GetRowFromKey(key);
  return IsPixelBad(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelSilentKey(Int_t key) const {
  // is this pixel silent (dead or inactive)?
  UInt_t eq = GetEqIdFromKey(key);
  UInt_t hs = GetHSFromKey(key);
  UInt_t chip = GetChipFromKey(key);
  UInt_t col = GetColFromKey(key);
  UInt_t row = GetRowFromKey(key);
  return IsPixelSilent(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDeadKey(Int_t key) const {
  // is this pixel dead?
  UInt_t eq = GetEqIdFromKey(key);
  UInt_t hs = GetHSFromKey(key);
  UInt_t chip = GetChipFromKey(key);
  UInt_t col = GetColFromKey(key);
  UInt_t row = GetRowFromKey(key);
  return IsPixelDead(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisyKey(Int_t key) const {
  // is this pixel noisy?
  UInt_t eq = GetEqIdFromKey(key);
  UInt_t hs = GetHSFromKey(key);
  UInt_t chip = GetChipFromKey(key);
  UInt_t col = GetColFromKey(key);
  UInt_t row = GetRowFromKey(key);
  return IsPixelNoisy(eq,hs,chip,col,row);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrBad() const {
  // returns the total nr of bad pixels (silent or noisy)
  UInt_t nrBad=0;
  nrBad+=GetNrSilent();
  UInt_t nrNoisy = GetNrNoisy();
  for (UInt_t i=0; i<nrNoisy; i++) {
    UInt_t eq   = GetNoisyEqIdAt(i);
    UInt_t hs   = GetNoisyHSAt(i);
    UInt_t chip = GetNoisyChipAt(i);
    UInt_t col  = GetNoisyColAt(i);
    UInt_t row  = GetNoisyRowAt(i);
    if (!IsPixelSilent(eq,hs,chip,col,row)) nrBad++;
  }
  return nrBad;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSilent() const {
  // returns the total nr of silent pixels (dead or inactive)
  UInt_t nrDead = 0;
  for (UInt_t eq=0; eq<20; eq++) {
    if (IsSilentEq(eq)) {
      nrDead+=81920*6;
      continue;
    }
    for (UInt_t hs=0; hs<6; hs++) {
      if (IsSilentHS(eq,hs)) {
	nrDead+=81920;
	continue;
      }
      for (UInt_t chip=0; chip<10; chip++) {
	if (IsSilentChip(eq,hs,chip)) {
	  nrDead+=8192;
	  continue;
	}
	else {
	  UInt_t gloChip = GetGloChip(eq,hs,chip);
	  nrDead+=fNrDead[gloChip];
	}
      }
    }
  }
  return nrDead;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDead() const {
  // returns the total nr of dead pixels
  UInt_t nrDead = 0;
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    nrDead+=fNrDead[gloChip];
  }
  return nrDead;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSparseDead() const {
  // returns the total nr of dead pixels
  UInt_t nrSparseDead = 0;
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    nrSparseDead+=fNrSparseDead[gloChip];
  }
  return nrSparseDead;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisy() const {
  // returns the total nr of noisy pixels
  UInt_t nrNoisy = 0;
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    nrNoisy+=fNrNoisy[gloChip];
  }
  return nrNoisy;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAt(UInt_t index) const {
  // get eq for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadEqIdAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAt(UInt_t index) const {
  // get eq for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyEqIdAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAt(UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadHSAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAt(UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyHSAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAt(UInt_t index) const {
  // get chip for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadChipAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAt(UInt_t index) const {
  // get chip for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyChipAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAt(UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadColAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAt(UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyColAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAt(UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadRowAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAt(UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyRowAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrBad(UInt_t module) const {
  // returns the number of bad pixels for a certain module (silent or noisy)
  UInt_t nrBad = 0;
  nrBad+=GetNrSilent(module);
  UInt_t nrNoisy = GetNrNoisy(module);
  for (UInt_t i=0; i<nrNoisy; i++) {
    UInt_t eq   = GetNoisyEqIdAt(module,i);
    UInt_t hs   = GetNoisyHSAt(module,i);
    UInt_t chip = GetNoisyChipAt(module,i);
    UInt_t col  = GetNoisyColAt(module,i);
    UInt_t row  = GetNoisyRowAt(module,i);
    if (!IsPixelSilent(eq,hs,chip,col,row)) nrBad++;
  }
  return nrBad;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSilent(UInt_t module) const {
  // returns the number of silent pixels for a certain module (dead or inactive)
  if (module>=240) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrSilent", "module nr (%d) out of bounds.",module);
    return 0;
  }
  UInt_t nrSilent = 0;
  UInt_t eq = GetEqIdFromOffline(module);
  if (IsSilentEq(eq)) return 160*256;
  UInt_t hs = GetHSFromOffline(module);
  if (IsSilentHS(eq,hs)) return 160*256;
  for (UInt_t ch=0; ch<5; ch++) {
    UInt_t chip = GetChipFromOffline(module,ch*32);
    if (IsSilentChip(eq,hs,chip)) {
      nrSilent+=8192;
    }
    else {
      UInt_t gloChip = GetGloChip(eq,hs,chip);
      nrSilent+=fNrDead[gloChip];
    }
  }
  return nrSilent;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDeadSingle(UInt_t module) const {
  // returns the number of single dead pixels (excluding the ones on silent chips) for a certain module
  if (module>=240) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrDeadSingle", "module nr (%d) out of bounds.",module);
    return 0;
  }
  UInt_t nrDead = 0;
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  for (UInt_t ch=0; ch<5; ch++) {
    UInt_t chip = GetChipFromOffline(module,ch*32);
    if (!IsSilentChip(eq,hs,chip)) {
      UInt_t gloChip = GetGloChip(eq,hs,chip);
      nrDead+=fNrDead[gloChip];
    }
  }
  return nrDead;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisySingle(UInt_t module) const {
  // returns the number of noisy pixels for a certain module
  return GetNrNoisy(module);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDead(UInt_t module) const {
  // returns the number of dead pixels for a certain module
  if (module>=240) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrDead", "module nr (%d) out of bounds.",module);
    return 0;
  }
  UInt_t nrDead = 0;
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  for (UInt_t ch=0; ch<5; ch++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,ch*32));
    nrDead+=fNrDead[gloChip];
  }
  return nrDead;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSparseDead(UInt_t module) const {
  // returns the number of sparse dead pixels for a certain module
  if (module>=240) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrSparseDead", "module nr (%d) out of bounds.",module);
    return 0;
  }
  UInt_t nrDead = 0;
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  for (UInt_t ch=0; ch<5; ch++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,ch*32));
    nrDead+=fNrSparseDead[gloChip];
  }
  return nrDead;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisy(UInt_t module) const {
  // returns the number of noisy pixels for a certain module
  if (module>=240) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrNoisy", "module nr (%d) out of bounds.",module);
    return 0;
  }
  UInt_t nrNoisy = 0;
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  for (UInt_t ch=0; ch<5; ch++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,ch*32));
    nrNoisy+=fNrNoisy[gloChip];
  }
  return nrNoisy;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAt(UInt_t module, UInt_t index) const {
  // get eq for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadEqIdAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAt(UInt_t module, UInt_t index) const {
  // get eq for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyEqIdAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAt(UInt_t module, UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadHSAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAt(UInt_t module, UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyHSAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAt(UInt_t module, UInt_t index) const {
  // get chip for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadChipAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAt(UInt_t module, UInt_t index) const {
  // get chip for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyChipAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAt(UInt_t module, UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadColAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAt(UInt_t module, UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyColAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAt(UInt_t module, UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadRowAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAt(UInt_t module, UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyRowAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrBadEq(UInt_t eq) const {
  // returns nr of bad for eq (silent or noisy)
  UInt_t nrBad = 0;
  nrBad+=GetNrSilentEq(eq);
  UInt_t nrNoisy = GetNrNoisy(eq);
  for (UInt_t i=0; i<nrNoisy; i++) {
    UInt_t hs   = GetNoisyHSAtEq(eq,i);
    UInt_t chip = GetNoisyChipAtEq(eq,i);
    UInt_t col  = GetNoisyColAtEq(eq,i);
    UInt_t row  = GetNoisyRowAtEq(eq,i);
    if (!IsPixelSilent(eq,hs,chip,col,row)) nrBad++;
  }
  return nrBad;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSilentEq(UInt_t eq) const {
  // returns nr of silent for eq (dead or inactive)
  UInt_t returnval=0;
  if (IsSilentEq(eq)) return 81920*6;
  for (UInt_t hs=0; hs<6; hs++) {
    if (IsSilentHS(eq,hs)) {
      returnval+=81920;
      continue;
    }
    for (UInt_t chip=0; chip<10; chip++) {
      if (IsSilentChip(eq,hs,chip)) {
	returnval+=8192;
	continue;
      }
      else { 
	returnval+=GetNrDeadC(eq,hs,chip);
      }
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDeadEq(UInt_t eq) const {
  // returns nr of dead for eq
  UInt_t returnval=0;
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      returnval+=GetNrDeadC(eq,hs,chip);
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSparseDeadEq(UInt_t eq) const {
  // returns nr of dead for eq
  UInt_t returnval=0;
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      returnval+=GetNrSparseDeadC(eq,hs,chip);
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisyEq(UInt_t eq) const {
  // returns nr of noisy for eq
  UInt_t returnval=0;
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      returnval+=GetNrNoisyC(eq,hs,chip);
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAtEq(UInt_t eq, UInt_t index) const {
  // get eq for the dead pixel at position index in list of dead
  if (eq<20 && index<GetNrDeadEq(eq)) {
    return eq;
  }
  else {
    return 20;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAtEq(UInt_t eq, UInt_t index) const {
  // get eq for the noisy pixel at position index in list of noisy
  if (eq<20 && index<GetNrNoisyEq(eq)) {
    return eq;
  }
  else {
    return 20;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAtEq(UInt_t eq, UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqDead(eq,index,gloChip,chipIndex);
  return GetDeadHSAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAtEq(UInt_t eq, UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqNoisy(eq,index,gloChip,chipIndex);
  return GetNoisyHSAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAtEq(UInt_t eq, UInt_t index) const {
  // get chip for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqDead(eq,index,gloChip,chipIndex);
  return GetDeadChipAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAtEq(UInt_t eq, UInt_t index) const {
  // get chip for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqNoisy(eq,index,gloChip,chipIndex);
  return GetNoisyChipAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAtEq(UInt_t eq, UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqDead(eq,index,gloChip,chipIndex);
  return GetDeadColAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAtEq(UInt_t eq, UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqNoisy(eq,index,gloChip,chipIndex);
  return GetNoisyColAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAtEq(UInt_t eq, UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqDead(eq,index,gloChip,chipIndex);
  return GetDeadRowAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAtEq(UInt_t eq, UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqNoisy(eq,index,gloChip,chipIndex);
  return GetNoisyRowAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrBadC(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // returns nr of bad for chip (silent or noisy)
  UInt_t nrBad = 0;
  nrBad+=GetNrSilentC(eq,hs,chip);
  UInt_t nrNoisy = GetNrNoisyC(eq,hs,chip);
  for (UInt_t i=0; i<nrNoisy; i++) {
    UInt_t col  = GetNoisyColAtC(eq,hs,chip,i);
    UInt_t row  = GetNoisyRowAtC(eq,hs,chip,i);
    if (!IsPixelSilent(eq,hs,chip,col,row)) nrBad++;
  }
  return nrBad;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSilentC(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // returns nr of silent for chip (dead or inactive)
  if (IsSilentChip(eq,hs,chip)) return 8192;
  else return GetNrDeadC(eq,hs,chip);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDeadC(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // returns nr of dead for chip
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNrDeadC2(gloChip);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSparseDeadC(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // returns nr of sparse dead for chip
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNrSparseDeadC2(gloChip);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisyC(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // returns nr of noisy for chip
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNrNoisyC2(gloChip);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadEqIdAtC2(gloChip,index);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNoisyEqIdAtC2(gloChip,index);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadHSAtC2(gloChip,index);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNoisyHSAtC2(gloChip,index);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadChipAtC2(gloChip,index);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNoisyChipAtC2(gloChip,index);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadColAtC2(gloChip,index);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNoisyColAtC2(gloChip,index);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadRowAtC2(gloChip,index);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNoisyRowAtC2(gloChip,index);
}
//____________________________________________________________________________________________
const Char_t* AliITSOnlineCalibrationSPDhandler::GetDeadPixelAsTextC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  // get a string of dead pixel info
  TString returnMess = "";
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadPixelAsTextC", "global chip nr (%d) out of bounds.",gloChip);
    return returnMess.Data();
  }
  if (index<fNrDead[gloChip]) {
    Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
    UInt_t col = GetColFromKey(key);
    UInt_t row = GetRowFromKey(key);
    UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
    UInt_t colM = AliITSRawStreamSPD::GetOfflineColFromOnline(eq,hs,chip,col);
    UInt_t rowM = AliITSRawStreamSPD::GetOfflineRowFromOnline(eq,hs,chip,row);
    returnMess = Form("%*d,%*d,%*d,%*d,%*d  |  %*d,%*d,%*d",2,eq,1,hs,1,chip,2,col,3,row,3,module,3,colM,3,rowM);
    return returnMess.Data();
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadPixelAsTextC", "Index %d out of bounds.", index);
    return returnMess.Data();
  }
}
//____________________________________________________________________________________________
const Char_t* AliITSOnlineCalibrationSPDhandler::GetNoisyPixelAsTextC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  // get a string of noisy pixel info
  TString returnMess = "";
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyPixelAsTextC", "global chip nr (%d) out of bounds.",gloChip);
    return returnMess.Data();
  }
  if (index<fNrNoisy[gloChip]) {
    Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(index);
    UInt_t col = GetColFromKey(key);
    UInt_t row = GetRowFromKey(key);    
    UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
    UInt_t colM = AliITSRawStreamSPD::GetOfflineColFromOnline(eq,hs,chip,col);
    UInt_t rowM = AliITSRawStreamSPD::GetOfflineRowFromOnline(eq,hs,chip,row);
    returnMess = Form("%*d,%*d,%*d,%*d,%*d  |  %*d,%*d,%*d",2,eq,1,hs,1,chip,2,col,3,row,3,module,3,colM,3,rowM);
    return returnMess.Data();
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyPixelAsTextC", "Index %d out of bounds.", index);
    return returnMess.Data();
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::AddSilentFrom(AliITSOnlineCalibrationSPDhandler* other) {
  // returns number of new silent pixels in this' list (dead or inactive)
  UInt_t tmpdead = GetNrSilent();

  for (UInt_t eq=0; eq<20; eq++) {
    if (!(other->IsActiveEq(eq))) ActivateEq(eq,kFALSE);
    if (other->IsDeadEq(eq))      SetDeadEq(eq,kTRUE);
    for (UInt_t hs=0; hs<6; hs++) {
      if (!(other->IsActiveHS(eq,hs))) ActivateHS(eq,hs,kFALSE);
      if (other->IsDeadHS(eq,hs))      SetDeadHS(eq,hs,kTRUE);
      for (UInt_t chip=0; chip<10; chip++) {
	if (!(other->IsActiveChip(eq,hs,chip))) ActivateChip(eq,hs,chip,kFALSE);
	if (other->IsDeadChip(eq,hs,chip))      SetDeadChip(eq,hs,chip,kTRUE);
      }
    }
  }

  AddDeadFrom(other);

  return GetNrSilent() - tmpdead;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::AddDeadFrom(AliITSOnlineCalibrationSPDhandler* other) {
  // returns number of new dead pixels in this' list
  UInt_t returnval=0;
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t ind1=0; ind1<other->fNrDead[gloChip]; ind1++) {
      UInt_t eq = other->GetDeadEqIdAtC2(gloChip,ind1);
      UInt_t hs   = other->GetDeadHSAtC2(gloChip,ind1);
      UInt_t chip = other->GetDeadChipAtC2(gloChip,ind1);
      UInt_t col  = other->GetDeadColAtC2(gloChip,ind1);
      UInt_t row  = other->GetDeadRowAtC2(gloChip,ind1);
      if (SetDeadPixel(eq,hs,chip,col,row)) {
	returnval++;
      }
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::AddNoisyFrom(AliITSOnlineCalibrationSPDhandler* other) {
  // returns number of new noisy pixels in this' list
  UInt_t returnval=0;
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t ind1=0; ind1<other->fNrNoisy[gloChip]; ind1++) {
      UInt_t eq = other->GetNoisyEqIdAtC2(gloChip,ind1);
      UInt_t hs   = other->GetNoisyHSAtC2(gloChip,ind1);
      UInt_t chip = other->GetNoisyChipAtC2(gloChip,ind1);
      UInt_t col  = other->GetNoisyColAtC2(gloChip,ind1);
      UInt_t row  = other->GetNoisyRowAtC2(gloChip,ind1);
      if (SetNoisyPixel(eq,hs,chip,col,row)) {
	returnval++;
      }
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns nr of dead/noisy in this' lists and not in other's lists (including inactive)
  return GetNrSilentDiff(other) + GetNrNoisyDiff(other);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSilentDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns nr of single silent pixels in this' lists and not in other's lists (dead or inactive)
  UInt_t returnval=0;
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	if (!IsActiveEq(eq) || IsDeadEq(eq) || !IsActiveHS(eq,hs) || IsDeadHS(eq,hs) || !IsActiveChip(eq,hs,chip) || IsDeadChip(eq,hs,chip)) {
	  if (other->IsActiveEq(eq) && !other->IsDeadEq(eq) && other->IsActiveHS(eq,hs) && !other->IsDeadHS(eq,hs) && other->IsActiveChip(eq,hs,chip) && !other->IsDeadChip(eq,hs,chip)) {
	    // if this is inactive and the other is active...
	    returnval+= 8192 - other->GetNrDeadC(eq,hs,chip);
	  }
	}
	else {
	  UInt_t gloChip = GetGloChip(eq,hs,chip);
	  for (UInt_t ind1=0; ind1<fNrDead[gloChip]; ind1++) {
	    Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(ind1);
	    UInt_t col = GetColFromKey(key);
	    UInt_t row = GetRowFromKey(key);
	    if (!(other->IsPixelSilent(eq,hs,chip,col,row))) {
	      returnval++;
	    }
	  }
	}
      }
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns nr of dead in this' lists and not in other's lists
  UInt_t returnval=0;
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t ind1=0; ind1<fNrDead[gloChip]; ind1++) {
      Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(ind1);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelDead(eq,hs,chip,col,row))) {
	returnval++;
      }
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns nr of noisy in this' lists and not in other's lists
  UInt_t returnval=0;
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t ind1=0; ind1<fNrNoisy[gloChip]; ind1++) {
      Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(ind1);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelNoisy(eq,hs,chip,col,row))) {
	returnval++;
      }
    }
  }
  return returnval;
}
//____________________________________________________________________________________________
AliITSOnlineCalibrationSPDhandler* AliITSOnlineCalibrationSPDhandler::GetDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns handler with active/dead/noisy in this' lists, removing those that are in other's lists
  AliITSOnlineCalibrationSPDhandler* newHandler = new AliITSOnlineCalibrationSPDhandler();

  for (UInt_t eq=0; eq<20; eq++) {
    if (!(IsActiveEq(eq))) {
      newHandler->ActivateEq(eq,kFALSE);
      if (!other->IsActiveEq(eq)) newHandler->ActivateEq(eq);
    }
    if (IsDeadEq(eq)) {
      newHandler->SetDeadEq(eq);
      if (other->IsDeadEq(eq)) newHandler->SetDeadEq(eq,kFALSE);
    }
    for (UInt_t hs=0; hs<6; hs++) {
      if (!(IsActiveHS(eq,hs))) {
	newHandler->ActivateHS(eq,hs,kFALSE);
	if (!other->IsActiveHS(eq,hs)) newHandler->ActivateHS(eq,hs);
      }
      if (IsDeadHS(eq,hs)) {
	newHandler->SetDeadHS(eq,hs);
	if (other->IsDeadHS(eq,hs)) newHandler->SetDeadHS(eq,hs,kFALSE);
      }
      for (UInt_t chip=0; chip<10; chip++) {
	if (!(IsActiveChip(eq,hs,chip))) {
	  newHandler->ActivateChip(eq,hs,chip,kFALSE);
	  if (!other->IsActiveChip(eq,hs,chip)) newHandler->ActivateHS(eq,hs,chip);
	}
	if (IsDeadChip(eq,hs,chip)) {
	  newHandler->SetDeadChip(eq,hs,chip);
	  if (other->IsDeadChip(eq,hs,chip)) newHandler->SetDeadChip(eq,hs,chip,kFALSE);
	}
      }
    }
  }

  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t ind1=0; ind1<fNrDead[gloChip]; ind1++) {
      Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(ind1);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelDead(eq,hs,chip,col,row))) {
	newHandler->SetDeadPixel(eq,hs,chip,col,row);
      }
    }
  }

  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t ind2=0; ind2<fNrNoisy[gloChip]; ind2++) {
      Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(ind2);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelNoisy(eq,hs,chip,col,row))) {
	newHandler->SetNoisyPixel(eq,hs,chip,col,row);
      }
    }
  }

  return newHandler;
}
//____________________________________________________________________________________________
AliITSOnlineCalibrationSPDhandler* AliITSOnlineCalibrationSPDhandler::GetSilentDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns handler with active/dead in this' lists, removing those that are in other's lists
  AliITSOnlineCalibrationSPDhandler* newHandler = new AliITSOnlineCalibrationSPDhandler();

  for (UInt_t eq=0; eq<20; eq++) {
    if (!(IsActiveEq(eq))) {
      newHandler->ActivateEq(eq,kFALSE);
      if (!other->IsActiveEq(eq)) newHandler->ActivateEq(eq);
    }
    if (IsDeadEq(eq)) {
      newHandler->SetDeadEq(eq);
      if (other->IsDeadEq(eq)) newHandler->SetDeadEq(eq,kFALSE);
    }
    for (UInt_t hs=0; hs<6; hs++) {
      if (!(IsActiveHS(eq,hs))) {
	newHandler->ActivateHS(eq,hs,kFALSE);
	if (!other->IsActiveHS(eq,hs)) newHandler->ActivateHS(eq,hs);
      }
      if (IsDeadHS(eq,hs)) {
	newHandler->SetDeadHS(eq,hs);
	if (other->IsDeadHS(eq,hs)) newHandler->SetDeadHS(eq,hs,kFALSE);
      }
      for (UInt_t chip=0; chip<10; chip++) {
	if (!(IsActiveChip(eq,hs,chip))) {
	  newHandler->ActivateChip(eq,hs,chip,kFALSE);
	  if (!other->IsActiveChip(eq,hs,chip)) newHandler->ActivateHS(eq,hs,chip);
	}
	if (IsDeadChip(eq,hs,chip)) {
	  newHandler->SetDeadChip(eq,hs,chip);
	  if (other->IsDeadChip(eq,hs,chip)) newHandler->SetDeadChip(eq,hs,chip,kFALSE);
	}
      }
    }
  }

  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t ind1=0; ind1<fNrDead[gloChip]; ind1++) {
      Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(ind1);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelDead(eq,hs,chip,col,row))) {
	newHandler->SetDeadPixel(eq,hs,chip,col,row);
      }
    }
  }

  return newHandler;
}
//____________________________________________________________________________________________
AliITSOnlineCalibrationSPDhandler* AliITSOnlineCalibrationSPDhandler::GetDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns handler with dead in this' lists, except for those in other's lists
  AliITSOnlineCalibrationSPDhandler* newHandler = new AliITSOnlineCalibrationSPDhandler();
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t ind1=0; ind1<fNrDead[gloChip]; ind1++) {
      Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(ind1);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelDead(eq,hs,chip,col,row))) {
	newHandler->SetDeadPixel(eq,hs,chip,col,row);
      }
    }
  }
  return newHandler;
}
//____________________________________________________________________________________________
AliITSOnlineCalibrationSPDhandler* AliITSOnlineCalibrationSPDhandler::GetNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns handler with noisy in this' lists, except for those in other's lists
  AliITSOnlineCalibrationSPDhandler* newHandler = new AliITSOnlineCalibrationSPDhandler();
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    for (UInt_t ind2=0; ind2<fNrNoisy[gloChip]; ind2++) {
      Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(ind2);
      UInt_t eq = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelNoisy(eq,hs,chip,col,row))) {
	newHandler->SetNoisyPixel(eq,hs,chip,col,row);
      }
    }
  }
  return newHandler;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexDead(UInt_t module, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const {
  // find gloChip and chipIndex from module and index
  if (index<GetNrDead(module)) {
    UInt_t eq = GetEqIdFromOffline(module);
    UInt_t hs = GetHSFromOffline(module);
    
    UInt_t glVal=0;
    for (UInt_t ch=0; ch<5; ch++) {
      gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,ch*32));
      if (glVal+fNrDead[gloChip]>index) {
	chipIndex = index-glVal;
	break;
      }
      else {
	glVal+=fNrDead[gloChip];
      }
    }

  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetChipAndIndexDead", "Index %d out of bounds.", index);
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexNoisy(UInt_t module, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const {
  // find gloChip and chipIndex from module and index
  if (index<GetNrNoisy(module)) {
    UInt_t eq = GetEqIdFromOffline(module);
    UInt_t hs = GetHSFromOffline(module);
    
    UInt_t glVal=0;
    for (UInt_t ch=0; ch<5; ch++) {
      gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,ch*32));
      if (glVal+fNrNoisy[gloChip]>index) {
	chipIndex = index-glVal;
	break;
      }
      else {
	glVal+=fNrNoisy[gloChip];
      }
    }

  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetChipAndIndexNoisy", "Index %d out of bounds.", index);
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexEqDead(UInt_t eq, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const {
  // find gloChip and chipIndex from module and index
  if (index<GetNrDeadEq(eq)) {

    UInt_t glVal=0;
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	gloChip = GetGloChip(eq,hs,chip);
	if (glVal+fNrDead[gloChip]>index) {
	  chipIndex = index-glVal;
	  break;
	}
	else {
	  glVal+=fNrDead[gloChip];
	}
      }
    }

  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetChipAndIndexEqDead", "Index %d out of bounds.", index);
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexEqNoisy(UInt_t eq, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const {
  // find gloChip and chipIndex from module and index
  if (index<GetNrNoisyEq(eq)) {

    UInt_t glVal=0;
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	gloChip = GetGloChip(eq,hs,chip);
	if (glVal+fNrNoisy[gloChip]>index) {
	  chipIndex = index-glVal;
	  break;
	}
	else {
	  glVal+=fNrNoisy[gloChip];
	}
      }
    }

  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetChipAndIndexEqNoisy", "Index %d out of bounds.", index);
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexTotDead(UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const {
  // find gloChip and chipIndex from global index
  if (index<GetNrDead()) {
    
    UInt_t glVal=0;
    for (gloChip=0; gloChip<1200; gloChip++) {
      if (glVal+fNrDead[gloChip]>index) {
	chipIndex = index-glVal;
	break;
      }
      else {
	glVal+=fNrDead[gloChip];
      }
    }

  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetChipAndIndexTotDead", "Index %d out of bounds.", index);
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexTotNoisy(UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const {
  // find gloChip and chipIndex from global index
  if (index<GetNrNoisy()) {
    
    UInt_t glVal=0;
    for (gloChip=0; gloChip<1200; gloChip++) {
      if (glVal+fNrNoisy[gloChip]>index) {
	chipIndex = index-glVal;
	break;
      }
      else {
	glVal+=fNrNoisy[gloChip];
      }
    }

  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetChipAndIndexTotNoisy", "Index %d out of bounds.", index);
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetEqIdFromOffline(UInt_t module) const {
  // module to eq mapping
  if (module>=240) {
    Error("AliITSOnlineCalibrationSPDhandler::GetEqIdFromOffline", "module nr (%d) out of bounds.",module);
    return 20;
  }
  return AliITSRawStreamSPD::GetOnlineEqIdFromOffline(module);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetHSFromOffline(UInt_t module) const {
  // module to hs mapping
  if (module>=240) {
    Error("AliITSOnlineCalibrationSPDhandler::GetHSFromOffline", "module nr (%d) out of bounds.",module);
    return 6;
  }
  return AliITSRawStreamSPD::GetOnlineHSFromOffline(module);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetChipFromOffline(UInt_t module, UInt_t colM) const {
  // module,colM to chip mapping
  if (module>=240 || colM>=160) {
    Error("AliITSOnlineCalibrationSPDhandler::GetChipFromOffline", "module,colM nrs (%d,%d) out of bounds.",module,colM);
    return 10;
  }
  return AliITSRawStreamSPD::GetOnlineChipFromOffline(module,colM);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetColFromOffline(UInt_t module, UInt_t colM) const {
  // colM to col mapping
  if (colM>=160) {
    Error("AliITSOnlineCalibrationSPDhandler::GetColFromOffline", "colM nr (%d) out of bounds.",colM);
    return 160;
  }
  return AliITSRawStreamSPD::GetOnlineColFromOffline(module,colM);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetRowFromOffline(UInt_t module, UInt_t rowM) const {
  // rowM to row mapping
  if (rowM>=256) {
    Error("AliITSOnlineCalibrationSPDhandler::GetRowFromOffline", "rowM nr (%d) out of bounds.",rowM);
    return 256;
  }
  return AliITSRawStreamSPD::GetOnlineRowFromOffline(module,rowM);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDeadC2(UInt_t gloChip) const {
  // returns nr of dead pixels on this chip
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrDeadC2", "global chip nr (%d) out of bounds.",gloChip);
    return 0;
  }
  return fNrDead[gloChip];
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrSparseDeadC2(UInt_t gloChip) const {
  // returns nr of dead pixels on this chip
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrSparseDeadC2", "global chip nr (%d) out of bounds.",gloChip);
    return 0;
  }
  return fNrSparseDead[gloChip];
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisyC2(UInt_t gloChip) const {
  // returns nr of noisy pixels on this chip
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrNoisyC2", "global chip nr (%d) out of bounds.",gloChip);
    return 0;
  }
  return fNrNoisy[gloChip];
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAtC2(UInt_t gloChip, UInt_t index) const {
  // get eq for the dead pixel at position index in list of dead
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAtC", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrDead[gloChip]) {
    Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
    return GetEqIdFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAtC", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAtC2(UInt_t gloChip, UInt_t index) const {
  // get eq for the noisy pixel at position index in list of noisy
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAtC", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrNoisy[gloChip]) {
    Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(index);
    return GetEqIdFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAtC", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAtC2(UInt_t gloChip, UInt_t index) const {
  // get hs for the dead pixel at position index in list of dead
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadHSAtC", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrDead[gloChip]) {
    Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
    return GetHSFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadHSAtC", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAtC2(UInt_t gloChip, UInt_t index) const {
  // get hs for the noisy pixel at position index in list of noisy
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyHSAtC", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrNoisy[gloChip]) {
    Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(index);
    return GetHSFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyHSAtC", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAtC2(UInt_t gloChip, UInt_t index) const {
  // get chip for the dead pixel at position index in list of dead
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadChipAtC", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrDead[gloChip]) {
    Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
    return GetChipFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadChipAtC", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAtC2(UInt_t gloChip, UInt_t index) const {
  // get chip for the noisy pixel at position index in list of noisy
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyChipAtC", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrNoisy[gloChip]) {
    Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(index);
    return GetChipFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyChipAtC", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAtC2(UInt_t gloChip, UInt_t index) const {
  // get col for the dead pixel at position index in list of dead
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadColAtC2", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrDead[gloChip]) {
    Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
    return GetColFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadColAtC2", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAtC2(UInt_t gloChip, UInt_t index) const {
  // get col for the noisy pixel at position index in list of noisy
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyColAtC", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrNoisy[gloChip]) {
    Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(index);
    return GetColFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyColAtC", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAtC2(UInt_t gloChip, UInt_t index) const {
  // get row for the dead pixel at position index in list of dead
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadRowAtC", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrDead[gloChip]) {
    Int_t key = fDeadPixelMap[gloChip]->GetKeyIndex(index);
    return GetRowFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadRowAtC", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAtC2(UInt_t gloChip, UInt_t index) const {
  // get row for the noisy pixel at position index in list of noisy
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyRowAtC", "global chip nr (%d) out of bounds.",gloChip);
    return 20;
  }
  if (index<fNrNoisy[gloChip]) {
    Int_t key = fNoisyPixelMap[gloChip]->GetKeyIndex(index);
    return GetRowFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyRowAtC", "Index %d out of bounds.", index);
    return 0;
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ActivateALL() {
  // activate all eq,hs,chips
  for (UInt_t eq=0; eq<20; eq++) {
    ActivateEq(eq);

    for (UInt_t hs=0; hs<6; hs++) {
      ActivateHS(eq,hs);
      for (UInt_t chip=0; chip<10; chip++) {
	ActivateChip(eq,hs,chip);
      }
    }
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ActivateEq(UInt_t eq, Bool_t setval) {
  // activate eq
  if (eq>=20) {
    Error("AliITSOnlineCalibrationSPDhandler::ActivateEq", "eq (%d) out of bounds.",eq);
    return;
  }
  fActiveEq[eq] = setval;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ActivateHS(UInt_t eq, UInt_t hs, Bool_t setval) {
  // activate hs
  if (eq>=20 || hs>=6) {
    Error("AliITSOnlineCalibrationSPDhandler::ActivateHS", "eq,hs (%d,%d) out of bounds.",eq,hs);
    return;
  }
  fActiveHS[eq][hs] = setval;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::ActivateChip(UInt_t eq, UInt_t hs, UInt_t chip, Bool_t setval) {
  // activate chip
  if (eq>=20 || hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPDhandler::ActivateChip", "eq,hs,chip (%d,%d,%d) out of bounds.",eq,hs,chip);
    return;
  }
  fActiveChip[eq][hs][chip] = setval;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsActiveEq(UInt_t eq) const {
  // Is eq active?
  if (eq>=20) {
    Error("AliITSOnlineCalibrationSPDhandler::IsActiveEq", "eq (%d) out of bounds.",eq);
    return kFALSE;
  }
  return fActiveEq[eq];
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsActiveHS(UInt_t eq, UInt_t hs) const {
  // Is hs active?
  if (eq>=20 || hs>=6) {
    Error("AliITSOnlineCalibrationSPDhandler::IsActiveHS", "eq,hs (%d,%d) out of bounds.",eq,hs);
    return kFALSE;
  }
  return fActiveHS[eq][hs];
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsActiveChip(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // Is chip active?
  if (eq>=20 || hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPDhandler::IsActiveChip", "eq,hs,chip (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  return fActiveChip[eq][hs][chip];
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::UnSetDeadALL() {
  // Clear all dead eq,hs,chips
  for (UInt_t eq=0; eq<20; eq++) {
    SetDeadEq(eq,kFALSE);
    for (UInt_t hs=0; hs<6; hs++) {
      SetDeadHS(eq,hs,kFALSE);
      for (UInt_t chip=0; chip<10; chip++) {
	SetDeadChip(eq,hs,chip,kFALSE);
      }
    }
  }
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::SetDeadEq(UInt_t eq, Bool_t setval) {
  // set eq dead
  if (eq>=20) {
    Error("AliITSOnlineCalibrationSPDhandler::SetDeadEq", "eq (%d) out of bounds.",eq);
    return;
  }
  fDeadEq[eq] = setval;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::SetDeadHS(UInt_t eq, UInt_t hs, Bool_t setval) {
  // set hs dead
  if (eq>=20 || hs>=6) {
    Error("AliITSOnlineCalibrationSPDhandler::SetDeadHS", "eq,hs (%d,%d) out of bounds.",eq,hs);
    return;
  }
  fDeadHS[eq][hs] = setval;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::SetDeadChip(UInt_t eq, UInt_t hs, UInt_t chip, Bool_t setval) {
  // set chip dead
  if (eq>=20 || hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPDhandler::SetDeadChip", "eq,hs,chip (%d,%d,%d) out of bounds.",eq,hs,chip);
    return;
  }
  fDeadChip[eq][hs][chip] = setval;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsDeadEq(UInt_t eq) const {
  // is eq dead?
  if (eq>=20) {
    Error("AliITSOnlineCalibrationSPDhandler::IsDeadEq", "eq (%d) out of bounds.",eq);
    return kFALSE;
  }
  return fDeadEq[eq];
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsDeadHS(UInt_t eq, UInt_t hs) const {
  // is hs dead?
  if (eq>=20 || hs>=6) {
    Error("AliITSOnlineCalibrationSPDhandler::IsDeadHS", "eq,hs (%d,%d) out of bounds.",eq,hs);
    return kFALSE;
  }
  return fDeadHS[eq][hs];
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsDeadChip(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // is chip dead?
  if (eq>=20 || hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPDhandler::IsDeadChip", "eq,hs,chip (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  return fDeadChip[eq][hs][chip];
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsSilentEq(UInt_t eq) const {
  // is eq silent?
  if (eq>=20) {
    Error("AliITSOnlineCalibrationSPDhandler::IsSilentEq", "eq (%d) out of bounds.",eq);
    return kFALSE;
  }
  return (!IsActiveEq(eq) || IsDeadEq(eq));
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsSilentHS(UInt_t eq, UInt_t hs) const {
  // is hs silent?
  if (eq>=20 || hs>=6) {
    Error("AliITSOnlineCalibrationSPDhandler::IsSilentHS", "eq,hs (%d,%d) out of bounds.",eq,hs);
    return kFALSE;
  }
  return (!IsActiveEq(eq) || IsDeadEq(eq) || !IsActiveHS(eq,hs) || IsDeadHS(eq,hs));
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsSilentChip(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // is chip silent?
  if (eq>=20 || hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPDhandler::IsSilentChip", "eq,hs,chip (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  return (!IsActiveEq(eq) || IsDeadEq(eq) || !IsActiveHS(eq,hs) || IsDeadHS(eq,hs) || !IsActiveChip(eq,hs,chip) || IsDeadChip(eq,hs,chip));
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsNoisyChip(UInt_t eq, UInt_t hs, UInt_t chip) const {
  // returns true if there is at least a noisy pixel in the chip

  if (eq>=20 || hs>=6 || chip>=10) {
    Error("AliITSOnlineCalibrationSPDhandler::IsNoisyChip", "eq,hs,chip (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  Bool_t isNoisy = kFALSE;

  UInt_t nrNoisy = GetNrNoisy();
  for (UInt_t i=0; i<nrNoisy; i++) {
    if(eq  == GetNoisyEqIdAt(i)){
      if(hs  == GetNoisyHSAt(i)){
        if(chip == GetNoisyChipAt(i)) {
            UInt_t col  = GetNoisyColAt(i);
            UInt_t row  = GetNoisyRowAt(i);
            if (IsPixelNoisy(eq,hs,chip,col,row)) isNoisy = kTRUE;
        }
      }
    }
  }
  return isNoisy;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::WritePITConditionsToDB(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage) {
  // writes noisy pixels to DB for given runNrs
  // overwrites any previous entries
  AliCDBManager* man = AliCDBManager::Instance();
  TString storageSTR = Form("%s",storage);
  if (storageSTR.CompareTo("default")==0) {
    if(!man->IsDefaultStorageSet()) {
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    }
  }
  else {
    storageSTR = Form("%s",storage);
    man->SetDefaultStorage(storageSTR.Data());
  }
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Annalisa Mastroserio");
  metaData->SetComment("Created by AliITSOnlineCalibrationSPDhandler.");
  AliCDBId idCalSPD("TRIGGER/SPD/PITConditions",runNrStart,runNrEnd);
  AliCDBEntry* cdbEntry = new AliCDBEntry((TObject*)fTriggerConditions,idCalSPD,metaData);
  man->Put(cdbEntry);
  delete cdbEntry;
  delete metaData;
  return kTRUE;
}

//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::SetInactiveChipInPITmask(UInt_t eq, UInt_t hs, UInt_t chip){
  //
  fTriggerConditions->SetInActiveChip(eq,hs,chip);
  return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetInactiveChipInPITmask(UInt_t eq, UInt_t hs, UInt_t chip){
  //
  fTriggerConditions->SetInActiveChip(eq,hs,chip);
  return kTRUE;
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::PrintDiffInDead(AliITSOnlineCalibrationSPDhandler *other) const {
  //
  // Printout of the differences between two ocdb files for SPD Dead channel map
  //
  UInt_t nrChipOk=0;
  UInt_t nrDeadChipOk=0; 
  UInt_t nrDeadHsOk=0;
  UInt_t nrDeadHs =0;
  UInt_t nrDeadChip=0;
  UInt_t nrDeadHsInOther =0;
  UInt_t nrDeadChipInOther=0;
  UInt_t nrMismatch =0;
  UInt_t nrMismatchInOther =0;
  printf("\n\n ****** loop over chips ***** \n");
  for(Int_t eq=0; eq<20; eq++){
   if(TMath::Abs((Int_t)GetNrBadEq(eq) - (Int_t)other->GetNrBadEq(eq)) >0) printf("-----> dead pixels differ in eq %i!   %i - %i in the other \n",eq,GetNrBadEq(eq),other->GetNrBadEq(eq));
   for(Int_t hs=0; hs<6; hs++){
    Short_t nchips =0;
    Short_t nchipsOther =0;
    Short_t nok=0;
    for(Int_t chip=0; chip<10; chip++){
      UInt_t chipkey = AliITSRawStreamSPD::GetOfflineChipKeyFromOnline(eq,hs,chip);
     // test if everything is coherent
     if(IsDeadChip(eq,hs,chip) && other->IsDeadChip(eq,hs,chip)) {
      nok++;
      nrChipOk++;
      nrDeadChipOk++;
     }
     if(!IsDeadChip(eq,hs,chip) && !other->IsDeadChip(eq,hs,chip)) nrChipOk++;
     // now testing if mismatches
     if(IsDeadChip(eq,hs,chip)) {
      nrDeadChip++;
      nchips++;
      if(!other->IsDeadChip(eq,hs,chip)) {
        nrMismatch++;
        printf("  mismatch -> eq %i  hs %i  chip %i is DEAD  - ALIVE in the other (chipkey %i)\n",eq,hs,chip,chipkey);
       }
      }
     if(other->IsDeadChip(eq,hs,chip)){
      nrDeadChipInOther++;
      nchipsOther++;
      if(!IsDeadChip(eq,hs,chip)) {
       nrMismatchInOther++;
       printf("  mismatch -> eq %i  hs %i  chip %i is ALIVE -  DEAD in the other (chipkey %i)\n",eq,hs,chip,chipkey);
      }
     }
    }
    if(nok==10) nrDeadHsOk++;
    if(nchips==10) nrDeadHs++;
    if(nchipsOther==10) nrDeadHsInOther++;
   }
  }

printf("\n\n\n*************SUMMARY****************\n");
printf(" BOTH have : %i Dead HS and %i Dead chips  with %i coherent chips \n",nrDeadHsOk,nrDeadChipOk,nrChipOk);
printf("\n_________MISMATCH RESULTS___________\n");
printf(" THIS  : Nr Dead HS %i - Nr Dead Chip %i \n",nrDeadHs,nrDeadChip);
printf(" OTHER : Nr Dead HS %i - Nr Dead Chip %i \n",nrDeadHsInOther,nrDeadChipInOther);
printf(" N Mismatches in Dead  chips (=ALIVE in the other) %i \n",nrMismatch);
printf(" N Mismatches in Alive chips (=DEAD  in the other) %i \n",nrMismatchInOther);
}
//____________________________________________________________________________________________
void AliITSOnlineCalibrationSPDhandler::PrintDiffInPITmask(AliITSOnlineCalibrationSPDhandler *other) const {
  //
  // Printout of the differences between two ocdb files for SPD Dead channel map
  //

Int_t nOk =0;
Int_t nMismatch =0;
Int_t nMismatchInOther =0;

printf("\n\n ****** loop over chips in PIT mask***** \n");
for(Int_t eq=0; eq<20; eq++){
  for(Int_t hs=0; hs<6; hs++){
   for(Int_t chip=0; chip<10; chip++){

  UInt_t chipkey = AliITSRawStreamSPD::GetOfflineChipKeyFromOnline(eq,hs,chip);

  if(fTriggerConditions->IsChipActive(eq,hs,chip) && (other->GetTriggerConditions())->IsChipActive(eq,hs,chip)) nOk++;
  if(fTriggerConditions->IsChipActive(eq,hs,chip) && !(other->GetTriggerConditions())->IsChipActive(eq,hs,chip)) {
   nMismatch++;
   printf("Mismatch -> eq %i  hs %i chip %i is ACTIVE - INACTIVE in the other (chipkey %i) \n",eq,hs,chip,chipkey);
  }
  if(!fTriggerConditions->IsChipActive(eq,hs,chip) && (other->GetTriggerConditions())->IsChipActive(eq,hs,chip)) {
   nMismatchInOther++;
   printf("Mismatch -> eq %i  hs %i chip %i is INACTIVE - ACTIVE in the other (chipkey %i) \n",eq,hs,chip,chipkey);
   }
  if(!fTriggerConditions->IsChipActive(eq,hs,chip) && !(other->GetTriggerConditions())->IsChipActive(eq,hs,chip)) nOk++;
  }
 }
}

printf("n Chips OK %i : ACTIVE mismatch %i  - INACTIVE mismatch in %i \n",nOk,nMismatch,nMismatchInOther);

}

