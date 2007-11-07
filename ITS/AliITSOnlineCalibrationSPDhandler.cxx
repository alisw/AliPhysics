//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// For easier handling of dead and noisy pixels they are kept in    //
// container maps (AliITSIntMap).                                   //
// The TArrayI objects that are put in the AliITSCalibrationSPD     //
// objects can be obtained from the methods GetDeadArray and        //
// GetNoisyArray.                                                   //
//////////////////////////////////////////////////////////////////////   

#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliITSOnlineCalibrationSPD.h"
#include "AliITSIntMap.h"
#include <TObjArray.h>
#include <TArrayI.h>
#include <TFile.h>
#include <TError.h>
#include <fstream>

#ifndef SPD_DA_OFF // you may want to exclude cdb classes, if using this class outside aliroot
#include "AliITSCalibrationSPD.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#endif

//____________________________________________________________________________________________
AliITSOnlineCalibrationSPDhandler::AliITSOnlineCalibrationSPDhandler():
  fFileLocation(".")
{
  // constructor
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    fNrDead[gloChip]=0;
    fNrNoisy[gloChip]=0;
    fDeadPixelMap[gloChip] = new AliITSIntMap();
    fNoisyPixelMap[gloChip] = new AliITSIntMap();
  }
}
AliITSOnlineCalibrationSPDhandler::AliITSOnlineCalibrationSPDhandler(const AliITSOnlineCalibrationSPDhandler& handle): 
  fFileLocation(".")
{
  // copy constructor
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    fNrDead[gloChip] = handle.fNrDead[gloChip];
    fNrNoisy[gloChip] = handle.fNrNoisy[gloChip];
    fDeadPixelMap[gloChip] = handle.fDeadPixelMap[gloChip]->Clone();
    fNoisyPixelMap[gloChip] = handle.fNoisyPixelMap[gloChip]->Clone();
  }
  fFileLocation = handle.fFileLocation;
}
AliITSOnlineCalibrationSPDhandler::~AliITSOnlineCalibrationSPDhandler() {
  //  ClearMaps();
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    delete fDeadPixelMap[gloChip];
    delete fNoisyPixelMap[gloChip];
  }
}
AliITSOnlineCalibrationSPDhandler& AliITSOnlineCalibrationSPDhandler::operator=(const AliITSOnlineCalibrationSPDhandler& handle) {
  // assignment operator
  if (this!=&handle) {
    this->ClearMaps();
    for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
      fNrDead[gloChip] = handle.fNrDead[gloChip];
      fNrNoisy[gloChip] = handle.fNrNoisy[gloChip];
      fDeadPixelMap[gloChip] = handle.fDeadPixelMap[gloChip]->Clone();
      fNoisyPixelMap[gloChip] = handle.fNoisyPixelMap[gloChip]->Clone();
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
}
void AliITSOnlineCalibrationSPDhandler::ResetDead() {
  // reset the dead pixel map
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    fNrDead[gloChip]=0;
    fDeadPixelMap[gloChip]->Clear();
  }
}
void AliITSOnlineCalibrationSPDhandler::ResetNoisy() {
  // clear the list of noisy pixels
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    fNrNoisy[gloChip]=0;
    fNoisyPixelMap[gloChip]->Clear();
  }
}
void AliITSOnlineCalibrationSPDhandler::ResetDeadForChip(UInt_t eq, UInt_t hs, UInt_t chip) {
  // clear the dead pixels for this chip
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      Int_t key = GetKey(eq,hs,chip,col,row);
      if (fDeadPixelMap[gloChip]->Remove(key)) {
	fNrDead[gloChip]--;
      }
    }
  }
}
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
  // read dead and noisy files from file location. returns true if at least one file found
  Bool_t b1 = ReadNoisyFromFiles();
  Bool_t b2 = ReadDeadFromFiles();
  return (b1 && b2);
}

Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromFiles() {
  // read dead files from file location. returns true if at least one file found
  Bool_t returnval=kFALSE;
  for (UInt_t eq=0; eq<20; eq++) {
    if (ReadDeadFromFile(eq)) {
      returnval=kTRUE;
    }
  }
  return returnval;
}
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
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromFile(UInt_t eq) {
  // read dead file for module from file location. 
  TString fileName = Form("%s/SPD_Dead_%d.root",fFileLocation.Data(),eq);
  return ReadDeadFromFileName(fileName.Data());
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromFile(UInt_t eq) {
  // read noisy file for module from file location. 
  TString fileName = Form("%s/SPD_Noisy_%d.root",fFileLocation.Data(),eq);
  return ReadNoisyFromFileName(fileName.Data());
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromFileName(const char *fileName) {
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
      }
    }
  }
  return kTRUE;
}
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
void AliITSOnlineCalibrationSPDhandler::WriteToFilesAlways() {
  // write the lists of dead and noisy to files
  for (UInt_t eq=0; eq<20; eq++) {
    WriteDeadToFile(eq);
    WriteNoisyToFile(eq);
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::WriteToFiles() {
  // write the lists of dead and noisy to files (only if there are >0 dead or noisy pixels) , returns nr of files produced
  return (WriteNoisyToFiles() + WriteDeadToFiles());
}
void AliITSOnlineCalibrationSPDhandler::WriteDeadToFilesAlways() {
  // write the lists of dead to files
  for (UInt_t eq=0; eq<20; eq++) {
      WriteDeadToFile(eq);
  }
}
void AliITSOnlineCalibrationSPDhandler::WriteNoisyToFilesAlways() {
  // write the lists of noisy to files
  for (UInt_t eq=0; eq<20; eq++) {
    WriteNoisyToFile(eq);
  }
}
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
void AliITSOnlineCalibrationSPDhandler::WriteDeadToFile(UInt_t eq) {
  // write the lists of dead and noisy for module to file
  AliITSOnlineCalibrationSPD* calib = new AliITSOnlineCalibrationSPD();
  calib->SetEqNr(eq);
  calib->SetBadList(GetDeadArrayOnline(eq));
  calib->SetNrBad(GetNrDeadEq(eq));
  TString fileName = Form("%s/SPD_Dead_%d.root",fFileLocation.Data(),eq);
  TFile file(fileName.Data(), "RECREATE");
  file.WriteTObject(calib, "AliITSOnlineCalibrationSPD");
  file.Close();
  delete calib;
}
void AliITSOnlineCalibrationSPDhandler::WriteNoisyToFile(UInt_t eq) {
  // write the lists of noisy and noisy for module to file
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
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadModuleFromDB(UInt_t module, Int_t runNr, Bool_t treeSerial) {
  // reads dead pixels from DB for given module and runNr
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) {
    man->SetDefaultStorage("local://$ALICE_ROOT");
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
  UInt_t nrDead = calibSPD->GetNrBad();
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
  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyModuleFromDB(UInt_t module, Int_t runNr, Bool_t treeSerial) {
  // reads noisy pixels from DB for given module and runNr
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) {
    man->SetDefaultStorage("local://$ALICE_ROOT");
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
  UInt_t nrNoisy = calibSPD->GetNrBad();
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
Bool_t AliITSOnlineCalibrationSPDhandler::ReadFromDB(Int_t runNr, Bool_t treeSerial) {
  // reads dead and noisy pixels from DB for given runNr
  // note that you may want to clear the lists (if they are not empty) before reading
  return (ReadNoisyFromDB(runNr,treeSerial) && ReadDeadFromDB(runNr,treeSerial));
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromDB(Int_t runNr, Bool_t treeSerial) {
  // reads dead pixels from DB for given runNr
  // note that you may want to clear the list (if it is not empty) before reading
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) {
    man->SetDefaultStorage("local://$ALICE_ROOT");
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
    //    printf("Reading module %d\n",module);
    calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    UInt_t nrDead = calibSPD->GetNrBad();
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
  }
  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromDB(Int_t runNr, Bool_t treeSerial) {
  // reads noisy pixels from DB for given runNr
  // note that you may want to clear the list (if it is not empty) before reading
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) {
    man->SetDefaultStorage("local://$ALICE_ROOT");
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
    //    printf("Reading module %d\n",module);
    calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    UInt_t nrNoisy = calibSPD->GetNrBad();
    if (nrNoisy>0) {    
      if (!treeSerial) {
	printf("*** mod %d nrnoisy=%d\n",module,nrNoisy);
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
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromCalibObj(TObjArray* calObj) {
  // reads dead pixels from calib object
  for (UInt_t module=0; module<240; module++) {
    for (Int_t i=0; i<((AliITSCalibrationSPD*)calObj->At(module))->GetNrBad(); i++) {
      SetDeadPixelM(module,
		    ((AliITSCalibrationSPD*)calObj->At(module))->GetBadColAt(i),
		    ((AliITSCalibrationSPD*)calObj->At(module))->GetBadRowAt(i));
    }
  }
  return kTRUE;
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromCalibObj(TObjArray* calObj) {
  // reads noisy pixels from calib object
  for (UInt_t module=0; module<240; module++) {
    for (Int_t i=0; i<((AliITSCalibrationSPD*)calObj->At(module))->GetNrBad(); i++) {
      SetNoisyPixelM(module,
		     ((AliITSCalibrationSPD*)calObj->At(module))->GetBadColAt(i),
		     ((AliITSCalibrationSPD*)calObj->At(module))->GetBadRowAt(i));
    }
  }
  return kTRUE;
}
Bool_t AliITSOnlineCalibrationSPDhandler::WriteToDB(Int_t runNrStart, Int_t runNrEnd) {
  // writes dead and noisy pixels to DB for given runNrs
  // overwrites any previous entries
  return (WriteNoisyToDB(runNrStart,runNrEnd) && WriteDeadToDB(runNrStart,runNrEnd));
}
Bool_t AliITSOnlineCalibrationSPDhandler::WriteDeadToDB(Int_t runNrStart, Int_t runNrEnd) {
  // writes dead pixels to DB for given runNrs
  // overwrites any previous entries
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) {
    man->SetDefaultStorage("local://$ALICE_ROOT");
  }
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Henrik Tydesjo");
  metaData->SetComment("Created by AliITSOnlineCalibrationSPDhandler.");
  AliCDBId idCalSPD("ITS/Calib/SPDDead",runNrStart,runNrEnd);
  TObjArray* spdEntry = new TObjArray(240);
  spdEntry->SetOwner(kTRUE);
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calObj = new AliITSCalibrationSPD();

    // *** this is temporarily hard coded here ********************
    // (later these parameters will be separated from the cal.obj.)
    calObj->SetThresholds(3000, 250);
    calObj->SetBiasVoltage(18.182);
    calObj->SetNoiseParam(0,0);
    calObj->SetCouplingParam(0.,0.055);
    // *** remove later...
    // ************************************************************

    spdEntry->Add(calObj);
  }
  for(UInt_t module=0; module<240; module++){
    ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNrBad( GetNrDead(module) );
    ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetBadList( GetDeadArray(module) );
  }
  AliCDBEntry* cdbEntry = new AliCDBEntry((TObject*)spdEntry,idCalSPD,metaData);
  man->Put(cdbEntry);
  delete spdEntry;
  delete cdbEntry;
  delete metaData;
  return kTRUE;
}
Bool_t AliITSOnlineCalibrationSPDhandler::WriteNoisyToDB(Int_t runNrStart, Int_t runNrEnd) {
  // writes noisy pixels to DB for given runNrs
  // overwrites any previous entries
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) {
    man->SetDefaultStorage("local://$ALICE_ROOT");
  }
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Henrik Tydesjo");
  metaData->SetComment("Created by AliITSOnlineCalibrationSPDhandler.");
  AliCDBId idCalSPD("ITS/Calib/SPDNoisy",runNrStart,runNrEnd);
  TObjArray* spdEntry = new TObjArray(240);
  spdEntry->SetOwner(kTRUE);
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calObj = new AliITSCalibrationSPD();

    // *** this is temporarily hard coded here ********************
    // (later these parameters will be separated from the cal.obj.)
    calObj->SetThresholds(3000, 250);
    calObj->SetBiasVoltage(18.182);
    calObj->SetNoiseParam(0,0);
    calObj->SetCouplingParam(0.,0.055);
    // *** remove later...
    // ************************************************************

    spdEntry->Add(calObj);
  }
  for(UInt_t module=0; module<240; module++){
    ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNrBad( GetNrNoisy(module) );
    ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetBadList( GetNoisyArray(module) );
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
void AliITSOnlineCalibrationSPDhandler::RecursiveInsertNoisy(AliITSCalibrationSPD* calibSPD, UInt_t module, Int_t lowInd, Int_t highInd) {
  // inserts noisy pixels recursively, used when reading from db
  printf("Rec mod %d : %d,%d\n",module,lowInd,highInd);
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
  dcsfile << "ActualDetCoonfiguration=" << "0,-1,-1\n"; // dummy values
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
TArrayI AliITSOnlineCalibrationSPDhandler::GetDeadArray(UInt_t module, Bool_t treeSerial) {
  // get a TArrayI of the dead pixels (format for the AliITSCalibrationSPD object)
  TArrayI returnArray;

  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t size=0;
  for (UInt_t chip=0; chip<5; chip++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,chip*32));
    size += fNrDead[gloChip];
  }
  returnArray.Set(size*2);
  UInt_t gloIndex=0;
  for (UInt_t chip=0; chip<5; chip++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,chip*32));
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
  return returnArray;
}
TArrayI AliITSOnlineCalibrationSPDhandler::GetNoisyArray(UInt_t module, Bool_t treeSerial) {
  // get a TArrayI of the noisy pixels (format for the AliITSCalibrationSPD object)
  TArrayI returnArray;

  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t size=0;
  for (UInt_t chip=0; chip<5; chip++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,chip*32));
    size += fNrNoisy[gloChip];
  }
  returnArray.Set(size*2);
  UInt_t gloIndex=0;
  for (UInt_t chip=0; chip<5; chip++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,chip*32));
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
TArrayI AliITSOnlineCalibrationSPDhandler::GetDeadArrayOnline(UInt_t eq) {
  // get a TArrayI of the dead pixels (format for the AliITSOnlineCalibrationSPD object)
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
TArrayI AliITSOnlineCalibrationSPDhandler::GetNoisyArrayOnline(UInt_t eq) {
  // get a TArrayI of the noisy pixels (format for the AliITSOnlineCalibrationSPD object)
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
    printf("Eq %*d: %*d dead , %*d noisy\n",2,eq,5,GetNrDeadEq(eq),5,GetNrNoisyEq(eq));
  }
}
void AliITSOnlineCalibrationSPDhandler::PrintDead() const {
  // print the dead pixels to screen
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
Bool_t AliITSOnlineCalibrationSPDhandler::SetDeadPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // set a dead pixel, returns false if pixel is already dead
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return SetDeadPixel(eq,hs,chip,col,row);
}
Bool_t AliITSOnlineCalibrationSPDhandler::SetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // set a noisy pixel, returns false if pixel is already noisy
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return SetNoisyPixel(eq,hs,chip,col,row);
}
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
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetDeadPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // unset a dead pixel, returns false if pixel is not dead
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return UnSetDeadPixel(eq,hs,chip,col,row);
}
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // unset a noisy pixel, returns false if pixel is not noisy
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return UnSetNoisyPixel(eq,hs,chip,col,row);
}
UInt_t AliITSOnlineCalibrationSPDhandler::SetDeadChip(UInt_t eq, UInt_t hs, UInt_t chip) {
  // set a full chip dead, returns nr of new dead pixels
  UInt_t nrNew = 0;
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      if (SetDeadPixel(eq,hs,chip,col,row)) {
	nrNew++;
      }
    }
  }
  return nrNew;
}
UInt_t AliITSOnlineCalibrationSPDhandler::SetNoisyChip(UInt_t eq, UInt_t hs, UInt_t chip) {
  // set a full chip noisy, returns nr of new noisy pixels
  UInt_t nrNew = 0;
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      if (SetNoisyPixel(eq,hs,chip,col,row)) {
	nrNew++;
      }
    }
  }
  return nrNew;
}
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetDeadChip(UInt_t eq, UInt_t hs, UInt_t chip) {
  // unset a full dead chip, returns false if it was not all dead before
  UInt_t nrUnset = 0;
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      if (UnSetDeadPixel(eq,hs,chip,col,row)) {
	nrUnset++;
      }
    }
  }
  if (nrUnset+GetNrNoisyC(eq,hs,chip)<8192) return kFALSE;
  else return kTRUE;
}
Bool_t AliITSOnlineCalibrationSPDhandler::UnSetNoisyChip(UInt_t eq, UInt_t hs, UInt_t chip) {
  // unset a full noisy chip, returns false if it was not all noisy before
  UInt_t nrUnset = 0;
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      if (UnSetNoisyPixel(eq,hs,chip,col,row)) {
	nrUnset++;
      }
    }
  }
  if (nrUnset<8192) return kFALSE;
  else return kTRUE;
}
//____________________________________________________________________________________________
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDead(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const {
  // is the pixel dead?
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::IsPixelDead", "eq,hs,chip nrs (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  UInt_t key = GetKey(eq,hs,chip,col,row);
  if ( fDeadPixelMap[gloChip]->Find(key) != NULL ) {
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisy(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const {
  // is the pixel noisy?
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::IsPixelNoisy", "eq,hs,chip nrs (%d,%d,%d) out of bounds.",eq,hs,chip);
    return kFALSE;
  }
  UInt_t key = GetKey(eq,hs,chip,col,row);
  if ( fNoisyPixelMap[gloChip]->Find(key) != NULL ) {
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDeadM(UInt_t module, UInt_t colM, UInt_t rowM) const {
  // is the pixel dead?
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return IsPixelDead(eq,hs,chip,col,row);
}
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisyM(UInt_t module, UInt_t colM, UInt_t rowM) const  {
  // is the pixel noisy?
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(module,colM);
  UInt_t row = GetRowFromOffline(module,rowM);
  return IsPixelNoisy(eq,hs,chip,col,row);
}
Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDeadKey(Int_t key) const {
  // is this pixel dead?
  UInt_t eq = GetEqIdFromKey(key);
  UInt_t hs = GetHSFromKey(key);
  UInt_t chip = GetChipFromKey(key);
  UInt_t col = GetColFromKey(key);
  UInt_t row = GetRowFromKey(key);
  return IsPixelDead(eq,hs,chip,col,row);
}
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
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDead() const {
  // returns the total nr of dead pixels
  UInt_t nrDead = 0;
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    nrDead+=fNrDead[gloChip];
  }
  return nrDead;
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisy() const {
  // returns the total nr of noisy pixels
  UInt_t nrNoisy = 0;
  for (UInt_t gloChip=0; gloChip<1200; gloChip++) {
    nrNoisy+=fNrNoisy[gloChip];
  }
  return nrNoisy;
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAt(UInt_t index) {
  // get eq for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadEqIdAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAt(UInt_t index) {
  // get eq for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyEqIdAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAt(UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadHSAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAt(UInt_t index) {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyHSAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAt(UInt_t index) {
  // get chip for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadChipAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAt(UInt_t index) {
  // get chip for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyChipAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAt(UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadColAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAt(UInt_t index) {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyColAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAt(UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotDead(index,gloChip,chipIndex);
  return GetDeadRowAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAt(UInt_t index) {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexTotNoisy(index,gloChip,chipIndex);
  return GetNoisyRowAtC2(gloChip,chipIndex);
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
  for (UInt_t chip=0; chip<5; chip++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,chip*32));
    nrDead+=fNrDead[gloChip];

  }
  return nrDead;
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisy(UInt_t module) const {
  // returns the number of noisy pixels for a certain module
  if (module>=240) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrNoisy", "module nr (%d) out of bounds.",module);
    return 0;
  }
  UInt_t nrNoisy = 0;
  UInt_t eq = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  for (UInt_t chip=0; chip<5; chip++) {
    UInt_t gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,chip*32));
    nrNoisy+=fNrNoisy[gloChip];

  }
  return nrNoisy;
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAt(UInt_t module, UInt_t index) {
  // get eq for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadEqIdAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAt(UInt_t module, UInt_t index) {
  // get eq for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyEqIdAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAt(UInt_t module, UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadHSAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAt(UInt_t module, UInt_t index) {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyHSAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAt(UInt_t module, UInt_t index) {
  // get chip for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadChipAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAt(UInt_t module, UInt_t index) {
  // get chip for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyChipAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAt(UInt_t module, UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadColAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAt(UInt_t module, UInt_t index) {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyColAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAt(UInt_t module, UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexDead(module,index,gloChip,chipIndex);
  return GetDeadRowAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAt(UInt_t module, UInt_t index) {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexNoisy(module,index,gloChip,chipIndex);
  return GetNoisyRowAtC2(gloChip,chipIndex);
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
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAtEq(UInt_t eq, UInt_t index) const {
  // get eq for the dead pixel at position index in list of dead
  if (eq<20 && index<GetNrDeadEq(eq)) {
    return eq;
  }
  else {
    return 20;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAtEq(UInt_t eq, UInt_t index) const {
  // get eq for the noisy pixel at position index in list of noisy
  if (eq<20 && index<GetNrNoisyEq(eq)) {
    return eq;
  }
  else {
    return 20;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAtEq(UInt_t eq, UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqDead(eq,index,gloChip,chipIndex);
  return GetDeadHSAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAtEq(UInt_t eq, UInt_t index) {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqNoisy(eq,index,gloChip,chipIndex);
  return GetNoisyHSAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAtEq(UInt_t eq, UInt_t index) {
  // get chip for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqDead(eq,index,gloChip,chipIndex);
  return GetDeadChipAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAtEq(UInt_t eq, UInt_t index) {
  // get chip for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqNoisy(eq,index,gloChip,chipIndex);
  return GetNoisyChipAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAtEq(UInt_t eq, UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqDead(eq,index,gloChip,chipIndex);
  return GetDeadColAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAtEq(UInt_t eq, UInt_t index) {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqNoisy(eq,index,gloChip,chipIndex);
  return GetNoisyColAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAtEq(UInt_t eq, UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqDead(eq,index,gloChip,chipIndex);
  return GetDeadRowAtC2(gloChip,chipIndex);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAtEq(UInt_t eq, UInt_t index) {
  // get hs for the noisy pixel at position index in list of noisy
  UInt_t gloChip;
  UInt_t chipIndex;
  GetChipAndIndexEqNoisy(eq,index,gloChip,chipIndex);
  return GetNoisyRowAtC2(gloChip,chipIndex);
}
//____________________________________________________________________________________________
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDeadC(UInt_t eq, UInt_t hs, UInt_t chip) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNrDeadC2(gloChip);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisyC(UInt_t eq, UInt_t hs, UInt_t chip) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNrNoisyC2(gloChip);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadEqIdAtC2(gloChip,index);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNoisyEqIdAtC2(gloChip,index);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadHSAtC2(gloChip,index);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNoisyHSAtC2(gloChip,index);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadChipAtC2(gloChip,index);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNoisyChipAtC2(gloChip,index);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadColAtC2(gloChip,index);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetNoisyColAtC2(gloChip,index);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const {
  UInt_t gloChip = GetGloChip(eq,hs,chip);
  return GetDeadRowAtC2(gloChip,index);
}
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
    UInt_t eq = GetEqIdFromKey(key);
    UInt_t hs = GetHSFromKey(key);
    UInt_t chip = GetChipFromKey(key);
    UInt_t col = GetColFromKey(key);
    UInt_t row = GetRowFromKey(key);    
    UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
    UInt_t colM = AliITSRawStreamSPD::GetOfflineColFromOnline(eq,hs,chip,col);
    UInt_t rowM = AliITSRawStreamSPD::GetOfflineRowFromOnline(eq,hs,chip,row);
    returnMess = Form("%*d,%*d,%*d,%*d,%*d  |  %*d,%*d,%*d\n",2,eq,1,hs,1,chip,2,col,3,row,3,module,3,colM,3,rowM);
    return returnMess.Data();
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadPixelAsTextC", "Index %d out of bounds.", index);
    return returnMess.Data();
  }
}
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
    UInt_t eq = GetEqIdFromKey(key);
    UInt_t hs = GetHSFromKey(key);
    UInt_t chip = GetChipFromKey(key);
    UInt_t col = GetColFromKey(key);
    UInt_t row = GetRowFromKey(key);    
    UInt_t module = AliITSRawStreamSPD::GetOfflineModuleFromOnline(eq,hs,chip);
    UInt_t colM = AliITSRawStreamSPD::GetOfflineColFromOnline(eq,hs,chip,col);
    UInt_t rowM = AliITSRawStreamSPD::GetOfflineRowFromOnline(eq,hs,chip,row);
    returnMess = Form("%*d,%*d,%*d,%*d,%*d  |  %*d,%*d,%*d\n",2,eq,1,hs,1,chip,2,col,3,row,3,module,3,colM,3,rowM);
    return returnMess.Data();
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyPixelAsTextC", "Index %d out of bounds.", index);
    return returnMess.Data();
  }
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
  // returns nr of dead/noisy in this' lists and not in other's lists
  return GetNrDeadDiff(other) + GetNrNoisyDiff(other);
}
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
AliITSOnlineCalibrationSPDhandler* AliITSOnlineCalibrationSPDhandler::GetDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns handler with dead/noisy in this' lists, except for those in other's lists
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
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexDead(UInt_t module, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) {
  // find gloChip and chipIndex from module and index
  if (index<GetNrDead(module)) {
    UInt_t eq = GetEqIdFromOffline(module);
    UInt_t hs = GetHSFromOffline(module);
    
    UInt_t glVal=0;
    for (UInt_t chip=0; chip<5; chip++) {
      gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,chip*32));
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
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexNoisy(UInt_t module, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) {
  // find gloChip and chipIndex from module and index
  if (index<GetNrNoisy(module)) {
    UInt_t eq = GetEqIdFromOffline(module);
    UInt_t hs = GetHSFromOffline(module);
    
    UInt_t glVal=0;
    for (UInt_t chip=0; chip<5; chip++) {
      gloChip = GetGloChip(eq,hs,GetChipFromOffline(module,chip*32));
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
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexEqDead(UInt_t eq, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) {
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
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexEqNoisy(UInt_t eq, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) {
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
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexTotDead(UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) {
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
void AliITSOnlineCalibrationSPDhandler::GetChipAndIndexTotNoisy(UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) {
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
UInt_t AliITSOnlineCalibrationSPDhandler::GetHSFromOffline(UInt_t module) const {
  // module to hs mapping
  if (module>=240) {
    Error("AliITSOnlineCalibrationSPDhandler::GetHSFromOffline", "module nr (%d) out of bounds.",module);
    return 6;
  }
  return AliITSRawStreamSPD::GetOnlineHSFromOffline(module);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetChipFromOffline(UInt_t module, UInt_t colM) const {
  // module,colM to chip mapping
  if (module>=240 || colM>=160) {
    Error("AliITSOnlineCalibrationSPDhandler::GetChipFromOffline", "module,colM nrs (%d,%d) out of bounds.",module,colM);
    return 10;
  }
  return AliITSRawStreamSPD::GetOnlineChipFromOffline(module,colM);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetColFromOffline(UInt_t module, UInt_t colM) const {
  // colM to col mapping
  if (colM>=160) {
    Error("AliITSOnlineCalibrationSPDhandler::GetColFromOffline", "colM nr (%d) out of bounds.",colM);
    return 160;
  }
  return AliITSRawStreamSPD::GetOnlineColFromOffline(module,colM);
}

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
UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisyC2(UInt_t gloChip) const {
  // returns nr of noisy pixels on this chip
  if (gloChip>=1200) {
    Error("AliITSOnlineCalibrationSPDhandler::GetNrNoisyC2", "global chip nr (%d) out of bounds.",gloChip);
    return 0;
  }
  return fNrNoisy[gloChip];
}
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



















