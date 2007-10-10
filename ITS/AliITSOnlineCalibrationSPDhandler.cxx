/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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


/* $Id$ */

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

AliITSOnlineCalibrationSPDhandler::AliITSOnlineCalibrationSPDhandler():
  fFileLocation("."),
  fModuleMapInited(kFALSE)
{
  // constructor
  for (UInt_t module=0; module<240; module++) {
    fNrDead[module]=0;
    fNrNoisy[module]=0;
    fDeadPixelMap[module] = new AliITSIntMap();
    fNoisyPixelMap[module] = new AliITSIntMap();
  }
}

AliITSOnlineCalibrationSPDhandler::AliITSOnlineCalibrationSPDhandler(const AliITSOnlineCalibrationSPDhandler& handle): 
  fFileLocation("."),
  fModuleMapInited(kFALSE)
{
  // copy constructor
  for (UInt_t module=0; module<240; module++) {
    fNrDead[module] = handle.fNrDead[module];
    fNrNoisy[module] = handle.fNrNoisy[module];
    fDeadPixelMap[module] = handle.fDeadPixelMap[module]->Clone();
    fNoisyPixelMap[module] = handle.fNoisyPixelMap[module]->Clone();
  }
  fFileLocation = handle.fFileLocation;
}

AliITSOnlineCalibrationSPDhandler::~AliITSOnlineCalibrationSPDhandler() {
  ClearMaps();
}

AliITSOnlineCalibrationSPDhandler& AliITSOnlineCalibrationSPDhandler::operator=(const AliITSOnlineCalibrationSPDhandler& handle) {
  // assignment operator
  if (this!=&handle) {
    this->ClearMaps();
    for (UInt_t module=0; module<240; module++) {
      fNrDead[module] = handle.fNrDead[module];
      fNrNoisy[module] = handle.fNrNoisy[module];
      fDeadPixelMap[module] = handle.fDeadPixelMap[module]->Clone();
      fNoisyPixelMap[module] = handle.fNoisyPixelMap[module]->Clone();
    }
    fFileLocation = handle.fFileLocation;
  }
  return *this;
}

void AliITSOnlineCalibrationSPDhandler::ClearMaps() {
  // clear the lists of dead and noisy
  ResetDead();
  ResetNoisy();
}

Bool_t AliITSOnlineCalibrationSPDhandler::ReadFromFiles() {
  // read dead and noisy files from file location. returns true if at least one file found
  Bool_t returnval=kFALSE;
  for (UInt_t module=0; module<240; module++) {
    if (ReadFromFile(module)) {
      returnval=kTRUE;
    }
  }
  return returnval;
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadFromFile(UInt_t module) {
  // read dead and noisy files for module from file location. 
  TString fileName = Form("%s/SPD_DeadNoisy_%d.root",fFileLocation.Data(),module);
  return ReadFromFileName(fileName.Data());
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadFromFileName(const char *fileName) {
  // read dead and noisy from file fileName
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
	UInt_t module = calib->GetModuleNr();
	Int_t nrDead=calib->GetNrDead();
	for (Int_t index=0; index<nrDead; index++) {
	  UInt_t colM=calib->GetDeadColAt(index);
	  UInt_t row=calib->GetDeadRowAt(index);
	  SetDeadPixelM(module,colM,row);
	}
	Int_t nrNoisy=calib->GetNrNoisy();
	for (Int_t index=0; index<nrNoisy; index++) {
	  UInt_t colM=calib->GetNoisyColAt(index);
	  UInt_t row=calib->GetNoisyRowAt(index);
	  SetNoisyPixelM(module,colM,row);
	}
      }
    }
  }
  return kTRUE;
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromFiles() {
  // read dead files from file location. returns true if at least one file found
  Bool_t returnval=kFALSE;
  for (UInt_t module=0; module<240; module++) {
    if (ReadDeadFromFile(module)) {
      returnval=kTRUE;
    }
  }
  return returnval;
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromFile(UInt_t module) {
  // read dead for module from file location. 
  TString fileName = Form("%s/SPD_DeadNoisy_%d.root",fFileLocation.Data(),module);
  return ReadDeadFromFileName(fileName.Data());
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
	UInt_t module = calib->GetModuleNr();
	Int_t nrDead=calib->GetNrDead();
	for (Int_t index=0; index<nrDead; index++) {
	  UInt_t colM=calib->GetDeadColAt(index);
	  UInt_t row=calib->GetDeadRowAt(index);
	  SetDeadPixelM(module,colM,row);
	}
      }
    }
  }
  return kTRUE;
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromFiles() {
  // read noisy files from file location. returns true if at least one file found
  Bool_t returnval=kFALSE;
  for (UInt_t module=0; module<240; module++) {
    if (ReadNoisyFromFile(module)) {
      returnval=kTRUE;
    }
  }
  return returnval;
}
Bool_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromFile(UInt_t module) {
  // read noisy for module from file location. 
  TString fileName = Form("%s/SPD_DeadNoisy_%d.root",fFileLocation.Data(),module);
  return ReadNoisyFromFileName(fileName.Data());
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
	UInt_t module = calib->GetModuleNr();
	Int_t nrNoisy=calib->GetNrNoisy();
	for (Int_t index=0; index<nrNoisy; index++) {
	  UInt_t colM=calib->GetNoisyColAt(index);
	  UInt_t row=calib->GetNoisyRowAt(index);
	  SetNoisyPixelM(module,colM,row);
	}
      }
    }
  }
  return kTRUE;
}


UInt_t AliITSOnlineCalibrationSPDhandler::ReadDeadFromText(const char *fileName) {
  // read dead from a text file (lines with eqId,hs,chip,col,row). returns nr of pixels added (not already there)
  UInt_t newNrDead=0;
  ifstream textFile;
  textFile.open(fileName, ifstream::in);
  if (textFile.fail()) {
    Error("AliITSOnlineCalibrationSPDhandler::ReadDeadFromText","No noisy text file (%s) present.",fileName);
  }
  else {
    while(1) {
      UInt_t eqId,hs,chip,col,row;
      textFile >> eqId; if (textFile.eof()) break;
      textFile >> hs; if (textFile.eof()) break;
      textFile >> chip; if (textFile.eof()) break;
      textFile >> col; if (textFile.eof()) break;
      textFile >> row; 
      if (SetDeadPixel(eqId,hs,chip,col,row)) {
	newNrDead++;
      }
      if (textFile.eof()) break;
    }
    textFile.close();
  }
  return newNrDead;
}

UInt_t AliITSOnlineCalibrationSPDhandler::ReadNoisyFromText(const char *fileName) {
  // read noisy from a text file (lines with eqId,hs,chip,col,row). returns nr of pixels added (not already here)
  UInt_t newNrNoisy=0;
  ifstream textFile;
  textFile.open(fileName, ifstream::in);
  if (textFile.fail()) {
    Error("AliITSOnlineCalibrationSPDhandler::ReadNoisyFromText","No noisy text file (%s) present.",fileName);
  }
  else {
    while(1) {
      UInt_t eqId,hs,chip,col,row;
      textFile >> eqId; if (textFile.eof()) break;
      textFile >> hs; if (textFile.eof()) break;
      textFile >> chip; if (textFile.eof()) break;
      textFile >> col; if (textFile.eof()) break;
      textFile >> row; 
      if (SetNoisyPixel(eqId,hs,chip,col,row)) {
	newNrNoisy++;
      }
      if (textFile.eof()) break;
    }
    textFile.close();
  }
  return newNrNoisy;
}


#ifndef SPD_DA_OFF
Bool_t AliITSOnlineCalibrationSPDhandler::ReadModuleFromDB(UInt_t module, Int_t runNr) {
  // reads dead and noisy pixels from DB for given module and runNr
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) {
    man->SetDefaultStorage("local://$ALICE_ROOT");
  }
  AliCDBEntry *cdbEntry = man->Get("ITS/Calib/CalibSPD", runNr);
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return kFALSE;
  }
  else {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadModuleFromDB","Calibration for run %d not found in database.",runNr);
    return kFALSE;
  }
  AliITSCalibrationSPD* calibSPD;
  calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
  UInt_t nrDead = calibSPD->GetNrDead();
  for (UInt_t index=0; index<nrDead; index++) {
    UInt_t colM = calibSPD->GetDeadColAt(index);
    UInt_t row  = calibSPD->GetDeadRowAt(index);
    SetDeadPixelM(module,colM,row);
  }
  UInt_t nrNoisy = calibSPD->GetNrNoisy();
  for (UInt_t index=0; index<nrNoisy; index++) {
    UInt_t colM = calibSPD->GetNoisyColAt(index);
    UInt_t row  = calibSPD->GetNoisyRowAt(index);
    SetNoisyPixelM(module,colM,row);
  }
  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::ReadFromDB(Int_t runNr) {
  // reads dead and noisy pixels from DB for given runNr
  // note that you may want to clear the lists (if they are not empty) before reading
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) {
    man->SetDefaultStorage("local://$ALICE_ROOT");
  }
  AliCDBEntry *cdbEntry = man->Get("ITS/Calib/CalibSPD", runNr);
  TObjArray* spdEntry;
  if(cdbEntry) {
    spdEntry = (TObjArray*)cdbEntry->GetObject();
    if(!spdEntry) return kFALSE;
  }
  else {
    Warning("AliITSOnlineCalibrationSPDhandler::ReadFromDB","Calibration for run %d not found in database.",runNr);
    return kFALSE;
  }
  AliITSCalibrationSPD* calibSPD;
  for (UInt_t module=0; module<240; module++) {
    //    printf("Reading module %d\n",module);
    calibSPD = (AliITSCalibrationSPD*) spdEntry->At(module);
    UInt_t nrDead = calibSPD->GetNrDead();
    for (UInt_t index=0; index<nrDead; index++) {
      UInt_t colM = calibSPD->GetDeadColAt(index);
      UInt_t row  = calibSPD->GetDeadRowAt(index);
      SetDeadPixelM(module,colM,row);
    }
    UInt_t nrNoisy = calibSPD->GetNrNoisy();
    for (UInt_t index=0; index<nrNoisy; index++) {
      UInt_t colM = calibSPD->GetNoisyColAt(index);
      UInt_t row  = calibSPD->GetNoisyRowAt(index);
      SetNoisyPixelM(module,colM,row);
    }
  }
  spdEntry->SetOwner(kTRUE);
  spdEntry->Clear();
  return kTRUE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::WriteToDB(Int_t runNrStart, Int_t runNrEnd) {
  // writes dead and noisy pixels to DB for given runNr
  // overwrites any previous entries
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) {
    man->SetDefaultStorage("local://$ALICE_ROOT");
  }
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("Henrik Tydesjo");
  metaData->SetComment("Created by AliITSOnlineCalibrationSPDhandler.");
  AliCDBId idCalSPD("ITS/Calib/CalibSPD",runNrStart,runNrEnd);
  TObjArray* spdEntry = new TObjArray(240);
  spdEntry->SetOwner(kTRUE);
  for(UInt_t module=0; module<240; module++){
    AliITSCalibrationSPD* calObj = new AliITSCalibrationSPD();

    // *** this is temporarily hard coded here ********************
    // (later these parameters will be separated from the cal.obj.)
    calObj->SetThresholds(3000, 250);
    calObj->SetBiasVoltage(18.182);
    calObj->SetNoiseParam(0,0);
    // CouplingRaw changed to 0.055 (fine tuning), was 0.047 in PDC06
    calObj->SetCouplingParam(0.,0.055);
    // *** remove later...
    // ************************************************************

    spdEntry->Add(calObj);
  }
  for(UInt_t module=0; module<240; module++){
    ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNrDead( GetNrDead(module) );
    ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetDeadList( GetDeadArray(module) );
    ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNrNoisy( GetNrNoisy(module) );
    ((AliITSCalibrationSPD*) spdEntry->At(module)) -> SetNoisyList( GetNoisyArray(module) );
  }
  AliCDBEntry* cdbEntry = new AliCDBEntry((TObject*)spdEntry,idCalSPD,metaData);
  man->Put(cdbEntry);
  delete spdEntry;
  delete cdbEntry;
  delete metaData;
  return kTRUE;
}
#endif

void AliITSOnlineCalibrationSPDhandler::WriteToFilesAlways() {
  // write the lists of dead and noisy to files (only if there are >0 dead or noisy pixels)
  for (UInt_t module=0; module<240; module++) {
    WriteToFile(module);
  }
}
void AliITSOnlineCalibrationSPDhandler::WriteToFiles() {
  // write the lists of dead and noisy to files (only if there are >0 dead or noisy pixels)
  for (UInt_t module=0; module<240; module++) {
    if (fNrDead[module]+fNrNoisy[module] > 0) {
      WriteToFile(module);
    }
  }
}
void AliITSOnlineCalibrationSPDhandler::WriteToFile(UInt_t module) {
  // write the lists of dead and noisy for module to file
  AliITSOnlineCalibrationSPD* calib = new AliITSOnlineCalibrationSPD();
  calib->SetModuleNr(module);
  calib->SetDeadList(GetDeadArray(module));
  calib->SetNoisyList(GetNoisyArray(module));
  calib->SetNrDead(GetNrDead(module));
  calib->SetNrNoisy(GetNrNoisy(module));
  TString fileName = Form("%s/SPD_DeadNoisy_%d.root",fFileLocation.Data(),module);
  TFile file(fileName.Data(), "RECREATE");
  file.WriteTObject(calib, "AliITSOnlineCalibrationSPD");
  file.Close();
  delete calib;
}
void AliITSOnlineCalibrationSPDhandler::WriteDeadToFiles() {
  // write the lists of dead to files (only if there are >0 dead pixels)
  for (UInt_t module=0; module<240; module++) {
    if (fNrDead[module] > 0) {
      WriteDeadToFile(module);
    }
  }
}
void AliITSOnlineCalibrationSPDhandler::WriteDeadToFile(UInt_t module) {
  // write the lists of dead for module to file
  AliITSOnlineCalibrationSPD* calib = new AliITSOnlineCalibrationSPD();
  calib->SetModuleNr(module);
  calib->SetDeadList(GetDeadArray(module));
  calib->SetNrDead(GetNrDead(module));
  TString fileName = Form("%s/SPD_DeadNoisy_%d.root",fFileLocation.Data(),module);
  TFile file(fileName.Data(), "RECREATE");
  file.WriteTObject(calib, "AliITSOnlineCalibrationSPD");
  file.Close();
  delete calib;
}
void AliITSOnlineCalibrationSPDhandler::WriteNoisyToFiles() {
  // write the lists of noisy to files (only if there are >0 dead pixels)
  for (UInt_t module=0; module<240; module++) {
    if (fNrNoisy[module] > 0) {
      printf("writing noisy to file for module %d\n",module);
      WriteNoisyToFile(module);
    }
  }
}
void AliITSOnlineCalibrationSPDhandler::WriteNoisyToFile(UInt_t module) {
  // write the lists of noisy for module to file
  AliITSOnlineCalibrationSPD* calib = new AliITSOnlineCalibrationSPD();
  calib->SetModuleNr(module);
  calib->SetNoisyList(GetNoisyArray(module));
  calib->SetNrNoisy(GetNrNoisy(module));
  TString fileName = Form("%s/SPD_DeadNoisy_%d.root",fFileLocation.Data(),module);
  TFile file(fileName.Data(), "RECREATE");
  file.WriteTObject(calib, "AliITSOnlineCalibrationSPD");
  file.Close();
  delete calib;
}

TArrayI AliITSOnlineCalibrationSPDhandler::GetDeadArray(UInt_t module) {
  // get a TArrayI of the dead pixels (format for the AliITSCalibrationSPD object)
  TArrayI returnArray;
  returnArray.Set(GetNrDead(module)*2);
  fDeadPixelMap[module]->PrepareSerialize(); // for tree ordered values
  for (UInt_t index=0; index<GetNrDead(module); index++) {
    Int_t key = fDeadPixelMap[module]->GetKeyIndex(index);
    Int_t colM = GetColMFromKey(key);
    Int_t row = GetRowFromKey(key);
    returnArray.AddAt(colM,index*2);
    returnArray.AddAt(row,index*2+1);
  }
  return returnArray;
}
TArrayI AliITSOnlineCalibrationSPDhandler::GetNoisyArray(UInt_t module) {
  // get a TArrayI of the noisy pixels (format for the AliITSCalibrationSPD object)
  TArrayI returnArray;
  returnArray.Set(GetNrNoisy(module)*2);
  fNoisyPixelMap[module]->PrepareSerialize(); // for tree ordered values
  for (UInt_t index=0; index<GetNrNoisy(module); index++) {
    Int_t key = fNoisyPixelMap[module]->GetKeyIndex(index);
    Int_t colM = GetColMFromKey(key);
    Int_t row = GetRowFromKey(key);
    returnArray.AddAt(colM,index*2);
    returnArray.AddAt(row,index*2+1);
  }
  return returnArray;
}

void AliITSOnlineCalibrationSPDhandler::ResetDead() {
  // reset the dead pixel map
  for (UInt_t module=0; module<240; module++) {
    fNrDead[module]=0;
    fDeadPixelMap[module]->Clear();
  }
}

void AliITSOnlineCalibrationSPDhandler::ResetDeadForChip(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // clear the noisy pixels for this chip
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      Int_t key = GetKey(eqId,hs,chip,col,row);
      if (fDeadPixelMap[module]->Remove(key)) {
	fNrDead[module]--;
      }
    }
  }
}

Bool_t AliITSOnlineCalibrationSPDhandler::SetDeadPixel(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // set a dead pixel, returns false if pixel is already dead
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  Int_t key = GetKey(eqId,hs,chip,col,row);
//!!!  // if noisy we dont want to add it...
//!!!  if (fNoisyPixelMap[module]->Find(key) != NULL) return kFALSE;
  if (fDeadPixelMap[module]->Insert(key,module)) {
    fNrDead[module]++;
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::SetDeadPixelM(UInt_t module, UInt_t colM, UInt_t row) {
  // set a dead pixel, returns false if pixel is already dead
  UInt_t eqId = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(colM);
  Int_t key = GetKey(eqId,hs,chip,col,row);
//!!!  // if noisy we dont want to add it...
//!!!  if (fNoisyPixelMap[module]->Find(key) != NULL) return kFALSE;
  if (fDeadPixelMap[module]->Insert(key,module)) {
    fNrDead[module]++;
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::UnSetDeadPixel(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // unset a dead pixel, returns false if pixel is not dead
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  Int_t key = GetKey(eqId,hs,chip,col,row);
  if (fDeadPixelMap[module]->Remove(key)) {
    fNrDead[module]--;
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::UnSetDeadPixelM(UInt_t module, UInt_t colM, UInt_t row) {
  // unset a dead pixel, returns false if pixel is not dead
  UInt_t eqId = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(colM);
  Int_t key = GetKey(eqId,hs,chip,col,row);
  if (fDeadPixelMap[module]->Remove(key)) {
    fNrDead[module]--;
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDeadKey(Int_t key) const {
  // is this pixel noisy?
  UInt_t eqId = GetEqIdFromKey(key);
  UInt_t hs = GetHSFromKey(key);
  UInt_t chip = GetChipFromKey(key);
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  return IsPixelDeadMKey(module,key);
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDeadMKey(UInt_t module, Int_t key) const {
  // is this pixel dead?
  if ( fDeadPixelMap[module]->Find(key) != NULL ) {
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDead(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const {
  // is the pixel dead?
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  Int_t key = GetKey(eqId,hs,chip,col,row);
  return IsPixelDeadMKey(module,key);
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDeadM(UInt_t module, UInt_t colM, UInt_t row) {
  // is the pixel dead?
  UInt_t eqId = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(colM);
  Int_t key = GetKey(eqId,hs,chip,col,row);
  return IsPixelDeadMKey(module,key);
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDead() const {
  // returns the total nr of dead pixels
  UInt_t nrDead = 0;
  for (UInt_t module=0; module<240; module++) {
    nrDead+=fNrDead[module];
  }
  return nrDead;
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDead(UInt_t module) const {
  // returns the number of dead pixels for a certain module
  if (module<240) {
    return fNrDead[module];
  }
  else {
    return 0;
  }
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAt(UInt_t module, UInt_t index) {
  // get eqId for the dead pixel at position index in list of dead
  if (index<GetNrDead(module)) {
    Int_t key = fDeadPixelMap[module]->GetKeyIndex(index);
    return GetEqIdFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadEqIdAt", "Index %d out of bounds.", index);
    return 0;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadHSAt(UInt_t module, UInt_t index) {
  // get hs for the dead pixel at position index in list of dead
  if (index<GetNrDead(module)) {
    Int_t key = fDeadPixelMap[module]->GetKeyIndex(index);
    return GetHSFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadHSAt", "Index %d out of bounds.", index);
    return 0;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadChipAt(UInt_t module, UInt_t index) {
  // get chip for the dead pixel at position index in list of dead
  if (index<GetNrDead(module)) {
    Int_t key = fDeadPixelMap[module]->GetKeyIndex(index);
    return GetChipFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadChipAt", "Index %d out of bounds.", index);
    return 0;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadColAt(UInt_t module, UInt_t index) {
  // get column for the dead pixel at position index in list of dead
  if (index<GetNrDead(module)) {
    Int_t key = fDeadPixelMap[module]->GetKeyIndex(index);
    return GetColFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadColAt", "Index %d out of bounds.", index);
    return 0;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAt(UInt_t module, UInt_t index) {
  // get row for the dead pixel at position index in list of dead
  if (index<GetNrDead(module)) {
    Int_t key = fDeadPixelMap[module]->GetKeyIndex(index);
    return GetRowFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetDeadRowAt", "Index %d out of bounds.", index);
    return 0;
  }
}

void AliITSOnlineCalibrationSPDhandler::ResetNoisy() {
  // clear the list of noisy pixels
  for (UInt_t module=0; module<240; module++) {
    fNrNoisy[module]=0;
    fNoisyPixelMap[module]->Clear();
  }
}
void AliITSOnlineCalibrationSPDhandler::ResetNoisyForChip(UInt_t eqId, UInt_t hs, UInt_t chip) {
  // clear the noisy pixels for this chip
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      Int_t key = GetKey(eqId,hs,chip,col,row);
      if (fNoisyPixelMap[module]->Remove(key)) {
	fNrNoisy[module]--;
      }
    }
  }
}


Bool_t AliITSOnlineCalibrationSPDhandler::SetNoisyPixel(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // set a noisy pixel, returns false if already there
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  Int_t key = GetKey(eqId,hs,chip,col,row);
//!!!  // if dead before - remove from the dead list 
//!!!  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
//!!!  if (fDeadPixelMap.Remove(key)) {
//!!!    fNrDead[module]--;
//!!!  }
  if (fNoisyPixelMap[module]->Insert(key,col)) {
    fNrNoisy[module]++;
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::SetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t row) {
  // set a noisy pixel, returns false if already there
  UInt_t eqId = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(colM);
  Int_t key = GetKey(eqId,hs,chip,col,row);
//!!!  // if dead before - remove from the dead list 
//!!!  if (fDeadPixelMap[module]->Remove(key)) {
//!!!    fNrDead[module]--;
//!!!  }
  if (fNoisyPixelMap[module]->Insert(key,col)) {
    fNrNoisy[module]++;
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::UnSetNoisyPixel(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // unset a noisy pixel, returns false if not there
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  Int_t key = GetKey(eqId,hs,chip,col,row);
  if (fNoisyPixelMap[module]->Remove(key)) {
    fNrNoisy[module]--;
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::UnSetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t row) {
  // unset a noisy pixel, returns false if not there
  UInt_t eqId = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(colM);
  Int_t key = GetKey(eqId,hs,chip,col,row);
  if (fNoisyPixelMap[module]->Remove(key)) {
    fNrNoisy[module]--;
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisyKey(Int_t key) const {
  // is this pixel noisy?
  UInt_t eqId = GetEqIdFromKey(key);
  UInt_t hs = GetHSFromKey(key);
  UInt_t chip = GetChipFromKey(key);
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  return IsPixelNoisyMKey(module,key);
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisyMKey(UInt_t module, Int_t key) const {
  // is this pixel noisy?
  if ( fNoisyPixelMap[module]->Find(key) != NULL ) {
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisy(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const {
  // is this pixel noisy?
  UInt_t module = AliITSRawStreamSPD::GetModuleNumber(eqId,hs,chip);
  Int_t key = GetKey(eqId,hs,chip,col,row);
  return IsPixelNoisyMKey(module,key);
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisyM(UInt_t module, UInt_t colM, UInt_t row) {
  // is this pixel noisy?
  UInt_t eqId = GetEqIdFromOffline(module);
  UInt_t hs = GetHSFromOffline(module);
  UInt_t chip = GetChipFromOffline(module,colM);
  UInt_t col = GetColFromOffline(colM);
  Int_t key = GetKey(eqId,hs,chip,col,row);
  return IsPixelNoisyMKey(module,key);
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisy() const {
  // returns the total nr of noisy pixels
  UInt_t nrNoisy = 0;
  for (UInt_t module=0; module<240; module++) {
    nrNoisy+=fNrNoisy[module];
  }
  return nrNoisy;
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisy(UInt_t module) const {
// returns the number of noisy pixels for a certain module
  if (module<240) {
    return fNrNoisy[module];
  }
  else {
    return 0;
  }
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAt(UInt_t module, UInt_t index) {
  // get chip for the noisy pixel at position index in list of noisy
  if (index<GetNrNoisy(module)) {
    Int_t key = fNoisyPixelMap[module]->GetKeyIndex(index);
    return GetEqIdFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyEqIdAt", "Index %d out of bounds.", index);
    return 0;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyHSAt(UInt_t module, UInt_t index) {
  // get chip for the noisy pixel at position index in list of noisy
  if (index<GetNrNoisy(module)) {
    Int_t key = fNoisyPixelMap[module]->GetKeyIndex(index);
    return GetHSFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyHSAt", "Index %d out of bounds.", index);
    return 0;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyChipAt(UInt_t module, UInt_t index) {
  // get chip for the noisy pixel at position index in list of noisy
  if (index<GetNrNoisy(module)) {
    Int_t key = fNoisyPixelMap[module]->GetKeyIndex(index);
    return GetChipFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyChipAt", "Index %d out of bounds.", index);
    return 0;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAt(UInt_t module, UInt_t index) {
  // get column for the noisy pixel at position index in list of noisy
  if (index<GetNrNoisy(module)) {
    Int_t key = fNoisyPixelMap[module]->GetKeyIndex(index);
    return GetColFromKey(key);
  }
  else {
    Warning("AliITSOnlineCalibrationSPDhandler::GetNoisyColAt", "Index %d out of bounds.", index);
    return 0;
  }
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAt(UInt_t module, UInt_t index) {
  // get row for the noisy pixel at position index in list of noisy
  if (index<GetNrNoisy(module)) {
    Int_t key = fNoisyPixelMap[module]->GetKeyIndex(index);
    return GetRowFromKey(key);
  }
  else {
    Error("AliITSOnlineCalibrationSPDhandler::GetNoisyRowAt", "Index %d out of bounds.", index);
    return 0;
  }
}

void AliITSOnlineCalibrationSPDhandler::PrintDead() const {
  // print the dead pixels to screen
  printf("-----------------------------------\n");
  printf("Dead Pixels: (eqId,hs,chip,col,row)\n");
  printf("-----------------------------------\n");
  for (UInt_t module=0; module<240; module++) {
    for (UInt_t index=0; index<GetNrDead(module); index++) {
      Int_t key = fDeadPixelMap[module]->GetKeyIndex(index);
      UInt_t eqId = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      printf("%d,%d,%d,%d,%d\n",eqId,hs,chip,col,row);
    }
  }
}

void AliITSOnlineCalibrationSPDhandler::PrintNoisy() const {
  // print the dead pixels to screen
  printf("-----------------------------------\n");
  printf("Noisy Pixels: (eqId,hs,chip,col,row)\n");
  printf("-----------------------------------\n");
  for (UInt_t module=0; module<240; module++) {
    for (UInt_t index=0; index<GetNrNoisy(module); index++) {
      Int_t key = fNoisyPixelMap[module]->GetKeyIndex(index);
      UInt_t eqId = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      printf("%d,%d,%d,%d,%d\n",eqId,hs,chip,col,row);
    }
  }
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns nr of dead in this' lists and not in other's lists
  UInt_t returnval=0;
  for (UInt_t module=0; module<240; module++) {
    for (UInt_t ind1=0; ind1<fNrDead[module]; ind1++) {
      Int_t key = fDeadPixelMap[module]->GetKeyIndex(ind1);
      UInt_t eqId = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelDead(eqId,hs,chip,col,row))) {
	returnval++;
      }
    }
  }
  return returnval;
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetNrNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns nr of noisy in this' lists and not in other's lists
  UInt_t returnval=0;
  for (UInt_t module=0; module<240; module++) {
    for (UInt_t ind1=0; ind1<fNrNoisy[module]; ind1++) {
      Int_t key = fNoisyPixelMap[module]->GetKeyIndex(ind1);
      UInt_t eqId = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelNoisy(eqId,hs,chip,col,row))) {
	returnval++;
      }
    }
  }
  return returnval;
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetNrDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns nr of dead/noisy in this' lists and not in other's lists
  return GetNrDeadDiff(other) + GetNrNoisyDiff(other);
}

AliITSOnlineCalibrationSPDhandler* AliITSOnlineCalibrationSPDhandler::GetDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns handler with dead/noisy in this' lists, except for those in other's lists
  AliITSOnlineCalibrationSPDhandler* newHandler = new AliITSOnlineCalibrationSPDhandler();
  for (UInt_t module=0; module<240; module++) {
    for (UInt_t ind1=0; ind1<fNrDead[module]; ind1++) {
      Int_t key = fDeadPixelMap[module]->GetKeyIndex(ind1);
      UInt_t eqId = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelDead(eqId,hs,chip,col,row))) {
	newHandler->SetDeadPixel(eqId,hs,chip,col,row);
      }
    }
    for (UInt_t ind2=0; ind2<fNrNoisy[module]; ind2++) {
      Int_t key = fNoisyPixelMap[module]->GetKeyIndex(ind2);
      UInt_t eqId = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelNoisy(eqId,hs,chip,col,row))) {
	newHandler->SetNoisyPixel(eqId,hs,chip,col,row);
      }
    }
  }
  return newHandler;
}

AliITSOnlineCalibrationSPDhandler* AliITSOnlineCalibrationSPDhandler::GetDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns handler with dead in this' lists, except for those in other's lists
  AliITSOnlineCalibrationSPDhandler* newHandler = new AliITSOnlineCalibrationSPDhandler();
  for (UInt_t module=0; module<240; module++) {
    for (UInt_t ind1=0; ind1<fNrDead[module]; ind1++) {
      Int_t key = fDeadPixelMap[module]->GetKeyIndex(ind1);
      UInt_t eqId = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelDead(eqId,hs,chip,col,row))) {
	newHandler->SetDeadPixel(eqId,hs,chip,col,row);
      }
    }
  }
  return newHandler;
}

AliITSOnlineCalibrationSPDhandler* AliITSOnlineCalibrationSPDhandler::GetNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const {
  // returns handler with noisy in this' lists, except for those in other's lists
  AliITSOnlineCalibrationSPDhandler* newHandler = new AliITSOnlineCalibrationSPDhandler();
  for (UInt_t module=0; module<240; module++) {
    for (UInt_t ind2=0; ind2<fNrNoisy[module]; ind2++) {
      Int_t key = fNoisyPixelMap[module]->GetKeyIndex(ind2);
      UInt_t eqId = GetEqIdFromKey(key);
      UInt_t hs = GetHSFromKey(key);
      UInt_t chip = GetChipFromKey(key);
      UInt_t col = GetColFromKey(key);
      UInt_t row = GetRowFromKey(key);
      if (!(other->IsPixelNoisy(eqId,hs,chip,col,row))) {
	newHandler->SetNoisyPixel(eqId,hs,chip,col,row);
      }
    }
  }
  return newHandler;
}

void AliITSOnlineCalibrationSPDhandler::InitModuleMaps() {
  // initializes the module mapping arrays needed for the methods below (GetEqIdFromOffline etc.)
  for (UInt_t iDDL=0; iDDL<20; iDDL++) {
    for (UInt_t iModule=0; iModule<12; iModule++) {
      UInt_t module = AliITSRawStreamSPD::GetModuleNumber(iDDL,iModule);
      fiDDL[module] = iDDL;
      fiModule[module] = iModule;
    }
  }
}

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
    for (UInt_t ind=0; ind<fNrNoisy[module]; ind++) {
      UInt_t newkey = GetNoisyEqIdAt(module,ind)*10*6 +
	GetNoisyHSAt(module,ind)*10 +
	GetNoisyChipAt(module,ind);
      if (newkey!=headkey) { // print eqId,hs,chip header
	headkey = newkey;
	dcsfile << "-" << newkey/(6*10) << "," << (newkey%(6*10))/10 << "," << (newkey%(6*10))%10 << "\n";
      }
      dcsfile << GetNoisyColAt(module,ind) << "," << GetNoisyRowAt(module,ind) << "\n";
    }
  }
  dcsfile.close();
}

UInt_t AliITSOnlineCalibrationSPDhandler::GetEqIdFromOffline(UInt_t module) {
  // module to eqId mapping
  if (!fModuleMapInited) InitModuleMaps();
  return fiDDL[module];
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetHSFromOffline(UInt_t module) {
  // module to hs mapping
  if (!fModuleMapInited) InitModuleMaps();
  return fiModule[module]/2;
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetChipFromOffline(UInt_t module, UInt_t colM) {
  // module,colM to chip mapping
  if (!fModuleMapInited) InitModuleMaps();
  return colM/32 + 5*(fiModule[module]%2);
}
UInt_t AliITSOnlineCalibrationSPDhandler::GetColFromOffline(UInt_t colM) const {
  // colM to col mapping
  return colM%32;
}
