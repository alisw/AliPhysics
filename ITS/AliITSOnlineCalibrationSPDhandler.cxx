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
#include <TArrayI.h>
#include <TFile.h>

AliITSOnlineCalibrationSPDhandler::AliITSOnlineCalibrationSPDhandler():
  fModuleNr(0),
  fDeadPixelMap(AliITSIntMap()),
  fNoisyPixelMap(AliITSIntMap())
{}

AliITSOnlineCalibrationSPDhandler::AliITSOnlineCalibrationSPDhandler(UInt_t module):
  fModuleNr(module),
  fDeadPixelMap(AliITSIntMap()),
  fNoisyPixelMap(AliITSIntMap())
{}

AliITSOnlineCalibrationSPDhandler::AliITSOnlineCalibrationSPDhandler(const AliITSOnlineCalibrationSPDhandler& handle): 
  fModuleNr(handle.fModuleNr),
  fDeadPixelMap(AliITSIntMap()),
  fNoisyPixelMap(AliITSIntMap())
{
  // copy constructor
  UInt_t nrDead = handle.GetNrDead();
  for (UInt_t index=0; index<nrDead; index++) {
    this->SetDeadPixel(handle.GetDeadColAt(index),handle.GetDeadRowAt(index));
  }
  UInt_t nrNoisy = handle.GetNrNoisy();
  for (UInt_t index=0; index<nrNoisy; index++) {
    this->SetNoisyPixel(handle.GetNoisyColAt(index),handle.GetNoisyRowAt(index));
  }
  sprintf(fFileLocation,"%s",handle.fFileLocation);
}

AliITSOnlineCalibrationSPDhandler::~AliITSOnlineCalibrationSPDhandler() {
  ClearMaps();
}

AliITSOnlineCalibrationSPDhandler& AliITSOnlineCalibrationSPDhandler::operator=(const AliITSOnlineCalibrationSPDhandler& handle) {
  // assignment operator
  if (this!=&handle) {
    this->ClearMaps();
    fModuleNr = handle.fModuleNr;
    UInt_t nrDead = handle.GetNrDead();
    for (UInt_t index=0; index<nrDead; index++) {
      this->SetDeadPixel(handle.GetDeadColAt(index),handle.GetDeadRowAt(index));
    }
    UInt_t nrNoisy = handle.GetNrNoisy();
    for (UInt_t index=0; index<nrNoisy; index++) {
      this->SetNoisyPixel(handle.GetNoisyColAt(index),handle.GetNoisyRowAt(index));
    }
    sprintf(fFileLocation,"%s",handle.fFileLocation);
  }
  return *this;
}

void AliITSOnlineCalibrationSPDhandler::ClearMaps() {
  // clear the lists of dead and noisy
  ResetDead();
  ResetNoisy();
}

void AliITSOnlineCalibrationSPDhandler::WriteToFile() {
  // write the lists of dead and noisy to default file
  Char_t fileName[200];
  sprintf(fileName,"%s/SPD_DeadNoisy_%d.root",fFileLocation,fModuleNr);
  WriteToFile(fileName);
}

void AliITSOnlineCalibrationSPDhandler::WriteToFile(Char_t* fileName) {
  // write the lists of dead and noisy to file fileName
  AliITSOnlineCalibrationSPD* calib = new AliITSOnlineCalibrationSPD();
  calib->SetModuleNr(GetModuleNr());
  calib->SetDeadList(GetDeadArray());
  calib->SetNoisyList(GetNoisyArray());
  calib->SetNrDead(GetNrDead());
  calib->SetNrNoisy(GetNrNoisy());

  TFile file(fileName, "RECREATE");
  file.WriteTObject(calib, "AliITSOnlineCalibrationSPD");
  file.Close();

  delete calib;
}

void AliITSOnlineCalibrationSPDhandler::ReadFromFile() {
  // read dead and noisy from default file
  Char_t fileName[200];
  sprintf(fileName,"%s/SPD_DeadNoisy_%d.root",fFileLocation,fModuleNr);
  ReadFromFile(fileName);
}
void AliITSOnlineCalibrationSPDhandler::ReadDeadFromFile() {
  // read dead from default file
  Char_t fileName[200];
  sprintf(fileName,"%s/SPD_DeadNoisy_%d.root",fFileLocation,fModuleNr);
  ReadDeadFromFile(fileName);
}
void AliITSOnlineCalibrationSPDhandler::ReadNoisyFromFile() {
  // read noisy from default file
  Char_t fileName[200];
  sprintf(fileName,"%s/SPD_DeadNoisy_%d.root",fFileLocation,fModuleNr);
  ReadNoisyFromFile(fileName);
}

void AliITSOnlineCalibrationSPDhandler::ReadFromFile(Char_t* fileName) {
  // read dead and noisy from file fileName (clear the previous list of dead and noisy pixels)
  ClearMaps();
  AliITSOnlineCalibrationSPD* calib;
  FILE* fp0 = fopen(fileName, "r");
  if (fp0 == NULL) {}
  else {
    fclose(fp0);
    TFile file(fileName, "READ");
    if (file.IsOpen()) {
      file.GetObject("AliITSOnlineCalibrationSPD", calib);
      file.Close();
      if (calib!=NULL) {
	SetModuleNr(calib->GetModuleNr());
	Int_t nrDead=calib->GetNrDead();
	for (Int_t index=0; index<nrDead; index++) {
	  Int_t col=calib->GetDeadColAt(index);
	  Int_t row=calib->GetDeadRowAt(index);
	  SetDeadPixel(col,row);
	}
	Int_t nrNoisy=calib->GetNrNoisy();
	for (Int_t index=0; index<nrNoisy; index++) {
	  Int_t col=calib->GetNoisyColAt(index);
	  Int_t row=calib->GetNoisyRowAt(index);
	  SetNoisyPixel(col,row);
	}
      }
    }
  }
}
void AliITSOnlineCalibrationSPDhandler::ReadDeadFromFile(Char_t* fileName) {
  // read dead from file fileName (clear the previous list of dead pixels)
  ResetDead();
  AliITSOnlineCalibrationSPD* calib;
  FILE* fp0 = fopen(fileName, "r");
  if (fp0 == NULL) {}
  else {
    fclose(fp0);
    TFile file(fileName, "READ");
    if (file.IsOpen()) {
      file.GetObject("AliITSOnlineCalibrationSPD", calib);
      file.Close();
      if (calib!=NULL) {
	SetModuleNr(calib->GetModuleNr());
	Int_t nrDead=calib->GetNrDead();
	for (Int_t index=0; index<nrDead; index++) {
	  Int_t col=calib->GetDeadColAt(index);
	  Int_t row=calib->GetDeadRowAt(index);
	  SetDeadPixel(col,row);
	}
      }
    }
  }
}
void AliITSOnlineCalibrationSPDhandler::ReadNoisyFromFile(Char_t* fileName) {
  // read noisy from file fileName (clear the previous list of noisy pixels)
  ResetNoisy();
  AliITSOnlineCalibrationSPD* calib;
  FILE* fp0 = fopen(fileName, "r");
  if (fp0 == NULL) {}
  else {
    fclose(fp0);
    TFile file(fileName, "READ");
    if (file.IsOpen()) {
      file.GetObject("AliITSOnlineCalibrationSPD", calib);
      file.Close();
      if (calib!=NULL) {
	SetModuleNr(calib->GetModuleNr());
	Int_t nrNoisy=calib->GetNrNoisy();
	for (Int_t index=0; index<nrNoisy; index++) {
	  Int_t col=calib->GetNoisyColAt(index);
	  Int_t row=calib->GetNoisyRowAt(index);
	  SetNoisyPixel(col,row);
	}
      }
    }
  }
}

TArrayI AliITSOnlineCalibrationSPDhandler::GetDeadArray() {
  // get a TArrayI of the dead pixels (format for the AliITSCalibrationSPD object)
  TArrayI returnArray;
  returnArray.Set(GetNrDead()*2);
  for (UInt_t index=0; index<fDeadPixelMap.GetNrEntries(); index++) {
    Int_t key = fDeadPixelMap.GetKey(index);
    Int_t col = GetColFromKey(key);
    Int_t row = GetRowFromKey(key);
    returnArray.AddAt(col,index*2);
    returnArray.AddAt(row,index*2+1);
  }
  return returnArray;
}

TArrayI AliITSOnlineCalibrationSPDhandler::GetNoisyArray() {
  // get a TArrayI of the noisy pixels (format for the AliITSCalibrationSPD object)
  TArrayI returnArray;
  returnArray.Set(GetNrNoisy()*2);
  for (UInt_t index=0; index<fNoisyPixelMap.GetNrEntries(); index++) {
    Int_t key = fNoisyPixelMap.GetKey(index);
    Int_t col = GetColFromKey(key);
    Int_t row = GetRowFromKey(key);
    returnArray.AddAt(col,index*2);
    returnArray.AddAt(row,index*2+1);
  }
  return returnArray;
}

void AliITSOnlineCalibrationSPDhandler::ResetDead() {
  fDeadPixelMap.Clear();
}


Bool_t AliITSOnlineCalibrationSPDhandler::SetDeadPixel(Int_t col, Int_t row) {
  // set a dead pixel, returns false if pixel is already dead or noisy
  Int_t key = GetKey(col,row);
  // if noisy we dont want to add it...
  if (fNoisyPixelMap.Find(key) != NULL) return kFALSE;
  return fDeadPixelMap.Insert(key,col);
}

Int_t AliITSOnlineCalibrationSPDhandler::GetDeadColAt(UInt_t index) const {
  // get column for the dead pixel at position index in list of dead
  if (index<fDeadPixelMap.GetNrEntries()) {
    Int_t key = fDeadPixelMap.GetKey(index);
    return GetColFromKey(key);
  }
  else return 0;
}

Int_t AliITSOnlineCalibrationSPDhandler::GetDeadRowAt(UInt_t index) const {
  // get row for the dead pixel at position index in list of dead
  if (index<fDeadPixelMap.GetNrEntries()) {
    Int_t key = fDeadPixelMap.GetKey(index);
    return GetRowFromKey(key);
  }
  else return 0;
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelDead(Int_t col, Int_t row) const {
  // is the pixel dead?
  Int_t key = GetKey(col,row);
  if ( fDeadPixelMap.Find(key) != NULL) {
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}

void AliITSOnlineCalibrationSPDhandler::ResetNoisy() {
  // clear the list of noisy pixels
  fNoisyPixelMap.Clear();
}

Bool_t AliITSOnlineCalibrationSPDhandler::SetNoisyPixel(Int_t col, Int_t row) {
  // set a noisy pixel, returns false if already there
  Int_t key = GetKey(col,row);
  // if dead before - remove from the dead list 
  fDeadPixelMap.Remove(key);
  return fNoisyPixelMap.Insert(key,col);
}

Int_t AliITSOnlineCalibrationSPDhandler::GetNoisyColAt(UInt_t index) const {
  // get column for the noisy pixel at position index in list of noisy
  if (index<fNoisyPixelMap.GetNrEntries()) {
    Int_t key = fNoisyPixelMap.GetKey(index);
    return GetColFromKey(key);
  }
  else return 0;
}

Int_t AliITSOnlineCalibrationSPDhandler::GetNoisyRowAt(UInt_t index) const {
  // get row for the noisy pixel at position index in list of noisy
  if (index<fNoisyPixelMap.GetNrEntries()) {
    Int_t key = fNoisyPixelMap.GetKey(index);
    return GetRowFromKey(key);
  }
  else return 0;
}

Bool_t AliITSOnlineCalibrationSPDhandler::IsPixelNoisy(Int_t col, Int_t row) const {
  // is this pixel noisy?
  Int_t key = GetKey(col,row);
  if ( fNoisyPixelMap.Find(key) != NULL ) {
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}

void AliITSOnlineCalibrationSPDhandler::PrintDead() const {
  // print the dead pixels to screen
  printf("-----------------------\n");
  printf("Dead Pixels Module %d:\n",fModuleNr);
  printf("-----------------------\n");
  for (UInt_t index=0; index<fDeadPixelMap.GetNrEntries(); index++) {
    Int_t key = fDeadPixelMap.GetKey(index);
    Int_t col = GetColFromKey(key);
    Int_t row = GetRowFromKey(key);
    printf("%d,%d\n",col,row);
  }
}

void AliITSOnlineCalibrationSPDhandler::PrintNoisy() const {
  // print the noisy pixels to screen
  printf("-----------------------\n");
  printf("Noisy Pixels Module %d:\n",fModuleNr);
  printf("-----------------------\n");
  for (UInt_t index=0; index<fNoisyPixelMap.GetNrEntries(); index++) {
    Int_t key = fNoisyPixelMap.GetKey(index);
    Int_t col = GetColFromKey(key);
    Int_t row = GetRowFromKey(key);
    printf("%d,%d\n",col,row);
  }
}
