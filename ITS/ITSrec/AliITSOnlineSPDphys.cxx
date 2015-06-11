////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online scan.   //
// Directly connected to a TFile with all containers.     //
// Handles reading and writing of this TFile. Hitmaps are //
// stored in this file (AliITSOnlineSPDHitArray).         //
// Also some general information is stored                //
// (AliITSOnlineSPDphysInfo).                             //
////////////////////////////////////////////////////////////

#include <math.h>
#include <TFile.h>
#include "AliITSOnlineSPDphys.h"
#include "AliITSOnlineSPDphysInfo.h"
#include "AliITSOnlineSPDHitArray.h"

AliITSOnlineSPDphys::AliITSOnlineSPDphys(const Char_t *fileName, Bool_t readFromGridFile) :
  fFile(NULL),
  fWrite(kFALSE),
  fModified(kFALSE),
  fInfoModified(kFALSE),
  fPhysInfo(NULL),
  fFileName(fileName)
{
  // constructor, open file for reading or writing
  // look for a previously saved info object 
  // (if file not found create a new one and return, else read)

  for(Int_t ihs =0; ihs<6; ihs++) fHitArray[ihs]=0x0;

  Bool_t bRead = readFromGridFile;

  if (!bRead) {
    FILE* fp0 = fopen(fFileName.Data(), "r");
    if (fp0 != NULL) {
      bRead=kTRUE;
      fclose(fp0);
    }
  }

  if (bRead) { // open file for reading
    fFile = TFile::Open(fFileName.Data(), "READ");
    if (fFile==NULL) { // grid file not found, create new local default file
      printf("ERROR: AliITSOnlineSPDphys: File %s not found! Creating 'test999.root' file instead\n",fFileName.Data());
      // create default empty file:
      fFileName = "test999.root";
      fPhysInfo = new AliITSOnlineSPDphysInfo();
      fInfoModified=kTRUE;
      fFile = new TFile(fFileName.Data(), "RECREATE");
      fWrite=kTRUE;
      InitHitmap();
    }
    else { // read from file (grid or local)
      fWrite=kFALSE;
      fFile->GetObject("AliITSOnlineSPDphysInfo", fPhysInfo);
      ReadHitmap();
    }
  }
  else { // create new local file
    fPhysInfo = new AliITSOnlineSPDphysInfo();
    fInfoModified=kTRUE;
    fFile = new TFile(fFileName.Data(), "RECREATE");
    fWrite=kTRUE;
    InitHitmap();
  }

}

AliITSOnlineSPDphys::AliITSOnlineSPDphys(const AliITSOnlineSPDphys& /*phys*/) :
  fFile(NULL),
  fWrite(kFALSE),
  fModified(kFALSE),
  fInfoModified(kFALSE),
  fPhysInfo(NULL),
  fFileName(".")
{
  for(Int_t i=0; i<6; i++) fHitArray[i]=0x0;
  printf("This object should not be copied!");
}

AliITSOnlineSPDphys::~AliITSOnlineSPDphys() {
  // destructor
  if (fModified) {
    SaveHitmap();
  }
  for (UInt_t hs=0; hs<6; hs++) {
    if (fHitArray[hs]!=NULL) {
      delete fHitArray[hs];
      fHitArray[hs]=NULL;
    }
  }
  if (fInfoModified) {
    if (!fWrite) {
      fFile->Close();
      delete fFile;
      fFile = new TFile(fFileName.Data(), "UPDATE");
      fWrite=kTRUE;
    }
    fFile->Delete("AliITSOnlineSPDphysInfo;*");
    fFile->WriteTObject(fPhysInfo, "AliITSOnlineSPDphysInfo");
  }
  if (fFile!=NULL) {
    delete fFile;
  }
}

AliITSOnlineSPDphys& AliITSOnlineSPDphys::operator=(const AliITSOnlineSPDphys& phys) {
  // assignment operator (should not be used)
  printf("This object should not be copied!");
  if (this!=&phys) {
    // still do nothing...
  }
  return *this;
}

void AliITSOnlineSPDphys::ClearThis() {
  // clear this phys, close file and open new
  for (UInt_t hs=0; hs<6; hs++) {
    if (fHitArray[hs]!=NULL) {
      delete fHitArray[hs];
    }
    fHitArray[hs] = NULL;
  }
  InitHitmap();
  fPhysInfo->ClearThis();
  fFile->Close();
  delete fFile;
  fFile = new TFile(fFileName.Data(), "RECREATE");
  fWrite=kTRUE;
  fFile->WriteTObject(fPhysInfo, "AliITSOnlineSPDphysInfo");
  fInfoModified=kTRUE;
}

void AliITSOnlineSPDphys::InitHitmap() {
  // init hit arrays and hit events
  for (UInt_t hs=0; hs<6; hs++) {
    fHitArray[hs] = new AliITSOnlineSPDHitArray();
  }
  fModified=kTRUE;
}

void AliITSOnlineSPDphys::AddPhys(AliITSOnlineSPDphys* phys2) {
  // add hitmap and info from another phys object
  if (phys2==NULL) return;
  if (GetEqNr()!=phys2->GetEqNr() && GetEqNr()!=999) {
    printf("AliITSOnlineSPDphys::AddPhys eqNr mismatch!\n");
    return;
  }
  if (GetEqNr()==999) {
    SetEqNr(phys2->GetEqNr());
  }
  UInt_t nrRuns = phys2->GetNrRuns();
  for (UInt_t i=0; i<nrRuns; i++) {
    AddRunNr(phys2->GetRunNr(i));
  }
  SetNrEvents(GetNrEvents() + phys2->GetNrEvents());
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {
      for (UInt_t col=0; col<32; col++) {
	for (UInt_t row=0; row<256; row++) {
	  SetHits(hs,chip,col,row,GetHits(hs,chip,col,row)+phys2->GetHits(hs,chip,col,row));
	}
      }
    }
  }
}

void AliITSOnlineSPDphys::ReadHitmap() {
  // read hitmap into memory
  for (UInt_t hs=0; hs<6; hs++) {
    TString hName = Form("HitArray_HS%d",hs);
    fFile->GetObject(hName.Data(), fHitArray[hs]);
  }
}

void AliITSOnlineSPDphys::SaveHitmap() {
  // save hitmap to file
  if (!fWrite) {
    fFile->Close();
    delete fFile;
    fFile = new TFile(fFileName.Data(), "UPDATE");
    fWrite=kTRUE;
  }
  for (UInt_t hs=0; hs<6; hs++) {
    TString hName = Form("HitArray_HS%d",hs);
    //TString hDelete = Form("%s;*",hName.Data());
    //fFile->Delete(hDelete.Data());
    fFile->WriteTObject(fHitArray[hs], hName.Data(),"overwrite");
  }
  fModified=kFALSE;
}

void AliITSOnlineSPDphys::SetHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi, UInt_t val) {
  // set nr of hits for pixel
  fHitArray[hs]->SetHits(chipi,coli,rowi,val);
  fModified=kTRUE;
}
void AliITSOnlineSPDphys::AddHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi, Int_t val) {
  // add val nr of hits for pixel (val could be negative)
  Int_t summedVal = fHitArray[hs]->GetHits(chipi,coli,rowi)+val;
  if (summedVal>0) {
    fHitArray[hs]->SetHits(chipi,coli,rowi,(UInt_t)summedVal);
  }
  else {
    fHitArray[hs]->SetHits(chipi,coli,rowi,0);
  }
  fModified=kTRUE;
}
void AliITSOnlineSPDphys::IncrementHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) {
  // increment nr of hits for pixel
  fHitArray[hs]->IncrementHits(chipi,coli,rowi);
  fModified=kTRUE;
}

UInt_t AliITSOnlineSPDphys::GetHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) {
  // get nr of hits for pixel
  return fHitArray[hs]->GetHits(chipi,coli,rowi);
}
Float_t AliITSOnlineSPDphys::GetHitsEfficiency(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) {
  // get the hit efficiency for pixel
  UInt_t ntr = GetNrEvents();
  if (ntr>0) {
    return ((Float_t)GetHits(hs,chipi,coli,rowi))/ntr;
  }
  else {
    return 0;
  }
}
Float_t AliITSOnlineSPDphys::GetHitsEfficiencyError(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) {
  // get error in hit efficiency for pixel
  Float_t hits = GetHits(hs,chipi,coli,rowi);
  UInt_t ntr = GetNrEvents();
  return sqrt(hits*(ntr-hits)/ntr)/ntr;
}
Float_t AliITSOnlineSPDphys::GetAverageMultiplicity(UInt_t hs, UInt_t chipi) {
  // get average multiplicity for a chip
  Float_t nrhits = 0;
  for (UInt_t chip=0;chip<10;chip++) {
    if (chipi==10 || chip==chipi) {
      for (Int_t col=0; col<32; col++) {
	for (Int_t row=0; row<256; row++) {
	  nrhits+=GetHits(hs,chip,col,row);
	}
      }
    }
  }
  UInt_t ntr = GetNrEvents();
  if (ntr>0) {
    return nrhits/ntr;
  }
  else {
    return 0;
  }
}
Float_t AliITSOnlineSPDphys::GetAverageMultiplicityTot(UInt_t hs) {
  // get average multiplicity for 10 chips
  return GetAverageMultiplicity(hs,10);
}

void AliITSOnlineSPDphys::AddRunNr(UInt_t val) {
  // add a run nr
  fPhysInfo->AddRunNr(val);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDphys::SetEqNr(UInt_t val) {
  // set router nr
  fPhysInfo->SetEqNr(val);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDphys::SetNrEvents(UInt_t val) {
  // set nr of events
  fPhysInfo->SetNrEvents(val);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDphys::AddNrEvents(Int_t val) {
  // add val nr of events (val could be negative)
  fPhysInfo->AddNrEvents(val);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDphys::IncrementNrEvents() {
  // increment nr of events
  fPhysInfo->IncrementNrEvents(); 
  fInfoModified=kTRUE;
}


UInt_t AliITSOnlineSPDphys::GetNrRuns() const {
  return fPhysInfo->GetNrRuns();
}
UInt_t AliITSOnlineSPDphys::GetRunNr(UInt_t posi) const {
  return fPhysInfo->GetRunNr(posi);
}
UInt_t AliITSOnlineSPDphys::GetEqNr() const {
  return fPhysInfo->GetEqNr();
}
UInt_t AliITSOnlineSPDphys::GetNrEvents() const {
  return fPhysInfo->GetNrEvents();
}

