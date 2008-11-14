////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online scan.   //
// Directly connected to a TFile with all containers.     //
// Handles reading and writing of this TFile.             //
// Hitmaps and information on nr of events with hits      //
// is stored in this file (AliITSOnlineSPDHitArray and    //
// AliITSOnlineSPDHitEvent). Also some general            //
// information is stored (AliITSOnlineSPDscanInfo).       //
// When switching between different steps of the scan,    //
// the previous step is automatically stored on the file. //
// With this scheme there is no risk of running out of    //
// memory.                                                //
////////////////////////////////////////////////////////////

#include <math.h>

#include <TFile.h>
#include "AliITSOnlineSPDscan.h"
#include "AliITSOnlineSPDscanInfo.h"
#include "AliITSOnlineSPDHitArray.h"
#include "AliITSOnlineSPDHitEvent.h"

AliITSOnlineSPDscan::AliITSOnlineSPDscan(const Char_t *fileName, Bool_t readFromGridFile) :
  fFile(NULL),
  fWrite(kFALSE),
  fCurrentStep(-1),
  fModified(kFALSE),
  fInfoModified(kFALSE),
  fScanInfo(NULL),
  fFileName(fileName)
{
  // constructor, open file for reading or writing
  // look for a previously saved info object 
  // (if file not found create a new one and return, else read)

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
      printf("ERROR: AliITSOnlineSPDscan: File %s not found! Creating 'test999.root' file instead\n",fFileName.Data());
      // create default empty file:
      fFileName = "test999.root";
      fScanInfo = new AliITSOnlineSPDscanInfo();
      fInfoModified=kTRUE;
      fFile = new TFile(fFileName.Data(), "RECREATE");
      fWrite=kTRUE;
    }
    else { // read from file (grid or local)
      fWrite=kFALSE;
      fFile->GetObject("AliITSOnlineSPDscanInfo", fScanInfo);
    }
  }
  else { // create new local file
    fScanInfo = new AliITSOnlineSPDscanInfo();
    fInfoModified=kTRUE;
    fFile = new TFile(fFileName.Data(), "RECREATE");
    fWrite=kTRUE;
  }

  Init();
}

AliITSOnlineSPDscan::AliITSOnlineSPDscan(const AliITSOnlineSPDscan& /*scan*/) :
  fFile(NULL),
  fWrite(kFALSE),
  fCurrentStep(-1),
  fModified(kFALSE),
  fInfoModified(kFALSE),
  fScanInfo(NULL),
  fFileName(".")
{
  printf("This object should not be copied!");
}

AliITSOnlineSPDscan::~AliITSOnlineSPDscan() {
  // destructor
  if (fModified) {
    SaveCurrentStep();
  }
  for (UInt_t hs=0; hs<6; hs++) {
    if (fCurrentHitArray[hs]!=NULL) {
      delete fCurrentHitArray[hs];
      fCurrentHitArray[hs]=NULL;
    }
    if (fCurrentHitEvent[hs]!=NULL) {
      delete fCurrentHitEvent[hs];
      fCurrentHitEvent[hs]=NULL;
    }
  }
  if (fInfoModified) {
    if (!fWrite) {
      fFile->Close();
      delete fFile;
      fFile = new TFile(fFileName.Data(), "UPDATE");
      fWrite=kTRUE;
    }
    fFile->Delete("AliITSOnlineSPDscanInfo;*");
    fFile->WriteTObject(fScanInfo, "AliITSOnlineSPDscanInfo");
  }
  if (fFile!=NULL) {
    delete fFile;
  }
}

AliITSOnlineSPDscan& AliITSOnlineSPDscan::operator=(const AliITSOnlineSPDscan& scan) {
  // assignment operator (should not be used)
  printf("This object should not be copied!");
  if (this!=&scan) {
    // still do nothing...
  }
  return *this;
}

void AliITSOnlineSPDscan::ClearThis() {
  // clear this scan, close file and open new
  for (UInt_t hs=0; hs<6; hs++) {
    if (fCurrentHitArray[hs]!=NULL) {
      delete fCurrentHitArray[hs];
    }
    fCurrentHitArray[hs] = NULL;
    if (fCurrentHitEvent[hs]!=NULL) {
      delete fCurrentHitEvent[hs];
    }
    fCurrentHitEvent[hs] = NULL;
  }
  fScanInfo->ClearThis();
  fFile->Close();
  delete fFile;
  fFile = new TFile(fFileName.Data(), "RECREATE");
  fWrite=kTRUE;
  fFile->WriteTObject(fScanInfo, "AliITSOnlineSPDscanInfo");
  fInfoModified=kTRUE;
}

void AliITSOnlineSPDscan::Init() {
  // init hit arrays and hit events
  for (UInt_t hs=0; hs<6; hs++) {
    fCurrentHitArray[hs]=NULL;
    fCurrentHitEvent[hs]=NULL;
  }

}

UInt_t AliITSOnlineSPDscan::AddScanStep() {
  // add a new scan step
  CreateNewStep();
  return fScanInfo->AddScanStep();
}

void AliITSOnlineSPDscan::CreateNewStep() {
  // create a new step
  // save current step to file (if modified)
  if (fModified) {
    SaveCurrentStep();
  }
  // create new step
  for (UInt_t hs=0; hs<6; hs++) {
    if (fCurrentHitArray[hs]!=NULL) {
      delete fCurrentHitArray[hs];
    }
    fCurrentHitArray[hs] = new AliITSOnlineSPDHitArray();
    if (fCurrentHitEvent[hs]!=NULL) {
      delete fCurrentHitEvent[hs];
    }
    fCurrentHitEvent[hs] = new AliITSOnlineSPDHitEvent();
  }
  fCurrentStep = fScanInfo->GetNSteps();
  fModified=kTRUE;
  fInfoModified=kTRUE;
}

void AliITSOnlineSPDscan::SwitchToStep(UInt_t nsi) {
  // switch to step nsi (save current step first if needed)
  if ((Int_t)nsi!=fCurrentStep) {
    if (fModified) {
      SaveCurrentStep();
    }
    for (UInt_t hs=0; hs<6; hs++) {
      if (fCurrentHitArray[hs]!=NULL) {
	delete fCurrentHitArray[hs];
	fCurrentHitArray[hs]=NULL;
      }
      if (fCurrentHitEvent[hs]!=NULL) {
	delete fCurrentHitEvent[hs];
	fCurrentHitEvent[hs]=NULL;
      }
    }
    if (nsi>=GetNSteps()) {
      FillGap(nsi); // makes fCurrentStep = nsi
    }
    else {
      fCurrentStep=nsi;
      ReadCurrentStep();
    }
  }
}

void AliITSOnlineSPDscan::FillGap(UInt_t nsi) {
  //create new steps until nsi is reached
  while (nsi>=GetNSteps()) {
    fCurrentStep = AddScanStep();
  }
}

void AliITSOnlineSPDscan::ReadCurrentStep() {
  // read current step index into memory
  for (UInt_t hs=0; hs<6; hs++) {
    TString stepName = Form("HitArray_HS%d_Step%d",hs,fCurrentStep);
    fFile->GetObject(stepName.Data(), fCurrentHitArray[hs]);
    TString stepName2 = Form("HitEvent_HS%d_Step%d",hs,fCurrentStep);
    fFile->GetObject(stepName2, fCurrentHitEvent[hs]);
  }
}

void AliITSOnlineSPDscan::SaveCurrentStep() {
  // save current step to file
  if (!fWrite) {
    fFile->Close();
    delete fFile;
    fFile = new TFile(fFileName.Data(), "UPDATE");
    fWrite=kTRUE;
  }
  for (UInt_t hs=0; hs<6; hs++) {
    TString stepName = Form("HitArray_HS%d_Step%d",hs,fCurrentStep);
    TString stepDelete = Form("%s;*",stepName.Data());
    fFile->Delete(stepDelete.Data());
    fFile->WriteTObject(fCurrentHitArray[hs], stepName.Data());
    TString stepName2 = Form("HitEvent_HS%d_Step%d",hs,fCurrentStep);
    TString stepDelete2 = Form("%s;*",stepName2.Data());
    fFile->Delete(stepDelete2.Data());
    fFile->WriteTObject(fCurrentHitEvent[hs], stepName2.Data());
  }
  fModified=kFALSE;
}

void AliITSOnlineSPDscan::SetHits(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi, UInt_t val) {
  // set nr of hits for pixel
  SwitchToStep(nsi);
  fCurrentHitArray[hs]->SetHits(chipi,coli,rowi,val);
  fModified=kTRUE;
}
void AliITSOnlineSPDscan::IncrementHits(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) {
  // increment nr of hits for pixel
  SwitchToStep(nsi);
  fCurrentHitArray[hs]->IncrementHits(chipi,coli,rowi);
  fModified=kTRUE;
}
void AliITSOnlineSPDscan::SetHitEvents(UInt_t nsi, UInt_t hs, UInt_t chipi, Int_t val) {
  // set nr of hit events for a chip
  SwitchToStep(nsi);
  fCurrentHitEvent[hs]->SetHitEvent(chipi,val);
  fModified=kTRUE;
}
void AliITSOnlineSPDscan::SetHitEventsTot(UInt_t nsi, UInt_t hs, Int_t val) {
  // set nr of hit events for 10 chips together
  SetHitEvents(nsi,hs,10,val);
}
void AliITSOnlineSPDscan::IncrementHitEvents(UInt_t nsi, UInt_t hs, UInt_t chipi) {
  // increment nr of hit events for a chip
  SwitchToStep(nsi);
  fCurrentHitEvent[hs]->IncrementHitEvent(chipi);
  fModified=kTRUE;
}
void AliITSOnlineSPDscan::IncrementHitEventsTot(UInt_t nsi, UInt_t hs) {
  // increment nr of hit events for 10 chips
  IncrementHitEvents(nsi,hs,10);
}


UInt_t AliITSOnlineSPDscan::GetHits(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) {
  // get nr of hits for pixel
  if (nsi<GetNSteps()) {
    SwitchToStep(nsi);
    return fCurrentHitArray[hs]->GetHits(chipi,coli,rowi);
  }
  else {
    return 0;
  }
}
Float_t AliITSOnlineSPDscan::GetHitsEfficiency(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) {
  // get the hit efficiency for pixel
  UInt_t ntr = GetTriggers(nsi);
  if (ntr>0) {
    return ((Float_t)GetHits(nsi,hs,chipi,coli,rowi))/ntr;
  }
  else {
    return 0;
  }
}
Float_t AliITSOnlineSPDscan::GetHitsEfficiencyError(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) {
  // get error in hit efficiency for pixel
  Float_t hits = GetHits(nsi,hs,chipi,coli,rowi);
  UInt_t ntr = GetTriggers(nsi);
  return sqrt(hits*(ntr-hits)/ntr)/ntr;
}
UInt_t AliITSOnlineSPDscan::GetHitEvents(UInt_t nsi, UInt_t hs, UInt_t chipi) {
  // get nr of hit events for a chip
  if (nsi<GetNSteps()) {
    SwitchToStep(nsi);
    return fCurrentHitEvent[hs]->GetHitEvent(chipi);
  }
  else {
    return 0;
  }
}
UInt_t AliITSOnlineSPDscan::GetHitEventsTot(UInt_t nsi, UInt_t hs) {
  // get nr of hit events for 10 chips
  return GetHitEvents(nsi,hs,10);
}
Float_t AliITSOnlineSPDscan::GetHitEventsEfficiency(UInt_t nsi, UInt_t hs, UInt_t chipi) {
  // get the hit events efficiency for a chip
  UInt_t ntr = GetTriggers(nsi);
  if (ntr>0) {
    return ((Float_t)GetHitEvents(nsi,hs,chipi))/ntr;
  }
  else {
    return 0;
  }
}
Float_t AliITSOnlineSPDscan::GetHitEventsTotEfficiency(UInt_t nsi, UInt_t hs) {
  // get the hit events efficiency for 10 chips
  return GetHitEventsEfficiency(nsi,hs,10);
}
Float_t AliITSOnlineSPDscan::GetHitEventsEfficiencyError(UInt_t nsi, UInt_t hs, UInt_t chipi) {
  // get error in hit events efficiency for a chip
  Float_t hitevents = (Float_t) GetHitEvents(nsi,hs,chipi);
  UInt_t ntr = GetTriggers(nsi);
  return sqrt(hitevents*(ntr-hitevents)/ntr)/ntr;
}
Float_t AliITSOnlineSPDscan::GetHitEventsTotEfficiencyError(UInt_t nsi, UInt_t hs) {
  // get error in hit events efficiency for 10 chips
  return GetHitEventsEfficiencyError(nsi,hs,10);
}
Float_t AliITSOnlineSPDscan::GetAverageMultiplicity(UInt_t nsi, UInt_t hs, UInt_t chipi) {
  // get average multiplicity for a chip
  Float_t nrhits = 0;
  for (UInt_t chip=0;chip<10;chip++) {
    if (chipi==10 || chip==chipi) {
      for (Int_t col=0; col<32; col++) {
	for (Int_t row=0; row<256; row++) {
	  nrhits+=GetHits(nsi,hs,chip,col,row);
	}
      }
    }
  }
  UInt_t ntr = GetTriggers(nsi);
  if (ntr>0) {
    return nrhits/ntr;
  }
  else {
    return 0;
  }
}
Float_t AliITSOnlineSPDscan::GetAverageMultiplicityTot(UInt_t nsi, UInt_t hs) {
  // get average multiplicity for 10 chips
  return GetAverageMultiplicity(nsi,hs,10);
}


void AliITSOnlineSPDscan::SetType(UInt_t val) {
  // set type
  fScanInfo->SetType(val);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetRunNr(UInt_t val) {
  // set run nr
  fScanInfo->SetRunNr(val);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetRouterNr(UInt_t val) {
  // set router nr
  fScanInfo->SetRouterNr(val);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetHalfStaveScanned(UInt_t val, Bool_t b) {
  // set half stave scanned
  fScanInfo->SetHalfStaveScanned(val,b);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetDataFormat(UInt_t val) {
  // set data format (0=normal 1=histogram)
  fScanInfo->SetDataFormat(val);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetTriggers(UInt_t nsi, UInt_t val) {
  // set nr of triggers
  SwitchToStep(nsi);
  fScanInfo->SetTriggers(nsi,val);
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetChipPresent(UInt_t hs, UInt_t chipi, Bool_t val){
  // set chip present
  fScanInfo->SetChipPresent(hs,chipi,val); 
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetRowStart(UInt_t val){
  // set row start
  fScanInfo->SetRowStart(val); 
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetRowEnd(UInt_t val){
  // set row end
  fScanInfo->SetRowEnd(val); 
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetDacStart(UInt_t val){
  // set dac start
  fScanInfo->SetDacStart(val); 
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetDacEnd(UInt_t val){
  // set dac end
  fScanInfo->SetDacEnd(val); 
  fInfoModified=kTRUE;
}  
void AliITSOnlineSPDscan::SetDacStep(UInt_t val){
  // set dac step
  fScanInfo->SetDacStep(val); 
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::SetDCSVersion(UInt_t val){
  // set dcs db version
  fScanInfo->SetDCSVersion(val); 
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscan::IncrementTriggers(UInt_t nsi) {
  // increment nr of triggers
  SwitchToStep(nsi);
  fScanInfo->IncrementTriggers(nsi); 
  fInfoModified=kTRUE;
}



UInt_t AliITSOnlineSPDscan::GetNSteps() const {
  return fScanInfo->GetNSteps();
}
UInt_t AliITSOnlineSPDscan::GetType() const {
  return fScanInfo->GetType();
}
UInt_t AliITSOnlineSPDscan::GetRunNr() const {
  return fScanInfo->GetRunNr();
}
UInt_t AliITSOnlineSPDscan::GetRouterNr() const {
  return fScanInfo->GetRouterNr();
}
Bool_t AliITSOnlineSPDscan::GetHalfStaveScanned(UInt_t val) const {
  return fScanInfo->GetHalfStaveScanned(val);
}
UInt_t AliITSOnlineSPDscan::GetDataFormat() const {
  return fScanInfo->GetDataFormat();
}
UInt_t AliITSOnlineSPDscan::GetTriggers(UInt_t nsi) const {
  return fScanInfo->GetTriggers(nsi);
}
Bool_t AliITSOnlineSPDscan::GetChipPresent(UInt_t hs, UInt_t chipi) const {
  return fScanInfo->GetChipPresent(hs,chipi);
}
UInt_t AliITSOnlineSPDscan::GetRowStart() const {
  return fScanInfo->GetRowStart();
}
UInt_t AliITSOnlineSPDscan::GetRowEnd() const {
  return fScanInfo->GetRowEnd();
}
UInt_t AliITSOnlineSPDscan::GetDacStart() const {
  return fScanInfo->GetDacStart();
}
UInt_t AliITSOnlineSPDscan::GetDacEnd() const {
  return fScanInfo->GetDacEnd();
}
UInt_t AliITSOnlineSPDscan::GetDacStep() const {
  return fScanInfo->GetDacStep();
}
UInt_t AliITSOnlineSPDscan::GetDCSVersion() const {
  return fScanInfo->GetDCSVersion();
}
