////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online scan    //
// with multiple steps.                                   //
////////////////////////////////////////////////////////////

#include <TFile.h>
#include "AliITSOnlineSPDscanMultiple.h"
#include "AliITSOnlineSPDscanInfoMultiple.h"

ClassImp(AliITSOnlineSPDscanMultiple)

AliITSOnlineSPDscanMultiple::AliITSOnlineSPDscanMultiple():AliITSOnlineSPDscan(){
// Default constructor
}
AliITSOnlineSPDscanMultiple::AliITSOnlineSPDscanMultiple(Char_t *fileName) {
  // constructor
  sprintf(fFileName,"%s",fileName);
  // look for a previously saved info object 
  // (if file not found create a new one and return, else read)
  FILE* fp0 = fopen(fFileName, "r");
  if (fp0 == NULL) {
    fScanInfo = new AliITSOnlineSPDscanInfoMultiple();
    fFile = new TFile(fFileName, "RECREATE");
    fWrite=kTRUE;
  }
  else {
    fclose(fp0);
    fFile = new TFile(fFileName, "READ");
    fWrite=kFALSE;
    fFile->GetObject("AliITSOnlineSPDscanInfo", fScanInfo);
  }
  Init();
}

AliITSOnlineSPDscanMultiple::AliITSOnlineSPDscanMultiple(const AliITSOnlineSPDscanMultiple& scan) :
  AliITSOnlineSPDscan(scan)
{}

AliITSOnlineSPDscanMultiple::~AliITSOnlineSPDscanMultiple() {}

AliITSOnlineSPDscanMultiple& AliITSOnlineSPDscanMultiple::operator=(const AliITSOnlineSPDscanMultiple& scan) {
  // Assignment operator, should not be called!!!
  printf("This object should not be copied!");
  if (this!=&scan) {
    // still do nothing...
  }
  return *this;
}

//void AliITSOnlineSPDscanMultiple::ReadFromTObjArray(TObjArray *arr) {
//  ClearThis();
//  Int_t nrEntries = arr->GetEntriesFast();
//  if (nrEntries>0 && nrEntries%2==1) {
//    fScanInfo = (AliITSOnlineSPDscanInfoMultiple*) arr->At(0);
//    fInfoModified=kTRUE;
//    FillFromTObjArray(arr,nrEntries);
//  }
//}

UInt_t AliITSOnlineSPDscanMultiple::AddScanStep() {
  CreateNewStep();
  return ((AliITSOnlineSPDscanInfoMultiple*)fScanInfo)->AddScanStep();
}

void AliITSOnlineSPDscanMultiple::SetDacId(Int_t val) {
  ((AliITSOnlineSPDscanInfoMultiple*)fScanInfo)->SetDacId(val); 
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscanMultiple::SetDacValue(UInt_t nsi, Int_t val) {
  // set dac value for step nsi
  SwitchToStep(nsi);
  ((AliITSOnlineSPDscanInfoMultiple*)fScanInfo)->SetDacValue(nsi,val); 
  fInfoModified=kTRUE;
}
Int_t AliITSOnlineSPDscanMultiple::GetDacId() {
  return ((AliITSOnlineSPDscanInfoMultiple*)fScanInfo)->GetDacId();
}
Int_t AliITSOnlineSPDscanMultiple::GetDacValue(UInt_t nsi) {
  return ((AliITSOnlineSPDscanInfoMultiple*)fScanInfo)->GetDacValue(nsi);
}
