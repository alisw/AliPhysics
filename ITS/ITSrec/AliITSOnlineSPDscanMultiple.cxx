////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online scan    //
// with multiple steps.                                   //
////////////////////////////////////////////////////////////

#include <TFile.h>
#include "AliITSOnlineSPDscanMultiple.h"
#include "AliITSOnlineSPDscanInfoMultiple.h"

AliITSOnlineSPDscanMultiple::AliITSOnlineSPDscanMultiple():AliITSOnlineSPDscan(){
// Default constructor
}
AliITSOnlineSPDscanMultiple::AliITSOnlineSPDscanMultiple(const Char_t *fileName, Bool_t readFromGridFile) {
  // constructor
  fFileName=fileName;
  fModified=kFALSE;
  fInfoModified=kFALSE;
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
      fScanInfo = new AliITSOnlineSPDscanInfoMultiple();
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
    fScanInfo = new AliITSOnlineSPDscanInfoMultiple();
    fInfoModified=kTRUE;
    fFile = new TFile(fFileName.Data(), "RECREATE");
    fWrite=kTRUE;
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
