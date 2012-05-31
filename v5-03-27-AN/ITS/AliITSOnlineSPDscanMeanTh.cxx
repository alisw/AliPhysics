////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online mean    //
// threshold scan.                                        //
////////////////////////////////////////////////////////////

#include <TFile.h>
#include "AliITSOnlineSPDscanMeanTh.h"
#include "AliITSOnlineSPDscanInfoMeanTh.h"

AliITSOnlineSPDscanMeanTh::AliITSOnlineSPDscanMeanTh(const Char_t *fileName, Bool_t readFromGridFile) {
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
      fScanInfo = new AliITSOnlineSPDscanInfoMeanTh();
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
    fScanInfo = new AliITSOnlineSPDscanInfoMeanTh();
    fInfoModified=kTRUE;
    fFile = new TFile(fFileName.Data(), "RECREATE");
    fWrite=kTRUE;
  }

  Init();
}

AliITSOnlineSPDscanMeanTh::AliITSOnlineSPDscanMeanTh(const AliITSOnlineSPDscanMeanTh& scan) :
  AliITSOnlineSPDscanMultiple(scan)
{}

AliITSOnlineSPDscanMeanTh::~AliITSOnlineSPDscanMeanTh() {}

AliITSOnlineSPDscanMeanTh& AliITSOnlineSPDscanMeanTh::operator=(const AliITSOnlineSPDscanMeanTh& scan) {
  // assignment operator (should not be used)
  printf("This object should not be copied!");
  if (this!=&scan) {
    // still do nothing...
  }
  return *this;
}

UInt_t AliITSOnlineSPDscanMeanTh::AddScanStep() {
  CreateNewStep();
  return ((AliITSOnlineSPDscanInfoMeanTh*)fScanInfo)->AddScanStep();
}

void AliITSOnlineSPDscanMeanTh::SetDacLow(UInt_t nsi, UInt_t hs, Int_t val) {
  // set dac low value for step nsi and half stave hs
  SwitchToStep(nsi);
  ((AliITSOnlineSPDscanInfoMeanTh*)fScanInfo)->SetDacLow(nsi,hs,val); 
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscanMeanTh::SetDacHigh(UInt_t nsi, UInt_t hs, Int_t val) {
  // set dac high value for step nsi and half stave hs
  SwitchToStep(nsi);
  ((AliITSOnlineSPDscanInfoMeanTh*)fScanInfo)->SetDacHigh(nsi,hs,val); 
  fInfoModified=kTRUE;
}
void AliITSOnlineSPDscanMeanTh::SetTPAmp(UInt_t nsi, UInt_t hs, Int_t val) {
  // set test pulse amplitude for step nsi and half stave hs
  SwitchToStep(nsi);
  ((AliITSOnlineSPDscanInfoMeanTh*)fScanInfo)->SetTPAmp(nsi,hs,val); 
  fInfoModified=kTRUE;
}

Int_t AliITSOnlineSPDscanMeanTh::GetDacLow(UInt_t nsi, UInt_t hs) {
  return ((AliITSOnlineSPDscanInfoMeanTh*)fScanInfo)->GetDacLow(nsi,hs);
}
Int_t AliITSOnlineSPDscanMeanTh::GetDacHigh(UInt_t nsi, UInt_t hs) {
  return ((AliITSOnlineSPDscanInfoMeanTh*)fScanInfo)->GetDacHigh(nsi,hs);
}
Int_t AliITSOnlineSPDscanMeanTh::GetTPAmp(UInt_t nsi, UInt_t hs) {
  return ((AliITSOnlineSPDscanInfoMeanTh*)fScanInfo)->GetTPAmp(nsi,hs);
}
