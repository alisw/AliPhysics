/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// It holds additional information needed for a scan with      //
// multiple steps. (dac scan, min thr. mean thr. etc.          //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscanMultiple class.                          //
/////////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDscanInfoMultiple.h"

ClassImp(AliITSOnlineSPDscanInfoMultiple)

AliITSOnlineSPDscanInfoMultiple::AliITSOnlineSPDscanInfoMultiple() :
  AliITSOnlineSPDscanInfo(), fDacId(-1), fDacValues(TArrayI())
{}

AliITSOnlineSPDscanInfoMultiple::~AliITSOnlineSPDscanInfoMultiple() {}

UInt_t AliITSOnlineSPDscanInfoMultiple::AddScanStep() {
  // add a new scan step, allocate space in the TArrayI
  UInt_t returnval = AliITSOnlineSPDscanInfo::AddScanStep();
  fDacValues.Set(fNSteps);
  fDacValues.AddAt(-1, fNSteps-1);
  return returnval;
}

void AliITSOnlineSPDscanInfoMultiple::SetDacValue(UInt_t nsi, Int_t val) {
  // set the dac value for step nsi
  if (nsi<fNSteps) {
    fDacValues.AddAt(val, nsi);
  }
}

Int_t AliITSOnlineSPDscanInfoMultiple::GetDacValue(UInt_t nsi) const {
  // get the dac value for step nsi
  if (nsi<fNSteps) {
    return fDacValues.At(nsi);
  }
  else return -1;
}
