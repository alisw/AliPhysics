/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// It holds additional information needed for a mean threshold //
// scan.                                                       //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscanMeanTh class.                            //
/////////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDscanInfoMeanTh.h"

ClassImp(AliITSOnlineSPDscanInfoMeanTh)

AliITSOnlineSPDscanInfoMeanTh::AliITSOnlineSPDscanInfoMeanTh() :
  AliITSOnlineSPDscanInfoMultiple() {}

AliITSOnlineSPDscanInfoMeanTh::~AliITSOnlineSPDscanInfoMeanTh() {}

UInt_t AliITSOnlineSPDscanInfoMeanTh::AddScanStep() {
  // add a new scan step, allocate space in the TArrays
  UInt_t returnval = AliITSOnlineSPDscanInfoMultiple::AddScanStep();
  for (UInt_t hs=0; hs<6; hs++) {
    fDacLow[hs].Set(fNSteps);
    fDacLow[hs].AddAt(-1, fNSteps-1);
    fDacHigh[hs].Set(fNSteps);
    fDacHigh[hs].AddAt(-1, fNSteps-1);
    fTPAmps[hs].Set(fNSteps);
    fTPAmps[hs].AddAt(-1, fNSteps-1);
  }
  return returnval;
}

void AliITSOnlineSPDscanInfoMeanTh::SetDacLow(UInt_t nsi, UInt_t hs, Int_t val) {
  if (nsi<fNSteps) fDacLow[hs].AddAt(val,nsi);
}
void AliITSOnlineSPDscanInfoMeanTh::SetDacHigh(UInt_t nsi, UInt_t hs, Int_t val) {
  if (nsi<fNSteps) fDacHigh[hs].AddAt(val,nsi);
}
void AliITSOnlineSPDscanInfoMeanTh::SetTPAmp(UInt_t nsi, UInt_t hs, Int_t val) {
  if (nsi<fNSteps) fTPAmps[hs].AddAt(val,nsi);
}
Int_t AliITSOnlineSPDscanInfoMeanTh::GetDacLow(UInt_t nsi, UInt_t hs) const {
  if (nsi<fNSteps) return fDacLow[hs].At(nsi);
  else return -1;
}
Int_t AliITSOnlineSPDscanInfoMeanTh::GetDacHigh(UInt_t nsi, UInt_t hs) const {
  if (nsi<fNSteps) return fDacHigh[hs].At(nsi);
  else return -1;
}
Int_t AliITSOnlineSPDscanInfoMeanTh::GetTPAmp(UInt_t nsi, UInt_t hs) const {
  if (nsi<fNSteps) return fTPAmps[hs].At(nsi);
  else return -1;
}
