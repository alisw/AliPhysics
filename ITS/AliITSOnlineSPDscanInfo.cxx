/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// It holds information needed for a scan.                     //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscan class.                                  //
/////////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDscanInfo.h"

ClassImp(AliITSOnlineSPDscanInfo)

AliITSOnlineSPDscanInfo::AliITSOnlineSPDscanInfo(): 
  fType(999),
  fRunNr(0),
  fRouterNr(999),
  fNSteps(0),
  fTriggers(0),
  fRowStart(0),
  fRowEnd(255),
  fDacStep(1),
  fDacStart(0),
  fDacEnd(255)
{
  ClearThis();
}

AliITSOnlineSPDscanInfo::~AliITSOnlineSPDscanInfo() {
}

void AliITSOnlineSPDscanInfo::ClearThis() {
  // reset all values for this object
  fNSteps=0;
  fTriggers=0;
  fType=999;
  fRunNr=0;
  fRouterNr=999;
  fRowStart=0;
  fRowEnd=255;
  for (Int_t i=0; i<10; i++) {
    fChipPresent[i]=kTRUE;
  }
}

UInt_t AliITSOnlineSPDscanInfo::AddScanStep() {
  // add a new scan step, allocate space for TArrayI
  fNSteps++;
  fTriggers.Set(fNSteps);
  fTriggers.AddAt(0, fNSteps-1);
  return fNSteps-1;
}

void AliITSOnlineSPDscanInfo::IncrementTriggers(UInt_t nsi) {
  // increment the nr of triggers for step nsi
  if (nsi<fNSteps) {
    fTriggers.AddAt(GetTriggers(nsi)+1,nsi);
  }
}
void AliITSOnlineSPDscanInfo::SetTriggers(UInt_t nsi, UInt_t val) {
  // set the nr of triggers for step nsi
  if (nsi<fNSteps) {
    fTriggers.AddAt(val,nsi);
  }
}

UInt_t AliITSOnlineSPDscanInfo::GetTriggers(UInt_t nsi) const {
  // get the nr of triggers for step nsi
  if (nsi<fNSteps) return fTriggers.At(nsi);
  else return 0;
}
