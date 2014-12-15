/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// It holds information needed for a scan.                     //
// This class should only be used through the interface of the //
// AliITSOnlineSPDphys class.                                  //
/////////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDphysInfo.h"

ClassImp(AliITSOnlineSPDphysInfo)

AliITSOnlineSPDphysInfo::AliITSOnlineSPDphysInfo(): 
  fNrRuns(0),
  fRunNrs(0),
  fEqNr(999),
  fNrEvents(0)
{}

AliITSOnlineSPDphysInfo::~AliITSOnlineSPDphysInfo() {}

void AliITSOnlineSPDphysInfo::ClearThis() {
  // reset all values for this object
  fNrRuns=0;
  fRunNrs=0;
  fEqNr=999;
  fNrEvents=0;
}

void AliITSOnlineSPDphysInfo::AddRunNr(UInt_t val) {
  // add a new run nr, allocate space for TArrayI
  fNrRuns++;
  fRunNrs.Set(fNrRuns);
  fRunNrs.AddAt(val, fNrRuns-1);
}

UInt_t AliITSOnlineSPDphysInfo::GetRunNr(UInt_t posi) const {
  // get run nr
  if (posi<fNrRuns) {
    return fRunNrs.At(posi);
  }
  else {
    return 0;
  }
}

void AliITSOnlineSPDphysInfo::AddNrEvents(Int_t val) {
  // add val nr of events (val could be negative)
  if (fNrEvents+val>0) {
    fNrEvents+=val;
  }
  else {
    fNrEvents=0;
  }
}
