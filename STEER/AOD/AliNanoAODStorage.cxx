#include "AliNanoAODStorage.h"
#include "AliNanoAODTrackMapping.h"
#include "AliLog.h"

ClassImp(AliNanoAODStorage)

void AliNanoAODStorage::AllocateInternalStorage(Int_t size) {
  // Creates the internal array
  if(size == 0){
    AliError("Zero size");
    return;
  }
  fNVars = size;
  fVars.clear();
  fVars.resize(size, 0);
  // if(fVars) {
  //   delete[] fVars;
  // }
  // fVars = new Double_t[fNVars];
  // for (Int_t ivar = 0; ivar<fNVars; ivar++) {
  //   fVars[ivar]=0;
  // }
}

AliNanoAODStorage& AliNanoAODStorage::operator=(const AliNanoAODStorage& sto)
{
  // Assignment operator
  AllocateInternalStorage(sto.fNVars);
  if(this!=&sto) {
    for (Int_t isize = 0; isize<sto.fNVars; isize++) {
      SetVar(isize, sto.GetVar(isize));    
    }
    
  }

  return *this;
}

void
AliNanoAODStorage::Complain(Int_t index) const {
  AliFatal(Form("Variable %d not included in this special aod", index));
}

