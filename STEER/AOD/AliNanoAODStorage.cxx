#include <TObjString.h>
#include "AliNanoAODStorage.h"
#include "AliLog.h"

ClassImp(AliNanoAODStorage)

void AliNanoAODStorage::AllocateInternalStorage(Int_t size) {
  AllocateInternalStorage(size, 0);
}

void AliNanoAODStorage::AllocateInternalStorage(Int_t size, Int_t sizeInt) {
  /// Creates the internal array

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

  if(sizeInt>0){
    fNVarsInt = sizeInt;
    fVarsInt.clear();
    fVarsInt.resize(sizeInt, 0);
  }
  
}

AliNanoAODStorage& AliNanoAODStorage::operator=(const AliNanoAODStorage& sto)
{
  /// Assignment operator

  AllocateInternalStorage(sto.fNVars, sto.fNVarsInt);
  if(this!=&sto) {
    for (Int_t isize = 0; isize<sto.fNVars; isize++) {
      SetVar(isize, sto.GetVar(isize));    
    }
    for (Int_t isize = 0; isize<sto.fNVarsInt; isize++) {
      SetVarInt(isize, sto.GetVarInt(isize));    
    }
    
  }

  return *this;
}

Int_t AliNanoAODStorage::GetIntParameters(const TString varListHeader){
  const TString stringVariables = "FiredTriggerClasses,BunchCrossNumber,OrbitNumber,PeriodNumber,NumberOfESDTracks,OfflineTrigger,RunNumber";//list of all possible string variables in AliNanoAODStorage

  TObjArray * vars = varListHeader.Tokenize(",");
  TIter it(vars);
  TObjString *token  = 0;
  Int_t stringVars=0;

  while ((token = (TObjString*) it.Next())) {
    TString var = token->GetString().Strip(TString::kBoth, ' ');
    if(stringVariables.Contains(var))
      stringVars++;
  }  
  
  return stringVars;
}

void
AliNanoAODStorage::Complain(Int_t index) const {
  AliFatal(Form("Variable %d not included in this special aod", index));
}

