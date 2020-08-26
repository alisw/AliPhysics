#include <TObjString.h>
#include "AliNanoAODHeader.h"
#include "AliLog.h"

ClassImp(AliNanoAODHeader)

Bool_t AliNanoAODHeader::fFatalMode = kFALSE;

AliNanoAODHeader::AliNanoAODHeader():
  AliVAODHeader(),
  AliNanoAODStorage(),
  fCentralityMethod("V0M"),
  fBunchCrossNumber(-1),
  fOrbitNumber(-1),
  fPeriodNumber(-1),
  fCentr(-1),
  fCentrTRK(-1),
  fCentrCL0(-1),
  fCentrCL1(-1),
  fFiredTriggerClasses(-1),
  fMagField(-1),
  fOfflineTrigger(-1),
  fRunNumber(-1),
  fNumberOfESDTracks(-1)
{

  // default ctor
  
  for (int i=0; i<4; i++)
    fT0Spread[i] = -1;
}

AliNanoAODHeader::AliNanoAODHeader(Int_t size):
  AliVAODHeader(),
  AliNanoAODStorage(),
  fCentralityMethod("V0M"),
  fBunchCrossNumber(-1),
  fOrbitNumber(-1),
  fPeriodNumber(-1),
  fCentr(-1),
  fCentrTRK(-1),
  fCentrCL0(-1),
  fCentrCL1(-1),
  fFiredTriggerClasses(-1),
  fMagField(-1),
  fOfflineTrigger(-1),
  fRunNumber(-1),
  fNumberOfESDTracks(-1)
{
  for (int i=0; i<4; i++)
    fT0Spread[i] = -1;

  AllocateInternalStorage(size, 0);
}

AliNanoAODHeader::AliNanoAODHeader(Int_t size, Int_t sizeInt):
  AliVAODHeader(),
  AliNanoAODStorage(),
  fCentralityMethod("V0M"),
  fBunchCrossNumber(-1),
  fOrbitNumber(-1),
  fPeriodNumber(-1),
  fCentr(-1),
  fCentrTRK(-1),
  fCentrCL0(-1),
  fCentrCL1(-1),
  fFiredTriggerClasses(-1),
  fMagField(-1),
  fOfflineTrigger(-1),
  fRunNumber(-1),
  fNumberOfESDTracks(-1)
{
  for (int i=0; i<4; i++)
    fT0Spread[i] = -1;
  
  AllocateInternalStorage(size, sizeInt);
}

AliNanoAODHeader& AliNanoAODHeader::operator=(const AliNanoAODHeader& evt) {

    AliVHeader::operator=(evt); // FIXME: ok?
    AliNanoAODStorage::operator=(evt);

    return *this;
}


void  AliNanoAODHeader::Clear(Option_t * /*opt*/) {
  /// empty storage

  fVars.clear();
  fVarsInt.clear();
  fNVars = 0;
  fNVarsInt = 0;
}

Double_t AliNanoAODHeader::GetCentrality () const {
    if (fCentralityMethod=="V0M") return GetCentr("V0M");
    if (fCentralityMethod=="TRK") return GetCentr("TRK");
    if (fCentralityMethod=="CL1") return GetCentr("CL1");
    if (fCentralityMethod=="CL0") return GetCentr("CL0");
    return -1;
}

Double_t AliNanoAODHeader::GetCentr (const char *x) const {
    TString method = x;
    if(method.CompareTo("V0M")==0)      return GetVar(fCentr);
    if(method.CompareTo("TRK")==0)      return GetVar(fCentrTRK);
    if(method.CompareTo("CL1")==0)      return GetVar(fCentrCL1);
    if(method.CompareTo("CL0")==0)      return GetVar(fCentrCL0);
    return -1;
} 

Int_t  AliNanoAODHeader::GetRunNumber() const { 
   if (fRunNumber>-1) return GetVarInt(fRunNumber);
   return 0;

} 


void AliNanoAODHeader::SetFiredTriggerClasses(TString varlist) {
  int firedTrigClasses = 0;

  TObjArray * vars = varlist.Tokenize("  ");
  TIter it(vars);
  TObjString *token  = 0;
  
  if(fMapFiredTriggerClasses.size()==0){
    AliFatal("fMapFiredTriggerClasses does not exist. Cannot set fired trigger classes.");
  }
  
  while ((token = (TObjString*) it.Next())) {
    TString var = token->GetString();
    std::map<TString,Int_t>::iterator it = fMapFiredTriggerClasses.find(var); // FIXME: do I need to delete "it"?
    if(it != fMapFiredTriggerClasses.end()) {
      //element found;
      firedTrigClasses |= 1 << it->second;
    }else{
      //ignore missing fired trigger classes
    } 
  }

  //create object to save these triggers
  fVarsInt[fFiredTriggerClasses] = firedTrigClasses;
}

TString  AliNanoAODHeader::GetFiredTriggerClasses() const {
  TString firedTrigClasses = "";

  if(fMapFiredTriggerClasses.size()==0){
    AliFatal("fMapFiredTriggerClasses does not exist. Cannot get fired trigger classes.");
  }

  for (std::map<TString, int>::const_iterator it = fMapFiredTriggerClasses.begin(); it != fMapFiredTriggerClasses.end(); it++){
    int bit = (fVarsInt[fFiredTriggerClasses] >> it->second) & 1;
    if(bit==1){
      if(firedTrigClasses.Length()>1)
	firedTrigClasses += "  ";
      firedTrigClasses += it->first;
    }
  }

  return firedTrigClasses;
  
}

void AliNanoAODHeader::SetMapFiredTriggerClasses (TString trigClasses){

  TObjArray * vars = trigClasses.Tokenize(",");
  TIter it(vars);
  TObjString *token  = 0;
  Int_t index=0;

  while ((token = (TObjString*) it.Next())) {
    TString var = token->GetString().Strip(TString::kBoth, ' ');
    fMapFiredTriggerClasses[var] = index;
    index++;
  }
}


Int_t AliNanoAODHeader::GetVarIndex(TString varName){
  std::map<TString,Int_t>::iterator it = fMapCstVar.find(varName); // FIXME: do I need to delete "it"?
  if(it != fMapCstVar.end()) {
    //element found;
    return it->second;
  }else{
    return -1;
  } 
} 

void AliNanoAODHeader::NotImplemented(void) const {
  if (fFatalMode)
    AliFatal("Not implemented");
  else
    AliError("Not implemented");
}
