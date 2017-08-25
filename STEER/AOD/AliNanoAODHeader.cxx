#include "AliNanoAODHeader.h"
#include "AliLog.h"

ClassImp(AliNanoAODHeader)

AliNanoAODHeader::AliNanoAODHeader():
  AliVAODHeader(),
  AliNanoAODStorage(),
  fCentralityMethod("V0M"),
  fCentr(-1),
  fCentrTRK(-1),
  fCentrCL0(-1),
  fCentrCL1(-1),
  fMagField(-1),
  fRunNumber(-1)
{

  // default ctor

}

AliNanoAODHeader::AliNanoAODHeader(Int_t size):
  AliVAODHeader(),
  AliNanoAODStorage(),
  fCentralityMethod("V0M"),
  fCentr(-1),
  fCentrTRK(-1),
  fCentrCL0(-1),
  fCentrCL1(-1),
  fMagField(-1),
  fRunNumber(-1)
{

AllocateInternalStorage(size);

}

AliNanoAODHeader& AliNanoAODHeader::operator=(const AliNanoAODHeader& evt) {

    AliVHeader::operator=(evt); // FIXME: ok?
    AliNanoAODStorage::operator=(evt);

    return *this;
}


void  AliNanoAODHeader::Clear(Option_t * /*opt*/) {
  // empty storage
  fVars.clear();
  fNVars = 0;
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
   if (fRunNumber>0) return Int_t(GetVar(fRunNumber));
   return 0;

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
  AliError("Not implemented");
}
