#include "AliESDHeader.h"
#include "AliESDRun.h"

#include "AliTRDeventInfo.h"

ClassImp(AliTRDeventInfo)

AliTRDeventInfo::AliTRDeventInfo():
  TObject()
  ,fHeader(0x0)
  ,fRun(0x0)
{
  //
  // Default Constructor
  // 
  SetBit(kOwner, 0);
}

AliTRDeventInfo::AliTRDeventInfo(AliESDHeader *header, AliESDRun *run):
  TObject()
  ,fHeader(header)
  ,fRun(run)
{
  //
  // Constructor with Arguments
  //
  SetBit(kOwner, 0);
}

AliTRDeventInfo::AliTRDeventInfo(const AliTRDeventInfo &info):
  TObject()
  ,fHeader(info.fHeader)
  ,fRun(info.fRun)
{
  //
  // Copy Constructor
  // Flat Copy
  // 
  SetBit(kOwner, 0);
}

AliTRDeventInfo& AliTRDeventInfo::operator=(const AliTRDeventInfo& info){
  //
  // Operator=
  // Flat Copy
  //
  this->fHeader = info.fHeader;
  this->fRun = info.fRun;
  SetBit(kOwner, 0);
  return *this;
}

AliTRDeventInfo::~AliTRDeventInfo(){
  //
  // Destructor
  // Delete the entries if it is the Owner
  //
  Delete("");
}

void AliTRDeventInfo::Delete(const Option_t *){
  //
  // Delete the Object
  // Delete the entries if it is the Owner
  // 
  if(IsOwner()){
    if(fHeader) delete fHeader;
    if(fRun) delete fRun;
  };
  fHeader = 0x0;
  fRun = 0x0;
}

void AliTRDeventInfo::SetOwner(){
  SetBit(kOwner, 1);
  // Do deep copy
  fHeader = new AliESDHeader(*fHeader);
  fRun = new AliESDRun(*fRun);
}
