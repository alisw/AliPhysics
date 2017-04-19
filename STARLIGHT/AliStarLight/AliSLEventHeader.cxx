// -*- C++ -*-

#include <TNamed.h>
#include <TVector2.h>
#include <TParameter.h>

#include "AliSLEventHeader.h"

ClassImp(AliSLEventHeader);

Double_t AliSLEventHeader::GetEventInfo(const char* name, const type2val<Double_t>&) const
{
  const TParameter<double>* p = dynamic_cast<const TParameter<double>* >(fEventInfo->GetValue(name));
  if (p) {
    return p->GetVal();
  }
  AliErrorF("key '%s' does not exist", name);
  return 0.0;
}
TVector2 AliSLEventHeader::GetEventInfo(const char* name, const type2val<TVector2>&) const
{
  const TVector2* v = dynamic_cast<const TVector2* >(fEventInfo->GetValue(name));
  if (v) {
    return TVector2(*v);
  }
  AliErrorF("key '%s' does not exist", name);
  return TVector2(0,0);
}

void AliSLEventHeader::Print(Option_t* opt) const
{
  AliGenEventHeader::Print(opt);
  Printf("STARLIGHT Event info: ");
  if (fEventInfo)
    fEventInfo->Print();
}
