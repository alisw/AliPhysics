// -*- C++ -*-

#include "AliSLEventHeader.h"

ClassImp(AliSLEventHeader);

Double_t AliSLEventHeader::GetEventInfo(const char* name) const
{
  const TParameter<double>* p = dynamic_cast<const TParameter<double>* >(fEventInfo->FindObject(name));
  if (p)
    return p->GetVal();

  AliErrorF("key '%s' does not exist", name);
  return 0.0;
}

void AliSLEventHeader::Print(Option_t* opt) const
{
  AliGenEventHeader::Print(opt);
  Printf("STARLIGHT Event info: ");
  if (fEventInfo)
    fEventInfo->Print();
}
