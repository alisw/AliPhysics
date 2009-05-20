
#include "AliESDv0.h"
#include "AliTRDv0Info.h"

ClassImp(AliTRDv0Info)

//_________________________________________________
AliTRDv0Info::AliTRDv0Info()
  : TObject()
  ,fStatus(0)
{
}

//_________________________________________________
AliTRDv0Info::AliTRDv0Info(AliESDv0 */*v0*/)
  : TObject()
  ,fStatus(0)
{
}


//_________________________________________________
void AliTRDv0Info::Print(Option_t */*opt*/) const
{
  printf("AliTRDv0Info::Print() : Nothing implemented so far.\n");
}
