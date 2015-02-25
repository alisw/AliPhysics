#include "AliVersion.h"
#include "ARVersion.h"
#include <TClass.h>

AliVersion* AliVersion::fgInstance = 0;

AliVersion* AliVersion::Instance()
{
  if (!fgInstance) fgInstance = new AliVersion;
  return (AliVersion*)fgInstance;
}

AliVersion::AliVersion()
  : TNamed("alirootVersion", "AliROOT Version"),
    fHash(ALIROOT_REVISION),
    fTag(ALIROOT_VERSION)
{   SetUniqueID(ALIROOT_SERIAL); }

Int_t AliVersion::Compare(const TObject* o) const
{
   if (!o->IsA()->InheritsFrom(AliVersion::Class()))
     Fatal("Compare", "Cannot compare an AliVersion object to a %s object",
           o->IsA()->GetName());
   const AliVersion* av = static_cast<const AliVersion*>(o);
   return (av->GetSerial() == GetSerial() ? 0 : 
           av->GetSerial() >  GetSerial() ? -1 : 1);
}

void AliVersion::Print(Option_t *) const
{
  // print aliroot version
  printf("AliRoot serial:\t%d\nAliRoot hash:\t%s\nAliRoot tag:\t%s\n",GetSerial(),fTag.Data(),fHash.Data());
}
