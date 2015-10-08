// We need the.cxx for the classimp and for CMake to compile the class
#include "AliVAODHeader.h"
#include "AliLog.h"

ClassImp(AliVAODHeader);

//_________________________________________
void AliVAODHeader::SetDAQAttributes(UInt_t)
{
  // warn
  AliWarning("not implmented");
}

//_________________________________________
UInt_t AliVAODHeader::GetDAQAttributes() const
{
  // warn
  AliWarning("not implmented");
  return 0;
}
