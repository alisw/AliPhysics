#include "AliITSUAux.h"

//_______________________________________________________________
void AliITSUAux::PrintBits(ULong64_t patt, Int_t maxBits)
{
  // print maxBits of the pattern
  maxBits = Min(64,maxBits);
  for (int i=0;i<maxBits;i++) printf("%c",((patt>>i)&0x1) ? '+':'-');
}
