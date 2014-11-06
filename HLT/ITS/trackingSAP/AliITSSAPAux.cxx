#include "AliITSSAPAux.h"
#include <stdio.h>

//_______________________________________________________________
void AliITSSAPAux::PrintBits(unsigned long long patt, int maxBits)
{
  // print maxBits of the pattern
  maxBits = maxBits>64 ? 64:maxBits;
  for (int i=0;i<maxBits;i++) printf("%c",((patt>>i)&0x1) ? '+':'-');
}
