#include "AliMC.h"

ClassImp(AliMC)

AliMC* AliMC::fgMC=0;

AliMC* gMC;

AliMC::AliMC()
{
}

AliMC::AliMC(const char *name, const char *title) : TNamed(name,title)
{
  if(fgMC) {
    printf("Cannot initialise twice MonteCarlo class\n");
  } else {
    fgMC=this;
    gMC=this;
  }
}

