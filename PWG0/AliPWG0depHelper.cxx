/* $Id$ */

#include <AliPWG0depHelper.h>

#include <AliHeader.h>

#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>

//____________________________________________________________________
ClassImp(AliPWG0depHelper)

//____________________________________________________________________
const Int_t AliPWG0depHelper::GetPythiaEventProcessType(AliHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());

  if (!pythiaGenHeader) {

    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(aHeader->GenEventHeader());
    if (!genCocktailHeader) {
      printf("AliPWG0depHelper::GetProcessType : Unknown header type (not Pythia or Cocktail). \n");
      return -1;
    }

    TList* headerList = genCocktailHeader->GetHeaders();
    if (!headerList) {
      return -1;
    }

    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }

    if (!pythiaGenHeader) {
      printf("AliPWG0depHelper::GetProcessType : Could not find Pythia header. \n");
      return -1;
    }
  }

  if (adebug) {
    printf("AliPWG0depHelper::GetProcessType : Pythia process type found: %d \n",pythiaGenHeader->ProcessType());
  }

  return pythiaGenHeader->ProcessType();
}
