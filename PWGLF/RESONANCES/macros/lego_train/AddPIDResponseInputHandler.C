#ifndef __CINT__
#include <AliPIDResponseInputHandler.h>
#endif

void AddPIDResponseInputHandler(AliMultiInputEventHandler *multiInputHandler,Bool_t isMC=kFALSE)
{
   if (multiInputHandler) {
      AliPIDResponseInputHandler *pidResponseIH = new AliPIDResponseInputHandler();
      pidResponseIH->SetIsMC(isMC);
      multiInputHandler->AddInputEventHandler(pidResponseIH);
   }
}
