#ifndef __CINT__
#include <AliPIDResponseInputHandler.h>
#endif

void AddPIDResponseInputHandler(AliMultiInputEventHandler *multiInputHandler)
{
   if (multiInputHandler) {
      AliPIDResponseInputHandler *pidResponseIH = new AliPIDResponseInputHandler();
      multiInputHandler->AddInputEventHandler(pidResponseIH);
   }
}
