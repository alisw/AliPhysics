#include <TApplication.h>
#include <TGClient.h>
#include "AliZMQMTviewerGUI.h"

int main(int argc, char **argv) {
   TApplication theApp("App", &argc, argv);
   AliZMQMTviewerGUI viewer(gClient->GetRoot(), 200, 200, argc, argv);
   if (viewer.GetInitStatus()<0) {
     printf("%s",viewer.fUSAGE);
     return 1;
   }
   theApp.Run();
   return 0;
}

