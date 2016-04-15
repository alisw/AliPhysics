#include <TApplication.h>
#include <TGClient.h>
#include "AliZMQMTviewerGUI.h"

int main(int argc, char **argv) {
   TApplication theApp("App", &argc, argv);
   AliZMQMTviewerGUI mymain (gClient->GetRoot(), 200, 200, argc, argv);
   theApp.Run();
   return 0;
}

