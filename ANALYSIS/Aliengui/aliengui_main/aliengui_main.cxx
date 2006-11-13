#include <TSystem.h>
#include <TApplication.h>
#include <Aliengui/AliAnalysisGUI.h>

int main(int argc, char **argv) {
  // main
  
  TApplication theApp("App",&argc,argv);
  new AliAnalysisGUI(gClient->GetRoot(),600,600);
  
  theApp.Run();
  
  return 0;
}
