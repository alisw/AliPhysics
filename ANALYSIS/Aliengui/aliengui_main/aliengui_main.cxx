#include <TSystem.h>
#include <TApplication.h>
#ifdef WITHXML
#include <Aliengui/AliAnalysisGUI.h>
#else
#include <Aliengui/AliAnalysisGUIdummy.h>
#endif
int main(int argc, char **argv) {
  // main
  
  TApplication theApp("App",&argc,argv);
  new AliAnalysisGUI(gClient->GetRoot(),600,600);
  
  theApp.Run();
  
  return 0;
}
