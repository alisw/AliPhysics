#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TInterpreter.h"
#include "MONITOR/AliMonitorProcess.h"
#include "MONITOR/AliMonitorControl.h"
#endif

void monitor(const char* alienDir = ".")
{
  if (!gSystem->Which(".", "galice.root")) {
    gInterpreter->ExecuteMacro("$ALICE_ROOT/MONITOR/galice.C");
  }
  if (!gSystem->Which(".", "Table0.dat")) {
    gSystem->Exec("cp $ALICE_ROOT/RAW/Table*.dat .");
  }
  AliMonitorProcess process(alienDir);
  process.Run();
//  new AliMonitorControl(&process);
}
