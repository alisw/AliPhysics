#include "AliTestChildProc.h"

#include <TRint.h>
#include <TSystem.h>

#include <TGFrame.h>

#include <cstdlib>

int main(int argc, char **argv)
{
  // Crash on SEGV.
  gSystem->IgnoreSignal(kSigSegmentationViolation, true);

  TRint app("App", &argc, argv);

  TGMainFrame* mf = new AliTestChildProc(atoi(argv[1]));
  mf->MapWindow();

  app.Run(kTRUE);
  return 0;
}
