#include "AliDimIntNotifier.h"
#include "AliOnlineReco.h"

#include <TRint.h>

int main(int argc, char **argv)
{
  AliDimIntNotifier::SetMainThreadId();

  Bool_t test = argc > 1 && strcmp("-test", argv[1]) == 0;

  TRint app("App", &argc, argv);

  AliOnlineReco *win = new AliOnlineReco;
  win->MapWindow();

  if (test)
  {
    win->SetTestMode();

    win->GetSOR()->infoHandlerTest(2214);
    win->GetSOR()->infoHandlerTest(2215);
    win->GetSOR()->infoHandlerTest(2224);
    win->GetSOR()->infoHandlerTest(2244);

    printf("win = (AliOnlineReco*) 0x%lx\n", (unsigned long)win);
  }

  app.Run(kTRUE);
  return 0;
}
