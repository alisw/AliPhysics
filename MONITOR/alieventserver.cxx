// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#include <TApplication.h>
#include "AliEventServerWindow.h"

int main(int argc, char **argv)
{
  TApplication app("AliEventServer", &argc, argv);

  AliEventServerWindow *win = new AliEventServerWindow;
  app.Run(kTRUE);
  
  if(win) delete win;
  win=0;
  
  return 0;
}
