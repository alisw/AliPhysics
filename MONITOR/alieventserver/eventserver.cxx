// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#include <TApplication.h>
#include "AliEventServerWindow.h"
#include "AliEventServer.h"
#include <iostream>
#include <string.h>

int main(int argc, char **argv)
{
  TApplication app("AliEventServer", &argc, argv);

  if(argc<2)
    {
 std::cout<<"Starting Event Server without GUI"<<std::endl;
      AliEventServer *server = new AliEventServer;
      app.Run(kTRUE);
      if(server){delete server;}
    }
  else if(strcmp(argv[1],"gui")==0)
    {
      std::cout<<"Starting Event Server in GUI mode"<<std::endl;
  AliEventServerWindow *win = new AliEventServerWindow;
  app.Run(kTRUE);
  
  if(win){delete win;}
    }
  else
    {
      std::cout<<"Call without parameters to run without GUI.\nCall with \"gui\" parameter to launch with GUI"<<std::endl;
    }

  return 0;
}
