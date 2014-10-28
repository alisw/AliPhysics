// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#include <TApplication.h>
#include <TSystem.h>
#include <TString.h>
#include <TROOT.h>
#include "AliEventServerWindow.h"
#include "AliEventServer.h"
#include <iostream>
#include <string.h>
#include <sstream>

int main(int argc, char **argv)
{
  // check if there is events server already running
  const char *pid = gSystem->GetFromPipe("pidof alieventserver").Data();
  int pidSize = gSystem->GetFromPipe("pidof alieventserver").Sizeof();
  std::string pidOfAll(pid,pidSize);
  std::stringstream pidStream(pidOfAll);

  int word_count=0; 
  std::string word;
  while( pidStream >> word ) ++word_count;

  if(word_count != 1)
    {
      std::cout<<"There are other servers. Cannot start multiple servers on the same machine. Quitting..."<<std::endl;
      return 0;
    }

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
