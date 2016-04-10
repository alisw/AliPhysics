#include "TApplication.h"
#include "signal.h"
#include "AliZMQhistViewer.h"
#include "TSystem.h"
#include "zmq.h"

const char* fUSAGE = 
"ZMQhstViewer: Draw() all ROOT drawables in a message\n"
"options: \n"
" -in : data in\n"
" -sleep : how long to sleep in between requests for data in s (if applicable)\n"
" -timeout : how long to wait for the server to reply (s)\n"
" -Verbose : be verbose\n"
" -select : only show selected histograms (by regexp)\n"
" -unselect : as select, only inverted\n"
" -drawoptions : what draw option to use\n"
" -log[xyz] : use log scale on [xyz] dimension\n"
" -histstats : histogram stat box options (default 0)\n"
" -AllowResetAtSOR : 0/1 to reset at change of run\n"
" -sort : 0/1 sort by title, by default on\n"
;

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  //process args
  AliZMQhistViewer viewer;
  int noptions = viewer.ProcessOptionString(argc,argv);
  if (noptions<=0) 
  {
    printf("%s",fUSAGE);
    return 1;
  }

  TApplication* gApp = new TApplication("viewer", &argc, argv); 
  gApp->SetReturnFromRun(true);
  //gApp->Run();

  gSystem->ResetSignal(kSigPipe);
  gSystem->ResetSignal(kSigQuit);
  gSystem->ResetSignal(kSigInterrupt);
  gSystem->ResetSignal(kSigTermination);

  viewer.Run(NULL);

  //destroy ZMQ sockets
  zmq_ctx_term(alizmq_context());
  return 0;
}

