/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// this program performs local monitoring on a GDC and checks the            //
// consistency of the data                                                   //
//                                                                           //
// If an argument is given, this is taken as the name of a date file which   //
// is used instead of the local node.                                        //
// The program can be stopped by pressing CTRL-C.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TError.h>
#include <TSysEvtHandler.h>
#include <cstdlib>
#ifdef ALI_DATE
#include <TROOT.h>
#include <TSystem.h>
#include <TDatime.h>
#include "AliRawReaderDate.h"
#include "event.h"
#include "monitor.h"
#endif

//_____________________________________________________________________________
class AliGDCInterruptHandler : public TSignalHandler {
public:
  AliGDCInterruptHandler();
  Bool_t Notify() {fStop = kTRUE; return kTRUE;};
  Bool_t Stop() const {return fStop;};
private:
  Bool_t fStop;  // CTRL-C pressed
};

//_____________________________________________________________________________
AliGDCInterruptHandler::AliGDCInterruptHandler() : 
  TSignalHandler(kSigInterrupt, kFALSE),
  fStop(kFALSE)  
{

}


//_____________________________________________________________________________
#ifdef ALI_DATE
int main(int argc, char** argv)
{
  // set ROOT in batch mode
  gROOT->SetBatch();   

  // get data from a file or online from this node
  Int_t status = 0;
  if (argc > 1) {
    status = monitorSetDataSource(argv[1]);
  } else {
    status = monitorSetDataSource(":");
  }
  if (status) ::Fatal("monitorSetDataSource", monitorDecodeError(status));

  // monitor only a sample of physics events
  char* table[] = {"Physics event", "yes", NULL};
  status = monitorDeclareTable(table);
  if (status) ::Fatal("monitorDeclareTable", monitorDecodeError(status));

  // declare this monitoring program to DATE
  status = monitorDeclareMp("data consistency check");
  if (status) ::Fatal("monitorDeclareMp", monitorDecodeError(status));

  // create the signal handler
  AliGDCInterruptHandler* handler = new AliGDCInterruptHandler;
  gSystem->AddSignalHandler(handler);

  // endless loop
  void* ptr = NULL;
  while (!handler->Stop()) {
    // get the next event
    status = monitorGetEventDynamic(&ptr);
    if (status == (Int_t)MON_ERR_EOF) break;
    if (status) ::Fatal("monitorGetEventDynamic", monitorDecodeError(status));

    // if no new event
    if (!ptr) {
      gSystem->Sleep(1000);   // sleep for 1 second
      continue;
    }

    // check the data
    AliRawReaderDate rawReader(ptr);
    Int_t errorCode = rawReader.CheckData();
    if ((errorCode != 0) && (errorCode != AliRawReader::kErrSize)) {
      TDatime time;
      printf("\n%s\n", time.AsString());
      const char* errMsg[6] = {"no error", "wrong magic word in sub event",
			       "no mini header", 
			       "wrong magic word in mini header",
			       "inconsistent size in sub event and mini header",
			       "access to data out of bounds"};
      printf("Error: %s\n", errMsg[errorCode]);
      printf("run: %d  event: %d %d\n", rawReader.GetRunNumber(), 
	     rawReader.GetEventId()[0], rawReader.GetEventId()[1]);
      printf("trigger: %08x %08x   detector: %08x\n",
	     rawReader.GetTriggerPattern()[0], 
	     rawReader.GetTriggerPattern()[1],
	     rawReader.GetDetectorPattern()[0]);
      printf("attributes: %08x %08x %08x\n", rawReader.GetAttributes()[0], 
	     rawReader.GetAttributes()[1], rawReader.GetAttributes()[2]);
    }

    free(ptr);
    gSystem->ProcessEvents();
  }

  gSystem->RemoveSignalHandler(handler);

  return 0;
}

#else
int main(int /*argc*/, char** /*argv*/)
{
  ::Fatal("main", "this program was compiled without DATE");

  return 1;
}
#endif
