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
// this program performs local monitoring on a GDC by running the HLT code   //
//                                                                           //
// If an argument is given, this is taken as the name of a date file which   //
// is used instead of the local node.                                        //
// The program can be stopped by pressing CTRL-C.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TError.h>
#include <TSysEvtHandler.h>
#ifdef DATE_SYS
#include <TROOT.h>
#include <TSystem.h>
#include <TDatime.h>
#include "AliRawReaderDate.h"
#include "event.h"
#include "monitor.h"
#include <AliLevel3.h>
#include <AliL3Transform.h>
#include <AliL3MemHandler.h>
#include <AliL3TrackArray.h>
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
  TSignalHandler(kSigInterrupt, kFALSE) 
{
  fStop = kFALSE;
};


//_____________________________________________________________________________
#ifdef DATE_SYS
int main(int argc, char** argv)
{
  // set ROOT in batch mode
  gROOT->SetBatch();   

  // open a log file
  FILE* file = fopen("monitorGDC.log", "w");
  TDatime time;

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
  status = monitorDeclareMp("GDC physics monitoring");
  if (status) ::Fatal("monitorDeclareMp", monitorDecodeError(status));

  // initialize HLT transformations
  if (!AliL3Transform::Init("./", kFALSE)) {
    ::Fatal("AliL3Transform::Init", "HLT initialization failed");
  }

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

    // read run, event, detector, DDL numbers and data size
    AliRawReaderDate rawReader(ptr);
    time.Set();
    printf("\n%s\n", time.AsString());
    printf("run: %d  event: %d %d\n", rawReader.GetRunNumber(), 
	   rawReader.GetEventId()[0], rawReader.GetEventId()[1]);
    while (rawReader.ReadHeader()) {
      printf(" detector: %d   DDL: %d  size: %d\n", 
	     rawReader.GetDetectorID(), rawReader.GetDDLID(), 
	     rawReader.GetDataSize());
    }

    if ((rawReader.GetAttributes()[0] & 0x02) != 0) {

      Int_t errorCode = rawReader.CheckData();
      if (errorCode && (errorCode != AliRawReader::kErrSize)) {
	time.Set();
	if (file) fprintf(file, "%s\n", time.AsString());
	if (file) fprintf(file, "run: %d  event: %d %d\n", 
			  rawReader.GetRunNumber(), 
			  rawReader.GetEventId()[0], 
			  rawReader.GetEventId()[1]);
	fprintf(file, "ERROR: %d\n\n", errorCode);

      } else {

	// run the HLT tracker
	AliLevel3* hlt = new AliLevel3((Char_t*)ptr);
	hlt->Init("./", AliLevel3::kDate, 1);

	hlt->SetClusterFinderParam(-1, -1, kTRUE);
  
	Int_t phiSegments = 50;
	Int_t etaSegments = 100;
	Int_t trackletlength = 3;
	Int_t tracklength = 20;//40 or 5
	Int_t rowscopetracklet = 2;
	Int_t rowscopetrack = 10;
	Double_t minPtFit = 0;
	Double_t maxangle = 0.1745;
	Double_t goodDist = 5;
	Double_t maxphi = 0.1;
	Double_t maxeta = 0.1;
	Double_t hitChi2Cut = 15;//100 or 15
	Double_t goodHitChi2 = 5;//20 or 5
	Double_t trackChi2Cut = 10;//50 or 10
	hlt->SetTrackerParam(phiSegments, etaSegments, 
			     trackletlength, tracklength,
			     rowscopetracklet, rowscopetrack,
			     minPtFit, maxangle, goodDist, hitChi2Cut,
			     goodHitChi2, trackChi2Cut, 50, maxphi, maxeta, 
			     kTRUE);
  
	gSystem->Exec("rm -rf hlt");
	gSystem->MakeDirectory("hlt");
	hlt->WriteFiles("./hlt/");
	hlt->ProcessEvent(0, 35, 0);

	time.Set();
	if (file) fprintf(file, "%s\n", time.AsString());
	if (file) fprintf(file, "run: %d  event: %d %d\n", 
			  rawReader.GetRunNumber(), 
			  rawReader.GetEventId()[0], 
			  rawReader.GetEventId()[1]);
	if (errorCode) fprintf(file, "ERROR: %d\n", errorCode);

	AliL3MemHandler memHandler;
	if (!memHandler.SetBinaryInput("hlt/tracks_0.raw")) {
	  if (file) fprintf(file, "no HLT tracks\n");
	  continue;
	}
	AliL3TrackArray* tracks = new AliL3TrackArray;
	memHandler.Binary2TrackArray(tracks);
	if (file) fprintf(file, "HLT found %d tracks\n", tracks->GetNTracks());
	delete tracks;
	memHandler.CloseBinaryInput();

	hlt->DoBench("hlt");
	if (file) {
	  FILE* bench = fopen("hlt.dat", "r");
	  while (bench && !feof(bench)) {
	    char buffer[256];
	    if (!fgets(buffer, 256, bench)) break;
	    fprintf(file, "%s", buffer);
	  }
	  fclose(bench);
	  fprintf(file, "\n");
	}

	gSystem->Exec("rm -rf hlt");
	delete hlt;
      }
    }

    gSystem->Sleep(100);   // sleep for 0.1 second
    free(ptr);

    gSystem->ProcessEvents();
    if (file) fflush(file);
  }

  gSystem->RemoveSignalHandler(handler);
  if (file) fclose(file);

  return 0;
}

#else
int main(int /*argc*/, char** /*argv*/)
{
  ::Fatal("main", "this program was compiled without DATE");

  return 1;
}
#endif
