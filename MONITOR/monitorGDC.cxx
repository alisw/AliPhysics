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
#include <AliL3StandardIncludes.h>
#include <AliL3Transform.h>
#include <AliL3MemHandler.h>
#include <AliL3TrackArray.h>
#include <AliL3HoughMaxFinder.h>
#include <AliL3HoughBaseTransformer.h>
#include <AliL3Hough.h>
#include <AliL3Benchmark.h>
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

    AliRawReaderDate rawReader(ptr);
    //    if ((rawReader.GetAttributes()[0] & 0x02) != 0) {

    //      Int_t errorCode = rawReader.CheckData();
    Int_t errorCode = 0;
      if (errorCode && (errorCode != AliRawReader::kErrSize)) {
	time.Set();
	if (file) fprintf(file, "%s\n", time.AsString());
	if (file) fprintf(file, "run: %d  event: %d %d\n", 
			  rawReader.GetRunNumber(), 
			  rawReader.GetEventId()[0], 
			  rawReader.GetEventId()[1]);
	fprintf(file, "ERROR: %d\n\n", errorCode);

      } else {

	AliL3Benchmark *fBenchmark = new AliL3Benchmark();
	fBenchmark->Start("Overall timing");

	Float_t ptmin = 0.1*AliL3Transform::GetSolenoidField();

	// Run the Hough Transformer
	fBenchmark->Start("Init");
	AliL3Hough *hough1 = new AliL3Hough();

	hough1->SetThreshold(4);
	hough1->CalcTransformerParams(ptmin);
	hough1->SetPeakThreshold(70,-1);
	//	printf("Pointer is %x\n",ptr);
	hough1->Init("./", kFALSE, 100, kFALSE,4,0,(Char_t*)ptr,0.0);
	hough1->SetAddHistograms();
	fBenchmark->Stop("Init");

	fBenchmark->Start("Init");
	AliL3Hough *hough2 = new AliL3Hough();

	hough2->SetThreshold(4);
	hough2->CalcTransformerParams(ptmin);
	hough2->SetPeakThreshold(70,-1);
	//	printf("Pointer is %x\n",ptr);
	hough2->Init("./", kFALSE, 100, kFALSE,4,0,(Char_t*)ptr,0.0);
	hough2->SetAddHistograms();
	fBenchmark->Stop("Init");

	Int_t n_global_tracks = 0;
	hough1->StartProcessInThread(0,17);
	hough2->StartProcessInThread(18,35);

	//	gSystem->Sleep(20000);
	if(hough1->WaitForThreadFinish())
	  ::Fatal("AliL3Hough::WaitForThreadFinish"," Can not join the required thread! ");
	if(hough2->WaitForThreadFinish())
	  ::Fatal("AliL3Hough::WaitForThreadFinish"," Can not join the required thread! ");

	gSystem->MakeDirectory("hough1");
	hough1->WriteTracks("./hough1");
	gSystem->MakeDirectory("hough2");
	hough2->WriteTracks("./hough2");

	/*
	for(int slice=0; slice<=35; slice++)
	{
	  //	  cout<<"Processing slice "<<slice<<endl;
	  fBenchmark->Start("ReadData");
	  hough->ReadData(slice,0);
	  fBenchmark->Stop("ReadData");
	  fBenchmark->Start("Transform");
	  hough->Transform();
	  fBenchmark->Stop("Transform");
	  hough->AddAllHistogramsRows();
	  hough->FindTrackCandidatesRow();
	  fBenchmark->Start("AddTracks");
	  hough->AddTracks();
	  fBenchmark->Stop("AddTracks");

	  AliL3TrackArray* tracks = (AliL3TrackArray*)hough->GetTracks(0);
	  n_global_tracks += tracks->GetNTracks();
	}
	*/
	fBenchmark->Stop("Overall timing");
	time.Set();
	if (file) fprintf(file, "%s\n", time.AsString());
	if (file) fprintf(file, "run: %d  event: %d %d\n", 
			  rawReader.GetRunNumber(), 
			  rawReader.GetEventId()[0], 
			  rawReader.GetEventId()[1]);
	if (errorCode) fprintf(file, "ERROR: %d\n", errorCode);

	if (file) fprintf(file, "Hough Transformer found %d tracks\n",n_global_tracks);

	hough1->DoBench("hough1");
	hough2->DoBench("hough2");
	fBenchmark->Analyze("overall");
	if (file) {
	  FILE* bench = fopen("hough1.dat", "r");
	  while (bench && !feof(bench)) {
	    char buffer[256];
	    if (!fgets(buffer, 256, bench)) break;
	    fprintf(file, "%s", buffer);
	  }
	  fclose(bench);
	}
	if (file) {
	  FILE* bench = fopen("hough2.dat", "r");
	  while (bench && !feof(bench)) {
	    char buffer[256];
	    if (!fgets(buffer, 256, bench)) break;
	    fprintf(file, "%s", buffer);
	  }
	  fclose(bench);
	}
	if (file) {
	  FILE* bench = fopen("overall.dat", "r");
	  while (bench && !feof(bench)) {
	    char buffer[256];
	    if (!fgets(buffer, 256, bench)) break;
	    fprintf(file, "%s", buffer);
	  }
	  fclose(bench);
	  fprintf(file, "\n\n");
	}

	delete hough1;
	delete hough2;
      }
    //    }

    /*
    // read run, event, detector, DDL numbers and data size
    AliRawReaderDate rawReader(ptr);
    time.Set();
    printf("\n%s\n", time.AsString());
    printf("run: %d  event: %d %d\n", rawReader.GetRunNumber(), 
	   rawReader.GetEventId()[0], rawReader.GetEventId()[1]);
    while (rawReader.ReadMiniHeader()) {
      printf(" detector: %d   DDL: %d  size: %d\n", 
	     rawReader.GetDetectorID(), rawReader.GetDDLID(), 
	     rawReader.GetDataSize());
    }

    */

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
