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

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data from a date monitoring libraries.
/// It supports two modes - event taken from shared memory via DATE monitoring
/// libs, or an emulation mode when the events are taken from a DATE file using
/// the same monitoring libs.
/// The constructor requires an argument:
///
/// : - events are taken from shared memory
///  or
/// <DATE_filename> - events are taken from date file
///
/// Cvetan Cheshkov 1/04/2008
///////////////////////////////////////////////////////////////////////////////
#include <TSystem.h>
#include <cstdlib>
#include "AliRawReaderDateOnline.h"
#include "AliLog.h"
#ifdef ALI_DATE
#include "event.h"
#include "monitor.h"
#endif

ClassImp(AliRawReaderDateOnline)

AliRawReaderDateOnline::AliRawReaderDateOnline(
#ifdef ALI_DATE
				   const char* filename
#else
				   const char* /* filename */
#endif
				   ) :
  AliRawReaderDate((void*)NULL),
  fStop(kFALSE)
{

// Constructor
// Initialize the DATE monitoring libs

#ifdef ALI_DATE


  //  Removal of the selection of physics events
  //  Requested by Filimon and FMD experts
  //  fSelectEventType = PHYSICS_EVENT;

  int status;

  /* define data source : this is argument 1 */  
  status=monitorSetDataSource( (char* )filename );
  if (status!=0) {
    AliFatal(Form("monitorSetDataSource() failed : %s",monitorDecodeError(status)));
  }

  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    AliFatal(Form("monitorDeclareMp() failed : %s",monitorDecodeError(status)));
  }

  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);

  const Char_t* table[]  = {"ALL", "few", "*", "*",
                            "EOR", "yes","*", "*",
                            NULL, NULL, NULL, NULL};
  monitorDeclareTableExtended(const_cast<char**>(table));

  // install SIGUSR1 handler to allow clean end-of-events loop
  gSystem->AddSignalHandler(new AliRawReaderDateIntHandler(this));

#else
  Fatal("AliRawReaderDateOnline", "this class was compiled without DATE");
#endif
}

Bool_t AliRawReaderDateOnline::NextEvent()
{
// wait and get the next event
// from shared memory

#ifdef ALI_DATE

  // Stop on SIGUSR1
  if (fStop) {
    AliInfo("Raw-data reading stopped by SIGUSR1");
    if (fEvent) free(fEvent);
    fEvent = NULL;
    return kFALSE;
  }

  // Event already loaded no need take a new one
  if (AliRawReaderDate::NextEvent()) return kTRUE;

  if (fEvent) free(fEvent);
  fEvent = NULL;

  while (1) {
    /* get next event (blocking call until timeout) */
    int status=monitorGetEventDynamic((void**)&fEvent);

    if (status==MON_ERR_EOF) {
      AliInfo("End of File detected");
      Reset();
      fEvent = NULL;
      return kFALSE; /* end of monitoring file has been reached */
    }
    
    if (status!=0) {
      AliError(Form("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status)));
      Reset();
      fEvent = NULL;
      return kFALSE;
    }
    
    /* retry if got no event */
    if (fEvent==NULL) {
      continue;
    }
    
    eventTypeType eventT=fEvent->eventType;
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      AliInfo("EOR event detected");
      Reset();
      free(fEvent);
      fEvent = NULL;
      return kFALSE;
    }
    
    if (!IsEventSelected()) {
      free(fEvent);
      fEvent = NULL;
      continue;
    }

    AliInfo(Form("Run #%lu, event size: %lu, BC:0x%x, Orbit:0x%x, Period:0x%x",
		 (unsigned long)fEvent->eventRunNb,
		 (unsigned long)fEvent->eventSize,
		 EVENT_ID_GET_BUNCH_CROSSING(fEvent->eventId),
		 EVENT_ID_GET_ORBIT(fEvent->eventId),
		 EVENT_ID_GET_PERIOD(fEvent->eventId)
		 ));
    break;
  }

  fEventNumber++;
  Reset();

  return kTRUE;

}

#else
  return kFALSE;
}
#endif

AliRawReaderDateOnline::~AliRawReaderDateOnline()
{
// Destructor
// Free the last event in shared memory

#ifdef ALI_DATE
  if (fEvent) free(fEvent);
#endif
}

void AliRawReaderDateOnline::SelectEvents(Int_t type,
  ULong64_t triggerMask,
  const char *triggerExpr)
{
  // Select event by using DATE monitoring
  // library
#ifdef ALI_DATE
  const Char_t* table[]  = {"ALL", "no", "*", "*",
			    "PHY", "yes","*", "*",
                            "EOR", "yes","*", "*",
 			    NULL, NULL, NULL, NULL};
  TString trSelection;
  for (Int_t i = 0; i < 50; i++) {
    if (triggerMask & (1ull << i)) {
	if (!trSelection.IsNull()) trSelection += "&";
	trSelection += Form("%d",i);
    }
  }
  table[7] = trSelection.Data();

  monitorLogout();
  monitorDeclareTableExtended(const_cast<char**>(table));
  
#endif
  AliRawReader::SelectEvents(type,triggerMask,triggerExpr);
}

//______________________________________________________________________________
void AliRawReaderDateOnline::Stop()
{
  // Stop the event loop (called on SIGUSR1)

  fStop = kTRUE; 
}
