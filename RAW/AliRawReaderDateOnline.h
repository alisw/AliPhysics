#ifndef ALIRAWREADERDATEONLINE_H
#define ALIRAWREADERDATEONLINE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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

#include <TSysEvtHandler.h>

#include "AliRawReaderDate.h"

class AliRawReaderDateOnline: public AliRawReaderDate {
  public :
    AliRawReaderDateOnline(const char* fileName, const Char_t** customTable = NULL);
    virtual ~AliRawReaderDateOnline();

    virtual Bool_t   NextEvent();
    //    virtual Bool_t   RewindEvents();

    // Method which can be used in order to force the auto-save on
    // ESD tree inside AliReconstruction. For the moment it will be
    // activated only for AliRawReaderDateOnline.
    virtual Bool_t   UseAutoSaveESD() const { return kTRUE; }

    // Method triggered by signal hanlder
    // Set fStop to false in which case
    // NextEvent() returns fFALSE and the
    // processing of raw data stops
    virtual void     Stop();

  protected:
    class AliRawReaderDateIntHandler : public TSignalHandler {
    public:
    AliRawReaderDateIntHandler(AliRawReaderDateOnline *rawReader):
      TSignalHandler(kSigUser1, kFALSE), fRawReader(rawReader) { }
      Bool_t Notify() {
	Info("Notify", "received a SIGUSR1 signal");
	fRawReader->Stop();
	return kTRUE;
      }
    private:
      AliRawReaderDateOnline *fRawReader;   // raw-reader to signal

      AliRawReaderDateIntHandler(const AliRawReaderDateIntHandler& handler); // Not implemented
      AliRawReaderDateIntHandler& operator=(const AliRawReaderDateIntHandler& handler); // Not implemented
    };

    virtual void     SelectEvents(Int_t type, ULong64_t triggerMask = 0, const char *triggerExpr = NULL);

  private:
    AliRawReaderDateOnline(const AliRawReaderDateOnline& rawReader);
    AliRawReaderDateOnline& operator = (const AliRawReaderDateOnline& rawReader);

    const Char_t**   fTable;// custom monitoring table
    Bool_t           fStop; // raw-reader signaled to stop

    ClassDef(AliRawReaderDateOnline, 0) // class for reading DATE raw data from shared memory
};

#endif
