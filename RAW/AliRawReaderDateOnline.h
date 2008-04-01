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

#include "AliRawReaderDate.h"

class AliRawReaderDateOnline: public AliRawReaderDate {
  public :
    AliRawReaderDateOnline(const char* fileName);
    virtual ~AliRawReaderDateOnline();

    virtual Bool_t   NextEvent();
    virtual Bool_t   RewindEvents();

  private:
    AliRawReaderDateOnline(const AliRawReaderDateOnline& rawReader);
    AliRawReaderDateOnline& operator = (const AliRawReaderDateOnline& rawReader);

    ClassDef(AliRawReaderDateOnline, 0) // class for reading DATE raw data from shared memory
};

#endif
