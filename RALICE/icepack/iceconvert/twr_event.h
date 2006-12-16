#ifndef TWR_EVENT_H
#define TWR_EVENT_H

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

/////////////////////////////////////////////////////////////////////////
// This header file originates from the file event.h
// of Wolfgang Wagner's TWR raw data reader.
// The types u_int32_t etc... have been replaced by the ROOT portable
// types UInt_t etc... to obtain identical performance on all ROOT
// supported platforms from within the Ralice/IcePack framework.
// In addition the low level C "read()" invokations have been replaced
// by the standard C/C++ "fread()" invokations to even more enhance
// portability. As a result several header files have become obsolete
// and as such are not included anymore.
//
// NvE 07-dec-2006 Utrecht University
/////////////////////////////////////////////////////////////////////////

#include "twr_gps.h"

typedef struct
{
    UShort_t value[2500];
} waveform_t;

typedef struct
{
    UInt_t timestamp;
    UInt_t n_wfm;
} twr_t;

typedef struct
{
    UInt_t eventcounter;
    UInt_t which_trigger;
    GPS_t gps;
    waveform_t wfm[N_OF_CHANNELS];
    short wfm_filled[N_OF_CHANNELS];
    UInt_t twr_id_of_om[N_OF_CHANNELS];
    twr_t twr[ (MAX_N_CRATES+1) * MAX_N_TWR_PER_CRATE];
} event_t;

#endif
