#ifndef TWR_GPS_H
#define TWR_GPS_H

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

/////////////////////////////////////////////////////////////////////////
// This header file originates from the file gps.h
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

/* GPS Configuration Registers        */
/* for the DAQ300                     */
/* Wolfgang Wagner - 22.1.02          */
/* according to new GPS Docu 13.12.01 */

/* Timestructure...*/

struct timebits
{
  UInt_t seconds:8,
    year:8,
    TQC:4, /* True time quality character 0 best, F worst */
    status:4,
    dummy:8;
};

union time
{
  UInt_t all;
  struct timebits bits;
};

typedef struct                          // GPS-event => 16 bytes
{
  UInt_t count_10MHz;
  UInt_t seconds;
  union time info;
  UInt_t counter;
}GPS_t;

#define GPS_CTRL_REG     0x0 // Write

#define GPS_CTRL_RESET   0x1
#define GPS_CTRL_WORK    0x0
#define GPS_CTRL_YR_SRC  0x2
// #define GPS_CTRL_CLKREF  0x4 // Clock Reference from GPS-Clock
#define GPS_CTRL_ENPULS  0x8
#define GPS_FIFO_NOTEMPT 0x10
#define GPS_FIFO_FULL    0x20
//+++++++++++++++++++++++++++

#define GPS_FREQ_DIV     0x4 // Only Read

//+++++++++++++++++++++++++++

#define GPS_YEAR         0x8
//+++++++++++++++++++++++++++

#define GPS_LED_REG      0xc
//+++++++++++++++++++++++++++

#define GPS_DATA_REG     0x10

#endif
