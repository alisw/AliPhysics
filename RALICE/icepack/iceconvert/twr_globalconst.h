#ifndef TWR_GLOBALCONST_H
#define TWR_GLOBALCONST_H

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

/////////////////////////////////////////////////////////////////////////
// This header file originates from the file globalconst.h
// of Wolfgang Wagner's TWR raw data reader.
//
// NvE 07-dec-2006 Utrecht University
/////////////////////////////////////////////////////////////////////////

/*
  TWR system
*/

#define VALUES_PER_WAVEFORM 1024
#define CHANNELS_PER_TWR 8
#define N_OF_EVENTS_PER_BANK 128
#define MEMORY_BANKS_PER_MODULE 2

#define MAX_N_CRATES 7
#define MAX_N_TWR_PER_CRATE 16

#define MAX_N_OF_FRAGS 20
#define MAX_N_OF_TDC_EDGES 16

#define N_OF_CHANNELS 700
#define N_OF_TRIGGERS 12

#define BASELINE_MEAN_MAGIC 3100

/*
  TWR Calibration
*/
#define NSECS_PER_TDC_BIN 1.041667 
#define NSECS_PER_TWR_BIN 10 
#define MV_PER_DIGIT (5000.0/4096.0)
#define BINxDIGIT2CHARGE ( 10.0 * 1e-9 *(5.0/4096.0) / 50.0 )
/*
  TWR analysis
*/

#define MAX_N_OF_PEAKS 300

#define THRESHOLDTAIL_PERCENTAGE 0.1
#define THRESHOLD_2ND_PEAK 10
#define THRESHOLD_2ND_PEAK_PERCENTAGE 0.5
#define THRESHOLD_OFFSET 20

#define MAX_SLOPE_LENGTH_ELEC 7

/*
  TWR Trigger stuff
*/

/* 2004 */

/* > 2005 */
#define RANDOM_TRIGGER_2005    690
#define MAIN_TRIG_OR_2005      689
#define MAIN_TRIG_LOGIC_2005   688
#define M12TRIGGER_2005        687
#define CAL_TRIG_LA_2005       686
#define CAL_TRIG_T0_2005       685
#define SPASE_TRI_2005         684
#define STRING_TRIG_2005       683
#define M18TRIGGER_2005        682
#define M24TRIGGER_2005        681

#endif
