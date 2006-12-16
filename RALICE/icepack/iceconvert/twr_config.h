#ifndef TWR_CONFIG_H
#define TWR_CONFIG_H

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

/////////////////////////////////////////////////////////////////////////
// This header file originates from the file config.h
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

#include "twr_globalconst.h"

#define LENGTH_OF_STRING 100
#define MAX_N_VLF_VETOS_PER_HOUR 16

/*
  Enhanced structures for Analysis
*/

typedef struct  /* of one module */
{
  UInt_t base;   /* Baseaddress in VME-BUS */
  UInt_t id;     /* identifier for TWR  */
  /* id = crate no * 0x10      (begining with 1) */
  /*     + TWR number in crate (beginnig with 0)  */
  
  UInt_t mod_id;    /* Module ID from internal register   0x4 */
  UInt_t ext_start; /* Extern Start Delay Register        0x14 */
  UInt_t ext_stop;  /* Extern Stop Delay Register         0x18 */
  
  UInt_t om_no[CHANNELS_PER_TWR]; 
  UInt_t om_is_optical[CHANNELS_PER_TWR];
  UInt_t baseline[CHANNELS_PER_TWR];
  UInt_t threshold[CHANNELS_PER_TWR];
} twr_config_t;


typedef struct
{
  UInt_t vme_base_bridge;   /* physical base local bridge in master crate */
  UInt_t vme_base_100MHz;   /* physical base local 100MHz module */
 
  char dsp_program_name[LENGTH_OF_STRING];
  UInt_t base_gps; /* physical base of GPSAMD in remote crate */
  UInt_t n_twr;  
  twr_config_t *twr[MAX_N_TWR_PER_CRATE];
} crate_config_t;

typedef struct
{
  Int_t pci_channel;
  char pci_driver_name[LENGTH_OF_STRING];
  char data_output_directory[LENGTH_OF_STRING];
  UInt_t n_events_per_file;
  UInt_t compression;

  /* DSP parameters */
  UInt_t preceeding_values_optical;
  UInt_t following_values_optical;
  UInt_t preceeding_values_electrical;
  UInt_t following_values_electrical;

  /* TWR */
  UInt_t clockdiv;  /* Clockpredivider for Time Stamp for TWRs  */ 
  UInt_t statreg;   /* Status Register                    0x0 */  
  UInt_t acq_ctrl;  /* Acquisition Control Register       0x10 */
  UInt_t evtconfig;
  UInt_t triggerflagclearcounter;

  /* GPS AMD / VLF stuff */
  UInt_t b_vlf_trigger_veto_on; /* Enable VLF trigger veto */
  UInt_t n_vlf_trigger_veto; /* no of trigger times per hour */
  UInt_t start_vlf_trigger_veto[MAX_N_VLF_VETOS_PER_HOUR];
  UInt_t duration_vlf_trigger_veto;
  UInt_t overlap_vlf_trigger_veto;

  UInt_t n_crates;
  UInt_t n_twr;

  crate_config_t *crate[MAX_N_CRATES+1];

  twr_config_t *twr_field[MAX_N_CRATES * MAX_N_TWR_PER_CRATE];

}sys_config_t;



#endif














