#ifndef TWR_READER_H
#define TWR_READER_H

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

/////////////////////////////////////////////////////////////////////////
// This header file originates from the files read_twr_binary_file.h
// and lib_wf2hit_new.h of Wolfgang Wagner's TWR raw data reader.
// The types u_int32_t etc... have been replaced by the ROOT portable
// types UInt_t etc... to obtain identical performance on all ROOT
// supported platforms from within the Ralice/IcePack framework.
// In addition the low level C "read()" invokations have been replaced
// by the standard C/C++ "fread()" invokations to even more enhance
// portability. As a result several header files have become obsolete
// and as such are not included anymore.
// The source code of read_twr_binary_file.cxx and lib_wf2hit_new.cxx
// has been implemented as memberfunctions of the class IceRawTWR,
// so the definitions of the function prototypes have been moved to
// the header file IceRawTWR.h.
//
// NvE 07-dec-2006 Utrecht University
/////////////////////////////////////////////////////////////////////////

#include "twr_globalconst.h"
#include "twr_config.h"
#include "twr_gps.h"
#include "twr_event.h"

/////////////////////////////////////
// Source of raw_twr_binary_file.h //
/////////////////////////////////////

/* Readout 2005 / 2006 */
/* offset for baseline */
// #define BASELINE_MEAN_MAGIC 3100 // Already defined in globalconst.h

/* Trigger defs */
#define M24_TRIG_CHA          681
#define M18_TRIG_CHA          682
#define STRING_TRIG_CHA       683
#define SPASE_TRIG_CHA        684
#define CAL_TRIG_T0_CHA       685
#define CAL_TRIG_LA_CHA       686
#define M12_TRIG_CHA          687
#define MAIN_TRIG_LOGIC_CHA   688
#define MAIN_TRIG_OR_CHA      689
#define RANDOM_TRIG_CHA       690

#define M20_FRAG_BIT       0x1       /* only software no hardware signal */
#define VOLUMEN_TRIG_BIT   0x2       /* only software no hardware signal */
#define M18_TRIG_BIT       0x4
#define M24_TRIG_BIT       0x8
#define STRING_TRIG_BIT    0x10
#define SPASE_TRIG_BIT     0x20
#define RANDOM_TRIG_BIT    0x40
#define CAL_TRIG_T0_BIT    0x80
#define CAL_TRIG_LA_BIT    0x100

#define MUON_VETO_MASK     0xffff0000

#define M24_TRIG_NUM         0
#define M18_TRIG_NUM         1
#define STRING_TRIG_NUM      2
#define SPASE_TRIG_NUM       3
#define CAL_TRIG_T0_NUM      4
#define CAL_TRIG_LA_NUM      5
#define M12_TRIG_NUM         6
#define MAIN_TRIG_LOGIC_NUM  7
#define MAIN_TRIG_OR_NUM     8
#define RANDOM_TRIG_NUM      9 
#define M20_FRAG_NUM        10 /* only software no hardware signal */
#define VOLUMEN_TRIG_NUM    11 /* only software no hardware signal */

// Software triggers			    
const UInt_t trigger_bits[N_OF_TRIGGERS]
= { M24_TRIG_BIT, M18_TRIG_BIT, STRING_TRIG_BIT, SPASE_TRIG_BIT,
    CAL_TRIG_T0_BIT, CAL_TRIG_LA_BIT, 0x0, 0x0, 
    0x0, RANDOM_TRIG_BIT, M20_FRAG_BIT, VOLUMEN_TRIG_BIT };

// Hardware triggers
const UInt_t trigger_channel[N_OF_TRIGGERS]
= { M24_TRIG_CHA, M18_TRIG_CHA, STRING_TRIG_CHA, SPASE_TRIG_CHA,
    CAL_TRIG_T0_CHA, CAL_TRIG_LA_CHA, M12_TRIG_CHA, MAIN_TRIG_LOGIC_CHA, 
    MAIN_TRIG_OR_CHA, RANDOM_TRIG_CHA, 0x0, 0x0 };

typedef struct{
  Int_t n_software_trigger;
  Int_t n_hardware_trigger;
  Int_t first_trigger;
  UInt_t first_trigger_time;
  UInt_t trigger_active[N_OF_TRIGGERS];
  UInt_t trigger_has_pulse[N_OF_TRIGGERS];
  UInt_t trigger_time[N_OF_TRIGGERS]; 
} trigger_hits_t;

/* Error codes */
#define ERROR_NOT_VALID_FILENAME 1
#define ERROR_TOO_MANY_CRATES 2
#define ERROR_TOO_MANY_TWRS   1

/* info from TWR filename */
typedef struct{
  Int_t day;
  Int_t year;
  Int_t run_no;
  Int_t file_no; 
  Int_t begin; 
  Int_t end;
} twr_raw_data_file_t;

////////////////////////////////
// Source of lib_wf2hit_new.h //
////////////////////////////////

typedef struct
{
  // output from analysis
  Int_t n_frag;

  Int_t frag_n_points[MAX_N_OF_FRAGS];
  Int_t frag_begin[MAX_N_OF_FRAGS];
  Int_t frag_end[MAX_N_OF_FRAGS];
  Float_t frag_begin_time[MAX_N_OF_FRAGS]; // in TDC frame
  //Float_t frag_t_end[MAX_N_OF_FRAGS];
  Float_t frag_mean[MAX_N_OF_FRAGS];


  Int_t n_peak;
  Int_t peak_begin[MAX_N_OF_PEAKS];  
  Int_t peak_end[MAX_N_OF_PEAKS];
  Int_t peak_max[MAX_N_OF_PEAKS];
  Int_t peak_local_minimum[MAX_N_OF_PEAKS];
  Int_t peak_in_fragment[MAX_N_OF_PEAKS];
  Int_t peak_TDC_edge[MAX_N_OF_PEAKS]; // holds the identified TDC edge
  Int_t peak_isolated[MAX_N_OF_PEAKS];
  Int_t crosstalk_charge_n_value[MAX_N_OF_PEAKS];
  
  Float_t peak_mean[MAX_N_OF_PEAKS];
  Float_t peak_m[MAX_N_OF_PEAKS];
  Float_t peak_b[MAX_N_OF_PEAKS];
  Float_t peak_t0[MAX_N_OF_PEAKS];
  Float_t peak_begin_time[MAX_N_OF_PEAKS];   
  Float_t peak_charge[MAX_N_OF_PEAKS];
  Float_t peak_height[MAX_N_OF_PEAKS];
  Float_t fitted_amplitude[MAX_N_OF_PEAKS];
  Float_t fitted_TOT[MAX_N_OF_PEAKS];

  Float_t crosstalk_charge[MAX_N_OF_PEAKS];
  Float_t crosstalk_slope[MAX_N_OF_PEAKS];
  Int_t n_point;
  Float_t wfm_x[1024];
  Float_t wfm_y[1024];
  Float_t wfm_min; 
  Float_t wfm_max; 
  Int_t b_out_of_range;

  // TDC data
  Int_t n_tdc_edges;
  Int_t identified_twr_hit[MAX_N_OF_TDC_EDGES];
  Double_t leading_edge[MAX_N_OF_TDC_EDGES];
  Double_t falling_edge[MAX_N_OF_TDC_EDGES];

} waveform_analyse_t;

#define ERROR_CALIBRATING_WFM      1

#define ERROR_FILL_OUT_OF_RANGE    2
#define ERROR_FILL_MAX_N_OF_FRAGS  3

#define ERROR_MISS_FRAG_STOP       4
#define ERROR_CORRUPTED_WF         5
#define ERROR_MAX_N_FRAGMENTS_EXCEEDED 6
		  		  
#endif
