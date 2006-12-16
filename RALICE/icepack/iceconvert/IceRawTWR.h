#ifndef IceRawTWR_h
#define IceRawTWR_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"

#include "AliJob.h"
#include "AliObjMatrix.h"

#include "IceAOM.h"
#include "IceEvent.h"

#include "twr_reader.h"

class IceRawTWR : public AliJob
{
 public :
  IceRawTWR(const char* name="IceRawTWR",const char* title=""); // Constructor
  virtual ~IceRawTWR();                                         // Destructor
  void SetMaxEvents(Int_t n);       // Set maximum number of events to be processed
  void SetPrintFreq(Int_t f);       // Set printfrequency to provide info every f events
  void SetSplitLevel(Int_t split);  // Set split level for the produced ROOT data file
  void SetBufferSize(Int_t bsize);  // Set buffersize for the produced ROO data file
  void AddInputFile(TString name);  // Add name of F2K input file to the list
  void SetOutputFile(TFile* ofile); // Set output file for the ROOT data structures       
  void SetOutputFile(TString name); // Create output file for the ROOT data structures
  TFile* GetOutputFile();           // Provide pointer to the ROOT output file
  virtual void Exec(Option_t* opt); // Perform the format conversion

 protected :
  Int_t fSplit;        // The split level of the produced ROOT data file
  Int_t fBsize;        // The buffersize of the produced ROOT data file
  Int_t fMaxevt;       // The maximum number of events to be processed
  Int_t fPrintfreq;    // The event info printing frequency
  TObjArray* fInfiles; // Names of all the raw data input files
  TFile* fOutfile;     // The ROOT output file
  void PutTrigger(Int_t year);   // Put the trigger info from the raw data event into the IcePack structure
  void PutWaveforms(Int_t year); // Put the waveforms from the raw data event into the IcePack structure

  FILE* fInput;              //! Pointer to the TWR raw data input file
  sys_config_t* fHeader;     //! Structure holding the raw configuration header info
  event_t fEvent;            //! Structure holding the actual raw event data
  trigger_hits_t fTrigger;   //! Structure holding the event trigger info
  waveform_analyse_t fWform; //! Waveform info for a certain OM from (merged) fragment(s)

  // Ralice/IcePack implementation of Wolfgang Wagner's original code
  Int_t extract_info_from_filename(char* fname,twr_raw_data_file_t* twr_file);
  Int_t clear_system(sys_config_t* sys);
  Int_t clear_event(event_t* event_ptr);
  Int_t read_header_from_file(FILE* fin,sys_config_t** system_ptr,UInt_t* header_length);
  Int_t update_system(sys_config_t* sys,Int_t run_number);
  Int_t read_event(FILE* fin,sys_config_t* sys,event_t* event_ptr);
  Int_t retrigger(event_t* ev,trigger_hits_t* trig);
  Int_t clear_waveform_analysis(waveform_analyse_t* wfm_om);
  Int_t restore_waveform(waveform_t f_wfm,waveform_analyse_t* wfm_om,Int_t year);

 ClassDef(IceRawTWR,1) // Job for conversion of TWR raw data into IceEvent data structures.
};
#endif
