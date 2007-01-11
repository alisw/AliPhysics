/*******************************************************************************
 * Copyright(c) 2003, IceCube Experiment at the South Pole. All rights reserved.
 *
 * Author: The IceCube RALICE-based Offline Project.
 * Contributors are mentioned in the code where appropriate.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation strictly for non-commercial purposes is hereby granted
 * without fee, provided that the above copyright notice appears in all
 * copies and that both the copyright notice and this permission notice
 * appear in the supporting documentation.
 * The authors make no claims about the suitability of this software for
 * any purpose. It is provided "as is" without express or implied warranty.
 *******************************************************************************/

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class IceRawTWR
// Conversion of Amanda raw TWR data into IceEvent data structures.
// The code to actually read the TWR raw data structures is an Ralice/IcePack
// implementation of Wolfgang Wagner's (Dortmund University, Germany)
// original read_twr_binary_file.cxx and wf2hit_new.cxx source code.
// The trigger information as encountered in the raw data, is available
// in the IceEvent structure via a device named "Trigger".
// The various triggers (and times) have been stored as different "hits"
// in this "Trigger" device, just like it was done in the IceF2k processor
// for the mu-daq F2K data.
// An indication of the active DAQ system is available in the IceEvent structure
// via a device named "Daq". Here the various daq systems (TWR, Muon, ...)
// from which the actual hits (ADC, LE, TOT) eventually will be composed
// are indicated as "signals" of the device itself. 
// This class is derived from AliJob providing a task-based processing
// structure on an event-by-event basis.
// The main object in the job environment is an IceEvent* pointer.
// In case the user has provided sub-tasks, these will be executed
// on an event-by-event basis after the IceEvent structure has been filled
// with the raw TWR data and before the final structures are written out.
// Note that the data structures are only written out if an outputfile has
// been specified via the SetOutputFile memberfunction.
// In case no outputfile has been specified, this class provides a facility
// to investigate/analyse raw TWR data using the Ralice/IcePack analysis tools.
//
// Usage example :
// ---------------
//
// gSystem->Load("ralice");
// gSystem->Load("icepack");
// gSystem->Load("iceconvert");
//
// IceRawTWR q("IceRawTWR","TWR raw data to IcePack data structure conversion");
//
// // Limit the number of entries for testing
// q.SetMaxEvents(10);
//
// // Print frequency to produce a short summary print every printfreq events
// q.SetPrintFreq(1);
//
// // The TWR raw data input filename(s)
// q.AddInputFile("twr_2005_101_009225_0983_57784_57850.dat.twr.to_tape_1");
//
// // Output file for the event structures
// q.SetOutputFile("events.root");
//
// ///////////////////////////////////////////////////////////////////
// // Here the user can specify his/her sub-tasks to be executed
// // on an event-by-event basis after the IceEvent structure
// // has been filled and before the data is written out.
// // Sub-tasks (i.e. a user classes derived from TTask) are entered
// // as follows :
// //
// //    MyXtalk task1("task1","Cross talk correction");
// //    MyClean task2("task2","Hit cleaning");
// //    q.Add(&task1);
// //    q.Add(&task2);
// //
// // The sub-tasks will be executed in the order as they are entered.
// ///////////////////////////////////////////////////////////////////
//
// // Perform the conversion and execute subtasks (if any)
// // on an event-by-event basis
// q.ExecuteJob();
//
//--- Author: Nick van Eijndhoven 12-dec-2006 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceRawTWR.h"
#include "Riostream.h"

ClassImp(IceRawTWR) // Class implementation to enable ROOT I/O

IceRawTWR::IceRawTWR(const char* name,const char* title) : AliJob(name,title)
{
// Default constructor.
// By default maxevent=-1, split=0, bsize=32000, printfreq=1.

 fSplit=0;
 fBsize=32000;
 fMaxevt=-1;
 fPrintfreq=1;
 fInfiles=0;
 fOutfile=0;
}
///////////////////////////////////////////////////////////////////////////
IceRawTWR::~IceRawTWR()
{
// Default destructor.

 if (fInfiles)
 {
  delete fInfiles;
  fInfiles=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::SetMaxEvents(Int_t n)
{
// Set the maximum number of events to be processed.
// n=-1 implies processing of the complete input file, which is the default
// initialisation in the constructor.
 fMaxevt=n;
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::SetPrintFreq(Int_t f)
{
// Set the printfrequency to produce info every f events.
// f=1 is the default initialisation in the constructor.
 if (f>=0) fPrintfreq=f;
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::SetSplitLevel(Int_t split)
{
// Set the split level for the ROOT data file.
// split=0 is the default initialisation in the constructor.
 if (split>=0) fSplit=split;
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::SetBufferSize(Int_t bsize)
{
// Set the buffer size for the ROOT data file.
// bsize=32000 is the default initialisation in the constructor.
 if (bsize>=0) fBsize=bsize;
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::AddInputFile(TString name)
{
// Add the name of this TWR raw data input file to the list to be processed.

 if (!fInfiles)
 {
  fInfiles=new TObjArray();
  fInfiles->SetOwner();
 }

 TObjString* s=new TObjString();
 s->SetString(name);
 fInfiles->Add(s);
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::SetOutputFile(TFile* ofile)
{
// Set the output file for the ROOT data.
 if (fOutfile) delete fOutfile;
 fOutfile=ofile;
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::SetOutputFile(TString name)
{
// Create the output file for the ROOT data.
 if (fOutfile) delete fOutfile;
 fOutfile=new TFile(name.Data(),"RECREATE","F2K data in IceEvent structure");
}
///////////////////////////////////////////////////////////////////////////
TFile* IceRawTWR::GetOutputFile()
{
// Provide pointer to the ROOT output file.
 return fOutfile;
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::Exec(Option_t* opt)
{
// Job to loop over the specified number of events and convert the 
// TWR raw data into the IceEvent structure.
// If maxevents<0 (default) all the entries of the input file
// will be processed.
// Every "printfreq" events a short event summary will be printed.
// The default value is printfreq=1.
// The output will be written on a standard output tree named "T".
//
// Notes :
// -------
// 1) This class is derived from AliJob, allowing a task based processing.
//    After the conversion of a raw data event into an IceEvent structure,
//    the processing of all available sub-tasks (if any) is invoked.
//    This provides an event-by-event (sub)task processing before the
//    final data structures are written out.
// 2) The main object in this job environment is an IceEvent* pointer.

 if (!fInfiles)
 {
  cout << " *IceRawTWR Exec* No data input file(s) specified." << endl;
  return;
 }

 Int_t ninfiles=fInfiles->GetEntries();
 if (!ninfiles)
 {
  cout << " *IceRawTWR Exec* No data input file(s) specified." << endl;
  return;
 }

 TTree* otree=0;
 if (fOutfile)
 {
  otree=new TTree("T","TWR raw data converted to IceEvent structures");
  otree->SetDirectory(fOutfile);
 }

 IceEvent* evt=new IceEvent();
 evt->SetTrackCopy(1);
 evt->SetDevCopy(1);

 // Branch in the tree for the event structure
 if (otree) otree->Branch("IceEvent","IceEvent",&evt,fBsize,fSplit); 

 // Initialise the job working environment
 SetMainObject(evt);
 if (fOutfile)
 {
  AddObject(fOutfile);
  AddObject(otree);
 }

 TString inputfile;

 cout << " ***" << endl;
 cout << " *** Start processing of job " << GetName() << " ***" << endl;
 cout << " ***" << endl;
 for (Int_t i=0; i<ninfiles; i++)
 {
  TObjString* sx=(TObjString*)fInfiles->At(i);
  if (!sx) continue;
  inputfile=sx->GetString(); 
  cout << " TWR raw data input file : " << inputfile.Data() << endl;
 }
 cout << " Maximum number of events to be processed : " << fMaxevt << endl;
 cout << " Print frequency : " << fPrintfreq << endl;
 if (fOutfile)
 {
  cout << " ROOT output file : " << fOutfile->GetName() << endl;
  cout << " Output characteristics : splitlevel = " << fSplit << " buffersize = " << fBsize << endl;
 }

 ListEnvironment();

 // Storage of the used parameters in the IceRawTWR device
 AliDevice params;
 params.SetNameTitle("IceRawTWR","IceRawTWR processor parameters");
 params.SetSlotName("Nchannels",1);
 params.SetSlotName("Ntriggers",2);
 params.SetSlotName("BaselineOffset",3);
 params.SetSignal(float(N_OF_CHANNELS),1);
 params.SetSignal(float(N_OF_TRIGGERS),2);
 params.SetSignal(float(BASELINE_MEAN_MAGIC),3);

 // Set DAQ device info
 AliDevice daq;
 daq.SetName("Daq");
 daq.SetSlotName("TWR",1);
 daq.SetSignal(1,1);

 twr_raw_data_file_t twr_file;
 Int_t year,runnum,evtnum;

 Int_t error;
 UInt_t nhead;

 GPS_t gps;
 UInt_t gpslow,gpshigh,gpssecs; // The GPS time information
 Int_t seconds,nsecs;           // Seconds and nanoseconds since start of the UT year

 Int_t nevt=0;
 fHeader=0;
 for (Int_t ifile=0; ifile<ninfiles; ifile++)
 {
  TObjString* sx=(TObjString*)fInfiles->At(ifile);
  if (!sx) continue;

  inputfile=sx->GetString(); 
  if (inputfile=="") continue;

  // Open the TWR raw data input file in binary mode
  fInput=fopen(inputfile.Data(),"rb");

  if (!fInput)
  {
   cout << " *IceRawTWR Exec* No input file found with name : " << inputfile.Data() << endl;
   continue;
  }

  // Extract info like run number, file number etc... from filename
  extract_info_from_filename((char*)inputfile.Data(),&twr_file);

  year=twr_file.year;
  runnum=twr_file.run_no;

  // Initialise the event structure 
  clear_event(&fEvent);

  // Read the file header information
  error=read_header_from_file(fInput,&fHeader,&nhead);

  if (error || !nhead) 
  {
   cout << " *IceRawTWR Exec* Error in header for input file : " << inputfile.Data() << endl;
   continue;
  }

  // Correct the mapping
  update_system(fHeader,runnum);
 
  while (!read_event(fInput,fHeader,&fEvent))
  {
   if (fMaxevt>-1 && nevt>=fMaxevt) break;

   evtnum=fEvent.eventcounter;

   // The GPS telegram info
   gps=fEvent.gps;
   gpslow=gps.seconds;            // The low 24 bits of the seconds count
   gpshigh=gps.info.bits.seconds; // The high 8 bits of the seconds count
   gpssecs=gpshigh<<24;
   gpssecs+=gpslow;

   // Seconds and nanoseconds since the start of the UT year
   seconds=gpssecs;
   nsecs=100*gps.count_10MHz;
   
   // Reset the complete Event structure
   evt->Reset();

   evt->SetRunNumber(runnum);
   evt->SetEventNumber(evtnum);
   evt->SetUT(year,0,seconds,nsecs);

   evt->AddDevice(params);
   evt->AddDevice(daq);

   PutTrigger(year);

   PutWaveforms(year);

   // Invoke all available sub-tasks (if any)
   CleanTasks();
   ExecuteTasks(opt);

   if (fPrintfreq)
   {
    if (!(nevt%fPrintfreq)) evt->HeaderData();
   }

   // Write the complete structure to the output Tree
   if (otree) otree->Fill();

   // Update event counter
   nevt++;

   // Reset the raw event structure 
   clear_event(&fEvent);
  } // End of event reading loop

  // Delete the file header structure
  clear_system(fHeader);

  if (fMaxevt>-1 && nevt>=fMaxevt) break;

 } // End of input file loop

 // Flush possible memory resident data to the output file
 if (fOutfile) fOutfile->Write();

 // Remove the IceEvent object from the environment
 // and delete it as well
 if (evt)
 {
  RemoveObject(evt);
  delete evt;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::PutWaveforms(Int_t year)
{
// Get the waveform info from the raw data event into the IcePack structure.

 IceEvent* evt=(IceEvent*)GetMainObject();
 if (!evt) return;

 // Loop over all the waveforms and add the histo(s) to the corresponding OM's
 TH1F histo;
 Int_t nbins=0;
 Float_t xlow=0;
 Float_t xup=0;
 TString hname;
 IceAOM om;
 IceAOM* omx=0;
 Int_t omid;
 Int_t omidmax=680;
 Int_t error;
 Float_t baseline;
 for (Int_t i=0; i<N_OF_CHANNELS; i++)
 {
  if (!fEvent.wfm_filled[i]) continue;

  omid=fEvent.twr_id_of_om[i];
  if (omid<=0 || omid>omidmax) continue; // Skip trigger channels

  // Get corresponding device from the current event structure  
  omx=(IceAOM*)evt->GetIdDevice(omid);
  if (!omx)
  {
   om.Reset(1);
   om.SetUniqueID(omid);
   evt->AddDevice(om);
   omx=(IceAOM*)evt->GetIdDevice(omid);
  }

  if (!omx) continue;

  clear_waveform_analysis(&fWform);
  error=restore_waveform(fEvent.wfm[i],&fWform,year);

  if (error) continue;

  baseline=fWform.frag_mean[0];

  hname="BASELINE-WF";
  hname+=omx->GetNwaveforms()+1;
  omx->AddNamedSlot(hname);
  omx->SetSignal(baseline,hname);

  // Fill the waveform histogram
  hname="OM";
  hname+=omid;
  hname+="-WF";
  hname+=omx->GetNwaveforms()+1;

  histo.Reset();
  histo.SetName(hname.Data());
  nbins=fWform.n_point;
  xlow=fWform.wfm_x[0];
  xup=fWform.wfm_x[nbins-1];
  histo.SetBins(nbins,xlow,xup);

  for (Int_t jbin=1; jbin<=nbins; jbin++)
  {
   histo.SetBinContent(jbin,baseline-fWform.wfm_y[jbin-1]);
  }

  omx->SetWaveform(&histo,omx->GetNwaveforms()+1);
 }
}
///////////////////////////////////////////////////////////////////////////
void IceRawTWR::PutTrigger(Int_t year)
{
// Get the trigger info from the raw data event into the IcePack structure.
// Currently only the trigger settings for the years 2005 and 2006 have been
// implemented.
// In addition to the hardware and software triggers as encountered in the
// raw data, an artificial "main" trigger has been introduced.
// This artificial "main" trigger is just an "or" of the standard hard and soft
// triggers (except calibration and random triggers) and serves only to
// provide a generic "main" trigger a la Amanda mu-daq so that the default
// "IceCleanHits" hit cleaning procedure will work correctly.
// The trigger time for the artificial "main" trigger is taken to be the
// time of the earliest hardware trigger pulse. In case there is no hardware
// trigger pulse available, the "main" trigger time is set to 0.
// For other years, only the artificial "main" trigger with a trigger time
// set to 0 will be stored in the IceEvent structure.

 // Fill the trigger structure
 Int_t error=retrigger(&fEvent,&fTrigger);
 if (error) return;

 IceEvent* evt=(IceEvent*)GetMainObject();
 if (!evt) return;

 AliDevice trig;
 trig.SetNameTitle("Trigger","Amanda/IceCube event triggers");
 AliSignal s;
 Float_t trigtime=0;

 if (year !=2005 && year != 2006)
 {
  s.SetName("main");
  s.SetUniqueID(0);
  s.SetSlotName("trig_pulse_le",1);
  s.SetSignal(trigtime,1);
  trig.AddHit(s);
  // Store the trigger data into the IceEvent structure
  evt->AddDevice(trig);
  return;
 }

 // Trigger settings for 2005 and 2006
 if (!fTrigger.n_software_trigger && !fTrigger.n_hardware_trigger) return;

 TString trignames[N_OF_TRIGGERS]={"m24","m18","string","spase","cal-t0","cal-la","m12",
                                   "main-logic","main-or","random","m20-frag","volume"};
 Int_t imain=0;
 for (Int_t i=0; i<N_OF_TRIGGERS; i++)
 {
  if (!fTrigger.trigger_active[i]) continue;

  s.Reset(1);
  s.SetName(trignames[i]);
  s.SetUniqueID(i);
  trigtime=0;
  if (fTrigger.trigger_has_pulse[i]) trigtime=fTrigger.trigger_time[i];
  s.SetSlotName("trig_pulse_le",1);
  s.SetSignal(trigtime,1);
  trig.AddHit(s);
  // Set flag to indicate creation of artificial "main" trigger
  if (i!=4 && i!=5 && i!=9) imain=1;
 }

 // Set the artificial "main" trigger
 if (imain)
 {
  s.Reset(1);
  s.SetName("main");
  s.SetUniqueID(N_OF_TRIGGERS);
  s.SetSlotName("trig_pulse_le",1);
  trigtime=0;
  if (fTrigger.first_trigger>=0) trigtime=fTrigger.first_trigger_time;
  s.SetSignal(trigtime,1);
  trig.AddHit(s);
 }

 // Store the trigger data into the IceEvent structure
 evt->AddDevice(trig);
}
///////////////////////////////////////////////////////////////////////////
Int_t IceRawTWR::extract_info_from_filename(char* fname,twr_raw_data_file_t* twr_file)
{
  char start_str[20],year_str[20],day_str[20],run_no_str[20], 
       file_no_str[20],begin_str[20],end_str[20];
  char* filename;

  filename = strstr(fname, "twr");
  if(filename == NULL)
  if(strncmp("twr_", start_str, 4)) 
    {
      printf("%s\n", filename);
      return(ERROR_NOT_VALID_FILENAME);
    }

  strncpy(start_str, filename, 4);
  if(strncmp("twr_", start_str, 4)) 
    {
      printf("%s %s\n", filename, start_str);
      return(ERROR_NOT_VALID_FILENAME);
    }
  strncpy(year_str, &filename[4], 4);
  twr_file->year = strtol(year_str, 0, 10);

  if(twr_file->year==2003)
    {
      strncpy(day_str, &filename[9], 3);
      day_str[3] = '\0';
      twr_file->day = strtol(day_str, 0, 10);
      
      strncpy(run_no_str, &filename[13], 4);
      run_no_str[4] = '\0';
      twr_file->run_no = strtol(run_no_str, 0, 10);
      
      strncpy(file_no_str, &filename[18], 4);
      file_no_str[4] = '\0';
      twr_file->file_no = strtol(file_no_str, 0, 10);
    }
  
  if(twr_file->year==2004)
    {
      strncpy(day_str, &filename[9], 3);
      day_str[3] = '\0';
      twr_file->day = strtol(day_str, 0, 10);
      
      strncpy(run_no_str, &filename[13], 4);
      run_no_str[4] = '\0';
      twr_file->run_no = strtol(run_no_str, 0, 10);
      
      strncpy(file_no_str, &filename[18], 4);
      file_no_str[4] = '\0';
      twr_file->file_no = strtol(file_no_str, 0, 10);
      
      strncpy(begin_str, &filename[23], 5);
      begin_str[5] = '\0';
      twr_file->begin = strtol(begin_str, 0, 10);
      
      strncpy(end_str, &filename[29], 5);
      end_str[5] = '\0';
      twr_file->end = strtol(end_str, 0, 10);
    }

  if(twr_file->year > 2004)      
    {
      strncpy(day_str, &filename[9], 3);
      day_str[3] = '\0';
      twr_file->day = strtol(day_str, 0, 10);
      
      strncpy(run_no_str, &filename[13], 6);
      run_no_str[6] = '\0';
      twr_file->run_no = strtol(run_no_str, 0, 10);
      
      strncpy(file_no_str, &filename[20], 4);
      file_no_str[4] = '\0';
      twr_file->file_no = strtol(file_no_str, 0, 10);
      
      strncpy(begin_str, &filename[25], 5);
      begin_str[5] = '\0';
      twr_file->begin = strtol(begin_str, 0, 10);
      
      strncpy(end_str, &filename[31], 5);
      end_str[5] = '\0';
      twr_file->end = strtol(end_str, 0, 10);
    }
  return(0);
}
///////////////////////////////////////////////////////////////////////////
Int_t IceRawTWR::clear_system(sys_config_t* sys)
{
// Deletion of the file header structure.

 if (!sys) return 0;

 for(Int_t icrate=0; icrate < int(sys->n_crates); icrate++)
 {
  if (!sys->crate[icrate]) continue;
  for(Int_t itwr=0; itwr < int(sys->crate[icrate]->n_twr); itwr++)
  {
   if (sys->crate[icrate]->twr[itwr]) delete sys->crate[icrate]->twr[itwr];
  }
  delete sys->crate[icrate];
 }
 delete sys;
 sys=0;
 return 0;
}
///////////////////////////////////////////////////////////////////////////
Int_t IceRawTWR::clear_event(event_t* event_ptr)
{
  Int_t i_value;
  Int_t *int_ptr = (int*) event_ptr;

  for(i_value=0; i_value < int(sizeof(event_t)/sizeof(Int_t)); i_value++)
    {
      *int_ptr++ = 0;
    }
  return(0);
}
///////////////////////////////////////////////////////////////////////////
Int_t IceRawTWR::read_header_from_file(FILE* fin,sys_config_t** system_ptr,UInt_t* header_length)
{
  Int_t i_crate, i_twr, i_channel;
  UInt_t count_twr_in_system = 0;
  UInt_t dummy;
  
  sys_config_t *sys;

  // allocating memory for sys_config structure
  sys = (sys_config_t*) malloc( sizeof(sys_config_t) );

  fread(&dummy,sizeof(UInt_t),1,fin); // Header Begin Mark

  fread(header_length,sizeof(UInt_t),1,fin);   // Length of header
  fread(&sys->clockdiv,sizeof(UInt_t),1,fin);  
  fread(&sys->n_crates,sizeof(UInt_t),1,fin);   

  if( (sys->n_crates > MAX_N_CRATES) || (sys->n_crates < 0) )
    return(ERROR_TOO_MANY_CRATES);

  for(i_crate=0; i_crate < int(sys->n_crates); i_crate++)
    {
      sys->crate[i_crate] = 
	(crate_config_t*) malloc( sizeof(crate_config_t) );

      fread(&sys->crate[i_crate]->vme_base_bridge,sizeof(UInt_t),1,fin); 
      fread(&sys->crate[i_crate]->vme_base_100MHz,sizeof(UInt_t),1,fin); 
      fread(&sys->crate[i_crate]->base_gps,sizeof(UInt_t),1,fin); 
      fread(&sys->crate[i_crate]->n_twr,sizeof(UInt_t),1,fin); 
      
      if( (sys->crate[i_crate]->n_twr > MAX_N_TWR_PER_CRATE) 
	  || (sys->crate[i_crate]->n_twr < 0) )
	return(ERROR_TOO_MANY_TWRS);

      for(i_twr=0; i_twr < int(sys->crate[i_crate]->n_twr); i_twr++)
	{
	  sys->crate[i_crate]->twr[i_twr] = 
	    (twr_config_t*) malloc( sizeof(twr_config_t) );
	  count_twr_in_system++;
	  fread(&sys->crate[i_crate]->twr[i_twr]->base, 
	       sizeof(UInt_t),1,fin);
	  fread(&sys->crate[i_crate]->twr[i_twr]->id, 
	       sizeof(UInt_t),1,fin);

	  sys->crate[i_crate]->twr[i_twr]->id 
	    = sys->crate[i_crate]->twr[i_twr]->id - 0x10; /* Correct */


	  fread(&dummy,sizeof(UInt_t),1,fin); /* stat_reg */
	  fread(&sys->crate[i_crate]->twr[i_twr]->mod_id, 
	       sizeof(UInt_t),1,fin);
	  fread(&dummy,sizeof(UInt_t),1,fin); /* acq_ctrl */
	  fread(&sys->crate[i_crate]->twr[i_twr]->ext_start, 
	       sizeof(UInt_t),1,fin);
	  fread(&sys->crate[i_crate]->twr[i_twr]->ext_stop, 
	       sizeof(UInt_t),1,fin);
	  fread(&dummy,sizeof(UInt_t),1,fin); /* evtconfig */

	  for(i_channel = 0; i_channel < CHANNELS_PER_TWR; i_channel++)
	    {
	      fread(&sys->crate[i_crate]->twr[i_twr]->om_no[i_channel], 
		   sizeof(UInt_t),1,fin);
	    }

	  for(i_channel = 0; i_channel < CHANNELS_PER_TWR; i_channel++)
	    {
	      fread(&sys->crate[i_crate]->twr[i_twr]->om_is_optical[i_channel], 
		   sizeof(UInt_t),1,fin);
	    }

	  for(i_channel = 0; i_channel < CHANNELS_PER_TWR; i_channel++)
	    {
	      fread(&sys->crate[i_crate]->twr[i_twr]->baseline[i_channel], 
		   sizeof(UInt_t),1,fin);
	    }

	  for(i_channel = 0; i_channel < CHANNELS_PER_TWR; i_channel++)
	    {
	      fread(&sys->crate[i_crate]->twr[i_twr]->threshold[i_channel],
		   sizeof(UInt_t),1,fin);
	    }

	  sys->twr_field[(i_crate * 0x10) + i_twr]
	      = sys->crate[i_crate]->twr[i_twr];
	  
	  /* Bug fix needed */
	  for(i_channel=0; i_channel < 8; i_channel++)
	  {
	     if( sys->crate[i_crate]->twr[i_twr]->om_no[i_channel] == 9000 )
		 sys->crate[i_crate]->twr[i_twr]->om_no[i_channel] 
		   = N_OF_CHANNELS - 1; 
	  }
	}
    }

  // Set number of TWRs in system
  sys->n_twr = count_twr_in_system;

  *system_ptr = sys;
  return(0);
}
///////////////////////////////////////////////////////////////////////////
Int_t IceRawTWR::update_system(sys_config_t* sys,Int_t run_number)
{
  Int_t i_crate, i_twr, i_channel;

  /* Data for bug fix 1 */
  UInt_t om_no_r1[CHANNELS_PER_TWR] 
    = {111, 112, 113, 114, 115, 116, 39, 118}; 
  UInt_t om_is_optical_r1[CHANNELS_PER_TWR] 
     = {0, 0, 0, 0, 0, 0, 0, 0};
  UInt_t threshold_r1[CHANNELS_PER_TWR]
    = {50, 50, 50, 50, 50, 50, 80, 50};

  UInt_t om_no_r2[CHANNELS_PER_TWR] 
    = {473, 484, 485, 486, 487, 475, 490, 491}; 
  UInt_t om_is_optical_r2[CHANNELS_PER_TWR] 
     = {1, 1, 1, 1, 1, 1, 1, 1};
  UInt_t threshold_r2[CHANNELS_PER_TWR]
    = {15, 50, 55, 40, 15, 23, 15, 15};

  
  /* Bugfix 1 Andreas Bug */ 

  /*
    By accident this TWR was counted twice in TWR.cnf
    as Crate 0 TWR 7 and Crate 4 TWR 7
    from run up to run
    TWR_OM          639     642     1       9       10      11      12      30
    OPTICAL         0       0       0       0       0       0       0       0
    TWR_BASELINE    110     120     110     140     150     160     170     180
    TWR_THRESHOLD   50      50      80      80      80      80      80      80

    Crate 4 TWR 7 should be replaced with this TWR
    TWR_OM          111     112     113     114     115     116     39      118
    OPTICAL         0       0       0       0       0       0       0       0
    TWR_BASELINE    110     120     130     140     150     160     170     180
    TWR_THRESHOLD   50      50      50      50      50      50      80      50
  */

  if( 
     (run_number >= 9153 )  /* Begin season 2005 13.2.05 */  
     && (run_number < 9800) /* Timo corrected TWR.cnf on after run ??? */
     /* Need to find exact date */
     )
    {
      i_crate = 4;
      i_twr = 7;
      for(i_channel = 0; i_channel < CHANNELS_PER_TWR; i_channel++)
	{
	  sys->crate[i_crate]->twr[i_twr]->om_no[i_channel]         
	    = om_no_r1[i_channel]; 
	  sys->crate[i_crate]->twr[i_twr]->om_is_optical[i_channel] 
	    = om_is_optical_r1[i_channel]; 
	  sys->crate[i_crate]->twr[i_twr]->threshold[i_channel]     
	    = threshold_r1[i_channel];
	}
    }

  /* Bugfix 2 Timos Bug */ 

  /*
    By accident this TWR was counted twice in TWR.cnf
    as Crate 0 TWR 1 and Crate 5 TWR b
    from run 9153 up to run 9188

    TWR_OM          492     493     495     496     497     499     500     501
    OPTICAL         1       1       1       1       1       1       1       1
    TWR_BASELINE    110     120     130     140     150     160     170     180
    TWR_THRESHOLD   16      45      25      42      35      46      15      15

    Crate 5 TWR b should be corrected to 
    TWR_OM          473     484     485     486     487     475     490     491
    OPTICAL         1       1       1       1       1       1       1       1
    TWR_BASELINE    4000    120     130     140     150     4000    170     180
    TWR_THRESHOLD   15      50      55      40      15      23      15      15
  */

  if( 
     (run_number >= 9153 )  /* Begin season 2005 = Feb 2nd 05 */  
     && (run_number < 9189) /* Timo corrected TWR.cnf on      */
     /* Mar 15th 05 = day 74 after run 9188                   */
     )
    {
      i_crate = 5;
      i_twr = 0xb;
      for(i_channel = 0; i_channel < CHANNELS_PER_TWR; i_channel++)
	{
	  sys->crate[i_crate]->twr[i_twr]->om_no[i_channel]         
	    = om_no_r2[i_channel]; 
	  sys->crate[i_crate]->twr[i_twr]->om_is_optical[i_channel] = 
	    om_is_optical_r2[i_channel]; 
	  sys->crate[i_crate]->twr[i_twr]->threshold[i_channel]     = 
	    threshold_r2[i_channel];
	}
    }
 return(0);
}
///////////////////////////////////////////////////////////////////////////
Int_t IceRawTWR::read_event(FILE* fin,sys_config_t* sys,event_t* event_ptr)
{
    Int_t i_wfm;
    UInt_t length_of_event_block;
    
    Int_t n_twr, n_of_waveforms_in_event, read_number;
    UInt_t length_wfm[CHANNELS_PER_TWR];
    UInt_t dummy, channel_no, om_no, twr_no;
    
    // Reset waveform filled register
    memset(&event_ptr->wfm_filled[0], 0, sizeof(UInt_t) * N_OF_CHANNELS);

    if( !fread(&dummy,sizeof(UInt_t),1,fin) ) return(1);

    if(dummy != 0xbbbbbbbb) 
    {
	printf("Wrong event begin mark %x\n", dummy); 
	while( (dummy !=0xbbbbbbbb) 
	       && (fread(&dummy,sizeof(UInt_t),1,fin) != 0) )
	  {;//printf("dummy:%x\n", dummy);
	  }
    }
    if( !fread(&length_of_event_block,sizeof(UInt_t),1,fin) ) return(1);
    if( !fread(&event_ptr->eventcounter,sizeof(UInt_t),1,fin) ) return(1);
    if( !fread(&event_ptr->which_trigger,sizeof(UInt_t),1,fin) ) return(1);
    if( !fread(&event_ptr->gps,sizeof(GPS_t),1,fin) ) return(1);
	
    // --reading waveforms from TWR blocks
    n_twr = 0;
    while(n_twr < int(sys->n_twr))
    {
	// --read TWR header
	if( !fread(&dummy,sizeof(UInt_t),1,fin) ) return(1);
	if(dummy != 0xffffffff) 
	{printf("Wrong twr begin mark %x\n", dummy); return(2);}
	if( !fread(&twr_no,sizeof(UInt_t),1,fin) ) return(1);

        // nur voruebergehend !!
	twr_no -= 0x10;

	if( !fread(&event_ptr->twr[twr_no].timestamp,sizeof(UInt_t),1,fin) ) 
	    return(1);
	if( !fread(&n_of_waveforms_in_event,sizeof(UInt_t),1,fin) ) 
	    return(1);
	event_ptr->twr[twr_no].n_wfm = n_of_waveforms_in_event;

	for(i_wfm=0; i_wfm < n_of_waveforms_in_event; i_wfm++)
	{
	    if( !fread(&length_wfm[i_wfm],sizeof(UInt_t),1,fin) ) return(1);
	}
	    
	// read waveforms
	for(i_wfm=0; i_wfm < n_of_waveforms_in_event; i_wfm++)
	{
	    if(length_wfm[i_wfm] != 0)
	    {
	      if( !fread(&channel_no,sizeof(UInt_t),1,fin) ) return(1);
	      if(sys->twr_field[twr_no]->om_no[channel_no] 
		   < N_OF_CHANNELS)
		    om_no = sys->twr_field[twr_no]->om_no[channel_no];
		else
		    om_no = N_OF_CHANNELS-1;

		/* Fix needed */

		event_ptr->twr_id_of_om[om_no] = twr_no;

		read_number = fread(&event_ptr->wfm[om_no], 
				   length_wfm[i_wfm]-sizeof(UInt_t),1,fin);
		event_ptr->wfm_filled[om_no] = 1;
		if( !read_number ) return(1);

                // read_number correction for usage of fread() instead of read()
                read_number*=length_wfm[i_wfm]-sizeof(UInt_t);

		if( read_number != int(length_wfm[i_wfm]-sizeof(UInt_t)) ) 
                {
                  cout << " read_number : " << read_number
                       << " length_wfm["<<i_wfm<<"] : " << length_wfm[i_wfm]
                       << " sizeof(UInt_t) : " << sizeof(UInt_t) << endl;
		  return(2);
                }
	    }	
	}
	n_twr++;
    } // end while n_twr
    return(0);
}
///////////////////////////////////////////////////////////////////////////
Int_t IceRawTWR::retrigger(event_t* ev,trigger_hits_t* trig)
{
// Returns the active trigger(s)

 // Initialise the trigger_hits_t structure with zeroes
 memset(trig, 0, sizeof(trigger_hits_t) );

 // Obtain the software trigger info
 trig->n_software_trigger=0;
 for(Int_t itrigger=0; itrigger<N_OF_TRIGGERS; itrigger++)
 {
  if(ev->which_trigger & trigger_bits[itrigger])
  {
   //printf("SetTrigger %i\n", i_trigger);
   trig->trigger_active[itrigger]=1;
   trig->n_software_trigger++;
  } 
  else
  {
   trig->trigger_active[itrigger]=0;
  }
 }

 // Obtain the hardware trigger info
 trig->n_hardware_trigger=0;
 trig->first_trigger_time=10000000;
 trig->first_trigger=-1;

 for(Int_t jtrigger=0; jtrigger<N_OF_TRIGGERS; jtrigger++)
 {
  if(!trigger_channel[jtrigger]) continue;

  if(ev->wfm_filled[trigger_channel[jtrigger]])
  {
   trig->trigger_active[jtrigger]=1;
   trig->trigger_time[jtrigger]=(ev->wfm[trigger_channel[jtrigger]].value[2] & 0xfff);
   trig->trigger_has_pulse[jtrigger]=1;
   if (trig->trigger_time[jtrigger] < trig->first_trigger_time)
   {
    trig->first_trigger_time=trig->trigger_time[jtrigger];
    trig->first_trigger=jtrigger;
   }
   trig->n_hardware_trigger++;
  }
 }
 return 0;
}
///////////////////////////////////////////////////////////////////////////
Int_t IceRawTWR::clear_waveform_analysis(waveform_analyse_t* wfm_om)
{
  Int_t i_value, i_frag, i_edge, i_peak;

  if(wfm_om == 0) return(1);

  // output from analysis
  wfm_om->n_frag = 0;
  for(i_frag=0; i_frag < MAX_N_OF_FRAGS; i_frag++)
    {
      wfm_om->frag_n_points[i_frag] = 0;
      wfm_om->frag_begin[i_frag] = 0;
      wfm_om->frag_end[i_frag] = 0;
      wfm_om->frag_mean[i_frag] = 0;
      wfm_om->frag_begin_time[i_frag] = 0;
    }
  
  wfm_om->n_peak = 0;
  for(i_peak=0; i_peak < MAX_N_OF_PEAKS; i_peak++)
    {
      wfm_om->peak_begin[i_peak] = 0;
      wfm_om->peak_end[i_peak] = 0;
      wfm_om->peak_max[i_peak] = 0;
      wfm_om->peak_TDC_edge[i_peak] = 0;
      wfm_om->peak_local_minimum[i_peak] = 0;
      wfm_om->crosstalk_charge_n_value[i_peak] = 0;
      wfm_om->peak_in_fragment[i_peak] = 0;
       
      wfm_om->peak_mean[i_peak] = 0.0;
    
      wfm_om->peak_m[i_peak] = 0.0;
      wfm_om->peak_b[i_peak] = 0.0;
      wfm_om->peak_t0[i_peak] = 0.0;
      wfm_om->peak_begin_time[i_peak] = 0.0; 
      wfm_om->peak_charge[i_peak] = 0.0;
      wfm_om->peak_height[i_peak] = 0.0;
      wfm_om->fitted_amplitude[i_peak] = 0.0;
      wfm_om->fitted_TOT[i_peak] = 0.0;
      wfm_om->crosstalk_charge[i_peak] = 0.0;
      wfm_om->crosstalk_slope[i_peak] = 0.0;
    }
  
  wfm_om->n_point = 0;
  wfm_om->wfm_min = 4095; 
  wfm_om->wfm_max = 0; 
  wfm_om->b_out_of_range = 0;
  
  for(i_value=0; i_value < 1024; i_value++)
    {
      wfm_om->wfm_x[i_value] = 0;
      wfm_om->wfm_y[i_value] = 0;
    }

  wfm_om->n_tdc_edges = 0;
  for(i_edge=0; i_edge < MAX_N_OF_TDC_EDGES; i_edge++)
    {
      wfm_om->leading_edge[i_edge] = 0.0;
      wfm_om->falling_edge[i_edge] = 0.0;
      wfm_om->identified_twr_hit[i_edge] = -1;
    }

 return(0);
}
///////////////////////////////////////////////////////////////////////////
Int_t IceRawTWR::restore_waveform(waveform_t f_wfm,waveform_analyse_t* wfm_om,Int_t year)
{
    UShort_t wfm_length, mean;
    static UShort_t tmp_wf[2000];

    Int_t debug = 0;
    Int_t fragment_start = 0; 
    Int_t frag_count = 0;     // position in current fragment
    Int_t n_position = 0;     // position in displayed waveform
    UInt_t  n_word = 2;      // position in featured waveform
    Int_t n_fragment = 0;     // actual fragment 
    Int_t b_wrong_value = 0;

    UShort_t assumed_frag_begin, last_value; /* bug in eventbuilder */

    wfm_om->wfm_min = 4095.0;
    wfm_om->wfm_max = 0.0;

    if( (f_wfm.value[0] & 0xf000) != 0xf000 ) return(1);
    wfm_length = (f_wfm.value[0] & 0xfff)/2;

    mean = f_wfm.value[1] + BASELINE_MEAN_MAGIC;    
    while( ((f_wfm.value[n_word] & 0xf000) == 0x4000) && 
	   (n_word < wfm_length) &&
	   (n_fragment < MAX_N_OF_FRAGS) )
    {
      fragment_start = f_wfm.value[n_word] & 0xfff;
      n_word++;
      wfm_om->frag_begin_time[n_fragment] 
	= fragment_start * NSECS_PER_TWR_BIN;
      wfm_om->frag_begin[n_fragment] = n_position;
      wfm_om->frag_mean[n_fragment] = mean;
      
      b_wrong_value = 0;
      frag_count = 0;
      
      while( ((f_wfm.value[n_word] & 0xf000) != 0x2000) &&
	     ((f_wfm.value[n_word] & 0xf000) != 0x4000) &&/*Reconstructable*/
	     !b_wrong_value && /* Buggy */
	     (n_word < wfm_length) )
	{
	  if(year > 2004)
	    {
	      /* 2005 2006 data */
	      if(frag_count == 0)
		{
		  tmp_wf[n_word] = f_wfm.value[n_word] + mean;
		  wfm_om->wfm_y[n_position] = (float) tmp_wf[n_word];
		  wfm_om->wfm_x[n_position] = (float) 
		    wfm_om->frag_begin_time[n_fragment]
		    + (frag_count * NSECS_PER_TWR_BIN);
		}
	      else if(frag_count == 1)
		{
		  tmp_wf[n_word] = f_wfm.value[n_word] + tmp_wf[n_word-1];
		  wfm_om->wfm_y[n_position] = (float) tmp_wf[n_word];
		  wfm_om->wfm_x[n_position] = (float) 
		    wfm_om->frag_begin_time[n_fragment]
		    + (frag_count * NSECS_PER_TWR_BIN);
		}
	      else 
		{
		  tmp_wf[n_word] = 
		    2*tmp_wf[n_word-1] + f_wfm.value[n_word];
		  tmp_wf[n_word] -= tmp_wf[n_word-2];
		  
		  wfm_om->wfm_y[n_position] = (float) tmp_wf[n_word];
		  wfm_om->wfm_x[n_position] = (float) 
		    wfm_om->frag_begin_time[n_fragment]
		    + (frag_count * NSECS_PER_TWR_BIN);
		}
	      
	      
	      /*
		Hack for wrongly merged overlapping fragments
	      */
	      if(tmp_wf[n_word] > 0x1fff)
		{
		  /* BUG FIXXXX                                         */
		  /* assume that fragment merge in eventbuilder caused  */
		  /* problem two fragments overlap in EXACTLY ONE point */
		  /* and are merged first point of the added part of    */
		  /* the fragment is encoded using the former fragment  */
		  /* start as a data point                              */
		  
		  last_value         = tmp_wf[n_word-1];
		  assumed_frag_begin = 0x4000 + fragment_start + frag_count;
		  tmp_wf[n_word] =  f_wfm.value[n_word] + 2 * last_value;
		  tmp_wf[n_word] -= assumed_frag_begin;
		  wfm_om->wfm_y[n_position] = (float) tmp_wf[n_word];
		  
		  /* Look if value is still buggy */
		  if(tmp_wf[n_word] > 0x1fff) b_wrong_value = 1;
		  
		  debug = ERROR_MISS_FRAG_STOP;
		}
	    } /* end year >= 2005 */
	  else
	    {
	      /* 2003 2004 data */
	      wfm_om->wfm_y[n_position] = (float) f_wfm.value[n_word];
	      wfm_om->wfm_x[n_position] = (float) 
		wfm_om->frag_begin_time[n_fragment]
		+ (frag_count * NSECS_PER_TWR_BIN);
	    } /* end year 2003 2004 */
	  
	  /* Set min and max Y */
	  
	  if(wfm_om->wfm_y[n_position] > wfm_om->wfm_max)
	    wfm_om->wfm_max = wfm_om->wfm_y[n_position];
	  if(wfm_om->wfm_y[n_position] < wfm_om->wfm_min)
	    wfm_om->wfm_min = wfm_om->wfm_y[n_position];
	  
	  n_position++;
	  n_word++;
	  frag_count++;
	}

	if((f_wfm.value[n_word] & 0xf000) == 0x2000) /* Normal wavf */
	  {

	    wfm_om->frag_end[n_fragment] = n_position - 1; 
	    wfm_om->frag_n_points[n_fragment] = 
	      wfm_om->frag_end[n_fragment] 
	      - wfm_om->frag_begin[n_fragment] + 1;
	    wfm_om->n_point += wfm_om->frag_n_points[n_fragment];
	    n_word++;
	  }
	else 
	  return(ERROR_CORRUPTED_WF);

	n_fragment++;
    } /* end while fragment */


    wfm_om->n_frag = n_fragment;
    if( !(n_word & 0x1) ) n_word++;

    if(n_fragment >= MAX_N_OF_FRAGS) return(ERROR_MAX_N_FRAGMENTS_EXCEEDED);
    

    // Hack to get rid of last value of waveform always set to 0
    if (wfm_om->wfm_y[wfm_om->n_point] == 0.0)
    {
     // erase last point of waveform
     wfm_om->n_point--;

     // Shorten last pulse if necessary
     // if( wfm_om.peak_end[wfm_om.n_peak-1] 
     //     == wfm_om.frag_end[wfm_om.n_frag-1] )
     //   wfm_om.peak_end[wfm_om.n_peak-1]--;
		  
     // Shorten last fragment
     wfm_om->frag_n_points[wfm_om->n_frag-1]--;
     wfm_om->frag_end[wfm_om->n_frag-1]--;

     wfm_om->wfm_min = 4095.0;
     wfm_om->wfm_max = 0.0;
     for (Int_t i_value=0; i_value < wfm_om->n_point; i_value++)
     {
      if (wfm_om->wfm_y[i_value] > wfm_om->wfm_max) wfm_om->wfm_max=wfm_om->wfm_y[i_value];
      if (wfm_om->wfm_y[i_value] < wfm_om->wfm_min) wfm_om->wfm_min=wfm_om->wfm_y[i_value];
     }
    }

 return(debug);
}
///////////////////////////////////////////////////////////////////////////
