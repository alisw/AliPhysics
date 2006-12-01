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
// Class IceF2k
// Conversion of Amanda F2K data into IceEvent physics event structures.
// This class is derived from AliJob providing a task-based processing
// structure on an event-by-event basis.
// The main object in the job environment is an IceEvent* pointer.
// In case the user has provided sub-tasks, these will be executed
// on an event-by-event basis after the IceEvent structure has been filled
// with the F2K data and before the final structures are written out.
// Note that the data structures are only written out if an outputfile has
// been specified via the SetOutputFile memberfunction.
// In case no outputfile has been specified, this class provides a facility
// to investigate/analyse F2K data using the Ralice/IcePack analysis tools.
//
// Note : Sometimes the filtering/reco process which produced the F2K file
//        may have introduced a shift (i.e. offset) in the hit times w.r.t.
//        the actual trigger time. The aim of this is to obtain the hit times
//        centered more or less around zero.
//        In case of real data, this is recorded in the F2K data itself and
//        as such will be taken automatically into account by this IceF2k
//        processor such that all times will be provided again unshifted.
//        In other words, all times will be w.r.t. the actual trigger time
//        as recorded in the trigger data device named "Trigger" in the IceEvent
//        structure.
//        In case of simulated data this shift is not available in the F2K data.
//        The offset denoted in the F2K record is related to the time of the
//        primary interaction to put it well ahead of the detector trigger.
//        This primary interaction time, however, is irrelevant for the
//        reconstruction of the recorded hit patterns. 
//        If a user had introduced a shift in producing the MC data,
//        very frequently (but not always) a value of -19000 is used.
//        For the IceF2k processing, the user can manually introduce a
//        time offset in case of MC data via the memberfunction SetMcToffset().
//        This user defined offset value will then be used to correct all
//        the hit times such that they will be provided again unshifted
//        w.r.t. the actual trigger time as recorded in the device named 
//        "Trigger" in the IceEvent structure.
//        By default the MC time offset is set to 0 in the constructor
//        of this class.
//
// Usage example :
// ---------------
//
// Note : This example creates automatically the ROOT output file, which
//        is the most user friendly way of running the conversion job. 
//        In the subdirectory /macros the example macro icef2k.cc provides
//        an example of how to create a ROOT output file yourself and passing
//        this file via a pointer to IceF2k. 
//
// gSystem->Load("ralice");
// gSystem->Load("icepack");
// gSystem->Load("iceconvert");
//
// IceF2k q("IceF2k","F2K to IcePack data structure conversion");
//
// // Limit the number of entries for testing
// q.SetMaxEvents(10);
//
// // Print frequency to produce a short summary print every printfreq events
// q.SetPrintFreq(1);
//
// // Split level for the output structures
// q.SetSplitLevel(0);
//
// // Buffer size for the output structures
// q.SetBufferSize(32000);
//
// // The F2K input filename(s)
// q.AddInputFile("run7825.f2k");
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
// // Select various objects to be added to the output file
//
// TFile* ofile=q.GetOutputFile();
// 
// if (ofile)
// {
//  ofile->cd(); // Switch to the output file directory
//
//  AliObjMatrix* omdb=q.GetOMdbase();
//  if (omdb) omdb->Write();
//
//  AliDevice* fitdefs=q.GetFitdefs();
//  if (fitdefs) fitdefs->Write();
//
//  TDatabasePDG* pdg=q.GetPDG();
//  if (pdg) pdg->Write();
//
//  // Flush additional objects to the output file.
//  // The output file is not explicitly closed here
//  // to allow interactive investigation of the data tree
//  // when this macro is run in an interactive ROOT/CINT session.
//  ofile->Write();
// }
//
//--- Author: Nick van Eijndhoven 11-mar-2005 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceF2k.h"
#include "Riostream.h"

ClassImp(IceF2k) // Class implementation to enable ROOT I/O

IceF2k::IceF2k(const char* name,const char* title) : AliJob(name,title)
{
// Default constructor.
// By default maxevent=-1, split=0, bsize=32000, printfreq=1.

 fSplit=0;
 fBsize=32000;
 fMaxevt=-1;
 fPrintfreq=1;
 fInfiles=0;
 fOutfile=0;

 fPdg=0;
 fOmdb=0;
 fFitdefs=0;
 fTrigdefs=0;
 fToffset=0;
 fMctoffset=0;
 fMctracks=3;
}
///////////////////////////////////////////////////////////////////////////
IceF2k::~IceF2k()
{
// Default destructor.

 if (fInfiles)
 {
  delete fInfiles;
  fInfiles=0;
 }

 if (fPdg)
 {
  delete fPdg;
  fPdg=0;
 }

 if (fOmdb)
 {
  delete fOmdb;
  fOmdb=0;
 }

 if (fFitdefs)
 {
  delete fFitdefs;
  fFitdefs=0;
 }

 if (fTrigdefs)
 {
  delete fTrigdefs;
  fTrigdefs=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetMaxEvents(Int_t n)
{
// Set the maximum number of events to be processed.
// n=-1 implies processing of the complete input file, which is the default
// initialisation in the constructor.
 fMaxevt=n;
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetPrintFreq(Int_t f)
{
// Set the printfrequency to produce info every f events.
// f=1 is the default initialisation in the constructor.
 if (f>=0) fPrintfreq=f;
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetSplitLevel(Int_t split)
{
// Set the split level for the ROOT data file.
// split=0 is the default initialisation in the constructor.
 if (split>=0) fSplit=split;
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetBufferSize(Int_t bsize)
{
// Set the buffer size for the ROOT data file.
// bsize=32000 is the default initialisation in the constructor.
 if (bsize>=0) fBsize=bsize;
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetMcToffset(Float_t toffset)
{
// Set a user defined time offset for Monte Carlo data.
// A very frequently (but not always) used value is -19000.
// See the introductory docs of this class for further details.
 fMctoffset=toffset;
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SelectMcTracks(Int_t mode)
{
// User selection of MC tracks to be stored in the event structure.
//
// mode = 0 : No MC tracks are stored
//        1 : Only muon and muon-neutrino MC tracks are stored
//        2 : All lepton MC tracks are stored
//        3 : All MC tracks (incl. brems, pairprod etc...) are stored
//
// By default mode=3 is set in the constructor of this class.

 if (mode<0 || mode >3) return;
 fMctracks=mode;
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetInputFile(TString name)
{
// Set the name of the F2K input file.
// This function has become obsolete but is kept for backward compatibility.
// The user is advised to use AddInputFile() instead, which allows processing
// of multiple F2K input files.
// This function will reset the list of all F2K input files and put the specified
// filename at the first position.
// Additional F2K input files can be specified via AddInputFile().

 if (fInfiles) delete fInfiles;

 fInfiles=new TObjArray();
 fInfiles->SetOwner();

 TObjString* s=new TObjString();
 s->SetString(name);
 fInfiles->Add(s);
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::AddInputFile(TString name)
{
// Add the name of this F2K input file to the list to be processed.

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
void IceF2k::SetOutputFile(TFile* ofile)
{
// Set the output file for the ROOT data.
 if (fOutfile) delete fOutfile;
 fOutfile=ofile;
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetOutputFile(TString name)
{
// Create the output file for the ROOT data.
 if (fOutfile) delete fOutfile;
 fOutfile=new TFile(name.Data(),"RECREATE","F2K data in IceEvent structure");
}
///////////////////////////////////////////////////////////////////////////
TFile* IceF2k::GetOutputFile()
{
// Provide pointer to the ROOT output file.
 return fOutfile;
}
///////////////////////////////////////////////////////////////////////////
TDatabasePDG* IceF2k::GetPDG()
{
// Provide pointer to the PDG database
 return fPdg;
}
///////////////////////////////////////////////////////////////////////////
AliObjMatrix* IceF2k::GetOMdbase()
{
// Provide pointer to the OM geometry, calib. etc... database
 return fOmdb;
}
///////////////////////////////////////////////////////////////////////////
AliDevice* IceF2k::GetFitdefs()
{
// Provide pointer to the fit definitions
 return fFitdefs;
}
///////////////////////////////////////////////////////////////////////////
AliDevice* IceF2k::GetTrigdefs()
{
// Provide pointer to the trigger definitions
 return fTrigdefs;
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::Exec(Option_t* opt)
{
// Job to loop over the specified number of events and convert the 
// F2K data into the IceEvent structure.
// If maxevents<0 (default) all the entries of the input file
// will be processed.
// Every "printfreq" events a short event summary will be printed.
// The default value is printfreq=1.
// The output will be written on a standard output tree named "T".
//
// Notes :
// -------
// 1) This class is derived from AliJob, allowing a task based processing.
//    After the conversion of an F2K event into an IceEvent structure,
//    the processing of all available sub-tasks (if any) is invoked.
//    This provides an event-by-event (sub)task processing before the
//    final data structures are written out.
// 2) The main object in this job environment is an IceEvent* pointer.

 if (!fInfiles)
 {
  cout << " *IceF2k Exec* No data input file(s) specified." << endl;
  return;
 }

 Int_t ninfiles=fInfiles->GetEntries();
 if (!ninfiles)
 {
  cout << " *IceF2k Exec* No data input file(s) specified." << endl;
  return;
 }

 TTree* otree=0;
 if (fOutfile)
 {
  otree=new TTree("T","F2K Data converted to IceEvent structures");
  otree->SetDirectory(fOutfile);
 }

 IceEvent* evt=new IceEvent();
 evt->SetTrackCopy(1);
 evt->SetDevCopy(1);

 // Branch in the tree for the event structure
 if (otree) otree->Branch("IceEvent","IceEvent",&evt,fBsize,fSplit); 

 // Create the particle database and extend it with some F2000 specific definitions
 if (!fPdg) fPdg=new TDatabasePDG();
 Double_t me=fPdg->GetParticle(11)->Mass();
 fPdg->AddParticle("brems"   ,"brems"   ,0,1,0,0,"none",10001001,0,0);
 fPdg->AddParticle("deltae"  ,"deltae"  ,me,1,0,-3,"Lepton",10001002,0,0);
 fPdg->AddParticle("pairprod","pairprod",0,1,0,0,"none",10001003,0,0);
 fPdg->AddParticle("nucl_int","nucl_Int",0,1,0,0,"none",10001004,0,0);
 fPdg->AddParticle("mu_pair" ,"mu_pair" ,0,1,0,0,"none",10001005,0,0);
 fPdg->AddParticle("hadrons" ,"hadrons" ,0,1,0,0,"none",10001006,0,0);
 fPdg->AddParticle("fiberlaser","fiberlaser",0,1,0,0,"none",10002100,0,0);
 fPdg->AddParticle("n2laser"   ,"n2laser"   ,0,1,0,0,"none",10002101,0,0);
 fPdg->AddParticle("yaglaser"  ,"yaglaser"  ,0,1,0,0,"none",10002201,0,0);
 fPdg->AddParticle("z_primary","z_primary",0,1,0,0,"none",10003000,0,0);
 fPdg->AddParticle("a_primary","a_primary",0,1,0,0,"none",10003500,0,0);

 // Storage of the used parameters in the IceF2k device
 AliSignal params;
 params.SetNameTitle("IceF2k","IceF2k processor parameters");
 params.SetSlotName("Toffset",1);
 params.SetSlotName("Mctoffset",2);
 params.SetSlotName("Mctracks",3);

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
  cout << " F2K input file : " << inputfile.Data() << endl;
 }
 cout << " Maximum number of events to be processed : " << fMaxevt << endl;
 cout << " Print frequency : " << fPrintfreq << endl;
 if (fOutfile)
 {
  cout << " ROOT output file : " << fOutfile->GetName() << endl;
  cout << " Output characteristics : splitlevel = " << fSplit << " buffersize = " << fBsize << endl;
 }

 ListEnvironment();

 Int_t nevt=0;
 for (Int_t ifile=0; ifile<ninfiles; ifile++)
 {
  TObjString* sx=(TObjString*)fInfiles->At(ifile);
  if (!sx) continue;

  inputfile=sx->GetString(); 
  if (inputfile=="") continue;

  // Open the input file in the default ascii format (autodetection) for reading 
  fInput=rdmc_mcopen(inputfile.Data(),"r",RDMC_DEFAULT_ASCII_F);

  if (!fInput)
  {
   cout << " *IceF2k Exec* No input file found with name : " << inputfile.Data() << endl;
   continue;
  }

  // Initialise the event structure 
  rdmc_init_mevt(&fEvent);

  // Read the file header information
  rdmc_rarr(fInput,&fHeader);

  // Fill the database with geometry, calib. etc... parameters
  // for all the devices
  FillOMdbase();

  // Set the fit definitions according to the F2000 header info
  SetFitdefs();

  // Set the trigger definitions according to the F2000 header info
  SetTrigdefs();
 
  while (!rdmc_revt(fInput,&fHeader,&fEvent))
  {
   if (fMaxevt>-1 && nevt>=fMaxevt) break;

   // Reset the complete Event structure
   evt->Reset();

   evt->SetRunNumber(fEvent.nrun);
   evt->SetEventNumber(fEvent.enr);
   evt->SetMJD(fEvent.mjd,fEvent.secs,fEvent.nsecs);

   // Take trigger offset into account which might have been
   // introduced during the filtering process.
   // For simulated data this will be treated separately in PutMcTracks().
   fToffset=fEvent.t_offset;

   PutTrigger();

   PutMcTracks();

   PutRecoTracks();

   PutHits();

   // Enter the IceF2k processor parameters into the event structure 
   params.SetSignal(fToffset,1);
   params.SetSignal(fMctoffset,2);
   params.SetSignal(fMctracks,3);
   evt->AddDevice(params);

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
  }
  if (fMaxevt>-1 && nevt>=fMaxevt) break;
 }

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
void IceF2k::FillOMdbase()
{
// Fill the database with geometry, calib. etc... parameters 
// for all the devices.

 if (fHeader.nch<=0)
 {
  if (fOmdb)
  {
   delete fOmdb;
   fOmdb=0;
  }
  return;
 }

 Int_t adccal=fHeader.is_calib.adc;
 Int_t tdccal=fHeader.is_calib.tdc;
 Int_t totcal=fHeader.is_calib.tot;

 TF1 fadccal("fadccal","(x-[1])*[0]");
 TF1 fadcdecal("fadcdecal","(x/[0])+[1]");
 fadccal.SetParName(0,"BETA-ADC");
 fadccal.SetParName(1,"PED-ADC");
 fadcdecal.SetParName(0,"BETA-ADC");
 fadcdecal.SetParName(1,"PED-ADC");

 TF1 ftdccal("ftdccal","(x*[0])-[1]-([0]-1.)*32767.-[2]/sqrt([3])");
 TF1 ftdcdecal("ftdcdecal","(x+([0]-1.)*32767.+[1]+[2]/sqrt([3]))/[0]");
 ftdccal.SetParName(0,"BETA-TDC");
 ftdccal.SetParName(1,"T0");
 ftdccal.SetParName(2,"ALPHA-TDC");
 ftdccal.SetParName(3,"ADC-SLEW");
 ftdcdecal.SetParName(0,"BETA-TDC");
 ftdcdecal.SetParName(1,"T0");
 ftdcdecal.SetParName(2,"ALPHA-TDC");
 ftdcdecal.SetParName(3,"ADC-SLEW");

 TF1 ftotcal("ftotcal","x*[0]");
 TF1 ftotdecal("ftotdecal","x/[0]");
 ftotcal.SetParName(0,"BETA-TOT");
 ftotdecal.SetParName(0,"BETA-TOT");

 if (fOmdb)
 {
  fOmdb->Reset();
 }
 else
 {
  fOmdb=new AliObjMatrix();
  fOmdb->SetNameTitle("OMDBASE","The OM geometry, calib. etc... database");
  fOmdb->SetOwner();
 }

 IceAOM* dev=0;
 Double_t pos[3]={0,0,0};
 for (Int_t i=0; i<fHeader.nch; i++)
 {
  dev=new IceAOM();
  dev->SetUniqueID(i+1);
  // Slots to hold the various (de)calibration functions  
  dev->SetSlotName("ADC",1);
  dev->SetSlotName("LE",2);
  dev->SetSlotName("TOT",3);
  // Slots to hold hardware parameters
  dev->SetSlotName("TYPE",4);
  dev->SetSlotName("ORIENT",5);
  dev->SetSlotName("THRESH",6);
  dev->SetSlotName("SENSIT",7);

  pos[0]=fHeader.x[i];
  pos[1]=fHeader.y[i];
  pos[2]=fHeader.z[i];
  dev->SetPosition(pos,"car");

  fadccal.SetParameter(0,fHeader.cal[i].beta_a);
  fadccal.SetParameter(1,fHeader.cal[i].ped);
  fadcdecal.SetParameter(0,fHeader.cal[i].beta_a);
  if (!fHeader.cal[i].beta_a) fadcdecal.SetParameter(0,1);
  fadcdecal.SetParameter(1,fHeader.cal[i].ped);

  ftdccal.SetParameter(0,fHeader.cal[i].beta_t);
  ftdccal.SetParameter(1,fHeader.cal[i].t_0);
  ftdccal.SetParameter(2,fHeader.cal[i].alpha_t);
  ftdccal.SetParameter(3,1.e20);
  ftdcdecal.SetParameter(0,fHeader.cal[i].beta_t);
  if (!fHeader.cal[i].beta_t) ftdcdecal.SetParameter(0,1);
  ftdcdecal.SetParameter(1,fHeader.cal[i].t_0);
  ftdcdecal.SetParameter(2,fHeader.cal[i].alpha_t);
  ftdcdecal.SetParameter(3,1.e20);

  ftotcal.SetParameter(0,fHeader.cal[i].beta_tot);
  ftotdecal.SetParameter(0,fHeader.cal[i].beta_tot);
  if (!fHeader.cal[i].beta_tot) ftotdecal.SetParameter(0,1);

  if (adccal)
  {
   dev->SetDecalFunction(&fadcdecal,1);
  }
  else
  {
   dev->SetCalFunction(&fadccal,1);
  }

  if (tdccal)
  {
   dev->SetDecalFunction(&ftdcdecal,2);
  }
  else
  {
   dev->SetCalFunction(&ftdccal,2);
  }

  if (totcal)
  {
   dev->SetDecalFunction(&ftotdecal,3);
  }
  else
  {
   dev->SetCalFunction(&ftotcal,3);
  }

  dev->SetSignal(fHeader.type[i],4);
  dev->SetSignal((Float_t)fHeader.costh[i],5);
  dev->SetSignal(fHeader.thresh[i],6);
  dev->SetSignal(fHeader.sensit[i],7);

  fOmdb->EnterObject(i+1,1,dev);
 }
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetFitdefs()
{
// Obtain the names of the variables for each fit procedure from the
// f2000 header. Each different fit procedure is then stored as a separate
// "hit" of an AliDevice object and the various fit variables are stored
// as separate signal slots of the corresponding "hit".
// Via the GetFitdefs() memberfunction this AliDevice object can be
// retrieved and stored in the ROOT output file if wanted.
// The name of the object is FitDefinitions and the stored data can be
// inspected via the AliDevice::Data() memberfunction and looks as follows :
//
//  *AliDevice::Data* Id :0 Name : FitDefinitions
//    Position Vector in car coordinates : 0 0 0
//    Err. in car coordinates : 0 0 0
//  The following 8 hits are registered :
//  *AliSignal::Data* Id :0
//    Position Vector in car coordinates : 0 0 0
//    Err. in car coordinates : 0 0 0
//    Owned by device : AliDevice Name : FitDefinitions
//    Slot : 1 Signal value : 0 name : id
//    Slot : 2 Signal value : 0 name : rchi2
//    Slot : 3 Signal value : 0 name : prob
//    Slot : 4 Signal value : 0 name : sigth
//    Slot : 5 Signal value : 0 name : covmin
//    Slot : 6 Signal value : 0 name : covmax
//    Slot : 7 Signal value : 0 name : cutflag
//    Slot : 8 Signal value : 0 name : chi2
//  *AliSignal::Data* Id :1
//    Position Vector in car coordinates : 0 0 0
//    Err. in car coordinates : 0 0 0
//    Owned by device : AliDevice Name : FitDefinitions
//    Slot : 1 Signal value : 0 name : id
//    Slot : 2 Signal value : 0 name : rchi2
//    Slot : 3 Signal value : 0 name : prob
// etc....  
//
// This memberfunction is based on the original idea/code by Adam Bouchta.

 if (fHeader.n_fit<=0)
 {
  if (fFitdefs)
  {
   delete fFitdefs;
   fFitdefs=0;
  }
  return;
 }

 if (fFitdefs)
 {
  fFitdefs->Reset(1);
 }
 else
 {
  fFitdefs=new AliDevice();
 }

 fFitdefs->SetName("FitDefinitions");
 fFitdefs->SetHitCopy (1);

 AliSignal s;
 s.Reset();

 for (Int_t i=0; i<fHeader.n_fit; i++)
 {
  s.SetUniqueID(fHeader.def_fit[i].id);
  s.SetName(TString(fHeader.def_fit[i].tag));

  for (Int_t j=0; j<fHeader.def_fit[i].nwords; j++)
  {
   s.SetSlotName(TString(fHeader.def_fit[i].words[j]),j+1);
   s.SetSignal(0,j+1);
  }

  fFitdefs->AddHit(s);
  s.Reset(1);
 }
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetTrigdefs()
{
// Obtain the names of the variables for each trigger procedure from the
// f2000 header. Each different trigger procedure is then stored as a separate
// "hit" of an AliDevice object and the various trigger variables are stored
// as separate signal slots of the corresponding "hit".
// Via the GetFitdefs() memberfunction this AliDevice object can be
// retrieved and stored in the ROOT output file if wanted.
// The name of the object is TrigDefinitions and the stored data can be
// inspected via the AliDevice::Data() memberfunction and looks as follows :
//
//  *AliDevice::Data* Id : 0 Name : TrigDefinitions
//    Position Vector in car (rad) coordinates : 0 0 0
//    Err. in car (rad) coordinates : 0 0 0
//  The following 9 hits are registered : 
//  *AliSignal::Data* Id : 1 Name : main
//    Position Vector in car (rad) coordinates : 0 0 0
//    Err. in car (rad) coordinates : 0 0 0
//    Owned by device : AliDevice Id : 0 Name : TrigDefinitions
//     Slot : 1 Signal value : 0 name : trig_pulse_le
//     Slot : 2 Signal value : 0 name : trig_pulse_tot
//     Slot : 3 Signal value : 0 name : regi_flag
//   *AliSignal::Data* Id : 2 Name : amaa
//     Position Vector in car (rad) coordinates : 0 0 0
//     Err. in car (rad) coordinates : 0 0 0
//     Owned by device : AliDevice Id : 0 Name : TrigDefinitions
//     Slot : 1 Signal value : 0 name : trig_pulse_le
//     Slot : 2 Signal value : 0 name : trig_pulse_tot
//     Slot : 3 Signal value : 0 name : regi_flag
//   *AliSignal::Data* Id : 3 Name : amab10
//     Position Vector in car (rad) coordinates : 0 0 0
//     Err. in car (rad) coordinates : 0 0 0
//     Owned by device : AliDevice Id : 0 Name : TrigDefinitions
//     Slot : 1 Signal value : 0 name : trig_pulse_le
//     Slot : 2 Signal value : 0 name : trig_pulse_tot
//     Slot : 3 Signal value : 0 name : regi_flag
// etc....  

 if (fHeader.n_trigger<=0)
 {
  if (fTrigdefs)
  {
   delete fTrigdefs;
   fTrigdefs=0;
  }
  return;
 }

 if (fTrigdefs)
 {
  fTrigdefs->Reset(1);
 }
 else
 {
  fTrigdefs=new AliDevice();
 }

 fTrigdefs->SetName("TrigDefinitions");
 fTrigdefs->SetHitCopy (1);

 AliSignal s;
 s.Reset();

 for (Int_t i=0; i<fHeader.n_trigger; i++)
 {
  s.SetUniqueID(fHeader.def_trig[i].id);
  s.SetName(TString(fHeader.def_trig[i].tag));

  for (Int_t j=0; j<fHeader.def_trig[i].nwords; j++)
  {
   s.SetSlotName(TString(fHeader.def_trig[i].words[j]),j+1);
   s.SetSignal(0,j+1);
  }

  fTrigdefs->AddHit(s);
  s.Reset(1);
 }
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::PutMcTracks()
{
// Get the MC tracks from the F2000 file into the IcePack structure.
// Note : MC tracks are given negative track id's in the event structure.
// This memberfunction is based on the original code by Adam Bouchta.

 IceEvent* evt=(IceEvent*)GetMainObject();
 if (!evt || fEvent.ntrack<=0) return;

 // User defined trigger offset in case of simulated data.
 // The offset in the F2K file is meant to put the primary interaction
 // well ahead of the detector trigger.
 // See the introductory docs of this IceF2k class for further details.
 fToffset=fMctoffset;

 if (!fMctracks) return;

 // Loop over all the tracks and add them to the current event
 AliTrack t;
 Double_t vec[3];
 AliPosition r;
 Ali3Vector p;
 Int_t tid=0;
 Int_t idpdg=0;
 Int_t idf2k=0;
 for (Int_t i=0; i<fEvent.ntrack; i++)
 {
  t.Reset ();

  // Beginpoint of the track
  vec[0]=fEvent.gen[i].x;
  vec[1]=fEvent.gen[i].y;
  vec[2]=fEvent.gen[i].z;
  r.SetPosition(vec,"car");
  t.SetBeginPoint(r);

  // Endpoint of the track
  vec[0]+=fEvent.gen[i].length*fEvent.gen[i].px;
  vec[1]+=fEvent.gen[i].length*fEvent.gen[i].py;
  vec[2]+=fEvent.gen[i].length*fEvent.gen[i].pz;
  r.SetPosition(vec,"car");
  t.SetEndPoint(r);

  // Momentum in GeV/c
  vec[0]=fEvent.gen[i].e*fEvent.gen[i].px*1e-3;
  vec[1]=fEvent.gen[i].e*fEvent.gen[i].py*1e-3;
  vec[2]=fEvent.gen[i].e*fEvent.gen[i].pz*1e-3;
  p.SetVector (vec,"car");
  t.Set3Momentum(p);

  // MC tracks are indicated by negative track id's
  tid=fEvent.gen[i].tag;
  t.SetId(-abs(tid));

  idf2k=fEvent.gen[i].id;
  idpdg=0;
  if (idf2k>1000)
  {
   idpdg=idf2k+10000000;
  }
  else if (idf2k <= 48)
  {
   idpdg=fPdg->ConvertGeant3ToPdg(idf2k);
  }
  else
  {
   if (idf2k==201) idpdg=12;
   if (idf2k==202) idpdg=14;
   if (idf2k==203) idpdg=16;
   if (idf2k==204) idpdg=-12;
   if (idf2k==205) idpdg=-14;
   if (idf2k==206) idpdg=-16;
  }

  // Check for the user selected MC track storage
  if (fMctracks==1) // Store only muon and muon-neutrino tracks
  {
   if (abs(idpdg)!=13 && abs(idpdg)!=14) continue;
  }
  else if (fMctracks==2) // Store all lepton tracks
  {
   if (abs(idpdg)<11 || abs(idpdg)>16) continue;
  }

  t.SetParticleCode(idpdg);
  t.SetName(fPdg->GetParticle(idpdg)->GetName());
  t.SetTitle("MC track");
  t.SetMass(fPdg->GetParticle(idpdg)->Mass());
  t.SetCharge(fPdg->GetParticle(idpdg)->Charge()/3.);

  evt->AddTrack(t);
 }

 // Create the pointers to each particle's parent particle.
 Int_t txid=0;
 Int_t parid=0;
 for (Int_t itk=1; itk<=evt->GetNtracks (); itk++)
 {
  AliTrack* tx=evt->GetTrack(itk);

  if (!tx) continue;

  txid=tx->GetId();

  parid=-1;
  for (Int_t j=0; j<fEvent.ntrack; j++)
  {
   tid=fEvent.gen[j].tag;
   if (-abs(tid) == txid) parid=fEvent.gen[j].parent;
  }

  if (parid<0) continue;

  AliTrack* tpar=evt->GetIdTrack(-abs(parid));

  if (!tpar) continue;

  tx->SetParentTrack(tpar);
 }
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::PutRecoTracks()
{
// Get the reconstructed tracks from the F2000 file into the IcePack structure.
// Note : Reco tracks are given positive track id's in the event structure.
// This memberfunction is based on the original code by Adam Bouchta.

 IceEvent* evt=(IceEvent*)GetMainObject();
 if (!evt || fEvent.nfit<=0) return;

 // Loop over all the tracks and add them to the current event
 AliTrack t;
 Double_t vec[3];
 AliPosition r;
 Ali3Vector p;
 Int_t tid=0;
 Int_t idpdg=0;
 Int_t idf2k=0;
 for (Int_t i=0; i<fEvent.nfit; i++)
 {
  t.Reset ();

  // Beginpoint of the track
  vec[0]=fEvent.rec[i].x;
  vec[1]=fEvent.rec[i].y;
  vec[2]=fEvent.rec[i].z;
  r.SetPosition(vec,"car");
  t.SetBeginPoint(r);

  // Endpoint of the track
  vec[0]+=fEvent.rec[i].length*fEvent.rec[i].px;
  vec[1]+=fEvent.rec[i].length*fEvent.rec[i].py;
  vec[2]+=fEvent.rec[i].length*fEvent.rec[i].pz;
  r.SetPosition(vec,"car");
  t.SetEndPoint(r);

  // Momentum in GeV/c
  if (fEvent.rec[i].e > 0)
  {
   vec[0]=fEvent.rec[i].e*fEvent.rec[i].px*1e-3;
   vec[1]=fEvent.rec[i].e*fEvent.rec[i].py*1e-3;
   vec[2]=fEvent.rec[i].e*fEvent.rec[i].pz*1e-3;
  }
  else // Give the track a nominal momentum of 1 GeV/c
  {
   vec[0]=fEvent.rec[i].px;
   vec[1]=fEvent.rec[i].py;
   vec[2]=fEvent.rec[i].pz;
  }
  p.SetVector (vec,"car");
  t.Set3Momentum(p);

  // Use the fit number as track id
  tid=fEvent.rec[i].tag;
  t.SetId(abs(tid));

  idf2k=fEvent.rec[i].id;
  idpdg=0;
  if (idf2k>1000)
  {
   idpdg=idf2k+10000000;
  }
  else if (idf2k <= 48)
  {
   idpdg=fPdg->ConvertGeant3ToPdg(idf2k);
  }
  else
  {
   if (idf2k==201) idpdg=12;
   if (idf2k==202) idpdg=14;
   if (idf2k==203) idpdg=16;
   if (idf2k==204) idpdg=-12;
   if (idf2k==205) idpdg=-14;
   if (idf2k==206) idpdg=-16;
  }

  t.SetParticleCode(idpdg);
  t.SetNameTitle("Sieglinde","RECO track");
  t.SetMass(fPdg->GetParticle(idpdg)->Mass());
  t.SetCharge(fPdg->GetParticle(idpdg)->Charge()/3.);

  // Retrieve the various fit parameters for this track
  AliSignal* fitdata=fFitdefs->GetIdHit(i);
  for (Int_t jval=0; jval<fEvent.fresult[i].nval; jval++)
  {
   fitdata->SetSignal(fEvent.fresult[i].val[jval],jval+1);
  }

  // Store the various fit parameters for this track
  t.SetFitDetails(fitdata);

  // Store the various reco tracks as track hypotheses.
  // A copy of the first reco track is entered as a new track instance
  // into the event and all reco tracks (incl. the first one) are
  // stored as hypotheses linked to this new reco track.
  if (i==0)
  {
   evt->AddTrack(t);
   AliTrack* tx=evt->GetTrack(evt->GetNtracks());
   Int_t nrec=evt->GetNtracks(1);
   tx->SetId(nrec+1);
  }
  AliTrack* tx=evt->GetTrack(evt->GetNtracks());
  if (tx) tx->AddTrackHypothesis(t);
 }
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::PutHits()
{
// Get the hit and waveform info from the F2000 file into the IcePack structure.
// This memberfunction is based on the original code by Adam Bouchta.

 IceEvent* evt=(IceEvent*)GetMainObject();
 if (!evt) return;

 // Loop over all the hits and add them to the current event
 IceAOM om;
 AliSignal s;
 s.SetSlotName("ADC",1);
 s.SetSlotName("LE",2);
 s.SetSlotName("TOT",3);
 Int_t chan=0;
 Int_t maxchan=800;
 if (fOmdb) maxchan=fHeader.nch;
 IceAOM* omx=0;
 AliSignal* sx=0;
 Int_t tid=0;
 AliTrack* tx=0;
 Float_t adc=0;
 Float_t adcfirst=0; // Adc value of the first hit of an OM
 for (Int_t i=0; i<fEvent.nhits; i++)
 {
  chan=fEvent.h[i].ch+1;
  if (chan>maxchan) continue; // Channels 9001, 9002 etc are trigger channels

  // Get corresponding device from the current event structure  
  omx=(IceAOM*)evt->GetIdDevice(chan);
  if (!omx)
  {
   if (fOmdb)
   {
    omx=(IceAOM*)fOmdb->GetObject(chan,1);
    evt->AddDevice(omx);
   }
   else
   {
    om.Reset(1);
    om.SetUniqueID(chan);
    evt->AddDevice(om);
   }
   omx=(IceAOM*)evt->GetIdDevice(chan);
  }

  if (!omx) continue;

  adc=fEvent.h[i].amp;

  // Multiple hits in the same OM with the same ADC value
  // are indicated by "*" in the F2K file.
  // This corresponds to a value of -2 in the data structure.
  if (int(adc) == -2)
  {
   adc=adcfirst;
  }
  else
  {
   adcfirst=adc;
  }
  s.Reset();
  s.SetUniqueID(fEvent.h[i].id);
  s.SetSignal(adc,1);
  s.SetSignal((fEvent.h[i].t-fToffset),2);
  s.SetSignal(fEvent.h[i].tot,3);

  omx->AddHit(s);

  sx=omx->GetHit(omx->GetNhits());
  if (!sx) continue;

  // ADC dependent TDC (de)calibration function for this hit
  TF1* fcal=omx->GetCalFunction("LE");
  TF1* fdecal=omx->GetDecalFunction("LE");
  if (fcal) sx->SetCalFunction(fcal,2);
  if (fdecal) sx->SetDecalFunction(fdecal,2);
  fcal=sx->GetCalFunction(2);
  fdecal=sx->GetDecalFunction(2);
  adc=sx->GetSignal(1,-4);
  if (adc>0)
  {
   if (fcal) fcal->SetParameter(3,adc);
   if (fdecal) fdecal->SetParameter(3,adc);
  }
  else
  {
   if (fcal) fcal->SetParameter(3,1.e20);
   if (fdecal) fdecal->SetParameter(3,1.e20);
  }

  // Bi-directional link between this hit and the track that caused the ADC value.
  // This F2K info is probably only present for MC tracks.
  tid=fEvent.h[i].ma;
  if (tid > 0)
  {
   tx=evt->GetIdTrack(tid); // Reco tracks
   if (!tx) tx=evt->GetIdTrack(-tid); // MC tracks
   if (tx) sx->AddTrack(*tx);
  }
  else
  {
   if (tid == -2) sx->SetNameTitle("N","Noise");
   if (tid == -3) sx->SetNameTitle("A","Afterpulse");
  }
 }

 // Loop over all the waveforms and add the histo(s) to the corresponding OM's
 TH1F histo;
 Int_t nbins=0;
 Float_t xlow=0;
 Float_t xup=0;
 TString hname;
 for (Int_t iwf=0; iwf<fEvent.nwf; iwf++)
 {
  chan=fEvent.wf[iwf].om;
  if (chan<=0 || chan>maxchan) continue; // Skip trigger channels

  // Get corresponding device from the current event structure  
  omx=(IceAOM*)evt->GetIdDevice(chan);
  if (!omx)
  {
   if (fOmdb)
   {
    omx=(IceAOM*)fOmdb->GetObject(chan,1);
    evt->AddDevice(omx);
   }
   else
   {
    om.Reset(1);
    om.SetUniqueID(chan);
    evt->AddDevice(om);
   }
   omx=(IceAOM*)evt->GetIdDevice(chan);
  }

  if (!omx) continue;

  hname="BASELINE-WF";
  hname+=omx->GetNwaveforms()+1;
  omx->AddNamedSlot(hname);
  omx->SetSignal(fEvent.wf[iwf].baseline,hname);

  // Fill the waveform histogram
  hname="OM";
  hname+=chan;
  hname+="-WF";
  hname+=omx->GetNwaveforms()+1;

  histo.Reset();
  histo.SetName(hname.Data());
  nbins=fEvent.wf[iwf].ndigi;
  xlow=fEvent.wf[iwf].t_start;
  xup=xlow+float(nbins)*fEvent.wf[iwf].t_bin;
  histo.SetBins(nbins,xlow,xup);

  for (Int_t jbin=1; jbin<=fEvent.wf[iwf].ndigi; jbin++)
  {
   histo.SetBinContent(jbin,fEvent.wf[iwf].baseline-fEvent.wf[iwf].digi[jbin-1]);
  }

  omx->SetWaveform(&histo,omx->GetNwaveforms()+1);
 }

 // Set bi-directional links between hits and reco track hypotheses.
 // Note : Reco tracks are recognised by their positive id.
 Int_t hid=0;
 TObjArray* rectracks=evt->GetTracks(1);
 for (Int_t jtk=0; jtk<rectracks->GetEntries(); jtk++)
 {
  tx=(AliTrack*)rectracks->At(jtk);
  if (!tx) continue;
  
  for (Int_t jhyp=1; jhyp<=tx->GetNhypotheses(); jhyp++)
  {
   AliTrack* hypx=tx->GetTrackHypothesis(jhyp);
   if (!hypx) continue;

   // Loop over all combinations of F2K fits and used OM hits
   for (Int_t k=0; k<fEvent.nfit_uses; k++)
   {
    if (fEvent.fit_uses[k].useid != hypx->GetId()) continue;
    hid=fEvent.fit_uses[k].hitid;
    sx=evt->GetIdHit(hid,"IceAOM");
    if (sx) sx->AddTrack(*hypx);
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::PutTrigger()
{
// Get the trigger info from the F2000 file into the IcePack structure.

 if (!fTrigdefs) return;

 IceEvent* evt=(IceEvent*)GetMainObject();
 if (!evt || fEvent.ntrig<=0) return;

 AliDevice trig;
 trig.SetNameTitle("Trigger","Amanda/IceCube event triggers");
 AliSignal s;
 TString trigname;
 TString slotname;
 Int_t id=0;
 Int_t nval=0;
 for (Int_t i=0; i<fEvent.ntrig; i++)
 {
  id=fEvent.ptrig[i].id;
  nval=fEvent.ptrig[i].nval;
  if (!nval) continue;
  AliSignal* tdef=fTrigdefs->GetIdHit(id+1);
  if (!tdef) continue;
  trigname=tdef->GetName();
  s.Reset(1);
  s.SetName(trigname);
  s.SetUniqueID(id+1);
  for (Int_t jval=0; jval<fEvent.ptrig[i].nval; jval++)
  {
   slotname=tdef->GetSlotName(jval+1);
   s.SetSlotName(slotname,jval+1);
   s.SetSignal(fEvent.ptrig[i].val[jval],jval+1);
  }
  trig.AddHit(s);
 }

 // Store the trigger data into the IceEvent structure
 evt->AddDevice(trig);
}
///////////////////////////////////////////////////////////////////////////
