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
// Class IceRootx
// Conversion of simple Root data into IceEvent data structures.
// This class reads data from the simple Root files as output by the
// original version of Martijn Duvoort's Walnut analyser.
// This class is only retained for backward compatibility. For output from
// newer versions of the Walnut analyser, use IceRoot.
// An indication of the active DAQ system is available in the IceEvent structure
// via a device named "Daq". Here the various daq systems (TWR, Muon, ...)
// from which the actual hits (ADC, LE, TOT) eventually will be composed
// are indicated as "signals" of the device itself. 
// This class is derived from AliJob providing a task-based processing
// structure on an event-by-event basis.
// The main object in the job environment is an IceEvent* pointer.
// In case the user has provided sub-tasks, these will be executed
// on an event-by-event basis after the IceEvent structure has been filled
// with the simple Root data and before the final structures are written out.
// Note that the data structures are only written out if an outputfile has
// been specified via the SetOutputFile memberfunction.
// In case no outputfile has been specified, this class provides a facility
// to investigate/analyse simple Root data using the Ralice/IcePack analysis tools.
//
// Usage example :
// ---------------
//
// gSystem->Load("ralice");
// gSystem->Load("icepack");
// gSystem->Load("iceconvert");
//
// IceRootx q("IceRootx","Simple Root data to IcePack data structure conversion");
//
// // Limit the number of entries for testing
// q.SetMaxEvents(10);
//
// // Print frequency to produce a short summary print every printfreq events
// q.SetPrintFreq(1);
//
// // The simple Root data input filename(s)
// q.AddInputFile("test-i3.root");
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
//--- Author: Garmt de Vries-Uiterweerd 13-Mar-2007 Utrecht University
//- Modified: GdV $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceRootx.h"
#include "Riostream.h"

ClassImp(IceRootx) // Class implementation to enable ROOT I/O

IceRootx::IceRootx(const char* name,const char* title) : AliJob(name,title)
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
IceRootx::~IceRootx()
{
// Default destructor.

 if (fInfiles)
 {
  delete fInfiles;
  fInfiles=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceRootx::SetMaxEvents(Int_t n)
{
// Set the maximum number of events to be processed.
// n=-1 implies processing of the complete input file, which is the default
// initialisation in the constructor.
 fMaxevt=n;
}
///////////////////////////////////////////////////////////////////////////
void IceRootx::SetPrintFreq(Int_t f)
{
// Set the printfrequency to produce info every f events.
// f=1 is the default initialisation in the constructor.
 if (f>=0) fPrintfreq=f;
}
///////////////////////////////////////////////////////////////////////////
void IceRootx::SetSplitLevel(Int_t split)
{
// Set the split level for the ROOT data file.
// split=0 is the default initialisation in the constructor.
 if (split>=0) fSplit=split;
}
///////////////////////////////////////////////////////////////////////////
void IceRootx::SetBufferSize(Int_t bsize)
{
// Set the buffer size for the ROOT data file.
// bsize=32000 is the default initialisation in the constructor.
 if (bsize>=0) fBsize=bsize;
}
///////////////////////////////////////////////////////////////////////////
void IceRootx::AddInputFile(TString name)
{
// Add the name of this simple Root data input file to the list to be processed.

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
void IceRootx::SetOutputFile(TFile* ofile)
{
// Set the output file for the ROOT data.
 if (fOutfile) delete fOutfile;
 fOutfile=ofile;
}
///////////////////////////////////////////////////////////////////////////
void IceRootx::SetOutputFile(TString name)
{
// Create the output file for the ROOT data.
 if (fOutfile) delete fOutfile;
 fOutfile=new TFile(name.Data(),"RECREATE","Simple Root data in IceEvent structure");
}
///////////////////////////////////////////////////////////////////////////
TFile* IceRootx::GetOutputFile()
{
// Provide pointer to the ROOT output file.
 return fOutfile;
}
///////////////////////////////////////////////////////////////////////////
void IceRootx::Exec(Option_t* opt)
{
// Job to loop over the specified number of events and convert the 
// simple Root data into the IceEvent structure.
// If maxevents<0 (default) all the entries of the input file
// will be processed.
// Every "printfreq" events a short event summary will be printed.
// The default value is printfreq=1.
// The output will be written on a standard output tree named "T".
//
// Notes :
// -------
// 1) This class is derived from AliJob, allowing a task based processing.
//    After the conversion of a simple Root data event into an IceEvent structure,
//    the processing of all available sub-tasks (if any) is invoked.
//    This provides an event-by-event (sub)task processing before the
//    final data structures are written out.
// 2) The main object in this job environment is an IceEvent* pointer.

 if (!fInfiles)
 {
  cout << " *IceRootx Exec* No data input file(s) specified." << endl;
  return;
 }

 Int_t ninfiles=fInfiles->GetEntries();
 if (!ninfiles)
 {
  cout << " *IceRootx Exec* No data input file(s) specified." << endl;
  return;
 }

 // Create output tree if necessary
 TTree* otree=0;
 if (fOutfile)
 {
  otree=new TTree("T","Simple Root data converted to IceEvent structures");
  otree->SetDirectory(fOutfile);
 }

 // Create IceEvent structure
 IceEvent* evt=new IceEvent();
 evt->SetTrackCopy(1);
 evt->SetDevCopy(1);

 // Branch in the tree for the event structure
 if (otree) otree->Branch("IceEvent","IceEvent",&evt,fBsize,fSplit); 

 // Initialise the job working environment
 SetMainObject(evt);

 // Some output for the user's convenience
 TString inputfile;
 cout << " ***" << endl;
 cout << " *** Start processing of job " << GetName() << " ***" << endl;
 cout << " ***" << endl;
 for (Int_t i=0; i<ninfiles; i++)
 {
  TObjString* sx=(TObjString*)fInfiles->At(i);
  if (!sx) continue;
  inputfile=sx->GetString(); 
  cout << " Simple Root data input file : " << inputfile.Data() << endl;
 }
 cout << " Maximum number of events to be processed : " << fMaxevt << endl;
 cout << " Print frequency : " << fPrintfreq << endl;
 if (fOutfile)
 {
  cout << " ROOT output file : " << fOutfile->GetName() << endl;
  cout << " Output characteristics : splitlevel = " << fSplit << " buffersize = " << fBsize << endl;
 }

 ListEnvironment();

 // Set DAQ device info
 AliDevice daq;
 daq.SetName("Daq");
 daq.SetSlotName("TWR",1);
 daq.SetSignal(1,1);

 // Trigger device
 AliDevice trig;
 trig.SetNameTitle("Trigger","Amanda/IceCube event triggers");
 AliSignal s;

 // Some variables
 Double_t triggertime=0;
 Int_t    eventtimemjd=0;
 Int_t    eventtimemjds=0;
 Int_t    eventtimemjdns=0;
 Int_t    eventid=0;
 Int_t    runid=0;
 Int_t    String=0;
 Int_t    OM=0;
 Float_t  baseline=0;
 Double_t binsize=10;
 Short_t  type=0;
 Int_t    numberbins=0;
 Float_t  starttime=0;
 Int_t    ntracks=0;

 Int_t nevt=0;
 Int_t lastevent=0;
 Int_t omid=0;

 TString hname;
 TH1F histo;
 IceAOM om;
 IceAOM* omx=0;
 AliTrack t;
 Double_t vec[3];
 AliPosition r;
 Ali3Vector p;

 Float_t pi=acos(-1.);

 Int_t firstomonstring[20];
 firstomonstring[0]=0;
 firstomonstring[1]=1;
 firstomonstring[2]=21;
 firstomonstring[3]=41;
 firstomonstring[4]=61;
 firstomonstring[5]=87;
 firstomonstring[6]=123;
 firstomonstring[7]=159;
 firstomonstring[8]=195;
 firstomonstring[9]=231;
 firstomonstring[10]=267;
 firstomonstring[11]=303;
 firstomonstring[12]=345;
 firstomonstring[13]=387;
 firstomonstring[14]=429;
 firstomonstring[15]=471;
 firstomonstring[16]=513;
 firstomonstring[17]=555;
 firstomonstring[18]=597;
 firstomonstring[19]=639;

 TFile* input=0;

 for (Int_t ifile=0; ifile<ninfiles; ifile++)
 {
  TObjString* sx=(TObjString*)fInfiles->At(ifile);
  if (!sx) continue;

  inputfile=sx->GetString(); 
  if (inputfile=="") continue;

  // Open the simple Root data input file
  input=TFile::Open(inputfile);

  if (!input)
  {
   cout << " *IceRootx Exec* No input file found with name : " << inputfile.Data() << endl;
   continue;
  }

  // Get simple Root tree
  fTree=(TTree*)input->Get("T");
  fTree->SetBranchAddress("triggertime",&triggertime);
  fTree->SetBranchAddress("eventtimemjd",&eventtimemjd);
  fTree->SetBranchAddress("eventtimemjds",&eventtimemjds);
  fTree->SetBranchAddress("eventtimemjdns",&eventtimemjdns);
  fTree->SetBranchAddress("eventid",&eventid);
  fTree->SetBranchAddress("runid",&runid);
  fTree->SetBranchAddress("String",&String);
  fTree->SetBranchAddress("OM",&OM);
  fTree->SetBranchAddress("baseline",&baseline);
  fTree->SetBranchAddress("type",&type);
  fTree->SetBranchAddress("numberbins",&numberbins);
  fTree->SetBranchAddress("starttime",&starttime);
  fTree->SetBranchAddress("ntracks",&ntracks);
  if(fTree->GetBranch("binsize")) fTree->SetBranchAddress("binsize",&binsize);

  Int_t nmaxbins=(Int_t)fTree->GetLeaf("numberbins")->GetMaximum();
  Float_t* wvform=new Float_t[nmaxbins+1];
  fTree->SetBranchAddress("wvform",wvform);

  Int_t nmaxtracks=(Int_t)fTree->GetLeaf("ntracks")->GetMaximum();
  Float_t* trackx=new Float_t[nmaxtracks];
  Float_t* tracky=new Float_t[nmaxtracks];
  Float_t* trackz=new Float_t[nmaxtracks];
  Double_t* trackenergy=new Double_t[nmaxtracks];
  Float_t* trackzenith=new Float_t[nmaxtracks];
  Float_t* trackazimuth=new Float_t[nmaxtracks];
  Int_t* tracktype=new Int_t[nmaxtracks];
  fTree->SetBranchAddress("trackx",trackx);
  fTree->SetBranchAddress("tracky",tracky);
  fTree->SetBranchAddress("trackz",trackz);
  fTree->SetBranchAddress("trackenergy",trackenergy);
  fTree->SetBranchAddress("trackzenith",trackzenith);
  fTree->SetBranchAddress("trackazimuth",trackazimuth);
  fTree->SetBranchAddress("tracktype",tracktype);

  // Prepare for loop over entries
  lastevent=0;

  // Loop over waveforms in tree
  for(Int_t ientry=0; ientry<fTree->GetEntries(); ientry++)
  {
   fTree->GetEntry(ientry);
   
   // If new event
   if(eventid!=lastevent)
   {
    // Write old event to tree (if it contains data)
    if(evt->GetDevices("IceAOM"))
    {
     // Invoke all available sub-tasks (if any) and write event to tree
     CleanTasks();
     ExecuteTasks(opt);
     if(otree) otree->Fill();
     if (fPrintfreq) { if (!(nevt%fPrintfreq)) evt->HeaderData(); }
     // Update event counter
     nevt++;
     if (fMaxevt>-1 && nevt>=fMaxevt) 
     {
      evt->Reset();
      break;
     }
    }

    // Start new event
    evt->Reset();
    evt->SetRunNumber(runid);
    evt->SetEventNumber(eventid);
    evt->SetMJD(eventtimemjd,eventtimemjds,eventtimemjdns,0);
    evt->AddDevice(daq);

    // Store trigger information
    trig.Reset(1);
    s.Reset(1);
    s.SetName("main");
    s.SetTitle("First trigger for TWR time reference");
    s.SetUniqueID(12);
    s.SetSlotName("trig_pulse_le",1);
    s.SetSignal(triggertime,1);
    trig.AddHit(s);
    //// TODO: store other triggers if available
    evt->AddDevice(trig);

    // Loop over all the tracks and add them to the current event
    Int_t nrecotracks=0;
    Int_t nmctracks=0;
    for (Int_t itrack=0; itrack<ntracks; itrack++)
    {
     t.Reset();

     // Beginpoint of the track
     vec[0]=trackx[itrack];
     vec[1]=tracky[itrack];
     vec[2]=trackz[itrack];
     r.SetPosition(vec,"car");
     t.SetBeginPoint(r);

     // Momentum in GeV/c
     vec[0]=trackenergy[itrack];
     vec[1]=pi-trackzenith[itrack];
     vec[2]=trackazimuth[itrack]+pi;
     if(vec[2]>=2*pi) vec[2]-=2*pi;
     p.SetVector(vec,"sph");
     t.Set3Momentum(p);

     // Check for unreconstructed tracks (NaN values) and reset track if necessary
     if(trackx[itrack]!=trackx[itrack] || tracky[itrack]!=tracky[itrack] || trackz[itrack]!=trackz[itrack] ||
        trackzenith[itrack]!=trackzenith[itrack] || trackazimuth[itrack]!=trackazimuth[itrack]){
      t.Reset();
     }

     // Track ID and name
     // Monte Carlo track
     if(tracktype[itrack]==7){
      nmctracks++;
      t.SetId(-nmctracks);
     }
     // Reco track
     else {
      nrecotracks++;
      t.SetId(nrecotracks);
     }
     if(tracktype[itrack]==1) { t.SetName("I3CFIRST"); t.SetTitle("IceTray CFIRST track"); }
     else if(tracktype[itrack]==2) { t.SetName("I3direct"); t.SetTitle("IceTray direct walk track"); }
     else if(tracktype[itrack]==3) { t.SetName("I3Jams_cluster"); t.SetTitle("IceTray Jams cluster"); }
     else if(tracktype[itrack]==4) { t.SetName("I3Jams_qual"); t.SetTitle("IceTray Jams quality"); }
     else if(tracktype[itrack]==5) { t.SetName("I3TOIFit"); t.SetTitle("IceTray TOI fit"); }
     else if(tracktype[itrack]==6) { t.SetName("I3line-direct"); t.SetTitle("IceTray line direct walk"); }
     else if(tracktype[itrack]==7) { t.SetName("I3MCTrack"); t.SetTitle("I3MCTrack"); }
     else { t.SetName("Unknown"); t.SetTitle("IceTray track of unknown type"); }

     // Add track to event
     evt->AddTrack(t);
    }

    // Remember event nr
    lastevent=eventid;
   }

   // Get OM from event strcture, or create and add it
   omid=firstomonstring[-String]+OM-1;
   if(omid==681) continue;  // Skip OM 681, which should never give data: avoid all risk of confusion
   omx=(IceAOM*)evt->GetIdDevice(omid);
   if (!omx)
   {
    om.Reset(1);
    om.SetUniqueID(omid);
    evt->AddDevice(om);
    omx=(IceAOM*)evt->GetIdDevice(omid);
   }
   if (!omx) continue;

   // Store baseline info
   hname="BASELINE-WF";
   hname+=omx->GetNwaveforms()+1;
   omx->AddNamedSlot(hname);
   omx->SetSignal(baseline,hname);

   // Store readout type
   omx->AddNamedSlot("READOUT");
   if(type==0) omx->SetSignal(1,"READOUT");       // Electrical
   else if(type==1) omx->SetSignal(2,"READOUT");  // Optical
   else if(type==2) omx->SetSignal(3,"READOUT");  // Digital 
   else omx->SetSignal(0,"READOUT");              // Unknown

   // Fill the waveform histogram with this fragment
   hname="OM";
   hname+=omid;
   hname+="-WF";
   hname+=omx->GetNwaveforms()+1;

   histo.Reset();
   histo.SetName(hname.Data());
   histo.SetBins(numberbins,triggertime+starttime,triggertime+starttime+numberbins*binsize);

   for (Int_t jbin=1; jbin<=numberbins; jbin++)
   {
    histo.SetBinContent(jbin,baseline-wvform[jbin-1]);
   }

   omx->SetWaveform(&histo,omx->GetNwaveforms()+1);

  }

  // Write last event to tree
  if(evt->GetDevices("IceAOM"))
  {
   CleanTasks();
   ExecuteTasks(opt);
   if (otree) otree->Fill();
   if (fPrintfreq) { if (!(nevt%fPrintfreq)) evt->HeaderData(); }
   // Update event counter
   nevt++;
   if (fMaxevt>-1 && nevt>=fMaxevt) break;
  }

  // Reset event
  evt->Reset();

  // Close input file
  input->Close();

  // Clean up
  delete[] wvform;
  delete[] trackx;
  delete[] tracky;
  delete[] trackz;
  delete[] trackenergy;
  delete[] trackzenith;
  delete[] trackazimuth;
  delete[] tracktype;

  // Stop looping over input files if max. nr. of events is reached
  if (fMaxevt>-1 && nevt>=fMaxevt)
  {
   break;
  }

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
