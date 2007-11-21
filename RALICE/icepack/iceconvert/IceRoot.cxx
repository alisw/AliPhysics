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
// Class IceRoot
// Conversion of simple Root data into IceEvent data structures.
// This class reads data from the simple Root files as output by Martijn
// Duvoort's Walnut analyser.
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
// IceRoot q("IceRoot","Simple Root data to IcePack data structure conversion");
//
// // Specify the relevant calibration file
// q.SetCalibFile("calib2006.root");
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
//- Modified: GdVU $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceRoot.h"
#include "Riostream.h"

ClassImp(IceRoot) // Class implementation to enable ROOT I/O

IceRoot::IceRoot(const char* name,const char* title) : AliJob(name,title)
{
// Default constructor.
// By default maxevent=-1, split=0, bsize=32000, printfreq=1.

 fSplit=0;
 fBsize=32000;
 fMaxevt=-1;
 fPrintfreq=1;
 fInfiles=0;
 fOutfile=0;
 fCalfile=0;
 fJEBTDaq=0;
}
///////////////////////////////////////////////////////////////////////////
IceRoot::~IceRoot()
{
// Default destructor.

 if (fInfiles)
 {
  delete fInfiles;
  fInfiles=0;
 }

 if (fCalfile)
 {
  delete fCalfile;
  fCalfile=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceRoot::SetMaxEvents(Int_t n)
{
// Set the maximum number of events to be processed.
// n=-1 implies processing of the complete input file, which is the default
// initialisation in the constructor.
 fMaxevt=n;
}
///////////////////////////////////////////////////////////////////////////
void IceRoot::SetPrintFreq(Int_t f)
{
// Set the printfrequency to produce info every f events.
// f=1 is the default initialisation in the constructor.
 if (f>=0) fPrintfreq=f;
}
///////////////////////////////////////////////////////////////////////////
void IceRoot::SetSplitLevel(Int_t split)
{
// Set the split level for the ROOT data file.
// split=0 is the default initialisation in the constructor.
 if (split>=0) fSplit=split;
}
///////////////////////////////////////////////////////////////////////////
void IceRoot::SetBufferSize(Int_t bsize)
{
// Set the buffer size for the ROOT data file.
// bsize=32000 is the default initialisation in the constructor.
 if (bsize>=0) fBsize=bsize;
}
///////////////////////////////////////////////////////////////////////////
void IceRoot::SetInputFile(TString name)
{
// Set the simple Root data input file.
 if (!fInfiles)
 {
  fInfiles=new TObjArray();
  fInfiles->SetOwner();
 }
 fInfiles->Clear();

 TObjString* s=new TObjString();
 s->SetString(name);
 fInfiles->Add(s);
}
///////////////////////////////////////////////////////////////////////////
void IceRoot::AddInputFile(TString name)
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
void IceRoot::SetOutputFile(TFile* ofile)
{
// Set the output file for the ROOT data.
 if (fOutfile) delete fOutfile;
 fOutfile=ofile;
}
///////////////////////////////////////////////////////////////////////////
void IceRoot::SetOutputFile(TString name)
{
// Create the output file for the ROOT data.
 if (fOutfile) delete fOutfile;
 fOutfile=new TFile(name.Data(),"RECREATE","Simple Root data in IceEvent structure");
}
///////////////////////////////////////////////////////////////////////////
TFile* IceRoot::GetOutputFile()
{
// Provide pointer to the ROOT output file.
 return fOutfile;
}
///////////////////////////////////////////////////////////////////////////
void IceRoot::SetCalibFile(TString name)
{
// Set the calibration ROOT file as created with IceCal2Root.
// Note: this will overrule a previously attached database. 
// In case no calibration file is specified, TWR waveforms cannot be processed.
 if (fCalfile) delete fCalfile;
 fCalfile=new TFile(name.Data());

 if(fCalfile){
  fJEBTDaq=(AliObjMatrix*)fCalfile->Get("JEBTDaq-OMDBASE");
 }
 if(!fJEBTDaq){
  cout << "*IceRoot* Warning: no calibration available for TWR. TWR waveforms cannot be processed." << endl;
 }

}
///////////////////////////////////////////////////////////////////////////
void IceRoot::SetOMdbase(AliObjMatrix* omdb)
{
// Set the calibration database as created with IceCal2Root.
// Note: this will overrule a previously attached database. 
// In case no calibration database is specified, TWR waveforms cannot be processed.
 fJEBTDaq=omdb;
 if(!fJEBTDaq){
  cout << "*IceRoot* Warning: no calibration available for TWR. TWR waveforms cannot be processed." << endl;
 }

}
///////////////////////////////////////////////////////////////////////////
void IceRoot::Exec(Option_t* opt)
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
  cout << " *IceRoot Exec* No data input file(s) specified." << endl;
  return;
 }

 Int_t ninfiles=fInfiles->GetEntries();
 if (!ninfiles)
 {
  cout << " *IceRoot Exec* No data input file(s) specified." << endl;
  return;
 }

 // Warning in case no calibration DB is specified
 if(!fJEBTDaq){
  cout << "*IceRoot* Warning: no calibration available for TWR. TWR waveforms cannot be processed." << endl;
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
 daq.SetNameTitle("Daq","Daq system used");

 // Trigger device
 AliDevice trig;
 trig.SetNameTitle("Trigger","Amanda/IceCube event triggers");
 AliSignal s;
 Char_t trigname[100];

 // Variables in the simple ROOT tree
 Int_t     pretrig=0;
 Int_t     trignr=0;
 Double_t* trigtime=0;
 Double_t* triglength=0;
 Char_t    trigsourceID[100][100];
 Char_t    trigtypeID[100][100];
 Char_t    trigsubtypeID[100][100];
 Int_t     eventtimemjd=0;
 Int_t     eventtimemjds=0;
 Int_t     eventtimemjdns=0;
 Int_t     eventid=0;
 Int_t     runid=0;
 Int_t     String=0;
 Int_t     OM=0;
 Short_t   waveformtype=0; // 1 = TWR, 2 = non-merged ATWD (not supported), 3 = InIce merged ATWD, 4 = InIce FADC, 5 = IceTop merged ATWD, 6 = IceTop FADC
 Short_t   twreventnr=0;
 UInt_t    baseline=0;
 Int_t     numberbinstwr=0;
 Int_t     startbin=0;
 UInt_t*   twrwaveform=0;
 Int_t     numberbinswaveform=0;
 Double_t  starttimewaveform=0;
 Double_t  binwidth=0;
 Short_t   source=0; // 0 = ATWD, 10 = FADC, 20 = TWR elec, 30 = TWR opt, 40 = etc
 Double_t* waveform=0;
 Int_t     ntracks=0;
 Double_t* trackx=0;
 Double_t* tracky=0;
 Double_t* trackz=0;
 Double_t* trackzenith=0;
 Double_t* trackazimuth=0;
 Double_t* trackenergy=0;
 Char_t    tracktype[100][100];

 // Some other variables
 Int_t nevt=0;
 Int_t lastevent=0;
 Int_t omid=0;
 Double_t starttime=0;
 Int_t extstop=0;
 Float_t twrtime=10240; // Number of bins in TWR
 Int_t twrbuffer=0;

 TString hname;
 TH1F histo;
 IceAOM aom;
 IceIDOM idom;
 IceTDOM tdom;
 IceGOM* om=0;
 IceGOM* calom=0;
 AliTrack t;
 Double_t vec[3];
 AliPosition r;
 Ali3Vector p;
 AliSample sample;

 Float_t pi=acos(-1.);

 TFile* input=0;

 Double_t* twreventtimes=0;
 Int_t nrtwrevents=0;
 Int_t* twreventorder=0;

 // Get global TWR time offsets from DB
 Float_t globalt0=0;
 Float_t twri3offset=0;
 if(fJEBTDaq){
  for(Int_t row=1; row<=fJEBTDaq->GetMaxRow(); row++){
   calom=(IceGOM*)fJEBTDaq->GetObject(row,1);
   if(calom){
    globalt0=calom->GetSignal("GLOBALT0");
    twri3offset=calom->GetSignal("TWRI3OFFSET");
    break;
   }
  }
 }

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
   cout << " *IceRoot Exec* No input file found with name : " << inputfile.Data() << endl;
   continue;
  }

  // Get simple Root tree
  fTree=(TTree*)input->Get("T");

  fTree->SetBranchAddress("pretrig",&pretrig);
  fTree->SetBranchAddress("trignr",&trignr);
  fTree->SetBranchAddress("eventtimemjd",&eventtimemjd);
  fTree->SetBranchAddress("eventtimemjds",&eventtimemjds);
  fTree->SetBranchAddress("eventtimemjdns",&eventtimemjdns);
  fTree->SetBranchAddress("eventid",&eventid);
  fTree->SetBranchAddress("runid",&runid);
  fTree->SetBranchAddress("String",&String);
  fTree->SetBranchAddress("OM",&OM);
  fTree->SetBranchAddress("waveformtype",&waveformtype);
  fTree->SetBranchAddress("twreventnr",&twreventnr);
  fTree->SetBranchAddress("baseline",&baseline);
  fTree->SetBranchAddress("numberbinstwr",&numberbinstwr);
  fTree->SetBranchAddress("startbin",&startbin);
  fTree->SetBranchAddress("numberbinswaveform",&numberbinswaveform);
  fTree->SetBranchAddress("starttimewaveform",&starttimewaveform);
  if(fTree->GetBranch("binwidth")) fTree->SetBranchAddress("binwidth",&binwidth);
  fTree->SetBranchAddress("source",&source);
  fTree->SetBranchAddress("ntracks",&ntracks);

  Int_t ntrig=(Int_t)fTree->GetLeaf("trignr")->GetMaximum();
  trigtime=new Double_t[ntrig+1];
  fTree->SetBranchAddress("trigtime",trigtime);
  triglength=new Double_t[ntrig+1];
  fTree->SetBranchAddress("triglength",triglength);
  fTree->SetBranchAddress("trigsourceID",trigsourceID);
  fTree->SetBranchAddress("trigtypeID",trigtypeID);
  fTree->SetBranchAddress("trigsubtypeID",trigsubtypeID);

  Int_t nmaxbinstwr=(Int_t)fTree->GetLeaf("numberbinstwr")->GetMaximum();
  twrwaveform=new UInt_t[nmaxbinstwr+1];
  fTree->SetBranchAddress("twrwaveform",twrwaveform);

  Int_t nmaxbinswaveform=(Int_t)fTree->GetLeaf("numberbinswaveform")->GetMaximum();
  waveform=new Double_t[nmaxbinswaveform+1];
  fTree->SetBranchAddress("waveform",waveform);

  Int_t nmaxtracks=(Int_t)fTree->GetLeaf("ntracks")->GetMaximum();
  trackx=new Double_t[nmaxtracks];
  tracky=new Double_t[nmaxtracks];
  trackz=new Double_t[nmaxtracks];
  trackzenith=new Double_t[nmaxtracks];
  trackazimuth=new Double_t[nmaxtracks];
  trackenergy=new Double_t[nmaxtracks];
  fTree->SetBranchAddress("trackx",trackx);
  fTree->SetBranchAddress("tracky",tracky);
  fTree->SetBranchAddress("trackz",trackz);
  fTree->SetBranchAddress("trackzenith",trackzenith);
  fTree->SetBranchAddress("trackazimuth",trackazimuth);
  fTree->SetBranchAddress("trackenergy",trackenergy);
  fTree->SetBranchAddress("tracktype",tracktype);

  // Prepare for loop over entries
  lastevent=0;

  // Loop over waveforms in input tree
  for(Int_t ientry=0; ientry<fTree->GetEntries(); ientry++)
  {
   fTree->GetEntry(ientry);
   // If new event
   if(eventid!=lastevent)
   {
    // Write old event to tree (if it contains data)
    if(evt->GetDevices("IceGOM"))
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

    // Daq: JEB for 2007 and later, TWR for 2006 and earlier
    daq.Reset(1);
    if(evt->GetDate()/10000 >= 2007){
     daq.AddNamedSlot("JEB");
     daq.SetSignal(1,"JEB");
    }
    else {
     daq.AddNamedSlot("TWR");
     daq.SetSignal(1,"TWR");
    }
    evt->AddDevice(daq);

    // Loop over all the tracks and add them to the current event
    Int_t nrecotracks=0;
    Int_t nmctracks=0;
    for (Int_t itrack=0; itrack<ntracks; itrack++)
    {
     t.Reset();

     TString tracktypestring(tracktype[itrack]);
     // Monte Carlo track
     if(tracktypestring.Index("MC")>=0){
      nmctracks++;
      t.SetId(-nmctracks);
     }
     // Reco track
     if(tracktypestring.Index("MC")<0){
      nrecotracks++;
      t.SetId(nrecotracks);
     }
     // Track ID and name
     t.SetName(tracktype[itrack]);
     t.SetTitle(tracktype[itrack]);

     // Beginpoint of the track
     vec[0]=trackx[itrack];
     vec[1]=tracky[itrack];
     vec[2]=trackz[itrack];
     r.SetPosition(vec,"car");
     t.SetBeginPoint(r);

     // Momentum in GeV/c
     vec[0]=trackenergy[itrack];
     if(vec[0]==0 || vec[0]!=vec[0]) vec[0]=1; // Energy unknown: set default value of 1 GeV
     vec[1]=pi-trackzenith[itrack];
     vec[2]=trackazimuth[itrack]+pi;
     if(vec[2]>=2*pi) vec[2]-=2*pi;
     p.SetVector(vec,"sph");
     t.Set3Momentum(p);

     // Check for unreconstructed tracks (NaN values) and reset track if necessary
     if(trackx[itrack]!=trackx[itrack] || tracky[itrack]!=tracky[itrack] || trackz[itrack]!=trackz[itrack] ||
        trackzenith[itrack]!=trackzenith[itrack] || trackazimuth[itrack]!=trackazimuth[itrack]){
      t.Reset();
     } else {
      // Add track to event
      evt->AddTrack(t);
     }
    }

    // Store trigger information
    if (twreventtimes) delete[] twreventtimes;
    twreventtimes=new Double_t[trignr];
    if (twreventorder) delete[] twreventorder;
    twreventorder=new Int_t[trignr];
    nrtwrevents=0;
    trig.Reset(1);
    // Loop over triggers
    for(Int_t itr=0; itr<trignr; itr++){
     // For TWR events:
     if(TString(trigsourceID[itr])=="AMANDA_TWR_DAQ"){
      // Add global time offset if there are no MC tracks (i.e. if the data are real data)
      if(!nmctracks) trigtime[itr]-=globalt0+twri3offset;
      // Check if this trigger belongs to a new event
      Int_t newevent=1;
      for(Int_t itwrevent=0; itwrevent<nrtwrevents; itwrevent++){
       if(fabs(trigtime[itr]-twreventtimes[itwrevent])<1500){
        newevent=0;
        break;
       }
      }
      // Add new TWR event time
      if(newevent){
       twreventtimes[nrtwrevents]=trigtime[itr];
       nrtwrevents++;
      }
     }

     // Add this trigger to trigger device
     s.Reset(1);
     sprintf(trigname,"%s/%s",trigsourceID[itr],trigtypeID[itr]);
     s.SetName(trigname);
     s.SetSlotName("trig_pulse_le",1);
     s.SetSignal(trigtime[itr],1);
     s.SetSlotName("trig_pulse_tot",2);
     s.SetSignal(triglength[itr],2);
     trig.AddHit(s);
    } // End of loop over triggers

    TMath::Sort(nrtwrevents,twreventtimes,twreventorder,0);

    // Add main trigger unless it is already present
    // TODO: select most appropriate trigger as artificial "main"
    // For now: select first trigger in event
    if(!trig.GetHit("main")){
     s.Reset(1);
     s.SetName("main");
     s.SetSlotName("trig_pulse_le",1);
     s.SetSignal(trigtime[0],1);
     s.SetSlotName("trig_pulse_tot",2);
     s.SetSignal(triglength[0],2);
     trig.AddHit(s);
    }
    evt->AddDevice(trig);

    // Remember event nr
    lastevent=eventid;
   }

   if(waveformtype==1){
    // Amanda module
    omid=aom.GetOMId(String,OM);
    if(omid==681) continue;  // Skip OM 681, which should never give data: avoid all risk of confusion
    om=(IceGOM*)evt->GetIdDevice(omid);
    if (!om)
    {
     aom.Reset(1);
     aom.SetUniqueID(omid);
     evt->AddDevice(aom);
     om=(IceGOM*)evt->GetIdDevice(omid);
    }
    if (!om) continue;
   }

   else if(waveformtype==2) {
    // Non-merged waveforms not supported
    cout << "IceRoot::Exec: Non-merged ATWD waveforms not supported in this module." << endl;
    continue;
   }

   else if(waveformtype==3 || waveformtype==4) {
    // InIce module
    omid=idom.GetOMId(String,OM);
    om=(IceGOM*)evt->GetIdDevice(omid);
    if (!om)
    {
     idom.Reset(1);
     idom.SetUniqueID(omid);
     evt->AddDevice(idom);
     om=(IceGOM*)evt->GetIdDevice(omid);
    }
    if (!om) continue;
   }

   else if(waveformtype==5 || waveformtype==6) {
    // IceTop module
    omid=tdom.GetOMId(String,OM);
    om=(IceGOM*)evt->GetIdDevice(omid);
    if (!om)
    {
     tdom.Reset(1);
     tdom.SetUniqueID(omid);
     evt->AddDevice(tdom);
     om=(IceGOM*)evt->GetIdDevice(omid);
    }
    if (!om) continue;
   }

   else {
    cout << "IceRoot::Exec: Unknown waveform type " << waveformtype << endl;
    continue;
   }

   // Fill the waveform histogram with this fragment
   // TWR waveform
   if(waveformtype==1){ 
    if(numberbinstwr>0 && twreventnr<nrtwrevents){
     histo.Reset();
     // Store baseline info
     hname="BASELINE-WF";
     hname+=om->GetNwaveforms()+1;
     om->AddNamedSlot(hname);
     om->SetSignal(baseline,hname);
     // Waveform name
     hname="OM";
     hname+=omid;
     hname+="-WF";
     hname+=om->GetNwaveforms()+1;
     hname+="-TWR";
     histo.SetName(hname.Data());
     // Add waveform
     calom=0;
     if(fJEBTDaq) calom=(IceGOM*)fJEBTDaq->GetObject(omid,1);
     if(!calom){
      cout << "IceRoot: No calibration info for OM " << omid << ", skipping this waveform" << endl;
      continue;
     }
     binwidth=calom->GetSignal("BINSIZE");
     if(binwidth<=0){
      cout << "IceRoot: Zero bin width for OM " << omid << ", skipping this waveform" << endl;
      continue;
     }
     twrbuffer=(Int_t)(twrtime/binwidth);
     extstop=(Int_t)calom->GetSignal("EXTSTOP");
     starttime=twreventtimes[twreventorder[twreventnr]]+binwidth*(startbin-twrbuffer+extstop);
     histo.SetBins(numberbinstwr,starttime,starttime+numberbinstwr*binwidth);
     for (Int_t jbin=1; jbin<=numberbinstwr; jbin++)
     {
      histo.SetBinContent(jbin,(Float_t)baseline-(Float_t)twrwaveform[jbin-1]);
     }
     om->SetWaveform(&histo,om->GetNwaveforms()+1);
    } else {
     cout << "IceRoot::Exec: waveformtype=1, but numberbinstwr=0 or nrtwrevents=0." << endl;
    }
   }

   // Separate waveforms from DOM launch
   else if(waveformtype==2){
    // Non-merged waveforms not supported
    cout << "IceRoot::Exec: Non-merged ATWD waveform found. You should never get to here anyway." << endl;
    continue;
   }

   // Merged ATWD waveform
   else if(waveformtype==3 || waveformtype==5){
    if(numberbinswaveform>0){
     histo.Reset();
     hname="OM";
     hname+=omid;
     hname+="-WF";
     hname+=om->GetNwaveforms()+1;
     hname+="-ATWD";
     histo.SetName(hname.Data());
     starttime=starttimewaveform;
     histo.SetBins(numberbinswaveform,starttime,starttime+numberbinswaveform*binwidth);
     sample.Reset();
     for (Int_t jbin=1; jbin<=numberbinswaveform; jbin++)
     {
      histo.SetBinContent(jbin,waveform[jbin-1]);
     }
     hname="BASELINE-WF";
     hname+=om->GetNwaveforms()+1;
     om->AddNamedSlot(hname);
     om->SetSignal(sample.GetMedian(&histo,2),hname);
     om->SetSignalError(sample.GetSpread(&histo,2),hname);
     om->SetWaveform(&histo,om->GetNwaveforms()+1);
    } else {
     cout << "IceRoot::Exec: waveformtype=" << waveformtype << ", but numberbinswaveform=0." << endl;
    }
   }

   // FADC waveform
   else if(waveformtype==4 || waveformtype==6){
    if(numberbinswaveform>0){
     histo.Reset();
     hname="OM";
     hname+=omid;
     hname+="-WF";
     hname+=om->GetNwaveforms()+1;
     hname+="-FADC";
     histo.SetName(hname.Data());
     starttime=starttimewaveform;
     histo.SetBins(numberbinswaveform,starttime,starttime+numberbinswaveform*binwidth);
     for (Int_t jbin=1; jbin<=numberbinswaveform; jbin++)
     {
      histo.SetBinContent(jbin,waveform[jbin-1]);
     }
     om->SetWaveform(&histo,om->GetNwaveforms()+1);
    } else {
     cout << "IceRoot::Exec: waveformtype=" << waveformtype << ", but numberbinswaveform=0." << endl;
    }
   }

   // Unknown waveform type
   else {
    cout << "IceRoot::Exec: Unknown waveform type " << waveformtype << endl;
    continue;
   }

  } // End of loop over waveforms in input tree

  // Write last event to tree
  if(evt->GetDevices("IceGOM"))
  {
   CleanTasks();
   ExecuteTasks(opt);
   if (otree) otree->Fill();
   if (fPrintfreq) { if (!(nevt%fPrintfreq)) evt->HeaderData(); }
   // Update event counter
   nevt++;
  }

  // Reset event
  evt->Reset();

  // Close input file
  input->Close();

  // Clean up
  if (trigtime) delete[] trigtime;
  if (triglength) delete[] triglength;
  if (twrwaveform) delete[] twrwaveform;
  if (waveform) delete[] waveform;
  if (trackx) delete[] trackx;
  if (tracky) delete[] tracky;
  if (trackz) delete[] trackz;
  if (trackzenith) delete[] trackzenith;
  if (trackazimuth) delete[] trackazimuth;
  if (trackenergy) delete[] trackenergy;
  if (twreventtimes)
  {
   delete[] twreventtimes;
   twreventtimes=0;
  }
  if (twreventorder)
  {
   delete[] twreventorder;
   twreventorder=0;
  }

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
