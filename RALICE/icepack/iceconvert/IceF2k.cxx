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
//
// Usage example :
// ---------------
//
// gSystem->Load("ralice");
// gSystem->Load("icepack");
// gSystem->Load("iceconvert");
//
// // Output file for the event structures
// TFile* ofile=new TFile("events.root","RECREATE","F2K data in IceEvent structure");
// TTree* otree=new TTree("T","Data of an Amanda run");
//
// // Limit the number of entries for testing
// Int_t nentries=300;
//
// // Print frequency to produce a short summary print every printfreq events
// Int_t printfreq=10;
//
// // Split level for the output structures
// Int_t split=2;
//
// // Buffer size for the output structures
// Int_t bsize=32000;
//
// IceF2k q("run8000.f2k",split,bsize);
// q.Loop(otree,nentries,printfreq);
//
// // Select various objects to be added to the output file
//
// AliObjMatrix* omdb=q.GetOMdbase();
// if (omdb) omdb->Write();
//
// AliDevice* fitdefs=q.GetFitdefs();
// if (fitdefs) fitdefs->Write();
//
// TDatabasePDG* pdg=q.GetPDG();
// if (pdg) pdg->Write();
//
// // Close output file
// ofile->Write();
// ofile->Close();
//
//--- Author: Nick van Eijndhoven 11-mar-2005 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceF2k.h"
#include "Riostream.h"

ClassImp(IceF2k) // Class implementation to enable ROOT I/O

IceF2k::IceF2k(char* fname,Int_t split,Int_t bsize)
{
// Default constructor.
// Initialise the input file and data structres to be converted.
// Also the required split level and buffer size of the output tree
// can be specified in this constructor.
// By default tree=0, split=0 and bsize=32000.

 fSplit=split;
 fBsize=bsize;

 fPdg=0;
 fOmdb=0;
 fFitdefs=0;

 if (!fname)
 {
  cout << " *IceF2k ctor* No data input file specified." << endl;
  return;
 }

 // Open the input file in the default ascii format (autodetection) for reading 
 fInput=rdmc_mcopen(fname,"r",RDMC_DEFAULT_ASCII_F);

 if (!fInput)
 {
  cout << " *IceF2k ctor* No input file found with name : " << fname << endl;
  return;
 }

 // Initialise the event structure 
 rdmc_init_mevt(&fEvent);

 // Read the file header information
 rdmc_rarr(fInput,&fHeader);
}
///////////////////////////////////////////////////////////////////////////
IceF2k::~IceF2k()
{
// Default destructor.
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
void IceF2k::Loop(TTree* otree,Int_t nentries,Int_t printfreq)
{
// Loop over the specified number of entries and convert the 
// F2K data into the IceEvent structure.
// The output will be written on the output tree specified as "otree".
// If otree=0, a default standard output tree will be created.
// If nentries<0 (default) all the entries of the input file
// will be processed.
// Every "printfreq" events a short event summary will be printed.
// The default value is printfreq=1.

 if (!fInput || fSplit<0) return;

 if (!otree) otree=new TTree("T","F2K Data");

 Double_t pi=acos(-1.);

 IceEvent* evt=new IceEvent();

 evt->SetTrackCopy(1);
 evt->SetDevCopy(1);

 // Branch in the tree for the event structure
 otree->Branch("IceEvent","IceEvent",&evt,fBsize,fSplit); 

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

 // Fill the database with geometry, calib. etc... parameters
 // for all the devices
 FillOMdbase();

 // Set the fit definitions according to the F2000 header info
 SetFitdefs();

/*************

 // The LEDA specific output data
 AliCalorimeter* ledaup=new AliCalorimeter(44,144);
 AliCalorimeter* ledalw=new AliCalorimeter(40,144);

 ledaup->SetName("LedaUp");
 ledalw->SetName("LedaDown");

 evt->InitLeda(ledaup);
 evt->InitLeda(ledalw);

 TDatime datim;
 Float_t pos[3],costh;
 AliSignal s;
 s.SetName("CPV signal ADC");
************************/

 for (Int_t jentry=0; jentry<nentries; jentry++)
 {
  if (rdmc_revt(fInput,&fHeader,&fEvent) != 0) break;

  // Reset the complete Event structure
  evt->Reset();

  evt->SetRunNumber(fEvent.nrun);
  evt->SetEventNumber(fEvent.enr);
  evt->SetMJD(fEvent.mjd,fEvent.secs,fEvent.nsecs);

  PutMcTracks(evt);

  PutRecoTracks(evt);

  PutHits(evt);

/*********************************
  datim.Set(Jdate,Jtime);
  evt->SetDayTime(datim);
  evt->SetProjectile(207,82,158);
  evt->SetTarget(207,82,0);
  evt->SetWeight(Jwscal);
  evt->SetTrig(Itword);
  evt->SetZdc(Zdc*1000.);
  evt->SetMiracE(1000.*Emir,Emire,Emirh);
  evt->SetMiracEt(Etm,Etme,Etmh);
 
  ledaup->Reset();
  ledalw->Reset();
  // Fill calorimeter with module data
  for (Int_t i=0; i<Nmod; i++)
  {
   if (Adcl[i] > 3) // Adc cut of 3 to remove noise
   {
    if (Irowl[i] > 0) ledaup->SetSignal(Irowl[i],Icoll[i],Adcl[i]);
    if (Irowl[i] < 0) ledalw->SetSignal(-Irowl[i],Icoll[i],Adcl[i]);
   }
  }

  // Store associated CPV signals
  for (Int_t j=0; j<Ncluv; j++)
  {
   s.Reset();
   s.SetSignal(Iadccv[j]);
   pos[1]=Thetacv[j]*pi/180.;
   pos[2]=Phicv[j]*pi/180.;
   costh=cos(pos[1]);
   pos[0]=0;
   if (costh) pos[0]=2103./costh;
   s.SetPosition(pos,"sph");
   pos[0]=0.4;
   pos[1]=2.2;
   pos[2]=0;
   s.SetPositionErrors(pos,"car");
   if (Phicv[j]>=0. && Phicv[j]<=180.)
   {
    ledaup->AddVetoSignal(s);
   }
   else
   {
    ledalw->AddVetoSignal(s);
   }
  }

  evt->AddDevice(ledaup);
  evt->AddDevice(ledalw);
************************************/

  if (!(jentry%printfreq))
  {
   evt->HeaderData();
  }

  // Write the complete structure to the output Tree
  otree->Fill();
 }

 if (evt) delete evt;
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::FillOMdbase()
{
// Fill the database with geometry, calib. etc... parameters 
// for all the devices.

 if (fHeader.nch<=0) return;

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
  dev->SetSlotName("TYPE",1);
  dev->SetSlotName("ORIENT",2);
  dev->SetSlotName("T0",3);
  dev->SetSlotName("ALPHA",4);
  dev->SetSlotName("KADC",5);
  dev->SetSlotName("KTOT",6);
  dev->SetSlotName("KTDC",7);

  pos[0]=fHeader.x[i];
  pos[1]=fHeader.y[i];
  pos[2]=fHeader.z[i];
  dev->SetPosition(pos,"car");
  dev->SetSignal(fHeader.type[i],1);
  dev->SetSignal((Float_t)fHeader.costh[i],2);
  dev->SetSignal(fHeader.cal[i].t_0,3);
  dev->SetSignal(fHeader.cal[i].alpha_t,4);
  dev->SetSignal(fHeader.cal[i].beta_a,5);
  dev->SetSignal(fHeader.cal[i].beta_tot,6);
  dev->SetSignal(fHeader.cal[i].beta_t,7);
  fOmdb->EnterObject(i+1,1,dev);
 }
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::SetFitdefs()
{
// Obtain the names of the variables for each fit procedure from the
// f2000 header. Each different fit procedure is then stored as a separate
// hit of an AliDevice object and the various fit variables are stored
// as separate signal slots of the corresponding hit.
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
//    Slot : 1 Signal value : 1 name : id
//    Slot : 2 Signal value : 2 name : rchi2
//    Slot : 3 Signal value : 3 name : prob
//    Slot : 4 Signal value : 4 name : sigth
//    Slot : 5 Signal value : 5 name : covmin
//    Slot : 6 Signal value : 6 name : covmax
//    Slot : 7 Signal value : 7 name : cutflag
//    Slot : 8 Signal value : 8 name : chi2
//  *AliSignal::Data* Id :1
//    Position Vector in car coordinates : 0 0 0
//    Err. in car coordinates : 0 0 0
//    Owned by device : AliDevice Name : FitDefinitions
//    Slot : 1 Signal value : 1 name : id
//    Slot : 2 Signal value : 2 name : rchi2
//    Slot : 3 Signal value : 3 name : prob
// etc....  
//
// This memberfunction is based on the original idea/code by Adam Bouchta.

 if (fHeader.n_fit<=0) return;

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

  for (Int_t j=0; j<fHeader.def_fit[i].nwords; j++)
  {
   s.SetSlotName(TString(fHeader.def_fit[i].words[j]),j+1);
   s.SetSignal(j+1,j+1);
  }

  fFitdefs->AddHit(s);
  s.Reset(1);
 }
}
///////////////////////////////////////////////////////////////////////////
void IceF2k::PutMcTracks(IceEvent* evt)
{
// Get the MC tracks from the F2000 file into the IcePack structure.
// Note : MC tracks are given negative track id's in the event structure.
// This memberfunction is based on the original code by Adam Bouchta.

 if (!evt || fEvent.ntrack<=0) return;

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
void IceF2k::PutRecoTracks(IceEvent* evt)
{
// Get the reconstructed tracks from the F2000 file into the IcePack structure.
// Note : Reco tracks are given positive track id's in the event structure.
// This memberfunction is based on the original code by Adam Bouchta.

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
  t.SetName(fPdg->GetParticle(idpdg)->GetName());
  t.SetTitle("RECO track");
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
void IceF2k::PutHits(IceEvent* evt)
{
// Get the hit and waveform info from the F2000 file into the IcePack structure.
// This memberfunction is based on the original code by Adam Bouchta.

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

  s.Reset();
  s.SetUniqueID(fEvent.h[i].id);
  s.SetSignal(fEvent.h[i].amp,1);
  s.SetSignal(fEvent.h[i].t,2);
  s.SetSignal(fEvent.h[i].tot,3);

  omx->AddHit(s);

  sx=omx->GetHit(omx->GetNhits());
  if (!sx) continue;

  // Bi-directional link between this hit and the track that caused the ADC value.
  // This F2K info is probably only present for MC tracks.
  tid=fEvent.h[i].ma;
  if (tid > 0)
  {
   tx=evt->GetIdTrack(tid); // Reco tracks
   if (!tx) tx=evt->GetIdTrack(-tid); // MC tracks
   if (tx) sx->AddLink(tx);
  }
  else
  {
   if (tid == -2) s.SetNameTitle("N","Noise");
   if (tid == -3) s.SetNameTitle("A","Afterpulse");
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

  omx->SetSlotName("BASELINE",omx->GetNnames()+1);
  omx->SetSignal(fEvent.wf[iwf].baseline,"BASELINE");

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
   histo.SetBinContent(jbin,fEvent.wf[iwf].digi[jbin]);
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
    if (sx) sx->AddLink(hypx);
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
