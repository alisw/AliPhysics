/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id: AliCollider.cxx,v 1.12 2004/05/04 15:33:04 nick Exp $

///////////////////////////////////////////////////////////////////////////
// Class AliCollider
// Pythia based universal physics event generator.
// This class is derived from TPythia6 and has some extensions to
// support also generation of nucleus-nucleus interactions and to allow
// investigation of the effect of detector resolving power.
// Furthermore, the produced event information is provided in a format
// using the AliEvent structure.
// For the produced AliTrack objects, the particle ID code is set to the
// Pythia KF value, which is compatible with the PDG identifier.
// This will allow a direct analysis of the produced data using the
// Ralice physics analysis tools.
//
// For further details concerning the produced output structure,
// see the docs of the memberfunctions SetVertexMode and SetResolution.
//
// Example job of minimum biased Pb+Pb interactions :
// --------------------------------------------------
// {
//  gSystem->Load("libEG");
//  gSystem->Load("libEGPythia6");
//  gSystem->Load("ralice");
//
//  AliCollider* gen=new AliCollider();
//
//  gen->SetOutputFile("test.root");
//  gen->SetVertexMode(3);    
//  gen->SetResolution(1e-6); // 1 micron vertex resolution
//
//  gen->SetRunNumber(1);
//
//  Int_t zp=82;
//  Int_t ap=208;
//  Int_t zt=82;
//  Int_t at=208;
//
//  gen->Init("fixt",zp,ap,zt,at,158);
//
//  gen->SetTitle("SPS Pb-Pb collision at 158A GeV/c beam energy");
//
//  Int_t nevents=5;
//
//  AliRandom rndm;
//  Float_t* rans=new Float_t[nevents];
//  rndm.Uniform(rans,nevents,2,ap+at);
//  Int_t npart;
//  for (Int_t i=0; i<nevents; i++)
//  {
//   npart=rans[i];
//   gen->MakeEvent(npart);
//
//   AliEvent* evt=gen->GetEvent();
//  
//   evt->List();
//  }
//
//  gen->EndRun();
// }
//
//
// Example job of a cosmic nu+p atmospheric interaction.
// -----------------------------------------------------
// {
//  gSystem->Load("libEG");
//  gSystem->Load("libEGPythia6");
//  gSystem->Load("ralice");
//
//  AliCollider* gen=new AliCollider();
//
//  gen->SetOutputFile("test.root");
//
//  gen->SetRunNumber(1);
//
//  gen->Init("fixt","nu_mu","p",1e11);
//
//  gen->SetTitle("Atmospheric nu_mu-p interaction at 1e20 eV");
//
//  Int_t nevents=10;
//
//  for (Int_t i=0; i<nevents; i++)
//  {
//   gen->MakeEvent(0,1);
//
//   AliEvent* evt=gen->GetEvent();
//  
//   evt->Data();
//  }
//
//  gen->EndRun();
// }
//
//
//--- Author: Nick van Eijndhoven 22-nov-2002 Utrecht University
//- Modified: NvE $Date: 2004/05/04 15:33:04 $ Utrecht University
///////////////////////////////////////////////////////////////////////////

#include "AliCollider.h"
#include "Riostream.h"
 
ClassImp(AliCollider) // Class implementation to enable ROOT I/O
 
AliCollider::AliCollider() : TPythia6()
{
// Default constructor.
// All variables initialised to default values.
//
// Some Pythia default MC parameters are automatically modified to provide
// more suitable running conditions for soft processes in view of
// nucleus-nucleus interactions and astrophysical processes.
// The user may initialise the generator with all the default Pythia
// parameters and obtain full user control to modify the settings by means
// of the SetUserControl memberfunction.
//
// Refer to the SetElastic memberfunction for the inclusion of elastic
// and diffractive processes.
// By default these processes are not included.

 fVertexmode=0;    // No vertex structure creation
 fResolution=1e-7; // Standard resolution is 0.1 micron
 fRunnum=0;
 fEventnum=0;
 fPrintfreq=1;
 fUserctrl=0; // Automatic optimisation of some MC parameters 
 fElastic=0;  // No elastic and diffractive processes

 fEvent=0;

 fSpecpmin=0;

 fFrame="none";
 fWin=0;

 fNucl=0;
 fZproj=0;
 fAproj=0;
 fZtarg=0;
 fAtarg=0;
 fFracpp=0;
 fFracnp=0;
 fFracpn=0;
 fFracnn=0;

 fOutFile=0;
 fOutTree=0;

 fSelections=0;
 fSelect=0;

 TString s=GetName();
 s+=" (AliCollider)";
 SetName(s.Data());
}
///////////////////////////////////////////////////////////////////////////
AliCollider::~AliCollider()
{
// Default destructor
 if (fEvent)
 {
  delete fEvent;
  fEvent=0;
 }
 if (fOutFile)
 {
  delete fOutFile;
  fOutFile=0;
 }
 if (fOutTree)
 {
  delete fOutTree;
  fOutTree=0;
 }
 if (fSelections)
 {
  delete fSelections;
  fSelections=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SetOutputFile(TString s)
{
// Create the output file containing all the data in ROOT output format.
 if (fOutFile)
 {
  delete fOutFile;
  fOutFile=0;
 }
 fOutFile=new TFile(s.Data(),"RECREATE","AliCollider data");

 if (fOutTree)
 {
  delete fOutTree;
  fOutTree=0;
 }
 fOutTree=new TTree("T","AliCollider event data");

 Int_t bsize=32000;
 Int_t split=0;
 fOutTree->Branch("Events","AliEvent",&fEvent,bsize,split);
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SetVertexMode(Int_t mode)
{
// Set the mode of the vertex structure creation.
//
// By default all generated tracks will only appear in the AliEvent
// structure without any primary (and secondary) vertex structure.
// The user can build the vertex structure if he/she wants by means
// of the beginpoint location of each AliTrack.
//
// However, one can also let AliCollider automatically create
// the primary (and secondary) vertex structure(s).
// In this case the primary vertex is given Id=1 and all sec. vertices
// are given Id's 2,3,4,....
// All vertices are created as standalone entities in the AliEvent structure
// without any linking between the various vertices.
// For this automated process, the user-selected resolution
// (see SetResolution) is used to decide whether or not certain vertex
// locations can be resolved.
// In case no vertex creation is selected (i.e. the default mode=0),
// the value of the resolution is totally irrelevant.
//
// The user can also let AliCollider automatically connect the sec. vertices
// to the primary vertex (i.e. mode=3). This process will also automatically
// generate the tracks connecting the vertices.
// Note that the result of the mode=3 operation may be very sensitive to
// the resolution parameter. Therefore, no attempt is made to distinguish
// between secondary, tertiary etc... vertices. All sec. vertices are
// linked to the primary one.
//  
// Irrespective of the selected mode, all generated tracks can be obtained
// directly from the AliEvent structure.
// In case (sec.) vertex creation is selected, all generated vertices can
// also be obtained directly from the AliEvent structure. 
// These (sec.) vertices contain only the corresponding pointers to the various
// tracks which are stored in the AliEvent structure.
//
// Overview of vertex creation modes :
// -----------------------------------
// mode = 0 ==> No vertex structure will be created
//        1 ==> Only primary vertex structure will be created
//        2 ==> Unconnected primary and secondary vertices will be created
//        3 ==> Primary and secondary vertices will be created where all the
//              sec. vertices will be connected to the primary vertex.
//              Also the vertex connecting tracks will be automatically
//              generated. 
//
 if (mode<0 || mode >3)
 {
  cout << " *AliCollider::SetVertexMode* Invalid argument mode : " << mode << endl;
  fVertexmode=0;
 }
 else
 {
  fVertexmode=mode;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCollider::GetVertexMode() const
{
// Provide the current mode for vertex structure creation.
 return fVertexmode;
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SetResolution(Double_t res)
{
// Set the resolution (in meter) for resolving (sec.) vertices.
// By default this resolution is set to 0.1 micron.
// Note : In case no vertex creation has been selected, the value of
//        the resolution is totally irrelevant.
 fResolution=fabs(res);
}
///////////////////////////////////////////////////////////////////////////
Double_t AliCollider::GetResolution() const
{
// Provide the current resolution (in meter) for resolving (sec.) vertices.
 return fResolution;
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SetRunNumber(Int_t run)
{
// Set the user defined run number.
// By default the run number is set to 0.
 fRunnum=run;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCollider::GetRunNumber() const
{
// Provide the user defined run number.
 return fRunnum;
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SetPrintFreq(Int_t n)
{
// Set the print frequency for every 'n' events.
// By default the printfrequency is set to 1 (i.e. every event).
 fPrintfreq=n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCollider::GetPrintFreq() const
{
// Provide the user selected print frequency.
 return fPrintfreq;
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SetUserControl(Int_t flag)
{
// Set the user control flag w.r.t. disabling automatic optimisation
// of some Pythia default MC parameters for soft interactions in view of
// nucleus-nucleus collisions and astrophysical processes.
// Flag = 0 : Limited user control (automatic optimisation enabled)
//        1 : Full user control (automatic optimisation disabled)
// By default the user control is set to 0 (i.e. automatic optimisation).
// See the Init() memberfunctions for further details w.r.t. the optimisations.
 fUserctrl=flag;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCollider::GetUserControl() const
{
// Provide the value of the user control flag.
 return fUserctrl;
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SetElastic(Int_t flag)
{
// Set the flag w.r.t. inclusion of elastic and diffractive processes.
// By default these processes are not included.
// Flag = 0 : Do not include elastic and diffractive processes
//        1 : Elastic and diffractive processes will be included
 fElastic=flag;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCollider::GetElastic() const
{
// Provide the value of the control flag for elastic and diffractive processes.
 return fElastic;
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::Init(char* frame,char* beam,char* target,Float_t win)
{
// Initialisation of the underlying Pythia generator package.
// The event number is reset to 0.
// This routine just invokes TPythia6::Initialize(...) and the arguments
// have the corresponding meaning.
// Some Pythia default MC parameters are automatically modified to provide
// more suitable running conditions for soft processes in view of
// astrophysical processes.
// The optimisations consist of : 
// * Usage of real photons for photon beams or targets
// * Minimum CMS energy of 3 GeV for the event
// * Activation of the default K factor values
//   with separate settings for ordinary and color annihilation graphs.
// The user may initialise the generator with all the default Pythia
// parameters and obtain full user control to modify the settings by means
// of invoking the SetUserControl memberfunction before this initialisation.
// Note that the inclusion of elastic and diffractive processes is controlled
// by invokation of the SetElastic memberfunction before this initialisation,
// irrespective of the UserControl selection.

 if (!fUserctrl) // Optimisation of some MC parameters
 {
  SetMSTP(14,10); // Real photons for photon beams or targets
  SetPARP(2,3.);  // Minimum CMS energy for the event
  SetMSTP(33,2);  // Activate K factor. Separate for ordinary and color annih. graphs
 }

 if (fElastic) SetMSEL(2); // Include low-Pt, elastic and diffractive events

 fEventnum=0;
 fNucl=0;
 fFrame=frame;
 fWin=win;
 Initialize(frame,beam,target,win);

 cout << endl;
 cout << " *AliCollider::Init* Standard Pythia initialisation." << endl;
 cout << " Beam particle : " << beam << " Target particle : " << target
      << " Frame = " << frame << " Energy = " << win
      << endl;
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::Init(char* frame,Int_t zp,Int_t ap,Int_t zt,Int_t at,Float_t win)
{
// Initialisation of the underlying Pythia generator package for the generation
// of nucleus-nucleus interactions.
// The event number is reset to 0.
// In addition to the Pythia standard arguments 'frame' and 'win', the user
// can specify here (Z,A) values of the projectile and target nuclei.
//
// Note : The 'win' value denotes either the cms energy per nucleon-nucleon collision
//        (i.e. frame="cms") or the momentum per nucleon in all other cases.
//
// Some Pythia default MC parameters are automatically modified to provide
// more suitable running conditions for soft processes in view of
// nucleus-nucleus interactions and astrophysical processes.
// The optimisations consist of : 
// * Minimum CMS energy of 3 GeV for the event
// * Activation of the default K factor values
//   with separate settings for ordinary and color annihilation graphs.
// The user may initialise the generator with all the default Pythia
// parameters and obtain full user control to modify the settings by means
// of invoking the SetUserControl memberfunction before this initialisation.
// Note that the inclusion of elastic and diffractive processes is controlled
// by invokation of the SetElastic memberfunction before this initialisation,
// irrespective of the UserControl selection.

 if (!fUserctrl) // Optimisation of some MC parameters
 {
  SetPARP(2,3.);  // Minimum CMS energy for the event
  SetMSTP(33,2);  // Activate K factor. Separate for ordinary and color annih. graphs
 }

 if (fElastic) SetMSEL(2); // Include low-Pt, elastic and diffractive events

 fEventnum=0;
 fNucl=1;
 fFrame=frame;
 fWin=win;
 fZproj=0;
 fAproj=0;
 fZtarg=0;
 fAtarg=0;
 fFracpp=0;
 fFracnp=0;
 fFracpn=0;
 fFracnn=0;

 if (ap<1 || at<1 || zp>ap || zt>at)
 {
  cout << endl;
  cout << " *AliCollider::Init* Invalid input value(s). Zproj = " << zp
       << " Aproj = " << ap << " Ztarg = " << zt << " Atarg = " << at << endl;
  return;
 }

 fZproj=zp;
 fAproj=ap;
 fZtarg=zt;
 fAtarg=at;

 cout << endl;
 cout << " *AliCollider::Init* Nucleus-Nucleus generator initialisation." << endl;
 cout << " Zproj = " << zp << " Aproj = " << ap << " Ztarg = " << zt << " Atarg = " << at
      << " Frame = " << frame << " Energy = " << win
      << endl;
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::GetFractions(Float_t zp,Float_t ap,Float_t zt,Float_t at)
{
// Determine the fractions for the various N-N collision processes.
// The various processes are : p+p, n+p, p+n and n+n.
 if (zp<0) zp=0;
 if (zt<0) zt=0;

 fFracpp=0;
 fFracnp=0;
 fFracpn=0;
 fFracnn=0;

 if (ap>0 && at>0)
 {
  fFracpp=(zp/ap)*(zt/at);
  fFracnp=(1.-zp/ap)*(zt/at);
  fFracpn=(zp/ap)*(1.-zt/at);
  fFracnn=(1.-zp/ap)*(1.-zt/at);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::MakeEvent(Int_t npt,Int_t mlist,Int_t medit)
{
// Generate one event.
// In case of a nucleus-nucleus interaction, the argument 'npt' denotes
// the number of participant nucleons.
// Normally also the spectator tracks will be stored into the event structure.
// The spectator tracks have a negative user Id to distinguish them from the
// ordinary generated tracks.
// In case the user has selected the creation of vertex structures, the spectator
// tracks will be linked to the primary vertex.
// However, specification of npt<0 will suppress the storage of spectator tracks.
// In the latter case abs(npt) will be taken as the number of participants.  
// In case of a standard Pythia run for 'elementary' particle interactions,
// the value of npt is totally irrelevant.
//
// The argument 'mlist' denotes the list mode used for Pylist().
// Note : mlist<0 suppresses the invokation of Pylist().
// By default, no listing is produced (i.e. mlist=-1).
//
// The argument 'medit' denotes the edit mode used for Pyedit().
// Note : medit<0 suppresses the invokation of Pyedit().
// By default, only 'stable' final particles are kept (i.e. medit=1). 
//
// In the case of a standard Pythia run concerning 'elementary' particle
// interactions, the projectile and target particle ID's for the created
// event structure are set to the corresponding Pythia KF codes.
// All the A and Z values are in that case set to zero.
// In case of a nucleus-nucleus interaction, the proper A and Z values for 
// the projectile and target particles are set in the event structure.
// However, in this case both particle ID's are set to zero.
//
// Note : Only in case an event passed the selection criteria as specified
//        via SelectEvent(), the event will appear on the output file.

 fEventnum++; 

 Int_t specmode=1;
 if (npt<0)
 {
  specmode=0;
  npt=abs(npt);
 }

 // Counters for the various (proj,targ) combinations : p+p, n+p, p+n and n+n
 Int_t ncols[4]={0,0,0,0};

 Int_t zp=0;
 Int_t ap=0;
 Int_t zt=0;
 Int_t at=0;

 Int_t ncol=1;
 if (fNucl)
 {
  if (npt<1 || npt>(fAproj+fAtarg))
  {
   cout << " *AliCollider::MakeEvent* Invalid input value. npt = " << npt
        << " Aproj = " << fAproj << " Atarg = " << fAtarg << endl;
   return;
  }

  // Determine the number of nucleon-nucleon collisions
  ncol=npt/2;
  if (npt%2 && fRan.Uniform()>0.5) ncol+=1;

  // Determine the number of the various types of N+N interactions
  zp=fZproj;
  ap=fAproj;
  zt=fZtarg;
  at=fAtarg;
  Int_t maxa=2; // Indicator whether proj (1) or target (2) has maximal A left
  if (ap>at) maxa=1;
  Float_t* rans=new Float_t[ncol];
  fRan.Uniform(rans,ncol);
  Float_t rndm=0;
  for (Int_t i=0; i<ncol; i++)
  {
   GetFractions(zp,ap,zt,at);
   rndm=rans[i];
   if (rndm<=fFracpp) // p+p interaction
   {
    ncols[0]++;
    if (maxa==2)
    {
     at--;
     zt--;
    } 
    else
    {
     ap--;
     zp--;
    }
   }
   if (rndm>fFracpp && rndm<=(fFracpp+fFracnp)) // n+p interaction
   {
    ncols[1]++;
    if (maxa==2)
    {
     at--;
     zt--;
    } 
    else
    {
     ap--;
    }
   }
   if (rndm>(fFracpp+fFracnp) && rndm<=(fFracpp+fFracnp+fFracpn)) // p+n interaction
   {
    ncols[2]++;
    if (maxa==2)
    {
     at--;
    } 
    else
    {
     ap--;
     zp--;
    }
   }
   if (rndm>(fFracpp+fFracnp+fFracpn)) // n+n interaction
   {
    ncols[3]++; 
    if (maxa==2)
    {
     at--;
    } 
    else
    {
     ap--;
    }
   }
  }
  delete [] rans;
 }

 if (!(fEventnum%fPrintfreq))
 {
  cout << " *AliCollider::MakeEvent* Run : " << fRunnum << " Event : " << fEventnum
       << endl;
  if (fNucl)
  {
   cout << " npart = " << npt << " ncol = " << ncol 
        << " ncolpp = " << ncols[0] << " ncolnp = " << ncols[1]
        << " ncolpn = " << ncols[2] << " ncolnn = " << ncols[3] << endl;
  }
 }

 if (!fEvent)
 {
  fEvent=new AliEvent();
  fEvent->SetOwner();
  fEvent->SetName(GetName());
  fEvent->SetTitle(GetTitle());
 }

 fEvent->Reset();
 fEvent->SetRunNumber(fRunnum);
 fEvent->SetEventNumber(fEventnum);

 AliTrack t;
 Ali3Vector p;
 AliPosition r,rx;
 Float_t v[3];
 AliVertex vert;
 Ali3Vector pproj,ptarg;

 if (fVertexmode)
 {
  // Make sure the primary vertex gets correct location and Id=1
  v[0]=0;
  v[1]=0;
  v[2]=0;
  r.SetPosition(v,"car");
  v[0]=fResolution;
  v[1]=fResolution;
  v[2]=fResolution;
  r.SetPositionErrors(v,"car");

  vert.SetId(1);
  vert.SetTrackCopy(0);
  vert.SetVertexCopy(0);
  vert.SetPosition(r);
  fEvent->AddVertex(vert,0);
 }

 Int_t kf=0;
 Float_t charge=0,mass=0;
 TString name;

 Int_t ntypes=4;

 // Singular settings for a normal Pythia elementary particle interation 
 if (!fNucl)
 {
  ntypes=1;
  ncols[0]=1;
 }

 // Generate all the various collisions
 fSelect=0;      // Flag to indicate whether the total event is selected or not 
 Int_t select=0; // Flag to indicate whether the sub-event is selected or not 
 Int_t first=1;  // Flag to indicate the first collision process
 Double_t pnucl;
 Int_t npart=0,ntk=0;
 Double_t dist=0;
 for (Int_t itype=0; itype<ntypes; itype++)
 {
  if (fNucl)
  {
   if (itype==0 && ncols[itype]) Initialize(fFrame,"p","p",fWin);
   if (itype==1 && ncols[itype]) Initialize(fFrame,"n","p",fWin);
   if (itype==2 && ncols[itype]) Initialize(fFrame,"p","n",fWin);
   if (itype==3 && ncols[itype]) Initialize(fFrame,"n","n",fWin);
  }
  for (Int_t jcol=0; jcol<ncols[itype]; jcol++)
  {
   GenerateEvent();

   select=IsSelected();
   if (select) fSelect=1;

   if (first) // Store generator parameter information in the event structure
   {
    // Enter generator parameters as a device in the event
    AliSignal params;
    params.SetNameTitle("AliCollider","AliCollider generator parameters");
    params.SetSlotName("Medit",1);
    params.SetSlotName("Vertexmode",2);
    params.SetSlotName("Resolution",3);
    params.SetSlotName("Userctrl",4);
    params.SetSlotName("Elastic",5);

    params.SetSignal(medit,1);
    params.SetSignal(fVertexmode,2);
    params.SetSignal(fResolution,3);
    params.SetSignal(fUserctrl,4);
    params.SetSignal(fElastic,5);

    // Store projectile and target information in the event structure
    if (fNucl)
    {
     v[0]=GetP(1,1);
     v[1]=GetP(1,2);
     v[2]=GetP(1,3);
     pproj.SetVector(v,"car");
     pnucl=pproj.GetNorm();
     fEvent->SetProjectile(fAproj,fZproj,pnucl);
     v[0]=GetP(2,1);
     v[1]=GetP(2,2);
     v[2]=GetP(2,3);
     ptarg.SetVector(v,"car");
     pnucl=ptarg.GetNorm();
     fEvent->SetTarget(fAtarg,fZtarg,pnucl);

     params.AddNamedSlot("specmode");
     params.AddNamedSlot("Specpmin");
     params.AddNamedSlot("npart");
     params.AddNamedSlot("ncolpp");
     params.AddNamedSlot("ncolnp");
     params.AddNamedSlot("ncolpn");
     params.AddNamedSlot("ncolnn");

     params.SetSignal(specmode,"specmode");
     params.SetSignal(fSpecpmin,"Specpmin");
     params.SetSignal(npt,"npart");
     params.SetSignal(ncols[0],"ncolpp");
     params.SetSignal(ncols[1],"ncolnp");
     params.SetSignal(ncols[2],"ncolpn");
     params.SetSignal(ncols[3],"ncolnn");
    }
    else
    {
     v[0]=GetP(1,1);
     v[1]=GetP(1,2);
     v[2]=GetP(1,3);
     pnucl=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
     kf=GetK(1,2);
     fEvent->SetProjectile(0,0,pnucl,kf);
     v[0]=GetP(2,1);
     v[1]=GetP(2,2);
     v[2]=GetP(2,3);
     pnucl=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
     kf=GetK(2,2);
     fEvent->SetTarget(0,0,pnucl,kf);
    }

    fEvent->AddDevice(params);

    first=0;
   }

   if (medit >= 0) Pyedit(medit); // Define which particles are to be kept

   if (mlist>=0 && select)
   {
    Pylist(mlist);
    cout << endl;
   }

   npart=GetN();
   for (Int_t jpart=1; jpart<=npart; jpart++)
   {
    kf=GetK(jpart,2);
    charge=Pychge(kf)/3.;
    mass=GetP(jpart,5);
    name=GetPyname(kf);

    // 3-momentum in GeV/c
    v[0]=GetP(jpart,1);
    v[1]=GetP(jpart,2);
    v[2]=GetP(jpart,3);
    p.SetVector(v,"car");

    // Production location in meter.
    v[0]=GetV(jpart,1)/1000.;
    v[1]=GetV(jpart,2)/1000.;
    v[2]=GetV(jpart,3)/1000.;
    r.SetPosition(v,"car");

    ntk++;

    t.Reset();
    t.SetId(ntk);
    t.SetParticleCode(kf);
    t.SetName(name.Data());
    t.SetCharge(charge);
    t.SetMass(mass);
    t.Set3Momentum(p);
    t.SetBeginPoint(r);

    fEvent->AddTrack(t);

    // Build the vertex structures if requested
    if (fVertexmode)
    {
     // Check if track belongs within the resolution to an existing vertex
     Int_t add=0;  
     for (Int_t jv=1; jv<=fEvent->GetNvertices(); jv++)
     {
      AliVertex* vx=fEvent->GetVertex(jv);
      if (vx)
      {
       rx=vx->GetPosition();
       dist=rx.GetDistance(r);
       if (dist < fResolution)
       {
        AliTrack* tx=fEvent->GetIdTrack(ntk);
        if (tx)
        { 
         vx->AddTrack(tx);
         add=1;
        }
       }
      }
      if (add) break; // No need to look further for vertex candidates
     }

     // If track was not close enough to an existing vertex
     // a new secondary vertex is created      
     if (!add && fVertexmode>1)
     {
      AliTrack* tx=fEvent->GetIdTrack(ntk);
      if (tx)
      {
       v[0]=fResolution;
       v[1]=fResolution;
       v[2]=fResolution;
       r.SetPositionErrors(v,"car");
       vert.Reset();
       vert.SetTrackCopy(0);
       vert.SetVertexCopy(0);
       vert.SetId((fEvent->GetNvertices())+1);
       vert.SetPosition(r);
       vert.AddTrack(tx);
       fEvent->AddVertex(vert,0);
      } 
     }
    }
   } // End of loop over the produced particles for each collision
  } // End of loop over number of collisions for each type
 } // End of loop over collision types

 // Link sec. vertices to primary if requested
 // Note that also the connecting tracks are automatically created
 if (fVertexmode>2)
 {
  AliVertex* vp=fEvent->GetIdVertex(1); // Primary vertex
  if (vp)
  {
   for (Int_t i=2; i<=fEvent->GetNvertices(); i++)
   {
    AliVertex* vx=fEvent->GetVertex(i);
    if (vx)
    {
     if (vx->GetId() != 1) vp->AddVertex(vx);
    }
   }
  }
 }

 // Include the spectator tracks in the event structure.
 if (fNucl && specmode)
 {
  v[0]=0;
  v[1]=0;
  v[2]=0;
  r.SetPosition(v,"car");

  zp=fZproj-(ncols[0]+ncols[2]);
  if (zp<0) zp=0;
  ap=fAproj-(ncols[0]+ncols[1]+ncols[2]+ncols[3]);
  if (ap<0) ap=0;
  zt=fZtarg-(ncols[0]+ncols[1]);
  if (zt<0) zt=0;
  at=fAtarg-(ncols[0]+ncols[1]+ncols[2]+ncols[3]);
  if (at<0) at=0;

  Int_t nspec=0;

  if (pproj.GetNorm() > fSpecpmin)
  {
   kf=2212; // Projectile spectator protons
   charge=Pychge(kf)/3.;
   mass=GetPMAS(Pycomp(kf),1);
   name=GetPyname(kf);
   for (Int_t iprojp=1; iprojp<=zp; iprojp++)
   {
    nspec++;
    t.Reset();
    t.SetId(-nspec);
    t.SetParticleCode(kf);
    t.SetName(name.Data());
    t.SetTitle("Projectile spectator proton");
    t.SetCharge(charge);
    t.SetMass(mass);
    t.Set3Momentum(pproj);
    t.SetBeginPoint(r);

    fEvent->AddTrack(t);
   }

   kf=2112; // Projectile spectator neutrons
   charge=Pychge(kf)/3.;
   mass=GetPMAS(Pycomp(kf),1);
   name=GetPyname(kf);
   for (Int_t iprojn=1; iprojn<=(ap-zp); iprojn++)
   {
    nspec++;
    t.Reset();
    t.SetId(-nspec);
    t.SetParticleCode(kf);
    t.SetName(name.Data());
    t.SetTitle("Projectile spectator neutron");
    t.SetCharge(charge);
    t.SetMass(mass);
    t.Set3Momentum(pproj);
    t.SetBeginPoint(r);

    fEvent->AddTrack(t);
   }
  }

  if (ptarg.GetNorm() > fSpecpmin)
  {
   kf=2212; // Target spectator protons
   charge=Pychge(kf)/3.;
   mass=GetPMAS(Pycomp(kf),1);
   name=GetPyname(kf);
   for (Int_t itargp=1; itargp<=zt; itargp++)
   {
    nspec++;
    t.Reset();
    t.SetId(-nspec);
    t.SetParticleCode(kf);
    t.SetName(name.Data());
    t.SetTitle("Target spectator proton");
    t.SetCharge(charge);
    t.SetMass(mass);
    t.Set3Momentum(ptarg);
    t.SetBeginPoint(r);

    fEvent->AddTrack(t);
   }

   kf=2112; // Target spectator neutrons
   charge=Pychge(kf)/3.;
   mass=GetPMAS(Pycomp(kf),1);
   name=GetPyname(kf);
   for (Int_t itargn=1; itargn<=(at-zt); itargn++)
   {
    nspec++;
    t.Reset();
    t.SetId(-nspec);
    t.SetParticleCode(kf);
    t.SetName(name.Data());
    t.SetTitle("Target spectator neutron");
    t.SetCharge(charge);
    t.SetMass(mass);
    t.Set3Momentum(ptarg);
    t.SetBeginPoint(r);

    fEvent->AddTrack(t);
   }
  }

 // Link the spectator tracks to the primary vertex.
 if (fVertexmode)
 {
  AliVertex* vp=fEvent->GetIdVertex(1);
  if (vp)
  {
   for (Int_t ispec=1; ispec<=nspec; ispec++)
   {
    AliTrack* tx=fEvent->GetIdTrack(-ispec);
    if (tx) vp->AddTrack(tx);
   }
  }
 }
}

 if (!(fEventnum%fPrintfreq) && (mlist || fEvent))
 {
  if (fEvent)
  {
   cout << " Number of tracks in the event structure : "
        << fEvent->GetNtracks() << endl;
  }
  cout << endl; // Create empty output line after the event
 }

 if (fOutTree && fSelect) fOutTree->Fill();
}
///////////////////////////////////////////////////////////////////////////
AliEvent* AliCollider::GetEvent(Int_t select) const
{
// Provide pointer to the generated event structure.
//
// select = 0 : Always return the pointer to the generated event.
//          1 : Only return the pointer to the generated event in case
//              the event passed the selection criteria as specified via
//              SelectEvent(). Otherwise the value 0 will be returned.
//
// By invoking GetEvent() the default of select=0 will be used.

 if (!select || fSelect)
 {
  return fEvent;
 }
 else
 {
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::EndRun()
{
// Properly close the output file (if needed).
 if (fOutFile)
 {
  fOutFile->Write();
  fOutFile->Close();
  cout << " *AliCollider::EndRun* Output file correctly closed." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SetStable(Int_t id,Int_t mode)
{
// Declare whether a particle must be regarded as stable or not.
// The parameter "id" indicates the Pythia KF particle code, which
// basically is the PDG particle identifier code.
// The parameter "mode" indicates the action to be taken.
//
// mode = 0 : Particle will be able to decay
//        1 : Particle will be regarded as stable.
//
// In case the user does NOT explicitly invoke this function, the standard
// Pythia settings for the decay tables are used.
//
// When this function is invoked without the "mode" argument, then the
// default of mode=1 will be used for the specified particle.
//
// Notes :
// -------
// 1) This function should be invoked after the initialisation call
//    to AliCollider::Init.
// 2) Due to the internals of Pythia, there is no need to specify particles
//    and their corresponding anti-particles separately as (un)stable.
//    Once a particle has been declared (un)stable, the corresponding 
//    anti-particle will be treated in the same way.

 if (mode==0 || mode==1)
 {
  Int_t kc=Pycomp(id);
  Int_t decay=1-mode;
  if (kc>0)
  {
   SetMDCY(kc,1,decay);
  }
  else 
  {
   cout << " *AliCollider::SetStable* Unknown particle code. id = " << id << endl;
  }
 }
 else
 {
  cout << " *AliCollider::SetStable* Invalid parameter. mode = " << mode << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SelectEvent(Int_t id)
{
// Add a particle to the event selection list.
// The parameter "id" indicates the Pythia KF particle code, which
// basically is the PDG particle identifier code.
// In case the user has built a selection list via this procedure, only the
// events in which one of the particles specified in the list was generated
// will be kept. 
// The investigation of the generated particles takes place when the complete
// event is in memory, including all (shortlived) mother particles and resonances.
// So, the settings of the various particle decay modes have no influence on
// the event selection described here.
//
// If no list has been specified, all events will be accepted.  
//
// Note : id=0 will delete the selection list.
//
// Be aware of the fact that severe selection criteria (i.e. selecting only
// rare events) may result in long runtimes before an event sample has been
// obtained.
//
 if (!id)
 {
  if (fSelections)
  {
   delete fSelections;
   fSelections=0;
  }
 }
 else
 {
  Int_t kc=Pycomp(id);
  if (!fSelections)
  {
   fSelections=new TArrayI(1);
   fSelections->AddAt(kc,0);
  }
  else
  {
   Int_t exist=0;
   Int_t size=fSelections->GetSize();
   for (Int_t i=0; i<size; i++)
   {
    if (kc==fSelections->At(i))
    {
     exist=1;
     break;
    }
   }
  
   if (!exist)
   {
    fSelections->Set(size+1);
    fSelections->AddAt(kc,size);
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCollider::GetSelectionFlag() const
{
// Return the value of the selection flag for the total event.
// When the event passed the selection criteria as specified via
// SelectEvent() the value 1 is returned, otherwise the value 0 is returned.
 return fSelect;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCollider::IsSelected()
{
// Check whether the generated (sub)event contains one of the particles
// specified in the selection list via SelectEvent().
// If this is the case or when no selection list is present, the value 1
// will be returned, indicating the event is selected to be kept.
// Otherwise the value 0 will be returned.

 if (!fSelections) return 1; 

 Int_t nsel=fSelections->GetSize();
 Int_t npart=GetN();
 Int_t kf,kc;

 Int_t select=0;
 for (Int_t jpart=1; jpart<=npart; jpart++)
 {
  kf=GetK(jpart,2);
  kc=Pycomp(kf);
  for (Int_t i=0; i<nsel; i++)
  {
   if (kc==fSelections->At(i))
   {
    select=1;
    break;
   }
  }
  if (select) break;
 }
 return select;
}
///////////////////////////////////////////////////////////////////////////
void AliCollider::SetSpectatorPmin(Float_t pmin)
{
// Set minimal momentum in GeV/c for spectator tracks to be stored.
// Spectator tracks with a momentum below this threshold will not be stored
// in the (output) event structure.
// This facility allows to minimise the output file size.
// Note that when the user wants to boost the event into another reference
// frame these spectator tracks might have got momenta above the threshold.
// However, when the spectator tracks were not stored in the event structure
// in the original frame, there is no way to retreive them anymore. 
 fSpecpmin=pmin;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCollider::GetSpectatorPmin() const
{
// Provide the minimal spectator momentum in GeV/c.
 return fSpecpmin;
}
///////////////////////////////////////////////////////////////////////////
TString AliCollider::GetPyname(Int_t kf)
{
// Provide the correctly truncated Pythia particle name for PGD code kf
//
// The TPythia6::Pyname returned name is copied into a TString and truncated
// at the first blank to prevent funny trailing characters due to incorrect
// stripping of empty characters in TPythia6::Pyname.
// The truncation at the first blank is allowed due to the Pythia convention
// that particle names never contain blanks.
 char name[16];
 TString sname;
 Pyname(kf,name);
 sname=name[0];
 for (Int_t i=1; i<16; i++)
 {
  if (name[i]==' ') break;
  sname=sname+name[i];
 }
 return sname;
}
///////////////////////////////////////////////////////////////////////////
