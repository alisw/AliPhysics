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

// $Id: AliEvent.cxx,v 1.11 2003/02/25 12:36:28 nick Exp $

///////////////////////////////////////////////////////////////////////////
// Class AliEvent
// Creation and investigation of an Alice physics event.
// An AliEvent can be constructed by adding AliTracks, Alivertices, AliJets
// and/or AliCalorimeters.
//
// The basic functionality of AliEvent is identical to the one of AliVertex.
// So, an AliEvent may be used as the primary vertex with some additional
// functionality compared to AliVertex.
//
// To provide maximal flexibility to the user, the two modes of track/jet/vertex
// storage as described in AliJet and AliVertex can be used.
// In addition an identical structure is provided for the storage of AliCalorimeter
// objects, which can be selected by means of the memberfunction SetCalCopy().
//
// a) SetCalCopy(0) (which is the default).
//    Only the pointers of the 'added' calorimeters are stored.
//    This mode is typically used by making cal. studies based on a fixed set
//    of calorimeters which stays under user control or is kept on an external
//    file/tree. 
//    In this way the AliEvent just represents a 'logical structure' for the
//    physics analysis.
//
//    Note :
//    Modifications made to the original calorimeters also affect the AliCalorimeter
//    objects which are stored in the AliEvent. 
//
// b) SetCalCopy(1).
//    Of every 'added' calorimeter a private copy will be made of which the pointer
//    will be stored.
//    In this way the AliEvent represents an entity on its own and modifications
//    made to the original calorimeters do not affect the AliCalorimeter objects
//    which are stored in the AliEvent. 
//    This mode will allow 'adding' many different AliCalorimeters into an AliEvent by
//    creating only one AliCalorimeter instance in the main programme and using the
//    AliCalorimeter::Reset() and AliCalorimeter parameter setting memberfunctions.
//
// See also the documentation provided for the memberfunction SetOwner(). 
//
// Coding example to make an event consisting of a primary vertex,
// 2 secondary vertices and a calorimeter.
// --------------------------------------------------------------
// vp contains the tracks 1,2,3 and 4 (primary vertex)
// v1 contains the tracks 5,6 and 7   (sec. vertex)
// v2 contains the jets 1 and 2       (sec. vertex)
//
//        AliEvent evt;
//
// Specify the event object as the repository of all objects
// for the event building and physics analysis.
// 
//        evt.SetCalCopy(1);
//        evt.SetTrackCopy(1);
//
// Fill the event structure with the basic objects
// 
//        AliCalorimeter emcal;
//         ...
//         ... // code to fill the calorimeter data
//         ...
//
//        evt.AddCalorimeter(emcal);
//
//        AliTrack* tx=new AliTrack();
//        for (Int_t i=0; i<10; i++)
//        {
//         ...
//         ... // code to fill the track data
//         ...
//         evt.AddTrack(tx);
//         tx->Reset(); 
//        }
//
//        if (tx)
//        {
//         delete tx;
//         tx=0;
//        }
//
// Build the event structure (vertices, jets, ...) for physics analysis
// based on the basic objects from the event repository.
//
//        AliJet j1,j2;
//        for (Int_t i=0; i<evt.GetNtracks(); i++)
//        {
//         tx=evt.GetTrack(i);
//         ...
//         ... // code to fill the jet data
//         ...
//        }
//
//        AliVertex vp;
//        tx=evt.GetTrack(1);
//        vp.AddTrack(tx);
//        tx=evt.GetTrack(2);
//        vp.AddTrack(tx);
//        tx=evt.GetTrack(3);
//        vp.AddTrack(tx);
//        tx=evt.GetTrack(4);
//        vp.AddTrack(tx);
//
//        Float_t rp[3]={2.4,0.1,-8.5};
//        vp.SetPosition(rp,"car");
//
//        AliVertex v1;
//        tx=evt.GetTrack(5);
//        v1.AddTrack(tx);
//        tx=evt.GetTrack(6);
//        v1.AddTrack(tx);
//        tx=evt.GetTrack(7);
//        v1.AddTrack(tx);
//
//        Float_t r1[3]={1.6,-3.2,5.7};
//        v1.SetPosition(r1,"car");
//
//
//        AliVertex v2;
//        v2.SetJetCopy(1);
//        v2.AddJet(j1);
//        v2.AddJet(j2);
//
//        Float_t r2[3]={6.2,4.8,1.3};
//        v2.SetPosition(r2,"car");
//
// Specify the vertices v1 and v2 as secondary vertices of the primary
//
//        vp.SetVertexCopy(1);
//        vp.AddVertex(v1);
//        vp.AddVertex(v2);
//
// Enter the physics structures into the event
//        evt.SetVertexCopy(1);
//        evt.AddVertex(vp,0);
//
// The jets j1 and j2 are already available via sec. vertex v2,
// but can be made available also from the event itself if desired.
//        AliJet* jx;
//        jx=v2.GetJet(1);
//        evt.AddJet(jx,0); 
//        jx=v2.GetJet(2);
//        evt.AddJet(jx,0); 
// 
//        evt.Data("sph");
//        v1.ListAll();
//        v2.List("cyl");
//
//        Float_t etot=evt.GetEnergy();
//        Ali3Vector ptot=evt.Get3Momentum();
//        Float_t loc[3];
//        evt.GetPosition(loc,"sph");
//        AliPosition r=v1.GetPosition();
//        r.Data(); 
//        Int_t nt=v2.GetNtracks();
//        AliTrack* tv=v2.GetTrack(1); // Access track number 1 of Vertex v2
//
//        evt.List();
//
//        Int_t nv=evt.GetNvtx();
//        AliVertex* vx=evt.GetVertex(1); // Access primary vertex
//        Float_t e=vx->GetEnergy();
//
//        Float_t M=evt.GetInvmass(); 
//
// Reconstruct the event from scratch
//
//        evt.Reset();
//        evt.SetNvmax(25); // Increase initial no. of sec. vertices
//        ...
//        ... // code to create tracks etc... 
//        ...
//
// Note : All quantities are in GeV, GeV/c or GeV/c**2
//
//--- Author: Nick van Eijndhoven 27-may-2001 UU-SAP Utrecht
//- Modified: NvE $Date: 2003/02/25 12:36:28 $ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliEvent.h"
#include "Riostream.h"
 
ClassImp(AliEvent) // Class implementation to enable ROOT I/O
 
AliEvent::AliEvent() : AliVertex()
{
// Default constructor.
// All variables initialised to default values.
 fDaytime.Set();
 fRun=0;
 fEvent=0;
 fAproj=0;
 fZproj=0;
 fPnucProj=0;
 fIdProj=0;
 fAtarg=0;
 fZtarg=0;
 fPnucTarg=0;
 fIdTarg=0;
 fNcals=0;
 fCalorimeters=0;
 fCalCopy=0;
}
///////////////////////////////////////////////////////////////////////////
AliEvent::AliEvent(Int_t n) : AliVertex(n)
{
// Create an event to hold initially a maximum of n tracks
// All variables initialised to default values
 if (n<=0)
 {
  cout << " *** This AliVertex initialisation was invoked via the AliEvent ctor." << endl;
 }
 fDaytime.Set();
 fRun=0;
 fEvent=0;
 fAproj=0;
 fZproj=0;
 fPnucProj=0;
 fIdProj=0;
 fAtarg=0;
 fZtarg=0;
 fPnucTarg=0;
 fIdTarg=0;
 fNcals=0;
 fCalorimeters=0;
 fCalCopy=0;
}
///////////////////////////////////////////////////////////////////////////
AliEvent::~AliEvent()
{
// Default destructor
 if (fCalorimeters)
 {
  delete fCalorimeters;
  fCalorimeters=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliEvent::AliEvent(AliEvent& evt) : AliVertex(evt)
{
// Copy constructor.
 fDaytime=evt.fDaytime;
 fRun=evt.fRun;
 fEvent=evt.fEvent;
 fAproj=evt.fAproj;
 fZproj=evt.fZproj;
 fPnucProj=evt.fPnucProj;
 fIdProj=evt.fIdProj;
 fAtarg=evt.fAtarg;
 fZtarg=evt.fZtarg;
 fPnucTarg=evt.fPnucTarg;
 fIdTarg=evt.fIdTarg;
 fNcals=evt.fNcals;
 fCalCopy=evt.fCalCopy;

 fCalorimeters=0;
 if (fNcals)
 {
  fCalorimeters=new TObjArray(fNcals);
  if (fCalCopy) fCalorimeters->SetOwner();
  for (Int_t i=1; i<=fNcals; i++)
  {
   AliCalorimeter* cal=evt.GetCalorimeter(i);
   if (cal)
   {
    if (fCalCopy)
    {
     fCalorimeters->Add(new AliCalorimeter(*cal));
    }
    else
    {
     fCalorimeters->Add(cal);
    }
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::Reset()
{
// Reset all variables to default values
// The max. number of tracks is set to the initial value again
// The max. number of vertices is set to the default value again
// Note : The CalCopy mode is maintained as it was set by the user before.

 AliVertex::Reset();

 fDaytime.Set();
 fRun=0;
 fEvent=0;
 fAproj=0;
 fZproj=0;
 fPnucProj=0;
 fIdProj=0;
 fAtarg=0;
 fZtarg=0;
 fPnucTarg=0;
 fIdTarg=0;

 fNcals=0;
 if (fCalorimeters)
 {
  delete fCalorimeters;
  fCalorimeters=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::SetOwner(Bool_t own)
{
// Set ownership of all added objects. 
// The default parameter is own=kTRUE.
//
// Invokation of this memberfunction also sets all the copy modes
// (e.g. TrackCopy & co.) according to the value of own.
//
// This function (with own=kTRUE) is particularly useful when reading data
// from a tree/file, since Reset() will then actually remove all the
// added objects from memory irrespective of the copy mode settings
// during the tree/file creation process. In this way it provides a nice way
// of preventing possible memory leaks in the reading/analysis process.
//
// In addition this memberfunction can also be used as a shortcut to set all
// copy modes in one go during a tree/file creation process.
// However, in this case the user has to take care to only set/change the
// ownership (and copy mode) for empty objects (e.g. newly created objects
// or after invokation of the Reset() memberfunction) otherwise it will
// very likely result in inconsistent destructor behaviour.

 Int_t mode=1;
 if (!own) mode=0;
 if (fCalorimeters) fCalorimeters->SetOwner(own);
 fCalCopy=mode;

 AliVertex::SetOwner(own);
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::SetDayTime(TDatime& stamp)
{
// Set the date and time stamp for this event
 fDaytime=stamp;
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::SetRunNumber(Int_t run)
{
// Set the run number for this event
 fRun=run;
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::SetEventNumber(Int_t evt)
{
// Set the event number for this event
 fEvent=evt;
}
///////////////////////////////////////////////////////////////////////////
TDatime AliEvent::GetDayTime()
{
// Provide the date and time stamp for this event
 return fDaytime;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetRunNumber()
{
// Provide the run number for this event
 return fRun;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetEventNumber()
{
// Provide the event number for this event
 return fEvent;
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::SetProjectile(Int_t a,Int_t z,Double_t pnuc,Int_t id)
{
// Set the projectile A, Z, momentum per nucleon and user defined particle ID.
// By default the particle ID is set to zero.
 fAproj=a;
 fZproj=z;
 fPnucProj=pnuc;
 fIdProj=id;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetProjectileA()
{
// Provide the projectile A value.
 return fAproj;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetProjectileZ()
{
// Provide the projectile Z value.
 return fZproj;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliEvent::GetProjectilePnuc()
{
// Provide the projectile momentum value per nucleon.
 return fPnucProj;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetProjectileId()
{
// Provide the user defined particle ID of the projectile.
 return fIdProj;
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::SetTarget(Int_t a,Int_t z,Double_t pnuc,Int_t id)
{
// Set the target A, Z, momentum per nucleon and user defined particle ID.
// By default the particle ID is set to zero.
 fAtarg=a;
 fZtarg=z;
 fPnucTarg=pnuc;
 fIdTarg=id;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetTargetA()
{
// Provide the target A value.
 return fAtarg;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetTargetZ()
{
// Provide the target Z value.
 return fZtarg;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliEvent::GetTargetPnuc()
{
// Provide the target momentum value per nucleon.
 return fPnucTarg;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetTargetId()
{
// Provide the user defined particle ID of the target.
 return fIdTarg;
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::HeaderData()
{
// Provide event header information
 cout << " *AliEvent::Data* Run : " << fRun << " Event : " << fEvent
      << " Date : " << fDaytime.AsString() << endl;

 ShowCalorimeters();
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::Data(TString f)
{
// Provide event information within the coordinate frame f
 HeaderData();
 AliVertex::Data(f);
} 
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetNcalorimeters()
{
// Provide the number of stored calorimeter systems
 return fNcals;
} 
///////////////////////////////////////////////////////////////////////////
void AliEvent::AddCalorimeter(AliCalorimeter& c)
{
// Add a calorimeter system to the event
 if (!fCalorimeters)
 {
  fCalorimeters=new TObjArray();
  if (fCalCopy) fCalorimeters->SetOwner();
 }
 
 // Add the calorimeter system to this event
 fNcals++;
 if (fCalCopy)
 {
  fCalorimeters->Add(new AliCalorimeter(c));
 }
 else
 {
  fCalorimeters->Add(&c);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::SetCalCopy(Int_t j)
{
// (De)activate the creation of private copies of the added calorimeters.
// j=0 ==> No private copies are made; pointers of original cals. are stored.
// j=1 ==> Private copies of the cals. are made and these pointers are stored.
//
// Note : Once the storage contains pointer(s) to AliCalorimeter(s) one cannot
//        change the CalCopy mode anymore.
//        To change the CalCopy mode for an existing AliEvent containing
//        calorimeters one first has to invoke Reset().
 if (!fCalorimeters)
 {
  if (j==0 || j==1)
  {
   fCalCopy=j;
  }
  else
  {
   cout << "*AliEvent::SetCalCopy* Invalid argument : " << j << endl;
  }
 }
 else
 {
  cout << "*AliEvent::SetCalCopy* Storage already contained calorimeters."
       << "  ==> CalCopy mode not changed." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliEvent::GetCalCopy()
{
// Provide value of the CalCopy mode.
// 0 ==> No private copies are made; pointers of original cals. are stored.
// 1 ==> Private copies of the cals. are made and these pointers are stored.
 return fCalCopy;
}
///////////////////////////////////////////////////////////////////////////
AliCalorimeter* AliEvent::GetCalorimeter(Int_t i)
{
// Return the i-th calorimeter of this event
 if (!fCalorimeters)
 {
  cout << " *AliEvent::GetCalorimeter* No calorimeters present." << endl;
  return 0;
 }
 else
 {
  if (i<=0 || i>fNcals)
  {
   cout << " *AliEvent::GetCalorimeter* Invalid argument i : " << i
        << " Ncals = " << fNcals << endl;
   return 0;
  }
  else
  {
   return (AliCalorimeter*)fCalorimeters->At(i-1);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
AliCalorimeter* AliEvent::GetCalorimeter(TString name)
{
// Return the calorimeter with name tag "name"
 if (!fCalorimeters)
 {
  cout << " *AliEvent::GetCalorimeter* No calorimeters present." << endl;
  return 0;
 }
 else
 {
  AliCalorimeter* cx;
  TString s;
  for (Int_t i=0; i<fNcals; i++)
  {
   cx=(AliCalorimeter*)fCalorimeters->At(i);
   if (cx)
   {
    s=cx->GetName();
    if (s == name) return cx;
   }
  }

  return 0; // No matching name found
 }
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::ShowCalorimeters()
{
// Provide an overview of the available calorimeter systems.
 if (fNcals>0)
 {
  cout << " The following " << fNcals << " calorimeter systems are available :" << endl; 
  for (Int_t i=1; i<=fNcals; i++)
  {
   AliCalorimeter* cal=GetCalorimeter(i);
   if (cal) cout << " Calorimeter number : " << i << " Name : " << (cal->GetName()).Data() << endl;
  }
 }
 else
 {
  cout << " No calorimeters present for this event." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////

