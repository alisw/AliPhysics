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

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class AliEvent
// Creation and investigation of an Alice physics event.
// An AliEvent can be constructed by adding AliTracks, Alivertices
// and/or AliJets.
//
// The basic functionality of AliEvent is identical to the one of AliVertex.
// So, an AliEvent may be regarded as the primary vertex with some
// additional functionality compared to AliVertex.
//
// Coding example to make an event consisting of a primary vertex
// and 2 secondary vertices.
// --------------------------------------------------------------
// v1 contains the tracks 1,2,3 and 4
// v2 contains the tracks 5,6 and 7
// v3 contains the jets 1 and 2
//
//        AliTrack t1,t2,t3,t4,t5,t6,t7;
//         ...
//         ... // code to fill the track data
//         ...
//
//        AliJet j1,j2;
//         ...
//         ... // code to fill the jet data
//         ...
//
//        AliEvent evt(5);
//
//        evt.AddTrack(t1);
//        evt.AddTrack(t2);
//        evt.AddTrack(t3);
//        evt.AddTrack(t4);
//
//        Float_t r0[3]={2.4,0.1,-8.5};
//        evt.SetPosition(r0,"car");
//
//        AliVertex v1(2);
//        v1.AddTrack(t5);
//        v1.AddTrack(t6);
//        v1.AddTrack(t7);
//
//        Float_t r1[3]={1.6,-3.2,5.7};
//        v1.SetPosition(r1,"car");
//
//        AliVertex v2;
//
//        v2.AddJet(j1);
//        v2.AddJet(j2);
//
//        Float_t r2[3]={6.2,4.8,1.3};
//        v2.SetPosition(r2,"car");
//
//        evt.Info("sph");
//        v1.ListAll();
//        v2.List("cyl");
//
//        Float_t etot=evt.GetEnergy();
//        Ali3Vector ptot=evt.Get3Momentum();
//        Float_t loc[3];
//        evt.GetPosition(loc,"sph");
//        AliPosition r=v1.GetPosition();
//        r.Info(); 
//        Int_t nt=v2.GetNtracks();
//        AliTrack* tv=v2.GetTrack(1); // Access track number 1 of Vertex v2
//
// Specify the vertices v2 and v3 as secondary vertices of the primary
//
//        evt.AddVertex(v2);
//        evt.AddVertex(v3);
//
//        evt.List();
//
//        Int_t nv=evt.GetNvtx();
//        AliVertex* vx=evt.GetVertex(1); // Access 1st secondary vertex
//        Float_t e=vx->GetEnergy();
//
//        Float_t M=evt.GetInvmass(); 
//
// Reconstruct the event from scratch
//
//        evt.Reset();
//        evt.SetNvmax(25); // Increase initial no. of sec. vertices
//        evt.AddTrack(t3);
//        evt.AddTrack(t7);
//        evt.AddJet(j2);
//        Float_t pos[3]={7,9,4};
//        evt.SetPosition(pos,"car");
//
// Note : All quantities are in GeV, GeV/c or GeV/c**2
//
//--- Author: Nick van Eijndhoven 27-may-2001 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliEvent.h"
 
ClassImp(AliEvent) // Class implementation to enable ROOT I/O
 
AliEvent::AliEvent()
{
// Default constructor.
// All variables initialised to default values.
 fDaytime.Set();
 fRun=0;
 fEvent=0;
 AliVertex::AliVertex();
}
///////////////////////////////////////////////////////////////////////////
AliEvent::AliEvent(Int_t n)
{
// Create an event to hold initially a maximum of n tracks
// All variables initialised to default values
 fDaytime.Set();
 fRun=0;
 fEvent=0;
 if (n > 0)
 {
  AliVertex::AliVertex(n);
 }
 else
 {
  cout << endl;
  cout << " *AliEvent* Initial max. number of tracks entered : " << n << endl;
  cout << " This is invalid. Default initial maximum will be used." << endl;
  cout << endl;
  AliVertex::AliVertex();
 }
}
///////////////////////////////////////////////////////////////////////////
AliEvent::~AliEvent()
{
// Default destructor
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::Reset()
{
// Reset all variables to default values
// The max. number of tracks is set to the initial value again
// The max. number of vertices is set to the default value again
 fDaytime.Set();
 fRun=0;
 fEvent=0;

 AliVertex::Reset();
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
void AliEvent::HeaderInfo()
{
// Provide event header information
 Int_t date=fDaytime.GetDate();
 Int_t time=fDaytime.GetTime();

 Int_t year=date/10000;
 Int_t month=(date%10000)/100;
 Int_t day=date%100;
 Int_t hh=time/10000;
 Int_t mm=(time%10000)/100;
 Int_t ss=time%100;

 char* c[12]={"jan","feb","mar","apr","may","jun",
              "jul","aug","sep","oct","nov","dec"};

 cout << " *AliEvent::Info* Run : " << fRun << " Event : " << fEvent;
 cout.fill('0');
 cout << " Date : " << setw(2) << day << "-" << c[month-1] << "-" << year
      << " Time : " << setw(2) << hh << ":" << setw(2) << mm << ":" << setw(2) << ss << endl;
 cout.fill(' ');
}
///////////////////////////////////////////////////////////////////////////
void AliEvent::Info(TString f)
{
// Provide event information within the coordinate frame f
 HeaderInfo();
 AliVertex::Info(f);
} 
///////////////////////////////////////////////////////////////////////////

