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

/*
$Log$
Revision 1.3  1999/09/29 09:24:28  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////
// Class AliVertex
// Creation and investigation of an AliVertex.
// An AliVertex can be constructed by adding AliTracks and/or AliJets.
//
// Note : Also (secondary) vertices can be added to a vertex.
//
// Coding example to make 3 vertices v1, v2 and v3.
// ------------------------------------------------
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
//        AliVertex v1(5);
//
//        v1.Add(t1);
//        v1.Add(t2);
//        v1.Add(t3);
//        v1.Add(t4);
//
//        Float_t r1[3]={2.4,0.1,-8.5};
//        v1.SetPosition(r1,"car");
//
//        AliVertex v2(2);
//        v2.Add(t5);
//        v2.Add(t6);
//        v2.Add(t7);
//
//        Float_t r2[3]={1.6,-3.2,5.7};
//        v2.SetPosition(r2,"car");
//
//        AliVertex v3;
//
//        v3.Add(j1);
//        v3.Add(j2);
//
//        Float_t r3[3]={6.2,4.8,1.3};
//        v3.SetPosition(r3,"car");
//
//        v1.Info("sph");
//        v2.ListAll();
//        v3.List("cyl");
//
//        Float_t e1=v1.GetEnergy();
//        Ali3Vector p1=v1.Get3Momentum();
//        Float_t loc[3];
//        v1.GetPosition(loc,"sph");
//        AliPosition r=v2.GetPosition();
//        r.Info(); 
//        Int_t nt=v2.GetNtracks();
//        AliTrack* tv=v2.GetTrack(1); // Access track number 1 of Vertex v2
//
// Specify the vertices v2 and v3 as secondary vertices of v1
//
//        v1.Add(v2);
//        v1.Add(v3);
//
//        v1.List();
//
//        Int_t nv=v1.GetNvtx();
//        AliVertex* vx=v1.GetVertex(1); // Access 1st secondary vertex of v1
//        Float_t e=vx->GetEnergy();
//
//        Float_t M=v1.GetInvmass(); 
//
// Reconstruct Vertex v1 from scratch
//
//        v1.Reset();
//        v1.SetNvmax(25); // Increase initial no. of sec. vertices
//        v1.Add(t3);
//        v1.Add(t7);
//        v1.Add(j2);
//        Float_t pos[3]={7,9,4};
//        v1.SetPosition(pos,"car");
//
// Note : All quantities are in GeV, GeV/c or GeV/c**2
//
//--- Author: Nick van Eijndhoven 04-apr-1998 UU-SAP Utrecht
//- Modified: NvE 08-apr-1999 UU-SAP Utrecht to inherit from AliJet
///////////////////////////////////////////////////////////////////////////

#include "AliVertex.h"
 
ClassImp(AliVertex) // Class implementation to enable ROOT I/O
 
AliVertex::AliVertex()
{
// Default constructor.
// All variables initialised to 0.
// Initial maximum number of tracks is set to the default value.
// Initial maximum number of sec. vertices is set to the default value.
 fNvmax=0;
 fVertices=0;
 Reset();
 SetNtinit();
 SetNvmax();
}
///////////////////////////////////////////////////////////////////////////
AliVertex::AliVertex(Int_t n)
{
// Create a vertex to hold initially a maximum of n tracks
// All variables initialised to 0
 fNvmax=0;
 fVertices=0;
 Reset();
 if (n > 0)
 {
  SetNtinit(n);
 }
 else
 {
  cout << endl;
  cout << " *AliVertex* Initial max. number of tracks entered : " << n << endl;
  cout << " This is invalid. Default initial maximum will be used." << endl;
  cout << endl;
  SetNtinit();
 }
 SetNvmax();
}
///////////////////////////////////////////////////////////////////////////
AliVertex::~AliVertex()
{
// Default destructor
 if (fVertices) delete fVertices;
 fVertices=0;
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::SetNvmax(Int_t n)
{
// Set the initial maximum number of (secondary) vertices
 if (n > 0)
 {
  fNvmax=n;
 }
 else
 {
  fNvmax=1;
 }
 if (fVertices) delete fVertices;
 fVertices=new TObjArray(fNvmax);
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::Reset()
{
// Reset all variables to 0
// The max. number of tracks is set to the initial value again
// The max. number of vertices is set to the default value again

 AliJet::Reset();

 fNvtx=0;
 if (fNvmax>0) SetNvmax(fNvmax);
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::Add(AliJet& j)
{
// Add the tracks of a jet to the vertex
 AliTrack* tj;
 for (Int_t i=1; i<=j.GetNtracks(); i++)
 {
  tj=j.GetTrack(i);
  AliJet::Add(tj);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::Add(AliVertex& v)
{
// Add a (secondary) vertex to the current vertex.
// In case the maximum number of (secondary) vertices has been reached,
// the array space will be extended automatically
//
// Note : The 4-momentum of the current (primary) vertex
//        is updated automatically, but the track connecting
//        both vertices has to be entered separately by the user.
//
 if (fNvtx == fNvmax) // Check if maximum vertex number is reached
 {
  fNvmax++;
  fVertices->Expand(fNvmax);
 }
 
 // Update 4-momentum for current vertex
 fNvtx++;
 fVertices->Add(&v);
 (Ali4Vector)(*this)+=v;
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::Info(TString f)
{
// Provide vertex information within the coordinate frame f
 cout << " *AliVertex::Info* Invmass : " << GetInvmass()
      << " Charge : " << GetCharge() << " Momentum : " << GetMomentum()
      << " Ntracks : " << GetNtracks() << " Nvertices : " << fNvtx << endl;
 cout << " ";
 Ali4Vector::Info(f);
 cout << "  Position";
 AliPosition::Info(f); 
} 
///////////////////////////////////////////////////////////////////////////
void AliVertex::List(TString f)
{
// Provide primary track and sec. vertex information within the coordinate frame f

 Info(f); // Information of the current vertex

 // The tracks of this vertex
 AliTrack* t; 
 for (Int_t it=1; it<=GetNtracks(); it++)
 {
  t=GetTrack(it);
  if (t)
  {
   cout << "  ---Track no. " << it << endl;
   cout << " ";
   t->Info(f); 
  }
  else
  {
   cout << " *AliVertex::List* Error : No track present." << endl; 
  }
 }

 // The secondary vertices of this vertex
 AliVertex* v; 
 for (Int_t iv=1; iv<=GetNvertices(); iv++)
 {
  v=GetVertex(iv);
  if (v)
  {
   cout << "  ---Level 1 sec. vertex no. " << iv << endl;
   cout << " ";
   v->Info(f); 
  }
  else
  {
   cout << " *AliVertex::List* Error : No sec. vertex present." << endl; 
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
void AliVertex::ListAll(TString f)
{
// Provide complete (sec) vertex and (decay) track info within the coordinate frame f

 Info(f); // Information of the current vertex

 // The tracks of this vertex
 AliTrack* t; 
 for (Int_t it=1; it<=GetNtracks(); it++)
 {
  t=GetTrack(it);
  if (t)
  {
   cout << "  ---Track no. " << it << endl;
   cout << " ";
   t->ListAll(f); 
  }
  else
  {
   cout << " *AliVertex::ListAll* Error : No track present." << endl; 
  }
 }

 AliVertex* v=this;
 Dump(v,1,f); // Information of all sec. vertices
}
//////////////////////////////////////////////////////////////////////////
void AliVertex::Dump(AliVertex* v,Int_t n,TString f)
{
// Recursively provide the info of all secondary vertices of this vertex
 AliVertex* vs; 
 for (Int_t iv=1; iv<=v->GetNvertices(); iv++)
 {
  vs=v->GetVertex(iv);
  if (vs)
  {
   cout << "  ---Level " << n << " sec. vertex no. " << iv << endl;
   cout << " ";
   vs->Info(f); 

   // The tracks of this vertex
   AliTrack* t; 
   for (Int_t it=1; it<=vs->GetNtracks(); it++)
   {
    t=vs->GetTrack(it);
    if (t)
    {
     cout << "  ---Track no. " << it << endl;
     cout << " ";
     t->ListAll(f); 
    }
    else
    {
     cout << " *AliVertex::Dump* Error : No track present." << endl; 
    }
   }

   // Go for next sec. vertex level of this sec. vertex recursively
   Dump(vs,n+1,f);
  }
  else
  {
   cout << " *AliVertex::Dump* Error : No sec. vertex present." << endl; 
  }
 }
} 
//////////////////////////////////////////////////////////////////////////
Int_t AliVertex::GetNvertices()
{
// Return the current number of (secondary) vertices
 return fNvtx;
}
///////////////////////////////////////////////////////////////////////////
AliVertex* AliVertex::GetVertex(Int_t i)
{
// Return the i-th (secondary) vertex of the current vertex
 return (AliVertex*)fVertices->At(i-1);
}
///////////////////////////////////////////////////////////////////////////
