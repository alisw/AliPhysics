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
// Class AliVertex
// Creation and investigation of an AliVertex.
// An AliVertex can be constructed by adding AliTracks and/or AliJets.
//
// Note : Also (secondary) vertices can be added to a vertex.
//
// To provide maximal flexibility to the user, two modes of vertex storage
// are provided by means of the memberfunction SetVertexCopy().
// The same holds for the storage of jets via SetJetCopy().
//
// a) SetVertexCopy(0) (which is the default).
//    Only the pointers of the 'added' vertices are stored.
//    This mode is typically used by making vertex studies based on a fixed list
//    of vertices which stays under user control or is contained for instance
//    in an AliEvent.  
//    In this way the AliVertex just represents a 'logical structure' for the
//    physics analysis which can be embedded in e.g. an AliEvent or AliVertex.
//
//    Note :
//    Modifications made to the original vertices also affect the AliVertex objects
//    which are stored.
// 
// b) SetVertexCopy(1).
//    Of every 'added' vertex a private copy will be made of which the pointer
//    will be stored.
//    In this way the AliVertex represents an entity on its own and modifications
//    made to the original vertices do not affect the AliVertex objects which are
//    stored. 
//    This mode will allow 'adding' many different AliVertex objects by
//    creating only one AliVertex instance in the main programme and using the
//    AliVertex::Reset, AliVertex::AddTrack and parameter setting memberfunctions.
//
// Coding example to make 3 vertices v1, v2 and v3.
// ------------------------------------------------
// v1 contains the tracks 1,2,3 and 4
// v2 contains many different tracks
// v3 contains the jets 1 and 2
//
//        AliTrack t1,t2,t3,t4;
//         ...
//         ... // code to fill the track data
//         ...
//
//        AliJet j1,j2;
//         ...
//         ... // code to fill the jet data
//         ...
//
//        AliVertex v1;
//        v1.SetVertexCopy(1);
//
//        v1.AddTrack(t1);
//        v1.AddTrack(t2);
//        v1.AddTrack(t3);
//        v1.AddTrack(t4);
//
//        Float_t r1[3]={2.4,0.1,-8.5};
//        v1.SetPosition(r1,"car");
//
//        AliVertex v2;
//        v2.SetTrackCopy(1);
//
//        AliTrack* tx=new AliTrack();
//        for (Int_t i=0; i<10; i++)
//        {
//         ...
//         ... // code to fill the track data
//         ...
//         v2.AddTrack(tx);
//         tx->Reset(); 
//        }
//
//        Float_t r2[3]={1.6,-3.2,5.7};
//        v2.SetPosition(r2,"car");
//
//        AliVertex v3;
//
//        v3.AddJet(j1);
//        v3.AddJet(j2);
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
//        v1.AddVertex(v2);
//        v1.AddVertex(v3);
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
//        v1.AddTrack(t3);
//        v1.AddTrack(t4);
//        v1.AddJet(j2);
//        Float_t pos[3]={7,9,4};
//        v1.SetPosition(pos,"car");
//
// Note : All quantities are in GeV, GeV/c or GeV/c**2
//
//--- Author: Nick van Eijndhoven 04-apr-1998 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
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
 fConnects=0;
 fVertexCopy=0;
 fNjmax=0;
 fJets=0;
 fJetCopy=0;
 Reset();
 SetNtinit();
 SetNvmax();
 SetNjmax();
}
///////////////////////////////////////////////////////////////////////////
AliVertex::AliVertex(Int_t n)
{
// Create a vertex to hold initially a maximum of n tracks
// All variables initialised to 0
 fNvmax=0;
 fVertices=0;
 fConnects=0;
 fVertexCopy=0;
 fNjmax=0;
 fJets=0;
 fJetCopy=0;
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
 SetNjmax();
}
///////////////////////////////////////////////////////////////////////////
AliVertex::~AliVertex()
{
// Default destructor
 if (fVertices)
 {
  if (fVertexCopy) fVertices->Delete();
  delete fVertices;
  fVertices=0;
 }
 if (fConnects)
 {
  fConnects->Delete();
  delete fConnects;
  fConnects=0;
 }
 if (fJets)
 {
  if (fJetCopy) fJets->Delete();
  delete fJets;
  fJets=0;
 }
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
 if (fVertices)
 {
  if (fVertexCopy) fVertices->Delete();
  delete fVertices;
  fVertices=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::SetNjmax(Int_t n)
{
// Set the initial maximum number of jets
 if (n > 0)
 {
  fNjmax=n;
 }
 else
 {
  fNjmax=1;
 }
 if (fJets)
 {
  if (fJetCopy) fJets->Delete();
  delete fJets;
  fJets=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::Reset()
{
// Reset all variables to 0 and reset all stored vertex and jet lists.
// The max. number of tracks is set to the initial value again
// The max. number of vertices is set to the default value again
// The max. number of jets is set to the default value again

 AliJet::Reset();

 Double_t a[3]={0,0,0};
 SetPosition(a,"sph");
 SetPositionErrors(a,"car");

 fNvtx=0;
 if (fNvmax>0) SetNvmax(fNvmax);
 if (fConnects)
 {
  fConnects->Delete();
  delete fConnects;
  fConnects=0;
 }

 fNjets=0;
 if (fNjmax>0) SetNjmax(fNjmax);
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::ResetVertices()
{
// Reset the stored vertex list and delete all connecting tracks which
// were generated automatically via connect=1 in AddVertex().
// The max. number of vertices is set to the default value again.
// All physics quantities are updated according to the removal of the
// connecting tracks.
 AliTrack* t;
 if (fConnects)
 {
  for (Int_t i=0; i<=fConnects->GetLast(); i++)
  {
   t=(AliTrack*)fConnects->At(i);
   AliTrack* test=(AliTrack*)fTracks->Remove(t);
   if (test)
   {
    fNtrk--;
    (Ali4Vector&)(*this)-=(Ali4Vector&)(*t);
    fQ-=t->GetCharge();
   }
  }
  fTracks->Compress();
 }

 fNvtx=0;
 if (fNvmax>0) SetNvmax(fNvmax);
 if (fConnects)
 {
  fConnects->Delete();
  delete fConnects;
  fConnects=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::AddJet(AliJet& j,Int_t tracks)
{
// Add a jet (and its tracks) to the vertex
// In case the maximum number of jets has been reached,
// the array space will be extended automatically
//
// Note : By default the tracks of the jet are added to the current (primary)
//        vertex.
//        The automatic addition of the tracks of the jet can be suppressed
//        by specifying tracks=0. In this case only the AliJet object will
//        be stored according to the mode specified by SetJetCopy().
//        The latter will enable jet studies based on a fixed list of tracks
//        as contained e.g. in an AliVertex or AliEvent. 
 if (!fJets) fJets=new TObjArray(fNjmax);
 if (fNjets == fNjmax) // Check if maximum jet number is reached
 {
  fNjmax++;
  fJets->Expand(fNjmax);
 }

 // Add the jet to the list 
 fNjets++;
 if (fJetCopy)
 {
  fJets->Add(j.Clone());
 }
 else
 {
  fJets->Add(&j);
 }

 // Add the tracks of the jet to this vertex
 if (tracks)
 {
  AliTrack* tj;
  for (Int_t i=1; i<=j.GetNtracks(); i++)
  {
   tj=j.GetTrack(i);
   AddTrack(tj);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::AddVertex(AliVertex& v,Int_t connect)
{
// Add a (secondary) vertex to the current vertex.
// In case the maximum number of (secondary) vertices has been reached,
// the array space will be extended automatically
//
// Note : By default the 4-momentum and charge of the current (primary) vertex
//        are updated by automatically creating the track connecting
//        both vertices. The track parameters are taken from the
//        4-momentum and charge of the secondary vertex.
//        The automatic creation of the connecting track and updating
//        of the (primary) vertex 4-momentum and charge can be suppressed
//        by specifying connect=0. In this case, however, the user
//        has to introduce the connecting track lateron by hand
//        explicitly in order to match the kinematics and charge.
//
 if (!fVertices) fVertices=new TObjArray(fNvmax);
 if (fNvtx == fNvmax) // Check if maximum vertex number is reached
 {
  fNvmax++;
  fVertices->Expand(fNvmax);
 }

 // Add the linked (secondary) vertex to the list 
 fNvtx++;
 if (fVertexCopy)
 {
  fVertices->Add(v.Clone());
 }
 else
 {
  fVertices->Add(&v);
 }

 // Create connecting track and update 4-momentum and charge for current vertex
 if (connect)
 {
  AliPosition r1=GetPosition();
  AliPosition r2=v.GetPosition();
  Float_t q=v.GetCharge();
  Ali3Vector p=v.Get3Momentum();
  Double_t v2=v.GetInvariant();
  Double_t dv2=v.Ali4Vector::GetResultError();

  AliTrack* t=new AliTrack;
  t->SetBeginPoint(r1);
  t->SetEndPoint(r2);
  t->SetCharge(q);
  t->Set3Momentum(p);
  t->SetInvariant(v2,dv2);

  AddTrack(t);

  if (!fConnects) fConnects=new TObjArray(fNvmax);
  fConnects->Add(t);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::Info(TString f)
{
// Provide vertex information within the coordinate frame f
 cout << " *AliVertex::Info* Invmass : " << GetInvmass()
      << " Charge : " << GetCharge() << " Momentum : " << GetMomentum()
      << " Ntracks : " << GetNtracks() << " Nvertices : " << fNvtx 
      << " Njets : " << fNjets << endl;
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
 if (!fVertices)
 {
  cout << " *AliVertex*::GetVertex* No (secondary) vertices present." << endl;
  return 0;
 }
 else
 {
  if (i<=0 || i>fNvtx)
  {
   cout << " *AliVertex*::GetVertex* Invalid argument i : " << i
        << " Nvtx = " << fNvtx << endl;
   return 0;
  }
  else
  {
   return (AliVertex*)fVertices->At(i-1);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::SetVertexCopy(Int_t j)
{
// (De)activate the creation of private copies of the added vertices.
// j=0 ==> No private copies are made; pointers of original vertices are stored.
// j=1 ==> Private copies of the vertices are made and these pointers are stored.
//
// Note : Once the storage contains pointer(s) to AliVertex objects one cannot
//        change the VertexCopy mode anymore.
//        To change the VertexCopy mode for an existing AliVertex containing
//        vertices one first has to invoke Reset().
 if (!fVertices)
 {
  if (j==0 || j==1)
  {
   fVertexCopy=j;
  }
  else
  {
   cout << "*AliVertex::SetVertexCopy* Invalid argument : " << j << endl;
  }
 }
 else
 {
  cout << "*AliVertex::SetVertexCopy* Storage already contained vertices."
       << "  ==> VertexCopy mode not changed." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliVertex::GetVertexCopy()
{
// Provide value of the VertexCopy mode.
// 0 ==> No private copies are made; pointers of original vertices are stored.
// 1 ==> Private copies of the vertices are made and these pointers are stored.
 return fVertexCopy;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliVertex::GetNjets()
{
// Return the current number of jets
 return fNjets;
}
///////////////////////////////////////////////////////////////////////////
AliJet* AliVertex::GetJet(Int_t i)
{
// Return the i-th jet of the current vertex
 if (!fJets)
 {
  cout << " *AliVertex*::GetJet* No jets present." << endl;
  return 0;
 }
 else
 {
  if (i<=0 || i>fNjets)
  {
   cout << " *AliVertex*::GetJet* Invalid argument i : " << i
        << " Njets = " << fNjets << endl;
   return 0;
  }
  else
  {
   return (AliJet*)fJets->At(i-1);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliVertex::SetJetCopy(Int_t j)
{
// (De)activate the creation of private copies of the added jets.
// j=0 ==> No private copies are made; pointers of original jets are stored.
// j=1 ==> Private copies of the jets are made and these pointers are stored.
//
// Note : Once the storage contains pointer(s) to AliJet objects one cannot
//        change the JetCopy mode anymore.
//        To change the JetCopy mode for an existing AliVertex containing
//        jets one first has to invoke Reset().
 if (!fJets)
 {
  if (j==0 || j==1)
  {
   fJetCopy=j;
  }
  else
  {
   cout << "*AliVertex::SetJetCopy* Invalid argument : " << j << endl;
  }
 }
 else
 {
  cout << "*AliVertex::SetJetCopy* Storage already contained jets."
       << "  ==> JetCopy mode not changed." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliVertex::GetJetCopy()
{
// Provide value of the JetCopy mode.
// 0 ==> No private copies are made; pointers of original jets are stored.
// 1 ==> Private copies of the jets are made and these pointers are stored.
 return fJetCopy;
}
///////////////////////////////////////////////////////////////////////////
