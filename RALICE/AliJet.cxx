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
Revision 1.2  1999/09/29 09:24:28  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////
// Class AliJet
// Creation and investigation of a jet of particle tracks.
// An AliJet can be constructed by adding AliTracks.
//
// Coding example to make 2 jets j1 and j2.
// ----------------------------------------
// j1 contains the AliTracks 1 and 2
// j2 contains the AliTracks 3 and 4
//
// AliTrack t1,t2,t3,t4;
//  ...
//  ... // code to fill the AliTrack data
//  ...
// AliJet j1(5);
// AliJet j2(12);
// j1.Add(t1);
// j1.Add(t2);
// j2.Add(t3);
// j2.Add(t4);
//
// j1.Info();
// j2.Info("sph");
//
// Float_t e1=j1.GetEnergy();
// Float_t pnorm=j1->GetMomentum();
// Ali3Vector p=j1->Get3Momentum();
// Float_t m=j1.GetInvmass();
// Int_t ntk=j1.GetNtracks();
// AliTrack* tj=j1.GetTrack(1);
//
// Note : All quantities are in GeV, GeV/c or GeV/c**2
//
//--- Author: Nick van Eijndhoven 10-jul-1997 UU-SAP Utrecht
//- Modified: NvE 06-apr-1999 UU-SAP Utrecht to inherit from Ali4Vector
///////////////////////////////////////////////////////////////////////////

#include "AliJet.h"
 
ClassImp(AliJet) // Class implementation to enable ROOT I/O
 
AliJet::AliJet()
{
// Default constructor
// All variables initialised to 0
// Initial maximum number of tracks is set to the default value
 fTracks=0;
 fNtinit=0;
 Reset();
 SetNtinit();
}
///////////////////////////////////////////////////////////////////////////
AliJet::AliJet(Int_t n)
{
// Create a jet to hold initially a maximum of n tracks
// All variables initialised to 0
 fTracks=0;
 fNtinit=0;
 Reset();
 if (n > 0)
 {
  SetNtinit(n);
 }
 else
 {
  cout << endl;
  cout << " *AliJet* Initial max. number of tracks entered : " << n << endl;
  cout << " This is invalid. Default initial maximum will be used." << endl;
  cout << endl;
  SetNtinit();
 }
}
///////////////////////////////////////////////////////////////////////////
AliJet::~AliJet()
{
// Default destructor
 if (fTracks) delete fTracks;
 fTracks=0;
}
///////////////////////////////////////////////////////////////////////////
void AliJet::SetNtinit(Int_t n)
{
// Set the initial maximum number of tracks for this jet
 fNtinit=n;
 fNtmax=n;
 if (fTracks) delete fTracks;
 fTracks=new TObjArray(fNtmax);
}
///////////////////////////////////////////////////////////////////////////
void AliJet::Reset()
{
// Reset all variables to 0
// The max. number of tracks is set to the initial value again
 fNtrk=0;
 fQ=0;
 Double_t a[4]={0,0,0,0};
 SetVector(a,"sph");
 if (fNtinit > 0) SetNtinit(fNtinit);
}
///////////////////////////////////////////////////////////////////////////
void AliJet::Add(AliTrack& t)
{
// Add a track to the jet
// In case the maximum number of tracks has been reached
// space will be extended to hold an additional amount of tracks as
// was initially reserved
 if (fNtrk == fNtmax) // Check if maximum track number is reached
 {
  fNtmax+=fNtinit;
  fTracks->Expand(fNtmax);
 }
 
 // Add the track to this jet
 fNtrk++;
 fTracks->Add(&t);
 (*this)+=(Ali4Vector&)t;
 fQ+=t.GetCharge();
}
///////////////////////////////////////////////////////////////////////////
void AliJet::Info(TString f)
{
// Provide jet information within the coordinate frame f
 cout << " *AliJet::Info* Invmass : " << GetInvmass() << " Charge : " << fQ
      << " Momentum : " << GetMomentum() << " Ntracks : " << fNtrk << endl;
 cout << " ";
 Ali4Vector::Info(f); 
} 
///////////////////////////////////////////////////////////////////////////
void AliJet::List(TString f)
{
// Provide jet and primary track information within the coordinate frame f

 Info(f); // Information of the current jet

 // The tracks of this jet
 AliTrack* t; 
 for (Int_t it=1; it<=fNtrk; it++)
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
   cout << " *AliJet::List* Error : No track present." << endl; 
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
void AliJet::ListAll(TString f)
{
// Provide jet and prim.+sec. track information within the coordinate frame f

 Info(f); // Information of the current jet

 // The tracks of this jet
 AliTrack* t; 
 for (Int_t it=1; it<=fNtrk; it++)
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
   cout << " *AliJet::List* Error : No track present." << endl; 
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
Int_t AliJet::GetNtracks()
{
// Return the current number of tracks of this jet
 return fNtrk;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetEnergy()
{
// Return the total energy of the jet
 return GetScalar();
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetMomentum()
{
// Return the value of the total jet 3-momentum
 Ali3Vector p=Get3Vector();
 Double_t p2=p.Dot(p);
 return sqrt(p2);
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector AliJet::Get3Momentum()
{
// Return the the total jet 3-momentum
 Ali3Vector p=Get3Vector();
 return p;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetInvmass()
{
// Return the invariant mass of the jet
 Double_t m2=Dot(*this);
 if (m2>0)
 {
  return sqrt(m2);
 }
 else
 {
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliJet::GetCharge()
{
// Return the total charge of the jet
 return fQ;
}
///////////////////////////////////////////////////////////////////////////
AliTrack* AliJet::GetTrack(Int_t i)
{
// Return the i-th track of this jet
 return (AliTrack*)fTracks->At(i-1);
}
///////////////////////////////////////////////////////////////////////////
