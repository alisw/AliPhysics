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
// Class AliJet
// Creation and investigation of a jet of particle tracks.
// An AliJet can be constructed by adding AliTracks.
//
// To provide maximal flexibility to the user, two modes of track storage
// are provided by means of the memberfunction SetTrackCopy().
//
// a) SetTrackCopy(0) (which is the default).
//    Only the pointers of the 'added' tracks are stored.
//    This mode is typically used by making jet studies based on a fixed list
//    of tracks which stays under user control or is contained for instance
//    in an AliEvent.  
//    In this way the AliJet just represents a 'logical structure' for the
//    physics analysis which can be embedded in e.g. an AliEvent or AliVertex.
//
//    Note :
//    Modifications made to the original tracks also affect the AliTrack objects
//    which are stored in the AliJet. 
//
// b) SetTrackCopy(1).
//    Of every 'added' track a private copy will be made of which the pointer
//    will be stored.
//    In this way the AliJet represents an entity on its own and modifications
//    made to the original tracks do not affect the AliTrack objects which are
//    stored in the AliJet. 
//    This mode will allow 'adding' many different AliTracks into an AliJet by
//    creating only one AliTrack instance in the main programme and using the
//    AliTrack::Reset() and AliTrack parameter setting memberfunctions.
//
// See also the documentation provided for the memberfunction SetOwner(). 
//
// Coding example to make 2 jets j1 and j2.
// ----------------------------------------
// j1 contains the AliTracks t1 and t2
// j2 contains 10 different AliTracks via tx
//
// AliTrack t1,t2;
//  ...
//  ... // code to fill the AliTrack data
//  ...
// AliJet j1();
// j1.AddTrack(t1);
// j1.AddTrack(t2);
//
// AliJet j2();
// j2.SetTrackCopy(1);
// AliTrack* tx=new AliTrack();
// for (Int_t i=0; i<10; i++)
// {
//  ...
//  ... // code to set momentum etc... of the track tx
//  ...
//  j2.AddTrack(tx);
//  tx->Reset();
// }
//
// j1.Data();
// j2.Data("sph");
//
// Float_t e1=j1.GetEnergy();
// Float_t pnorm=j1->GetMomentum();
// Ali3Vector p=j1->Get3Momentum();
// Float_t m=j1.GetInvmass();
// Int_t ntk=j1.GetNtracks();
// AliTrack* tj=j1.GetTrack(1);
//
// delete tx;
//
// Note : All quantities are in GeV, GeV/c or GeV/c**2
//
//--- Author: Nick van Eijndhoven 10-jul-1997 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliJet.h"
#include "Riostream.h"
 
ClassImp(AliJet) // Class implementation to enable ROOT I/O
 
AliJet::AliJet() : TNamed(),Ali4Vector()
{
// Default constructor
// All variables initialised to 0
// Initial maximum number of tracks is set to the default value
 Init();
 Reset();
 SetNtinit();
}
///////////////////////////////////////////////////////////////////////////
void AliJet::Init()
{
// Initialisation of pointers etc...
 fTracks=0;
 fNtinit=0;
 fTrackCopy=0;
}
///////////////////////////////////////////////////////////////////////////
AliJet::AliJet(Int_t n) : TNamed(),Ali4Vector()
{
// Create a jet to hold initially a maximum of n tracks
// All variables initialised to 0
 Init();
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
 if (fTracks)
 {
  delete fTracks;
  fTracks=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliJet::SetOwner(Bool_t own)
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
 if (fTracks) fTracks->SetOwner(own);
 fTrackCopy=mode;
}
///////////////////////////////////////////////////////////////////////////
AliJet::AliJet(const AliJet& j) : TNamed(j),Ali4Vector(j)
{
// Copy constructor
 fNtinit=j.fNtinit;
 fNtmax=j.fNtmax;
 fQ=j.fQ;
 fNtrk=j.fNtrk;
 fTrackCopy=j.fTrackCopy;
 fUserId=j.fUserId;

 fTracks=0;
 if (fNtrk)
 {
  fTracks=new TObjArray(fNtmax);
  if (fTrackCopy) fTracks->SetOwner();
 }

 for (Int_t i=1; i<=fNtrk; i++)
 {
  AliTrack* tx=j.GetTrack(i);
  if (fTrackCopy)
  {
   fTracks->Add(tx->Clone());
  }
  else
  {
   fTracks->Add(tx);
  }
 } 
}
///////////////////////////////////////////////////////////////////////////
void AliJet::SetNtinit(Int_t n)
{
// Set the initial maximum number of tracks for this jet
 fNtinit=n;
 fNtmax=n;
 if (fTracks)
 {
  delete fTracks;
  fTracks=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliJet::Reset()
{
// Reset all variables to 0
// The max. number of tracks is set to the initial value again
 fNtrk=0;
 fQ=0;
 fUserId=0;
 Double_t a[4]={0,0,0,0};
 SetVector(a,"sph");
 if (fNtinit > 0) SetNtinit(fNtinit);
}
///////////////////////////////////////////////////////////////////////////
void AliJet::AddTrack(AliTrack& t)
{
// Add a track to the jet.
// In case the maximum number of tracks has been reached
// space will be extended to hold an additional amount of tracks as
// was initially reserved.
// See SetTrackCopy() to tailor the functionality of the stored structures.
//
// Note :
// In case a private copy is made, this is performed via the Clone() memberfunction.
// All AliTrack and derived classes have the default TObject::Clone() memberfunction.
// However, derived classes generally contain an internal data structure which may
// include pointers to other objects. Therefore it is recommended to provide
// for all derived classes a specific copy constructor and override the default Clone()
// memberfunction using this copy constructor.
// An example for this may be seen from AliTrack.   

 AddTrack(t,1);
}
///////////////////////////////////////////////////////////////////////////
void AliJet::AddTrack(AliTrack& t,Int_t copy)
{
// Internal memberfunction to actually add a track to the jet.
// In case the maximum number of tracks has been reached
// space will be extended to hold an additional amount of tracks as
// was initially reserved.
//
// If copy=0 NO copy of the track will be made, irrespective of the setting
// of the TrackCopy flag.
// This allows a proper treatment of automatically generated connecting
// tracks between vertices.
//
// Note :
// In case a private copy is made, this is performed via the Clone() memberfunction.

 if (!fTracks)
 {
  fTracks=new TObjArray(fNtmax);
  if (fTrackCopy) fTracks->SetOwner();
 }
 if (fNtrk == fNtmax) // Check if maximum track number is reached
 {
  fNtmax+=fNtinit;
  fTracks->Expand(fNtmax);
 }
 
 // Add the track to this jet
 fNtrk++;
 if (fTrackCopy && copy)
 {
  fTracks->Add(t.Clone());
 }
 else
 {
  fTracks->Add(&t);
 }

 (*this)+=(Ali4Vector&)t;
 fQ+=t.GetCharge();

}
///////////////////////////////////////////////////////////////////////////
void AliJet::Data(TString f)
{
// Provide jet information within the coordinate frame f
 const char* name=GetName();
 const char* title=GetTitle();

 cout << " *AliJet::Data*";
 if (strlen(name))  cout << " Name : " << GetName();
 if (strlen(title)) cout << " Title : " << GetTitle();
 cout << endl;
 cout << " Id : " << fUserId << " Invmass : " << GetInvmass() << " Charge : " << fQ
      << " Momentum : " << GetMomentum() << " Ntracks : " << fNtrk << endl;

 Ali4Vector::Data(f); 
} 
///////////////////////////////////////////////////////////////////////////
void AliJet::List(TString f)
{
// Provide jet and primary track information within the coordinate frame f

 Data(f); // Information of the current jet

 // The tracks of this jet
 AliTrack* t; 
 for (Int_t it=1; it<=fNtrk; it++)
 {
  t=GetTrack(it);
  if (t)
  {
   cout << "  ---Track no. " << it << endl;
   cout << " ";
   t->Data(f); 
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

 Data(f); // Information of the current jet

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
Int_t AliJet::GetNtracks() const
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
// The error can be obtained by invoking GetResultError() after
// invokation of GetMomentum().
 Double_t norm=fV.GetNorm();
 fDresult=fV.GetResultError();
 return norm;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector AliJet::Get3Momentum() const
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
Float_t AliJet::GetCharge() const
{
// Return the total charge of the jet
 return fQ;
}
///////////////////////////////////////////////////////////////////////////
AliTrack* AliJet::GetTrack(Int_t i) const
{
// Return the i-th track of this jet
 if (!fTracks)
 {
  cout << " *AliJet*::GetTrack* No tracks present." << endl;
  return 0;
 }
 else
 {
  if (i<=0 || i>fNtrk)
  {
   cout << " *AliJet*::GetTrack* Invalid argument i : " << i
        << " Ntrk = " << fNtrk << endl;
   return 0;
  }
  else
  {
   return (AliTrack*)fTracks->At(i-1);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
AliTrack* AliJet::GetIdTrack(Int_t id) const
{
// Return the track with user identifier "id" of this jet
 AliTrack* tx=0;
 AliTrack* t=0;
 if (!fTracks)
 {
  cout << " *AliJet*::GetIdTrack* No tracks present." << endl;
  return 0;
 }
 else
 {
  for (Int_t i=0; i<fNtrk; i++)
  {
   tx=(AliTrack*)fTracks->At(i);
   if (id == tx->GetId()) t=tx;
  }
  return t;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetPt()
{
// Provide trans. momentum value w.r.t. z-axis.
// The error on the value can be obtained by GetResultError()
// after invokation of GetPt().
 Ali3Vector v;
 v=GetVecTrans();
 Double_t norm=v.GetNorm();
 fDresult=v.GetResultError();

 return norm;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetPl()
{
// Provide long. momentum value w.r.t. z-axis.
// Note : the returned value can also be negative.
// The error on the value can be obtained by GetResultError()
// after invokation of GetPl().
 Ali3Vector v;
 v=GetVecLong();

 Double_t pl=v.GetNorm();
 fDresult=v.GetResultError();

 Double_t a[3];
 v.GetVector(a,"sph");
 if (cos(a[1])<0) pl=-pl;

 return pl;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetEt()
{
// Provide trans. energy value w.r.t. z-axis.
// The error on the value can be obtained by GetResultError()
// after invokation of GetEt().
 Double_t et=GetScaTrans();

 return et;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetEl()
{
// Provide long. energy value w.r.t. z-axis.
// Note : the returned value can also be negative.
// The error on the value can be obtained by GetResultError()
// after invokation of GetEl().
 Double_t el=GetScaLong();

 return el;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetMt()
{
// Provide transverse mass value w.r.t. z-axis.
// The error on the value can be obtained by GetResultError()
// after invokation of GetMt().
 Double_t pt=GetPt();
 Double_t dpt=GetResultError();
 Double_t m=GetInvmass();
 Double_t dm=GetResultError();

 Double_t mt=sqrt(pt*pt+m*m);
 Double_t dmt2=0;
 if (mt) dmt2=(pow((pt*dpt),2)+pow((m*dm),2))/(mt*mt);

 fDresult=sqrt(dmt2);
 return mt;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetRapidity()
{
// Provide rapidity value w.r.t. z-axis.
// The error on the value can be obtained by GetResultError()
// after invokation of GetRapidity().
// Note : Also GetPseudoRapidity() is available since this class is
//        derived from Ali4Vector.
 Double_t e=GetEnergy();
 Double_t de=GetResultError();
 Double_t pl=GetPl();
 Double_t dpl=GetResultError();
 Double_t sum=e+pl;
 Double_t dif=e-pl;

 Double_t y=9999,dy2=0;
 if (sum && dif) y=0.5*log(sum/dif);

 if (sum*dif) dy2=(1./(sum*dif))*(pow((pl*de),2)+pow((e*dpl),2));

 fDresult=sqrt(dy2);
 return y;
}
///////////////////////////////////////////////////////////////////////////
void AliJet::SetTrackCopy(Int_t j)
{
// (De)activate the creation of private copies of the added tracks.
// j=0 ==> No private copies are made; pointers of original tracks are stored.
// j=1 ==> Private copies of the tracks are made and these pointers are stored.
//
// Note : Once the storage contains pointer(s) to AliTrack(s) one cannot
//        change the TrackCopy mode anymore.
//        To change the TrackCopy mode for an existing AliJet containing
//        tracks one first has to invoke Reset().
 if (!fTracks)
 {
  if (j==0 || j==1)
  {
   fTrackCopy=j;
  }
  else
  {
   cout << "*AliJet::SetTrackCopy* Invalid argument : " << j << endl;
  }
 }
 else
 {
  cout << "*AliJet::SetTrackCopy* Storage already contained tracks."
       << "  ==> TrackCopy mode not changed." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliJet::GetTrackCopy() const
{
// Provide value of the TrackCopy mode.
// 0 ==> No private copies are made; pointers of original tracks are stored.
// 1 ==> Private copies of the tracks are made and these pointers are stored.
 return fTrackCopy;
}
///////////////////////////////////////////////////////////////////////////
void AliJet::SetId(Int_t id)
{
// Set a user defined identifier for this jet.
 fUserId=id;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliJet::GetId() const
{
// Provide the user defined identifier of this jet.
 return fUserId;
}
///////////////////////////////////////////////////////////////////////////
TObject* AliJet::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.
// This feature allows e.g. AliVertex to store either AliJet objects or
// objects derived from AliJet via the AddJet memberfunction, provided
// these derived classes also have a proper Clone memberfunction. 

 AliJet* jet=new AliJet(*this);
 if (name)
 {
  if (strlen(name)) jet->SetName(name);
 }
 return jet;
}
///////////////////////////////////////////////////////////////////////////
