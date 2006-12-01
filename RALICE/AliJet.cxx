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
// Note : By default all quantities are in GeV, GeV/c or GeV/c**2
//        but the user can indicate the usage of a different scale
//        for the energy-momentum units via the SetEscale() memberfunction.
//        The actual energy-momentum unit scale can be obtained via the
//        GetEscale() memberfunction.
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
 fRef=0;
 fSelected=0;
 fEscale=1;
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
 if (fRef)
 {
  delete fRef;
  fRef=0;
 }
 if (fSelected)
 {
  delete fSelected;
  fSelected=0;
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
 fEscale=j.fEscale;
 fNtrk=j.fNtrk;
 fTrackCopy=j.fTrackCopy;
 fUserId=j.fUserId;
 if (j.fRef) fRef=new AliPositionObj(*(j.fRef));

 fSelected=0;

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
 if (fRef)
 {
  delete fRef;
  fRef=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliJet::Reset()
{
// Reset all variables to 0
// The max. number of tracks is set to the initial value again
// Note : The scale for the energy/momentum units will not be changed.
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
// Notes :
// -------
// In case a private copy is made, this is performed via the Clone() memberfunction.
// All AliTrack and derived classes have the default TObject::Clone() memberfunction.
// However, derived classes generally contain an internal data structure which may
// include pointers to other objects. Therefore it is recommended to provide
// for all derived classes a specific copy constructor and override the default Clone()
// memberfunction using this copy constructor.
// An example for this may be seen from AliTrack.   
//
// In case NO private copy is made, a check will be performed if this
// specific track is already present in the jet.
// If this is the case, no action is performed to prevent multiple
// additions of the same track.


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
// In case NO copy of the track is made, a check will be performed if this
// specific track is already present in the jet.
// If this is the case, no action is performed to prevent multiple
// additions of the same track.
//
// Note :
// In case a private copy is made, this is performed via the Clone() memberfunction.

 if (!fTracks)
 {
  fTracks=new TObjArray(fNtmax);
  if (fTrackCopy) fTracks->SetOwner();
 }
 else if (!fTrackCopy || !copy) // Check if this track is already present
 {
  for (Int_t i=0; i<fNtrk; i++)
  {
   AliTrack* tx=(AliTrack*)fTracks->At(i);
   if (tx == &t) return;
  }
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

 fQ+=t.GetCharge();

 Ali4Vector p4=(Ali4Vector)t;
 Float_t tscale=t.GetEscale();
 if ((tscale/fEscale > 1.1) || (fEscale/tscale > 1.1)) p4=p4*(tscale/fEscale);
 (*this)+=p4;

}
///////////////////////////////////////////////////////////////////////////
void AliJet::Data(TString f,TString u)
{
// Provide jet information within the coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 const char* name=GetName();
 const char* title=GetTitle();

 cout << " *AliJet::Data*";
 if (strlen(name))  cout << " Name : " << GetName();
 if (strlen(title)) cout << " Title : " << GetTitle();
 cout << endl;
 cout << " Id : " << fUserId << " Invmass : " << GetInvmass() << " Charge : " << fQ
      << " Momentum : " << GetMomentum() << " Energy scale : " << fEscale << " GeV" << endl;

 ShowTracks(0);

 Ali4Vector::Data(f,u); 
} 
///////////////////////////////////////////////////////////////////////////
void AliJet::List(TString f,TString u)
{
// Provide jet and primary track information within the coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 Data(f,u); // Information of the current jet
 if (fRef)   { cout << " Ref-point   :"; fRef->Data(f,u); }

 // The tracks of this jet
 AliTrack* t; 
 for (Int_t it=1; it<=fNtrk; it++)
 {
  t=GetTrack(it);
  if (t)
  {
   cout << "  ---Track no. " << it << endl;
   cout << " ";
   t->Data(f,u); 
  }
  else
  {
   cout << " *AliJet::List* Error : No track present." << endl; 
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
void AliJet::ListAll(TString f,TString u)
{
// Provide jet and prim.+sec. track information within the coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 Data(f,u); // Information of the current jet
 if (fRef)   { cout << " Ref-point   :"; fRef->Data(f,u); }

 // The tracks of this jet
 AliTrack* t; 
 for (Int_t it=1; it<=fNtrk; it++)
 {
  t=GetTrack(it);
  if (t)
  {
   cout << "  ---Track no. " << it << endl;
   cout << " ";
   t->ListAll(f,u); 
  }
  else
  {
   cout << " *AliJet::List* Error : No track present." << endl; 
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
Int_t AliJet::GetNtracks(Int_t idmode,Int_t chmode,Int_t pcode)
{
// Provide the number of user selected tracks in this jet based on the
// idmode, chmode and pcode selections as specified by the user.
// For specification of the selection parameters see GetTracks().
// The default parameters correspond to no selection, which implies
// that invokation of GetNtracks() just returns the total number of
// tracks registered in this jet.
//
// Note : In case certain selections are specified, this function
//        invokes GetTracks(idmode,chmode,pcode) to determine the
//        number of tracks corresponding to the selections.
//        When the jet contains a large number of tracks, invokation
//        of GetTracks(idmode,chmode,pcode) and subsequently invoking
//        GetEntries() for the resulting TObjArray* might be slightly
//        faster.

 Int_t n=0;
 if (idmode==0 && chmode==2 && pcode==0)
 {
  return fNtrk;
 }
 else
 {
  TObjArray* arr=GetTracks(idmode,chmode,pcode);
  n=arr->GetEntries();
  return n;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetEnergy(Float_t scale)
{
// Return the total energy of the jet.
// By default the energy is returned in the units as it was stored in the jet
// structure. However, the user can select a different energy unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV, so specification
// of scale=0.001 will provide the energy in MeV.
// The error can be obtained by invoking GetResultError() after
// invokation of GetEnergy().
 Double_t E=GetScalar();
 if (E>0)
 {
  if (scale>0)
  {
   E*=fEscale/scale;
   fDresult*=fEscale/scale;
  }
  return E;
 }
 else
 {
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetMomentum(Float_t scale)
{
// Return the value of the total jet 3-momentum
// By default the momentum is returned in the units as it was stored in the jet
// structure. However, the user can select a different momentum unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c, so specification
// of scale=0.001 will provide the momentum in MeV/c.
// The error can be obtained by invoking GetResultError() after
// invokation of GetMomentum().

 Double_t norm=fV.GetNorm();
 fDresult=fV.GetResultError();
 if (scale>0)
 {
  norm*=fEscale/scale;
  fDresult*=fEscale/scale;
 }
 return norm;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector AliJet::Get3Momentum(Float_t scale) const
{
// Return the the total jet 3-momentum
// By default the components of the 3-momentum are returned in the units
// as they were stored in the jet structure.
// However, the user can select a different momentum unit scale for the
// components by specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c, so specification
// of scale=0.001 will provide the 3-momentum in MeV/c.

 Ali3Vector p=Get3Vector();
 if (scale>0) p*=fEscale/scale;
 return p;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetInvmass(Float_t scale)
{
// Return the invariant mass of the jet.
// By default the mass is returned in the units as it was stored in the jet
// structure. However, the user can select a different mass unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c**2, so specification
// of scale=0.001 will provide the invariant mass in MeV/c**2.
// The error can be obtained by invoking GetResultError() after
// invokation of GetInvmass().

 Double_t inv=Dot(*this);
 Double_t dinv=GetResultError();
 Double_t dm=0;
 if (inv >= 0)
 {
 Double_t m=sqrt(inv);
 if (m) dm=dinv/(2.*m);
 if (scale>0)
 {
  m*=fEscale/scale;
  dm*=fEscale/scale;
 }
 fDresult=dm;
 return m;
 }
 else
 {
  fDresult=dm;
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

 if (!fTracks) return 0;

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
///////////////////////////////////////////////////////////////////////////
AliTrack* AliJet::GetIdTrack(Int_t id) const
{
// Return the track with user identifier "id" of this jet
 if (!fTracks) return 0;

 AliTrack* tx=0;
 for (Int_t i=0; i<fNtrk; i++)
 {
  tx=(AliTrack*)fTracks->At(i);
  if (id == tx->GetId()) return tx;
 }
 return 0; // No matching id found
}
///////////////////////////////////////////////////////////////////////////
TObjArray* AliJet::GetTracks(Int_t idmode,Int_t chmode,Int_t pcode)
{
// Provide references to user selected tracks based on the idmode, chmode
// and pcode selections as specified by the user.
//
// The following selection combinations are available :
// ----------------------------------------------------
// idmode = -1 ==> Select tracks with negative user identifier "id"
//           0 ==> No selection on user identifier
//           1 ==> Select tracks with positive user identifier "id"
//
// chmode = -1 ==> Select tracks with negative charge
//           0 ==> Select neutral tracks
//           1 ==> Select tracks with positive charge
//           2 ==> No selection on charge
//           3 ==> Select all charged tracks
//
// pcode  =  0 ==> No selection on particle code
//           X ==> Select tracks with particle code +X or -X
//                 This allows selection of both particles and anti-particles
//                 in case of PDG particle codes.
//                 Selection of either particles or anti-particles can be
//                 obtained in combination with the "chmode" selector.
//
// Examples :
// ----------
// idmode=-1 chmode=0 pcode=0   : Selection of all neutral tracks with negative id.
// idmode=0  chmode=2 pcode=211 : Selection of all charged pions (PDG convention).
// idmode=0  chmode=1 pcode=321 : Selection of all positive kaons (PDG convention).
//
// The default values are idmode=0 chmode=2 pcode=0 (i.e. no selections applied).
//
// Notes :
// -------
// 1) In case the user has labeled simulated tracks with negative id and
//    reconstructed tracks with positive id, this memberfunction provides
//    easy access to either all simulated or reconstructed tracks.
// 2) Subsequent invokations of this memberfunction with e.g. chmode=-1 and chmode=1
//    provides a convenient way to investigate particle pairs with opposite charge
//    (e.g. for invariant mass analysis).
// 3) The selected track pointers are returned via a multi-purpose array,
//    which will be overwritten by subsequent selections.
//    In case the selected track list is to be used amongst other selections,
//    the user is advised to store the selected track pointers in a local
//    TObjArray or TRefArray.  

 if (fSelected)
 {
  fSelected->Clear();
 }
 else
 {
  fSelected=new TObjArray();
 }

 if (!fTracks) return fSelected;

 AliTrack* tx=0;
 Int_t code=0;
 Int_t id=0;
 Float_t q=0;
 for (Int_t i=0; i<fNtrk; i++)
 {
  tx=(AliTrack*)fTracks->At(i);
  if (!tx) continue;

  code=tx->GetParticleCode();
  if (pcode && abs(pcode)!=abs(code)) continue;

  id=tx->GetId();
  if (idmode==-1 && id>=0) continue;
  if (idmode==1 && id<=0) continue;

  q=tx->GetCharge();
  if (chmode==-1 && q>=0) continue;
  if (chmode==0 && fabs(q)>1e-10) continue;
  if (chmode==1 && q<=0) continue;
  if (chmode==3 && fabs(q)<1e-10) continue;

  fSelected->Add(tx);
 }

 return fSelected;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* AliJet::GetTracks(TString name)
{
// Provide references to all tracks with the specified name.
//
// Notes :
// -------
// 1) In case the user has labeled reconstructed tracks with the name of
//    the applied reconstruction algorithm, this memberfunction provides
//    easy access to all tracks reconstructed by a certain method.
// 2) The selected track pointers are returned via a multi-purpose array,
//    which will be overwritten by subsequent selections.
//    In case the selected track list is to be used amongst other selections,
//    the user is advised to store the selected track pointers in a local
//    TObjArray or TRefArray.  

 if (fSelected)
 {
  fSelected->Clear();
 }
 else
 {
  fSelected=new TObjArray();
 }

 if (!fTracks) return fSelected;

 AliTrack* tx=0;
 TString s;
 for (Int_t i=0; i<fNtrk; i++)
 {
  tx=(AliTrack*)fTracks->At(i);
  if (!tx) continue;

  s=tx->GetName();
  if (s == name) fSelected->Add(tx);
 }

 return fSelected;
}
///////////////////////////////////////////////////////////////////////////
void AliJet::RemoveTracks(TString name)
{
// Remove all tracks with the specified name.
// If name="*" all tracks will be removed.
//
// Note :
// ------
// In case the user has labeled reconstructed tracks with the name of
// the applied reconstruction algorithm, this memberfunction provides
// easy removal of all tracks reconstructed by a certain method.

 if (!fTracks) return;

 AliTrack* tx=0;
 TString s;
 TObject* obj=0;
 for (Int_t i=0; i<fNtrk; i++)
 {
  tx=(AliTrack*)fTracks->At(i);
  if (!tx) continue;

  s=tx->GetName();
  if (s==name || name=="*")
  {
   obj=fTracks->Remove(tx);
   if (obj && fTracks->IsOwner()) delete tx;
  }
 }
 fTracks->Compress();
 fNtrk=fTracks->GetEntries();
}
///////////////////////////////////////////////////////////////////////////
void AliJet::RemoveTracks(Int_t idmode,Int_t chmode,Int_t pcode)
{
// Remove user selected tracks based on the idmode, chmode and pcode
// selections as specified by the user.
// For defintions of these selections see the corresponding GetTracks()
// memberfunction.

 if (!fTracks) return;

 TObjArray* arr=GetTracks(idmode,chmode,pcode);
 if (!arr) return;
 
 Int_t ntk=arr->GetEntries();
 if (!ntk) return;

 AliTrack* tx=0;
 TObject* obj=0;
 for (Int_t i=0; i<ntk; i++)
 {
  tx=(AliTrack*)arr->At(i);
  if (!tx) continue;

  obj=fTracks->Remove(tx);
  if (obj && fTracks->IsOwner()) delete tx;
 }
 fTracks->Compress();
 fNtrk=fTracks->GetEntries();
 arr->Clear();
}
///////////////////////////////////////////////////////////////////////////
void AliJet::ShowTracks(Int_t mode)
{
// Provide an overview of the available tracks.
// The argument mode determines the amount of information as follows :
// mode = 0 ==> Only printout of the number of tracks
//        1 ==> Provide a listing with 1 line of info for each track
//
// The default is mode=1.
//
 Int_t ntk=GetNtracks();
 if (ntk)
 {
  if (!mode)
  {
   cout << " There are " << ntk << " tracks available." << endl; 
  }
  else
  {
   cout << " The following " << ntk << " tracks are available :" << endl; 
   for (Int_t i=1; i<=ntk; i++)
   {
    AliTrack* tx=GetTrack(i);
    if (tx)
    {
     const char* name=tx->GetName();
     const char* title=tx->GetTitle();
     cout << " Track : " << i;
     cout << " Id : " << tx->GetId();
     cout << " Q : " << tx->GetCharge() << " m : " << tx->GetMass() << " p : " << tx->GetMomentum();
     if (strlen(name)) cout << " Name : " << name;
     if (strlen(title)) cout << " Title : " << title;
     cout << endl;
    }
   }
  }
 }
 else
 {
  cout << " No tracks are present." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetPt(Float_t scale)
{
// Provide the transverse momentum value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the jet
// structure. However, the user can select a different momentum unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c, so specification
// of scale=0.001 will provide the transverse momentum in MeV/c.
// The error on the value can be obtained by GetResultError()
// after invokation of GetPt().
 Ali3Vector v;
 v=GetVecTrans();
 Double_t norm=v.GetNorm();
 fDresult=v.GetResultError();
 if (scale>0)
 {
  norm*=fEscale/scale;
  fDresult*=fEscale/scale;
 }

 return norm;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetPl(Float_t scale)
{
// Provide the longitudinal momentum value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the jet
// structure. However, the user can select a different momentum unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c, so specification
// of scale=0.001 will provide the longitudinal momentum in MeV/c.
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
 if (scale>0)
 {
  pl*=fEscale/scale;
  fDresult*=fEscale/scale;
 }

 return pl;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetEt(Float_t scale)
{
// Provide transverse energy value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the jet
// structure. However, the user can select a different energy unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV, so specification
// of scale=0.001 will provide the transverse energy in MeV.
// The error on the value can be obtained by GetResultError()
// after invokation of GetEt().

 Double_t et=GetScaTrans();
 if (scale>0)
 {
  et*=fEscale/scale;
  fDresult*=fEscale/scale;
 }

 return et;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetEl(Float_t scale)
{
// Provide longitudinal energy value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the jet
// structure. However, the user can select a different energy unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV, so specification
// of scale=0.001 will provide the longitudinal energy in MeV.
// Note : the returned value can also be negative.
// The error on the value can be obtained by GetResultError()
// after invokation of GetEl().

 Double_t el=GetScaLong();
 if (scale>0)
 {
  el*=fEscale/scale;
  fDresult*=fEscale/scale;
 }

 return el;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetMt(Float_t scale)
{
// Provide transverse mass value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the jet
// structure. However, the user can select a different energy unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV, so specification
// of scale=0.001 will provide the transverse mass in MeV.
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
 if (scale>0)
 {
  mt*=fEscale/scale;
  fDresult*=fEscale/scale;
 }
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
void AliJet::SetReferencePoint(AliPosition& p)
{
// Store the position of the jet reference-point.
// The reference-point of a jet provides a means to define a generic
// space-time location for the jet as a whole.
// This doesn't have to be necessarily the location where all the constituent
// tracks originate (e.g. a bundle of parallel tracks doesn't have such
// a location). As such the meaning of this reference-point is different from
// a normal vertex position and allows to provide complimentary information. 
// This reference point is the preferable point to start e.g. extrapolations
// and investigate coincidences in space and/or time.
 if (fRef) delete fRef;
 fRef=new AliPositionObj(p);
}
///////////////////////////////////////////////////////////////////////////
AliPosition* AliJet::GetReferencePoint()
{
// Provide the position of the jet reference-point.
// The reference-point of a jet provides a means to define a generic
// space-time location for the jet as a whole.
// This doesn't have to be necessarily the location where all the constituent
// tracks originate (e.g. a bundle of parallel tracks doesn't have such
// a location). As such the meaning of this reference-point is different from
// a normal vertex position and allows to provide complimentary information. 
// This reference point is the preferable point to start e.g. extrapolations
// and investigate coincidences in space and/or time.
 return fRef;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* AliJet::SortTracks(Int_t mode,TObjArray* tracks)
{
// Order the references to an array of tracks by looping over the input array "tracks"
// and checking the value of a certain observable.
// The ordered array is returned as a TObjArray.
// In case tracks=0 (default), the registered tracks of the current jet are used. 
// Note that the original track array is not modified.
// Via the "mode" argument the user can specify the observable to be checked upon
// and specify whether sorting should be performed in decreasing order (mode<0)
// or in increasing order (mode>0).
//
// The convention for the observable selection is the following :
// mode : 1 ==> Number of signals associated to the track
//        2 ==> Track energy
//        3 ==> Track momentum
//        4 ==> Mass of the track
//        5 ==> Transverse momentum of the track
//        6 ==> Longitudinal momentum of the track
//        7 ==> Transverse energy of the track
//        8 ==> Longitudinal energy of the track
//        9 ==> Transverse mass of the track
//       10 ==> Track rapidity
//       11 ==> Pseudo-rapidity of the track
//       12 ==> Charge of the track
//       13 ==> Probability of the track hypothesis
//
// The default is mode=-1.
//
// Note : This sorting routine uses a common area in memory, which is used
//        by various other sorting facilities as well.
//        This means that the resulting sorted TObjArray may be overwritten
//        when another sorting is invoked.
//        To retain the sorted list of pointers, the user is advised to copy
//        the pointers contained in the returned TObjArray into a private
//        TObjArray instance.

 if (fSelected)
 {
  delete fSelected;
  fSelected=0;
 }

 if (!tracks) tracks=fTracks;
 
 if (!mode || abs(mode)>13 || !tracks) return fSelected;

 Int_t ntracks=tracks->GetEntries();
 if (!ntracks)
 {
  return fSelected;
 }
 else
 {
  fSelected=new TObjArray(ntracks);
 }

 Double_t val1,val2; // Values of the observable to be tested upon
 
 Int_t nord=0;
 for (Int_t i=0; i<ntracks; i++) // Loop over all tracks of the array
 {
  AliTrack* tx=(AliTrack*)tracks->At(i);

  if (!tx) continue;
 
  if (nord == 0) // store the first track at the first ordered position
  {
   nord++;
   fSelected->AddAt(tx,nord-1);
   continue;
  }
 
  for (Int_t j=0; j<=nord; j++) // put track in the right ordered position
  {
   if (j == nord) // track has smallest (mode<0) or largest (mode>0) observable value seen so far
   {
    nord++;
    fSelected->AddAt(tx,j); // add track at the end
    break; // go for next track
   }

   val1=0;
   val2=0;

   switch (abs(mode))
   {
    case 1:
     val1=tx->GetNsignals();
     val2=((AliTrack*)fSelected->At(j))->GetNsignals();
     break;
    case 2:
     val1=tx->GetEnergy(1);
     val2=((AliTrack*)fSelected->At(j))->GetEnergy(1);
     break;
    case 3:
     val1=tx->GetMomentum(1);
     val2=((AliTrack*)fSelected->At(j))->GetMomentum(1);
     break;
    case 4:
     val1=tx->GetMass(1);
     val2=((AliTrack*)fSelected->At(j))->GetMass(1);
     break;
    case 5:
     val1=tx->GetPt(1);
     val2=((AliTrack*)fSelected->At(j))->GetPt(1);
     break;
    case 6:
     val1=tx->GetPl(1);
     val2=((AliTrack*)fSelected->At(j))->GetPl(1);
     break;
    case 7:
     val1=tx->GetEt(1);
     val2=((AliTrack*)fSelected->At(j))->GetEt(1);
     break;
    case 8:
     val1=tx->GetEl(1);
     val2=((AliTrack*)fSelected->At(j))->GetEl(1);
     break;
    case 9:
     val1=tx->GetMt(1);
     val2=((AliTrack*)fSelected->At(j))->GetMt(1);
     break;
    case 10:
     val1=tx->GetRapidity();
     val2=((AliTrack*)fSelected->At(j))->GetRapidity();
     break;
    case 11:
     val1=tx->GetPseudoRapidity();
     val2=((AliTrack*)fSelected->At(j))->GetPseudoRapidity();
     break;
    case 12:
     val1=tx->GetCharge();
     val2=((AliTrack*)fSelected->At(j))->GetCharge();
     break;
    case 13:
     val1=tx->GetProb();
     val2=((AliTrack*)fSelected->At(j))->GetProb();
     break;
   }

   if (mode<0 && val1 <= val2) continue;
   if (mode>0 && val1 >= val2) continue;
 
   nord++;
   for (Int_t k=nord-1; k>j; k--) // create empty position
   {
    fSelected->AddAt(fSelected->At(k-1),k);
   }
   fSelected->AddAt(tx,j); // put track at empty position
   break; // go for next track
  }
 }
 return fSelected;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetDistance(AliPosition* p,Float_t scale)
{
// Provide distance of the current jet to the position p.
// The error on the result can be obtained as usual by invoking
// GetResultError() afterwards. 
//
// By default the distance will be provided in the metric unit scale of
// the AliPosition p.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the distance in cm.
// As such it is possible to obtain a correctly computed distance even in case
// the jet parameters have a different unit scale.
// However, it is recommended to work always with one single unit scale.
//
// Note : In case of incomplete information, a distance value of -1 is
//        returned.
 
 Double_t dist=-1.;
 fDresult=0.;

 if (!p) return dist;

 // Obtain a defined position on this jet
 AliPosition* rx=fRef;

 if (!rx) return dist;

 Ali3Vector pj=Get3Momentum();

 if (pj.GetNorm() <= 0.) return dist;

 AliTrack tj;
 tj.Set3Momentum(pj);
 tj.SetReferencePoint(*rx);
 dist=tj.GetDistance(p,scale);
 fDresult=tj.GetResultError();
 return dist;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetDistance(AliTrack* t,Float_t scale)
{
// Provide distance of the current jet to the track t.
// The error on the result can be obtained as usual by invoking
// GetResultError() afterwards. 
//
// By default the distance will be provided in the metric unit scale of
// the current jet.
// However, the user can specify a required metric unit scale by specification
// of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the distance in cm.
// As such it is possible to obtain a correctly computed distance even in case
// the jet and track parameters have a different unit scale.
// However, it is recommended to work always with one single unit scale.
//
// Note : In case of incomplete information, a distance value of -1 is
//        returned.
 
 Double_t dist=-1.;
 fDresult=0.;

 if (!t) return dist;

 // Obtain a defined position on this jet
 AliPosition* rx=fRef;

 if (!rx) return dist;

 Ali3Vector pj=Get3Momentum();

 if (pj.GetNorm() <= 0.) return dist;

 AliTrack tj;
 tj.Set3Momentum(pj);
 tj.SetReferencePoint(*rx);
 dist=tj.GetDistance(t,scale);
 fDresult=tj.GetResultError();
 return dist;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliJet::GetDistance(AliJet* j,Float_t scale)
{
// Provide distance of the current jet to the jet j.
// The error on the result can be obtained as usual by invoking
// GetResultError() afterwards. 
//
// By default the distance will be provided in the metric unit scale of
// the current jet.
// This implies that the results of j1.GetDistance(j2) and j2.GetDistance(j1)
// may be numerically different in case j1 and j2 have different metric units.
// However, the user can specify a required metric unit scale by specification
// of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the distance in cm.
// As such it is possible to obtain a correctly computed distance even in case
// the jet parameters have a different unit scale.
// However, it is recommended to work always with one single unit scale.
//
// Note : In case of incomplete information, a distance value of -1 is
//        returned.
 
 Double_t dist=-1.;
 fDresult=0.;

 if (!j) return dist;

 // Obtain a defined position on jet j
 AliPosition* rx=j->GetReferencePoint();

 if (!rx) return dist;

 Ali3Vector pj=j->Get3Momentum();

 if (pj.GetNorm() <= 0.) return dist;

 AliTrack tj;
 tj.Set3Momentum(pj);
 tj.SetReferencePoint(*rx);
 dist=GetDistance(tj,scale);
 return dist;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliJet::GetNsignals() const
{
// Provide the number of signals associated to the jet tracks.
// Note : Multiple occurrences of the same signal are only counted once. 

 if (fNtrk<1) return 0;

 TObjArray arr;

 Int_t n=0;
 AliTrack* tx=0;
 Int_t exists=0;
 for (Int_t i=1; i<=fNtrk; i++)
 {
  tx=GetTrack(i);
  for (Int_t j=1; j<=tx->GetNsignals(); j++)
  {
   AliSignal* sx=tx->GetSignal(j);
   if (!sx) continue;
   exists=0;
   for (Int_t k=0; k<arr.GetEntries(); k++)
   {
    if (sx==(AliSignal*)arr.At(k))
    {
     exists=1;
     break;
    } 
   }
   if (!exists) arr.Add(sx);
  } 
 }
 n=arr.GetEntries();
 return n;
}
///////////////////////////////////////////////////////////////////////////
void AliJet::SetEscale(Float_t scale)
{
// Indicate the energy/momentum scale as used by the user.
// The convention is that scale=1 indicates values in units
// of GeV, GeV/c or GeV/c**2.
// So, in case one decides to store values in units of MeV, MeV/c or MeV/c**2
// the scale indicator should be set to scale=0.001.
//
// By default scale=1 is set in the constructor.

 if (scale>0)
 {
  fEscale=scale;
 }
 else
 {
  cout << " *AliJet::SetEscale* Invalid scale value : " << scale << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliJet::GetEscale() const
{
// Provide the energy/momentum scale as used by the user.
// The convention is that scale=1 indicates values in units
// of GeV, GeV/c or GeV/c**2.
// So, a value of scale=0.001 indicates that energy/momentum values are
// stored in units of MeV, MeV/c or MeV/c**2.
 return fEscale;
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
