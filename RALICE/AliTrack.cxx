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
// Class AliTrack
// Handling of the attributes of a reconstructed particle track.
//
// Coding example :
// ----------------
//
// Float_t a[4]={195.,1.2,-0.04,8.5};
// Ali4Vector pmu;
// pmu.SetVector(a,"car");
// AliTrack t1;
// t1.Set4Momentum(pmu);
//
// Float_t b[3]={1.2,-0.04,8.5};
// Ali3Vector p;
// p.SetVector(b,"car");
// AliTrack t2;
// t2.Set3Momentum(p);
// t2.SetCharge(0);
// t2.SetMass(1.115);
//
// t1.Data();
// t2.Data();
//
// Float_t pi=acos(-1.);
// Float_t thcms=0.2*pi; // decay theta angle in cms
// Float_t phicms=pi/4.; // decay theta angle in cms
// Float_t m1=0.938;
// Float_t m2=0.140;
// t2.Decay(m1,m2,thcms,phicms); // Track t2 decay : Lambda -> proton + pion
//
// t2.List();
//
// Int_t ndec=t2.GetNdecay();
// AliTrack* d1=t2.GetDecayTrack(1); // Access to decay track number 1
// AliTrack* d2=t2.GetDecayTrack(2); // Access to decay track number 2
//
// AliSignal s1,s2,s3,s4;
//
// .... // Code (e.g. detector readout) to fill AliSignal data
//
// AliTrack trec; // Track which will be reconstructed from signals
// trec.AddSignal(s1);
// trec.AddSignal(s3);
// trec.AddSignal(s4);
//
// Ali3Vector P;
// Float_t Q,M;
//
// ... // Code which accesses signals from trec and reconstructs
//        3-momentum P, charge Q, mass M etc...
//
// trec.Set3Momentum(P);
// trec.SetCharge(Q);
// trec.SetMass(M);
//
// Float_t r1[3]={1.6,-3.8,25.7};
// Float_t er1[3]={0.2,0.5,1.8};
// Float_t r2[3]={8.6,23.8,-6.7};
// Float_t er2[3]={0.93,1.78,0.8};
// AliPosition begin,end;
// begin.SetPosition(r1,"car");
// begin.SetPositionErrors(er1,"car");
// end.SetPosition(r2,"car");
// end.SetPositionErrors(er2,"car");
// trec.SetBeginPoint(begin);
// trec.SetEndPoint(end);
// 
// Note : All quantities are in GeV, GeV/c or GeV/c**2
//
//--- Author: Nick van Eijndhoven 10-jul-1997 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliTrack.h"
#include "Riostream.h"
 
ClassImp(AliTrack) // Class implementation to enable ROOT I/O
 
AliTrack::AliTrack() : TNamed(),Ali4Vector()
{
// Default constructor
// All variables initialised to 0
 Init();
 Reset();
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Init()
{
// Initialisation of pointers etc...
 fDecays=0;
 fSignals=0;
 fHypotheses=0;
 fBegin=0;
 fEnd=0;
 fImpactXY=0;
 fImpactXZ=0;
 fImpactYZ=0;
 fClosest=0;
 fParent=0;
}
///////////////////////////////////////////////////////////////////////////
AliTrack::~AliTrack()
{
// Destructor to delete memory allocated for decay tracks array.
// This destructor automatically cleares the pointer of this AliTrack
// from all the link slots of the related AliSignal objects.

 Int_t nsig=GetNsignals();
 for (Int_t i=1; i<=nsig; i++)
 {
  AliSignal* s=GetSignal(i);
  if (s) s->ResetLinks(this);
 }
 
 if (fDecays)
 {
  delete fDecays;
  fDecays=0;
 }
 if (fSignals)
 {
  fSignals->Clear();
  delete fSignals;
  fSignals=0;
 }
 if (fHypotheses)
 {
  delete fHypotheses;
  fHypotheses=0;
 }
 if (fBegin)
 {
  delete fBegin;
  fBegin=0;
 }
 if (fEnd)
 {
  delete fEnd;
  fEnd=0;
 }
 if (fImpactXY)
 {
  delete fImpactXY;
  fImpactXY=0;
 }
 if (fImpactXZ)
 {
  delete fImpactXZ;
  fImpactXZ=0;
 }
 if (fImpactYZ)
 {
  delete fImpactYZ;
  fImpactYZ=0;
 }
 if (fClosest)
 {
  delete fClosest;
  fClosest=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliTrack::AliTrack(const AliTrack& t) : TNamed(t),Ali4Vector(t)
{
// Copy constructor
 Init();

 fQ=t.fQ;
 fProb=t.fProb;
 if (t.fBegin) fBegin=new AliPositionObj(*(t.fBegin));
 if (t.fEnd) fEnd=new AliPositionObj(*(t.fEnd));
 if (t.fImpactXY) fImpactXY=new AliPositionObj(*(t.fImpactXY));
 if (t.fImpactXZ) fImpactXZ=new AliPositionObj(*(t.fImpactXZ));
 if (t.fImpactYZ) fImpactYZ=new AliPositionObj(*(t.fImpactYZ));
 if (t.fClosest) fClosest=new AliPositionObj(*(t.fClosest));
 fUserId=t.fUserId;
 fChi2=t.fChi2;
 fNdf=t.fNdf;
 fCode=t.fCode;
 fParent=t.fParent;

 Int_t ndec=t.GetNdecay();
 if (ndec)
 {
  fDecays=new TObjArray(ndec);
  fDecays->SetOwner();
  for (Int_t it=1; it<=ndec; it++)
  {
   AliTrack* tx=t.GetDecayTrack(it);
   fDecays->Add(new AliTrack(*tx));
  }
 }

 Int_t nsig=t.GetNsignals();
 if (nsig)
 {
  fSignals=new TObjArray(nsig);
  for (Int_t is=1; is<=nsig; is++)
  {
   AliSignal* sx=t.GetSignal(is);
   fSignals->Add(sx);
  }
 }

 Int_t nhyp=t.GetNhypotheses();
 if (nhyp)
 {
  fHypotheses=new TObjArray(nhyp);
  fHypotheses->SetOwner();
  for (Int_t ih=1; ih<=nhyp; ih++)
  {
   AliTrack* tx=t.GetTrackHypothesis(ih);
   fHypotheses->Add(new AliTrack(*tx));
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Reset()
{
// Reset all variables to 0 and delete all auto-generated decay tracks.
 fQ=0;
 fChi2=0;
 fNdf=0;
 fUserId=0;
 fCode=0;
 fProb=0;
 Double_t a[4]={0,0,0,0};
 SetVector(a,"sph");
 fParent=0;
 if (fDecays)
 {
  delete fDecays;
  fDecays=0;
 }
 if (fSignals)
 {
  fSignals->Clear();
  delete fSignals;
  fSignals=0;
 }
 if (fHypotheses)
 {
  delete fHypotheses;
  fHypotheses=0;
 }
 if (fBegin)
 {
  delete fBegin;
  fBegin=0;
 }
 if (fEnd)
 {
  delete fEnd;
  fEnd=0;
 }
 if (fImpactXY)
 {
  delete fImpactXY;
  fImpactXY=0;
 }
 if (fImpactXZ)
 {
  delete fImpactXZ;
  fImpactXZ=0;
 }
 if (fImpactYZ)
 {
  delete fImpactYZ;
  fImpactYZ=0;
 }
 if (fClosest)
 {
  delete fClosest;
  fClosest=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Set3Momentum(Ali3Vector& p)
{
// Set the track parameters according to the 3-momentum p.
// In case the mass was not yet set, the energy is set to correspond to m=0. 
 Set3Vector(p);
 Double_t inv=GetInvariant();
 if (inv<0) SetMass(0.);
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Set4Momentum(Ali4Vector& p)
{
// Set the track parameters according to the 4-momentum p
 Double_t E=p.GetScalar();
 Double_t dE=p.GetResultError();
 Ali3Vector pv=p.Get3Vector();
 SetVector(E,pv);
 SetScalarError(dE);
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetMass(Double_t m,Double_t dm)
{
// Set the particle mass
// The default value for the error dm is 0.
 Double_t inv=pow(m,2);
 Double_t dinv=fabs(2.*m*dm);
 SetInvariant(inv,dinv);
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetCharge(Float_t q)
{
// Set the particle charge
 fQ=q;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Data(TString f)
{
// Provide track information within the coordinate frame f
 Double_t m=GetMass();
 Double_t dm=GetResultError();
 const char* name=GetName();
 const char* title=GetTitle();

 cout << " *" << ClassName() << "::Data*";
 if (strlen(name))  cout << " Name : " << GetName();
 if (strlen(title)) cout << " Title : " << GetTitle();
 cout << endl;
 cout << " Id : " << fUserId << " Code : " << fCode
      << " m : " << m << " dm : " << dm << " Charge : " << fQ
      << " p : " << GetMomentum() << endl;
 cout << " Nhypotheses : " << GetNhypotheses() << " Ndecay-tracks : " << GetNdecay()
      << " Nsignals : " << GetNsignals() << endl;
 Ali4Vector::Data(f); 
} 
///////////////////////////////////////////////////////////////////////////
void AliTrack::List(TString f)
{
// Provide current track and decay level 1 information within coordinate frame f

 Data(f); // Information of the current track

 // Decay products of this track
 AliTrack* td; 
 for (Int_t id=1; id<=GetNdecay(); id++)
 {
  td=GetDecayTrack(id);
  if (td)
  {
   cout << "  ---Level 1 sec. track no. " << id << endl;
   td->Data(f); 
  }
  else
  {
   cout << " *AliTrack::List* Error : Empty decay track slot." << endl; 
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
void AliTrack::ListAll(TString f)
{
// Provide complete track and decay information within the coordinate frame f

 Data(f); // Information of the current track
 if (fBegin) { cout << " Begin-point :"; fBegin->Data(f); }
 if (fEnd)   { cout << " End-point   :"; fEnd->Data(f); }

 Int_t nhyp=GetNhypotheses();
 if (nhyp)
 {
  cout << " List of the " << nhyp << " track hypotheses : " << endl;
  for (Int_t ih=1; ih<=nhyp; ih++)
  {
   AliTrack* tx=GetTrackHypothesis(ih);
   if (tx) tx->Data(f);
  }
 }

 Int_t nsig=GetNsignals();
 if (nsig)
 {
  cout << " List of the corresponding slots for the " << nsig
       << " related signals : " << endl;
  AliPosition r;
  Int_t nrefs,jslot;
  TArrayI slotarr;
  for (Int_t is=1; is<=nsig; is++)
  {
   AliSignal* sx=GetSignal(is);
   if (sx)
   {
    nrefs=sx->GetIndices(this,slotarr,0);
    for (Int_t jref=0; jref<nrefs; jref++)
    {
     jslot=slotarr.At(jref);
     sx->List(jslot);
    }
    r=sx->GetPosition();
    cout << "   Position";
    r.Data(f);
   }
  }
 }

 AliTrack* t=this;
 Dumps(t,1,f); // Information of all decay products
}
//////////////////////////////////////////////////////////////////////////
void AliTrack::Dumps(AliTrack* t,Int_t n,TString f)
{
// Recursively provide the info of all decay levels of this track
 AliTrack* td; 
 for (Int_t id=1; id<=t->GetNdecay(); id++)
 {
  td=t->GetDecayTrack(id);
  if (td)
  {
   cout << "  ---Level " << n << " sec. track no. " << id << endl;
   td->Data(f); 

   Int_t nhyp=td->GetNhypotheses();
   if (nhyp)
   {
    cout << " List of the " << nhyp << " track hypotheses : " << endl;
    for (Int_t ih=1; ih<=nhyp; ih++)
    {
     AliTrack* tx=td->GetTrackHypothesis(ih);
     if (tx) tx->Data(f);
    }
   }

   Int_t nsig=td->GetNsignals();
   if (nsig)
   {
    cout << " List of the " << nsig << " related signals : " << endl;
    for (Int_t is=1; is<=nsig; is++)
    {
     AliSignal* sx=td->GetSignal(is);
     if (sx) sx->Data(f);
    }
   }

   // Go for next decay level of this decay track recursively
   Dumps(td,n+1,f);
  }
  else
  {
   cout << " *AliTrack::Dumps* Error : Empty decay track slot." << endl; 
  }
 }
} 
//////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetMomentum()
{
// Provide the value of the track 3-momentum.
// The error can be obtained by invoking GetResultError() after
// invokation of GetMomentum().
 Double_t norm=fV.GetNorm();
 fDresult=fV.GetResultError();
 return norm;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector AliTrack::Get3Momentum() const
{
// Provide the track 3-momentum
 return (Ali3Vector)Get3Vector();
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetMass()
{
// Provide the particle mass.
// The error can be obtained by invoking GetResultError() after
// invokation of GetMass().
 Double_t inv=GetInvariant();
 Double_t dinv=GetResultError();
 Double_t dm=0;
 if (inv >= 0)
 {
 Double_t m=sqrt(inv);
 if (m) dm=dinv/(2.*m);
 fDresult=dm;
 return m;
 }
 else
 {
  cout << "*AliTrack::GetMass* Unphysical situation m**2 = " << inv << endl;
  cout << " Value 0 will be returned." << endl;
  fDresult=dm;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliTrack::GetCharge() const
{
// Provide the particle charge
 return fQ;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetEnergy()
{
// Provide the particle's energy.
// The error can be obtained by invoking GetResultError() after
// invokation of GetEnergy().
 Double_t E=GetScalar();
 if (E>0)
 {
  return E;
 }
 else
 {
  cout << "*AliTrack::GetEnergy* Unphysical situation E = " << E << endl;
  cout << " Value 0 will be returned." << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Decay(Double_t m1,Double_t m2,Double_t thcms,Double_t phicms)
{
// Perform 2-body decay of current track
// m1     : mass of decay product 1
// m2     : mass of decay product 2
// thcms  : cms theta decay angle (in rad.) of m1
// phicms : cms phi decay angle (in rad.) of m1
 
 Double_t M=GetMass();
 
// Compute the 4-momenta of the decay products in the cms
// Note : p2=p1=pnorm for a 2-body decay
 Double_t e1=0;
 if (M) e1=((M*M)+(m1*m1)-(m2*m2))/(2.*M);
 Double_t e2=0;
 if (M) e2=((M*M)+(m2*m2)-(m1*m1))/(2.*M);
 Double_t pnorm=(e1*e1)-(m1*m1);
 if (pnorm>0.)
 {
  pnorm=sqrt(pnorm);
 }
 else
 {
  pnorm=0;
 }
 
 Double_t a[3];
 a[0]=pnorm;
 a[1]=thcms;
 a[2]=phicms;
 Ali3Vector p;
 p.SetVector(a,"sph");

 Ali4Vector pprim1;
 pprim1.SetVector(e1,p);
 pprim1.SetInvariant(m1*m1);

 Ali4Vector pprim2;
 p*=-1;
 pprim2.SetVector(e2,p);
 pprim2.SetInvariant(m2*m2);

 // Determine boost parameters from the parent particle
 Double_t E=GetEnergy();
 p=Get3Vector();
 Ali4Vector pmu;
 pmu.SetVector(E,p);

 AliBoost q;
 q.Set4Momentum(pmu);
 
 Ali4Vector p1=q.Inverse(pprim1); // Boost decay product 1
 Ali4Vector p2=q.Inverse(pprim2); // Boost decay product 2
 
 // Enter the boosted data into the decay tracks array
 if (fDecays)
 {
  delete fDecays;
  fDecays=0;
 }
 fDecays=new TObjArray(2);
 fDecays->SetOwner();

 fDecays->Add(new AliTrack);
 ((AliTrack*)fDecays->At(0))->Set4Momentum(p1);
 ((AliTrack*)fDecays->At(0))->SetMass(m1);
 fDecays->Add(new AliTrack);
 ((AliTrack*)fDecays->At(1))->Set4Momentum(p2);
 ((AliTrack*)fDecays->At(1))->SetMass(m2);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetNdecay() const
{
// Provide the number of decay produced tracks
 Int_t ndec=0;
 if (fDecays) ndec=fDecays->GetEntries();
 return ndec;
}
///////////////////////////////////////////////////////////////////////////
AliTrack* AliTrack::GetDecayTrack(Int_t j) const
{
// Provide decay produced track number j
// Note : j=1 denotes the first decay track
 if (!fDecays)
 {
  cout << " *AliTrack::GetDecayTrack* No tracks present." << endl;
  return 0;
 }
 else
 {
  if ((j >= 1) && (j <= GetNdecay()))
  {
   return (AliTrack*)fDecays->At(j-1);
  }
  else
  {
   cout << " *AliTrack* decay track number : " << j << " out of range."
        << " Ndec = " << GetNdecay() << endl;
   return 0;  
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::RemoveDecays()
{
// Remove all decay tracks from this track.
 if (fDecays)
 {
  delete fDecays;
  fDecays=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::AddSignal(AliSignal& s)
{
// Relate an AliSignal object to this track.
 if (!fSignals) fSignals=new TObjArray(1);

 // Check if this signal is already stored for this track
 Int_t nsig=GetNsignals();
 for (Int_t i=0; i<nsig; i++)
 {
  if (&s==fSignals->At(i)) return; 
 }

 fSignals->Add(&s);
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::RemoveSignal(AliSignal& s)
{
// Remove related AliSignal object from this track.
 if (fSignals)
 {
  AliSignal* test=(AliSignal*)fSignals->Remove(&s);
  if (test) fSignals->Compress();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::RemoveSignals()
{
// Remove all related AliSignal objects from this track.
 if (fSignals)
 {
  fSignals->Clear();
  delete fSignals;
  fSignals=0;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetNsignals() const
{
// Provide the number of related AliSignals.
 Int_t nsig=0;
 if (fSignals) nsig=fSignals->GetEntries();
 return nsig;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliTrack::GetSignal(Int_t j) const
{
// Provide the related AliSignal number j.
// Note : j=1 denotes the first signal.
 if (!fSignals)
 {
  cout << " *AliTrack::GetSignal* No signals present." << endl;
  return 0;
 }
 else
 {
  if ((j >= 1) && (j <= GetNsignals()))
  {
   return (AliSignal*)fSignals->At(j-1);
  }
  else
  {
   cout << " *AliTrack* signal number : " << j << " out of range."
        << " Nsig = " << GetNsignals() << endl;
   return 0;
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::AddTrackHypothesis(AliTrack& t)
{
// Relate a track hypothesis to this track.
// Note : a private copy of the input track will be made via the Clone()
//        facility.
 if (!fHypotheses)
 {
  fHypotheses=new TObjArray(1);
  fHypotheses->SetOwner();
 }
 fHypotheses->Add(t.Clone());
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::AddTrackHypothesis(Double_t prob,Double_t m,Double_t dm)
{
// Add a track hypothesis by explicitly setting the mass and probability.
// This will affect e.g. the hypothesis track's energy, since the momentum
// and all other attributes will be copied from the current track.
//
// Input arguments :
// ----------------- 
// prob=probalility  m=mass value  dm=error on the mass value.
// The default value for the mass error dm is 0.

 AliTrack t(*this);
 t.RemoveDecays();
 t.RemoveTrackHypotheses();
 t.RemoveSignals();
 t.SetTitle("Mass hypothesis");
 t.SetMass(m,dm);
 t.SetProb(prob);
 AddTrackHypothesis(t);
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::RemoveTrackHypothesis(AliTrack& t)
{
// Remove the specified track hypothesis from this track.
 if (fHypotheses)
 {
  AliTrack* test=(AliTrack*)fHypotheses->Remove(&t);
  if (test) fHypotheses->Compress();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::RemoveTrackHypotheses()
{
// Remove all track hypotheses from this track.
 if (fHypotheses)
 {
  delete fHypotheses;
  fHypotheses=0;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetNhypotheses() const
{
// Provide the number of track hypotheses.
 Int_t nhyp=0;
 if (fHypotheses) nhyp=fHypotheses->GetEntries();
 return nhyp;
}
///////////////////////////////////////////////////////////////////////////
AliTrack* AliTrack::GetTrackHypothesis(Int_t j) const
{
// Provide the j-th track hypothesis.
// Note : j=1 denotes the first hypothesis.
// Default : j=0 ==> Hypothesis with highest probability.

 if (!fHypotheses) return 0;

 Int_t nhyp=GetNhypotheses();

 // Check validity of index j
 if (j<0 || j>nhyp)
 {
   cout << " *AliTrack* hypothesis number : " << j << " out of range."
        << " Nhyp = " << nhyp << endl;
   return 0;
 } 

 AliTrack* t=0;

 if (j==0) // Provide track hypothesis with highest probability
 {
  Float_t prob=0;   
  t=(AliTrack*)fHypotheses->At(0);
  if (t) prob=t->GetProb();
  Float_t probx=0;
  for (Int_t ih=1; ih<nhyp; ih++)
  {
   AliTrack* tx=(AliTrack*)fHypotheses->At(ih);
   if (tx)
   {
    probx=tx->GetProb();
    if (probx > prob) t=tx; 
   }
  }
  return t;
 }
 else // Provide requested j-th track hypothesis
 {
  return (AliTrack*)fHypotheses->At(j-1);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetBeginPoint(AliPosition& p)
{
// Store the position of the track begin-point.
 if (!fBegin)
 {
  fBegin=new AliPositionObj(p);
 }
 else
 {
  fBegin->Load(p);
 }
}
///////////////////////////////////////////////////////////////////////////
AliPosition* AliTrack::GetBeginPoint()
{
// Provide the position of the track begin-point.
 return fBegin;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetEndPoint(AliPosition& p)
{
// Store the position of the track end-point.
 if (!fEnd)
 {
  fEnd=new AliPositionObj(p);
 }
 else
 {
  fEnd->Load(p);
 }
}
///////////////////////////////////////////////////////////////////////////
AliPosition* AliTrack::GetEndPoint()
{
// Provide the position of the track end-point.
 return fEnd;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetMass()
{
// Set the mass and error to the value of the hypothesis with highest prob.

 Double_t m=0,dm=0;

 // Select mass hypothesis with highest probability
 AliTrack* t=GetTrackHypothesis(0);
 if (t) 
 {
  m=t->GetMass();
  dm=t->GetResultError();
  SetMass(m,dm);
 }
 else
 {
  cout << " *AliTrack::SetMass()* No hypothesis present => No action." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetPt()
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
Double_t AliTrack::GetPl()
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
Double_t AliTrack::GetEt()
{
// Provide trans. energy value w.r.t. z-axis.
// The error on the value can be obtained by GetResultError()
// after invokation of GetEt().
 Double_t et=GetScaTrans();

 return et;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetEl()
{
// Provide long. energy value w.r.t. z-axis.
// Note : the returned value can also be negative.
// The error on the value can be obtained by GetResultError()
// after invokation of GetEl().
 Double_t el=GetScaLong();

 return el;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetMt()
{
// Provide transverse mass value w.r.t. z-axis.
// The error on the value can be obtained by GetResultError()
// after invokation of GetMt().
 Double_t pt=GetPt();
 Double_t dpt=GetResultError();
 Double_t m=GetMass();
 Double_t dm=GetResultError();

 Double_t mt=sqrt(pt*pt+m*m);
 Double_t dmt2=0;
 if (mt) dmt2=(pow((pt*dpt),2)+pow((m*dm),2))/(mt*mt);

 fDresult=sqrt(dmt2);
 return mt;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetRapidity()
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
void AliTrack::SetImpactPoint(AliPosition& p,TString q)
{
// Store the position of the impact-point in the plane "q=0".
// Here q denotes one of the axes X, Y or Z.
// Note : The character to denote the axis may be entered in lower or
//        in uppercase.
 Int_t axis=0;
 if (q=="x" || q=="X") axis=1;
 if (q=="y" || q=="Y") axis=2;
 if (q=="z" || q=="Z") axis=3;

 switch (axis)
 {
  case 1: // Impact-point in the plane X=0
   if (!fImpactYZ)
   {
    fImpactYZ=new AliPositionObj(p);
   }
   else
   {
    fImpactYZ->Load(p);
   }
   break;

  case 2: // Impact-point in the plane Y=0
   if (!fImpactXZ)
   {
    fImpactXZ=new AliPositionObj(p);
   }
   else
   {
    fImpactXZ->Load(p);
   }
   break;

  case 3: // Impact-point in the plane Z=0
   if (!fImpactXY)
   {
    fImpactXY=new AliPositionObj(p);
   }
   else
   {
    fImpactXY->Load(p);
   }
   break;

  default: // Unsupported axis
   cout << "*AliTrack::SetImpactPoint* Unsupported axis : " << q << endl
        << " Possible axes are 'X', 'Y' and 'Z'." << endl; 
   break;
 }
}
///////////////////////////////////////////////////////////////////////////
AliPosition* AliTrack::GetImpactPoint(TString q)
{
// Provide the position of the impact-point in the plane "q=0".
// Here q denotes one of the axes X, Y or Z.
// Note : The character to denote the axis may be entered in lower or
//        in uppercase.
 Int_t axis=0;
 if (q=="x" || q=="X") axis=1;
 if (q=="y" || q=="Y") axis=2;
 if (q=="z" || q=="Z") axis=3;

 switch (axis)
 {
  case 1: // Impact-point in the plane X=0
   return fImpactYZ;

  case 2: // Impact-point in the plane Y=0
   return fImpactXZ;

  case 3: // Impact-point in the plane Z=0
   return fImpactXY;

  default: // Unsupported axis
   cout << "*AliTrack::GetImpactPoint* Unsupported axis : " << q << endl
        << " Possible axes are 'X', 'Y' and 'Z'." << endl; 
   return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetId(Int_t id)
{
// Set a user defined unique identifier for this track.
 fUserId=id;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetId() const
{
// Provide the user defined unique identifier of this track.
 return fUserId;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetClosestPoint(AliPosition& p)
{
// Set position p as the point of closest approach w.r.t. some reference
 if (!fClosest)
 {
  fClosest=new AliPositionObj(p);
 }
 else
 {
  fClosest->Load(p);
 }
}
///////////////////////////////////////////////////////////////////////////
AliPosition* AliTrack::GetClosestPoint()
{
// Provide the point of closest approach w.r.t. some reference
 return fClosest;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetChi2(Float_t chi2)
{
// Set the chi-squared value of the track fit.
 if (chi2<0)
 {
  cout << " *AliTrack::SetChi2* Invalid chi2 value : " << chi2 << endl;
 }
 else
 {
  fChi2=chi2;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetNdf(Int_t ndf)
{
// Set the number of degrees of freedom for the track fit.
 if (ndf<0)
 {
  cout << " *AliTrack::SetNdf* Invalid ndf value : " << ndf << endl;
 }
 else
 {
  fNdf=ndf;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliTrack::GetChi2() const
{
// Provide the chi-squared value of the track fit.
 return fChi2;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetNdf() const
{
// Provide the number of degrees of freedom for the track fit.
 return fNdf;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetParticleCode(Int_t code)
{
// Set the user defined particle id code (e.g. the PDF convention).
 fCode=code;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetParticleCode() const
{
// Provide the user defined particle id code.
 return fCode;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetParentTrack(AliTrack* t)
{
// Set pointer to the parent track.
 fParent=t;
}
///////////////////////////////////////////////////////////////////////////
AliTrack* AliTrack::GetParentTrack()
{
// Provide pointer to the parent track.
 return fParent;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetProb(Double_t prob)
{
// Set hypothesis probability for this track.
 fProb=prob;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliTrack::GetProb() const
{
// Provide the hypothesis probability for this track.
 return fProb;
}
///////////////////////////////////////////////////////////////////////////
TObject* AliTrack::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.
// This feature allows e.g. AliJet to store either AliTrack objects or
// objects derived from AliTrack via the AddTrack memberfunction, provided
// these derived classes also have a proper Clone memberfunction. 

 AliTrack* trk=new AliTrack(*this);
 if (name)
 {
  if (strlen(name)) trk->SetName(name);
 }
 return trk;
}
///////////////////////////////////////////////////////////////////////////
