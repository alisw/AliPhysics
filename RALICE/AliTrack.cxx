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
// t1.Info();
// t2.Info();
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
//- Modified: NvE 29-oct-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliTrack.h"
 
ClassImp(AliTrack) // Class implementation to enable ROOT I/O
 
AliTrack::AliTrack()
{
// Default constructor
// All variables initialised to 0
 fDecays=0;
 fSignals=0;
 Reset();
}
///////////////////////////////////////////////////////////////////////////
AliTrack::~AliTrack()
{
// Destructor to delete memory allocated for decay tracks array
 if (fDecays)
 {
  fDecays->Delete();
  delete fDecays;
  fDecays=0;
 }
 if (fSignals)
 {
  fSignals->Clear();
  delete fSignals;
  fSignals=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Reset()
{
// Reset all variables to 0
 fQ=0;
 fNdec=0;
 fNsig=0;
 Double_t a[4]={0,0,0,0};
 SetVector(a,"sph");
 if (fDecays)
 {
  fDecays->Delete();
  delete fDecays;
  fDecays=0;
 }
 if (fSignals)
 {
  fSignals->Clear();
  delete fSignals;
  fSignals=0;
 }
 Double_t b[3]={0,0,0};
 fBegin.SetPosition(b,"sph");
 fEnd.SetPosition(b,"sph");
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Set3Momentum(Ali3Vector& p)
{
// Set the track parameters according to the 3-momentum p
 Set3Vector(p);
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
void AliTrack::Info(TString f)
{
// Provide track information within the coordinate frame f
 Double_t m=GetMass();
 Double_t dm=GetResultError();
 cout << " *AliTrack::Info* Mass : " << m
      << " error : " << dm << " Charge : " << fQ
      << " Momentum : " << GetMomentum() << " Ntracks : " << fNdec
      << " Nsignals : " << fNsig << endl;
 Ali4Vector::Info(f); 
} 
///////////////////////////////////////////////////////////////////////////
void AliTrack::List(TString f)
{
// Provide current track and decay level 1 information within coordinate frame f

 Info(f); // Information of the current track

 // Decay products of this track
 AliTrack* td; 
 for (Int_t id=1; id<=fNdec; id++)
 {
  td=GetDecayTrack(id);
  if (td)
  {
   cout << "  ---Level 1 sec. track no. " << id << endl;
   td->Info(f); 
  }
  else
  {
   cout << " *AliTrack::List* Error : No decay track present." << endl; 
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
void AliTrack::ListAll(TString f)
{
// Provide complete track and decay information within the coordinate frame f

 Info(f); // Information of the current track
 cout << " Begin-point :"; fBegin.Info(f);
 cout << " End-point   :"; fEnd.Info(f);
 for (Int_t is=1; is<=GetNsignals(); is++)
 {
  ((AliSignal*)GetSignal(is))->Info(f);
 }

 AliTrack* t=this;
 Dump(t,1,f); // Information of all decay products
}
//////////////////////////////////////////////////////////////////////////
void AliTrack::Dump(AliTrack* t,Int_t n,TString f)
{
// Recursively provide the info of all decay levels of this track
 AliTrack* td; 
 for (Int_t id=1; id<=t->GetNdecay(); id++)
 {
  td=t->GetDecayTrack(id);
  if (td)
  {
   cout << "  ---Level " << n << " sec. track no. " << id << endl;
   td->Info(f); 
   for (Int_t is=1; is<=td->GetNsignals(); is++)
   {
    ((AliSignal*)td->GetSignal(is))->Info(f);
   }

   // Go for next decay level of this decay track recursively
   Dump(td,n+1,f);
  }
  else
  {
   cout << " *AliTrack::Dump* Error : No decay track present." << endl; 
  }
 }
} 
//////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetMomentum()
{
// Provide the value of the track 3-momentum.
// The error can be obtained by invoking GetResultError() after
// invokation of GetMomentum().

// Ali3Vector p=Get3Vector();
// return sqrt(p.Dot(p));
 Double_t norm=fV.GetNorm();
 return norm;
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector AliTrack::Get3Momentum()
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
Float_t AliTrack::GetCharge()
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
 
 fNdec=2; // it's a 2-body decay

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
  fDecays->Delete();
  delete fDecays;
 }
 fDecays=new TObjArray();

 fDecays->Add(new AliTrack);
 ((AliTrack*)fDecays->At(0))->Set4Momentum(p1);
 ((AliTrack*)fDecays->At(0))->SetMass(m1);
 fDecays->Add(new AliTrack);
 ((AliTrack*)fDecays->At(1))->Set4Momentum(p2);
 ((AliTrack*)fDecays->At(1))->SetMass(m2);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetNdecay()
{
// Provide the number of decay produced tracks
 return fNdec;
}
///////////////////////////////////////////////////////////////////////////
AliTrack* AliTrack::GetDecayTrack(Int_t j)
{
// Provide decay produced track number j
// Note : j=1 denotes the first decay track
 if ((j >= 1) && (j <= fNdec))
 {
  return (AliTrack*)fDecays->At(j-1);
 }
 else
 {
  cout << " *AliTrack* decay track number : " << j << " out of range." << endl;
  cout << " -- Decay track number 1 (if any) returned." << endl;
  return (AliTrack*)fDecays->At(0);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::AddSignal(AliSignal& s)
{
// Relate an AliSignal object to this track.
 if (!fSignals) fSignals=new TObjArray();
 fNsig++;
 fSignals->Add(&s);
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::RemoveSignal(AliSignal& s)
{
// Remove related AliSignal object to this track.
 if (fSignals)
 {
  AliSignal* test=(AliSignal*)fSignals->Remove(&s);
  if (test)
  {
   fNsig--;
   fSignals->Compress();
  }
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetNsignals()
{
// Provide the number of related AliSignals.
 return fNsig;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliTrack::GetSignal(Int_t j)
{
// Provide the related AliSignal number j.
// Note : j=1 denotes the first signal.
 if ((j >= 1) && (j <= fNsig))
 {
  return (AliSignal*)fSignals->At(j-1);
 }
 else
 {
  cout << " *AliTrack* signal number : " << j << " out of range." << endl;
  cout << " -- Signal number 1 (if any) returned." << endl;
  return (AliSignal*)fDecays->At(0);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetBeginPoint(AliPosition p)
{
// Store the position of the track begin-point.
 fBegin=p;
}
///////////////////////////////////////////////////////////////////////////
AliPosition AliTrack::GetBeginPoint()
{
// Provide the position of the track begin-point.
 return fBegin;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetEndPoint(AliPosition p)
{
// Store the position of the track end-point.
 fEnd=p;
}
///////////////////////////////////////////////////////////////////////////
AliPosition AliTrack::GetEndPoint()
{
// Provide the position of the track end-point.
 return fEnd;
}
///////////////////////////////////////////////////////////////////////////
