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
 
AliTrack::AliTrack() : TObject(),Ali4Vector()
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
 fMasses=0;
 fDmasses=0;
 fPmasses=0;
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
// Destructor to delete memory allocated for decay tracks array
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
 if (fMasses)
 {
  delete fMasses;
  fMasses=0;
 }
 if (fDmasses)
 {
  delete fDmasses;
  fDmasses=0;
 }
 if (fPmasses)
 {
  delete fPmasses;
  fPmasses=0;
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
AliTrack::AliTrack(AliTrack& t) : TObject(t),Ali4Vector(t)
{
// Copy constructor
 Init();

 fQ=t.fQ;
 fNdec=t.fNdec;
 fNsig=t.fNsig;
 fNmasses=t.fNmasses;
 if (fNmasses)
 {
  fMasses=new TArrayD(*(t.fMasses));
  fDmasses=new TArrayD(*(t.fDmasses));
  fPmasses=new TArrayD(*(t.fPmasses));
 }
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

 if (fNdec)
 {
  fDecays=new TObjArray(fNdec);
  fDecays->SetOwner();
  for (Int_t it=1; it<=fNdec; it++)
  {
   AliTrack* tx=t.GetDecayTrack(it);
   fDecays->Add(new AliTrack(*tx));
  }
 }

 if (fNsig)
 {
  fSignals=new TObjArray(fNsig);
  for (Int_t is=1; is<=fNsig; is++)
  {
   AliSignal* sx=t.GetSignal(is);
   fSignals->Add(sx);
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
 fNdec=0;
 fNsig=0;
 fNmasses=0;
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
 if (fMasses)
 {
  delete fMasses;
  fMasses=0;
 }
 if (fDmasses)
 {
  delete fDmasses;
  fDmasses=0;
 }
 if (fPmasses)
 {
  delete fPmasses;
  fPmasses=0;
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
void AliTrack::Data(TString f)
{
// Provide track information within the coordinate frame f
 Double_t m=GetMass();
 Double_t dm=GetResultError();
 cout << " *AliTrack::Data* Id : " << fUserId << " Code : " << fCode
      << " Mass : " << m << " error : " << dm << " Charge : " << fQ
      << " Momentum : " << GetMomentum() << " Nmass hyp. : " << fNmasses
      << " Ntracks : " << fNdec << " Nsignals : " << fNsig << endl;
 for (Int_t i=0; i<fNmasses; i++)
 {
  cout << " Mass hypothesis " << (i+1) << " Mass : " << fMasses->At(i)
       << " error : " << fDmasses->At(i) << " prob. : " << fPmasses->At(i)
       << endl;
 }
 Ali4Vector::Data(f); 
} 
///////////////////////////////////////////////////////////////////////////
void AliTrack::List(TString f)
{
// Provide current track and decay level 1 information within coordinate frame f

 Data(f); // Information of the current track

 // Decay products of this track
 AliTrack* td; 
 for (Int_t id=1; id<=fNdec; id++)
 {
  td=GetDecayTrack(id);
  if (td)
  {
   cout << "  ---Level 1 sec. track no. " << id << endl;
   td->Data(f); 
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

 Data(f); // Information of the current track
 if (fBegin) { cout << " Begin-point :"; fBegin->Data(f); }
 if (fEnd)   { cout << " End-point   :"; fEnd->Data(f); }
 for (Int_t is=1; is<=GetNsignals(); is++)
 {
  ((AliSignal*)GetSignal(is))->Data(f);
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
   td->Data(f); 
   for (Int_t is=1; is<=td->GetNsignals(); is++)
   {
    ((AliSignal*)td->GetSignal(is))->Data(f);
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
 Double_t norm=fV.GetNorm();
 fDresult=fV.GetResultError();
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
  delete fDecays;
  fDecays=0;
 }
 fDecays=new TObjArray();
 fDecays->SetOwner();

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
 if (!fDecays)
 {
  cout << " *AliTrack::GetDecayTrack* No tracks present." << endl;
  return 0;
 }
 else
 {
  if ((j >= 1) && (j <= fNdec))
  {
   return (AliTrack*)fDecays->At(j-1);
  }
  else
  {
   cout << " *AliTrack* decay track number : " << j << " out of range."
        << " Ndec = " << fNdec << endl;
   return 0;  
  }
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
 if (!fSignals)
 {
  cout << " *AliTrack::GetSignal* No signals present." << endl;
  return 0;
 }
 else
 {
  if ((j >= 1) && (j <= fNsig))
  {
   return (AliSignal*)fSignals->At(j-1);
  }
  else
  {
   cout << " *AliTrack* signal number : " << j << " out of range."
        << " Nsig = " << fNsig << endl;
   return 0;
  }
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
void AliTrack::AddMassHypothesis(Double_t prob,Double_t m,Double_t dm)
{
// Add a mass hypothesis for this current track.
// prob=probalility  m=mass value  dm=error on the mass value.
// The default value for the mass error dm is 0.
 if (!fMasses) fMasses=new TArrayD();
 if (!fDmasses) fDmasses=new TArrayD();
 if (!fPmasses) fPmasses=new TArrayD();

 fNmasses++;
 fMasses->Set(fNmasses);
 fDmasses->Set(fNmasses);
 fPmasses->Set(fNmasses);

 fMasses->AddAt(m,fNmasses-1);
 fDmasses->AddAt(dm,fNmasses-1);
 fPmasses->AddAt(prob,fNmasses-1);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetNMassHypotheses()
{
// Provide the number of mass hypotheses for this track.
 return fNmasses;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetMassHypothesis(Int_t j)
{
// Provide the mass of the jth hypothesis for this track.
// Note : the first hypothesis is indicated by j=1.
// Default : j=0 ==> Hypothesis with highest probability.
// The error on the mass can be obtained by invoking GetResultError()
// after invokation of GetMassHypothesis(j).

 Double_t m=0,dm=0,prob=0;

 // Check validity of index j
 if (j<0 || j>fNmasses)
 {
  cout << " *AliTrack::GetMassHypothesis* Invalid index j : " << j
       << " Number of mass hypotheses : " << fNmasses << endl;
  fDresult=0;
  return 0;
 }

 // Select mass hypothesis with highest probability
 if (j==0) 
 {
  if (fNmasses) 
  {
   m=fMasses->At(0);
   dm=fDmasses->At(0);
   prob=fPmasses->At(0);
   for (Int_t i=1; i<fNmasses; i++)
   {
    if (fPmasses->At(i)>prob)
    {
     m=fMasses->At(i);
     dm=fDmasses->At(i);
    }
   }
  }
  fDresult=dm;
  return m;  
 }

 // Provide data of requested mass hypothesis
 m=fMasses->At(j-1);
 fDresult=fDmasses->At(j-1);
 return m;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliTrack::GetMassHypothesisProb(Int_t j)
{
// Provide the probability of the jth hypothesis for this track.
// Note : the first hypothesis is indicated by j=1.
// Default : j=0 ==> Hypothesis with highest probability.

 Double_t prob=0;

 // Check validity of index j
 if (j<0 || j>fNmasses)
 {
  cout << " *AliTrack::GetMassHypothesisProb* Invalid index j : " << j
       << " Number of mass hypotheses : " << fNmasses << endl;
  return 0;
 }

 // Select mass hypothesis with highest probability
 if (j==0) 
 {
  if (fNmasses) 
  {
   prob=fPmasses->At(0);
   for (Int_t i=1; i<fNmasses; i++)
   {
    if (fPmasses->At(i)>prob) prob=fPmasses->At(i);
   }
  }
  return prob;  
 }

 // Provide probability of requested mass hypothesis
 prob=fPmasses->At(j-1);
 return prob;
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetMass()
{
// Set the mass and error to the value of the hypothesis with highest prob.

 Double_t m=0,dm=0,prob=0;

 // Select mass hypothesis with highest probability
 if (fNmasses) 
 {
  m=fMasses->At(0);
  dm=fDmasses->At(0);
  prob=fPmasses->At(0);
  for (Int_t i=1; i<fNmasses; i++)
  {
   if (fPmasses->At(i)>prob)
   {
    m=fMasses->At(i);
    dm=fDmasses->At(i);
   }
  }
  SetMass(m,dm);
 }
 else
 {
  cout << " *AliTrack::SetMass()* No hypothesis present => No action." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::RemoveMassHypothesis(Int_t j)
{
// Remove the jth mass hypothesis for this track.
// Note : the first hypothesis is indicated by j=1.

 if (j<=0 || j>fNmasses) // Check validity of index j
 {
  cout << " *AliTrack::RemoveMassHypothesis* Invalid index j : " << j
       << " Number of mass hypotheses : " << fNmasses << endl;
 }
 else
 {
  if (j != fNmasses)
  {
   fMasses->AddAt(fMasses->At(fNmasses-1),j-1);
   fDmasses->AddAt(fDmasses->At(fNmasses-1),j-1);
   fPmasses->AddAt(fPmasses->At(fNmasses-1),j-1);
  }
  fMasses->AddAt(0,fNmasses-1);
  fDmasses->AddAt(0,fNmasses-1);
  fPmasses->AddAt(0,fNmasses-1);
  fNmasses--;
  fMasses->Set(fNmasses);
  fDmasses->Set(fNmasses);
  fPmasses->Set(fNmasses);
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
Double_t AliTrack::GetMt(Int_t j)
{
// Provide transverse mass value w.r.t. z-axis and jth mass hypothesis.
// Note : the first hypothesis is indicated by j=1.
//        j=0 ==> Hypothesis with highest probability.
// The error on the value can be obtained by GetResultError()
// after invokation of GetMt(j).
 Double_t pt=GetPt();
 Double_t dpt=GetResultError();
 Double_t m=GetMassHypothesis(j);
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
Int_t AliTrack::GetId()
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
Float_t AliTrack::GetChi2()
{
// Provide the chi-squared value of the track fit.
 return fChi2;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliTrack::GetNdf()
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
Int_t AliTrack::GetParticleCode()
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
