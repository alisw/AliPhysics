#include "AliTrack.h"
 
ClassImp(AliTrack) // Class implementation to enable ROOT I/O
 
AliTrack::AliTrack()
{
// Default constructor
// All variables initialised to 0
 fDecays=0;
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
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Reset()
{
// Reset all variables to 0
 fM=0;
 fQ=0;
 fNdec=0;
 Double_t a[4]={0,0,0,0};
 SetVector(a,"sph");
 if (fDecays)
 {
  fDecays->Delete();
  delete fDecays;
  fDecays=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Set3Momentum(Ali3Vector& p)
{
// Set the track parameters according to the 3-momentum p
 Double_t E=sqrt(p.Dot(p)+fM*fM);
 SetVector(E,p);
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::Set4Momentum(Ali4Vector& p)
{
// Set the track parameters according to the 4-momentum p
 Double_t E=p.GetScalar();
 Ali3Vector pv=p.Get3Vector();
 SetVector(E,pv);

 Double_t m2=p.Dot(p);
 fM=0;
 if (m2 > 0.) fM=sqrt(m2);
}
///////////////////////////////////////////////////////////////////////////
void AliTrack::SetMass(Double_t m)
{
// Set the particle mass
 fM=m;
 Ali3Vector p=Get3Vector();
 Double_t E=sqrt(p.Dot(p)+fM*fM);
 SetVector(E,p);
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
 cout << " *AliTrack::Info* Mass : " << fM << " Charge : " << fQ
      << " Momentum : " << GetMomentum() << " Ntracks : " << fNdec << endl;
 cout << " ";
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
   cout << " ";
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
   cout << " ";
   td->Info(f); 

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
// Provide the value of the track 3-momentum
 Ali3Vector p=Get3Vector();
 return sqrt(p.Dot(p));
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
// Provide the particle mass
 return fM;
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
// Provide the particle's energy
 return GetScalar();
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
 
// Compute the 4-momenta of the decay products in the cms
// Note : p2=p1=pnorm for a 2-body decay
 Double_t e1=((fM*fM)+(m1*m1)-(m2*m2))/(2.*fM);
 Double_t e2=((fM*fM)+(m2*m2)-(m1*m1))/(2.*fM);
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

 Ali4Vector pprim2;
 p*=-1;
 pprim2.SetVector(e2,p);

 // Determine boost parameters from the parent particle
 Double_t E=GetScalar();
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
 fDecays->Add(new AliTrack);
 ((AliTrack*)fDecays->At(1))->Set4Momentum(p2);
 
// Set the mass values to m1 and m2 to omit roundoff errors
 ((AliTrack*)fDecays->At(0))->SetMass(m1);
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
