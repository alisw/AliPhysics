#include "AliVertex.h"
 
ClassImp(AliVertex) // Class implementation to enable ROOT I/O
 
AliVertex::AliVertex()
{
// Default constructor
// All variables initialised to 0
// Initial maximum number of tracks is set to the default value
// Initial maximum number of sec. vertices is set to the default value
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
