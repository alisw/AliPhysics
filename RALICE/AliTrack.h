#ifndef ALITRACK_H
#define ALITRACK_H
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
// Note : All quantities are in GeV, GeV/c or GeV/c**2
//
//--- NvE 10-jul-1997 UU-SAP Utrecht
//--- Modified : NvE 06-apr-1999 UU-SAP Utrecht to inherit from Ali4Vector
///////////////////////////////////////////////////////////////////////////
 
#include "TObject.h"
#include "TObjArray.h"
 
#include "AliBoost.h"
 
class AliTrack : public TObject,public Ali4Vector
{
 public:
  AliTrack();                       // Default constructor
  ~AliTrack();                      // Destructor
  void Reset();                     // Reset all values to 0
  void Set4Momentum(Ali4Vector& p); // Set track 4-momentum
  void Set3Momentum(Ali3Vector& p); // Set track 3-momentum
  void SetMass(Double_t m);         // Set particle mass
  void SetCharge(Float_t q);        // Set particle charge
  void Info(TString f="car");       // Print track information for coord. frame f
  void List(TString f="car");       // Print track and decay level 1 information for coord. frame f
  void ListAll(TString f="car");    // Print track and all decay level information for coord. frame f
  Ali3Vector Get3Momentum();        // Provide track 3-momentum
  Double_t GetMomentum();           // Provide value of track 3-momentum
  Double_t GetMass();               // Provide particle mass
  Float_t GetCharge();              // Provide particle charge
  Double_t GetEnergy();             // Provide particle total energy
  void Decay(Double_t m1,Double_t m2,Double_t thcms,Double_t phicms); // Perform 2-body decay
  Int_t GetNdecay();                // Provide number of decay products
  AliTrack* GetDecayTrack(Int_t j); // Access to decay produced track number j
 
 protected:
  Double_t fM;        // The mass of the particle
  Float_t fQ;         // The charge of the particle
  Int_t fNdec;        // The number of decay products
  TObjArray* fDecays; // The array of decay produced tracks for output

 private:
  void Dump(AliTrack* t,Int_t n,TString f); // Recursively print all decay levels
 
 ClassDef(AliTrack,1) // Class definition to enable ROOT I/O
};
#endif
