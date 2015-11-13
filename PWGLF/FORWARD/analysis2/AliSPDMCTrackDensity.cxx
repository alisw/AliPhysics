#include "AliSPDMCTrackDensity.h"
#include "AliMCEvent.h"
#include "AliTrackReference.h"
#include "AliForwardUtil.h"
#include <TMath.h>
#include <AliLog.h>
#include <TROOT.h>
#include <TH2D.h>
#include <iostream>

//____________________________________________________________________
AliSPDMCTrackDensity::AliSPDMCTrackDensity()
  : AliBaseMCTrackDensity(), 
    fMinR(3.5), 
    fMaxR(4.5),
    fMinZ(-15), // -14.1), 
    fMaxZ(+15), // +14.1)
    fStored(0), 
    fOutput(0)
{
  // Default constructor 
}

//____________________________________________________________________
AliSPDMCTrackDensity::AliSPDMCTrackDensity(const char*)
  : AliBaseMCTrackDensity("spdMCTrackDensity"), 
    fMinR(3.5), 
    fMaxR(4.5),
    fMinZ(-14.1), 
    fMaxZ(+14.1),
    fStored(0), 
    fOutput(0)
{
  // Normal constructor constructor 
}

//____________________________________________________________________
AliSPDMCTrackDensity::AliSPDMCTrackDensity(const AliSPDMCTrackDensity& o)
  : AliBaseMCTrackDensity(o),
    fMinR(o.fMinR), 
    fMaxR(o.fMaxR),
    fMinZ(o.fMinZ), 
    fMaxZ(o.fMaxZ), 
    fStored(o.fStored), 
    fOutput(o.fOutput)
{
  // Normal constructor constructor 
}

//____________________________________________________________________
AliSPDMCTrackDensity&
AliSPDMCTrackDensity::operator=(const AliSPDMCTrackDensity& o)
{
  // Assignment operator 
  if (&o == this) return *this;
  AliBaseMCTrackDensity::operator=(o);
  fMinR             = o.fMinR;
  fMaxR             = o.fMaxR;
  fMinZ             = o.fMinZ;
  fMaxZ             = o.fMaxZ;
  fStored           = o.fStored;
  fOutput           = o.fOutput;

  return *this;
}

//____________________________________________________________________
Int_t
AliSPDMCTrackDensity::GetDetectorId() const
{
  return AliTrackReference::kITS;
}


//____________________________________________________________________
void
AliSPDMCTrackDensity::BeginTrackRefs()
{
  fStored = 0;
}

//____________________________________________________________________
Bool_t
AliSPDMCTrackDensity::CheckTrackRef(AliTrackReference*   ref) const
{
  // Get radius and z where the track reference was made 
  Double_t r = ref->R();
  Double_t z = ref->Z();
  if (r > fMaxR || r < fMinR) return false;
  if (z > fMaxZ || z < fMinZ) return false;

  return true;
}
//____________________________________________________________________
AliTrackReference*
AliSPDMCTrackDensity::ProcessRef(AliMCParticle*       /*particle*/,
				 const AliMCParticle* /*mother*/,
				 AliTrackReference*   ref)
{
  if (fStored) return 0;

  return fStored = ref;
}

//____________________________________________________________________
Double_t
AliSPDMCTrackDensity::StoreParticle(AliMCParticle* particle, 
				    const AliMCParticle* mother, 
				    AliTrackReference*   ref) const
{
  Double_t w = AliBaseMCTrackDensity::StoreParticle(particle, mother, ref);
  // Double_t r = ref->R();
  Double_t x = ref->X()-fIP.X();
  Double_t y = ref->Y()-fIP.Y();
  Double_t z = ref->Z()-fIP.Z();
  Double_t r = TMath::Sqrt(x*x+y*y);
  
  Double_t zr = z; // -fVz;
  Double_t th = TMath::ATan2(r,zr);
  if (th < 0) th += 2*TMath::Pi();
  Double_t et = -TMath::Log(TMath::Tan(th/2));
  Double_t ph = TMath::ATan2(y,x);
  if (ph < 0) ph += 2*TMath::Pi();
  fOutput->Fill(et,ph,w);

  return w;
}


//____________________________________________________________________
Bool_t
AliSPDMCTrackDensity::Calculate(const AliMCEvent& event, 
				const TVector3&   ip, 
				TH2D&             output, 
				TH2D*             primary)
{
  // 
  // Filter the input kinematics and track references, using 
  // some of the ESD information
  // 
  // Parameters:
  //    input   Input ESD event
  //    event   Input MC event
  //    vz      Vertex position 
  //    output  Output ESD-like object
  //    primary Per-event histogram of primaries 
  //
  // Return:
  //    True on succes, false otherwise 
  //
  fOutput = &output;

  return ProcessTracks(event, ip, primary);
}
#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)

#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)
//____________________________________________________________________
void
AliSPDMCTrackDensity::Print(Option_t* option) const 
{
  AliBaseMCTrackDensity::Print(option);
  gROOT->IncreaseDirLevel();
  PF("R range", "[%f,%f]", fMinR, fMaxR);
  PF("Z range", "[%f,%f]", fMinZ, fMaxZ);
  gROOT->DecreaseDirLevel();
}

//____________________________________________________________________
//
// EOF
//
