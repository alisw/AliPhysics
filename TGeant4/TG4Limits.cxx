// $Id$
// Category: global
//
// See the class description in the header file.

#include "TG4Limits.h"
#include "TG4G3CutVector.h"
#include "TG4G3ControlVector.h"

TG4Limits::TG4Limits()
  : G4UserLimits(),              
    // default values of G4UserLimits data members are set: 
    // fMaxStep (DBL_MAX), fMaxTrack(DBL_MAX),fMaxTime(DBL_MAX),
    // fMinEkine(0.), fMinRange(0.)
    fIsCut(false),
    fIsControl(false) 
{
//
  fCutVector = new TG4G3CutVector();
  fControlVector = new TG4G3ControlVector();
}

TG4Limits::TG4Limits(const TG4Limits& right)
  : G4UserLimits(right) 
{
//    
  fCutVector = new TG4G3CutVector(*right.fCutVector);
  fControlVector = new TG4G3ControlVector(*right.fControlVector);
}  

TG4Limits::~TG4Limits() {
//
  delete fCutVector;
  delete fControlVector;
}

// operators

TG4Limits& TG4Limits::operator=(const TG4Limits& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  G4UserLimits::operator=(right);

  *fCutVector  = *right.fCutVector;
  *fControlVector = *right.fControlVector;

  return *this;  
}    
          
// private methods

G4double TG4Limits::GetUserMinEkine(const G4Track& track)
{
// Returns the kinetic energy cut for the particle
// associated with the specified track.
// ---

  if (fIsCut)
    return fCutVector->GetMinEkine(track);
  else 
    return fMinEkine;
}

// public methods

void TG4Limits::SetG3Cut(TG4G3Cut cut, G4double cutValue)
{
// Sets the cut value for the specified cut.
// ---

  fCutVector->SetG3Cut(cut, cutValue);
  fIsCut = true;
}

void TG4Limits::SetG3Control(TG4G3Control control, G4double flagValue)
{
// Sets the process control value for the specified flag.
// ---

  fControlVector->SetG3Control(control, flagValue);
  if (flagValue - kUnset > 0.01) fIsControl = true;
}

void TG4Limits::SetG3DefaultCuts()
{
// Sets the G3 default cut values for all cuts.
// ---

  fCutVector->SetG3Defaults();
  fIsCut = true;  
}

void TG4Limits::SetG3DefaultControls()
{
// Sets the G3 default process control values for all flags.
// ---

  fControlVector->SetG3Defaults();
  fIsControl = true;
}

G4double TG4Limits::GetMinEkineForGamma(const G4Track& track) const
{
// Returns the cut value for gamma.
// ---

  if (fIsCut)
    return fCutVector->GetMinEkine(track);
  else 
    return fMinEkine;
}

G4double TG4Limits::GetMinEkineForElectron(const G4Track& track) const
{
// Returns the cut value for e-.
// ---

  if (fIsCut)
    return fCutVector->GetMinEkineForElectron(track);
  else 
    return fMinEkine;
}

G4double TG4Limits::GetMinEkineForHadron(const G4Track& track) const
{
// Returns the cut value for charged hadron.
// ---

  if (fIsCut)
    return fCutVector->GetMinEkineForHadron(track);
  else 
    return fMinEkine;
}

G4double TG4Limits::GetMinEkineForNeutralHadron(const G4Track& track) const
{
// Returns the cut value for neutral hadron.
// ---

  if (fIsCut)
    return fCutVector->GetMinEkineForNeutralHadron(track);
  else 
    return fMinEkine;
}

G4double TG4Limits::GetMinEkineForMuon(const G4Track& track) const
{
// Returns the cut value for neutral muon.
// ---

  if (fIsCut)
    return fCutVector->GetMinEkineForMuon(track);
  else 
    return fMinEkine;
}

G4double TG4Limits::GetMinEkineForOther(const G4Track& track) const
{
 // Returns 0.
// ---

 return fMinEkine;
}

G4int TG4Limits::GetControl(G4VProcess* process) const 
{
// Returns the flag value for the particle associated with
// the specified process.
// ---

  if (fIsControl)
    return fControlVector->GetControl(process);
  else 
    return kUnset;
}
