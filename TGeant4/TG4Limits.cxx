// $Id$
// Category: global
//
// See the class description in the header file.

#include "TG4Limits.h"
#include "TG4CutVector.h"
#include "TG4FlagVector.h"

TG4Limits::TG4Limits()
  : G4UserLimits(),              
    // default values of G4UserLimits data members are set: 
    // fMaxStep (DBL_MAX), fMaxTrack(DBL_MAX),fMaxTime(DBL_MAX),
    // fMinEkine(0.), fMinRange(0.)
    fIsCut(false),
    fIsFlag(false) 
{
//
  fCutVector = new TG4CutVector();
  fFlagVector = new TG4FlagVector();
}

TG4Limits::TG4Limits(const TG4Limits& right)
  : G4UserLimits(right) 
{
//    
  fCutVector = new TG4CutVector(*right.fCutVector);
  fFlagVector = new TG4FlagVector(*right.fFlagVector);
}  

TG4Limits::~TG4Limits() {
//
  delete fCutVector;
  delete fFlagVector;
}

// operators

TG4Limits& TG4Limits::operator=(const TG4Limits& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  G4UserLimits::operator=(right);

  *fCutVector  = *right.fCutVector;
  *fFlagVector = *right.fFlagVector;

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

void TG4Limits::SetG3Cut(TG3Cut g3Cut, G4double cutValue)
{
// Sets the cut value for the specified cut.
// ---

  fCutVector->SetG3Cut(g3Cut, cutValue);
  fIsCut = true;
}

void TG4Limits::SetG3Flag(TG3Flag g3Flag, G4double flagValue)
{
// Sets the process control value for the specified flag.
// ---

  fFlagVector->SetG3Flag(g3Flag, flagValue);
  if (flagValue - kUnset > 0.01) fIsFlag = true;
}

void TG4Limits::SetG3DefaultCuts()
{
// Sets the G3 default cut values for all cuts.
// ---

  fCutVector->SetG3Defaults();
  fIsCut = true;  
}

void TG4Limits::SetG3DefaultFlags()
{
// Sets the G3 default process control values for all flags.
// ---

  fFlagVector->SetG3Defaults();
  fIsFlag = true;
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

G4int TG4Limits::GetFlag(G4VProcess* process) const 
{
// Returns the flag value for the particle associated with
// the specified process.
// ---

  if (fIsFlag)
    return fFlagVector->GetFlag(process);
  else 
    return kUnset;
}
