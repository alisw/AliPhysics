// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4FlagVector.h"
#include "TG4CutVector.h"
#include "TG4G3Defaults.h"
#include "TG4Globals.h"

#include <G4ParticleDefinition.hh>
#include <G4VProcess.hh>

#include <math.h>

TG4FlagVector::TG4FlagVector()
{
  // initialize fFlagVector 
  fFlagVector = new TG3FlagVector;
  for (G4int i=0; i<kNoG3Flags; i++) fFlagVector->insert(kUnset); 
}

TG4FlagVector::TG4FlagVector(const TG4FlagVector& right)
{
  // copy fFlagVector 
  fFlagVector = new TG3FlagVector;
  for (G4int i=0; i<kNoG3Flags; i++) {
    fFlagVector->insert((*right.fFlagVector)[i]);
  }   
}

TG4FlagVector::~TG4FlagVector() {
//
  delete fFlagVector;
}

// operators

TG4FlagVector& TG4FlagVector::operator=(const TG4FlagVector& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // initialize fFlagVector 
  fFlagVector->clear();
  for (G4int i=0; i<kNoG3Flags; i++) {
    fFlagVector->insert((*right.fFlagVector)[i]);
  }
  
  return *this;   
}  

G4double TG4FlagVector::operator[](G4int index) const
{
//
  if (index < kNoG3Flags)
    return (*fFlagVector)[index];
  else {
    TG4Globals::Exception(
      "TG4FlagVector::operator[]: index out of the vector scope");
    return 0.;  
  }    
}  

// public methods

void TG4FlagVector::SetG3Flag(TG3Flag g3Flag, G4double flagValue)
{
// Sets the flagValue for the specified flag.
// ---

  if (g3Flag<kNoG3Flags) {
    // conversion G4double -> G3FlagValue
    if (abs(flagValue - kUnset) < 0.01) {
        (*fFlagVector)[g3Flag] = kUnset ;
     }	 
     else if (abs(flagValue - kInActivate) < 0.01) {
        (*fFlagVector)[g3Flag] = kInActivate; 
     }
     else if (abs(flagValue - kActivate) < 0.01) {
        (*fFlagVector)[g3Flag] = kActivate; 
     }	
     else if (abs(flagValue - kActivate2) < 0.01) {
        (*fFlagVector)[g3Flag] = kActivate2; 
     }	     	  
     else {
      G4String text = "TG4FlagVector::SetG3Flag:\n ";
      text = text + "Inconsistent/Not-yet-implemented flag has been ignored.";
      TG4Globals::Warning(text);
     }	
  }
}

void TG4FlagVector::SetG3Defaults()
{
// Sets G3 default values for all flags.
// ---

  for (G4int i=0; i<kNoG3Flags; i++) {
   (*fFlagVector)[i] = TG4G3Defaults::FlagValue(i);
  } 
}

G4int TG4FlagVector::GetFlag(G4VProcess* process) const 
{
// Returns the flag value for the particle associated with
// the specified process.
// ---

  G4String name = process->GetProcessName();
  if       (name == "conv")    return (*fFlagVector)(kPAIR);
  else if  (name == "compt")   return (*fFlagVector)(kCOMP);
  else if  (name == "phot")    return (*fFlagVector)(kPHOT);
  // else if (name == "??")  return (*fFlagVector)(kPFIS);
  else if ((name == "eIoni") || 
           (name == "IeIoni") || 
	   (name == "eIoni+") ||
           (name == "MuIoni") || 
	   (name == "IMuIonisation") ||
	   (name == "hIoni") || 
	   (name == "IhIoni"))    
                               return (*fFlagVector)(kDRAY); 
  else if  (name == "annihil") return (*fFlagVector)(kANNI);
  else if ((name == "eBrem") || 
           (name == "eBrem+") || 
	   (name == "IeBrems") || 
           (name == "MuBrems") || 
	   (name == "IMuBremsstrahlung"))   
                              return (*fFlagVector)(kBREM);
  // else if (name == "??")  return (*fFlagVector)(kHADR);
  else if (name == "MuNucl")  return (*fFlagVector)(kMUNU);
  else if (name == "Decay")   return (*fFlagVector)(kDCAY);
  // else if (name == "??")  return (*fFlagVector)(kLOSS);
     // !!! not yet implemented 
  else if ((name == "msc") || 
           (name == "Imsc"))  return (*fFlagVector)(kMULS);
  else return kUnset;
}
