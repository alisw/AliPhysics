// $Id$
// Category: global
//
// See the class description in the header file.

#include "TG4G3ControlVector.h"
#include "TG4G3Defaults.h"
#include "TG4Globals.h"

#include <G4VProcess.hh>

#include <math.h>

TG4G3ControlVector::TG4G3ControlVector()
{
  // initialize fControlVector 
  fControlVector = new TG4ControlValueVector;
  for (G4int i=0; i<kNoG3Controls; i++) fControlVector->insert(kUnset); 
}

TG4G3ControlVector::TG4G3ControlVector(const TG4G3ControlVector& right)
{
  // copy fControlVector 
  fControlVector = new TG4ControlValueVector;
  for (G4int i=0; i<kNoG3Controls; i++) {
    fControlVector->insert((*right.fControlVector)[i]);
  }   
}

TG4G3ControlVector::~TG4G3ControlVector() {
//
  delete fControlVector;
}

// operators

TG4G3ControlVector& TG4G3ControlVector::operator=(
                                          const TG4G3ControlVector& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // initialize fControlVector 
  fControlVector->clear();
  for (G4int i=0; i<kNoG3Controls; i++) {
    fControlVector->insert((*right.fControlVector)[i]);
  }
  
  return *this;   
}  

G4double TG4G3ControlVector::operator[](G4int index) const
{
//
  if (index < kNoG3Controls)
    return (*fControlVector)[index];
  else {
    TG4Globals::Exception(
      "TG4G3ControlVector::operator[]: index out of the vector scope");
    return 0.;  
  }    
}  

// public methods

void TG4G3ControlVector::SetG3Control(TG4G3Control control, 
                                      G4double controlValue)
{
// Sets the controlValue for the specified process control.
// ---

  if (control<kNoG3Controls) {
    // conversion G4double -> G3ControlValue
    if (abs(controlValue - kUnset) < 0.01) {
        (*fControlVector)[control] = kUnset ;
     }	 
     else if (abs(controlValue - kInActivate) < 0.01) {
        (*fControlVector)[control] = kInActivate; 
     }
     else if (abs(controlValue - kActivate) < 0.01) {
        (*fControlVector)[control] = kActivate; 
     }	
     else if (abs(controlValue - kActivate2) < 0.01) {
        (*fControlVector)[control] = kActivate2; 
     }	     	  
     else {
      G4String text = "TG4G3ControlVector::SetG3Control:\n ";
      text = text + "Inconsistent/Not-yet-implemented control has been ignored.";
      TG4Globals::Warning(text);
     }	
  }
}

void TG4G3ControlVector::SetG3Defaults()
{
// Sets G3 default values for all controls.
// ---

  for (G4int i=0; i<kNoG3Controls; i++) {
   (*fControlVector)[i] = TG4G3Defaults::ControlValue(i);
  } 
}

G4int TG4G3ControlVector::GetControl(G4VProcess* process) const 
{
// Returns the control value for the particle associated with
// the specified process.
// ---

  G4String name = process->GetProcessName();
  if       (name == "conv")    return (*fControlVector)(kPAIR);
  else if  (name == "compt")   return (*fControlVector)(kCOMP);
  else if  (name == "phot")    return (*fControlVector)(kPHOT);
  // else if (name == "??")  return (*fControlVector)(kPFIS);
  else if ((name == "eIoni") || 
           (name == "IeIoni") || 
	   (name == "eIoni+") ||
           (name == "MuIoni") || 
	   (name == "IMuIonisation") ||
	   (name == "hIoni") || 
	   (name == "IhIoni"))    
                               return (*fControlVector)(kDRAY); 
  else if  (name == "annihil") return (*fControlVector)(kANNI);
  else if ((name == "eBrem") || 
           (name == "eBrem+") || 
	   (name == "IeBrems") || 
           (name == "MuBrems") || 
	   (name == "IMuBremsstrahlung"))   
                              return (*fControlVector)(kBREM);
  // else if (name == "??")  return (*fControlVector)(kHADR);
  else if (name == "MuNucl")  return (*fControlVector)(kMUNU);
  else if (name == "Decay")   return (*fControlVector)(kDCAY);
  // else if (name == "??")  return (*fControlVector)(kLOSS);
     // !!! not yet implemented 
  else if ((name == "msc") || 
           (name == "Imsc"))  return (*fControlVector)(kMULS);
  else return kUnset;
}
