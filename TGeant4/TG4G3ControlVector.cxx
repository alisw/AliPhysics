// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4G3ControlVector 
// ------------------------
// See the class description in the header file.

#include "TG4G3ControlVector.h"
#include "TG4G3CutVector.h"
#include "TG4ProcessControlMap.h"
#include "TG4G3Defaults.h"
#include "TG4Globals.h"

#include <G4VProcess.hh>
#include <g4std/strstream>

#include <math.h>

TG4StringVector TG4G3ControlVector::fgControlNameVector;

//_____________________________________________________________________________
TG4G3ControlVector::TG4G3ControlVector()
{
  // initialize fControlVector 
  for (G4int i=0; i<=kNoG3Controls; i++) fControlVector.push_back(kUnset); 
  
  // fill name vector
  if (fgControlNameVector.size() == 0) FillControlNameVector(); 
}

//_____________________________________________________________________________
TG4G3ControlVector::TG4G3ControlVector(const TG4G3ControlVector& right)
  :  fControlVector(right.fControlVector.size())
{
  // copy stuff
  *this = right;  
}

//_____________________________________________________________________________
TG4G3ControlVector::~TG4G3ControlVector() {
//
}

// operators

//_____________________________________________________________________________
TG4G3ControlVector& TG4G3ControlVector::operator=(
                                          const TG4G3ControlVector& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // initialize fControlVector 
  for (G4int i=0; i<=kNoG3Controls; i++) 
    fControlVector[i] = right.fControlVector[i];
  
  return *this;   
}  

//_____________________________________________________________________________
TG4G3ControlValue TG4G3ControlVector::operator[](G4int index) const
{
//
  if (index <= kNoG3Controls)
    return fControlVector[index];
  else {
    TG4Globals::Exception(
      "TG4G3ControlVector::operator[]: index out of the vector scope");
    return kUnset;  
  }    
}  

// private methods

//_____________________________________________________________________________
void TG4G3ControlVector::FillControlNameVector() 
{
// Defines fControlNameVector.
// ---

  fgControlNameVector.push_back("PAIR");
  fgControlNameVector.push_back("COMP");
  fgControlNameVector.push_back("PHOT");
  fgControlNameVector.push_back("PFIS");
  fgControlNameVector.push_back("DRAY");
  fgControlNameVector.push_back("ANNI");
  fgControlNameVector.push_back("BREM");
  fgControlNameVector.push_back("HADR");
  fgControlNameVector.push_back("MUNU");
  fgControlNameVector.push_back("DCAY");
  fgControlNameVector.push_back("LOSS");
  fgControlNameVector.push_back("MULS");
  fgControlNameVector.push_back("CKOV");
  fgControlNameVector.push_back("RAYL");
  fgControlNameVector.push_back("LABS");
  fgControlNameVector.push_back("SYNC");
  fgControlNameVector.push_back("NONE");
}

// public methods

//_____________________________________________________________________________
TG4G3Control TG4G3ControlVector::GetControl(const G4String& controlName)
{
// Retrieves corresponding TG4G3Control constant from the controlName.
// ---

  if      (controlName == fgControlNameVector[kPAIR]) return kPAIR;
  else if (controlName == fgControlNameVector[kCOMP]) return kCOMP;
  else if (controlName == fgControlNameVector[kPHOT]) return kPHOT;
  else if (controlName == fgControlNameVector[kPFIS]) return kPFIS;
  else if (controlName == fgControlNameVector[kDRAY]) return kDRAY;
  else if (controlName == fgControlNameVector[kANNI]) return kANNI;
  else if (controlName == fgControlNameVector[kBREM]) return kBREM;
  else if (controlName == fgControlNameVector[kHADR]) return kHADR;
  else if (controlName == fgControlNameVector[kMUNU]) return kMUNU;
  else if (controlName == fgControlNameVector[kDCAY]) return kDCAY;
  else if (controlName == fgControlNameVector[kLOSS]) return kLOSS;
  else if (controlName == fgControlNameVector[kMULS]) return kMULS;
  else return kNoG3Controls;
}

//_____________________________________________________________________________
const G4String& TG4G3ControlVector::GetControlName(TG4G3Control control)
{
// Returns name of a specified cut.
// ---

  // fill name vector
  if (fgControlNameVector.size() == 0) 
    TG4G3ControlVector::FillControlNameVector(); 

  return fgControlNameVector[control];
}  

//_____________________________________________________________________________
TG4G3ControlValue TG4G3ControlVector::GetControlValue(G4int value,
                                                      TG4G3Control control)
{
// Conversion G4int -> G3ControlValue,
// special treatment for LOSS values 3,4,5.
// ---

  switch (value) {
    case kInActivate: 
      return kInActivate;
      ;;
    case kActivate:
      return kActivate;
      ;;
    case kActivate2:
      return kActivate2;
      ;;
    case 3: case 4: case 5:
      if (control == kLOSS) 
        return kActivate;
      else
        return kUnset;
      ;;	      	  
  }    
  return kUnset;
}    

//_____________________________________________________________________________
TG4G3ControlValue TG4G3ControlVector::GetControlValue(G4double value, 
                                                      TG4G3Control control)
{
// Conversion G4double -> G3ControlValue
// ---

  return TG4G3ControlVector::GetControlValue((G4int)value, control);
}    


//_____________________________________________________________________________
G4bool TG4G3ControlVector::SetControl(TG4G3Control control, 
                                      TG4G3ControlValue controlValue,
				      TG4G3CutVector& cuts)
{
// Sets the controlValue for the specified process control.
// Modifies cuts if necessary.
// Returns true if the control value was set.
// ---

  if (control == kDRAY)
    if (controlValue == kActivate &&
        GetControlValue(kLOSS) == kActivate2) {
      TG4Globals::Warning(
        "TG4Limits::SetG3Control: Cannot set DRAY=1 when LOSS=2.");    
      return false;
    }
    else 
      cuts.SetDeltaRaysOn(true);           

  if (control == kLOSS && controlValue == kActivate2) {
    SetControl(kDRAY, kInActivate, cuts);
    cuts.SetDeltaRaysOn(false);  
  }	

  fControlVector[control] = controlValue;
  return true;
}

//_____________________________________________________________________________
void TG4G3ControlVector::SetG3Defaults()
{
// Sets G3 default values for all controls.
// ---

  for (G4int i=0; i<=kNoG3Controls; i++) 
    fControlVector[i] = TG4G3Defaults::Instance()->ControlValue(i);
}

//_____________________________________________________________________________
G4bool TG4G3ControlVector::Update(const TG4G3ControlVector& vector)
{
// Unset value of DRAY (this information was passed to cut vector.)
// Resets value of LOSS (the special controls process operates only with
// activate/inactivate options.)
// Returns true if some value was modified.
// ---

  G4bool result = false;

  if (fControlVector[kDRAY] != kUnset ) {
      fControlVector[kDRAY] = kUnset;
       result = true;
  }
  
  // if both kLOSS values will have the same effect
  // unset this control

  TG4G3ControlValue passed  = vector[kLOSS];
  TG4G3ControlValue current = fControlVector[kLOSS];

  if (passed  == kActivate2) passed = kActivate;
  if (current == kActivate2) current = kActivate;
           // there is no need to distinguish 
	   // kActivate, kActivate2 after Init phase

  if (current == passed) current = kUnset;
           // if both kLOSS values will have the same effect
           // unset this control

  if (current != fControlVector[kLOSS]) {
     fControlVector[kLOSS] = current;
     result = true;
  }
  return result;     
}

//_____________________________________________________________________________
G4String TG4G3ControlVector::Format() const
{
// Formats the output into a string.
// ---

  strstream tmpStream;

  tmpStream << "  G3 control vector:" << G4endl; 
  for (G4int i=0; i<kNoG3Controls; i++) 
    //if (i != kDRAY) {
      tmpStream << "    " << fgControlNameVector[i] 
                << " control value: " << fControlVector[i] << G4endl; 
    //}	     
    
  return tmpStream.str();  
}	   

//_____________________________________________________________________________
void TG4G3ControlVector::Print() const
{
// Prints the controls.
// ---

  G4cout << Format();	     
}	   

//_____________________________________________________________________________
TG4G3ControlValue 
TG4G3ControlVector::GetControlValue(G4VProcess* process) const 
{
// Returns the control value for the particle associated with
// the specified process.
// ---

  TG4G3Control control 
    = TG4ProcessControlMap::Instance()->GetControl(process);
    
  return fControlVector[control];
}

//_____________________________________________________________________________
TG4G3ControlValue 
TG4G3ControlVector::GetControlValue(TG4G3Control control) const 
{
// Returns the control value for the particle associated with
// the specified process.
// ---

  return fControlVector[control];
}

//_____________________________________________________________________________
G4bool TG4G3ControlVector::IsControl() const
{
// Returns true if any of controls is set.
// ---

  for (G4int i=0; i<kNoG3Controls; i++) 
    if (fControlVector[i] != kUnset) return true;
    
  return false;  
}  
