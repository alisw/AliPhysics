// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4G3Defaults.h"
#include "TG4Globals.h"

#include <math.h>

// static const data members

// precision tolerance
const G4double TG4G3Defaults::fgkCutTolerance = 1. * keV;

// kinetic energy cuts
const G4double TG4G3Defaults::fgkCUTGAM = 0.001 * GeV;
const G4double TG4G3Defaults::fgkCUTELE = 0.001 * GeV;
const G4double TG4G3Defaults::fgkCUTNEU = 0.01 * GeV;
const G4double TG4G3Defaults::fgkCUTHAD = 0.01 * GeV;
const G4double TG4G3Defaults::fgkCUTMUO = 0.01 * GeV;
const G4double TG4G3Defaults::fgkBCUTE  = fgkCUTGAM;
const G4double TG4G3Defaults::fgkBCUTM  = fgkCUTGAM;
const G4double TG4G3Defaults::fgkDCUTE  = 10. * TeV;
const G4double TG4G3Defaults::fgkDCUTM  = 10. * TeV;
const G4double TG4G3Defaults::fgkPPCUTM = 0.01 * GeV;

// physics processes
const TG3FlagValue TG4G3Defaults::fgkPAIR = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgkCOMP = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgkPHOT = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgkPFIS = kInActivate; // 0
const TG3FlagValue TG4G3Defaults::fgkDRAY = kActivate2;  // 2
const TG3FlagValue TG4G3Defaults::fgkANNI = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgkBREM = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgkHADR = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgkMUNU = kInActivate; // 0
const TG3FlagValue TG4G3Defaults::fgkDCAY = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgkLOSS = kActivate2;  // 2
const TG3FlagValue TG4G3Defaults::fgkMULS = kActivate;   // 1

TG4G3Defaults::TG4G3Defaults() {
//
}
  
TG4G3Defaults::~TG4G3Defaults() {
//
}

G4double TG4G3Defaults::CutValue(G4int g3Cut)
{
// Returns the G3 default value for the specified cut.
// ---

  switch (g3Cut) {
    case kCUTGAM:  
      return fgkCUTGAM;
    case kCUTELE:  
      return fgkCUTELE;
    case kCUTNEU:  
      return fgkCUTNEU;
    case kCUTHAD:  
      return fgkCUTHAD;
    case kCUTMUO:  
      return fgkCUTMUO;
    case kBCUTE:   
      return fgkBCUTE; 
    case kBCUTM:   
      return fgkBCUTM; 
    case kDCUTE:   
      return fgkDCUTE; 
    case kDCUTM:   
      return fgkDCUTM; 
    case kPPCUTM:  
      return fgkPPCUTM;
    default:
      TG4Globals::Warning("TG4G3Defaults::CutValue: Inconsistent cut.");
      return 0.;      
  }
}          

TG3FlagValue TG4G3Defaults::FlagValue(G4int g3Flag)
{
// Returns the G3 default value for the specified flag.
// ---

  switch (g3Flag) {
    case kPAIR:
      return fgkPAIR;
    case kCOMP:
      return fgkCOMP;
    case kPHOT:
      return fgkPHOT;
    case kPFIS:
      return fgkPFIS;
    case kDRAY:
      return fgkDRAY;
    case kANNI:
      return fgkANNI;
    case kBREM:
      return fgkBREM;
    case kHADR:
      return fgkHADR;
    case kMUNU:
      return fgkMUNU;
    case kDCAY:
      return fgkDCAY;
    case kLOSS:
      return fgkLOSS;
    case kMULS:
      return fgkMULS;
    default:
      TG4Globals::Warning("TG4G3Defaults::FlagValue: Inconsistent flag.");
      return kUnset;      
  }
}          

G4bool TG4G3Defaults::IsDefaultCut(TG3Cut g3Cut, G4double value)
{
// Tests if the parameter value is equal to the G3 default value.
// ---

  if (abs(value*GeV - CutValue(g3Cut)) > fgkCutTolerance) 
    return false;
  else  
    return true;
}

G4bool TG4G3Defaults::IsDefaultFlag(TG3Flag g3Flag, G4double value)
{
// Tests if the parameter value is equal to the G3 default value.
// ---

  if (value == FlagValue(g3Flag)) 
    return true;
  else  
    return false;
}
