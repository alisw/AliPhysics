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
      return fgkCUTGAM; break;
    case kCUTELE:  
      return fgkCUTELE; break;
    case kCUTNEU:  
      return fgkCUTNEU; break;
    case kCUTHAD:  
      return fgkCUTHAD; break;
    case kCUTMUO:  
      return fgkCUTMUO; break;
    case kBCUTE:   
      return fgkBCUTE;  break;
    case kBCUTM:   
      return fgkBCUTM;  break; 
    case kDCUTE:   
      return fgkDCUTE;  break;
    case kDCUTM:   
      return fgkDCUTM;  break;
    case kPPCUTM:  
      return fgkPPCUTM; break;
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
      return fgkPAIR; break;
    case kCOMP:
      return fgkCOMP; break;
    case kPHOT:
      return fgkPHOT; break;
    case kPFIS:
      return fgkPFIS; break;
    case kDRAY:
      return fgkDRAY; break;
    case kANNI:
      return fgkANNI; break;
    case kBREM:
      return fgkBREM; break;
    case kHADR:
      return fgkHADR; break;
    case kMUNU:
      return fgkMUNU; break;
    case kDCAY:
      return fgkDCAY; break;
    case kLOSS:
      return fgkLOSS; break;
    case kMULS:
      return fgkMULS; break;
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
