// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4G3Defaults.h"
#include "TG4Globals.h"

#include <math.h>

// static const data members

// precision tolerance
const G4double TG4G3Defaults::fgCutTolerance = 1. * keV;

// kinetic energy cuts
const G4double TG4G3Defaults::fgCUTGAM = 0.001 * GeV;
const G4double TG4G3Defaults::fgCUTELE = 0.001 * GeV;
const G4double TG4G3Defaults::fgCUTNEU = 0.01 * GeV;
const G4double TG4G3Defaults::fgCUTHAD = 0.01 * GeV;
const G4double TG4G3Defaults::fgCUTMUO = 0.01 * GeV;
const G4double TG4G3Defaults::fgBCUTE  = fgCUTGAM;
const G4double TG4G3Defaults::fgBCUTM  = fgCUTGAM;
const G4double TG4G3Defaults::fgDCUTE  = 10. * TeV;
const G4double TG4G3Defaults::fgDCUTM  = 10. * TeV;
const G4double TG4G3Defaults::fgPPCUTM = 0.01 * GeV;

// physics processes
const TG3FlagValue TG4G3Defaults::fgPAIR = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgCOMP = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgPHOT = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgPFIS = kInActivate; // 0
const TG3FlagValue TG4G3Defaults::fgDRAY = kActivate2;  // 2
const TG3FlagValue TG4G3Defaults::fgANNI = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgBREM = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgHADR = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgMUNU = kInActivate; // 0
const TG3FlagValue TG4G3Defaults::fgDCAY = kActivate;   // 1
const TG3FlagValue TG4G3Defaults::fgLOSS = kActivate2;  // 2
const TG3FlagValue TG4G3Defaults::fgMULS = kActivate;   // 1

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
      return fgCUTGAM; break;
    case kCUTELE:  
      return fgCUTELE; break;
    case kCUTNEU:  
      return fgCUTNEU; break;
    case kCUTHAD:  
      return fgCUTHAD; break;
    case kCUTMUO:  
      return fgCUTMUO; break;
    case kBCUTE:   
      return fgBCUTE;  break;
    case kBCUTM:   
      return fgBCUTM;  break; 
    case kDCUTE:   
      return fgDCUTE;  break;
    case kDCUTM:   
      return fgDCUTM;  break;
    case kPPCUTM:  
      return fgPPCUTM; break;
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
      return fgPAIR; break;
    case kCOMP:
      return fgCOMP; break;
    case kPHOT:
      return fgPHOT; break;
    case kPFIS:
      return fgPFIS; break;
    case kDRAY:
      return fgDRAY; break;
    case kANNI:
      return fgANNI; break;
    case kBREM:
      return fgBREM; break;
    case kHADR:
      return fgHADR; break;
    case kMUNU:
      return fgMUNU; break;
    case kDCAY:
      return fgDCAY; break;
    case kLOSS:
      return fgLOSS; break;
    case kMULS:
      return fgMULS; break;
    default:
      TG4Globals::Warning("TG4G3Defaults::FlagValue: Inconsistent flag.");
      return kUnset;      
  }
}          

G4bool TG4G3Defaults::IsDefaultCut(TG3Cut g3Cut, G4double value)
{
// Tests if the parameter value is equal to the G3 default value.
// ---

  if (abs(value*GeV - CutValue(g3Cut)) > fgCutTolerance) 
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
