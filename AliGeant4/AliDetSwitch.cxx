// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliDetSwitch.h"
#include "AliGlobals.h"

AliDetSwitch::AliDetSwitch(G4String detName, G4int nofVersions, 
                 G4int defaultVersion, G4int pprVersion,
		 AliModuleType modType, G4bool isStandalone)
  : fDetName(detName),
    fNofVersions(nofVersions),
    fDefaultVersion(defaultVersion),
    fPPRVersion(pprVersion),
    fIsStandalone(isStandalone),
    fType(modType),
    fSwitchedVersion(-1)
{
//
}

AliDetSwitch::AliDetSwitch(const AliDetSwitch& right) {
//
  fDetName = right.fDetName;
  fNofVersions = right.fNofVersions;
  fDefaultVersion = right.fDefaultVersion;
  fPPRVersion = right.fPPRVersion;
  fIsStandalone = right.fIsStandalone;
  fType = right.fType;
  fSwitchedVersion = right.fSwitchedVersion;
}

AliDetSwitch::~AliDetSwitch(){
//
}

// operators
AliDetSwitch& AliDetSwitch::operator=(const AliDetSwitch& right)
{    
  // check assignement to self
  if (this == &right) return *this;

  fDetName = right.fDetName;
  fNofVersions = right.fNofVersions;
  fDefaultVersion = right.fDefaultVersion;
  fPPRVersion = right.fPPRVersion;
  fSwitchedVersion = right.fSwitchedVersion;
  fType = right.fType;
  fIsStandalone = right.fIsStandalone;
  
  return *this;
}

G4int AliDetSwitch::operator==(const AliDetSwitch& right) const
{    
//
  G4int returnValue = 0;
  if (fDetName == right.fDetName )
     returnValue = 1;

  return returnValue;  
}

G4int AliDetSwitch::operator!=(const AliDetSwitch& right) const
{ 
//   
  G4int returnValue = 1;
  if (*this == right) returnValue = 0; 
  
  return returnValue;
}
  
// public methods

void AliDetSwitch::SwitchOn(G4int iVersion)
{
// Switchs on the iVersion version.
// ---

  if ((iVersion < 0) || (iVersion >= fNofVersions)) {
    G4String text = "Wrong version number for ";
    text = text + fDetName + ".";
    AliGlobals::Exception(text);
  }  
   
  fSwitchedVersion = iVersion;
}

void AliDetSwitch::SwitchOnDefault()
{
// Switchs on the default version.
// ---

  fSwitchedVersion = fDefaultVersion;
}

void AliDetSwitch::SwitchOnPPR()
{
// Switchs on the default version.
// ---

  fSwitchedVersion = fPPRVersion;
}

void AliDetSwitch::SwitchOff()
{
// No version is switched on.
// ---

  fSwitchedVersion = -1;
}
