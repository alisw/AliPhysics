// $Id$
// Category: geometry
//
// Data type class that stores available detector option.
// Used in interactive detector setup.

#ifndef ALI_DET_SWITCH_H
#define ALI_DET_SWITCH_H

#include "AliModuleType.h"

#include <globals.hh>

class AliDetSwitch
{
  public:
    AliDetSwitch(G4String detName, G4int nofVersions, G4int defaultVersion,
                 AliModuleType modType = kDetector, G4bool isStandalone = true);
    AliDetSwitch(const AliDetSwitch& right);
    virtual ~AliDetSwitch();

    //operators
    AliDetSwitch& operator=(const AliDetSwitch& right);
    G4int operator==(const AliDetSwitch& right) const;
    G4int operator!=(const AliDetSwitch& right) const;
    
    // methods
    void SwitchOn(G4int version); 
    void SwitchOnDefault(); 
    void SwitchOff(); 

    // get methods
    G4String GetDetName() const;
    G4int GetNofVersions() const;
    G4int GetDefaultVersion() const;
    G4bool IsStandalone() const;
    AliModuleType GetType() const;
    G4int GetSwitchedVersion() const;

  private:
    // data members
    G4String       fDetName;         //module name
    G4int          fNofVersions;     //number of versions
    G4int          fDefaultVersion;  //default version
    G4bool         fIsStandalone;    //true if module can be built standalone
    AliModuleType  fType;            //type of module (detector or structure)
    G4int          fSwitchedVersion; //current selected version
};
    
// inline methods    
    
inline G4String AliDetSwitch::GetDetName() const
{ return fDetName; }

inline G4int AliDetSwitch::GetNofVersions() const
{ return fNofVersions; }

inline G4int AliDetSwitch::GetDefaultVersion() const
{ return fDefaultVersion; }

inline G4int AliDetSwitch::GetSwitchedVersion() const
{ return fSwitchedVersion; }

inline AliModuleType AliDetSwitch::GetType() const
{ return fType; }

inline G4bool AliDetSwitch::IsStandalone() const
{ return fIsStandalone; }

#endif //ALI_DET_SWITCH_H
