// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliDetSwitchVector
// ---------------------------
// The class contains a vector of detector switches
// and provides methods for their interactive setting.

#ifndef ALI_DET_SWITCH_VECTOR_H
#define ALI_DET_SWITCH_VECTOR_H

#include "AliDetSwitchVectorMessenger.h"

#include <globals.hh>
#include <g4std/vector>

class AliDetSwitch;

class AliDetSwitchVector
{
  typedef G4std::vector<AliDetSwitch*>    DetSwitchVector;
  typedef DetSwitchVector::iterator       DetSwitchIterator;
  typedef DetSwitchVector::const_iterator DetSwitchConstIterator;

  public:
    AliDetSwitchVector();
    // --> protected
    // AliDetSwitchVector(const AliDetSwitchVector& right);
    virtual ~AliDetSwitchVector();

    // methods
    void Add(AliDetSwitch* detSwitch);
    void UpdateMessenger();
    void SwitchDetOn(const G4String& moduleNameVer);
    void SwitchDetOn(const G4String& moduleName, G4int version);
    void SwitchDetOnDefault(const G4String& moduleName);
    void SwitchDetOff(const G4String& moduleName);
    void PrintSwitchedDets() const;
    void PrintAvailableDets() const;

    // get methods
    G4int GetSize() const;
    AliDetSwitch* GetDetSwitch(G4int i) const;
    AliDetSwitch* GetDetSwitch(const G4String& moduleName) const;
    G4String GetSwitchedDetsList() const;
    G4String GetAvailableDetsList() const;
    G4String GetAvailableDetsListWithCommas() const;
    G4String GetDetNamesList() const;
    G4String GetDetNamesListWithCommas() const;
    
  protected:
    AliDetSwitchVector(const AliDetSwitchVector& right);

    // operators
    AliDetSwitchVector& operator=(const AliDetSwitchVector& right);
    
  private:    
    // data members
    AliDetSwitchVectorMessenger  fMessenger;       //messenger
    DetSwitchVector              fDetSwitchVector; //vector of AliDetSwitch
};

// inline methods

inline void AliDetSwitchVector::UpdateMessenger()
{ fMessenger.Update(); }

inline G4int AliDetSwitchVector::GetSize() const
{ return fDetSwitchVector.size(); }

inline AliDetSwitch* AliDetSwitchVector::GetDetSwitch(G4int i) const
{ return fDetSwitchVector[i]; }

#endif //ALI_DET_SWITCH_VECTOR_H

