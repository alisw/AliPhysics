// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4ProcessControlMap
// --------------------------
// Singleton map container for associated pairs
// G4 process name and TG4G3Control.

#ifndef TG4_PROCESS_CONTROL_MAP_H
#define TG4_PROCESS_CONTROL_MAP_H

#include <g4std/map>
#include <globals.hh>

#include "TG4G3Control.h"

class G4VProcess;

class TG4ProcessControlMap
{
  typedef G4std::map<G4String, TG4G3Control, G4std::less<G4String> >  Map;
  typedef Map::iterator       MapIterator;
  typedef Map::const_iterator MapConstIterator;

  public:
    TG4ProcessControlMap();
    // --> protected
    // TG4ProcessControlMap(const TG4ProcessControlMap& right);
    virtual ~TG4ProcessControlMap();

    // static access method
    static TG4ProcessControlMap* Instance();

    // methods
    G4bool Add(G4VProcess* process,  TG4G3Control second);  
    G4bool Add(G4String processName, TG4G3Control second);  
    void PrintAll() const;
    void Clear();

    // get methods
    TG4G3Control    GetControl(const G4VProcess* process);
    TG4G3Control    GetControl(const G4String& processName);
    const G4String& GetControlName(const G4VProcess* process);
    const G4String& GetControlName(const G4String& processName);

  protected:
    TG4ProcessControlMap(const TG4ProcessControlMap& right);

    // operators
    TG4ProcessControlMap& operator=(const TG4ProcessControlMap& right);
  
  private:
    // methods
    G4bool IsDefined(const G4String& processName);

    // static data members
    static TG4ProcessControlMap*  fgInstance; //this instance

    // data members
    Map  fMap; //map container
};

// inline methods

inline TG4ProcessControlMap* TG4ProcessControlMap::Instance() 
{ return fgInstance; }

#endif //TG4_PROCESS_CONTROL_MAP_H
