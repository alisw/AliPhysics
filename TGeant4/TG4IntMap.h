// $Id$
// Category: global
//
// The map container for integer numbers associated with names.

#ifndef TG4_INT_MAP_H
#define TG4_INT_MAP_H

#include "g4std/map"
#include "globals.hh"

class TG4IntMap
{
  public:
    TG4IntMap();
    // --> protected
    // TG4IntMap(const TG4IntMap& right);
    virtual ~TG4IntMap();

    // methods
    G4bool Add(const G4String& first, G4int second);  
    G4int GetSecond(const G4String& name);
    void PrintAll();
    void Clear();

  protected:
    TG4IntMap(const TG4IntMap& right);

    // operators
    TG4IntMap& operator=(const TG4IntMap& right);
  
  private:
    G4bool IsDefined(const G4String& first);
  
    // data members
    G4std::map<G4String, G4int, G4std::less<G4String> > fMap;
                                  //map container
};

#endif //TG4_NAME_MAP_H
