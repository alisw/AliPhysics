// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4IntMap
// ---------------
// The map container for integer numbers associated with names.

#ifndef TG4_INT_MAP_H
#define TG4_INT_MAP_H

#include <g4std/map>
#include <globals.hh>

class TG4IntMap
{
  typedef G4std::map<G4String, G4int, G4std::less<G4String> >  Map;
  typedef Map:: iterator       MapIterator;
  typedef Map:: const_iterator MapConstIterator;

  public:
    TG4IntMap();
    // --> protected
    // TG4IntMap(const TG4IntMap& right);
    virtual ~TG4IntMap();

    // methods
    G4bool Add(const G4String& first, G4int second);  
    G4int GetSecond(const G4String& name, G4bool warn = true);
    void PrintAll() const;
    void Clear();

  protected:
    TG4IntMap(const TG4IntMap& right);

    // operators
    TG4IntMap& operator=(const TG4IntMap& right);
  
  private:
    G4bool IsDefined(const G4String& first);
  
    // data members
    Map  fMap; //map container
};

#endif //TG4_NAME_MAP_H
