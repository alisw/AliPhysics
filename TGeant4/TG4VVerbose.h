// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4VVerbose
// -----------------
// Class defines the verbose level:
// 0 - no output
// 1 - minimal output (default)
// 2 and more - more detailed output

#ifndef TG4_V_VERBOSE_H
#define TG4_V_VERBOSE_H

#include <globals.hh>

class G4UImessenger;

class TG4VVerbose
{
  public:
    TG4VVerbose();
    TG4VVerbose(G4int verboseLevel);      
    virtual ~TG4VVerbose();

    // set methods
    virtual void  VerboseLevel(G4int level);

    // get methods
    virtual G4int VerboseLevel() const;

  private:
    // methods
    virtual G4UImessenger* CreateMessenger() = 0;    

    // static data members
    static const G4int  fgkDefaultVerboseLevel; // default verbose level

    // data members
    G4int  fVerboseLevel; // verbose level
};     

// inline methods

inline void TG4VVerbose::VerboseLevel(G4int level)
{ fVerboseLevel =  level; }

inline G4int TG4VVerbose::VerboseLevel() const
{ return fVerboseLevel; }

#endif //TG4_V_VERBOSE_H
