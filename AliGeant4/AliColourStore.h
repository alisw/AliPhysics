// $Id$
// Category: visualization
//
// Singleton data type class - store for the predefined colours.

#ifndef ALI_COLOUR_STORE_H
#define ALI_COLOUR_STORE_H

#include "AliColour.h"

#include <g4rw/tvordvec.h>

class AliColourStore 
{
  typedef G4RWTValOrderedVector<AliColour>  AliColourVector;

  public:
    // --> protected
    // AliColourStore();
    // AliColourStore(const AliColourStore& right);
    virtual ~AliColourStore();
    
    // static methods
    static AliColourStore* Instance();

    // get methods
    G4Colour GetColour(G4String name) const;
    G4String GetColoursList() const;
    G4String GetColoursListWithCommas() const;
    
  protected:
    AliColourStore();  
    AliColourStore(const AliColourStore& right);

    // operators
    AliColourStore& operator=(const AliColourStore& right);

  private:
    // static data members
    static AliColourStore*  fgInstance; //this instance

    // data members
    AliColourVector  fColours;          //vector of AliColour
};   

#endif //ALI_COLOUR_STORE_H
