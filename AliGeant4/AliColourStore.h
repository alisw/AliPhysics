// $Id$
// Category: visualization
//
// Author: I. Hrivnacova
//
// Class AliColourStore
// --------------------
// Singleton data type class - store for the predefined colours.

#ifndef ALI_COLOUR_STORE_H
#define ALI_COLOUR_STORE_H

#include <G4Colour.hh>
 
class AliColour;

#include <g4std/vector>

class AliColourStore 
{
  typedef G4std::vector<AliColour>     ColourVector;
  typedef ColourVector::iterator       ColourIterator;
  typedef ColourVector::const_iterator ColourConstIterator;

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
    ColourVector  fColours; //vector of AliColour
};   

#endif //ALI_COLOUR_STORE_H
