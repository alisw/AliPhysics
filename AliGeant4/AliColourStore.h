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
#include <globals.hh> 
#include <g4std/vector>

#include <TColor.h>

class AliColourStore 
{
  typedef G4std::vector<TColor>        ColourVector;
  typedef ColourVector::iterator       ColourIterator;
  typedef ColourVector::const_iterator ColourConstIterator;

  public:
    // --> protected
    // AliColourStore();
    // AliColourStore(const AliColourStore& right);
    virtual ~AliColourStore();
    
    // static methods
    static AliColourStore* Instance();

    // modifiers
    G4Colour AddColour(const G4String& name, 
                       G4double red, G4double blue, G4double green);

    // get methods
    G4Colour GetColour(const G4String& name) const;
    G4Colour GetColour(const TColor& color) const;
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
