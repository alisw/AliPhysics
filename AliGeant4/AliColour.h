// $Id$
// Category: visualization
//
// Data type class that defines colours with names.

#ifndef ALI_COLOUR_H
#define ALI_COLOUR_H

#include <G4Colour.hh>
#include <globals.hh>

class AliColour
{
  public:
    AliColour();
    AliColour(G4String name, G4double red, G4double blue, G4double green);
    AliColour(const AliColour& right);
    virtual ~AliColour();

    // operators
    AliColour& operator=(const AliColour& right);
    G4int operator==(const AliColour& right) const;
    G4int operator!=(const AliColour& right) const;

    // get methods
    G4Colour GetColour() const;
    G4String GetName() const;
  
  private:
    G4String  fName;  //colour name
    G4double  fRed;   //red component
    G4double  fBlue;  //blue component
    G4double  fGreen; //green component
};

// inline methods

inline G4Colour AliColour::GetColour() const
{ return G4Colour(fRed, fBlue, fGreen); }

inline G4String AliColour::GetName() const
{ return fName; }

#endif //ALCOLOUR_H

