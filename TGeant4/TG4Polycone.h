// $Id$
// Category: geometry
// by I. Hrivnacova, 12.10.2000 
//
// G4Polycone class extented for public access method
// to the original parameters.

#ifndef TG4_POLYCONE_H
#define TG4_POLYCONE_H

#include <G4Polycone.hh>
#include <globals.hh>

class TG4Polycone : public G4Polycone
{
  public:
    TG4Polycone(const G4Polycone& rhs);
    // TG4Polycone(); --> moved to private
    virtual ~TG4Polycone();

    // methods 
   G4int GetNofZPlanes();
   G4double* GetRmin();
   G4double* GetRmax();
   G4double* GetZ();

  private:
    TG4Polycone();
    
    void CheckOrigin();
};

#endif //TG4_POLYCONE_H

