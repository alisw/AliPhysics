// $Id$
// Category: global
//
// Vector of kinetic energy cut values with
// convenient set/get methods

#ifndef TG4_G3_CUT_VECTOR_H
#define TG4_G3_CUT_VECTOR_H

#include "TG4Globals.h"
#include "TG4G3Cut.h"

class G4Track;

class TG4G3CutVector
{
  public:
    TG4G3CutVector();
    TG4G3CutVector(const TG4G3CutVector& right);
    virtual ~TG4G3CutVector();
    
    // operators
    TG4G3CutVector& operator=(const TG4G3CutVector& right);
    G4double operator[](G4int index) const;
        
    // set methods
    void SetG3Cut(TG4G3Cut cut, G4double cutValue);
    void SetG3Defaults();
    
    // get methods
    G4double GetMinEkine(const G4Track& track) const;
    G4double GetMinEkineForGamma(const G4Track& track) const;
    G4double GetMinEkineForElectron(const G4Track& track) const;
    G4double GetMinEkineForHadron(const G4Track& track) const;
    G4double GetMinEkineForNeutralHadron(const G4Track& track) const;
    G4double GetMinEkineForMuon(const G4Track& track) const;
    G4double GetMinEkineForOther(const G4Track& track) const;

  private:
    // data members
    TG4doubleVector*  fCutVector; //vector of kinetic energy cut values 
};

#endif //TG4_CUT_VECTOR_H



