// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4G3CutVector
// --------------------
// Vector of kinetic energy cut values with
// convenient set/get methods.

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
    
    // static methods
    static G4double Tolerance();
    static TG4G3Cut GetCut(const G4String& cutName);
    static const G4String& GetCutName(TG4G3Cut cut);
    
    // set methods
    void SetCut(TG4G3Cut cut, G4double cutValue);
    void SetG3Defaults();
    G4bool Update(const TG4G3CutVector& vector);
    void SetDeltaRaysOn(G4bool value);
    
    // methods
    void Print() const;

    // get methods
    G4double GetMinEkineForGamma(const G4Track& track) const;
    G4double GetMinEkineForElectron(const G4Track& track) const;
    G4double GetMinEkineForEplus(const G4Track& track) const;
    G4double GetMinEkineForChargedHadron(const G4Track& track) const;
    G4double GetMinEkineForNeutralHadron(const G4Track& track) const;
    G4double GetMinEkineForMuon(const G4Track& track) const;
    G4double GetMinEkineForOther(const G4Track& track) const;
    G4bool IsCut() const;

  private:
    // static methods 
    static void FillCutNameVector();
  
    // static data members
    static const G4double   fgkDCUTEOff;  //cut for delta rays by e- (if off)
    static const G4double   fgkDCUTMOff;  //cut for delta rays by mu (if off)
    static const G4double   fgkTolerance; //tolerance for comparing cuts
    static TG4StringVector  fgCutNameVector; //vector of cut parameters names
  
    // data members
    TG4doubleVector  fCutVector;  //vector of kinetic energy cut values 
    G4bool           fDeltaRaysOn;//delta rays process control
};

// inline methods

inline G4double TG4G3CutVector::Tolerance()
{ return fgkTolerance; }

inline void TG4G3CutVector::SetDeltaRaysOn(G4bool value)
{ fDeltaRaysOn = value; }

#endif //TG4_CUT_VECTOR_H



