// $Id$
// Category: global
//
// G4UserLimits derived class extended with
// vectors of kinetic energy cuts and control process flags
// data members

#ifndef TG4_LIMITS_H
#define TG4_LIMITS_H

#include "TG4Globals.h"
#include "TG4G3Cut.h"
#include "TG4G3Control.h"

#include <G4UserLimits.hh>

class TG4G3CutVector;
class TG4G3ControlVector;

class G4VProcess;

class TG4Limits: public G4UserLimits
{
  public:
    TG4Limits();
    TG4Limits(const TG4Limits& right);
    virtual ~TG4Limits();
    
    // operators
    TG4Limits& operator=(const TG4Limits& right);

    // set methods
    void SetG3Cut(TG4G3Cut cut, G4double cutValue);
    void SetG3Control(TG4G3Control control, G4double controlValue);
    void SetG3DefaultCuts();
    void SetG3DefaultControls();
    
    // get methods
    G4bool IsCut() const;
    G4bool IsControl() const;
    virtual G4double GetUserMinEkine(const G4Track& track);
    G4double GetMinEkineForGamma(const G4Track& track) const;
    G4double GetMinEkineForElectron(const G4Track& track) const;
    G4double GetMinEkineForHadron(const G4Track& track) const;
    G4double GetMinEkineForNeutralHadron(const G4Track& track) const;
    G4double GetMinEkineForMuon(const G4Track& track) const;
    G4double GetMinEkineForOther(const G4Track& track) const;
    G4int GetControl(G4VProcess* process) const; 

  private:
    // data members
    G4bool              fIsCut;    //true if any cut value is set
    G4bool              fIsControl;//true if any flag value is set
    TG4G3CutVector*     fCutVector;    //TG4CutVector
    TG4G3ControlVector* fControlVector;//TG4ControlVector
};

// inline methods

inline G4bool TG4Limits::IsCut() const  { return fIsCut; }
inline G4bool TG4Limits::IsControl() const { return fIsControl; }

#endif //TG4_USER_LIMITS_H



