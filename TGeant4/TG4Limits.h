// $Id$
// Category: global
//
// G4UserLimits derived class extended with
// vectors of kinetic energy cuts and control process flags
// data members

#ifndef TG4_LIMITS_H
#define TG4_LIMITS_H

#include "TG4Globals.h"
#include "TG3Cut.h"
#include "TG3Flag.h"

#include <G4UserLimits.hh>

class TG4CutVector;
class TG4FlagVector;

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
    void SetG3Cut(TG3Cut g3Cut, G4double cutValue);
    void SetG3Flag(TG3Flag g3Flag, G4double flagValue);
    void SetG3DefaultCuts();
    void SetG3DefaultFlags();
    
    // get methods
    G4bool IsCut() const;
    G4bool IsFlag() const;
    virtual G4double GetUserMinEkine(const G4Track& track);
    G4double GetMinEkineForGamma(const G4Track& track) const;
    G4double GetMinEkineForElectron(const G4Track& track) const;
    G4double GetMinEkineForHadron(const G4Track& track) const;
    G4double GetMinEkineForNeutralHadron(const G4Track& track) const;
    G4double GetMinEkineForMuon(const G4Track& track) const;
    G4double GetMinEkineForOther(const G4Track& track) const;
    G4int GetFlag(G4VProcess* process) const; 

  private:
    // data members
    G4bool            fIsCut;      //true if any cut value is set
    G4bool            fIsFlag;     //true if any flag value is set
    TG4CutVector*     fCutVector;  //TG4CutVector
    TG4FlagVector*    fFlagVector; //TG4FlagVector
};

// inline methods

inline G4bool TG4Limits::IsCut() const  { return fIsCut; }
inline G4bool TG4Limits::IsFlag() const { return fIsFlag; }

#endif //TG4_USER_LIMITS_H



