// $Id$
// Category: global
//
// Class defines the G3 default units of physical quantities; 
// all physical quantities returned by MC are expressed in these units.

#ifndef TG3_UNITS_H
#define TG3_UNITS_H

#include <globals.hh>

class TG3Units
{
  public:
    // --> protected
    // TG3Units();  
    virtual ~TG3Units();

    // static get methods
    static G4double Length(); 
    static G4double Time(); 
    static G4double Charge(); 
    static G4double Energy(); 
    static G4double Mass(); 
    static G4double MassDensity(); 
    static G4double AtomicWeight();     
    static G4double Field(); 
      
  protected:
    TG3Units();      
        // only static data members and methods

  private:
    // static data members  
    static const G4double fgLength;       //G3 length unit 
    static const G4double fgTime;         //G3 time unit 
    static const G4double fgCharge;       //G3 charge unit  
    static const G4double fgEnergy;       //G3 energy unit  
    static const G4double fgMass;         //G3 mass unit
    static const G4double fgMassDensity;  //G3 mass density unit 
    static const G4double fgAtomicWeight; //G3 atomic weight unit  
    static const G4double fgField;        //G3 magnetic field unit 
};     

// inline methods

inline G4double TG3Units::Length() { return fgLength; }
inline G4double TG3Units::Time()   { return fgTime; }
inline G4double TG3Units::Charge() { return fgCharge; }
inline G4double TG3Units::Energy() { return fgEnergy; }
inline G4double TG3Units::Mass()   { return fgMass; }
inline G4double TG3Units::MassDensity()  { return fgMassDensity; }
inline G4double TG3Units::AtomicWeight() { return fgAtomicWeight; }
inline G4double TG3Units::Field()  { return fgField; }

#endif //TG3_UNITS_H
