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
    static const G4double fgkLength;       //G3 length unit 
    static const G4double fgkTime;         //G3 time unit 
    static const G4double fgkCharge;       //G3 charge unit  
    static const G4double fgkEnergy;       //G3 energy unit  
    static const G4double fgkMass;         //G3 mass unit
    static const G4double fgkMassDensity;  //G3 mass density unit 
    static const G4double fgkAtomicWeight; //G3 atomic weight unit  
    static const G4double fgkField;        //G3 magnetic field unit 
};     

// inline methods

inline G4double TG3Units::Length() { return fgkLength; }
inline G4double TG3Units::Time()   { return fgkTime; }
inline G4double TG3Units::Charge() { return fgkCharge; }
inline G4double TG3Units::Energy() { return fgkEnergy; }
inline G4double TG3Units::Mass()   { return fgkMass; }
inline G4double TG3Units::MassDensity()  { return fgkMassDensity; }
inline G4double TG3Units::AtomicWeight() { return fgkAtomicWeight; }
inline G4double TG3Units::Field()  { return fgkField; }

#endif //TG3_UNITS_H
