// $Id$
// Category: geometry
//
// This class adds integer identifier data member to G4VSensitiveDetector

#ifndef TG4V_SENSITIVE_DETECTOR_H
#define TG4V_SENSITIVE_DETECTOR_H

#include <G4VSensitiveDetector.hh>
#include <globals.hh>

class TG4VSensitiveDetector : public G4VSensitiveDetector
{
  public:
    TG4VSensitiveDetector(G4String sdName);
    TG4VSensitiveDetector(G4String sdName, G4int id);
    TG4VSensitiveDetector(const TG4VSensitiveDetector& right);
    // --> protected
    // TG4VSensitiveDetector();
    virtual ~TG4VSensitiveDetector();

    // operators
    TG4VSensitiveDetector& operator=(const TG4VSensitiveDetector &right);

    // static get method
    static G4int GetTotalNofSensitiveDetectors();

    // get methods
    G4int GetID() const;
    
  protected:
    TG4VSensitiveDetector(); 

    // data members
    G4int  fID;                //sensitive detector ID
    
  private:          
    // data members
    static G4int  fgSDCounter; //sensitive detector counter
};

// inline methods

inline G4int TG4VSensitiveDetector::GetTotalNofSensitiveDetectors()
{ return fgSDCounter; }

inline G4int TG4VSensitiveDetector::GetID() const
{ return fID; }

#endif //TG4V_SENSITIVE_DETECTOR_H


