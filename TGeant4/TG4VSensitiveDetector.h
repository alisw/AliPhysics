// $Id$
// Category: geometry
//
// This class adds integer identifier data member to G4VSensitiveDetector.
// It also takes care of setting step status (kBoundary, kNormalStep)
// and passing G4Step to TG4StepManager before calling user derived
// sensitive detector class. 

#ifndef TG4V_SENSITIVE_DETECTOR_H
#define TG4V_SENSITIVE_DETECTOR_H

#include <G4VSensitiveDetector.hh>
#include <globals.hh>

class TG4StepManager;

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

    // methods
    virtual void UserProcessHits(const G4Track* track, const G4Step* step) = 0;
                   // the following methods should not
		   // be overwritten in a derived class
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual G4bool ProcessHitsOnBoundary(G4Step* step);
 
    // static get method
    static G4int GetTotalNofSensitiveDetectors();

    // get methods
    G4int GetID() const;
    
  protected:
    TG4VSensitiveDetector(); 

    // data members
    G4int  fID;                    //sensitive detector ID
    TG4StepManager*  fStepManager; //TG4StepManager
    
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


