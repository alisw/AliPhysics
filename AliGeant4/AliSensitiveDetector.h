// $Id$
// Category: geometry
//
// Sensitive detector class,
// that implements G4VSensitiveDetector::ProcessHits()
// with AliDetector:: StepManager().

#ifndef ALI_SENSITIVE_DETECTOR_H
#define ALI_SENSITIVE_DETECTOR_H

#include "TG4VSensitiveDetector.h"

#include <globals.hh>

class AliModule;
class TG4StepManager;

class G4HCofThisEvent;
class G4Step;

class AliSensitiveDetector : public TG4VSensitiveDetector
{
  public:
    AliSensitiveDetector(G4String sdName, AliModule* module);
    AliSensitiveDetector(G4String sdName, AliModule* module, G4int id);
    AliSensitiveDetector(const AliSensitiveDetector& right);
    // --> protected
    // AliSensitiveDetector();
    virtual ~AliSensitiveDetector();

    // operators
    AliSensitiveDetector& operator=(const AliSensitiveDetector& right);
  
    // methods
    virtual void Initialize(G4HCofThisEvent* hc);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void EndOfEvent(G4HCofThisEvent* hce);
    //virtual void clear();
    virtual void PrintAll();
    virtual void DrawAll();
    
  protected:  
    AliSensitiveDetector();
    
  private:
    // data members
    AliModule*       fModule;      //AliModule
    TG4StepManager*  fStepManager; //TG4StepManager
};

#endif //ALI_SENSITIVE_DETECTOR_H


