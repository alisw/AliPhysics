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
class AliMCQA;

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
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual void UserProcessHits(const G4Track* track, const G4Step* step);
    
  protected:  
    AliSensitiveDetector();
    
  private:
    // data members    
    AliModule*       fModule;      //AliModule
    G4int            fModuleID;    //AliModule index in AliRun::fModules
    AliMCQA*         fMCQA;        //AliMCQA    
};

#endif //ALI_SENSITIVE_DETECTOR_H


