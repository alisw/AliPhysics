// $Id$
// Category: geometry
//
// Special sensitive detector class for lego run.
// It implements G4VSensitiveDetector::ProcessHits() 
// with AliLego:: StepManager().

#ifndef ALI_LEGO_SENSITIVE_DETECTOR_H
#define ALI_LEGO_SENSITIVE_DETECTOR_H

#include "TG4VSensitiveDetector.h"

#include <globals.hh>

class AliLego;
class TG4StepManager;

class G4HCofThisEvent;
class G4Step;

class AliLegoSensitiveDetector : public TG4VSensitiveDetector
{
  public:
    AliLegoSensitiveDetector(G4String name, AliLego* lego,
                             G4VSensitiveDetector* standardSD);
    AliLegoSensitiveDetector(const AliLegoSensitiveDetector& right);
    // --> protected
    // AliLegoSensitiveDetector();
    virtual ~AliLegoSensitiveDetector();

    // operators
    AliLegoSensitiveDetector& operator=(const AliLegoSensitiveDetector& right);

    // methods
    virtual void Initialize(G4HCofThisEvent* hc);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void EndOfEvent(G4HCofThisEvent* hce);
    //virtual void clear();
    virtual void PrintAll();
    virtual void DrawAll();
    
    // get methods
    G4VSensitiveDetector* GetStandardSD() const;
    
  protected:  
    AliLegoSensitiveDetector();
    
  private:
    // data members
    AliLego*               fLego;        //lego from AliRoot
    G4VSensitiveDetector*  fStandardSD;  //standard sensitive detector
    TG4StepManager*        fStepManager; //TG4StepManager
};

// inline methods

inline G4VSensitiveDetector* AliLegoSensitiveDetector::GetStandardSD() const
{ return fStandardSD; }

#endif //ALI_LEGO_SENSITIVE_DETECTOR_H


