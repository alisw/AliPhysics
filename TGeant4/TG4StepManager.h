// $Id$
// Category: event
//
// Geant4 implementation of the MonteCarlo interface methods                    
// for access to Geant4 at step level
//
// The public methods that do not implement AliMC methods
// are commented as G4 specific

#ifndef TG4_STEP_MANAGER_H
#define TG4_STEP_MANAGER_H

#include "TG4StepStatus.h"

#include <G4Step.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

#include <Rtypes.h>

class G4Track;
class G4SteppingManager;
class G4VPhysicalVolume;

class TLorentzVector;

class TG4StepManager
{
  public:
    TG4StepManager();
    // --> protected
    // TG4StepManager(const TG4StepManager& right);
    virtual ~TG4StepManager();

    // static access method
    static TG4StepManager* Instance();
        
    // methods
    void StopTrack(); //new
    void StopEvent(); //new
    void Rndm(Float_t* array, const Int_t size) const;
    
    // set methods
    void SetStep(G4Step* step, TG4StepStatus status);    // G4 specific
    void SetStep(G4Track* track, TG4StepStatus status);  // G4 specific
    void SetSteppingManager(G4SteppingManager* manager); // G4 specific
    void SetMaxStep(Float_t step);
    void SetMaxNStep(Int_t maxNofSteps);  //??
    void SetUserDecay(Int_t pdg);  //NEW
    
    // get methods
    G4Track* GetTrack() const;                            // G4 specific
    G4Step*  GetStep() const;                             // G4 specific
    TG4StepStatus GetStepStatus() const;                  // G4 specific
        
        // tracking volume(s) 
    G4VPhysicalVolume* GetCurrentPhysicalVolume() const;  // G4 specific
    Int_t CurrentVolID(Int_t& copyNo) const;
    Int_t CurrentVolOffID(Int_t off, Int_t& copyNo) const;
    const char* CurrentVolName() const;
    const char* CurrentVolOffName(Int_t off) const;
    Int_t CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens, 
                    Float_t &radl, Float_t &absl) const;
    void Gmtod(Float_t* xm, Float_t* xd, Int_t iflag);  //new
    void Gdtom(Float_t* xd, Float_t* xm, Int_t iflag); //new
    Float_t MaxStep() const;
    Int_t GetMaxNStep() const;  //??                       
    Int_t GetMedium() const;  //??

        // tracking particle 
        // dynamic properties
    void TrackPosition(TLorentzVector& position) const;
    void TrackMomentum(TLorentzVector& momentum) const;
    void TrackVertexPosition(TLorentzVector& position) const;
    void TrackVertexMomentum(TLorentzVector& momentum) const;
    Float_t TrackStep() const;  
    Float_t TrackLength() const;   
    Float_t TrackTime() const;  
    Float_t Edep() const;
        // static properties
    Int_t TrackPid() const;
    Float_t TrackCharge() const;
    Float_t TrackMass() const;
    Float_t Etot() const;

        // track status
    Bool_t IsTrackInside() const;
    Bool_t IsTrackEntering() const;
    Bool_t IsTrackExiting() const;
    Bool_t IsTrackOut() const;
    Bool_t IsTrackDisappeared() const;
    Bool_t IsTrackStop() const;
    Bool_t IsTrackAlive() const;
    Bool_t IsNewTrack() const;

        // secondaries
    Int_t NSecondaries() const;
    void GetSecondary(Int_t isec, Int_t& particleId,
                      TLorentzVector& position, TLorentzVector& momentum);      
    const char* ProdProcess() const; 

  protected:
    TG4StepManager(const TG4StepManager& right);

    // operators
    TG4StepManager& operator=(const TG4StepManager& right);

  private:
    // methods
    void CheckTrack() const;
    void CheckStep() const;
    void CheckSteppingManager() const;
    void SetTLorentzVector(G4ThreeVector xyz, G4double t, 
                           TLorentzVector& lv) const;    
    G4VPhysicalVolume* GetCurrentOffPhysicalVolume(G4int off) const;
    G4int GetVolumeID(G4VPhysicalVolume* volume) const;

    // static data members
    static TG4StepManager*  fgInstance;   //this instance
    
    // data members
    G4Track*            fTrack;           //current track
    G4Step*             fStep;            //current step
    TG4StepStatus       fStepStatus;      //step status that decides whether
                                          //track properties will be returned
					  //from PreStepPoint or PostStepPoint
    G4SteppingManager*  fSteppingManager; //G4SteppingManager
};

// inline methods

inline TG4StepManager* TG4StepManager::Instance() 
{ return fgInstance; }

inline void TG4StepManager::SetStep(G4Step* step, TG4StepStatus status)
{ fTrack = step->GetTrack(); fStep = step; fStepStatus = status; }

inline void TG4StepManager::SetStep(G4Track* track, TG4StepStatus status)
{ fTrack = track; fStep = 0; fStepStatus = status; }

inline void TG4StepManager::SetSteppingManager(G4SteppingManager* manager)
{ fSteppingManager = manager; }

inline G4Track* TG4StepManager::GetTrack() const
{ return fTrack; }

inline G4Step* TG4StepManager::GetStep() const
{ return fStep; }

inline TG4StepStatus TG4StepManager::GetStepStatus() const
{ return fStepStatus; }

#endif //TG4_STEP_MANAGER_H

