// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsManager
// -----------------------
// Geant4 implementation of the MonteCarlo interface methods                    
// for building Geant4 physics and access to it.

#ifndef TG4_PHYSICS_MANAGER_H
#define TG4_PHYSICS_MANAGER_H

#include "TG4Verbose.h"
#include "TG4PhysicsMessenger.h"
#include "TG4ProcessControlMap.h"
#include "TG4ProcessMCMap.h"
#include "TG4NameMap.h"
#include "TG4G3Cut.h"
#include "TG4G3Control.h"
#include "TG4Globals.h"

#include <Rtypes.h>
#include "AliMCProcess.h"

#include <globals.hh>

class AliDecayer;
class TG4ParticlesManager;
class TG4G3PhysicsManager;
class TG4G3ProcessMap;

class G4ParticleDefinition;
class G4VProcess;
class TG4ModularPhysicsList;

class TG4PhysicsManager : public TG4Verbose
{
  public:
    TG4PhysicsManager(TG4ModularPhysicsList* physicsList);
    // --> protected
    // TG4PhysicsManager();
    // TG4PhysicsManager(const TG4PhysicsManager& right);
    virtual ~TG4PhysicsManager();

    // static access method
    static TG4PhysicsManager* Instance();
        
    // methods
    void Gstpar(Int_t itmed, const char *param, Float_t parval); 

    // set methods
    void SetCut(const char* cutName, Float_t cutValue);
    void SetProcess(const char* controlName, Int_t controlValue);
    Float_t Xsec(char* reac, Float_t energy, Int_t part, Int_t mate);
    void SetExternalDecayer(AliDecayer* decayer);
     
    // get methods
    AliDecayer* Decayer() const;

        // particle table usage         
    Int_t IdFromPDG(Int_t pdgID) const;
    Int_t PDGFromId(Int_t mcID) const;
    void  DefineParticles();      
    
    //
    // methods for Geant4 only 
    //

    void CreatePhysicsConstructors();
    void SetProcessActivation();  
    AliMCProcess GetMCProcess(const G4VProcess* process);
    AliMCProcess GetOpBoundaryStatus(const G4VProcess* process);

    // set methods
    void SetPhysicsList(TG4ModularPhysicsList* physicsList);
    void SetEMPhysics(G4bool value);
    void SetMuonPhysics(G4bool value);
    void SetHadronPhysics(G4bool value);
    void SetOpticalPhysics(G4bool value);
    void SetSpecialCutsPhysics(G4bool value);
    void SetSpecialControlsPhysics(G4bool value);
    
    // get methods
    TG4ModularPhysicsList* GetPhysicsList() const; 
   
  protected:
    TG4PhysicsManager();
    TG4PhysicsManager(const TG4PhysicsManager& right);

    // operators
    TG4PhysicsManager& operator=(const TG4PhysicsManager& right);

  private:
    // methods
    void FillProcessMap();
    void GstparCut(G4int itmed, TG4G3Cut par, G4double parval);
    void GstparControl(G4int itmed, TG4G3Control control, 
                       TG4G3ControlValue parval);

    // static data members
    static TG4PhysicsManager*  fgInstance; //this instance
    
    // data members
    TG4PhysicsMessenger    fMessenger;        //messenger
    TG4ParticlesManager*   fParticlesManager; //particles manager
    TG4G3PhysicsManager*   fG3PhysicsManager; //G3 physics manager
    TG4ModularPhysicsList* fPhysicsList; //physics list
    AliDecayer*            fDecayer;     //external decayer
    TG4ProcessMCMap        fProcessMCMap;//the mapping between G4 process names
                                         //and AliMCProcess codes
    TG4ProcessControlMap   fProcessControlMap; //the mapping between G4 processes
                                         //and G3 process controls
    G4bool  fSetEMPhysics;          //electromagnetic physics control
    G4bool  fSetMuonPhysics;        //muon physics control
    G4bool  fSetHadronPhysics;      //hadron physics control
    G4bool  fSetOpticalPhysics;     //optical physics control
    G4bool  fSetSpecialCutsPhysics; //special cuts process control 
                                    //(under development)                  
    G4bool  fSetSpecialControlsPhysics;//special controls process control
                                    //(under development)
};

// inline methods

inline TG4PhysicsManager* TG4PhysicsManager::Instance() 
{ return fgInstance; }

inline void TG4PhysicsManager::SetExternalDecayer(AliDecayer* decayer) 
{ fDecayer = decayer; }

inline AliDecayer* TG4PhysicsManager::Decayer() const
{ return fDecayer; }

inline void TG4PhysicsManager::SetPhysicsList(TG4ModularPhysicsList* physicsList)
{ fPhysicsList = physicsList; }

inline void TG4PhysicsManager::SetEMPhysics(G4bool value)
{ fSetEMPhysics = value; }

inline void TG4PhysicsManager::SetMuonPhysics(G4bool value)
{ fSetMuonPhysics = value; }

inline void TG4PhysicsManager::SetHadronPhysics(G4bool value)
{ fSetHadronPhysics = value; }

inline void TG4PhysicsManager::SetOpticalPhysics(G4bool value)
{ fSetOpticalPhysics = value; }

inline void TG4PhysicsManager::SetSpecialCutsPhysics(G4bool value)
{ fSetSpecialCutsPhysics = value; }

inline void TG4PhysicsManager::SetSpecialControlsPhysics(G4bool value)
{ fSetSpecialControlsPhysics = value; }

inline TG4ModularPhysicsList* TG4PhysicsManager::GetPhysicsList() const
{ return fPhysicsList; }

#endif //TG4_PHYSICS_MANAGER_H

