// $Id$
// Category: physics
//
// Geant4 implementation of the MonteCarlo interface methods                    
// for building Geant4 physics and access to it

#ifndef TG4_PHYSICS_MANAGER_H
#define TG4_PHYSICS_MANAGER_H

#include "TG4Globals.h"
#include "TG4NameMap.h"
#include "TG4IntMap.h"

#include <Rtypes.h>
#include "AliMCProcess.h"

#include <globals.hh>

class AliDecayer;
class TG4ParticlesManager;
class TG4G3PhysicsManager;

class G4ParticleDefinition;
class G4VProcess;
class G4VModularPhysicsList;

class TG4PhysicsManager
{
  public:
    TG4PhysicsManager(G4VModularPhysicsList* physicsList);
    // --> protected
    // TG4PhysicsManager();
    // TG4PhysicsManager(const TG4PhysicsManager& right);
    virtual ~TG4PhysicsManager();

    // static access method
    static TG4PhysicsManager* Instance();
        
    // methods
    void BuildPhysics();

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

    // set methods
    void SetPhysicsList(G4VModularPhysicsList* physicsList);
    void SetEMPhysics(G4bool value);
    void SetOpticalPhysics(G4bool value);
    void SetHadronPhysics(G4bool value);
    void SetSpecialCutsPhysics(G4bool value);
    void SetSpecialControlsPhysics(G4bool value);
    
    // get methods
    G4VModularPhysicsList* GetPhysicsList() const; 
   
  protected:
    TG4PhysicsManager();
    TG4PhysicsManager(const TG4PhysicsManager& right);

    // operators
    TG4PhysicsManager& operator=(const TG4PhysicsManager& right);

  private:
    // methods
    void FillProcessMap();

    // static data members
    static TG4PhysicsManager*  fgInstance; //this instance
    
    // data members
    TG4ParticlesManager*   fParticlesManager; //particles manager
    TG4G3PhysicsManager*   fG3PhysicsManager; //G3 physics manager
    G4VModularPhysicsList* fPhysicsList; //physics list
    AliDecayer*            fDecayer;     //external decayer
    TG4IntMap              fProcessMap;  //the mapping between G4 process names
                                         //and AliMCProcess codes
    G4bool  fSetEMPhysics;          //electromagnetic physics control
    G4bool  fSetOpticalPhysics;     //optical physics control
    G4bool  fSetHadronPhysics;      //hadron physics control
    G4bool  fSetSpecialCutsPhysics; //special cuts process control 
                                    //(under development)                  
    G4bool  fSetSpecialControlsPhysics;//special controls process control
                                    //(under development)
};

// inline methods

inline TG4PhysicsManager* TG4PhysicsManager::Instance() 
{ return fgInstance; }

inline void TG4PhysicsManager::SetPhysicsList(G4VModularPhysicsList* physicsList)
{ fPhysicsList = physicsList; }

inline void TG4PhysicsManager::SetEMPhysics(G4bool value)
{ fSetEMPhysics = value; }

inline void TG4PhysicsManager::SetOpticalPhysics(G4bool value)
{ fSetOpticalPhysics = value; }

inline void TG4PhysicsManager::SetHadronPhysics(G4bool value)
{ fSetHadronPhysics = value; }

inline void TG4PhysicsManager::SetSpecialCutsPhysics(G4bool value)
{ fSetSpecialCutsPhysics = value; }

inline void TG4PhysicsManager::SetSpecialControlsPhysics(G4bool value)
{ fSetSpecialControlsPhysics = value; }

inline G4VModularPhysicsList* TG4PhysicsManager::GetPhysicsList() const
{ return fPhysicsList; }

#endif //TG4_PHYSICS_MANAGER_H

