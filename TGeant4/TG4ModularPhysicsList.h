// $Id$
// Category: physics
//
// The modular physics list.
// Creates all standard particles.
// The physics processes have to be created
// using the G4VPhysicsCreator derived classes
// and registered to this physics list.
// Only G4Decay is created directly in this modular
// physics list.

#ifndef TG4_PHYSICS_LIST_H
#define TG4_PHYSICS_LIST_H

#include <G4VModularPhysicsList.hh>
#include <globals.hh>

class G4VProcess;

class TG4ModularPhysicsList: public G4VModularPhysicsList
{
  public:
    TG4ModularPhysicsList();
    virtual ~TG4ModularPhysicsList();
  
    // methods
    virtual void SetCuts();
    void SetProcessActivation();
    void PrintAllProcesses() const;
    
  protected:
    // methods
    virtual void ConstructParticle();
    virtual void ConstructProcess();

         // construct all particles in each category
    void ConstructAllBosons();
    void ConstructAllLeptons();
    void ConstructAllMesons();
    void ConstructAllBaryons();
    void ConstructAllIons();
    void ConstructAllShortLiveds();
    
        // construct general processes
    void ConstructGeneral();	
};

#endif //TG4_MODULAR_PHYSICS_LIST_H

