// $Id$
// Category: physics
//
// Empty physics list for geometry tests.

#ifndef ALI_EMPTY_PHYSICS_LIST_H
#define ALI_EMPTY_PHYSICS_LIST_H

#include <G4VUserPhysicsList.hh>
#include <globals.hh>

class AliEmptyPhysicsList: public G4VUserPhysicsList
{
  public:
    AliEmptyPhysicsList();
    virtual ~AliEmptyPhysicsList();

  protected:
    // methods
    // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
 
    virtual void SetCuts();
    
    // construct particles methods
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBarions();

    // construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();    
};

#endif //ALI_EMPTY_PHYSICS_LIST_H



