// $Id$
// Category: physics
//
// Physics list (mandatory) class.
// According to:
// ExN04PhysicsList.hh, GEANT4 tag Name: geant4-01-01

#ifndef TG4_PHYSICS_LIST_H
#define TG4_PHYSICS_LIST_H

#include <G4VUserPhysicsList.hh>
#include <globals.hh>

class TG4PhysicsListMessenger;

class G4VProcess;

class TG4PhysicsList: public G4VUserPhysicsList
{
  public:
    TG4PhysicsList();
    TG4PhysicsList(const TG4PhysicsList& right);
    virtual ~TG4PhysicsList();
  
    // operators
    TG4PhysicsList& operator=(const TG4PhysicsList& right);
    
    // methods
    void PrintAllProcesses() const;
    void SetProcessActivation();
    
    // set methods
    void SetOptical(G4bool optical);
    void SetHadron(G4bool hadron);
    void SetSpecialCuts(G4bool specialCuts);
    void SetSpecialFlags(G4bool specialFlags);

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
  
         // construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    void ConstructHad();
    void ConstructOp();
    void ConstructNewSpecialCuts();
    void ConstructSpecialCuts();
    void ConstructSpecialFlags();

         // construct all particles in each category
    void ConstructAllBosons();
    void ConstructAllLeptons();
    void ConstructAllMesons();
    void ConstructAllBaryons();
    void ConstructAllIons();
    void ConstructAllShortLiveds();

    // set methods
    virtual void SetCuts();
   
  private:
    // methods
    G4VProcess* FindProcess(G4String processName) const;

	 // only for tests - to be removed  
    void InActivateEM();
    void InActivateProcess(G4String processName, 
           G4ParticleDefinition* particle);  
	   
 
    // data members
    G4bool  fSetOptical;     //optical processes control
    G4bool  fSetHadron;      //hadron processes control
    G4bool  fSetSpecialCuts; //special cuts process control (under development)                  
    G4bool  fSetSpecialFlags;//special flags process control (under development)

    TG4PhysicsListMessenger*  fMessenger;   //messenger
};

// inline methods

inline void TG4PhysicsList::SetOptical(G4bool optical)
{ fSetOptical = optical; }

inline void TG4PhysicsList::SetHadron(G4bool hadron)
{ fSetHadron = hadron; }

inline void TG4PhysicsList::SetSpecialCuts(G4bool specialCuts)
{ fSetSpecialCuts = specialCuts; }

inline void TG4PhysicsList::SetSpecialFlags(G4bool specialFlags)
{ fSetSpecialFlags = specialFlags; }

#endif //TG4_PHYSICS_LIST_H

