// $Id$
// Category: physics
//
// This class provides a Geant3 way control
// to Geant4 physics. 
// The G3 cuts and process controls are
// stored in fCutVector and fControlVector vectors.
// These special cuts/controls are activated 
// by registering their physics constructors
// (TG4PhysicsConstructorSpecialCuts, G4PhysicsConstructorSpecialControl)
// to the modular physics list (TG4ModularPhysicsList)
// by physics manager (TG4PhysicsManager).

#ifndef TG4_G3_PHYSICS_MANAGER_H
#define TG4_G3_PHYSICS_MANAGER_H

#include "TG4Globals.h"
#include "TG4NameMap.h"
#include "TG4IntMap.h"
#include "TG4G3Cut.h"
#include "TG4G3Control.h"
#include "TG4G3ParticleWSP.h"

#include <Rtypes.h>
#include "AliMCProcess.h"

#include <globals.hh>

class TG4G3CutVector;
class TG4G3ControlVector;

class G4ParticleDefinition;
class G4VProcess;

class TG4G3PhysicsManager
{
  public:
    TG4G3PhysicsManager();
    // --> protected
    // TG4G3PhysicsManager(const TG4G3PhysicsManager& right);
    virtual ~TG4G3PhysicsManager();
        
    // static access method
    static TG4G3PhysicsManager* Instance();

    // methods
    void Lock();     
    void CheckLock();     
    G4VProcess*  FindProcess(G4String processName) const;
    G4bool CheckCutWithTheVector(
             G4String name, G4double value, TG4G3Cut& cut);   
    G4bool CheckControlWithTheVector(
             G4String name, G4double value, TG4G3Control& control); 
    G4bool CheckCutWithG3Defaults(
             G4String name, G4double value, TG4G3Cut& cut); 
    G4bool CheckControlWithG3Defaults(
             G4String name, G4double value, TG4G3Control& control); 

    // set methods
    void SetCut(TG4G3Cut cut, G4double cutValue);
    void SetProcess(TG4G3Control control, G4int controlValue);
    void SetG3DefaultCuts();                    
    void SetG3DefaultControls();               
    
    // get methods
    G4bool IsSpecialCuts() const;
    G4bool IsSpecialControls() const;
    TG4G3CutVector* GetCutVector() const;    
    TG4G3ControlVector* GetControlVector() const;   
    TG4boolVector* GetIsCutVector() const;
    TG4boolVector* GetIsControlVector() const;
    TG4G3Cut  GetG3Cut(G4String cutName);
    TG4G3Control GetG3Control(G4String controlName);
    TG4G3ParticleWSP GetG3ParticleWSP(G4ParticleDefinition* particle) const;
    void GetG3ParticleWSPName(G4int particleWSP, G4String& name) const;
    
  protected:
    TG4G3PhysicsManager(const TG4G3PhysicsManager& right);

    // operators
    TG4G3PhysicsManager& operator=(const TG4G3PhysicsManager& right);

  private:
    // methods
    void FillG3CutNameVector();
    void FillG3ControlNameVector();

    // set methods
    void SwitchIsCutVector(TG4G3Cut cut);
    void SwitchIsControlVector(TG4G3Control control);

    // static data members
    static TG4G3PhysicsManager*  fgInstance; //this instance

    // data members
    TG4G3CutVector*      fCutVector;    //TG4CutVector  
    TG4G3ControlVector*  fControlVector;//TG4ControlVector
    TG4boolVector*    fIsCutVector;  //vector of booleans which cuts are set
    TG4boolVector*    fIsControlVector; //vector of booleans which controls are set
    TG4StringVector   fG3CutNameVector; //vector of cut parameters names     
    TG4StringVector   fG3ControlNameVector; //vector of process control 
                                        // parameters names   
    G4bool            fLock;  //if true: cut/control vectors cannot be modified
};

// inline methods

inline TG4G3PhysicsManager* TG4G3PhysicsManager::Instance() 
{ return fgInstance; }

inline void TG4G3PhysicsManager::Lock() 
{ fLock = true; }

inline TG4G3CutVector* TG4G3PhysicsManager::GetCutVector() const
{ return fCutVector; }

inline TG4G3ControlVector* TG4G3PhysicsManager::GetControlVector() const
{ return fControlVector; }

inline TG4boolVector* TG4G3PhysicsManager::GetIsCutVector() const
{ return fIsCutVector; }

inline TG4boolVector* TG4G3PhysicsManager::GetIsControlVector() const
{ return fIsControlVector; }

#endif //TG4_PHYSICS_MANAGER_H

