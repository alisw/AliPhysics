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
#include "TG3Cut.h"
#include "TG3Flag.h"
#include "TG3ParticleWSP.h"

#include <Rtypes.h>
#include "AliMCProcess.h"

#include <globals.hh>


class TG4CutVector;
class TG4FlagVector;
class TG4PhysicsList;

class G4ParticleDefinition;

class TG4PhysicsManager
{
  public:
    TG4PhysicsManager();
    // --> protected
    // TG4PhysicsManager(const TG4PhysicsManager& right);
    virtual ~TG4PhysicsManager();

    // static access method
    static TG4PhysicsManager* Instance();
        
    // methods
    void BuildPhysics();

    // set methods
    void SetCut(const char* cutName, Float_t cutValue);
    void SetProcess(const char* flagName, Int_t flagValue);
    Float_t Xsec(char* reac, Float_t energy, Int_t part, Int_t mate);
 
        // particle table usage         
    Int_t IdFromPDG(Int_t pdgID) const;
    Int_t PDGFromId(Int_t mcID) const;
    void  DefineParticles();      
    
    //
    // methods for Geant4 only 
    //

    void Lock();     
    void SetProcessActivation();  
    G4int GetPDGEncodingFast(G4ParticleDefinition* particle);
    AliMCProcess GetMCProcess(const G4String& g4ProcessName);
    G4bool CheckCutWithCutVector(
             G4String name, G4double value, TG3Cut& cut);   
    G4bool CheckFlagWithFlagVector(
             G4String name, G4double value, TG3Flag& flag); 
    G4bool CheckCutWithG3Defaults(
             G4String name, G4double value, TG3Cut& cut); 
    G4bool CheckFlagWithG3Defaults(
             G4String name, G4double value, TG3Flag& flag); 

    // set methods
    void SetPhysicsList(TG4PhysicsList* physicsList);
    void SetG3DefaultCuts();                    
    void SetG3DefaultProcesses();               
    
    // get methods
    G4bool IsLock() const;   
    G4bool IsSpecialCuts() const;
    G4bool IsSpecialFlags() const;
    TG4CutVector*  GetCutVector() const;    
    TG4FlagVector* GetFlagVector() const;   
    TG4boolVector* GetIsCutVector() const;
    TG4boolVector* GetIsFlagVector() const;
    TG3ParticleWSP GetG3ParticleWSP(G4ParticleDefinition* particle) const;
    void GetG3ParticleWSPName(G4int particleWSP, G4String& name) const;
    
  protected:
    TG4PhysicsManager(const TG4PhysicsManager& right);

    // operators
    TG4PhysicsManager& operator=(const TG4PhysicsManager& right);

  private:
    // methods
    void LockException() const;
    void FillG3CutNameVector();
    void FillG3FlagNameVector();
    void FillProcessMap();
    G4int GetPDGEncoding(G4ParticleDefinition* particle);
    G4int GetPDGEncoding(G4String particleName);

    // set methods
    void SetCut(TG3Cut cut, G4double cutValue);
    void SetProcess(TG3Flag flag, G4int flagValue);
    void SwitchIsCutVector(TG3Cut cut);
    void SwitchIsFlagVector(TG3Flag flag);

    // get methods
    TG3Cut  GetG3Cut(G4String cutName);
    TG3Flag GetG3Flag(G4String flagName);

    // static data members
    static TG4PhysicsManager*  fgInstance; //this instance
    
    // data members
    G4bool           fLock;         //if true: cut/flag vectors cannot be modified
    TG4PhysicsList*  fPhysicsList;  //physics list
    TG4CutVector*    fCutVector;    //TG4CutVector  
    TG4FlagVector*   fFlagVector;   //TG4FlagVector
    TG4boolVector*   fIsCutVector;  //vector of booleans which cuts are set
    TG4boolVector*   fIsFlagVector; //vector of booleans which flags are set
    TG4StringVector  fG3CutNameVector;  //vector of cut parameters names     
    TG4StringVector  fG3FlagNameVector; //vector of process control parameters
                                        //names   
    TG4NameMap       fParticleNameMap;  //the mapping between G4 particle names
                                        //and TDatabasePDG names 
    TG4IntMap        fParticlePDGMap;   //the mapping between G4 particle names
                                        //and TDatabasePDG codes
    TG4IntMap        fProcessMap;       //the mapping between G4 process names
                                        //and AliMCProcess codes
};

// inline methods

inline TG4PhysicsManager* TG4PhysicsManager::Instance() 
{ return fgInstance; }

inline void TG4PhysicsManager::Lock() 
{ fLock = true; }

inline void TG4PhysicsManager::SetPhysicsList(TG4PhysicsList* physicsList)
{ fPhysicsList = physicsList; }

inline G4bool TG4PhysicsManager::IsLock() const
{ return fLock; }

inline TG4CutVector* TG4PhysicsManager::GetCutVector() const
{ return fCutVector; }

inline TG4FlagVector* TG4PhysicsManager::GetFlagVector() const
{ return fFlagVector; }

inline TG4boolVector* TG4PhysicsManager::GetIsCutVector() const
{ return fIsCutVector; }

inline TG4boolVector* TG4PhysicsManager::GetIsFlagVector() const
{ return fIsFlagVector; }

#endif //TG4_PHYSICS_MANAGER_H

