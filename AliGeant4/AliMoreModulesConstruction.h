// $Id$
// Category: geometry
//
// Class for geometry construction of a set of dependent
// modules (AliModule). 


#ifndef ALI_MORE_MODULES_CONSTRUCTION_H
#define ALI_MORE_MODULES_CONSTRUCTION_H

#include "AliModuleType.h"

#include <globals.hh>

#include <g4std/vector>

class AliSingleModuleConstruction;
class AliSDManager;
class AliFiles;

class AliMoreModulesConstruction
{  
  typedef G4std::vector<AliSingleModuleConstruction*>  
                                   AliSingleModuleConstructionPtrVector;

  public:
    AliMoreModulesConstruction();
    AliMoreModulesConstruction(const AliMoreModulesConstruction& right);
    virtual ~AliMoreModulesConstruction();

    // operators
    AliMoreModulesConstruction& operator=(
                               const AliMoreModulesConstruction& right);

    // methods
    void AddModule(G4String moduleName, G4int version, 
                   AliModuleType moduleType);
    void Configure(const AliFiles& files);
    void Construct();
    
    // get methods
    G4int GetNofModules() const;
    AliSingleModuleConstruction* GetModuleConstruction(G4int i) const;
        
  private:    
    // methods
    void CreateSensitiveDetectors(G4bool allLVSensitive);
    void CreateSensitiveDetectors1();
    void CreateSensitiveDetectors2();

    // data members
    AliSingleModuleConstructionPtrVector  fModuleConstructionVector; //..
                                //vector of AliSingleModuleConstruction
    AliSDManager*  fSDManager;  //AliSDManager						   
};						   

// inline methods

inline G4int AliMoreModulesConstruction::GetNofModules() const
{ return fModuleConstructionVector.size(); }

inline AliSingleModuleConstruction* 
  AliMoreModulesConstruction::GetModuleConstruction(G4int i) const 
{ return fModuleConstructionVector[i]; }

#endif //ALI_MORE_MODULES_CONSTRUCTION_H
