// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliModuleConstruction
// ---------------------------
// Composite of AliDetector and additional parameters
// used for optional writing/reading geometry ASCII
// files, interactive detector setup and processing 
// of module related G4 macros. 


#ifndef ALI_MODULE_CONSTRUCTION_H
#define ALI_MODULE_CONSTRUCTION_H

#include "AliModuleType.h"

#include <globals.hh>

class AliModule;
class AliFiles;

class AliModuleConstruction
{
  public:
    AliModuleConstruction(const G4String& moduleName, 
                          G4int version, 
			  AliModuleType moduleType = kDetector);
    AliModuleConstruction(const AliModuleConstruction& right);
    // --> protected
    // AliModuleConstruction();
    virtual ~AliModuleConstruction();

    // operators
    AliModuleConstruction& operator=(const AliModuleConstruction &right);
    G4int operator==(const AliModuleConstruction& right) const;
    G4int operator!=(const AliModuleConstruction& right) const;

    // methods
    void Configure();    

    // set methods
    void SetModuleType(AliModuleType type);
    void SetProcessConfig(G4bool processConfig);
    void SetReadGeometry(G4bool readGeometry);
    void SetWriteGeometry(G4bool writeGeometry);

    // get methods
    AliModule* GetAliModule() const;
    AliModuleType GetType() const;
    G4String   GetDetName() const;
    G4int      GetVersion() const;
    G4bool     GetProcessConfig() const;
    G4bool     GetReadGeometry() const;
    G4bool     GetWriteGeometry() const;
    G4String   GetDataFilePath() const;

  private:
    AliModuleConstruction(); 

    // data members
    AliModule*      fAliModule;     //AliModule
    G4String        fModuleName;    //module name
    AliModuleType   fType;          //module type (detector/structure)
    G4int           fVersion;       //module version
    G4bool          fProcessConfig; //control for processing Config.C
    G4bool          fReadGeometry;  //if true: geometry is read from file
    G4bool          fWriteGeometry; //if true: geometry is written to file
    G4String        fDataFilePath;  //path to geometry data file

    // data members
};

// inline methods

inline void AliModuleConstruction::SetModuleType(AliModuleType type)
{ fType = type; }

inline void AliModuleConstruction::SetProcessConfig(G4bool processConfig)
{ fProcessConfig = processConfig; }

inline void AliModuleConstruction::SetReadGeometry(G4bool readGeometry)
{ fReadGeometry = readGeometry; }  

inline void AliModuleConstruction::SetWriteGeometry(G4bool writeGeometry)
{ fWriteGeometry = writeGeometry; }  

inline AliModuleType AliModuleConstruction::GetType() const
{ return fType; }

inline AliModule* AliModuleConstruction::GetAliModule() const
{ return fAliModule; }

inline G4String AliModuleConstruction::GetDetName() const
{ return fModuleName; }

inline G4int AliModuleConstruction::GetVersion() const
{ return fVersion; }

inline G4bool AliModuleConstruction::GetProcessConfig() const
{ return fProcessConfig; }

inline G4bool AliModuleConstruction::GetReadGeometry() const
{ return fReadGeometry; }

inline G4bool AliModuleConstruction::GetWriteGeometry() const
{ return fWriteGeometry; }

inline G4String AliModuleConstruction::GetDataFilePath() const
{ return fDataFilePath; }

#endif //ALI_MODULE_CONSTRUCTION_H

