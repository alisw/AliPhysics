// $Id$
// Category: global
//
// See the class description in the header file.

#include "AliFiles.h"
#include "AliGlobals.h"

#include <stdlib.h>

// static data members

AliFiles* AliFiles::fgInstance = 0;
const G4String AliFiles::fgkTop = getenv("AG4_INSTALL");    
const G4String AliFiles::fgkDefaultMacroName = "Config";
const G4String AliFiles::fgkDefaultG3CallsName = "g3calls";
const G4String AliFiles::fgkRootMacroExtension = ".C";
const G4String AliFiles::fgkG4MacroExtension = ".in";
const G4String AliFiles::fgkG3CallsExtension = ".dat";
const G4String AliFiles::fgkXMLFileExtension = ".xml";   

AliFiles::AliFiles()
  : fMacroName(fgkDefaultMacroName),
    fG3CallsName(fgkDefaultG3CallsName)
{
//    
  if (fgInstance) {
    AliGlobals::Exception(
      "AliFiles: attempt to create two instances of singleton.");
  };
      
  fgInstance = this;      
}
  
AliFiles::AliFiles(const G4String& config)  
  : fMacroName(config),
    fG3CallsName(fgkDefaultG3CallsName)
{
//
  if (fgInstance) {
    AliGlobals::Exception(
      "AliFiles: attempt to create two instances of singleton.");
  };
      
  fgInstance = this;      
}
  
AliFiles::AliFiles(const G4String& config, const G4String& g3calls)  
  : fMacroName(config),
    fG3CallsName(g3calls)
{
//
  if (fgInstance) {
    AliGlobals::Exception(
      "AliFiles: attempt to create two instances of singleton.");
  };
      
  fgInstance = this;      
}
  
AliFiles::~AliFiles() {
//
}

// private methods 

G4String AliFiles::GetMacroPath(const G4String& macroName,
                                const G4String& moduleName,
                                G4bool isStructure) const
{
// Returns the filepath to Config.C/in with filename extension:
// $AG4_INSTALL/macro/XXX/Config
// $AG4_INSTALL/macro/STRUCT/XXXConfig
// ---
    
  G4String name = fgkTop + "/macro/";
  
  if (!isStructure) 
    name = name + moduleName + "/" + macroName;
  else  			      
    name = name + "STRUCT/"+ fgkDefaultMacroName + moduleName;
			      
  return name;
}  			      
			      
// public methods

G4String AliFiles::GetRootMacroPath() const
{
// Returns the filepath:
// $ALICE_ROOT/macros/Config.C
// ---
    
  G4String name 
    = fgkTop + "/../macros/" + fMacroName + fgkRootMacroExtension;
			      
  return name;
}  			      
			      
G4String AliFiles::GetRootMacroPath(const G4String& moduleName,
                                    G4bool isStructure) const
{
// Returns the filepath:
// $AG4_INSTALL/macro/XXX/Config.C
// $AG4_INSTALL/macro/STRUCT/XXXConfig.C
// ---
    
  G4String name = GetMacroPath(fMacroName, moduleName, isStructure);
  name = name + fgkRootMacroExtension;
			      
  return name;
}  			      
			      
G4String AliFiles::GetG4MacroPath(const G4String& moduleName, 
                                  G4bool isStructure) const
{
// Returns the filepath:
// $AG4_INSTALL/macro/XXX/Config.in
// $AG4_INSTALL/macro/STRUCT/XXXConfig.in
// ---
    
  G4String name = GetMacroPath(fgkDefaultMacroName, moduleName, isStructure); 
  name = name + fgkG4MacroExtension;
			      
  return name;
}  			      
			      				  				  
G4String AliFiles::GetG3CallsDatPath(const G4String& moduleName, 
                              G4int moduleVersion, G4bool isStructure) const
{
// Returns the filepath:
// $AG4_INSTALL/macro/XXX/g3calls_vN.dat
// $AG4_INSTALL/macro/STRUCT/g3calls_XXXvN.dat
// ---
    
  G4String version("v");
  AliGlobals::AppendNumberToString(version, moduleVersion);

  G4String name = fgkTop + "/macro/";
  
  if (!isStructure) 
    name = name + moduleName + "/" + fG3CallsName + "_";
  else  			      
    name = name + "STRUCT/" + fG3CallsName + "_" + moduleName;

  name = name + version + fgkG3CallsExtension;
			      
  return name;
}  			      
			      
G4String AliFiles::GetXMLFilePath(const G4String& moduleName,
                                  G4int moduleVersion) const
{
// Returns the filepath:
// $AG4_INSTALL/xml/XXXvN.xml
// ---
    
  G4String version = "v";
  AliGlobals::AppendNumberToString(version, moduleVersion); 

  G4String name 
    = fgkTop + "/xml/" + moduleName + version + fgkXMLFileExtension;
  
  return name;
}
