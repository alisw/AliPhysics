// $Id$
// Category: global
//
// Class for generating file names and paths.
// The input files:
// Config.C    - the basic AliRoot configuration file (Root macro)
//               When detector setup is defined interactively,
//               the alternative Config.C files per detector
//               are used.
// Config.in   - G4 specific configuration macros per detector
// The output files:
// g3calls.dat - the ASCII file with G3 geometry calls;
//               generation is switch on with /aliDet/writeGeometry command
//               in PreInit phase
// DetvN.xml   - XML geometry description file;                
//               generation is switch on with /aliDet/generateXML command
//               in Init phase

#ifndef ALI_FILES_H
#define ALI_FILES_H

#include <globals.hh>

class AliFiles
{
  public:
    AliFiles();
    AliFiles(const G4String& config);
    AliFiles(const G4String& config, const G4String& g3calls);
    virtual ~AliFiles();
    
    // static access method
    static AliFiles* Instance();

    // methods
    G4String GetRootMacroPath() const;
    G4String GetRootMacroPath(const G4String& moduleName,
                              G4bool isStructure) const;
    G4String GetG4MacroPath(const G4String& moduleName, 
                              G4bool isStructure) const;
    G4String GetG3CallsDatPath(const G4String& moduleName, 
                              G4int moduleVersion, G4bool isStructure) const;
    G4String GetXMLFilePath(const G4String& moduleName, 
                              G4int moduleVersion) const;
			      
    // set methods
    void SetMacroName(const G4String& name);			      
    void SetG3CallsName(const G4String& name);			      
 
    // get methods
    G4String GetMacroName() const;			      
    G4String GetG3CallsName() const;			      
    G4String GetDefaultMacroName() const;			      
    G4String GetDefaultG3CallsName() const;			      
       
  private: 
    // methods
    G4String GetMacroPath(const G4String& macroName,
                          const G4String& moduleName,
                          G4bool isStructure) const;    
        
    // static data members  
    static AliFiles*       fgInstance; //this instance
    static const G4String  fgkTop;     //top directory
    static const G4String  fgkDefaultMacroName;   // default config. macro name
    static const G4String  fgkDefaultG3CallsName; // default g3calls name
    static const G4String  fgkRootMacroExtension; //".C"  Root macro extension
    static const G4String  fgkG4MacroExtension;   //".in" G4 macro extension
    static const G4String  fgkG3CallsExtension;   //".dat"   
    static const G4String  fgkXMLFileExtension;   //".xml"   

    // data members  
    G4String        fMacroName;      //configuration macro name
    G4String        fG3CallsName;        //g3calls data file name  
};  

// inline methods

inline AliFiles* AliFiles::Instance() 
{ return fgInstance; }

inline void AliFiles::SetMacroName(const G4String& name)
{ fMacroName = name; }
		      
inline void AliFiles::SetG3CallsName(const G4String& name)
{ fG3CallsName = name; }

inline G4String AliFiles::GetMacroName() const
{ return fMacroName; }

inline G4String AliFiles::GetG3CallsName() const
{ return fG3CallsName; }

inline G4String AliFiles::GetDefaultMacroName() const
{ return fgkDefaultMacroName; }
			      
inline G4String AliFiles::GetDefaultG3CallsName() const			      
{ return fgkDefaultG3CallsName; }

#endif //ALI_FILES_H
