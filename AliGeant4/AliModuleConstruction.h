// $Id$
// Category: geometry
//
// Abstract base class for modular construction of geometry,
// providing methods for browsing geometry (list volumes trees, 
// visualization).

#ifndef ALI_MODULE_CONSTRUCTION_H
#define ALI_MODULE_CONSTRUCTION_H

#include <globals.hh>

class AliLVStructure;
class AliModuleConstructionMessenger;
class AliModule;

class G4VPhysicalVolume;
class G4LogicalVolume;
#ifdef ALICE_VISUALIZE
class G4Colour;
#endif

class AliModuleConstruction
{
  public:
    AliModuleConstruction(G4String moduleName);
    AliModuleConstruction(const AliModuleConstruction& right);
    // --> protected
    // AliModuleConstruction();
    virtual ~AliModuleConstruction();

    // operators
    AliModuleConstruction& operator=(const AliModuleConstruction &right);
    G4int operator==(const AliModuleConstruction& right) const;
    G4int operator!=(const AliModuleConstruction& right) const;

    // methods
    virtual void Construct() = 0;
    void ListAllLVTree();
    void ListAllLVTreeLong();
    void ListLVTree(G4String lvName);
    void ListLVTreeLong(G4String lvName);
    G4LogicalVolume* FindLogicalVolume(G4String name, 
                                       G4bool silent = false) const;

    // set methods
    void SetDetFrame(G4bool warn = true);
    void SetDetFrame(G4String frameName, G4bool warn = true);
    void SetReadGeometry(G4bool readGeometry);
    void SetWriteGeometry(G4bool writeGeometry);
#ifdef ALICE_VISUALIZE
    void SetDetVisibility(G4bool visibility);
    void SetLVTreeVisibility(G4LogicalVolume* lv, G4bool visibility);
    void SetVolumeVisibility(G4LogicalVolume* lv, G4bool visibility);
    void SetDetColour(G4String colName);
    void SetLVTreeColour(G4LogicalVolume* lv, G4String colName);
    void SetVolumeColour(G4LogicalVolume* lv, G4String colName);     
#endif

    // get methods
    G4String GetDetName() const;
    G4LogicalVolume* GetDetFrame() const;
    AliModule* GetAliModule() const;
    G4bool GetReadGeometry() const;
    G4bool GetWriteGeometry() const;
    G4String GetDataFilePath() const;

  protected:
    AliModuleConstruction(); 

    // data members
    G4String            fModuleName;      //module name
    G4String            fModuleFrameName; //module frame name
                                          //(used for retrieving the frame LV)
    G4LogicalVolume*    fModuleFrameLV;   //module frame logical volume

    // to be moved to AliSingleModuleConstruction
    // in order to make AliModuleConstruction independent on
    // AliRoot
    AliModule*          fAliModule;       //AliModule
    G4bool              fReadGeometry;    //if true: geometry is read from file
    G4bool              fWriteGeometry;   //if true: geometry is written to file
    G4String            fDataFilePath;    //path to geometry data file

  private:
    // methods
    void RegisterLogicalVolume(G4LogicalVolume* lv, G4String path, 
           AliLVStructure& lvStructure);

    // data members
    AliModuleConstructionMessenger*  fMessenger; //messenger     
};

// inline methods

inline void AliModuleConstruction::SetReadGeometry(G4bool readGeometry)
{ fReadGeometry = readGeometry; }  

inline void AliModuleConstruction::SetWriteGeometry(G4bool writeGeometry)
{ fWriteGeometry = writeGeometry; }  

inline G4String AliModuleConstruction::GetDetName() const
{ return fModuleName; }

inline G4LogicalVolume* AliModuleConstruction::GetDetFrame() const
{ return fModuleFrameLV; }

inline AliModule* AliModuleConstruction::GetAliModule() const
{ return fAliModule; }

inline G4bool AliModuleConstruction::GetReadGeometry() const
{ return fReadGeometry; }

inline G4bool AliModuleConstruction::GetWriteGeometry() const
{ return fWriteGeometry; }

inline G4String AliModuleConstruction::GetDataFilePath() const
{ return fDataFilePath; }

#endif //ALI_MODULE_CONSTRUCTION_H

