// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliModuleConstruction.h"
#include "AliGlobals.h"
#include "AliLVStructure.h"

#ifdef ALICE_VISUALIZE
#include "AliColourStore.h"

#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#endif //ALICE_VISUALIZE
#include <G4LogicalVolumeStore.hh>
#include <G4LogicalVolume.hh>

#include <fstream.h>

AliModuleConstruction::AliModuleConstruction(G4String moduleName) 
  : fModuleName(moduleName), 
    fModuleFrameName(moduleName),
    fModuleFrameLV(0),
    fAliModule(0),
    fReadGeometry(false),
    fWriteGeometry(false),
    fDataFilePath("")    
{
//
  moduleName.toLower();
  fMessenger = new AliModuleConstructionMessenger(this, moduleName);
}

AliModuleConstruction::AliModuleConstruction(const AliModuleConstruction& right)
{
//
  fModuleName = right.fModuleName; 
  fModuleFrameName = right.fModuleFrameName;
  fModuleFrameLV = right.fModuleFrameLV;
  fAliModule = right.fAliModule;
  fReadGeometry = right.fReadGeometry;
  fWriteGeometry = right.fWriteGeometry;
  fDataFilePath = right.fDataFilePath;

  // new messenger
  G4String moduleName = fModuleName;
  moduleName.toLower();
  fMessenger = new AliModuleConstructionMessenger(this, moduleName);
}

AliModuleConstruction::AliModuleConstruction()
  : fModuleName(""), 
    fModuleFrameName(""),
    fModuleFrameLV(0),
    fMessenger(0),
    fAliModule(0),
    fReadGeometry(false),
    fWriteGeometry(false),
    fDataFilePath("")    
{
//
}

AliModuleConstruction::~AliModuleConstruction()
{
//
  delete fMessenger;
  delete fAliModule;
}

// operators

AliModuleConstruction& 
AliModuleConstruction::operator=(const AliModuleConstruction& right)
{    
  // check assignement to self
  if (this == &right) return *this;
  
  fModuleName = right.fModuleName; 
  fModuleFrameName = right.fModuleFrameName;
  fModuleFrameLV = right.fModuleFrameLV;
  fAliModule = right.fAliModule;
  fReadGeometry = right.fReadGeometry;
  fWriteGeometry = right.fWriteGeometry;
  fDataFilePath = right.fDataFilePath;

  return *this;
}

G4int 
AliModuleConstruction::operator==(const AliModuleConstruction& right) const
{
//    
  return 0;
}

G4int 
AliModuleConstruction::operator!=(const AliModuleConstruction& right) const
{
//    
  G4int returnValue = 1;
  if (*this == right) returnValue = 0; 
  
  return returnValue;
}

// protected methods

void AliModuleConstruction::RegisterLogicalVolume(G4LogicalVolume* lv,
       G4String path, AliLVStructure& lvStructure)
{
// Registers logical volume lv in the structure.
// ---        

  G4String lvName = lv->GetName();
  lvStructure.AddNewVolume(lv, path);
  
  // register daughters
  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters>0) {
    G4String previousName = "";
    for (G4int i=0; i<nofDaughters; i++) {
      G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume();
      G4String currentName = lvd->GetName();
      if (currentName != lvName && currentName != previousName) { 
        G4String newPath = path + lvName +"/";
        RegisterLogicalVolume(lvd, newPath, lvStructure);
	previousName = currentName;
      }
    }
  }     
}          

// public methods

void AliModuleConstruction::SetDetFrame(G4bool warn)
{ 
// The logical volume with name identical with
// fModuleName is retrieved from G4LogicalVolumeStore.
// ---

  fModuleFrameLV = FindLogicalVolume(fModuleFrameName, true);
  
  if (fModuleFrameLV == 0 && warn) {
    G4String text = "AliModuleConstruction: Detector frame for ";
    text = text + fModuleFrameName + " has not been found.";
    AliGlobals::Warning(text); 
  }  
}

void AliModuleConstruction::SetDetFrame(G4String frameName, G4bool warn)
{ 
// The logical volume with frameName
// is retrieved from G4LogicalVolumeStore.
// ---

  fModuleFrameName = frameName;
  SetDetFrame(warn);
}

void AliModuleConstruction::ListAllLVTree()
{
// Lists all logical volumes tree if the frame logical volume 
// is defined.
// ---- 

  if (fModuleFrameLV) 
    ListLVTree(fModuleFrameLV->GetName());
  else {
    G4String text = "AliModuleConstruction::ListAllLVTree:\n";
    text = text + "    Detector frame is not defined.";    
    AliGlobals::Warning(text);
  }   
}

void AliModuleConstruction::ListAllLVTreeLong()
{
// Lists all logical volume tree if the frame logical volume 
// is defined with numbers of daughters (physical volumes).
// ---- 

  if (fModuleFrameLV) 
    ListLVTreeLong(fModuleFrameLV->GetName());
  else {
    G4String text = "AliModuleConstruction::ListAllLVTreeLong:\n";
    text = text + "    Detector frame is not defined.";    
    AliGlobals::Warning(text);
  }  
}

void AliModuleConstruction::ListLVTree(G4String lvName)
{
// Lists logical volumes tree (daughters) of the logical volume 
// with specified lvName.
// ---- 

  G4LogicalVolume* lv = FindLogicalVolume(lvName);
  if (lv)
  {
    G4String path = "";
    AliLVStructure lvStructure(path);
    RegisterLogicalVolume(lv, path, lvStructure);
    lvStructure.ListTree();
  }
}

void AliModuleConstruction::ListLVTreeLong(G4String lvName)
{
// Lists logical volumes tree (daughters) of the logical volume 
// with specified lvName with numbers of daughters (physical volumes).
// ---- 

  G4LogicalVolume* lv = FindLogicalVolume(lvName);
  if (lv) {
    G4String path = "";
    AliLVStructure lvStructure(path);
    RegisterLogicalVolume(lv, path, lvStructure);
    lvStructure.ListTreeLong();
  }
}

G4LogicalVolume* AliModuleConstruction::FindLogicalVolume(
                                          G4String name, G4bool silent) const
{
// Finds logical volume with specified name in G4LogicalVolumeStore.
// (To do: use this method only for retrieving detector frame;
//  add method FindLogicalVolumeInDet - that will search only
//  in the detector frame LVTree.)
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  
  for (G4int i=0; i<pLVStore->entries(); i++) {
    G4LogicalVolume* lv = pLVStore->at(i);
    if (lv->GetName() == name) return lv;
  }
  
  G4String text = "AliModuleConstruction: Logical volume "; 
  text = text + name + " does not exist.";
  if (!silent) AliGlobals::Warning(text);
  return 0;	       	         
}  

#ifdef ALICE_VISUALIZE

void AliModuleConstruction::SetDetVisibility(G4bool visibility)
{
// Sets visibility to all detector logical volumes if
// frame logical volume is defined.
// ---

  if (fModuleFrameLV) 
    SetLVTreeVisibility(fModuleFrameLV, visibility);  
  else  {
    G4String text = "AliModuleConstruction::SetDetVisibility:\n";
    text = text + "    Detector frame is not defined.";    
    AliGlobals::Warning(text);
  }  
}


void AliModuleConstruction::SetLVTreeVisibility(G4LogicalVolume* lv, 
                             G4bool visibility)
{ 
// Sets visibility to the logical volumes tree (daughters) of 
// the logical volume lv.
// ---

  if (lv) {
    G4String path = "";
    AliLVStructure lvStructure(path);
    RegisterLogicalVolume(lv, path, lvStructure);
    lvStructure.SetTreeVisibility(visibility);
  }
}

void AliModuleConstruction::SetVolumeVisibility(G4LogicalVolume* lv, 
                             G4bool visibility)
{ 
// Sets visibility to the specified logical volume.
// ---

  if (lv) {
    const G4VisAttributes* kpVisAttributes = lv->GetVisAttributes ();
    G4VisAttributes* newVisAttributes = new G4VisAttributes(kpVisAttributes); 
    delete kpVisAttributes;

    newVisAttributes->SetVisibility(visibility); 

    lv->SetVisAttributes(newVisAttributes);
  }
}

void AliModuleConstruction::SetDetColour(G4String colName)
{
// Sets colour to all detector logical volumes if
// frame logical volume is defined.
// ---

  if (fModuleFrameLV) 
    SetLVTreeColour(fModuleFrameLV, colName);  
  else { 
    G4String text = "AliModuleConstruction::SetDetColour:\n";
    text = text + "    Detector frame is not defined.";    
    AliGlobals::Warning(text);
  }  
}

void AliModuleConstruction::SetLVTreeColour(G4LogicalVolume* lv, 
                             G4String colName)
{ 
// Sets colour to the logical volumes tree (daughters) of 
// the logical volume lv.
// ---

  if (lv) {
    G4String path = "";
    AliLVStructure lvStructure(path);
    RegisterLogicalVolume(lv, path, lvStructure);
    lvStructure.SetTreeVisibility(true);
    lvStructure.SetTreeColour(colName);
  }
}

void AliModuleConstruction::SetVolumeColour(G4LogicalVolume* lv,
                             G4String colName)
{
// Sets colour to the specified logical volume.
// ---

  if (lv) {
    const G4VisAttributes* kpVisAttributes = lv->GetVisAttributes ();
    G4VisAttributes* newVisAttributes = new G4VisAttributes(kpVisAttributes); 
    delete kpVisAttributes;

    AliColourStore* pColours = AliColourStore::Instance();
    const G4Colour kColour = pColours->GetColour(colName);
    newVisAttributes->SetVisibility(true); 
    newVisAttributes->SetColour(kColour);

    lv->SetVisAttributes(newVisAttributes);
  }      
}

#endif //ALICE_VISUALIZE

