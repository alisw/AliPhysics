// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliModuleConstruction
// ---------------------------
// See the class description in the header file.

#include "AliModuleConstruction.h"
#include "AliGlobals.h"
#include "AliLVStructure.h"
#include "AliModule.h"
#ifdef ALICE_VISUALIZE
#include "AliColourStore.h"
#endif

#include "TG4GeometryServices.h"

#include <G4LogicalVolumeStore.hh>
#include <G4LogicalVolume.hh>
#ifdef ALICE_VISUALIZE
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#endif

//_____________________________________________________________________________
AliModuleConstruction::AliModuleConstruction(const G4String& moduleName) 
  : fModuleName(moduleName), 
    fModuleFrameName(moduleName),
    fModuleFrameLV(0),
    fAliModule(0),
    fReadGeometry(false),
    fWriteGeometry(false),
    fDataFilePath(""),    
    fMessenger(this, moduleName) {
//
}

//_____________________________________________________________________________
AliModuleConstruction::AliModuleConstruction(const AliModuleConstruction& right)
  : fMessenger(this, right.fModuleName)
{
//
  // copy stuff
  *this = right;
}

//_____________________________________________________________________________
AliModuleConstruction::AliModuleConstruction()
  : fModuleName(""), 
    fModuleFrameName(""),
    fModuleFrameLV(0),
    fAliModule(0),
    fReadGeometry(false),
    fWriteGeometry(false),
    fDataFilePath(""),    
    fMessenger(this, "") {
//
}

//_____________________________________________________________________________
AliModuleConstruction::~AliModuleConstruction()
{
//
  delete fAliModule;
}

// operators

//_____________________________________________________________________________
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
  //fMessenger = right.fMessenger;

  return *this;
}

//_____________________________________________________________________________
G4int 
AliModuleConstruction::operator==(const AliModuleConstruction& right) const
{
//    
  return 0;
}

//_____________________________________________________________________________
G4int 
AliModuleConstruction::operator!=(const AliModuleConstruction& right) const
{
//    
  G4int returnValue = 1;
  if (*this == right) returnValue = 0; 
  
  return returnValue;
}

// protected methods

//_____________________________________________________________________________
void AliModuleConstruction::RegisterLogicalVolume(
                                    G4LogicalVolume* lv, const G4String& path,
				    AliLVStructure& lvStructure) const
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

//_____________________________________________________________________________
void AliModuleConstruction::SetDetFrame(G4bool warn)
{ 
// The logical volume with name identical with
// fModuleName is retrieved from G4LogicalVolumeStore.
// ---

  fModuleFrameLV 
    = TG4GeometryServices::Instance()->FindLogicalVolume(fModuleFrameName, true);
  
  if (fModuleFrameLV == 0 && warn) {
    G4String text = "AliModuleConstruction: Detector frame for ";
    text = text + fModuleFrameName + " has not been found.";
    AliGlobals::Warning(text); 
  }  
}

//_____________________________________________________________________________
void AliModuleConstruction::SetDetFrame(const G4String& frameName, G4bool warn)
{ 
// The logical volume with frameName
// is retrieved from G4LogicalVolumeStore.
// ---

  fModuleFrameName = frameName;
  SetDetFrame(warn);
}

//_____________________________________________________________________________
void AliModuleConstruction::ListAllLVTree() const
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

//_____________________________________________________________________________
void AliModuleConstruction::ListAllLVTreeLong() const
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

//_____________________________________________________________________________
void AliModuleConstruction::ListLVTree(const G4String& lvName) const
{
// Lists logical volumes tree (daughters) of the logical volume 
// with specified lvName.
// ---- 

  G4LogicalVolume* lv 
    = TG4GeometryServices::Instance()->FindLogicalVolume(lvName);

  if (lv)
  {
    G4String path = "";
    AliLVStructure lvStructure(path);
    RegisterLogicalVolume(lv, path, lvStructure);
    lvStructure.ListTree();
  }
}

//_____________________________________________________________________________
void AliModuleConstruction::ListLVTreeLong(const G4String& lvName) const
{
// Lists logical volumes tree (daughters) of the logical volume 
// with specified lvName with numbers of daughters (physical volumes).
// ---- 

  G4LogicalVolume* lv 
    = TG4GeometryServices::Instance()->FindLogicalVolume(lvName);

  if (lv) {
    G4String path = "";
    AliLVStructure lvStructure(path);
    RegisterLogicalVolume(lv, path, lvStructure);
    lvStructure.ListTreeLong();
  }
}

#ifdef ALICE_VISUALIZE

//_____________________________________________________________________________
void AliModuleConstruction::SetDetVisibility(G4bool visibility) const
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


//_____________________________________________________________________________
void AliModuleConstruction::SetLVTreeVisibility(G4LogicalVolume* lv, 
                                                G4bool visibility) const
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

//_____________________________________________________________________________
void AliModuleConstruction::SetVolumeVisibility(G4LogicalVolume* lv, 
                                                G4bool visibility) const
{ 
// Sets visibility to the specified logical volume.
// ---

  if (lv) {
    const G4VisAttributes* kpVisAttributes = lv->GetVisAttributes ();
    G4VisAttributes* newVisAttributes; 
    if (kpVisAttributes) {
      G4Colour oldColour   = kpVisAttributes->GetColour();
      newVisAttributes = new G4VisAttributes(oldColour); 
    }  
    else
      newVisAttributes = new G4VisAttributes();
    delete kpVisAttributes;

    newVisAttributes->SetVisibility(visibility); 

    lv->SetVisAttributes(newVisAttributes);
  }
}

//_____________________________________________________________________________
void AliModuleConstruction::SetDetColour(G4String colName) const
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

//_____________________________________________________________________________
void AliModuleConstruction::SetLVTreeColour(G4LogicalVolume* lv, 
                                            const G4String& colName) const
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

//_____________________________________________________________________________
void AliModuleConstruction::SetVolumeColour(G4LogicalVolume* lv,
                                            const G4String& colName) const
{
// Sets colour to the specified logical volume.
// ---

  if (lv) {
    const G4VisAttributes* kpVisAttributes = lv->GetVisAttributes ();
    delete kpVisAttributes;

    G4VisAttributes* newVisAttributes = new G4VisAttributes(); 

    AliColourStore* pColours = AliColourStore::Instance();
    const G4Colour kColour = pColours->GetColour(colName);
    newVisAttributes->SetVisibility(true); 
    newVisAttributes->SetColour(kColour);

    lv->SetVisAttributes(newVisAttributes);
  }      
}

#endif //ALICE_VISUALIZE

