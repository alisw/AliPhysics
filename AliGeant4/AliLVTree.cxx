// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliLVTree
// ---------------------------
// See the class description in the header file.

#include "AliLVTree.h"
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

AliLVTree* AliLVTree::fgInstance = 0;

//_____________________________________________________________________________
AliLVTree* AliLVTree::Instance() 
{
// Static singleton access method.
// ---

  if (!fgInstance) new AliLVTree();
  
  return fgInstance;
}  

//_____________________________________________________________________________
AliLVTree::AliLVTree()
  : fMessenger(this) 
{
// Protected singleton constructor.
// ---

  fgInstance = this;
}

//_____________________________________________________________________________
AliLVTree::AliLVTree(const AliLVTree& right)
  : fMessenger(this)
{
// Protected singleton copy constructor.
// ---
//
  AliGlobals::Exception(
    "Attempt to copy AliLVTree singleton.");
}

//_____________________________________________________________________________
AliLVTree::~AliLVTree() {
//
}

// operators

//_____________________________________________________________________________
AliLVTree& AliLVTree::operator=(const AliLVTree& right)
{    
  // check assignement to self
  if (this == &right) return *this;

  AliGlobals::Exception(
    "Attempt to assign AliLVTree singleton.");
    
  return *this;  

}

// protected methods

//_____________________________________________________________________________
void AliLVTree::RegisterLogicalVolume(G4LogicalVolume* lv, const G4String& path,
				      AliLVStructure& lvStructure) const
{
// Registers logical volume lv in the structure.
// ---        

  lvStructure.AddNewVolume(lv, path);
  
  // register daughters
  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters>0) {
    G4String previousName = "";
    for (G4int i=0; i<nofDaughters; i++) {
      G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume();
      G4String currentName = lvd->GetName();
      if (currentName != lv->GetName() && currentName != previousName) { 
        G4String newPath = path + lv->GetName() +"/";
        RegisterLogicalVolume(lvd, newPath, lvStructure);
	previousName = currentName;
      }
    }
  }     
}          

//_____________________________________________________________________________
void AliLVTree::Warn(const G4String& where, const G4String& lvName) const			       
{
// Prints warning "volume not found".
// ---
  
   G4String text("AliLVTree::" + where + ": " + lvName + " volume not found.");
   AliGlobals::Warning(text);
}

//_____________________________________________________________________________
void AliLVTree::Warn(const G4String& where) const
{			       
// Prints warning "volume not specified".
// ---
  
   G4String text("AliLVTree::" + where + ": " + " volume not specified.");
   AliGlobals::Warning(text);
}

// public methods

//_____________________________________________________________________________
void AliLVTree::List(const G4String& lvName) const
{
// Lists logical volumes tree (daughters) of the logical volume 
// with specified lvName.
// ---- 

  G4LogicalVolume* lv 
    = TG4GeometryServices::Instance()->FindLogicalVolume(lvName);

  if (lv) {
    G4String path = "";
    AliLVStructure lvStructure(path);
    RegisterLogicalVolume(lv, path, lvStructure);
    lvStructure.ListTree();
  }
  else 
    Warn("List", lvName);
}    

//_____________________________________________________________________________
void AliLVTree::List(G4LogicalVolume* lv) const
{
// Lists logical volumes tree of the given lv.
// ---- 

  if (!lv) {
    Warn("List");
    return; 
  }  
  
  List(lv->GetName());
}

//_____________________________________________________________________________
void AliLVTree::ListLong(const G4String& lvName) const
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
  else 
    Warn("ListLong", lvName);
}

//_____________________________________________________________________________
void AliLVTree::ListLong(G4LogicalVolume* lv) const
{
// Lists logical volumes tree with numbers of daughters 
// (physical volumes) of the given lv.
// ---- 

  if (!lv) {
    Warn("ListLong");
    return; 
  }  
  
  ListLong(lv->GetName());
}

#ifdef G4VIS_USE

//_____________________________________________________________________________
void AliLVTree::SetLVTreeVisibility(G4LogicalVolume* lv, 
                                    G4bool visibility) const
{ 
// Sets visibility to the logical volumes tree (daughters) of 
// the logical volume lv.
// ---

  if (!lv) {
    Warn("SetLVTreeVisibility");
    return; 
  }  
  
  G4String path = "";
  AliLVStructure lvStructure(path);
  RegisterLogicalVolume(lv, path, lvStructure);
  lvStructure.SetTreeVisibility(visibility);
}

//_____________________________________________________________________________
void AliLVTree::SetVolumeVisibility(G4LogicalVolume* lv, 
                                    G4bool visibility) const
{ 
// Sets visibility to the specified logical volume.
// ---

  if (!lv) {
    Warn("SetVolumeVisibility");
    return; 
  }  
  
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

//_____________________________________________________________________________
void AliLVTree::SetLVTreeColour(G4LogicalVolume* lv, 
                                const G4String& colName) const
{ 
// Sets colour to the logical volumes tree (daughters) of 
// the logical volume lv.
// ---

  if (!lv) {
    Warn("SetLVTreeColour");
    return; 
  }  
  
  G4String path = "";
  AliLVStructure lvStructure(path);
  RegisterLogicalVolume(lv, path, lvStructure);
  lvStructure.SetTreeVisibility(true);
  lvStructure.SetTreeColour(colName);
}

//_____________________________________________________________________________
void AliLVTree::SetVolumeColour(G4LogicalVolume* lv,
                                const G4String& colName) const
{
// Sets colour to the specified logical volume.
// ---

  if (!lv) {
    Warn("SetVolumeColour");
    return; 
  }  
  
  const G4VisAttributes* kpVisAttributes = lv->GetVisAttributes ();
  delete kpVisAttributes;

  G4VisAttributes* newVisAttributes = new G4VisAttributes(); 

  AliColourStore* pColours = AliColourStore::Instance();
  const G4Colour kColour = pColours->GetColour(colName);
  newVisAttributes->SetVisibility(true); 
  newVisAttributes->SetColour(kColour);

  lv->SetVisAttributes(newVisAttributes);
}

#endif //G4VIS_USE

