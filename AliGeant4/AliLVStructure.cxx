// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "AliLVStructure.h"
#include "AliGlobals.h"

#ifdef ALICE_VISUALIZE
#include "AliColourStore.h"

#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#endif //ALICE_VISUALIZE
#include <G4LogicalVolume.hh>

AliLVStructure::AliLVStructure(G4String path)
  : fPathName(path),
    fDirName(path),
    fVerboseLevel(0)
{
//
  G4int i = fDirName.length();
  if (i > 1) {
    fDirName.remove(i-1);
    G4int isl = fDirName.last('/');
    fDirName.remove(0,isl+1);
    fDirName += "/";
  }
}

AliLVStructure::AliLVStructure(const AliLVStructure& right)
{
  // copy vector of structures
  fStructures.clearAndDestroy();
  G4int i;
  for (i=0; i<right.fStructures.entries(); i++) {
    // new full structure tree has to be created
    AliLVStructure* rhsStructure = right.fStructures[i];
    fStructures.insert(new AliLVStructure(*rhsStructure)); 
  }  
  
  // copy vector of logical volumes
  fLogicalVolumes.clear();
  for (i=0; i<right.fLogicalVolumes.entries(); i++) {
    G4LogicalVolume* rhsLV = right.fLogicalVolumes[i];
    fLogicalVolumes.insert(rhsLV); 
  }  
  
  fPathName = right.fPathName;
  fDirName = right.fPathName;
  fVerboseLevel = right.fVerboseLevel;
}

AliLVStructure::AliLVStructure() {
//
}

AliLVStructure::~AliLVStructure() {
//
  fStructures.clearAndDestroy();
  fLogicalVolumes.clear();
}

// operators

AliLVStructure& AliLVStructure::operator=(const AliLVStructure &right)
{
  // check assignement to self
  if (this == &right) return *this;

  // copy vector of structures
  fStructures.clearAndDestroy();
  G4int i;
  for (i=0; i<right.fStructures.entries(); i++) {
    // new full structure tree has to be created
    AliLVStructure* rhsStructure = right.fStructures[i];
    fStructures.insert(new AliLVStructure(*rhsStructure)); 
  }  
  
  // copy vector of logical volumes
  fLogicalVolumes.clear();
  for (i=0; i<right.fLogicalVolumes.entries(); i++) {
    G4LogicalVolume* rhsLV = right.fLogicalVolumes[i];
    fLogicalVolumes.insert(rhsLV); 
  }  
  
  fPathName = right.fPathName;
  fDirName = right.fPathName;
  fVerboseLevel = right.fVerboseLevel;

  return *this;
}

G4int AliLVStructure::operator==(const AliLVStructure &right) const
{
  // check == to self
  if (this == &right) return true;

  return false;
}

// private methods

AliLVStructure* AliLVStructure::FindSubDirectory(G4String subDir)
{
// Finds the subdirectory.
// ---

  for( G4int i=0; i<fStructures.entries(); i++ ) {
    if (subDir == fStructures(i)->fDirName) return fStructures(i);
  } 
  return 0;
}

G4String AliLVStructure::ExtractDirName(G4String name)
{
// Extracts the directory name from the path.
// ---

  G4String subDir = name;
  G4int i = name.first('/');
  if (i != G4std::string::npos) subDir.remove(i+1);
  return subDir;
}

// public methods

void AliLVStructure::AddNewVolume(G4LogicalVolume* lv, 
                      G4String treeStructure)
{
// Adds new logical volume to the structure.
// ---

  G4String remainingPath = treeStructure;
  remainingPath.remove(0, fPathName.length());  
  if (!remainingPath.isNull()) { 
    // The lv should be kept in subdirectory.
    // First, check if the subdirectoy exists.
    G4String subDir = ExtractDirName( remainingPath );
    AliLVStructure* targetLVS = FindSubDirectory(subDir);
    if (targetLVS == 0) { 
      // Subdirectory not found. Create a new directory.
      subDir.prepend(fPathName);
      targetLVS = new AliLVStructure(subDir);
      fStructures.insert( targetLVS );
    }
    targetLVS->AddNewVolume(lv, treeStructure);
  }
  else { 
    // the logical volumes should be kept in this directory.
    G4LogicalVolume* targetLV = GetVolume(lv->GetName());
    if (targetLV != 0) {
      // G4cout << lv->GetName() << " had already stored in "
      //        << fPathName << G4endl;
    }
    else {
      fLogicalVolumes.insert(lv);
    }
  }
}

G4LogicalVolume* AliLVStructure::GetVolume(G4String lvName)
{
// Returns logical volume of lvName if present in the structure,
// returns 0 otherwise.
// ---

  for (G4int i=0; i<fLogicalVolumes.entries(); i++) {
    G4LogicalVolume* targetLV = fLogicalVolumes(i);
    if (lvName == targetLV->GetName()) return targetLV;
  }
  return 0;
}

G4LogicalVolume* AliLVStructure::FindVolume(G4String name)
{
// Finds logical volume of given name in all structure tree.
// ---

  G4String path = name;
  path.remove(0, fPathName.length());
  if (path.first('/') != G4std::string::npos) { 
    // SD exists in sub-directory
    G4String subDir = ExtractDirName(path);
    AliLVStructure* targetLVS = FindSubDirectory(subDir);
    if (targetLVS == 0) {  
      // The subdirectory is not found
      G4String text = subDir + " is not found in " + fPathName;
      AliGlobals:: Warning(text);
      return 0;
    }
    else { 
      return targetLVS->FindVolume(name); 
    }
  }
  else { 
    // LV must exist in this directory
    G4LogicalVolume* targetLV = GetVolume(path);
    if (targetLV == 0) {  
      // The fLogicalVolumes is not found.
      G4String text = path + " is not found in " + fPathName;
      AliGlobals::Warning(text);
    }
    return targetLV;
  }
}

void AliLVStructure::ListTree() const
{
// Prints LV tree structure.
// ---

  for (G4int i=0; i<fLogicalVolumes.entries(); i++) {
    G4LogicalVolume* lv = fLogicalVolumes(i);
    G4cout << fPathName << lv->GetName() << G4endl;
  }
  for (G4int j=0; j<fStructures.entries(); j++) { 
    fStructures(j)->ListTree(); 
  }
}
        
void AliLVStructure::ListTreeLong() const
{
// Prints LV tree structure with number of
// daughters (physical volume)
// ---

  for (G4int i=0; i<fLogicalVolumes.entries(); i++) {
    G4LogicalVolume* lv = fLogicalVolumes(i);
    G4cout << fPathName << lv->GetName() 
           << " (" << lv->GetNoDaughters() << ")" << G4endl;
  }
  for (G4int j=0; j<fStructures.entries(); j++) { 
    fStructures(j)->ListTreeLong(); 
  }
}
        
void AliLVStructure::SetVerboseLevel(G4int verbose) 
{
// Sets verbose level.
// ---

  fVerboseLevel = verbose;  
  for (G4int i=0; i<fStructures.entries(); i++) { 
    fStructures(i)->SetVerboseLevel(verbose); 
  }
}

#ifdef ALICE_VISUALIZE
void AliLVStructure::SetTreeVisibility(G4bool visibility)       
{
// Sets visibility to all logical volumes in the structure 
// tree.
// ---

  for (G4int i=0; i<fLogicalVolumes.entries(); i++) {
    G4LogicalVolume* lv = fLogicalVolumes(i);

    const G4VisAttributes* kpVisAttributes = lv->GetVisAttributes();
    G4VisAttributes* newVisAttributes; 
    if (kpVisAttributes) {
      G4Colour colour   = kpVisAttributes->GetColour();
      newVisAttributes = new G4VisAttributes(colour); 
    }
    else
      newVisAttributes = new G4VisAttributes();
    delete kpVisAttributes;

    newVisAttributes->SetVisibility(visibility); 

    lv->SetVisAttributes(newVisAttributes);
  }
  for (G4int j=0; j<fStructures.entries(); j++) { 
    fStructures(j)->SetTreeVisibility(visibility); 
  }
}

void AliLVStructure::SetTreeColour(G4String colName)
{
// Sets colour specified  by name to all logical volumes
// in the structure tree.
// ---

  for (G4int i=0; i<fLogicalVolumes.entries(); i++) {
    G4LogicalVolume* lv = fLogicalVolumes(i);

    const G4VisAttributes* kpVisAttributes = lv->GetVisAttributes ();
    G4VisAttributes* newVisAttributes; 
    if (kpVisAttributes) {
      G4bool oldVisibility = kpVisAttributes->IsVisible();
      newVisAttributes = new G4VisAttributes(oldVisibility); 
    }
    else
      newVisAttributes = new G4VisAttributes();
    delete kpVisAttributes;

    AliColourStore* pColours = AliColourStore::Instance();
    G4Colour colour = pColours->GetColour(colName);
    newVisAttributes->SetColour(colour);

    lv->SetVisAttributes(newVisAttributes);
  }
  for (G4int j=0; j<fStructures.entries(); j++) { 
    fStructures(j)->SetTreeColour(colName); 
  }
}
#endif             


