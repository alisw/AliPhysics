// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliLVStructure
// --------------------
// See the class description in the header file.

#include "AliLVStructure.h"
#include "AliGlobals.h"

#ifdef G4VIS_USE
#include "AliColourStore.h"

#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#endif //G4VIS_USE
#include <G4LogicalVolume.hh>
#include <G4BooleanSolid.hh>
#include <g4std/set>

//_____________________________________________________________________________
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

//_____________________________________________________________________________
AliLVStructure::AliLVStructure(const AliLVStructure& right)
{
  // copy stuff
  *this = right;
}

//_____________________________________________________________________________
AliLVStructure::AliLVStructure() {
//
}

//_____________________________________________________________________________
AliLVStructure::~AliLVStructure() {
//


  ClearAndDestroy(&fStructures);
  
  
  fLogicalVolumes.resize(0);
}

// operators

//_____________________________________________________________________________
AliLVStructure& AliLVStructure::operator=(const AliLVStructure &right)
{
  // check assignement to self
  if (this == &right) return *this;

  // copy vector of structures
  ClearAndDestroy(&fStructures);
  G4int i;
  for (i=0; i<right.fStructures.size(); i++) {
    // new full structure tree has to be created
    AliLVStructure* rhsStructure = right.fStructures[i];
    fStructures.push_back(new AliLVStructure(*rhsStructure)); 
  }  
  
  // copy vector of logical volumes
  fLogicalVolumes.resize(0);
  for (i=0; i<right.fLogicalVolumes.size(); i++) {
    G4LogicalVolume* rhsLV = right.fLogicalVolumes[i];
    fLogicalVolumes.push_back(rhsLV); 
  }  
  
  fPathName = right.fPathName;
  fDirName = right.fPathName;
  fVerboseLevel = right.fVerboseLevel;

  return *this;
}

//_____________________________________________________________________________
G4int AliLVStructure::operator==(const AliLVStructure &right) const
{
  // check == to self
  if (this == &right) return true;

  return false;
}

// private methods

//_____________________________________________________________________________
AliLVStructure* AliLVStructure::FindSubDirectory(const G4String& subDir) const
{
// Finds the subdirectory.
// ---

  for( G4int i=0; i<fStructures.size(); i++ ) {
    if (subDir == fStructures[i]->fDirName) return fStructures[i];
  } 
  return 0;
}

//_____________________________________________________________________________
G4String AliLVStructure::ExtractDirName(const G4String& name) const
{
// Extracts the directory name from the path.
// ---

  G4String subDir = name;
  G4int i = name.first('/');
  if (i != G4std::string::npos) subDir.remove(i+1);
  return subDir;
}

void AliLVStructure::ClearAndDestroy(LVStructuresVector* structures) 
{
// Clears the structures vectors and deletes all contained 
// elements.
// (According to geant4.3.2/source/global/STLInterface/g4rw/tpordvec.icc.)
// ---

  G4std::set<AliLVStructure*,G4std::greater<AliLVStructure*> > tmp;
  for (size_t sz=0; sz<structures->size(); sz++)
    {
      AliLVStructure* current;
      current=(*structures)[sz];
      if (current)
         tmp.insert(current);
    }
    
  G4std::set<AliLVStructure*, G4std::greater<AliLVStructure*> >::iterator it;
  for (it=tmp.begin(); it!=tmp.end(); it++)
    {
      delete *it;
    }
  // std_pvector::erase(std_pvector::begin(), std_pvector::end());
  structures->resize(0);
}

// public methods

//_____________________________________________________________________________
void AliLVStructure::AddNewVolume(G4LogicalVolume* lv, 
                                  const G4String& treeStructure)
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
      fStructures.push_back( targetLVS );
    }
    targetLVS->AddNewVolume(lv, treeStructure);
  }
  else { 
    // the logical volumes should be kept in this directory.
    G4LogicalVolume* targetLV = GetVolume(lv->GetName());
    if (targetLV != 0) {
      // G4cout << lv->GetName() << " has been already stored in "
      //        << fPathName << G4endl;
    }
    else {
      fLogicalVolumes.push_back(lv);
    }
  }
}

//_____________________________________________________________________________
G4LogicalVolume* AliLVStructure::GetVolume(const G4String& lvName) const
{
// Returns logical volume of lvName if present in the structure,
// returns 0 otherwise.
// ---

  for (G4int i=0; i<fLogicalVolumes.size(); i++) {
    G4LogicalVolume* targetLV = fLogicalVolumes[i];
    if (lvName == targetLV->GetName()) return targetLV;
  }
  return 0;
}

//_____________________________________________________________________________
G4LogicalVolume* AliLVStructure::FindVolume(const G4String& name) const
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

//_____________________________________________________________________________
void AliLVStructure::ListTree() const
{
// Prints LV tree structure.
// ---

  for (G4int i=0; i<fLogicalVolumes.size(); i++) {
    G4LogicalVolume* lv = fLogicalVolumes[i];
    G4cout << fPathName << lv->GetName() << G4endl;
  }
  for (G4int j=0; j<fStructures.size(); j++) { 
    fStructures[j]->ListTree(); 
  }
}
        
//_____________________________________________________________________________
void AliLVStructure::ListTreeLong() const
{
// Prints LV tree structure with number of
// daughters (physical volume), indicates Boolean solid.
// ---

  for (G4int i=0; i<fLogicalVolumes.size(); i++) {
    G4LogicalVolume* lv = fLogicalVolumes[i];

    G4cout << fPathName << lv->GetName() << " (" << lv->GetNoDaughters();
	    
    if (dynamic_cast<G4BooleanSolid*>(lv->GetSolid()))
      G4cout << ", B";

    G4cout << ")" << G4endl;
  }
  for (G4int j=0; j<fStructures.size(); j++) { 
    fStructures[j]->ListTreeLong(); 
  }
}
        
//_____________________________________________________________________________
void AliLVStructure::SetVerboseLevel(G4int verbose) 
{
// Sets verbose level.
// ---

  fVerboseLevel = verbose;  
  for (G4int i=0; i<fStructures.size(); i++) { 
    fStructures[i]->SetVerboseLevel(verbose); 
  }
}

#ifdef G4VIS_USE
//_____________________________________________________________________________
void AliLVStructure::SetTreeVisibility(G4bool visibility)       
{
// Sets visibility to all logical volumes in the structure 
// tree.
// ---

  for (G4int i=0; i<fLogicalVolumes.size(); i++) {
    G4LogicalVolume* lv = fLogicalVolumes[i];

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
  for (G4int j=0; j<fStructures.size(); j++) { 
    fStructures[j]->SetTreeVisibility(visibility); 
  }
}

//_____________________________________________________________________________
void AliLVStructure::SetTreeColour(const G4String& colName)
{
// Sets colour specified  by name to all logical volumes
// in the structure tree.
// ---

  for (G4int i=0; i<fLogicalVolumes.size(); i++) {
    G4LogicalVolume* lv = fLogicalVolumes[i];

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
  for (G4int j=0; j<fStructures.size(); j++) { 
    fStructures[j]->SetTreeColour(colName); 
  }
}
#endif             


