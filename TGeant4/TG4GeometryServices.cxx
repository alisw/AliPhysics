// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "TG4GeometryServices.h"
#include "TG4VSensitiveDetector.h"
#include "TG4Globals.h"
#include "TG3Units.h"

#include <G4VSensitiveDetector.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G3toG4.hh> 

TG4GeometryServices* TG4GeometryServices::fgInstance = 0;

TG4GeometryServices::TG4GeometryServices(TG4intVector* mediumIdVector, 
                                         TG4NameMap* nameMap) 
{
//
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4GeometryServices: attempt to create two instances of singleton.");
  }

  fMediumIdVector = mediumIdVector;
  fNameMap = nameMap;

  fgInstance = this;
}

TG4GeometryServices::TG4GeometryServices() {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4GeometryServices singleton.");
}


TG4GeometryServices::TG4GeometryServices(const TG4GeometryServices& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4GeometryServices singleton.");
}


TG4GeometryServices::~TG4GeometryServices() {
//
}

//=============================================================================
//
// operators
//
//=============================================================================

TG4GeometryServices& 
TG4GeometryServices::operator=(const TG4GeometryServices& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4GeometryServices singleton.");
    
  return *this;  
}    
          

//=============================================================================
//
// public methods
//
//=============================================================================


G4double* TG4GeometryServices::CreateG4doubleArray(Float_t* array, 
               G4int size) const
{
// Converts Float_t* array to G4double*,
// !! The new array has to be deleted by user.
// ---

  G4double* doubleArray;
  if (size>0) {
    doubleArray = new G4double[size]; 
    for (G4int i=0; i<size; i++) doubleArray[i] = array[i];
  }
  else {
    doubleArray = 0; 
  }  
  return doubleArray;
}


G4String TG4GeometryServices::CutName(const char* name) const
{
// Removes spaces after the name if present.
// ---

  G4String cutName = name;
  G4int i = cutName.length();
  while (cutName(--i) == ' ') cutName = cutName(0,i);

  return cutName;
}  


void TG4GeometryServices::G4ToG3VolumeName(G4String& name) const
{
// Cuts _copyNo extension added to logical volume name in case 
// the logical volume was created by Gsposp method.
// ---

  if (name.contains(gSeparator)) 
  name = name(0,name.first(gSeparator));
}


G4Material* TG4GeometryServices::MixMaterials(G4String name, G4double density, 
        TG4StringVector* matNames, TG4doubleVector* matWeights)
{
// Creates a mixture of selected materials
// ---

  // number of materials to be mixed  
  G4int nofMaterials = matNames->entries();
  if (nofMaterials != matWeights->entries()) {
    G4String text = "TG4GeometryServices::MixMaterials: ";
    text = text +  "different number of material names and weigths.";
    TG4Globals::Exception(text);
  }    
  // add verbose
  // G4cout << "Nof of materials to be mixed: " << nofMaterials << G4endl;

  // fill vector of materials
  TG4MaterialVector matVector;  
  G4int im;
  for (im=0; im< nofMaterials; im++) {
    // material
    G4Material* material = G4Material::GetMaterial((*matNames)[im]);
    matVector.insert(material);
  } 

  // create the mixed material
  G4Material* mixture = new G4Material(name, density, nofMaterials);
  for (im=0; im< nofMaterials; im++) {
    G4Material* material = matVector[im];
    G4double fraction = (*matWeights)[im];
    mixture->AddMaterial(material, fraction);
  }

  return mixture;
}  


 Int_t TG4GeometryServices::NofG3Volumes() const
{
// Returns the total number of logical volumes corresponding
// to G3 volumes. (
// The logical volume that were created by Gsposp method 
// with a generic name (name_copyNo) are NOT included.
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();

  G4int counter = 0;  
  for (G4int i=0; i<pLVStore->entries(); i++) {
    G4LogicalVolume* lv = (*pLVStore)[i];
    if (IsG3Volume(lv->GetName())) counter++;
  }
  
  return counter;  
}


Int_t TG4GeometryServices::NofG4LogicalVolumes() const
{
// Returns the total number of logical volumes in the geometry.
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  return pLVStore->entries();
}


Int_t TG4GeometryServices::NofG4PhysicalVolumes() const
{
// Returns the total number of physical volumes in the geometry.
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();

  G4int counter = 0;
  for (G4int i=0; i<pLVStore->entries(); i++) {
    counter += ((*pLVStore)[i])->GetNoDaughters();
  }
  
  return counter;  
}


Int_t TG4GeometryServices::NofSensitiveDetectors() const
{
// Returns the total number of sensitive detectors.
// ---

  return TG4VSensitiveDetector::GetTotalNofSensitiveDetectors();
}

 
G4int TG4GeometryServices::GetVolumeID(const G4String& volName) const
{ 
// Returns the sensitive detector identifier.
// !! Gives exception in case logical volume is not associated with 
// a sensitive detector.
// ---

  G4String g4VolName = CutName(volName);
  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  
  for (G4int i=0; i<pLVStore->entries(); i++) {
    G4LogicalVolume* lv = pLVStore->at(i);
    G4VSensitiveDetector* sd = lv->GetSensitiveDetector();
  
    if ((sd) && (sd->GetName()==g4VolName)) {
      TG4VSensitiveDetector* tsd = dynamic_cast<TG4VSensitiveDetector*>(sd);
      if (tsd)
        return tsd->GetID();
      else {
        TG4Globals::Exception(
          "TG4GeometryServices::GetVolumeID: Unknown sensitive detector type");
        return 0;
      }
    }   
  }

  G4String text = "TG4GeometryServices::VolId: Sensitive detector ";
  text = text + g4VolName;
  text = text + " is not defined.\n"; 
  TG4Globals::Warning(text);
  return 0;
}


G4int TG4GeometryServices::GetVolumeID(G4LogicalVolume* logicalVolume) const 
{
// Returns the sensitive detector ID of the specified
// logical volume.
// ---
 
  // sensitive detector ID
  G4VSensitiveDetector* sd
    = logicalVolume->GetSensitiveDetector();
  if (sd) {
    TG4VSensitiveDetector* tsd = dynamic_cast<TG4VSensitiveDetector*>(sd);
    if (tsd)
      return tsd->GetID();
    else {
      TG4Globals::Exception(
        "TG4GeometryServices::GetVolumeID: Unknown sensitive detector type");
      return 0;
    }   	
  }  
  else {
    G4String text = "TG4GeometryServices::GetVolumeID: \n";
    text = text + "    Volume " + logicalVolume->GetName();
    text = text + " has not a sensitive detector.";
    //TG4Globals::Exception(text);
    TG4Globals::Warning(text);
    return 0;
  }      	
} 


G4String TG4GeometryServices::GetVolumeName(G4int volumeId) const
{
// Returns the name of the sensitive detector with the given identifier.
// Gives a warning in case logical volume is not associated with 
// a sensitive detector.
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  
  for (G4int i=0; i<pLVStore->entries(); i++) {
    G4LogicalVolume* lv = pLVStore->at(i);
    G4VSensitiveDetector* sd = lv->GetSensitiveDetector();
    
    if (sd) {
      G4int sdID;
      TG4VSensitiveDetector* tsd = dynamic_cast<TG4VSensitiveDetector*>(sd);
      if (tsd)
        sdID = tsd->GetID();
      else {
        TG4Globals::Exception(
          "TG4GeometryServices::VolId: Unknown sensitive detector type");
        return "";
      }
      if (sdID == volumeId) return sd->GetName();
    }  
  }

  G4String text = "TG4GeometryServices::VolName:\n";
  text = text + "    Sensitive detector with given id is not defined. \n";
  TG4Globals::Warning(text);
  return "";	       	         
}


G4LogicalVolume* TG4GeometryServices::GetLogicalVolume(G4int volumeId) const 
{
// Finds the first logical volume with specified volumeId 
// (sensitive detector ID) in G4LogicalVolumeStore.
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  
  for (G4int i=0; i<pLVStore->entries(); i++) {
    G4LogicalVolume* lv = pLVStore->at(i);
    if (GetVolumeID(lv) == volumeId) return lv;
  }
  
  G4String text = "TG4GeometryServices::GetLogicalVolume: \n";
  text = text + "    Logical volume with given ID does not exist.";
  return 0;	       	         
}  


G4bool TG4GeometryServices::IsG3Volume(const G4String& lvName) const
{
// Returns true if the logical volume of given volumeName
// was not created by Gsposp method with a generic name 
// (name_copyNo).
// ---

  if (lvName.contains(gSeparator))
    return false;  
  else
    return true;   
}


const G4String& TG4GeometryServices::GetMapSecond(const G4String& name)
{ 
// Returns the second string associated with the name in
// the name map.
// ---

  return fNameMap->GetSecond(name); 
}


G4int TG4GeometryServices::GetMediumId(G4Material* material) const
{
// Returns the second index for materials (having its origin in
// G4 tracking media concept)
// ---

  return (*fMediumIdVector)[material->GetIndex()];
}  


Int_t TG4GeometryServices::GetMediumId(G4int volumeId)  const
{
// Return the material number for a given volume id
// ---

  G4LogicalVolume* logicalVolume = GetLogicalVolume(volumeId);
    
  G4Material* material = logicalVolume->GetMaterial();  

  return GetMediumId(material);	       	         
}
 

G4double TG4GeometryServices::GetEffA(G4Material* material) const
{
// Returns A or effective A=sum(pi*Ai) (if compound/mixture)
// of given material.
// ---

  G4double a = 0.;
  G4int nofElements = material->GetNumberOfElements();
  if (nofElements > 1) {
    G4String text = "Effective A for material mixture (";
    text = text + material->GetName();
    text = text + ") is used.";
    //TG4Globals::Warning(text);

    for (G4int i=0; i<nofElements; i++) {
      G4double aOfElement = material->GetElement(i)->GetA();
      G4double massFraction = material->GetFractionVector()[i];      
      a += aOfElement*massFraction /(TG3Units::AtomicWeight());
    }
  }
  else { 
    a = material->GetA();
    a /= TG3Units::AtomicWeight();
  }
  return a;
}


G4double TG4GeometryServices::GetEffZ(G4Material* material) const
{
// Returns Z or effective Z=sum(pi*Zi) (if compound/mixture)
// of given material.
// ---

  G4double z = 0.;
  G4int nofElements = material->GetNumberOfElements();
  if (nofElements > 1) {
    G4String text = "Effective Z for material mixture (";
    text = text + material->GetName();
    text = text + ") is used.";
    //TG4Globals::Warning(text);

    for (G4int i=0; i<nofElements; i++) {
      G4double zOfElement = material->GetElement(i)->GetZ();
      G4double massFraction = material->GetFractionVector()[i];
      z += zOfElement*massFraction;
    }
  }
  else { 
    z = material->GetZ(); 
  }  
  return z;
}
