// $Id$
// Category: geometry
// by I. Hrivnacova, 27.07.2000 
//
// See the class description in the header file.

#include "TG4XMLGeometryGenerator.h"
#include "TG4XMLConvertor.h"
#include "TG4Globals.h"

#include <G4Material.hh>
#include <G4VSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>

#include <g4std/iomanip>
#include <g4std/vector>

TG4XMLGeometryGenerator* TG4XMLGeometryGenerator::fgInstance = 0;

TG4XMLGeometryGenerator::TG4XMLGeometryGenerator() 
{
//
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4XMLGeometryGenerator: attempt to create two instances of singleton.");
  }

  fConvertor = new TG4XMLConvertor(fOutFile);  
}

TG4XMLGeometryGenerator::TG4XMLGeometryGenerator(const TG4XMLGeometryGenerator& right) 
{
// 
  TG4Globals::Exception(
    "TG4XMLGeometryGenerator: attempt to copy singleton.");
}


TG4XMLGeometryGenerator::~TG4XMLGeometryGenerator() {
//
}

// operators

TG4XMLGeometryGenerator& 
TG4XMLGeometryGenerator::operator=(const TG4XMLGeometryGenerator& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4XMLGeometryGenerator singleton.");
    
  return *this;  
}    
          

// private methods

void TG4XMLGeometryGenerator::CutName(G4String& name) const
{
// Removes spaces after the name if present.
// ---

  G4int i = name.length();
  while (name(--i) == ' ') name = name(0,i);
}  

void TG4XMLGeometryGenerator::ProcessSolids(G4LogicalVolume* lv) 
{
// Writes all solids of given logical volume.
// ---

  G4VSolid* solid = lv->GetSolid();
  G4String lvName = lv->GetName();
  G4String material = lv->GetMaterial()->GetName();
  fConvertor->WriteSolid(lvName, solid, material);
  
  // store the name of logical volume in the set
  fVolumeNames.insert(fVolumeNames.begin(), lvName); 

  // process daughters
  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters>0) 
    for (G4int i=0; i<nofDaughters; i++) {
      //G4cout << "processing " << i << "th daughter of " 
      //       << lv->GetName() << G4endl;
      G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume();
      G4String lvdName = lvd->GetName();

      if (fVolumeNames.find(lvdName) == fVolumeNames.end()) {
        // process lvd only if it was not yet processed
        ProcessSolids(lvd);
      }	
    }
}  

void TG4XMLGeometryGenerator::ProcessMaterials(G4LogicalVolume* lv) 
{
// Writes all materials of given logical volume.
// ---

  G4Material* material = lv->GetMaterial();
  
  // check if this material was already written
  G4bool written = false;
  G4String name = material->GetName();
  CutName(name);
  if (fMaterialNames.find(name) != fMaterialNames.end()) written = true;
  
  if (!written) {
    fConvertor->WriteMaterial(material);
    fMaterialNames.insert(fMaterialNames.begin(), name); 
  }  
  
  // store the name of logical volume in the set
  G4String lvName = lv->GetName();
  fVolumeNames.insert(fVolumeNames.begin(), lvName); 

  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters>0) 
    for (G4int i=0; i<nofDaughters; i++) {
      G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume();
      G4String lvdName = lvd->GetName();

      if (fVolumeNames.find(lvdName) == fVolumeNames.end()) {
        // process lvd only if it was not yet processed
        ProcessMaterials(lvd);
      }	
    }
}  

void TG4XMLGeometryGenerator::ProcessRotations(G4LogicalVolume* lv) 
{
// Writes all rotation matrices of given logical volume.
// ---

  G4String lvName = lv->GetName();

  // store the name of logical volume in the set
  fVolumeNames.insert(fVolumeNames.begin(), lvName); 
  
  G4int nofDaughters = lv->GetNoDaughters();

  if (nofDaughters>0) {
    G4int i; 
    for (i=0; i<nofDaughters; i++) {
      
      G4VPhysicalVolume* pvd = lv->GetDaughter(i);
      const G4RotationMatrix* kRotation = pvd->GetRotation();
      if (kRotation) 
        fConvertor->WriteRotation(kRotation);

      G4LogicalVolume* lvd = pvd->GetLogicalVolume();
      G4String lvdName = lvd->GetName();

      if (fVolumeNames.find(lvdName) == fVolumeNames.end()) {
        // process lvd only if it was not yet processed
        ProcessRotations(lvd);
      }	
    }
  }  
}  

void TG4XMLGeometryGenerator::ProcessLogicalVolume(G4LogicalVolume* lv) 
{
// Writes logical volume tree.
// ---
  
  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters == 0) return;
  
  // open composition
  G4String lvName = lv->GetName();
  G4String name = lvName;
  name.append("_comp");
  fConvertor->OpenComposition(name);
  
  // write positions  
  G4int i;
  for (i=0; i<nofDaughters; i++) {
   // G4cout << "processing " << i << "th daughter of " 
   //        << lv->GetName() << G4endl;
   
    G4VPhysicalVolume* vpvd = lv->GetDaughter(i);
    G4LogicalVolume* lvd = vpvd->GetLogicalVolume();
      
    // get parameters
    G4String lvName = lvd->GetName();
    G4String compName = lvd->GetName();
    compName.append("_comp");      
    G4int nd = lvd->GetNoDaughters(); 

    G4PVPlacement* pvd = dynamic_cast<G4PVPlacement*>(vpvd);
    if (pvd) {
      // placement
      G4ThreeVector  position = vpvd->GetTranslation();
      const G4RotationMatrix* kMatrix = vpvd->GetObjectRotation();      

      if (!kMatrix) {
  	fConvertor->WritePosition(lvName, position);
        // if volume is not leaf node place its logical volume
        if (nd>0) 
    	  fConvertor->WritePosition(compName, position);
      }	  
      else {  
  	fConvertor->WritePositionWithRotation(lvName, position, kMatrix);
        if (nd>0) 
      	   fConvertor->WritePositionWithRotation(compName, position, kMatrix);
      }
    }
    else {
      G4PVReplica* pvr = dynamic_cast<G4PVReplica*>(vpvd);
      if (pvr) {
        // replica
    	fConvertor->WriteReplica(lvName, pvr);
        // if volume is not leaf node place its logical volume
        if (nd>0) 
      	  fConvertor->WriteReplica(compName, pvr);
      }
      else {
        G4String text = "TG4XMLGeometryGenerator::ProcessLogicalVolume: \n";
        text = text + "    Limitation: \n";
        text = text + "    Other physical volumes than PVPlacement and PVReplica";
        text = text + " are not implemented.";
        TG4Globals::Exception(text);
      }
    }  
  }  

  // close composition
  fConvertor->CloseComposition();	
  fConvertor->WriteEmptyLine();

  // store the name of logical volume in the set
  fVolumeNames.insert(fVolumeNames.begin(), lvName); 

  // process daughters
  for (i=0; i<nofDaughters; i++) {
    G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume();
    G4String lvdName = lvd->GetName();

    if (fVolumeNames.find(lvdName) == fVolumeNames.end()) {
      // process lvd only if it was not yet processed
      ProcessLogicalVolume(lvd);
    }
  }    
}  

void TG4XMLGeometryGenerator::ClearMaterialNames() 
{
// Clears the set of material names.
// ---

  fMaterialNames.erase(fMaterialNames.begin(), fMaterialNames.end());
}  

void TG4XMLGeometryGenerator::ClearVolumeNames() 
{
// Clears the set of volume names.
// ---

  fVolumeNames.erase(fVolumeNames.begin(), fVolumeNames.end());
}  

// public methods

void TG4XMLGeometryGenerator::GenerateMaterials( 
                        const G4String& version, const G4String& date,
		        const G4String& author,  const G4String dtdVersion,
			G4LogicalVolume* lv)
{
// Generates the XML material element containing
// all materials present in given logical volume.
// ---

  // create section
  fConvertor->OpenMaterials(version, date, author, dtdVersion);  
  fConvertor->WriteEmptyLine();
  
  // process materials
  ProcessMaterials(lv);
  fConvertor->WriteEmptyLine();
  ClearMaterialNames();
  ClearVolumeNames();

  // close section
  fConvertor->CloseMaterials();
  fConvertor->WriteEmptyLine();
}   

void TG4XMLGeometryGenerator::GenerateSection(const G4String& name, 
                        const G4String& version, const G4String& date,
		        const G4String& author, const G4String& topVolume,
                        G4LogicalVolume* lv)
{
// Generates the XML section element containing
// all geometry objects defined in given logical volume:
// rotation matrices, solids and volumes hierarchy.
// ---

  // create section
  fConvertor->OpenSection(name, version, date, author, topVolume);  
  fConvertor->WriteEmptyLine();
  
  // process rotations
  //ProcessRotations(lv);
  //fConvertor->WriteEmptyLine();
  //ClearRotations();
  //ClearVolumeNames();
    
  // process solids
  ProcessSolids(lv);
  fConvertor->WriteEmptyLine();
  ClearVolumeNames();
    
  // process geometry tree
  ProcessLogicalVolume(lv);
  fConvertor->WriteEmptyLine();
  ClearVolumeNames();
  
  // close section
  fConvertor->CloseSection();
}   

void TG4XMLGeometryGenerator::OpenFile(G4String filePath)
{ 
// Opens output file.
// ---

  fOutFile.open(filePath, G4std::ios::out); 
  
  if (!fOutFile) {
    G4String text = "Cannot open ";
    text = text + filePath;
    TG4Globals::Warning(text);  
  }
  
  // use FORTRAN compatibility output
  fOutFile.setf(G4std::ios::fixed, G4std::ios::floatfield);
}


void TG4XMLGeometryGenerator::CloseFile()
{ 
// Closes output file.
// ---

  fOutFile.close(); 
}
