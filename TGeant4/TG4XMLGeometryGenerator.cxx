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

void TG4XMLGeometryGenerator::ProcessSolids(G4LogicalVolume* lv) 
{
// Writes all solids of given logical volume.
// ---

  G4VSolid* solid = lv->GetSolid();
  G4String material = lv->GetMaterial()->GetName();
  fConvertor->WriteSolid(solid, material);
  
  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters>0) 
    for (G4int i=0; i<nofDaughters; i++) {
      G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume();
      ProcessSolids(lvd);
    }
}  

void TG4XMLGeometryGenerator::ProcessMaterials(G4LogicalVolume* lv) 
{
// Writes all materials of given logical volume.
// ---

  G4Material* material = lv->GetMaterial();
  fConvertor->WriteMaterial(material);
  
  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters>0) 
    for (G4int i=0; i<nofDaughters; i++) {
      G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume();
      ProcessMaterials(lvd);
    }
}  

void TG4XMLGeometryGenerator::ProcessRotations(G4LogicalVolume* lv) 
{
// Writes all rotation matrices of given logical volume.
// ---

  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters>0) 
    for (G4int i=0; i<nofDaughters; i++) {
      G4VPhysicalVolume* pvd = lv->GetDaughter(i);
      const G4RotationMatrix* kRotation = pvd->GetRotation();
      if (kRotation) fConvertor->WriteRotation(kRotation);
      G4LogicalVolume* lvd = pvd->GetLogicalVolume();
      ProcessRotations(lvd);
    }
}  

void TG4XMLGeometryGenerator::ProcessLogicalVolume(G4LogicalVolume* lv) 
{
// Writes logical volume tree.
// ---
  
  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters == 0) return;
  
  // open composition
  G4String name = lv->GetName();
  name.append("_comp");
  fConvertor->OpenComposition(name);
      
  // write positions  
  for (G4int i=0; i<nofDaughters; i++) {
    G4VPhysicalVolume* vpvd = lv->GetDaughter(i);
    G4LogicalVolume* lvd = vpvd->GetLogicalVolume();
      
    // only placements are processed
    G4PVPlacement* pvd = dynamic_cast<G4PVPlacement*>(vpvd);
    if (pvd) {
      G4String solidName = lvd->GetSolid()->GetName();
      G4String compName = lvd->GetName();
      compName.append("_comp");      
      G4int nd = lvd->GetNoDaughters(); 
      G4ThreeVector  position = vpvd->GetFrameTranslation();
      const G4RotationMatrix* kMatrix = vpvd->GetFrameRotation();      

      if (!kMatrix) {
  	fConvertor->WritePosition(solidName, position);
        // if volume is not leaf node place its logical volume
        if (nd>0) 
    	  fConvertor->WritePosition(compName, position);
      }	  
      else {  
  	fConvertor->WritePositionWithRotation(solidName, position, kMatrix);
        if (nd>0) 
      	   fConvertor->WritePositionWithRotation(compName, position, kMatrix);
      }   
	   
    }
    else {
      G4String text = "TG4XMLGeometryGenerator::ProcessLogicalVolume: \n";
      text = text + "    Limitation: \n";
      text = text + "    Other physical volumes than PVPlacement";
      text = text + " are not implemented.";
      TG4Globals::Warning(text);
    }
  }  

  // close composition
  fConvertor->CloseComposition();	
  fConvertor->WriteEmptyLine();

  // make a vector of contained logical volumes
  // with daughters
  // -> change to a global map of names of written compositions 
  //    and test against this map 
  G4std::vector<G4LogicalVolume*> vect;
  for (G4int i=0; i<nofDaughters; i++) {
    G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume();
    G4bool store = true;
    for (G4int j=0; j<vect.size(); j++) 
      if (vect[j] == lvd || lvd->GetNoDaughters()==0) store = false;
    if (store) vect.push_back(lvd);
  }   

  // process contained logical volumes with daughters  
  fConvertor->IncreaseIndention();
  for (G4int j=0; j<vect.size(); j++) 
    ProcessLogicalVolume(vect[j]);
  fConvertor->DecreaseIndention();
}  

/*
void TG4XMLGeometryGenerator::ProcessLogicalVolumeOld(G4LogicalVolume* lv) 
{
// Writes logical volume tree.
// ---
  
  G4int nofDaughters = lv->GetNoDaughters();
  if (nofDaughters == 0) return;
  
  // make a vector of contained logical volumes
  G4std::vector<G4LogicalVolume*> vect;
  for (G4int i=0; i<nofDaughters; i++) {
    G4LogicalVolume* lvd = lv->GetDaughter(i)->GetLogicalVolume();
    G4bool store = true;
    for (G4int j=0; j<vect.size(); j++) 
      if (vect[j] == lvd) store = false;
    if (store) vect.push_back(lvd);
  }   

  // loop over contained logical volumes
  for (G4int j=0; j<vect.size(); j++) {
    G4LogicalVolume* lvd = vect[j];
 
    // open composition
    if(lvd->GetNoDaughters()>0) {
      G4String name = lvd->GetName();
      name.append("_lv");
      fConvertor->OpenComposition(name);
    }	
      
    // write positions  
    for (G4int i=0; i<nofDaughters; i++) {
      G4VPhysicalVolume* vpvd = lv->GetDaughter(i);
      G4LogicalVolume* lvdi = vpvd->GetLogicalVolume();
      
      if (lvdi == lvd) {
        // only placements are processed
        G4PVPlacement* pvd = dynamic_cast<G4PVPlacement*>(vpvd);
        if (pvd) {
          G4String solidName = lvd->GetSolid()->GetName();
          G4ThreeVector  position = vpvd->GetFrameTranslation();
          const G4RotationMatrix* kMatrix = vpvd->GetFrameRotation();      
	  if (!kMatrix)
  	    fConvertor->WritePosition(solidName, position);
	  else  
  	    fConvertor->WritePositionWithRotation(solidName, position, kMatrix);
        }
        else {
          G4String text = "TG4XMLGeometryGenerator::ProcessLogicalVolume: \n";
          text = text + "    Limitation: \n";
	  text = text + "    Other physical volumes than PVPlacement";
	  text = text + " are not implemented.";
          TG4Globals::Warning(text);
        }
      } 
    }

    if(lvd->GetNoDaughters()>0) {
      // process daughters recursively
      ProcessLogicalVolume(lvd);
      fConvertor->CloseComposition();	
      fConvertor->WriteEmptyLine();
    }
  }
}  
*/

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
  ProcessRotations(lv);
  fConvertor->WriteEmptyLine();
    
  // process solids
  ProcessSolids(lv);
  fConvertor->WriteEmptyLine();
    
  // process geometry tree
  ProcessLogicalVolume(lv);
  fConvertor->WriteEmptyLine();
  
  // close section
  fConvertor->CloseSection();
}   

void TG4XMLGeometryGenerator::OpenFile(G4String filePath)
{ 
// Opens output file.
// ---

  G4cout << "TG4XMLGeometryGenerator::OpenFile: " << filePath << G4endl;
  
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
