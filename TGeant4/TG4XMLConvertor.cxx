// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "TG4XMLConvertor.h"
#include "TG3Units.h"

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VSolid.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Trd.hh>
#include <globals.hh>

#include <g4std/iostream>
#include <g4std/iomanip>

const G4int TG4XMLConvertor::fgkMaxVolumeNameLength   = 10;
const G4int TG4XMLConvertor::fgkMaxMaterialNameLength = 20;

TG4XMLConvertor::TG4XMLConvertor(G4std::ofstream& outFile) 
  : fOutFile(outFile),
    fIndention("")
{
//
}

TG4XMLConvertor::~TG4XMLConvertor() {
//
}

// private methods

void TG4XMLConvertor::CutName(G4String& name) const
{
// Removes spaces after the name if present.
// ---

  G4int i = name.length();
  while (name(--i) == ' ') name = name(0,i);
}  

void TG4XMLConvertor::CutName(G4String& name, G4int size) const
{
// Cuts name to given size.
// ---

  if (name.length() > size) name = name(0, size);
}  

void TG4XMLConvertor::PutName(G4String& element, G4String name, 
                              G4String templ) const
{
// Replaces given template in string element with a give name.
// ---

  CutName(name);
  // make better
  if (templ == "#") 
    CutName(name, fgkMaxVolumeNameLength);
  else if (templ == "!")  
    CutName(name, fgkMaxMaterialNameLength);
  
  element.replace(element.find(templ), name.size(), name);
  element.replace(element.find(templ), 1, "\"");
  while (element.contains(templ)) element.replace(element.find(templ), 1 , " ");
}    
  
void TG4XMLConvertor::WriteBox(const G4Box* box, G4String materialName)
{
// Writes G4box solid.
// ---

  // get parameters
  G4String solidName = box->GetName();
  G4double x = box->GetXHalfLength()/TG3Units::Length();
  G4double y = box->GetYHalfLength()/TG3Units::Length();
  G4double z = box->GetZHalfLength()/TG3Units::Length();

  // compose element string template
  G4String element1 
    = "<box   name=\"###########   material=\"!!!!!!!!!!!!!!!!!!!!!   X_Y_Z=\"";
  G4String element2 = "\" />";
  
  // put solid and material names
  PutName(element1, solidName, "#");
  PutName(element1, materialName, "!");
  
  // write element
  //fOutFile.setf(G4std::ios::fixed, G4std::ios::floatfield);
  fOutFile << element1
           << G4std::setw(7) << G4std::setprecision(2) << x << "  "
           << G4std::setw(7) << G4std::setprecision(2) << y << "  "
           << G4std::setw(7) << G4std::setprecision(2) << z
	   << element2
	   << G4endl;
}
 
void TG4XMLConvertor::WriteTubs(const G4Tubs* tubs, G4String materialName)
{
// Writes G4tubs solid.
// ---

  // get parameters
  G4String solidName = tubs->GetName();
  G4double rmin = tubs->GetInnerRadius()/TG3Units::Length();
  G4double rmax = tubs->GetOuterRadius()/TG3Units::Length();
  G4double hz   = tubs->GetZHalfLength()/TG3Units::Length();
  G4double sphi = tubs->GetStartPhiAngle()/TG3Units::Angle();
  G4double dphi = tubs->GetDeltaPhiAngle()/TG3Units::Angle();

  // compose element string template
  G4String element1 
    = "<tubs  name=\"###########   material=\"!!!!!!!!!!!!!!!!!!!!!   Rio_Z=\"";
  G4String element2 = "\"  profile=\"";
  G4String element3 = "\" />";
  
  // put solid and material names
  PutName(element1, solidName, "#");
  PutName(element1, materialName, "!");
  
  // write element
  fOutFile << element1
           << G4std::setw(7) << G4std::setprecision(2) << rmin << "  "
           << G4std::setw(7) << G4std::setprecision(2) << rmax << "  "
           << G4std::setw(7) << G4std::setprecision(2) << hz
	   << element2
           << G4std::setw(7) << G4std::setprecision(2) << sphi << "  "
           << G4std::setw(7) << G4std::setprecision(2) << dphi
	   << element3
	   << G4endl;
}  


void TG4XMLConvertor::WriteTrd(const G4Trd* trd, G4String materialName)
{
// Writes G4Trd solid.
// ---

  // get parameters
  G4String solidName = trd->GetName();
  G4double x1 = trd->GetXHalfLength1()/TG3Units::Length();
  G4double x2 = trd->GetXHalfLength2()/TG3Units::Length();
  G4double y1 = trd->GetYHalfLength1()/TG3Units::Length();
  G4double y2 = trd->GetYHalfLength2()/TG3Units::Length();
  G4double hz = trd->GetZHalfLength()/TG3Units::Length();

  // compose element string template
  G4String element1 
    = "<trd   name=\"###########   material=\"!!!!!!!!!!!!!!!!!!!!!   Xmp_Ymp_Z=\"";
  G4String element2 = "\" />";
  
  // put solid and material names
  // put solid and material names
  PutName(element1, solidName, "#");
  PutName(element1, materialName, "!");
  
  // write element
  fOutFile << element1
           << G4std::setw(7) << G4std::setprecision(2) << x1 << "  "
           << G4std::setw(7) << G4std::setprecision(2) << x2 << "  "
           << G4std::setw(7) << G4std::setprecision(2) << y1 << "  "
           << G4std::setw(7) << G4std::setprecision(2) << y2 << "  "
           << G4std::setw(7) << G4std::setprecision(2) << hz
	   << element2
	   << G4endl;
}  


// public methods

void TG4XMLConvertor::OpenSection(const G4String& name, const G4String& version,
 	                 const G4String& date, const G4String& author,
                         const G4String& topVolume)
{
// Writes section opening.
// ---
			 
  G4String element1 = "<section name       = \"";
  G4String element2 = "         version    = \"";
  G4String element3 = "         date       = \"";
  G4String element4 = "         author     = \"";
  G4String element5 = "         topVolume  = \"";
  G4String element6 = "  >";
  G4String quota = "\"";   
  
  // write element
  fOutFile << element1 << name    << quota << G4endl
           << element2 << version << quota << G4endl
           << element3 << date    << quota << G4endl
           << element4 << author  << quota << G4endl
           << element5 << topVolume << quota
           << element6 << G4endl;
}  

void TG4XMLConvertor::OpenComposition(const G4String& name)
{
// Writes composition opening.
// ---
			 
  G4String element = "<composition name=\"";
  element.append(name);
  element.append(">");

  // write element
  fOutFile << fIndention
           << element
	   << G4endl;

  // increase indention
  fIndention.append("   ");	   
}  

void TG4XMLConvertor::CloseSection()
{
// Writes section closing.
// ---

  // define element
  G4String element = "</section>";

  // write element
  fOutFile << element
	   << G4endl;
}  

void TG4XMLConvertor::CloseComposition()
{
// Writes composition closing.
// ---

  // decrease indention
  fIndention.replace(fIndention.find("   "), 3 , "");

  // define element
  G4String element = "</composition>";

  // write element
  fOutFile << fIndention
           << element
	   << G4endl;
}  

void TG4XMLConvertor::WriteMaterial(const G4Material* material) 
{
// Writes G4Material. 
// Not yet implemented, only XML comment element is written.
// ---

  G4String name = material->GetName();
  CutName(name);

  // return if material of this name was already written
  if (fMaterialNames.find(name) != fMaterialNames.end()) return;
    
  fMaterialNames.insert(fMaterialNames.begin(), name); 

  // only comment line
  G4String element1 = "<!-- material = \""; 
  G4String element2 = "\" -->";
  
  // write element
  fOutFile << element1 << name
	   << element2
           << G4endl;
}  

void TG4XMLConvertor::WriteSolid(const G4VSolid* solid, G4String materialName) 
{
// Finds G4Solid concrete type and calls writing function. 
// For not yet implemented solids, only XML comment element is written.
// ---

  G4String name = solid->GetName();

  // return if solid of this name was already written
  if (fSolidNames.find(name) != fSolidNames.end()) return;
    
  fSolidNames.insert(fSolidNames.begin(), name); 

  // find concrete solid type and write it

  const G4Box* box = dynamic_cast<const G4Box*>(solid);
  if (box) { 
    WriteBox(box, materialName); 
    return;
  }
  
  const G4Tubs* tubs = dynamic_cast<const G4Tubs*>(solid);
  if (tubs) { 
    WriteTubs(tubs, materialName); 
    return;
  }
  
  const G4Trd* trd = dynamic_cast<const G4Trd*>(solid);
  if (trd) { 
    WriteTrd(trd, materialName); 
    return;
  }
  
  // write comment line in case of unsupported
  // shape

  // only comment line
  G4String element1 = "<!-- unsupported shape   name= \""; 
  G4String element2 = "\" -->";
  
  // write element
  fOutFile << element1 << name
	   << element2
           << G4endl;
}  

void TG4XMLConvertor::WriteRotation(const G4RotationMatrix* rotation)
{
// Writes G4RotationMatrix. 
// Not yet implemented, only XML comment element is written.
// ---

  // only comment line
  G4String element = "<!-- rotation matrix-->";
  
  // write element
  fOutFile << element 
           << G4endl;
}  

void TG4XMLConvertor::WritePosition(G4String solidName, G4ThreeVector position) 
{
// Writes position without rotation with a given solid name. 
// ---

  // get parameters
  G4double x = position.x()/TG3Units::Length();
  G4double y = position.y()/TG3Units::Length();
  G4double z = position.z()/TG3Units::Length();

  // compose element string template
  G4String element1 = "<posXYZ   volume=\"###########   X_Y_Z=\"";
  G4String element2 = "\" />";
  
  // put solid name
  PutName(element1, solidName, "#");
  
  // write element
  fOutFile << fIndention
           << element1
           << G4std::setw(7) << G4std::setprecision(2) << x << "  "
           << G4std::setw(7) << G4std::setprecision(2) << y << "  "
           << G4std::setw(7) << G4std::setprecision(2) << z
	   << element2
	   << G4endl;
}  

void TG4XMLConvertor::WritePositionWithRotation(
                           G4String solidName, G4ThreeVector position, 
			   const G4RotationMatrix* rotation)
{
// Writes position with rotation with a given solid name. 
// Not yet implemented, only XML comment element is written.
// ---

  // only comment line
  G4String element = "<!-- position with rotation -->"; 
  
  // write element
  fOutFile << fIndention
           << element 
           << G4endl;
}  

void TG4XMLConvertor::WriteEmptyLine()
{
// Writes empty line.
// ---

  fOutFile << G4endl;
}  


