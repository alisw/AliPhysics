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
#include <G4Trap.hh>
#include <globals.hh>

#include <g4std/iostream>
#include <g4std/iomanip>

const G4int TG4XMLConvertor::fgkMaxVolumeNameLength   = 10;
const G4int TG4XMLConvertor::fgkMaxMaterialNameLength = 20;

TG4XMLConvertor::TG4XMLConvertor(G4std::ofstream& outFile) 
  : fOutFile(outFile),
    fBasicIndention("   "),
    fIndention(fBasicIndention),
    fRotationCounter(0)
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
  
void TG4XMLConvertor::WriteBox(G4String lvName, const G4Box* box, 
                               G4String materialName)
{
// Writes G4box solid.
// ---

  // get parameters
  G4double x = box->GetXHalfLength()/TG3Units::Length();
  G4double y = box->GetYHalfLength()/TG3Units::Length();
  G4double z = box->GetZHalfLength()/TG3Units::Length();

  // compose element string template
  G4String element1 
    = "<box   name=\"###########   material=\"!!!!!!!!!!!!!!!!!!!!!   X_Y_Z=\"";
  G4String element2 = "\" />";
  
  // put solid and material names
  PutName(element1, lvName, "#");
  PutName(element1, materialName, "!");
  
  // write element
  fOutFile << fBasicIndention
           << element1
           << G4std::setw(7) << G4std::setprecision(2) << x*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << y*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << z*2.
	   << element2
	   << G4endl;
}
 
void TG4XMLConvertor::WriteTubs(G4String lvName, const G4Tubs* tubs, 
                                G4String materialName)
{
// Writes G4tubs solid.
// ---

  // get parameters
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
  PutName(element1, lvName, "#");
  PutName(element1, materialName, "!");
  
  // write element
  fOutFile << fBasicIndention
           << element1
           << G4std::setw(7) << G4std::setprecision(2) << rmin << "  "
           << G4std::setw(7) << G4std::setprecision(2) << rmax << "  "
           << G4std::setw(7) << G4std::setprecision(2) << hz*2.
	   << element2
           << G4std::setw(7) << G4std::setprecision(2) << sphi << "  "
           << G4std::setw(7) << G4std::setprecision(2) << sphi+dphi
	   << element3
	   << G4endl;
}  


void TG4XMLConvertor::WriteTrd(G4String lvName, const G4Trd* trd, 
                               G4String materialName)
{
// Writes G4Trd solid.
// ---

  // get parameters
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
  PutName(element1, lvName, "#");
  PutName(element1, materialName, "!");
  
  // write element
  fOutFile << fBasicIndention
           << element1
           << G4std::setw(7) << G4std::setprecision(2) << x1*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << x2*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << y1*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << y2*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << hz*2.
	   << element2
	   << G4endl;
}  


void TG4XMLConvertor::WriteTrap(G4String lvName, const G4Trap* trap, 
                                G4String materialName)
{
// Writes G4Trap solid.
// ---

  // get parameters
  G4double dz = trap->GetZHalfLength()/TG3Units::Length();
  G4ThreeVector symAxis = trap->GetSymAxis();
  G4double y1 = trap->GetYHalfLength1()/TG3Units::Length();
  G4double x1 = trap->GetXHalfLength1()/TG3Units::Length();
  G4double x2 = trap->GetXHalfLength2()/TG3Units::Length();
  G4double tanAlpha1 = trap->GetTanAlpha1();
  G4double y2 = trap->GetYHalfLength2()/TG3Units::Length();
  G4double x3 = trap->GetXHalfLength3()/TG3Units::Length();
  G4double x4 = trap->GetXHalfLength4()/TG3Units::Length();
  G4double tanAlpha2 = trap->GetTanAlpha2();

  // ordering of parameters in XML element
  // Xmumdpupd_Ymp_Z: 2x2 2x1 2x4 2x3 2y2 2y1 2dz
  // inclination: atan(symAxis.x/symAxis.z), atan(symAxis.y/symAxis.z)
  // declination: alpha1, alpha2

  // get angles
  G4double inc1 = atan(symAxis.x()/symAxis.z()) / deg;
  G4double inc2 = atan(symAxis.y()/symAxis.z()) / deg;
  G4double alpha1 = atan(tanAlpha1) / deg;
  G4double alpha2 = atan(tanAlpha2) / deg;

  // compose element string template
  G4String element1 
    = "<trap  name=\"###########   material=\"!!!!!!!!!!!!!!!!!!!!!   Xmumdpupd_Ymp_Z=\"";
  G4String element2
    = "                                                              inclination=\""; 
  G4String element3 = "   declination=\""; 
  G4String element4 = "\" />";
  G4String quota = "\"";   
  
  // put solid and material names
  // put solid and material names
  PutName(element1, lvName, "#");
  PutName(element1, materialName, "!");
  
  // write element
  fOutFile << fBasicIndention
           << element1
           << G4std::setw(7) << G4std::setprecision(2) << x2*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << x1*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << x4*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << x3*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << y2*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << y1*2. << "  "
           << G4std::setw(7) << G4std::setprecision(2) << dz*2. << quota << G4endl
           << fBasicIndention
	   << element2
           << G4std::setw(7) << G4std::setprecision(2) << inc1 << "  "
           << G4std::setw(7) << G4std::setprecision(2) << inc2 << quota
	   << element3
           << G4std::setw(7) << G4std::setprecision(2) << alpha1 << "  "
           << G4std::setw(7) << G4std::setprecision(2) << alpha2 
	   << element4
	   << G4endl;
}  

// public methods

void TG4XMLConvertor::OpenMaterials(const G4String& version, 
 	                 const G4String& date, const G4String& author,
                         const G4String dtdVersion)
{
// Writes section opening.
// ---
			 
  G4String element1 = "<materials  version = \"";
  G4String element2 = "            date    = \"";
  G4String element3 = "            author  = \"";
  G4String element4 = "            DTD_version=\"";
  G4String element5 = "  >";
  G4String quota = "\"";   
  
  // write element
  fOutFile << element1 << version << quota << G4endl
           << element2 << date    << quota << G4endl
           << element3 << author  << quota << G4endl
           << element4 << dtdVersion << quota
           << element5 << G4endl;
}  

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
  element.append("\">");

  // write element
  fOutFile << fIndention
           << element
	   << G4endl;

  // increase indention
  IncreaseIndention();	   
}  

void TG4XMLConvertor::CloseMaterials()
{
// Writes materials closing.
// ---

  // define element
  G4String element = "</materials>";

  // write element
  fOutFile << element
	   << G4endl;
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
  DecreaseIndention();

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

  // only comment line
  G4String element1 = "<!-- material = \""; 
  G4String element2 = "\" -->";
  
  // write element
  fOutFile << fBasicIndention
           << element1 << name
	   << element2
           << G4endl;
}  

void TG4XMLConvertor::WriteSolid(G4String lvName, const G4VSolid* solid, 
                                 G4String materialName) 
{
// Finds G4Solid concrete type and calls writing function. 
// For not yet implemented solids, only XML comment element is written.
// ---

  // to be removed when materials are supported
  materialName = "Hydrogen";
  
  const G4Box* box = dynamic_cast<const G4Box*>(solid);
  if (box) { 
    WriteBox(lvName, box, materialName); 
    return;
  }
  
  const G4Tubs* tubs = dynamic_cast<const G4Tubs*>(solid);
  if (tubs) { 
    WriteTubs(lvName, tubs, materialName); 
    return;
  }
  
  const G4Trd* trd = dynamic_cast<const G4Trd*>(solid);
  if (trd) { 
    WriteTrd(lvName, trd, materialName); 
    return;
  }
  
  const G4Trap* trap = dynamic_cast<const G4Trap*>(solid);
  if (trap) { 
    WriteTrap(lvName, trap, materialName); 
    return;
  }
  
  // write comment line in case of unsupported
  // shape

  // only comment line
  G4String element1 = "<!-- unsupported shape   name= \""; 
  G4String element2 = "\" -->";
  
  // write element
  fOutFile << fBasicIndention
           << element1 << lvName
	   << element2
           << G4endl;
}  

void TG4XMLConvertor::WriteRotation(const G4RotationMatrix* rotation)
{
// Writes G4RotationMatrix. 
// Not yet implemented, only XML comment element is written.
// ---

  // return if this rotation was already written
  G4int nofRotations = fRotations.size();
  if (nofRotations>0)
    for (G4int i=0; i<nofRotations; i++) 
      if (fRotations[i] == rotation) return;

  fRotations.push_back(rotation);  


  // get parameters
  G4double xx = rotation->xx();
  G4double xy = rotation->xy();
  G4double xz = rotation->xz();
  G4double yx = rotation->yx();
  G4double yy = rotation->yy();
  G4double yz = rotation->yz();
  G4double zx = rotation->zx();
  G4double zy = rotation->zy();
  G4double zz = rotation->zz();
  G4String id = "RM";
  TG4Globals::AppendNumberToString(id, ++fRotationCounter);
 
  // compose element string template
  G4String quota = "\"\n";
  G4String element1 = "<rot_matrix   id=\"#######  XX_XY_XZ=\"";
  G4String element2 = "                           YX_YY_YZ=\"";
  G4String element3 = "                           ZX_ZY_ZZ=\"";
  G4String element4 = "\" />";
  
  // put identifier
  PutName(element1, id, "#");

  // write element
  fOutFile << fBasicIndention
           << element1
	   << G4std::setw(8) << G4std::setprecision(5) << xx << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << xy << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << xz << quota
           << fBasicIndention
           << element2
	   << G4std::setw(8) << G4std::setprecision(5) << yx << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << yy << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << yz << quota
	   << fBasicIndention
           << element3
	   << G4std::setw(8) << G4std::setprecision(5) << zx << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << zy << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << zz 
	   << element4 	   
           << G4endl;
}  

void TG4XMLConvertor::WritePosition(G4String lvName, G4ThreeVector position) 
{
// Writes position without rotation with a given solid name. 
// ---

  // get parameters
  G4double x = position.x()/TG3Units::Length();
  G4double y = position.y()/TG3Units::Length();
  G4double z = position.z()/TG3Units::Length();

  // compose element string template
  G4String element1 = "<posXYZ      volume=\"###########   X_Y_Z=\"";
  G4String element2 = "\" />";
  
  // put solid name
  PutName(element1, lvName, "#");
  
  // write element
  fOutFile << fIndention
           << element1
           << G4std::setw(8) << G4std::setprecision(2) << x << "  "
           << G4std::setw(8) << G4std::setprecision(2) << y << "  "
           << G4std::setw(8) << G4std::setprecision(2) << z
	   << element2
	   << G4endl;
}  

void TG4XMLConvertor::WritePositionWithRotation(
                           G4String lvName, G4ThreeVector position, 
			   const G4RotationMatrix* rotation)
{
// Writes position with rotation with a given solid name. 
// Not yet implemented, only XML comment element is written.
// ---

  // get parameters
  G4double x = position.x()/TG3Units::Length();
  G4double y = position.y()/TG3Units::Length();
  G4double z = position.z()/TG3Units::Length();
  G4double xx = rotation->xx();
  G4double xy = rotation->xy();
  G4double xz = rotation->xz();
  G4double yx = rotation->yx();
  G4double yy = rotation->yy();
  G4double yz = rotation->yz();
  G4double zx = rotation->zx();
  G4double zy = rotation->zy();
  G4double zz = rotation->zz();
  
/*
  // find rotation
  G4int i=0;
  while (i<fRotations.size() && fRotations[i] != rotation) i++; 
  if (i==fRotations.size()) {
    G4String text = "TG4XMLConvertor::WritePositionWithRotation: ";
    text = text + "    Unknown rotation - fatal error.";    
    TG4Globals::Exception(text);
  }  
  G4String id = "RM";
  TG4Globals::AppendNumberToString(id, i); 
*/  

  // compose element string template
  G4String quota = "\"\n";
  G4String element1 = "<transform   volume=\"###########     pos=\"";
  G4String element2 = "                                     rot=\"";
  G4String element3 = "                                          ";
  G4String element4 = "\" />";
  
  // put solid name
  PutName(element1, lvName, "#");
  
  // write element
  fOutFile << fIndention
           << element1
           << G4std::setw(8) << G4std::setprecision(2) << x << "  "
           << G4std::setw(8) << G4std::setprecision(2) << y << "  "
           << G4std::setw(8) << G4std::setprecision(2) << z << quota
	   << fIndention
	   << element2 
	   << G4std::setw(8) << G4std::setprecision(5) << xx << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << xy << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << xz << G4endl
           << fIndention
           << element3
	   << G4std::setw(8) << G4std::setprecision(5) << yx << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << yy << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << yz << G4endl
	   << fIndention
           << element3
	   << G4std::setw(8) << G4std::setprecision(5) << zx << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << zy << "  "  
	   << G4std::setw(8) << G4std::setprecision(5) << zz 
	   << element4
	   << G4endl;
}  

void TG4XMLConvertor::WriteEmptyLine()
{
// Writes empty line.
// ---

  fOutFile << G4endl;
}  

void TG4XMLConvertor::IncreaseIndention()
{
  // increase indention
  fIndention.append(fBasicIndention);	   
}

void TG4XMLConvertor::DecreaseIndention()
{
  // decrease indention
  fIndention.replace(fIndention.find(fBasicIndention), 3 , "");
}
