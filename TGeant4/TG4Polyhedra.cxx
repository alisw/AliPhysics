// $Id$
// Category: geometry
//
// Author: I. Hrivnacova, 12.10.2000 
//
// Class TG4Polyhedra
// ------------------
// See the class description in the header file.

#include "TG4Polyhedra.h"
#include "TG4Globals.h"

//_____________________________________________________________________________
TG4Polyhedra::TG4Polyhedra(const G4Polyhedra& rhs) 
  : G4Polyhedra(rhs) {
//
}

//_____________________________________________________________________________
TG4Polyhedra::~TG4Polyhedra() {
//
}

// private methods

//_____________________________________________________________________________
void TG4Polyhedra::CheckOrigin() 
{
// Checks if polycone was created in a "historical way"
// and give exception otherwise.
// ---

  if (!original_parameters) {
    G4String text = "TG4Polyhedra::CheckOrigin: \n";
    text = text + "    Polycone has not defined original parameters.";
    TG4Globals::Exception(text);
  }  
}


// public methods

//_____________________________________________________________________________
G4int TG4Polyhedra::GetNofZPlanes()
{
// Returns nof z planes.
// ----
  
  CheckOrigin();

  return original_parameters->Num_z_planes;
}  


//_____________________________________________________________________________
G4double* TG4Polyhedra::GetRmin()
{
// Returns array of rmin parameters of the planes.
// ----
  
  CheckOrigin();

  return original_parameters->Rmin;
}  

//_____________________________________________________________________________
G4double* TG4Polyhedra::GetRmax()
{
// Returns array of rmax parameters of the planes.
// ----

  CheckOrigin();
  
  return original_parameters->Rmax;
}  

//_____________________________________________________________________________
G4double* TG4Polyhedra::GetZ()
{
// Returns array of z parameters of the planes.
// ----

  CheckOrigin();
  
  return original_parameters->Z_values;
}  
