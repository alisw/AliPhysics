// $Id$
// Category: interfaces
//
// Author: D. Adamova
//
//==================================================================
//
//----------------TG4GuiVolume.cxx-------------------------------//
//----Creating link for Logical Volume Tree in AG4 Geometry----//
//
//===================================================================
 
 
 
#include "TG4GuiVolume.h"
#include <G4LogicalVolume.hh>
 

ClassImp(TG4GuiVolume)

TG4GuiVolume::TG4GuiVolume(const char* name, G4LogicalVolume* lvolume)
{
// Constructor
    fkName   = name;
    fItem   = 0;
    fLogicalVolume = lvolume; 
}

G4LogicalVolume* TG4GuiVolume::GetLogicalVolume() const
{
// Returns logical volume
  return fLogicalVolume;
}
