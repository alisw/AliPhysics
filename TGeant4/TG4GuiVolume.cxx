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
#include "TG4Globals.h"
#include <G4LogicalVolume.hh> 

ClassImp(TG4GuiVolume)

TG4GuiVolume::TG4GuiVolume(const char* name, G4LogicalVolume* lvolume)
{
// Constructor
    fItem   = 0;
    fLogicalVolume = lvolume; 
    
    G4String lName = fLogicalVolume->GetName();
    
    if ( lName != name ) TG4Globals::Exception(
       "A wrong name assigned to the guiVolume in the ctor" );
}

TG4GuiVolume::TG4GuiVolume(const TG4GuiVolume& gv) 
{
// Dummy copy constructor 
  TG4Globals::Exception(
    "Attempt to use TG4GuiVolume copy constructor.");
}

TG4GuiVolume& TG4GuiVolume::operator=(const TG4GuiVolume& gv)
{
  // check assignement to self
  if (this == &gv) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4GuiVolume singleton.");
    
  return *this;  
}

const char* TG4GuiVolume::GetName() const
{
// Returns the gui/logical volume name
   
  G4String lName = fLogicalVolume->GetName();
  return lName; 
}    

TGListTreeItem* TG4GuiVolume::GetItem() const
{
// Returns ListTree item

    return fItem;
    
}

G4LogicalVolume* TG4GuiVolume::GetLogicalVolume() const
{
// Returns logical volume

  return fLogicalVolume;

}

