// $Id$
// Category: digits+hits
//
// See the class description in the header file.

#include "TG4SDManager.h"
#include "TG4VSDConstruction.h"
#include "TG4SDServices.h"
#include "TG4Globals.h"

TG4SDManager* TG4SDManager::fgInstance = 0;

//_____________________________________________________________________________
TG4SDManager::TG4SDManager(TG4VSDConstruction* sdConstruction)
  : fSDConstruction(sdConstruction) {
//
  if (fgInstance)
    TG4Globals::Exception(
      "TG4SDManager: attempt to create two instances of singleton.");
      
  fgInstance = this; 
  
  fSDServices = new TG4SDServices();
}

//_____________________________________________________________________________
TG4SDManager::TG4SDManager(const TG4SDManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4SDManager singleton.");
}


//_____________________________________________________________________________
TG4SDManager::TG4SDManager() {
// 
  TG4Globals::Exception(
    "TG4SDManager has protected the default constructor.");
}


//_____________________________________________________________________________
TG4SDManager::~TG4SDManager(){
//

  delete fSDServices;
}

// operators

//_____________________________________________________________________________
TG4SDManager& TG4SDManager::operator=(const TG4SDManager& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4SDManager singleton.");
    
  return *this;  
}    
          

// public methods 

//_____________________________________________________________________________
void TG4SDManager::Initialize() 
{
// Creates sensitive detectors,
// sets second indexes for materials (corresponding to G3 tracking 
// media) and clears remaing G3 tables.
// ---

  fSDConstruction->Construct();
}  
  

//_____________________________________________________________________________
Int_t TG4SDManager::VolId(const Text_t* volName) const
{ 
// Returns the sensitive detector identifier.
// ---

  return fSDServices->GetVolumeID(volName);
}


//_____________________________________________________________________________
const char* TG4SDManager::VolName(Int_t id) const
{
// Returns the name of the sensitive detector with the given identifier.
// ---

  return fSDServices->GetVolumeName(id);
}


//_____________________________________________________________________________
Int_t TG4SDManager::NofVolumes() const
{
// Returns the total number of sensitive detectors.
// ---

  return fSDServices->NofSensitiveDetectors();
}


//_____________________________________________________________________________
Int_t TG4SDManager::VolId2Mate(Int_t volumeId)  const
{
// Return the material number for a given volume id
// ---

  return fSDServices->GetMediumId(volumeId);	       	         
}
 
