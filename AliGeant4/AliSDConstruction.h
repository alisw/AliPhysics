// $Id$
// Category: digits+hits
//
// Author: I. Hrivnacova
//
// Class AliSDConstruction
// -----------------------
// Mandatory (TGeant4) class providing method for creating 
// sensitive detectors.
// It provides also methods for switching between lego 
// sensitive detectors and the standard ones.

#ifndef ALI_SD_CONSTRUCTION_H
#define ALI_SD_CONSTRUCTION_H

#include "TG4VSDConstruction.h"

#include <globals.hh>

class AliLego;
class AliModule;

class G4VPhysicalVolume;
class G4LogicalVolume;

class AliSDConstruction : public TG4VSDConstruction 
{  
  public:
    AliSDConstruction();
    virtual ~AliSDConstruction();
    
    // methods    
    virtual void Construct();

    void SetLego(AliLego* lego) const;
    void UnsetLego() const;

  private:
    // methods
    void InitializeModules();
    AliModule* FindAliModule(G4LogicalVolume* lv) const;
    void CreateSD(G4LogicalVolume* lv, AliModule* module) const;
    void CreateLegoSD(G4LogicalVolume* lv, AliLego* lego) const;
    void UnsetLegoSD(G4LogicalVolume* lv) const;    
};             

#endif //ALI_SD_CONSTRUCTION_H
