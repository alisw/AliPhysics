// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliDetConstruction
// ------------------------
// Detector construction class with interactive setting of detectors setup 
// available.  
// In case the detector setup is not defined in Root macro Config.C
// the commands /alDet/switchOn/Off can be used either interactively or
// in Geant4 macro.

#ifndef ALI_DET_CONSTRUCTION_H
#define ALI_DET_CONSTRUCTION_H

#include "AliModulesComposition.h"
#include "AliDetSwitchVector.h"

#include <globals.hh>

class G4VPhysicalVolume;

class AliDetConstruction : public AliModulesComposition
{
  public:
    AliDetConstruction();
    // --> protected
    // AliDetConstruction(const AliDetConstruction& right);
    virtual ~AliDetConstruction();

    // methods
    virtual G4VPhysicalVolume* Construct();
    virtual void GenerateXMLGeometry() const;
    
  protected:
    AliDetConstruction(const AliDetConstruction& right);

    // operators
    AliDetConstruction& operator=(const AliDetConstruction& right);

  private:  
    // methods
    void BuildDetectors();
    void CreateDetectors();
    void CheckDetDependencies();
  
    // data members
    AliDetSwitchVector  fDetSwitchVector; //vector of AliDetSwitch
};

#endif //ALI_DET_CONSTRUCTION_H

