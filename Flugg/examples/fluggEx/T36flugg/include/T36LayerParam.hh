// $Id$
// Flugg tag $Name$

//
//  Parameterisation for t36 EM and HAD layers
//

#ifndef T36LayerParam_H
#define T36LayerParam_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4Box;

class T36LayerParam : public G4VPVParameterisation
{ 
  public:
    T36LayerParam(G4double LayerThickness, G4double NbOfLayers); 
    ~T36LayerParam();
    void ComputeTransformation
    (const G4int copyNo,G4VPhysicalVolume *physVol) const;

  private:

    G4double fLayerThickness;  
    G4double fNbOfLayers;      
};

#endif


