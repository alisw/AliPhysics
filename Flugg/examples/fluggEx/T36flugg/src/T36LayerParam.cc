// $Id$
// Flugg tag $Name$

//
// T36LayerParam.cc, 3/III/99, Sara Vanini
// parametrisation for t36 layers
//

#include "T36LayerParam.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

T36LayerParam::T36LayerParam(G4double LayerThickness, G4double NbOfLayers)
{
  fLayerThickness = LayerThickness; 
  fNbOfLayers = NbOfLayers;
}

T36LayerParam::~T36LayerParam()
{}

void T36LayerParam::ComputeTransformation
(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  static int counter = 0;
  counter +=1;
  G4double Xposition= - fLayerThickness/2*(fNbOfLayers-1) + copyNo*fLayerThickness;
  G4ThreeVector origin(Xposition,0,0);

  physVol->SetTranslation(origin);
}
