#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#include "TFlukaGeo.h"
#include "TGeoMaterial.h"
#include "TFlukaCerenkov.h"
#include "TGeoManager.h"


#ifndef WIN32
# define abscff abscff_
#else
# define abscff ABSCFF
#endif
extern "C" {
Double_t abscff(Double_t& wvlngt, Double_t& /*omgpho*/, Int_t& mmat)
{
//
//  Return absorption length  for given photon energy and material
//
    TFluka* fluka =  (TFluka*) gMC;
    TGeoMaterial*    material =  (TGeoMaterial*) (gGeoManager->GetListOfMaterials())->At(fluka->GetMaterialIndex(mmat));
    TFlukaCerenkov*  cerenkov = dynamic_cast<TFlukaCerenkov*> (material->GetCerenkovProperties());
    Double_t y = (cerenkov->GetAbsorptionCoefficientByWaveLength(wvlngt));
    return (y);
    
}
}


