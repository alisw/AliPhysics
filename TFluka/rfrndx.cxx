#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#include "TFluka.h"
#include "TGeoMaterial.h"
#include "TFlukaCerenkov.h"

#ifndef WIN32
# define rfrndx rfrndx_
#else
# define rfrndx RFRNDX
#endif
extern "C" {
Double_t rfrndx(Double_t& wvlngt, Double_t& /*omgpho*/, Int_t& mmat)
{
//
//  Return refraction index for given photon energy and material
//
    TFluka* fluka =  (TFluka*) gMC;
    TGeoMaterial*    material =  (TGeoMaterial*) (fluka->GetFlukaMaterials())->At(fluka->GetMaterialIndex(mmat));
    TFlukaCerenkov*  cerenkov = dynamic_cast<TFlukaCerenkov*> (material->GetCerenkovProperties());
    Double_t y =  (cerenkov->GetRefractionIndexByWaveLength(wvlngt));
    return (y);
}
}
