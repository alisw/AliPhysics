#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#include "Fopphcm.h"  //(OPPHCM) fluka common
#include "TGeoMaterial.h"
#include "TFlukaCerenkov.h"
#include "TFlukaGeo.h"
#include "TGeoManager.h"

#ifndef WIN32
# define queffc queffc_
#else
# define queffc QUEFFC
#endif
extern "C" {
Double_t queffc(Double_t& wvlngt, Double_t& omgpho, Int_t& mmat)
{
    printf("queffc called  %e %e %d \n", wvlngt, omgpho, mmat);
    TFluka* fluka =  (TFluka*) gMC;
    TGeoMaterial*    material =  (TGeoMaterial*) (gGeoManager->GetListOfMaterials())->At(fluka->GetMaterialIndex(mmat));
    TFlukaCerenkov*  cerenkov = dynamic_cast<TFlukaCerenkov*> (material->GetCerenkovProperties());
    return (cerenkov->GetQuantumEfficiencyByWaveLength(wvlngt));
}
}
