#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#include "TFluka.h"
#include "TGeoMaterial.h"
#include "TFlukaCerenkov.h"


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
//
//  Check if stopping has been required by user
//
    if (fluka->GetStoppingCondition()) {
	fluka->ResetStoppingCondition();
	return (1.e15);
    }
//
//  Get absorption coefficient for current material
//    
    TGeoMaterial*    material =  (TGeoMaterial*) (fluka->GetFlukaMaterials())->At(fluka->GetMaterialIndex(mmat));
    TFlukaCerenkov*  cerenkov = dynamic_cast<TFlukaCerenkov*> (material->GetCerenkovProperties());
    Double_t y = (cerenkov->GetAbsorptionCoefficientByWaveLength(wvlngt));
    return (y);
    
}
}


