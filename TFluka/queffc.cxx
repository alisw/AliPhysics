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
    Double_t queffc(Double_t& wvlngt, Double_t& /*omgpho*/)
    {
	TGeoMaterial* material = (gGeoManager->GetCurrentVolume())->GetMaterial();
	Int_t nmat = material->GetIndex();
	TFlukaCerenkov*  cerenkov = dynamic_cast<TFlukaCerenkov*> (material->GetCerenkovProperties());
	Double_t y = 1.;
	if (cerenkov->IsSensitive()) 
	    y = (cerenkov->GetQuantumEfficiencyByWaveLength(wvlngt));
	
//	printf("queff: %e %d %e\n", wvlngt, nmat, y);
	return (y);
    }
}

