#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#include "Fopphcm.h"  //(OPPHCM) fluka common
#include "TGeoMaterial.h"
#include "TFlukaCerenkov.h"
#include "TFluka.h"
#include "TGeoManager.h"

#ifndef WIN32
# define queffc queffc_
#else
# define queffc QUEFFC
#endif
extern "C" {
    Double_t queffc(Double_t& /*wvlngt*/, Double_t& /*omgpho*/)
    {
	Double_t eff = TFlukaCerenkov::GetGlobalMaximumEfficiency();
	return (eff);
    }
}

