#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#include "Fopphcm.h"  //(OPPHCM) fluka common
#ifndef WIN32
# define queffc queffc_
#else
# define queffc QUEFFC
#endif
extern "C" {
Double_t queffc(Double_t& wvlngt, Double_t& omgpho, Int_t& mmat)
{
    printf("queffc called  %e %e %d \n", wvlngt, omgpho, mmat);
    if (wvlngt  > OPPHCM.wvmxsn || wvlngt < OPPHCM.wvmnsn) {
	return (0.);
    } else {
	return (OPPHCM.opsnmx);
    }
}
}
