#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#ifndef WIN32
# define rfrndx rfrndx_
#else
# define rfrndx RFRNDX
#endif
extern "C" {
Double_t rfrndx(Double_t& wvlngt, Double_t& omgpho, Int_t& mmat)
{
    printf("rfrndx called  %e %e %d \n", wvlngt, omgpho, mmat);
    return (0.);
}
}

