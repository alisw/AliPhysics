#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#ifndef WIN32
# define rflctv rflctv_
#else
# define rflctv RFLCTV
#endif
extern "C" {
Double_t rflctv(Double_t& wvlngt, Double_t& omgpho, Int_t& mmat)
{
    printf("rflctv called  %e %e %d \n", wvlngt, omgpho, mmat);
    return (0.);
}
}
