#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#ifndef WIN32
# define abscff abscff_
#else
# define abscff ABSCFF
#endif
extern "C" {
Double_t abscff(Double_t& wvlngt, Double_t& omgpho, Int_t& mmat)
{
    printf("abscff called  %e %e %d \n", wvlngt, omgpho, mmat);
    return (0.);
}
}


