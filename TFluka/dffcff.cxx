#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fiounit.h"  //(IOUNIT) fluka common
#ifndef WIN32
# define dffcff dffcff_
#else
# define dffcff DFFCFF
#endif
extern "C" {
Double_t dffcff(Double_t& /*wvlngt*/, Double_t& /*omgpho*/, Int_t& /*mmat*/)
{
//    printf("dffcff called  %e %e %d \n", wvlngt, omgpho, mmat);
    return (0.);
}
}

    
