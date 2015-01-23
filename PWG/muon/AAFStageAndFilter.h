#ifndef AAFSTAGEANDFILTER_H
#define AAFSTAGEANDFILTER_H

#ifndef TROOT_TString
#  include "TString.h"
#endif

namespace AAF {
  
  Int_t StageAndFilter(TString from, const TString& to,
                       TString filterName);
}


#endif
