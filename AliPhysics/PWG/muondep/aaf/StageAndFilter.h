#ifndef AAFSTAGEANDFILTER_H
#define AAFSTAGEANDFILTER_H

#ifndef TROOT_TString
#  include "TString.h"
#endif

namespace AAF {
  
  Int_t StageAndFilter(const TString& from, const TString& to,
                       	              const TString& filterName, Int_t verboseLevel);
}


#endif
