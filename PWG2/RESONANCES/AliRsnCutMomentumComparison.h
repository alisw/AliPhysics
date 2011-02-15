//
// Class AliRsnCutRange
//
// General implementation of cuts which check a value inside a range.
// This range can be defined by two integers or two doubles.
// A user-friendly enumeration allows to define what is checked.
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTMOMENTUMCOMPARISON_H
#define ALIRSNCUTMOMENTUMCOMPARISON_H

#include "AliRsnCut.h"

class AliRsnCutMomentumComparison : public AliRsnCut {
public:

   enum EMode {
      kFirstLargerP,
      kFirstSmallerP,
      kFirstLargerPt,
      kFirstSmallerPt
   };

   AliRsnCutMomentumComparison(const char *name = "cutMomComparison", EMode mode = kFirstLargerPt);
   AliRsnCutMomentumComparison(const AliRsnCutMomentumComparison& copy);
   AliRsnCutMomentumComparison& operator=(const AliRsnCutMomentumComparison& copy);
   virtual ~AliRsnCutMomentumComparison() {;};

   virtual Bool_t IsSelected(TObject *object);

protected:

   EMode fMode;     // comparison mode

   ClassDef(AliRsnCutMomentumComparison, 1)
};

#endif
