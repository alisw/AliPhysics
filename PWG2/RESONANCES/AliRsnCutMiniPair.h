#ifndef ALIRSNCUTMINIPAIR_H
#define ALIRSNCUTMINIPAIR_H

//
// This cut definition works on mini-pairs for the 'mini' subpackage.
// Since cuts on mini-pairs can usually be just kinematic selections
// or kinematic comparisons between 4-momenta of daughters, they are all
// implemented in this class, by the use of an enumeration which allows
// the user to choose between all possibilities.
//

#include "AliRsnCut.h"

class AliRsnCutMiniPair : public AliRsnCut {

public:

   enum EType {
      kRapidityRange,
      kRapidityRangeMC,
      kMomentumComparison,
      kTypes
   };

   AliRsnCutMiniPair(const char *name = "cut", EType type = kTypes);
   virtual ~AliRsnCutMiniPair() { }
   
   virtual Bool_t IsSelected(TObject *obj);

private:

   EType fType;    // cut type

   ClassDef(AliRsnCutMiniPair,1)

};

#endif
