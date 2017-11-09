#ifndef ALIPICOBASE_H
#define ALIPICOBASE_H

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace AliPicoBase {

  enum {
    kPP = 0,  // pp    collisions
    kPA = 1,  // p-Pb  collisions
    kAP = 2,  // Pb-p  collisions
    kAA = 3   // Pb-Pb collisions
  };

  enum {
    kKshort     = BIT(0),  //
    kLambda     = BIT(1),  //
    kAntiLambda = BIT(2),  //
    kXiPos      = BIT(3),  //
    kXiNeg      = BIT(4),  //
    kDzero      = BIT(5),  //
    kDstar      = BIT(6)   //
  };

  enum {
    kEventCheck   = BIT(0),  //
    kEventMult    = BIT(1),  //
    kEventTrigger = BIT(2),  //
    kEventVertex  = BIT(3),  //
  };

  enum {
    kPrimary                = BIT(0), //
    kPhysicalPrimary        = BIT(1), //
    kSecondaryFromWeakDecay = BIT(2), //
    kSecondaryFromMaterial  = BIT(3)  //
  };

  inline Double_t MassPion()   { return 0.13957;  }
  inline Double_t MassKshort() { return 0.497614; }
  inline Double_t MassProton() { return 0.938272; }
  inline Double_t MassLambda() { return 1.11568;  }
}

#endif
