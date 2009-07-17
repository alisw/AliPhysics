//
// Class AliRsnCut
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()]
// - a value equal to a given reference     [--> IsEqual()  ]
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUT_H
#define ALIRSNCUT_H

#include "TNamed.h"

class AliRsnDaughter;
class AliRsnPairParticle;
class AliRsnEvent;

class AliRsnCut : public TNamed
{
  public:

    // possible targets for a cut
    enum ETarget {
      kParticle = 0,
      kPair,
      kEvent,
      kMixEvent,
      kLastCutTarget
    };

    // data type for check
    enum EVarType {
      kInt = 0,
      kULong,
      kDouble
    };

    AliRsnCut();
    AliRsnCut(const char *name, Int_t    min, Int_t    max = 0);
    AliRsnCut(const char *name, ULong_t  min, ULong_t  max = 0);
    AliRsnCut(const char *name, Double_t min, Double_t max = 0);
    virtual ~AliRsnCut() {;};

    void             SetRange(Int_t    min, Int_t    max) {fMinI = min; fMaxI = max; fVarType = kInt;}
    void             SetRange(ULong_t  min, ULong_t  max) {fMinU = min; fMaxU = max; fVarType = kULong;}
    void             SetRange(Double_t min, Double_t max) {fMinD = min; fMaxD = max; fVarType = kDouble;}

    void             SetValue(Int_t value)    {fMinI = value; fVarType = kInt;}
    void             SetValue(ULong_t value)  {fMinU = value; fVarType = kULong;}
    void             SetValue(Double_t value) {fMinD = value; fVarType = kDouble;}

    virtual Bool_t   IsSelected(ETarget tgt, AliRsnDaughter *daughter);
    virtual Bool_t   IsSelected(ETarget tgt, AliRsnPairParticle *pair);
    virtual Bool_t   IsSelected(ETarget tgt, AliRsnEvent *event);
    virtual Bool_t   IsSelected(ETarget tgt, AliRsnEvent *ev1, AliRsnEvent *ev2);

  protected:

    Bool_t  OkValue();
    Bool_t  OkRange();

    EVarType  fVarType;    // type of checked variable

    Int_t     fMinI;       // lower edge of INT range or ref. value for INT CUT
    Int_t     fMaxI;       // upper edge of INT range (not used for value cuts)
    ULong_t   fMinU;       // lower edge of ULONG range or ref. value for INT CUT
    ULong_t   fMaxU;       // upper edge of ULONG range (not used for value cuts)
    Double_t  fMinD;       // lower edge of DOUBLE range or ref. value for INT CUT
    Double_t  fMaxD;       // upper edge of DOUBLE range (not used for value cuts)

    Int_t     fCutValueI;  // cut value
    ULong_t   fCutValueU;  // cut value
    Double_t  fCutValueD;  // cut value
    Bool_t    fCutResult;  // tells if the cut is passed or not

    ClassDef(AliRsnCut, 1)
};

#endif
