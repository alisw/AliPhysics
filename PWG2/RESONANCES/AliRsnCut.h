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

class AliRsnEvent;

class AliRsnCut : public TNamed
{
  public:

    // possible targets for a cut
    enum ETarget 
    {
      kDaughter = 0,
      kMother,
      kEvent,
      kMixEvent,
      kLastCutTarget
    };

    // data type for check
    enum EVarType 
    {
      kNoVar = 0,
      kInt,
      kDouble
    };

    AliRsnCut(ETarget target = kLastCutTarget);
    AliRsnCut(const AliRsnCut& copy);
    AliRsnCut& operator=(const AliRsnCut& copy);
    AliRsnCut(const char *name, ETarget target, Int_t    min, Int_t    max = 0 );
    AliRsnCut(const char *name, ETarget target, Double_t min, Double_t max = 0.);
    virtual ~AliRsnCut() { /*nothing*/ };
    
    EVarType         GetVarTypeEnum() {return fVarType;}
    Char_t           GetVarTypeChar() {if (fVarType == kInt) return 'I'; else if (fVarType == kDouble) return 'D'; else return 'X';}
    ETarget          GetTargetEnum()  {return fTarget;}
    Char_t           GetTargetChar()  {if (fTarget == kDaughter) return 'D'; else if (fTarget == kMother) return 'M'; else if (fTarget == kEvent) return 'E'; else return 'X';}
    Bool_t           IsTarget(ETarget target) {return (fTarget == target);}
    Bool_t           TargetOK(TObject *obj1, TObject *obj2 = 0x0);
    Int_t            GetMinI()        {return fMinI;}
    Int_t            GetMaxI()        {return fMaxI;}
    Double_t         GetMinD()        {return fMinD;}
    Double_t         GetMaxD()        {return fMaxD;}
    Int_t            GetCutValueI()   {return fCutValueI;}
    Double_t         GetCutValueD()   {return fCutValueD;}
    Bool_t           GetCutResult()   {return fCutResult;}
    
    void             SetRange(Int_t    min, Int_t    max) {fMinI = min; fMaxI = max; fVarType = kInt;}
    void             SetRange(Double_t min, Double_t max) {fMinD = min; fMaxD = max; fVarType = kDouble;}

    void             SetValue(Int_t value)    {fMinI = value; fVarType = kInt;}
    void             SetValue(Double_t value) {fMinD = value; fVarType = kDouble;}
    
    Bool_t           OkValue();
    Bool_t           OkRange();
    Bool_t           OkValueI();
    Bool_t           OkRangeI();
    Bool_t           OkValueD();
    Bool_t           OkRangeD();

    virtual void     SetEvent(AliRsnEvent *event);
    AliRsnEvent*     GetEvent() {return fEvent;}
    
    virtual Bool_t   IsSelected(TObject *obj1, TObject *obj2 = 0x0);
    virtual void     Print(Option_t *opt = "") const;

  protected:

    EVarType     fVarType;    // type of checked variable
    ETarget      fTarget;     // type of object on which the cut is checked

    Int_t        fMinI;       // lower edge of INT range or ref. value    for INT CUT
    Int_t        fMaxI;       // upper edge of INT range (not used    for value cuts)
    Double_t     fMinD;       // lower edge of DOUBLE range or ref. value    for INT CUT
    Double_t     fMaxD;       // upper edge of DOUBLE range (not used    for value cuts)

    Int_t        fCutValueI;  // cut value INT
    Double_t     fCutValueD;  // cut value DOUBLE
    
    Bool_t       fCutResult;  // tells if the cut is passed or not

    AliRsnEvent *fEvent;      //! pointer to current event (can be needed sometimes, but never streamed)

    ClassDef(AliRsnCut, 1)
};

#endif
