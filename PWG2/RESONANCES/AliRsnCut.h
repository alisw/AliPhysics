/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//=========================================================================
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
//=========================================================================

#ifndef ALIRSNCUT_H
#define ALIRSNCUT_H

#include "TNamed.h"

class AliRsnDaughter;
class AliRsnPairParticle;
class AliRsnPairDef;

class AliRsnCut : public TNamed
{
public:
    enum ERsnCutType
    {
        kMomentum = 0,
        kTransMomentum,
        kEta,
        kRadialImpactParam,
        kMomentumMC,
        kTransMomentumMC,
        kEtaMC,
        kRestMomentum,
        kNSigma,
        kNSigmaCalculate,
        kStatus,
        kIsPdgEqual,
        kIsLabelEqual,
        kIsTruePair,
        kChargePos,
        kChargeNeg,
        kPIDType,
        kLastCutType
    };

    enum ERsnCutVarType
    {
        kDouble_t = 0,
        kInt_t,
        kUInt_t
    };

    enum ECutSetType
    {
        kParticle = 0,
        kPair,
        kMixEventFinderCut,
        kLastCutSetIndex
    };

    AliRsnCut();
    AliRsnCut (const char *name, const char *title, ERsnCutType type);
    AliRsnCut (const char *name, const char *title, ERsnCutType type, Double_t min, Double_t max = 1e-100);
    AliRsnCut (const char *name, const char *title, ERsnCutType type, Int_t min, Int_t max = 32767);
    AliRsnCut (const char *name, const char *title, ERsnCutType type, UInt_t min, UInt_t max = 65534);

    ~AliRsnCut();

    void      SetCutValues (ERsnCutType type, const Double_t& theValue, const Double_t& theValue2);
    void      SetCutValues (ERsnCutType type, const Int_t& theValue, const Int_t& theValue2);
    void      SetCutValues (ERsnCutType type, const UInt_t& theValue, const UInt_t& theValue2);

    Bool_t    IsSelected (ECutSetType type,  AliRsnDaughter *daughter);
    Bool_t    IsSelected (ECutSetType type,  AliRsnPairParticle *pair);

    void      PrintAllValues();

    Bool_t    IsBetween (const Double_t &theValue);
    Bool_t    IsEqual (const Int_t &theValue);
    Bool_t    IsEqual (const UInt_t &theValue);
    Bool_t    IsEqual (const Double_t &theValue);

private:

    Bool_t CheckRestMomentum(AliRsnPairParticle *pair);

    Double_t        fDMin;          // min. double value
    Double_t        fDMax;          // max. double value
    Int_t           fIMin;          // min. int value
    Int_t           fIMax;          // max. int value
    UInt_t          fUIMin;         // min. uint value
    UInt_t          fUIMax;         // max. uint value

    ERsnCutType     fRsnCutType;    // cut type
    ERsnCutVarType  fRsnCutVarType; // variable type

    static const Double_t fgkDSmallNumber;  // small double value
    static const Double_t fgkDBigNumber;    // big double value
    static const Int_t    fgkIBigNumber;    // big int value

    ClassDef (AliRsnCut, 1)
};

#endif
