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
class AliRsnPairDef;
class AliRsnEvent;

class AliRsnCut : public TNamed
{
 public:

  // available cut types
  // some ones work both for pairs and single tracks
  enum EType {
    kMomentum = 0,
    kTransMomentum,
    kEta,
    kRadialImpactParam,
    kRadialImpactParamMC,
    kMomentumMC,
    kTransMomentumMC,
    kEtaMC,
    kNSigma,
    kNSigmaCalculate,
    kStatus,
    kIsLabelEqual,
    kIsTruePair,
    kIsPrimary,
    kIsKinkDaughter,
    kIsKinkMother,
    kMCTracked,
    kChargePos,
    kChargeNeg,
    kPIDType,
    kPIDProb,
    kPIDProbForSpecies,
    kAssignedPID,
    kTruePID,
    kVz,
    kTrueMultiplicity,
    kMultiplicity,
    kMultiplicityDifference,
    kMultiplicityRatio,
    kPhiMeanDifference,
    kVzDifference,
    kLastCutType
  };

  // types of cut variables
  enum EVarType {
    kDouble_t = 0,
    kInt_t,
    kUInt_t
  };

  // possible targets for a cut
  enum ETarget {
    kParticle = 0,
    kPair,
    kEvent,
    kMixEvent,
    kLastCutTarget
  };

  AliRsnCut();
  AliRsnCut(const char *name, const char *title, EType type);
  AliRsnCut(const char *name, const char *title, EType type, Double_t min, Double_t max = 1e-100);
  AliRsnCut(const char *name, const char *title, EType type, Int_t min, Int_t max = 32767);
  AliRsnCut(const char *name, const char *title, EType type, UInt_t min, UInt_t max = 65534);
  AliRsnCut(const char *name, const char *title, EType type, ULong_t min, ULong_t max = 65534);

  ~AliRsnCut();

  void      SetCutValues(EType type, const Double_t& theValue, const Double_t& theValue2);
  void      SetCutValues(EType type, const Int_t& theValue, const Int_t& theValue2);
  void      SetCutValues(EType type, const UInt_t& theValue, const UInt_t& theValue2);
  void      SetCutValues(EType type, const ULong_t& theValue, const ULong_t& theValue2);

  Bool_t    IsSelected(ETarget tgt,  AliRsnDaughter *daughter);
  Bool_t    IsSelected(ETarget tgt,  AliRsnPairParticle *pair);
  Bool_t    IsSelected(ETarget tgt,  AliRsnEvent *event);
  Bool_t    IsSelected(ETarget tgt,  AliRsnEvent *ev1, AliRsnEvent *ev2);

  void      PrintAllValues();

  Bool_t    IsBetween(const Double_t &theValue);
  Bool_t    IsBetween(const Int_t &theValue);
  Bool_t    MatchesValue(const Int_t &theValue);
  Bool_t    MatchesValue(const UInt_t &theValue);
  Bool_t    MatchesValue(const ULong_t &theValue);
  Bool_t    MatchesValue(const Double_t &theValue);

 private:

  Double_t        fDMin;          // min. double value
  Double_t        fDMax;          // max. double value
  Int_t           fIMin;          // min. int value
  Int_t           fIMax;          // max. int value
  UInt_t          fUIMin;         // min. uint value
  UInt_t          fUIMax;         // max. uint value
  ULong_t         fULMin;         // min. ulong value
  ULong_t         fULMax;         // max. ulong value

  EType           fType;          // cut type
  EVarType        fVarType;       // variable type

  static const Double_t fgkDSmallNumber;  // small double value
  static const Double_t fgkDBigNumber;    // big double value
  static const Int_t    fgkIBigNumber;    // big int value

  ClassDef(AliRsnCut, 1)
};

#endif
