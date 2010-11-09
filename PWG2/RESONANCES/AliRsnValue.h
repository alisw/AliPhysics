//
// Class AliRsnValue
//
// Definition of a single value which can be computed
// from any of the defined input objects implemented
// in the resonance package.
//

#ifndef ALIRSNVALUE_H
#define ALIRSNVALUE_H

#include "TNamed.h"
#include "TArrayD.h"
#include "AliESDtrackCuts.h"

class AliRsnPairDef;
class AliRsnMother;

class AliRsnValue : public TNamed
{
  public:

    enum EValueType
    {
      kTrack1P,
      kTrack2P,
      kTrack1Pt,
      kTrack2Pt,
      kTrack1Px,
      kTrack1Py,
      kTrack1Pz,
      kTrack2Px,
      kTrack2Py,
      kTrack2Pz,
      kPairInvMass,
      kPairInvMassMC,
      kPairInvMassRes,
      kPairPt,
      kPairEta,
      kPairMt,
      kPairY,
      kPairPhi,
      kPairPhiMC,
      kPairPtRatio,
      kPairDipAngle,
      kPairCosThetaStar,
      kAngleToLeading,
      kLeadingPt,
      kQInv,
      kEventMult,
      kEventMultESDCuts,
      kValueTypes
    };

    AliRsnValue();
    AliRsnValue(const char *name, EValueType type, Int_t n = 0, Double_t min = 0.0, Double_t max = 0.0);
    AliRsnValue(const char *name, EValueType type, Double_t min, Double_t max, Double_t step);
    AliRsnValue(const char *name, EValueType type, Int_t n, Double_t *array);
    AliRsnValue(const AliRsnValue& copy);
    AliRsnValue& operator=(const AliRsnValue& copy);
    virtual ~AliRsnValue() { }
    
    TArrayD                GetArray() const {return fArray;}
    Double_t               GetValue() const {return fValue;}
    EValueType             GetValueType() const {return fType;}
    const AliESDtrackCuts* GetCuts() const {return &fESDCuts;}

    void        SetValueType(EValueType type) {fType = type;}
    void        SetBins(Int_t n, Double_t min, Double_t max);
    void        SetBins(Double_t min, Double_t max, Double_t step);
    void        SetBins(Int_t n, Double_t *array);
    void        Set(EValueType type, Int_t n = 0, Double_t min = 0.0, Double_t max = 0.0) {fType = type; SetBins(n, min, max);}
    void        Set(EValueType type, Double_t min, Double_t max, Double_t step) {fType = type; SetBins(min, max, step);}
    void        Set(EValueType type, Int_t n, Double_t *array) {fType = type; SetBins(n, array);}
    
    virtual Bool_t  Eval(AliRsnMother *mother, AliRsnPairDef *pairDef, AliRsnEvent *event);
    virtual Bool_t  Eval(AliRsnDaughter *daughter, AliRsnEvent *event);
    virtual void    Print(Option_t *option = "") const;

  protected:
  
    Double_t        fValue;   // computed value
    EValueType      fType;    // value type
    TArrayD         fArray;   // array of bins (when necessary)
    AliESDtrackCuts fESDCuts; // ESD track cuts used for a way to compute multiplicity
    
    // ROOT dictionary
    ClassDef(AliRsnValue, 1)
};

#endif
