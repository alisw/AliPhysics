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
      kPairInvMass,
      kPairInvMassMC,
      kPairInvMassRes,
      kPairPt,
      kPairEta,
      kPairMt,
      kPairY,
      kPairCosThetaStar,
      kPairCosThetaStar1,
      kPairCosThetaStar2,
      kPairCosThetaStarMC1,
      kPairCosThetaStarMC2,
      kEventMult,
      kValueTypes
    };

    AliRsnValue();
    AliRsnValue(const char *name, EValueType type, Int_t n, Double_t min, Double_t max);
    AliRsnValue(const char *name, EValueType type, Double_t min, Double_t max, Double_t step);
    AliRsnValue(const AliRsnValue& copy) : TNamed(copy),fType(copy.fType),fNBins(copy.fNBins),fMin(copy.fMin),fMax(copy.fMax),fValue(copy.fValue) {}
    AliRsnValue& operator=(const AliRsnValue& copy) {SetName(copy.GetName());fType=copy.fType;fNBins=copy.fNBins;fMin=copy.fMin;fMax=copy.fMax;fValue=copy.fValue;return (*this);}
    virtual ~AliRsnValue() { }
    
    //const char* GetName() const;
    Int_t       GetNBins() const {return fNBins;}
    Double_t    GetMin() const {return fMin;}
    Double_t    GetMax() const {return fMax;}
    TArrayD     GetArray() const;
    Double_t    GetValue() const {return fValue;}
    EValueType  GetValueType() {return fType;}

    void        SetValueType(EValueType type) {fType = type;}
    void        SetBins(Int_t n, Double_t min, Double_t max);
    void        SetBins(Double_t min, Double_t max, Double_t step);
    
    Bool_t      Eval(AliRsnMother * const mother, AliRsnPairDef * const pairDef, AliRsnEvent * const event);

  private:

    EValueType fType;    // value type
    Int_t      fNBins;   // number of bins (when applicable)
    Double_t   fMin;     // lower edge
    Double_t   fMax;     // upper edge
    Double_t   fValue;   // computed value
    
    // ROOT dictionary
    ClassDef(AliRsnValue, 1)
};

#endif
