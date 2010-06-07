//
// Class AliRsnValue
//
// Definition of a single value which can be computed
// from any of the defined input objects implemented
// in the resonance package.
//

#ifndef ALIRSNVALUE_H
#define ALIRSNVALUE_H

class AliRsnPairDef;

class AliRsnValue : public TObject
{
  public:

    enum EAxisType
    {
      kTrackPt,
      kTrackEta,
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
      kEventMult,
      kAxisTypes
    };

    enum EAxisObject
    {
      kParticle,
      kPair,
      kEvent,
      kNone
    };

    AliRsnValue(EAxisType type = kAxisTypes);
    virtual ~AliRsnValue() { }

    virtual const char* GetName() const;
    EAxisObject         GetAxisObject() const;
    void                SetType(EAxisType type) {fType = type;}
    EAxisType           GetAxisType() {return fType;}
    Double_t            Eval(TObject * const obj, const AliRsnPairDef *pairDef = 0x0) const;

  private:

    EAxisType fType;    // value type

    // ROOT dictionary
    ClassDef(AliRsnValue, 1)
};

#endif
