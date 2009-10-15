//
// Class AliRsnFunctionAxis
//
// Histogram definition.
// Contains required informations to create a histogram
// with fixed bin size: number of bins, minimum and maximum.
// Variable bin sizes are not considered because they are
// not used as typical output of analysis in this package.
//

#ifndef ALIRSNFUNCTIONAXIS_H
#define ALIRSNFUNCTIONAXIS_H


class TArrayD;
class AliRsnPairDef;

class AliRsnFunctionAxis : public TObject
{
  public:

    enum EAxisType {
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

    enum EAxisObject {
      kParticle,
      kPair,
      kEvent,
      kNone
    };

    AliRsnFunctionAxis();
    AliRsnFunctionAxis(EAxisType type, Int_t n, Double_t min, Double_t max);
    AliRsnFunctionAxis(EAxisType type, Double_t min, Double_t max, Double_t step);
    virtual ~AliRsnFunctionAxis() { }

    virtual const char* GetName() const;

    Int_t       GetNBins() const {return fNBins;}
    Double_t    GetMin() const {return fMin;}
    Double_t    GetMax() const {return fMax;}
    TArrayD     GetArray() const;
    EAxisObject GetAxisObject() const;

    void     SetType(EAxisType type) {fType = type;}
    void     SetBins(Int_t n, Double_t min, Double_t max);
    void     SetBins(Double_t min, Double_t max, Double_t step);
    void     SetMass(Double_t mass) {fMass = mass;}

    Double_t Eval(AliRsnDaughter *daughter) const;
    Double_t Eval(AliRsnPairParticle*const pair, AliRsnPairDef*const pairDef) const;
    Double_t Eval(AliRsnEvent *const event) const;

  private:

    EAxisType fType;    // binning type

    Int_t     fNBins;   // number of bins
    Double_t  fMin;     // lower edge
    Double_t  fMax;     // upper edge
    Double_t  fMass;    // reference mass for Y and Mt bins

    // ROOT dictionary
    ClassDef(AliRsnFunctionAxis, 1)
};

#endif
