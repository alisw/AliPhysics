//
// Class AliRsnValue
//
// This class implements all the computations which could be useful
// during the analysis, both for cuts and for output histograms.
// 
// It inherits from the AliRsnTarget base class since it can operate
// on tracks, pairs and events, and the kind of expected object to
// be processed depends on the kind of requested computation.
//
// Since this class is used to produce the outputs, it contains the
// facilities to define a binning in an output histogram.
//

#ifndef ALIRSNVALUE_H
#define ALIRSNVALUE_H

#include "TArrayD.h"
#include "AliRsnTarget.h"

class AliRsnValue : public AliRsnTarget
{
  public:
  
    // this enumeration lists all available computations
    // any user feedback proposing new ones is welcome
    enum EValueType
    {
      kTrackP,             // single track total momentum
      kTrackPt,            // single track transverse momentum
      kTrackEta,           // single track pseudo-rapidity
      kPairP1,             // total momentum of 1st daughter of a pair
      kPairP2,             // total momentum of 2nd daughter of a pair
      kPairP1t,            // total momentum of 1st daughter of a pair
      kPairP2t,            // total momentum of 2nd daughter of a pair
      kPairP1z,            // total momentum of 1st daughter of a pair
      kPairP2z,            // total momentum of 2nd daughter of a pair
      kPairInvMass,        // pair invariant mass (with reconstructed momenta)
      kPairInvMassMC,      // pair invariant mass (with MC momenta)
      kPairInvMassRes,     // pair invariant mass resolution
      kPairPt,             // pair transverse momentum
      kPairPz,             // pair longitudinal momentum
      kPairEta,            // pair pseudo-rapidity
      kPairMt,             // pair transverse mass (need a reference mass)
      kPairY,              // pair rapidity (need a reference mass)
      kPairPhi,            // pair azimuthal angle (with reconstructed momenta)
      kPairPhiMC,          // pair azimuthal angle (with MC momenta)
      kPairPtRatio,        // ratio |pt1 - pt2|/(pt1 + pt2) of daughter transverse momenta
      kPairDipAngle,       // inverse cosine of the angle between daughter vector momenta
      kPairCosThetaStar,   // polarization angle
      kPairQInv,           // invariant relative momentum of the two daughters
      kPairAngleToLeading, // angle between the pair momentum and that of the event leading particle
      kEventLeadingPt,     // transverse momentum of the event leading particle
      kEventMult,          // multiplicity computed as the number of tracks
      kEventMultESDCuts,   // multiplicity computed as the number of track passing an ESD quality cut (need this cut defined)
      kEventVz,            // Z position of event primary vertex
      
      kValueTypes          // last value is used to have a meaningless enum value for initializations
    };

    AliRsnValue();
    AliRsnValue(const char *name, EValueType type, Int_t nbins = 0, Double_t min = 0.0, Double_t max = 0.0);
    AliRsnValue(const char *name, EValueType type, Double_t min, Double_t max, Double_t step);
    AliRsnValue(const char *name, EValueType type, Int_t nbins, Double_t *array);
    AliRsnValue(const AliRsnValue& copy);
    AliRsnValue& operator=(const AliRsnValue& copy);
    virtual ~AliRsnValue() { /*does nothing, since pointers are not owned by this object*/ }
    
    TArrayD     GetArray() const               {return fBinArray;}
    Double_t    GetComputedValue() const       {return fComputedValue;}
    EValueType  GetValueType() const           {return fValueType;}
    const char* GetValueTypeName() const;
    TObject*    GetSupportObject()             {return fSupportObject;}
    void        SetSupportObject(TObject *obj) {fSupportObject = obj;}
    void        SetValueType(EValueType type)  {fValueType = type;}
    void        AssignTarget();
    
    void        SetBins(Int_t n, Double_t min, Double_t max);
    void        SetBins(Int_t n, Double_t *array);
    void        SetBins(Double_t min, Double_t max, Double_t step);
    
    void        Set(EValueType type, Int_t n, Double_t min, Double_t max)       {fValueType = type; AssignTarget(); SetBins(n, min, max);}
    void        Set(EValueType type, Int_t n, Double_t *array)                  {fValueType = type; AssignTarget(); SetBins(n, array);}
    void        Set(EValueType type, Double_t min, Double_t max, Double_t step) {fValueType = type; AssignTarget(); SetBins(min, max, step);}
    
    
    virtual Bool_t  Eval(TObject *object, Bool_t useMC = kFALSE);
    virtual void    Print(Option_t *option = "") const;

  protected:
  
    Double_t        fComputedValue;  // computed value
    EValueType      fValueType;      // value type
    TArrayD         fBinArray;       // array of bins (when used for a histogram axis)
    TObject        *fSupportObject;  // support object needed for computing some of the values
    
    // ROOT dictionary
    ClassDef(AliRsnValue, 2)
};

#endif
