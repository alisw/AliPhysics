//
// Class AliRsnCutRange
//
// General implementation of cuts which check a value inside a range.
// This range can be defined by two integers or two doubles.
// A user-friendly enumeration allows to define what is checked.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTSTD_H
#define ALIRSNCUTSTD_H

#include "AliRsnCut.h"

class AliRsnDaughter;
class AliRsnMother;
class AliRsnEvent;

class AliRsnCutStd : public AliRsnCut
{
  public:

    // available cut types
    // some ones work both for pairs and single tracks
    enum EType 
    {
      kP = 0,
      kPt,
      kEta,
      kY,
      kThetaDeg,
      kMult,
      kPtLeading,
      kDipAngle,
      // cut without reference values
      kCharge,
      kSameLabel,
      kPhysPrimary,
      // last
      kLastType
    };

    AliRsnCutStd();
    AliRsnCutStd(const char *name, ETarget target, EType type, Int_t    val1, Int_t    val2 = 0 , Bool_t useMC = kFALSE);
    AliRsnCutStd(const char *name, ETarget target, EType type, Double_t val1, Double_t val2 = 0., Bool_t useMC = kFALSE);
    virtual ~AliRsnCutStd() { }
    
    void           SetMass(Double_t mass) {fMass = mass;}
    EVarType       CheckType();
    
    virtual Bool_t IsSelected(TObject *obj1, TObject *obj2 = 0x0);

  protected:
  
    virtual Bool_t IsDaughterSelected(AliRsnDaughter *daughter);
    virtual Bool_t IsMotherSelected(AliRsnMother *mother);
    virtual Bool_t IsEventSelected(AliRsnEvent *event);

    EType     fType;       // cut type
    Bool_t    fUseMC;      // use or not MC values (when applicable)
    Double_t  fMass;       // mass hypothesis (used for Y and Mt)

    ClassDef(AliRsnCutStd, 1)
};

#endif
