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

class AliRsnCutStd : public AliRsnCut
{
 public:

  // available cut types
  // some ones work both for pairs and single tracks
  enum EType {
    kP = 0,
    kPt,
    kEta,
    kThetaDeg,
    kDr,
    kDz,
    kTPCsignal,
    kMult,
    kMultDiff,
    kMultDiffRel,
    kVzDiff,
    // value cuts
    kStatus,
    kKink,
    kKinkMother,
    kAssignedPID,
    kTruePID,
    // cut without reference values
    kCharge,
    kSameLabel,
    kTruePair,
    // last
    kLastType
  };

  AliRsnCutStd();
  AliRsnCutStd(const char *name, EType type, Int_t val1, Int_t val2 = 0, Bool_t useMC = kFALSE);
  AliRsnCutStd(const char *name, EType type, ULong_t val1, ULong_t val2 = 0, Bool_t useMC = kFALSE);
  AliRsnCutStd(const char *name, EType type, Double_t val1, Double_t val2 = 0.0, Bool_t useMC = kFALSE);
  virtual ~AliRsnCutStd() { }

  virtual Bool_t IsSelected(ETarget tgt, AliRsnDaughter *daughter);
  virtual Bool_t IsSelected(ETarget tgt, AliRsnPairParticle *pair);
  virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *event);
  virtual Bool_t IsSelected(ETarget tgt, AliRsnEvent *ev1, AliRsnEvent *ev2);

protected:

  EType     fType;       // cut type
  Bool_t    fUseMC;      // use or not MC values (when applicable)

  ClassDef(AliRsnCutStd, 1)
};

#endif
