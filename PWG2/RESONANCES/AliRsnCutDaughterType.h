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

#ifndef ALIRSNCUTDAUGHTERTYPE_H
#define ALIRSNCUTDAUGHTERTYPE_H

#include "AliRsnCut.h"

class AliRsnCutDaughterType : public AliRsnCut
{
  public:
  
    enum EType
    {
      kTrackTPC,
      kTrackITSSA,
      kV0,
      kTypes
    };

    AliRsnCutDaughterType();
    AliRsnCutDaughterType(const char *name, EType type);
    virtual ~AliRsnCutDaughterType() {;};

    virtual Bool_t IsSelected(TObject *object);

  protected:
  
    EType fRefType;   // type to which the track format is compared

    ClassDef(AliRsnCutDaughterType, 1)
};

#endif
