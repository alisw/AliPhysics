//
// Class AliRsnPIDIndex
//
// It sorts the indexes of all tracks in an AliRsnEvent
// for a fast retrieval of them according to charge and PID.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNPIDINDEX_H
#define ALIRSNPIDINDEX_H

#include <TArrayI.h>

#include "AliRsnDaughter.h"
#include "AliRsnPIDDefESD.h"

class AliAODTrack;
class AliRsnEvent;

class AliRsnPIDIndex : public TObject
{
  public:

    AliRsnPIDIndex(Int_t num = 1000);
    AliRsnPIDIndex(const AliRsnPIDIndex &copy);
    AliRsnPIDIndex& operator= (const AliRsnPIDIndex& copy);

    virtual  ~AliRsnPIDIndex();

    void      Print(Option_t *option = "") const;
    void      ResetAll(Int_t num = 1000);

    void      FillFromEvent(AliRsnEvent*const event = 0);

    void      AddIndex(const Int_t index, AliRsnDaughter::EPIDMethod meth, Char_t sign, AliPID::EParticleType type);
    void      SetCorrectIndexSize();

    TArrayI*  GetTracksArray(AliRsnDaughter::EPIDMethod meth, Char_t sign, AliPID::EParticleType type);
    TArrayI*  GetCharged(Char_t sign);

  private:

    Int_t     ChargeIndex(Char_t sign) const;
    Char_t    IndexCharge(Short_t sign) const;

    TArrayI   fIndex[AliRsnDaughter::kMethods][2][AliPID::kSPECIES + 1];       // index arrays of pos/neg particles of each PID
    Int_t     fNumOfIndex[AliRsnDaughter::kMethods][2][AliPID::kSPECIES + 1];  //! array size

    ClassDef(AliRsnPIDIndex, 1);
};

#endif
