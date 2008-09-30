//
// Class AliRsnPIDIndex
//
// It sorts the indexes of all tracks in an AliRsnEvent
// for a fast retrieval of them according to charge and PID.
//
// author: M. Vala (email: martin.vala@cern.ch)
//

#ifndef AliRsnPIDIndex_h
#define AliRsnPIDIndex_h

#include <TArrayI.h>
#include "AliRsnPID.h"

class AliRsnPIDIndex : public TObject
{
  public:

    AliRsnPIDIndex(Int_t num = 1000);
    AliRsnPIDIndex(const AliRsnPIDIndex &copy);
    AliRsnPIDIndex& operator= (const AliRsnPIDIndex& copy);

    virtual ~AliRsnPIDIndex();

    void     Print(Option_t *option = "") const;

    void     AddIndex(const Int_t index, Char_t sign, AliRsnPID::EType type);
    void     AddIndex(const Int_t index, Short_t sign, Int_t type);
    void     SetCorrectIndexSize();

    TArrayI* GetTracksArray(Char_t sign, AliRsnPID::EType type);
    TArrayI* GetCharged(Char_t sign);

  private:

    Int_t    ChargeIndex(Char_t sign) const;

    TArrayI  fIndex[2][AliRsnPID::kSpecies+1];       // index arrays of pos/neg particles of each PID
    Int_t    fNumOfIndex[2][AliRsnPID::kSpecies+1];  //! array size

    ClassDef(AliRsnPIDIndex, 1);
};

#endif
