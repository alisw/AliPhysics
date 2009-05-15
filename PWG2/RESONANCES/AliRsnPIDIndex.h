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

#include "AliRsnEvent.h"
#include "AliAODTrack.h"
#include "AliRsnDaughter.h"

class AliESDtrackCuts;

class AliRsnPIDIndex : public TObject
{
  public:

    AliRsnPIDIndex(Int_t num = 1000);
    AliRsnPIDIndex(const AliRsnPIDIndex &copy);
    AliRsnPIDIndex& operator= (const AliRsnPIDIndex& copy);

    virtual  ~AliRsnPIDIndex();

    void      Print(Option_t *option = "") const;
    void      ResetAll(Int_t num=1000);
    
    void      FillFromEvent(AliRsnEvent *event = 0, AliESDtrackCuts *cuts = 0);

    void      AddIndex(const Int_t index, AliRsnDaughter::EPIDMethod meth, Char_t sign, AliPID::EParticleType type);
    void      SetCorrectIndexSize();

    TArrayI*  GetTracksArray(AliRsnDaughter::EPIDMethod meth, Char_t sign, AliPID::EParticleType type);
    TArrayI*  GetCharged(Char_t sign);
    
    // Prior probs
    void            SetPriorProbability(AliPID::EParticleType type, Double_t p);
    void            SetPriorProbability(Double_t *out);
    void            DumpPriors();
    void            GetPriorProbability(Double_t *out);

  private:

    Int_t     ChargeIndex(Char_t sign) const;
    Char_t    IndexCharge(Short_t sign) const;

    TArrayI   fIndex[AliRsnDaughter::kMethods][2][AliPID::kSPECIES + 1];       // index arrays of pos/neg particles of each PID
    Int_t     fNumOfIndex[AliRsnDaughter::kMethods][2][AliPID::kSPECIES + 1];  //! array size
    
    Double_t  fPrior[AliPID::kSPECIES]; // prior probabilities

    ClassDef(AliRsnPIDIndex, 1);
};

#endif
