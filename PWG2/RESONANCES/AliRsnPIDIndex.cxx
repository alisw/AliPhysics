//
// Class AliRsnPIDIndex
//
// It sorts the indexes of all tracks in an AliRsnEvent
// for a fast retrieval of them according to charge and PID.
//
// author: M. Vala (email: martin.vala@cern.ch)
//

#include <TObject.h>

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliRsnPIDIndex.h"

ClassImp(AliRsnPIDIndex)

//_____________________________________________________________________________
AliRsnPIDIndex::AliRsnPIDIndex(Int_t num)
{
//
// Default constructor
//
  Int_t i, j, k, l;

  for (l = 0; l < AliRsnDaughter::kMethods; l++) {
    for (i = 0; i < 2; i++) {
      for (j = 0; j <= AliPID::kSPECIES; j++) {
        fNumOfIndex[l][i][j] = 0;
        fIndex[l][i][j].Set(num);
        for (k = 0; k < num; k++) fIndex[l][i][j].AddAt(-1, k);
      }
    }
  }
}

//_____________________________________________________________________________
AliRsnPIDIndex::AliRsnPIDIndex(const AliRsnPIDIndex & copy)
    : TObject(copy)
{
//
// Copy constructor.
// Creates new instances of all collections
// to store a copy of all objects.
//

  Int_t l, i, j, k, size;

  for (l = 0; l < AliRsnDaughter::kMethods; l++) {
    for (i = 0; i < 2; i++) {
      for (j = 0; j <= AliPID::kSPECIES; j++) {
        fNumOfIndex[l][i][j] = copy.fNumOfIndex[l][i][j];
        size = copy.fIndex[l][i][j].GetSize();
        fIndex[l][i][j].Set(size);
        for (k = 0; k < size; k++) {
          fIndex[l][i][j].AddAt(copy.fIndex[l][i][j].At(k), k);
        }
      }
    }
  }
}

//_____________________________________________________________________________
AliRsnPIDIndex& AliRsnPIDIndex::operator= (const AliRsnPIDIndex & copy)
{
//
// Assignment operator.
// Creates new instances of all collections
// to store a copy of all objects.
//

  Int_t l, i, j, k, size;

  for (l = 0; l < AliRsnDaughter::kMethods; l++) {
    for (i = 0; i < 2; i++) {
      for (j = 0; j < AliPID::kSPECIES; j++) {
        fNumOfIndex[l][i][j] = copy.fNumOfIndex[l][i][j];
        size = copy.fIndex[l][i][j].GetSize();
        fIndex[l][i][j].Set(size);
        for (k = 0; k < size; k++) {
          fIndex[l][i][j].AddAt(copy.fIndex[l][i][j].At(k), k);
        }
      }
    }
  }

  // return this object
  return (*this);
}

AliRsnPIDIndex::~AliRsnPIDIndex()
{
//
// Destructor.
// Does nothing.
//
}
//_____________________________________________________________________________
void AliRsnPIDIndex::ResetAll(Int_t num)
{
//
// Resets All
//

  Int_t i, j, k, l;

  for (l = 0; l < AliRsnDaughter::kMethods; l++) {
    for (i = 0; i < 2; i++) {
      for (j = 0; j <= AliPID::kSPECIES; j++) {
        fNumOfIndex[l][i][j] = 0;
        fIndex[l][i][j].Set(num);
        for (k = 0; k < num; k++) fIndex[l][i][j].AddAt(-1, k);
      }
    }
  }
}

//_____________________________________________________________________________
void AliRsnPIDIndex::Print(Option_t* /*option*/) const
{
//
// Prints AliRsnPIDIndex info
//
  Int_t i, j, l;

  for (l = 0; l < AliRsnDaughter::kMethods; l++) {
    for (i = 0; i < 2; i++) {
      for (j = 0; j <= AliPID::kSPECIES; j++) {
        AliInfo(Form(" [%d][%d][%d] %d %d",l, i, j, fIndex[l][i][j].GetSize(), fNumOfIndex[i][j]));
      }
    }
  }
}

//_____________________________________________________________________________
void AliRsnPIDIndex::AddIndex(const Int_t index, AliRsnDaughter::EPIDMethod meth, Char_t sign, AliPID::EParticleType  type)
{
//
// Adds index to corresponding TArrayI
//

  Int_t iMethod = (Int_t)meth;
  Int_t iCharge = ChargeIndex(sign);
  Int_t iType = (Int_t)type;
  // in AliPID kUnknown = 10 and kSPECIES = 5!
  if (type == AliPID::kUnknown) iType = (Int_t)AliPID::kSPECIES;
  AliDebug(AliLog::kDebug+1,Form("Adding index=%d method=%d sign='%c', pidType=%d",index,meth,sign,type));
  fIndex[iMethod][iCharge][iType].AddAt(index, fNumOfIndex[iMethod][iCharge][iType]);
  fNumOfIndex[iMethod][iCharge][iType]++;

}

//_____________________________________________________________________________
void AliRsnPIDIndex::SetCorrectIndexSize()
{
//
// Sets Correct sizes to all TArrayI
//

  Int_t i, j, l;
  for (l = 0; l < AliRsnDaughter::kMethods; l++) {
    for (i = 0; i < 2; i++) {
      for (j = 0; j <= AliPID::kSPECIES; j++) {
        fIndex[l][i][j].Set(fNumOfIndex[l][i][j]);
      }
    }
  }
}

//_____________________________________________________________________________
TArrayI* AliRsnPIDIndex::GetTracksArray(AliRsnDaughter::EPIDMethod meth, Char_t sign, AliPID::EParticleType type)
{
//
// Returns the array of indexes of tracks whose charge
// and PID correspond to the passed arguments:
//   1) sign of particle ('+' or '-')
//   2) PID of particle (from AliRsnPID::EType)
// Otherwise returns null pointer.
//

  Int_t iMethod = (Int_t)meth;
  Int_t icharge = ChargeIndex(sign);
  Int_t itype   = 0;
  if (icharge < 0) return (TArrayI *) 0x0;
  if (type == AliPID::kUnknown) itype = (Int_t)AliPID::kSPECIES; else itype = (Int_t)type;
  if (itype < 0 || itype > (Int_t)AliPID::kSPECIES) {
    AliError(Form("Index %d out of range", itype));
    return (TArrayI *) 0x0;
  }

  return &fIndex[iMethod][icharge][itype];
}

//_____________________________________________________________________________
TArrayI* AliRsnPIDIndex::GetCharged(Char_t sign)
{
//
// Returns the array of indexes of tracks whose charge
// corresponds to the passed argument
// Otherwise returns a null pointer.
//
  return GetTracksArray(AliRsnDaughter::kNoPID, sign, AliPID::kUnknown);

}

//_____________________________________________________________________________
Int_t AliRsnPIDIndex::ChargeIndex(Char_t sign) const
{
//
// Returns the array index corresponding to charge
//

  if (sign == '+') return 0;
  else if (sign == '-') return 1;
  else {
    AliError(Form("Character '%c' not recognized as charge sign", sign));
    return -1;
  }
}

//_____________________________________________________________________________
Char_t AliRsnPIDIndex::IndexCharge(Short_t sign) const
{
  //
  // Returns the array index corresponding to charge
  //

  if (sign == 1) return '+';
  else if (sign == -1) return '-';
  else {
    AliError(Form("Charge is different the 0(+) and 1(-) value is '%d'...", sign));
    return '+';
  }
}

//_____________________________________________________________________________
void AliRsnPIDIndex::FillFromEvent(AliRsnEvent* event, AliESDtrackCuts *cuts)
{
//
// Fills indexes from event
//

  Int_t numOfTracks = event->GetMultiplicity();
  AliRsnDaughter daughter;
  Int_t i;
  for (i=0;i<numOfTracks;i++) {
    daughter = event->GetDaughter(i);
    // if ESD track cuts are specified, 
    // and if the reference is an ESD track
    // skip all tracks not passing the cut
    if (cuts) {
      AliESDtrack *track = daughter.GetRefESD();
      if (track) if (!cuts->IsSelected(track)) continue;
    }
    daughter.CombineWithPriors(fPrior);
//     daughter.Print("ALL");

    AddIndex(i,AliRsnDaughter::kNoPID,IndexCharge(daughter.Charge()),AliPID::kUnknown);
    AddIndex(i,AliRsnDaughter::kRealistic,IndexCharge(daughter.Charge()),(AliPID::EParticleType)daughter.RealisticPID());
    if (daughter.GetParticle())
      AddIndex(i,AliRsnDaughter::kPerfect,IndexCharge(daughter.Charge()),(AliPID::EParticleType)daughter.PerfectPID());
  }
}


//_____________________________________________________________________________
void AliRsnPIDIndex::SetPriorProbability(AliPID::EParticleType type, Double_t p)
{
  //
  // Sets the prior probability for Realistic PID, for a
  // given particle species.
  //

  if (type >= 0 && type < (Int_t)AliPID::kSPECIES) {
    fPrior[type] = p;
  }

}

//_____________________________________________________________________________
void AliRsnPIDIndex::DumpPriors()
{
  //
  // Print all prior probabilities
  //

  Int_t i;
  for (i = 0; i < AliPID::kSPECIES; i++) {
    AliInfo(Form("Prior probability for %10s = %3.5f", AliPID::ParticleName((AliPID::EParticleType)i), fPrior[i]));
  }
}

//_____________________________________________________________________________
void AliRsnPIDIndex::GetPriorProbability(Double_t *out)
{

  Int_t i;
  for (i=0; i<AliPID::kSPECIES; i++) {
    out[i] = fPrior[i];
  }

}

//_____________________________________________________________________________
void AliRsnPIDIndex::SetPriorProbability(Double_t *out)
{

  Int_t i;
  for (i=0;i<AliPID::kSPECIES;i++) {
    fPrior[i] = out[i];
  }

}

