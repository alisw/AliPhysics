//
// Class AliRsnPIDIndex
//
// Sorts the indexes of all tracks in an AliRsnEvent
// for a fast retrieval of them according to charge and PID.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliLog.h"

#include "AliRsnEvent.h"
#include "AliRsnPIDIndex.h"

ClassImp(AliRsnPIDIndex)

//_____________________________________________________________________________
AliRsnPIDIndex::AliRsnPIDIndex(Int_t num)
{
//
// Default constructor.
// Reset all arrays to the size specified in the argument (default = 1000).
//

  Int_t imethod, icharge, ipid, i;

  for (imethod = 0; imethod < AliRsnDaughter::kMethods; imethod++)
  {
    for (icharge = 0; icharge < 2; icharge++)
    {
      for (ipid = 0; ipid <= AliPID::kSPECIES; ipid++)
      {
        fNumOfIndex[imethod][icharge][ipid] = 0;
        fIndex[imethod][icharge][ipid].Set(num);
        for (i = 0; i < num; i++) fIndex[imethod][icharge][ipid].AddAt(-1, i);
      }
    }
  }
}

//_____________________________________________________________________________
AliRsnPIDIndex::AliRsnPIDIndex(const AliRsnPIDIndex & copy) :
  TObject(copy)
{
//
// Copy constructor.
// Duplicates all arrays.
//

  Int_t imethod, icharge, ipid, i, size;

  for (imethod = 0; imethod < AliRsnDaughter::kMethods; imethod++)
  {
    for (icharge = 0; icharge < 2; icharge++)
    {
      for (ipid = 0; ipid <= AliPID::kSPECIES; ipid++)
      {
        fNumOfIndex[imethod][icharge][ipid] = copy.fNumOfIndex[imethod][icharge][ipid];
        size = copy.fIndex[imethod][icharge][ipid].GetSize();
        fIndex[imethod][icharge][ipid].Set(size);
        for (i = 0; i < size; i++) {
          fIndex[imethod][icharge][ipid].AddAt(copy.fIndex[imethod][icharge][ipid].At(i), i);
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
// Duplicates all arrays.
//

  Int_t imethod, icharge, ipid, k, size;

  for (imethod = 0; imethod < AliRsnDaughter::kMethods; imethod++)
  {
    for (icharge = 0; icharge < 2; icharge++)
    {
      for (ipid = 0; ipid <= AliPID::kSPECIES; ipid++)
      {
        fNumOfIndex[imethod][icharge][ipid] = copy.fNumOfIndex[imethod][icharge][ipid];
        size = copy.fIndex[imethod][icharge][ipid].GetSize();
        fIndex[imethod][icharge][ipid].Set(size);
        for (k = 0; k < size; k++) {
          fIndex[imethod][icharge][ipid].AddAt(copy.fIndex[imethod][icharge][ipid].At(k), k);
        }
      }
    }
  }

  // return this object
  return (*this);
}

//_____________________________________________________________________________
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
// Resets all arrays to a given size, and storing '-1' in all values
//

  Int_t imethod, icharge, ipid, i;

  for (imethod = 0; imethod < AliRsnDaughter::kMethods; imethod++)
  {
    for (icharge = 0; icharge < 2; icharge++)
    {
      for (ipid = 0; ipid <= AliPID::kSPECIES; ipid++)
      {
        fNumOfIndex[imethod][icharge][ipid] = 0;
        fIndex[imethod][icharge][ipid].Set(num);
        for (i = 0; i < num; i++) fIndex[imethod][icharge][ipid].AddAt(-1, i);
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
void AliRsnPIDIndex::AddIndex
(const Int_t index, AliRsnDaughter::EPIDMethod meth, Char_t sign, AliPID::EParticleType type)
{
//
// Adds one index (1st arg) to the corresponding array defined by:
// - PID method (2nd arg)
// - charge     (3rd arg)
// - pid type   (4th arg)
//

  Int_t iMethod = (Int_t)meth;
  Int_t iCharge = ChargeIndex(sign);
  Int_t iType   = (Int_t)type;

  // since AliPID kUnknown = 10 and kSPECIES = 5, the PID is "unknown"
  // a correction in the storage index must be done
  if (type == AliPID::kUnknown) iType = (Int_t)AliPID::kSPECIES;

  // debug message (can this be removed?)
  AliDebug(AliLog::kDebug+1,Form("Adding index=%d method=%d sign='%c', pidType=%d",index,meth,sign,type));

  // check indexes before storing
  if (iMethod < 0 || iMethod >= AliRsnDaughter::kMethods) return;
  if (iCharge < 0 || iCharge > 1) return;
  if (iType < 0 || iType > AliPID::kSPECIES) return;

  // debug to know if adding was successful
  AliDebug(AliLog::kDebug+1,Form("Adding succeeded for index=%d method=%d sign='%c', pidType=%d",index,meth,sign,type));

  // insert the index in the array and increment the counter of occupied slots
  fIndex[iMethod][iCharge][iType].AddAt(index, fNumOfIndex[iMethod][iCharge][iType]);
  fNumOfIndex[iMethod][iCharge][iType]++;
}

//_____________________________________________________________________________
void AliRsnPIDIndex::SetCorrectIndexSize()
{
//
// Removes in every array the unused slots, to compress its size.
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
TArrayI* AliRsnPIDIndex::GetTracksArray
(AliRsnDaughter::EPIDMethod meth, Char_t sign, AliPID::EParticleType type)
{
//
// Returns the array of indexes of tracks whose charge
// and PID correspond to the passed arguments:
//   1) PID method
//   2) sign of particle ('+' or '-')
//   3) PID of particle (from AliRsnPID::EType)
// Otherwise returns null pointer.
//

  Int_t iMethod = (Int_t)meth;
  Int_t iCharge = ChargeIndex(sign);
  Int_t iType   = (Int_t)type;

  // since AliPID kUnknown = 10 and kSPECIES = 5, the PID is "unknown"
  // a correction in the storage index must be done
  if (type == AliPID::kUnknown) iType = (Int_t)AliPID::kSPECIES;

  // check indexes
  if (iMethod < 0 || iMethod >= AliRsnDaughter::kMethods)
  {
    AliError("Unrecognized PID method");
    return 0x0;
  }
  if (iCharge < 0 || iCharge > 1)
  {
    AliError(Form("Charge index (%d) out of range", iCharge));
    return 0x0;
  }
  if (iType < 0 || iType > AliPID::kSPECIES)
  {
    AliError(Form("PID index(%d) out of range", iType));
    return 0x0;
  }

  return &fIndex[iMethod][iCharge][iType];
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
// Returns the array index corresponding to charge.
// The inverse operation (returning '+' or '-' depending on charge)
// is done by 'AliRsnDaughter::GetChargeC()' method.
//

  if (sign == '+') return 0;
  else if (sign == '-') return 1;
  else if (sign == '0') return 2;
  else {
    AliError(Form("Character '%c' not recognized as charge sign", sign));
    return -1;
  }
}

//_____________________________________________________________________________
void AliRsnPIDIndex::FillFromEvent(AliRsnEvent* const event)
{
//
// Scans a whole event and fills the arrays depending
// on the identification of all tracks with all methods.
//

  Int_t numOfTracks = event->GetMultiplicity();

  // first, reset all arrays to the multiplicity,
  // in order to avoid 'out of range' errors
  ResetAll(numOfTracks);

  // now loop on tracks using a unique AliRsnDaughter object
  // which will point to one track at a time, and allows
  // to recovery the informations to classify it
  Int_t i;
  AliRsnDaughter daughter;
  for (i = 0; i < numOfTracks; i++)
  {
    // link the cursor to current track
    // this method recoveries the true PID and
    // combines weights with prior probabilities (defined in AliRsnEvent)
    // to obtain realistic PID
    event->SetDaughter(daughter, i);

    // assign current index to right slot for NoPID (essentially the charge)
    AddIndex(i,AliRsnDaughter::kNoPID,daughter.ChargeC(),AliPID::kUnknown);

    // assign current index to right slot for realistic PID
    AddIndex(i, AliRsnDaughter::kRealistic, daughter.ChargeC(), daughter.RealisticPID());

    // if MC info is present, assign current index to right slot for true PID
    if (daughter.GetParticle())
      AddIndex(i, AliRsnDaughter::kPerfect, daughter.ChargeC(), daughter.PerfectPID());
  }

  // at the end, compress all arrays to the correct size
  SetCorrectIndexSize();
}
