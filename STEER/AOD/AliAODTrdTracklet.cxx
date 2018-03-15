#include "AliAODTrdTracklet.h"
AliAODTrdTracklet::AliAODTrdTracklet() :
  AliVTrdTracklet(),
  fHCId(-1),
  fTrackletWord(0),
  fLabel(-1)
{
  /// default constructor

}

AliAODTrdTracklet::AliAODTrdTracklet(const AliVTrdTracklet &rhs) :
  AliVTrdTracklet(rhs),
  fHCId(rhs.GetHCId()),
  fTrackletWord(rhs.GetTrackletWord()),
  fLabel(rhs.GetLabel())
{
  /// default constructor

}

AliAODTrdTracklet::AliAODTrdTracklet(UInt_t trackletWord, Short_t hcid, Int_t label) :
  AliVTrdTracklet(),
  fHCId(hcid),
  fTrackletWord(trackletWord),
  fLabel(label)
{
  /// constructor

}

AliAODTrdTracklet::AliAODTrdTracklet(const AliAODTrdTracklet& rhs) :
  AliVTrdTracklet(rhs),
  fHCId(rhs.fHCId),
  fTrackletWord(rhs.fTrackletWord),
  fLabel(rhs.fLabel)
{
  /// copy constructor

}

AliAODTrdTracklet& AliAODTrdTracklet::operator=(const AliAODTrdTracklet& rhs)
{
  /// assignment operator

  if (&rhs != this) {
    AliVTrdTracklet::operator=(rhs);

    fHCId = rhs.fHCId;
    fTrackletWord = rhs.fTrackletWord;
    fLabel = rhs.fLabel;
  }

  return *this;
}

void AliAODTrdTracklet::Copy(TObject &rhs) const
{
  /// copy

  AliVTrdTracklet::Copy(rhs);
}

Int_t AliAODTrdTracklet::GetBinY() const
{
  /// returns (signed) value of Y

  if (fTrackletWord & 0x1000) {
    return -((~(fTrackletWord-1)) & 0x1fff);
  }
  else {
    return (fTrackletWord & 0x1fff);
  }
}

Int_t AliAODTrdTracklet::GetBinDy() const
{
  /// returns (signed) value of the deflection length

  if (fTrackletWord & (1 << 19)) {
    return -((~((fTrackletWord >> 13) - 1)) & 0x7f);
  }
  else {
    return ((fTrackletWord >> 13) & 0x7f);
  }
}

// Float_t AliAODTrdTracklet::GetDyDx() const
// {
//   // returns the deflection over 3 cm drift length

//   return GetBinDy() * fgkBinWidthDy/fgkDriftLength;
// }
