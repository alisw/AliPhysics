//
// Class AliMixEventCutObj
//
// AliMixEventCutObj object contains information about one cut on for event mixing
// used by AliMixEventPool class
//
// authors: 
//          Martin Vala (martin.vala@cern.ch)
//

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMultiplicity.h"

#include "AliMixEventCutObj.h"

ClassImp(AliMixEventCutObj)

//_________________________________________________________________________________________________
AliMixEventCutObj::AliMixEventCutObj(EEPAxis_t type, Float_t min, Float_t max, Float_t step) : TObject(),
    fCutType((Short_t)type),
    fCutMin(min),
    fCutMax(max),
    fCutStep(step),
    fCutSmallVal((max - min) / (step*1e6)),
    fCurrentVal(min),
    fNoMore(kFALSE) {
  //
  // Default constructor.
  //

  AliDebug(AliLog::kDebug + 2, "<-");

  AliDebug(AliLog::kDebug + 2, "->");
}
//_________________________________________________________________________________________________
AliMixEventCutObj::AliMixEventCutObj(const AliMixEventCutObj& obj) : TObject(obj),
    fCutType(obj.fCutType),
    fCutMin(obj.fCutMin),
    fCutMax(obj.fCutMax),
    fCutStep(obj.fCutStep),
    fCutSmallVal(obj.fCutSmallVal),
    fCurrentVal(obj.fCurrentVal),
    fNoMore(obj.fNoMore) {
  //
  // Copy constructor.
  //
}

//_________________________________________________________________________________________________
void AliMixEventCutObj::Reset() {
  //
  // Reset cut.
  //

  AliDebug(AliLog::kDebug + 2, "<-");
  fCurrentVal = fCutMin - fCutStep;
  fNoMore = kFALSE;
  AliDebug(AliLog::kDebug + 2, "->");
}
//_________________________________________________________________________________________________
Bool_t AliMixEventCutObj::HasMore() const {
  //
  // Return kTRUE when fCurrentVal is in interval of cut range
  //
  
  return ((fCurrentVal + fCutStep) < fCutMax);
}

//_________________________________________________________________________________________________
void AliMixEventCutObj::AddStep() {
  //
  // Adds step
  //
  
  fCurrentVal += fCutStep;
}

//_________________________________________________________________________________________________
void AliMixEventCutObj::Print(const Option_t *) const {
  //
  // Prints cut information
  //
  
  AliInfo(Form("%d %f %f %f", fCutType, fCutMin, fCutMax, fCutStep));
}
//_________________________________________________________________________________________________
void AliMixEventCutObj::PrintCurrentInterval() {
  //
  // Prints current cut interval information
  //
  
  AliInfo(Form("%s <%f,%f>", GetNameOfCut(fCutType), GetMin(), GetMax()));
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetNumberOfBins() const {
  //
  // Returns number of bins
  //
  
  if (!fCutStep) return -1;
  return (Int_t)((fCutMax -fCutMin) / fCutStep);
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetBinNumber(Float_t num) const {
  //
  // Returns bin (index) number in current cut.
  // Returns -1 in case of out of range
  //
  
  Int_t binNum = 0;
  for (Float_t i = fCutMin;i <= fCutMax - fCutSmallVal;i += fCutStep) {
    binNum++;
    if ((num >= i) && (num <= i + fCutStep - fCutSmallVal)) return binNum;

  }
  return -1;
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetIndex(AliVEvent* ev) {
  //
  // Finds bin (index) in current cut from event information.
  //
  
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(ev);
  if (esd) return GetIndex(esd);
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(ev);
  if (aod) return GetIndex(aod);
  return -1;
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetIndex(AliESDEvent* ev) {
  //
  // Finds bin (index) in current cut from ESD event information.
  //
  switch (fCutType) {
  case kMultiplicity:
    return GetBinNumber((Float_t)ev->GetNumberOfTracks());
  case kZVertex:
    return GetBinNumber(ev->GetVertex()->GetZ());
  case kNumberV0s:
    return GetBinNumber(ev->GetNumberOfV0s());
  case kNumberTracklets:
    const AliMultiplicity* multESD = ev->GetMultiplicity();
    return GetBinNumber(multESD->GetNumberOfTracklets());
  }
  return -1;
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetIndex(AliAODEvent* ev) {
  //
  // Finds bin (index) in current cut from AOD event information.
  //
  switch (fCutType) {
  case kMultiplicity:
    return GetBinNumber((Float_t)ev->GetNumberOfTracks());
  case kZVertex:
    return GetBinNumber(ev->GetVertex(0)->GetZ());
  case kNumberV0s:
    return GetBinNumber(ev->GetNumberOfV0s());
  }
  return -1;
}
//_________________________________________________________________________________________________
const char* AliMixEventCutObj::GetNameOfCut(Short_t index) const
{
  //
  // Retruns name of cut
  //
  switch ((Int_t)index) {
    case kMultiplicity:
      return "Multiplicity";
    case kZVertex:
      return "ZVertex";
    case kNumberV0s:
      return "NumberV0s";
    case kNumberTracklets:
      return "NumberTracklets";
  }
  return "";
}

