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
   fCutType((Int_t)type),
   fCutMin(min),
   fCutMax(max),
   fCutStep(step),
   fCutSmallVal(0),
   fCurrentVal(min)
{
   //
   // Default constructor
   //
   AliDebug(AliLog::kDebug + 5, "<-");
   if (fCutStep < 1e-5) AliError("fCutStep is too small !!! This cut will not work !!!");
   AliDebug(AliLog::kDebug + 5, "->");
}

//_________________________________________________________________________________________________
AliMixEventCutObj::AliMixEventCutObj(const AliMixEventCutObj &obj) : TObject(obj),
   fCutType(obj.fCutType),
   fCutMin(obj.fCutMin),
   fCutMax(obj.fCutMax),
   fCutStep(obj.fCutStep),
   fCutSmallVal(obj.fCutSmallVal),
   fCurrentVal(obj.fCurrentVal)
{
   //
   // Copy constructor
   //
   AliDebug(AliLog::kDebug + 5, "<-");
   AliDebug(AliLog::kDebug + 5, "->");
}

//_________________________________________________________________________________________________
AliMixEventCutObj& AliMixEventCutObj::operator=(const AliMixEventCutObj& obj)
{
   //
   // Assigned operator
   //
   if (&obj != this) {
      TObject::operator=(obj);
      fCutType = obj.fCutType;
      fCutMin = obj.fCutMin;
      fCutMax = obj.fCutMax;
      fCutStep = obj.fCutStep;
      fCutSmallVal = obj.fCutSmallVal;
      fCurrentVal = obj.fCurrentVal;
//       fNoMore = obj.fNoMore;
   }
   return *this;
}


//_________________________________________________________________________________________________
void AliMixEventCutObj::Reset()
{
   //
   // Reset cut
   //
   AliDebug(AliLog::kDebug + 5, "<-");
   fCurrentVal = fCutMin - fCutStep;
   AliDebug(AliLog::kDebug + 5, "->");
}
//_________________________________________________________________________________________________
Bool_t AliMixEventCutObj::HasMore() const
{
   //
   // Return kTRUE when fCurrentVal is in interval of cut range
   //
   return ((fCurrentVal + fCutStep) < fCutMax);
}

//_________________________________________________________________________________________________
void AliMixEventCutObj::AddStep()
{
   //
   // Adds step
   //
   fCurrentVal += fCutStep;
}

//_________________________________________________________________________________________________
void AliMixEventCutObj::Print(const Option_t *) const
{
   //
   // Prints cut information
   //
   AliInfo(Form("%s %f %f %f", GetCutName(fCutType), fCutMin, fCutMax, fCutStep));
}
//_________________________________________________________________________________________________
void AliMixEventCutObj::PrintCurrentInterval()
{
   //
   // Prints current cut interval information
   //
   AliDebug(AliLog::kDebug, Form("%s <%f,%f>", GetCutName(fCutType), GetMin(), GetMax()));
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetNumberOfBins() const
{
   //
   // Returns number of bins
   //
   if (fCutStep < 1e-5) return -1;
   return (Int_t)((fCutMax - fCutMin) / fCutStep);
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetBinNumber(Float_t num) const
{
   //
   // Returns bin (index) number in current cut.
   // Returns -1 in case of out of range
   //
   Int_t binNum = 0;
   for (Float_t i = fCutMin; i <= fCutMax - fCutSmallVal; i += fCutStep) {
      binNum++;
      if ((num >= i) && (num <= i + fCutStep - fCutSmallVal)) return binNum;
   }
   return -1;
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetIndex(AliVEvent *ev)
{
   //
   // Finds bin (index) in current cut from event information.
   //
   AliESDEvent *esd = dynamic_cast<AliESDEvent *>(ev);
   if (esd) return GetIndex(esd);
   AliAODEvent *aod = dynamic_cast<AliAODEvent *>(ev);
   if (aod) return GetIndex(aod);
   return -1;
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetIndex(AliESDEvent *ev)
{
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
         const AliMultiplicity *multESD = ev->GetMultiplicity();
         return GetBinNumber(multESD->GetNumberOfTracklets());
   }
   return -1;
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetIndex(AliAODEvent *ev)
{
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
const char *AliMixEventCutObj::GetCutName(Int_t index) const
{
   //
   // Retruns name of cut
   //

   if (index < 0) index = fCutType;
   switch (index) {
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

//_________________________________________________________________________________________________
void AliMixEventCutObj::SetCurrentValueToIndex(Int_t index)
{
   //
   // Sets current value to index
   //
   for (Int_t i = 0; i < index; i++) AddStep();
}
