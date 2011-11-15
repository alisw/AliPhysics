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
#include "AliAODVertex.h"

#include "AliMixEventCutObj.h"

ClassImp(AliMixEventCutObj)

//_________________________________________________________________________________________________
AliMixEventCutObj::AliMixEventCutObj(AliMixEventCutObj::EEPAxis_t type, Float_t min, Float_t max, Float_t step, const char *opt) : TObject(),
   fCutType((Int_t)type),
   fCutOpt(opt),
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
   fCutOpt(obj.fCutOpt),
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
AliMixEventCutObj &AliMixEventCutObj::operator=(const AliMixEventCutObj &obj)
{
   //
   // Assigned operator
   //
   if (&obj != this) {
      TObject::operator=(obj);
      fCutType = obj.fCutType;
      fCutOpt = obj.fCutOpt;
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
   for (Float_t iCurrent = fCutMin; iCurrent < fCutMax; iCurrent += fCutStep) {
      binNum++;
      if ((num >= iCurrent) && (num < iCurrent + fCutStep - fCutSmallVal)) {
         return binNum;
      }
   }
   return -1;
}

//_________________________________________________________________________________________________
Int_t AliMixEventCutObj::GetIndex(AliVEvent *ev)
{
   //
   // Finds bin (index) in current cut from event information.
   //
   return GetBinNumber(GetValue(ev));
}

//_________________________________________________________________________________________________
Double_t AliMixEventCutObj::GetValue(AliVEvent *ev)
{
   //
   // Returns value from event
   //

   AliESDEvent *esd = dynamic_cast<AliESDEvent *>(ev);
   if (esd) return GetValue(esd);
   AliAODEvent *aod = dynamic_cast<AliAODEvent *>(ev);
   if (aod) return GetValue(aod);

   AliFatal("Event is not supported in Event Mixing cuts!!!!");
   return -99999;
}

//_________________________________________________________________________________________________
Double_t AliMixEventCutObj::GetValue(AliESDEvent *ev)
{
   //
   // Returns value from esd event
   //

   const AliMultiplicity *multESD = 0;

   switch (fCutType) {
      case kMultiplicity:
         return (Double_t)ev->GetNumberOfTracks();
      case kZVertex:
         return ev->GetVertex()->GetZ();
      case kNumberV0s:
         return ev->GetNumberOfV0s();
      case kNumberTracklets:
         multESD = ev->GetMultiplicity();
         if (multESD) return multESD->GetNumberOfTracklets();
         else AliFatal("esd->GetMultiplicity() is null");
         break;
      case kCentrality:
         AliCentrality *c = ev->GetCentrality();
         if (!c) AliFatal("esd->GetCentrality() is null");
         if (fCutOpt.IsNull()) AliFatal("fCutOpt is null");
         return c->GetCentralityPercentile(fCutOpt.Data());

   }

   AliFatal("Mixing Cut TYPE is not supported for ESD");
   return -99999;

}

//_________________________________________________________________________________________________
Double_t AliMixEventCutObj::GetValue(AliAODEvent *ev)
{
   //
   // Returns value from aod event
   //

   AliAODVertex *v=0;
   switch (fCutType) {
      case kMultiplicity:
         return (Double_t) ev->GetNumberOfTracks();
      case kZVertex:
         v = ev->GetVertex(0);
         if (!v)  return -99999;
         return ev->GetVertex(0)->GetZ();
         // if verttex is null return -9999
         return -99999;
      case kNumberV0s:
         return ev->GetNumberOfV0s();
      case kCentrality:
         AliCentrality *c = ev->GetCentrality();
         if (!c) AliFatal("aod->GetCentrality() is null");
         if (fCutOpt.IsNull()) AliFatal("fCutOpt is null");
         return c->GetCentralityPercentile(fCutOpt.Data());
   }

   AliFatal("Mixing Cut TYPE is not supported for AOD");
   return -99999;
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
      case kCentrality:
         return Form("kCentrality[%s]", fCutOpt.Data());
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

//_________________________________________________________________________________________________
void AliMixEventCutObj::PrintValues(AliVEvent *main, AliVEvent *mix)
{
   //
   // Prints values of both events for current type
   //
   AliInfo(Form("name=%s main=%f mix=%f", GetCutName(), GetValue(main), GetValue(mix)));
}

