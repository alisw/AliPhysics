//
// Class AliMixInfo
//
// AliMixInfo object contains information about one cut on for event mixing
// available for users containing mixing information
//
// authors:
//          Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIMIXINFO_H
#define ALIMIXINFO_H

#include <TNamed.h>
#include <Rtypes.h>

class AliMixEventPool;
class TH1I;
class TList;
class TCollection;
class AliMixInfo : public TNamed {
public:
   enum EInfoHistorgramType { kMainEvents = 0, kMixedEvents = 1, kNumTypes };

   AliMixInfo(const char *name = "mix", const char *title = "MixInfo");
   AliMixInfo(const AliMixInfo &obj);
   virtual ~AliMixInfo();

   void Reset();
   virtual void Print(Option_t *option = "") const;
   virtual void Draw(Option_t *option = "");
   virtual Long64_t Merge(TCollection *list);

   void Add(AliMixInfo *mi);

   void SetOutputList(TList * const list) { fHistogramList = list; }
   void CreateHistogram(EInfoHistorgramType type, Int_t nbins, Int_t min, Int_t max);
   void FillHistogram(AliMixInfo::EInfoHistorgramType type, Int_t value);
   const char *GetNameHistogramByType(Int_t index) const;
   const char *GetTitleHistogramByType(Int_t index) const;
   TH1I  *GetHistogramByType(Int_t index) const;

   void    SetEventPool(AliMixEventPool *evPool);
   AliMixEventPool *GetEventPool(const char *name);


   static void DynamicExec(AliMixInfo*const mixInfo);
private:

   TList     *fHistogramList;  // histogram list

   AliMixInfo &operator=(const AliMixInfo &) { return *this; }

   ClassDef(AliMixInfo, 1)

};

#endif
