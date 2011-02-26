//
// Class AliMixEventCutObj
//
// AliMixEventCutObj object contains information about one cut on for event mixing
// used by AliMixEventPool class
//
// authors:
//          Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIMIXEVENTCUTOBJ_H
#define ALIMIXEVENTCUTOBJ_H

#include <TObject.h>
class AliVEvent;
class AliAODEvent;
class AliESDEvent;
class AliMixEventCutObj : public TObject {
public:
   enum EEPAxis_t {kMultiplicity = 0, kZVertex = 1, kNumberV0s = 2, kNumberTracklets = 3, kAllEventAxis = 4};

   AliMixEventCutObj(EEPAxis_t type = kMultiplicity, Float_t min = 0.0, Float_t max = 0.0, Float_t step = 1.0);
   AliMixEventCutObj(const AliMixEventCutObj &obj);
   AliMixEventCutObj &operator=(const AliMixEventCutObj &obj);

   virtual void Print(const Option_t *) const;
   void PrintCurrentInterval();
   void PrintValues(AliVEvent *main, AliVEvent* mix);
   void Reset();
   void AddStep();

   Bool_t      HasMore() const;

   Int_t       GetNumberOfBins() const;
   Float_t     GetMin() const { return fCurrentVal; }
   Float_t     GetMax() const { return fCurrentVal + fCutStep - fCutSmallVal; }
   Float_t     GetStep() const { return fCutStep; }
   Short_t     GetType() const { return fCutType; }
   Int_t       GetBinNumber(Float_t num) const;
   Int_t       GetIndex(AliVEvent *ev);
   Double_t    GetValue(AliVEvent *ev);
   Double_t    GetValue(AliESDEvent *ev);
   Double_t    GetValue(AliAODEvent *ev);

   const char *GetCutName(Int_t index = -1) const;

   void        SetCurrentValueToIndex(Int_t index);

private:
   Int_t       fCutType;       // cut type
   Float_t     fCutMin;        // cut min
   Float_t     fCutMax;        // cut max
   Float_t     fCutStep;       // cut step
   Float_t     fCutSmallVal;   // small value

   Float_t     fCurrentVal;    // current value

   ClassDef(AliMixEventCutObj, 1)
};

#endif
