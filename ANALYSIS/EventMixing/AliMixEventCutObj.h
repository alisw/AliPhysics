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
class AliMixEventCutObj : public TObject
{
public:
  enum EEPAxis_t {kMultiplicity=0, kZVertex=1, kNumberV0s=2, kNumberTracklets=3, kAllEventAxis=4};
    
    AliMixEventCutObj(EEPAxis_t type=kMultiplicity,Float_t min=0,Float_t max=0, Float_t step=0);
    AliMixEventCutObj(const AliMixEventCutObj& obj);
    virtual ~AliMixEventCutObj() {;}

    virtual void Print(const Option_t *) const;
    void PrintCurrentInterval();
    void Reset();
    void AddStep();

    Bool_t      HasMore() const;

    Int_t       GetNumberOfBins() const;
    Float_t     GetMin() const { return fCurrentVal;}
    Float_t     GetMax() const { return fCurrentVal+fCutStep-fCutSmallVal;}
    Float_t     GetStep() const { return fCutStep;}
    Short_t     GetType() const {return fCutType;}
    Int_t       GetBinNumber(Float_t num) const;
    Int_t       GetIndex(AliVEvent*ev);
    Int_t       GetIndex(AliESDEvent*ev);
    Int_t       GetIndex(AliAODEvent*ev);
    const char *GetNameOfCut(Short_t index) const;
    
private:
    Short_t     fCutType;       // cut type
    Float_t     fCutMin;        // cut min
    Float_t     fCutMax;        // cut max
    Float_t     fCutStep;       // cut step
    Float_t     fCutSmallVal;   // small value

    Float_t     fCurrentVal;    // current value
    Bool_t      fNoMore;        // flag for no more bins

    AliMixEventCutObj& operator=(const AliMixEventCutObj& ) { return *this;}

    ClassDef(AliMixEventCutObj, 1)
};

#endif
