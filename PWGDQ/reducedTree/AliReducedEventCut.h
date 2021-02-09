// Class for cutting on ALICE Var manager and other event specific information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
//   07/09/2016

#ifndef ALIREDUCEDEVENTCUT_H
#define ALIREDUCEDEVENTCUT_H

#include <vector>
#include "AliReducedVarCut.h"
#include "AliReducedVarManager.h"

//_________________________________________________________________________
class AliReducedEventCut : public AliReducedVarCut {

 public:
  AliReducedEventCut();
  AliReducedEventCut(const Char_t* name, const Char_t* title);
  virtual ~AliReducedEventCut();

  void AddEventTagFilter(ULong64_t filter, ULong64_t exclude = 0) {
    fEventFilter |= filter;
    fEventTagFilterEnabled = kTRUE;
    fEventTagFilterExclude |= exclude;
  };
  void AddEventTagFilterBit(UShort_t bit, Bool_t exclude = kFALSE)  {
    if(bit<64) {
      fEventFilter |= (ULong64_t(1)<<bit);
      fEventTagFilterEnabled = kTRUE;
      if(exclude) fEventTagFilterExclude |= (ULong64_t(1)<<bit);
    }
  };
  void EnableVertexDistanceCut() {
    fEventFilter |= (ULong64_t(1)<<13);
    fEventTagFilterEnabled = kTRUE;
  };
  
  void AddEventTriggerFilter(ULong64_t filter) {fEventTriggerMask |= filter; fEventTriggerMaskEnabled = kTRUE;};
  void AddEventTriggerFilterBit(UShort_t bit)  {if(bit>63) return; fEventTriggerMask |= (ULong64_t(1)<<bit); fEventTriggerMaskEnabled = kTRUE;};
  void AddEventTriggerClassFilter(TString triggerClass) {fEventTriggerClass.push_back(triggerClass); fEventTriggerClassEnabled = kTRUE;};
  void EventTriggerClassFilterUseOr() {fEventTriggerClassLogicalOr = kTRUE;};
  
  void AddEventL1InputFilter(UInt_t filter) {fEventL1Mask |= filter; fEventL1MaskEnabled = kTRUE;};
  void AddEventL1InputFilterBit(UShort_t bit)  {if(bit>31) return; fEventL1Mask |= (UInt_t(1)<<bit); fEventL1MaskEnabled = kTRUE;};
  void AddEventL0InputFilter(UInt_t filter) {fEventL0Mask |= filter; fEventL0MaskEnabled = kTRUE;};
  void AddEventL0InputFilterBit(UShort_t bit)  {if(bit>31) return; fEventL0Mask |= (UInt_t(1)<<bit); fEventL0MaskEnabled = kTRUE;};
  
  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(TObject* obj, Float_t* values);
  
 protected: 
      
  // Cuts on event specific quantities
  Bool_t                fEventTagFilterEnabled;       // if true apply the filter to the event tag
  ULong64_t             fEventTagFilterExclude;       // exclusion map for filter for the event tag map
  ULong64_t             fEventFilter;                 // filter for the event tag map
  Bool_t                fEventTriggerMaskEnabled;     // if true apply the filter on the event trigger mask
  ULong64_t             fEventTriggerMask;            // filter for the event trigger mask
  Bool_t                fEventTriggerClassEnabled;    // if true apply the filter on the event trigger class
  Bool_t                fEventTriggerClassLogicalOr;  // if true and >1 event trigger classes, use logical OR, default is AND
  std::vector<TString>  fEventTriggerClass;           // filter for the event trigger class
  Bool_t                fEventL1MaskEnabled;          // if true apply the filter on the event L1 trigger inputs mask
  ULong64_t             fEventL1Mask;                 // filter for the event L1 trigger inputs mask
  Bool_t                fEventL0MaskEnabled;          // if true apply the filter on the event L0 trigger inputs mask
  ULong64_t             fEventL0Mask;                 // filter for the event L0 trigger inputs mask
   
  AliReducedEventCut(const AliReducedEventCut &c);
  AliReducedEventCut& operator= (const AliReducedEventCut &c);
  
  ClassDef(AliReducedEventCut,6);
};

#endif
