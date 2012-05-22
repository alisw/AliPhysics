#ifndef ALIEMCALPHYSICSSELECTION_H
#define ALIEMCALPHYSICSSELECTION_H

// $Id$

#include "AliPhysicsSelection.h"

class AliEmcalPhysicsSelection: public AliPhysicsSelection
{
 public:
  enum EOfflineEmcalTypes { 
    kEmcalHC = BIT(28), //=true when EMCAL cell above given Et found
    kEmcalHT = BIT(29), //=true when EMCAL cluster above given Et found
    kEmcalOk = BIT(31), //=true when EMCAL good event criteria are found
  };

  AliEmcalPhysicsSelection();
  virtual ~AliEmcalPhysicsSelection() {;}

  virtual UInt_t GetSelectionMask(const TObject* obj);

  void           SetCellMinE(Double_t e)       { fCellMinE     = e; }
  void           SetClusMinE(Double_t e)       { fClusMinE     = e; }
  void           SetMarkFastOnly(Bool_t b)     { fMarkFastOnly = b; }
  void           SetMarkLedEvent(Bool_t b)     { fMarkLedEvent = b; }
  void           SetSkipFastOnly(Bool_t b)     { fSkipFastOnly = b; }
  void           SetSkipLedEvent(Bool_t b)     { fSkipLedEvent = b; }
  void           SetTrackMinPt(Double_t p)     { fTrackMinPt   = p; }
  void           SetTriggers(UInt_t t)         { fTriggers     = t; }

  Double_t       GetCellMaxE()   const { return fCellMaxE;    }
  Double_t       GetClusMaxE()   const { return fClusMaxE;    }
  Double_t       GetTrackMaxPt() const { return fTrackMaxPt;  }
  Bool_t         IsFastOnly()    const { return fIsFastOnly;  }
  Bool_t         IsLedEvent()    const { return fIsLedEvent;  }
  Bool_t         IsGoodEvent()   const { return fIsGoodEvent; }

 protected:
  Bool_t         fMarkFastOnly;     //=true then mark FastOnly events (only for LHC11a)
  Bool_t         fMarkLedEvent;     //=true then mark Led events (only for LHC11a)
  Bool_t         fSkipFastOnly;     //=true then skip FastOnly events (only for LHC11a)
  Bool_t         fSkipLedEvent;     //=true then skip Led events (only for LHC11a)
  Double_t       fCellMinE;         //minimum cell energy (<0 -> do not compute)
  Double_t       fClusMinE;         //minimum clus energy (<0 -> do not compute)
  Double_t       fTrackMinPt;       //minimum track pt    (<0 -> do not compute)
  UInt_t         fTriggers;         //if not zero only process given trigges
  Bool_t         fIsFastOnly;       //!=true if FASTONLY event is found
  Bool_t         fIsLedEvent;       //!=true if LED event is found
  Bool_t         fIsGoodEvent;      //!=true if good EMCAL event
  Double_t       fCellMaxE;         //!maximum cell energy in event
  Double_t       fClusMaxE;         //!maximum clus energy in event
  Double_t       fTrackMaxPt;       //!maximum track pt in event
  
  ClassDef(AliEmcalPhysicsSelection, 3); // Emcal physics selection
};
#endif
