/* $Id: AliOfflineTrigger.h 35782 2009-10-22 11:54:31Z jgrosseo $ */

#ifndef ALIOFFLINETRIGGER_H
#define ALIOFFLINETRIGGER_H

#include <TObject.h>
#include <AliPWG0Helper.h>

class AliESDEvent;

class AliOfflineTrigger : public TObject
{
  public:
    enum AliceSide { kASide = 1, kCSide, kCentralBarrel };
    
    AliOfflineTrigger();
    virtual ~AliOfflineTrigger() {}

    Bool_t IsEventTriggered(const AliESDEvent* aEsd, AliPWG0Helper::Trigger trigger) const;
    
    void SetSPDGFOThreshhold(Int_t t) { fSPDGFOThreshold = t; }
    void SetV0Threshhold(Int_t aSide, Int_t cSide) { fV0AThreshold = aSide; fV0CThreshold = cSide; }
    
    Int_t GetSPDGFOThreshhold() const { return fSPDGFOThreshold; }
    Int_t GetV0AThreshold() const { return fV0AThreshold; }
    Int_t GetV0CThreshold() const { return fV0CThreshold; }

  protected:
    Bool_t SPDGFOTrigger(const AliESDEvent* aEsd) const;
    Bool_t V0Trigger(const AliESDEvent* aEsd, AliceSide side) const;
    Bool_t ZDCTrigger(const AliESDEvent* aEsd, AliceSide side) const;
    Bool_t FMDTrigger(const AliESDEvent* aEsd) const;

    Int_t fSPDGFOThreshold;         // number of chips to accept a SPD GF0 trigger
    Int_t fV0AThreshold;            // threshold for number of BB triggers in V0A
    Int_t fV0CThreshold;            // threshold for number of BB triggers in V0C
 
    ClassDef(AliOfflineTrigger, 0)
    
  private:
    AliOfflineTrigger(const AliOfflineTrigger&);
    AliOfflineTrigger& operator=(const AliOfflineTrigger&);
};

#endif
