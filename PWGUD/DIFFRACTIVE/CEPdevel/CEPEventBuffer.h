#ifndef CEPEVENTBUFFER
#define CEPEVENTBUFFER

#include "TObject.h"
#include "TList.h"
#include "AliPBBase.h"
#include "CEPTrackBuffer.h"

class CEPEventBuffer : public TObject {

  private:
    Int_t fRunNumber;
    Int_t fEventNumber;
    Int_t fnumTracks;
    Int_t fnumSoftTracks;
    Int_t fnumResiduals;
    
    // see AliPBBase.h for the definition of fGapCondition
    Int_t fGapCondition;
    
    TList   *fTracks;
    
  public:
    CEPEventBuffer();
    ~CEPEventBuffer();
    
    // Modifiers
    void Reset();
    void SetRunNumber(Int_t runnum)     { fRunNumber = runnum; }
    void SetEventNumber(Int_t evnum)    { fEventNumber = evnum; }
    void SetnumResiduals(Int_t nres)    { fnumResiduals = nres; }
    void SetGapCondition(Int_t gapcond) { fGapCondition = gapcond; }
    void AddTrack(CEPTrackBuffer* trk);
    
    // Accessors
    Int_t GetRunNumber()     const { return fRunNumber; }
    Int_t GetEventNumber()   const { return fEventNumber; }
    Int_t GetnumTracks()     const { return fnumTracks; }
    Int_t GetnumSoftTracks() const { return fnumSoftTracks; }

    // readout gap condition
    Int_t  GetGapCondition() const { return fGapCondition; }
    Bool_t isTPCA() const { return fGapCondition & AliPBBase::kBitTPCA; }
    Bool_t isTPCC() const { return fGapCondition & AliPBBase::kBitTPCC; }
    Bool_t isSPDA() const { return fGapCondition & AliPBBase::kBitSPDA; }
    Bool_t isSPDC() const { return fGapCondition & AliPBBase::kBitSPDC; }
    Bool_t isFMDA() const { return fGapCondition & AliPBBase::kBitFMDA; }
    Bool_t isFMDC() const { return fGapCondition & AliPBBase::kBitFMDC; }
    Bool_t isV0A()  const { return fGapCondition & AliPBBase::kBitV0A; }
    Bool_t isV0C()  const { return fGapCondition & AliPBBase::kBitV0C; }
    Bool_t isADA()  const { return fGapCondition & AliPBBase::kBitADA; }
    Bool_t isADC()  const { return fGapCondition & AliPBBase::kBitADC; }
    Bool_t isZDCA() const { return fGapCondition & AliPBBase::kBitZDCA; }
    Bool_t isZDCC() const { return fGapCondition & AliPBBase::kBitZDCC; }
    
    CEPTrackBuffer* GetTrack(Int_t ind);

  ClassDef(CEPEventBuffer,1)     // CEP event buffer

};

#endif
