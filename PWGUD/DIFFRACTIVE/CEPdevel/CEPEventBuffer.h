#ifndef CEPEVENTBUFFER
#define CEPEVENTBUFFER

#include "TObject.h"
#include "TList.h"
#include "AliPBBase.h"
#include "TClonesArray.h"
#include "CEPTrackBuffer.h"

class CEPEventBuffer : public TObject {

  private:
    Int_t fRunNumber;
    Int_t fEventNumber;
    Bool_t fPhysel;
    Int_t fnumTracks;
    Int_t fnumSoftTracks;
    Int_t fnumResiduals;
    Bool_t fisPileup;
    TVector3 fVertexPos;
    
    // Monte Carlo information
    TString fMCGenerator;
    Int_t fMCProcessType; 
    TVector3 fMCVertexPos;
    
    // see AliPBBase.h for the definition of fGapCondition
    Int_t fGapCondition;
    
    // list of tracks as TClonesArray
    TClonesArray *fCEPTracks;
    TClonesArray &ftrb;
    
  public:
    CEPEventBuffer();
    ~CEPEventBuffer();
    
    // Modifiers
    void Reset();
    void SetRunNumber(Int_t runnum)     { fRunNumber = runnum; }
    void SetEventNumber(Int_t evnum)    { fEventNumber = evnum; }
    void SetPhyssel(Bool_t physsel)     { fPhysel = physsel; }
    void SetnumResiduals(Int_t nres)    { fnumResiduals = nres; }
    void SetGapCondition(Int_t gapcond) { fGapCondition = gapcond; }
    void SetPileup(Bool_t isPileup)     { fisPileup = isPileup; }
    void SetVertexPos(Double_t xp,Double_t yp,Double_t zp)
      { fVertexPos.SetXYZ(xp,yp,zp); }
    void AddTrack(CEPTrackBuffer* trk);

    void SetMCGenerator(TString MCGenerator) { fMCGenerator = MCGenerator; }
    void SetMCProcessType(Int_t MCProcess) { fMCProcessType = MCProcess; }
    void SetMCVertexPos(Double_t xp,Double_t yp,Double_t zp)
      { fMCVertexPos.SetXYZ(xp,yp,zp); }

    // Accessors
    Int_t GetRunNumber()     const { return fRunNumber; }
    Int_t GetEventNumber()   const { return fEventNumber; }
    Bool_t GetPhyssel()      const { return fPhysel; }
    Int_t GetnumTracks()     const { return fnumTracks; }
    Int_t GetnumSoftTracks() const { return fnumSoftTracks; }
    TVector3 GetVertexPos()  const { return fVertexPos; }

    TString GetMCGenerator() const { return fMCGenerator; }
    Int_t GetMCProcessType() const { return fMCProcessType; }
    TVector3 GetMCVertexPos()const { return fMCVertexPos; }
  
    // readout gap condition
    Int_t  GetGapCondition() const { return fGapCondition; }
    Bool_t isTPCA() const { return fGapCondition & AliPBBase::kBitTPCA; }
    Bool_t isTPCC() const { return fGapCondition & AliPBBase::kBitTPCC; }
    Bool_t isSPDA() const { return fGapCondition & AliPBBase::kBitSPDA; }
    Bool_t isSPDC() const { return fGapCondition & AliPBBase::kBitSPDC; }
    Bool_t isFMDA() const { return fGapCondition & AliPBBase::kBitFMDA; }
    Bool_t isFMDC() const { return fGapCondition & AliPBBase::kBitFMDC; }
    Bool_t isV0A()  const { return fGapCondition & AliPBBase::kBitV0A;  }
    Bool_t isV0C()  const { return fGapCondition & AliPBBase::kBitV0C;  }
    Bool_t isADA()  const { return fGapCondition & AliPBBase::kBitADA;  }
    Bool_t isADC()  const { return fGapCondition & AliPBBase::kBitADC;  }
    Bool_t isZDCA() const { return fGapCondition & AliPBBase::kBitZDCA; }
    Bool_t isZDCC() const { return fGapCondition & AliPBBase::kBitZDCC; }
    Bool_t isPileup() const { return fisPileup; }
    
    CEPTrackBuffer* GetTrack(Int_t ind);

  ClassDef(CEPEventBuffer,2)     // CEP event buffer

};

#endif
