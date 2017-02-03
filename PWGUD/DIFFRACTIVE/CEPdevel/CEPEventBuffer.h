#ifndef CEPEVENTBUFFER
#define CEPEVENTBUFFER

#include "TObject.h"
#include "TList.h"
#include "TArrayI.h"
#include "AliCEPBase.h"
#include "TClonesArray.h"
#include "CEPTrackBuffer.h"

class CEPEventBuffer : public TObject {

  private:
    // event header information
    Int_t fRunNumber;
    Int_t fEventNumber;
    Int_t fPeriodNumber;
    Int_t fOrbitNumber;
    Int_t fBunchCrossNumber;
    Int_t fCollissionType;
    Double_t fMagnField;
    
    // general event features
    Bool_t fPhysel;
    Bool_t fEvHandlersel;
    Bool_t fisPileup;
    Bool_t fisClusterCut;
    Bool_t fisMBOR;
    Bool_t fisMBAND;

    // see AliCEPBase.h for the definition of fGapCondition
    Int_t fGapCondition;

    // summary track information
    Int_t fnTracklets;      // total number of tracklets
    Int_t fnTracks;         // total number of tracks
    Int_t fnTracksCombined; // ITS+TPC
    Int_t fnTracksITSpure;  // ITS only
    Int_t fnResiduals;      // tracklets without associated track
    Int_t fnMSelection;     // tarcks after Martin's selection
                            // negative value: full event is rejected
    Int_t fnV0;             // number of V0s

    Int_t fVtxType;         // see AliCEPBase.h for a definition of the bits
    TVector3 fVtxPos;
    
    // Monte Carlo information
    TString  fMCGenerator;
    Int_t    fMCProcessType; 
    TVector3 fMCVtxPos;
    
    // list of tracks
    TList *fCEPTracks;
    
  public:
    CEPEventBuffer();
    ~CEPEventBuffer();
    
    // Modifiers
    void Reset();
    
    void SetRunNumber(Int_t runnum)     { fRunNumber = runnum; }
    void SetEventNumber(Int_t evnum)    { fEventNumber = evnum; }
    void SetPeriodNumber(Int_t pnum)    { fPeriodNumber = pnum; }
    void SetOrbitNumber(Int_t onum)     { fOrbitNumber = onum; }
    void SetBunchCrossNumber(Int_t bcnum) { fBunchCrossNumber = bcnum; }
    void SetCollissionType(Int_t coltype) { fCollissionType = coltype; }
    void SetMagnField(Double_t magnf)   { fMagnField = magnf; }

    void SetPhyssel(Bool_t physsel)     { fPhysel = physsel; }
    void SetEvHandlersel(Bool_t evhsel) { fEvHandlersel = evhsel; }
    void SetPileup(Bool_t isPileup)     { fisPileup = isPileup; }
    void SetClusterCut(Bool_t isClusterCut) { fisClusterCut = isClusterCut; }

    void SetGapCondition(Int_t gapcond) { fGapCondition = gapcond; }
 
    // fnTracks, fnTracksCombined, and fnTracksITSpure are incremented
    // automatically when tracks are added with the method AddTrack
    void AddTrack(CEPTrackBuffer* trk);
    
    // the number of tracklets and residuals has to be set separately
    void SetnTracklets(Int_t ntrklts)   { fnTracklets = ntrklts; }
    void SetnResiduals(Int_t nres)      { fnResiduals = nres; }
    
    void SetnV0(Int_t nV0)              { fnV0 = nV0; }
    void SetVtxType(Int_t vtxtype)      { fVtxType = vtxtype; }
    void SetVtxPos(Double_t xp,Double_t yp,Double_t zp)
                                        { fVtxPos.SetXYZ(xp,yp,zp); }
      
    void SetMCGenerator(TString MCGenerator) { fMCGenerator = MCGenerator; }
    void SetMCProcessType(Int_t MCProcess)   { fMCProcessType = MCProcess; }
    void SetMCVtxPos(Double_t xp,Double_t yp,Double_t zp)
                                             { fMCVtxPos.SetXYZ(xp,yp,zp); }

    // Accessors
    Int_t GetRunNumber()     const { return fRunNumber; }
    Int_t GetEventNumber()   const { return fEventNumber; }
    Int_t GetPeriodNumber()  const { return fPeriodNumber; }
    Int_t GetOrbitNumber()   const { return fOrbitNumber; }
    Int_t GetBunchCrossNumber() const { return fBunchCrossNumber; }
    Int_t GetCollissionType()const { return fCollissionType; }
    Double_t GetMagnField()  const { return fMagnField; }

    Bool_t GetPhyssel()      const { return fPhysel; }
    Bool_t GetEvHandlersel() const { return fEvHandlersel; }

    // different ways of retrieving number of tracks
    Int_t GetnTracks()       const { return fnTracks; }
    Int_t GetnTracks(Int_t mask, Int_t pattern);
    Int_t GetnTracks(Int_t mask, Int_t pattern, TArrayI *indices);
    Int_t GetnTracks(TArrayI *masks, TArrayI *patterns);
    Int_t GetnTracks(TArrayI *masks, TArrayI *patterns, TArrayI *indices);

    Int_t GetnTracklets()    const { return fnTracklets; }
    Int_t GetnTracksCombined() const { return fnTracksCombined; }
    Int_t GetnITSpureTracks()const { return fnTracksITSpure; }
    Int_t GetnResiduals()    const { return fnResiduals; }
    Int_t GetnMSelection()   const { return fnMSelection; }

    Int_t GetnTracksisSet(Int_t TTest);
    Int_t GetnTracksisEqual(Int_t TTest);


    Int_t GetVtxType()   const { return fVtxType; }
    TVector3 GetVtxPos() const { return fVtxPos; }

    TString GetMCGenerator() const { return fMCGenerator; }
    Int_t GetMCProcessType() const { return fMCProcessType; }

    TVector3 GetMCVtxPos() const { return fMCVtxPos; }
  
    // readout gap condition
    Int_t  GetGapCondition() const { return fGapCondition; }
    Bool_t isTPCA() const { return fGapCondition & AliCEPBase::kBitTPCA; }
    Bool_t isTPCC() const { return fGapCondition & AliCEPBase::kBitTPCC; }
    Bool_t isSPD()  const { return isSPDA() || isSPDC(); }
    Bool_t isSPDA() const { return fGapCondition & AliCEPBase::kBitSPDA; }
    Bool_t isSPDC() const { return fGapCondition & AliCEPBase::kBitSPDC; }
    Bool_t isFMD()  const { return isFMDA() || isFMDC(); }
    Bool_t isFMDA() const { return fGapCondition & AliCEPBase::kBitFMDA; }
    Bool_t isFMDC() const { return fGapCondition & AliCEPBase::kBitFMDC; }
    Bool_t isV0()   const { return isV0A() || isV0C(); }
    Bool_t isV0A()  const { return fGapCondition & AliCEPBase::kBitV0A;  }
    Bool_t isV0C()  const { return fGapCondition & AliCEPBase::kBitV0C;  }
    Bool_t isAD()   const { return isADA() || isADC();  }
    Bool_t isADA()  const { return fGapCondition & AliCEPBase::kBitADA;  }
    Bool_t isADC()  const { return fGapCondition & AliCEPBase::kBitADC;  }
    Bool_t isZDC()  const { return isZDCA() || isZDCC(); }
    Bool_t isZDCA() const { return fGapCondition & AliCEPBase::kBitZDCA; }
    Bool_t isZDCC() const { return fGapCondition & AliCEPBase::kBitZDCC; }
    Bool_t isZDN()  const { return isZDNA() || isZDNC(); }
    Bool_t isZDNA() const { return fGapCondition & AliCEPBase::kBitZDNA; }
    Bool_t isZDNC() const { return fGapCondition & AliCEPBase::kBitZDNC; }
    Bool_t isPileup()     const { return fisPileup; }
    Bool_t isClusterCut() const { return fisClusterCut; }
    Bool_t isMBOR()       const { return fisMBOR; }
    Bool_t isMBAND()      const { return fisMBAND; }
    
    CEPTrackBuffer* GetTrack(Int_t ind);

  ClassDef(CEPEventBuffer, 3)     // CEP event buffer

};

#endif
