#ifndef CEPEVENTBUFFER
#define CEPEVENTBUFFER

#include "TObject.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "AliVEvent.h"
#include "AliVVZERO.h"
#include "AliESDAD.h"
#include "AliCEPBase.h"
#include "CEPTrackBuffer.h"

class CEPEventBuffer : public TObject {

  private:
    Int_t fRunNumber;
    Int_t fEventNumber;
    Int_t fPeriodNumber;
    Int_t fOrbitNumber;
    Int_t fBunchCrossNumber;
    Int_t fCollissionType;
    Double_t fMagnField;
    TString fFiredTriggerClasses;
    
    Bool_t fPFBBFlagV0[21];
    Bool_t fPFBGFlagV0[21];
    Bool_t fPFBBFlagAD[21];
    Bool_t fPFBGFlagAD[21];

    // general event features
    Bool_t fEventCutsel;
    Bool_t fPhysel;
    Bool_t fisPileup;
    Bool_t fisClusterCut;
    Bool_t fisDGTrigger;

    // see AliCEPBase.h for the definition of fEventCondition
    UInt_t fEventCondition;

    // summary track information
    Int_t fnTracklets;      // total number of tracklets
    Int_t fnTracksTotal;    // total number of tracks
    Int_t fnTracks;         // number of tracks in fCEPTracks
    Int_t fnTracksCombined; // ITS+TPC
    Int_t fnTracksITSpure;  // ITS only
    Int_t fnResiduals;      // tracklets without associated track
    Int_t fnMSelection;     // tracks after Martin's selection
                            // negative value: full event is rejected
    Int_t fnV0;             // number of V0s

    UInt_t   fVtxType;      // see AliCEPBase.h for a definition of the bits
    TVector3 fVtxPos;       // vertex position
    
    // Monte Carlo information
    TString  fMCGenerator;
    UInt_t   fMCProcessType; 
    TVector3 fMCVtxPos;
    
    // list of tracks
    TObjArray *fCEPTracks;
    
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
    void SetFiredTriggerClasses(TString ftc)   { fFiredTriggerClasses = ftc; }
    void SetPFFlags(AliVEvent *Event);

    void SetEventCondition(UInt_t evcond) { fEventCondition = evcond; }
 
    // fnTracks, fnTracksCombined, and fnTracksITSpure are incremented
    // automatically when tracks are added with the method AddTrack
    void AddTrack(CEPTrackBuffer* trk);
    
    // the number of tracklets and residuals, as well as the enumber
    // of tracks passing Martin's selection have to be set separately
    void SetnTracksTotal(Int_t ntrks) { fnTracksTotal = ntrks; }
    void SetnTracklets(Int_t ntrklts) { fnTracklets = ntrklts; }
    void SetnResiduals(Int_t nres)    { fnResiduals = nres; }
    void SetnMSelection(Int_t nMsel)  { fnMSelection = nMsel; }
    
    void SetnV0(Int_t nV0)              { fnV0 = nV0; }
    void SetVtxType(UInt_t vtxtype)     { fVtxType = vtxtype; }
    void SetVtxPos(Double_t xp,Double_t yp,Double_t zp)
                                        { fVtxPos.SetXYZ(xp,yp,zp); }
    void SetVtxPos(TVector3 vtxpos)     { fVtxPos = TVector3(vtxpos); }
    void SetMCGenerator(TString MCGenerator) { fMCGenerator = MCGenerator; }
    void SetMCProcessType(UInt_t MCProcess)  { fMCProcessType = MCProcess; }
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
    TString GetFiredTriggerClasses() const { return fFiredTriggerClasses; }
    Bool_t* GetPFBBFlagV0() { return fPFBBFlagV0; }
    Bool_t* GetPFBGFlagV0() { return fPFBGFlagV0; }
    Bool_t* GetPFBBFlagAD() { return fPFBBFlagAD; }
    Bool_t* GetPFBGFlagAD() { return fPFBGFlagAD; }

    // different ways of retrieving number of tracks
    Int_t GetnTracksTotal()  const { return fnTracksTotal; }
    Int_t GetnTracks()       const { return fnTracks; }
    Int_t GetnTracks(UInt_t mask, UInt_t pattern);
    Int_t GetnTracks(UInt_t mask, UInt_t pattern, TArrayI *indices);
    Int_t GetnTracks(TArrayI *masks, TArrayI *patterns);
    Int_t GetnTracks(TArrayI *masks, TArrayI *patterns, TArrayI *indices);

    Int_t GetnTracksCombined() const { return fnTracksCombined; }
    Int_t GetnITSpureTracks()  const { return fnTracksITSpure; }
    Int_t GetnTracklets()      const { return fnTracklets; }
    Int_t GetnResiduals()      const { return fnResiduals; }
    Int_t GetnMSelection()     const { return fnMSelection; }

    Int_t GetnV0()       const { return fnV0; }
    UInt_t GetVtxType()  const { return fVtxType; }
    TVector3 GetVtxPos() const { return fVtxPos; }

    TString GetMCGenerator() const { return fMCGenerator; }
    UInt_t GetMCProcessType()const { return fMCProcessType; }
    TVector3 GetMCVtxPos()   const { return fMCVtxPos; }
  
    // readout gap condition
    UInt_t  GetEventCondition() const { return fEventCondition; }
    
    Bool_t isEventCut()   const { return fEventCondition & AliCEPBase::kBitEventCut; }
    Bool_t isPhyssel()    const { return fEventCondition & AliCEPBase::kBitPhyssel;   }
    Bool_t isPileup()     const { return fEventCondition & AliCEPBase::kBitPileup;    }
    Bool_t isClusterCut() const { return fEventCondition & AliCEPBase::kBitClusterCut;}
    Bool_t isDGTrigger()  const { return fEventCondition & AliCEPBase::kBitDGTrigger; }

    Bool_t isMBOR() const { return fEventCondition & AliCEPBase::kBitMBOR; }
    Bool_t isMBAND()const { return fEventCondition & AliCEPBase::kBitMBAND; }
    Bool_t isTPC()  const { return isTPCA() || isTPCC(); }
    Bool_t isTPCA() const { return fEventCondition & AliCEPBase::kBitTPCA; }
    Bool_t isTPCC() const { return fEventCondition & AliCEPBase::kBitTPCC; }
    Bool_t isSPD()  const { return isSPDA() || isSPDC(); }
    Bool_t isSPDA() const { return fEventCondition & AliCEPBase::kBitSPDA; }
    Bool_t isSPDC() const { return fEventCondition & AliCEPBase::kBitSPDC; }
    Bool_t isFMD()  const { return isFMDA() || isFMDC(); }
    Bool_t isFMDA() const { return fEventCondition & AliCEPBase::kBitFMDA; }
    Bool_t isFMDC() const { return fEventCondition & AliCEPBase::kBitFMDC; }
    Bool_t isV0()   const { return isV0A() || isV0C(); }
    Bool_t isV0A()  const { return fEventCondition & AliCEPBase::kBitV0A;  }
    Bool_t isV0C()  const { return fEventCondition & AliCEPBase::kBitV0C;  }
    Bool_t isAD()   const { return isADA() || isADC();  }
    Bool_t isADA()  const { return fEventCondition & AliCEPBase::kBitADA;  }
    Bool_t isADC()  const { return fEventCondition & AliCEPBase::kBitADC;  }
    Bool_t isZDC()  const { return isZDCA() || isZDCC(); }
    Bool_t isZDCA() const { return fEventCondition & AliCEPBase::kBitZDCA; }
    Bool_t isZDCC() const { return fEventCondition & AliCEPBase::kBitZDCC; }
    Bool_t isZDN()  const { return isZDNA() || isZDNC(); }
    Bool_t isZDNA() const { return fEventCondition & AliCEPBase::kBitZDNA; }
    Bool_t isZDNC() const { return fEventCondition & AliCEPBase::kBitZDNC; }

    CEPTrackBuffer* GetTrack(Int_t ind);
    Bool_t RemoveTrack(Int_t ind);

  ClassDef(CEPEventBuffer, 3)     // CEP event buffer

};

#endif
