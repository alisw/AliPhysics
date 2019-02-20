#ifndef CEPEVENTBUFFER
#define CEPEVENTBUFFER

#include "TObject.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "AliVEvent.h"
#include "AliVVZERO.h"
#include "AliESDAD.h"
#include "AliCEPBase.h"
#include "TParticle.h"
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
    Short_t fnITSCluster[6];   // Number of ITS clusters on all 6 layers
                               // and number of offline fired chips
    Short_t fFiredChips[4];    // Number of FastOR chips in the two SPD layers
                               // and number of offline fired chips
    TBits  fSPD_0STG_Online;   // using FastOrMap    (online)
    TBits  fSPD_0STG_Offline;  // using FiredChipMap (offline)


    Bool_t fisDGTrigger;
    UInt_t fisSTGTriggerFired;
    Int_t  fnTOFmaxipads;
    AliTOFTriggerMask *fTOFTriggerMask;
    
    // see AliCEPBase.h for the definition of fEventCondition
    UInt_t fEventCondition;

    // summary track information
    Int_t fnTracklets;      // total number of tracklets
    Int_t fnSingles;        // number of clusters on SPD layer 1 and 2
                            // not associated with a tracklet on other SPD 
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
    
    Int_t fnCaloCluster[2]; // number of EMCal and PHOS cluster
    Double_t fCaloEnergy[2];// total energy in EMCal/PHOS
    Double_t fdPhiEtaMinMax;// max distance of EMC cluster and charged track hit
    
    // Monte Carlo information
    TString  fMCGenerator;
    UInt_t   fMCProcessType; 
    TLorentzVector fMCIniSystem, fMCParticle;
    TVector3 fMCVtxPos;
    Int_t fnMCParticles[6];
    
    // list of tracks
    TObjArray *fCEPTracks;
    // tracklet - track associations
    TObjArray *fTrl2Tr;
    
    
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
    void SetnITSCluster(Short_t nclus[6]) {
      for (Int_t ii=0; ii<6; ii++) fnITSCluster[ii]=nclus[ii];
    }
    void SetnFiredChips(Short_t chips[4]) {
      for (Int_t ii=0; ii<4; ii++) fFiredChips[ii]=chips[ii];
    }
    // every 100ns, 1200 FastOR signals
    void SetSPMapOnline(TBits map)     { fSPD_0STG_Online=map; }
    void SetSPMapOffline(TBits map)    { fSPD_0STG_Offline=map; }

    void SetisSTGTriggerFired(UInt_t stgtrig) { fisSTGTriggerFired = stgtrig; }
    void SetnTOFmaxipads(Int_t nmaxipads) { fnTOFmaxipads = nmaxipads; }
    void SetTOFTriggerMask(AliTOFTriggerMask *triggermask)
      { fTOFTriggerMask = triggermask; }

    void SetEventCondition(UInt_t evcond) { fEventCondition = evcond; }
 
    // fnTracks, fnTracksCombined, and fnTracksITSpure are incremented
    // automatically when tracks are added with the method AddTrack
    void AddTrack(CEPTrackBuffer* trk);
    
    // the number of tracklets and residuals, as well as the enumber
    // of tracks passing Martin's selection have to be set separately
    void SetnTracksTotal(Int_t ntrks)   { fnTracksTotal = ntrks; }
    
    void SetnTracklets(Int_t ntrklts)   { fnTracklets = ntrklts; }
    void SetnSingles(Int_t nsingle)     { fnSingles = nsingle; }
    void AddTrl2Tr(TVector2 *vec, Int_t pos)        {
      fTrl2Tr->AddAt((TObject*)vec,pos);
    }
      
    void SetnResiduals(Int_t nres)      { fnResiduals = nres; }
    void SetnMSelection(Int_t nMsel)    { fnMSelection = nMsel; }
    
    void SetnV0(Int_t nV0)              { fnV0 = nV0; }
    void SetVtxType(UInt_t vtxtype)     { fVtxType = vtxtype; }
    void SetVtxPos(Double_t xp,Double_t yp,Double_t zp)
                                        { fVtxPos.SetXYZ(xp,yp,zp); }
    void SetnCaloCluster(Int_t nclus, Int_t ind)
                                        { if (ind==0 || ind==1) fnCaloCluster[ind] = nclus; }
    void SetCaloEnergy(Double_t ene, Int_t ind)
                                        { if (ind==0 || ind==1) fCaloEnergy[ind] = ene; }
    void SetdPhiEtaMinMax(Double_t dminmax)
                                        { fdPhiEtaMinMax = dminmax; }
    
    void SetVtxPos(TVector3 vtxpos)     { fVtxPos = TVector3(vtxpos); }
    
    void SetMCGenerator(TString MCGenerator) { fMCGenerator = MCGenerator; }
    void SetMCProcessType(UInt_t MCProcess)  { fMCProcessType = MCProcess; }
    void SetMCVtxPos(Double_t xp,Double_t yp,Double_t zp)
                                             { fMCVtxPos.SetXYZ(xp,yp,zp); }
    void SetMCIniSystem(TLorentzVector part) { fMCIniSystem = part; }
    void SetMCParticle(TLorentzVector part)  { fMCParticle = part; }
    void SetMCnParticles(Int_t *nparts)      {
      for (Int_t ii=0; ii<6; ii++) fnMCParticles[ii]=nparts[ii];
    }

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
    Short_t GetnITSCluster(Int_t layer) {
      return (layer>=0 && layer<=5) ? fnITSCluster[layer] : -1; }
    Short_t GetnFiredChips(Int_t layer) {
      return (layer>=0 && layer<=3) ? fFiredChips[layer] : -1; }
    TBits GetSPMapOnline()  { return fSPD_0STG_Online;  }
    TBits GetSPMapOffline() { return fSPD_0STG_Offline; }

    UInt_t  GetisSTGTriggerFired() { return fisSTGTriggerFired; }
    Int_t   GetnTOFmaxipads() { return fnTOFmaxipads; }
    AliTOFTriggerMask *GetTOFTriggerMask() { return fTOFTriggerMask; }

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
    Int_t GetnSingles()        const { return fnSingles; }
    TObjArray* GetTrl2Tr()           { return fTrl2Tr; }
    Int_t GetnResiduals()      const { return fnResiduals; }
    Int_t GetnMSelection()     const { return fnMSelection; }

    Int_t GetnV0()       const { return fnV0; }
    UInt_t GetVtxType()  const { return fVtxType; }
    TVector3 GetVtxPos() const { return fVtxPos; }
    Double_t GetdPhiEtaMinMax() const { return fdPhiEtaMinMax; }
    
    Int_t GetnCaloCluster(Int_t ind) const { return (ind==0 || ind==1) ? fnCaloCluster[ind] : 0; }
    Double_t GetCaloEnergy(Int_t ind) const { return (ind==0 || ind==1) ? fCaloEnergy[ind] : 0; }

    TString GetMCGenerator()        const { return fMCGenerator; }
    UInt_t GetMCProcessType()       const { return fMCProcessType; }
    TVector3 GetMCVtxPos()          const { return fMCVtxPos; }
    TLorentzVector GetMCIniSystem() const { return fMCIniSystem; }
    TLorentzVector GetMCParticle()  const { return fMCParticle; }
    Int_t* GetnMCParticles()              { return fnMCParticles; }
     
    // readout gap condition
    UInt_t  GetEventCondition() const { return fEventCondition; }
    
    Bool_t isEventCut()   const { return fEventCondition & AliCEPBase::kETEventCut; }
    Bool_t isPhyssel()    const { return fEventCondition & AliCEPBase::kETPhyssel;   }
    Bool_t isPileup()     const { return fEventCondition & AliCEPBase::kETPileup;    }
    Bool_t isClusterCut() const { return fEventCondition & AliCEPBase::kETClusterCut;}
    Bool_t isDGTrigger()  const { return fEventCondition & AliCEPBase::kETDGTrigger; }

    Bool_t isMBOR() const { return fEventCondition & AliCEPBase::kETMBOR; }
    Bool_t isMBAND()const { return fEventCondition & AliCEPBase::kETMBAND; }
    Bool_t isTPC()  const { return isTPCA() || isTPCC(); }
    Bool_t isTPCA() const { return fEventCondition & AliCEPBase::kETTPCA; }
    Bool_t isTPCC() const { return fEventCondition & AliCEPBase::kETTPCC; }
    Bool_t isSPD()  const { return isSPDA() || isSPDC(); }
    Bool_t isSPDA() const { return fEventCondition & AliCEPBase::kETSPDA; }
    Bool_t isSPDC() const { return fEventCondition & AliCEPBase::kETSPDC; }
    Bool_t isFMD()  const { return isFMDA() || isFMDC(); }
    Bool_t isFMDA() const { return fEventCondition & AliCEPBase::kETFMDA; }
    Bool_t isFMDC() const { return fEventCondition & AliCEPBase::kETFMDC; }
    Bool_t isV0()   const { return isV0A() || isV0C(); }
    Bool_t isV0A()  const { return fEventCondition & AliCEPBase::kETV0A;  }
    Bool_t isV0C()  const { return fEventCondition & AliCEPBase::kETV0C;  }
    Bool_t isAD()   const { return isADA() || isADC();  }
    Bool_t isADA()  const { return fEventCondition & AliCEPBase::kETADA;  }
    Bool_t isADC()  const { return fEventCondition & AliCEPBase::kETADC;  }
    Bool_t isZDC()  const { return isZDCA() || isZDCC(); }
    Bool_t isZDCA() const { return fEventCondition & AliCEPBase::kETZDCA; }
    Bool_t isZDCC() const { return fEventCondition & AliCEPBase::kETZDCC; }
    Bool_t isZDN()  const { return isZDNA() || isZDNC(); }
    Bool_t isZDNA() const { return fEventCondition & AliCEPBase::kETZDNA; }
    Bool_t isZDNC() const { return fEventCondition & AliCEPBase::kETZDNC; }

    CEPTrackBuffer* GetTrack(Int_t ind);
    Bool_t RemoveTrack(Int_t ind);

    ClassDef(CEPEventBuffer, 5)     // CEP event buffer

};

#endif
