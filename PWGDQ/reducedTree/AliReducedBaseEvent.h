// Class for reduced event information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 

#ifndef ALIREDUCEDBASEEVENT_H
#define ALIREDUCEDBASEEVENT_H

#include <TObject.h>
#include <TClonesArray.h>

class AliReducedBaseTrack;
class AliReducedPairInfo;

//_________________________________________________________________________
class AliReducedBaseEvent : public TObject {

  friend class AliAnalysisTaskReducedTreeMaker;     // friend analysis task which fills the object
  friend class AliReducedAnalysisFilterTrees;
  
 public:
  enum ETrackOption {
     kNoInit=0,
     kUseBaseTracks,            // use AliReducedBaseTrack for the track array
     kUseReducedTracks            // use AliReducedTrackInfo for the track array
  };
  
  // bits toggled in the fEventTag data member
  enum EventTagBits {      
     kAnaUtils2013pPb=0,   // 0 - 2013 p-Pb event selection
     kAnaUtilPileupMV,     // 1 - multi-vertexer (MV) pileup 
     kAnaUtilPileupMV2,    // 2 - MV pileup without bunch-crossing check
     kAnaUtilPileupMV3,    // 3 - MV pileup with min weighted distance 10 (instead of 15)
     kAnaUtilPileupMV4,    // 4 - MV pileup with min weighted distance 5 (instead of 15)
     kIsPileupFromSPD1,    // 5 - event->IsPileupFromSPD(3,0.6,3.,2.,5.)
     kIsPileupFromSPD2,    // 6 - event->IsPileupFromSPD(4,0.6,3.,2.,5.)
     kIsPileupFromSPD3,    // 7 - event->IsPileupFromSPD(5,0.6,3.,2.,5.)
     kIsPileupFromSPD4,    // 8 - event->IsPileupFromSPD(6,0.6,3.,2.,5.)
     kIsPileupFromSPD5,    // 9 - event->IsPileupFromSPD(3,0.8,3.,2.,5.)
     kIsPileupFromSPD6,    // 10 - event->IsPileupFromSPD(4,0.8,3.,2.,5.)
     kIsPileupFromSPD7,    // 11 - event->IsPileupFromSPD(5,0.8,3.,2.,5.)
     kIsPileupFromSPD8,    // 12 - event->IsPileupFromSPD(6,0.8,3.,2.,5.)
     kVtxDistanceSelected, // 13 - Improved cut on the distance between SPD and track vertices 
     kUnbiasedEvent,       // 14 - event selected for writing in the trees on a random basis 
     kTimeRange,           // 15 - selected by AliTimeRangeCut (to be rejected)
     kNEventTagBits
  };
  
 public:
  AliReducedBaseEvent();
  AliReducedBaseEvent(const Char_t* name, Int_t trackOption=kNoInit, Int_t track2Option=kNoInit);
  virtual ~AliReducedBaseEvent();
  
  virtual void CopyEventHeader(const AliReducedBaseEvent* other);

  // getters
  ULong64_t EventTag()                        const {return fEventTag;}
  Bool_t    EventTag(UShort_t bit)            const {return (bit<8*sizeof(ULong64_t) ? (fEventTag&(ULong64_t(1)<<bit)) : kFALSE);}
  Int_t     RunNo()                           const {return fRunNo;}
  Float_t   Vertex(Int_t axis)                const {return (axis>=0 && axis<=2 ? fVtx[axis] : 0);}
  Int_t     VertexNContributors()             const {return fNVtxContributors;}
  Float_t   CentralityVZERO()                 const {return fCentrality[0];}
  Float_t   CentralitySPD()                   const {return fCentrality[1];}
  Float_t   CentralityTPC()                   const {return fCentrality[2];}
  Float_t   CentralityZEMvsZDC()              const {return fCentrality[3];}
  Float_t   CentralityVZEROA()                const {return fCentrality[4];}
  Float_t   CentralityVZEROC()                const {return fCentrality[5];}
  Float_t   CentralityZNA()                   const {return fCentrality[6];}
  Int_t     CentralityQuality()               const {return fCentQuality;}
  Int_t     NTracksTotal()                    const {return fNtracks[0];}
  Int_t     NTracks()                         const {return fNtracks[1];}
  Int_t     NTracks1()                       const {return (fTracks ? fTracks->GetEntries() : 0);}
  Int_t     NTracks2()                       const {return (fTracks2 ? fTracks2->GetEntries() : 0);}
  Int_t     NV0CandidatesTotal()              const {return fNV0candidates[0];}
  Int_t     NV0Candidates()                   const {return fNV0candidates[1];}
  Int_t     NPairs()                   const {return fCandidates->GetEntries();}
  
  AliReducedBaseTrack* GetTrack(Int_t i) const {return (fTracks && i>=0 && i<fTracks->GetEntries() ? (AliReducedBaseTrack*)fTracks->At(i) : 0x0);}
  AliReducedBaseTrack* GetTrack2(Int_t i) const {return (fTracks2 && i>=0 && i<fTracks2->GetEntries() ? (AliReducedBaseTrack*)fTracks2->At(i) : 0x0);}
  TClonesArray* GetTracks()          const {return fTracks;}
  TClonesArray* GetTracks2()        const {return fTracks2;}
  
  AliReducedPairInfo* GetV0Pair(Int_t i)         const 
  {return (i>=0 && i<fNV0candidates[1] ? (AliReducedPairInfo*)fCandidates->At(i) : 0x0);}
  AliReducedPairInfo* GetPair(Int_t i)         const 
  {return (i>=0 && i<fCandidates->GetEntries() ? (AliReducedPairInfo*)fCandidates->At(i) : 0x0);}
  TClonesArray* GetPairs()                       const {return fCandidates;}
  
  Bool_t    TestEventTag(UShort_t iflag) const {return (iflag<8*sizeof(ULong64_t) ? fEventTag&(ULong64_t(1)<<iflag) : kFALSE);}
  Bool_t    SetEventTag(UShort_t iflag)        {if (iflag>=8*sizeof(ULong64_t)) return kFALSE; fEventTag|=(ULong64_t(1)<<iflag); return kTRUE;}
  
  virtual void ClearEvent();
  
 protected:
  ULong64_t fEventTag;        // Event tags to be used either during analysis or to filter events
  Int_t     fRunNo;                 // run number
  Float_t   fVtx[3];                // global event vertex vector in cm
  Int_t     fNVtxContributors;      // global event vertex contributors
  Float_t   fCentrality[7];         // centrality; 0-V0M, 1-CL1, 2-TRK, 3-ZEMvsZDC, 4-V0A, 5-V0C, 6-ZNA
  Int_t     fCentQuality;           // quality flag for the centrality 
  Int_t     fNtracks[2];            // number of tracks, [0]-total, [1]-selected for the tree
  Int_t     fNV0candidates[2];      // number of V0 candidates, [0]-total, [1]-selected for the tree
    
  TClonesArray* fTracks;            //->   array containing particles
  static TClonesArray* fgTracks;    //       global tracks

  TClonesArray* fTracks2;               //->   array containing additional particles
  static TClonesArray* fgTracks2;    //       global tracks
  
  TClonesArray* fCandidates;        //->   array containing pair candidates
  static TClonesArray* fgCandidates;  // pair candidates
  
  AliReducedBaseEvent& operator= (const AliReducedBaseEvent &c);
  AliReducedBaseEvent(const AliReducedBaseEvent &c);

  ClassDef(AliReducedBaseEvent, 3);
};

#endif
