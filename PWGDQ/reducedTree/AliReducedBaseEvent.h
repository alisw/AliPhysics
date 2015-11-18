// Class for reduced event information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 

#ifndef ALIREDUCEDBASEEVENT_H
#define ALIREDUCEDBASEEVENT_H

#include <TObject.h>
#include <TClonesArray.h>

class AliReducedBaseTrack;

//_________________________________________________________________________
class AliReducedBaseEvent : public TObject {

  friend class AliAnalysisTaskReducedTreeMaker;     // friend analysis task which fills the object

 public:
  AliReducedBaseEvent();
  AliReducedBaseEvent(const Char_t* name);
  virtual ~AliReducedBaseEvent();

  // getters
  ULong64_t EventTag()                        const {return fEventTag;}
  Bool_t    EventTag(UShort_t bit)            const {return (bit<8*sizeof(ULong64_t) ? (fEventTag&(ULong64_t(1)<<bit)) : kFALSE);}
  Int_t     RunNo()                           const {return fRunNo;}
  Float_t   Vertex(Int_t axis)                const {return (axis>=0 && axis<=2 ? fVtx[axis] : 0);}
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
  
  AliReducedBaseTrack* GetTrack(Int_t i) const {return (i<fNtracks[1] ? (AliReducedBaseTrack*)fTracks->At(i) : 0x0);}
  TClonesArray* GetTracks()          const {return fTracks;}
  
  Bool_t    TestEventTag(UShort_t iflag) const {return (iflag<8*sizeof(ULong64_t) ? fEventTag&(ULong64_t(1)<<iflag) : kFALSE);}
  Bool_t    SetEventTag(UShort_t iflag)        {if (iflag>=8*sizeof(ULong64_t)) return kFALSE; fEventTag|=(ULong64_t(1)<<iflag); return kTRUE;}
  
 protected:
  ULong64_t fEventTag;              // Event tags to be used either during analysis or to filter events
  Int_t     fRunNo;                 // run number
  Float_t   fVtx[3];                // global event vertex vector in cm
  Float_t   fCentrality[7];         // centrality; 0-V0M, 1-CL1, 2-TRK, 3-ZEMvsZDC, 4-V0A, 5-V0C, 6-ZNA
  Int_t     fCentQuality;           // quality flag for the centrality 
  Int_t     fNtracks[2];            // number of tracks, [0]-total, [1]-selected for the tree
    
  TClonesArray* fTracks;            //->   array containing particles
  static TClonesArray* fgTracks;    //       global tracks
  
  void ClearEvent();
  AliReducedBaseEvent(const AliReducedBaseEvent &c);
  AliReducedBaseEvent& operator= (const AliReducedBaseEvent &c);

  ClassDef(AliReducedBaseEvent, 1);
};

#endif
