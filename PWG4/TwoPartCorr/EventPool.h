// $Id$

#ifndef EventPool_h
#define EventPool_h

#include <vector>
#include <deque>
#include <Rtypes.h>
#include <Riostream.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>
#include "TreeClasses.h"

class EventPool : public TObject
{
public:
  EventPool(Int_t d) : fTracks(0), fMixDepth(d), fMultMin(-999), fMultMax(+999), 
                       fZvtxMin(-999), fZvtxMax(+999), fMult(0), fZvtx(0), 
                       fWasUpdated(0), fMultBinIndex(0), fZvtxBinIndex(0), fDebug(0) {;}
  EventPool(Int_t d, Int_t multMin, Int_t multMax, 
	    Double_t zvtxMin, Double_t zvtxMax) : fTracks(0), fMixDepth(d), fMultMin(multMin), fMultMax(multMax), 
                                                  fZvtxMin(zvtxMin), fZvtxMax(zvtxMax), fMult(0), fZvtx(0), 
                       fWasUpdated(0), fMultBinIndex(0), fZvtxBinIndex(0), fDebug(0) {;}
  ~EventPool() {;}

  Bool_t   EventMatchesBin(Int_t mult, Short_t zvtx) const;
  Bool_t   IsPoolReady()                const { return (Int_t)fEvents.size()==fMixDepth; }
  Int_t    Depth()                      const { return fEvents.size();                   }
  MyPart  *GetRandomTrack()             const;
  Int_t    MultBinIndex()               const { return fMultBinIndex;                    }
  Int_t    NTracksInEvent(Int_t iEvent) const;
  Int_t    NTracksInCurrentEvent()      const { return fNTracksInEvent.back();           }
  void     PrintInfo()                  const;
  Int_t    TracksInPool()               const;
  Bool_t   WasUpdated()                 const { return fWasUpdated;   }
  Int_t    ZvtxBinIndex()               const { return fZvtxBinIndex; }
  void     SetDebug(Bool_t b)                 { fDebug = b;           }
  Int_t    SetEventMultRange(Int_t multMin, Int_t multMax);
  Int_t    SetEventZvtxRange(Int_t zvtxMin, Int_t zvtxMax);
  void     SetMultBinIndex(Int_t iM) { fMultBinIndex = iM; }
  void     SetZvtxBinIndex(Int_t iZ) { fZvtxBinIndex = iZ; }
  Int_t    UpdatePool(int iEvent, const MyHeader *ev, TClonesArray* trk);
  
protected:
  TClonesArray         *fTracks;              // Copy of trk array. Refreshes each event.
  deque<TClonesArray*>  fEvents;              // The guts of the class. Holds TObjArrays of MyParts.
  deque<int>            fNTracksInEvent;      // Tracks in event
  deque<int>            fEventIndex;          // Original event index
  Int_t                 fMixDepth;            // Number of evts. to mix with
  Int_t                 fMultMin, fMultMax;   // Track multiplicity bin range
  Double_t              fZvtxMin, fZvtxMax;   // Event z-vertex bin range
  Int_t                 fMult;                // Tracks in current event
  Short_t               fZvtx;                // Current z-vertex
  Bool_t                fWasUpdated;          // Evt. succesfully passed selection?
  Int_t                 fMultBinIndex;        // Multiplicity bin
  Int_t                 fZvtxBinIndex;        // Zvertex bin
  Int_t                 fDebug;               // if 1 then debug on

  ClassDef(EventPool,2) // Event Pool class
};
#endif
