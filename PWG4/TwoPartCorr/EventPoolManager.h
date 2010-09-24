// $Id$

#ifndef EventPoolManager_h
#define EventPoolManager_h

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

class GenericEventPool : public TObject
{
public:
  GenericEventPool(Int_t d) : 
    fMixDepth(d), fMultMin(-999), fMultMax(+999), 
    fZvtxMin(-999), fZvtxMax(+999), fWasUpdated(0), fMultBinIndex(0), 
    fZvtxBinIndex(0), fDebug(0) {;}
  GenericEventPool(Int_t d, Int_t multMin, Int_t multMax, 
                   Double_t zvtxMin, Double_t zvtxMax) : 
    fMixDepth(d), fMultMin(multMin), fMultMax(multMax), 
    fZvtxMin(zvtxMin), fZvtxMax(zvtxMax), fWasUpdated(0), fMultBinIndex(0), 
    fZvtxBinIndex(0), fDebug(0) {;}
  ~GenericEventPool() {;}

  Bool_t      EventMatchesBin(Int_t mult, Double_t zvtx) const;
  Bool_t      IsReady()                    const { return (Int_t)fEvents.size()==fMixDepth; }
  Int_t       Depth()                      const { return fEvents.size();                   }
  TObject    *GetRandomTrack()             const;
  TObjArray  *GetRandomEvent()             const;
  Int_t       MultBinIndex()               const { return fMultBinIndex;                    }
  Int_t       NTracksInEvent(Int_t iEvent) const;
  Int_t       NTracksInCurrentEvent()      const { return fNTracksInEvent.back();           }
  void        PrintInfo()                  const;
  Int_t       TracksInPool()               const;
  Bool_t      WasUpdated()                 const { return fWasUpdated;   }
  Int_t       ZvtxBinIndex()               const { return fZvtxBinIndex; }
  void        SetDebug(Bool_t b)                 { fDebug = b;           }
  Int_t       SetEventMultRange(Int_t multMin, Int_t multMax);
  Int_t       SetEventZvtxRange(Double_t zvtxMin, Double_t zvtxMax);
  void        SetMultBinIndex(Int_t iM) { fMultBinIndex = iM; }
  void        SetZvtxBinIndex(Int_t iZ) { fZvtxBinIndex = iZ; }
  Int_t       UpdatePool(int iEvent, const MyHeader *ev, TObjArray *trk);
  
protected:
  deque<TObjArray*>     fEvents;              //Holds TObjArrays of MyTracklets
  deque<int>            fNTracksInEvent;      //Tracks in event
  deque<int>            fEventIndex;          //Original event index
  Int_t                 fMixDepth;            //Number of evts. to mix with
  Int_t                 fMultMin, fMultMax;   //Track multiplicity bin range
  Double_t              fZvtxMin, fZvtxMax;   //Event z-vertex bin range
  Bool_t                fWasUpdated;          //Evt. succesfully passed selection?
  Int_t                 fMultBinIndex;        //Multiplicity bin
  Int_t                 fZvtxBinIndex;        //Zvertex bin
  Int_t                 fDebug;               //If 1 then debug on

  ClassDef(GenericEventPool,1) // Event pool class
};

class EventPoolManager : public TObject
{
public:
  EventPoolManager() : fNMultBins(0), fNZvtxBins(0) {;}
  ~EventPoolManager() {;}
  
  GenericEventPool       *GetEventPool(Int_t iMult, Int_t iZvtx) const;
  Int_t                   InitEventPools(Int_t depth, 
                                         Int_t nmultbins, Int_t *multbins, 
                                         Int_t nzvtxbins, Double_t *zvtxbins);
  Int_t                   UpdatePools(Int_t iEvent, const MyHeader *ev, TObjArray *trk);

protected:
  Int_t      fNMultBins;                                // number mult bins
  Int_t      fNZvtxBins;                                // number vertex bins
  std::vector<std::vector<GenericEventPool*> > fEvPool; // pool in bins of [fMultBin][fZvtxBin]

  ClassDef(EventPoolManager,1)
};
#endif
