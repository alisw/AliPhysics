// $Id$

#ifndef AutoCorr_h
#define AutoCorr_h

#include <Rtypes.h>
#include <Riostream.h>
//#include <TTree.h>
//#include <TChain.h>
//#include <TCanvas.h>
//#include <TClonesArray.h>
//#include <TH1F.h>
//#include <TH2F.h>
#include <TMath.h>
//#include "TFile.h"
//#include "TSystem.h"
//#include "TROOT.h"
#include <vector>
//#include <deque>
#include <TRandom.h>
#include "TreeClasses.h"
#include "EventPool.h"

class AutoCorr : public TObject
{
public:
  AutoCorr(){;}
  ~AutoCorr(){;}
  
  
  Double_t   DeltaPhi(const MyPart &t1, const MyPart &t2, 
                      Double_t rangeMin = -TMath::Pi()/2, 
                      Double_t rangeMax = 3*TMath::Pi()/2)      const;
  Double_t   DeltaEta(const MyPart &t1, const MyPart &t2)       const;  
  Bool_t     IsMixedPairOk(const MyPart &t1, const MyPart &t2)  const;
  Bool_t     IsPairOk(const MyPart &t1, const MyPart &t2)       const;
  Bool_t     IsTrackOk(const MyPart &t, 
                       Double_t etaMin, Double_t etaMax)        const;
  Bool_t     IsTrackOk(const MyPart &t, 
                     Double_t etaMin, Double_t etaMax,
                       Double_t ptMin, Double_t ptMax)          const;
  Double_t   InvMass(const MyPart &p1, const MyPart &p2)        const;
  Bool_t     IsEventOk(const MyHeader &ev, Int_t minVc, 
                       Int_t maxNTracklets, Double_t zMin,
                       Double_t zMax)                           const;
  Bool_t     InBounds(Double_t val, Double_t min, Double_t max) const;
  Bool_t     InBounds(Int_t val, Int_t min, Int_t max)          const;

  EventPool* GetEventPool(Int_t iMult, Int_t iZvtx)             const;
  Int_t      InitEventPools(Int_t depth, 
                            Int_t nmultbins, Double_t multbins[], 
                            Int_t nzvtxbins, Double_t zvtxbins[]);
  Int_t      UpdatePools(int iEvent, MyHeader* ev, TClonesArray* trk);

protected:
  Int_t      fNMultBins;                         // mult bins
  Int_t      fNZvtxBins;                         // vertex bins
  std::vector<std::vector<EventPool*> > fEvPool; // pool in bins of [fMultBin][fZvtxBin]

  ClassDef(AutoCorr,1)
};
#endif
