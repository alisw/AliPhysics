// $Id:$

#ifndef AutoCorr_h
#define AutoCorr_h

#include <Riostream.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include "Rtypes.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include <vector>
#include <deque>
#include "TRandom.h"
#include "TreeClasses.h"
#include "EventPool.h"

class AutoCorr : public TObject
{
public:
  AutoCorr(){;}
  ~AutoCorr(){;}
  
  // RangeMin,Max specify periodic boundaries
  Double_t DeltaPhi(MyPart* t1, MyPart* t2, 
		    double rangeMin = -TMath::Pi()/2, 
		    double rangeMax = 3*TMath::Pi()/2);
  Double_t DeltaEta(MyPart* t1, MyPart* t2);  
  Bool_t IsTrackOk(MyPart* t);
  Bool_t IsPairOk(MyPart* t1, MyPart* t2);
  Bool_t IsMixedPairOk(MyPart* t1, MyPart* t2);
//  Bool_t IsEventOk(MyHeader *h);

  Int_t InitEventPools(Int_t depth, 
		       Int_t nmultbins, Double_t multbins[], 
		       Int_t nzvtxbins, Double_t zvtxbins[]);

  Int_t UpdatePools(int iEvent, MyHeader* ev, TClonesArray* trk);
  EventPool* GetEventPool(int iMult, int iZvtx);// { return fEvPool.at(iMult).at(iZvtx); }
protected:
  Int_t fNMultBins;
  Int_t fNZvtxBins;

  std::vector<std::vector<EventPool*> > fEvPool; // [fMultBin][fZvtxBin]

  ClassDef(AutoCorr,1)
};
#endif
