// $Id:$

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
#include "AutoCorr.h"

Int_t AutoCorr::InitEventPools(Int_t depth, 
			       Int_t nMultBins, Double_t multbin[], 
			       Int_t nZvtxBins, Double_t zvtxbin[])
{
  // First assign AutoCorr members
  fNMultBins = nMultBins;
  fNZvtxBins = nZvtxBins;

  for (int iM=0; iM<nMultBins; iM++) {
    std::vector<EventPool*> evp;
    for (int iZ=0; iZ<nZvtxBins; iZ++) {
      evp.push_back(new EventPool(depth, 
				  multbin[iM], multbin[iM+1], 
				  zvtxbin[iZ], zvtxbin[iZ+1] ));
    }
    fEvPool.push_back(evp);
  }
  
  for (int iM=0; iM<nMultBins; iM++) {
    for (int iZ=0; iZ<nZvtxBins; iZ++) {
      fEvPool.at(iM).at(iZ)->SetMultBinIndex(iM);
      fEvPool.at(iM).at(iZ)->SetZvtxBinIndex(iZ);
    }
  }
  
  bool print_this = false;
  if (print_this) {
    cout << "fEvPool outer size: " << fEvPool.size() << endl;
    for (int iM=0; iM<nMultBins; iM++) {
      for (int iZ=0; iZ<nZvtxBins; iZ++) {
	if(fEvPool.at(iM).at(iZ)) {
	  cout << "multiplicity bin: " << iM;
	  cout << ", z-vertex bin: " << iZ;
	  fEvPool.at(iM).at(iZ)->PrintInfo();
	}
      }
    }
  }
  
  return 0;
}

EventPool* AutoCorr::GetEventPool(int iMult, int iZvtx)
{
  if (iMult < 0 || iMult >= fNMultBins) return 0x0;
  if (iZvtx < 0 || iZvtx >= fNZvtxBins) return 0x0;
  return fEvPool.at(iMult).at(iZvtx);
}

Int_t AutoCorr::UpdatePools(int iEvent, MyHeader* ev, TClonesArray* trk)
{
  for (int iM=0; iM<fNMultBins; iM++) {
    for (int iZ=0; iZ<fNZvtxBins; iZ++) {
      fEvPool.at(iM).at(iZ)->UpdatePool(iEvent, ev, trk);
    }
  }  
  return 0;
}


Double_t AutoCorr::DeltaPhi(MyPart* t1, MyPart* t2,
			    double rangeMin, double rangeMax)
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();
  
  if (!t1 || !t2) return -99.;
  double phia = t1->Phi();  
  double phib = t2->Phi();  
  
  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += 2*pi;
  else if (dphi > rangeMax) dphi -= 2*pi;
  
  return dphi;
}

Double_t AutoCorr::DeltaEta(MyPart* t1, MyPart* t2)
{
  if (!t1 || !t2) return -99.;
  return t1->Eta() - t2->Eta();
}

Bool_t AutoCorr::IsTrackOk(MyPart* t)
{
  if (!t) return false;
  return fabs(t->Eta()) <= 1.2;    
}

Bool_t AutoCorr::IsPairOk(MyPart* t1, MyPart* t2)
{
  if (!t1 || !t2) return false;
  if (!IsTrackOk(t1) || !IsTrackOk(t2)) return false;

  double deta = DeltaEta(t1, t2);
  double dphi = DeltaPhi(t1, t2);
  double dpmax = 0.03;
  double demax = 0.01;

  double dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
  return (dr > 1);
}

// So far same as IsPairOk()
Bool_t AutoCorr::IsMixedPairOk(MyPart* t1, MyPart* t2)
{
  if (!t1 || !t2) return false;
  if (!IsTrackOk(t1) || !IsTrackOk(t2)) return false;

  double deta = DeltaEta(t1, t2);
  double dphi = DeltaPhi(t1, t2);
  double dpmax = 0.03;
  double demax = 0.01;

  double dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
  return (dr > 1);
}

