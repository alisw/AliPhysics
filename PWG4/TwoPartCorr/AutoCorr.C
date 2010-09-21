// $Id$

#include "AutoCorr.h"

ClassImp(AutoCorr)

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

EventPool* AutoCorr::GetEventPool(int iMult, int iZvtx) const
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

Double_t AutoCorr::DeltaPhi(const MyPart &t1, const MyPart &t2,
			    Double_t rangeMin, Double_t rangeMax) const
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();
  Double_t phia = t1.fPhi;  
  Double_t phib = t2.fPhi;  
  
  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += 2*pi;
  else if (dphi > rangeMax) dphi -= 2*pi;
  
  return dphi;
}

Double_t AutoCorr::DeltaEta(const MyPart &t1, const MyPart &t2) const
{
  return t1.fEta - t2.fEta;
}

Bool_t AutoCorr::InBounds(Double_t val, Double_t min, Double_t max) const
{
  if (val<min)
    return 0;
  if (val>max)
    return 0;
  return 1;
}

Bool_t AutoCorr::InBounds(Int_t val, Int_t min, Int_t max) const
{
  if (val<min)
    return 0;
  if (val>max)
    return 0;
  return 1;
}

Bool_t AutoCorr::IsEventOk(const MyHeader &ev, Int_t minVc, 
			   Int_t maxNTracklets, Double_t zMin, Double_t zMax) const
{
  Bool_t VcOk = ev.fVc >= minVc;
  Bool_t NTrackletsOK = ev.fNTracklets <= maxNTracklets;
  Bool_t zOk = InBounds(ev.fVz, zMin, zMax);
  return (!ev.fIsPileupSPD && VcOk && NTrackletsOK && zOk);
}

Bool_t AutoCorr::IsTrackOk(const MyPart &t, Double_t etaMin, Double_t etaMax) const
{
  return InBounds(t.fEta, etaMin, etaMax);    
}

Bool_t AutoCorr::IsTrackOk(const MyPart &t, Double_t etaMin, Double_t etaMax,
			   Double_t ptMin, Double_t ptMax) const
{
  Bool_t etaOk = InBounds(t.fEta, etaMin, etaMax);
  Bool_t ptOk  = InBounds(t.fPt, ptMin, ptMax);    
  return  etaOk && ptOk;
}

Bool_t AutoCorr::IsPairOk(const MyPart &t1, const MyPart &t2) const
{
  Double_t deta = DeltaEta(t1, t2);
  Double_t dphi = DeltaPhi(t1, t2);
  Double_t dpmax = 0.03;
  Double_t demax = 0.01;
  Double_t dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
  return (dr > 1);
}

Bool_t AutoCorr::IsMixedPairOk(const MyPart &t1, const MyPart &t2) const
{
  Double_t deta = DeltaEta(t1, t2);
  Double_t dphi = DeltaPhi(t1, t2);
  Double_t dpmax = 0.04;
  Double_t demax = 0.04;
  Double_t dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
  return (dr > 1);
}

Double_t AutoCorr::InvMass(const MyPart &p1, const MyPart &p2) const
{
  Double_t px1 = p1.Px();
  Double_t py1 = p1.Py();
  Double_t pz1 = p1.Pz();
  Double_t px2 = p2.Px();
  Double_t py2 = p2.Py();
  Double_t pz2 = p2.Pz();
  Double_t pm1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
  Double_t pm2 = TMath::Sqrt(px2*px2+py2*py2+pz1*pz2);
  Double_t p12 = px1*px2+py1*py2+pz1*pz2;
  Double_t m = TMath::Sqrt(pm1*pm2-p12);
  return m;
}
