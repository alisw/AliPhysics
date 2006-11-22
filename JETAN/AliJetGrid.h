#ifndef ALIJETGRID_H
#define ALIJETGRID_H
/* Copyright(c) 2001-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class description :
//
// Author : Magali Estienne, IPHC Strasbourg - e-mail: magali.estienne@ires.in2p3.fr
//
// --- Standard library ---
#include <Riostream.h>
//-- Root headers ---
#include <TNamed.h>
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TArrayI.h>
//-------------------

class AliJetGrid : public TNamed {

 public:
  
  AliJetGrid();
  AliJetGrid(Int_t nphi,Int_t neta,Double_t phiMin,Double_t etaMin,Double_t phiMax,Double_t etaMax);
  AliJetGrid(const AliJetGrid& grid);
  virtual ~AliJetGrid();

  // Getter
  TArrayD*  GetArrayEta();
  TArrayD*  GetArrayPhi();
  TMatrixD* GetIndexObject();
  void      InitParams(Double_t phiMinCal,Double_t phiMaxCal,Double_t etaMinCal,Double_t etaMaxCal);
  Int_t     GetIndexFromEtaPhi(Double_t phi,Double_t eta) const;
  void      GetEtaPhiFromIndex(Int_t index,Float_t &eta,Float_t &phi);
  Int_t     GetBinsEta() const {return fNeta+1;}
  Int_t     GetBinsPhi() const {return fNphi+1;}
  Int_t     GetTotBins() const {return (fNeta+1)*(fNphi+1);}
  Int_t     GetPointsEta() const {return fNeta;}
  Int_t     GetPointsPhi() const {return fNphi;}
  Int_t     GetTotPoints() const {return fNeta*fNphi;}
  Int_t     GetMatrixIndex(Int_t iphi,Int_t ieta) {return static_cast<Int_t>((*fIndex)(iphi,ieta));}
  Int_t     GetIndexIFromPhi(Double_t phi);
  Int_t     GetIndexJFromEta(Double_t eta);
  Int_t     GetIndex(Double_t phi,Double_t eta);
  void      GetAccParam(Int_t &nphi, Int_t &neta, Float_t &minphi, 
			Float_t &maxphi, Float_t &mineta, Float_t &maxeta);
  void      GetBinParam(Int_t &phibintpc, Int_t &etabintpc, 
			Int_t &phibinemc, Int_t &etabinemc, Int_t &nbinphi);
  void      GetIJFromIndex(Int_t index, Int_t i, Int_t j);
  void      GetEtaPhiFromIndex2(Int_t index, Float_t &phi, Float_t &eta);
  Int_t     GetNEntries();
  Int_t     GetNEntries2();
  Int_t     GetDeta() {return static_cast<Int_t>((fEtaMax-fEtaMin)/fNeta); 
    if(fDebug>21) cout << "static_cast<Int_t>((fEtaMax-fEtaMin)/fNeta) : " << 
      static_cast<Int_t>((fEtaMax-fEtaMin)/fNeta);}
  Int_t     GetDphi() {return static_cast<Int_t>((fPhiMax-fPhiMin)/fNphi); 
    if(fDebug>21) cout << "static_cast<Int_t>((fPhiMax-fPhiMin)/fNphi) : " << 
      static_cast<Int_t>((fPhiMax-fPhiMin)/fNphi);}
  Int_t     GetGridType() {return fGrid;}

  // Setter
  void      SetEtaRange(Double_t etaMin, Double_t etaMax) {fEtaMin = etaMin; fEtaMax = etaMax;}
  void      SetPhiRange(Double_t phiMin, Double_t phiMax) {fPhiMin = phiMin; fPhiMax = phiMax;}
  void      SetNeta(Int_t neta) {fNeta = neta;}
  void      SetNphi(Int_t nphi) {fNphi = nphi;}
  void      SetMatrixIndexes();
  void      SetMatrixIndex(Int_t i,Double_t par);
  void      SetMatrixIndex(Int_t iphi,Int_t ieta,Double_t par) {
                                 (*fIndex)(iphi,ieta)=par; return; }
  void      SetGridType(Int_t type) {fGrid = type;}
  void      SetIndexIJ();

 private:
  Int_t     fGrid;              // Close the type of grid you want to fill
                                // 0 = grid in tpc acceptance, 1 = grid in (tpc-emcal) acceptance  
  Int_t     fNphi;              // number of points in the grid in phi
  Int_t     fNeta;              //               "                 eta
  TArrayD*  fPhi;               // grid points in phi
  TArrayD*  fEta;               // grid points in eta
  TMatrixD* fIndex;             // matrix of indexes in the grid points  
  TArrayI*  fIndexI;            // grid points in phi
  TArrayI*  fIndexJ;            // grid points in eta
  Double_t  fPhiMin;            // phi acceptance min
  Double_t  fPhiMax;            // phi acceptance max
  Double_t  fEtaMin;            // eta acceptance min 
  Double_t  fEtaMax;            // eta acceptance max
  Int_t     fEtaBinInTPCAcc;    // number of points in TPC acceptance in eta
  Int_t     fPhiBinInTPCAcc;    // number of points in TPC acceptance in phi
  Int_t     fEtaBinInEMCalAcc;  // number of points in EMCal acceptance in eta
  Int_t     fPhiBinInEMCalAcc;  // number of points in EMCal acceptance in phi
  Int_t     fNbinEta;
  Int_t     fNbinPhi;
  Double_t  fMaxPhi;
  Double_t  fMinPhi;
  Double_t  fMaxEta;
  Double_t  fMinEta;
  Int_t     fDebug;

  ClassDef(AliJetGrid,1) // Parameters used by AliTPCtrackerParam 
};


#endif




