#ifndef ALIJETFILLUNITARRAYTRACKS_H
#define ALIJETFILLUNITARRAYTRACKS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet Fill Unit Array 
// Called by ESD Reader for jet analysis
// Author: Magali Estienne (magali.estienne@subatech.in2p3.fr)
//---------------------------------------------------------------------

#ifndef ROOT_TTask
#include "TTask.h"
#endif

#include <TMatrixD.h>
#include <TArrayD.h>

class AliEMCALGeometry;
class AliJetHadronCorrection;
class AliJetReader;
class AliJetESDReader;
class TClonesArray;
class AliJetUnitArray;
//class AliJetReaderHeader;
class AliJetReader;
class AliJetGrid;

class AliJetFillUnitArrayTracks : public TTask
{
 public: 
  AliJetFillUnitArrayTracks();
  virtual ~AliJetFillUnitArrayTracks();
  
  // Setter
  void SetReaderHeader(AliJetReaderHeader *readerHeader) {fReaderHeader = readerHeader;}
  void SetMomentumArray(TClonesArray *momentumArray) {fMomentumArray = momentumArray;}
  void SetUnitArray(AliJetUnitArray *unitArray) {fUnitArray = unitArray;}
  void SetHadCorrection(Int_t flag = 1) {fHCorrection = flag;}
  void SetHadCorrector(AliJetHadronCorrectionv1* corr) {fHadCorr = corr;}
  void SetTPCGrid(AliJetGrid *grid) {fTPCGrid = grid;}
  void SetEMCalGrid(AliJetGrid *grid) {fEMCalGrid = grid;}
  void SetGrid(Double_t phiMin,Double_t phiMax,Double_t etaMin,Double_t etaMax);

  // Getter
  AliJetUnitArray* GetUnitArray() {return fUnitArray;}
  //  Int_t GetIndexFromEtaPhi(Double_t eta,Double_t phi) const;
  void  GetEtaPhiFromIndex(Int_t index,Float_t &eta,Float_t &phi);
  Int_t GetNeta() {return fNeta;}
  Int_t GetNphi() {return fNphi;}

  void Exec(Option_t*);

 protected:
  Int_t   fNumUnits;      // Number of units in the unit object array (same as num towers in EMCAL)
  Float_t fEtaMinCal;     // Define EMCal acceptance in Eta
  Float_t fEtaMaxCal;     // Define EMCal acceptance in Eta
  Float_t fPhiMinCal;     // Define EMCal acceptance in Phi
  Float_t fPhiMaxCal;     // Define EMCal acceptance in Phi
  AliJetHadronCorrectionv1   *fHadCorr;         // Pointer to Hadron Correction Object
  Int_t                       fHCorrection;     //  Hadron correction flag
  Int_t                       fNIn;             // Number of Array filled in UnitArray
  Int_t                       fOpt;             // Detector to be used for jet reconstruction
  Int_t                       fDebug;           // Debug option

  AliJetReaderHeader          *fReaderHeader;   // ReaderHeader
  TClonesArray                *fMomentumArray;  // MomentumArray
  AliJetUnitArray             *fUnitArray;      // UnitArray
  AliJetGrid                  *fTPCGrid;        // Define filled grid
  AliJetGrid                  *fEMCalGrid;      // Define filled grid

  // geometry info
  static AliEMCALGeometry *fGeom;     //!

  Int_t     fNphi;                    // number of points in the grid:   phi
  Int_t     fNeta;                    //               "                 eta
  TArrayD*  fPhi2;                    // grid points in phi
  TArrayD*  fEta2;                    // grid points in eta
  TArrayD*  fPhi;                     // grid points in phi
  TArrayD*  fEta;                     // grid points in eta
  TMatrixD* fIndex;                   // grid points in (phi,eta) 
  TMatrixD* fParams;                  // matrix of parameters in the grid points  
  Int_t     fGrid;                    // Select the grid acceptance you want to fill
                                      // 0 = TPC acceptance, 1 = TPC-EMCal acceptance
  Float_t   fPhiMin;
  Float_t   fPhiMax;
  Float_t   fEtaMin;
  Float_t   fEtaMax;
  Int_t     fEtaBinInTPCAcc;
  Int_t     fPhiBinInTPCAcc;
  Int_t     fEtaBinInEMCalAcc;
  Int_t     fPhiBinInEMCalAcc;
  Int_t     fNbinPhi;

 private:
  void SetEMCALGeometry();
  void InitParameters();

  //  Int_t fEvent;
  
  ClassDef(AliJetFillUnitArrayTracks,1) // Fill Unit Array with tpc and/or emcal information
};

#endif
