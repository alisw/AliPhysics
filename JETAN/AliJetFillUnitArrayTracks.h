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

class AliJetHadronCorrection;
class AliJetReader;
class AliJetESDReader;
class TClonesArray;
class TRefArray;
class AliJetGrid;
class AliJetDummyGeo;
class AliESD;
class AliESDEvent;

class AliJetFillUnitArrayTracks : public TTask
{
 public: 
  AliJetFillUnitArrayTracks();
  AliJetFillUnitArrayTracks(AliESD *fESD);
  AliJetFillUnitArrayTracks(AliESDEvent *fESD);
  virtual ~AliJetFillUnitArrayTracks();
  
  // Setter
  void SetReaderHeader(AliJetReaderHeader *readerHeader) {fReaderHeader = readerHeader;}
  void SetGeom(AliJetDummyGeo *geom) {fGeom = geom;}
  void SetMomentumArray(TClonesArray *momentumArray) {fMomentumArray = momentumArray;}
  void SetUnitArray(TClonesArray *unitArray) {fUnitArray = unitArray;}
  void SetRefArray(TRefArray *refArray) {fRefArray = refArray;}
  void SetHadCorrection(Int_t flag = 1) {fHCorrection = flag;}
  void SetHadCorrector(AliJetHadronCorrectionv1* corr) {fHadCorr = corr;}
  void SetTPCGrid(AliJetGrid *grid) {fTPCGrid = grid;}
  void SetEMCalGrid(AliJetGrid *grid) {fEMCalGrid = grid;}
  void SetGrid(Double_t phiMin,Double_t phiMax,Double_t etaMin,Double_t etaMax);
  //  void SetESD(AliESD *esd) {fESD = esd;}
  void SetESD(AliESDEvent *esd) {fESD = esd;}
  void SetGrid0(AliJetGrid *grid0){fGrid0 = grid0;}
  void SetGrid1(AliJetGrid *grid1){fGrid1 = grid1;}
  void SetGrid2(AliJetGrid *grid2){fGrid2 = grid2;}
  void SetGrid3(AliJetGrid *grid3){fGrid3 = grid3;}
  void SetGrid4(AliJetGrid *grid4){fGrid4 = grid4;}

  // Getter
  TClonesArray* GetUnitArray() {return fUnitArray;}
  TRefArray*    GetRefArray() {return fRefArray;}
  //  Int_t         GetIndexFromEtaPhi(Double_t eta,Double_t phi) const;
  void          GetEtaPhiFromIndex(Int_t index,Float_t &eta,Float_t &phi);
  Int_t         GetNeta()          const {return fNeta;}
  Int_t         GetNphi()          const {return fNphi;}
  Int_t         GetHadCorrection() const {return fHCorrection;}
  Int_t         GetMult()          const {return fNTracks;}
  Int_t         GetMultCut()       const {return fNTracksCut;}
  void          Exec(Option_t*);

 protected:
  Int_t   fNumUnits;      // Number of units in the unit object array (same as num towers in EMCAL)
  Float_t fEtaMinCal;     // Define EMCal acceptance in Eta
  Float_t fEtaMaxCal;     // Define EMCal acceptance in Eta
  Float_t fPhiMinCal;     // Define EMCal acceptance in Phi
  Float_t fPhiMaxCal;     // Define EMCal acceptance in Phi
  AliJetHadronCorrectionv1   *fHadCorr;         // Pointer to Hadron Correction Object
  Int_t                       fHCorrection;     //  Hadron correction flag
  Int_t                       fNTracks;         // Number of tracks stored in UnitArray
  Int_t                       fNTracksCut;      // Number of tracks stored in UnitArray with a pt cut 
  Int_t                       fOpt;             // Detector to be used for jet reconstruction
  Bool_t                      fDZ;              // Use or not dead zones
  Int_t                       fDebug;           // Debug option

  AliJetReaderHeader          *fReaderHeader;   // ReaderHeader
  TClonesArray                *fMomentumArray;  // MomentumArray
  TClonesArray                *fUnitArray;      // UnitArray
  TRefArray                   *fRefArray;       // UnitArray
  AliJetGrid                  *fTPCGrid;        // Define filled grid
  AliJetGrid                  *fEMCalGrid;      // Define filled grid
  AliJetDummyGeo              *fGeom;           // Define EMCal geometry
  AliESDEvent                 *fESD;            // ESD
  AliJetGrid                  *fGrid0;          // Pointer to Grid 1
  AliJetGrid                  *fGrid1;          // Pointer to Grid 2
  AliJetGrid                  *fGrid2;          // Pointer to Grid 3
  AliJetGrid                  *fGrid3;          // Pointer to Grid 4
  AliJetGrid                  *fGrid4;          // Pointer to Grid 5  

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
  Float_t   fPhiMin;                  // minimum phi
  Float_t   fPhiMax;                  // maximum phi
  Float_t   fEtaMin;                  // minimum eta
  Float_t   fEtaMax;                  // maximum eta   
  Int_t     fEtaBinInTPCAcc;          // eta bins in tpc acceptance
  Int_t     fPhiBinInTPCAcc;          // phi bins in tpc acceptance 
  Int_t     fEtaBinInEMCalAcc;        // eta bins in emcal acceptance
  Int_t     fPhiBinInEMCalAcc;        // phi bins in emcal acceptance
  Int_t     fNbinPhi;                 // number of phi bins

 private:
  //  void SetEMCALGeometry();
  void InitParameters();

  ClassDef(AliJetFillUnitArrayTracks,1) // Fill Unit Array with tpc and/or emcal information
};

#endif
