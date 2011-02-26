#ifndef ALIJETAODREADER_H
#define ALIJETAODREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet AOD Reader
// AOD reader for jet analysis
// Author: Davide Perrino (davide.perrino@cern.ch)
//---------------------------------------------------------------------

#include "AliJetReader.h"
#include "AliJetUnitArray.h"
#include "AliJetGrid.h"
class AliJetUnitArray;
class AliEMCALGeoUtils;
class AliJetHadronCorrection;
class AliJetAODReaderHeader;
class AliJetReaderHeader;
class AliAODMCParticle;
class AliAODEvent;
class TRefArray;

class AliJetAODReader : public AliJetReader
{
 public: 
  AliJetAODReader();
  virtual ~AliJetAODReader();

  TRefArray*   GetReferences() const {return fRef;}

  Bool_t FillMomentumArray(); 
  static Bool_t AcceptAODMCParticle(AliAODMCParticle *mcP,Short_t flag);

  void   OpenInputFiles();
  void   ConnectTree(TTree* tree, TObject* data);
  void   InitUnitArray();
  void   CreateTasks(TChain* tree);
  Bool_t ExecTasks(Bool_t procid, TRefArray* refArray);

  void SetInputEvent(const TObject* /*esd*/, const TObject* aod, const TObject* /*mc*/) {fAOD = (AliAODEvent*) aod;}
  void SetTPCGrid(AliJetGrid *grid)   {fTpcGrid = grid;}
  void SetEMCalGrid(AliJetGrid *grid) {fEmcalGrid = grid;}
  // Correction of hadronic energy
  void SetHadronCorrection(const Int_t flag = 1) {fHCorrection = flag;}
  void SetHadronCorrector(AliJetHadronCorrection* corr) {fHadCorr = corr;}
  void SetApplyElectronCorrection(const Int_t flag = 1) {fECorrection = flag; fEFlag=kTRUE;}
  void SetApplyMIPCorrection(const Bool_t val);
  void SetApplyFractionHadronicCorrection(Bool_t val);
  void SetFractionHadronicCorrection(Double_t val);


 private:
  Bool_t SetEMCALGeometry();
  Bool_t FillMomentumArrayMC();


  void InitParameters();
  AliJetAODReader(const AliJetAODReader &det);
  AliJetAODReader &operator=(const AliJetAODReader &det);

 private:
  AliAODEvent                *fAOD;    //! pointer to aod
  TRefArray                  *fRef;    //! pointer to array of references to tracks
  Int_t                       fDebug;  // Debug option
  Int_t                       fOpt;    // Detector to be used for jet reconstruction
  AliJetHadronCorrection     *fHadCorr;          //! Pointer to Hadron Correction Object
  AliJetGrid                 *fTpcGrid;          //! Pointer to grid object
  AliJetGrid                 *fEmcalGrid;        //! Pointer to grid object
  AliJetGrid                 *fGrid0;            // Pointer to grid object
  AliJetGrid                 *fGrid1;            // Pointer to grid object
  AliJetGrid                 *fGrid2;            // Pointer to grid object
  AliJetGrid                 *fGrid3;            // Pointer to grid object
  AliJetGrid                 *fGrid4;            // Pointer to grid object
  Int_t                       fApplyElectronCorrection;      // Electron correction flag
  Bool_t                      fApplyMIPCorrection; // Apply MIP or not ? Exclusive with fApplyFractionHadronicCorrection
  Bool_t                      fApplyFractionHadronicCorrection; // Another type of charged particle energy deposition in EMC
  Double_t                    fFractionHadronicCorrection; // Fraction of momentum of the TPC track to be subtracted from EMC tower
  Int_t                       fNumUnits;         // Number of units in the unit object array
                                                 // (same as num towers in EMCAL)
  Float_t                     fMass;             // Particle mass
  Int_t                       fSign;             // Particle sign
  Int_t                       fNIn;              // Number of Array filled in UnitArray
  Bool_t                      fDZ;               // Use or not dead zones
  Int_t                       fNeta;             // Number of bins in eta of tpc grid
  Int_t                       fNphi;             // Number of bins in phi of tpc grid
  TRefArray                  *fRefArray;         // array of digit position and energy
  Bool_t                      fProcId;           // Bool_t for TProcessID synchronization
  ClassDef(AliJetAODReader,2)
};
 
#endif
