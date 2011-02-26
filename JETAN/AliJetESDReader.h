#ifndef ALIJETESDREADER_H
#define ALIJETESDREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet ESD Reader 
// ESD reader for jet analysis
// Author: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)
//=========================================================================
// Modified in order to use a fUnitArray object instead of a fMomentumArray
// Includes EMCal Geometry, fUnitArray, grid objects and tools for Hadron correction
// Author : magali.estienne@subatech.in2p3.fr
//---------------------------------------------------------------------

#include "AliJetReader.h"
#include "AliJetUnitArray.h"
#include "AliJetGrid.h"

class TRefArray;
class AliJetESDReaderHeader;
class AliEMCALGeometry;
class AliEMCALGeoUtils;
class AliJetHadronCorrection;
class AliJetUnitArray;
class AliJetReaderHeader;
class AliESDEvent;

class AliJetESDReader : public AliJetReader
{
 public: 
  AliJetESDReader();
  virtual ~AliJetESDReader();

  Bool_t     FillMomentumArray(); 
  void       OpenInputFiles();
  void       InitUnitArray();
  void       CreateTasks(TChain* tree);
  Bool_t     ExecTasks(const Bool_t procid, TRefArray* refArray);

  // Getters
  Float_t    GetTrackMass() const {return fMass;}  // returns mass of the track
  Int_t      GetTrackSign() const {return fSign;}  // returns sign of the track
 
  // Setters
  void       SetInputEvent(const TObject* esd, const TObject* aod, const TObject* mc);
  void       SetTPCGrid(AliJetGrid *grid)   {fTpcGrid = grid;}
  void       SetEMCalGrid(AliJetGrid *grid) {fEmcalGrid = grid;}
  // Correction of hadronic energy
  void       SetHadronCorrector(AliJetHadronCorrection* corr) {fHadCorr = corr;}
  void       SetApplyElectronCorrection(Int_t flag = 1) {fECorrection = flag; fEFlag=kTRUE;}
  void       SetApplyMIPCorrection(Bool_t val);
  void       SetApplyFractionHadronicCorrection(Bool_t val);
  void       SetFractionHadronicCorrection(Double_t val);

 protected:
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
  Int_t                       fDebug;            //! Debug option
  Float_t                     fMass;             // Particle mass
  Int_t                       fSign;             // Particle sign
  Int_t                       fNIn;              // Number of Array filled in UnitArray
  Int_t                       fOpt;              // Detector to be used for jet reconstruction
  Bool_t                      fDZ;               // Use or not dead zones
  Int_t                       fNeta;             // Number of bins in eta of tpc grid
  Int_t                       fNphi;             // Number of bins in phi of tpc grid
  TRefArray                  *fRefArray;         // array of digit position and energy
  Bool_t                      fProcId;           // Bool_t for TProcessID synchronization
  void InitParameters();
  AliJetESDReader(const AliJetESDReader &det);
  AliJetESDReader &operator=(const AliJetESDReader &det);



  ClassDef(AliJetESDReader,2)
};
 
#endif
