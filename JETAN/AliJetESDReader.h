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
// Author : magali.estienne@ires.in2p3.fr
//---------------------------------------------------------------------

#include "AliJetReader.h"
#include "AliJetUnitArray.h"
#include "AliJetGrid.h"
class AliJetESDReaderHeader;
class AliJetDummyGeo;
class AliJetHadronCorrection;
class AliJetUnitArray;
class AliJetReaderHeader;
class AliESDEvent;

class AliJetESDReader : public AliJetReader
{
 public: 
  AliJetESDReader();
  virtual ~AliJetESDReader();

  // Getters
  Float_t GetTrackMass() const {return fMass;}  // returns mass of the track
  Int_t   GetTrackSign() const {return fSign;}  // returns sign of the track

  // Setters
  Bool_t FillMomentumArray(Int_t event); 
  void   OpenInputFiles();
  void   InitUnitArray();
  void   CreateTasks();
  //  void   ExecTasks(Int_t event);
  Bool_t   ExecTasks(Int_t event);
  void   SetInputEvent(TObject* esd, TObject* aod, TObject* mc);
  virtual void SetTPCGrid(AliJetGrid *grid)   {fTpcGrid = grid;}
  virtual void SetEMCalGrid(AliJetGrid *grid) {fEmcalGrid = grid;}
  // Correction of hadronic energy
  virtual void SetHadronCorrection(Int_t flag = 1) {fHCorrection = flag;}
  virtual void SetHadronCorrector(AliJetHadronCorrectionv1* corr) {fHadCorr = corr;}
 private:
  void SetEMCALGeometry();
  void InitParameters();
 protected:
  AliJetDummyGeo             *fGeom;             //! EMCAL Geometry 
  TChain                     *fChain;            //! chain for reconstructed tracks
  AliESDEvent                *fESD;              //! pointer to esd
  AliJetHadronCorrectionv1   *fHadCorr;          //! Pointer to Hadron Correction Object 
  AliJetGrid                 *fTpcGrid;          //! Pointer to grid object
  AliJetGrid                 *fEmcalGrid;        //! Pointer to grid object
  AliJetGrid                 *fGrid0;            // Pointer to grid object
  AliJetGrid                 *fGrid1;            // Pointer to grid object
  AliJetGrid                 *fGrid2;            // Pointer to grid object
  AliJetGrid                 *fGrid3;            // Pointer to grid object
  AliJetGrid                 *fGrid4;            // Pointer to grid object
  Float_t                     fPtCut;            // Pt cut for tracks to minimise background contribution
  Int_t                       fHCorrection;      // Hadron correction flag
  Int_t                       fNumUnits;         // Number of units in the unit object array
                                                 // (same as num towers in EMCAL)
  Int_t                       fDebug;            // Debug option
  Float_t                     fMass;    // Particle mass
  Int_t                       fSign;    // Particle sign
  Int_t                       fNIn;              // Number of Array filled in UnitArray
  Int_t                       fOpt;              // Detector to be used for jet reconstruction
  Bool_t                      fDZ;               // Use or not dead zones
  Int_t                       fNeta;             // Number of bins in eta of tpc grid
  Int_t                       fNphi;             // Number of bins in phi of tpc grid
  Bool_t                      fArrayInitialised; // To check that array of units is initialised
  
  ClassDef(AliJetESDReader,1)
};
 
#endif
