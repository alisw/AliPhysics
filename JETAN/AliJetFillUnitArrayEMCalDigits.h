#ifndef ALIJETFILLUNITARRAYEMCALDIGITS_H
#define ALIJETFILLUNITARRAYEMCALDIGITS_H
 
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

class AliEMCALGeometry;
class AliJetDummyGeo;
class AliESDCaloCluster;
class AliJetReader;
class AliJetESDReader;
class TClonesArray;
class TRefArray;
class AliJetUnitArray;
class AliESDEvent;
class AliJetGrid;
class AliEMCALCalibData ;

class AliJetFillUnitArrayEMCalDigits : public TTask
{
 public: 
  AliJetFillUnitArrayEMCalDigits();
  AliJetFillUnitArrayEMCalDigits(AliESDEvent *fESD);
  virtual ~AliJetFillUnitArrayEMCalDigits();
  
  // Setter
  void SetReaderHeader(AliJetReaderHeader *readerHeader) {fReaderHeader = readerHeader;}
  void SetGeom(AliJetDummyGeo *geom) {fGeom = geom;}
  void SetMomentumArray(TClonesArray *momentumArray) {fMomentumArray = momentumArray;}
  void SetUnitArray(TClonesArray *unitArray) {fUnitArray = unitArray;}
  void SetRefArray(TRefArray *refArray) {fRefArray = refArray;}
  void SetTPCGrid(AliJetGrid *grid) {fTPCGrid = grid;}
  void SetEMCalGrid(AliJetGrid *grid) {fEMCalGrid = grid;}
  void SetEleCorrection(Int_t flag = 1) {fECorrection = flag;}
  void SetESD(AliESDEvent *esd) {fESD = esd;}
  void SetInitMult(Int_t mult) {fNDigitEmcal = mult;}
  void SetInitMultCut(Int_t multcut) {fNDigitEmcalCut = multcut;}
  void SetProcId(Bool_t id) {fProcId = id;}

  // Getter
  TClonesArray* GetUnitArray() {return fUnitArray;}
  TRefArray*    GetRefArray()  {return fRefArray;}
  Int_t         GetMult()      {return fNDigitEmcal;}
  Int_t         GetMultCut()   {return fNDigitEmcalCut;}
  // For calibration => Under construction
  //  virtual Float_t Calibrate(Int_t amp, Int_t cellId) ;  // Tranforms Amp to energy

  // Other
  void          Exec(Option_t*);
  Float_t       EtaToTheta(Float_t arg);

 protected:
  AliESDEvent        *fESD;            // ESD
  Int_t               fNumUnits;       // Number of units in the unit object array (same as num towers in EMCAL)
  Float_t             fEtaMinCal;      // Define EMCAL acceptance in Eta
  Float_t             fEtaMaxCal;      // Define EMCAL acceptance in Eta
  Float_t             fPhiMinCal;      // Define EMCAL acceptance in Phi
  Float_t             fPhiMaxCal;      // Define EMCAL acceptance in Phi
  Int_t               fNIn;            // Number of Array filled in UnitArray
  Int_t               fOpt;            // Detector to be used for jet reconstruction
  Int_t               fCluster;        // Use all cells or cells in clusters for jet finding 
  Int_t               fDebug;          // Debug option
  Int_t               fNCEMCAL;        // Number of clusters in EMCAL
  Int_t               fNCPHOS;         // Number of clusters in PHOS
  Int_t               fNCCalo;         // Number of cluster in EMCAL + PHOS calorimeters

  AliJetGrid         *fTPCGrid;        // Define filled grid
  AliJetGrid         *fEMCalGrid;      // Define filled grid
  Int_t               fECorrection;    // Electron correction flag

  AliJetReaderHeader *fReaderHeader;   // ReaderHeader
  TClonesArray       *fMomentumArray;  // MomentumArray
  TClonesArray       *fUnitArray;      // UnitArray
  TRefArray          *fRefArray;       // UnitArray
  Bool_t              fProcId;         // Bool_t for TProcessID synchronization
  AliJetDummyGeo     *fGeom;           // Set EMCal geometry
  
  AliESDCaloCluster  *fClus;           //! 
  Int_t               fNDigitEmcal;    //!
  Int_t               fNDigitEmcalCut; //!
  //Calibration parameters... to be replaced by database
  AliEMCALCalibData  *fCalibData;      //! Calibration database if aval
  Float_t             fADCchannelECA;  // width of one ADC channel for EC section (GeV)
  Float_t             fADCpedestalECA; // pedestal of ADC for EC section (GeV)

 private:
  AliJetFillUnitArrayEMCalDigits(const AliJetFillUnitArrayEMCalDigits &det);
  AliJetFillUnitArrayEMCalDigits &operator=(const AliJetFillUnitArrayEMCalDigits &det);

  void InitParameters();
  // Under construction
  //  void    GetCalibrationParameters(void) ;
  
  ClassDef(AliJetFillUnitArrayEMCalDigits,1) // Fill Unit Array with tpc and/or emcal information
};

#endif
