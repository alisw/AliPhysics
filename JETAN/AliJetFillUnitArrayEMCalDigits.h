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

class AliJetDummyGeo;
class AliESDCaloCluster;
class AliEMCALCalibData;
class AliJetReader;
class AliJetESDReader;
class TClonesArray;
class TRefArray;
class AliJetUnitArray;
//class AliESD;
class AliESDEvent;
class AliJetGrid;

class AliJetFillUnitArrayEMCalDigits : public TTask
{
 public: 
  AliJetFillUnitArrayEMCalDigits();
  AliJetFillUnitArrayEMCalDigits(Int_t event);
  AliJetFillUnitArrayEMCalDigits(AliESD *fESD);
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
  //  void SetESD(AliESD *esd) {fESD = esd;}
  void SetESD(AliESDEvent *esd) {fESD = esd;}
  void SetInitMult(Int_t mult) {fNDigitEmcal = mult;}
  void SetInitMultCut(Int_t multcut) {fNDigitEmcalCut = multcut;}

  // Getter
  TClonesArray* GetUnitArray() {return fUnitArray;}
  TRefArray*    GetRefArray()  {return fRefArray;}
  Int_t         GetMult()      {return fNDigitEmcal;}
  Int_t         GetMultCut()   {return fNDigitEmcalCut;}

  // Other
  void          Exec(Option_t*);
  Float_t       EtaToTheta(Float_t arg);
 private:
  void InitParameters();

 protected:
  AliESDEvent  *fESD; // ESD
  Int_t   fNumUnits;  // Number of units in the unit object array (same as num towers in EMCAL)
  Float_t fEtaMinCal; // Define EMCAL acceptance in Eta
  Float_t fEtaMaxCal; // Define EMCAL acceptance in Eta
  Float_t fPhiMinCal; // Define EMCAL acceptance in Phi
  Float_t fPhiMaxCal; // Define EMCAL acceptance in Phi
  Int_t   fNIn;       // Number of Array filled in UnitArray
  Int_t   fOpt;       // Detector to be used for jet reconstruction
  Int_t   fDebug;     // Debug option
  Int_t   fNCEMCAL;   // Number of clusters in EMCAL
  Int_t   fNCPHOS;    // Number of clusters in PHOS
  Int_t   fNCCalo;    // Number of cluster in EMCAL + PHOS calorimeters

  AliJetGrid         *fTPCGrid;           // Define filled grid
  AliJetGrid         *fEMCalGrid;         // Define filled grid

  AliJetReaderHeader *fReaderHeader;      // ReaderHeader
  TClonesArray       *fMomentumArray;     // MomentumArray
  TClonesArray       *fUnitArray;         // UnitArray
  TRefArray          *fRefArray;          // UnitArray
  AliJetDummyGeo     *fGeom;              // Set EMCal geometry
  
  AliESDCaloCluster *fClus;               //! 
  Int_t fNDigitEmcal;                     //!
  Int_t fNDigitEmcalCut;                  //!

  

  ClassDef(AliJetFillUnitArrayEMCalDigits,1) // Fill Unit Array with tpc and/or emcal information
};

#endif
