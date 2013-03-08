#ifndef ALITRDTENDERSUPPLY_H
#define ALITRDTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  TRD tender, reapply pid on the fly                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////



#include <AliTenderSupply.h>

class AliTRDCalDet;
class AliESDEvent;
class AliOADBContainer;
class AliTRDonlineTrackMatching;

class AliTRDTenderSupply: public AliTenderSupply {
  
public:
  enum{
    kNNpid = 0,
    k1DLQpid = 1,
    k2DLQpid = 2
  };
  AliTRDTenderSupply();
  AliTRDTenderSupply(const char *name, const AliTender *tender=NULL);
  virtual ~AliTRDTenderSupply();

  void SetRunByRunCorrection(const char *filename) { fNameRunByRunCorrection = filename; }
//  void SetLoadReferencesFromCDB() { fLoadReferences = kTRUE; fLoadReferencesFromCDB = kTRUE; }
//  void SetLoadReferencesFromFile() { fLoadReferences = kTRUE; fLoadReferencesFromCDB = kFALSE; }
  void SetLoadDeadChambersFromCDB(){ fLoadDeadChambers = kTRUE;} ;
  void SetPIDmethod(Int_t pidMethod) { fPIDmethod = pidMethod; }
  void SetNormalizationFactor(Double_t norm, Int_t runMin, Int_t runMax);
  void SetNormalizationFactor(Double_t norm) { fNormalizationFactor = norm; }
  void SetCalibLowpThreshold(Double_t pmin) { fPthreshold = pmin; };
  void SetGeoFile(const char *filename) { fGeoFile = filename; }
  void SetDebugMode() { fDebugMode = kTRUE; }
  void SetRedoTRDMatching(Bool_t redo = kTRUE) {fRedoTrdMatching = redo;}

  virtual void              Init();
  virtual void              ProcessEvent();
  
  void SwitchOnGainCorrection() { fGainCorrection = kTRUE; }
  void SwitchOffGainCorrection() { fGainCorrection = kFALSE; }
  void SetSlicesForPID(UInt_t min, UInt_t max) { fSlicesForPID[0] = min; fSlicesForPID[1] = max;}
  void AddBadChamber(Int_t chamberID){fBadChamberID[fNBadChambers++] = chamberID;};
  
private:
  enum{
    kNPlanes = 6,
    kNStacks = 5,
    kNChambers = 540
  };

  Bool_t GetTRDchamberID(AliESDtrack * const track, Int_t *detectors);
  void SetChamberGain();
  void ApplyGainCorrection(AliESDtrack *track, const Int_t * const detectors);
  void ApplyRunByRunCorrection(AliESDtrack *const track);
  void MaskChambers(AliESDtrack * const track, const Int_t * const detectors); 
  void LoadReferences();
  void LoadDeadChambersFromCDB();
  void LoadRunByRunCorrection(const char *filename);
  Bool_t IsBadChamber(Int_t chamberID);
  Double_t GetNormalizationFactor(Int_t runnumber);
  
  AliESDEvent           *fESD;       //! the ESD Event
  AliESDpid             *fESDpid;    //! ESD PID object
  AliTRDonlineTrackMatching  *fTrdOnlineTrackMatcher;   //! TRD online track matcher

  AliTRDCalDet *fChamberGainOld;     // TRD Chamber Gain Factor used for producing the ESD
  AliTRDCalDet *fChamberGainNew;     // New TRD Chamber Gain Factor
  AliTRDCalDet *fChamberVdriftOld;   // Old drift velocity calibration
  AliTRDCalDet *fChamberVdriftNew;   // New drift velocity calibration
  AliOADBContainer *fRunByRunCorrection;    // Run by run gain correction

  Int_t fPIDmethod;                  // PID method
  Double_t fNormalizationFactor;     // dE/dx Normalization Factor 
  Double_t fPthreshold;              // Low Momentum threshold for calibration
  Int_t fBadChamberID[kNChambers];   // List of Bad Chambers
  UInt_t fSlicesForPID[2];           // Select range of slices used in the PID response
  UInt_t fNBadChambers;              // Number of bad chambers
  const char *fGeoFile;              // File with geometry.root
  Bool_t fGainCorrection;            // Apply gain correction 
//  Bool_t fLoadReferences;            // Tender Load references
//  Bool_t fLoadReferencesFromCDB;     // Load References from CDB
  Bool_t fLoadDeadChambers;          // Load dead chambers
  Bool_t fHasReferences;             // has references loaded
  Bool_t fHasNewCalibration;         // has new calibration
  Bool_t fDebugMode;                 // Run in debug mode
  Bool_t fRedoTrdMatching;           // Redo Track Matching
  TString fNameRunByRunCorrection;   // filename with the run-by-run gain correction
  TObjArray *fNormalizationFactorArray; // Array with normalisation Factors
  
  AliTRDTenderSupply(const AliTRDTenderSupply&c);
  AliTRDTenderSupply& operator= (const AliTRDTenderSupply&c);
  
  ClassDef(AliTRDTenderSupply, 1);  // TRD tender task
};
#endif

