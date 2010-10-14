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
class AliESDtrack;
class AliESDEvent;

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

  void SetPIDmethod(Int_t pidMethod) { fPIDmethod = pidMethod; }
  void SetCalibLowpThreshold(Double_t pmin) { fPthreshold = pmin; };
  
  virtual void              Init();
  virtual void              ProcessEvent();
  
  void SwitchOnGainCorrection() { fGainCorrection = kTRUE; }
  void SwitchOffGainCorrection() { fGainCorrection = kFALSE; }
  void AddBadChamber(Int_t chamberID){fBadChamberID[fNBadChambers++] = chamberID;};
  
private:
  enum{
    kNPlanes = 6,
    kNStacks = 5,
    kNChambers = 540
  };

  Bool_t GetTRDchamberID(AliESDtrack * const track, Int_t *detectors);
  void SetChamberGain();
  void ApplyGainCorrection(AliESDtrack *track);
  
  AliESDEvent           *fESD;       //! the ESD Event
  AliESDpid             *fESDpid;    //! ESD PID object

  AliTRDCalDet *fChamberGainOld;     // TRD Chamber Gain Factor used for producing the ESD
  AliTRDCalDet *fChamberGainNew;     // New TRD Chamber Gain Factor
  AliTRDCalDet *fChamberVdriftOld;   // Old drift velocity calibration
  AliTRDCalDet *fChamberVdriftNew;   // New drift velocity calibration

  Int_t fPIDmethod;                  // PID method
  Double_t fPthreshold;              // Low Momentum threshold for calibration
  Int_t fBadChamberID[kNChambers];   // List of Bad Chambers
  UInt_t fNBadChambers;              // Number of bad chambers
  Bool_t fGainCorrection;            // Apply gain correction 
  
  AliTRDTenderSupply(const AliTRDTenderSupply&c);
  AliTRDTenderSupply& operator= (const AliTRDTenderSupply&c);
  
  ClassDef(AliTRDTenderSupply, 1);  // TRD tender task
};
#endif

