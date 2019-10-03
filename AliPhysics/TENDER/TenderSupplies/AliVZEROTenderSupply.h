#ifndef ALIVZEROTENDERSUPPLY_H
#define ALIVZEROTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Recalculate VZERO timing and decision using the tender            //
//  (in case the LHC phase drift is updated in OCDB)                  //
//                                                                    //
////////////////////////////////////////////////////////////////////////



#include <AliTenderSupply.h>

class TF1;
class AliVZEROCalibData;
class AliVZERORecoParam;

class AliVZEROTenderSupply: public AliTenderSupply {
  
public:
  AliVZEROTenderSupply();
  AliVZEROTenderSupply(const char *name, const AliTender *tender=NULL);
  
  virtual ~AliVZEROTenderSupply(){;}

  virtual void              Init();
  virtual void              ProcessEvent();
  
  void GetPhaseCorrection();

  void SetDebug(Bool_t flag) { fDebug = flag; }

private:
  AliVZEROCalibData* fCalibData;      //! calibration data
  TF1*               fTimeSlewing;    //! Function for time slewing correction
  AliVZERORecoParam* fRecoParam;      //! pointer to reco-param object
  Float_t            fLHCClockPhase;  //! the correction to the LHC-clock phase
  Bool_t             fDebug;          //  debug on/off
  
  AliVZEROTenderSupply(const AliVZEROTenderSupply&c);
  AliVZEROTenderSupply& operator= (const AliVZEROTenderSupply&c);
  
  ClassDef(AliVZEROTenderSupply, 2)  // VZERO tender task
};


#endif

