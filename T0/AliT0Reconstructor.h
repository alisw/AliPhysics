#ifndef ALIT0RECONSTRUCTOR_H
#define ALIT0RECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/******************************************************************** 
 * header class T0 reconstruction 
 * Alla Maevskaya INR RAS alla@inr.ru      *
 * Alla.Maevskaya@cern.ch
 *******************************************************************/

#include "AliReconstructor.h"
#include "AliT0Parameters.h"
#include "AliT0Calibrator.h"
#include "AliT0RecoParam.h"
#include "AliESDTZEROfriend.h"
#include "AliESDTZERO.h"

class AliT0Reconstructor: public AliReconstructor {
 public:
  AliT0Reconstructor();
  virtual ~AliT0Reconstructor() {};

  virtual  void   Reconstruct(TTree* fdigits, TTree * frecpoints) const;
  virtual  void   Reconstruct(AliRawReader*rawReader , TTree* recTree) const;
  
  virtual void     FillESD( AliRawReader*/*rawReader*/,  TTree*clustersTree, AliESDEvent*esd ) const
  {FillESD((TTree*)NULL,clustersTree,esd);}
  virtual void     FillESD( TTree* digitsTree,  TTree*clustersTree, AliESDEvent*esd ) const;

  virtual Bool_t   HasDigitConversion() const {return kFALSE;}
  static const AliT0RecoParam* GetRecoParam()
    { return dynamic_cast<const AliT0RecoParam*>(AliReconstructor::GetRecoParam(11)); } // getting RecoParam obj

  //!!!!!!!!!!!!!!!!!!!!!!!!!!
  Bool_t  PileupFlag() const;
  Bool_t  BackgroundFlag() const;
  Bool_t  SatelliteFlag() const;
  //!!!!!!!!!!!!!!!!!!!!!!!
 
 protected:
  Float_t             fdZonA;             // Zideal - Zreal side A 
  Float_t             fdZonC;             // Zideal - Zreal side C
  Float_t             fZposition;        // vertex position
  Float_t             fTime0vertex[24];  // time position if Zvertex=0
  AliT0Parameters     *fParam;           //pointer to T0 parameters class     
  TObjArray           fAmpLEDrec;        // amp LED-CFD 
  TObjArray           fQTC;        // QTC vs #MIPs
  TObjArray           fAmpLED;        // LED-CFD vs #MIPs
  AliT0Calibrator     *fCalib;           //pointer to T0 Calibrator     
  Float_t fLatencyHPTDC;  //latency HPTDC
  Float_t fLatencyL1;     //  latency for (T0A+T0C)/2
  Float_t fLatencyL1A;    // latency for T0A
  Float_t fLatencyL1C;    //latency for T0C
  Float_t fGRPdelays;    //latency for T0C
  Float_t *fTimeMeanShift;
  Float_t *fTimeSigmaShift;

  AliESDTZEROfriend*  fESDTZEROfriend; // ESD friend object 
  AliESDTZERO*        fESDTZERO;       // ESD output object  
 
 private:
  AliT0Reconstructor( const AliT0Reconstructor&r ); //Not implemented
  AliT0Reconstructor& operator=(const AliT0Reconstructor&r); //Not implemented

  ClassDef(AliT0Reconstructor, 8)   // class for the T0 reconstruction

};

typedef AliT0Reconstructor AliSTARTReconstructor; // for backward compatibility

#endif
