#ifndef ALITOFSDIGITIZER_H
#define ALITOFSDIGITIZER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//__________________________________________//
//                                          //
//       Class for making SDigits in TOF    // 
//                                          //
//-- Authors: F. Pierella, A. De Caro       //
//                                          //
//__________________________________________//

/* $Id$ */

#include "TNamed.h"

class TF1;
class TString;

class AliLoader;
class AliRunLoader;

class AliTOFcalib;

class AliTOFSDigitizer: public TNamed {

public:
  AliTOFSDigitizer() ;          // ctor
  //AliTOFSDigitizer(const char* HeaderFile) ; // par ctor
  AliTOFSDigitizer(const char* HeaderFile, Int_t evNumber1=-1, Int_t nEvents=0) ; // par ctor

  AliTOFSDigitizer(const AliTOFSDigitizer &source); // copy constructor
  AliTOFSDigitizer& operator=(const AliTOFSDigitizer &/*source*/); // ass. op.

  virtual ~AliTOFSDigitizer() ; // dtor

  //static Float_t WidthTdcBin() {return fgkTdcBin;};

  virtual void  Digitize(Option_t *verboseOption); 
  void SetSDigitsFile(char * /*file*/ ) const {;}
  
  void InitParameters();
  virtual void PrintParameters() const ;
  virtual void SimulateDetectorResponse(Float_t z0, Float_t x0, Float_t geantTime, Int_t& nActivatedPads, Int_t& nFiredPads, Bool_t* isFired, Int_t* nPlace, Float_t* qInduced, Float_t* tofTime, Float_t& averageTime);
  virtual void SimulateDetectorResponseOLD(Float_t z0, Float_t x0, Float_t geantTime, Int_t& nActivatedPads, Int_t& nFiredPads, Bool_t* isFired, Int_t* nPlace, Float_t* qInduced, Float_t* tofTime, Float_t& averageTime);
  virtual void Print(Option_t* opt) const ;
  void  SetFirstEvent(Int_t event1)      {fEvent1 = event1;}
  void  SetSecondEvent(Int_t event2)     {fEvent2 = event2;}
  Int_t GetFirstEvent()  const {return fEvent1;}
  Int_t GetSecondEvent() const {return fEvent2;}
  Int_t GetNEvents() const {return (fEvent2-fEvent1);}
  void  SelectSectorAndPlate(Int_t sector, Int_t plate);

  // setters and getters for detector simulation
  // it summarizes all it is known about TOF strip 
  void  SetPadefficiency(Float_t padefficiency)      {fpadefficiency=padefficiency;}
  void  SetEdgeEffect(Int_t   edgeEffect)            {fEdgeEffect=edgeEffect;}
  void  SetEdgeTails(Int_t   edgeTails)              {fEdgeTails=edgeTails;}
  void  SetHparameter(Float_t hparameter)            {fHparameter=hparameter;}
  void  SetH2parameter(Float_t h2parameter)          {fH2parameter=h2parameter;}
  void  SetKparameter(Float_t kparameter)            {fKparameter=kparameter;}
  void  SetK2parameter(Float_t k2parameter)          {fK2parameter=k2parameter;}
  void  SetEffCenter(Float_t effCenter)              {fEffCenter=effCenter;}
  void  SetEffBoundary(Float_t effBoundary)          {fEffBoundary=effBoundary;}
  void  SetEff2Boundary(Float_t eff2Boundary)        {fEff2Boundary=eff2Boundary;}
  void  SetEff3Boundary(Float_t eff3Boundary)        {fEff3Boundary=eff3Boundary;}
  void  SetAddTRes(Float_t addTRes)                  {fAddTRes=addTRes;}
  void  SetResCenter (Float_t resCenter)             {fResCenter=resCenter;}
  void  SetResBoundary(Float_t resBoundary)          {fResBoundary=resBoundary;}
  void  SetResSlope(Float_t resSlope)                {fResSlope=resSlope;}
  void  SetTimeWalkCenter(Float_t timeWalkCenter)    {fTimeWalkCenter=timeWalkCenter;}
  void  SetTimeWalkBoundary(Float_t timeWalkBoundary){fTimeWalkBoundary=timeWalkBoundary;}
  void  SetTimeWalkSlope(Float_t timeWalkSlope)      {fTimeWalkSlope=timeWalkSlope;}

  void  SetTimeDelayFlag(Int_t timeDelayFlag)        {fTimeDelayFlag=timeDelayFlag;}
  void  SetPulseHeightSlope(Float_t pulseHeightSlope){fPulseHeightSlope=pulseHeightSlope;}
  void  SetTimeDelaySlope(Float_t timeDelaySlope)    {fTimeDelaySlope=timeDelaySlope;}
  void  SetMinimumCharge(Float_t minimumCharge)      {fMinimumCharge=minimumCharge;}
  void  SetChargeSmearing(Float_t chargeSmearing)    {fChargeSmearing=chargeSmearing;}
  void  SetLogChargeSmearing(Float_t logChargeSmearing){fLogChargeSmearing=logChargeSmearing;}
  void  SetTimeSmearing(Float_t timeSmearing)        {fTimeSmearing=timeSmearing;}
  void  SetAverageTimeFlag(Int_t averageTimeFlag)    {fAverageTimeFlag=averageTimeFlag;}

  void  SetAdcBin(Float_t adcBin)                    {fAdcBin=adcBin;}
  void  SetAdcMean(Float_t adcMean)                  {fAdcMean=adcMean;}
  void  SetAdcRms(Float_t adcRms)                    {fAdcRms=adcRms;}

  void SetTimeResolution(Float_t time) {fTimeResolution=time;}

  Float_t  GetPadefficiency()    const {return fpadefficiency;}
  Int_t    GetEdgeEffect()       const {return fEdgeEffect;}
  Int_t    GetEdgeTails()        const {return fEdgeTails;}
  Float_t  GetHparameter()       const {return fHparameter;}
  Float_t  GetH2parameter()      const {return fH2parameter;}
  Float_t  GetKparameter()       const {return fKparameter;}
  Float_t  GetK2parameter()      const {return fK2parameter;}
  Float_t  GetEffCenter()        const {return fEffCenter;}
  Float_t  GetEffBoundary()      const {return fEffBoundary;}
  Float_t  GetEff2Boundary()     const {return fEff2Boundary;}
  Float_t  GetEff3Boundary()     const {return fEff3Boundary;}
  Float_t  GetAddTRes ()         const {return fAddTRes;}
  Float_t  GetResCenter ()       const {return fResCenter;}
  Float_t  GetResBoundary()      const {return fResBoundary;}
  Float_t  GetResSlope()         const {return fResSlope;}
  Float_t  GetTimeWalkCenter()   const {return fTimeWalkCenter;}
  Float_t  GetTimeWalkBoundary() const {return fTimeWalkBoundary;}
  Float_t  GetTimeWalkSlope()    const {return fTimeWalkSlope;}
  Int_t    GetTimeDelayFlag()    const {return fTimeDelayFlag;}
  Float_t  GetPulseHeightSlope() const {return fPulseHeightSlope;}
  Float_t  GetTimeDelaySlope()   const {return fTimeDelaySlope;}
  Float_t  GetMinimumCharge()    const {return fMinimumCharge;}
  Float_t  GetChargeSmearing()   const {return fChargeSmearing;}
  Float_t  GetLogChargeSmearing()const {return fLogChargeSmearing;}
  Float_t  GetTimeSmearing()     const {return fTimeSmearing;}
  Int_t    GetAverageTimeFlag()  const {return fAverageTimeFlag;}

  Float_t  GetAdcBin()           const {return fAdcBin;}
  Float_t  GetAdcMean()          const {return fAdcMean;}
  Float_t  GetAdcRms()           const {return fAdcRms;}
  
  Float_t  GetTimeResolution()  const {return fTimeResolution;}


protected:


private:
  Int_t   fEvent1;          // lower bound for events to sdigitize
  Int_t   fEvent2;          // upper bound for events to sdigitize
  TF1     *ftail;           // pointer to formula for time with tail
  TString fHeadersFile;     // input file
  AliRunLoader* fRunLoader; //! Run Loader
  AliLoader* fTOFLoader;    //! Loader

  Int_t fSelectedSector;    // sector number for sdigitization
  Int_t fSelectedPlate ;    // plate  number for sdigitization

  // detector response simulation
  // Intrisic MRPC time resolution and pad (edge effect) parameters
  Float_t fTimeResolution;   // time resolution (ps)
  Float_t fpadefficiency;   // intrinsic pad efficiency, used if fEdgeEffect==0
  Int_t   fEdgeEffect;      // edge effects option
  Int_t   fEdgeTails;       // edge tails option
  Float_t fHparameter;      // sensitive edge (to produce hits on the neighbouring pads)
                            //                 0.7 cm (old); 0.4 cm (new)
  Float_t fH2parameter;     // parameter to fit the efficiency
  Float_t fKparameter;      // sensitive edge (going ahead towards the center
                            // no delay effects are suffered) 1.0 cm (old); 0.5 cm (new)
  Float_t fK2parameter;     // parameter to fit the efficiency
  // Pad Efficiency and Resolution parameters
  Float_t fEffCenter;       // efficiency in the central region of the pad
  Float_t fEffBoundary;     // efficiency at the boundary of the pad
  Float_t fEff2Boundary;    // efficiency value at H2parameter
  Float_t fEff3Boundary;    // efficiency value at K2parameter
  Float_t fAddTRes;         // additional contribution to 
                            // the intrinsic MRPC time resolution (ps)
  Float_t fResCenter;       // resolution (ps) in the central region of the pad
  Float_t fResBoundary;     // resolution (ps)  at the boundary of the pad
  Float_t fResSlope;        // slope (ps/K) for neighbouring pad
  // Time Walk parameters
  Float_t fTimeWalkCenter;  // time walk (ps) in the central region of the pad
  Float_t fTimeWalkBoundary;// time walk (ps) at the boundary of the pad
  Float_t fTimeWalkSlope;   // slope (ps/K) for neighbouring pad
  Int_t   fTimeDelayFlag;   // flag for delay due to the PulseHeightEffect
  Float_t fPulseHeightSlope;// It determines the charge amount induced
  // due to edge effect, using the formula
  // qInduced=exp(-PulseHeightSlope*x)
  Float_t fTimeDelaySlope;  // It determines the time delay. This is the slope
  // in the T1-T2 vs log(q1/q2) plot
  // ADC-TDC correlation parameters
  Float_t fMinimumCharge;   // Minimum charge amount which could be induced
  Float_t fChargeSmearing;  // Smearing in charge in (q1/q2) vs x plot
  Float_t fLogChargeSmearing;// Smearing in log of charge ratio
  Float_t fTimeSmearing;    // Smearing in time in time vs log(q1/q2) plot
  Int_t   fAverageTimeFlag; // flag (see the setter for details)

  Float_t fAdcBin;      // charge-window for the ADC bins [pC]
  Float_t fAdcMean;     // mean value for the ADC spectrum [bins]
  Float_t fAdcRms;      // rms value for the ADC spectrum [bins]

  AliTOFcalib * fCalib; //! calibration object

  ClassDef(AliTOFSDigitizer,5)  // creates TOF SDigits

};

#endif // AliTOFSDIGITIZER_H
