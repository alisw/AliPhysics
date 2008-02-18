#ifndef ALIITSDCSANALYZERSDD_H
#define ALIITSDCSANALYZERSDD_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////
// Class for SDD dcs data analysis called by AliITSPreprocessorSDD   //
// Origin: F.Prino, Torino, prino@to.infn.it                         //
//         V.Pospisil, CTU Prague, gdermog@seznam.cz                 //
///////////////////////////////////////////////////////////////////////

#include <TMap.h>
#include <TObjArray.h>
#include <TFile.h>
#include "AliITSDCSDataSDD.h"
#include "AliITSgeomTGeo.h"

class AliITSDCSAnalyzerSDD : public TObject 
{

 public:
  AliITSDCSAnalyzerSDD();
  ~AliITSDCSAnalyzerSDD();


  void SetVoltageDelays( Int_t HVDelay, Int_t MVDelay )     { fHVDelay = HVDelay; fMVDelay = MVDelay; }
  void SetTemperatureDelays( Int_t TLDelay, Int_t TRDelay ) { fTLDelay = TLDelay; fTRDelay = TRDelay; }
  void SetStatusDelays( Int_t StTLDelay, Int_t StTRDelay, Int_t OKDelay ) 
                                                            { fStTLDelay = StTLDelay; fStTRDelay = StTRDelay; fOKDelay = OKDelay; }
                        // There is some delay between variable readout and setting up the time stamp. Delays differs
                        //  in voltage and temperature readouts. So it is necessary to substract some value from time stamps
                        //  during the data processing 

  void SetHVThreshold( Float_t percent=0.1 ) { fHVThresholdFrac = percent / 100.0; };
  void SetMVThreshold( Float_t percent=0.1 ) { fMVThresholdFrac = percent / 100.0; };
  void SetTLThreshold( Float_t percent=0.1 ) { fTLThresholdFrac = percent / 100.0; };
  void SetTRThreshold( Float_t percent=0.1 ) { fTRThresholdFrac = percent / 100.0; };
                        // It is possible to keep value arrays condensed - new value is stored
                        //  only if it crosses a treshold. This threshold is determined from
                        //  last stored value by
                        //
                        //  low edge = value * (100 - percent) %
                        //  upper edge = value * (100 + percent) %
                        //
                        // These four Setters allows to tune tresholds for each DCS variable
                        //   independently. By default all tresholds are set to 0 % .


  void AnalyzeData( TMap* dcsMap );
                        // Processes the data

  void PrintDCSDPNames( FILE *output = stdout );
                        // Prints module identifications in text mode

  AliITSDCSDataSDD* GetDCSData( Int_t iModule ) const { return fDCSData[iModule]; }
                        // Returns data for module specified by its index in range 0..259

  AliITSDCSDataSDD* GetDCSData( Int_t iLayer, Int_t iLadder, Int_t iModule ) const
                                  { return fDCSData[AliITSgeomTGeo::GetModuleIndex( iLayer, iLadder, iModule ) - 240]; }
                        // Returns data for module specified by layer[3..4], ladder[1..22] and module number[1..8]

  void Export( char *outputDCSFileName );
                        // Exports all stored AliITSDCSDataSDD type object into specified root file. Objects are named as
                        //
                        // DCSDataSDD_module<number>
                        //
                        // where <number> is in range 0..256 and it is obtained by calling
                        //
                        // AliITSgeomTGeo::GetModuleIndex( layer, ladder, moduleInLadder ) - 240

 protected:
  AliITSDCSAnalyzerSDD(const AliITSDCSAnalyzerSDD& /* dcsa  */);
  AliITSDCSAnalyzerSDD& operator=(const AliITSDCSAnalyzerSDD& /* dcsa */);
                        // Copy constructor and assignment operator not allowed.
                        // They are protected to avoid misuse

  void Init();          // Creates module text identifications

 private:

  enum { kNmodules=260,
         kNladders3=14,
         kNladders4=22,
         kNmodLad3=6,
         kNmodLad4=8 };                   // Basic SDD geometry

  TString fHVDPNames[kNmodules];          // DCS DP names for High Voltage  
  TString fMVDPNames[kNmodules];          // DCS DP names for Medium Voltage
  TString fOKDPNames[kNmodules];          // DCS DP names for Medium Voltage
  TString fTLDPNames[kNmodules];          // DCS DP names for Temperature Left
  TString fTRDPNames[kNmodules];          // DCS DP names for Temperature Right
  TString fTLStDPNames[kNmodules];        // DCS DP names for status of Temperature Left
  TString fTRStDPNames[kNmodules];        // DCS DP names for status of Temperature Right
  AliITSDCSDataSDD *fDCSData[kNmodules];  // values of DCS data points

  Int_t fHVDelay;                         // There is some delay between variable readout
  Int_t fMVDelay;                         // and setting up the time stamp. Delays differs
  Int_t fTLDelay;                         // in voltage and temperature readouts. So it is
  Int_t fTRDelay;                         // necessary to substract some value from time stamp
  Int_t fStTLDelay;                       // during the data processing. 
  Int_t fStTRDelay;                       // Here are the values of delays stored
  Int_t fOKDelay;                         // for the 7 variables.

  Float_t fHVThresholdFrac;                // Value of HV won't be stored until it won't exceeds this fraction of last stored value
  Float_t fMVThresholdFrac;                // Same for MV
  Float_t fTLThresholdFrac;                // Same for TL
  Float_t fTRThresholdFrac;                // Same for TR

  ClassDef(AliITSDCSAnalyzerSDD, 3);

}; /*class AliITSDCSAnalyzerSDD*/

#endif
