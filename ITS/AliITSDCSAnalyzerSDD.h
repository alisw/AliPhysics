#ifndef ALIITSDCSANALYZERSDD_H
#define ALIITSDCSANALYZERSDD_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for SDD dcs data analysis                               //
//  called by AliITSPreprocessorSDD                              //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////



#include <TMap.h>
#include "AliITSDCSDataSDD.h"
#include "AliLog.h"

class AliITSDCSAnalyzerSDD : public TObject { 

 public:
  AliITSDCSAnalyzerSDD();
  ~AliITSDCSAnalyzerSDD();


  void AnalyzeData(TMap* dcsMap);
  void PrintDCSDPNames();
  AliITSDCSDataSDD* GetDCSData(Int_t iModule) const {return fDCSData[iModule];}

 protected:
  // Copy constructor and assignment operator not allowed.
  // They are protected to avoid misuse
  AliITSDCSAnalyzerSDD(const AliITSDCSAnalyzerSDD& /* dcsa  */);
  AliITSDCSAnalyzerSDD& operator=(const AliITSDCSAnalyzerSDD& /* dcsa */);
  void Init();

 private:

  enum {kNmodules=260, 
	kNladders3=14, 
	kNladders4=22, 
	kNmodLad3=6, 
	kNmodLad4=8};

  static const Int_t fgkNcathodes; // number of SDD cathodes
  static const Float_t fgkCathodePitch; // cathode pitch cm
  static const Int_t fgkTemperatureStatusOK; // max. Drift Field variations

  TString fHVDPNames[kNmodules];   // DCS DP names for High Voltage  
  TString fMVDPNames[kNmodules];   // DCS DP names for Medium Voltage
  TString fTLDPNames[kNmodules];   // DCS DP names for Temperature Left
  TString fTRDPNames[kNmodules];   // DCS DP names for Temperature Right
  TString fTLStDPNames[kNmodules];   // DCS DP names for status of Temperature Left
  TString fTRStDPNames[kNmodules];   // DCS DP names for status of Temperature Right
  AliITSDCSDataSDD **fDCSData;  // values of DCS data points

  ClassDef(AliITSDCSAnalyzerSDD, 1);
};

#endif
