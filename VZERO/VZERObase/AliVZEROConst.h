#ifndef ALIVZEROCONST_H
#define ALIVZEROCONST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */

const Float_t kIntTimeRes = 0.39; // intrinsic time resolution of the scintillator
const Float_t kIntTimeResRing[8] = { 0.27, 0.34, 0.41, 0.49, 0.26, 0.33, 0.37, 0.55 }; // intrinsic time resolution ring by ring
const Float_t kV0CDelayCables = 8.1; // delay cables on the C side (in ns)
const Float_t kV0Offset = -1338.6; // general V0 offset between the TDCs and the trigger
const Int_t   kNClocks = 21; // Number of ADC clocks that are read out
const Float_t kChargePerADC = 0.6e-12; // Charge per ADC
const Int_t   kMinTDCWidth = 7; // minimum signal width measured by TDC
const Int_t   kMaxTDCWidth = 127; // maximum signal width measured by TDC
const Float_t kPMRespTime = 6.0; // PM response time (corresponds to 1.9 ns rise time)
const Float_t kPMTransparency = 0.25; // Transparency of the first dynode of the PM
const Float_t kPMNbOfSecElec = 6.0;   // Number of secondary electrons emitted from first dynode (per ph.e.)

#endif

