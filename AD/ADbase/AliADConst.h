#ifndef ALIADCONST_H
#define ALIADCONST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */

#include "TROOT.h"

const Float_t kADIntTimeRes = 0.39; // intrinsic time resolution of the scintillator
const Float_t kADOffset = 962; // general AD offset between the TDCs and the trigger
const Int_t   kADNClocks = 21; // Number of ADC clocks that are read out
const Float_t kADChargePerADC = 0.6e-12; // Charge per ADC
const Int_t   kNPhotonsPerMIP = 49000;// Number of photons per MIP
const Int_t   kADMinTDCWidth = 13; // minimum signal width measured by TDC
const Int_t   kADMaxTDCWidth = 128; // maximum signal width measured by TDC
const Float_t kADPMRespTime = 6.0; // PM response time (corresponds to 1.9 ns rise time)
const Float_t kADPMTransparency = 0.25; // Transparency of the first dynode of the PM
const Float_t kADPMNbOfSecElec = 6.0;   // Number of secondary electrons emitted from first dynode (per ph.e.)
const Float_t kPhotoCathodeEfficiency = 0.18; // Photocathode efficiency
const Int_t   kNCIUBoards = 2; //Number of CIU boards
/*				    |------------Cside------------|----------Aside-------|   */	
const Int_t   kOfflineChannel[16] = {15, 11, 14, 10, 13, 9, 12, 8, 7, 3, 6, 2, 5, 1, 4, 0};
/*	      Online->Offline								     */
const Int_t   kOnlineChannel[16] =  {15, 13, 11, 9, 14, 12, 10, 8, 7, 5, 3, 1, 6, 4, 2, 0};
/*	      Offline->Online								     */

#endif

