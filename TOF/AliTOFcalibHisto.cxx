/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// *
// *
// * AliTOFcalibHisto - class to handle TOF calibration histograms,
// *                    map histograms and more
// *
// *
// * autors:   Roberto Preghenella (R+)
// * concacts: preghenella@bo.infn.it
// *
// *

#include "AliTOFcalibHisto.h"
#include "AliLog.h"
#include "TH1D.h"
#include "TFile.h"
#include "AliTOFRawStream.h"
#include "AliTOFCableLengthMap.h"
#include "AliESDtrack.h"

#define SLEW_TOTMIN 10.
#define SLEW_TOTMAX 16.

ClassImp(AliTOFcalibHisto)

//__________________________________________________________________________

TFile *AliTOFcalibHisto::fgCalibHistoFile = NULL;
TFile *AliTOFcalibHisto::fgCalibParFile = NULL;
TFile *AliTOFcalibHisto::fgCalibStatFile = NULL;

//__________________________________________________________________________

TString AliTOFcalibHisto::fgCalibHistoFileName = "$ALICE_ROOT/TOF/data/AliTOFcalibHisto.root";
TString AliTOFcalibHisto::fgCalibParFileName = "$ALICE_ROOT/TOF/data/AliTOFcalibPar.root";
TString AliTOFcalibHisto::fgCalibStatFileName = "$ALICE_ROOT/TOF/data/AliTOFcalibStat.root";

//__________________________________________________________________________

const TString AliTOFcalibHisto::fgkCalibConstName[kNcalibConsts] = {
  "LHCperiod",
  "AmphenolCableDelay",
  "FlatCableDelay",
  "InterfaceCardDelay"
};

//__________________________________________________________________________

const TString AliTOFcalibHisto::fgkCalibMapName[kNcalibMaps] = {
  /* main index */
  "Index",
  /* EO indices */
  "DDL",
  "TRM", 
  "Chain", 
  "TDC", 
  "Channel", 
  /* DO indices */
  "Sector", 
  "Plate", 
  "Strip", 
  "SectorStrip", 
  "PadZ", 
  "PadX", 
  "Pad",
  "InterfaceCardIndex",
  /* calib constants */
  "DDLBCshift",
  "FlatCableLength",
  "InterfaceCardLength",
  "AmphenolCableLength"
};

//__________________________________________________________________________

const TString AliTOFcalibHisto::fgkCalibParName[kNcalibPars] = {
  "hDDLDelay",
  "hHPTDCDelay",
  "hLeftFEAchDelay",
  "hRightFEAchDelay",
  "hFEADelay",
  "hICDelay",
  "hTRMDelay",
  "hStripDelay",
  "hIndexDelay",
  "hSlewing"
};

//__________________________________________________________________________

const TString AliTOFcalibHisto::fgkCalibStatName[kNcalibStats] = {
  "hStripStat"
};

//__________________________________________________________________________

/* LHC clock period [ns] */
const Double_t AliTOFcalibHisto::fgkLHCperiod = (24.4e-3 * 1024); /* ns */

//__________________________________________________________________________

/* Amphenol cable delay [ns/cm] */
const Double_t AliTOFcalibHisto::fgkAmphenolCableDelay = 5.13e-2; /* from measurement */

//__________________________________________________________________________

/* flat cable delay [ns/cm] */
//const Double_t AliTOFcalibHisto::fgkFlatCableDelay = 5.3e-2; /* from Amphenol 132-2829 series data-sheet */
const Double_t AliTOFcalibHisto::fgkFlatCableDelay = 5.124e-2; /* from LHC08d calibration */

//__________________________________________________________________________

/* interface card delay [ns/cm] */
//const Double_t AliTOFcalibHisto::fgkInterfaceCardDelay = 6.9e-2; /* from HyperLinx simulation */
//const Double_t AliTOFcalibHisto::fgkInterfaceCardDelay = 5.7898e-2; /* from LHC08d calibration */
const Double_t AliTOFcalibHisto::fgkInterfaceCardDelay = 6.31360207815420404e-02; /* from LHC09c calibration */

//__________________________________________________________________________

/* number of readout channels (DO/EO) */
const Int_t AliTOFcalibHisto::fgkNchannels = 157248;
const Int_t AliTOFcalibHisto::fgkNchannelsEO = 172800;

//__________________________________________________________________________

/* DDL BC shifts due to TTC fibers [LHCperiod] */
const Int_t AliTOFcalibHisto::fgkDDLBCshift[72] = {
  2, 2, -1, -1,
  2, 2, 0, 0,
  2, 2, 0, 0,
  2, 2, 0, 0,
  2, 2, 0, 0,
  2, 2, 0, 0,
  2, 2, 0, 0,
  2, 2, 0, 0,
  2, 2, 0, 0,
  2, 2, 0, 0,
  2, 2, -1, -1,
  2, 2, -1, -1,
  2, 2, -2, -2,
  2, 2, -2, -2,
  2, 2, -2, -2,
  2, 2, -1, -1,
  2, 2, -1, -1,
  2, 2, -1, -1
};

//__________________________________________________________________________

/* strip flat-cable length (preliminary) [cm] */
const Double_t AliTOFcalibHisto::fgkFlatCableLength[91] = {
  18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 17.,
  21., 21., 21., 21., 21., 17., 17., 21., 21., 17., 21., 21., 21., 17., 21., 21., 17., 21., 23.,
  17., 19., 17., 19., 17., 19., 17., 19., 17., 19., 17., 19., 17., 19., 17.,
  23., 21., 17., 21., 21., 17., 21., 21., 21., 17., 21., 21., 17., 17., 21., 21., 21., 21., 21.,
  17., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18.
};

//__________________________________________________________________________

/* interface card lenght (preliminary) [cm] */
const Double_t AliTOFcalibHisto::fgkInterfaceCardLength[48] = {
  13.97, 12.57, 14.52, 13.10, 15.44, 13.60, 10.58, 9.14, 
  11.21, 9.76, 12.11, 10.76, 8.67, 7.58, 9.32, 8.09,
  10.24, 8.4, 5.51, 4.31, 6.54, 5.23, 7.48, 6.28,
  10.43, 8.76, 11.05, 9.43, 11.72, 10.14, 7.2, 5.69,
  7.71, 6.26, 8.36, 7.19, 4.85, 4.09, 5.57, 4.35, 
  6.59, 5.12, 2.49, 2.96, 2.70, 2.76, 2.91, 2.55
};

//__________________________________________________________________________

Bool_t AliTOFcalibHisto::fgCableCorrectionFlag[kNcorrections] = {
  kFALSE, // kDDLBCcorr
  kTRUE, // kAmphenolCableCorr
  kTRUE, // kFlatCableCorr
  kTRUE, // kInterfaceCardCorr
  kFALSE, // kDDLdelayCorr
  kFALSE, // kHPTDCdelayCorr
  kFALSE, // kFEAchDelayCorr
  kFALSE, // kFEAdelayCorr
  kFALSE, // kTRMdelayCorr
  kFALSE, // kICdelayCorr
  kFALSE, // kStripDelayCorr
  kFALSE, // kIndexDelayCorr
  kFALSE, // kTimeSlewingCorr
};

//__________________________________________________________________________

Bool_t AliTOFcalibHisto::fgFullCorrectionFlag[kNcorrections] = {
  kFALSE, // kDDLBCcorr
  kTRUE, // kAmphenolCableCorr
  kTRUE, // kFlatCableCorr
  kTRUE, // kInterfaceCardCorr
  kTRUE, // kDDLdelayCorr
  kTRUE, // kHPTDCdelayCorr
  kTRUE, // kFEAchDelayCorr
  kTRUE, // kFEAdelayCorr
  kTRUE, // kTRMdelayCorr
  kFALSE, // kICdelayCorr
  kTRUE, // kStripDelayCorr
  kTRUE, // kIndexDelayCorr
  kTRUE, // kTimeSlewingCorr
};

//__________________________________________________________________________

const Int_t AliTOFcalibHisto::fgkStripStat[18][91] = {

  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S00 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S01 */

  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S02 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S03 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S04 */

  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S05 */

  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S06 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S07 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S08 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S09 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S10 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S11 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S12 */

  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S13 */

  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S14 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S15 */

  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, /* S16 */
  
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}  /* S17 */

};

//__________________________________________________________________________

AliTOFcalibHisto::AliTOFcalibHisto() :
  TObject(),
  fCalibConst(),
  fCalibMap(),
  fCalibPar(),
  fCalibStat()
{
  /* default constructor */
}

//__________________________________________________________________________

AliTOFcalibHisto::~AliTOFcalibHisto()
{
  /* default destructor */
}

//__________________________________________________________________________

void 
AliTOFcalibHisto::LoadHisto(TFile * const file, TH1D **histo, const Char_t *name) 
{
  /* load histo */
  *histo = (TH1D *)file->Get(name);
  if (!*histo)
    AliWarning(Form("error while getting %s histo", name));
}

//__________________________________________________________________________

void 
AliTOFcalibHisto::CreateHisto(TH1D **histo, const Char_t *name, Int_t size) 
{
  /* create histo */
  *histo = new TH1D(name, Form(";index;%s", name), size, 0, size);
  if (!*histo)
    AliWarning(Form("error while creating %s histo", name));
}

//__________________________________________________________________________

void 
AliTOFcalibHisto::WriteHisto(TFile *file, TH1D *histo) 
{
  /* write histo */
  if (!file || !file->IsOpen() || !histo)
    return;
  file->cd(); 
  histo->Write();
}

//__________________________________________________________________________

void
AliTOFcalibHisto::SetHisto(TH1D *histo, Int_t index, Double_t value)
{
  /* set histo */
  if (!histo)
    return;
  histo->SetBinContent(index + 1, value);
}

//__________________________________________________________________________

Double_t
AliTOFcalibHisto::GetHisto(TH1D *histo, Int_t index)
{
  /* get histo */
  if (!histo) {
    AliWarning("cannot get histo");
    return 0.;
  }
  return histo->GetBinContent(index + 1);
}

//__________________________________________________________________________

void
AliTOFcalibHisto::LoadCalibHisto()
{
  /* load calib histo */

  if (fgCalibHistoFile && fgCalibHistoFile->IsOpen())
    AliWarning("calib histo file already open: reloading"); 

  /* open input file */
  TFile *fileIn = TFile::Open(GetCalibHistoFileName());
  if (!fileIn || !fileIn->IsOpen()) {
    AliFatal(Form("cannot open input file %s", GetCalibHistoFileName()));
    return;
  }

  /* set calib histo file */
  fgCalibHistoFile = fileIn;

  /* load consts */
  for (Int_t iConst = 0; iConst < kNcalibConsts; iConst++)
    LoadHisto(fileIn, &fCalibConst[iConst], fgkCalibConstName[iConst].Data());
  /* load maps */
  for (Int_t iMap = 0; iMap < kNcalibMaps; iMap++)
    LoadHisto(fileIn, &fCalibMap[iMap], fgkCalibMapName[iMap].Data());
}

//__________________________________________________________________________

void
AliTOFcalibHisto::LoadCalibPar()
{
  /* load calib par */

  if (fgCalibParFile && fgCalibParFile->IsOpen())
    AliWarning("calib par file already open: reloading"); 

  /* load calib histo */
  LoadCalibHisto();

  /* open input file */
  TFile *fileIn = TFile::Open(GetCalibParFileName());
  if (!fileIn || !fileIn->IsOpen()) {
    AliError(Form("cannot open input file %s", GetCalibParFileName()));
    return;
  }

  /* set calib par file */
  fgCalibParFile = fileIn;

  /* load pars */
  for (Int_t i = 0; i < kNcalibPars; i++)
    LoadHisto(fileIn, &fCalibPar[i], fgkCalibParName[i].Data());

}

//__________________________________________________________________________

void
AliTOFcalibHisto::LoadCalibStat()
{
  /* load calib stat */

  if (fgCalibStatFile && fgCalibStatFile->IsOpen())
    AliWarning("calib par file already open: reloading"); 

  /* load calib histo */
  LoadCalibHisto();

  /* open input file */
  TFile *fileIn = TFile::Open(GetCalibStatFileName());
  if (!fileIn || !fileIn->IsOpen())
    AliError(Form("cannot open input file %s", GetCalibStatFileName()));

  /* set calib par file */
  fgCalibStatFile = fileIn;

  /* load pars */
  for (Int_t i = 0; i < kNcalibStats; i++)
    LoadHisto(fileIn, &fCalibStat[i], fgkCalibStatName[i].Data());

}

//__________________________________________________________________________

void
AliTOFcalibHisto::WriteCalibHisto()
{
  /* write calib histo */

  /* open output file */
  TFile *fileOut = TFile::Open(GetCalibHistoFileName(), "RECREATE");
  if (!fileOut || !fileOut->IsOpen()) {
    AliFatal(Form("cannot open output file %s", GetCalibHistoFileName()));
    return;
  }

  /* create consts */
  for (Int_t iConst = 0; iConst < kNcalibConsts; iConst++)
    CreateHisto(&fCalibConst[iConst], fgkCalibConstName[iConst].Data(), 1);
  /* create maps */
  for (Int_t iMap = 0; iMap < kNcalibMaps; iMap++)
    if (iMap == kIndex)
      CreateHisto(&fCalibMap[iMap], fgkCalibMapName[iMap].Data(), fgkNchannelsEO);
    else
      CreateHisto(&fCalibMap[iMap], fgkCalibMapName[iMap].Data(), fgkNchannels);

  /*** SETUP CONSTS ***/

  SetHisto(fCalibConst[kLHCperiod], 0, fgkLHCperiod);
  SetHisto(fCalibConst[kAmphenolCableDelay], 0, fgkAmphenolCableDelay);
  SetHisto(fCalibConst[kFlatCableDelay], 0, fgkFlatCableDelay);
  SetHisto(fCalibConst[kInterfaceCardDelay], 0, fgkInterfaceCardDelay);
  
  /***  SETUP MAPS  ***/

  AliTOFRawStream rawStream;
  Int_t indexEO, det[5], dummy, index, sector, plate, strip, sectorStrip, padz, padx, pad, icIndex;

  /* temporarly disable warnings */
  AliLog::EType_t logLevel = (AliLog::EType_t)AliLog::GetGlobalLogLevel();
  AliLog::SetGlobalLogLevel(AliLog::kError);

  /* loop over electronics oriented (EO) indices */
  for (Int_t ddl = 0; ddl < 72; ddl++)
    for (Int_t trm = 0; trm < 10; trm++)
      for (Int_t chain = 0; chain < 2; chain++)
	for (Int_t tdc = 0; tdc < 15; tdc++)
	  for (Int_t channel = 0; channel < 8; channel++) {
	    
	    /* compute index EO */
	    indexEO = GetIndexEO(ddl, trm, chain, tdc, channel);

	    /* convert EO indices into detector oriented (DO) indices
	       (this call causes some warnings because the loop includes
	       EO indices which are not connected to physical channels) */
	    rawStream.EquipmentId2VolumeId(ddl, trm + 3, chain, tdc, channel, det);
	    
	    /* swap det[3] and det[4] */
	    dummy = det[3]; det[3] = det[4]; det[4] = dummy;
	    
	    /* check detector indices */
	    if (det[0] < 0 || det[0] > 17 ||
		det[1] < 0 || det[1] > 4 ||
		det[2] < 0 || det[2] > 18 ||
		det[3] < 0 || det[3] > 1 ||
		det[4] < 0 || det[4] > 47) {
	      SetHisto(fCalibMap[kIndex], indexEO, -1);
	      continue;
	    }
	    
	    /* setup information */
	    index = AliTOFGeometry::GetIndex(det);
	    sector = det[0];
	    plate = det[1];
	    strip = det[2];
	    sectorStrip = plate < 3 ? plate * 19 + strip : plate * 19 - 4 + strip;
	    padz = det[3];
	    padx = det[4];
	    pad = padz + 2 * padx;
	    icIndex = pad < 48 ? pad : 95 - pad;

	    /* set maps */

	    /* main index */
	    SetHisto(fCalibMap[kIndex], indexEO, index);
	    /* EO indices */
	    SetHisto(fCalibMap[kDDL], index, ddl);
	    SetHisto(fCalibMap[kTRM], index, trm);
	    SetHisto(fCalibMap[kChain], index, chain);
	    SetHisto(fCalibMap[kTDC], index, tdc);
	    SetHisto(fCalibMap[kChannel], index, channel);
	    /* DO indices */
	    SetHisto(fCalibMap[kSector], index, sector);
	    SetHisto(fCalibMap[kPlate], index, plate);
	    SetHisto(fCalibMap[kStrip], index, strip);
	    SetHisto(fCalibMap[kSectorStrip], index, sectorStrip);
	    SetHisto(fCalibMap[kPadZ], index, padz);
	    SetHisto(fCalibMap[kPadX], index, padx);
	    SetHisto(fCalibMap[kPad], index, pad);
	    SetHisto(fCalibMap[kInterfaceCardIndex], index, icIndex);
	    /* calib constants */
	    SetHisto(fCalibMap[kDDLBCshift], index, fgkDDLBCshift[ddl]);
	    SetHisto(fCalibMap[kFlatCableLength], index, fgkFlatCableLength[sectorStrip]);
	    SetHisto(fCalibMap[kInterfaceCardLength], index, fgkInterfaceCardLength[icIndex]);
	    SetHisto(fCalibMap[kAmphenolCableLength], index, AliTOFCableLengthMap::GetCableLength(ddl, trm + 3, chain, tdc));
	    
	  } /* loop over electronics oriented (EO) indices */

  /* re-enable warnings */
  AliLog::SetGlobalLogLevel(logLevel);

  /* write consts */
  for (Int_t iConst = 0; iConst < kNcalibConsts; iConst++)
    WriteHisto(fileOut, fCalibConst[iConst]);
  /* write maps */
  for (Int_t iMap = 0; iMap < kNcalibMaps; iMap++)
    WriteHisto(fileOut, fCalibMap[iMap]);

  /* close output file */
  fileOut->Close();
}

//__________________________________________________________________________

void
AliTOFcalibHisto::WriteCalibStat()
{
  /* write calib stat */

  /* open output file */
  TFile *fileOut = TFile::Open(GetCalibStatFileName(), "RECREATE");
  if (!fileOut || !fileOut->IsOpen()) {
    AliFatal(Form("cannot open output file %s", GetCalibStatFileName()));
    return;
  }

  /* create stats */
  for (Int_t iStat = 0; iStat < kNcalibStats; iStat++)
    CreateHisto(&fCalibStat[iStat], fgkCalibStatName[iStat].Data(), fgkNchannels);

  /***  SETUP STATS  ***/

  Int_t sector, sectorStrip;

  /* load calib histo */
  LoadCalibHisto();

  /* loop over channels */
  for (Int_t index = 0; index < fgkNchannels; index++) {
    sector = (Int_t)GetCalibMap(kSector, index);
    sectorStrip = (Int_t)GetCalibMap(kSectorStrip, index);
    /* strip stat */
    SetHisto(fCalibStat[kStripStat], index, (Double_t)fgkStripStat[sector][sectorStrip]);
  } /* loop over channels */

  /* write stats */
  for (Int_t iStat = 0; iStat < kNcalibStats; iStat++)
    WriteHisto(fileOut, fCalibStat[iStat]);

  /* close output file */
  fileOut->Close();
}

//__________________________________________________________________________

Double_t
AliTOFcalibHisto::GetCorrection(Int_t corr, Int_t index, Double_t tot)
{
  /* apply correction */

  Int_t ddl, chain, tdc, channel, hptdc, pbCh, feaIndex, sector, plate, padx, trm, icIndex, sectorStrip;
  Double_t slewing;
  
  switch (corr) {
  case kDDLBCcorr:
    return -GetCalibConst(kLHCperiod) * GetCalibMap(kDDLBCshift, index);
  case kAmphenolCableCorr:
    return GetCalibConst(kAmphenolCableDelay) * GetCalibMap(kAmphenolCableLength, index);
  case kFlatCableCorr:
    return GetCalibConst(kFlatCableDelay) * GetCalibMap(kFlatCableLength, index);
  case kInterfaceCardCorr:
    return GetCalibConst(kInterfaceCardDelay) * GetCalibMap(kInterfaceCardLength, index);
  case kDDLdelayCorr:
    ddl = (Int_t)GetCalibMap(kDDL, index);
    return GetCalibPar(kDDLdelayPar, ddl);
  case kHPTDCdelayCorr:
    chain = (Int_t)GetCalibMap(kChain, index);
    tdc = (Int_t)GetCalibMap(kTDC, index);
    hptdc = tdc + 15 * chain;
    return GetCalibPar(kHPTDCdelayPar, hptdc);
  case kFEAchDelayCorr:
    ddl = (Int_t)GetCalibMap(kDDL, index);
    tdc = (Int_t)GetCalibMap(kTDC, index);
    channel = (Int_t)GetCalibMap(kChannel, index);
    pbCh = channel + 8 * (tdc % 3);
    if (ddl % 2 == 0)
      return GetCalibPar(kRightFEAchDelayPar, pbCh);
    else
      return GetCalibPar(kLeftFEAchDelayPar, pbCh);
  case kFEAdelayCorr:
    sector = (Int_t)GetCalibMap(kSector, index);
    plate = (Int_t)GetCalibMap(kPlate, index);
    sectorStrip = (Int_t)GetCalibMap(kSectorStrip, index);
    padx = (Int_t)GetCalibMap(kPadX, index);
    feaIndex = padx / 12 + 4 * sectorStrip + 364 * sector;      
    return GetCalibPar(kFEAdelayPar, feaIndex);
  case kTRMdelayCorr:
    ddl = (Int_t)GetCalibMap(kDDL, index);
    trm = (Int_t)GetCalibMap(kTRM, index);
    return GetCalibPar(kTRMdelayPar, trm + 10 * ddl);
  case kICdelayCorr:
    icIndex = (Int_t)GetCalibMap(kInterfaceCardIndex, index);
    return GetCalibPar(kICdelayPar, icIndex);
  case kStripDelayCorr:
    sector = (Int_t)GetCalibMap(kSector, index);
    sectorStrip = (Int_t)GetCalibMap(kSectorStrip, index);
    return GetCalibPar(kStripDelayPar, sectorStrip + 91 * sector);
  case kIndexDelayCorr:
    return GetCalibPar(kIndexDelayPar, index);
  case kTimeSlewingCorr:
    tot = tot < SLEW_TOTMIN ? SLEW_TOTMIN : tot;
    tot = tot > SLEW_TOTMAX ? SLEW_TOTMAX : tot;
    slewing = 0.;
    for (Int_t i = 0; i < fCalibPar[kTimeSlewingPar]->GetNbinsX(); i++)
      slewing += GetCalibPar(kTimeSlewingPar, i) * TMath::Power(tot, i);
    return slewing;
  default:
    AliWarning(Form("unknown correction flag (%d)", corr));
    return 0.;
  }
}

//__________________________________________________________________________

Double_t
AliTOFcalibHisto::GetNominalCorrection(Int_t index)
{
  /* get nominal correction */
  Double_t corr = 0;
  for (Int_t iCorr = 0; iCorr < kNcorrections; iCorr++)
    corr += GetCorrection(iCorr, index);
  return corr;
}

//__________________________________________________________________________

void
AliTOFcalibHisto::ApplyNominalCorrection(AliESDtrack *track)
{
  /* apply nominal correction */
  
  Double_t rawTime = track->GetTOFsignalRaw();
  Int_t index = track->GetTOFCalChannel();
  Double_t time = rawTime - 1.e3 * GetNominalCorrection(index);
  track->SetTOFsignal(time);
}

//__________________________________________________________________________

Double_t
AliTOFcalibHisto::GetCableCorrection(Int_t index)
{
  /* get cable correction */
  Double_t corr = 0;
  for (Int_t iCorr = 0; iCorr < kNcorrections; iCorr++)
    if (fgCableCorrectionFlag[iCorr])
      corr += GetCorrection(iCorr, index);
  return corr;
}

//__________________________________________________________________________

Double_t
AliTOFcalibHisto::GetFullCorrection(Int_t index, Double_t tot)
{
  /* get full correction */
  Double_t corr = 0;
  for (Int_t iCorr = 0; iCorr < kNcorrections; iCorr++)
    if (fgFullCorrectionFlag[iCorr]) {
      corr += GetCorrection(iCorr, index, tot);
    }
  return corr;
}

//__________________________________________________________________________

Bool_t
AliTOFcalibHisto::GetStatus(Int_t index)
{
  /* get status */

  Bool_t status = kTRUE;

  for (Int_t iStat = 0; iStat < kNcalibStats; iStat++)
    status &= GetCalibStat(iStat, index);

  return status;
}
