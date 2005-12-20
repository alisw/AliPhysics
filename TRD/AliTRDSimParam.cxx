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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class containing constant simulation parameters                           //
//                                                                           //
// Request an instance with AliTRDSimParam::Instance()                 //
// Then request the needed values                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>

#include "AliRun.h"

#include "AliTRDSimParam.h"

ClassImp(AliTRDSimParam)

AliTRDSimParam* AliTRDSimParam::fgInstance = 0;
Bool_t AliTRDSimParam::fgTerminated = kFALSE;

//_ singleton implementation __________________________________________________
AliTRDSimParam* AliTRDSimParam::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  // 
  
  if (fgTerminated != kFALSE)
    return 0;

  if (fgInstance == 0)
    fgInstance = new AliTRDSimParam();
  
  return fgInstance;
}

void AliTRDSimParam::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag, instances cannot be requested anymore
  // This function can be called several times.
  //
  
  fgTerminated = kTRUE;
  
  if (fgInstance != 0)
  {
    delete fgInstance;
    fgInstance = 0;
  }
}

//_____________________________________________________________________________
AliTRDSimParam::AliTRDSimParam()
{
  //
  // default constructor
  //
  
  fGasGain            = 0.0;
  fNoise              = 0.0;
  fChipGain           = 0.0;
  
  fADCoutRange        = 0.0;
  fADCinRange         = 0.0;
  fADCthreshold       = 0;
  fADCbaseline        = 0;        
  
  fDiffusionOn        = kFALSE;
  
  fElAttachOn         = kFALSE;
  fElAttachProp       = 0.0;
  
  fTRFOn              = kFALSE;
  fTRFsmp             = 0;
  fTRFbin             = 0;
  fTRFlo              = 0.0;
  fTRFhi              = 0.0;
  fTRFwid             = 0.0;
  
  fCTOn               = kFALSE;
  fCTsmp              = 0;
  
  fTCOn               = kFALSE;
  
  fTCnexp             = 0;
  
  fAnodeWireOffset    = 0.0;
  fPadCoupling        = 0.0;
  fTimeCoupling       = 0.0;
  fTimeStructOn       = kFALSE;
  
  Init();
}

//_____________________________________________________________________________
void AliTRDSimParam::Init()
{
  // 
  // default initializiation
  //
  
  // The default parameter for the digitization
  fGasGain        = 4000.;
  fChipGain       = 12.4;
  fNoise          = 1000.;
  fADCoutRange    = 1023.;          // 10-bit ADC
  fADCinRange     = 2000.;          // 2V input range
  fADCthreshold   = 1;
  fADCbaseline    = 0;

  // Diffusion on
  fDiffusionOn    = kTRUE;
  
  // Propability for electron attachment
  fElAttachOn     = kFALSE;
  fElAttachProp   = 0.0;

  // The time response function
  fTRFOn          = kTRUE;

  // The cross talk
  fCTOn           = kTRUE;

  // The tail cancelation
  fTCOn           = kTRUE;
  
  // The number of exponentials
  fTCnexp         = 1;

  // The pad coupling factor
  //fPadCoupling    = 0.3;
  // Use 0.46 instead which reproduces better the test beam
  // data, even tough it is not understood why.
  fPadCoupling    = 0.46;

  // The time coupling factor (same number as for the TPC)
  fTimeCoupling   = 0.4;

  // Distance of first Anode wire from first pad edge
  fAnodeWireOffset = 0.25;

  // Use drift time maps
  fTimeStructOn = kTRUE;
  
  ReInit();
}

//_____________________________________________________________________________
AliTRDSimParam::~AliTRDSimParam() 
{
  //
  // destructor
  //
  
  if (fTRFsmp) {
    delete [] fTRFsmp;
    fTRFsmp = 0;
  }

  if (fCTsmp) {
    delete [] fCTsmp;
    fCTsmp  = 0;
  }
}

//_____________________________________________________________________________
AliTRDSimParam::AliTRDSimParam(const AliTRDSimParam &p):TObject(p)
{
  //
  // copy constructor
  //

  ((AliTRDSimParam &) p).Copy(*this);
}


//_____________________________________________________________________________
AliTRDSimParam &AliTRDSimParam::operator=(const AliTRDSimParam &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDSimParam &) p).Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliTRDSimParam::Copy(TObject &p) const
{
  //
  // Copy function
  //
  
  AliTRDSimParam* target = dynamic_cast<AliTRDSimParam*> (&p);
  if (!target)
    return;

  target->fGasGain            = fGasGain;
  //target->fField              = fField;
  target->fNoise              = fNoise;
  target->fChipGain           = fChipGain;
  
  target->fADCoutRange        = fADCoutRange;
  target->fADCinRange         = fADCinRange;
  target->fADCthreshold       = fADCthreshold;
  target->fADCbaseline        = fADCbaseline; 
  
  target->fDiffusionOn        = fDiffusionOn; 
  
  target->fElAttachOn         = fElAttachOn;
  target->fElAttachProp       = fElAttachProp;
  
  target->fTRFOn              = fTRFOn;
  if (target->fTRFsmp) 
    delete[] target->fTRFsmp;
  target->fTRFsmp = new Float_t[fTRFbin];
  for (Int_t iBin = 0; iBin < fTRFbin; iBin++) {
    target->fTRFsmp[iBin] = fTRFsmp[iBin];
  }
  target->fTRFbin             = fTRFbin;
  target->fTRFlo              = fTRFlo;
  target->fTRFhi              = fTRFhi;
  target->fTRFwid             = fTRFwid;
  
  target->fCTOn               = fCTOn;
  if (target->fCTsmp) 
    delete[] target->fCTsmp;
  target->fCTsmp  = new Float_t[fTRFbin];
  for (Int_t iBin = 0; iBin < fTRFbin; iBin++) {
    target->fCTsmp[iBin]  = fCTsmp[iBin];
  }
  
  target->fTCOn               = fTCOn;
  target->fTCnexp             = fTCnexp;

  target->fAnodeWireOffset    = fAnodeWireOffset;
  target->fPadCoupling        = fPadCoupling;
  target->fTimeCoupling       = fTimeCoupling;
}

//_____________________________________________________________________________
void AliTRDSimParam::ReInit()
{
  //
  // Reinitializes the parameter class after a change
  //

  // The range and the binwidth for the sampled TRF 
  fTRFbin = 100;
  // Start 0.2 mus before the signal
  fTRFlo  = -0.2;
  // End the maximum drift time after the signal 
  fTRFhi  = 2.2;
  // 
  fTRFwid = (fTRFhi - fTRFlo) / ((Float_t) fTRFbin);

  // Create the sampled TRF
  SampleTRF();
}

//_____________________________________________________________________________
void AliTRDSimParam::SampleTRF()
{
  //
  // Samples the time response function
  //
  // New TRF from Venelin Angelov, simulated with CADENCE
  // Pad-ground capacitance = 25 pF
  // Pad-pad cross talk capacitance = 6 pF   
  //

  Int_t   ipos1;
  Int_t   ipos2;
  Float_t diff;

  const Int_t kNpasa     = 252;

  Float_t time[kNpasa]   = { -0.220000, -0.210000, -0.200000, -0.190000 
                           , -0.180000, -0.170000, -0.160000, -0.150000 
                           , -0.140000, -0.130000, -0.120000, -0.110000 
                           , -0.100000, -0.090000, -0.080000, -0.070000 
                           , -0.060000, -0.050000, -0.040000, -0.030000 
                           , -0.020000, -0.010000, -0.000000,  0.010000 
                           ,  0.020000,  0.030000,  0.040000,  0.050000 
                           ,  0.060000,  0.070000,  0.080000,  0.090000 
                           ,  0.100000,  0.110000,  0.120000,  0.130000 
                           ,  0.140000,  0.150000,  0.160000,  0.170000 
                           ,  0.180000,  0.190000,  0.200000,  0.210000 
                           ,  0.220000,  0.230000,  0.240000,  0.250000 
                           ,  0.260000,  0.270000,  0.280000,  0.290000 
                           ,  0.300000,  0.310000,  0.320000,  0.330000 
                           ,  0.340000,  0.350000,  0.360000,  0.370000 
                           ,  0.380000,  0.390000,  0.400000,  0.410000 
                           ,  0.420000,  0.430000,  0.440000,  0.450000 
                           ,  0.460000,  0.470000,  0.480000,  0.490000 
                           ,  0.500000,  0.510000,  0.520000,  0.530000 
                           ,  0.540000,  0.550000,  0.560000,  0.570000 
                           ,  0.580000,  0.590000,  0.600000,  0.610000 
                           ,  0.620000,  0.630000,  0.640000,  0.650000 
                           ,  0.660000,  0.670000,  0.680000,  0.690000 
                           ,  0.700000,  0.710000,  0.720000,  0.730000 
                           ,  0.740000,  0.750000,  0.760000,  0.770000 
                           ,  0.780000,  0.790000,  0.800000,  0.810000 
                           ,  0.820000,  0.830000,  0.840000,  0.850000 
                           ,  0.860000,  0.870000,  0.880000,  0.890000 
                           ,  0.900000,  0.910000,  0.920000,  0.930000 
                           ,  0.940000,  0.950000,  0.960000,  0.970000 
                           ,  0.980000,  0.990000,  1.000000,  1.010000 
                           ,  1.020000,  1.030000,  1.040000,  1.050000 
                           ,  1.060000,  1.070000,  1.080000,  1.090000 
                           ,  1.100000,  1.110000,  1.120000,  1.130000 
                           ,  1.140000,  1.150000,  1.160000,  1.170000 
                           ,  1.180000,  1.190000,  1.200000,  1.210000 
                           ,  1.220000,  1.230000,  1.240000,  1.250000 
                           ,  1.260000,  1.270000,  1.280000,  1.290000 
                           ,  1.300000,  1.310000,  1.320000,  1.330000 
                           ,  1.340000,  1.350000,  1.360000,  1.370000 
                           ,  1.380000,  1.390000,  1.400000,  1.410000 
                           ,  1.420000,  1.430000,  1.440000,  1.450000 
                           ,  1.460000,  1.470000,  1.480000,  1.490000 
                           ,  1.500000,  1.510000,  1.520000,  1.530000 
                           ,  1.540000,  1.550000,  1.560000,  1.570000 
                           ,  1.580000,  1.590000,  1.600000,  1.610000 
                           ,  1.620000,  1.630000,  1.640000,  1.650000 
                           ,  1.660000,  1.670000,  1.680000,  1.690000 
                           ,  1.700000,  1.710000,  1.720000,  1.730000 
                           ,  1.740000,  1.750000,  1.760000,  1.770000 
                           ,  1.780000,  1.790000,  1.800000,  1.810000 
                           ,  1.820000,  1.830000,  1.840000,  1.850000 
                           ,  1.860000,  1.870000,  1.880000,  1.890000 
                           ,  1.900000,  1.910000,  1.920000,  1.930000 
                           ,  1.940000,  1.950000,  1.960000,  1.970000 
                           ,  1.980000,  1.990000,  2.000000,  2.010000 
                           ,  2.020000,  2.030000,  2.040000,  2.050000 
                           ,  2.060000,  2.070000,  2.080000,  2.090000 
                           ,  2.100000,  2.110000,  2.120000,  2.130000 
                           ,  2.140000,  2.150000,  2.160000,  2.170000 
                           ,  2.180000,  2.190000,  2.200000,  2.210000 
                           ,  2.220000,  2.230000,  2.240000,  2.250000 
                           ,  2.260000,  2.270000,  2.280000,  2.290000 };

  Float_t signal[kNpasa] = {  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000396 
                           ,  0.005096,  0.022877,  0.061891,  0.126614 
                           ,  0.215798,  0.324406,  0.444507,  0.566817 
                           ,  0.683465,  0.787089,  0.873159,  0.937146 
                           ,  0.979049,  0.999434,  1.000000,  0.983579 
                           ,  0.954134,  0.913364,  0.866365,  0.813703 
                           ,  0.759910,  0.706116,  0.653454,  0.603624 
                           ,  0.556625,  0.514156,  0.475085,  0.439977 
                           ,  0.408834,  0.380578,  0.355549,  0.333352 
                           ,  0.313647,  0.296093,  0.280351,  0.266195 
                           ,  0.253397,  0.241789,  0.231257,  0.221574 
                           ,  0.212627,  0.204417,  0.196772,  0.189581 
                           ,  0.182956,  0.176784,  0.171008,  0.165515 
                           ,  0.160419,  0.155606,  0.151076,  0.146716 
                           ,  0.142639,  0.138845,  0.135221,  0.131767 
                           ,  0.128482,  0.125368,  0.122424,  0.119592 
                           ,  0.116931,  0.114326,  0.111891,  0.109513 
                           ,  0.107248,  0.105096,  0.103058,  0.101019 
                           ,  0.099151,  0.097282,  0.095527,  0.093715 
                           ,  0.092129,  0.090544,  0.088958,  0.087429 
                           ,  0.086014,  0.084598,  0.083239,  0.081880 
                           ,  0.080634,  0.079388,  0.078143,  0.077010 
                           ,  0.075878,  0.074745,  0.073669,  0.072593 
                           ,  0.071574,  0.070612,  0.069649,  0.068686 
                           ,  0.067780,  0.066874,  0.066025,  0.065176 
                           ,  0.064326,  0.063533,  0.062684,  0.061948 
                           ,  0.061212,  0.060419,  0.059740,  0.059003 
                           ,  0.058324,  0.057644,  0.057022,  0.056342 
                           ,  0.055663,  0.055096,  0.054473,  0.053851 
                           ,  0.053284,  0.052718,  0.052152,  0.051585 
                           ,  0.051019,  0.050566,  0.050000,  0.049490 
                           ,  0.048981,  0.048528,  0.048018,  0.047508 
                           ,  0.047055,  0.046602,  0.046149,  0.045696 
                           ,  0.045300,  0.044904,  0.044451,  0.044054 
                           ,  0.043658,  0.043205,  0.042865,  0.042469 
                           ,  0.042072,  0.041733,  0.041336,  0.040997 
                           ,  0.040657,  0.040260,  0.039921,  0.039581 
                           ,  0.039241,  0.038958,  0.038618,  0.038335 
                           ,  0.037995,  0.037656,  0.037373,  0.037089 
                           ,  0.036806,  0.036467,  0.036183,  0.035900 
                           ,  0.035617,  0.035334,  0.035108,  0.034824 
                           ,  0.034541,  0.034315,  0.034032,  0.033805 
                           ,  0.033522,  0.033296,  0.033069,  0.032786 
                           ,  0.032559,  0.032333,  0.032106,  0.031880 
                           ,  0.031653,  0.031427,  0.031200,  0.030974 
                           ,  0.030804,  0.030578,  0.030351,  0.030125 
                           ,  0.029955,  0.029785,  0.029558,  0.029332 
                           ,  0.029162,  0.028992,  0.028766,  0.028596 
                           ,  0.028426,  0.028199,  0.028086,  0.027860 
                           ,  0.027746,  0.027633,  0.027463,  0.027293 
                           ,  0.027180,  0.027067,  0.026954,  0.026954 
                           ,  0.026840,  0.026727,  0.026727,  0.026614 
                           ,  0.026614,  0.026614,  0.026557,  0.026501 
                           ,  0.026501,  0.026501,  0.026501,  0.026501 
                           ,  0.026501,  0.026501,  0.026501,  0.026387 
                           ,  0.026387,  0.026387,  0.026387,  0.026387 
                           ,  0.026387,  0.026387,  0.026387,  0.026387 
                           ,  0.026387,  0.026387,  0.026387,  0.026387 
                           ,  0.026387,  0.026274,  0.026274,  0.026274 
                           ,  0.026274,  0.026274,  0.026274,  0.026274 
                           ,  0.026274,  0.026274,  0.026274,  0.026274 
                           ,  0.026274,  0.026274,  0.026274,  0.026161 };

  Float_t xtalk[kNpasa]  = {  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000113 
                           ,  0.000793,  0.003058,  0.007305,  0.013194 
                           ,  0.019706,  0.025821,  0.030634,  0.033465 
                           ,  0.034145,  0.032729,  0.029615,  0.025198 
                           ,  0.019989,  0.014496,  0.009003,  0.003964 
                           , -0.000510, -0.004190, -0.007191, -0.009400 
                           , -0.010872, -0.011835, -0.012288, -0.012288 
                           , -0.012005, -0.011495, -0.010872, -0.010136 
                           , -0.009343, -0.008607, -0.007871, -0.007191 
                           , -0.006512, -0.005946, -0.005379, -0.004926 
                           , -0.004473, -0.004077, -0.003737, -0.003398 
                           , -0.003114, -0.002831, -0.002605, -0.002378 
                           , -0.002208, -0.002039, -0.001869, -0.001699 
                           , -0.001585, -0.001472, -0.001359, -0.001246 
                           , -0.001132, -0.001019, -0.001019, -0.000906 
                           , -0.000906, -0.000793, -0.000793, -0.000680 
                           , -0.000680, -0.000680, -0.000566, -0.000566 
                           , -0.000566, -0.000566, -0.000453, -0.000453 
                           , -0.000453, -0.000453, -0.000453, -0.000453 
                           , -0.000340, -0.000340, -0.000340, -0.000340 
                           , -0.000340, -0.000340, -0.000340, -0.000340 
                           , -0.000340, -0.000340, -0.000340, -0.000340 
                           , -0.000340, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000227 
                           , -0.000227, -0.000227, -0.000227, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113, -0.000113 
                           , -0.000113, -0.000113, -0.000113,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 
                           ,  0.000000,  0.000000,  0.000000,  0.000000 };

  // increase CrossTalk to measurements
  for (Int_t ipasa = 0; ipasa < kNpasa; ipasa++) {
    xtalk[ipasa] *= 1.75;
  }

  if (fTRFsmp) delete [] fTRFsmp;
  fTRFsmp = new Float_t[fTRFbin];
  if (fCTsmp)  delete [] fCTsmp;
  fCTsmp  = new Float_t[fTRFbin];

  Float_t loTRF    = TMath::Max(fTRFlo, time[0]);
  Float_t hiTRF    = TMath::Min(fTRFhi, time[kNpasa-1]);
  Float_t binWidth = (hiTRF - loTRF) / ((Float_t) fTRFbin);

  // Take the linear interpolation
  for (Int_t iBin = 0; iBin < fTRFbin; iBin++) {

    Float_t bin = (((Float_t) iBin) + 0.5) * binWidth + loTRF;
    ipos1 = ipos2 = 0;
    diff  = 0;
    do {
      diff = bin - time[ipos2++];
    } while (diff > 0);
    ipos2--;
    if (ipos2 >= kNpasa) ipos2 = kNpasa - 1;
    ipos1 = ipos2 - 1;

    fTRFsmp[iBin] = signal[ipos2] 
                  + diff * (signal[ipos2] - signal[ipos1]) 
                         / (  time[ipos2] -   time[ipos1]);

    fCTsmp[iBin]  = xtalk[ipos2] 
                  + diff * (xtalk[ipos2]  -  xtalk[ipos1]) 
                         / (  time[ipos2] -   time[ipos1]);

  }

}

//_____________________________________________________________________________
Double_t AliTRDSimParam::TimeResponse(Double_t time) const
{
  //
  // Applies the preamp shaper time response
  //

  Int_t iBin = ((Int_t) ((time) / fTRFwid));
  //Int_t iBin = ((Int_t) ((time - fTRFlo) / fTRFwid)); 
  if ((iBin >= 0) && (iBin < fTRFbin)) {
    return fTRFsmp[iBin];
  }
  else {
    return 0.0;
  }    

}

//_____________________________________________________________________________
Double_t AliTRDSimParam::CrossTalk(Double_t time) const
{
  //
  // Applies the pad-pad capacitive cross talk
  //

  Int_t iBin = ((Int_t) ((time - fTRFlo) / fTRFwid)); 
  if ((iBin >= 0) && (iBin < fTRFbin)) {
    return fCTsmp[iBin];
  }
  else {
    return 0.0;
  }    

}

