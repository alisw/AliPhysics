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

/*
$Log$
Revision 1.23  2001/05/07 08:04:48  cblume
New TRF and PRF. Speedup of the code. Digits from amplification region included

Revision 1.22  2001/03/30 14:40:14  cblume
Update of the digitization parameter

Revision 1.21  2001/03/13 09:30:35  cblume
Update of digitization. Moved digit branch definition to AliTRD

Revision 1.20  2001/02/25 20:19:00  hristov
Minor correction: loop variable declared only once for HP, Sun

Revision 1.19  2001/02/14 18:22:26  cblume
Change in the geometry of the padplane

Revision 1.18  2001/01/26 19:56:57  hristov
Major upgrade of AliRoot code

Revision 1.17  2000/12/08 12:53:27  cblume
Change in Copy() function for HP-compiler

Revision 1.16  2000/12/07 12:20:46  cblume
Go back to array compression. Use sampled PRF to speed up digitization

Revision 1.15  2000/11/23 14:34:08  cblume
Fixed bug in expansion routine of arrays (initialize buffers properly)

Revision 1.14  2000/11/20 08:54:44  cblume
Switch off compression as default

Revision 1.13  2000/11/10 14:57:52  cblume
Changes in the geometry constants for the DEC compiler

Revision 1.12  2000/11/01 14:53:20  cblume
Merge with TRD-develop

Revision 1.1.4.9  2000/10/26 17:00:22  cblume
Fixed bug in CheckDetector()

Revision 1.1.4.8  2000/10/23 13:41:35  cblume
Added protection against Log(0) in the gas gain calulation

Revision 1.1.4.7  2000/10/17 02:27:34  cblume
Get rid of global constants

Revision 1.1.4.6  2000/10/16 01:16:53  cblume
Changed timebin 0 to be the one closest to the readout

Revision 1.1.4.5  2000/10/15 23:34:29  cblume
Faster version of the digitizer

Revision 1.1.4.4  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.4.3  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.1.4.2  2000/09/22 14:41:10  cblume
Bug fix in PRF. Included time response. New structure

Revision 1.10  2000/10/05 07:27:53  cblume
Changes in the header-files by FCA

Revision 1.9  2000/10/02 21:28:19  fca
Removal of useless dependecies via forward declarations

Revision 1.8  2000/06/09 11:10:07  cblume
Compiler warnings and coding conventions, next round

Revision 1.7  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.6  2000/06/07 16:27:32  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.5  2000/05/09 16:38:57  cblume
Removed PadResponse(). Merge problem

Revision 1.4  2000/05/08 15:53:45  cblume
Resolved merge conflict

Revision 1.3  2000/04/28 14:49:27  cblume
Only one declaration of iDict in MakeDigits()

Revision 1.1.4.1  2000/05/08 14:42:04  cblume
Introduced AliTRDdigitsManager

Revision 1.1  2000/02/28 19:00:13  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Creates and handles digits from TRD hits                                 //
//                                                                           //
//  The following effects are included:                                      //
//      - Diffusion                                                          //
//      - ExB effects                                                        //
//      - Gas gain including fluctuations                                    //
//      - Pad-response (simple Gaussian approximation)                       //
//      - Electronics noise                                                  //
//      - Electronics gain                                                   //
//      - Digitization                                                       //
//      - ADC threshold                                                      //
//  The corresponding parameter can be adjusted via the various              //
//  Set-functions. If these parameters are not explicitly set, default       //
//  values are used (see Init-function).                                     //
//  To produce digits from a root-file with TRD-hits use the                 //
//  slowDigitsCreate.C macro.                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TMath.h>
#include <TVector.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

#include "AliRun.h"
#include "AliMagF.h"

#include "AliTRD.h"
#include "AliTRDhit.h"
#include "AliTRDdigitizer.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdataArrayF.h"
#include "AliTRDsegmentArray.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDdigitizer)

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer():TNamed()
{
  //
  // AliTRDdigitizer default constructor
  //

  fInputFile      = NULL;
  fDigits         = NULL;
  fTRD            = NULL;
  fGeo            = NULL;
  fPRFsmp         = NULL;
  fTRFsmp         = NULL;

  fEvent          = 0;
  fGasGain        = 0.0;
  fNoise          = 0.0;
  fChipGain       = 0.0;
  fSinRange       = 0.0;
  fSoutRange      = 0.0;
  fADCoutRange    = 0.0;
  fADCinRange     = 0.0;
  fADCthreshold   = 0;
  fDiffusionOn    = 0;
  fDiffusionT     = 0.0;
  fDiffusionL     = 0.0;
  fElAttachOn     = 0;
  fElAttachProp   = 0.0;
  fExBOn          = 0;
  fOmegaTau       = 0.0;
  fPRFOn          = 0;
  fTRFOn          = 0;
  fDriftVelocity  = 0.0;
  fPadCoupling    = 0.0;
  fTimeCoupling   = 0.0;
  fTimeBinWidth   = 0.0;
  fField          = 0.0;

  fPRFbin         = 0;
  fPRFlo          = 0.0;
  fPRFhi          = 0.0;
  fPRFwid         = 0.0;
  fPRFpad         = 0;
  fTRFbin         = 0;
  fTRFlo          = 0.0;
  fTRFhi          = 0.0;
  fTRFwid         = 0.0;

  fCompress       = kTRUE;
  fVerbose        = 1;
  fSDigits        = kFALSE;

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDdigitizer default constructor
  //

  fInputFile     = NULL;
  fDigits        = NULL;
  fTRD           = NULL;
  fGeo           = NULL;
  fPRFsmp        = NULL;
  fTRFsmp        = NULL;

  fEvent         = 0;

  fCompress      = kTRUE;
  fVerbose       = 1;
  fSDigits       = kFALSE;

  Init();

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const AliTRDdigitizer &d)
{
  //
  // AliTRDdigitizer copy constructor
  //

  ((AliTRDdigitizer &) d).Copy(*this);

}

//_____________________________________________________________________________
AliTRDdigitizer::~AliTRDdigitizer()
{
  //
  // AliTRDdigitizer destructor
  //

  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
  }

  if (fDigits) {
    delete fDigits;
  }

}

//_____________________________________________________________________________
AliTRDdigitizer &AliTRDdigitizer::operator=(const AliTRDdigitizer &d)
{
  //
  // Assignment operator
  //

  if (this != &d) ((AliTRDdigitizer &) d).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDdigitizer::Copy(TObject &d)
{
  //
  // Copy function
  //

  Int_t iBin;

  ((AliTRDdigitizer &) d).fInputFile      = NULL;
  ((AliTRDdigitizer &) d).fDigits         = NULL;
  ((AliTRDdigitizer &) d).fTRD            = NULL;
  ((AliTRDdigitizer &) d).fGeo            = NULL;

  ((AliTRDdigitizer &) d).fEvent          = 0;

  ((AliTRDdigitizer &) d).fGasGain        = fGasGain;
  ((AliTRDdigitizer &) d).fNoise          = fNoise;
  ((AliTRDdigitizer &) d).fChipGain       = fChipGain;
  ((AliTRDdigitizer &) d).fSoutRange      = fSoutRange;
  ((AliTRDdigitizer &) d).fSinRange       = fSinRange;
  ((AliTRDdigitizer &) d).fADCoutRange    = fADCoutRange;
  ((AliTRDdigitizer &) d).fADCinRange     = fADCinRange;
  ((AliTRDdigitizer &) d).fADCthreshold   = fADCthreshold;
  ((AliTRDdigitizer &) d).fDiffusionOn    = fDiffusionOn; 
  ((AliTRDdigitizer &) d).fDiffusionT     = fDiffusionT;
  ((AliTRDdigitizer &) d).fDiffusionL     = fDiffusionL;
  ((AliTRDdigitizer &) d).fElAttachOn     = fElAttachOn;
  ((AliTRDdigitizer &) d).fElAttachProp   = fElAttachProp;
  ((AliTRDdigitizer &) d).fExBOn          = fExBOn;
  ((AliTRDdigitizer &) d).fOmegaTau       = fOmegaTau;
  ((AliTRDdigitizer &) d).fLorentzFactor  = fLorentzFactor;
  ((AliTRDdigitizer &) d).fDriftVelocity  = fDriftVelocity;
  ((AliTRDdigitizer &) d).fPadCoupling    = fPadCoupling;
  ((AliTRDdigitizer &) d).fTimeCoupling   = fTimeCoupling;
  ((AliTRDdigitizer &) d).fTimeBinWidth   = fTimeBinWidth;
  ((AliTRDdigitizer &) d).fField          = fField;
  ((AliTRDdigitizer &) d).fPRFOn          = fPRFOn;
  ((AliTRDdigitizer &) d).fTRFOn          = fTRFOn;

  ((AliTRDdigitizer &) d).fCompress       = fCompress;
  ((AliTRDdigitizer &) d).fVerbose        = fVerbose;
  ((AliTRDdigitizer &) d).fSDigits        = fSDigits;

  ((AliTRDdigitizer &) d).fPRFbin         = fPRFbin;
  ((AliTRDdigitizer &) d).fPRFlo          = fPRFlo;
  ((AliTRDdigitizer &) d).fPRFhi          = fPRFhi;
  ((AliTRDdigitizer &) d).fPRFwid         = fPRFwid;
  ((AliTRDdigitizer &) d).fPRFpad         = fPRFpad;
  if (((AliTRDdigitizer &) d).fPRFsmp) delete ((AliTRDdigitizer &) d).fPRFsmp;
  ((AliTRDdigitizer &) d).fPRFsmp = new Float_t[fPRFbin];
  for (iBin = 0; iBin < fPRFbin; iBin++) {
    ((AliTRDdigitizer &) d).fPRFsmp[iBin] = fPRFsmp[iBin];
  }                                                                             
  ((AliTRDdigitizer &) d).fTRFbin         = fTRFbin;
  ((AliTRDdigitizer &) d).fTRFlo          = fTRFlo;
  ((AliTRDdigitizer &) d).fTRFhi          = fTRFhi;
  ((AliTRDdigitizer &) d).fTRFwid         = fTRFwid;
  if (((AliTRDdigitizer &) d).fTRFsmp) delete ((AliTRDdigitizer &) d).fTRFsmp;
  ((AliTRDdigitizer &) d).fTRFsmp = new Float_t[fTRFbin];
  for (iBin = 0; iBin < fTRFbin; iBin++) {
    ((AliTRDdigitizer &) d).fTRFsmp[iBin] = fTRFsmp[iBin];
  }                                      
                                       
}

//_____________________________________________________________________________
Int_t AliTRDdigitizer::Diffusion(Float_t driftlength, Float_t *xyz)
{
  //
  // Applies the diffusion smearing to the position of a single electron
  //

  Float_t driftSqrt = TMath::Sqrt(driftlength);
  Float_t sigmaT = driftSqrt * fDiffusionT;
  Float_t sigmaL = driftSqrt * fDiffusionL;
  xyz[0] = gRandom->Gaus(xyz[0], sigmaL * fLorentzFactor);
  xyz[1] = gRandom->Gaus(xyz[1], sigmaT * fLorentzFactor);
  xyz[2] = gRandom->Gaus(xyz[2], sigmaT);

  return 1;

}

//_____________________________________________________________________________
Int_t AliTRDdigitizer::ExB(Float_t driftlength, Float_t *xyz)
{
  //
  // Applies E x B effects to the position of a single electron
  //

  xyz[0] = xyz[0];
  xyz[1] = xyz[1] + fOmegaTau * driftlength;
  xyz[2] = xyz[2];

  return 1;

}

//_____________________________________________________________________________
Int_t AliTRDdigitizer::PadResponse(Float_t signal, Float_t dist, Float_t *pad)
{
  //
  // Applies the pad response
  //

  Int_t iBin =  ((Int_t) (( - dist - fPRFlo) / fPRFwid));

  Int_t iBin0 = iBin - fPRFpad;
  Int_t iBin1 = iBin;
  Int_t iBin2 = iBin + fPRFpad;

  if ((iBin0 >= 0) && (iBin2 < fPRFbin)) {

    pad[0] = signal * fPRFsmp[iBin0];
    pad[1] = signal * fPRFsmp[iBin1];
    pad[2] = signal * fPRFsmp[iBin2];

    return 1;

  }
  else {

    return 0;

  }

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::TimeResponse(Float_t time)
{
  //
  // Applies the preamp shaper time response
  //

  Int_t iBin = ((Int_t) ((time - fTRFlo) / fTRFwid)); 
  if ((iBin >= 0) && (iBin < fTRFbin)) {
    return fTRFsmp[iBin];
  }
  else {
    return 0.0;
  }    

}

//_____________________________________________________________________________
void AliTRDdigitizer::Init()
{
  //
  // Initializes the digitization procedure with standard values
  //

  // The default parameter for the digitization
  fGasGain        = 2800.;
  fChipGain       = 6.1;
  fNoise          = 1000.;
  fADCoutRange    = 1023.;          // 10-bit ADC
  fADCinRange     = 1000.;          // 1V input range
  fADCthreshold   = 1;

  // For the summable digits
  fSinRange       = 1000000.;
  fSoutRange      = 1000000.;

  // The drift velocity (cm / mus)
  fDriftVelocity  = 1.5;

  // The magnetic field strength in Tesla
  fField          = 0.2 * gAlice->Field()->Factor();

  // Diffusion on
  fDiffusionOn    = 1;

  // E x B effects
  fExBOn          = 0;

  // Propability for electron attachment
  fElAttachOn     = 0;
  fElAttachProp   = 0.0;

  // The pad response function
  fPRFOn          =  1;

  // The time response function
  fTRFOn          =  1;

  // The pad coupling factor (same number as for the TPC)
  fPadCoupling    = 0.5;

  // The time coupling factor (same number as for the TPC)
  fTimeCoupling   = 0.4;

}

//_____________________________________________________________________________
void AliTRDdigitizer::ReInit()
{
  //
  // Reinitializes the digitization procedure after a change in the parameter
  //

  if (!fGeo) {
    printf("AliTRDdigitizer::ReInit -- ");
    printf("No geometry defined. Run InitDetector() first\n");
    exit(1);
  }

  // Calculate the time bin width in ns
  fTimeBinWidth   = fGeo->GetTimeBinSize() / fDriftVelocity * 1000.0;

  // The range and the binwidth for the sampled TRF 
  fTRFbin = 100;
  // Start 0.2 mus before the signal
  fTRFlo  = -0.2 * fDriftVelocity;
  // End the maximum driftlength after the signal 
  fTRFhi  = AliTRDgeometry::DrThick() 
          + fGeo->GetTimeAfter() * fGeo->GetTimeBinSize();
  fTRFwid = (fTRFhi - fTRFlo) / ((Float_t) fTRFbin);

  // Transverse and longitudinal diffusion coefficients (Xe/CO2)
  fDiffusionT     = GetDiffusionT(fDriftVelocity,fField);
  fDiffusionL     = GetDiffusionL(fDriftVelocity,fField);

  // omega * tau.= tan(Lorentz-angle)
  fOmegaTau       = GetOmegaTau(fDriftVelocity,fField);

  // The Lorentz factor
  if (fExBOn) {
    fLorentzFactor = 1.0 / (1.0 + fOmegaTau*fOmegaTau);
  }
  else {
    fLorentzFactor = 1.0;
  }

}

//_____________________________________________________________________________
void AliTRDdigitizer::SampleTRF()
{
  //
  // Samples the time response function
  // It is defined according to Vasiles simulation of the preamp shaper
  // output and includes the effect of the ion tail (based on Tariqs 
  // Garfield simulation) and a shaping time of 125 ns FWHM
  //

  Int_t   ipos1;
  Int_t   ipos2;
  Float_t diff;

  const Int_t kNpasa = 200;
  Float_t time[kNpasa]   = {  -0.280000, -0.270000, -0.260000, -0.250000
                            , -0.240000, -0.230000, -0.220000, -0.210000
                            , -0.200000, -0.190000, -0.180000, -0.170000
                            , -0.160000, -0.150000, -0.140000, -0.130000
                            , -0.120000, -0.110000, -0.100000, -0.090000
                            , -0.080000, -0.070000, -0.060000, -0.050000
                            , -0.040000, -0.030000, -0.020000, -0.010000
                            , -0.000000,  0.010000,  0.020000,  0.030000
                            ,  0.040000,  0.050000,  0.060000,  0.070000
                            ,  0.080000,  0.090000,  0.100000,  0.110000
                            ,  0.120000,  0.130000,  0.140000,  0.150000
                            ,  0.160000,  0.170000,  0.180000,  0.190000
                            ,  0.200000,  0.210000,  0.220000,  0.230000
                            ,  0.240000,  0.250000,  0.260000,  0.270000
                            ,  0.280000,  0.290000,  0.300000,  0.310000
                            ,  0.320000,  0.330000,  0.340000,  0.350000
                            ,  0.360000,  0.370000,  0.380000,  0.390000
                            ,  0.400000,  0.410000,  0.420000,  0.430000
                            ,  0.440000,  0.450000,  0.460000,  0.470000
                            ,  0.480000,  0.490000,  0.500000,  0.510000
                            ,  0.520000,  0.530000,  0.540000,  0.550000
                            ,  0.560000,  0.570000,  0.580000,  0.590000
                            ,  0.600000,  0.610000,  0.620000,  0.630000
                            ,  0.640000,  0.650000,  0.660000,  0.670000
                            ,  0.680000,  0.690000,  0.700000,  0.710000
                            ,  0.720000,  0.730000,  0.740000,  0.750000
                            ,  0.760000,  0.770000,  0.780000,  0.790000
                            ,  0.800000,  0.810000,  0.820000,  0.830000
                            ,  0.840000,  0.850000,  0.860000,  0.870000
                            ,  0.880000,  0.890000,  0.900000,  0.910000
                            ,  0.920000,  0.930000,  0.940000,  0.950000
                            ,  0.960000,  0.970000,  0.980000,  0.990000
                            ,  1.000000,  1.010000,  1.020000,  1.030000
                            ,  1.040000,  1.050000,  1.060000,  1.070000
                            ,  1.080000,  1.090000,  1.100000,  1.110000
                            ,  1.120000,  1.130000,  1.140000,  1.150000
                            ,  1.160000,  1.170000,  1.180000,  1.190000
                            ,  1.200000,  1.210000,  1.220000,  1.230000
                            ,  1.240000,  1.250000,  1.260000,  1.270000
                            ,  1.280000,  1.290000,  1.300000,  1.310000
                            ,  1.320000,  1.330000,  1.340000,  1.350000
                            ,  1.360000,  1.370000,  1.380000,  1.390000
                            ,  1.400000,  1.410000,  1.420000,  1.430000
                            ,  1.440000,  1.450000,  1.460000,  1.470000
                            ,  1.480000,  1.490000,  1.500000,  1.510000
                            ,  1.520000,  1.530000,  1.540000,  1.550000
                            ,  1.560000,  1.570000,  1.580000,  1.590000
                            ,  1.600000,  1.610000,  1.620000,  1.630000
                            ,  1.640000,  1.650000,  1.660000,  1.670000
			    ,  1.680000,  1.690000,  1.700000,  1.710000 };

  Float_t signal[kNpasa] = {   0.000000,  0.000000,  0.000000,  0.000000 
                            ,  0.000000,  0.000000,  0.000000,  0.000000
                            ,  0.000000,  0.000000,  0.000000,  0.000000
                            ,  0.000000,  0.000000,  0.000000,  0.000098
                            ,  0.003071,  0.020056,  0.066053,  0.148346
                            ,  0.263120,  0.398496,  0.540226,  0.674436
                            ,  0.790977,  0.883083,  0.947744,  0.985714
                            ,  0.999248,  0.992105,  0.967669,  0.930827
                            ,  0.884586,  0.833083,  0.778571,  0.723684
                            ,  0.669173,  0.617293,  0.567669,  0.521805
                            ,  0.479699,  0.440977,  0.405639,  0.373985
                            ,  0.345526,  0.320038,  0.297256,  0.276917
                            ,  0.258797,  0.242632,  0.228195,  0.215301
                            ,  0.203759,  0.193383,  0.184023,  0.175564
                            ,  0.167895,  0.160940,  0.154549,  0.148722
                            ,  0.143308,  0.138346,  0.133722,  0.129398
                            ,  0.125376,  0.121617,  0.118045,  0.114699
                            ,  0.111541,  0.108571,  0.105714,  0.103008
                            ,  0.100414,  0.097970,  0.095602,  0.093346
                            ,  0.091165,  0.089060,  0.087068,  0.085150
                            ,  0.083308,  0.081541,  0.079812,  0.078158
                            ,  0.076541,  0.075000,  0.073496,  0.072068
                            ,  0.070677,  0.069286,  0.068008,  0.066729
                            ,  0.065489,  0.064286,  0.063120,  0.061992
                            ,  0.060902,  0.059850,  0.058797,  0.057820
                            ,  0.056842,  0.055902,  0.054962,  0.054060
                            ,  0.053158,  0.052293,  0.051466,  0.050639
                            ,  0.049850,  0.049060,  0.048308,  0.047556
                            ,  0.046842,  0.046128,  0.045451,  0.044774
                            ,  0.044098,  0.043459,  0.042820,  0.042218
                            ,  0.041617,  0.041015,  0.040451,  0.039887
                            ,  0.039323,  0.038797,  0.038271,  0.037744
                            ,  0.037237,  0.036744,  0.036259,  0.035786
                            ,  0.035323,  0.034872,  0.034429,  0.033996
                            ,  0.033575,  0.033162,  0.032756,  0.032361
                            ,  0.031974,  0.031594,  0.031222,  0.030857
                            ,  0.030496,  0.030143,  0.029793,  0.029451
                            ,  0.029109,  0.028774,  0.028444,  0.028113
                            ,  0.027793,  0.027477,  0.027165,  0.026861
                            ,  0.026564,  0.026271,  0.025981,  0.025699
                            ,  0.025421,  0.025147,  0.024880,  0.024613
                            ,  0.024353,  0.024094,  0.023842,  0.023590
                            ,  0.023346,  0.023102,  0.022865,  0.022628
                            ,  0.022398,  0.022173,  0.021951,  0.021733
                            ,  0.021519,  0.021308,  0.021098,  0.020891
                            ,  0.020688,  0.020485,  0.020286,  0.020090
                            ,  0.019895,  0.019707,  0.019519,  0.019335
                            ,  0.019150,  0.018974,  0.018797,  0.018624
                            ,  0.018451,  0.018282,  0.018113,  0.017947
			    ,  0.017782,  0.017617,  0.017455,  0.017297 };

  //for (Int_t ipasa = 0; ipasa < kNpasa; ipasa++) {
  //  time[ipasa] += 0.13; 
  //}

  if (fTRFsmp) delete fTRFsmp;
  fTRFsmp = new Float_t[fTRFbin];

  Float_t loTRF    = TMath::Max(fTRFlo / fDriftVelocity,time[0]);
  Float_t hiTRF    = TMath::Min(fTRFhi / fDriftVelocity,time[kNpasa-1]);
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
    if (ipos2 > kNpasa) ipos2 = kNpasa - 1;
    ipos1 = ipos2 - 1;

    fTRFsmp[iBin] = signal[ipos2] 
                  + diff * (signal[ipos2] - signal[ipos1]) 
                         / (  time[ipos2] -   time[ipos1]);

  }

}

//_____________________________________________________________________________
void AliTRDdigitizer::SamplePRF()
{
  //
  // Samples the pad response function
  //

  const Int_t kPRFbin = 61;
  Float_t prf[kPRFbin] = { 0.002340, 0.003380, 0.004900, 0.007080, 0.010220
                         , 0.014740, 0.021160, 0.030230, 0.042800, 0.059830
                         , 0.082030, 0.109700, 0.142550, 0.179840, 0.220610
                         , 0.263980, 0.309180, 0.355610, 0.402790, 0.450350
                         , 0.497930, 0.545190, 0.591740, 0.637100, 0.680610
                         , 0.721430, 0.758400, 0.790090, 0.814720, 0.830480
                         , 0.835930, 0.830480, 0.814710, 0.790070, 0.758390
                         , 0.721410, 0.680590, 0.637080, 0.591730, 0.545180
                         , 0.497920, 0.450340, 0.402790, 0.355610, 0.309190
                         , 0.263990, 0.220630, 0.179850, 0.142570, 0.109720
                         , 0.082040, 0.059830, 0.042820, 0.030230, 0.021170
                         , 0.014740, 0.010230, 0.007080, 0.004900, 0.003380
			 , 0.002340 };

  fPRFlo  = -1.5;
  fPRFhi  =  1.5;
  fPRFbin = kPRFbin;
  fPRFwid = (fPRFhi - fPRFlo) / ((Float_t) fPRFbin);
  fPRFpad = ((Int_t) (1.0 / fPRFwid));

  if (fPRFsmp) delete fPRFsmp;
  fPRFsmp = new Float_t[fPRFbin];
  for (Int_t iBin = 0; iBin < fPRFbin; iBin++) {
    fPRFsmp[iBin] = prf[iBin];
  }

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::Open(const Char_t *name, Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the hit-tree
  //

  // Connect the AliRoot file containing Geometry, Kine, and Hits
  fInputFile = (TFile*) gROOT->GetListOfFiles()->FindObject(name);
  if (!fInputFile) {
    printf("AliTRDdigitizer::Open -- ");
    printf("Open the ALIROOT-file %s.\n",name);
    fInputFile = new TFile(name,"UPDATE");
  }
  else {
    printf("AliTRDdigitizer::Open -- ");
    printf("%s is already open.\n",name);
  }

  gAlice = (AliRun*) fInputFile->Get("gAlice");
  if (gAlice) {
    printf("AliTRDdigitizer::Open -- ");
    printf("AliRun object found on file.\n");
  }
  else {
    printf("AliTRDdigitizer::Open -- ");
    printf("Could not find AliRun object.\n");
    return kFALSE;
  }

  fEvent = nEvent;

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(fEvent);
  if (nparticles <= 0) {
    printf("AliTRDdigitizer::Open -- ");
    printf("No entries in the trees for event %d.\n",fEvent);
    return kFALSE;
  }

  return InitDetector();

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::InitDetector()
{
  //
  // Sets the pointer to the TRD detector and the geometry
  //

  // Get the pointer to the detector class and check for version 1
  fTRD = (AliTRD*) gAlice->GetDetector("TRD");
  if (fTRD->IsVersion() != 1) {
    printf("AliTRDdigitizer::InitDetector -- ");
    printf("TRD must be version 1 (slow simulator).\n");
    exit(1);
  }

  // Get the geometry
  fGeo = fTRD->GetGeometry();
  printf("AliTRDdigitizer::InitDetector -- ");
  printf("Geometry version %d\n",fGeo->IsVersion());

  ReInit();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::SumSDigits()
{
  //
  // Sums up the summable digits and creates final digits
  // Not yet implemented
  //

  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::MakeDigits()
{
  //
  // Creates digits.
  //

  ///////////////////////////////////////////////////////////////
  // Parameter 
  ///////////////////////////////////////////////////////////////

  // Converts number of electrons to fC
  const Double_t kEl2fC  = 1.602E-19 * 1.0E15; 

  ///////////////////////////////////////////////////////////////

  // Number of pads included in the pad response
  const Int_t kNpad  = 3;

  // Number of track dictionary arrays
  const Int_t kNDict = AliTRDdigitsManager::kNDict;

  // Half the width of the amplification region
  const Float_t kAmWidth = AliTRDgeometry::AmThick() / 2.;

  Int_t   iRow, iCol, iTime, iPad;
  Int_t   iDict  = 0;
  Int_t   nBytes = 0;

  Int_t   totalSizeDigits = 0;
  Int_t   totalSizeDict0  = 0;
  Int_t   totalSizeDict1  = 0;
  Int_t   totalSizeDict2  = 0;

  Int_t   timeTRDbeg = 0;
  Int_t   timeTRDend = 1;

  Float_t pos[3];
  Float_t rot[3];
  Float_t xyz[3];
  Float_t padSignal[kNpad];
  Float_t signalOld[kNpad];

  AliTRDdataArrayF *signals = 0;
  AliTRDdataArrayI *digits  = 0;
  AliTRDdataArrayI *dictionary[kNDict];

  // Create a digits manager
  fDigits = new AliTRDdigitsManager();

  // Create a container for the amplitudes
  AliTRDsegmentArray *signalsArray 
                     = new AliTRDsegmentArray("AliTRDdataArrayF",AliTRDgeometry::Ndet());

  if (fTRFOn) {
    timeTRDbeg = ((Int_t) (-fTRFlo / fGeo->GetTimeBinSize())) - 1;
    timeTRDend = ((Int_t) ( fTRFhi / fGeo->GetTimeBinSize())) - 1;
    printf("AliTRDdigitizer::MakeDigits -- ");
    printf("Sample the TRF between -%d and %d\n",timeTRDbeg,timeTRDend);
  }

  Float_t elAttachProp = fElAttachProp / 100.; 

  // Create the sampled PRF
  SamplePRF();

  // Create the sampled TRF
  SampleTRF();

  if (!fGeo) {
    printf("AliTRDdigitizer::MakeDigits -- ");
    printf("No geometry defined\n");
    return kFALSE;
  }

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Start creating digits.\n");
  if (fVerbose > 0) this->Dump();

  // Get the pointer to the hit tree
  TTree *HitTree = gAlice->TreeH();

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) HitTree->GetEntries();
  if (fVerbose > 0) {
    printf("AliTRDdigitizer::MakeDigits -- ");
    printf("Found %d primary particles\n",nTrack);
  } 

  Int_t detectorOld = -1;
  Int_t countHits   =  0; 

  // Loop through all entries in the tree
  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

    gAlice->ResetHits();
    nBytes += HitTree->GetEvent(iTrack);

    // Get the number of hits in the TRD created by this particle
    Int_t nHit = fTRD->Hits()->GetEntriesFast();
    if (fVerbose > 0) {
      printf("AliTRDdigitizer::MakeDigits -- ");
      printf("Found %d hits for primary particle %d\n",nHit,iTrack);
    }

    // Loop through the TRD hits  
    for (Int_t iHit = 0; iHit < nHit; iHit++) {

      countHits++;

      AliTRDhit *hit = (AliTRDhit *) fTRD->Hits()->UncheckedAt(iHit);
              pos[0]   = hit->X();
              pos[1]   = hit->Y();
              pos[2]   = hit->Z();
      Float_t q        = hit->GetCharge();
      Int_t   track    = hit->Track();
      Int_t   detector = hit->GetDetector();
      Int_t   plane    = fGeo->GetPlane(detector);
      Int_t   sector   = fGeo->GetSector(detector);
      Int_t   chamber  = fGeo->GetChamber(detector);

      if (!(CheckDetector(plane,chamber,sector))) continue;

      Int_t   nRowMax     = fGeo->GetRowMax(plane,chamber,sector);
      Int_t   nColMax     = fGeo->GetColMax(plane);
      Int_t   nTimeMax    = fGeo->GetTimeMax();
      Int_t   nTimeBefore = fGeo->GetTimeBefore();
      Int_t   nTimeAfter  = fGeo->GetTimeAfter();
      Int_t   nTimeTotal  = fGeo->GetTimeTotal();
      Float_t row0        = fGeo->GetRow0(plane,chamber,sector);
      Float_t col0        = fGeo->GetCol0(plane);
      Float_t time0       = fGeo->GetTime0(plane);
      Float_t rowPadSize  = fGeo->GetRowPadSize(plane,chamber,sector);
      Float_t colPadSize  = fGeo->GetColPadSize(plane);
      Float_t timeBinSize = fGeo->GetTimeBinSize();
      Float_t divideRow   = 1.0 / rowPadSize;
      Float_t divideCol   = 1.0 / colPadSize;
      Float_t divideTime  = 1.0 / timeBinSize;

      if (fVerbose > 1) {
        printf("Analyze hit no. %d ",iHit);
        printf("-----------------------------------------------------------\n");
        hit->Dump();
        printf("plane = %d, sector = %d, chamber = %d\n"
              ,plane,sector,chamber);
        printf("nRowMax = %d, nColMax = %d, nTimeMax = %d\n" 
              ,nRowMax,nColMax,nTimeMax);
        printf("nTimeBefore = %d, nTimeAfter = %d, nTimeTotal = %d\n"
	      ,nTimeBefore,nTimeAfter,nTimeTotal);
        printf("row0 = %f, col0 = %f, time0 = %f\n"
              ,row0,col0,time0);
        printf("rowPadSize = %f, colPadSize = %f, timeBinSize = %f\n"
	       ,rowPadSize,colPadSize,timeBinSize); 
      }
       
      // Don't analyze test hits
      if (hit->FromTest()) continue;

      if (detector != detectorOld) {

        if (fVerbose > 1) {
          printf("AliTRDdigitizer::MakeDigits -- ");
          printf("Get new container. New det = %d, Old det = %d\n"
                ,detector,detectorOld);
	}
        // Compress the old one if enabled
        if ((fCompress) && (detectorOld > -1)) {
          if (fVerbose > 1) {
            printf("AliTRDdigitizer::MakeDigits -- ");
            printf("Compress the old container ...");
	  }
          signals->Compress(1,0);
          for (iDict = 0; iDict < kNDict; iDict++) {
            dictionary[iDict]->Compress(1,0);
	  }
          if (fVerbose > 1) printf("done\n");
	}
	// Get the new container
        signals = (AliTRDdataArrayF *) signalsArray->At(detector);
        if (signals->GetNtime() == 0) {
          // Allocate a new one if not yet existing
          if (fVerbose > 1) {
            printf("AliTRDdigitizer::MakeDigits -- ");
            printf("Allocate a new container ... ");
	  }
          signals->Allocate(nRowMax,nColMax,nTimeTotal);
	}
        else {
	  // Expand an existing one
          if (fCompress) {
            if (fVerbose > 1) {
              printf("AliTRDdigitizer::MakeDigits -- ");
              printf("Expand an existing container ... ");
	    }
            signals->Expand();
	  }
	}
	// The same for the dictionary
        for (iDict = 0; iDict < kNDict; iDict++) {       
          dictionary[iDict] = fDigits->GetDictionary(detector,iDict);
          if (dictionary[iDict]->GetNtime() == 0) {
            dictionary[iDict]->Allocate(nRowMax,nColMax,nTimeTotal);
	  }
          else {
            if (fCompress) dictionary[iDict]->Expand();
	  }
        }      
        if (fVerbose > 1) printf("done\n");
        detectorOld = detector;
      }

      // Rotate the sectors on top of each other
      fGeo->Rotate(detector,pos,rot);

      // The driftlength. It is negative if the hit is in the 
      // amplification region.
      Float_t driftlength = time0 - rot[0];

      // Take also the drift in the amplification region into account
      // The drift length is at the moment still the same, regardless of
      // the position relativ to the wire. This non-isochronity needs still
      // to be implemented.
      Float_t driftlengthL = TMath::Abs(driftlength + kAmWidth);
      if (fExBOn) driftlengthL /= TMath::Sqrt(fLorentzFactor);

      // Loop over all electrons of this hit
      // TR photons produce hits with negative charge
      Int_t nEl = ((Int_t) TMath::Abs(q));
      for (Int_t iEl = 0; iEl < nEl; iEl++) {

        xyz[0] = rot[0];
        xyz[1] = rot[1];
        xyz[2] = rot[2];

        // Electron attachment
        if (fElAttachOn) {
          if (gRandom->Rndm() < (driftlengthL * elAttachProp)) 
            continue;
        }

        // Apply the diffusion smearing
        if (fDiffusionOn) {
          if (!(Diffusion(driftlengthL,xyz))) continue;
	}

        // Apply E x B effects (depends on drift direction)
        if (fExBOn) { 
          if (!(ExB(driftlength+kAmWidth,xyz))) continue;   
	}

        // The electron position after diffusion and ExB in pad coordinates 
        // The pad row (z-direction)
        Int_t  rowE = ((Int_t) ((xyz[2] -  row0) * divideRow));
        if ((rowE < 0) || (rowE >= nRowMax)) continue;   

        // The pad column (rphi-direction)
        Int_t  colE = ((Int_t) ((xyz[1] -  col0) * divideCol));
        if ((colE < 0) || (colE >= nColMax)) continue;   
  
        // The time bin (negative for hits in the amplification region)
	// In the amplification region the electrons drift from both sides
	// to the middle (anode wire plane)
        Float_t timeDist   = time0 - xyz[0];
        Float_t timeOffset = 0;
        Int_t   timeE      = 0;
        if (timeDist > 0) {
	  // The time bin
          timeE      = ((Int_t) (timeDist * divideTime));
          // The distance of the position to the middle of the timebin
          timeOffset = ((((Float_t) timeE) + 0.5) * timeBinSize) - timeDist;
	}
        else {
	  // Difference between half of the amplification gap width and
	  // the distance to the anode wire
          Float_t anodeDist = kAmWidth - TMath::Abs(timeDist + kAmWidth);
          // The time bin
          timeE      = -1 * (((Int_t ) (anodeDist * divideTime)) + 1);
          // The distance of the position to the middle of the timebin
          timeOffset = ((((Float_t) timeE) + 0.5) * timeBinSize) + anodeDist;
	}
 
        // Apply the gas gain including fluctuations
        Float_t ggRndm = 0.0;
        do {
          ggRndm = gRandom->Rndm();
	} while (ggRndm <= 0);
        Int_t signal = (Int_t) (-fGasGain * TMath::Log(ggRndm));

        // Apply the pad response 
        if (fPRFOn) {
  	  // The distance of the electron to the center of the pad 
	  // in units of pad width
          Float_t dist = (xyz[1] - col0 - (colE + 0.5) * colPadSize) 
                       * divideCol;
          if (!(PadResponse(signal,dist,padSignal))) continue;
	}
	else {
          padSignal[0] = 0.0;
          padSignal[1] = signal;
          padSignal[2] = 0.0;
	}

	// Sample the time response inside the drift region
	// + additional time bins before and after.
        // The sampling is done always in the middle of the time bin
        for (Int_t iTimeBin = TMath::Max(timeE-timeTRDbeg,        -nTimeBefore) 
  	          ;iTimeBin < TMath::Min(timeE+timeTRDend,nTimeMax+nTimeAfter ) 
        	  ;iTimeBin++) {

   	  // Apply the time response
          Float_t timeResponse = 1.0;
          if (fTRFOn) {
            Float_t time = (iTimeBin - timeE) * timeBinSize + timeOffset;
            timeResponse = TimeResponse(time);
	  }

          signalOld[0] = 0.0;
          signalOld[1] = 0.0;
          signalOld[2] = 0.0;

          for (iPad = 0; iPad < kNpad; iPad++) {

            Int_t colPos = colE + iPad - 1;
            if (colPos <        0) continue;
            if (colPos >= nColMax) break;

            // Add the signals
	    // Note: The time bin number is shifted by nTimeBefore to avoid negative
	    // time bins. This has to be subtracted later.
            Int_t iCurrentTimeBin = iTimeBin + nTimeBefore;
            signalOld[iPad]  = signals->GetDataUnchecked(rowE,colPos,iCurrentTimeBin);
            signalOld[iPad] += padSignal[iPad] * timeResponse;
            signals->SetDataUnchecked(rowE,colPos,iCurrentTimeBin,signalOld[iPad]);

            // Store the track index in the dictionary
            // Note: We store index+1 in order to allow the array to be compressed
            if (signalOld[iPad] > 0) {
              for (iDict = 0; iDict < kNDict; iDict++) {
                Int_t oldTrack = dictionary[iDict]->GetDataUnchecked(rowE
                                                                    ,colPos
                                                                    ,iCurrentTimeBin);
                if (oldTrack == track+1) break;
                if (oldTrack ==       0) {
                  dictionary[iDict]->SetDataUnchecked(rowE,colPos,iCurrentTimeBin,track+1);
                  break;
                }
              }
            }

	  }

	}

      }

    }

  } // All hits finished

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Finished analyzing %d hits\n",countHits);

  // The total conversion factor
  Float_t convert = kEl2fC * fPadCoupling * fTimeCoupling * fChipGain;

  // Loop through all chambers to finalize the digits
  for (Int_t iDet = 0; iDet < AliTRDgeometry::Ndet(); iDet++) {

    Int_t plane       = fGeo->GetPlane(iDet);
    Int_t sector      = fGeo->GetSector(iDet);
    Int_t chamber     = fGeo->GetChamber(iDet);
    Int_t nRowMax     = fGeo->GetRowMax(plane,chamber,sector);
    Int_t nColMax     = fGeo->GetColMax(plane);
    Int_t nTimeMax    = fGeo->GetTimeMax();
    Int_t nTimeTotal  = fGeo->GetTimeTotal();

    if (fVerbose > 0) {
      printf("AliTRDdigitizer::MakeDigits -- ");
      printf("Digitization for chamber %d\n",iDet);
    }

    // Add a container for the digits of this detector
    digits = fDigits->GetDigits(iDet);        
    // Allocate memory space for the digits buffer
    digits->Allocate(nRowMax,nColMax,nTimeTotal);

    // Get the signal container
    signals = (AliTRDdataArrayF *) signalsArray->At(iDet);
    if (signals->GetNtime() == 0) {
      // Create missing containers
      signals->Allocate(nRowMax,nColMax,nTimeTotal);      
    }
    else {
      // Expand the container if neccessary
      if (fCompress) signals->Expand();
    }
    // Create the missing dictionary containers
    for (iDict = 0; iDict < kNDict; iDict++) {       
      dictionary[iDict] = fDigits->GetDictionary(iDet,iDict);
      if (dictionary[iDict]->GetNtime() == 0) {
        dictionary[iDict]->Allocate(nRowMax,nColMax,nTimeTotal);
      }
    }

    Int_t nDigits = 0;

    // Don't create noise in detectors that are switched off
    if (CheckDetector(plane,chamber,sector)) {

      // Create the digits for this chamber
      for (iRow  = 0; iRow  <  nRowMax;   iRow++ ) {
        for (iCol  = 0; iCol  <  nColMax;   iCol++ ) {
          for (iTime = 0; iTime < nTimeTotal; iTime++) {         

	    // Create summable digits
            if (fSDigits) {

              Float_t signalAmp = signals->GetDataUnchecked(iRow,iCol,iTime);
              Int_t adc  = 0;
              if (signalAmp >= fSinRange) {
                adc = ((Int_t) fSoutRange);
	      }
              else {
                adc = ((Int_t) (signalAmp * (fSoutRange / fSinRange)));
	      }
              nDigits++;
              digits->SetDataUnchecked(iRow,iCol,iTime,adc);

	    }
	    // Create normal digits
            else {

              Float_t signalAmp = signals->GetDataUnchecked(iRow,iCol,iTime);

              // Add the noise
              signalAmp  = TMath::Max((Double_t) gRandom->Gaus(signalAmp,fNoise),0.0);
              // Convert to mV
              signalAmp *= convert;
	      // Convert to ADC counts. Set the overflow-bit fADCoutRange if the 
	      // signal is larger than fADCinRange
              Int_t adc  = 0;
              if (signalAmp >= fADCinRange) {
                adc = ((Int_t) fADCoutRange);
	      }
              else {
                adc = ((Int_t) (signalAmp * (fADCoutRange / fADCinRange)));
	      }

              // Store the amplitude of the digit if above threshold
              if (adc > fADCthreshold) {
                if (fVerbose > 2) {
                  printf("  iRow = %d, iCol = %d, iTime = %d\n"
                        ,iRow,iCol,iTime);
                  printf("  signal = %f, adc = %d\n",signalAmp,adc);
	        }
                nDigits++;
                digits->SetDataUnchecked(iRow,iCol,iTime,adc);
  	      }

	    }

	  }
        }
      }

    }

    // Compress the arrays
    digits->Compress(1,0);
    for (iDict = 0; iDict < kNDict; iDict++) {
      dictionary[iDict]->Compress(1,0);
    }

    totalSizeDigits += digits->GetSize();
    totalSizeDict0  += dictionary[0]->GetSize();
    totalSizeDict1  += dictionary[1]->GetSize();
    totalSizeDict2  += dictionary[2]->GetSize();

    Float_t nPixel = nRowMax * nColMax * nTimeMax;
    printf("AliTRDdigitizer::MakeDigits -- ");
    printf("Found %d digits in detector %d (%3.0f).\n"
          ,nDigits,iDet
          ,100.0 * ((Float_t) nDigits) / nPixel);
 
    if (fCompress) signals->Compress(1,0);

  }

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Total number of analyzed hits = %d\n",countHits);

  printf("AliTRDdigitizer::MakeDigits -- ");
  printf("Total digits data size = %d, %d, %d, %d\n",totalSizeDigits
                                                    ,totalSizeDict0
                                                    ,totalSizeDict1
                                                    ,totalSizeDict2);        

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::CheckDetector(Int_t plane, Int_t chamber, Int_t sector)
{
  //
  // Checks whether a detector is enabled
  //

  if ((fTRD->GetSensChamber() >=       0) &&
      (fTRD->GetSensChamber() != chamber)) return kFALSE;
  if ((fTRD->GetSensPlane()   >=       0) &&
      (fTRD->GetSensPlane()   !=   plane)) return kFALSE;
  if ( fTRD->GetSensSector()  >=       0) {
    Int_t sens1 = fTRD->GetSensSector();
    Int_t sens2 = sens1 + fTRD->GetSensSectorRange();
    sens2 -= ((Int_t) (sens2 / AliTRDgeometry::Nsect())) 
           * AliTRDgeometry::Nsect();
    if (sens1 < sens2) {
      if ((sector < sens1) || (sector >= sens2)) return kFALSE;
    }
    else {
      if ((sector < sens1) && (sector >= sens2)) return kFALSE;
    }
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::WriteDigits()
{
  //
  // Writes out the TRD-digits and the dictionaries
  //

  // Create the branches
  if (!(gAlice->TreeD()->GetBranch("TRDdigits"))) { 
    return kFALSE;
  }

  // Store the digits and the dictionary in the tree
  fDigits->WriteDigits();

  // Write the new tree into the input file (use overwrite option)
  Char_t treeName[15];
  sprintf(treeName,"TreeD%d",fEvent);
  printf("AliTRDdigitizer::WriteDigits -- ");
  printf("Write the digits tree %s for event %d.\n"
        ,treeName,fEvent);
  gAlice->TreeD()->Write(treeName,TObject::kOverwrite);
 
  return kTRUE;

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::GetDiffusionL(Float_t vd, Float_t b)
{
  //
  // Returns the longitudinal diffusion coefficient for a given drift 
  // velocity <vd> and a B-field <b> for Xe/CO2 (15%).
  // The values are according to a GARFIELD simulation.
  //

  const Int_t kNb = 5;
  Float_t p0[kNb] = {  0.007440,  0.007493,  0.007513,  0.007672,  0.007831 };
  Float_t p1[kNb] = {  0.019252,  0.018912,  0.018636,  0.018012,  0.017343 };
  Float_t p2[kNb] = { -0.005042, -0.004926, -0.004867, -0.004650, -0.004424 };
  Float_t p3[kNb] = {  0.000195,  0.000189,  0.000195,  0.000182,  0.000169 };

  Int_t ib = ((Int_t) (10 * (b - 0.15)));
  ib       = TMath::Max(  0,ib);
  ib       = TMath::Min(kNb,ib);

  Float_t diff = p0[ib] 
               + p1[ib] * vd
               + p2[ib] * vd*vd
               + p3[ib] * vd*vd*vd;

  return diff;

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::GetDiffusionT(Float_t vd, Float_t b)
{
  //
  // Returns the transverse diffusion coefficient for a given drift 
  // velocity <vd> and a B-field <b> for Xe/CO2 (15%).
  // The values are according to a GARFIELD simulation.
  //

  const Int_t kNb = 5;
  Float_t p0[kNb] = {  0.009550,  0.009599,  0.009674,  0.009757,  0.009850 };
  Float_t p1[kNb] = {  0.006667,  0.006539,  0.006359,  0.006153,  0.005925 };
  Float_t p2[kNb] = { -0.000853, -0.000798, -0.000721, -0.000635, -0.000541 };
  Float_t p3[kNb] = {  0.000131,  0.000122,  0.000111,  0.000098,  0.000085 };

  Int_t ib = ((Int_t) (10 * (b - 0.15)));
  ib       = TMath::Max(  0,ib);
  ib       = TMath::Min(kNb,ib);

  Float_t diff = p0[ib] 
               + p1[ib] * vd
               + p2[ib] * vd*vd
               + p3[ib] * vd*vd*vd;

  return diff;

}

//_____________________________________________________________________________
Float_t AliTRDdigitizer::GetOmegaTau(Float_t vd, Float_t b)
{
  //
  // Returns omega*tau (tan(Lorentz-angle)) for a given drift velocity <vd> 
  // and a B-field <b> for Xe/CO2 (15%).
  // The values are according to a GARFIELD simulation.
  //

  const Int_t kNb = 5;
  Float_t p0[kNb] = {  0.004810,  0.007412,  0.010252,  0.013409,  0.016888 };
  Float_t p1[kNb] = {  0.054875,  0.081534,  0.107333,  0.131983,  0.155455 };
  Float_t p2[kNb] = { -0.008682, -0.012896, -0.016987, -0.020880, -0.024623 };
  Float_t p3[kNb] = {  0.000155,  0.000238,  0.000330,  0.000428,  0.000541 };

  Int_t ib = ((Int_t) (10 * (b - 0.15)));
  ib       = TMath::Max(  0,ib);
  ib       = TMath::Min(kNb,ib);

  Float_t alphaL = p0[ib] 
                 + p1[ib] * vd
                 + p2[ib] * vd*vd
                 + p3[ib] * vd*vd*vd;

  return TMath::Tan(alphaL);

}

