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
//  TRD parameter class                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TRandom.h>

#include "AliRun.h"
#include "AliMagF.h"

#include "AliTRDparameter.h"
#include "AliTRDgeometryFull.h"

ClassImp(AliTRDparameter)

//_____________________________________________________________________________
AliTRDparameter::AliTRDparameter():TNamed()
{
  //
  // AliTRDparameter default constructor
  //

  fGeo                = 0;
  fPRFsmp             = 0;
  fTRFsmp             = 0;
  fCTsmp              = 0;
  fGasGain            = 0.0;
  fNoise              = 0.0;
  fChipGain           = 0.0;
  fADCoutRange        = 0.0;
  fADCinRange         = 0.0;
  fADCthreshold       = 0;
  fADCbaseline        = 0;        
  fDiffusionOn        = 0;
  fDiffusionT         = 0.0;
  fDiffusionL         = 0.0;
  fElAttachOn         = 0;
  fElAttachProp       = 0.0;
  fExBOn              = 0;
  fOmegaTau           = 0.0;
  fPRFOn              = 0;
  fTRFOn              = 0;
  fCTOn               = 0;
  fTCOn               = 0;
  fDriftVelocity      = 0.0;
  fPadCoupling        = 0.0;
  fTimeCoupling       = 0.0;
  fTimeBinWidth       = 0.0;
  fField              = 0.0;
  fTiltingAngle       = 0.0;
  fPRFbin             = 0;
  fPRFlo              = 0.0;
  fPRFhi              = 0.0;
  fPRFwid             = 0.0;
  fPRFpad             = 0;
  fTRFbin             = 0;
  fTRFlo              = 0.0;
  fTRFhi              = 0.0;
  fTRFwid             = 0.0;
  fTCnexp             = 0;

  fLUTOn              = 0;  
  fLUT                = 0;
  fClusMaxThresh      = 0;
  fClusSigThresh      = 0;

}

//_____________________________________________________________________________
AliTRDparameter::AliTRDparameter(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDparameter constructor
  //

  fGeo                = new AliTRDgeometryFull();
  fPRFsmp             = 0;
  fTRFsmp             = 0;
  fCTsmp              = 0;
  fGasGain            = 0.0;
  fNoise              = 0.0;
  fChipGain           = 0.0;
  fADCoutRange        = 0.0;
  fADCinRange         = 0.0;
  fADCthreshold       = 0;
  fADCbaseline        = 0;        
  fDiffusionOn        = 0;
  fDiffusionT         = 0.0;
  fDiffusionL         = 0.0;
  fElAttachOn         = 0;
  fElAttachProp       = 0.0;
  fExBOn              = 0;
  fOmegaTau           = 0.0;
  fPRFOn              = 0;
  fTRFOn              = 0;
  fCTOn               = 0;
  fTCOn               = 0;
  fDriftVelocity      = 0.0;
  fPadCoupling        = 0.0;
  fTimeCoupling       = 0.0;
  fTimeBinWidth       = 0.0;
  fField              = 0.0;
  fTiltingAngle       = 0.0;
  fPRFbin             = 0;
  fPRFlo              = 0.0;
  fPRFhi              = 0.0;
  fPRFwid             = 0.0;
  fPRFpad             = 0;
  fTRFbin             = 0;
  fTRFlo              = 0.0;
  fTRFhi              = 0.0;
  fTRFwid             = 0.0;
  fTCnexp             = 0;

  fLUTOn              = 0;  
  fLUT                = 0;
  fClusMaxThresh      = 0;
  fClusSigThresh      = 0;

  Init();

}


//_____________________________________________________________________________
AliTRDparameter::AliTRDparameter(const AliTRDparameter &p):TNamed(p)
{
  //
  // AliTRDparameter copy constructor
  //

  ((AliTRDparameter &) p).Copy(*this);

}

///_____________________________________________________________________________
AliTRDparameter::~AliTRDparameter()
{
  //
  // AliTRDparameter destructor
  //

  if (fTRFsmp) {
    delete [] fTRFsmp;
    fTRFsmp = 0;
  }

  if (fPRFsmp) {
    delete [] fPRFsmp;
    fPRFsmp = 0;
  }

  if (fCTsmp) {
    delete [] fCTsmp;
    fCTsmp  = 0;
  }

  if (fLUT) {
    delete [] fLUT;
    fLUT    = 0;
  }

  if (fGeo) {
    delete fGeo;
    fGeo    = 0;
  }

}

//_____________________________________________________________________________
AliTRDparameter &AliTRDparameter::operator=(const AliTRDparameter &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDparameter &) p).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDparameter::Copy(TObject &p)
{
  //
  // Copy function
  //

  Int_t iBin;

  ((AliTRDparameter &) p).fGasGain            = fGasGain;
  ((AliTRDparameter &) p).fNoise              = fNoise;
  ((AliTRDparameter &) p).fChipGain           = fChipGain;
  ((AliTRDparameter &) p).fADCoutRange        = fADCoutRange;
  ((AliTRDparameter &) p).fADCinRange         = fADCinRange;
  ((AliTRDparameter &) p).fADCthreshold       = fADCthreshold;
  ((AliTRDparameter &) p).fADCbaseline        = fADCbaseline; 
  ((AliTRDparameter &) p).fDiffusionOn        = fDiffusionOn; 
  ((AliTRDparameter &) p).fDiffusionT         = fDiffusionT;
  ((AliTRDparameter &) p).fDiffusionL         = fDiffusionL;
  ((AliTRDparameter &) p).fElAttachOn         = fElAttachOn;
  ((AliTRDparameter &) p).fElAttachProp       = fElAttachProp;
  ((AliTRDparameter &) p).fExBOn              = fExBOn;
  ((AliTRDparameter &) p).fOmegaTau           = fOmegaTau;
  ((AliTRDparameter &) p).fLorentzFactor      = fLorentzFactor;
  ((AliTRDparameter &) p).fDriftVelocity      = fDriftVelocity;
  ((AliTRDparameter &) p).fPadCoupling        = fPadCoupling;
  ((AliTRDparameter &) p).fTimeCoupling       = fTimeCoupling;
  ((AliTRDparameter &) p).fTimeBinWidth       = fTimeBinWidth;
  ((AliTRDparameter &) p).fField              = fField;
  ((AliTRDparameter &) p).fPRFOn              = fPRFOn;
  ((AliTRDparameter &) p).fTRFOn              = fTRFOn;
  ((AliTRDparameter &) p).fCTOn               = fCTOn;
  ((AliTRDparameter &) p).fTCOn               = fTCOn;
  ((AliTRDparameter &) p).fTiltingAngle       = fTiltingAngle;
  ((AliTRDparameter &) p).fPRFbin             = fPRFbin;
  ((AliTRDparameter &) p).fPRFlo              = fPRFlo;
  ((AliTRDparameter &) p).fPRFhi              = fPRFhi;
  ((AliTRDparameter &) p).fPRFwid             = fPRFwid;
  ((AliTRDparameter &) p).fPRFpad             = fPRFpad;
  if (((AliTRDparameter &) p).fPRFsmp) delete [] ((AliTRDparameter &) p).fPRFsmp;
  ((AliTRDparameter &) p).fPRFsmp = new Float_t[fPRFbin];
  for (iBin = 0; iBin < fPRFbin; iBin++) {
    ((AliTRDparameter &) p).fPRFsmp[iBin] = fPRFsmp[iBin];
  }                                                                             
  ((AliTRDparameter &) p).fTRFbin             = fTRFbin;
  ((AliTRDparameter &) p).fTRFlo              = fTRFlo;
  ((AliTRDparameter &) p).fTRFhi              = fTRFhi;
  ((AliTRDparameter &) p).fTRFwid             = fTRFwid;
  if (((AliTRDparameter &) p).fTRFsmp) delete [] ((AliTRDparameter &) p).fTRFsmp;
  ((AliTRDparameter &) p).fTRFsmp = new Float_t[fTRFbin];
  for (iBin = 0; iBin < fTRFbin; iBin++) {
    ((AliTRDparameter &) p).fTRFsmp[iBin] = fTRFsmp[iBin];
  }                                      
  if (((AliTRDparameter &) p).fCTsmp)  delete [] ((AliTRDparameter &) p).fCTsmp;
  ((AliTRDparameter &) p).fCTsmp  = new Float_t[fTRFbin];
  for (iBin = 0; iBin < fTRFbin; iBin++) {
    ((AliTRDparameter &) p).fCTsmp[iBin]  = fCTsmp[iBin];
  }                                      
  ((AliTRDparameter &) p).fTCnexp             = fTCnexp;

  ((AliTRDparameter &) p).fLUTOn              = fLUTOn;
  ((AliTRDparameter &) p).fLUTbin             = fLUTbin;
  if (((AliTRDparameter &) p).fLUT)    delete [] ((AliTRDparameter &) p).fLUT;
  ((AliTRDparameter &) p).fLUT  = new Float_t[fLUTbin];
  for (iBin = 0; iBin < fLUTbin; iBin++) {
    ((AliTRDparameter &) p).fLUT[iBin]  = fLUT[iBin];
  }                                      
  ((AliTRDparameter &) p).fClusMaxThresh      = fClusMaxThresh;
  ((AliTRDparameter &) p).fClusSigThresh      = fClusSigThresh;

}

//_____________________________________________________________________________
void AliTRDparameter::Init()
{
  //
  // Initializes the parameter
  //
  // The maximum number of pads
  // and the position of pad 0,0,0
  //
  // chambers seen from the top:
  //     +----------------------------+
  //     |                            |
  //     |                            |      ^
  //     |                            |  rphi|
  //     |                            |      |
  //     |0                           |      |
  //     +----------------------------+      +------>
  //                                             z
  // chambers seen from the side:            ^
  //     +----------------------------+ drift|
  //     |0                           |      |
  //     |                            |      |
  //     +----------------------------+      +------>
  //                                             z
  //
  // IMPORTANT: time bin 0 is now the first one in the drift region
  // closest to the readout !!!
  //

  //
  // ----------------------------------------------------------------------------
  // The pad definition
  // ----------------------------------------------------------------------------
  //

  // The pad size in column direction (rphi-direction)
  SetColPadSize(0,0.65);
  SetColPadSize(1,0.68);
  SetColPadSize(2,0.71);
  SetColPadSize(3,0.74);
  SetColPadSize(4,0.77);
  SetColPadSize(5,0.80);

  // New pad size? Needs to be checked!
  //SetColPadSize(0,0.664);
  //SetColPadSize(1,0.695);
  //SetColPadSize(2,0.726);
  //SetColPadSize(3,0.756);
  //SetColPadSize(4,0.788);
  //SetColPadSize(5,0.818);

  // The pad row (z-direction)
  SetNRowPad();

  // The number of time bins. Default is 100 ns timbin size
  SetNTimeBin(15);

  // Additional time bins before and after the drift region.
  // Default is to only sample the drift region
  SetExpandTimeBin(0,0);

  //
  // ----------------------------------------------------------------------------
  // The digitization parameter
  // ----------------------------------------------------------------------------
  //

  // The default parameter for the digitization
  fGasGain        = 4500.;
  fChipGain       = 12.4;
  fNoise          = 1000.;
  fADCoutRange    = 1023.;          // 10-bit ADC
  fADCinRange     = 2000.;          // 2V input range
  fADCthreshold   = 1;
  fADCbaseline    = 0;

  // The drift velocity (cm / mus)
  fDriftVelocity  = 1.5;

  // Diffusion on
  fDiffusionOn    = 1;

  // E x B effects
  fExBOn          = 1;

  // Propability for electron attachment
  fElAttachOn     = 0;
  fElAttachProp   = 0.0;

  // The pad response function
  fPRFOn          = 1;

  // The time response function
  fTRFOn          = 1;

  // The cross talk
  fCTOn           = 1;

  // The tail cancelation
  fTCOn           = 1;
  
  // The number of exponentials
  fTCnexp         = 2;

  // The pad coupling factor (same number as for the TPC)
  fPadCoupling    = 0.5;

  // The time coupling factor (same number as for the TPC)
  fTimeCoupling   = 0.4;

  // The tilting angle for the readout pads
  SetTiltingAngle(2.0);

  // The magnetic field strength in Tesla
  fField           = 0.4;

  //
  // ----------------------------------------------------------------------------
  // The clusterization parameter
  // ----------------------------------------------------------------------------
  //

  // The default parameter for the clustering
  fClusMaxThresh = 3;
  fClusSigThresh = 1;

  // Use the LUT
  fLUTOn         = 1;  

  ReInit();

}

//_____________________________________________________________________________
void AliTRDparameter::ReInit()
{
  //
  // Reinitializes the parameter class after a change
  //

  // Calculate the time bin width in ns
  fTimeBinWidth   = fTimeBinSize / fDriftVelocity * 1000.0;

  // The range and the binwidth for the sampled TRF 
  fTRFbin = 100;
  // Start 0.2 mus before the signal
  fTRFlo  = -0.2 * fDriftVelocity;
  // End the maximum driftlength after the signal 
  fTRFhi  = AliTRDgeometry::DrThick() 
          + fTimeAfter * fTimeBinSize;
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

  // Create the sampled PRF
  SamplePRF();

  // Create the sampled TRF
  SampleTRF();

  // Create the LUT
  FillLUT();

}

//_____________________________________________________________________________
void AliTRDparameter::SetNRowPad(Int_t p, Int_t c, Int_t npad)
{
  //
  // Redefines the number of pads in raw direction for
  // a given plane and chamber number
  //

  for (Int_t isect = 0; isect < AliTRDgeometry::Nsect(); isect++) {

    fRowMax[p][c][isect] = npad;

    fRowPadSize[p][c][isect] = (fGeo->GetChamberLength(p,c) 
				- 2.* AliTRDgeometry::RpadW())
                             / ((Float_t) npad);

  }

}

//_____________________________________________________________________________
void AliTRDparameter::SetNRowPad()
{
  //
  // Defines the number of pads in row direction
  //

  Int_t isect;
  Int_t icham;
  Int_t iplan;

  Int_t rowMax[kNplan][kNcham] = { { 16, 16, 12, 16, 16 }
                                 , { 16, 16, 12, 16, 16 }
                                 , { 16, 16, 12, 16, 16 }
                                 , { 16, 16, 12, 16, 16 }
                                 , { 16, 16, 12, 16, 16 }
                                 , { 16, 16, 12, 16, 16 } };

  Float_t rpadW = AliTRDgeometry::RpadW();

  for (isect = 0; isect < kNsect; isect++) {
    for (icham = 0; icham < kNcham; icham++) {
      for (iplan = 0; iplan < kNplan; iplan++) {

        fRowMax[iplan][icham][isect]     = rowMax[iplan][icham];

        fRowPadSize[iplan][icham][isect] = (fGeo->GetChamberLength(iplan,icham) 
                                            - 2.*rpadW)
                                         / ((Float_t) rowMax[iplan][icham]);

        Float_t row0 = rpadW - fGeo->GetChamberLength(iplan,0)
  	                     - fGeo->GetChamberLength(iplan,1)
                             - fGeo->GetChamberLength(iplan,2) / 2.;
        for (Int_t ic = 0; ic < icham; ic++) {
          row0 += fGeo->GetChamberLength(iplan,ic);
        }
        fRow0[iplan][icham][isect]       = row0;

      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDparameter::SetColPadSize(Int_t p, Float_t s)
{
  //
  // Redefines the pad size in column direction
  //

  Float_t cpadW  = AliTRDgeometry::CpadW();

  fColPadSize[p] = s;
  fCol0[p]       = - fGeo->GetChamberWidth(p)/2. + cpadW;
  fColMax[p]     = ((Int_t) ((fGeo->GetChamberWidth(p) - 2.*cpadW) / s));

}

//_____________________________________________________________________________
void AliTRDparameter::SetNTimeBin(Int_t nbin)
{
  //
  // Redefines the number of time bins in the drift region.
  // The time bin width is defined by the length of the
  // drift region divided by <nbin>.
  //

  fTimeMax     = nbin;
  fTimeBinSize = AliTRDgeometry::DrThick() / ((Float_t) fTimeMax);
  for (Int_t iplan = 0; iplan < AliTRDgeometry::Nplan(); iplan++) {
    fTime0[iplan] = AliTRDgeometry::Rmin()
                  + AliTRDgeometry::CraHght()
                  + AliTRDgeometry::CdrHght()
                  + iplan * (AliTRDgeometry::Cheight() 
                           + AliTRDgeometry::Cspace());
  }

}

//_____________________________________________________________________________
Float_t AliTRDparameter::CrossTalk(Float_t time) const
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

//_____________________________________________________________________________
Int_t AliTRDparameter::Diffusion(Float_t driftlength, Float_t *xyz)
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
Int_t AliTRDparameter::ExB(Float_t driftlength, Float_t *xyz) const
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
Int_t AliTRDparameter::PadResponse(Float_t signal, Float_t dist
                                 , Int_t plane, Float_t *pad) const
{
  //
  // Applies the pad response
  //

  const Int_t kNplan = AliTRDgeometry::kNplan;

  Int_t iBin  = ((Int_t) (( - dist - fPRFlo) / fPRFwid));
  Int_t iOff  = plane * fPRFbin;

  Int_t iBin0 = iBin - fPRFpad + iOff;
  Int_t iBin1 = iBin           + iOff;
  Int_t iBin2 = iBin + fPRFpad + iOff;

  pad[0] = 0.0;
  pad[1] = 0.0;
  pad[2] = 0.0;
  if ((iBin1 >= 0) && (iBin1 < (fPRFbin*kNplan))) {

    if (iBin0 >= 0) {
      pad[0] = signal * fPRFsmp[iBin0];
    }
    pad[1] = signal * fPRFsmp[iBin1];
    if (iBin2 < (fPRFbin*kNplan)) {
      pad[2] = signal * fPRFsmp[iBin2];
    }

    return 1;

  }
  else {

    return 0;

  }

}

//_____________________________________________________________________________
Float_t AliTRDparameter::TimeResponse(Float_t time) const
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
Float_t AliTRDparameter::Col0Tilted(Float_t col0, Float_t rowOffset
                                  , Int_t plane)
{
  //
  // Calculates col0 for tilted pads
  //

  Float_t diff = fTiltingAngle * rowOffset;
  return (col0 + TMath::Power(-1.0,plane) * diff);

}

//_____________________________________________________________________________
void AliTRDparameter::SampleTRF()
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
void AliTRDparameter::SamplePRF()
{
  //
  // Samples the pad response function
  //

  const Int_t kNplan  = AliTRDgeometry::kNplan;
  const Int_t kPRFbin = 61;

  Float_t prf[kNplan][kPRFbin] = { { 0.018570, 0.022270, 0.026710, 0.032010
                                   , 0.038350, 0.045920, 0.054930, 0.065650
                                   , 0.078370, 0.093420, 0.111150, 0.131940
                                   , 0.156160, 0.184160, 0.216220, 0.252470
                                   , 0.292860, 0.337030, 0.384330, 0.433750
                                   , 0.484010, 0.533630, 0.581150, 0.625200
                                   , 0.664710, 0.698860, 0.727130, 0.749230
                                   , 0.765050, 0.774540, 0.777700, 0.774540
                                   , 0.765050, 0.749230, 0.727130, 0.698860
                                   , 0.664710, 0.625200, 0.581150, 0.533630
                                   , 0.484010, 0.433750, 0.384330, 0.337030
                                   , 0.292860, 0.252470, 0.216220, 0.184160
                                   , 0.156160, 0.131940, 0.111150, 0.093420
                                   , 0.078370, 0.065650, 0.054930, 0.045920
                                   , 0.038350, 0.032010, 0.026710, 0.022270
				   , 0.018570                               }
                                 , { 0.015730, 0.019040, 0.023030, 0.027840
                                   , 0.033650, 0.040650, 0.049060, 0.059160
                                   , 0.071260, 0.085710, 0.102910, 0.123270
                                   , 0.147240, 0.175220, 0.207590, 0.244540
                                   , 0.286090, 0.331920, 0.381350, 0.433290
                                   , 0.486290, 0.538710, 0.588870, 0.635280
                                   , 0.676760, 0.712460, 0.741890, 0.764810
                                   , 0.781150, 0.790930, 0.794180, 0.790930
                                   , 0.781150, 0.764810, 0.741890, 0.712460
                                   , 0.676760, 0.635280, 0.588870, 0.538710
                                   , 0.486290, 0.433290, 0.381350, 0.331920
                                   , 0.286090, 0.244540, 0.207590, 0.175220
                                   , 0.147240, 0.123270, 0.102910, 0.085710
                                   , 0.071260, 0.059160, 0.049060, 0.040650
                                   , 0.033650, 0.027840, 0.023030, 0.019040
				   , 0.015730                               }
				 , { 0.013330, 0.016270, 0.019850, 0.024210
                                   , 0.029510, 0.035960, 0.043790, 0.053280
                                   , 0.064740, 0.078580, 0.095190, 0.115070
                                   , 0.138700, 0.166570, 0.199120, 0.236660
                                   , 0.279260, 0.326660, 0.378140, 0.432540
                                   , 0.488260, 0.543440, 0.596200, 0.644900
                                   , 0.688240, 0.725380, 0.755840, 0.779470
                                   , 0.796260, 0.806280, 0.809610, 0.806280
                                   , 0.796260, 0.779470, 0.755840, 0.725380
                                   , 0.688240, 0.644900, 0.596200, 0.543440
                                   , 0.488260, 0.432540, 0.378140, 0.326660
                                   , 0.279260, 0.236660, 0.199120, 0.166570
                                   , 0.138700, 0.115070, 0.095190, 0.078580
                                   , 0.064740, 0.053280, 0.043790, 0.035960
                                   , 0.029510, 0.024210, 0.019850, 0.016270
				   , 0.013330                               }
                                 , { 0.011280, 0.013890, 0.017090, 0.021030
                                   , 0.025870, 0.031800, 0.039060, 0.047940
                                   , 0.058790, 0.071980, 0.087990, 0.107330
                                   , 0.130550, 0.158220, 0.190850, 0.228870
                                   , 0.272410, 0.321270, 0.374740, 0.431560
                                   , 0.489960, 0.547870, 0.603180, 0.654080
                                   , 0.699190, 0.737640, 0.769030, 0.793260
                                   , 0.810410, 0.820620, 0.824010, 0.820620
                                   , 0.810410, 0.793260, 0.769030, 0.737640
                                   , 0.699190, 0.654080, 0.603180, 0.547870
                                   , 0.489960, 0.431560, 0.374740, 0.321270
                                   , 0.272410, 0.228870, 0.190850, 0.158220
                                   , 0.130550, 0.107330, 0.087990, 0.071980
                                   , 0.058790, 0.047940, 0.039060, 0.031800
                                   , 0.025870, 0.021030, 0.017090, 0.013890
				   , 0.011280                               }
                                 , { 0.009550, 0.011860, 0.014720, 0.018270
                                   , 0.022660, 0.028100, 0.034820, 0.043120
                                   , 0.053340, 0.065900, 0.081280, 0.100040
                                   , 0.122800, 0.150180, 0.182800, 0.221170
                                   , 0.265550, 0.315790, 0.371180, 0.430370
                                   , 0.491430, 0.552030, 0.609840, 0.662860
                                   , 0.709630, 0.749290, 0.781490, 0.806220
                                   , 0.823650, 0.834000, 0.837430, 0.834000
                                   , 0.823650, 0.806220, 0.781490, 0.749290
                                   , 0.709630, 0.662860, 0.609840, 0.552030
                                   , 0.491430, 0.430370, 0.371180, 0.315790
                                   , 0.265550, 0.221170, 0.182800, 0.150180
                                   , 0.122800, 0.100040, 0.081280, 0.065900
                                   , 0.053340, 0.043120, 0.034820, 0.028100
                                   , 0.022660, 0.018270, 0.014720, 0.011860
				   , 0.009550                               }
                                 , { 0.008080, 0.010120, 0.012670, 0.015860
                                   , 0.019840, 0.024820, 0.031030, 0.038760
                                   , 0.048370, 0.060300, 0.075040, 0.093200
                                   , 0.115430, 0.142450, 0.174980, 0.213610
                                   , 0.258720, 0.310250, 0.367480, 0.429010
                                   , 0.492690, 0.555950, 0.616210, 0.671280
                                   , 0.719600, 0.760350, 0.793250, 0.818380
                                   , 0.836020, 0.846460, 0.849920, 0.846460
                                   , 0.836020, 0.818380, 0.793250, 0.760350
                                   , 0.719600, 0.671280, 0.616210, 0.555950
                                   , 0.492690, 0.429010, 0.367480, 0.310250
                                   , 0.258720, 0.213610, 0.174980, 0.142450
                                   , 0.115430, 0.093200, 0.075040, 0.060300
                                   , 0.048370, 0.038760, 0.031030, 0.024820
                                   , 0.019840, 0.015860, 0.012670, 0.010120
				   , 0.008080                               } };

  // More sampling precision with linear interpolation
  fPRFlo  = -1.5;
  fPRFhi  =  1.5;
  Float_t pad[kPRFbin];
  Int_t   sPRFbin = kPRFbin;  
  Float_t sPRFwid = (fPRFhi - fPRFlo) / ((Float_t) sPRFbin);
  for (Int_t iPad = 0; iPad < sPRFbin; iPad++) {
    pad[iPad] = ((Float_t) iPad + 0.5) * sPRFwid + fPRFlo;
  }
  fPRFbin = 500;  
  fPRFwid = (fPRFhi - fPRFlo) / ((Float_t) fPRFbin);
  fPRFpad = ((Int_t) (1.0 / fPRFwid));

  if (fPRFsmp) delete [] fPRFsmp;
  fPRFsmp = new Float_t[kNplan*fPRFbin];

  Int_t   ipos1;
  Int_t   ipos2;
  Float_t diff;

  for (Int_t iPla = 0; iPla < kNplan; iPla++) {

    for (Int_t iBin = 0; iBin < fPRFbin; iBin++) {

      Float_t bin = (((Float_t) iBin) + 0.5) * fPRFwid + fPRFlo;
      ipos1 = ipos2 = 0;
      diff  = 0;
      do {
        diff = bin - pad[ipos2++];
      } while ((diff > 0) && (ipos2 < kPRFbin));
      if      (ipos2 == kPRFbin) {
        fPRFsmp[iPla*fPRFbin+iBin] = prf[iPla][ipos2-1];
      }
      else if (ipos2 == 1) {
        fPRFsmp[iPla*fPRFbin+iBin] = prf[iPla][ipos2-1];
      }
      else {
        ipos2--;
        if (ipos2 >= kPRFbin) ipos2 = kPRFbin - 1;
        ipos1 = ipos2 - 1;
        fPRFsmp[iPla*fPRFbin+iBin] = prf[iPla][ipos2] 
                                   + diff * (prf[iPla][ipos2] - prf[iPla][ipos1]) 
                                          / sPRFwid;
      }

    }
  } 

}

//_____________________________________________________________________________
void AliTRDparameter::FillLUT()
{
  //
  // Create the LUT
  //

  const Int_t kNplan = AliTRDgeometry::kNplan;
  const Int_t kNlut  = 128;

  fLUTbin = kNplan * kNlut;

  // The lookup table from Bogdan
  Float_t lut[kNplan][kNlut] = {  
    {
      0.0070, 0.0150, 0.0224, 0.0298, 0.0374, 0.0454, 0.0533, 0.0611, 
      0.0684, 0.0755, 0.0827, 0.0900, 0.0975, 0.1049, 0.1120, 0.1187, 
      0.1253, 0.1318, 0.1385, 0.1453, 0.1519, 0.1584, 0.1646, 0.1704, 
      0.1762, 0.1821, 0.1879, 0.1938, 0.1996, 0.2053, 0.2108, 0.2160, 
      0.2210, 0.2260, 0.2310, 0.2361, 0.2411, 0.2461, 0.2509, 0.2557, 
      0.2602, 0.2646, 0.2689, 0.2732, 0.2774, 0.2816, 0.2859, 0.2901, 
      0.2942, 0.2983, 0.3022, 0.3061, 0.3099, 0.3136, 0.3172, 0.3207, 
      0.3242, 0.3278, 0.3312, 0.3347, 0.3382, 0.3416, 0.3450, 0.3483, 
      0.3515, 0.3547, 0.3579, 0.3609, 0.3639, 0.3669, 0.3698, 0.3727, 
      0.3756, 0.3785, 0.3813, 0.3842, 0.3870, 0.3898, 0.3926, 0.3952, 
      0.3979, 0.4005, 0.4032, 0.4057, 0.4082, 0.4108, 0.4132, 0.4157, 
      0.4181, 0.4205, 0.4228, 0.4252, 0.4275, 0.4299, 0.4322, 0.4345, 
      0.4367, 0.4390, 0.4412, 0.4434, 0.4456, 0.4478, 0.4499, 0.4520, 
      0.4541, 0.4562, 0.4583, 0.4603, 0.4623, 0.4643, 0.4663, 0.4683, 
      0.4702, 0.4722, 0.4741, 0.4758, 0.4774, 0.4790, 0.4805, 0.4824, 
      0.4844, 0.4863, 0.4883, 0.4902, 0.4921, 0.4940, 0.4959, 0.4978 
    },
    {
      0.0072, 0.0156, 0.0235, 0.0313, 0.0394, 0.0478, 0.0561, 0.0642, 
      0.0718, 0.0792, 0.0868, 0.0947, 0.1025, 0.1101, 0.1172, 0.1241, 
      0.1309, 0.1378, 0.1449, 0.1518, 0.1586, 0.1650, 0.1710, 0.1770, 
      0.1830, 0.1891, 0.1952, 0.2011, 0.2070, 0.2125, 0.2177, 0.2229, 
      0.2280, 0.2332, 0.2383, 0.2435, 0.2484, 0.2533, 0.2581, 0.2627, 
      0.2670, 0.2714, 0.2757, 0.2799, 0.2842, 0.2884, 0.2927, 0.2968, 
      0.3008, 0.3048, 0.3086, 0.3123, 0.3159, 0.3195, 0.3231, 0.3266, 
      0.3301, 0.3335, 0.3370, 0.3404, 0.3438, 0.3471, 0.3504, 0.3536, 
      0.3567, 0.3598, 0.3628, 0.3657, 0.3686, 0.3715, 0.3744, 0.3772, 
      0.3800, 0.3828, 0.3856, 0.3884, 0.3911, 0.3938, 0.3965, 0.3991, 
      0.4016, 0.4042, 0.4067, 0.4092, 0.4116, 0.4140, 0.4164, 0.4187, 
      0.4211, 0.4234, 0.4257, 0.4280, 0.4302, 0.4325, 0.4347, 0.4369, 
      0.4391, 0.4413, 0.4434, 0.4456, 0.4477, 0.4497, 0.4518, 0.4538, 
      0.4558, 0.4578, 0.4598, 0.4618, 0.4637, 0.4656, 0.4675, 0.4694, 
      0.4713, 0.4732, 0.4750, 0.4766, 0.4781, 0.4797, 0.4813, 0.4832, 
      0.4851, 0.4870, 0.4888, 0.4906, 0.4925, 0.4942, 0.4960, 0.4978
    },
    {
      0.0075, 0.0163, 0.0246, 0.0328, 0.0415, 0.0504, 0.0592, 0.0674, 
      0.0753, 0.0832, 0.0914, 0.0996, 0.1077, 0.1154, 0.1225, 0.1296, 
      0.1369, 0.1442, 0.1515, 0.1585, 0.1652, 0.1714, 0.1776, 0.1839, 
      0.1902, 0.1965, 0.2025, 0.2085, 0.2141, 0.2194, 0.2247, 0.2299, 
      0.2352, 0.2405, 0.2457, 0.2507, 0.2557, 0.2604, 0.2649, 0.2693, 
      0.2737, 0.2780, 0.2823, 0.2867, 0.2909, 0.2951, 0.2992, 0.3033, 
      0.3072, 0.3110, 0.3146, 0.3182, 0.3218, 0.3253, 0.3288, 0.3323, 
      0.3357, 0.3392, 0.3426, 0.3459, 0.3492, 0.3524, 0.3555, 0.3586, 
      0.3616, 0.3645, 0.3674, 0.3703, 0.3731, 0.3759, 0.3787, 0.3815, 
      0.3843, 0.3870, 0.3897, 0.3925, 0.3950, 0.3976, 0.4002, 0.4027, 
      0.4052, 0.4076, 0.4101, 0.4124, 0.4148, 0.4171, 0.4194, 0.4217, 
      0.4239, 0.4262, 0.4284, 0.4306, 0.4328, 0.4350, 0.4371, 0.4393, 
      0.4414, 0.4435, 0.4455, 0.4476, 0.4496, 0.4516, 0.4536, 0.4555, 
      0.4575, 0.4594, 0.4613, 0.4632, 0.4650, 0.4669, 0.4687, 0.4705, 
      0.4723, 0.4741, 0.4758, 0.4773, 0.4789, 0.4804, 0.4821, 0.4839, 
      0.4857, 0.4875, 0.4893, 0.4910, 0.4928, 0.4945, 0.4961, 0.4978
    },
    {
      0.0078, 0.0171, 0.0258, 0.0345, 0.0438, 0.0532, 0.0624, 0.0708, 
      0.0791, 0.0875, 0.0962, 0.1048, 0.1130, 0.1206, 0.1281, 0.1356, 
      0.1432, 0.1508, 0.1582, 0.1651, 0.1716, 0.1780, 0.1845, 0.1910, 
      0.1975, 0.2038, 0.2099, 0.2155, 0.2210, 0.2263, 0.2317, 0.2371, 
      0.2425, 0.2477, 0.2528, 0.2578, 0.2626, 0.2671, 0.2715, 0.2759, 
      0.2803, 0.2846, 0.2890, 0.2933, 0.2975, 0.3016, 0.3056, 0.3095, 
      0.3132, 0.3168, 0.3204, 0.3239, 0.3274, 0.3309, 0.3344, 0.3378, 
      0.3412, 0.3446, 0.3479, 0.3511, 0.3543, 0.3574, 0.3603, 0.3633, 
      0.3662, 0.3690, 0.3718, 0.3747, 0.3774, 0.3802, 0.3829, 0.3857, 
      0.3883, 0.3910, 0.3936, 0.3962, 0.3987, 0.4012, 0.4037, 0.4061, 
      0.4085, 0.4109, 0.4132, 0.4155, 0.4177, 0.4200, 0.4222, 0.4244, 
      0.4266, 0.4288, 0.4309, 0.4331, 0.4352, 0.4373, 0.4394, 0.4414, 
      0.4435, 0.4455, 0.4475, 0.4494, 0.4514, 0.4533, 0.4552, 0.4571, 
      0.4590, 0.4608, 0.4626, 0.4645, 0.4662, 0.4680, 0.4698, 0.4715, 
      0.4733, 0.4750, 0.4766, 0.4781, 0.4796, 0.4812, 0.4829, 0.4846, 
      0.4863, 0.4880, 0.4897, 0.4914, 0.4930, 0.4946, 0.4963, 0.4979
    },
    {
      0.0081, 0.0178, 0.0270, 0.0364, 0.0463, 0.0562, 0.0656, 0.0744, 
      0.0831, 0.0921, 0.1013, 0.1102, 0.1183, 0.1261, 0.1339, 0.1419, 
      0.1499, 0.1576, 0.1648, 0.1715, 0.1782, 0.1849, 0.1917, 0.1984, 
      0.2048, 0.2110, 0.2167, 0.2223, 0.2278, 0.2333, 0.2389, 0.2444, 
      0.2497, 0.2548, 0.2598, 0.2645, 0.2691, 0.2735, 0.2780, 0.2824, 
      0.2868, 0.2912, 0.2955, 0.2997, 0.3038, 0.3078, 0.3116, 0.3152, 
      0.3188, 0.3224, 0.3259, 0.3294, 0.3329, 0.3364, 0.3398, 0.3432, 
      0.3465, 0.3497, 0.3529, 0.3561, 0.3591, 0.3620, 0.3649, 0.3677, 
      0.3705, 0.3733, 0.3761, 0.3788, 0.3816, 0.3843, 0.3869, 0.3896, 
      0.3922, 0.3948, 0.3973, 0.3998, 0.4022, 0.4047, 0.4070, 0.4094, 
      0.4117, 0.4139, 0.4162, 0.4184, 0.4206, 0.4227, 0.4249, 0.4270, 
      0.4291, 0.4313, 0.4334, 0.4354, 0.4375, 0.4395, 0.4415, 0.4435, 
      0.4455, 0.4474, 0.4493, 0.4512, 0.4531, 0.4550, 0.4568, 0.4586, 
      0.4604, 0.4622, 0.4639, 0.4657, 0.4674, 0.4691, 0.4708, 0.4725, 
      0.4742, 0.4758, 0.4773, 0.4788, 0.4803, 0.4819, 0.4836, 0.4852, 
      0.4869, 0.4885, 0.4901, 0.4917, 0.4933, 0.4948, 0.4964, 0.4979
    },
    {
      0.0085, 0.0189, 0.0288, 0.0389, 0.0497, 0.0603, 0.0699, 0.0792, 
      0.0887, 0.0985, 0.1082, 0.1170, 0.1253, 0.1336, 0.1421, 0.1505, 
      0.1587, 0.1662, 0.1733, 0.1803, 0.1874, 0.1945, 0.2014, 0.2081, 
      0.2143, 0.2201, 0.2259, 0.2316, 0.2374, 0.2431, 0.2487, 0.2541, 
      0.2593, 0.2642, 0.2689, 0.2735, 0.2781, 0.2826, 0.2872, 0.2917, 
      0.2961, 0.3003, 0.3045, 0.3086, 0.3125, 0.3162, 0.3198, 0.3235, 
      0.3270, 0.3306, 0.3342, 0.3377, 0.3411, 0.3446, 0.3479, 0.3511, 
      0.3543, 0.3575, 0.3605, 0.3634, 0.3663, 0.3691, 0.3720, 0.3748, 
      0.3775, 0.3803, 0.3830, 0.3857, 0.3884, 0.3911, 0.3937, 0.3962, 
      0.3987, 0.4012, 0.4036, 0.4060, 0.4084, 0.4107, 0.4129, 0.4152, 
      0.4174, 0.4196, 0.4218, 0.4239, 0.4261, 0.4282, 0.4303, 0.4324, 
      0.4344, 0.4365, 0.4385, 0.4405, 0.4425, 0.4445, 0.4464, 0.4483, 
      0.4502, 0.4521, 0.4539, 0.4558, 0.4576, 0.4593, 0.4611, 0.4629, 
      0.4646, 0.4663, 0.4680, 0.4697, 0.4714, 0.4730, 0.4747, 0.4759, 
      0.4769, 0.4780, 0.4790, 0.4800, 0.4811, 0.4827, 0.4843, 0.4859, 
      0.4874, 0.4889, 0.4905, 0.4920, 0.4935, 0.4950, 0.4965, 0.4979
    }
  }; 

  if (fLUT) delete [] fLUT;
  fLUT = new Float_t[fLUTbin];

  for (Int_t iplan = 0; iplan < kNplan; iplan++) {
    for (Int_t ilut  = 0; ilut  < kNlut; ilut++) {
      fLUT[iplan*kNlut+ilut] = lut[iplan][ilut];
    }
  }

}

//_____________________________________________________________________________
void AliTRDparameter::SetTiltingAngle(Float_t v)
{
  //
  // Set the tilting angle for the readout pads
  //

  fTiltingAngle = TMath::Tan(TMath::Pi()/180.0 * v);

}

//_____________________________________________________________________________
Float_t AliTRDparameter::GetTiltingAngle() const
{
  //
  // Get the tilting angle for the readout pads
  //

  return 180.0 / TMath::Pi() * TMath::ATan(fTiltingAngle);

}

//_____________________________________________________________________________
Float_t AliTRDparameter::GetDiffusionL(Float_t vd, Float_t b)
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
Float_t AliTRDparameter::GetDiffusionT(Float_t vd, Float_t b)
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
Float_t AliTRDparameter::GetOmegaTau(Float_t vd, Float_t b)
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

//_____________________________________________________________________________
Double_t AliTRDparameter::LUTposition(Int_t iplane, Double_t ampL 
                                                  , Double_t ampC
                                                  , Double_t ampR) const
{
  //
  // Calculates the cluster position using the lookup table.
  // Method provided by Bogdan Vulpescu.
  //

  const Int_t kNplan = AliTRDgeometry::kNplan;
  const Int_t kNlut  = 128;

  Double_t pos;
  Double_t x = 0.0;
  Double_t xmin;
  Double_t xmax;
  Double_t xwid;

  Int_t    side = 0;
  Int_t    ix;

  Double_t xMin[kNplan] = { 0.006492, 0.006377, 0.006258, 0.006144, 0.006030, 0.005980 };
  Double_t xMax[kNplan] = { 0.960351, 0.965870, 0.970445, 0.974352, 0.977667, 0.996101 };

  if      (ampL > ampR) {
    x    = (ampL - ampR) / ampC;
    side = -1;
  } 
  else if (ampL < ampR) {
    x    = (ampR - ampL) / ampC;
    side = +1;
  }

  if (ampL != ampR) {

    xmin = xMin[iplane] + 0.000005;
    xmax = xMax[iplane] - 0.000005;
    xwid = (xmax - xmin) / 127.0;

    if      (x < xmin) {
      pos = 0.0000;
    } 
    else if (x > xmax) {
      pos = side * 0.5000;
    } 
    else {
      ix  = (Int_t) ((x - xmin) / xwid);
      pos = side * fLUT[iplane*kNlut+ix];
    }
       
  } 
  else {

    pos = 0.0;

  }

  return pos;

}
