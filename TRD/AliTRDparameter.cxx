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

  fTimeStructOn       = 0;
  fTimeStruct         = 0;
  fStaggeringOn       = 0;
  fAnodeWireOffset    = 0.0;

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

  fTimeStructOn       = 0;
  fTimeStruct         = 0;
  fStaggeringOn       = 0;
  fAnodeWireOffset    = 0.0;

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

  if (fTimeStruct) {
    delete [] fTimeStruct;
    fTimeStruct = 0;
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

  ((AliTRDparameter &) p).fAnodeWireOffset    = fAnodeWireOffset;
  ((AliTRDparameter &) p).fStaggeringOn       = fStaggeringOn;
  ((AliTRDparameter &) p).fTimeStructOn       = fTimeStructOn;
  if (((AliTRDparameter &) p).fTimeStruct) 
    delete [] ((AliTRDparameter &) p).fTimeStruct;
  ((AliTRDparameter &) p).fTimeStruct = new Float_t[38*26];
  for (Int_t i = 0; i < 38*26; i++) {
    ((AliTRDparameter &) p).fTimeStruct[i] = fTimeStruct[i];
  }                                      

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
  //SetColPadSize(0,0.65);
  //SetColPadSize(1,0.68);
  //SetColPadSize(2,0.71);
  //SetColPadSize(3,0.74);
  //SetColPadSize(4,0.77);
  //SetColPadSize(5,0.80);

  SetColPadSize(0,0.664);
  SetColPadSize(1,0.695);
  SetColPadSize(2,0.726);
  SetColPadSize(3,0.756);
  SetColPadSize(4,0.788);
  SetColPadSize(5,0.818);

  // The pad row (z-direction)
  SetNRowPad();

  // The number of time bins. Default is 100 ns timbin size
  SetNTimeBin(20);

  // Additional time bins before and after the drift region.
  // Default is to only sample the drift region
  SetExpandTimeBin(0,0);

  //
  // ----------------------------------------------------------------------------
  // The digitization parameter
  // ----------------------------------------------------------------------------
  //

  // The default parameter for the digitization
  fGasGain        = 4000.;
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

  // The pad coupling factor
  //fPadCoupling    = 0.3;
  // Use 0.49 instead which reproduces better the test beam
  // data, even tough it is not understood why.
  fPadCoupling    = 0.49;

  // The time coupling factor (same number as for the TPC)
  fTimeCoupling   = 0.4;

  // Drift time non-isochronity on
  fTimeStructOn   = 1;

  // Wire planes staggered
  fStaggeringOn   = 1;

  // Distance of first Anode wire from first pad edge
  fAnodeWireOffset = 0.25;

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

  // Create the sampled timing structure
  SampleTimeStruct();

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
Int_t AliTRDparameter::TimeStruct(Float_t time, Float_t z, Float_t *xyz) const
{
  //
  // Applies the time structure of the drift cells (by C.Lippmann).
  // The drift time of electrons to the anode wires depends on the
  // distance to the wire (z) and on the position in the drift region.
  //
  // input : time (radial distance from pad cathode) in [cm] and z (distance
  // from anode wire parallel to cathode planes) in [cm]
  //

  if ( z > 0.24 )                   z = 0.24;
  if ( !fStaggeringOn && z < 0.01 ) z = 0.01;

  Int_t r1 = (Int_t)(10*time);
  Int_t r2 = r1+1;
  Int_t z1 = (Int_t)(100*z);
  Int_t z2 = z1+1;

  if (r1<0 || r1>37 || z1<0 || z1>25) {
    printf("<AliTRDparameter::TimeStruct> Warning. Indices out of range: ");
    printf("time=%.2f, z=%.2f, r1=%d, z1=%d\n",time,z,r1,z1);
    return kFALSE;
    if (r1<0)  r1 = 0;
    if (r1>37) r1 = 37;
    if (z1<0)  z1 = 0;
    if (z1>25) z1 = 25;
  }

  // 2D Interpolation:
  Float_t y11, y12, y21, y22, y1, y2;
  y11 = fTimeStruct[r1+38*z1];
  y22 = (r2 <= 37 && z2 <= 25) ? fTimeStruct[r2+38*z2] : fTimeStruct[37+38*25];
  y12 = (z2 <= 25)             ? fTimeStruct[r1+38*z2] : fTimeStruct[r1+38*25];
  y21 = (r2 <= 37)             ? fTimeStruct[r2+38*z1] : fTimeStruct[37+38*z1];

  y1  = (y21-y11)*10*time + y11 - (y21-y11)*r1;
  y2  = (y22-y12)*10*time + y12 - (y22-y12)*r1;

  Float_t AmTh   = AliTRDgeometry::AmThick() / 2.;
  Float_t tdrift = 0.0;                       // drift time (GARFIELD)
  if (TMath::Abs(time-AmTh)>0.005 || z>0.005)
    tdrift = (y2-y1)*100*z + y1 - (y2-y1)*z1;
  if (time < AmTh) tdrift *= -1.;
  Float_t xdrift = tdrift*(AmTh + AliTRDgeometry::DrThick())/fTimeStruct[37] + AmTh;
  Float_t offset = time - xdrift;

  xyz[0] = xyz[0] + offset;
  xyz[1] = xyz[1];
  xyz[2] = xyz[2];

  return 1;

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
  return (col0 + TMath::Power(-1.0,(plane+1)) * diff);

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
void AliTRDparameter::SampleTimeStruct()
{
  //
  // Samples the timing structure of a drift cell
  // Drift Time data calculated with Garfield (by C.Lippmann)
  //

  const Int_t ktimebin  = 38;
  const Int_t kZbin = 26;

  Float_t stagg[ktimebin][kZbin] =
    {{9.1221e-02, 9.1269e-02, 9.1426e-02, 9.1678e-02, 9.2041e-02, 9.2500e-02,
      9.3081e-02, 9.3760e-02, 9.4540e-02, 9.5434e-02, 9.6436e-02, 9.7554e-02,
      9.8796e-02, 1.0016e-01, 1.0167e-01, 1.0336e-01, 1.0528e-01, 1.0748e-01,
      1.1012e-01, 1.1348e-01, 1.1790e-01, 1.2395e-01, 1.3317e-01, 1.5036e-01,
      1.9348e-01, 5.7191e-01},
     {6.5708e-02, 6.5787e-02, 6.5954e-02, 6.6237e-02, 6.6646e-02, 6.7170e-02,
      6.7804e-02, 6.8554e-02, 6.9424e-02, 7.0392e-02, 7.1492e-02, 7.2693e-02,
      7.4006e-02, 7.5459e-02, 7.7042e-02, 7.8795e-02, 8.0786e-02, 8.3081e-02,
      8.5798e-02, 8.9290e-02, 9.3841e-02, 1.0007e-01, 1.0959e-01, 1.2741e-01,
      1.7152e-01, 5.4720e-01},
     {3.9602e-02, 3.9699e-02, 3.9969e-02, 4.0417e-02, 4.1015e-02, 4.1802e-02,
      4.2726e-02, 4.3804e-02, 4.4992e-02, 4.6309e-02, 4.7720e-02, 4.9244e-02,
      5.0874e-02, 5.2592e-02, 5.4437e-02, 5.6461e-02, 5.8704e-02, 6.1340e-02,
      6.4533e-02, 6.8509e-02, 7.3620e-02, 8.0737e-02, 9.1920e-02, 1.1276e-01,
      1.6078e-01, 5.2455e-01},
     {1.1236e-02, 1.1514e-02, 1.2355e-02, 1.3666e-02, 1.5354e-02, 1.7292e-02,
      1.9469e-02, 2.1756e-02, 2.4118e-02, 2.6505e-02, 2.8902e-02, 3.1273e-02,
      3.3618e-02, 3.5960e-02, 3.8367e-02, 4.0957e-02, 4.3859e-02, 4.7286e-02,
      5.1516e-02, 5.7160e-02, 6.5086e-02, 7.6740e-02, 9.4905e-02, 1.2426e-01,
      1.8267e-01, 4.8444e-01},
     {1.1285e-02, 1.1581e-02, 1.2430e-02, 1.3739e-02, 1.5383e-02, 1.7329e-02,
      1.9477e-02, 2.1732e-02, 2.4063e-02, 2.6412e-02, 2.8746e-02, 3.1073e-02,
      3.3329e-02, 3.5733e-02, 3.8146e-02, 4.0847e-02, 4.3908e-02, 4.7650e-02,
      5.2475e-02, 5.9069e-02, 6.8420e-02, 8.2181e-02, 1.0199e-01, 1.3525e-01,
      1.9150e-01, 4.6381e-01},
     {3.9345e-02, 3.9430e-02, 3.9698e-02, 4.0124e-02, 4.0714e-02, 4.1462e-02,
      4.2325e-02, 4.3327e-02, 4.4476e-02, 4.5727e-02, 4.7075e-02, 4.8532e-02,
      5.0118e-02, 5.1838e-02, 5.3758e-02, 5.5958e-02, 5.8528e-02, 6.1589e-02,
      6.5285e-02, 6.9903e-02, 7.6022e-02, 8.4799e-02, 9.8707e-02, 1.2316e-01,
      1.7495e-01, 5.2371e-01},
     {6.4837e-02, 6.4889e-02, 6.5062e-02, 6.5329e-02, 6.5688e-02, 6.6204e-02,
      6.6770e-02, 6.7494e-02, 6.8306e-02, 6.9234e-02, 7.0297e-02, 7.1484e-02,
      7.2842e-02, 7.4318e-02, 7.6041e-02, 7.8040e-02, 8.0428e-02, 8.3315e-02,
      8.6827e-02, 9.1171e-02, 9.6759e-02, 1.0445e-01, 1.1651e-01, 1.3841e-01,
      1.8746e-01, 5.4644e-01},
     {9.0789e-02, 9.0824e-02, 9.0908e-02, 9.1065e-02, 9.1268e-02, 9.1522e-02,
      9.1838e-02, 9.2207e-02, 9.2643e-02, 9.3230e-02, 9.4073e-02, 9.5341e-02,
      9.7229e-02, 1.0633e-01, 1.0703e-01, 1.0868e-01, 1.1124e-01, 1.1451e-01,
      1.1853e-01, 1.2342e-01, 1.2981e-01, 1.3868e-01, 1.5245e-01, 1.7661e-01,
      2.2796e-01, 5.6931e-01},
     {1.2735e-01, 1.2744e-01, 1.2775e-01, 1.2828e-01, 1.2909e-01, 1.3021e-01,
      1.3172e-01, 1.3377e-01, 1.3658e-01, 1.4055e-01, 1.4624e-01, 1.5549e-01,
      1.7490e-01, 2.1347e-01, 1.8445e-01, 1.7581e-01, 1.7263e-01, 1.7216e-01,
      1.7377e-01, 1.7731e-01, 1.8337e-01, 1.9306e-01, 2.0875e-01, 2.3610e-01,
      2.9100e-01, 6.0131e-01},
     {1.7654e-01, 1.7663e-01, 1.7711e-01, 1.7789e-01, 1.7901e-01, 1.8055e-01,
      1.8252e-01, 1.8511e-01, 1.8844e-01, 1.9288e-01, 1.9889e-01, 2.0816e-01,
      2.2673e-01, 2.7253e-01, 2.3949e-01, 2.3043e-01, 2.2695e-01, 2.2620e-01,
      2.2751e-01, 2.3091e-01, 2.3691e-01, 2.4671e-01, 2.6258e-01, 2.9036e-01,
      3.4575e-01, 6.4208e-01},
     {2.2845e-01, 2.2860e-01, 2.2907e-01, 2.2987e-01, 2.3103e-01, 2.3262e-01,
      2.3466e-01, 2.3727e-01, 2.4063e-01, 2.4510e-01, 2.5116e-01, 2.6032e-01,
      2.7853e-01, 3.2713e-01, 2.9240e-01, 2.8314e-01, 2.7962e-01, 2.7878e-01,
      2.8008e-01, 2.8341e-01, 2.8944e-01, 2.9918e-01, 3.1503e-01, 3.4275e-01,
      3.9786e-01, 6.8482e-01},
     {2.8077e-01, 2.8094e-01, 2.8142e-01, 2.8222e-01, 2.8341e-01, 2.8496e-01,
      2.8702e-01, 2.8962e-01, 2.9301e-01, 2.9745e-01, 3.0350e-01, 3.1263e-01,
      3.3069e-01, 3.8049e-01, 3.4500e-01, 3.3562e-01, 3.3206e-01, 3.3126e-01,
      3.3253e-01, 3.3588e-01, 3.4183e-01, 3.5155e-01, 3.6729e-01, 3.9494e-01,
      4.4955e-01, 7.2458e-01},
     {3.3321e-01, 3.3336e-01, 3.3385e-01, 3.3466e-01, 3.3584e-01, 3.3744e-01,
      3.3946e-01, 3.4213e-01, 3.4544e-01, 3.4987e-01, 3.5592e-01, 3.6502e-01,
      3.8298e-01, 4.3361e-01, 3.9758e-01, 3.8814e-01, 3.8456e-01, 3.8373e-01,
      3.8498e-01, 3.8831e-01, 3.9424e-01, 4.0392e-01, 4.1961e-01, 4.4703e-01,
      5.0121e-01, 7.6594e-01},
     {3.8572e-01, 3.8588e-01, 3.8634e-01, 3.8716e-01, 3.8834e-01, 3.8993e-01,
      3.9196e-01, 3.9456e-01, 3.9792e-01, 4.0236e-01, 4.0848e-01, 4.1745e-01,
      4.3532e-01, 4.8671e-01, 4.5019e-01, 4.4067e-01, 4.3708e-01, 4.3625e-01,
      4.3746e-01, 4.4078e-01, 4.4663e-01, 4.5626e-01, 4.7192e-01, 4.9913e-01,
      5.5270e-01, 8.0632e-01},
     {4.3847e-01, 4.3855e-01, 4.3892e-01, 4.3976e-01, 4.4089e-01, 4.4247e-01,
      4.4451e-01, 4.4711e-01, 4.5048e-01, 4.5489e-01, 4.6092e-01, 4.6997e-01,
      4.8773e-01, 5.3996e-01, 5.0285e-01, 4.9328e-01, 4.8966e-01, 4.8878e-01,
      4.9006e-01, 4.9327e-01, 4.9912e-01, 5.0870e-01, 5.2424e-01, 5.5126e-01,
      6.0420e-01, 8.4589e-01},
     {4.9089e-01, 4.9105e-01, 4.9153e-01, 4.9234e-01, 4.9351e-01, 4.9509e-01,
      4.9712e-01, 4.9970e-01, 5.0306e-01, 5.0749e-01, 5.1348e-01, 5.2250e-01,
      5.4016e-01, 5.9329e-01, 5.5559e-01, 5.4594e-01, 5.4229e-01, 5.4138e-01,
      5.4257e-01, 5.4584e-01, 5.5162e-01, 5.6114e-01, 5.7660e-01, 6.0338e-01,
      6.5560e-01, 8.8395e-01},
     {5.4358e-01, 5.4372e-01, 5.4422e-01, 5.4500e-01, 5.4618e-01, 5.4776e-01,
      5.4978e-01, 5.5238e-01, 5.5572e-01, 5.6014e-01, 5.6613e-01, 5.7508e-01,
      5.9263e-01, 6.4680e-01, 6.0840e-01, 5.9867e-01, 5.9497e-01, 5.9416e-01,
      5.9522e-01, 5.9846e-01, 6.0424e-01, 6.1369e-01, 6.2902e-01, 6.5556e-01,
      7.0707e-01, 9.1967e-01},
     {5.9631e-01, 5.9648e-01, 5.9694e-01, 5.9776e-01, 5.9892e-01, 6.0049e-01,
      6.0252e-01, 6.0511e-01, 6.0845e-01, 6.1287e-01, 6.1881e-01, 6.2775e-01,
      6.4519e-01, 7.0051e-01, 6.6128e-01, 6.5147e-01, 6.4773e-01, 6.4680e-01,
      6.4793e-01, 6.5114e-01, 6.5690e-01, 6.6624e-01, 6.8147e-01, 7.0775e-01,
      7.5854e-01, 9.5826e-01},
     {6.4916e-01, 6.4928e-01, 6.4976e-01, 6.5057e-01, 6.5173e-01, 6.5329e-01,
      6.5533e-01, 6.5792e-01, 6.6124e-01, 6.6563e-01, 6.7159e-01, 6.8051e-01,
      6.9779e-01, 7.5441e-01, 7.1424e-01, 7.0434e-01, 7.0057e-01, 6.9962e-01,
      7.0071e-01, 7.0392e-01, 7.0959e-01, 7.1891e-01, 7.3396e-01, 7.6003e-01,
      8.1008e-01, 9.9861e-01},
     {7.0207e-01, 7.0221e-01, 7.0264e-01, 7.0345e-01, 7.0461e-01, 7.0621e-01,
      7.0820e-01, 7.1079e-01, 7.1410e-01, 7.1849e-01, 7.2442e-01, 7.3329e-01,
      7.5047e-01, 8.0863e-01, 7.6728e-01, 7.5729e-01, 7.5348e-01, 7.5250e-01,
      7.5358e-01, 7.5674e-01, 7.6239e-01, 7.7162e-01, 7.8659e-01, 8.1236e-01,
      8.6164e-01, 1.0402},
     {7.5498e-01, 7.5573e-01, 7.5560e-01, 7.5642e-01, 7.5759e-01, 7.5913e-01,
      7.6116e-01, 7.6373e-01, 7.6704e-01, 7.7142e-01, 7.7734e-01, 7.8617e-01,
      8.0322e-01, 8.6277e-01, 8.2039e-01, 8.1030e-01, 8.0646e-01, 8.0546e-01,
      8.0650e-01, 8.0964e-01, 8.1524e-01, 8.2443e-01, 8.3926e-01, 8.6476e-01,
      9.1328e-01, 1.0826},
     {8.0804e-01, 8.0817e-01, 8.0865e-01, 8.0945e-01, 8.1060e-01, 8.1219e-01,
      8.1420e-01, 8.1674e-01, 8.2006e-01, 8.2444e-01, 8.3032e-01, 8.3912e-01,
      8.5604e-01, 9.1704e-01, 8.7360e-01, 8.6340e-01, 8.5952e-01, 8.5850e-01,
      8.5951e-01, 8.6262e-01, 8.6815e-01, 8.7727e-01, 8.9200e-01, 9.1725e-01,
      9.6506e-01, 1.1263},
     {8.6112e-01, 8.6129e-01, 8.6175e-01, 8.6255e-01, 8.6370e-01, 8.6526e-01,
      8.6729e-01, 8.6986e-01, 8.7316e-01, 8.7749e-01, 8.8338e-01, 8.9214e-01,
      9.0895e-01, 9.7166e-01, 9.2687e-01, 9.1657e-01, 9.1266e-01, 9.1157e-01,
      9.1260e-01, 9.1569e-01, 9.2125e-01, 9.3021e-01, 9.4481e-01, 9.6985e-01,
      1.0169, 1.1709},
     {9.1435e-01, 9.1448e-01, 9.1492e-01, 9.1575e-01, 9.1689e-01, 9.1846e-01,
      9.2046e-01, 9.2304e-01, 9.2633e-01, 9.3068e-01, 9.3652e-01, 9.4523e-01,
      9.6191e-01, 1.0265, 9.8022e-01, 9.6983e-01, 9.6588e-01, 9.6480e-01,
      9.6576e-01, 9.6877e-01, 9.7426e-01, 9.8323e-01, 9.9774e-01, 1.0225,
      1.0689, 1.2167},
     {9.6759e-01, 9.6773e-01, 9.6823e-01, 9.6900e-01, 9.7017e-01, 9.7170e-01,
      9.7372e-01, 9.7628e-01, 9.7957e-01, 9.8389e-01, 9.8975e-01, 9.9845e-01,
      1.0150, 1.0814, 1.0336, 1.0232, 1.0192, 1.0181, 1.0190, 1.0220, 1.0274,
      1.0363, 1.0507, 1.0752, 1.1210, 1.2635},
     {1.0219, 1.0211, 1.0215, 1.0223, 1.0235, 1.0250, 1.0270, 1.0296, 1.0329,
      1.0372, 1.0432, 1.0517, 1.0681, 1.1366, 1.0872, 1.0766, 1.0725, 1.0714,
      1.0723, 1.0753, 1.0807, 1.0895, 1.1038, 1.1281, 1.1732, 1.3110},
     {1.0743, 1.0748, 1.0749, 1.0757, 1.0770, 1.0784, 1.0804, 1.0830, 1.0863,
      1.0906, 1.0964, 1.1050, 1.1214, 1.1919, 1.1407, 1.1300, 1.1260, 1.1248,
      1.1257, 1.1287, 1.1340, 1.1428, 1.1570, 1.1810, 1.2256, 1.3591},
     {1.1278, 1.1279, 1.1284, 1.1292, 1.1304, 1.1319, 1.1339, 1.1364, 1.1397,
      1.1440, 1.1498, 1.1584, 1.1747, 1.2472, 1.1943, 1.1836, 1.1795, 1.1783,
      1.1792, 1.1821, 1.1874, 1.1961, 1.2102, 1.2340, 1.2781, 1.4079},
     {1.1912, 1.1815, 1.1820, 1.1827, 1.1839, 1.1854, 1.1875, 1.1900, 1.1933,
      1.1975, 1.2033, 1.2119, 1.2281, 1.3027, 1.2480, 1.2371, 1.2330, 1.2318,
      1.2327, 1.2356, 1.2409, 1.2495, 1.2635, 1.2872, 1.3308, 1.4571},
     {1.2349, 1.2351, 1.2356, 1.2363, 1.2375, 1.2390, 1.2410, 1.2436, 1.2468,
      1.2511, 1.2569, 1.2654, 1.2815, 1.3583, 1.3017, 1.2908, 1.2866, 1.2854,
      1.2863, 1.2892, 1.2944, 1.3030, 1.3170, 1.3404, 1.3836, 1.5071},
     {1.3170, 1.2887, 1.2892, 1.2900, 1.2911, 1.2927, 1.2947, 1.2972, 1.3005,
      1.3047, 1.3105, 1.3190, 1.3350, 1.4139, 1.3555, 1.3445, 1.3403, 1.3391,
      1.3399, 1.3428, 1.3480, 1.3566, 1.3704, 1.3937, 1.4366, 1.5575},
     {1.3423, 1.3424, 1.3429, 1.3437, 1.3449, 1.3464, 1.3484, 1.3509, 1.3541,
      1.3584, 1.3642, 1.3727, 1.3886, 1.4695, 1.4093, 1.3982, 1.3940, 1.3928,
      1.3936, 1.3966, 1.4017, 1.4103, 1.4240, 1.4472, 1.4896, 1.6085},
     {1.3960, 1.3962, 1.3966, 1.3974, 1.3985, 1.4001, 1.4021, 1.4046, 1.4079,
      1.4122, 1.4179, 1.4264, 1.4423, 1.5251, 1.4632, 1.4520, 1.4477, 1.4465,
      1.4474, 1.4502, 1.4554, 1.4639, 1.4776, 1.5007, 1.5429, 1.6600},
     {1.4498, 1.4499, 1.4504, 1.4512, 1.4524, 1.4539, 1.4559, 1.4584, 1.4616,
      1.4659, 1.4717, 1.4801, 1.4959, 1.5806, 1.5170, 1.5058, 1.5016, 1.5003,
      1.5011, 1.5040, 1.5091, 1.5177, 1.5313, 1.5542, 1.5962, 1.7119},
     {1.5036, 1.5037, 1.5042, 1.5050, 1.5061, 1.5077, 1.5097, 1.5122, 1.5154,
      1.5197, 1.5254, 1.5339, 1.5497, 1.6359, 1.5708, 1.5596, 1.5554, 1.5541,
      1.5549, 1.5577, 1.5629, 1.5714, 1.5850, 1.6079, 1.6497, 1.7643},
     {1.5574, 1.5575, 1.5580, 1.5588, 1.5600, 1.5615, 1.5635, 1.5660, 1.5693,
      1.5735, 1.5792, 1.5877, 1.6035, 1.6909, 1.6247, 1.6135, 1.6092, 1.6080,
      1.6088, 1.6117, 1.6167, 1.6252, 1.6388, 1.6616, 1.7033, 1.8171},
     {1.6173, 1.6114, 1.6119, 1.6127, 1.6138, 1.6153, 1.6173, 1.6198, 1.6231,
      1.6274, 1.6331, 1.6415, 1.6573, 1.7454, 1.6786, 1.6674, 1.6631, 1.6618,
      1.6626, 1.6654, 1.6706, 1.6791, 1.6926, 1.7154, 1.7570, 1.8704},
     {1.6651, 1.6652, 1.6657, 1.6665, 1.6676, 1.6692, 1.6712, 1.6737, 1.6769,
      1.6812, 1.6869, 1.6953, 1.7111, 1.7995, 1.7324, 1.7212, 1.7169, 1.7157,
      1.7164, 1.7193, 1.7244, 1.7328, 1.7465, 1.7693, 1.8108, 1.9241}};

  Float_t nonstagg[ktimebin][kZbin] =
    {{0.0912, 0.0913, 0.0914, 0.0917, 0.0920, 0.0925, 0.0931, 0.0938, 0.0946,
      0.0954, 0.0964, 0.0975, 0.0988, 0.1001, 0.1017, 0.1034, 0.1053, 0.1075,
      0.1101, 0.1135, 0.1179, 0.1239, 0.1332, 0.1504, 0.1935, 0.5719},
     {0.0657, 0.0658, 0.0660, 0.0662, 0.0666, 0.0672, 0.0678, 0.0685, 0.0694,
      0.0704, 0.0715, 0.0727, 0.0740, 0.0755, 0.0770, 0.0788, 0.0808, 0.0831,
      0.0858, 0.0893, 0.0938, 0.1001, 0.1096, 0.1274, 0.1715, 0.5475},
     {0.0396, 0.0397, 0.0400, 0.0404, 0.0410, 0.0418, 0.0427, 0.0438, 0.0450,
      0.0463, 0.0477, 0.0492, 0.0509, 0.0526, 0.0544, 0.0565, 0.0587, 0.0613,
      0.0645, 0.0685, 0.0736, 0.0807, 0.0919, 0.1128, 0.1608, 0.5233},
     {0.0112, 0.0115, 0.0124, 0.0137, 0.0154, 0.0173, 0.0195, 0.0218, 0.0241,
      0.0265, 0.0289, 0.0313, 0.0336, 0.0360, 0.0384, 0.0410, 0.0439, 0.0473,
      0.0515, 0.0571, 0.0651, 0.0767, 0.0949, 0.1243, 0.1826, 0.4858},
     {0.0113, 0.0116, 0.0124, 0.0137, 0.0154, 0.0173, 0.0195, 0.0217, 0.0241,
      0.0264, 0.0287, 0.0311, 0.0334, 0.0357, 0.0381, 0.0408, 0.0439, 0.0476,
      0.0525, 0.0591, 0.0684, 0.0821, 0.1019, 0.1352, 0.1914, 0.4601},
     {0.0394, 0.0394, 0.0397, 0.0401, 0.0407, 0.0415, 0.0423, 0.0434, 0.0445,
      0.0457, 0.0471, 0.0485, 0.0501, 0.0518, 0.0538, 0.0560, 0.0585, 0.0616,
      0.0652, 0.0698, 0.0759, 0.0846, 0.0983, 0.1226, 0.1742, 0.5225},
     {0.0650, 0.0650, 0.0652, 0.0655, 0.0658, 0.0663, 0.0669, 0.0676, 0.0685,
      0.0694, 0.0704, 0.0716, 0.0729, 0.0744, 0.0761, 0.0780, 0.0802, 0.0828,
      0.0860, 0.0899, 0.0950, 0.1020, 0.1130, 0.1334, 0.1809, 0.5444},
     {0.0929, 0.0929, 0.0922, 0.0921, 0.0923, 0.0928, 0.0934, 0.0941, 0.0950,
      0.0958, 0.0968, 0.0979, 0.0990, 0.1002, 0.1015, 0.1029, 0.1044, 0.1062,
      0.1080, 0.1101, 0.1126, 0.1155, 0.1192, 0.1241, 0.1312, 0.0000},
     {0.1688, 0.1657, 0.1518, 0.1450, 0.1409, 0.1382, 0.1366, 0.1357, 0.1352,
      0.1352, 0.1354, 0.1359, 0.1367, 0.1377, 0.1389, 0.1403, 0.1420, 0.1441,
      0.1465, 0.1495, 0.1532, 0.1580, 0.1643, 0.1736, 0.1902, 0.2919},
     {0.3413, 0.2191, 0.2052, 0.1981, 0.1936, 0.1905, 0.1883, 0.1869, 0.1860,
      0.1856, 0.1856, 0.1859, 0.1865, 0.1875, 0.1888, 0.1904, 0.1924, 0.1948,
      0.1977, 0.2011, 0.2054, 0.2107, 0.2175, 0.2272, 0.2440, 0.3413},
     {0.4016, 0.2716, 0.2576, 0.2505, 0.2459, 0.2428, 0.2406, 0.2392, 0.2383,
      0.2378, 0.2378, 0.2381, 0.2387, 0.2396, 0.2409, 0.2426, 0.2446, 0.2470,
      0.2500, 0.2535, 0.2578, 0.2631, 0.2700, 0.2796, 0.2964, 0.3877},
     {0.4388, 0.3240, 0.3100, 0.3029, 0.2983, 0.2952, 0.2931, 0.2916, 0.2907,
      0.2902, 0.2901, 0.2904, 0.2911, 0.2920, 0.2933, 0.2950, 0.2970, 0.2994,
      0.3024, 0.3059, 0.3102, 0.3155, 0.3224, 0.3319, 0.3486, 0.4320},
     {0.3580, 0.3764, 0.3625, 0.3554, 0.3509, 0.3477, 0.3455, 0.3440, 0.3431,
      0.3427, 0.3426, 0.3429, 0.3435, 0.3445, 0.3458, 0.3474, 0.3494, 0.3519,
      0.3548, 0.3583, 0.3626, 0.3679, 0.3747, 0.3843, 0.4008, 0.4768},
     {0.4023, 0.4289, 0.4150, 0.4079, 0.4033, 0.4002, 0.3980, 0.3965, 0.3956,
      0.3952, 0.3951, 0.3954, 0.3960, 0.3970, 0.3983, 0.3999, 0.4019, 0.4043,
      0.4072, 0.4107, 0.4150, 0.4203, 0.4271, 0.4367, 0.4530, 0.5228},
     {0.4334, 0.4815, 0.4676, 0.4604, 0.4559, 0.4527, 0.4506, 0.4491, 0.4482,
      0.4477, 0.4476, 0.4479, 0.4486, 0.4495, 0.4508, 0.4524, 0.4544, 0.4568,
      0.4597, 0.4632, 0.4675, 0.4728, 0.4796, 0.4890, 0.5052, 0.5699},
     {0.4695, 0.5341, 0.5202, 0.5130, 0.5085, 0.5053, 0.5032, 0.5017, 0.5008,
      0.5003, 0.5002, 0.5005, 0.5012, 0.5021, 0.5034, 0.5050, 0.5070, 0.5095,
      0.5123, 0.5158, 0.5200, 0.5253, 0.5321, 0.5415, 0.5576, 0.6179},
     {0.5905, 0.5868, 0.5729, 0.5657, 0.5612, 0.5580, 0.5559, 0.5544, 0.5535,
      0.5530, 0.5529, 0.5532, 0.5538, 0.5548, 0.5560, 0.5577, 0.5597, 0.5621,
      0.5649, 0.5684, 0.5726, 0.5779, 0.5846, 0.5940, 0.6099, 0.6666},
     {0.5940, 0.6307, 0.6254, 0.6185, 0.6141, 0.6108, 0.6086, 0.6071, 0.6062,
      0.6057, 0.6057, 0.6059, 0.6066, 0.6075, 0.6088, 0.6104, 0.6124, 0.6148,
      0.6176, 0.6211, 0.6253, 0.6306, 0.6373, 0.6465, 0.6623, 0.7156},
     {0.6986, 0.6924, 0.6784, 0.6713, 0.6667, 0.6636, 0.6614, 0.6599, 0.6590,
      0.6585, 0.6585, 0.6587, 0.6594, 0.6603, 0.6616, 0.6632, 0.6651, 0.6675,
      0.6704, 0.6739, 0.6781, 0.6833, 0.6900, 0.6992, 0.7147, 0.7651},
     {0.7449, 0.7453, 0.7312, 0.7242, 0.7196, 0.7165, 0.7143, 0.7128, 0.7119,
      0.7117, 0.7113, 0.7116, 0.7122, 0.7132, 0.7144, 0.7160, 0.7180, 0.7204,
      0.7233, 0.7267, 0.7309, 0.7361, 0.7427, 0.7519, 0.7672, 0.8151},
     {0.8011, 0.7982, 0.7843, 0.7772, 0.7726, 0.7694, 0.7673, 0.7658, 0.7649,
      0.7644, 0.7643, 0.7646, 0.7652, 0.7661, 0.7674, 0.7690, 0.7709, 0.7733,
      0.7763, 0.7796, 0.7838, 0.7890, 0.7956, 0.8046, 0.8198, 0.8656},
     {0.8261, 0.8513, 0.8374, 0.8302, 0.8256, 0.8225, 0.8203, 0.8188, 0.8179,
      0.8174, 0.8173, 0.8176, 0.8182, 0.8191, 0.8204, 0.8220, 0.8240, 0.8263,
      0.8292, 0.8326, 0.8367, 0.8419, 0.8485, 0.8575, 0.8725, 0.9163},
     {0.9145, 0.9043, 0.8905, 0.8834, 0.8788, 0.8756, 0.8734, 0.8720, 0.8710,
      0.8706, 0.8705, 0.8707, 0.8713, 0.8722, 0.8735, 0.8751, 0.8770, 0.8794,
      0.8822, 0.8856, 0.8898, 0.8949, 0.9015, 0.9104, 0.9252, 0.9674},
     {0.8915, 0.9576, 0.9437, 0.9365, 0.9320, 0.9288, 0.9267, 0.9251, 0.9242,
      0.9237, 0.9236, 0.9239, 0.9245, 0.9254, 0.9267, 0.9282, 0.9302, 0.9326,
      0.9354, 0.9388, 0.9429, 0.9480, 0.9545, 0.9635, 0.9781, 1.0188},
     {1.0290, 1.0109, 0.9970, 0.9898, 0.9852, 0.9821, 0.9799, 0.9784, 0.9775,
      0.9770, 0.9769, 0.9772, 0.9778, 0.9787, 0.9800, 0.9815, 0.9835, 0.9858,
      0.9886, 0.9920, 0.9961, 1.0012, 1.0077, 1.0165, 1.0310, 1.0705},
     {1.0470, 1.1170, 1.0503, 1.0431, 1.0386, 1.0354, 1.0333, 1.0318, 1.0308,
      1.0303, 1.0303, 1.0305, 1.0311, 1.0320, 1.0333, 1.0348, 1.0368, 1.0391,
      1.0419, 1.0453, 1.0494, 1.0545, 1.0609, 1.0697, 1.0840, 1.1223},
     {1.2284, 1.1176, 1.1037, 1.0966, 1.0922, 1.0889, 1.0867, 1.0852, 1.0842,
      1.0838, 1.0837, 1.0839, 1.0845, 1.0854, 1.0867, 1.0882, 1.0902, 1.0925,
      1.0953, 1.0987, 1.1028, 1.1078, 1.1143, 1.1230, 1.1371, 1.1744},
     {1.1548, 1.1712, 1.1572, 1.1501, 1.1455, 1.1424, 1.1402, 1.1387, 1.1377,
      1.1373, 1.1371, 1.1374, 1.1380, 1.1389, 1.1401, 1.1417, 1.1436, 1.1460,
      1.1488, 1.1521, 1.1562, 1.1613, 1.1676, 1.1763, 1.1903, 1.2267},
     {1.2050, 1.2248, 1.2108, 1.2037, 1.1991, 1.1959, 1.1937, 1.1922, 1.1913,
      1.1908, 1.1907, 1.1909, 1.1915, 1.1925, 1.1937, 1.1952, 1.1972, 1.1995,
      1.2023, 1.2056, 1.2097, 1.2147, 1.2211, 1.2297, 1.2436, 1.2791},
     {1.2799, 1.2784, 1.2644, 1.2573, 1.2527, 1.2495, 1.2473, 1.2458, 1.2449,
      1.2444, 1.2443, 1.2446, 1.2451, 1.2460, 1.2473, 1.2488, 1.2508, 1.2531,
      1.2559, 1.2592, 1.2633, 1.2683, 1.2746, 1.2833, 1.2970, 1.3317},
     {1.2857, 1.3320, 1.3181, 1.3109, 1.3063, 1.3032, 1.3010, 1.2995, 1.2986,
      1.2981, 1.2980, 1.2982, 1.2988, 1.2997, 1.3009, 1.3025, 1.3044, 1.3067,
      1.3095, 1.3129, 1.3169, 1.3219, 1.3283, 1.3368, 1.3505, 1.3845},
     {1.3924, 1.3856, 1.3718, 1.3646, 1.3601, 1.3569, 1.3547, 1.3532, 1.3523,
      1.3518, 1.3517, 1.3519, 1.3525, 1.3534, 1.3546, 1.3562, 1.3581, 1.3604,
      1.3632, 1.3665, 1.3706, 1.3756, 1.3819, 1.3904, 1.4041, 1.4375},
     {1.5732, 1.4395, 1.4255, 1.4184, 1.4138, 1.4107, 1.4085, 1.4070, 1.4060,
      1.4055, 1.4054, 1.4057, 1.4063, 1.4071, 1.4086, 1.4099, 1.4119, 1.4142,
      1.4170, 1.4203, 1.4243, 1.4293, 1.4356, 1.4441, 1.4576, 1.4906},
     {1.4352, 1.4933, 1.4793, 1.4722, 1.4676, 1.4644, 1.4623, 1.4607, 1.4598,
      1.4593, 1.4592, 1.4595, 1.4600, 1.4609, 1.4622, 1.4637, 1.4656, 1.4680,
      1.4707, 1.4741, 1.4781, 1.4831, 1.4894, 1.4979, 1.5113, 1.5439},
     {1.4933, 1.5471, 1.5332, 1.5260, 1.5214, 1.5183, 1.5161, 1.5146, 1.5136,
      1.5131, 1.5131, 1.5133, 1.5139, 1.5148, 1.5160, 1.5175, 1.5195, 1.5218,
      1.5247, 1.5279, 1.5319, 1.5369, 1.5432, 1.5516, 1.5650, 1.5973},
     {1.5709, 1.6010, 1.5870, 1.5798, 1.5753, 1.5721, 1.5699, 1.5684, 1.5675,
      1.5670, 1.5669, 1.5671, 1.5677, 1.5686, 1.5698, 1.5714, 1.5733, 1.5756,
      1.5784, 1.5817, 1.5857, 1.5907, 1.5970, 1.6054, 1.6188, 1.6509},
     {1.6675, 1.6548, 1.6408, 1.6337, 1.6291, 1.6259, 1.6238, 1.6223, 1.6213,
      1.6208, 1.6207, 1.6210, 1.6216, 1.6225, 1.6237, 1.6252, 1.6272, 1.6295,
      1.6322, 1.6356, 1.6396, 1.6445, 1.6508, 1.6593, 1.6727, 1.7046},
     {1.7042, 1.7087, 1.6947, 1.6876, 1.6830, 1.6798, 1.6776, 1.6761, 1.6752,
      1.6747, 1.6746, 1.6749, 1.6754, 1.6763, 1.6775, 1.6791, 1.6810, 1.6833,
      1.6861, 1.6894, 1.6935, 1.6984, 1.7047, 1.7132, 1.7265, 1.7584}};

  if (fTimeStruct)  delete [] fTimeStruct;
  fTimeStruct  = new Float_t[ktimebin*kZbin];

  for (Int_t ctrt = 0; ctrt<ktimebin; ctrt++)
    for (Int_t ctrz = 0; ctrz<kZbin; ctrz++) {
      if (fStaggeringOn) fTimeStruct[ctrt+ctrz*ktimebin] = stagg[ctrt][ctrz];
      else               fTimeStruct[ctrt+ctrz*ktimebin] = nonstagg[ctrt][ctrz];
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

  Float_t prf[kNplan][kPRFbin] = { {0.0267681, 0.031056,  0.036111,  0.0419945, 0.0489578,
				    0.0571508, 0.0666984, 0.078064,  0.0912508, 0.106579,
				    0.124429,  0.144803,  0.168187,  0.194645,  0.224143,
				    0.25717,   0.293146,  0.331943,  0.373128,  0.415798,
				    0.459149,  0.502028,  0.54349,   0.582109,  0.617335,
				    0.648306,  0.674189,  0.695082,  0.710222,  0.719711,
				    0.723751,  0.719711,  0.710222,  0.695082,  0.674189,
				    0.648306,  0.617335,  0.582109,  0.54349,   0.502028,
				    0.459149,  0.415798,  0.373128,  0.331943,  0.293146,
				    0.25717,   0.224143,  0.194645,  0.168187,  0.144803,
				    0.124429,  0.106579,  0.0912508, 0.078064,  0.0666984,
				    0.0571508, 0.0489578, 0.0419945, 0.036111,  0.031056,
				    0.0267681},
				   {0.0240831, 0.0281427, 0.0328987, 0.0384801, 0.0451249,
				    0.0530119, 0.0623139, 0.0733268, 0.0864095, 0.101704,
				    0.119473,  0.140308,  0.164191,  0.191333,  0.222141,
				    0.256551,  0.29433,   0.335294,  0.37892,   0.424206,
				    0.470202,  0.515621,  0.559255,  0.600068,  0.636821,
				    0.668836,  0.695767,  0.717076,  0.73233,   0.741818,
				    0.745491,  0.741818,  0.73233,   0.717076,  0.695767,
				    0.668836,  0.636821,  0.600068,  0.559255,  0.515621,
				    0.470202,  0.424206,  0.37892,   0.335294,  0.29433,
				    0.256551,  0.222141,  0.191333,  0.164191,  0.140308,
				    0.119473,  0.101704,  0.0864095, 0.0733268, 0.0623139,
				    0.0530119, 0.0451249, 0.0384801, 0.0328987, 0.0281427,
				    0.0240831},
				   {0.0206855, 0.0243011, 0.0285603, 0.0335902, 0.0396015,
				    0.0467843, 0.0553424, 0.0655316, 0.0776688, 0.0921782,
				    0.109269,  0.129259,  0.152459,  0.179369,  0.210091,
				    0.244681,  0.283099,  0.325176,  0.370396,  0.417732,
				    0.466086,  0.514192,  0.560518,  0.603805,  0.643028,
				    0.67738,   0.706095,  0.728728,  0.745307,  0.755731,
				    0.759954,  0.755731,  0.745307,  0.728728,  0.706095,
				    0.67738,   0.643028,  0.603805,  0.560518,  0.514192,
				    0.466086,  0.417732,  0.370396,  0.325176,  0.283099,
				    0.244681,  0.210091,  0.179369,  0.152459,  0.129259,
				    0.109269,  0.0921782, 0.0776688, 0.0655316, 0.0553424,
				    0.0467843, 0.0396015, 0.0335902, 0.0285603, 0.0243011,
				    0.0206855},
				   {0.0186168, 0.0219999, 0.0260102, 0.0307769, 0.0364947,
				    0.0433655, 0.0516213, 0.0615466, 0.0734611, 0.0877121,
				    0.104666,  0.124855,  0.14853,   0.176033,  0.207639,
				    0.243511,  0.283633,  0.327786,  0.37537,   0.425281,
				    0.476227,  0.526727,  0.575268,  0.620462,  0.66101,
				    0.69611,   0.725313,  0.748317,  0.764969,  0.775206,
				    0.779006,  0.775206,  0.764969,  0.748317,  0.725313,
				    0.69611,   0.66101,   0.620462,  0.575268,  0.526727,
				    0.476227,  0.425281,  0.37537,   0.327786,  0.283633,
				    0.243511,  0.207639,  0.176033,  0.14853,   0.124855,
				    0.104666,  0.0877121, 0.0734611, 0.0615466, 0.0516213,
				    0.0433655, 0.0364947, 0.0307769, 0.0260102, 0.0219999,
				    0.0186168, },
 				   {0.0159737, 0.0189921, 0.0225916, 0.0268927, 0.0320634,
				    0.0382995, 0.0458393, 0.0549765, 0.0660512, 0.0794439,
				    0.095565,  0.114844,  0.137713,  0.164586,  0.195824,
				    0.231681,  0.272223,  0.31727,   0.366248,  0.418078,
				    0.471358,  0.524425,  0.575561,  0.623193,  0.666055,
				    0.703243,  0.734192,  0.7586,    0.776331,  0.787347,
				    0.791646,  0.787347,  0.776331,  0.7586,    0.734192,
				    0.703243,  0.666055,  0.623193,  0.575561,  0.524425,
				    0.471358,  0.418078,  0.366248,  0.31727,   0.272223,
				    0.231681,  0.195824,  0.164586,  0.137713,  0.114844,
				    0.095565,  0.0794439, 0.0660512, 0.0549765, 0.0458393,
				    0.0382995, 0.0320634, 0.0268927, 0.0225916, 0.0189921,
				    0.0159737},
				   {0.0143532, 0.0171745, 0.0205622, 0.024635,  0.0295461,
				    0.0354915, 0.0427198, 0.0515368, 0.0623058, 0.0754416,
				    0.0913994, 0.110662,  0.133721,  0.161057,  0.193093,
				    0.230138,  0.27229,   0.319327,  0.370597,  0.424961,
				    0.480832,  0.536333,  0.589555,  0.638805,  0.682781,
				    0.720621,  0.751845,  0.776243,  0.793764,  0.804426,
				    0.808259,  0.804426,  0.793764,  0.776243,  0.751845,
				    0.720621,  0.682781,  0.638805,  0.589555,  0.536333,
				    0.480832,  0.424961,  0.370597,  0.319327,  0.27229,
				    0.230138,  0.193093,  0.161057,  0.133721,  0.110662,
				    0.0913994, 0.0754416, 0.0623058, 0.0515368, 0.0427198,
				    0.0354915, 0.0295461, 0.024635,  0.0205622, 0.0171745,
				    0.0143532} };

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
