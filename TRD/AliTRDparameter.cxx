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
void AliTRDparameter::Copy(TObject &p) const
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
  ((AliTRDparameter &) p).fTimeStructOn       = fTimeStructOn;
  if (((AliTRDparameter &) p).fTimeStruct) 
    delete [] ((AliTRDparameter &) p).fTimeStruct;
  ((AliTRDparameter &) p).fTimeStruct = new Float_t[38*11];
  for (Int_t i = 0; i < 38*11; i++) {
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
  SetColPadSize(0,0.664);
  SetColPadSize(1,0.695);
  SetColPadSize(2,0.726);
  SetColPadSize(3,0.756);
  SetColPadSize(4,0.788);
  SetColPadSize(5,0.818);

  // The pad row (z-direction)
  SetNRowPad();

  // The number of time bins. Default is 100 ns timbin size
  SetNTimeBin(18);

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
  fDriftVelocity  = 1.62;

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
  fTCnexp         = 1;

  // The pad coupling factor
  //fPadCoupling    = 0.3;
  // Use 0.46 instead which reproduces better the test beam
  // data, even tough it is not understood why.
  fPadCoupling    = 0.46;

  // The time coupling factor (same number as for the TPC)
  fTimeCoupling   = 0.4;

  // Drift time non-isochronity on
  fTimeStructOn   = 1;

  // Distance of first Anode wire from first pad edge
  fAnodeWireOffset = 0.25;

  // The tilting angle for the readout pads
  SetTiltingAngle(2.0);

  // The magnetic field strength in Tesla
  Double_t x[3] = { 0.0, 0.0, 0.0 };
  Double_t b[3];
  gAlice->Field(x,b);  // b[] is in kilo Gauss
  fField = b[2] * 0.1; // Tesla

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

        fRow0[iplan][icham][isect]          = row0;
	// For new chamber ordering
        //fRow0[iplan][kNcham-icham-1][isect] = row0;

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
Int_t AliTRDparameter::TimeStruct(Float_t dist, Float_t z, Float_t *xyz) const
{
  //
  // Applies the time structure of the drift cells (by C.Lippmann).
  // The drift time of electrons to the anode wires depends on the
  // distance to the wire (z) and on the position in the drift region.
  // Here we omnly reposition each electron, since the drift velocity
  // in the digitisation is still constant. If the electrons at a given
  // position has a long drift time, it is just moved further away from
  // the anode wire plane. 
  // 
  // input :
  // dist = radial distance from anode wire plane [cm]
  // z    = distance from anode wire (parallel to cathode planes) [cm]
  //


  // indices:
  Int_t r1 = (Int_t)(10*dist);
  Int_t r2 = r1+1;
  Int_t z1 = (Int_t)(100*z/2.5);
  Int_t z2 = z1+1;

  if (r1<0 || r1>37 || z1<0 || z1>10) {
    printf("<AliTRDparameter::TimeStruct> Warning. Indices out of range: ");
    printf("dist=%.2f, z=%.2f, r1=%d, z1=%d\n",dist,z,r1,z1);
    return kFALSE;
  }

  // 2D Interpolation:
  Float_t y11 = fTimeStruct[r1+38*z1];
  Float_t y22 = (r2 <= 37 && z2 <= 10) ? fTimeStruct[r2+38*z2] : fTimeStruct[37+38*10];
  Float_t y12 = (z2 <= 10)             ? fTimeStruct[r1+38*z2] : fTimeStruct[r1+38*10];
  Float_t y21 = (r2 <= 37)             ? fTimeStruct[r2+38*z1] : fTimeStruct[37+38*z1];

  Float_t y1  = (y21-y11)*10*dist + y11 - (y21-y11)*r1;
  Float_t y2  = (y22-y12)*10*dist + y12 - (y22-y12)*r1;

  Float_t AmTh = AliTRDgeometry::AmThick()/2.0;

  // dist now is the drift distance to anode wires (negative if electrons are
  // between anode wire plane and cathode pad plane)
  dist -= AmTh;

  // Get the drift time from the interpolation:
  Float_t tdrift =
    ( TMath::Abs(dist)>0.005 || z>0.005 ) ? tdrift = (y2-y1)*100*z/2.5+y1-(y2-y1)*z1 : 0.0;

  // We move electrons further away from the anode wire plane, if the drift time is
  // larger than the one expected for that distance for a constant drift velocity. We
  // increase (decrease) the drift distance by a factor
  //
  //                           c = t(x0)*vd/x0 ,
  //
  // where x0 is the distance from the anode wire plane (dist), vd is the constant
  // drift velocity (fDriftVelocity) and t(x0) is the real drift time (tdrift, as
  // calculated with GARFIELD). 
  //  
  // The factor is negative for electrons between anode wire plane and cathode pads.
  //
  // The new position of the electron is then given by: t(x0)*vd:

  Float_t zdrift = tdrift * fDriftVelocity;
  if (dist < 0.) zdrift *= -1.; 

  xyz[0] = xyz[0] + dist - zdrift;
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
  const Int_t kZbin = 11;
  
  Float_t time[ktimebin][kZbin] =
    {{0.10176, 0.10197, 0.10356, 0.10479, 0.10497, 0.10500, 0.10500,
      0.10698, 0.12213, 0.16118, 0.22154},
     {0.07500, 0.07500, 0.07500, 0.07500, 0.07500, 0.07500, 0.07509,
      0.08010, 0.10176, 0.13244, 0.22848},
     {0.04500, 0.04500, 0.04500, 0.04500, 0.04500, 0.04503, 0.05130,
      0.07398, 0.07557, 0.11064, 0.26833},
     {0.01500, 0.01500, 0.01500, 0.01500, 0.03219, 0.04500, 0.04500,
      0.04500, 0.07482, 0.11628, 0.40846},
     {0.01500, 0.01500, 0.01500, 0.01500, 0.02886, 0.04500, 0.04500,
      0.04512, 0.07503, 0.13230, 0.43602},
     {0.04500, 0.04500, 0.04500, 0.04500, 0.04500, 0.04500, 0.04872,
      0.07437, 0.07635, 0.11925, 0.30372},
     {0.07500, 0.07497, 0.07500, 0.07500, 0.07500, 0.07500, 0.07500,
      0.08232, 0.10458, 0.14100, 0.24799},
     {0.10098, 0.10125, 0.10218, 0.10299, 0.10467, 0.15000, 0.10959,
      0.12207, 0.14241, 0.19485, 0.25688},
     {0.13500, 0.13503, 0.13509, 0.13782, 0.16182, 0.25076, 0.19593,
      0.19254, 0.21060, 0.25561, 0.29469},
     {0.19305, 0.19374, 0.19512, 0.19974, 0.22020, 0.29420, 0.25992,
      0.25461, 0.27515, 0.31938, 0.35238},
     {0.25143, 0.25320, 0.25476, 0.26166, 0.27993, 0.34470, 0.32160,
      0.31623, 0.33687, 0.37420, 0.40797},
     {0.31044, 0.31197, 0.31491, 0.32280, 0.34125, 0.39702, 0.38401,
      0.37794, 0.39545, 0.43494, 0.46956},
     {0.37128, 0.37146, 0.37668, 0.38424, 0.40326, 0.45271, 0.44643,
      0.43783, 0.45882, 0.49726, 0.52348},
     {0.43065, 0.43242, 0.43722, 0.44379, 0.46551, 0.51167, 0.50685,
      0.49924, 0.52016, 0.55609, 0.57829},
     {0.49161, 0.49329, 0.49725, 0.50629, 0.52803, 0.56777, 0.56755,
      0.56057, 0.58062, 0.61293, 0.63839},
     {0.55284, 0.55311, 0.55743, 0.56703, 0.58824, 0.62747, 0.62833,
      0.62083, 0.63899, 0.67412, 0.70006},
     {0.61323, 0.61494, 0.61782, 0.62802, 0.65025, 0.68675, 0.68847,
      0.68248, 0.70177, 0.73588, 0.76314},
     {0.67293, 0.67482, 0.68034, 0.68814, 0.71169, 0.74695, 0.75024,
      0.74462, 0.76083, 0.79194, 0.81463},
     {0.73470, 0.73665, 0.74013, 0.75054, 0.77107, 0.80588, 0.81071,
      0.80508, 0.82337, 0.85380, 0.87045},
     {0.79509, 0.79677, 0.80115, 0.81150, 0.83400, 0.86360, 0.87159,
      0.86608, 0.88323, 0.91756, 0.93721},
     {0.85695, 0.85800, 0.86307, 0.87276, 0.89599, 0.92639, 0.93266,
      0.92962, 0.94291, 0.97445, 0.99502},
     {0.91680, 0.91875, 0.92400, 0.93219, 0.95555, 0.98729, 0.99281,
      0.99145, 1.00683, 1.03343, 1.05299},
     {0.97833, 0.98058, 0.98421, 0.99357, 1.01602, 1.04438, 1.05299,
      1.05168, 1.06528, 1.09644, 1.11557},
     {1.03938, 1.04061, 1.04568, 1.05617, 1.07923, 1.10648, 1.11466,
      1.11385, 1.12606, 1.15817, 1.17352},
     {1.10166, 1.10268, 1.10766, 1.11785, 1.14014, 1.16717, 1.17496,
      1.17644, 1.19115, 1.21622, 1.23605},
     {1.16334, 1.16364, 1.16985, 1.17996, 1.19935, 1.22911, 1.23561,
      1.23526, 1.25144, 1.27936, 1.29018},
     {1.22403, 1.22598, 1.23000, 1.24092, 1.26173, 1.28817, 1.29991,
      1.29698, 1.31238, 1.33793, 1.35595},
     {1.28556, 1.28898, 1.29300, 1.30356, 1.32617, 1.34912, 1.35925,
      1.35835, 1.37312, 1.40111, 1.41497},
     {1.34787, 1.34877, 1.35468, 1.36449, 1.38654, 1.41077, 1.41975,
      1.42136, 1.43843, 1.45897, 1.48025},
     {1.41045, 1.41171, 1.41555, 1.42656, 1.44698, 1.47026, 1.48324,
      1.48479, 1.50180, 1.52354, 1.53392},
     {1.47090, 1.47354, 1.47647, 1.48905, 1.51211, 1.53038, 1.54349,
      1.54803, 1.56136, 1.58662, 1.60480},
     {1.53357, 1.53429, 1.54017, 1.55116, 1.57247, 1.59498, 1.60576,
      1.60746, 1.62389, 1.64239, 1.66245},
     {1.59522, 1.59636, 1.60247, 1.61385, 1.63641, 1.65659, 1.66619,
      1.67044, 1.68213, 1.70379, 1.72372},
     {1.65867, 1.65711, 1.66320, 1.67434, 1.69603, 1.71788, 1.72927,
      1.73214, 1.74713, 1.76485, 1.77775},
     {1.71825, 1.72074, 1.72497, 1.73712, 1.75690, 1.77755, 1.78711,
      1.79849, 1.80587, 1.82576, 1.84050},
     {1.78194, 1.78320, 1.78875, 1.80047, 1.81716, 1.84044, 1.85328,
      1.85511, 1.86911, 1.88919, 1.89926},
     {1.84335, 1.84431, 1.84938, 1.86063, 1.88022, 1.90234, 1.91177,
      1.92173, 1.93202, 1.95097, 1.97105},
     {1.90491, 1.90665, 1.91303, 1.92575, 1.94432, 1.96400, 1.97574,
      1.98378, 1.99368, 2.01964, 2.02605}};

  if (fTimeStruct)  delete [] fTimeStruct;
  fTimeStruct  = new Float_t[ktimebin*kZbin];

  for (Int_t ctrt = 0; ctrt<ktimebin; ctrt++)
    for (Int_t ctrz = 0; ctrz<kZbin; ctrz++) {
      fTimeStruct[ctrt+ctrz*ktimebin] = time[ctrt][ctrz];
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
