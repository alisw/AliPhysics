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
#include "AliTRDpadPlane.h"

ClassImp(AliTRDparameter)

//_____________________________________________________________________________
AliTRDparameter::AliTRDparameter():TNamed()
{
  //
  // AliTRDparameter default constructor
  //

  fGeo                = 0;
  fPadPlaneArray      = 0;
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
  fDiffusionT         = 0.0;
  fDiffusionL         = 0.0;
  fElAttachProp       = 0.0;
  fOmegaTau           = 0.0;
  fDriftVelocity      = 0.0;
  fSamplingFrequency  = 0.0;
  fPadCoupling        = 0.0;
  fTimeCoupling       = 0.0;
  fField              = 0.0;
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

  fLUT                = 0;
  fClusMaxThresh      = 0;
  fClusSigThresh      = 0;

  fTimeStruct1        = 0;
  fTimeStruct2        = 0;
  fAnodeWireOffset    = 0.0;

  fDiffusionOn        = kFALSE;
  fElAttachOn         = kFALSE;
  fExBOn              = kFALSE;
  fPRFOn              = kFALSE;
  fTRFOn              = kFALSE;
  fCTOn               = kFALSE;
  fTCOn               = kFALSE;
  fLUTOn              = kFALSE;  
  fTimeStructOn       = kFALSE;

}

//_____________________________________________________________________________
AliTRDparameter::AliTRDparameter(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDparameter constructor
  //

  fGeo                = new AliTRDgeometryFull();
  fPadPlaneArray      = new TObjArray(AliTRDgeometry::Nplan() 
                                    * AliTRDgeometry::Ncham());
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
  fDiffusionT         = 0.0;
  fDiffusionL         = 0.0;
  fElAttachProp       = 0.0;
  fExBOn              = 0;
  fOmegaTau           = 0.0;
  fDriftVelocity      = 0.0;
  fSamplingFrequency  = 0.0;
  fPadCoupling        = 0.0;
  fTimeCoupling       = 0.0;
  fField              = 0.0;
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

  fLUT                = 0;
  fClusMaxThresh      = 0;
  fClusSigThresh      = 0;

  fTimeStruct1        = 0;
  fTimeStruct2        = 0;
  fAnodeWireOffset    = 0.0;

  fDiffusionOn        = kFALSE;
  fElAttachOn         = kFALSE;
  fExBOn              = kFALSE;
  fPRFOn              = kFALSE;
  fTRFOn              = kFALSE;
  fCTOn               = kFALSE;
  fTCOn               = kFALSE;
  fLUTOn              = kFALSE;  
  fTimeStructOn       = kFALSE;

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

  if (fPadPlaneArray) {
    fPadPlaneArray->Delete();
    delete fPadPlaneArray;
    fPadPlaneArray = 0;
  }

  if (fTimeStruct1) {
    delete [] fTimeStruct1;
    fTimeStruct1 = 0;
  }

  if (fTimeStruct2) {
    delete [] fTimeStruct2;
    fTimeStruct2 = 0;
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
  ((AliTRDparameter &) p).fSamplingFrequency  = fSamplingFrequency;
  ((AliTRDparameter &) p).fPadCoupling        = fPadCoupling;
  ((AliTRDparameter &) p).fTimeCoupling       = fTimeCoupling;
  ((AliTRDparameter &) p).fField              = fField;
  ((AliTRDparameter &) p).fPRFOn              = fPRFOn;
  ((AliTRDparameter &) p).fTRFOn              = fTRFOn;
  ((AliTRDparameter &) p).fCTOn               = fCTOn;
  ((AliTRDparameter &) p).fTCOn               = fTCOn;
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
  if (((AliTRDparameter &) p).fTimeStruct1) 
    delete [] ((AliTRDparameter &) p).fTimeStruct1;
  ((AliTRDparameter &) p).fTimeStruct1 = new Float_t[38*11];
  for (Int_t i = 0; i < 38*11; i++) {
    ((AliTRDparameter &) p).fTimeStruct1[i] = fTimeStruct1[i];
  }
  if (((AliTRDparameter &) p).fTimeStruct2) 
    delete [] ((AliTRDparameter &) p).fTimeStruct2;
  ((AliTRDparameter &) p).fTimeStruct2 = new Float_t[38*11];
  for (Int_t i = 0; i < 38*11; i++) {
    ((AliTRDparameter &) p).fTimeStruct2[i] = fTimeStruct2[i];
  }

}

//_____________________________________________________________________________
void AliTRDparameter::Init()
{
  //
  // Initializes the parameter
  //

  Int_t iplan;
  Int_t icham;

  //
  // ----------------------------------------------------------------------------
  // The pad planes
  // ----------------------------------------------------------------------------
  //
  for (iplan = 0; iplan < AliTRDgeometry::Nplan(); iplan++) {
    for (icham = 0; icham < AliTRDgeometry::Ncham(); icham++) {
      Int_t ipp = iplan + icham * AliTRDgeometry::Nplan();
      fPadPlaneArray->AddAt(new AliTRDpadPlane(iplan,icham),ipp);
    }
  }

  // !CHANGED!
  // Time position of anode wire plane
  for (iplan = 0; iplan < AliTRDgeometry::Nplan(); iplan++) {
    fTime0[iplan] = AliTRDgeometry::Rmin() 
                  + AliTRDgeometry::CraHght() 
                  + AliTRDgeometry::CdrHght()
	          + AliTRDgeometry::CamHght()/2. 
                  + iplan * (AliTRDgeometry::Cheight() + AliTRDgeometry::Cspace());
  }

  //
  // ----------------------------------------------------------------------------
  // The digitization parameters
  // ----------------------------------------------------------------------------
  //

  // Use drift time maps
  fTimeStructOn = kTRUE;

  // The Sampling Frequency. Default is 10MHz
  fSamplingFrequency = 10.0;

  // The drift velocity (in drift region) and the time structure of the
  // drift cells. Default is 1.5 cm/mus
  SetDriftVelocity(1.5);

  // Additional time bins before and after the drift region.
  // Default is to only sample the drift region
  SetExpandTimeBin(0,0);

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

  // E x B effects
  fExBOn          = kTRUE;

  // Propability for electron attachment
  fElAttachOn     = kFALSE;
  fElAttachProp   = 0.0;

  // The pad response function
  fPRFOn          = kTRUE;

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
  fLUTOn         = kTRUE;

  ReInit();

}

//_____________________________________________________________________________
void AliTRDparameter::ReInit()
{
  //
  // Reinitializes the parameter class after a change
  //

  Float_t maxdrift =  AliTRDgeometry::DrThick() + AliTRDgeometry::AmThick()/2.0;

  // Number of time bins:
  fTimeMax = (Int_t) ( maxdrift / fDriftVelocity * fSamplingFrequency );

  // The range and the binwidth for the sampled TRF 
  fTRFbin = 100;
  // Start 0.2 mus before the signal
  fTRFlo  = -0.2;
  // End the maximum drift time after the signal 
  fTRFhi  = 2.2;
  // 
  fTRFwid = (fTRFhi - fTRFlo) / ((Float_t) fTRFbin);

  // Transverse and longitudinal diffusion coefficients (Xe/CO2)
  fDiffusionT     = GetDiffusionT(fDriftVelocity,fField);
  fDiffusionL     = GetDiffusionL(fDriftVelocity,fField);

  // omega * tau = tan(Lorentz-angle)
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
AliTRDpadPlane *AliTRDparameter::GetPadPlane(Int_t p, Int_t c) const
{
  //
  // Returns the pad plane for a given plane <p> and chamber <c> number
  //

  Int_t ipp = p + c * AliTRDgeometry::Nplan();
  return ((AliTRDpadPlane *) fPadPlaneArray->At(ipp));

}

//_____________________________________________________________________________
Int_t AliTRDparameter::GetRowMax(Int_t p, Int_t c, Int_t /*s*/) const
{
  //
  // Returns the number of rows on the pad plane
  //

  return GetPadPlane(p,c)->GetNrows();

}

//_____________________________________________________________________________
Int_t AliTRDparameter::GetColMax(Int_t p) const
{
  //
  // Returns the number of rows on the pad plane
  //

  return GetPadPlane(p,0)->GetNcols();

}

//_____________________________________________________________________________
Double_t AliTRDparameter::GetRow0(Int_t p, Int_t c, Int_t /*s*/) const
{
  //
  // Returns the position of the border of the first pad in a row
  //

  return GetPadPlane(p,c)->GetRow0();

}

//_____________________________________________________________________________
Double_t AliTRDparameter::GetCol0(Int_t p) const
{
  //
  // Returns the position of the border of the first pad in a column
  //

  return GetPadPlane(p,0)->GetCol0();

}

//_____________________________________________________________________________
Double_t AliTRDparameter::CrossTalk(Double_t time) const
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
void AliTRDparameter::PrintDriftVelocity()
{
  //
  // Prints the used drift velocity
  //

  printf("<AliTRDparameter::PrintDriftVelocity> Driftvelocity = %.3f\n"
        ,fDriftVelocity);

}

//_____________________________________________________________________________
Double_t AliTRDparameter::TimeStruct(Double_t dist, Double_t z) const
{
  //
  // Applies the time structure of the drift cells (by C.Lippmann).
  // The drift time of electrons to the anode wires depends on the
  // distance to the wire (z) and on the position in the drift region.
  // 
  // input :
  // dist = radial distance from (cathode) pad plane [cm]
  // z    = distance from anode wire (parallel to cathode planes) [cm]
  //
  // output :
  // tdrift = the drift time of an electron at the given position
  //
  // We interpolate between the drift time values at the two drift
  // velocities fVDlo and fVDhi, being smaller and larger than
  // fDriftVelocity. We use the two stored drift time maps fTimeStruct1
  // and fTimeStruct2, calculated for the two mentioned drift velocities.
  //

  // indices:
  Int_t r1 = (Int_t)(10*dist);
  Int_t r2 = r1+1;
  if (r1<0)  r1 = 0;
  if (r1>37) r1 = 37;
  const Int_t z1 = (Int_t)(100*z/2.5);
  const Int_t z2 = z1+1;

  if (r1<0 || r1>37 || z1<0 || z1>10) {
    printf("<AliTRDparameter::TimeStruct> Warning. Indices out of range: ");
    printf("dist=%.2f, z=%.2f, r1=%d, z1=%d\n",dist,z,r1,z1);
  }

  const Float_t y111 = fTimeStruct1[r1+38*z1];
  const Float_t y221 = (r2 <= 37 && z2 <= 10) ? fTimeStruct1[r2+38*z2] : fTimeStruct1[37+38*10];
  const Float_t y121 = (z2 <= 10)             ? fTimeStruct1[r1+38*z2] : fTimeStruct1[r1+38*10];
  const Float_t y211 = (r2 <= 37)             ? fTimeStruct1[r2+38*z1] : fTimeStruct1[37+38*z1];

  // 2D Interpolation, lower drift time map
  const Float_t y11  = (y211-y111)*10*dist + y111 - (y211-y111)*r1;
  const Float_t y21  = (y221-y121)*10*dist + y121 - (y221-y121)*r1;

  const Float_t y112 = fTimeStruct2[r1+38*z1];
  const Float_t y222 = (r2 <= 37 && z2 <= 10) ? fTimeStruct2[r2+38*z2] : fTimeStruct2[37+38*10];
  const Float_t y122 = (z2 <= 10)             ? fTimeStruct2[r1+38*z2] : fTimeStruct2[r1+38*10];
  const Float_t y212 = (r2 <= 37)             ? fTimeStruct2[r2+38*z1] : fTimeStruct2[37+38*z1];

  // 2D Interpolation, larger drift time map
  const Float_t y12  = (y212-y112)*10*dist + y112 - (y212-y112)*r1;
  const Float_t y22  = (y222-y122)*10*dist + y122 - (y222-y122)*r1;

  // dist now is the drift distance to the anode wires (negative if electrons are
  // between anode wire plane and cathode pad plane)
  dist -= AliTRDgeometry::AmThick()/2.0;

  // Get the drift times for the drift velocities fVDlo and fVDhi
  const Float_t tdrift1 =
    ( TMath::Abs(dist)>0.005 || z>0.005 ) ? (y21-y11)*100*z/2.5+y11-(y21-y11)*z1 : 0.0;
  const Float_t tdrift2 =
    ( TMath::Abs(dist)>0.005 || z>0.005 ) ? (y22-y12)*100*z/2.5+y12-(y22-y12)*z1 : 0.0;

  // 1D Interpolation between the values at fVDlo and fVDhi
  Float_t a = (tdrift2 - tdrift1) / (fVDhi - fVDlo);
  Float_t b = tdrift2 - a * fVDhi;


  //printf("(%.2f, %.2f): %f, %f -> %f\n",
  //	 dist+AliTRDgeometry::AmThick()/2.0, z, tdrift1, tdrift2, a*fDriftVelocity+b);

  return a * fDriftVelocity + b;

}

//_____________________________________________________________________________
Int_t AliTRDparameter::Diffusion(Double_t driftlength, Double_t *xyz)
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
Int_t AliTRDparameter::ExB(Double_t driftlength, Double_t *xyz) const
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
Int_t AliTRDparameter::PadResponse(Double_t signal, Double_t dist
                                 , Int_t plane, Double_t *pad) const
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
Double_t AliTRDparameter::TimeResponse(Double_t time) const
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
void AliTRDparameter::SampleTimeStruct()
{

  //
  // Samples the timing structure of a drift cell
  // Drift Time data calculated with Garfield (by C.Lippmann)
  //

  // drift time maps are saved for some drift velocity values (in drift region):
  Float_t  fVDsmp[8];
  fVDsmp[0] = 1.032;
  fVDsmp[1] = 1.158;
  fVDsmp[2] = 1.299;
  fVDsmp[3] = 1.450;
  fVDsmp[4] = 1.610;
  fVDsmp[5] = 1.783;
  fVDsmp[6] = 1.959;
  fVDsmp[7] = 2.134;

  if ( fDriftVelocity < fVDsmp[0] ) {
    printf("<AliTRDparameter::SampleTimeStruct> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("<AliTRDparameter::SampleTimeStruct> Drift Velocity too small (%.3f<%.3f)\n"
          , fDriftVelocity, fVDsmp[0]);
    printf("<AliTRDparameter::SampleTimeStruct> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fDriftVelocity = fVDsmp[0];
  } else if ( fDriftVelocity > fVDsmp[7] ) {
    printf("<AliTRDparameter::SampleTimeStruct> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("<AliTRDparameter::SampleTimeStruct> Drift Velocity too large (%.3f>%.3f)\n"
          , fDriftVelocity,fVDsmp[6]);
    printf("<AliTRDparameter::SampleTimeStruct> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fDriftVelocity = fVDsmp[7];
  }

  const Int_t ktimebin  = 38;
  const Int_t kZbin = 11;

  // Garfield simulation at UD = -1500V; vd = 1.032cm/microsec, <driftfield> = 525V/cm
  Float_t time1500[ktimebin][kZbin] =
      {{0.09170, 0.09205, 0.09306, 0.09475, 0.09716, 0.10035,
	0.10445, 0.10993, 0.11838, 0.13986, 0.37858},
       {0.06588, 0.06626, 0.06739, 0.06926, 0.07186, 0.07524,
	0.07951, 0.08515, 0.09381, 0.11601, 0.35673},
       {0.03946, 0.04003, 0.04171, 0.04435, 0.04780, 0.05193,
	0.05680, 0.06306, 0.07290, 0.09873, 0.34748},
       {0.01151, 0.01283, 0.01718, 0.02282, 0.02880, 0.03479,
	0.04098, 0.04910, 0.06413, 0.10567, 0.36897},
       {0.01116, 0.01290, 0.01721, 0.02299, 0.02863, 0.03447,
	0.04074, 0.04984, 0.06839, 0.11625, 0.37277},
       {0.03919, 0.03974, 0.04131, 0.04380, 0.04703, 0.05102,
	0.05602, 0.06309, 0.07651, 0.10938, 0.36838},
       {0.06493, 0.06560, 0.06640, 0.06802, 0.07051, 0.07392,
	0.07853, 0.08510, 0.09690, 0.12621, 0.38058},
       {0.09174, 0.09186, 0.09225, 0.09303, 0.09477, 0.00000,
	0.11205, 0.11952, 0.13461, 0.16984, 0.43017},
       {0.14356, 0.14494, 0.14959, 0.16002, 0.18328, 0.27981,
	0.22785, 0.21240, 0.21948, 0.25965, 0.52392},
       {0.23120, 0.23366, 0.24046, 0.25422, 0.28071, 0.36914,
	0.32999, 0.31208, 0.31772, 0.35804, 0.62249},
       {0.32686, 0.32916, 0.33646, 0.35053, 0.37710, 0.46292,
	0.42773, 0.40948, 0.41497, 0.45527, 0.71955},
       {0.42353, 0.42583, 0.43317, 0.44727, 0.47380, 0.55884,
	0.52479, 0.50650, 0.51194, 0.55225, 0.81658},
       {0.52038, 0.52271, 0.53000, 0.54415, 0.57064, 0.65545,
	0.62172, 0.60341, 0.60885, 0.64915, 0.91339},
       {0.61724, 0.61953, 0.62694, 0.64098, 0.66756, 0.75226,
	0.71862, 0.70030, 0.70575, 0.74604, 1.01035},
       {0.71460, 0.71678, 0.72376, 0.73786, 0.76447, 0.84913,
	0.81551, 0.79720, 0.80264, 0.84292, 1.10723},
       {0.81101, 0.81334, 0.82066, 0.83475, 0.86127, 0.94599,
	0.91240, 0.89408, 0.89952, 0.93981, 1.20409},
       {0.90788, 0.91023, 0.91752, 0.93163, 0.95815, 1.04293,
	1.00929, 0.99096, 0.99640, 1.03669, 1.30106},
       {1.00477, 1.00707, 1.01449, 1.02852, 1.05504, 1.13976,
	1.10617, 1.08784, 1.09329, 1.13358, 1.39796},
       {1.10166, 1.10397, 1.11130, 1.12541, 1.15257, 1.23672,
	1.20307, 1.18472, 1.19018, 1.23046, 1.49486},
       {1.19854, 1.20084, 1.20818, 1.22235, 1.24885, 1.33355,
	1.29992, 1.28156, 1.28707, 1.32735, 1.59177},
       {1.29544, 1.29780, 1.30507, 1.31917, 1.34575, 1.43073,
	1.39681, 1.37851, 1.38396, 1.42377, 1.68854},
       {1.39236, 1.39462, 1.40205, 1.41607, 1.44259, 1.52745,
	1.49368, 1.47541, 1.48083, 1.52112, 1.78546},
       {1.49314, 1.49149, 1.49885, 1.51297, 1.53949, 1.62420,
	1.59016, 1.57231, 1.57772, 1.61800, 1.88048},
       {1.58610, 1.58839, 1.59572, 1.60983, 1.63635, 1.72109,
	1.68651, 1.66921, 1.67463, 1.71489, 1.97918},
       {1.68400, 1.68529, 1.69261, 1.70671, 1.73331, 1.81830,
	1.78341, 1.76605, 1.77150, 1.81179, 2.07608},
       {1.77991, 1.78215, 1.78952, 1.80385, 1.83014, 1.91486,
	1.88128, 1.86215, 1.86837, 1.90865, 2.17304},
       {1.87674, 1.87904, 1.88647, 1.90052, 1.92712, 2.01173,
	1.97812, 1.95905, 1.96527, 2.00710, 2.26979},
       {1.97369, 1.97594, 1.98326, 1.99869, 2.02442, 2.10865,
	2.07501, 2.05666, 2.06214, 2.10243, 2.36669},
       {2.07052, 2.07281, 2.08016, 2.09425, 2.12132, 2.20555,
	2.17182, 2.15341, 2.15904, 2.19933, 2.46363},
       {2.16742, 2.16971, 2.17707, 2.19114, 2.21766, 2.30240,
	2.26877, 2.25015, 2.25573, 2.29586, 2.56060},
       {2.26423, 2.26659, 2.27396, 2.28803, 2.31456, 2.40828,
	2.36567, 2.34705, 2.35282, 2.39765, 2.65744},
       {2.36153, 2.36349, 2.37330, 2.38501, 2.41159, 2.49940,
	2.46257, 2.44420, 2.44843, 2.48987, 2.75431},
       {2.46558, 2.46035, 2.46822, 2.48181, 2.50849, 2.59630,
	2.55947, 2.54112, 2.54513, 2.58677, 2.85094},
       {2.56248, 2.55723, 2.56486, 2.57871, 2.60520, 2.68998,
	2.65626, 2.63790, 2.64316, 2.68360, 2.94813},
       {2.65178, 2.65441, 2.66153, 2.67556, 2.70210, 2.78687,
	2.75319, 2.73463, 2.74032, 2.78060, 3.04503},
       {2.74868, 2.75131, 2.75870, 2.77245, 2.79385, 2.88700,
	2.85009, 2.83177, 2.83723, 2.87750, 3.14193},
       {2.84574, 2.84789, 2.85560, 2.86935, 2.89075, 2.98060,
	2.94576, 2.92868, 2.93356, 2.97436, 3.23868},
       {2.94239, 2.94469, 2.95221, 2.96625, 2.99345, 3.07747,
	3.04266, 3.02545, 3.03051, 3.07118, 3.33555}};
  
  // Garfield simulation at UD = -1600V; vd = 1.158cm/microsec, <driftfield> = 558V/cm
  Float_t time1600[ktimebin][kZbin] =
      {{0.09169, 0.09203, 0.09305, 0.09473, 0.09714, 0.10032,
	0.10441, 0.10990, 0.11835, 0.13986, 0.37845},
       {0.06589, 0.06626, 0.06738, 0.06924, 0.07184, 0.07521,
	0.07947, 0.08512, 0.09379, 0.11603, 0.35648},
       {0.03947, 0.04003, 0.04171, 0.04434, 0.04778, 0.05190,
	0.05678, 0.06304, 0.07292, 0.09876, 0.34759},
       {0.01111, 0.01281, 0.01718, 0.02281, 0.02879, 0.03477,
	0.04097, 0.04910, 0.06415, 0.10573, 0.36896},
       {0.01120, 0.01311, 0.01721, 0.02279, 0.02862, 0.03446,
	0.04074, 0.04981, 0.06825, 0.11595, 0.37255},
       {0.03919, 0.03980, 0.04132, 0.04380, 0.04704, 0.05102,
	0.05602, 0.06302, 0.07633, 0.10896, 0.36743},
       {0.06531, 0.06528, 0.06631, 0.06805, 0.07053, 0.07392,
	0.07853, 0.08505, 0.09669, 0.12578, 0.37967},
       {0.09157, 0.09171, 0.09216, 0.09301, 0.09475, 0.00000,
	0.11152, 0.11879, 0.13352, 0.16802, 0.42750},
       {0.13977, 0.14095, 0.14509, 0.15433, 0.17534, 0.26406,
	0.21660, 0.20345, 0.21113, 0.25067, 0.51434},
       {0.21816, 0.22041, 0.22631, 0.23850, 0.26208, 0.34340,
	0.30755, 0.29237, 0.29878, 0.33863, 0.60258},
       {0.30344, 0.30547, 0.31241, 0.32444, 0.34809, 0.42696,
	0.39464, 0.37919, 0.38546, 0.42530, 0.68926},
       {0.38969, 0.39164, 0.39810, 0.41059, 0.43441, 0.51246,
	0.48112, 0.46562, 0.47191, 0.51172, 0.77558},
       {0.47592, 0.47799, 0.48442, 0.49689, 0.52061, 0.59855,
	0.56752, 0.55201, 0.55826, 0.59808, 0.86202},
       {0.56226, 0.56428, 0.57074, 0.58324, 0.60696, 0.68483,
	0.65388, 0.63837, 0.64461, 0.68445, 0.94830},
       {0.64881, 0.65063, 0.65709, 0.66958, 0.69331, 0.77117,
	0.74023, 0.72472, 0.73098, 0.77079, 1.03486},
       {0.73506, 0.73698, 0.74344, 0.75596, 0.77964, 0.85751,
	0.82658, 0.81107, 0.81731, 0.85712, 1.12106},
       {0.82132, 0.82333, 0.82979, 0.84228, 0.86608, 0.94386,
	0.91293, 0.89742, 0.90367, 0.94335, 1.20737},
       {0.90767, 0.90968, 0.91614, 0.92863, 0.95236, 1.03021,
	0.99928, 0.98377, 0.99001, 1.02984, 1.29371},
       {0.99410, 0.99602, 1.00257, 1.01498, 1.03869, 1.11720,
	1.08563, 1.07011, 1.07637, 1.11621, 1.37873},
       {1.08036, 1.08240, 1.08884, 1.10138, 1.12504, 1.20301,
	1.17198, 1.15647, 1.16272, 1.20255, 1.46651},
       {1.16670, 1.16872, 1.17525, 1.18783, 1.21139, 1.28934,
	1.25833, 1.24281, 1.24909, 1.28889, 1.55275},
       {1.25307, 1.25510, 1.26153, 1.27404, 1.29773, 1.37584,
	1.34469, 1.32916, 1.33536, 1.37524, 1.63915},
       {1.33942, 1.34146, 1.34788, 1.36040, 1.38410, 1.46438,
	1.43105, 1.41537, 1.42176, 1.46158, 1.72538},
       {1.42579, 1.42782, 1.43458, 1.44674, 1.47043, 1.55085,
	1.51675, 1.50168, 1.50810, 1.54793, 1.81174},
       {1.51207, 1.51454, 1.52060, 1.53307, 1.55684, 1.63478,
	1.60336, 1.58820, 1.59446, 1.63414, 1.89814},
       {1.59856, 1.60047, 1.60693, 1.61942, 1.64317, 1.72257,
	1.69008, 1.67454, 1.68080, 1.72063, 1.98433},
       {1.68481, 1.68682, 1.69330, 1.70584, 1.72949, 1.80752,
	1.77643, 1.76089, 1.76716, 1.80692, 2.07069},
       {1.77117, 1.77319, 1.77969, 1.79260, 1.81583, 1.89376,
	1.86226, 1.84720, 1.85355, 1.89256, 2.15343},
       {1.85748, 1.85967, 1.86605, 1.87848, 1.90222, 1.98010,
	1.94913, 1.93271, 1.93981, 1.97968, 2.24355},
       {1.94386, 1.94587, 1.95233, 1.96484, 1.98854, 2.06646,
	2.03542, 2.01755, 2.02617, 2.06604, 2.32993},
       {2.03022, 2.03230, 2.03868, 2.05134, 2.07488, 2.15367,
	2.12178, 2.10391, 2.11252, 2.15432, 2.41623},
       {2.11656, 2.11857, 2.12505, 2.13772, 2.16147, 2.23919,
	2.20817, 2.19265, 2.20744, 2.23872, 2.49996},
       {2.20291, 2.20611, 2.21137, 2.22387, 2.24758, 2.32563,
	2.29450, 2.27901, 2.28525, 2.32507, 2.58897},
       {2.28922, 2.29172, 2.29774, 2.31345, 2.33400, 2.41287,
	2.38086, 2.36535, 2.37160, 2.40869, 2.67113},
       {2.37572, 2.37764, 2.38410, 2.39803, 2.42046, 2.49817,
	2.46721, 2.45171, 2.45794, 2.49505, 2.76061},
       {2.46190, 2.46396, 2.47043, 2.48340, 2.50665, 2.58453,
	2.55357, 2.53728, 2.54430, 2.58407, 2.84816},
       {2.54833, 2.55032, 2.55679, 2.56976, 2.59312, 2.67103,
	2.63993, 2.62364, 2.63062, 2.67040, 2.93444},
       {2.63456, 2.63660, 2.64304, 2.65555, 2.67938, 2.75739,
	2.72629, 2.71064, 2.71688, 2.75671, 3.01886}};
  
  // Garfield simulation at UD = -1700V; vd = 1.299cm/microsec, <driftfield> = 590V/cm
  Float_t time1700[ktimebin][kZbin] =
      {{0.09167, 0.09201, 0.09302, 0.09471, 0.09712, 0.10029,
	0.10438, 0.10986, 0.11832, 0.13986, 0.37824},
       {0.06591, 0.06626, 0.06736, 0.06923, 0.07183, 0.07519,
	0.07944, 0.08511, 0.09378, 0.11603, 0.35625},
       {0.03946, 0.04003, 0.04170, 0.04433, 0.04777, 0.05189,
	0.05676, 0.06301, 0.07291, 0.09880, 0.34724},
       {0.01110, 0.01281, 0.01718, 0.02280, 0.02879, 0.03476,
	0.04096, 0.04910, 0.06417, 0.10582, 0.36861},
       {0.01115, 0.01294, 0.01721, 0.02276, 0.02862, 0.03447,
	0.04074, 0.04980, 0.06812, 0.11565, 0.37231},
       {0.03920, 0.03974, 0.04133, 0.04381, 0.04706, 0.05105,
	0.05603, 0.06299, 0.07618, 0.10860, 0.36646},
       {0.06498, 0.06529, 0.06634, 0.06808, 0.07055, 0.07395,
	0.07852, 0.08500, 0.09650, 0.12532, 0.37850},
       {0.09143, 0.09159, 0.09207, 0.09297, 0.09473, 0.00000,
	0.11102, 0.11809, 0.13245, 0.16627, 0.42496},
       {0.13646, 0.13750, 0.14112, 0.14926, 0.16806, 0.24960,
	0.20627, 0.19536, 0.20366, 0.24256, 0.50557},
       {0.20678, 0.20848, 0.21384, 0.22450, 0.24552, 0.32018,
	0.28740, 0.27477, 0.28196, 0.32128, 0.58475},
       {0.28287, 0.28461, 0.29020, 0.30108, 0.32224, 0.39467,
	0.36500, 0.35217, 0.35926, 0.39860, 0.66194},
       {0.35972, 0.36145, 0.36713, 0.37797, 0.39912, 0.47091,
	0.44212, 0.42925, 0.43632, 0.47563, 0.73892},
       {0.43667, 0.43841, 0.44413, 0.45494, 0.47607, 0.54780,
	0.51912, 0.50627, 0.51334, 0.55254, 0.81595},
       {0.51365, 0.51540, 0.52101, 0.53193, 0.55305, 0.62463,
	0.59617, 0.58328, 0.59035, 0.62965, 0.89303},
       {0.59064, 0.59240, 0.59801, 0.60893, 0.63009, 0.70169,
	0.67317, 0.66028, 0.66735, 0.70666, 0.96995},
       {0.66765, 0.66939, 0.67501, 0.68592, 0.70724, 0.77863,
	0.75016, 0.73728, 0.74435, 0.78366, 1.04696},
       {0.74464, 0.74636, 0.75200, 0.76293, 0.78405, 0.85561,
	0.82716, 0.81427, 0.82133, 0.86064, 1.12396},
       {0.82165, 0.82340, 0.82902, 0.83991, 0.86104, 0.93266,
	0.90414, 0.89128, 0.89834, 0.93763, 1.20100},
       {0.89863, 0.90042, 0.90659, 0.91705, 0.93805, 1.00960,
	0.98115, 0.96825, 0.97533, 1.01462, 1.27801},
       {0.97563, 0.97740, 0.98310, 0.99391, 1.01504, 1.08659,
	1.05814, 1.04526, 1.05233, 1.09163, 1.35503},
       {1.05276, 1.05451, 1.06002, 1.07090, 1.09099, 1.16357,
	1.13516, 1.12225, 1.12933, 1.16863, 1.43195},
       {1.12977, 1.13138, 1.13700, 1.14792, 1.16797, 1.24061,
	1.21212, 1.19926, 1.20626, 1.24554, 1.50900},
       {1.20664, 1.20839, 1.21400, 1.22490, 1.24606, 1.31772,
	1.28914, 1.27382, 1.28329, 1.32262, 1.58550},
       {1.28367, 1.28541, 1.29099, 1.30189, 1.32312, 1.39460,
	1.36612, 1.34924, 1.36030, 1.39961, 1.66310},
       {1.36064, 1.36249, 1.36799, 1.37896, 1.40004, 1.48030,
	1.44314, 1.43032, 1.43731, 1.47659, 1.73442},
       {1.43762, 1.43937, 1.44497, 1.45618, 1.47704, 1.54932,
	1.52012, 1.50725, 1.51430, 1.55357, 1.81708},
       {1.51462, 1.51937, 1.52203, 1.53316, 1.55403, 1.62572,
	1.59713, 1.58424, 1.59128, 1.63061, 1.89406},
       {1.59162, 1.59338, 1.59947, 1.60989, 1.63103, 1.70270,
	1.67411, 1.66124, 1.66799, 1.70759, 1.97103},
       {1.66874, 1.67037, 1.67597, 1.68687, 1.70814, 1.77969,
	1.75112, 1.73806, 1.74530, 1.78457, 2.04794},
       {1.74693, 1.74749, 1.75297, 1.76476, 1.78500, 1.85667,
	1.82811, 1.81504, 1.82101, 1.86161, 2.12492},
       {1.82260, 1.82437, 1.82995, 1.84174, 1.86202, 1.93372,
	1.90509, 1.89202, 1.89930, 1.93859, 2.20189},
       {1.89964, 1.90135, 1.90693, 1.91789, 1.93952, 2.01080,
	1.98207, 1.96921, 1.97628, 2.01384, 2.27887},
       {1.97662, 1.97917, 1.98611, 1.99487, 2.01601, 2.08778,
	2.05846, 2.04623, 2.05330, 2.09244, 2.35585},
       {2.05359, 2.05615, 2.06309, 2.07187, 2.09867, 2.16459,
	2.13610, 2.12322, 2.13029, 2.16942, 2.43199},
       {2.13063, 2.13233, 2.13795, 2.14886, 2.17008, 2.24199,
	2.21310, 2.20020, 2.20727, 2.24659, 2.50983},
       {2.20761, 2.20931, 2.21955, 2.22624, 2.24708, 2.32147,
	2.29009, 2.27725, 2.28276, 2.32359, 2.58680},
       {2.28459, 2.29108, 2.29202, 2.30286, 2.32007, 2.39559,
	2.36683, 2.35422, 2.36119, 2.40058, 2.66081},
       {2.36153, 2.36806, 2.36889, 2.37985, 2.40092, 2.47828,
	2.44381, 2.43099, 2.43819, 2.47750, 2.73779}};
  
  // Garfield simulation at UD = -1800V; vd = 1.450cm/microsec, <driftfield> = 623V/cm
  Float_t time1800[ktimebin][kZbin] =
      {{0.09166, 0.09199, 0.09300, 0.09470, 0.09709, 0.10026,
	0.10434, 0.10983, 0.11831, 0.13987, 0.37802},
       {0.06585, 0.06623, 0.06735, 0.06921, 0.07180, 0.07520,
	0.07941, 0.08507, 0.09376, 0.11604, 0.35624},
       {0.03945, 0.04004, 0.04169, 0.04432, 0.04776, 0.05187,
	0.05674, 0.06300, 0.07290, 0.09884, 0.34704},
       {0.01108, 0.01287, 0.01717, 0.02280, 0.02880, 0.03476,
	0.04095, 0.04909, 0.06419, 0.10589, 0.36843},
       {0.01115, 0.01287, 0.01720, 0.02276, 0.02862, 0.03448,
	0.04073, 0.04973, 0.06799, 0.11535, 0.37224},
       {0.03918, 0.03975, 0.04134, 0.04382, 0.04707, 0.05105,
	0.05603, 0.06296, 0.07605, 0.10822, 0.36560},
       {0.06498, 0.06532, 0.06635, 0.06809, 0.07058, 0.07395,
	0.07855, 0.08495, 0.09632, 0.12488, 0.37730},
       {0.09130, 0.09160, 0.09199, 0.09300, 0.09472, 0.00000,
	0.11059, 0.11747, 0.13146, 0.16462, 0.42233},
       {0.13364, 0.13449, 0.13767, 0.14481, 0.16147, 0.23635,
	0.19706, 0.18812, 0.19704, 0.23520, 0.49749},
       {0.19698, 0.19844, 0.20311, 0.21236, 0.23082, 0.29920,
	0.26936, 0.25927, 0.26732, 0.30601, 0.56871},
       {0.26518, 0.26670, 0.27160, 0.28099, 0.29955, 0.36597,
	0.33885, 0.32858, 0.33653, 0.37524, 0.63801},
       {0.33441, 0.33553, 0.34040, 0.34987, 0.36841, 0.43428,
	0.40797, 0.39763, 0.40556, 0.44425, 0.70698},
       {0.40296, 0.40447, 0.40934, 0.41881, 0.43737, 0.50306,
	0.47695, 0.46662, 0.47455, 0.51329, 0.77600},
       {0.47296, 0.47344, 0.47830, 0.48779, 0.50632, 0.57200,
	0.54593, 0.53559, 0.54351, 0.58222, 0.84489},
       {0.54089, 0.54264, 0.54727, 0.55673, 0.57529, 0.64094,
	0.61490, 0.60457, 0.61249, 0.65118, 0.91394},
       {0.60987, 0.61138, 0.61624, 0.62573, 0.64428, 0.70989,
	0.68397, 0.67354, 0.68147, 0.72015, 0.98291},
       {0.67883, 0.68035, 0.68521, 0.69469, 0.71324, 0.77896,
	0.75287, 0.74251, 0.75043, 0.78912, 1.04458},
       {0.74780, 0.74932, 0.75421, 0.76367, 0.78221, 0.84785,
	0.82185, 0.81148, 0.81941, 0.85811, 1.12085},
       {0.81690, 0.81830, 0.82316, 0.83263, 0.85120, 0.91683,
	0.89077, 0.88045, 0.88837, 0.92707, 1.18976},
       {0.88574, 0.88726, 0.89228, 0.90198, 0.92017, 0.98578,
	0.95974, 0.94947, 0.95734, 0.99604, 1.25873},
       {0.95493, 0.95624, 0.96110, 0.97058, 0.98913, 1.05481,
	1.02873, 1.01839, 1.02631, 1.06503, 1.32772},
       {1.02392, 1.02524, 1.03008, 1.03955, 1.05810, 1.12378,
	1.09757, 1.08605, 1.09530, 1.13399, 1.39669},
       {1.09270, 1.09418, 1.09911, 1.10854, 1.12714, 1.19281,
	1.16502, 1.15633, 1.16427, 1.20271, 1.46574},
       {1.16168, 1.16323, 1.16801, 1.17772, 1.19604, 1.26190,
	1.23399, 1.22531, 1.23323, 1.27194, 1.53475},
       {1.23061, 1.23214, 1.23698, 1.24669, 1.26503, 1.33073,
	1.30461, 1.29428, 1.30220, 1.34091, 1.60372},
       {1.29960, 1.30110, 1.30596, 1.31544, 1.33398, 1.39962,
	1.37228, 1.36323, 1.37121, 1.40988, 1.67273},
       {1.36851, 1.37007, 1.37512, 1.38441, 1.40297, 1.46865,
	1.44256, 1.43222, 1.44017, 1.47878, 1.74155},
       {1.43752, 1.43904, 1.44773, 1.45338, 1.47220, 1.53759,
	1.51136, 1.50119, 1.50914, 1.54775, 1.81050},
       {1.50646, 1.50802, 1.51288, 1.52237, 1.54097, 1.60697,
	1.58049, 1.57018, 1.57811, 1.61678, 1.87947},
       {1.57545, 1.57720, 1.58185, 1.59134, 1.60996, 1.67787,
	1.64929, 1.63914, 1.64707, 1.68570, 1.94851},
       {1.64442, 1.64617, 1.65081, 1.66035, 1.67893, 1.74684,
	1.71826, 1.70745, 1.71604, 1.75310, 2.01748},
       {1.71337, 1.71513, 1.71978, 1.72932, 1.74645, 1.81346,
	1.78739, 1.77642, 1.78501, 1.82151, 2.08644},
       {1.78238, 1.78410, 1.78876, 1.79824, 1.81678, 1.88332,
	1.85639, 1.84262, 1.85397, 1.89270, 2.15533},
       {1.85135, 1.85306, 1.85778, 1.86728, 1.88580, 1.95615,
	1.92536, 1.91171, 1.92283, 1.96165, 2.22428},
       {1.92774, 1.92184, 1.92672, 1.93618, 1.95477, 2.02048,
	1.99427, 1.98068, 1.99192, 2.03062, 2.29338},
       {1.98929, 1.99081, 1.99567, 2.00515, 2.02373, 2.08987,
	2.06332, 2.05249, 2.05939, 2.09928, 2.36227},
       {2.05829, 2.05978, 2.06464, 2.07414, 2.09272, 2.15850,
	2.12928, 2.12194, 2.12987, 2.16825, 2.43083},
       {2.12726, 2.12869, 2.13360, 2.14425, 2.16160, 2.22872,
	2.20118, 2.19078, 2.19876, 2.23416, 2.50007}};
  
  // Garfield simulation at UD = -1900V; vd = 1.610cm/microsec, <driftfield> = 655V/cm
  Float_t time1900[ktimebin][kZbin] =
      {{0.09166, 0.09198, 0.09298, 0.09467, 0.09707, 0.10023,
	0.10431, 0.10980, 0.11828, 0.13988, 0.37789},
       {0.06584, 0.06622, 0.06735, 0.06920, 0.07179, 0.07514,
	0.07938, 0.08505, 0.09374, 0.11606, 0.35599},
       {0.03945, 0.04002, 0.04169, 0.04432, 0.04775, 0.05185,
	0.05672, 0.06298, 0.07290, 0.09888, 0.34695},
       {0.01109, 0.01281, 0.01717, 0.02279, 0.02878, 0.03476,
	0.04094, 0.04909, 0.06421, 0.10597, 0.36823},
       {0.01115, 0.01287, 0.01720, 0.02294, 0.02862, 0.03448,
	0.04074, 0.04973, 0.06783, 0.11506, 0.37206},
       {0.03940, 0.03975, 0.04135, 0.04386, 0.04708, 0.05106,
	0.05604, 0.06293, 0.07592, 0.10787, 0.36484},
       {0.06500, 0.06534, 0.06638, 0.06811, 0.07060, 0.07413,
	0.07852, 0.08491, 0.09614, 0.12446, 0.37626},
       {0.09119, 0.09140, 0.09194, 0.09293, 0.09471, 0.00000,
	0.11013, 0.11685, 0.13050, 0.16302, 0.41991},
       {0.13111, 0.13190, 0.13466, 0.14091, 0.15554, 0.22409,
	0.18846, 0.18167, 0.19113, 0.22854, 0.48995},
       {0.18849, 0.18975, 0.19380, 0.20185, 0.21797, 0.28050,
	0.25368, 0.24575, 0.25446, 0.29249, 0.55442},
       {0.24995, 0.25125, 0.25563, 0.26366, 0.27986, 0.34065,
	0.31605, 0.30815, 0.31680, 0.35482, 0.61697},
       {0.31187, 0.31324, 0.31745, 0.32580, 0.34205, 0.40217,
	0.37825, 0.37031, 0.37897, 0.41696, 0.67890},
       {0.37401, 0.37531, 0.37955, 0.38777, 0.40395, 0.46408,
	0.44037, 0.43242, 0.44108, 0.47906, 0.74122},
       {0.43610, 0.43741, 0.44161, 0.44986, 0.46604, 0.52614,
	0.50248, 0.49452, 0.50316, 0.54116, 0.80326},
       {0.49820, 0.49988, 0.50372, 0.51196, 0.52814, 0.58822,
	0.56459, 0.55661, 0.56527, 0.60326, 0.86526},
       {0.56032, 0.56161, 0.56582, 0.57408, 0.59024, 0.65032,
	0.62670, 0.61872, 0.62737, 0.66537, 0.92749},
       {0.62240, 0.62371, 0.62792, 0.63624, 0.65236, 0.71241,
	0.68881, 0.68081, 0.68947, 0.72750, 0.98941},
       {0.68449, 0.68581, 0.69002, 0.69828, 0.71444, 0.77452,
	0.75089, 0.74295, 0.75157, 0.78957, 1.05157},
       {0.74660, 0.74790, 0.75212, 0.76036, 0.77654, 0.83748,
	0.81299, 0.80501, 0.81193, 0.85168, 1.11375},
       {0.80870, 0.81017, 0.81423, 0.82246, 0.83867, 0.89908,
	0.87509, 0.86660, 0.87577, 0.91376, 1.17586},
       {0.87080, 0.87233, 0.87632, 0.88458, 0.90074, 0.96083,
	0.93718, 0.92922, 0.93787, 0.97588, 1.23794},
       {0.93291, 0.93422, 0.93844, 0.94667, 0.96293, 1.02295,
	0.99929, 0.99127, 0.99997, 1.03797, 1.29995},
       {0.99500, 0.99645, 1.00308, 1.00877, 1.02497, 1.08504,
	1.06140, 1.05343, 1.06203, 1.10006, 1.36216},
       {1.05712, 1.05926, 1.06262, 1.07092, 1.08706, 1.14754,
	1.12350, 1.11550, 1.12417, 1.16218, 1.42427},
       {1.11921, 1.12059, 1.12473, 1.13297, 1.14916, 1.21140,
	1.18560, 1.17284, 1.18625, 1.22414, 1.48629},
       {1.18140, 1.18262, 1.18690, 1.19508, 1.21125, 1.27139,
	1.24164, 1.23495, 1.24838, 1.28634, 1.54852},
       {1.24340, 1.24473, 1.24901, 1.25732, 1.27336, 1.33358,
	1.30793, 1.30179, 1.31047, 1.34848, 1.61066},
       {1.30551, 1.30684, 1.31104, 1.32056, 1.33553, 1.39609,
	1.37004, 1.36392, 1.37045, 1.41057, 1.67259},
       {1.36755, 1.36892, 1.37315, 1.39148, 1.39755, 1.45820,
	1.43215, 1.42602, 1.43467, 1.47268, 1.73477},
       {1.42966, 1.43101, 1.43549, 1.45359, 1.45976, 1.52031,
	1.49601, 1.48811, 1.49677, 1.53477, 1.79691},
       {1.49180, 1.49321, 1.49760, 1.51570, 1.52175, 1.58185,
	1.55771, 1.55023, 1.55888, 1.59672, 1.85501},
       {1.55391, 1.55527, 1.55943, 1.57782, 1.58391, 1.64395,
	1.62008, 1.61233, 1.62085, 1.65883, 1.92091},
       {1.61599, 1.61732, 1.62154, 1.63993, 1.64612, 1.70608,
	1.68237, 1.67108, 1.68301, 1.72110, 1.97931},
       {1.67808, 1.67948, 1.68364, 1.70204, 1.70823, 1.76858,
	1.74404, 1.73539, 1.74512, 1.78321, 2.04522},
       {1.74019, 1.74152, 1.74573, 1.76415, 1.77015, 1.83040,
	1.80615, 1.79366, 1.80723, 1.84509, 2.10742},
       {1.80235, 1.80362, 1.80783, 1.82626, 1.83227, 1.89246,
	1.86795, 1.85405, 1.86938, 1.90720, 2.16953},
       {1.86442, 1.86572, 1.86994, 1.88837, 1.89438, 1.95445,
	1.93006, 1.92283, 1.93148, 1.96931, 2.23147},
       {1.92700, 1.92783, 1.93197, 1.95049, 1.95649, 2.01668,
	1.99217, 1.98486, 1.99352, 2.03143, 2.29358}};

  // Garfield simulation at UD = -2000V; vd = 1.783cm/microsec, <driftfield> = 688V/cm
  Float_t time2000[ktimebin][kZbin] =
      {{0.09176, 0.09196, 0.09296, 0.09465, 0.09704, 0.10020,
	0.10427, 0.10977, 0.11825, 0.13988, 0.37774},
       {0.06583, 0.06620, 0.06735, 0.06918, 0.07177, 0.07513,
	0.07936, 0.08503, 0.09372, 0.11606, 0.35586},
       {0.03944, 0.04001, 0.04170, 0.04431, 0.04774, 0.05184,
	0.05670, 0.06296, 0.07291, 0.09893, 0.34680},
       {0.01108, 0.01281, 0.01719, 0.02279, 0.02879, 0.03474,
	0.04093, 0.04908, 0.06422, 0.10605, 0.36800},
       {0.01114, 0.01287, 0.01720, 0.02276, 0.02863, 0.03449,
	0.04073, 0.04970, 0.06774, 0.11478, 0.37179},
       {0.03925, 0.03977, 0.04135, 0.04386, 0.04711, 0.05108,
	0.05604, 0.06290, 0.07580, 0.10748, 0.36386},
       {0.06501, 0.06536, 0.06640, 0.06814, 0.07062, 0.07398,
	0.07852, 0.08487, 0.09598, 0.12405, 0.37519},
       {0.09109, 0.09127, 0.09188, 0.09292, 0.09472, 0.00000,
	0.10964, 0.11630, 0.12960, 0.16150, 0.41765},
       {0.12898, 0.12968, 0.13209, 0.13749, 0.15034, 0.21286,
	0.18088, 0.17590, 0.18591, 0.22254, 0.48315},
       {0.18122, 0.18227, 0.18574, 0.19263, 0.20674, 0.26376,
	0.23960, 0.23375, 0.24316, 0.28047, 0.54179},
       {0.23674, 0.23784, 0.24142, 0.24847, 0.26264, 0.31810,
	0.29602, 0.29008, 0.29944, 0.33675, 0.59795},
       {0.29279, 0.29382, 0.29742, 0.30448, 0.31865, 0.37364,
	0.35215, 0.34629, 0.35555, 0.39286, 0.65411},
       {0.34875, 0.34987, 0.35346, 0.36054, 0.37472, 0.42956,
	0.40825, 0.40229, 0.41167, 0.44894, 0.71033},
       {0.40484, 0.40594, 0.40954, 0.41660, 0.43077, 0.48560,
	0.46433, 0.45840, 0.46772, 0.50500, 0.76632},
       {0.46090, 0.46201, 0.46560, 0.47267, 0.48684, 0.54167,
	0.52041, 0.51449, 0.52382, 0.56108, 0.82227},
       {0.51698, 0.51809, 0.52173, 0.52874, 0.54291, 0.59776,
	0.57646, 0.57052, 0.57986, 0.61717, 0.87836},
       {0.57306, 0.57418, 0.57782, 0.58513, 0.59899, 0.65380,
	0.63255, 0.62661, 0.63594, 0.67325, 0.93460},
       {0.62912, 0.63024, 0.63383, 0.64103, 0.65506, 0.70988,
	0.68484, 0.68267, 0.69202, 0.72878, 0.99046},
       {0.68521, 0.68633, 0.68990, 0.69699, 0.71115, 0.76595,
	0.74468, 0.73872, 0.74814, 0.78538, 1.04674},
       {0.74127, 0.74239, 0.74605, 0.75303, 0.77022, 0.82204,
	0.80078, 0.79484, 0.80416, 0.84147, 1.10261},
       {0.79736, 0.79846, 0.80206, 0.80947, 0.82330, 0.87813,
	0.85688, 0.85087, 0.86023, 0.89752, 1.15874},
       {0.85342, 0.85454, 0.85815, 0.86519, 0.87936, 0.93417,
	0.91293, 0.90428, 0.91631, 0.95360, 1.20760},
       {0.90949, 0.91061, 0.91423, 0.92128, 0.93544, 0.99026,
	0.96807, 0.96305, 0.97239, 1.00967, 1.27078},
       {0.96556, 0.96669, 0.97111, 0.97734, 0.99151, 1.04664,
	1.02508, 1.01879, 1.02846, 1.06167, 1.32695},
       {1.02167, 1.02279, 1.02656, 1.03341, 1.04759, 1.10242,
	1.08115, 1.07003, 1.08453, 1.12184, 1.38304},
       {1.07776, 1.07883, 1.08242, 1.08950, 1.10384, 1.16422,
	1.13725, 1.13133, 1.14061, 1.17793, 1.43910},
       {1.13379, 1.13492, 1.13864, 1.14567, 1.15973, 1.21455,
	1.19323, 1.18734, 1.19668, 1.23401, 1.49528},
       {1.18988, 1.19098, 1.19457, 1.20164, 1.21582, 1.27064,
	1.24937, 1.24044, 1.25275, 1.29004, 1.55137},
       {1.24592, 1.24706, 1.25087, 1.25773, 1.27188, 1.32670,
	1.30544, 1.29953, 1.30883, 1.34613, 1.60743},
       {1.30202, 1.30313, 1.30673, 1.31381, 1.32797, 1.38278,
	1.36151, 1.35167, 1.36490, 1.40221, 1.66306},
       {1.35809, 1.35921, 1.36282, 1.36986, 1.38403, 1.43888,
	1.41760, 1.41174, 1.42083, 1.45830, 1.71915},
       {1.41419, 1.41528, 1.41890, 1.42595, 1.44011, 1.49496,
	1.47368, 1.46769, 1.47706, 1.51436, 1.77523},
       {1.47131, 1.47141, 1.47494, 1.48850, 1.49620, 1.55137,
	1.52977, 1.51820, 1.53315, 1.57042, 1.83158},
       {1.52635, 1.52750, 1.53103, 1.53814, 1.55228, 1.60736,
	1.58503, 1.57986, 1.58920, 1.62649, 1.88767},
       {1.58418, 1.58355, 1.58711, 1.59526, 1.60833, 1.66316,
	1.63345, 1.63261, 1.64556, 1.68204, 1.94359},
       {1.64027, 1.63958, 1.64489, 1.65024, 1.66443, 1.71925,
	1.69794, 1.69201, 1.70143, 1.73865, 1.99968},
       {1.69450, 1.69566, 1.69940, 1.70697, 1.71841, 1.77819,
	1.75396, 1.74814, 1.75743, 1.79083, 2.05427},
       {1.75054, 1.75221, 1.75527, 1.76306, 1.77662, 1.83428,
	1.81006, 1.81173, 1.81345, 1.85076, 2.10289}};

  // Garfield simulation at UD = -2100V; vd = 1.959cm/microsec, <driftfield> = 720V/cm
  Float_t time2100[ktimebin][kZbin] =
      {{0.09160, 0.09194, 0.09294, 0.09462, 0.09701, 0.10017,
	0.10424, 0.10974, 0.11823, 0.13988, 0.37762},
       {0.06585, 0.06619, 0.06731, 0.06916, 0.07174, 0.07509,
	0.07933, 0.08500, 0.09370, 0.11609, 0.35565},
       {0.03960, 0.04001, 0.04171, 0.04430, 0.04774, 0.05182,
	0.05668, 0.06294, 0.07291, 0.09896, 0.34676},
       {0.01109, 0.01280, 0.01716, 0.02279, 0.02876, 0.03474,
	0.04096, 0.04908, 0.06424, 0.10612, 0.36790},
       {0.01114, 0.01285, 0.01719, 0.02287, 0.02863, 0.03449,
	0.04073, 0.04964, 0.06759, 0.11446, 0.37162},
       {0.03922, 0.03977, 0.04146, 0.04386, 0.04711, 0.05109,
	0.05605, 0.06287, 0.07575, 0.10713, 0.36298},
       {0.06504, 0.06538, 0.06641, 0.06818, 0.07064, 0.07426,
	0.07852, 0.08483, 0.09581, 0.12363, 0.37424},
       {0.09103, 0.09129, 0.09186, 0.09291, 0.09476, 0.00000,
	0.10923, 0.11578, 0.12873, 0.16005, 0.41525},
       {0.12723, 0.12777, 0.12988, 0.13458, 0.14579, 0.20264,
	0.17421, 0.17078, 0.18132, 0.21708, 0.47699},
       {0.17508, 0.17601, 0.17897, 0.18487, 0.19698, 0.24881,
	0.22737, 0.22337, 0.23348, 0.27000, 0.53032},
       {0.22571, 0.22663, 0.22969, 0.23570, 0.24787, 0.29826,
	0.27871, 0.27462, 0.28471, 0.32122, 0.58166},
       {0.27664, 0.27759, 0.28067, 0.28669, 0.29891, 0.34898,
	0.32982, 0.32570, 0.33576, 0.37229, 0.63268},
       {0.32766, 0.32862, 0.33170, 0.33778, 0.34988, 0.39973,
	0.38088, 0.37675, 0.38680, 0.42333, 0.68159},
       {0.37872, 0.37966, 0.38275, 0.38875, 0.40093, 0.45073,
	0.43192, 0.42780, 0.43786, 0.47438, 0.73480},
       {0.42974, 0.43070, 0.43378, 0.43982, 0.45196, 0.50177,
	0.48297, 0.47884, 0.48889, 0.52544, 0.78581},
       {0.48081, 0.48175, 0.48482, 0.49084, 0.50302, 0.55290,
	0.53398, 0.52988, 0.53994, 0.57647, 0.83687},
       {0.53645, 0.53295, 0.53586, 0.54188, 0.55408, 0.60398,
	0.58504, 0.58092, 0.59100, 0.62768, 0.88773},
       {0.58345, 0.58409, 0.58690, 0.59292, 0.60510, 0.65562,
	0.63609, 0.63197, 0.64203, 0.67856, 0.93907},
       {0.63397, 0.63490, 0.63795, 0.64403, 0.65613, 0.70612,
	0.68714, 0.68301, 0.69294, 0.72955, 0.99000},
       {0.68496, 0.68592, 0.68899, 0.69504, 0.70733, 0.75708,
	0.73818, 0.73405, 0.74412, 0.78064, 1.04100},
       {0.73600, 0.73696, 0.74003, 0.74624, 0.75828, 0.80805,
	0.78904, 0.78512, 0.79517, 0.83152, 1.09205},
       {0.78709, 0.78801, 0.79108, 0.79709, 0.80931, 0.85906,
	0.84027, 0.83614, 0.84621, 0.88269, 1.14058},
       {0.83808, 0.83905, 0.84215, 0.84816, 0.86031, 0.91011,
	0.89139, 0.88718, 0.89725, 0.93377, 1.19413},
       {0.88916, 0.89010, 0.89320, 0.89920, 0.91136, 0.96117,
	0.94235, 0.93822, 0.94828, 0.98480, 1.24538},
       {0.94036, 0.94113, 0.94422, 0.95023, 0.96241, 1.01220,
	0.99310, 0.98927, 0.99933, 1.03585, 1.29629},
       {0.99139, 0.99220, 0.99525, 1.00127, 1.01344, 1.06324,
	1.04451, 1.04033, 1.04836, 1.08690, 1.34727},
       {1.04261, 1.04325, 1.04628, 1.05232, 1.06448, 1.12090,
	1.09546, 1.09136, 1.10142, 1.13795, 1.39831},
       {1.09331, 1.09429, 1.09742, 1.10336, 1.11557, 1.16547,
	1.14658, 1.13642, 1.15247, 1.18898, 1.44936},
       {1.14436, 1.14539, 1.14847, 1.15443, 1.16662, 1.21794,
	1.19763, 1.19329, 1.20349, 1.23956, 1.50043},
       {1.19533, 1.19651, 1.19943, 1.20548, 1.21666, 1.26753,
	1.24862, 1.24434, 1.25455, 1.29106, 1.55142},
       {1.24638, 1.24756, 1.25046, 1.25648, 1.26764, 1.31858,
	1.29967, 1.29538, 1.30499, 1.34211, 1.60250},
       {1.29747, 1.29847, 1.30175, 1.30753, 1.31869, 1.36969,
	1.35069, 1.34656, 1.35663, 1.39316, 1.64644},
       {1.35537, 1.34952, 1.35255, 1.35869, 1.36973, 1.41387,
	1.40173, 1.39761, 1.40768, 1.44396, 1.70238},
       {1.39956, 1.40056, 1.40380, 1.40961, 1.42178, 1.46492,
	1.45278, 1.45423, 1.45872, 1.49522, 1.75557},
       {1.45080, 1.45159, 1.45463, 1.46109, 1.47287, 1.52263,
	1.50382, 1.50050, 1.50977, 1.54502, 1.80670},
       {1.50165, 1.50264, 1.50570, 1.51214, 1.52233, 1.57370,
	1.55484, 1.55155, 1.56080, 1.59731, 1.85778},
       {1.55269, 1.55364, 1.55675, 1.56274, 1.57491, 1.62598,
	1.60590, 1.60259, 1.61185, 1.64836, 1.90883},
       {1.60368, 1.60469, 1.60779, 1.61373, 1.62596, 1.67738,
	1.65651, 1.65249, 1.66290, 1.69936, 1.95959}};

  // Garfield simulation at UD = -2200V; vd = 2.134cm/microsec, <driftfield> = 753V/cm
  Float_t time2200[ktimebin][kZbin] =
      {{0.09162, 0.09194, 0.09292, 0.09460, 0.09702, 0.10014,
	0.10421, 0.10971, 0.11820, 0.13990, 0.37745},
       {0.06581, 0.06618, 0.06730, 0.06915, 0.07173, 0.07507,
	0.07931, 0.08497, 0.09368, 0.11609, 0.35560},
       {0.03947, 0.04001, 0.04167, 0.04429, 0.04772, 0.05183,
	0.05667, 0.06293, 0.07292, 0.09900, 0.34673},
       {0.01111, 0.01280, 0.01716, 0.02279, 0.02876, 0.03473,
	0.04091, 0.04907, 0.06426, 0.10620, 0.36766},
       {0.01113, 0.01285, 0.01719, 0.02276, 0.02863, 0.03452,
	0.04076, 0.04960, 0.06745, 0.11419, 0.37139},
       {0.03923, 0.03978, 0.04137, 0.04387, 0.04713, 0.05110,
	0.05605, 0.06284, 0.07551, 0.10677, 0.36210},
       {0.06505, 0.06540, 0.06644, 0.06820, 0.07069, 0.07401,
	0.07852, 0.08479, 0.09565, 0.12325, 0.37313},
       {0.09107, 0.09127, 0.09181, 0.09291, 0.09474, 0.00000,
	0.10883, 0.11528, 0.12789, 0.15865, 0.41313},
       {0.12559, 0.12622, 0.12800, 0.13206, 0.14166, 0.19331,
	0.16832, 0.16632, 0.17724, 0.21218, 0.47098},
       {0.16992, 0.17070, 0.17325, 0.17831, 0.18871, 0.23557,
	0.21690, 0.21451, 0.22514, 0.26082, 0.52034},
       {0.21640, 0.21722, 0.21987, 0.22500, 0.23540, 0.28097,
	0.26399, 0.26154, 0.27214, 0.30784, 0.56734},
       {0.26318, 0.26400, 0.26679, 0.27181, 0.28220, 0.32739,
	0.31090, 0.30842, 0.31902, 0.35474, 0.61415},
       {0.31001, 0.31085, 0.31348, 0.31866, 0.32903, 0.37412,
	0.35777, 0.35546, 0.36588, 0.40159, 0.66103},
       {0.35687, 0.35769, 0.36033, 0.36556, 0.37588, 0.42094,
	0.40523, 0.40214, 0.41273, 0.44841, 0.70785},
       {0.40372, 0.40457, 0.40723, 0.41234, 0.42273, 0.46778,
	0.45148, 0.44903, 0.45961, 0.49526, 0.75486},
       {0.45062, 0.45139, 0.45404, 0.45966, 0.46958, 0.51470,
	0.49833, 0.49584, 0.50644, 0.54211, 0.80160},
       {0.49742, 0.49825, 0.50088, 0.50605, 0.51644, 0.56148,
	0.54518, 0.54270, 0.55330, 0.58897, 0.84854},
       {0.54427, 0.54510, 0.54774, 0.55290, 0.56329, 0.60846,
	0.59203, 0.58955, 0.60014, 0.63578, 0.89528},
       {0.59119, 0.59199, 0.59471, 0.59975, 0.61014, 0.65533,
	0.63889, 0.63636, 0.64699, 0.68269, 0.94197},
       {0.63866, 0.63880, 0.64145, 0.64664, 0.65701, 0.70639,
	0.68574, 0.68325, 0.69385, 0.72949, 0.98900},
       {0.68483, 0.68566, 0.68831, 0.69347, 0.70386, 0.74890,
	0.73260, 0.73010, 0.74069, 0.77638, 1.03320},
       {0.73168, 0.73255, 0.73515, 0.74031, 0.75072, 0.79576,
	0.77117, 0.77501, 0.78755, 0.82318, 1.08006},
       {0.77854, 0.78310, 0.78200, 0.79525, 0.79756, 0.84281,
	0.81803, 0.82393, 0.83441, 0.87008, 1.12692},
       {0.82541, 0.82642, 0.82916, 0.83528, 0.84442, 0.89086,
	0.87315, 0.87079, 0.88125, 0.91694, 1.17648},
       {0.87226, 0.87308, 0.87602, 0.88086, 0.89128, 0.93772,
	0.92001, 0.91751, 0.92811, 0.95587, 1.22328},
       {0.91921, 0.91994, 0.92256, 0.92772, 0.94713, 0.98566,
	0.96690, 0.96436, 0.97495, 1.01064, 1.26882},
       {0.96790, 0.96679, 0.96941, 0.97463, 0.99399, 1.03001,
	1.01376, 1.01112, 1.02181, 1.05749, 1.31568},
       {1.01278, 1.01390, 1.01674, 1.02147, 1.03182, 1.07820,
	1.06056, 1.05798, 1.06867, 1.10433, 1.36390},
       {1.05964, 1.06076, 1.06331, 1.06833, 1.07870, 1.13297,
	1.10742, 1.10520, 1.11543, 1.15120, 1.41069},
       {1.10664, 1.10762, 1.10997, 1.11519, 1.12556, 1.17531,
	1.15427, 1.14620, 1.16229, 1.19805, 1.45758},
       {1.15352, 1.15421, 1.15683, 1.16218, 1.17242, 1.21910,
	1.20035, 1.19863, 1.20579, 1.24473, 1.50412},
       {1.20019, 1.20115, 1.20369, 1.20892, 1.21928, 1.26913,
	1.24721, 1.24549, 1.25605, 1.29159, 1.54920},
       {1.24707, 1.24846, 1.25052, 1.25602, 1.26608, 1.31558,
	1.29448, 1.29232, 1.30293, 1.33675, 1.59798},
       {1.29391, 1.29475, 1.29738, 1.30255, 1.31294, 1.36244,
	1.34167, 1.33918, 1.34979, 1.38229, 1.64496},
       {1.34078, 1.34304, 1.34424, 1.35565, 1.35980, 1.40930,
	1.38853, 1.38229, 1.39664, 1.42863, 1.69162},
       {1.38762, 1.38847, 1.39110, 1.39627, 1.40666, 1.45183,
	1.43539, 1.43289, 1.44348, 1.47549, 1.73876},
       {1.43524, 1.43533, 1.43796, 1.44310, 1.45371, 1.49305,
	1.48224, 1.47941, 1.49034, 1.52601, 1.78552},
       {1.48122, 1.48219, 1.48482, 1.48991, 1.50030, 1.53991,
	1.52898, 1.52653, 1.53653, 1.57282, 1.82386}};

  if (fTimeStruct1)  delete [] fTimeStruct1;
  fTimeStruct1  = new Float_t[ktimebin*kZbin];

  if (fTimeStruct2)  delete [] fTimeStruct2;
  fTimeStruct2  = new Float_t[ktimebin*kZbin];

  for (Int_t ctrt = 0; ctrt<ktimebin; ctrt++) {
    for (Int_t ctrz = 0; ctrz<kZbin; ctrz++) {
      if ( fDriftVelocity > fVDsmp[6] ) {
	fTimeStruct1[ctrt+ctrz*ktimebin] = time2100[ctrt][ctrz];
	fTimeStruct2[ctrt+ctrz*ktimebin] = time2200[ctrt][ctrz];	    
	fVDlo = fVDsmp[6];
	fVDhi = fVDsmp[7];
      } else if ( fDriftVelocity > fVDsmp[5] ) {
	fTimeStruct1[ctrt+ctrz*ktimebin] = time2000[ctrt][ctrz];
	fTimeStruct2[ctrt+ctrz*ktimebin] = time2100[ctrt][ctrz];	    
	fVDlo = fVDsmp[5];
	fVDhi = fVDsmp[6];
      } else if ( fDriftVelocity > fVDsmp[4] ) {
	fTimeStruct1[ctrt+ctrz*ktimebin] = time1900[ctrt][ctrz];
	fTimeStruct2[ctrt+ctrz*ktimebin] = time2000[ctrt][ctrz];	    
	fVDlo = fVDsmp[4];
	fVDhi = fVDsmp[5];
      } else if ( fDriftVelocity > fVDsmp[3] ) {
	fTimeStruct1[ctrt+ctrz*ktimebin] = time1800[ctrt][ctrz];
	fTimeStruct2[ctrt+ctrz*ktimebin] = time1900[ctrt][ctrz];	    
	fVDlo = fVDsmp[3];
	fVDhi = fVDsmp[4];
      } else if ( fDriftVelocity > fVDsmp[2] ) {
	fTimeStruct1[ctrt+ctrz*ktimebin] = time1700[ctrt][ctrz];
	fTimeStruct2[ctrt+ctrz*ktimebin] = time1800[ctrt][ctrz];	    
	fVDlo = fVDsmp[2];
	fVDhi = fVDsmp[3];
      } else if ( fDriftVelocity > fVDsmp[1] ) {
	fTimeStruct1[ctrt+ctrz*ktimebin] = time1600[ctrt][ctrz];
	fTimeStruct2[ctrt+ctrz*ktimebin] = time1700[ctrt][ctrz];	    
	fVDlo = fVDsmp[1];
	fVDhi = fVDsmp[2];
      } else if ( fDriftVelocity > (fVDsmp[0] - 1.e-5) ) {
	fTimeStruct1[ctrt+ctrz*ktimebin] = time1500[ctrt][ctrz];
	fTimeStruct2[ctrt+ctrz*ktimebin] = time1600[ctrt][ctrz];	    
	fVDlo = fVDsmp[0];
	fVDhi = fVDsmp[1];
      }
    }
  }

  //printf("fVDlo = %f, fVDhi = %f\n", fVDlo, fVDhi);

  ReInit();

}

//_____________________________________________________________________________
void AliTRDparameter::SamplePRF()
{
  //
  // Samples the pad response function
  //

  const Int_t kNplan  = AliTRDgeometry::kNplan;
  const Int_t kPRFbin = 61;

  Float_t prf[kNplan][kPRFbin] = { {2.9037e-02, 3.3608e-02, 3.9020e-02, 4.5292e-02,
                    5.2694e-02, 6.1362e-02, 7.1461e-02, 8.3362e-02,
                    9.7063e-02, 1.1307e-01, 1.3140e-01, 1.5235e-01,
                    1.7623e-01, 2.0290e-01, 2.3294e-01, 2.6586e-01,
                    3.0177e-01, 3.4028e-01, 3.8077e-01, 4.2267e-01,
                    4.6493e-01, 5.0657e-01, 5.4655e-01, 5.8397e-01,
                    6.1767e-01, 6.4744e-01, 6.7212e-01, 6.9188e-01,
                    7.0627e-01, 7.1499e-01, 7.1851e-01, 7.1499e-01,
                    7.0627e-01, 6.9188e-01, 6.7212e-01, 6.4744e-01,
                    6.1767e-01, 5.8397e-01, 5.4655e-01, 5.0657e-01,
                    4.6493e-01, 4.2267e-01, 3.8077e-01, 3.4028e-01,
                    3.0177e-01, 2.6586e-01, 2.3294e-01, 2.0290e-01,
                    1.7623e-01, 1.5235e-01, 1.3140e-01, 1.1307e-01,
                    9.7063e-02, 8.3362e-02, 7.1461e-02, 6.1362e-02,
                    5.2694e-02, 4.5292e-02, 3.9020e-02, 3.3608e-02,
                    2.9037e-02},
                   {2.5478e-02, 2.9695e-02, 3.4655e-02, 4.0454e-02,
                    4.7342e-02, 5.5487e-02, 6.5038e-02, 7.6378e-02,
                    8.9696e-02, 1.0516e-01, 1.2327e-01, 1.4415e-01,
                    1.6794e-01, 1.9516e-01, 2.2573e-01, 2.5959e-01,
                    2.9694e-01, 3.3719e-01, 3.7978e-01, 4.2407e-01,
                    4.6889e-01, 5.1322e-01, 5.5569e-01, 5.9535e-01,
                    6.3141e-01, 6.6259e-01, 6.8882e-01, 7.0983e-01,
                    7.2471e-01, 7.3398e-01, 7.3761e-01, 7.3398e-01,
                    7.2471e-01, 7.0983e-01, 6.8882e-01, 6.6259e-01,
                    6.3141e-01, 5.9535e-01, 5.5569e-01, 5.1322e-01,
                    4.6889e-01, 4.2407e-01, 3.7978e-01, 3.3719e-01,
                    2.9694e-01, 2.5959e-01, 2.2573e-01, 1.9516e-01,
                    1.6794e-01, 1.4415e-01, 1.2327e-01, 1.0516e-01,
                    8.9696e-02, 7.6378e-02, 6.5038e-02, 5.5487e-02,
                    4.7342e-02, 4.0454e-02, 3.4655e-02, 2.9695e-02,
                    2.5478e-02},
                   {2.2363e-02, 2.6233e-02, 3.0782e-02, 3.6140e-02,
                    4.2535e-02, 5.0157e-02, 5.9197e-02, 6.9900e-02,
                    8.2707e-02, 9.7811e-02, 1.1548e-01, 1.3601e-01,
                    1.5998e-01, 1.8739e-01, 2.1840e-01, 2.5318e-01,
                    2.9182e-01, 3.3373e-01, 3.7837e-01, 4.2498e-01,
                    4.7235e-01, 5.1918e-01, 5.6426e-01, 6.0621e-01,
                    6.4399e-01, 6.7700e-01, 7.0472e-01, 7.2637e-01,
                    7.4206e-01, 7.5179e-01, 7.5551e-01, 7.5179e-01,
                    7.4206e-01, 7.2637e-01, 7.0472e-01, 6.7700e-01,
                    6.4399e-01, 6.0621e-01, 5.6426e-01, 5.1918e-01,
                    4.7235e-01, 4.2498e-01, 3.7837e-01, 3.3373e-01,
                    2.9182e-01, 2.5318e-01, 2.1840e-01, 1.8739e-01,
                    1.5998e-01, 1.3601e-01, 1.1548e-01, 9.7811e-02,
                    8.2707e-02, 6.9900e-02, 5.9197e-02, 5.0157e-02,
                    4.2535e-02, 3.6140e-02, 3.0782e-02, 2.6233e-02,
                    2.2363e-02},
                   {1.9635e-02, 2.3167e-02, 2.7343e-02, 3.2293e-02,
                    3.8224e-02, 4.5335e-02, 5.3849e-02, 6.4039e-02,
                    7.6210e-02, 9.0739e-02, 1.0805e-01, 1.2841e-01,
                    1.5216e-01, 1.7960e-01, 2.1099e-01, 2.4671e-01,
                    2.8647e-01, 3.2996e-01, 3.7660e-01, 4.2547e-01,
                    4.7536e-01, 5.2473e-01, 5.7215e-01, 6.1632e-01,
                    6.5616e-01, 6.9075e-01, 7.1939e-01, 7.4199e-01,
                    7.5838e-01, 7.6848e-01, 7.7227e-01, 7.6848e-01,
                    7.5838e-01, 7.4199e-01, 7.1939e-01, 6.9075e-01,
                    6.5616e-01, 6.1632e-01, 5.7215e-01, 5.2473e-01,
                    4.7536e-01, 4.2547e-01, 3.7660e-01, 3.2996e-01,
                    2.8647e-01, 2.4671e-01, 2.1099e-01, 1.7960e-01,
                    1.5216e-01, 1.2841e-01, 1.0805e-01, 9.0739e-02,
                    7.6210e-02, 6.4039e-02, 5.3849e-02, 4.5335e-02,
                    3.8224e-02, 3.2293e-02, 2.7343e-02, 2.3167e-02,
                    1.9635e-02},
                   {1.7224e-02, 2.0450e-02, 2.4286e-02, 2.8860e-02,
                    3.4357e-02, 4.0979e-02, 4.8966e-02, 5.8612e-02,
                    7.0253e-02, 8.4257e-02, 1.0102e-01, 1.2094e-01,
                    1.4442e-01, 1.7196e-01, 2.0381e-01, 2.4013e-01,
                    2.8093e-01, 3.2594e-01, 3.7450e-01, 4.2563e-01,
                    4.7796e-01, 5.2991e-01, 5.7974e-01, 6.2599e-01,
                    6.6750e-01, 7.0344e-01, 7.3329e-01, 7.5676e-01,
                    7.7371e-01, 7.8410e-01, 7.8793e-01, 7.8410e-01,
                    7.7371e-01, 7.5676e-01, 7.3329e-01, 7.0344e-01,
                    6.6750e-01, 6.2599e-01, 5.7974e-01, 5.2991e-01,
                    4.7796e-01, 4.2563e-01, 3.7450e-01, 3.2594e-01,
                    2.8093e-01, 2.4013e-01, 2.0381e-01, 1.7196e-01,
                    1.4442e-01, 1.2094e-01, 1.0102e-01, 8.4257e-02,
                    7.0253e-02, 5.8612e-02, 4.8966e-02, 4.0979e-02,
                    3.4357e-02, 2.8860e-02, 2.4286e-02, 2.0450e-02,
                    1.7224e-02},
                   {1.5096e-02, 1.8041e-02, 2.1566e-02, 2.5793e-02,
                    3.0886e-02, 3.7044e-02, 4.4515e-02, 5.3604e-02,
                    6.4668e-02, 7.8109e-02, 9.4364e-02, 1.1389e-01,
                    1.3716e-01, 1.6461e-01, 1.9663e-01, 2.3350e-01,
                    2.7527e-01, 3.2170e-01, 3.7214e-01, 4.2549e-01,
                    4.8024e-01, 5.3460e-01, 5.8677e-01, 6.3512e-01,
                    6.7838e-01, 7.1569e-01, 7.4655e-01, 7.7071e-01,
                    7.8810e-01, 7.9871e-01, 8.0255e-01, 7.9871e-01,
                    7.8810e-01, 7.7071e-01, 7.4655e-01, 7.1569e-01,
                    6.7838e-01, 6.3512e-01, 5.8677e-01, 5.3460e-01,
                    4.8024e-01, 4.2549e-01, 3.7214e-01, 3.2170e-01,
                    2.7527e-01, 2.3350e-01, 1.9663e-01, 1.6461e-01,
                    1.3716e-01, 1.1389e-01, 9.4364e-02, 7.8109e-02,
                    6.4668e-02, 5.3604e-02, 4.4515e-02, 3.7044e-02,
                    3.0886e-02, 2.5793e-02, 2.1566e-02, 1.8041e-02,
                    1.5096e-02}};

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
