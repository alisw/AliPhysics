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
//__________________________________________________________________
////////////////////////////////////////////////////////////////////
//
// class AliMevSimConfig
//
// Class containing configuation inforamtion for MeVSim generator
// --------------------------------------------
// --------------------------------------------
// --------------------------------------------
// author: radomski@if.pw.edu.pl
//
//    (3) model_type - equals 1,2,3,4,5 or 6 so far.  See comments in
//                     Function dNdpty to see what is calculated.
//                     The models included are:
//                   = 1, Factorized mt exponential, Gaussian rapidity model
//                   = 2, Pratt non-expanding, spherical thermal source model
//                   = 3, Bertsch non-expanding spherical thermal source model
//                   = 4, Pratt spherically expanding, thermally equilibrated
//                        source model.
//                   = 5, Factorized pt and eta distributions input bin-by-bin.
//                   = 6, Fully 2D pt,eta distributions input bin-by-bin.
//                        NOTE: model_type = 1-4 are functions of (pt,y)
//                              model_type = 5,6 are functions of (pt,eta)
//    (4) reac_plane_cntrl - Can be either 1,2,3 or 4 where:
//                         = 1 to ignore reaction plane and anisotropic flow,
//                             all distributions will be azimuthally symm.
//                         = 2 to use a fixed reaction plane angle for all
//                             events in the run.
//                         = 3 to assume a randomly varying reaction plane
//                             angle for each event as determined by a
//                             Gaussian distribution.
//                         = 4 to assume a randomly varying reaction plane
//                             for each event in the run as determined by
//                             a uniform distribution from 0 to 360 deg.
//    (5) PSIr_mean, PSIr_stdev - Reaction plane angle mean and Gaussian
//                                std.dev. (both are in degrees) for cases
//                                with reac_plane_cntrl = 2 (use mean value)
//                                and 3.  Note: these are read in regardless
//                                of the value of reac_plane_cntrl.
//    (6) MultFac_mean, MultFac_stdev - Overall multiplicity scaling factor
//                                      for all PID types; mean and std.dev.;
//                                      for trigger fluctuations event-to-evt.
//    (7) pt_cut_min,pt_cut_max - Range of transverse momentum in GeV/c.
//    (8) eta_cut_min,eta_cut_max - Pseudorapidity range
//    (9) phi_cut_min,phi_cut_max - Azimuthal angular range in degrees.
//   (10) n_stdev_mult - Number of standard deviations about the mean value
//                       of multiplicity to include in the random event-to-
//                       event selection process.  The maximum number of
//                       steps that can be covered is determined by
//                       parameter n_mult_max_steps in the accompanying
//                       include file 'Parameter_values.inc' which is
//                       presently set at 1000, but the true upper limit for
//                       this is n_mult_max_steps - 1 = 999.
//   (11) n_stdev_temp - Same, except for the "Temperature" parameter.
//   (12) n_stdev_sigma- Same, except for the rapidity width parameter.
//   (13) n_stdev_expvel - Same, except for the expansion velocity parameter.
//   (14) n_stdev_PSIr   - Same, except for the reaction plane angle
//   (15) n_stdev_Vn     - Same, except for the anisotropy coefficients, Vn.
//   (16) n_stdev_MultFac - Same, except for the multiplicity scaling factor.
//   (17) n_integ_pts - Number of mesh points to use in the random model
//                      parameter selection process.  The upper limit is
//                      set by parameter nmax_integ in the accompanying
//                      include file 'Parameter_values.inc' which is presently
//                      set at 100, but the true upper limit for n_integ_pts
//                      is nmax_integ - 1 = 99.
//   (18) n_scan_pts  - Number of mesh points to use to scan the (pt,y)
//                      dependence of the model distributions looking for
//                      the maximum value.  The 2-D grid has
//                      n_scan_pts * n_scan_pts points; no limit to size of
//                      n_scan_pts.
//
////////////////////////////////////////////////////////////////////

#include "AliMevSimConfig.h"

ClassImp(AliMevSimConfig)


//////////////////////////////////////////////////////////////////////////////////////////////////

AliMevSimConfig::AliMevSimConfig() {
//def ctor
  Init();
}

//////////////////////////////////////////////////////////////////////////////////////////////////

AliMevSimConfig::AliMevSimConfig(Int_t modelType) {
//ctor
  Init();
  SetModelType(modelType);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

AliMevSimConfig::~AliMevSimConfig() {
//dtor
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::Init() {
  // Default Values

  fModelType = 1;
  fReacPlaneCntrl = 4;
  fPsiRMean = fPsiRStDev = 0;

  fMultFacMean  = 1.0;
  fMultFacStDev = 0.0;

  fNStDevMult = fNStDevTemp = fNStDevSigma = 3.0;
  fNStDevExpVel = fNStdDevPSIr = fNStDevVn = fNStDevMultFac = 3.0;
  
  fNIntegPts = fNScanPts = 100;

}

//////////////////////////////////////////////////////////////////////////////////////////////////
 
void AliMevSimConfig::SetModelType(Int_t modelType) {
//Sets type of the model
  if (modelType < 0 || modelType > fgkMAX_MODEL)
    Error("SetModelType","Wrog Model Type indentifier (%d)",modelType);

  fModelType = modelType;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::SetRectPlane(Int_t ctrl, Float_t psiRMean, Float_t psiRStDev) {
//Sets reaction plane parameters
  if (ctrl < 0 || ctrl > fgkMAX_CTRL)
    Error("SetReactPlane","Wrong Control Parameter (%d)", ctrl);

  fReacPlaneCntrl = ctrl;
  fPsiRMean = psiRMean;
  fPsiRStDev = psiRStDev;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::SetMultFac(Float_t mean, Float_t stDev) {
  //Sets multiplicity mean and variance
  fMultFacMean = mean;
  fMultFacStDev = stDev;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::SetStDev(Float_t mult, Float_t temp, Float_t sigma,
  Float_t expVel, Float_t psiR, Float_t Vn, Float_t multFac) {
//sets Dev parameters (whatever Dev is)
  fNStDevMult = mult;
  fNStDevTemp = temp;
  fNStDevSigma = sigma;
  fNStDevExpVel = expVel;
  fNStdDevPSIr = psiR;
  fNStDevVn = Vn;
  fNStDevMultFac =multFac;

}
void AliMevSimConfig::GetStDev(Float_t& mult, Float_t& temp, Float_t& sigma,
                Float_t& expVel, Float_t& psiR, Float_t& Vn, Float_t& multFac) const
{
 //returns dev parameters
   mult  = fNStDevMult;
   temp  = fNStDevTemp;
   sigma  = fNStDevSigma;
   expVel  = fNStDevExpVel;
   psiR  = fNStdDevPSIr;
   Vn  = fNStDevVn;
   multFac  = fNStDevMultFac;
 
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::SetGrid(Int_t integr, Int_t scan) {
//Sets grid 
  fNIntegPts = integr;
  fNScanPts = scan;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
 
