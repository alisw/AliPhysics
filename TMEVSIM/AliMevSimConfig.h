#ifndef ALIMEVSIMCONFIG_H
#define ALIMEVSIMCONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

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
//Comments taken from MevSim fortran source file (by Lanny Ray):
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

#include "TObject.h"

class AliGenMevSim;

class AliMevSimConfig : public TObject {

 protected:

  static const Int_t fgkMAX_MODEL = 4; //Maximum number of available models
  static const Int_t fgkMAX_CTRL = 4;//Maximum number of available controls

  Int_t fModelType;  //current type of model

  Int_t fReacPlaneCntrl; //reaction plane simulation model
  Float_t fPsiRMean; //fPsiRMean mean psi
  Float_t fPsiRStDev; //fPsiRStDev  psi variance

  Float_t fMultFacMean;//fMultFacMean Mean multiplicity
  Float_t fMultFacStDev;//fMultFacStDev multiplicity variance

  Float_t fNStDevMult;//see (10) n_stdev_mult
  Float_t fNStDevTemp;//(11) n_stdev_temp
  Float_t fNStDevSigma;//see (12) n_stdev_sigma
  Float_t fNStDevExpVel;//see (13) n_stdev_expvel
  Float_t fNStdDevPSIr;//see (14) n_stdev_PSIr
  Float_t fNStDevVn;//see (15) n_stdev_Vn
  Float_t fNStDevMultFac;//see (16) n_stdev_MultFac

  Int_t fNIntegPts;//see (17) n_integ_pts
  Int_t fNScanPts;//see (18) n_scan_pts

  void Init();

 public:

  AliMevSimConfig();
  AliMevSimConfig(Int_t modelType);

  ~AliMevSimConfig();

  void SetModelType(Int_t modelType);
  Int_t  GetModelType() const {return fModelType;}

  void SetRectPlane(Int_t ctrl, Float_t psiRMean = 0, Float_t psiRStDev = 0);
  void GetRectPlane(Int_t& ctrl, Float_t& psiRMean, Float_t& psiRStDev ) const
   {ctrl  = fReacPlaneCntrl; psiRMean = fPsiRMean; psiRStDev = fPsiRStDev;}
  
  void SetMultFac(Float_t mean, Float_t stDev);
  void GetMultFac(Float_t& mean, Float_t& stDev) const {mean = fMultFacMean ;stDev = fMultFacStDev;}

  void SetStDev(Float_t mult, Float_t temp, Float_t sigma,
                Float_t expVel, Float_t psiR, Float_t Vn, Float_t multFac);
  void GetStDev(Float_t& mult, Float_t& temp, Float_t& sigma,
                Float_t& expVel, Float_t& psiR, Float_t& Vn, Float_t& multFac) const;
  void SetGrid(Int_t integr, Int_t scan);
  void GetGrid(Int_t& integr, Int_t& scan) const {scan=fNScanPts;integr=fNIntegPts;}

  ClassDef(AliMevSimConfig,1)

};


#endif
