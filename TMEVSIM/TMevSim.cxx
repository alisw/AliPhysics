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

////////////////////////////////////////////////////////////////////////////
// 
// TMevSim
//
// TMevSim is an interface class between the event generator MEVSIM and 
// the ROOT system. The current implementation is based on the 6.11.2000
// version provided by Lanny Ray on /afs/cern.ch/user/y/yiota/ray/mult_gen.
// 
// Authors of MEVSIM:
// For The STAR Collaboration
//
//  Lanny Ray      
//     Dept. of Physics
//     The University of Texas at Austin
//     Austin, Texas 78712
//     (512) 471-6107
//     ray@physics.utexas.edu
//
//  Ron Longacre   email:
//
//
////////////////////////////////////////////////////////////////////////////
//
// I.    OVERVIEW
//
//    This code is intended to provide a quick means of producing
//    uncorrelated simulated events for event-by-event studies,
//    detector acceptance and efficiency studies, etc.  The
//    user selects the number of events, the one-particle distribution
//    model, the Geant particles to include, the ranges in transverse
//    momentum, pseudorapidity and azimuthal angle, the mean
//    multiplicity for each particle type for the event run, the
//    mean temperature, Rapidity width, etc., and the standard deviations
//    for the event-to-event variation in the model parameters.
//    Note that these events are produced in the c.m. frame only.
//    
//    Anisotropic flow may also be simulated by introducing explicit
//    phi-dependence (azimuthal angle) in the particle distributions.  
//    The assumed model is taken from Poskanzer and Voloshin, Phys. Rev.
//    C58, 1671 (1998), Eq.(1), where we use,
//
//         E d^3N/dp^3 = (1/2*pi*pt)*[d^2N/dpt*dy]
//            * [1 + SUM(n=1,nflowterms){2*Vn*cos[n(phi-PSIr)]}]
//
//    with up to 'nflowterms' (currently set to 6, see file
//    Parameter_values.inc) Fourier components allowed.  Vn are
//    coefficients and PSIr is the reaction plane angle.
//    The algebraic signs of the Vn terms for n=odd are reversed
//    from their input values for particles with rapidity (y) < 0
//    as suggested in Poskanzer and Voloshin.
//    The flow parameters can depend on pt and rapidity (y) according
//    to the model suggested by Art Poskanzer (Feb. 2000) and as
//    defined in the Function Vn_pt_y.
//
//    The user may select either to have the same multiplicity per
//    particle type for each event or to let the multiplicity vary
//    randomly according to a Poisson distribution. In addition, an
//    overall multiplicative scale factor can be applied to each
//    particle ID's multiplicity (same factor applied to each PID).
//    This scaling can vary randomly according to a Gaussian from
//    event-to-event.  This is to simulate trigger acceptance
//    fluctuations.  Similarly the
//    parameters of the one-particle distribution models may either
//    be fixed to the same value for each event or allowed to randomly
//    vary about a specified mean with a specified standard deviation
//    and assuming a Gaussian distribution.
//
//    With respect to the reaction plane and anisotropic flow simulation,
//    the user may select among four options:
//       (1) ignore reaction plane and anisotropic flow effects; all
//           distributions will be azimuthally invariant, on average.
//       (2) assume a fixed reaction plane angle, PSIr, for all events
//           in the run.
//       (3) assume a Gaussian distribution with specified mean and
//           standard deviation for the reaction plane angles for the
//           events in the run.  PSIr is randomly determined for each
//           event.
//       (4) assume uniformly distributed, random values for the reaction  
//           plane angles from 0 to 360 deg., for each event in the run.
//
//    The user may also select the anisotropic flow parameters, Vn,
//    to either be fixed for each event, or to randomly vary from event
//    to event according to a Gaussian distribution where the user must
//    specify the mean and std. dev.  For both cases the input file must
//    list the 'nflowterms' (e.g. 6) values of the mean and Std.dev. for
//    the Vn parameters for all particle ID types included in the run.
//
//    The available list of particles has been increased to permit a
//    number of meson and baryon resonances.  For those with broad widths
//    the code samples the mass distribution for the resonance and outputs
//    the resonance mass for each instance in a special kinematic file
//    (see file unit=9, filename = 'mult_gen.kin').  The resonance shapes
//    are approximately Breit-Wigner and are specific for each resonance 
//    case.  The additional particle/resonances include: rho(+,-,0),
//    omega(0), eta', phi, J/Psi, Delta(-,0,+,++) and K*(+,-,0).  Masses
//    are sampled for the rho, omega, phi, Deltas and D*s. 
//    Refer to SUBR: Particle_prop and Particle_mass for the explicit
//    parameters, resonance shape models, and sampling ranges.
//
//    The input is from a file, named 'mult_gen.in'.  The output is
//    loaded into a file named 'mult_gen.out' which includes run
//    header information, event header information and the EVENT: and
//    TRACK: formats as in the new STAR TEXT Format for event generator
//    input to GSTAR.  A log file, 'mult_gen.log' is also written which
//    may contain error messages.  Normally this file should be empty
//    after a successful run.  These filenames can easily be changed
//    to more suitable names by the script that runs the program or
//    by hand.
//
//
// II.   ALGORITHM
//
//
//
//    The method for generating random multiplicities and model parameter
//    values involves the following steps:
//       (1) The Poisson or Gaussian distributions are computed and
//           loaded into function f().
//       (2) The distribution f(x') is integrated from xmin to x
//           and saved from x = xmin to x = xmax.  The range and mesh
//           spaces are specified by the user.
//       (3) The integral of f is normalized to unity where 
//           integral[f(x')](at x = xmin) = 0.0
//           integral[f(x')](at x = xmax) = 1.0
//       (4) A random number generator is called which delivers values
//           between 0.0 and 1.0.  
//       (5) We consider the coordinate x (from xmin to xmax) to be
//           dependent on the integral[f].  Using the random number
//           for the selected value of integral[f] the value of x
//           is obtained by interpolation.
//
//    An interpolation subroutine from Rubin Landau, Oregon State Univ.,
//    is used to do this interpolation; it involves uneven mesh point 
//    spacing.
//
//    The method for generating the particle momenta uses the
//    standard random elimination method and involves the following
//    steps:
//
//    For model_type = 1,2,3,4 which are functions of pt,y (see following):
//       (1) The y range is computed using the pseudorapidity (eta)
//           range and includes ample cushioning around the sides
//           along the eta acceptance edges.
//       (2) The transverse momentum (pt) and rapidity (y) are
//           randomly chosen within the specified ranges.
//       (3) The pseudorapidity is computed for this (pt,y) value
//           (and the mass for each pid) and checked against the
//           pseudorapidity acceptance range.
//       (4) If the pseudorapidity is within range then the one-particle
//           model distribution is calculated at this point and its ratio
//           to the maximum value throughout (pt,eta) acceptance region
//           is calculated.
//       (5) Another random number is called and if less than the ratio
//           from step#4 the particle momentum is used; if not, then 
//           another trial value of (pt,y) is obtained.
//       (6) This continues until the required multiplicity for the
//           specific event and particle type has been satisfied.
//       (7) This process is repeated for the requested number of particle
//           types and events.
//
//    For model_type = 5,6 (see following) which are input bin-by-bin
//    in pt,eta:
//       (1) The transverse momentum (pt) and pseudorapidity (eta) are 
//           randomly chosen within the specified ranges.
//       (2) The one-particle model distribution is calculated at this
//           point and its ratio to the maximum value throughout the
//           (pt,eta) region is calculated.
//       (3) Another random number is called and if less than the ratio
//           from step(2) the particle momentum is used; if not then
//           another trial value of (pt,eta) is obtained.
//       (4) This continues until the required multiplicity for the 
//           specific event and particle type has been satisfied.
//       (5) This process is repeated for the requested number of particle
//           types and events. 
//
//    Problematic parameter values are tested, bad input values are checked
//    and in some cases may be changed so that the program will not crash.
//    In some cases the code execution is stopped.
//    Some distributions and/or unusual model parameter values may cause the
//    code to hang up due to the poor performance of the "elimination"
//    method for very strongly peaked distributions.  These are tested for
//    certain problematic values and if necessary these events are aborted.
//    A message, "*** Event No.    2903 ABORTED:" for example is printed
//    in the 'mult_gen.out' file.  Temperatures .le. 0.01 GeV and rapidity
//    width parameters .le. 0.01 will cause the event to abort.
//
//
//
// III.  DESCRIPTION OF THE INPUT:
//
//
//    The input is described below in the 'read' statements and also in
//    the sample input file.  Some additional comments are as follows:
//
//    (1) n_events - Selected number of events in run. Can be anything
//                   .ge. 1.
//    (2) n_pid_type - Number of particle ID types to include in the
//                     particle list. e.g. pi(+) and pi(-) are counted
//                     separately.  The limit is set by parameter npid
//                     in the accompanying include file 'Parameter_values.inc'
//                     and is presently set at 20.
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
//   (19) irand       - Starting random number seed.
//
//**************************************************************************
//   FOR MODEL_TYPE = 1,2,3 or 4:
//   Input the following 7 lines for each particle type; repeat these
//   set of lines n_pid_type times:
//
//        (a) gpid - Geant Particle ID code number
//        (b) mult_mean,mult_variance_control - Mean multiplicity and
//                                              variance control where:
//            mult_variance_control = 0 for no variance in multiplicity 
//            mult_variance_control = 1 to allow Poisson distribution for
//                                      particle multiplicities for all events.
//            Note that a hard limit exists for the maximum possible
//            multiplicity for a given particle type per event.  This is
//            determined by parameter factorial_max in accompanying include
//            file 'common_facfac.inc' and is presently set at 10000.
//        (c) Temp_mean, Temp_stdev - Temperature parameter mean (in GeV)
//            and standard deviation (Gaussian distribution assumed).
//        (d) sigma_mean, sigma_stdev - Rapidity distribution width (sigma)
//            parameter mean and standard deviation (Gaussian distribution
//            assumed).
//        (e) expvel_mean, expvel_stdev - S. Pratt expansion velocity
//            (in units of c) mean and standard deviation (Gaussian 
//            distribution assumed).
//        (f) Vn_mean(k);  k=1,4  - Anisotropic flow parameters, mean values
//                                  for Fourier component n=1.
//        (g) Vn_stdev(k); k=1,4  - Anisotropic flow parameters, std.dev.
//                                  values for Fourier component n=1.
//
//            Repeat the last two lines of input for remaining Fourier
//            components n=2,3...6.  Include all 6 sets of parameters
//            even if these are not used by the model for Vn(pt,y) (set
//            unused parameter means and std.dev. to 0.0).  List 4 values
//            on every line, even though for n=even the 4th quantity is
//            not used.
//
//**************************************************************************
//   FOR MODEL_TYPE = 5 input the following set of lines for each particle
//                      type; repeat these n_pid_type times.
//
//        (a) gpid - Geant Particle ID code number
//        (b) mult_mean,mult_variance_control - Mean multiplicity and
//                                              variance control where:
//            mult_variance_control = 0 for no variance in multiplicity
//            mult_variance_control = 1 to allow Poisson distribution for
//                                      particle multiplicities for all events.
//        (c) pt_start, eta_start - minimum starting values for pt, eta 
//                                  input for the bin-by-bin distributions.
//        (d) n_pt_bins, n_eta_bins - # input pt and eta bins.
//        (e) delta_pt, pt_bin - pt bin size and function value, repeat for
//                               each pt bin.
//        (f) delta_eta, eta_bin - eta bin size and function value, repeat
//                                 for each eta bin.
//        (g) Vn_mean(k);  k=1,4  - Anisotropic flow parameters, mean values
//                                  for Fourier component n=1.
//        (h) Vn_stdev(k); k=1,4  - Anisotropic flow parameters, std.dev.
//                                  values for Fourier component n=1.
//
//            Repeat the last two lines of input for remaining Fourier
//            components n=2,3...6.  Include all 6 sets of parameters
//            even if these are not used by the model for Vn(pt,y) (set
//            unused parameter means and std.dev. to 0.0).  List 4 values
//            on every line, even though for n=even the 4th quantity is
//            not used.
//
//        NOTE: The pt, eta ranges must fully include the requested ranges
//              in input #4 and 5 above; else the code execution will stop.
//
//        Also, variable bin sizes are permitted for the input distributions.
//
//        Also, this input distribution is used for all events in the run;
//        no fluctuations in this "parent" distribution are allowed from 
//        event-to-event.
//
//**************************************************************************
//   FOR MODEL_TYPE = 6 input the following set of lines for each particle
//                      type; repeat these n_pid_type times.
//
//        (a) gpid - Geant Particle ID code number
//        (b) mult_mean,mult_variance_control - Mean multiplicity and
//                                              variance control where:
//            mult_variance_control = 0 for no variance in multiplicity
//            mult_variance_control = 1 to allow Poisson distribution for
//                                      particle multiplicities for all events.
//        (c) pt_start, eta_start - minimum starting values for pt, eta
//                                  input for the bin-by-bin distributions.
//        (d) n_pt_bins, n_eta_bins - # input pt and eta bins.
//        (e) delta_pt - pt bin size, repeat for each pt bin. 
//        (f) delta_eta - eta bin size, repeat for each eta bin.
//        (g) i,j,pt_eta_bin(i,j) - read pt (index = i) and eta (index = j)
//                                  bin numbers and bin value for full 2D space
//        (h) Vn_mean(k);  k=1,4  - Anisotropic flow parameters, mean values
//                                  for Fourier component n=1.
//        (i) Vn_stdev(k); k=1,4  - Anisotropic flow parameters, std.dev.
//                                  values for Fourier component n=1.
//
//            Repeat the last two lines of input for remaining Fourier
//            components n=2,3...6.  Include all 6 sets of parameters
//            even if these are not used by the model for Vn(pt,y) (set
//            unused parameter means and std.dev. to 0.0).  List 4 values
//            on every line, even though for n=even the 4th quantity is
//            not used.
//
//        NOTE: The pt, eta ranges must fully include the requested ranges
//              in input #4 and 5 above; else the code execution will stop.
//
//        Also, variable bin sizes are permitted for the input distributions.
//
//        Also, this input distribution is used for all events in the run;
//        no fluctuations in this "parent" distribution are allowed from
//        event-to-event.
//
///////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <TParticle.h>
#include <TClonesArray.h>

#include "TMevSim.h"
#include "TMevSimPartTypeParams.h"

#ifndef WIN32
# define multgen multgen_
# define type_of_call
#else
# define multgen MULTGEN
# define type_of_call _stdcall
#endif


ClassImp(TMevSim)


extern "C" void type_of_call multgen();

//______________________________________________________________________________
TMevSim::TMevSim(Int_t nEvents, Int_t modelType, Int_t reacPlaneCntrl,
	   Float_t psiRMean, Float_t psiRStDev, Float_t multFacMean, Float_t multFacStDev,
	   Float_t ptCutMin, Float_t ptCutMax, Float_t etaCutMin, Float_t etaCutMax, 
	   Float_t phiCutMin, Float_t phiCutMax, Int_t irand) : TGenerator("MevSim", "MevSim")
{
// TMevSim constructor: initializes all the event-wide variables of MevSim with
// user supplied values, or with the default ones (declared in the header file).
// It also allocates space for the array which will store parameters specific to
// each particle species.
// Caution: Setting nEvents > 1 will have no effect, since only the last generated 
// event will be stored in POUT COMMON, and therefore only one event can be 
// accessible at a time. 
 
   fNEvents = nEvents;
   fModelType = modelType;
   fReacPlaneCntrl = reacPlaneCntrl;
   fPsiRMean = psiRMean;
   fPsiRStDev = psiRStDev;
   fMultFacMean = multFacMean;
   fMultFacStDev = multFacStDev;
   fPtCutMin = ptCutMin;
   fPtCutMax = ptCutMax;
   fEtaCutMin = etaCutMin;
   fEtaCutMax = etaCutMax;
   fPhiCutMin = phiCutMin;
   fPhiCutMax = phiCutMax;
   fNStDevMult = fNStDevTemp = fNStDevSigma = fNStDevExpVel = fNStdDevPSIr = fNStDevVn = fNStDevMultFac = 3.0;
   fNIntegPts = 100;
   fNScanPts = 100;
   firand = irand;
   fParticleTypeParameters = new TClonesArray("TMevSimPartTypeParams",10);
   fNPDGCodes = 0;
   DefineParticles();
}
//______________________________________________________________________________
TMevSim::~TMevSim()
{
// TMevSim destructor: destroys the object and all the particle information stored
// in the list.
   
  if (fParticleTypeParameters) {
    fParticleTypeParameters->Clear();
    delete fParticleTypeParameters;
    fParticleTypeParameters = 0;
  }
}
//______________________________________________________________________________
TMevSim::TMevSim(TMevSim& mevsim) : TGenerator(mevsim) {
// The copy constructor
   
   *this = mevsim;
}
//______________________________________________________________________________

TMevSim& TMevSim::operator=(const TMevSim& mevsim) {
// An assignment operator: initializes all the event-wide variables of MevSim with
// the ones from a copied object. It also copies the parameters specific to
// each particle species.
   
   fNEvents = mevsim.GetNEvents();
   fModelType = mevsim.GetModelType();
   fReacPlaneCntrl = mevsim.GetReacPlaneCntrl();
   fPsiRMean = mevsim.GetPsiRMean();
   fPsiRStDev = mevsim.GetPsiRStDev();
   fMultFacMean = mevsim.GetMultFacMean();
   fMultFacStDev = mevsim.GetMultFacStDev();
   fPtCutMin = mevsim.GetPtCutMin();
   fPtCutMax = mevsim.GetPtCutMax();
   fEtaCutMin = mevsim.GetEtaCutMin();
   fEtaCutMax = mevsim.GetEtaCutMax();
   fPhiCutMin = mevsim.GetPhiCutMin();
   fPhiCutMax = mevsim.GetPhiCutMax();
   fNStDevMult = mevsim.GetNStDevMult();
   fNStDevTemp = mevsim.GetNStDevTemp();
   fNStDevSigma =GetNStDevSigma();
   fNStDevExpVel = mevsim.GetNStDevExpVel();
   fNStdDevPSIr = mevsim.GetNStDevPSIr(); 
   fNStDevVn =  mevsim.GetNStDevVn();
   fNStDevMultFac =  mevsim.GetNStDevMultFac();
   fNIntegPts = mevsim.GetNintegPts();
   fNScanPts = mevsim.GetNScanPts();
   firand = mevsim.firand;
   fParticleTypeParameters = new TClonesArray("TMevSimPartTypeParams",mevsim.GetNPidTypes());
   for (int i=0; i< mevsim.GetNPidTypes(); i++) 
     {	
	TMevSimPartTypeParams *temp = 0;
	mevsim.GetPartTypeParamsByIndex(i,temp);
	fParticleTypeParameters->AddLast(temp);
     }
   DefineParticles();
   return (*this);
}
//______________________________________________________________________________
void        TMevSim::Initialize() {
// TMevSim initialization: creates an input file for the FORTRAN
// program MevSim. Converts all the event-wide information and particle
// specific information to the format readable by MevSim and writes it
// to disk in current directory.
// Caution: At least one TMevSimPartTypeParams object must be created and 
// added to the collection before event generation can start.
   
  TMevSimPartTypeParams * params = 0;
  

   ofstream *file = new ofstream("mult_gen.in",ios::out | ios::trunc);
   // Write out the parameters to the pramameter file
   *file << "  " << fNEvents << "  ! Number of Events \n";
   *file << "  " << GetNPidTypes() << " \n";
   *file << "  " << fModelType << " \n";
   *file << "  " << fReacPlaneCntrl << " \n";
   file->setf(ios::showpoint);
   *file << "  " << fPsiRMean << "   " << fPsiRStDev << " \n";
   *file << "  " << fMultFacMean << "   " << fMultFacStDev << " \n";
   *file << "  " << fPtCutMin << "   " << fPtCutMax << " \n";
   *file << "  " << fEtaCutMin << "   " << fEtaCutMax << " \n";
   *file << "  " << fPhiCutMin << "   " << fPhiCutMax << " \n";
   *file << "  " << fNStDevMult << " \n";
   *file << "  " << fNStDevTemp << " \n";
   *file << "  " << fNStDevSigma << " \n";
   *file << "  " << fNStDevExpVel << " \n";
   *file << "  " << fNStdDevPSIr << " \n";
   *file << "  " << fNStDevVn << " \n";
   *file << "  " << fNStDevMultFac << " \n";
   *file << "  " << fNIntegPts << " \n";
   *file << "  " << fNScanPts << " \n";
   *file << "  " << firand << " \n";
   // Write out particle specific information
   for (Int_t i=0; i< (fParticleTypeParameters->GetLast() + 1); i++) {

     params = (TMevSimPartTypeParams *) ((*fParticleTypeParameters)[i]);
     
     *file << "  " << params->GetGPid() << "       ! Particle GEANT Pid \n";
     *file << "  " << params->GetMultMean() << "  " << params->GetMultVarianceControl() << "  \n";
     *file << "  " << params->GetTempMean() << "  " << params->GetTempStDev() << "  \n";
     *file << "  " << params->GetSigmaMean() << "  " << params->GetSigmaStDev() << "  \n";
     *file << "  " << params->GetExpVelMean() << "  " << params->GetExpVelStDev() << "  \n";

     for (Int_t cnt1 = 0; cnt1 < NFLOWTERMS; cnt1++) {
       *file << "  ";
       Int_t cnt2;
       for (cnt2 = 0; cnt2 < 4; cnt2++) *file << params->GetVnMeanComponent(cnt1, cnt2) << "  ";
       *file << "  \n  ";
       for (cnt2 = 0; cnt2 < 4; cnt2++) *file << params->GetVnStDevComponent(cnt1, cnt2) << "  ";
       *file << "  \n";
     }
   }
   file->close();

}
//______________________________________________________________________________
void TMevSim::GenerateEvent() {
// Generates one MevSim event. TMevSim::Initialize() must be called prior
// to calling this function.
   
   Info("GenerateEvent","Calling FORTRAN multgen()");
   multgen();
}

//______________________________________________________________________________
Int_t TMevSim::ImportParticles(TClonesArray *particles, Option_t */*option*/)
{
// Read in particles created by MevSim into the TClonesArray(). The Initialize()
// and GenrateEvent() functions must be called prior to calling this funtion.
// The particles are read from the COMMON POUT. Right now the only provided 
// information is Geant PID, 3 momentum components and the energy of the particle.
   
   if (particles == 0) return 0;
   TClonesArray &aParticles = *particles;
   aParticles.Clear();
   
   Int_t totpart = 0;
   for (Int_t nrpidtype=0; nrpidtype < (fParticleTypeParameters->GetLast() + 1); nrpidtype++) {
      Int_t nrpart = 0;
      Int_t pidcode = ((TMevSimPartTypeParams *) (*fParticleTypeParameters)[nrpidtype])->GetGPid();
      while ((TRACK.pout[(4*nrpart+3)*NPID+nrpidtype] > 0.0) || (TRACK.pout[(4*nrpart)*NPID+nrpidtype] != 0.0)) {
	 int poffset = 4*nrpart*NPID+nrpidtype;
	 Float_t px = TRACK.pout[poffset];
	 poffset += NPID;
	 Float_t py = TRACK.pout[poffset];
	 poffset += NPID;
	 Float_t pz = TRACK.pout[poffset];
	 poffset += NPID;
	 Float_t mass = TRACK.pout[poffset];
	 new(aParticles[totpart+nrpart]) TParticle(
					  PDGFromId(pidcode),  // Get the PDG ID from GEANT ID
					  0,
					  0,
					  0,
					  0,
					  0,
					  px,
					  py,
					  pz,
					  sqrt(mass*mass+px*px+py*py+pz*pz),                               
					  0,
					  0,
					  0,
					  0);
	 nrpart++;
      }
      totpart += nrpart;
   }
   return totpart;
}
//______________________________________________________________________________
TObjArray * TMevSim::ImportParticles(Option_t */*option*/)
{
// Read in particles created by MevSim into the TClonesArray(). The Initialize()
// and GenrateEvent() functions must be called prior to calling this funtion.
// The particles are read from the COMMON POUT. Right now the only provided 
// information is Geant PID, 3 momentum components and the energy of the particle.
   
   fParticles->Clear();
   
   for (Int_t nrpidtype=0; nrpidtype < (fParticleTypeParameters->GetLast() + 1); nrpidtype++) {
      Int_t nrpart = 0;
      Int_t pidcode = ((TMevSimPartTypeParams *) (*fParticleTypeParameters)[nrpidtype])->GetGPid();
      while ((TRACK.pout[(4*nrpart+3)*NPID+nrpidtype] > 0.0) || (TRACK.pout[(4*nrpart)*NPID+nrpidtype] != 0.0)) {
	 int poffset = 4*nrpart*NPID+nrpidtype;
	 Float_t px = TRACK.pout[poffset];
	 poffset += NPID;
	 Float_t py = TRACK.pout[poffset];
	 poffset += NPID;
	 Float_t pz = TRACK.pout[poffset];
	 poffset += NPID;
	 Float_t mass = TRACK.pout[poffset];
	 TParticle * p = new TParticle(
					  PDGFromId(pidcode),  // Get the PDG ID from GEANT ID
					  0,
					  0,
					  0,
					  0,
					  0,
					  px,
					  py,
					  pz,
					  sqrt(mass*mass+px*px+py*py+pz*pz),                               
					  0,
					  0,
					  0,
					  0);
	 fParticles->Add(p);
	 nrpart++;
      }
   }
   return fParticles;
}
//______________________________________________________________________________
void        TMevSim::SetNEvents(Int_t nEvents ) { 
// Sets the number of events to be generated by MevSim.
// Caution: Setting nEvents > 1 will have no effect, since only the last generated 
// event will be stored in POUT COMMON, and therefore only one event can be 
// accessible at a time. 

   fNEvents = nEvents;
}
//______________________________________________________________________________
Int_t       TMevSim::GetNEvents() const {
   return fNEvents;
}
//______________________________________________________________________________
Int_t       TMevSim::GetNPidTypes() const {
   return fParticleTypeParameters->GetLast()+1;
}
//______________________________________________________________________________
void        TMevSim::SetModelType(Int_t modelType) {
   fModelType  = modelType;
}
//______________________________________________________________________________
Int_t       TMevSim::GetModelType() const {
   return fModelType;
}
//______________________________________________________________________________
void        TMevSim::SetReacPlaneCntrl(Int_t reacPlaneCntrl) {
   fReacPlaneCntrl = reacPlaneCntrl; 
}
//______________________________________________________________________________
Int_t       TMevSim::GetReacPlaneCntrl() const {
   return fReacPlaneCntrl;
}
//______________________________________________________________________________
void        TMevSim::SetPsiRParams(Float_t psiRMean,  Float_t psiRStDev) {
   fPsiRMean = psiRMean;
   fPsiRStDev = psiRStDev;
}
//______________________________________________________________________________
Float_t     TMevSim::GetPsiRMean() const { 
   return fPsiRMean;
}
//______________________________________________________________________________
Float_t     TMevSim::GetPsiRStDev() const { 
   return fPsiRStDev;
}
//______________________________________________________________________________
void        TMevSim::SetMultFacParams(Float_t multFacMean,  Float_t multFacStDev) {
   fMultFacMean = multFacMean;
   fMultFacStDev = multFacStDev;
}
//______________________________________________________________________________
Float_t     TMevSim::GetMultFacMean() const { 
   return fMultFacMean;     
}
//______________________________________________________________________________
Float_t     TMevSim::GetMultFacStDev() const { 
   return fMultFacStDev;
}
//______________________________________________________________________________
void        TMevSim::SetPtCutRange(Float_t ptCutMin,  Float_t ptCutMax) { 
   fPtCutMin = ptCutMin;
   fPtCutMax = ptCutMax;
}
//______________________________________________________________________________
Float_t     TMevSim::GetPtCutMin() const { 
   return fPtCutMin;
}
//______________________________________________________________________________
Float_t     TMevSim::GetPtCutMax() const { 
   return fPtCutMax;
}
//______________________________________________________________________________
void        TMevSim::SetEtaCutRange(Float_t etaCutMin,  Float_t etaCutMax) {     fEtaCutMin = etaCutMin;
   fEtaCutMax = etaCutMax;
}

//______________________________________________________________________________
   Float_t     TMevSim::GetEtaCutMin() const { 
     return fEtaCutMin;
}
//______________________________________________________________________________
   Float_t     TMevSim::GetEtaCutMax() const {
     return fEtaCutMax;
}
//______________________________________________________________________________
void        TMevSim::SetPhiCutRange(Float_t phiCutMin,  Float_t phiCutMax) {
   fPhiCutMin = phiCutMin;
   fPhiCutMax = phiCutMax;
}
//______________________________________________________________________________
Float_t     TMevSim::GetPhiCutMin() const { 
   return fPhiCutMin;
}
//______________________________________________________________________________
Float_t     TMevSim::GetPhiCutMax() const { 
   return fPhiCutMax;
}
//______________________________________________________________________________
void        TMevSim::SetNStDevMult(Float_t nStDevMult) { 
   fNStDevMult = nStDevMult;
}
//______________________________________________________________________________
Float_t       TMevSim::GetNStDevMult() const { 
   return  fNStDevMult;
}
//______________________________________________________________________________
void        TMevSim::SetNStDevTemp(Float_t nStDevTemp) { 
   fNStDevTemp = nStDevTemp;
}
//______________________________________________________________________________
Float_t       TMevSim::GetNStDevTemp() const { 
   return fNStDevTemp;
}
//______________________________________________________________________________
void        TMevSim::SetNStDevSigma(Float_t nStDevSigma) { 
   fNStDevSigma = nStDevSigma;
}
//______________________________________________________________________________
Float_t       TMevSim::GetNStDevSigma() const { 
   return fNStDevSigma;  
}
//______________________________________________________________________________
void        TMevSim::SetNStDevExpVel(Float_t nStDevExpVel) { 
   fNStDevExpVel = nStDevExpVel;
}
//______________________________________________________________________________
Float_t       TMevSim::GetNStDevExpVel() const { 
   return fNStDevExpVel;
}
//______________________________________________________________________________
void        TMevSim::SetNStDevPSIr(Float_t nStDevPSIr) { 
   fNStdDevPSIr = nStDevPSIr;
}
//______________________________________________________________________________
Float_t       TMevSim::GetNStDevPSIr() const { 
   return fNStdDevPSIr;
}
//______________________________________________________________________________
void        TMevSim::SetNStDevVn(Float_t nStDevVn) { 
   fNStDevVn = nStDevVn;
}
//______________________________________________________________________________
Float_t       TMevSim::GetNStDevVn() const { 
   return fNStDevVn; 
}
//______________________________________________________________________________
void        TMevSim::SetNStDevMultFac(Float_t nStDevMultFac) { 
   fNStDevMultFac = nStDevMultFac;
}
//______________________________________________________________________________
Float_t       TMevSim::GetNStDevMultFac() const { 
   return fNStDevMultFac;  
}
//______________________________________________________________________________
void        TMevSim::SetNIntegPts(Int_t nIntegPts) { 
   fNIntegPts = nIntegPts;
}
//______________________________________________________________________________
Int_t       TMevSim::GetNintegPts() const { 
   return fNIntegPts;
}
//______________________________________________________________________________
void        TMevSim::SetNScanPts(Int_t nScanPts) { 
   fNScanPts = nScanPts;
}
//______________________________________________________________________________
Int_t       TMevSim::GetNScanPts() const { 
   return fNScanPts;  
}
//______________________________________________________________________________
void        TMevSim::AddPartTypeParams(TMevSimPartTypeParams *params) {
// Add the particle specied parameters and the end of the list.
   
  //cout << params << " " <<  fParticleTypeParameters <<  endl;

  //fParticleTypeParameters->Dump(); 
  params->Dump();
  
  Int_t last = fParticleTypeParameters->GetLast();
  new ((*fParticleTypeParameters)[last+1]) TMevSimPartTypeParams(*params);
}
//______________________________________________________________________________
void        TMevSim::SetPartTypeParams(Int_t index, TMevSimPartTypeParams *params) 
{
// Create the new copy particle species parameters provided by params, and store
// them in the position designated by index. 
   
   *((TMevSimPartTypeParams *) ((*fParticleTypeParameters)[index])) = *params;
}
//______________________________________________________________________________
void TMevSim::GetPartTypeParamsByIndex(Int_t index, TMevSimPartTypeParams *params) const
{
// Return the particle parameters stored in the list at the postion index.   
// Returns NULL if index is out of bounds.
   
   if ((index < fParticleTypeParameters->GetLast()) && (index >= 0))
	params = (TMevSimPartTypeParams *) (*fParticleTypeParameters)[index];
   else
     params = NULL;
}
//______________________________________________________________________________
void TMevSim::GetPartTypeParamsByGPid(Int_t gpid, TMevSimPartTypeParams *params) const
{
// Return the particle parameters for the particle with Geant PID gpid.
// Returns NULL if the parameters for such particle do not exist in the list.   
   
   Int_t i = -1;
   
   while (++i <= fParticleTypeParameters->GetLast())
     {
	if (((TMevSimPartTypeParams *) (*fParticleTypeParameters)[i])->GetGPid() == gpid)
	  {
	     params = (TMevSimPartTypeParams *) (*fParticleTypeParameters)[i];
	     return;
	  }
     }
   params = NULL;
   return;
}
//_____________________________________________________________________________
Int_t TMevSim::PDGFromId(Int_t id) const 
{
  //
  // Return PDG code and pseudo ENDF code from Geant3 code
  //
  if(id>0 && id<fNPDGCodes) return fPDGCode[id];
  else return -1;
}
//_____________________________________________________________________________
void TMevSim::DefineParticles() 
{
  //
  // Load standard numbers for GEANT particles and PDG conversion
  fPDGCode[fNPDGCodes++]=-99;   //  0 = unused location
  fPDGCode[fNPDGCodes++]=22;    //  1 = photon
  fPDGCode[fNPDGCodes++]=-11;   //  2 = positron
  fPDGCode[fNPDGCodes++]=11;    //  3 = electron
  fPDGCode[fNPDGCodes++]=12;    //  4 = neutrino e
  fPDGCode[fNPDGCodes++]=-13;   //  5 = muon +
  fPDGCode[fNPDGCodes++]=13;    //  6 = muon -
  fPDGCode[fNPDGCodes++]=111;   //  7 = pi0
  fPDGCode[fNPDGCodes++]=211;   //  8 = pi+
  fPDGCode[fNPDGCodes++]=-211;  //  9 = pi-
  fPDGCode[fNPDGCodes++]=130;   // 10 = Kaon Long
  fPDGCode[fNPDGCodes++]=321;   // 11 = Kaon +
  fPDGCode[fNPDGCodes++]=-321;  // 12 = Kaon -
  fPDGCode[fNPDGCodes++]=2112;  // 13 = Neutron
  fPDGCode[fNPDGCodes++]=2212;  // 14 = Proton
  fPDGCode[fNPDGCodes++]=-2212; // 15 = Anti Proton
  fPDGCode[fNPDGCodes++]=310;   // 16 = Kaon Short
  fPDGCode[fNPDGCodes++]=221;   // 17 = Eta
  fPDGCode[fNPDGCodes++]=3122;  // 18 = Lambda
  fPDGCode[fNPDGCodes++]=3222;  // 19 = Sigma +
  fPDGCode[fNPDGCodes++]=3212;  // 20 = Sigma 0
  fPDGCode[fNPDGCodes++]=3112;  // 21 = Sigma -
  fPDGCode[fNPDGCodes++]=3322;  // 22 = Xi0
  fPDGCode[fNPDGCodes++]=3312;  // 23 = Xi-
  fPDGCode[fNPDGCodes++]=3334;  // 24 = Omega-
  fPDGCode[fNPDGCodes++]=-2112; // 25 = Anti Proton
  fPDGCode[fNPDGCodes++]=-3122; // 26 = Anti Proton
  fPDGCode[fNPDGCodes++]=-3222; // 27 = Anti Sigma -
  fPDGCode[fNPDGCodes++]=-3212; // 28 = Anti Sigma 0
  fPDGCode[fNPDGCodes++]=-3112; // 29 = Anti Sigma 0
  fPDGCode[fNPDGCodes++]=-3322; // 30 = Anti Xi 0
  fPDGCode[fNPDGCodes++]=-3312; // 31 = Anti Xi +
  fPDGCode[fNPDGCodes++]=-3334; // 32 = Anti Omega +
}

