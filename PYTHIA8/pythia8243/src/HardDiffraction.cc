// HardDiffraction.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Christine O. Rasmussen.

// Function definitions (not found in the header) for the
// HardDiffraction class.

#include "Pythia8/HardDiffraction.h"
namespace Pythia8 {

//==========================================================================

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Lower limit on PDF value in order to avoid division by zero.
const double HardDiffraction::TINYPDF = 1e-10;

// Ficticious Pomeron mass to leave room for beam remnant
const double HardDiffraction::POMERONMASS = 1.;
const double HardDiffraction::RHOMASS     = 0.77549;
const double HardDiffraction::PROTONMASS  = 0.93827;

// Safetymargin for diffractive masses
const double HardDiffraction::DIFFMASSMARGIN = 0.2;
//--------------------------------------------------------------------------

void HardDiffraction::init(Info* infoPtrIn, Settings& settingsPtrIn,
  Rndm* rndmPtrIn, BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  BeamParticle* beamPomAPtrIn, BeamParticle* beamPomBPtrIn,
  SigmaTotal* sigTotPtrIn) {

  // Store pointers.
  infoPtr     = infoPtrIn;
  settings    = settingsPtrIn;
  rndmPtr     = rndmPtrIn;
  beamAPtr    = beamAPtrIn;
  beamBPtr    = beamBPtrIn;
  beamPomAPtr = beamPomAPtrIn;
  beamPomBPtr = beamPomBPtrIn;
  sigTotPtr   = sigTotPtrIn;

  // Set diffraction parameters.
  pomFlux     = settings.mode("SigmaDiffractive:PomFlux");

  // Read out some properties of beams to allow shorthand.
  idA = (beamAPtr != 0) ? beamAPtr->id() : 0;
  idB = (beamBPtr != 0) ? beamBPtr->id() : 0;
  mA  = (beamAPtr != 0) ? beamAPtr->m()  : 0.;
  mB  = (beamBPtr != 0) ? beamBPtr->m()  : 0.;
  isGammaA = (beamAPtr != 0) ? beamAPtr->isGamma() : false;
  isGammaB = (beamBPtr != 0) ? beamBPtr->isGamma() : false;
  isGammaGamma = (isGammaA && isGammaB);

  // Set up Pomeron flux constants.
  rescale = settings.parm("Diffraction:PomFluxRescale");
  a0      = 1. + settings.parm("SigmaDiffractive:PomFluxEpsilon");
  ap      = settings.parm("SigmaDiffractive:PomFluxAlphaPrime");

  if (pomFlux == 1) {
    double sigmaRefPomP = settings.parm("Diffraction:sigmaRefPomP");
    normPom = pow2(sigmaRefPomP) * 0.02;
    b0      = 2.3;
  } else if (pomFlux == 2) {
    normPom = 1/2.3;
    A1      = 6.38;
    A2      = 0.424;
    a1      = 8.;
    a2      = 3.;
  } else if (pomFlux == 3) {
    double beta = 10.;
    normPom = pow2(beta)/(16.*M_PI);
    a1      = 4.7;
  } else if (pomFlux == 4) {
    double beta = 1.8;
    normPom = 9. * pow2(beta) / (4. * pow2(M_PI));
    A1      = 0.27;
    a1      = 8.38;
    A2      = 0.56;
    a2      = 3.78;
    A3      = 0.18;
    a3      = 1.36;
  } else if (pomFlux == 5) {
    A1      = 0.9;
    a1      = 4.6;
    A2      = 0.1;
    a2      = 0.6;
    a0      = 1. + settings.parm("SigmaDiffractive:MBRepsilon");
    ap      = settings.parm("SigmaDiffractive:MBRalpha");
    bool renormalize   = settings.flag("Diffraction:useMBRrenormalization");
    double cflux       = 0.858;
    double m2min       = settings.parm("SigmaDiffractive:MBRm2Min");
    double dyminSDflux = settings.parm("SigmaDiffractive:MBRdyminSDflux");
    double dymaxSD     = log(infoPtr->eCM()*infoPtr->eCM() / m2min);
    double nGap        = 0.;
    if (renormalize){
      double step        = (dymaxSD - dyminSDflux) / 1000.;
      // Calculate the integral of the flux
      // to renormalize the gap:
      for (int i = 0; i < 1000; ++i) {
        double dy = dyminSDflux + (i + 0.5) * step;
        double f  = exp(2.*(a0 - 1.)*dy) * ( (A1/(a1 + 2.*ap*dy))
                     + (A2/(a2 + 2. * ap*dy)) );
        nGap      += step * cflux * f;
      }
    }
    if (nGap < 1.) nGap = 1.;
    normPom = cflux/nGap;
  } else if (pomFlux == 6 || pomFlux == 7) {
    // Has fixed values of eps and alpha' to get normalisation correct
    ap = 0.06;
    b0 = 5.5;
    if (pomFlux == 6) a0 = 1.1182;
    else a0 = 1.1110;
    double xNorm = 0.003;
    double b     = b0 + 2. * ap * log(1./xNorm);
    double mMin  = (isGammaA || isGammaB) ? RHOMASS : PROTONMASS;
    double tmin  = -pow(mMin * xNorm, 2.)/(1. - xNorm);
    double tcut  = -1.;
    double fNorm = exp(log(1./xNorm) * ( 2.*a0 - 2.));
    fNorm       *= (exp(b*tmin) - exp(b*tcut))/b;
    normPom      = 1./fNorm;
  }

  // Initialise Pomeron values to zero.
  xPomA = tPomA = thetaPomA = 0.;
  xPomB = tPomB = thetaPomB = 0.;

  // Calculate rescaling factor for Pomeron flux in photons.
  sigTotRatio = 1.;
  if (isGammaA || isGammaB) {
    sigTotPtr->calc(22, 2212, infoPtr->eCM());
    double sigGamP = sigTotPtr->sigmaTot();
    sigTotPtr->calc(2212, 2212, infoPtr->eCM());
    double sigPP = sigTotPtr->sigmaTot();
    sigTotRatio = sigGamP / sigPP;
  }

  // Done.
}

//--------------------------------------------------------------------------

bool HardDiffraction::isDiffractive( int iBeamIn, int partonIn,
  double xIn, double Q2In, double xfIncIn) {

  // iBeam = 1 means A B -> A' + (Pom + B) -> A' X
  // iBeam = 2 means A B -> (A + Pom) + B' -> X B'

  // Store incoming values.
  iBeam          = iBeamIn;
  int parton     = partonIn;
  double x       = xIn;
  double Q2      = Q2In;
  double xfInc   = xfIncIn;
  tmpPomPtr      = (iBeam == 1) ? beamPomAPtr : beamPomBPtr;
  usePomInPhoton = ((iBeam == 1 && isGammaA) || (iBeam == 2 && isGammaB))
                 ? true : false;

  // Return false if value of inclusive PDF is zero.
  if (xfInc < TINYPDF) {
    infoPtr->errorMsg("Warning in HardDiffraction::isDiffractive: "
      "inclusive PDF is zero");
    return false;
  }

  // Generate an xNow = x_P according to dx_P / x_P.
  double xNow = pow(x, rndmPtr->flat());

  // Find estimate of diffractive PDF based on x_P choice above.
  // f_i(xP) = int_x^1 d(xP)/xP * xP f_{P/p}(xP) * x/xP f_{i/P}(x/xP, Q^2)
  //         = ln (1/x) < xP f_{P/p}(xP) * x/xP f_{i/P}(x/xP, Q^2) >
  double xfEst = log(1./x) * xfPom(xNow) * tmpPomPtr->xf(parton, x/xNow, Q2);

  // Warn if the estimated function exceeds the inclusive PDF.
  if (xfEst > xfInc) {
    ostringstream msg;
    msg << ", id = " << parton;
    infoPtr->errorMsg("Warning in HardDiffraction::isDiffractive: "
      "weight above unity", msg.str());
  }
  // Discard if estimate/inclusive PDF is less than random number.
  if (xfEst < rndmPtr->flat() * xfInc) return false;

  // Make sure there is momentum left for beam remnant.
  // Make Pomeron massless and let mMin be the mass of the subcollision
  // particle: if Pp then mMin = mp, if Pgamma then mMin = mrho.
  double mMin    = (usePomInPhoton) ? RHOMASS : PROTONMASS;
  double m2Diff  = xNow * pow2( infoPtr->eCM());
  double mDiff   = sqrt(m2Diff);
  double mDiffA  = (iBeam == 1) ? 0. : mMin;
  double mDiffB  = (iBeam == 2) ? 0. : mMin;
  double m2DiffA = mDiffA * mDiffA;
  double m2DiffB = mDiffB * mDiffB;
  double eDiff   = (iBeam == 1)
    ? 0.5 * (m2Diff + m2DiffA - m2DiffB) / mDiff
    : 0.5 * (m2Diff + m2DiffB - m2DiffA) / mDiff;
  if ( 1. - x / xNow < POMERONMASS / eDiff) {
    infoPtr->errorMsg("Warning in HardDiffraction::isDiffractive: "
    "No momentum left for beam remnant.");
    return false;
  }

  // Make sure that the diffractive mass is not too high.
  double m3 = ((iBeam == 1 && isGammaA) || (iBeam == 2 && isGammaB))
            ? RHOMASS : PROTONMASS;
  double m4 = mDiff;
  if (m3 + m4 + DIFFMASSMARGIN >= infoPtr->eCM()) {
    infoPtr->errorMsg("Warning in HardDiffraction::isDiffractive: "
    "Too high diffractive mass.");
    return false;
  }

  // The chosen xNow is accepted, now find t and theta.
  double tNow = pickTNow(xNow);
  double thetaNow = getThetaNow(xNow, tNow);

  // Set the chosen diffractive values, to be able to use them later.
  if (iBeam == 1) {
    xPomA     = xNow;
    tPomA     = tNow;
    thetaPomA = thetaNow;
  } else {
    xPomB     = xNow;
    tPomB     = tNow;
    thetaPomB = thetaNow;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Return x*f_P/p( x), ie. Pomeron flux inside proton, integrated over t.

double HardDiffraction::xfPom(double xIn) {

  // Setup t range.
  pair<double, double> tLim = tRange(xIn);
  double tMin  = tLim.first;
  double tMax  = tLim.second;
  double x     = xIn;
  double xFlux = 0.;

  // Schuler-Sjostrand Pomeron flux, see Phys. Rev. D.49 (1994) 2259.
  // flux = normPom * 1/x * exp(2t(2.3 + 0.25 * log(1/x)))
  // => x * flux = normPom * exp(2t(2.3 + 0.25*log(1/x)))
  if (pomFlux == 1) {
    double b = b0 + ap * log(1./x);
    xFlux    = normPom/(2.*b) * ( exp(2.*b*tMax) - exp(2.*b*tMin));
  }

  // Bruni-Ingelman Pomeron flux, see Phys. Lett. B311 (1993) 317.
  // flux = normPom * (1/x) * (6.38 *exp(8*t)+ 0.424 * exp(3*t))
  // => x * flux = normPom * (6.38 *exp(8*t)+ 0.424 * exp(3*t))
  else if (pomFlux == 2) {
    xFlux = normPom * (A1/a1 * (exp(a1*tMax) - exp(a1*tMin))
                     + A2/a2 * (exp(a2*tMax) - exp(a2*tMin)));
  }

  // Streng-Berger Pomeron flux, see Comp. Phys. Comm. 86 (1995) 147.
  // flux = normPom * x^(1 - 2*alpha(t)) * exp(-R_N^2 * t)
  // => x * flux = normPom * x^(2 - 2*alpha(t)) * exp(-R_N^2 * t)
  else if (pomFlux == 3) {
    double b = (a1 + 2. * ap * log(1./x));
    xFlux    = normPom * exp(log(1./x) * (2.*a0 - 2.));
    xFlux   *= (exp(b*tMax) - exp(b*tMin))/b;
  }

  // Donnachie-Landshoff Pomeron flux, see Phys. Lett. B 191 (1987) 309.
  // flux = 9 beta^2(0)/(4 pi^2) x^(1 - 2*alpha(t)) F_1(t)^2 with
  // F_1(t)^2 = 0.27 exp(8.38 t) + 0.56 exp(3.78 t) + 0.18 exp(1.36 t)
  //          = (4m_p^2-2.8t)^2/(4m_p^2-t)^2*(1/(1-t/0.7))^4
  // => x * flux = 9 beta^2(0)/(4 pi^2) * x^(2 - 2*\alpha(t)) F_1(t)^2
  else if (pomFlux == 4) {
    double Q = 2. * ap * log(1./x);
    xFlux    = normPom * exp(log(1./x) * (2.*a0 - 2.));
    xFlux   *= (A1/(Q + a1) * (exp((Q + a1)*tMax) - exp((Q + a1)*tMin))
              + A2/(Q + a2) * (exp((Q + a2)*tMax) - exp((Q + a2)*tMin))
              + A3/(Q + a3) * (exp((Q + a3)*tMax) - exp((Q + a3)*tMin)));
  }

  // MBR Pomeron flux, see arXiv:1205.1446v2 [hep-ph] 2012.
  // flux = normPom * F_1(t)^2 * exp((2 * alpha(t) - 1)*log(1/x))
  // F_1(t)^2 = 0.9 exp(4.6 t) + 0.1 exp(0.6 t)
  //          = (4m_p^2-2.8t)^2/(4m_p^2-t)^2*(1/(1-t/0.7))^4
  // => x * flux = normPom * F_1(t)^2 * exp( 2*(alpha(t) -1)*log(1/x))
  else if (pomFlux == 5) {
    double Q = 2. * ap * log(1./x);
    xFlux    = normPom * exp(log(1./x) * ( 2.*a0 - 2.));
    xFlux   *= (A1/(Q + a1) * (exp((Q + a1)*tMax) - exp((Q + a1)*tMin))
              + A2/(Q + a2) * (exp((Q + a2)*tMax) - exp((Q + a2)*tMin)));
  }

  // H1 Fit A, B Pomeron flux, see Eur. Phys. J. C48 (2006) 715, ibid. 749
  // flux = normPom * exp(B_Pom*t)/x^(2*\alpha(t)-1)
  // => x * flux = normPom * exp(B_Pom * t) / x^(2*\alpha(t)-2)
  else if (pomFlux == 6 || pomFlux == 7) {
    double b = b0 + 2. * ap * log(1./x);
    xFlux    = normPom * exp(log(1./x) * ( 2.*a0 - 2.));
    xFlux   *= (exp(b*tMax) - exp(b*tMin))/b;
  }

  // Done
  if (usePomInPhoton) return xFlux * rescale * sigTotRatio;
  else return xFlux * rescale;

}

//--------------------------------------------------------------------------

// Pick a t value according to different Pomeron flux parametrizations.

double HardDiffraction::pickTNow(double xIn) {

  // Get kinematical limits for t. initial values.
  pair<double, double> tLim = HardDiffraction::tRange(xIn);
  double tMin = tLim.first;
  double tMax = tLim.second;
  double tTmp = 0.;
  double rndm = rndmPtr->flat();

  // Schuler-Sjostrand Pomeron flux, see Phys. Rev. D.49 (1994) 2259.
  if (pomFlux == 1) {
    double b = b0 + ap * log(1./xIn);
    tTmp     = log( rndm*exp(2.*b*tMin) + (1. - rndm)*exp(2.*b*tMax))/(2.*b);
  }

  // Bruni-Ingelman Pomeron flux, see Phys. Lett. B311 (1993) 317.
  else if (pomFlux == 2) {
    double prob1 = A1/a1 * (exp(a1*tMax) - exp(a1*tMin));
    double prob2 = A2/a2 * (exp(a2*tMax) - exp(a2*tMin));
    prob1       /= (prob1 + prob2);
    tTmp         = (prob1 > rndmPtr->flat())
      ? log( rndm * exp(a1*tMin) + (1. - rndm) * exp(a1*tMax))/a1
      : log( rndm * exp(a2*tMin) + (1. - rndm) * exp(a2*tMax))/a2;
  }

  // Streng-Berger Pomeron flux, see Comp. Phys. Comm. 86 (1995) 147.
  else if (pomFlux == 3) {
    double b = (2. * ap * log(1./xIn) + a1);
    tTmp     = log( rndm * exp(b*tMin) + (1. - rndm) * exp(b*tMax))/b;
  }

  // Donnachie-Landshoff Pomeron flux, see Phys. Lett. B 191 (1987) 309.
  else if (pomFlux == 4) {
    double b1       = 2. * ap * log(1./xIn) + a1;
    double b2       = 2. * ap * log(1./xIn) + a2;
    double b3       = 2. * ap * log(1./xIn) + a3;
    double prob1    = A1/b1 * ( exp(b1*tMax) - exp(b1*tMin));
    double prob2    = A2/b2 * ( exp(b2*tMax) - exp(b2*tMin));
    double prob3    = A3/b3 * ( exp(b3*tMax) - exp(b3*tMin));
    double rndmProb = (prob1 + prob2 + prob3) * rndmPtr->flat() ;
    if (rndmProb < prob1)
      tTmp = log( rndm * exp(b1*tMin) + (1. - rndm) * exp(b1*tMax))/b1;
    else if ( rndmProb < prob1 + prob2)
      tTmp = log( rndm * exp(b2*tMin) + (1. - rndm) * exp(b2*tMax))/b2;
    else
      tTmp = log( rndm * exp(b3*tMin) + (1. - rndm) * exp(b3*tMax))/b3;
  }

  // MBR Pomeron flux, see arXiv:1205.1446v2 [hep-ph] 2012.
  else if (pomFlux == 5) {
    double b1    = a1 + 2. * ap * log(1./xIn);
    double b2    = a2 + 2. * ap * log(1./xIn);
    double prob1 = A1/b1 * (exp(b1*tMax) - exp(b1*tMin));
    double prob2 = A2/b2 * (exp(b2*tMax) - exp(b2*tMin));
    prob1       /= (prob1 + prob2);
    tTmp         = (prob1 > rndmPtr->flat())
      ? log( rndm * exp(b1*tMin) + (1. - rndm) * exp(b1*tMax))/b1
      : log( rndm * exp(b2*tMin) + (1. - rndm) * exp(b2*tMax))/b2;
  }

  // H1 Pomeron flux, see Eur. Phys. J. C48 (2006) 715, ibid. 749
  else if (pomFlux == 6 || pomFlux == 7){
    double b = b0 + 2. * ap * log(1./xIn);
    tTmp     = log( rndm * exp(b*tMin) + (1. - rndm) * exp(b*tMax))/b;
  }

  // Done.
  return tTmp;

}

//--------------------------------------------------------------------------

// Return x*f_P/p( x, t), ie. Pomeron flux inside proton differential in t.

double HardDiffraction::xfPomWithT(double xIn, double tIn) {

  // Initial values.
  double x        = xIn;
  double t        = tIn;
  double xFlux    = 0.;

  // Schuler-Sjostrand Pomeron flux, see Phys. Rev. D.49 (1994) 2259.
  if (pomFlux == 1) {
    double b = b0 + ap * log(1./x);
    xFlux    = normPom * exp( 2.*b*t);
  }

  // Bruni-Ingelman Pomeron flux, see Phys. Lett. B311 (1993) 317.
  else if (pomFlux == 2)
    xFlux = normPom * (A1 * exp(a1*t) + A2 * exp(a2*t));

  // Streng-Berger Pomeron flux, see Comp. Phys. Comm. 86 (1995) 147.
  else if (pomFlux == 3) {
    xFlux = normPom * exp(log(1./x) * (2.*a0 - 2.))
                    * exp(t * (a1 + 2.*ap*log(1./x)));
  }

  // Donnachie-Landshoff Pomeron flux, see Phys. Lett. B 191 (1987) 309.
  else if (pomFlux == 4){
    double sqrF1 = A1 * exp(a1*t) + A2 * exp(a2*t) + A3 * exp(a3*t);
    xFlux        = normPom * pow(x, 2. +  2. * (a0 + ap*t)) * sqrF1;
  }

  // MBR Pomeron flux, see arXiv:1205.1446v2 [hep-ph] 2012.
  else if (pomFlux == 5) {
    double sqrF1 = A1 * exp(a1*t) + A2 * exp(a2*t);
    xFlux        = normPom * sqrF1 * exp(log(1./x) * (-2. + a0 + ap*t));
  }

  // H1 Pomeron flux, see Eur. Phys. J. C48 (2006) 715, ibid. 749
  else if (pomFlux == 6 || pomFlux == 7)
    xFlux = normPom * exp(b0*t)/pow(x, 2. * (a0 + ap*t) - 2.);

  // Done
  if (usePomInPhoton) return xFlux * rescale * sigTotRatio;
  else return xFlux * rescale;

}

//--------------------------------------------------------------------------

// Set up t range. See p. 113 of 6.4 manual.

pair<double, double> HardDiffraction::tRange(double xIn) {

  // Set up diffractive masses.
  // s1 = mA^2, s2 = mB^2,
  // s3 = M^2 (= mA^2 if A scatteres elastically)
  // s4 = M^2 (= mB^2 if B scatteres elastically)
  double eCM = infoPtr->eCM();
  s          = eCM * eCM;
  double M2  = xIn * s;
  s1         = pow2(mA);
  s2         = pow2(mB);
  s3         = (iBeam == 1) ? s1 : M2;
  s4         = (iBeam == 2) ? s2 : M2;

  // Calculate kinematics.
  double lambda12 = sqrtpos(pow2(s - s1 - s2) - 4. * s1 * s2);
  double lambda34 = sqrtpos(pow2(s - s3 - s4) - 4. * s3 * s4);
  double tmp1     = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
  double tmp2     = lambda12 * lambda34 / s;
  double tmp3     = (s1 + s4 - s2 - s3) * (s1 * s4 - s2 * s3) / s
                  + (s3 - s1) * (s4 - s2);
  double tMin     = -0.5 * (tmp1 + tmp2);
  double tMax     = tmp3 / tMin;

  // Done.
  return make_pair(tMin, tMax);

}

//--------------------------------------------------------------------------

// Get the scattering angle from the chosen t.

double HardDiffraction::getThetaNow( double xIn, double tIn) {

  // Set up diffractive masses.
  // s1 = mA^2, s2 = mB^2,
  // s3 = M^2 (= mA^2 if A scatteres elastically)
  // s4 = M^2 (= mB^2 if B scatteres elastically)
  double eCM = infoPtr->eCM();
  s          = eCM * eCM;
  double M2  = xIn * s;
  s1         = pow2(mA);
  s2         = pow2(mB);
  s3         = (iBeam == 1) ? s1 : M2;
  s4         = (iBeam == 2) ? s2 : M2;

  // Find theta from the chosen t.
  double lambda12 = sqrtpos(pow2(s - s1 - s2) - 4. * s1 * s2);
  double lambda34 = sqrtpos(pow2(s - s3 - s4) - 4. * s3 * s4);
  double tmp1     = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4)/s;
  double tmp2     = lambda12 * lambda34 / s;
  double tmp3     = (s1 + s4 - s2 - s3) * (s1 * s4 - s2 * s3) / s
                  + (s3 - s1) * (s4 - s2);
  double cosTheta = min(1., max(-1., (tmp1 + 2. * tIn) / tmp2));
  double sinTheta = 2. * sqrtpos( -(tmp3 + tmp1 * tIn + tIn * tIn) ) / tmp2;
  double theta    = asin( min(1., sinTheta));

  if (cosTheta < 0.) theta = M_PI - theta;

  // Done.
  return theta;

}

//==========================================================================

} // end namespace Pythia8
