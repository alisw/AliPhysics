// HardDiffraction.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
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

// Proton mass
const double HardDiffraction::PROTONMASS = 0.938;

//--------------------------------------------------------------------------

void HardDiffraction::init(Info* infoPtrIn, Settings& settingsPtrIn,
  Rndm* rndmPtrIn, BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  BeamParticle* beamPomAPtrIn, BeamParticle* beamPomBPtrIn) {

  // Store pointers.
  infoPtr     = infoPtrIn;
  settings    = settingsPtrIn;
  rndmPtr     = rndmPtrIn;
  beamAPtr    = beamAPtrIn;
  beamBPtr    = beamBPtrIn;
  beamPomAPtr = beamPomAPtrIn;
  beamPomBPtr = beamPomBPtrIn;

  // Set diffraction parameters.
  pomSet  = settings.mode("PDF:PomSet");
  pomFlux = settings.mode("Diffraction:PomFlux");

  // Read out some properties of beams to allow shorthand.
  idA = (beamAPtr != 0) ? beamAPtr->id() : 0;
  idB = (beamBPtr != 0) ? beamBPtr->id() : 0;
  mA  = (beamAPtr != 0) ? beamAPtr->m()  : 0.;
  mB  = (beamBPtr != 0) ? beamBPtr->m()  : 0.;

  // Set up Pomeron flux constants.
  if (pomFlux == 1) {
    double sigmaRefPomP = settings.parm("Diffraction:sigmaRefPomP");
    normPom = pow2(sigmaRefPomP) * 0.02;
    b0      = 2.3;
    ap      = 0.25;
  } else if (pomFlux == 2) {
    normPom = 1/2.3;
    A1      = 6.38;
    A2      = 0.424;
    a1      = 8.;
    a2      = 3.;
  } else if (pomFlux == 3) {
    normPom = 1.99;
    a0      = 1.085;
    ap      = 0.25;
    a1      = 4.7;
  } else if (pomFlux == 4) {
    normPom = 0.74;
    ap      = 0.06;
    a0      = 1.1182;
    A1      = 0.27;
    a1      = 8.38;
    A2      = 0.56;
    a2      = 3.78;
    A3      = 0.18;
    a3      = 1.36;
  } else if (pomFlux == 5) {
    normPom = 0.858;
    A1      = 0.9;
    a1      = 4.6;
    A2      = 0.1;
    a2      = 0.6;
    a0      = 1.104;
    ap      = 0.25;
  } else if (pomFlux == 6 || pomFlux == 7) {
    normPom = 1.57285;
    ap      = 0.06;
    b0      = 5.5;
    if (pomFlux == 6) a0 = 1.1182;
    else a0 = 1.1110;
  }

  // Initialise Pomeron values to zero.
  xPomA = tPomA = thetaPomA = 0.;
  xPomB = tPomB = thetaPomB = 0.;

  // Done.
}

//--------------------------------------------------------------------------

bool HardDiffraction::isDiffractive( int iBeamIn, int partonIn, double xIn,
  double Q2In, double xfIncIn) {

  // Store incoming values.
  iBeam        = iBeamIn;
  int parton   = partonIn;
  double x     = xIn;
  double Q2    = Q2In;
  double xfInc = xfIncIn;
  tmpPDFPtr    = (iBeam == 1) ? beamPomAPtr : beamPomBPtr;

  // Return false if value of inclusive PDF is zero.
  if (xfInc < TINYPDF) {
    infoPtr->errorMsg("Warning in HardDiffraction::isDiffractive: "
      "inclusive PDF is zero");
    return false;
  }

  // Generate an xNow according to 1 / x.
  double xNow = pow(xIn, rndmPtr->flat());

  // Overestimated function:
  // g(xP) = c/xP * x*f_{P/p}(x) * x*f_{i/P}(x, Q^2_{max})
  // G(x)  = int_x^1 dxP g(xP)
  //       = c * x*f_{P/p}(x) * x*f_{i/P}(x, Q^2_{max}) * log( 1/x )
  // True function:
  // f(xP) = 1/xP * xP *f_{P/p}(xP) * x/xP * f_{i/P}(x/xP, Q^2)
  // F(x)  = int_x^1 dxP f(xP)
  // We have diffraction if G(x)/xfInc > R and if f(xP)/g(xP) > R =>
  // R < (G(x) * f(xP)) / (g(xP) * xfInc) = (log(1/x) * f(xP)) / (xP * xfInc)
  double over = log(1./x) * xfPom(xNow) * tmpPDFPtr->xf(parton, x/xNow, Q2);
  if (over > xfInc) {
    stringstream msg;
    msg << " Weight above unity with parton " << parton;
    infoPtr->errorMsg("Warning in HardDiffraction::isDiffractive:", msg.str());
  }

  // Discard if overestimate/inclusive PDF is less than random number.
  if (over/xfInc < rndmPtr->flat()) return false;

  // Make sure there is momentum left for beam remnant
  double m2Diff  = xNow * pow2( infoPtr->eCM());
  double mDiff   = sqrt(m2Diff);
  double mDiffA  = (iBeam == 1) ? 0. : PROTONMASS;
  double mDiffB  = (iBeam == 2) ? 0. : PROTONMASS;
  double m2DiffA = mDiffA * mDiffA;
  double m2DiffB = mDiffB * mDiffB;
  double eDiff   = (iBeam == 1) ? 0.5 * (m2Diff + m2DiffA - m2DiffB) / mDiff :
    0.5 * (m2Diff + m2DiffB - m2DiffA) / mDiff;
  if ( 1. - x / xNow < POMERONMASS / eDiff) {
    infoPtr->errorMsg("Warning in HardDiffraction::isDiffractive: "
    "No momentum left for beam remnant.");
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

  // Schuler-Sjöstrand Pomeron flux, see Phys. Rev. D.49 (1994) 2259.
  // flux = normPom * 1/x * exp(2t(2.3 + 0.25 * log(1/x)))
  // => x * flux = normPom * exp(2t(2.3 + 0.25*log(1/x)))
  if (pomFlux == 1) {
    double b = b0 + ap * log(1./x);
    xFlux = normPom/(2.*b) * ( exp(2.*b*tMax) - exp(2.*b*tMin));
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
    xFlux = normPom * exp(log(1./x) * (2.*a0 - 2.));
    xFlux *= (exp(b*tMax) - exp(b*tMin))/b;
  }

  // Donnachie-Landshoff Pomeron flux, see Phys. Lett. B 191 (1987) 309.
  // flux = beta^2(0)/(16 pi) x^(1 - 2*alpha(t)) F_1(t)^2 with
  // F_1(t)^2 = 0.27 exp(8.38 t) + 0.56 exp(3.78 t) + 0.18 exp(1.36 t)
  //          = (4m_p^2-2.8t)^2/(4m_p^2-t)^2*(1/(1-t/0.7))^4
  // => x * flux = beta^2(0)/(16 pi) * x^(2 - 2*\alpha(t)) F_1(t)^2
  else if (pomFlux == 4) {
    double Q = 2. * ap * log(1./x);
    xFlux = normPom * exp(log(1./x) * (2.*a0 - 2.));
    xFlux *= (A1/(Q + a1) * (exp((Q + a1)*tMax) - exp((Q + a1)*tMin))
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
    xFlux = normPom * exp(log(1./x) * ( 2.*a0 - 2.));
    xFlux *= (A1/(Q + a1) * (exp((Q + a1)*tMax) - exp((Q + a1)*tMin))
            + A2/(Q + a2) * (exp((Q + a2)*tMax) - exp((Q + a2)*tMin)));
  }

  // H1 Pomeron flux, see Eur. Phys. J. C48 (2006) 715, ibid. 749
  // flux = normPom * exp(B_Pom*t)/x^(2*\alpha(t)-1)
  // => x * flux = normPom * exp(B_Pom * t) / x^(2*\alpha(t)-2)
  else if (pomFlux == 6 || pomFlux == 7) {
    double b = b0 + 2. * ap * log(1./x);
    xFlux = normPom * exp(log(1./x) * ( 2.*a0 - 2.));
    xFlux *= (exp(b*tMax) - exp(b*tMin))/b;
  }

  // Done
  return xFlux;
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

  // Schuler-Sjöstrndm Pomeron flux, see Phys. Rev. D.49 (1994) 2259.
  if (pomFlux == 1) {
    double b = b0 + ap * log(1./xIn);
    tTmp = log( rndm*exp(2.*b*tMin) + (1. - rndm)*exp(2.*b*tMax))/(2.*b);
  }

  // Bruni-Ingelman Pomeron flux, see Phys. Lett. B311 (1993) 317.
  else if (pomFlux == 2) {
    double prob1 = A1/a1 * (exp(a1*tMax) - exp(a1*tMin));
    double prob2 = A2/a2 * (exp(a2*tMax) - exp(a2*tMin));
    prob1 /= (prob1 + prob2);
    tTmp  = (prob1 > rndmPtr->flat())
      ? log( rndm * exp(a1*tMin) + (1. - rndm) * exp(a1*tMax))/a1
      : log( rndm * exp(a2*tMin) + (1. - rndm) * exp(a2*tMax))/a2;
  }

  // Streng-Berger Pomeron flux, see Comp. Phys. Comm. 86 (1995) 147.
  else if (pomFlux == 3) {
    double b = (2. * ap * log(1./xIn) + a1);
    tTmp = log( rndm * exp(b*tMin) + (1. - rndm) * exp(b*tMax))/b;
  }

  // Donnachie-Landshoff Pomeron flux, see Phys. Lett. B 191 (1987) 309.
  else if (pomFlux == 4) {
    double b1 = 2. * ap * log(1./xIn) + a1;
    double b2 = 2. * ap * log(1./xIn) + a2;
    double b3 = 2. * ap * log(1./xIn) + a3;
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
    double b1 = a1 + 2. * ap * log(1./xIn);
    double b2 = a2 + 2. * ap * log(1./xIn);
    double prob1 = A1/b1 * (exp(b1*tMax) - exp(b1*tMin));
    double prob2 = A2/b2 * (exp(b2*tMax) - exp(b2*tMin));
    prob1 /= (prob1 + prob2);
    tTmp  = (prob1 > rndmPtr->flat())
      ? log( rndm * exp(b1*tMin) + (1. - rndm) * exp(b1*tMax))/b1
      : log( rndm * exp(b2*tMin) + (1. - rndm) * exp(b2*tMax))/b2;
  }

  // H1 Pomeron flux, see Eur. Phys. J. C48 (2006) 715, ibid. 749
  else if (pomFlux == 6 || pomFlux == 7){
    double b = b0 + 2. * ap * log(1./xIn);
    tTmp = log( rndm * exp(b*tMin) + (1. - rndm) * exp(b*tMax))/b;
  }

  // Done.
  return tTmp;
}

//--------------------------------------------------------------------------

// Return x*f_P/p( x, t), ie. Pomeron flux inside proton differential in t.

double HardDiffraction::xfPomWithT(double xIn, double tIn) {

  // Initial values.
  double x     = xIn;
  double t     = tIn;
  double xFlux = 0.;

  // Schuler-Sjöstrand Pomeron flux, see Phys. Rev. D.49 (1994) 2259.
  if (pomFlux == 1) {
    double b = b0 + ap * log(1./x);
    xFlux = normPom * exp( 2.*b*t);
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
    xFlux = normPom * pow(x, 2. +  2. * (a0 + ap*t)) * sqrF1;
  }

  // MBR Pomeron flux, see arXiv:1205.1446v2 [hep-ph] 2012.
  else if (pomFlux == 5) {
    double sqrF1 = A1 * exp(a1*t) + A2 * exp(a2*t);
    xFlux = normPom * sqrF1 * exp(log(1./x) * (-2. + a0 + ap*t));
  }

  // H1 Pomeron flux, see Eur. Phys. J. C48 (2006) 715, ibid. 749
  else if (pomFlux == 6 || pomFlux == 7)
    xFlux = normPom * exp(b0*t)/pow(x, 2. * (a0 + ap*t) - 2.);

  // Done
  return xFlux;
}

//--------------------------------------------------------------------------

// Set up t range. See p. 113 of 6.4 manual.

pair<double, double> HardDiffraction::tRange(double xIn) {

  // Set up diffractive masses.
  double eCM = infoPtr->eCM();
  s          = eCM * eCM;
  s1         = pow2(mA);
  s2         = pow2(mB);
  s3         = (iBeam == 1) ? xIn * s : s1;
  s4         = (iBeam == 2) ? xIn * s : s2;

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
  double eCM = infoPtr->eCM();
  s          = eCM * eCM;
  s1         = pow2( mA);
  s2         = pow2( mB);
  s3         = (iBeam == 1) ? xIn * s : s1;
  s4         = (iBeam == 2) ? xIn * s : s2;

  // Find theta from the chosen t.
  double lambda12 = sqrtpos(pow2(s - s1 - s2) - 4. * s1 * s2);
  double lambda34 = sqrtpos(pow2(s - s3 - s4) - 4. * s3 * s4);
  double tmp1     = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4)/s;
  double tmp2     = lambda12 * lambda34 / s;
  double tmp3     = (s1 + s4 - s2 - s3) * (s1 * s4 - s2 * s3) / s
    + (s3 - s1) * (s4 - s2);
  double cosTheta = min(1., max(-1., (tmp1 + 2. * tIn) / tmp2));
  double sinTheta = 2. * sqrtpos( -(tmp3 + tmp1 * tIn + tIn * tIn) ) / tmp2;
  double theta = asin( min(1., sinTheta));
  if (cosTheta < 0.) theta = M_PI - theta;

  // Done.
  return theta;
}

//==========================================================================

} // end namespace Pythia8
