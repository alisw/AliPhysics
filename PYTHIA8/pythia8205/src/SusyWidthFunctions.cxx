// SusyWidthFunctions.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand
// Authors: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// SUSY Resonance three-body decay width classes.

#include "Pythia8/SusyResonanceWidths.h"
#include "Pythia8/SusyWidthFunctions.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/PythiaComplex.h"

namespace Pythia8 {

//==========================================================================

// The WidthFunctions class.
// Functions to be integrated for calculating the 3-body decay widths.

//--------------------------------------------------------------------------

void WidthFunction::setPointers( ParticleData* particleDataPtrIn,
  CoupSUSY* coupSUSYPtrIn, Info* infoPtrIn) {

  particleDataPtr = particleDataPtrIn;
  coupSUSYPtr = coupSUSYPtrIn;
  infoPtr = infoPtrIn;
}

//--------------------------------------------------------------------------

double WidthFunction::function(double) {

  cout << "Warning using dummy width function" << endl;
  return 0.;
}

//--------------------------------------------------------------------------

// Adapted from the CERNLIB DGAUSS routine by K.S. Kolbig.

double WidthFunction::integrateGauss(double xlo, double xhi, double tol) {

  // 8-point unweighted.
  static double x8[4]={0.96028985649753623,
                       0.79666647741362674,
                       0.52553240991632899,
                       0.18343464249564980};
  static double w8[4]={0.10122853629037626,
                       0.22238103445337447,
                       0.31370664587788729,
                       0.36268378337836198};
  // 16-point unweighted.
  static double x16[8]={0.98940093499164993,
                       0.94457502307323258,
                       0.86563120238783174,
                       0.75540440835500303,
                       0.61787624440264375,
                       0.45801677765722739,
                       0.28160355077925891,
                       0.09501250983763744};
  static double w16[8]={0.027152459411754095,
                       0.062253523938647893,
                       0.095158511682492785,
                       0.12462897125553387,
                       0.14959598881657673,
                       0.16915651939500254,
                       0.18260341504492359,
                       0.18945061045506850};
  // Boundary checks.
  if (xlo == xhi) {
    cerr<<"xlo = xhi"<<endl;
    return 0.0;
  }
  if ( xlo > xhi ) {
    cerr<<" (integrateGauss:) -> xhi < xlo"<<endl;
    return 0.0;
  }
  // Initialize.
  double sum = 0.0;
  double c = 0.001/abs(xhi-xlo);
  double zlo = xlo;
  double zhi = xhi;
    
  bool nextbin = true;
  
  while ( nextbin ) {
    
    double zmi = 0.5*(zhi+zlo); // midpoint
    double zmr = 0.5*(zhi-zlo); // midpoint, relative to zlo
    
    // Calculate 8-point and 16-point quadratures.
    double s8=0.0;
    for (int i=0;i<4;i++) {
      double dz = zmr * x8[i];
      s8 += w8[i]*(function(zmi+dz) + function(zmi-dz));
    }
    s8 *= zmr;
    double s16=0.0;
    for (int i=0;i<8;i++) {
      double dz = zmr * x16[i];
      s16 += w16[i]*(function(zmi+dz) + function(zmi-dz));
    }
    s16 *= zmr;
    if (abs(s16-s8) < tol*(1+abs(s16))) {
      // Precision in this bin OK, add to cumulative and go to next.
      nextbin=true;
      sum += s16;
      // Next bin: LO = end of current, HI = end of integration region.
      zlo=zhi;
      zhi=xhi;
      if ( zlo == zhi ) nextbin = false;
    } else {
      // Precision in this bin not OK, subdivide.
      if (1.0 + c*abs(zmr) == 1.0) {
        cerr << " (integrateGauss:) too high accuracy required"<<endl;
        sum = 0.0 ;
        break;
      }
      zhi=zmi;
      nextbin=true;
    }
  }
  return sum;
}

//==========================================================================

// The StauWidths class.
// Width functions for 3-body stau decays.

//--------------------------------------------------------------------------

double StauWidths::getWidth(int idResIn, int idIn){

  setChannel(idResIn, idIn);
    
  // Calculate integration limits and return integrated width.
  if (idResIn % 2 == 0) return 0.0; 
  double width = integrateGauss(0.0,1.0,1.0e-3);
  return width;

}

//--------------------------------------------------------------------------

void StauWidths::setChannel(int idResIn, int idIn) {

  // Common masses.
  idRes = abs(idResIn);
  idIn = abs(idIn);
  mRes = particleDataPtr->m0(idRes); 
  m1 = particleDataPtr->m0(1000022);
  m2 = particleDataPtr->m0(idIn);

  mInt = particleDataPtr->m0(15); 
  gammaInt = particleDataPtr->mWidth(15);

  // Couplings etc.
  f0 = 92.4; //MeV
  gf =   coupSUSYPtr->GF();
  delm = mRes - m1;
  cons = pow2(f0)*pow2(gf)*(pow2(delm) - pow2(m2)) 
       * coupSUSYPtr->V2CKMid(1,1) / (128.0 * pow(mRes*M_PI,3));
  
  if (idIn == 900111) wparam = 1.16;
  else if (idIn == 113) wparam = 0.808;
  else wparam = 1.0;

  double g = coupSUSYPtr->alphaEM(mRes * mRes);
  int ksusy = 1000000;
  int isl = (abs(idRes)/ksusy == 2) ? (abs(idRes)%10+1)/2 + 3
                                    : (abs(idRes)%10+1)/2;

  gL = g * coupSUSYPtr->LsllX[isl][3][1] / ( sqrt(2.0) * coupSUSYPtr->cosb);
  gR = g * coupSUSYPtr->RsllX[isl][3][1] / ( sqrt(2.0) * coupSUSYPtr->cosb);

  // Set function switch and internal propagators depending on decay product.
  if (idIn == 111) fnSwitch = 1;
  else if (idIn == 900111 || idIn == 113)  fnSwitch = 2;
  else if (idIn == 14 || idIn == 12) {
    m2 = particleDataPtr->m0(idIn-1);
    fnSwitch = 3;
  } 
  else { 
    stringstream mess;
    mess <<  " unknown decay channel idIn = " << idIn;
    infoPtr->errorMsg("Warning in StauWidths::setChannel:", mess.str() );
  }
  
  return;
}

//------------------------------------------------------------------------

double StauWidths::function(double x){

  // Decay width functions documented in arXiv:1212.2886 Citron et. al.
  double value = 0.0;
  double qf2 = pow2(delm) - (pow2(delm) - pow2(m2)) * x;
  double fac = 1.0 / pow3(mRes);
  double term3 = (norm(gL) * qf2 + norm(gR) * mInt * mInt) 
               * (delm * delm + 2.0 * m1 * delm - qf2);
  double term4 = -2.0 * real(gL * conj(gR)) * m2 * mInt * qf2;

  //  ~tau -> pi0 nu_tau ~chi0.
  if (fnSwitch == 1 ) {
    fac *= pow2(delm) - pow2(m2);
    double term1 = sqrt((pow2(delm) - qf2) * (pow2(delm + 2 * m1) - qf2));
    double term2 = pow2(qf2 - pow2(m2)) / qf2 / (pow2(qf2 - pow2(mInt)) 
      + pow2(mInt*gammaInt));
    value = fac * (term1 * term2 * (term3 + term4));
  }

  else if (fnSwitch == 2 ) {
    double term1 = sqrt((pow2(delm) - qf2) * (pow2(delm + 2 * m1) - qf2));
    double term2 = pow2(qf2 - pow2(m2)) * (qf2 + pow2(m2)) 
      / (qf2 * qf2 * (pow2(qf2 - pow2(mInt)) + pow2(mInt*gammaInt)));
    value = fac * (term1 * term2 * (term3 + term4)); 
  }

  else if (fnSwitch == 3 ) {
    double qf4 = qf2 * qf2;
    double m24 = pow2(m2*m2);
    double term1 = sqrt((pow2(delm) - qf2) * (pow2(delm + 2 * m1) - qf2));
    double term2a = 1.0 / (pow2(qf2 - pow2(mInt)) + pow2(mInt*gammaInt)) / qf4;
    double term2b = 12 * m24 * qf4 * log(qf2/pow2(m2)) 
      + (qf4 - m24) * (qf4 - 8 * m2 * m2 * qf2 + m24);
    value = fac * (term1 * term2a * term2b * (term3 + term4)); 
  }

  else {
    stringstream mess;
    mess <<  " unknown decay channel fnSwitch = " << fnSwitch;
    infoPtr->errorMsg("Warning in StauWidths::function:", mess.str() );
  }

  return value;

}

//==========================================================================

} //end namespace Pythia8

