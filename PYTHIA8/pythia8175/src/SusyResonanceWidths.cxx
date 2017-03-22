// SusyResonanceWidths.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand
// Authors: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// SUSY Resonance widths classes. 

#include "SusyResonanceWidths.h"
#include "SusyCouplings.h"
#include "ParticleData.h"
#include "PythiaComplex.h"

namespace Pythia8 {

//==========================================================================

// WidthFunctions Class 

// Contains functions to be integrated for calculating the 3-body
// decay widths

//--------------------------------------------------------------------------

void WidthFunction::init( ParticleData* particleDataPtrIn, 
  CoupSUSY* coupSUSYPtrIn) {

  particleDataPtr = particleDataPtrIn;
  coupSUSYPtr = coupSUSYPtrIn;

}

//--------------------------------------------------------------------------

void WidthFunction::setInternal2(int idResIn, int id1In, int id2In, 
  int id3In, int idIntIn) {

  // Res -> 1,2,3
  idRes = idResIn;
  id1 = id1In; 
  id2 = id2In; 
  id3 = id3In; 
  idInt = idIntIn; 

  mRes = particleDataPtr->m0(idRes);
  m1 = particleDataPtr->m0(id1);
  m2 = particleDataPtr->m0(id2);
  m3 = particleDataPtr->m0(id3);

  // Internal propagator
  mInt = particleDataPtr->m0(idInt);
  gammaInt = particleDataPtr->mWidth(idInt);

  return;
}

//--------------------------------------------------------------------------

double WidthFunction::function(double) {

  cout<<"Warning using dummy width function"<<endl;
  return 0.;
}

//--------------------------------------------------------------------------

double WidthFunction::function(double,double) {

  cout<<"Warning using dummy width function"<<endl;
  return 0.;
}

//==========================================================================

// Psi, Upsilon and Phi classes.

//--------------------------------------------------------------------------

void Psi::setInternal (int idResIn, int id1In, int id2In, int id3In, 
  int idIntIn, int) {

  setInternal2(idResIn, id1In, id2In, id3In, idIntIn);

  mInt = particleDataPtr->m0(idInt);
  gammaInt = particleDataPtr->mWidth(idInt);
  iX = coupSUSYPtr->typeNeut(idRes);  
  iQ = (id3+1)/2;
  iSq = (idInt>1000000)? 3 + (idInt%10+1)/2 :  (idInt%10+1)/2;
  isSqDown = (idInt % 2 == 1)? true : false;
  m1 = particleDataPtr->m0(id1);
  m2 = particleDataPtr->m0(id2);
  m3 = particleDataPtr->m0(id3);
  return;
}

//--------------------------------------------------------------------------

void Upsilon::setInternal (int idResIn, int id1In, int id2In, int id3In, 
  int idIntIn, int idInt2In) {

  setInternal2(idResIn, id1In, id2In, id3In, idIntIn);

  idInt2 = idInt2In;
  mInt = particleDataPtr->m0(idInt);
  gammaInt = particleDataPtr->mWidth(idInt);
  mInt2 = particleDataPtr->m0(idInt2);
  gammaInt2 = particleDataPtr->mWidth(idInt2);

  iX = coupSUSYPtr->typeNeut(idRes);  
  iQ = (id3+1)/2;
  iSq  = (idInt>1000000)? 3 + (idInt%10+1)/2 :  (idInt%10+1)/2;
  iSq2 = (idInt2>1000000)? 3 + (idInt2%10+1)/2 :  (idInt2%10+1)/2;
  isSqDown = (idIntIn % 2 == 1)? true : false;
  m1 = particleDataPtr->m0(id1);
  m2 = particleDataPtr->m0(id2);
  m3 = particleDataPtr->m0(id3);
  return;
}

//--------------------------------------------------------------------------

void Phi::setInternal (int idResIn, int id1In, int id2In, int id3In, 
  int idIntIn, int idInt2In) {

  setInternal2(idResIn, id1In, id2In, id3In, idIntIn);

  idInt2 = idInt2In;
  mInt = particleDataPtr->m0(idInt);
  gammaInt = particleDataPtr->mWidth(idInt);
  mInt2 = particleDataPtr->m0(idInt2);
  gammaInt2 = particleDataPtr->mWidth(idInt2);

  iX = coupSUSYPtr->typeNeut(idRes);  
  iQ = (id3+1)/2;
  iSq  = (idInt>1000000)? 3 + (idInt%10+1)/2 :  (idInt%10+1)/2;
  iSq2 = (idInt2>1000000)? 3 + (idInt2%10+1)/2 :  (idInt2%10+1)/2;
  isSqDown = (idIntIn % 2 == 1)? true : false;
  m1 = particleDataPtr->m0(id1);
  m2 = particleDataPtr->m0(id2);
  m3 = particleDataPtr->m0(id3);
  return;
}

//--------------------------------------------------------------------------

double Psi::function(double m12sq) {

  double R, factor1, factor2, value;

  // Check that the propagators are offshell
  if (m12sq > pow2(mInt) || abs(m12sq - pow2(mInt)) < gammaInt) return 0; 

  R = 1.0/(pow2(m12sq-pow2(mInt)) + pow2(gammaInt*mInt));
  if (isSqDown){
    factor1 = (norm(coupSUSYPtr->LsddX[iSq][iQ][iX])
	       + norm(coupSUSYPtr->RsddX[iSq][iQ][iX]))*
      (mRes*mRes + m3*m3 - m12sq);
    factor2 = 4.0 * real(coupSUSYPtr->LsddX[iSq][iQ][iX] 
			 * conj(coupSUSYPtr->RsddX[iSq][iQ][iX])) * m3 * mRes;

  } else {
    factor1 = (norm(coupSUSYPtr->LsuuX[iSq][iQ][iX]) 
	       + norm(coupSUSYPtr->RsuuX[iSq][iQ][iX]))*
      (mRes*mRes + m3*m3 - m12sq);
    factor2 = 4.0 * real(coupSUSYPtr->LsuuX[iSq][iQ][iX] 
			 * conj(coupSUSYPtr->RsuuX[iSq][iQ][iX])) * m3 * mRes;
  }

  value = R * (m12sq - m1*m1 - m2*m2) * (factor1+factor2);

  return value;
}


//--------------------------------------------------------------------------

double Upsilon::function(double m12sq) {

  double R1,R2, S, factor1, factor2, value;

  // Check that the propagators are offshell
  if (m12sq > pow2(mInt) || abs(m12sq - pow2(mInt)) < gammaInt) return 0; 
  if (m12sq > pow2(mInt2) || abs(m12sq - pow2(mInt2)) < gammaInt2) return 0; 
  
  R1 = 1.0/(pow2(m12sq-pow2(mInt)) + pow2(gammaInt*mInt));
  R2 =  1.0/(pow2(m12sq-pow2(mInt2)) + pow2(gammaInt2*mInt2));
  S = R1 * R2 * ((m12sq - pow2(mInt)) * (m12sq - pow2(mInt2)) 
		 + gammaInt * mInt * gammaInt2 * mInt2);

  if (isSqDown){
    factor1 = real(coupSUSYPtr->LsddX[iSq][iQ][iX] 
		   * conj(coupSUSYPtr->LsddX[iSq2][iQ][iX])) 
      + real(coupSUSYPtr->RsddX[iSq][iQ][iX] 
	     * conj(coupSUSYPtr->RsddX[iSq2][iQ][iX]));
    factor2 = real(coupSUSYPtr->LsddX[iSq][iQ][iX] 
		   * conj(coupSUSYPtr->RsddX[iSq2][iQ][iX])) 
      + real(coupSUSYPtr->RsddX[iSq][iQ][iX] 
	     * conj(coupSUSYPtr->LsddX[iSq2][iQ][iX]));
  }else{
    factor1 = real(coupSUSYPtr->LsuuX[iSq][iQ][iX]
		   * conj(coupSUSYPtr->LsuuX[iSq2][iQ][iX])) 
      + real(coupSUSYPtr->RsuuX[iSq][iQ][iX]
	     * conj(coupSUSYPtr->RsuuX[iSq2][iQ][iX]));
    factor2 = real(coupSUSYPtr->LsuuX[iSq][iQ][iX]
		   * conj(coupSUSYPtr->RsuuX[iSq2][iQ][iX])) 
      + real(coupSUSYPtr->RsuuX[iSq][iQ][iX]
	     * conj(coupSUSYPtr->LsuuX[iSq2][iQ][iX]));
  }

  value = S * (m12sq - pow2(m1) - pow2(m2)) *
    ( factor1 * (pow2(mRes) + pow2(m3) - m12sq) + 2.0 * factor2 * m3 * mRes);

  //  cout<<"I1: "<<idInt<<" I2:"<<idInt2<<" factor1: "<<factor1
  //        <<" factor2:"<<factor2<<" value:"<<value<<endl;

  return value;

}

//--------------------------------------------------------------------------

double Phi::function(double m12sqIn) {

  m12sq =  m12sqIn;
  // Check that the propagators are offshell
  if (m12sq > pow2(mInt) || abs(m12sq - pow2(mInt)) < gammaInt) return 0; 

  double m23max, m23min, E2, E3;

  E2 = (m12sq - pow2(m1) - pow2(m2))/(2.0 * sqrt(m12sq));
  E3 = (pow2(mRes) - m12sq - pow2(m3))/(2.0 * sqrt(m12sq));
  m23max = pow2(E2+E3) - (sqrt(E2*E2 - m2*m2) - sqrt(E3*E3 - m3*m3)) ;
  m23min = pow2(E2+E3) - (sqrt(E2*E2 - m2*m2) + sqrt(E3*E3 - m3*m3)) ;

  if (E2 < m2 || E3 < m3){
    cout <<"Error in Phi:function"<<endl;
    return 0.;
  }

  double integral2 = integrateGauss(m23min,m23max,1.0e-4);
  return integral2;
}

//--------------------------------------------------------------------------

double Phi::function2(double m23sq) {

  // Check that the propagators are offshell
  if (m23sq > pow2(mInt2) || abs(m23sq - pow2(mInt2)) < gammaInt2) return 0; 

  double R1, R2, S, factor1, factor2, factor3, factor4, value, fac;
  int iQ2;

  R1 = 1.0/(pow2(m12sq-pow2(mInt)) + pow2(gammaInt*mInt));
  R2 = 1.0/(pow2(m12sq-pow2(mInt2)) + pow2(gammaInt2*mInt2));
  S = R1 * R2 * ((m12sq - pow2(mInt)) * (m12sq - pow2(mInt2)) 
    + gammaInt * mInt * gammaInt2 * mInt2);

  fac = 1.0;

  if (isSqDown){
    // Only factor is when both d_i and d_j are near.
    iQ2 = (id1%2 == 1)? (id1+1)/2 : (id2+1)/2;

    if (mRes > mInt2 + particleDataPtr->m0(iQ2)) fac = 0.;

    factor1 = m1 * m3 * real(coupSUSYPtr->LsddX[iSq][iQ][iX] 
			     * conj(coupSUSYPtr->LsddX[iSq2][iQ2][iX]))
      * (m12sq + m23sq - pow2(m1) - pow2(m3));
    factor2 = m1 * mRes * real(coupSUSYPtr->RsddX[iSq][iQ][iX] 
			       * conj(coupSUSYPtr->LsddX[iSq2][iQ2][iX]))
      * (m23sq - pow2(m2) - pow2(m3));
    factor3 = m3 * mRes * real(coupSUSYPtr->LsddX[iSq][iQ][iX] 
			       * conj(coupSUSYPtr->RsddX[iSq2][iQ2][iX]))
      * (m12sq - pow2(m1) - pow2(m2));
    factor4 = m3 * mRes * real(coupSUSYPtr->RsddX[iSq][iQ][iX] 
			       * conj(coupSUSYPtr->RsddX[iSq2][iQ2][iX]))
      * (m12sq * m23sq - pow2(m1 * m3) - pow2(m2 * mRes));

  }else{
    // Factor A: u and d_1
    iQ2 = (id1+1)/2;

    if (mRes > mInt2 + particleDataPtr->m0(iQ2)) fac = 0.;

    factor1 = m1 * m3 * real(coupSUSYPtr->LsuuX[iSq][iQ][iX] 
			     * conj(coupSUSYPtr->LsddX[iSq2][iQ2][iX]))
      * (m12sq + m23sq - pow2(m1) - pow2(m3));
    factor2 = m1 * mRes * real(coupSUSYPtr->RsuuX[iSq][iQ][iX] 
			       * conj(coupSUSYPtr->LsddX[iSq2][iQ2][iX]))
      * (m23sq - pow2(m2) - pow2(m3));
    factor3 = m3 * mRes * real(coupSUSYPtr->LsuuX[iSq][iQ][iX] 
			       * conj(coupSUSYPtr->RsddX[iSq2][iQ2][iX]))
      * (m12sq - pow2(m1) - pow2(m2));
    factor4 = m3 * mRes * real(coupSUSYPtr->RsuuX[iSq][iQ][iX] 
			       * conj(coupSUSYPtr->RsddX[iSq2][iQ2][iX]))
      * (m12sq * m23sq - pow2(m1 * m3) - pow2(m2 * mRes));


    // Factor B: u and d_2; change 1 <=> 2
    iQ2 = (id2+1)/2;

    if (mRes > mInt2 + particleDataPtr->m0(iQ2)) fac = 0.;

    factor1 += m2 * m3 * real(coupSUSYPtr->LsuuX[iSq][iQ][iX] 
			      * conj(coupSUSYPtr->LsddX[iSq2][iQ2][iX]))
      * (m12sq + m23sq - pow2(m2) - pow2(m3));
    factor2 += m2 * mRes * real(coupSUSYPtr->RsuuX[iSq][iQ][iX] 
				* conj(coupSUSYPtr->LsddX[iSq2][iQ2][iX]))
      * (m23sq - pow2(m1) - pow2(m3));
    factor3 += m3 * mRes * real(coupSUSYPtr->LsuuX[iSq][iQ][iX] 
				* conj(coupSUSYPtr->RsddX[iSq2][iQ2][iX]))
      * (m12sq - pow2(m2) - pow2(m1));
    factor4 += m3 * mRes * real(coupSUSYPtr->RsuuX[iSq][iQ][iX] 
				* conj(coupSUSYPtr->RsddX[iSq2][iQ2][iX]))
      * (m12sq * m23sq - pow2(m2 * m3) - pow2(m1 * mRes));
  }

  value = S * (factor1 + factor2 + factor3 + factor4);

  //  cout<<"I1: "<<idInt<<" I2:"<<idInt2<<" factor1: "<<factor1
  //      <<" factor2:"<<factor2<<" value:"<<value<<endl;

  return (fac * value);
}

//--------------------------------------------------------------------------

double Phi::integrateGauss(double xlo, double xhi, double tol) {

  //8-point unweighted
  static double x8[4]={0.96028985649753623,
		       0.79666647741362674,
		       0.52553240991632899, 
		       0.18343464249564980};
  static double w8[4]={0.10122853629037626,
		       0.22238103445337447,
		       0.31370664587788729,
		       0.36268378337836198};
  //16-point unweighted
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
  //boundary checks
  if (xlo == xhi) {
    cerr<<"xlo = xhi"<<endl;
    return 0.0;
  }
  if ( xlo > xhi ) {
    cerr<<" (integrateGauss:) -> xhi < xlo"<<endl;
    return 0.0;
  }
  //initialize
  double sum = 0.0;
  double c = 0.001/abs(xhi-xlo);
  double zlo = xlo;
  double zhi = xhi;
    
  bool nextbin = true;
  
  while ( nextbin ) {
    
    double zmi = 0.5*(zhi+zlo); // midpoint
    double zmr = 0.5*(zhi-zlo); // midpoint, relative to zlo
    
    //calculate 8-point and 16-point quadratures
    double s8=0.0;
    for (int i=0;i<4;i++) {
      double dz = zmr * x8[i];
      s8 += w8[i]*(function2(zmi+dz) + function2(zmi-dz));
    }
    s8 *= zmr;
    double s16=0.0;
    for (int i=0;i<8;i++) {
      double dz = zmr * x16[i];
      s16 += w16[i]*(function2(zmi+dz) + function2(zmi-dz)); 
    }
    s16 *= zmr;
    if (abs(s16-s8) < tol*(1+abs(s16))) { 
      //precision in this bin OK, add to cumulative and go to next
      nextbin=true;
      sum += s16;
      //next bin: LO = end of current, HI = end of integration region.
      zlo=zhi;
      zhi=xhi;
      if ( zlo == zhi ) nextbin = false;
    } else {
      //precision in this bin not OK, subdivide.
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

// The SUSYResonanceWidths Class
// Derived class for SUSY resonances

const bool SUSYResonanceWidths::DBSUSY = false;

//--------------------------------------------------------------------------

bool SUSYResonanceWidths::initBSM(){

  if (couplingsPtr->isSUSY) {
    coupSUSYPtr     = (CoupSUSY *) couplingsPtr;
    return true;
  }

  return false;
}


//--------------------------------------------------------------------------

bool SUSYResonanceWidths::allowCalc(){

  // Check if decay calculations at all possible
  if ( !couplingsPtr->isSUSY ) return false;

  // Next check if SLHA decay tables are ignored (= always do calculation)
  if ( !settingsPtr->flag("SLHA:useDecayTable") ) return true;

  // Next check if decay table was read in via SLHA and takes precedence
  for ( int iDec = 1; iDec < int((coupSUSYPtr->slhaPtr)->decays.size()); ++iDec)
    if ( (coupSUSYPtr->slhaPtr)->decays[iDec].getId() == abs(idRes) ) {  
      if (DBSUSY) cout<<"Using external decay table for:"<<idRes<<endl;
      return false;
    }
  
  // Else we should do the calculation
  return true;
}

//==========================================================================

// The ResonanceSquark class
// Derived class for Squark resonances

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceSquark::initConstants() {

  // Locally stored properties and couplings.
  alpS  = coupSUSYPtr->alphaS(mHat * mHat );
  alpEM = coupSUSYPtr->alphaEM(mHat * mHat);
  s2W   = coupSUSYPtr->sin2W;
}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceSquark::calcPreFac(bool) {

  // Common coupling factors.
  preFac = 1.0 / (s2W * pow(mHat,3));

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceSquark::calcWidth(bool) {

  // Squark type -- in u_i/d_i and generation
  int ksusy = 1000000;
  bool idown = (abs(idRes)%2 == 0 ? false : true);
  int isq = (abs(idRes)/ksusy == 2) ? 
    (abs(idRes)%10+1)/2 + 3: (abs(idRes)%10+1)/2;
  // int isqgen = (abs(idRes)%10 + 1)/2;

  // Check that mass is above threshold.
  if (ps == 0.) return;
  else{
    // Two-body decays
    kinFac = (mHat * mHat - mf1 * mf1 - mf2 * mf2);
    
    double fac = 0.0 , wid = 0.0;
  
    //RPV decays
    //Case 1a:  UDD-type 
    if (id1Abs < 7 && id2Abs < 7){
      
      // Quark generations
      int iq1 = (id1Abs + 1)/2;
      int iq2 = (id2Abs + 1)/2;
      
      // Check for RPV UDD couplings
      if (!coupSUSYPtr->isUDD) {widNow = 0; return;}
      
      // ~q -> q_i + q_j
      
      fac = 2.0 * kinFac / (16.0 * M_PI * pow(mHat,3)); 
      wid = 0.0;
      if (idown) {
	if ((id1Abs+id2Abs)%2 == 1){
	  if (id1Abs%2==1)
	    for (int isq2 = 1; isq2 < 4; isq2++)
	      wid += norm(coupSUSYPtr->rvUDD[iq2][iq1][isq2] 
                   * coupSUSYPtr->Rdsq[isq][isq2+3]);
	  else
	    for (int isq2 = 1; isq2 < 4; isq2++)
	      wid += norm(coupSUSYPtr->rvUDD[iq1][iq2][isq2] 
                   * coupSUSYPtr->Rdsq[isq][isq2+3]);
	}
      }
      else {
	if ((id1Abs+id2Abs)%2 != 0) widNow = 0.0;
	else
	  for (int isq2 = 1; isq2 < 4; isq2++)
	    wid += norm(coupSUSYPtr->rvUDD[isq2][iq1][iq2] 
                 * coupSUSYPtr->Rusq[isq][isq2+3]);
      }
  }
    
    //Case 1b:  LQD-type 
    else if (id1Abs < 17 && id2Abs < 7){
      if (!coupSUSYPtr->isLQD) {widNow = 0; return;}
      
      int ilep = (id1Abs - 9)/2;
      int iq = (id2Abs + 1)/2;
    
      fac = kinFac / (16.0 * M_PI * pow(mHat,3)); 
      wid = 0.0;
      if (idown){
	if (iq%2 == 0){
	  // q is up-type; ~q is right-handed down type
	  for (int isq2=1; isq2<3; isq2++)
	    wid += norm(coupSUSYPtr->Rdsq[isq][isq2+3] 
                 * coupSUSYPtr->rvLQD[ilep][iq][isq2]);
	}else{
	  //q is down type; ~q left-handed down-type
	  for (int isq2=1; isq2<3; isq2++)
	    wid += norm(coupSUSYPtr->Rdsq[isq][isq2] 
                 * coupSUSYPtr->rvLQD[ilep][isq2][isq2]);
	}
      }
      else{
	if (iq%2 == 0) {widNow = 0.0; return;}
	// q is down type; ~q is left-handed up-type
	for (int isq2=1; isq2<3; isq2++)
	  wid += norm(coupSUSYPtr->Rusq[isq][isq2] 
               * coupSUSYPtr->rvLQD[ilep][isq2][iq]);
      }
    }
    
    //Case 2: quark + gaugino 
    else if (id1Abs > ksusy && id2Abs < 7) {
      
      int iq = (id2Abs + 1)/2;
      
      // ~q -> ~g + q
      if (id1Abs == 1000021 && idRes%10 == id2Abs) {
	// Removed factor of s2W in denominator: strong process -- no EW
	fac = 2.0 * alpS / (3.0 * pow3(mHat));
	if (idown)
	  wid = kinFac * (norm(coupSUSYPtr->LsddG[isq][iq]) 
              + norm(coupSUSYPtr->RsddG[isq][iq]))
	      - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsddG[isq][iq] 
              * conj(coupSUSYPtr->RsddG[isq][iq]));
	else
	  wid = kinFac * (norm(coupSUSYPtr->LsuuG[isq][iq]) 
              + norm(coupSUSYPtr->RsuuG[isq][iq]))
	      - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsuuG[isq][iq] 
              * conj(coupSUSYPtr->RsuuG[isq][iq]));
      } 
      else 
	for (int i=1; i<6 ; i++){
	  // ~q -> ~chi0 + q
	  if (coupSUSYPtr->idNeut(i)==id1Abs && idRes%2 == id2Abs%2){
	    fac = alpEM *  preFac / (2.0 * (1 - s2W));
	    if (idown)
	      wid = kinFac * (norm(coupSUSYPtr->LsddX[isq][iq][i]) 
                  + norm(coupSUSYPtr->RsddX[isq][iq][i]))
		  - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsddX[isq][iq][i] 
                  * conj(coupSUSYPtr->RsddX[isq][iq][i]));
	    else
	      wid = kinFac * (norm(coupSUSYPtr->LsuuX[isq][iq][i]) 
                  + norm(coupSUSYPtr->RsuuX[isq][iq][i]))
		  - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsuuX[isq][iq][i] 
                  * conj(coupSUSYPtr->RsuuX[isq][iq][i]));
	  }
	  
	  // ~q -> chi- + q
	  else if (i < 3 && coupSUSYPtr->idChar(i)==id1Abs 
            && idRes%2 != id2Abs%2){
	    
	    fac = alpEM *  preFac / (4.0 * (1 - s2W));
	    if (idown)
	      wid = kinFac * (norm(coupSUSYPtr->LsduX[isq][iq][i]) 
                  + norm(coupSUSYPtr->RsduX[isq][iq][i]))
		  - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsduX[isq][iq][i] 
                  * conj(coupSUSYPtr->RsduX[isq][iq][i]));
	    else
	      wid = kinFac * (norm(coupSUSYPtr->LsudX[isq][iq][i]) 
                  + norm(coupSUSYPtr->RsudX[isq][iq][i]))
		  - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsudX[isq][iq][i] 
                  * conj(coupSUSYPtr->RsudX[isq][iq][i]));
	  }
	}
    }
    
    //Case 3: ~q_i -> ~q_j + Z/W
    else if (id1Abs > ksusy && id1Abs%100 < 7 
      && (id2Abs == 23 || id2Abs == 24)){
      
      // factor of lambda^(3/2) = ps^(3/2) ; 
      fac = alpEM * preFac/(16.0 * pow2(particleDataPtr->m0(id2Abs)) 
          * (1.0 - s2W)) * pow2(ps) ;
      
      int isq2 = (id1Abs/ksusy == 2) ? (id1Abs%10+1)/2 + 3: (id1Abs%10+1)/2;
      
      if (id2Abs == 23 && id1Abs%2 == idRes%2){
	if (idown)
	  wid = norm(coupSUSYPtr->LsdsdZ[isq][isq2] 
		     + coupSUSYPtr->RsdsdZ[isq][isq2]);
	else
	  wid = norm(coupSUSYPtr->LsusuZ[isq][isq2] 
		     + coupSUSYPtr->RsusuZ[isq][isq2]);
      }
      else if (id2Abs == 24 && id1Abs%2 != idRes%2){
	if (idown)
	  wid = norm(coupSUSYPtr->LsusdW[isq2][isq]);
	else
	  wid = norm(coupSUSYPtr->LsusdW[isq][isq2]);
      }
    }
    
    // TODO: Case ~q_i -> ~q_j + h/H
    widNow = fac * wid * ps;
    if (DBSUSY) cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs
		  <<" Width: "<<widNow<<endl;
    return;
  }
	
}

//==========================================================================

// The ResonanceGluino class
// Derived class for Gluino resonances

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceGluino::initConstants() {

  // Locally stored properties and couplings.
  alpS  = coupSUSYPtr->alphaS(mHat * mHat);
  return;
}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceGluino::calcPreFac(bool) {
  // Common coupling factors.
  preFac = alpS /( 8.0 * pow(mHat,3));
  return;
}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceGluino::calcWidth(bool) {


  widNow = 0.0;
  if (ps == 0.) return;
  kinFac = (mHat * mHat - mf1 * mf1 + mf2 * mf2);

  if (id1Abs > 1000000 && (id1Abs % 100) < 7 && id2Abs < 7) {

    int isq = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                         : (abs(id1Abs)%10+1)/2;
    bool idown = id2Abs%2;
    int iq = (id2Abs + 1)/2;

    // ~g -> ~q + q
    if (idown){
      widNow = kinFac * (norm(coupSUSYPtr->LsddG[isq][iq]) 
             + norm(coupSUSYPtr->RsddG[isq][iq]))
	     + 4.0 * mHat * mf2 * real(coupSUSYPtr->LsddG[isq][iq] 
             * conj(coupSUSYPtr->RsddG[isq][iq]));
    }
    else{
      widNow = kinFac * (norm(coupSUSYPtr->LsuuG[isq][iq]) 
             + norm(coupSUSYPtr->RsuuG[isq][iq]))
	     + 4.0 * mHat * mf2 * real(coupSUSYPtr->LsuuG[isq][iq] 
             * conj(coupSUSYPtr->RsuuG[isq][iq]));
    }
    widNow = widNow * preFac * ps;
    if (DBSUSY) {
      cout<<"Gluino:: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: ";
      cout<<scientific<<widNow<<endl;
    }
    return;
  }
}

//==========================================================================

//  Class ResonanceNeut
//  Derived class for Neutralino Resonances
//
//--------------------------------------------------------------------------


void ResonanceNeut::initConstants(){
  
  alpEM = coupSUSYPtr->alphaEM(mHat * mHat);
  s2W   = coupSUSYPtr->sin2W;

  // Initialize functions for calculating 3-body widths
  psi.init(particleDataPtr,coupSUSYPtr);
  phi.init(particleDataPtr,coupSUSYPtr);
  upsil.init(particleDataPtr,coupSUSYPtr);
}
 
//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.
void  ResonanceNeut::calcPreFac(bool){

  // Common coupling factors.
  preFac = alpEM / (8.0 * s2W * pow(mHat,3));
  return;
}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.
void  ResonanceNeut::calcWidth(bool){

  widNow = 0.0;

  if (ps ==0.) return;
  else if (mult ==2){
    // Two-body decays
    
    kinFac = mHat * mHat - mf1 * mf1 + mf2 * mf2;
    kinFac2 = pow(mHat,4) + pow(mf1,4) - 2.0 * pow(mf2,4) 
            + pow2(mHat) * pow2(mf2) + pow2(mf1) * pow2(mf2) 
            - 2.0 * pow2(mHat) * pow2(mf1);
    
    // Stable lightest neutralino
    if (idRes == 1000022) return;
    
    double fac = 0.0;
    int iNeut1 = coupSUSYPtr->typeNeut(idRes);
    int iNeut2 = coupSUSYPtr->typeNeut(id1Abs);
    int iChar1 = coupSUSYPtr->typeChar(id1Abs);
    
    if (iNeut2>0 && id2Abs == 23){
      // ~chi0_i -> chi0_j + Z
      fac = kinFac2 * (norm(coupSUSYPtr->OLpp[iNeut1][iNeut2]) 
          + norm(coupSUSYPtr->ORpp[iNeut1][iNeut2]));
      fac -= 12.0 * mHat * mf1 * pow2(mf2) 
           * real(coupSUSYPtr->OLpp[iNeut1][iNeut2] 
	   * conj(coupSUSYPtr->ORpp[iNeut1][iNeut2]));
      fac /= pow2(mf2) * (1.0 - s2W);
    }
    else if (iChar1>0 && id2Abs==24){
      // ~chi0_i -> chi+_j + W- (or c.c.)
      
      fac = kinFac2 * (norm(coupSUSYPtr->OL[iNeut1][iChar1]) 
          + norm(coupSUSYPtr->OR[iNeut1][iChar1]));
      fac -= 12.0 * mHat * mf1 * pow2(mf2) 
           * real(coupSUSYPtr->OL[iNeut1][iChar1] 
	   * conj(coupSUSYPtr->OR[iNeut1][iChar1]));
      fac /= pow2(mf2);
    }
    else if (id1Abs > 1000000 && id1Abs%100 < 7 && id2Abs < 7){
      // ~chi0_k -> ~q + q
      bool idown = (id1Abs%2 == 1);
      int iq = (id2Abs + 1 )/ 2;
      int isq = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                           : (abs(id1Abs)%10+1)/2;
      
      if (idown){
	fac  = kinFac * (norm(coupSUSYPtr->LsddX[isq][iq][iNeut1]) 
             + norm(coupSUSYPtr->RsddX[isq][iq][iNeut1]));
	fac += 4.0 * mHat * mf2 * real(coupSUSYPtr->LsddX[isq][iq][iNeut1] 
	     * conj(coupSUSYPtr->RsddX[isq][iq][iNeut1]));
      }
      else{
	fac = kinFac * (norm(coupSUSYPtr->LsuuX[isq][iq][iNeut1]) 
            + norm(coupSUSYPtr->RsuuX[isq][iq][iNeut1]));
	fac += 4.0 * mHat * mf2 * real(coupSUSYPtr->LsuuX[isq][iq][iNeut1] 
	     * conj(coupSUSYPtr->RsuuX[isq][iq][iNeut1]));
      }
      // Extra multiplicative factor of 3 over sleptons
      fac *= 6.0/(1 - s2W);
    }
    else if (id1Abs > 1000000 && id1Abs%100 > 10 && id1Abs%100 < 17 
      && id2Abs < 17){
      // ~chi0_k -> ~l + l
      bool idown = id2Abs%2;
      int il = (id2Abs - 9)/ 2;
      int isl = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                           : (abs(id1Abs)%10+1)/2;
      
      if (idown){
	fac  = kinFac * (norm(coupSUSYPtr->LsllX[isl][il][iNeut1]) 
             + norm(coupSUSYPtr->RsllX[isl][il][iNeut1]));
	fac += 4.0 * mHat * mf2 * real(coupSUSYPtr->LsllX[isl][il][iNeut1] 
	     * conj(coupSUSYPtr->RsllX[isl][il][iNeut1]));
      }
      else{
	fac = kinFac * (norm(coupSUSYPtr->LsvvX[isl][il][iNeut1]));
      }
      fac *= 2.0/(1 - s2W);
    }
    // TODO: Decays in higgs
    // Final width for 2-body decays
    widNow = fac * preFac * ps ;
    if (DBSUSY) {
      cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: ";
      cout<<scientific<<widNow<<endl;
    }
  }
  else {
    //RPV 3-body decays
    //Case: UDD-type

    if (!coupSUSYPtr->isUDD) return;

    if (id1Abs < 7 && id2Abs < 7 && id3Abs < 7){

      // Check that quarks compatible with UDD structure
      if ((id1Abs+id2Abs+id3Abs)%2 == 1) return; 
      double rvfac,m12min,m12max,integral;
      int idInt;
      
      // Loop over mass eigenstate in actual propagator
      for (int idIntRes = 1; idIntRes <= 6; idIntRes++){ 
	// Right handed field in the UDD coupling
	for (int iSq = 1; iSq <= 3; iSq++){ 
	  double m1, m2, m3, mixfac1(0.), mixfac2(0.), mixfac3(0.);
	  int itemp1,itemp2,itemp3;
	  // up/down-type internal lines
	  for (int itype = 1; itype<=3; itype++){ 
	    //itype = 1: up
	    //itype = 2: down1
	    //itype = 3: down2
	    if (itype ==1 ) idInt = coupSUSYPtr->idSup(idIntRes);
	    else idInt = coupSUSYPtr->idSdown(idIntRes);
	    if (id1Abs%2 == 0){
	      if (itype == 1){
		itemp3 = id1Abs;
		itemp1 = id2Abs;
		itemp2 = id3Abs;
		rvfac = pow2(
                  coupSUSYPtr->rvUDD[iSq][(id2Abs+1)/2][(id3Abs+1)/2]);
	      }else if (itype ==2){
		itemp3 = id2Abs;
		itemp1 = id1Abs;
		itemp2 = id3Abs;
		rvfac = pow2(
                  coupSUSYPtr->rvUDD[(id1Abs+1)/2][iSq][(id3Abs+1)/2]);
	      } else{
		itemp3 = id3Abs;
		itemp1 = id1Abs;
		itemp2 = id2Abs;
		rvfac = pow2(
                  coupSUSYPtr->rvUDD[(id1Abs+1)/2][(id2Abs+1)/2][iSq]);
	      }
	    }else if (id2Abs%2 == 0){
	      if (itype==1){
		itemp3 = id2Abs;
		itemp1 = id1Abs;
		itemp2 = id3Abs;
		rvfac = pow2(
                  coupSUSYPtr->rvUDD[iSq][(id1Abs+1)/2][(id3Abs+1)/2]);
	      }else if (itype ==2){
		itemp3 = id1Abs;
		itemp1 = id2Abs;
		itemp2 = id3Abs;
		rvfac = pow2(
                  coupSUSYPtr->rvUDD[(id2Abs+1)/2][iSq][(id3Abs+1)/2]);
	      } else{
		itemp3 = id3Abs;
		itemp1 = id2Abs;
		itemp2 = id1Abs;
		rvfac = pow2(
                  coupSUSYPtr->rvUDD[(id2Abs+1)/2][(id2Abs+1)/2][iSq]);
	      }
	    }else{
	      if (itype==1){
		itemp3 = id3Abs;
		itemp1 = id1Abs;
		itemp2 = id2Abs;
		rvfac = pow2(
                  coupSUSYPtr->rvUDD[iSq][(id1Abs+1)/2][(id2Abs+1)/2]);
	      }else if (itype ==2){
		itemp3 = id2Abs;
		itemp1 = id1Abs;
		itemp2 = id3Abs;
		rvfac = pow2(
                  coupSUSYPtr->rvUDD[(id3Abs+1)/2][iSq][(id2Abs+1)/2]);
	      } else{
		itemp3 = id3Abs;
		itemp1 = id1Abs;
		itemp2 = id2Abs;
		rvfac = pow2(
                  coupSUSYPtr->rvUDD[(id3Abs+1)/2][(id1Abs+1)/2][iSq]);
	      }
	      
	    }
	    
	    m1 = particleDataPtr->m0(itemp1);
	    m2 = particleDataPtr->m0(itemp2);
	    m3 = particleDataPtr->m0(itemp3);

	    m12min = pow2(m1+m2);
	    m12max = pow2(mHat-m3);

	    // Ignore mode when 2-body decay is possible
	    if (mRes > particleDataPtr->m0(idInt) + particleDataPtr->m0(itemp3))
	      continue;
	    
	    // Single diagram squared terms
	    psi.setInternal(idRes, itemp1, itemp2, itemp3, idInt, 0);
	    // Mixing with R-states
	    if (itype == 1)
	      mixfac1 = norm(coupSUSYPtr->Rusq[idIntRes][iSq+3]); 
	    else
	      mixfac1 = norm(coupSUSYPtr->Rdsq[idIntRes][iSq+3]); 
	    
	    if (abs(rvfac * mixfac1) > 1.0e-8) {
	      integral =  integrateGauss(&psi,m12min,m12max,1.0e-4);
	      widNow += rvfac * mixfac1 * integral;
	    //   if (DBSUSY || idRes == 1000023)
	    // 	cout << scientific << setw(10) <<"Psi: intRes: "<<idInt
	    // 	     <<" integral:"<<integral<<" mixfac:"<<mixfac1
	    // 	     <<" widNow:"<<widNow<<endl;
	    }
	    
	    // Mixing of diagrams with different internal squarks 
            // of same isospin
	    for (int idIntRes2 = 1; idIntRes2 <= 6; idIntRes2++){
	      if (idIntRes2 == idIntRes) continue;
	      int idInt2;
	      if (itype == 1 ){
	        idInt2 = coupSUSYPtr->idSup(idIntRes2);
		mixfac2 = 2.0 * real(coupSUSYPtr->Rusq[idIntRes][iSq+3] 
			* conj(coupSUSYPtr->Rusq[idIntRes2][iSq+3]));
	      } else {
		idInt2 = coupSUSYPtr->idSdown(idIntRes2);
		mixfac2 = 2.0 * real(coupSUSYPtr->Rdsq[idIntRes][iSq+3]
			* conj(coupSUSYPtr->Rdsq[idIntRes2][iSq+3]));
	      }

	      // Ignore mode when 2-body decay is possible
	      if (mRes > particleDataPtr->m0(idInt2) 
                      + particleDataPtr->m0(itemp3)) continue;

	      upsil.setInternal(idRes,itemp1, itemp2,itemp3,idInt,idInt2);
	      if (abs(rvfac * mixfac2) > 0.0) {
		integral =  integrateGauss(&upsil,m12min,m12max,1.0e-4);
		widNow += rvfac * mixfac2 * integral;
		// if (DBSUSY || idRes == 1000023)
		//   cout << scientific << setw(10) <<"Upsilon: intRes: "
		//        <<idInt<<" intRes2:"<<idInt2<<" integral:"<<integral
		//        <<" mixfac:"<<mixfac2<<" widNow:"<<widNow<<endl;
	      }
	    }
	    
	    // Interference between two diagrams with quarks 
            // of different isospin

	    for (int idIntRes2 = 1; idIntRes2 <= 6; idIntRes2++){
	      if (itype != 1 && idIntRes2 == idIntRes) continue;
	      int idInt2;

	      for (int iSq2 = 1; iSq2 <= 3; iSq2++){
		if (itype == 1 ){
		  idInt2 = coupSUSYPtr->idSdown(idIntRes2);
		  mixfac3 = 2.0 * real(coupSUSYPtr->Rusq[idIntRes][iSq+3]  
			  * conj(coupSUSYPtr->Rdsq[idIntRes2][iSq2+3]));
		} else {
		  idInt2 = coupSUSYPtr->idSdown(idIntRes2);
		  mixfac3 = 2.0 * real(coupSUSYPtr->Rdsq[idIntRes][iSq+3]  
			  * conj(coupSUSYPtr->Rdsq[idIntRes2][iSq2+3]));
		}

		if (abs(rvfac * mixfac3) > 0.0) {
		  phi.setInternal(idRes,itemp1, itemp2,itemp3,idInt,idInt2);
		  //integral = 0;
		  //if (idIntRes == 2 && iSq2 ==4)
		  integral =  integrateGauss(&phi,m12min,m12max,1.0e-4);
		  widNow -= rvfac * mixfac2 * integral;
		}
	      }
	    }
	  }
	}
      }
    }
    // Normalisation.  Extra factor of 2 to account for the fact that
    // d_i, d_j will always be ordered in ascending order. N_c! = 6
    widNow *= 12.0 /(pow3(2.0 * M_PI * mHat) * 32.0); 
  }

  return;
}

//==========================================================================

//  Class ResonanceChar
//  Derived class for Neutralino Resonances
//  Decays into higgses/sleptons not yet implemented

//--------------------------------------------------------------------------

void ResonanceChar::initConstants(){

  alpEM = coupSUSYPtr->alphaEM(mHat * mHat);
  s2W   = coupSUSYPtr->sin2W;
  return;
}
 
//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.
void  ResonanceChar::calcPreFac(bool){

  preFac = alpEM / (8.0 * s2W * pow(mHat,3));
  return;
}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.
void  ResonanceChar::calcWidth(bool){

  widNow = 0.0;
  if (ps == 0.) return;

  if (mult ==2){
    double fac = 0.0;
    kinFac = mHat * mHat - mf1 * mf1 + mf2 * mf2;
    kinFac2 = pow(mHat,4) + pow(mf1,4) - 2.0 * pow(mf2,4) 
      + pow2(mHat) * pow2(mf2) + pow2(mf1) 
      * pow2(mf2) - 2.0 * pow2(mHat) * pow2(mf1);
    
    int idChar1 = coupSUSYPtr->typeChar(idRes);
    int idChar2 = coupSUSYPtr->typeChar(id1Abs);
    int idNeut1 = coupSUSYPtr->typeNeut(id1Abs);
    
    if (idChar2>0 && id2Abs == 23){
      // ~chi_i -> chi_j + Z
      fac = kinFac2 * (norm(coupSUSYPtr->OLp[idChar1][idChar2]) 
          + norm(coupSUSYPtr->ORp[idChar1][idChar2]));
      fac -= 12.0 * mHat * mf1 * pow2(mf2) 
           * real(coupSUSYPtr->OLp[idChar1][idChar2] 
	   * conj(coupSUSYPtr->ORp[idChar1][idChar2]));
      fac /= pow2(mf2) * (1.0 - s2W);
    }
    else if (idNeut1>0 && id2Abs==24){
      // ~chi_i -> chi0_j + W- (or c.c.)
      
      fac  = kinFac2 * (norm(coupSUSYPtr->OL[idNeut1][idChar1]) 
           + norm(coupSUSYPtr->OR[idNeut1][idChar1]));
      fac -= 12.0 * mHat * mf1 * pow2(mf2) 
           * real(coupSUSYPtr->OL[idNeut1][idChar1] 
	   * conj(coupSUSYPtr->OR[idNeut1][idChar1]));
      fac /= pow2(mf2);
    }
    else if (id1Abs > 1000000 && id1Abs%100 < 7 && id2Abs < 7){
      // ~chi0_k -> ~q + q
      bool idown = (id1Abs%2 == 1);
      int iq = (id2Abs + 1 )/ 2;
      int isq = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                           : (abs(id1Abs)%10+1)/2;
      
      if (idown){
	fac  = kinFac * (norm(coupSUSYPtr->LsduX[isq][iq][idChar1]) 
	     + norm(coupSUSYPtr->RsduX[isq][iq][idChar1]));
	fac += 4.0 * mHat * mf2  
	     * real(coupSUSYPtr->LsduX[isq][iq][idChar1] 
	     * conj(coupSUSYPtr->RsduX[isq][iq][idChar1]));
      }
      else{
	fac  = kinFac * (norm(coupSUSYPtr->LsudX[isq][iq][idChar1]) 
	     + norm(coupSUSYPtr->RsudX[isq][iq][idChar1]));
	fac += 4.0 * mHat * mf2  
	     * real(coupSUSYPtr->LsudX[isq][iq][idChar1] 
	     * conj(coupSUSYPtr->RsudX[isq][iq][idChar1]));
      }
      fac *= 6.0/(1 - s2W);
    }
    else if (id1Abs > 1000000 && id1Abs%100 > 10 && id1Abs%100 < 17 
      && id2Abs < 17){
      // ~chi+_k -> ~l + l
      bool idown = id2Abs%2;
      int il = (id2Abs - 9)/ 2;
      int isl = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                           : (abs(id1Abs)%10+1)/2;
      
      if (idown){
	fac  = kinFac * (norm(coupSUSYPtr->LslvX[isl][il][idChar1]) 
             + norm(coupSUSYPtr->RslvX[isl][il][idChar1]));
	fac += 4.0 * mHat * mf2 * real(coupSUSYPtr->LslvX[isl][il][idChar1] 
	     * conj(coupSUSYPtr->RslvX[isl][il][idChar1]));
      }
      else{
	fac = kinFac * (norm(coupSUSYPtr->LsvlX[isl][il][idChar1]));
      }
      fac *= 2.0/(1 - s2W);
    }

    // TODO: Decays in higgs
    widNow = fac * preFac * ps ;
    if (DBSUSY) {
      cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: ";
      cout<<scientific<<widNow<<endl;
    }
  }else{
    //TODO: Implement Chargino 3-body decays
  }
  return;
}


//==========================================================================
// The ResonanceSlepton class
// Derived class for Slepton (and sneutrino) resonances

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceSlepton::initConstants() {

  // Locally stored properties and couplings.
  alpEM = coupSUSYPtr->alphaEM(mHat * mHat);
  s2W   = coupSUSYPtr->sin2W;
}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceSlepton::calcPreFac(bool) {

  // Common coupling factors.
  preFac = 1.0 / (s2W * pow(mHat,3));

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceSlepton::calcWidth(bool) {

  // Slepton type -- in u_i/d_i and generation
  int ksusy = 1000000;
  int isl = (abs(idRes)/ksusy == 2) ? (abs(idRes)%10+1)/2 + 3
                                    : (abs(idRes)%10+1)/2;
  int il = (id2Abs-9)/2;
  bool islep = idRes%2;

  // Check that mass is above threshold.
  if (ps == 0.) return;
  else{
    // Two-body decays
    kinFac = (mHat * mHat - mf1 * mf1 - mf2 * mf2);
    
    double fac = 0.0 , wid = 0.0;

    //Case 1: RPV: To be implemented
    //Case 2: slepton + gaugino 

    if (id1Abs > ksusy && id2Abs > 10 && id2Abs < 17) {
      for (int i=1; i<6 ; i++){
	// ~ell/~nu -> ~chi0 + ell/nu
	if (coupSUSYPtr->idNeut(i)==id1Abs && idRes%2 == id2Abs%2){
	  fac = alpEM *  preFac / (2.0 * (1 - s2W));
	  if (islep)
	    wid = kinFac * (norm(coupSUSYPtr->LsllX[isl][il][i]) 
                + norm(coupSUSYPtr->RsllX[isl][il][i]))
	        - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsllX[isl][il][i] 
		* conj(coupSUSYPtr->RsllX[isl][il][i]));
	  else
	    wid = kinFac * (norm(coupSUSYPtr->LsvvX[isl][il][i]) 
                + norm(coupSUSYPtr->RsvvX[isl][il][i]))
	        - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsvvX[isl][il][i] 
                * conj(coupSUSYPtr->RsvvX[isl][il][i]));
	}
	
	// ~ell/~nu -> ~chi- + nu/ell
	else if (i < 3 && coupSUSYPtr->idChar(i)==id1Abs 
          && idRes%2 != id2Abs%2){
	  
	  fac = alpEM *  preFac / (4.0 * (1 - s2W));
	  if (islep)
	    wid = kinFac * (norm(coupSUSYPtr->LslvX[isl][il][i]) 
                + norm(coupSUSYPtr->RslvX[isl][il][i]))
	        - 4.0 * mHat * mf2 * real(coupSUSYPtr->LslvX[isl][il][i] 
                * conj(coupSUSYPtr->RslvX[isl][il][i]));
	  else
	    wid = kinFac * (norm(coupSUSYPtr->LslvX[isl][il][i]) 
                + norm(coupSUSYPtr->RslvX[isl][il][i]))
	        - 4.0 * mHat * mf2 * real(coupSUSYPtr->LslvX[isl][il][i] 
                * conj(coupSUSYPtr->RslvX[isl][il][i]));
	}
      }
    }
    
    //Case 3: ~l_i -> ~l_j + Z/W
    else if (id1Abs > ksusy+10 && id1Abs%100 < 17 
      && (id2Abs == 23 || id2Abs == 24)){
      
      // factor of lambda^(3/2) = ps^3; 
      fac = alpEM * preFac/(16.0 * pow2(mf2) * (1.0 - s2W)) * pow2(ps) ;
      
      int isl2 = (id1Abs/ksusy == 2) ? (id1Abs%10+1)/2 + 3: (id1Abs%10+1)/2;
      
      if (id2Abs == 23 && id1Abs%2 == idRes%2){
	if (islep)
	  wid = norm(coupSUSYPtr->LslslZ[isl][isl2] 
              + coupSUSYPtr->RslslZ[isl][isl2]);
	else
	  wid = norm(coupSUSYPtr->LsvsvZ[isl][isl2] 
              + coupSUSYPtr->RsvsvZ[isl][isl2]);
      }
      else if (id2Abs == 24 && id1Abs%2 != idRes%2){
	if (islep)
	  wid = norm(coupSUSYPtr->LslsvW[isl2][isl]);
	else
	  wid = norm(coupSUSYPtr->LslsvW[isl][isl2]);
      }
    }
    
    // TODO: Case ~l_i -> ~l_j + h/H
    
    
    widNow = fac * wid * ps;
    if (DBSUSY) cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs
		  <<" Width: "<<widNow<<endl;
    return;
  }
	
}

//==========================================================================

// Gaussian Integrator for 3-body decay widths

double SUSYResonanceWidths::integrateGauss(WidthFunction* widthFn, 
  double xlo, double xhi, double tol) {

  //8-point unweighted
  static double x8[4]={0.96028985649753623,
		       0.79666647741362674,
		       0.52553240991632899, 
		       0.18343464249564980};
  static double w8[4]={0.10122853629037626,
		       0.22238103445337447,
		       0.31370664587788729,
		       0.36268378337836198};
  //16-point unweighted
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
  //boundary checks
  if (xlo == xhi) {
    cerr<<"xlo = xhi"<<endl;
    return 0.0;
  }
  if (xlo > xhi) {
    cerr<<" (integrateGauss:) -> xhi < xlo"<<endl;
    return 0.0;
  }
  //initialize
  double sum = 0.0;
  double c = 0.001/abs(xhi-xlo);
  double zlo = xlo;
  double zhi = xhi;
    
  bool nextbin = true;
  
  while ( nextbin ) {
    
    double zmi = 0.5*(zhi+zlo); // midpoint
    double zmr = 0.5*(zhi-zlo); // midpoint, relative to zlo
    
    //calculate 8-point and 16-point quadratures
    double s8=0.0;
    for (int i=0;i<4;i++) {
      double dz = zmr * x8[i];
      s8 += w8[i]*(widthFn->function(zmi+dz) + widthFn->function(zmi-dz));
    }
    s8 *= zmr;
    double s16=0.0;
    for (int i=0;i<8;i++) {
      double dz = zmr * x16[i];
      s16 += w16[i]*(widthFn->function(zmi+dz) + widthFn->function(zmi-dz)); 
    }
    s16 *= zmr;
    if (abs(s16-s8) < tol*(1+abs(s16))) { 
      //precision in this bin OK, add to cumulative and go to next
      nextbin=true;
      sum += s16;
      //next bin: LO = end of current, HI = end of integration region.
      zlo=zhi;
      zhi=xhi;
      if ( zlo == zhi ) nextbin = false;
    } else {
      //precision in this bin not OK, subdivide.
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

//======================================================================

} //end namespace Pythia8


