// SigmaSUSY.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// supersymmetry simulation classes. 

#include "SigmaSUSY.h"

namespace Pythia8 {
  
//==========================================================================

// Sigma2qqbar2chi0chi0 
// Cross section for gaugino pair production: neutralino pair

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qqbar2chi0chi0::initProc() {

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  // Construct name of process. 
  nameSave = "q qbar' -> " + particleDataPtr->name(id3) + " " 
    + particleDataPtr->name(id4);

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(id3, id4);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qqbar2chi0chi0::sigmaKin() {

  // Common flavour-independent factor.
  sigma0 = M_PI /3.0/ sH2 / pow2(coupSUSYPtr->sin2W) * pow2(alpEM) 
    * openFracPair; 

  // Auxiliary factors for use below
  ui       = uH - s3;
  uj       = uH - s4;
  ti       = tH - s3;
  tj       = tH - s4;
  double sV= sH - pow2(coupSUSYPtr->mZpole);
  double d = pow2(sV) + pow2(coupSUSYPtr->mZpole * coupSUSYPtr->wZpole);
  propZ    = complex( sV / d, coupSUSYPtr->mZpole * coupSUSYPtr->wZpole / d);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qqbar2chi0chi0::sigmaHat() {

  // Only allow quark-antiquark incoming states
  if (id1*id2 >= 0) {
    return 0.0;    
  }
  
  // Only allow incoming states with sum(charge) = 0
  if ((id1+id2) % 2 != 0) {
    return 0.0;    
  }

  if(id1<0) swapTU = true;
  
  // Shorthands
  int idAbs1    = abs(id1);  
  int idAbs2    = abs(id2);

  // Flavour-dependent kinematics-dependent couplings.
  complex QuLL(0.0),QtLL(0.0),QuRR(0.0),QtRR(0.0);
  complex QuLR(0.0),QtLR(0.0),QuRL(0.0),QtRL(0.0);

  // s-channel Z couplings
  if (idAbs1 == idAbs2) {
    QuLL = coupSUSYPtr->LqqZ[idAbs1] * coupSUSYPtr->OLpp[id3chi][id4chi] 
         * propZ / 2.0;
    QtLL = coupSUSYPtr->LqqZ[idAbs1] * coupSUSYPtr->ORpp[id3chi][id4chi] 
         * propZ / 2.0;
    QuRR = coupSUSYPtr->RqqZ[idAbs1] * coupSUSYPtr->ORpp[id3chi][id4chi] 
         * propZ / 2.0;
    QtRR = coupSUSYPtr->RqqZ[idAbs1] * coupSUSYPtr->OLpp[id3chi][id4chi] 
         * propZ / 2.0;
  }

  // Flavour indices
  int ifl1 = (idAbs1+1) / 2;
  int ifl2 = (idAbs2+1) / 2;

  // Add t-channel squark flavour sums to QmXY couplings
  for (int ksq=1; ksq<=6; ksq++) {    

    // squark id and squark-subtracted u and t
    int idsq=((ksq+2)/3)*1000000 + 2*((ksq-1) % 3) + (idAbs1+1) % 2 + 1;
    double msq2    = pow(particleDataPtr->m0(idsq),2);
    double usq     = uH - msq2;
    double tsq     = tH - msq2;
    
    // Couplings
    complex Lsqq1X3 = coupSUSYPtr->LsuuX[ksq][ifl1][id3chi];
    complex Lsqq1X4 = coupSUSYPtr->LsuuX[ksq][ifl1][id4chi];
    complex Lsqq2X3 = coupSUSYPtr->LsuuX[ksq][ifl2][id3chi];
    complex Lsqq2X4 = coupSUSYPtr->LsuuX[ksq][ifl2][id4chi];
    complex Rsqq1X3 = coupSUSYPtr->RsuuX[ksq][ifl1][id3chi];
    complex Rsqq1X4 = coupSUSYPtr->RsuuX[ksq][ifl1][id4chi];
    complex Rsqq2X3 = coupSUSYPtr->RsuuX[ksq][ifl2][id3chi];
    complex Rsqq2X4 = coupSUSYPtr->RsuuX[ksq][ifl2][id4chi];
    if (idAbs1 % 2 != 0) {
      Lsqq1X3 = coupSUSYPtr->LsddX[ksq][ifl1][id3chi];
      Lsqq1X4 = coupSUSYPtr->LsddX[ksq][ifl1][id4chi];
      Lsqq2X3 = coupSUSYPtr->LsddX[ksq][ifl2][id3chi];
      Lsqq2X4 = coupSUSYPtr->LsddX[ksq][ifl2][id4chi];
      Rsqq1X3 = coupSUSYPtr->RsddX[ksq][ifl1][id3chi];
      Rsqq1X4 = coupSUSYPtr->RsddX[ksq][ifl1][id4chi];
      Rsqq2X3 = coupSUSYPtr->RsddX[ksq][ifl2][id3chi];
      Rsqq2X4 = coupSUSYPtr->RsddX[ksq][ifl2][id4chi];      
    }

    // QuXY
    QuLL += conj(Lsqq1X4)*Lsqq2X3/usq;
    QuRR += conj(Rsqq1X4)*Rsqq2X3/usq;
    QuLR += conj(Lsqq1X4)*Rsqq2X3/usq;
    QuRL += conj(Rsqq1X4)*Lsqq2X3/usq;
    
    
    // QtXY
    QtLL -= conj(Lsqq1X3)*Lsqq2X4/tsq;
    QtRR -= conj(Rsqq1X3)*Rsqq2X4/tsq;
    QtLR += conj(Lsqq1X3)*Rsqq2X4/tsq;
    QtRL += conj(Rsqq1X3)*Lsqq2X4/tsq;

  }

  // Overall factor multiplying each coupling; multiplied at the end as fac^2
  double fac = (1.0-coupSUSYPtr->sin2W);
  if(abs(id3)==abs(id4)) fac *= sqrt(2.); // for identical final particles

  // Compute matrix element weight
  double weight = 0;
  double facLR = uH*tH - s3*s4;
  double facMS = m3*m4*sH;

  // Average over separate helicity contributions
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += norm(QuLL) * ui * uj + norm(QtLL) * ti * tj
    + 2 * real(conj(QuLL) * QtLL) * facMS;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += norm(QtRR) * ti * tj + norm(QuRR) * ui * uj  
    + 2 * real(conj(QuRR) * QtRR) * facMS;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += norm(QuRL) * ui * uj + norm(QtRL) * ti * tj
    + real(conj(QuRL) * QtRL) * facLR;
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += norm(QuLR) * ui * uj + norm(QtLR) * ti * tj
    + real(conj(QuLR) * QtLR) * facLR;

  // Cross section, including colour factor.
  double sigma = sigma0 * weight / pow2(fac);

  // Answer.
  return sigma;    

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2chi0chi0::setIdColAcol() {

  // Set flavours.
  setId( id1, id2, id3, id4);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2charchi0
// Cross section for gaugino pair production: neutralino-chargino

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qqbar2charchi0::sigmaKin() {

  // Common flavour-independent factor.
  
  sigma0 = M_PI / sH2 / 3.0 / pow2(coupSUSYPtr->sin2W) * pow2(alpEM) ; 
  sigma0 /= 2.0 * (1 - coupSUSYPtr->sin2W) ; 

  // Auxiliary factors for use below
  ui        = uH - s3;
  uj        = uH - s4;
  ti        = tH - s3;
  tj        = tH - s4;
  double sW = sH - pow2(coupSUSYPtr->mWpole);
  double d  = pow2(sW) + pow2(coupSUSYPtr->mWpole * coupSUSYPtr->wWpole);
  propW     = complex( sW / d, coupSUSYPtr->mWpole * coupSUSYPtr->wWpole / d);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qqbar2charchi0::sigmaHat() {

  // Only allow particle-antiparticle incoming states
  if (id1*id2 >= 0) {
    return 0.0;    
  }
  
  // Only allow incoming states with sum(charge) = final state
  if (abs(id1) % 2 == abs(id2) % 2) return 0.0;
  int isPos  = (id3chi > 0 ? 1 : 0);
  if (id1 < 0 && id1 > -10 && abs(id1) % 2 == 1-isPos ) return 0.0;
  else if (id1 > 0 && id1 < 10 && abs(id1) % 2 == isPos ) return 0.0;

  // Flavour-dependent kinematics-dependent couplings.
  int idAbs1  = abs(id1);  
  int iChar = abs(id3chi);
  int iNeut = abs(id4chi); 
  
  complex QuLL(0.0),QtLL(0.0),QuRR(0.0),QtRR(0.0);
  complex QuLR(0.0),QtLR(0.0),QuRL(0.0),QtRL(0.0);
  
  // Calculate everything from udbar -> ~chi+ ~chi0 template process
  
  // u dbar , ubar d : do nothing
  // dbar u , d ubar : swap 1<->2 and t<->u

  int iGu = abs(id1)/2;
  int iGd = (abs(id2)+1)/2;

  if (idAbs1 % 2 != 0) {
    swapTU = true;
    iGu = abs(id2)/2;
    iGd = (abs(id1)+1)/2;
  }

  // s-channel W contribution
  QuLL = conj(coupSUSYPtr->LudW[iGu][iGd]) 
    * conj(coupSUSYPtr->OL[iNeut][iChar])
    * propW / sqrt(2.0);
  QtLL = conj(coupSUSYPtr->LudW[iGu][iGd]) 
    * conj(coupSUSYPtr->OR[iNeut][iChar])
    * propW / sqrt(2.0);

  // Add t-channel squark flavour sums to QmXY couplings
  for (int jsq=1; jsq<=6; jsq++) {    

    int idsu=((jsq+2)/3)*1000000 + 2*((jsq-1) % 3) + 2;
    int idsd=((jsq+2)/3)*1000000 + 2*((jsq-1) % 3) + 1;
    double msd2 = pow(particleDataPtr->m0(idsd),2);
    double msu2 = pow(particleDataPtr->m0(idsu),2);
    double tsq  = tH - msd2;
    double usq  = uH - msu2;

    QuLL += conj(coupSUSYPtr->LsuuX[jsq][iGu][iNeut])
      * conj(coupSUSYPtr->LsudX[jsq][iGd][iChar])/usq;
    QuLR += conj(coupSUSYPtr->LsuuX[jsq][iGu][iNeut])
      * conj(coupSUSYPtr->RsudX[jsq][iGd][iChar])/usq;
    QuRR += conj(coupSUSYPtr->RsuuX[jsq][iGu][iNeut])
      * conj(coupSUSYPtr->RsudX[jsq][iGd][iChar])/usq;
    QuRL += conj(coupSUSYPtr->RsuuX[jsq][iGu][iNeut])
      * conj(coupSUSYPtr->LsudX[jsq][iGd][iChar])/usq;

    QtLL -= conj(coupSUSYPtr->LsduX[jsq][iGu][iChar])
      * coupSUSYPtr->LsddX[jsq][iGd][iNeut]/tsq;
    QtRR -= conj(coupSUSYPtr->RsduX[jsq][iGu][iChar])
      * coupSUSYPtr->RsddX[jsq][iGd][iNeut]/tsq;
    QtLR += conj(coupSUSYPtr->LsduX[jsq][iGu][iChar])
      * coupSUSYPtr->RsddX[jsq][iGd][iNeut]/tsq;
    QtRL += conj(coupSUSYPtr->RsduX[jsq][iGu][iChar])
      * coupSUSYPtr->LsddX[jsq][iGd][iNeut]/tsq;
  }

  // Compute matrix element weight
  double weight = 0;

  // Average over separate helicity contributions
  // (if swapped, swap ha, hb if computing polarized cross sections)
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += norm(QuLL) * ui * uj + norm(QtLL) * ti * tj
    + 2 * real(conj(QuLL) * QtLL) * m3 * m4 * sH;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += norm(QtRR) * ti * tj + norm(QuRR) * ui * uj  
    + 2 * real(conj(QuRR) * QtRR) * m3 * m4 * sH;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += norm(QuRL) * ui * uj + norm(QtRL) * ti * tj
    + real(conj(QuRL) * QtRL) * (uH * tH - s3 * s4);
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += norm(QuLR) * ui * uj + norm(QtLR) * ti * tj
    + real(conj(QuLR) * QtLR) * (uH * tH - s3 * s4);

  // Cross section, including colour factor.
  double sigma = sigma0 * weight;

  // Answer.
  return sigma;    

}

//==========================================================================

// Sigma2qqbar2charchar
// Cross section for gaugino pair production: chargino-chargino

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qqbar2charchar::sigmaKin() {

  // Common flavour-independent factor.
  sigma0 = M_PI / 3.0 / sH2 / pow2(coupSUSYPtr->sin2W) * pow2(alpEM) ; 

  // Auxiliary factors for use below
  ui       = uH - s3;
  uj       = uH - s4;
  ti       = tH - s3;
  tj       = tH - s4;
  double sV= sH - pow2(coupSUSYPtr->mZpole);
  double d = pow2(sV) + pow2(coupSUSYPtr->mZpole * coupSUSYPtr->wZpole);
  propZ    = complex( sV / d, coupSUSYPtr->mZpole * coupSUSYPtr->wZpole / d);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qqbar2charchar::sigmaHat() { 

  // Only allow quark-antiquark incoming states
  if (id1*id2 >= 0) {
    return 0.0;
  }
  
  // Only allow incoming states with sum(charge) = 0
  if ((id1+id2) % 2 != 0) {
    return 0.0;    
  }
  
  //if (id1 > 0 || id1==-1 || id1==-3 || id1==-5) return 0.0;
  //if (id1 < 0 || id1==1 || id1==3 || id1==5) return 0.0;
  
  swapTU = (id1 < 0 ? true : false);

  // Flavour-dependent kinematics-dependent couplings.
  int idAbs1    = abs(id1);  
  int idAbs2    = abs(id2);  
  int i3        = abs(id3chi);
  int i4        = abs(id4chi);
  
  // Flavour-dependent kinematics-dependent couplings.
  complex QuLL(0.0),QtLL(0.0),QuRR(0.0),QtRR(0.0);
  complex QuLR(0.0),QtLR(0.0),QuRL(0.0),QtRL(0.0);

  // Add Z/gamma* for same-flavour in-quarks
  if (idAbs1 == idAbs2) {
    
    QuLL = -coupSUSYPtr->LqqZ[idAbs1]*conj(coupSUSYPtr->ORp[i3][i4]);
    QtLL = -coupSUSYPtr->LqqZ[idAbs1]*conj(coupSUSYPtr->OLp[i3][i4]);
    QuRR = -coupSUSYPtr->RqqZ[idAbs1]*conj(coupSUSYPtr->OLp[i3][i4]);
    QtRR = -coupSUSYPtr->RqqZ[idAbs1]*conj(coupSUSYPtr->ORp[i3][i4]);

    QuLL *= propZ / 2.0 / (1.0-coupSUSYPtr->sin2W);
    QtLL *= propZ / 2.0 / (1.0-coupSUSYPtr->sin2W);
    QuRR *= propZ / 2.0 / (1.0-coupSUSYPtr->sin2W);
    QtRR *= propZ / 2.0 / (1.0-coupSUSYPtr->sin2W);  
  
    // s-channel gamma* (only for same-type charginos)
    if (i3 == i4) {
    
      // Charge of in-particles
      double q = 2.0/3.0;
      if (idAbs1 % 2 == 1) q = -1.0/3.0;      
      QuLL += q * coupSUSYPtr->sin2W / sH;
      QuRR += q * coupSUSYPtr->sin2W / sH;
      QtLL += q * coupSUSYPtr->sin2W / sH;
      QtRR += q * coupSUSYPtr->sin2W / sH;
    
    }
  }

  int iG1    = (abs(id1)+1)/2;
  int iG2    = (abs(id2)+1)/2;

  // Add t- or u-channel squark flavour sums to QmXY couplings
  for (int ksq=1; ksq<=6; ksq++) {    
    
    if(id1 % 2 == 0) { 

      // u-channel diagrams only
      // up-type incoming; u-channel ~d 
      
      int idsd    = ((ksq+2)/3)*1000000 + 2*((ksq-1) % 3) + 1;
      double msq  = particleDataPtr->m0(idsd);
      double ufac = 2.0 * (uH - pow2(msq));

      //u-ubar -> chi-chi+
      QuLL += coupSUSYPtr->LsduX[ksq][iG2][i3] 
            * conj(coupSUSYPtr->LsduX[ksq][iG1][i4]) / ufac;
      QuRR += coupSUSYPtr->RsduX[ksq][iG2][i3] 
            * conj(coupSUSYPtr->RsduX[ksq][iG1][i4]) / ufac;
      QuLR += coupSUSYPtr->RsduX[ksq][iG2][i3] 
            * conj(coupSUSYPtr->LsduX[ksq][iG1][i4]) / ufac;
      QuRL += coupSUSYPtr->LsduX[ksq][iG2][i3] 
            * conj(coupSUSYPtr->RsduX[ksq][iG1][i4]) / ufac;


    }else{
      // t-channel diagrams only;
      // down-type incoming; t-channel ~u 

      int idsu    = ((ksq+2)/3)*1000000 + 2*((ksq-1) % 3) + 2;
      double msq  = particleDataPtr->m0(idsu);
      double tfac = 2.0 * (tH - pow2(msq));

      //d-dbar -> chi-chi+
      QtLL -= coupSUSYPtr->LsudX[ksq][iG1][i3] 
            * conj(coupSUSYPtr->LsudX[ksq][iG2][i4]) / tfac;
      QtRR -= coupSUSYPtr->RsudX[ksq][iG1][i3] 
            * conj(coupSUSYPtr->RsudX[ksq][iG2][i4]) / tfac;
      QtLR += coupSUSYPtr->LsudX[ksq][iG1][i3] 
            * conj(coupSUSYPtr->RsudX[ksq][iG2][i4]) / tfac;
      QtRL += coupSUSYPtr->RsudX[ksq][iG1][i3] 
            * conj(coupSUSYPtr->LsudX[ksq][iG2][i4]) / tfac;

    }
  }
   // Compute matrix element weight
   double weight = 0;

  // Average over separate helicity contributions
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += norm(QuLL) * ui * uj + norm(QtLL) * ti * tj
    + 2 * real(conj(QuLL) * QtLL) * m3 * m4 * sH;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += norm(QtRR) * ti * tj + norm(QuRR) * ui * uj  
    + 2 * real(conj(QuRR) * QtRR) * m3 * m4 * sH;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += norm(QuRL) * ui * uj + norm(QtRL) * ti * tj
    + real(conj(QuRL) * QtRL) * (uH * tH - s3 * s4);
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += norm(QuLR) * ui * uj + norm(QtLR) * ti * tj
    + real(conj(QuLR) * QtLR) * (uH * tH - s3 * s4);

  // Cross section, including colour factor.
  double sigma = sigma0 * weight;

  // Answer.
  return sigma;    

}

//==========================================================================

// Sigma2qgchi0squark 
// Cross section for gaugino-squark production: neutralino-squark

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qg2chi0squark::initProc() {

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  // Construct name of process. 
  if (id4 % 2 == 0) {
    nameSave = "q g -> " + particleDataPtr->name(id3) + " " 
      + particleDataPtr->name(id4) + " + c.c. (q=u,c)";
  } 
  else {
    nameSave = "q g -> " + particleDataPtr->name(id3) + " " 
      + particleDataPtr->name(id4) + " + c.c. (q=d,s,b)";
  }

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(id3, id4);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qg2chi0squark::sigmaKin() {

  // Common flavour-independent factor.
  // tmp: alphaS = 0.1 for counter-checks
  sigma0 = M_PI / sH2 / coupSUSYPtr->sin2W * alpEM * alpS
    * openFracPair;

  // Auxiliary factors for use below
  ui       = uH - s3;
  uj       = uH - s4;
  ti       = tH - s3;
  tj       = tH - s4;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qg2chi0squark::sigmaHat() {

  // Antiquark -> antisquark
  int idq = id1;
  if (id1 == 21 || id1 == 22) idq = id2;
  if (idq < 0) {
    id4 = -abs(id4);
  } else {
    id4 = abs(id4);
  }

  // tmp: only allow incoming quarks on side 1
  //  if (id1 < 0 || id1 == 21) return 0.0;

  // Generation index
  int iGq = (abs(idq)+1)/2;

  // Only accept u(bar) -> ~u(bar) and d(bar) -> ~d(bar)
  if (particleDataPtr->chargeType(idq) != particleDataPtr->chargeType(id4))
    return 0.0;
  
  // Couplings
  complex LsqqX, RsqqX;
  if (idq % 2 == 0) {
    LsqqX = coupSUSYPtr->LsuuX[id4sq][iGq][id3chi];
    RsqqX = coupSUSYPtr->RsuuX[id4sq][iGq][id3chi];
  }
  else { 
    LsqqX = coupSUSYPtr->LsddX[id4sq][iGq][id3chi];
    RsqqX = coupSUSYPtr->RsddX[id4sq][iGq][id3chi];
  }  

  // Prefactors : swap u and t if gq instead of qg
  double fac1, fac2;
  if (idq == id1) {
    fac1 = -ui/sH + 2.0 * ( uH*tH - s4*s3 )/sH/tj;
    fac2 = ti/tj * ( (tH + s4)/tj + (ti - uj)/sH );
  } else {
    fac1 = -ti/sH + 2.0 * ( uH*tH - s4*s3 )/sH/uj;
    fac2 = ui/uj * ( (uH + s4)/uj + (ui - tj)/sH );
  }

  // Compute matrix element weight
  double weight = 0.0;

  // Average over separate helicity contributions
  // (for qbar g : ha -> -ha )
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += fac2 * norm(LsqqX) / 2.0;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += fac2 * norm(RsqqX) / 2.0;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += fac2 * norm(RsqqX) / 2.0 + fac1 * norm(RsqqX);
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += fac2 * norm(LsqqX) / 2.0 + fac1 * norm(LsqqX);

  double sigma = sigma0 * weight;
  if (abs(idq) < 9) sigma /= 3.;

  // Answer.
  return sigma;    

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2chi0squark::setIdColAcol() {

  // Set flavours.
  setId( id1, id2, id3, (id1*id2 > 0 ? abs(id4) : -abs(id4)));

  // Colour flow topology. Swap if first is gluon, or when antiquark.
  if (id1 != 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else setColAcol( 1, 2, 2, 0, 0, 0, 1, 0);
  if (id1*id2 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qg2charsquark
// Cross section for gaugino-squark production: chargino-squark

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qg2charsquark::initProc() {

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  // Construct name of process. 
  if (id4 % 2 == 0) {
    nameSave = "q g -> " + particleDataPtr->name(id3) + " " 
      + particleDataPtr->name(id4) + " + c.c. (q=d,s,b)";
  } 
  else {
    nameSave = "q g -> " + particleDataPtr->name(id3) + " " 
      + particleDataPtr->name(id4) + " + c.c. (q=u,c)";
  }

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(id3Sav, id4Sav);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qg2charsquark::sigmaHat() {

  // Antiquark -> antisquark
  int idq = (id1 == 21) ? id2 : id1;
  if (idq > 0) {
    id3 = id3Sav;
    id4 = id4Sav;
  } else {
    id3 = -id3Sav;
    id4 = -id4Sav;
  }

  // Only accept u(bar) -> ~d(bar) and d(bar) -> ~u(bar)
  if (particleDataPtr->chargeType(idq) == particleDataPtr->chargeType(id4))
    return 0.0;
  
  // Generation index
  int iGq = (abs(idq)+1)/2;

  // Couplings
  complex LsqqX, RsqqX;
  if (idq % 2 == 0) {
    LsqqX = coupSUSYPtr->LsduX[id4sq][iGq][id3chi];
    RsqqX = coupSUSYPtr->RsduX[id4sq][iGq][id3chi];
  }
  else { 
    LsqqX = coupSUSYPtr->LsduX[id4sq][iGq][id3chi];
    RsqqX = coupSUSYPtr->RsduX[id4sq][iGq][id3chi];
  }  

  // Prefactors : swap u and t if gq instead of qg
  double fac1, fac2;
  if (idq == id1) {
    fac1 = -ui/sH + 2.0 * ( uH*tH - s4*s3 )/sH/tj;
    fac2 = ti/tj * ( (tH + s4)/tj + (ti - uj)/sH );
  } else {
    fac1 = -ti/sH + 2.0 * ( uH*tH - s4*s3 )/sH/uj;
    fac2 = ui/uj * ( (uH + s4)/uj + (ui - tj)/sH );
  }

  // Compute matrix element weight
  double weight = 0.0;

  // Average over separate helicity contributions
  // (a, b refers to qg configuration)
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += fac2 * norm(LsqqX) / 2.0;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += fac2 * norm(RsqqX) / 2.0;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += fac2 * norm(RsqqX) / 2.0 + fac1 * norm(RsqqX);
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += fac2 * norm(LsqqX) / 2.0 + fac1 * norm(LsqqX);

  double sigma = sigma0 * weight;
  if (abs(idq) < 9) sigma /= 3.;

  // Answer.
  return sigma * openFracPair;    

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2charsquark::setIdColAcol() {

  // Set flavours.
  if (id1 > 0 && id2 > 0) {
    setId( id1, id2, id3Sav, id4Sav);
  } else {
    setId( id1, id2,-id3Sav,-id4Sav);
  }

  // Colour flow topology. Swap if first is gluon, or when antiquark.
  if (id1 != 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else setColAcol( 1, 2, 2, 0, 0, 0, 1, 0);
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qq2squarksquark
// Cross section for squark-squark production

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qq2squarksquark::initProc() {

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  // Extract mass-ordering indices
  iGen3 = 3*(abs(id3Sav)/2000000) + (abs(id3Sav)%10+1)/2;
  iGen4 = 3*(abs(id4Sav)/2000000) + (abs(id4Sav)%10+1)/2;

  // Is this a ~u_i ~d_j fial state or ~d_i ~d_j, ~u_i ~u_j
  if (abs(id3Sav) % 2 == abs(id4Sav) % 2) isUD = false;
  else isUD = true;

  // Derive name
  nameSave = "q q' -> "+particleDataPtr->name(abs(id3Sav))+" "
    +particleDataPtr->name(abs(id4Sav))+" + c.c.";

  // Count 5 neutralinos in NMSSM
  nNeut = (coupSUSYPtr->isNMSSM ? 5 : 4);

  // Store mass squares of all possible internal propagator lines
  m2Glu     = pow2(particleDataPtr->m0(1000021));
  m2Neut.resize(nNeut+1);
  for (int iNeut=1;iNeut<=nNeut;iNeut++) {
    m2Neut[iNeut] = pow2(particleDataPtr->m0(coupSUSYPtr->idNeut(iNeut)));
  }
  m2Char.resize(3);
  m2Char[1] = pow2(particleDataPtr->m0(coupSUSYPtr->idChar(1)));
  m2Char[2] = pow2(particleDataPtr->m0(coupSUSYPtr->idChar(2)));
  
  // Set sizes of some arrays to be used below
  tNeut.resize(nNeut+1);
  uNeut.resize(nNeut+1);
  tChar.resize(3);
  uChar.resize(3);

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(id3Sav, id4Sav);
}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qq2squarksquark::sigmaKin() {

  // Weak mixing
  double xW = coupSUSYPtr->sin2W;

  // pi/sH2
  double comFacHat = M_PI/sH2 * openFracPair;

  // Channel-dependent but flavor-independent pre-factors
  sigmaNeut     = comFacHat * pow2(alpEM) / pow2(xW) / pow2(1-xW);
  sigmaGlu      = comFacHat * 2.0 * pow2(alpS) / 9.0;
  if (isUD) {
    sigmaChar     = comFacHat * pow2(alpEM) / 4.0 / pow2(xW);
    sigmaCharNeut = comFacHat * pow2(alpEM) / 3.0 / pow2(xW) / (1-xW); 
    sigmaCharGlu  = comFacHat * 4.0 * alpEM * alpS / 9.0 / xW;
    sigmaNeutGlu  = 0.0;
  } else {
    sigmaChar     = 0.0;
    sigmaCharNeut = 0.0;
    sigmaCharGlu  = 0.0;
    sigmaNeutGlu  = comFacHat * 8.0 * alpEM * alpS / 9.0 / xW/(1-xW);
  }

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qq2squarksquark::sigmaHat() {

  // In-pair must be same-sign
  if (id1 * id2 < 0) return 0.0;
  
  // Check correct charge sum
  if (isUD && abs(id1) %2 == abs(id2) % 2) return 0.0;
  if (!isUD && abs(id1) % 2 != abs(id2) % 2) return 0.0;
  if (!isUD && abs(id1) % 2 != abs(id3Sav) % 2) return 0.0;

  // Coded sigma is for ud -> ~q~q'. Swap t and u for du -> ~q~q'.
  swapTU = (isUD && abs(id1) % 2 == 0); 
  int    idIn1A = (swapTU) ? abs(id2) : abs(id1);
  int    idIn2A = (swapTU) ? abs(id1) : abs(id2);

  // Auxiliary factors for use below
  tGlu     = tH - m2Glu;
  uGlu     = uH - m2Glu;
  for (int i=1; i<= nNeut; i++) {
    tNeut[i] = tH - m2Neut[i];
    uNeut[i] = uH - m2Neut[i];
    if (isUD && i <= 2) {
      tChar[i] = tH - m2Char[i];
      uChar[i] = uH - m2Char[i];
    }
  }

  // Generation indices of incoming particles
  int iGen1 = (abs(idIn1A)+1)/2;
  int iGen2 = (abs(idIn2A)+1)/2;

  // Initial values for pieces used for color-flow selection below
  sumCt = 0.0;
  sumCu = 0.0;
  sumNt = 0.0;
  sumNu = 0.0;
  sumGt = 0.0;
  sumGu = 0.0;
  sumInterference = 0.0;

  // Common factor for LR and RL contributions
  double facTU =  uH*tH-s3*s4;

  // Case A) Opposite-isospin: qq' -> ~d~u 
  if ( isUD ) {

    // t-channel Charginos
    for (int k=1;k<=2;k++) {

      // Skip if only including gluinos
      if (settingsPtr->flag("SUSY:qq2squarksquark:onlyQCD")) continue;
      
      for (int l=1;l<=2;l++) {

	// kl-dependent factor for LL and RR contributions
	double facMS = sH*sqrt(m2Char[k]*m2Char[l]);

	// Note: Ckl defined as in [Boz07] with sigmaChar factored out
	// [1][1] = LL, [1][2] = LR, [2][1] = RL, [2][2] = RR
	complex Ckl[3][3];
	Ckl[1][1] = facMS
	  * coupSUSYPtr->LsudX[iGen4][iGen2][k]
	  * conj(coupSUSYPtr->LsduX[iGen3][iGen1][k])
	  * conj(coupSUSYPtr->LsudX[iGen4][iGen2][l])
	  * coupSUSYPtr->LsduX[iGen3][iGen1][l];
	Ckl[1][2] = facTU
	  * coupSUSYPtr->RsudX[iGen4][iGen2][k]
	  * conj(coupSUSYPtr->LsduX[iGen3][iGen1][k])
	  * conj(coupSUSYPtr->RsudX[iGen4][iGen2][l])
	  * coupSUSYPtr->LsduX[iGen3][iGen1][l];
	Ckl[2][1] = facTU 
	  * coupSUSYPtr->LsudX[iGen4][iGen2][k]
	  * conj(coupSUSYPtr->RsduX[iGen3][iGen1][k])
	  * conj(coupSUSYPtr->LsudX[iGen4][iGen2][l])
	  * coupSUSYPtr->RsduX[iGen3][iGen1][l];
	Ckl[2][2] = facMS  
	  * coupSUSYPtr->RsudX[iGen4][iGen2][k]
	  * conj(coupSUSYPtr->RsduX[iGen3][iGen1][k])
	  * conj(coupSUSYPtr->RsudX[iGen4][iGen2][l])
	  * coupSUSYPtr->RsduX[iGen3][iGen1][l];

	// Add to sum of t-channel charginos
	sumCt += sigmaChar * real(Ckl[1][1] + Ckl[1][2] + Ckl[2][1] 
               + Ckl[2][2]) / tChar[k] / tChar[l];

      }
    }
    
    // u-channel Neutralinos
    for (int k=1;k<=nNeut;k++) {

      // Skip if only including gluinos
      if (settingsPtr->flag("SUSY:qq2squarksquark:onlyQCD")) continue;
      
      for (int l=1;l<=nNeut;l++) {

	// kl-dependent factor for LL, RR contributions
	double facMS = sH*sqrt(m2Neut[k]*m2Neut[l]);

	// Note: Nkl defined as in [Boz07] with sigmaNeut factored out
	// [1][1] = LL, [1][2] = LR, [2][1] = RL, [2][2] = RR
	complex Nkl[3][3];
	Nkl[1][1] = facMS
	  * conj(coupSUSYPtr->LsuuX[iGen4][iGen1][k])
	  * conj(coupSUSYPtr->LsddX[iGen3][iGen2][k])
	  * coupSUSYPtr->LsuuX[iGen4][iGen1][l]
	  * coupSUSYPtr->LsddX[iGen3][iGen2][l];
	Nkl[1][2] = facTU 
	  * conj(coupSUSYPtr->LsuuX[iGen4][iGen1][k])
	  * conj(coupSUSYPtr->RsddX[iGen3][iGen2][k])
	  * coupSUSYPtr->LsuuX[iGen4][iGen1][l]
	  * coupSUSYPtr->RsddX[iGen3][iGen2][l];
	Nkl[2][1] =  facTU 
	  * conj(coupSUSYPtr->RsuuX[iGen4][iGen1][k])
	  * conj(coupSUSYPtr->LsddX[iGen3][iGen2][k])
	  * coupSUSYPtr->RsuuX[iGen4][iGen1][l]
	  * coupSUSYPtr->LsddX[iGen3][iGen2][l];
	Nkl[2][2] =  facMS 
	  * conj(coupSUSYPtr->RsuuX[iGen4][iGen1][k])
	  * conj(coupSUSYPtr->RsddX[iGen3][iGen2][k])
	  * coupSUSYPtr->RsuuX[iGen4][iGen1][l]
	  * coupSUSYPtr->RsddX[iGen3][iGen2][l];

	// Add to sum of u-channel neutralinos	
	sumNu += sigmaNeut / uNeut[k] / uNeut[l] 
	  * real(Nkl[1][1] + Nkl[1][2] + Nkl[2][1] + Nkl[2][2]);
	  	
      }
    }

    // u-channel gluino
    // Note: Gij defined as in [Boz07] with sigmaGlu factored out
    // [1][1] = LL, [1][2] = LR, [2][1] = RL, [2][2] = RR
    double Gij[3][3];
    Gij[1][1] = norm(coupSUSYPtr->LsuuG[iGen4][iGen1]
                * coupSUSYPtr->LsddG[iGen3][iGen2]);
    Gij[1][2] = norm(coupSUSYPtr->LsuuG[iGen4][iGen1]
                * coupSUSYPtr->RsddG[iGen3][iGen2]);
    Gij[2][1] = norm(coupSUSYPtr->RsuuG[iGen4][iGen1]
                * coupSUSYPtr->LsddG[iGen3][iGen2]);
    Gij[2][2] = norm(coupSUSYPtr->RsuuG[iGen4][iGen1]
                * coupSUSYPtr->RsddG[iGen3][iGen2]);
    Gij[1][1] *= sH*m2Glu;
    Gij[1][2] *= facTU;
    Gij[2][1] *= facTU;
    Gij[2][2] *= sH*m2Glu;
    // Sum over polarizations
    sumGu += sigmaGlu * (Gij[1][1] + Gij[1][2] + Gij[2][1] + Gij[2][2]) 
      / pow2(uGlu); 

    // chargino-neutralino interference
    for (int k=1;k<=2;k++) {
      for (int l=1;l<=nNeut;l++) {

	// Skip if only including gluinos
	if (settingsPtr->flag("SUSY:qq2squarksquark:onlyQCD")) continue;

	// Note: CNkl defined as in [Boz07] with pi/sH2 factored out
	// [1][1] = LL, [1][2] = LR, [2][1] = RL, [2][2] = RR
	double CNkl[3][3];
	CNkl[1][1] = real(coupSUSYPtr->LsudX[iGen4][iGen2][k] 
			  * conj(coupSUSYPtr->LsduX[iGen3][iGen1][k]) 
			  * coupSUSYPtr->LsuuX[iGen4][iGen1][l] 
			  * coupSUSYPtr->LsddX[iGen3][iGen2][l]);
	CNkl[1][2] = real(coupSUSYPtr->RsudX[iGen4][iGen2][k] 
			  * conj(coupSUSYPtr->LsduX[iGen3][iGen1][k]) 
			  * coupSUSYPtr->LsuuX[iGen4][iGen1][l] 
			  * coupSUSYPtr->RsddX[iGen3][iGen2][l]);
	CNkl[2][1] = real(coupSUSYPtr->LsudX[iGen4][iGen2][k] 
			  * conj(coupSUSYPtr->RsduX[iGen3][iGen1][k]) 
			  * coupSUSYPtr->RsuuX[iGen4][iGen1][l] 
			  * coupSUSYPtr->LsddX[iGen3][iGen2][l]);
	CNkl[2][2] = real(coupSUSYPtr->RsudX[iGen4][iGen2][k] 
			  * conj(coupSUSYPtr->RsduX[iGen3][iGen1][k]) 
			  * coupSUSYPtr->RsuuX[iGen4][iGen1][l] 
			  * coupSUSYPtr->RsddX[iGen3][iGen2][l]);
	CNkl[1][1] *= sH*sqrt(m2Char[k]*m2Neut[l]);
	CNkl[1][2] *= uH*tH-s3*s4;
	CNkl[2][1] *= uH*tH-s3*s4;
	CNkl[2][2] *= sH*sqrt(m2Char[k]*m2Neut[l]);	
	// Sum over polarizations
	sumInterference += sigmaCharNeut * (CNkl[1][1] + CNkl[1][2] 
                         + CNkl[2][1] + CNkl[2][2]) / tChar[k] / uNeut[l];
      }
    }

    // chargino-gluino interference
    for (int k=1;k<=2;k++) {

      // Skip if only including gluinos
      if (settingsPtr->flag("SUSY:qq2squarksquark:onlyQCD")) continue;
	
      // Note: CGk defined as in [Boz07] with sigmaCharGlu factored out
      // [1][1] = LL, [1][2] = LR, [2][1] = RL, [2][2] = RR
      double CGk[3][3];
      CGk[1][1] = real(coupSUSYPtr->LsudX[iGen4][iGen2][k]
		       * conj(coupSUSYPtr->LsduX[iGen3][iGen1][k])
		       * conj(coupSUSYPtr->LsuuG[iGen4][iGen1])
		       * conj(coupSUSYPtr->LsddG[iGen3][iGen2]));
      CGk[1][2] = real(coupSUSYPtr->RsudX[iGen4][iGen2][k]
		       * conj(coupSUSYPtr->LsduX[iGen3][iGen1][k])
		       * conj(coupSUSYPtr->LsuuG[iGen4][iGen1])
		       * conj(coupSUSYPtr->RsddG[iGen3][iGen2]));
      CGk[2][1] = real(coupSUSYPtr->LsudX[iGen4][iGen2][k]
		       * conj(coupSUSYPtr->RsduX[iGen3][iGen1][k])
		       * conj(coupSUSYPtr->RsuuG[iGen4][iGen1])
		       * conj(coupSUSYPtr->LsddG[iGen3][iGen2]));
      CGk[2][2] = real(coupSUSYPtr->RsudX[iGen4][iGen2][k]
		       * conj(coupSUSYPtr->RsduX[iGen3][iGen1][k])
		       * conj(coupSUSYPtr->RsuuG[iGen4][iGen1])
		       * conj(coupSUSYPtr->RsddG[iGen3][iGen2]));
      CGk[1][1] *= sH*sqrt(m2Glu*m2Char[k]);
      CGk[1][2] *= uH*tH-s3*s4;
      CGk[2][1] *= uH*tH-s3*s4;
      CGk[2][2] *= sH*sqrt(m2Glu*m2Char[k]);
      // Sum over polarizations
      sumInterference += sigmaGlu * (CGk[1][1] + CGk[1][2] + CGk[2][1] 
        + CGk[2][2]) / uGlu / tChar[k]; 
    }    
  }

  // Case B) Same-isospin: qq' -> ~d~d , ~u~u
  else {    

    // t-channel + u-channel Neutralinos + t/u interference
    for (int k=1;k<=nNeut;k++) {

      // Skip if only including gluinos
      if (settingsPtr->flag("SUSY:qq2squarksquark:onlyQCD")) continue;

      for (int l=1;l<=nNeut;l++) {

	// kl-dependent factor for LL and RR contributions
	double facMS = sH * particleDataPtr->m0(coupSUSYPtr->idNeut(k))
	  * particleDataPtr->m0(coupSUSYPtr->idNeut(l));

	// Note: Nxkl defined as in [Boz07] with sigmaNeut factored out
	// [1][1] = LL, [1][2] = LR, [2][1] = RL, [2][2] = RR
	complex NTkl[3][3], NUkl[3][3], NTUkl[3][3];
	NTkl[1][1] = facMS 
	  * conj(coupSUSYPtr->getLsqqX(iGen4,idIn2A,k))
	  * conj(coupSUSYPtr->getLsqqX(iGen3,idIn1A,k))
	  * coupSUSYPtr->getLsqqX(iGen4,idIn2A,l)
	  * coupSUSYPtr->getLsqqX(iGen3,idIn1A,l);
	NTkl[1][2] = facTU 
	  * conj(coupSUSYPtr->getRsqqX(iGen4,idIn2A,k))
	  * conj(coupSUSYPtr->getLsqqX(iGen3,idIn1A,k))
	  * coupSUSYPtr->getRsqqX(iGen4,idIn2A,l)
	  * coupSUSYPtr->getLsqqX(iGen3,idIn1A,l);
	NTkl[2][1] = facTU 
	  * conj(coupSUSYPtr->getLsqqX(iGen4,idIn2A,k))
	  * conj(coupSUSYPtr->getRsqqX(iGen3,idIn1A,k))
	  * coupSUSYPtr->getLsqqX(iGen4,idIn2A,l)
	  * coupSUSYPtr->getRsqqX(iGen3,idIn1A,l);
	NTkl[2][2] = facMS  
	  * conj(coupSUSYPtr->getRsqqX(iGen4,idIn2A,k))
	  * conj(coupSUSYPtr->getRsqqX(iGen3,idIn1A,k))
	  * coupSUSYPtr->getRsqqX(iGen4,idIn2A,l)
	  * coupSUSYPtr->getRsqqX(iGen3,idIn1A,l);
	NUkl[1][1] = facMS
	  * conj(coupSUSYPtr->getLsqqX(iGen3,idIn2A,k))
	  * conj(coupSUSYPtr->getLsqqX(iGen4,idIn1A,k))
	  * coupSUSYPtr->getLsqqX(iGen3,idIn2A,l)
	  * coupSUSYPtr->getLsqqX(iGen4,idIn1A,l);
	NUkl[1][2] = facTU 
	  * conj(coupSUSYPtr->getRsqqX(iGen3,idIn2A,k))
	  * conj(coupSUSYPtr->getLsqqX(iGen4,idIn1A,k))
	  * coupSUSYPtr->getRsqqX(iGen3,idIn2A,l)
	  * coupSUSYPtr->getLsqqX(iGen4,idIn1A,l);
	NUkl[2][1] = facTU
	  * conj(coupSUSYPtr->getLsqqX(iGen3,idIn2A,k))
	  * conj(coupSUSYPtr->getRsqqX(iGen4,idIn1A,k))
	  * coupSUSYPtr->getLsqqX(iGen3,idIn2A,l)
	  * coupSUSYPtr->getRsqqX(iGen4,idIn1A,l);
	NUkl[2][2] = facMS
	  * conj(coupSUSYPtr->getRsqqX(iGen3,idIn2A,k))
	  * conj(coupSUSYPtr->getRsqqX(iGen4,idIn1A,k))
	  * coupSUSYPtr->getRsqqX(iGen3,idIn2A,l)
	  * coupSUSYPtr->getRsqqX(iGen4,idIn1A,l);
	NTUkl[1][1] = facMS  
	  * real( conj(coupSUSYPtr->getLsqqX(iGen4,idIn2A,k))
		  * conj(coupSUSYPtr->getLsqqX(iGen3,idIn1A,k))
		  * coupSUSYPtr->getLsqqX(iGen3,idIn2A,l)
		  * coupSUSYPtr->getLsqqX(iGen4,idIn1A,l) );
	NTUkl[1][2] = facTU 
	  * real( conj(coupSUSYPtr->getRsqqX(iGen4,idIn2A,k))
		  * conj(coupSUSYPtr->getLsqqX(iGen3,idIn1A,k))
		  * coupSUSYPtr->getRsqqX(iGen3,idIn2A,l)
		  * coupSUSYPtr->getLsqqX(iGen4,idIn1A,l) );
	NTUkl[2][1] = facTU 
	  * real( conj(coupSUSYPtr->getLsqqX(iGen4,idIn2A,k))
		  * conj(coupSUSYPtr->getRsqqX(iGen3,idIn1A,k))
		  * coupSUSYPtr->getLsqqX(iGen3,idIn2A,l)
		  * coupSUSYPtr->getRsqqX(iGen4,idIn1A,l) );
	NTUkl[2][2] = facMS  
	  * real( conj(coupSUSYPtr->getRsqqX(iGen4,idIn2A,k))
		  * conj(coupSUSYPtr->getRsqqX(iGen3,idIn1A,k))
		  * coupSUSYPtr->getRsqqX(iGen3,idIn2A,l)
		  * coupSUSYPtr->getRsqqX(iGen4,idIn1A,l) );
	  
	// Add to sums
	sumNt += sigmaNeut / tNeut[k] / tNeut[l] 
	  * real(NTkl[1][1] + NTkl[1][2] + NTkl[2][1] + NTkl[2][2]);
	sumNu += sigmaNeut / uNeut[k] / uNeut[l]
	  * real(NUkl[1][1] + NUkl[1][2] + NUkl[2][1] + NUkl[2][2]);
	sumInterference += 2.0 / 3.0 * sigmaNeut 
	  * real(NTUkl[1][1] + NTUkl[1][2] + NTUkl[2][1] + NTUkl[2][2])
	  / tNeut[k] / uNeut[l];
      }     

      // Neutralino / Gluino interference

      // k-dependent factor for LL and RR contributions
      double facMS = sH * particleDataPtr->m0(coupSUSYPtr->idNeut(k))
	* particleDataPtr->m0(1000021);
      
      // Note: Nxkl defined as in [Boz07] with sigmaNeutGlu factored out
      // [1][1] = LL, [1][2] = LR, [2][1] = RL, [2][2] = RR
      complex NGA[3][3], NGB[3][3];
      NGA[1][1] = facMS
	* real( conj(coupSUSYPtr->getLsqqX(iGen4,idIn2A,k))
		* conj(coupSUSYPtr->getLsqqX(iGen3,idIn1A,k))
		* conj(coupSUSYPtr->getLsqqG(iGen3,idIn2A))
		* conj(coupSUSYPtr->getLsqqG(iGen4,idIn1A)) );
      NGA[1][2] = facTU
	* real( conj(coupSUSYPtr->getRsqqX(iGen4,idIn2A,k))
		* conj(coupSUSYPtr->getLsqqX(iGen3,idIn1A,k))
		* conj(coupSUSYPtr->getLsqqG(iGen3,idIn2A))
		* conj(coupSUSYPtr->getRsqqG(iGen4,idIn1A)) );
      NGA[2][1] = facTU
	* real( conj(coupSUSYPtr->getLsqqX(iGen4,idIn2A,k))
		* conj(coupSUSYPtr->getRsqqX(iGen3,idIn1A,k))
		* conj(coupSUSYPtr->getRsqqG(iGen3,idIn2A))
		* conj(coupSUSYPtr->getLsqqG(iGen4,idIn1A)) );
      NGA[2][2] = facMS
	* real( conj(coupSUSYPtr->getRsqqX(iGen4,idIn2A,k))
		* conj(coupSUSYPtr->getRsqqX(iGen3,idIn1A,k))
		* conj(coupSUSYPtr->getRsqqG(iGen3,idIn2A))
		* conj(coupSUSYPtr->getRsqqG(iGen4,idIn1A)) );
      NGB[1][1] = facMS
	* real( conj(coupSUSYPtr->getLsqqX(iGen3,idIn2A,k))
		* conj(coupSUSYPtr->getLsqqX(iGen4,idIn1A,k))
		* conj(coupSUSYPtr->getLsqqG(iGen4,idIn2A))
		* conj(coupSUSYPtr->getLsqqG(iGen3,idIn1A)) );
      NGB[1][2] = facMS
	* real( conj(coupSUSYPtr->getRsqqX(iGen3,idIn2A,k))
		* conj(coupSUSYPtr->getLsqqX(iGen4,idIn1A,k))
		* conj(coupSUSYPtr->getRsqqG(iGen4,idIn2A))
		* conj(coupSUSYPtr->getLsqqG(iGen3,idIn1A)) );
      NGB[2][1] = facMS
	* real( conj(coupSUSYPtr->getLsqqX(iGen3,idIn2A,k))
		* conj(coupSUSYPtr->getRsqqX(iGen4,idIn1A,k))
		* conj(coupSUSYPtr->getLsqqG(iGen4,idIn2A))
		* conj(coupSUSYPtr->getRsqqG(iGen3,idIn1A)) );
      NGB[2][2] = facMS
	* real( conj(coupSUSYPtr->getRsqqX(iGen3,idIn2A,k))
		* conj(coupSUSYPtr->getRsqqX(iGen4,idIn1A,k))
		* conj(coupSUSYPtr->getRsqqG(iGen4,idIn2A))
		* conj(coupSUSYPtr->getRsqqG(iGen3,idIn1A)) );

      // Add to sums
      sumInterference += sigmaNeutGlu * 
	( real(NGA[1][1] + NGA[1][2] + NGA[2][1] + NGA[2][2]) 
        / tNeut[k] / uGlu
	+ real(NGB[1][1] + NGB[1][2] + NGB[2][1] + NGB[2][2]) 
        / uNeut[k] / tGlu );
    }
    
    // t-channel + u-channel Gluinos + t/u interference

    // factor for LL and RR contributions
    double facMS = sH * m2Glu;

    // Note: GT, GU defined as in [Boz07] with sigmaGlu factored out
    // [1][1] = LL, [1][2] = LR, [2][1] = RL, [2][2] = RR
    complex GT[3][3], GU[3][3], GTU[3][3];
    GT[1][1] = facMS 
      * norm(coupSUSYPtr->getLsqqG(iGen4,idIn2A) 
      * coupSUSYPtr->getLsqqG(iGen3,idIn1A));
    GT[1][2] = facTU
      * norm(coupSUSYPtr->getRsqqG(iGen4,idIn2A) 
      * coupSUSYPtr->getLsqqG(iGen3,idIn1A));
    GT[2][1] = facTU
      * norm(coupSUSYPtr->getLsqqG(iGen4,idIn2A) 
      * coupSUSYPtr->getRsqqG(iGen3,idIn1A));
    GT[2][2] = facMS
      * norm(coupSUSYPtr->getRsqqG(iGen4,idIn2A) 
      * coupSUSYPtr->getRsqqG(iGen3,idIn1A));
    GU[1][1] = facMS 
      * norm(coupSUSYPtr->getLsqqG(iGen3,idIn2A) 
      * coupSUSYPtr->getLsqqG(iGen4,idIn1A));
    GU[1][2] = facTU
      * norm(coupSUSYPtr->getLsqqG(iGen3,idIn2A) 
      * coupSUSYPtr->getRsqqG(iGen4,idIn1A));
    GU[2][1] = facTU
      * norm(coupSUSYPtr->getRsqqG(iGen3,idIn2A) 
      * coupSUSYPtr->getLsqqG(iGen4,idIn1A));
    GU[2][2] = facMS
      * norm(coupSUSYPtr->getRsqqG(iGen3,idIn2A) 
      * coupSUSYPtr->getRsqqG(iGen4,idIn1A));

    GTU[1][1] = facMS 
      * real(coupSUSYPtr->getLsqqG(iGen3,idIn1A) 
	     * coupSUSYPtr->getLsqqG(iGen4,idIn2A)
	     * conj(coupSUSYPtr->getLsqqG(iGen3,idIn2A))
	     * conj(coupSUSYPtr->getLsqqG(iGen4,idIn1A)) );

    GTU[1][2] = facTU 
      * real(coupSUSYPtr->getLsqqG(iGen3,idIn1A) 
	     * coupSUSYPtr->getRsqqG(iGen4,idIn2A)
	     * conj(coupSUSYPtr->getRsqqG(iGen3,idIn2A))
	     * conj(coupSUSYPtr->getLsqqG(iGen4,idIn1A)) );

    GTU[2][1] = facTU 
      * real(coupSUSYPtr->getRsqqG(iGen3,idIn1A) 
	     * coupSUSYPtr->getLsqqG(iGen4,idIn2A)
	     * conj(coupSUSYPtr->getLsqqG(iGen3,idIn2A))
	     * conj(coupSUSYPtr->getRsqqG(iGen4,idIn1A)) );

    GTU[2][2] = facMS 
      * real(coupSUSYPtr->getRsqqG(iGen3,idIn1A) 
	     * coupSUSYPtr->getRsqqG(iGen4,idIn2A)
	     * conj(coupSUSYPtr->getRsqqG(iGen3,idIn2A))
	     * conj(coupSUSYPtr->getRsqqG(iGen4,idIn1A)) );
 
    // Add to sums
    sumGt += sigmaGlu * real(GT[1][1] + GT[1][2] + GT[2][1] + GT[2][2])
      / pow2(tGlu) ;
    sumGu += sigmaGlu * real(GU[1][1] + GU[1][2] + GU[2][1] + GU[2][2])
      / pow2(uGlu) ;
    sumInterference += - 2.0 / 3.0 * sigmaGlu
      * real(GTU[1][1] + GTU[1][2] + GTU[2][1] + GTU[2][2])
      / tGlu / uGlu;    

  }

  // Cross section
  double sigma = sumNt + sumNu + sumCt + sumCu + sumGt + sumGu 
    + sumInterference;

  // Identical particles?
  if (id3Sav == id4Sav) sigma /= 2.0;

  // Return answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qq2squarksquark::setIdColAcol() {

  // Set flavours.
  if (id1 > 0 && id2 > 0) {
    setId( id1, id2, id3Sav, id4Sav);
  } else {
    // 1,2 -> -3,-4
    setId( id1, id2,-id3Sav,-id4Sav);
  }

  // Coded sigma is for ud -> ~q~q'. Swap t and u for du -> ~q~q'.
  swapTU = (isUD && abs(id1) % 2 == 0); 

  // Select colour flow topology 
  // A: t-channel neutralino, t-channel chargino, or u-channel gluino
  double fracA = sumNt + sumCt + sumGu 
    / (sumNt + sumNu + sumCt + sumCu + sumGt + sumGu);
  if (swapTU) fracA = 1.0 - fracA;
  setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  // B: t-channel gluino or u-channel neutralino 
  if (rndmPtr->flat() > fracA) setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);

  // Switch to anti-colors if antiquarks
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2squarkantisquark
// Cross section for qqbar-initiated squark-antisquark production

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qqbar2squarkantisquark::initProc() {

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  // Extract isospin and mass-ordering indices
  iGen3 = 3*(abs(id3Sav)/2000000) + (abs(id3Sav)%10+1)/2;
  iGen4 = 3*(abs(id4Sav)/2000000) + (abs(id4Sav)%10+1)/2;

  // Is this a ~u_i ~d*_j, ~d_i ~u*_j final state or ~d_i ~d*_j, ~u_i ~u*_j
  if (abs(id3Sav) % 2 == abs(id4Sav) % 2) isUD = false;
  else isUD = true;

  // Derive name
  nameSave = "q qbar' -> "+particleDataPtr->name(abs(id3Sav))+" "+
    particleDataPtr->name(-abs(id4Sav));
  if (isUD && abs(id3Sav) != abs(id4Sav)) nameSave +=" + c.c.";

  // Count 5 neutralinos in NMSSM
  nNeut = (coupSUSYPtr->isNMSSM ? 5 : 4);

  // Store mass squares of all possible internal propagator lines
  m2Glu     = pow2(particleDataPtr->m0(1000021));
  m2Neut.resize(nNeut+1);
  for (int iNeut=1;iNeut<=nNeut;iNeut++) 
    m2Neut[iNeut] = pow2(particleDataPtr->m0(coupSUSYPtr->idNeut(iNeut)));
  
  // Set sizes of some arrays to be used below
  tNeut.resize(nNeut+1);
  uNeut.resize(nNeut+1);

  // Shorthand for Weak mixing
  xW = coupSUSYPtr->sin2W;

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(id3Sav, id4Sav);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qqbar2squarkantisquark::sigmaKin() {

  // Z/W propagator
  if (! isUD) {
    double sV= sH - pow2(coupSUSYPtr->mZpole);
    double d = pow2(sV) + pow2(coupSUSYPtr->mZpole * coupSUSYPtr->wZpole);
    propZW   = complex( sV / d, coupSUSYPtr->mZpole * coupSUSYPtr->wZpole / d);
  } else {
    double sV= sH - pow2(coupSUSYPtr->mWpole);
    double d = pow2(sV) + pow2(coupSUSYPtr->mWpole * coupSUSYPtr->wWpole);
    propZW   = complex( sV / d, coupSUSYPtr->mWpole * coupSUSYPtr->wWpole / d);
  }

  // Flavor-independent pre-factors
  double comFacHat = M_PI/sH2 * openFracPair;

  sigmaEW       = comFacHat * pow2(alpEM);
  sigmaGlu      = comFacHat * 2.0 * pow2(alpS) / 9.0;
  sigmaEWG      = comFacHat * 8.0 * alpEM * alpS / 9.0;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qqbar2squarkantisquark::sigmaHat() {

  // In-pair must be opposite-sign
  if (id1 * id2 > 0) return 0.0;
  
  // Check correct charge sum
  if (isUD && abs(id1) %2 == abs(id2) % 2) return 0.0;
  if (!isUD && abs(id1) % 2 != abs(id2) % 2) return 0.0;

  // Check if using QCD diagrams only
  bool onlyQCD = settingsPtr->flag("SUSY:qqbar2squarkantisquark:onlyQCD");

  // Coded UD sigma is for udbar -> ~u~d'*. Swap t<->u for dbaru -> ~u~d'*.
  swapTU = (isUD && abs(id1) % 2 != 0); 

  // Coded QQ sigma is for qqbar -> ~q~q*. Swap t<->u for qbarq -> ~q~q*.
  if (!isUD && id1 < 0) swapTU = true;

  // Generation indices of incoming particles
  int idIn1A = (swapTU) ? abs(id2) : abs(id1);
  int idIn2A = (swapTU) ? abs(id1) : abs(id2);
  int iGen1  = (idIn1A+1)/2;
  int iGen2  = (idIn2A+1)/2;  

  // Auxiliary factors for use below  
  tGlu     = tH - m2Glu;
  uGlu     = uH - m2Glu;
  for (int i=1; i<= nNeut; i++) {
    tNeut[i] = tH - m2Neut[i];
    uNeut[i] = uH - m2Neut[i];
  }

  // Initial values for pieces used for color-flow selection below
  sumColS   = 0.0;
  sumColT   = 0.0;
  sumInterference = 0.0;

  // Common factor for LR and RL contributions
  double facTU =  uH*tH-s3*s4;

  // Case A) Opposite-isospin: udbar -> ~u~d* 
  if ( isUD ) {

    // s-channel W contribution (only contributes to LL helicities)
    if ( !onlyQCD ) {
      sumColS += sigmaEW / 16.0 / pow2(xW) / pow2(1.0-xW)
	* norm(conj(coupSUSYPtr->LudW[iGen1][iGen2])
	       * coupSUSYPtr->LsusdW[iGen3][iGen4]) * facTU  
	* norm(propZW);
    }

    // t-channel gluino contributions
    double GT[3][3];
    double facLR = m2Glu * sH;
    // LL, LR, RL, RR
    GT[1][1] = facTU * norm(conj(coupSUSYPtr->LsddG[iGen4][iGen2])
			    *coupSUSYPtr->LsuuG[iGen3][iGen1]);
    GT[1][2] = facLR * norm(conj(coupSUSYPtr->RsddG[iGen4][iGen2])
			    *coupSUSYPtr->LsuuG[iGen3][iGen1]);
    GT[2][1] = facLR * norm(conj(coupSUSYPtr->LsddG[iGen4][iGen2])
			    *coupSUSYPtr->RsuuG[iGen3][iGen1]);
    GT[2][2] = facTU * norm(conj(coupSUSYPtr->RsddG[iGen4][iGen2])
			    *coupSUSYPtr->RsuuG[iGen3][iGen1]);
    // leading color flow for t-channel gluino is annihilation-like
    sumColS += sigmaGlu / pow2(tGlu)
      * (GT[1][1] + GT[1][2] + GT[2][1] + GT[2][2]);
      
    // W-Gluino interference (only contributes to LL helicities)
    if ( !onlyQCD ) {
      sumColS += sigmaEWG / 4.0 / xW / (1-xW) 
	* real(conj(coupSUSYPtr->LsuuG[iGen3][iGen1])
	       * coupSUSYPtr->LsddG[iGen4][iGen2]
	       * conj(coupSUSYPtr->LudW[iGen1][iGen2])
	       * coupSUSYPtr->LsusdW[iGen3][iGen4]) * facTU 
	/ tGlu * sqrt(norm(propZW));
    }

    // t-channel neutralinos
    // NOT YET IMPLEMENTED !

  }

  // Case B) Same-isospin: qqbar -> ~d~d* , ~u~u*
  else {    
    
    double eQ  = (idIn1A % 2 == 0) ? 2./3. : 1./3. ;
    double eSq = (abs(id3Sav) % 2 == 0) ? 2./3. : 1./3. ;

    // s-channel gluon (strictly flavor-diagonal)
    if (abs(id3Sav) == abs(id4Sav) && abs(id1) == abs(id2)) {
      // Factor 2 since contributes to both ha != hb helicities
      sumColT += 2. * sigmaGlu * facTU / pow2(sH);      
    }

    // t-channel gluino (only for in-isospin = out-isospin). 
    if (eQ == eSq) {
      // Sum over helicities.     
      double GT[3][3];
      double facLR = sH * m2Glu;
      GT[1][1] = facTU * norm(coupSUSYPtr->getLsqqG(iGen3,idIn1A)
			      * conj(coupSUSYPtr->getLsqqG(iGen4,idIn2A)));
      GT[1][2] = facLR * norm(coupSUSYPtr->getLsqqG(iGen3,idIn1A)
			      * conj(coupSUSYPtr->getRsqqG(iGen4,idIn2A)));
      GT[2][1] = facLR * norm(coupSUSYPtr->getRsqqG(iGen3,idIn1A)
			      * conj(coupSUSYPtr->getLsqqG(iGen4,idIn2A)));
      GT[2][2] = facTU * norm(coupSUSYPtr->getRsqqG(iGen3,idIn1A)
			      * conj(coupSUSYPtr->getRsqqG(iGen4,idIn2A)));    
      // Add contribution to color topology: S
      sumColS += sigmaGlu / pow2(tGlu)
	* ( GT[1][1] + GT[1][2] + GT[2][1] + GT[2][2]);
      
      // gluon-gluino interference (strictly flavor-diagonal)
      if (abs(id3Sav) == abs(id4Sav) && abs(id1) == abs(id2)) {
	double GG11, GG22;
	GG11 = - facTU * 2./3. 
             * real( conj(coupSUSYPtr->getLsqqG(iGen3,idIn1A))
	     * coupSUSYPtr->getLsqqG(iGen4,idIn2A));
	GG22 = - facTU * 2./3. 
             * real( conj(coupSUSYPtr->getRsqqG(iGen3,idIn1A))
	     * coupSUSYPtr->getRsqqG(iGen4,idIn2A)); 
	// Sum over two contributing helicities
	sumInterference += sigmaGlu / sH / tGlu
	  * ( GG11 + GG22 );
      }

    }

    // Skip the rest if only including QCD diagrams
    if (onlyQCD) return sumColT+sumColS+sumInterference;

    // s-channel photon (strictly flavor-diagonal) and Z/gamma interference
    if (abs(id3Sav) == abs(id4Sav) && abs(id1) == abs(id2)) {

      // gamma
      // Factor 2 since contributes to both ha != hb helicities
      sumColS += 2. * pow2(eQ) * pow2(eSq) * sigmaEW * facTU / pow2(sH);

      // Z/gamma interference
      double CsqZ = real(coupSUSYPtr->LsusuZ[iGen3][iGen4] 
			 + coupSUSYPtr->RsusuZ[iGen3][iGen4]);
      if (abs(id3Sav)%2 != 0) CsqZ = real(coupSUSYPtr->LsdsdZ[iGen3][iGen4] 
					  + coupSUSYPtr->RsdsdZ[iGen3][iGen4]);
      sumColS += eQ * eSq * sigmaEW * facTU / 2.0 / xW / (1.-xW) 
	* sqrt(norm(propZW)) / sH * CsqZ
	* (coupSUSYPtr->LqqZ[idIn1A] + coupSUSYPtr->LqqZ[idIn2A]);
      
      // Gluino/gamma interference (only for same-isospin)
      if (eQ == eSq) {
	double CsqG11 = real(conj(coupSUSYPtr->LsuuG[iGen3][iGen1]) 
			     *coupSUSYPtr->LsuuG[iGen4][iGen2]);
	double CsqG22 = real(conj(coupSUSYPtr->RsuuG[iGen3][iGen1]) 
			     *coupSUSYPtr->RsuuG[iGen4][iGen2]);
	if (id3Sav%2 != 0) {
	  CsqG11 = real(conj(coupSUSYPtr->LsddG[iGen3][iGen1]) 
			*coupSUSYPtr->LsddG[iGen4][iGen2]);
	  CsqG22 = real(conj(coupSUSYPtr->RsddG[iGen3][iGen1]) 
			*coupSUSYPtr->RsddG[iGen4][iGen2]);
	}
	sumColS += eQ * eSq * sigmaEWG * facTU
	  * (CsqG11 + CsqG22) / sH / tGlu; 
      }
    }
    
    // s-channel Z (only for q flavor = qbar flavor)
    if (abs(id1) == abs(id2)) {
      double CsqZ = norm(coupSUSYPtr->LsusuZ[iGen3][iGen4] 
			 + coupSUSYPtr->RsusuZ[iGen3][iGen4]);
      if (abs(id3Sav)%2 != 0) CsqZ = norm(coupSUSYPtr->LsdsdZ[iGen3][iGen4] 
					  + coupSUSYPtr->RsdsdZ[iGen3][iGen4]);
      sumColS += sigmaEW * facTU / 16.0 / pow2(xW) / pow2(1.0-xW)
	* norm(propZW) * CsqZ * ( pow2(coupSUSYPtr->LqqZ[idIn1A]) 
        + pow2(coupSUSYPtr->RqqZ[idIn1A]) );

      // Z/gluino interference (only for in-isospin = out-isospin)
      if (eQ == eSq) {
	double GZ11 = real(conj(coupSUSYPtr->getLsqqG(iGen3,idIn1A))
			   *coupSUSYPtr->getLsqqG(iGen4,idIn2A)	         
			   *(coupSUSYPtr->getLsqsqZ(id3Sav,id4Sav)
			     +coupSUSYPtr->getRsqsqZ(id3Sav,id4Sav)))
	  *coupSUSYPtr->LqqZ[idIn1A];
	double GZ22 = real(conj(coupSUSYPtr->getRsqqG(iGen3,idIn1A))
			   *coupSUSYPtr->getRsqqG(iGen4,idIn2A)
			   *(coupSUSYPtr->getLsqsqZ(id3Sav,id4Sav)
			     +coupSUSYPtr->getRsqsqZ(id3Sav,id4Sav)))
	  *coupSUSYPtr->RqqZ[idIn1A];
	sumColS += sigmaEWG * facTU / 4.0 / xW / (1.-xW) 
	  * ( GZ11 + GZ22 ) * sqrt(norm(propZW)) / tGlu;	
      }
    }
    
    // t-channel neutralinos
    // NOT YET IMPLEMENTED !
    
  }

  // Cross section
  double sigma = sumColS + sumColT + sumInterference;

  // Return answer.
  return sigma;
  
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2squarkantisquark::setIdColAcol() {

  // Check if charge conjugate final state?
  isCC = false;
  if (isUD && ( (id1-1)%2 < 0 || (id2-1)%2 < 0 )) isCC = true;
  
  //check if charge conjugate
  id3 = (isCC) ? -id3Sav : id3Sav;
  id4 = (isCC) ? -id4Sav : id4Sav;

  // Set flavours.
  setId( id1, id2, id3, id4);                 

  // Coded UD sigma is for udbar -> ~u~d'*. Swap t<->u for dbaru -> ~u~d'*.
  // Coded QQ sigma is for qqbar -> ~q~q*. Swap t<->u for qbarq -> ~q~q*.
  if (isUD) {
    swapTU = (abs(id1) % 2 != 0); 
  } else {
    swapTU = (id1 < 0);
  }

  // Select colour flow topology 
  double R = rndmPtr->flat();
  double fracS = sumColS / (sumColS + sumColT) ;
  // S: color flow as in S-channel singlet
  if (R < fracS) {
    setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
    if (swapTU) setColAcol( 0, 1, 1, 0, 2, 0, 0, 2);
  } 
  // T: color flow as in T-channel singlet
  else {
    setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
    if (swapTU) setColAcol( 0, 1, 2, 0, 2, 0, 0, 1);
  }

  if (isCC) swapColAcol();

}

//==========================================================================

// Sigma2gg2squarkantisquark
// Cross section for gg-initiated squark-antisquark production

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2gg2squarkantisquark::initProc() {

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  // Process Name
  nameSave = "g g -> "+particleDataPtr->name(abs(id3Sav))+" "
    +particleDataPtr->name(-abs(id4Sav)); 

  // Squark pole mass
  m2Sq = pow2(particleDataPtr->m0(id3Sav));
  
  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(id3Sav, id4Sav);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2gg2squarkantisquark::sigmaKin() {

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  // tHSq = tHat - m_squark^2; uHSq = uHat - m_squark^2. 
  //  double s34Avg = 0.5 * (s3 + s4) - 0.25 * pow2(s3 - s4) / sH; 
  double tHSq    = -0.5 * (sH - tH + uH);
  double uHSq    = -0.5 * (sH + tH - uH); 
  // ! (NEED TO CHECK THAT THESE APPLIED CORRECTLY BELOW)   ! 
  // ! (PRELIMINARY CROSS-CHECKS WITH PYTHIA 6 COME OUT OK) !

  // Helicity-independent prefactor
  double comFacHat = M_PI/sH2 * pow2(alpS) / 128.0
    * ( 24.0 * (1.0 - 2*tHSq*uHSq/sH2) - 8.0/3.0 );

  // Helicity-dependent factors
  sigma = 0.0;
  for (int ha=-1;ha<=1;ha += 2) {
    for (int hb=-1;hb<=1;hb += 2) {
      // Divide by 4 for helicity average      
      sigma += comFacHat / 4.0 
	* ( (1.0-ha*hb) 
	    - 2.0 * sH*m2Sq/tHSq/uHSq 
	    * ( 1.0 - ha*hb - sH*m2Sq/tHSq/uHSq));
    }    
  }  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2squarkantisquark::setIdColAcol() {

  // Set flavours.
  setId( id1, id2, id3Sav, id4Sav);                 

  // Set color flow (random for now)
  double R = rndmPtr->flat();
  if (R < 0.5) setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else setColAcol( 1, 2, 3, 1, 3, 0, 0, 2);

}

//==========================================================================

// Sigma2qg2squarkgluino
// Cross section for squark-gluino production

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qg2squarkgluino::initProc() {

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  // Derive name
  nameSave = "q g -> "+particleDataPtr->name(abs(id3Sav))+" gluino + c.c.";

  // Final-state mass squares
  m2Glu     = pow2(particleDataPtr->m0(1000021));
  m2Sq      = pow2(particleDataPtr->m0(id3Sav));

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(id3Sav, 1000021);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qg2squarkgluino::sigmaKin() {
  
  // Common pre-factor
  comFacHat = (M_PI / sH2) * pow2(alpS) * 0.5 * openFracPair;

  // Invariants (still with Pythia 6 sign convention)
  double tGlu = m2Glu-tH;
  double uGlu = m2Glu-uH;
  double tSq  = m2Sq-tH;
  double uSq  = m2Sq-uH;

  // Color flow A: quark color annihilates with anticolor of g
  sigmaA = 0.5*4./9.* tGlu/sH + (tGlu*sH+2.*m2Glu*tSq)/pow2(tGlu) -
    ( (sH-m2Sq+m2Glu)*(-tSq)-sH*m2Glu )/sH/(-tGlu) 
    + 0.5*1./2.*( tSq*(tH+2.*uH+m2Glu)-tGlu*(sH-2.*tSq) 
		  + (-uGlu)*(tH+m2Glu+2.*m2Sq) )/2./tGlu/uSq;
  // Color flow B: quark and gluon colors iterchanged
  sigmaB =     4./9.*(-uGlu)*(uH+m2Sq)/pow2(uSq) 
    + 1./18.* (sH*(uH+m2Glu) + 2.*(m2Sq-m2Glu)*uGlu)/sH/(-uSq) 
    + 0.5*4./9.*tGlu/sH 
    + 0.5*1./2.*(tSq*(tH+2.*uH+m2Glu)-tGlu*(sH-2.*tSq)
		 + (-uGlu)*(tH+m2Glu+2.*m2Sq))/2./tGlu/uSq;

}

double Sigma2qg2squarkgluino::sigmaHat() {
  
  // Check whether right incoming flavor
  int idQA = (id1 == 21) ? abs(id2) : abs(id1);
  int idSqSM = id3Sav%1000000;
  if (idQA != idSqSM) return 0.0;
  // NOTE: ONLY WORKS FOR SLHA1 ENUMERATION !!!
  // (should replace this by squark mixing matrix squares)

  return comFacHat * (sigmaA + sigmaB);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2squarkgluino::setIdColAcol() {

  // Check if charge conjugate final state?
  int idQ = (id1 == 21) ? id2 : id1;
  id3 = (idQ > 0) ? id3Sav : -id3Sav;
  id4 = 1000021;
  
  // Set flavors
  setId( id1, id2, id3, id4);                 

  // Select color flow A or B (see above)
  double R = rndmPtr->flat()*(sigmaA+sigmaB);  
  if (idQ == id1) {
    setColAcol(1,0,2,1,3,0,2,3);
    if (R > sigmaA) setColAcol(1,0,2,3,2,0,1,3);
  } else {
    setColAcol(2,1,1,0,3,0,2,3);
    if (R > sigmaB) setColAcol(2,3,1,0,2,0,1,3);    
  }
  if (idQ < 0) swapColAcol();

  // Use reflected kinematics if gq initial state
  if (id1 == 21) swapTU = true;

}

//==========================================================================

// Sigma2gg2gluinogluino
// Cross section for gluino pair production from gg initial states
// (validated against Pythia 6)

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2gg2gluinogluino::initProc() {

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(1000021, 1000021);
  
} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence. 

void Sigma2gg2gluinogluino::sigmaKin() { 

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  // tHG = tHat - m_gluino^2; uHG = uHat - m_gluino^2. 
  double s34Avg = 0.5 * (s3 + s4) - 0.25 * pow2(s3 - s4) / sH; 
  double tHG    = -0.5 * (sH - tH + uH);
  double uHG    = -0.5 * (sH + tH - uH); 
  double tHG2   = tHG * tHG;
  double uHG2   = uHG * uHG;

  // Calculate kinematics dependence.
  sigTS  = (tHG * uHG - 2. * s34Avg * (tHG + 2. * s34Avg)) / tHG2
         + (tHG * uHG + s34Avg * (uHG - tHG)) / (sH * tHG);  
  sigUS  = (tHG * uHG - 2. * s34Avg * (uHG + 2. * s34Avg)) / uHG2
         + (tHG * uHG + s34Avg * (tHG - uHG)) / (sH * uHG);
  sigTU  = 2. * tHG * uHG / sH2 + s34Avg * (sH - 4. * s34Avg) 
         / (tHG * uHG);
  sigSum = sigTS + sigUS + sigTU;
    
  // Answer contains factor 1/2 from identical gluinos.
  sigma  = (M_PI / sH2) * pow2(alpS) * (9./4.) * 0.5 * sigSum 
         * openFracPair;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2gluinogluino::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, 1000021, 1000021);

  // Three colour flow topologies, each with two orientations.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 4, 4, 3);
  else if (sigRand < sigTS + sigUS) 
                       setColAcol( 1, 2, 3, 1, 3, 4, 4, 2);
  else                 setColAcol( 1, 2, 3, 4, 1, 4, 3, 2); 
  if (rndmPtr->flat() > 0.5) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2gluinogluino
// Cross section for gluino pair production from qqbar initial states
// (validated against Pythia 6)

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qqbar2gluinogluino::initProc() {

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(1000021, 1000021);
  
} 

//--------------------------------------------------------------------------

// Begin evaluate d(sigmaHat)/d(tHat); flavour-independent part. 

void Sigma2qqbar2gluinogluino::sigmaKin() { 

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  // tHG = tHat - m_gluino^2; uHG = uHat - m_gluino^2. 
  // (Note: tHG and uHG defined with opposite sign wrt Pythia 6)
  s34Avg = 0.5 * (s3 + s4) - 0.25 * pow2(s3 - s4) / sH; 
  tHG    = -0.5 * (sH - tH + uH);
  uHG    = -0.5 * (sH + tH - uH); 
  tHG2   = tHG * tHG;
  uHG2   = uHG * uHG;

  // s-channel contribution.
  sigS   = (tHG2 + uHG2 + 2. * s34Avg * sH) / sH2; 

}

//--------------------------------------------------------------------------


// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qqbar2gluinogluino::sigmaHat() {  

  // Squarks (L/R or 1/2) can contribute in t or u channel.
  // Assume identical CKM matrices in quark and squark sector. 
  // (Note: tHQL,R and uHQL,R defined with opposite sign wrt Pythia 6. 
  //  This affects the sign of the interference term below)
  double sQL    = pow2( particleDataPtr->m0(1000000 + abs(id1)) );
  double tHQL   = tHG + s34Avg - sQL; 
  double uHQL   = uHG + s34Avg - sQL; 
  double sQR    = pow2( particleDataPtr->m0(2000000 + abs(id1)) );
  double tHQR   = tHG + s34Avg - sQR; 
  double uHQR   = uHG + s34Avg - sQR; 
 
  // Calculate kinematics dependence.
  double sigQL  = (4./9.) * (tHG2 / pow2(tHQL) + uHG2 / pow2(uHQL)) 
                + (1./9.) * s34Avg * sH / (tHQL * uHQL)
                + (tHG2 + sH * s34Avg) /(sH * tHQL)   
                + (uHG2 + sH * s34Avg) /(sH * uHQL);   
  double sigQR  = (4./9.) * (tHG2 / pow2(tHQR) + uHG2 / pow2(uHQR)) 
                + (1./9.) * s34Avg * sH / (tHQR * uHQR)
                + (tHG2 + sH * s34Avg) /(sH * tHQR)   
                + (uHG2 + sH * s34Avg) /(sH * uHQR);
  double sigSum = sigS + 0.5 * (sigQL + sigQR);     
    
  // Answer contains factor 1/2 from identical gluinos.
  double sigma  = (M_PI / sH2) * pow2(alpS) * (8./3.) * 0.5 * sigSum 
                * openFracPair;  
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2gluinogluino::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, 1000021, 1000021);

  // Two colour flow topologies. Swap if first is antiquark.
  if (rndmPtr->flat() < 0.5) setColAcol( 1, 0, 0, 2, 1, 3, 3, 2);
  else                    setColAcol( 1, 0, 0, 2, 3, 2, 1, 3); 
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma1qq2antisquark
// R-parity violating squark production

//--------------------------------------------------------------------------

// Initialise process

void Sigma1qq2antisquark::initProc(){

  //Typecast to the correct couplings
  coupSUSYPtr = (CoupSUSY*) couplingsPtr;

  //Construct name of the process from lambda'' couplings

  nameSave = "q q' -> " + coupSUSYPtr->getName(idRes)+"* + c.c";
  codeSave = 2000 + 10*abs(idRes)/1000000 + abs(idRes)%10;
}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma1qq2antisquark::sigmaKin() {

  // Check if at least one RPV coupling non-zero
  if(!coupSUSYPtr->isUDD) {
    sigBW = 0.0;
    return;
  }

  mRes = particleDataPtr->m0(abs(idRes));
  GammaRes = particleDataPtr->mWidth(abs(idRes));
  m2Res = pow2(mRes);
    
  sigBW        = sH * GammaRes/ ( pow2(sH - m2Res) + pow2(mRes * GammaRes) );
  sigBW       *= 2.0/3.0/mRes;

  // Width out only includes open channels. 
  widthOut     = GammaRes * particleDataPtr->resOpenFrac(id3);
}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma1qq2antisquark::sigmaHat() {

  // Only allow (anti)quark-(anti)quark incoming states
  if (id1*id2 <= 0) return 0.0;    

  // Generation indices
  int iA = (abs(id1)+1)/2;
  int iB = (abs(id2)+1)/2;

  //Covert from pdg-code to ~u_i/~d_i basis
  bool idown = (abs(idRes)%2 == 1) ? true : false;
  int iC = (abs(idRes)/1000000 == 2) 
         ? (abs(idRes)%10+1)/2 + 3: (abs(idRes)%10+1)/2; 

  // UDD structure
  if (abs(id1)%2 == 0 && abs(id2)%2 == 0) return 0.0;  
  if (abs(id1)%2 == 1 && abs(id2)%2 == 1 && idown) return 0.0;
  if ((abs(id1) + abs(id2))%2 == 1 && !idown) return 0.0;

  double sigma = 0.0;

  if(!idown){
   // d_i d_j --> ~u*_k
    for(int isq = 1; isq <=3; isq++){
      // Loop over R-type squark contributions
      sigma += pow2(coupSUSYPtr->rvUDD[isq][iA][iB]) 
	* norm(coupSUSYPtr->Rusq[iC][isq+3]);
    }
  }else{
    // u_i d_j --> d*_k
    // Pick the correct coupling for d-u in-state
    if(abs(id1)%2==1){
      iA = (abs(id2)+1)/2;
      iB = (abs(id1)+1)/2;
    }
    for(int isq = 1; isq <= 3; isq++){
      // Loop over R-type squark contributions
      sigma += pow2(coupSUSYPtr->rvUDD[iA][iB][isq]) 
	* norm(coupSUSYPtr->Rdsq[iC][isq+3]);
    }
  }

  sigma *= sigBW;
  return sigma;    

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1qq2antisquark::setIdColAcol() {

  // Set flavours.
  if(id1 < 0 && id2 < 0 ) setId( id1, id2, idRes);
  else setId( id1, id2, -idRes);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 2, 0, 0, 3);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}


//==========================================================================

} // end namespace Pythia8

