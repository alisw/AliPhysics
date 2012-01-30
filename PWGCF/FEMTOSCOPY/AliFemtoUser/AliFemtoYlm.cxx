////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoYlm - the class to calculate varous components of spherical        //
//  harmonics                                                                 //
//                                                                            //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoYlm.h"
#include <TMath.h>
#include <iostream>

double *AliFemtoYlm::fgPrefactors = 0x0;
int    *AliFemtoYlm::fgPrefshift = 0x0;
int    *AliFemtoYlm::fgPlmshift = 0x0;

AliFemtoYlm::AliFemtoYlm() {
  InitializeYlms();
}

AliFemtoYlm::~AliFemtoYlm() {}


AliFemtoYlm::AliFemtoYlm(const AliFemtoYlm& aYlm){
  fgPrefshift = aYlm.fgPrefshift;
  InitializeYlms();
}

AliFemtoYlm& AliFemtoYlm::operator=(const AliFemtoYlm& aYlm){
  if (this == &aYlm)
    return *this;

  InitializeYlms();

  return *this;
}

std::complex<double> AliFemtoYlm::Ceiphi(double phi){
  return std::complex<double>(cos(phi),sin(phi));
}

double AliFemtoYlm::Legendre(int ell, int em, double ctheta){
  // Calculate a single Legendre value
  // *** Warning - NOT optimal - calculated all Plms up to L !!!
  double lbuf[36];
  AliFemtoYlm::LegendreUpToYlm(ell, ctheta, lbuf);

  return lbuf[fgPlmshift[ell]-abs(em)];
}

std::complex<double> AliFemtoYlm::Ylm(int ell,int m,double theta,double phi){
  // Calculate Ylm spherical input
  double ctheta;
  std::complex<double> answer;
  std::complex<double> ci(0.0,1.0);
  ctheta=cos(theta);
  answer = (fgPrefactors[fgPrefshift[ell]+m]*Legendre(ell,m,ctheta))*Ceiphi(m*phi);

  return answer;
}

std::complex<double> AliFemtoYlm::Ylm(int ell, int m, double x, double y, double z){
  // Calculate Ylm cartesian input
  std::complex<double> answer; 
  double ctheta,phi;
  double r = sqrt(x*x+y*y+z*z);
  if ( r < 1e-10 || fabs(z) < 1e-10 ) ctheta = 0.0;
  else ctheta=z/r;
  phi=atan2(y,x);
  answer = (fgPrefactors[fgPrefshift[ell]+m]*Legendre(ell,m,ctheta))*Ceiphi(m*phi);

  return answer;	
}

double AliFemtoYlm::ReYlm(int ell, int m, double theta, double phi){
	return real(AliFemtoYlm::Ylm(ell,m,theta,phi));
}

double AliFemtoYlm::ImYlm(int ell, int m, double theta, double phi){
	return imag(AliFemtoYlm::Ylm(ell,m,theta,phi));
}

double AliFemtoYlm::ReYlm(int ell, int m, double x,double y,double z){
	return real(AliFemtoYlm::Ylm(ell,m,x,y,z));
}

double AliFemtoYlm::ImYlm(int ell, int m, double x,double y,double z){
	return imag(AliFemtoYlm::Ylm(ell,m,x,y,z));
}

void AliFemtoYlm::InitializeYlms()
{
  // Calculate prefactors for fast Ylm calculation

  double oneoversqrtpi = 1.0/TMath::Sqrt(TMath::Pi());

  fgPrefactors = (double *) malloc(sizeof(double) * 36);
  fgPrefshift  = (int *) malloc(sizeof(int) * 6);
  fgPlmshift   = (int *) malloc(sizeof(int) * 6);

  // l=0 prefactors
  fgPrefactors[0] = 0.5*oneoversqrtpi;

  // l=1 prefactors
  fgPrefactors[1] = 0.5*sqrt(3.0/2.0)*oneoversqrtpi;
  fgPrefactors[2] = 0.5*sqrt(3.0)*oneoversqrtpi;
  fgPrefactors[3] = -fgPrefactors[1];

  // l=2 prefactors
  fgPrefactors[4] = 0.25*sqrt(15.0/2.0)*oneoversqrtpi;
  fgPrefactors[5] = 0.5*sqrt(15.0/2.0)*oneoversqrtpi;
  fgPrefactors[6] = 0.25*sqrt(5.0)*oneoversqrtpi;
  fgPrefactors[7] = -fgPrefactors[5];
  fgPrefactors[8] = fgPrefactors[4];

  // l=3 prefactors
  fgPrefactors[9]  = 0.125*sqrt(35.0)*oneoversqrtpi;
  fgPrefactors[10] = 0.25*sqrt(105.0/2.0)*oneoversqrtpi;
  fgPrefactors[11] = 0.125*sqrt(21.0)*oneoversqrtpi;
  fgPrefactors[12] = 0.25*sqrt(7.0)*oneoversqrtpi;
  fgPrefactors[13] = -fgPrefactors[11];
  fgPrefactors[14] = fgPrefactors[10];
  fgPrefactors[15] = -fgPrefactors[9];

  // l=4 prefactors
  fgPrefactors[16] = 3.0/16.0*sqrt(35.0/2.0)*oneoversqrtpi;
  fgPrefactors[17] = 3.0/8.0*sqrt(35.0)*oneoversqrtpi;
  fgPrefactors[18] = 3.0/8.0*sqrt(5.0/2.0)*oneoversqrtpi;
  fgPrefactors[19] = 3.0/8.0*sqrt(5.0)*oneoversqrtpi;
  fgPrefactors[20] = 3.0/16.0*oneoversqrtpi;
  fgPrefactors[21] = -fgPrefactors[19];
  fgPrefactors[22] = fgPrefactors[18];
  fgPrefactors[23] = -fgPrefactors[17];
  fgPrefactors[24] = fgPrefactors[16];

  // l=5 prefactors
  fgPrefactors[25] = 3.0/32.0*sqrt(77.0)*oneoversqrtpi;
  fgPrefactors[26] = 3.0/16.0*sqrt(385.0/2.0)*oneoversqrtpi;
  fgPrefactors[27] = 1.0/32.0*sqrt(385.0)*oneoversqrtpi;
  fgPrefactors[28] = 1.0/8.0*sqrt(1155.0/2.0)*oneoversqrtpi;
  fgPrefactors[29] = 1.0/16.0*sqrt(165.0/2.0)*oneoversqrtpi;
  fgPrefactors[30] = 1.0/16.0*sqrt(11.0)*oneoversqrtpi;
  fgPrefactors[31] = -fgPrefactors[29];
  fgPrefactors[32] = fgPrefactors[28];
  fgPrefactors[33] = -fgPrefactors[27];
  fgPrefactors[34] = fgPrefactors[26];
  fgPrefactors[35] = -fgPrefactors[25];

  fgPrefshift[0] = 0;
  fgPrefshift[1] = 2;
  fgPrefshift[2] = 6;
  fgPrefshift[3] = 12;
  fgPrefshift[4] = 20;
  fgPrefshift[5] = 30;

  fgPlmshift[0] = 0;
  fgPlmshift[1] = 2;
  fgPlmshift[2] = 5;
  fgPlmshift[3] = 9;
  fgPlmshift[4] = 14;
  fgPlmshift[5] = 20;
}

void AliFemtoYlm::LegendreUpToYlm(int lmax, double ctheta, double *lbuf)
{
  // Calculate a set of legendre polynomials up to a given l
  // with spherical input
  double sins[6];
  double coss[6];
  sins[0] = 0.0;
  coss[0] = 1.0;
  sins[1] = sqrt(1-ctheta*ctheta);
  coss[1] = ctheta;
  for (int iter=2; iter<6; iter++) {
    sins[iter] = sins[iter-1]*sins[1];
    coss[iter] = coss[iter-1]*coss[1];
  }

  // Legendre polynomials l=0
  lbuf[0] = 1.0;

  // Legendre polynomials l=1
  if (lmax>0) {
    lbuf[1] = sins[1];
    lbuf[2] = coss[1];
  }

  // Legendre polynomials l=2
  if (lmax>1) {
    lbuf[3] = sins[2];
    lbuf[4] = sins[1]*coss[1];
    lbuf[5] = 3*coss[2]-1;
  }

  // Legendre polynomials l=3
  if (lmax>2) {
    lbuf[6] = sins[3];
    lbuf[7] = sins[2]*coss[1];
    lbuf[8] = (5*coss[2]-1)*sins[1];
    lbuf[9] = 5*coss[3]-3*coss[1];
  }

  // Legendre polynomials l=4
  if (lmax>3) {
    lbuf[10] = sins[4];
    lbuf[11] = sins[3]*coss[1];
    lbuf[12] = (7*coss[2]-1)*sins[2];
    lbuf[13] = (7*coss[3]-3*coss[1])*sins[1];
    lbuf[14] = 35*coss[4]-30*coss[2]+3;
  }

  // Legendre polynomials l=5
  if (lmax>4) {
    lbuf[15] = sins[5];
    lbuf[16] = sins[4]*coss[1];
    lbuf[17] = (9*coss[2]-1)*sins[3];
    lbuf[18] = (3*coss[3]-1*coss[1])*sins[2];
    lbuf[19] = (21*coss[4]-14*coss[2]+1)*sins[1];
    lbuf[20] = 63*coss[5]-70*coss[3]+15*coss[1];
  }
}

void AliFemtoYlm::YlmUpToL(int lmax, double x, double y, double z, std::complex<double> *ylms)
{
  // Calculate a set of Ylms up to a given l
  // with cartesian input
  double ctheta,phi;

  double r = sqrt(x*x+y*y+z*z);
  if ( r < 1e-10 || fabs(z) < 1e-10 ) ctheta = 0.0;
  else ctheta=z/r;
  phi=atan2(y,x);
  
  YlmUpToL(lmax, ctheta, phi, ylms);

}

void AliFemtoYlm::YlmUpToL(int lmax, double ctheta, double phi, std::complex<double> *ylms)
{
  // Calculate a set of Ylms up to a given l
  // with spherical input
  int lcur = 0;  
  double lpol;
  
  double coss[6];
  double sins[6];

  double lbuf[36];
  LegendreUpToYlm(lmax, ctheta, lbuf);

  for (int iter=1; iter<=lmax; iter++) {
    coss[iter-1] = cos(iter*phi);
    sins[iter-1] = sin(iter*phi);
  }
  ylms[lcur++] = fgPrefactors[0]*lbuf[0] * std::complex<double>(1,0);
  
  for (int il = 1; il<=lmax; il++) {
    // First im = 0
    ylms[lcur+il] = fgPrefactors[fgPrefshift[il]]*lbuf[fgPlmshift[il]]*std::complex<double>(1.0,0.0);
    // Im != 0
    for (int im=1; im<=il; im++) {
      lpol = lbuf[fgPlmshift[il]-im];
      ylms[lcur+il-im] = fgPrefactors[fgPrefshift[il]-im]*lpol*std::complex<double>(coss[im-1],-sins[im-1]);
      ylms[lcur+il+im] = fgPrefactors[fgPrefshift[il]+im]*lpol*std::complex<double>(coss[im-1],sins[im-1]);
    }
    lcur += 2*il + 1;
  }
}
