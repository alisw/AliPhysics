////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoYlm - the class to calculate varous components of spherical        //
//  harmonics                                                                 //
//                                                                            //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOYLM_H
#define ALIFEMTOYLM_H
#include <cstdlib>
#include <cmath>
#include <complex>
#include <TMath.h>

class AliFemtoYlm {
 public:
  AliFemtoYlm();
  ~AliFemtoYlm();

  AliFemtoYlm(const AliFemtoYlm& aYlm);
  AliFemtoYlm& operator=(const AliFemtoYlm& aYlm);

  static double Legendre(int ell, int emm, double ctheta);
  static void   LegendreUpToYlm(int lmax, double ctheta, double *lbuf);

  static std::complex<double> Ylm(int ell,int m,double theta,double phi);
  static std::complex<double> Ylm(int ell, int m, double x, double y, double z);

  static void YlmUpToL(int lmax, double x, double y, double z, std::complex<double> *ylms);
  static void YlmUpToL(int lmax, double ctheta, double phi, std::complex<double> *ylms);

  static double ReYlm(int ell, int m, double theta, double phi);
  static double ReYlm(int ell, int m, double x, double y, double z);
  static double ImYlm(int ell, int m, double theta, double phi);
  static double ImYlm(int ell, int m, double x, double y, double z);

  static void InitializeYlms();
  
 private:
  static std::complex<double> Ceiphi(double phi);

  static double *fgPrefactors;
  static int    *fgPrefshift;
  static int    *fgPlmshift;
};

#endif
