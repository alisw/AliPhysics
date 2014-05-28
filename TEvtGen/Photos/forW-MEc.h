#ifndef _forW_MEc_h_included_
#define _forW_MEc_h_included_
#include <complex>
using std::complex;

namespace Photospp
{

class PhotosMEforW
{
public:
  static void   PHOBWnlo(double *WT);

private:
  static double WDecayEikonalSqrKS_1ph(double p3[4],double p1[4],double p2[4],double k[4]);
  static double WDecayBornAmpSqrKS_1ph(double p3[4],double p1[4],double p2[4]);
  static double WDecayAmplitudeSqrKS_1ph(double p3[4],double p1[4],double p2[4],double k[4]);
  static double SANC_WT(double PW[4],double PNE[4],double PMU[4],double PPHOT[4],double B_PW[4],double B_PNE[4],double B_PMU[4]);
  static void   SANC_INIT1(double QB0,double QF20,double MF10,double MF20,double MB0);
  static void   SANC_INIT(double ALPHA,int PHLUN);

private:
  static complex<double> InProd_zero(double p1[4],int l1,double p2[4],int l2);
  static double          InSqrt(double p[4],double q[4]);
  static complex<double> InProd_mass(double p1[4],double m1,int l1,double p2[4],double m2,int l2);
  static complex<double> BsFactor(int s,double k[4],double p[4],double m);
  static complex<double> WDecayEikonalKS_1ph(double p3[4],double p1[4],double p2[4],double k[4],int s);
  static complex<double> SoftFactor(int s,double k[4],double p1[4],double m1,double p2[4],double m2,double Gmass2);
  static complex<double> TrMatrix_zero(double p1[4],double m1,int l1,double k[4],int s,double p2[4],double m2,int l2);
  static complex<double> TrMatrix_mass(double p1[4],double m1,int l1,double k[4],double m,int s,double p2[4],double m2,int l2);
  static complex<double> WDecayBornAmpKS_1ph(double p3[4],int l3,double p1[4],int l1,double p2[4],int l2);
  static complex<double> WDecayAmplitudeKS_1ph(double p3[4],int l3,double p1[4],int l1,double p2[4],int l2,double k[4],int s);

private:
  // COMMON /Kleiss_Stirling/spV,bet
  static double spV[4],bet[4];

  // COMMON /mc_parameters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN
  static double pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN;
};

} // namespace Photospp
#endif
