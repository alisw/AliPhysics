//---------------------------------------------------------------------------------
//
// Wilson coeficients according to A.J.Buras and M.Munz, Phys.Rev. D52, 186. (1995)
// Thanks to N. Nikitine for example code for Pythia
// Coefficient C8eff and C2 correction to C7eff taken from:
//   A.J.Buras, M.Misiak, M.Munz, S.Pokorski, Nucl.Phys. B424, 374 (1994)
//
// Used constants come from PDG 2004
//
// P. Reznicek 18.02.2005
//
//  04/03/2005  PR   Added h-function                                                                                                    
//
//---------------------------------------------------------------------------------

#ifndef EVTWILSONCOEFICIENTS_HH
#define EVTWILSONCOEFICIENTS_HH

#include "EvtGenBase/EvtComplex.hh"

class EvtWilsonCoefficients {

public:

  EvtWilsonCoefficients();
  //~EvtWilsonCoefficients() {};

  // calculate strong coupling constant for n_f flavours and scale mu
  double alphaS(double mu,int n_f,double Lambda);
  // calculate Lambda matching alphaS using simple iterative method
  double Lambda(double alpha,int n_f,double mu,double epsilon,int maxstep);
  // eta-function: ratio of strong coupling constants
  double eta(double mu,int n_f,double Lambda,double M_W);

  // Wilson coeficients C1-C6
  EvtComplex C1(double mu,int n_f,double Lambda,double M_W);
  EvtComplex C2(double mu,int n_f,double Lambda,double M_W);
  EvtComplex C3(double mu,int n_f,double Lambda,double M_W);
  EvtComplex C4(double mu,int n_f,double Lambda,double M_W);
  EvtComplex C5(double mu,int n_f,double Lambda,double M_W);
  EvtComplex C6(double mu,int n_f,double Lambda,double M_W);
  // Wilson coeficietns C7,C8 => C7eff
  EvtComplex C7(double M_t,double M_W);
  EvtComplex C8(double M_t,double M_W);
  EvtComplex C7eff0(double mu,int n_f,double Lambda,double M_t,double M_W);
  EvtComplex C8eff0(double mu,int n_f,double Lambda,double M_t,double M_W);
  // Wilson coeficient C10
  EvtComplex C10tilda(double sin2W,double M_t,double M_W);
  EvtComplex C10(double sin2W,double M_t,double M_W,double ialpha);
  // Wilson coeficient C9
  double PE(double mu,int n_f,double Lambda,double M_W);
  EvtComplex P0(int ksi,double mu,int n_f,double Lambda,double M_W);
  EvtComplex C9tilda(int ksi,double mu,int n_f,double Lambda,double sin2W,double M_t,double M_W);
  EvtComplex C9(int ksi,double mu,int n_f,double Lambda,double sin2W,double M_t,double M_W,double ialpha);

  // Intermediate functions A-F,Y,Z
  double A(double x);
  double B(double x);
  double C(double x);
  double D(double x);
  double E(double x);
  double F(double x);
  double Y(double x);
  double Z(double x);

  // Mode decay specific functions
  EvtComplex hzs(double z,double shat,double mu,double M_b);
  double     fz(double z);
  double     kappa(double z,double alpha_S);
  double     etatilda(double shat,double alpha_S);
  double     omega(double shat);
  EvtComplex C9efftilda(double z,double shat,double alpha_S,EvtComplex c1,EvtComplex c2,EvtComplex c3,EvtComplex c4,EvtComplex c5,EvtComplex c6,EvtComplex c9tilda,int ksi);
  EvtComplex C7b2sg(double alpha_S,double et,EvtComplex c2,double M_t,double M_W);
  EvtComplex Yld(double q2,double ki[],double Gi[],double Mi[],int ni,EvtComplex c1,EvtComplex c2,EvtComplex c3,EvtComplex c4,EvtComplex c5,EvtComplex c6,double ialpha);

  // User function
  void CalculateAllCoefficients();
  // Set parameters
  void SetLambda(double lambda) { m_Lambda=lambda; }
  void CalculateLambda(double epsilon,int maxstep) { m_Lambda=Lambda(m_alphaMZ,m_n_f,m_mu,epsilon,maxstep); }
  void SetStrongCouplingAtZMass(double alphaMZ) { m_alphaMZ=alphaMZ; }
  void SetScale(double mu) { m_mu=mu; }
  void SetNumberOfFlavours(int n_f) { m_n_f=n_f; }
  void SetZMass(double M_Z) { m_M_Z=M_Z; }
  void SetWMass(double M_W) { m_M_W=M_W; }
  void SetTopMass(double M_t) { m_M_t=M_t; }
  void SetSin2WeinbergAngle(double sin2W) { m_sin2W=sin2W;}
  void SetInvElMagCoupling(double ialpha) { m_ialpha=ialpha; }
  void SetRenormalizationScheme(std::string scheme);
  // Get parameters
  double GetLambda() { return m_Lambda; }
  double GetStrongCouplingAtZMass() { return m_alphaMZ; }
  double GetStrongCouplingConst() { return m_alphaS; }
  double GetScale() { return m_mu; }
  int GetNumberOfFlavours() { return m_n_f; }
  int GetRenormSchemePar() { return m_ksi; }
  double GetZMass() { return m_M_Z; }
  double GetWMass() { return m_M_W; }
  double GetTopMass() { return m_M_t; }
  double GetSin2WeinbergAngle() { return m_sin2W;}
  double GetInvElMagCoupling() { return m_ialpha; }
  double GetEta() { return m_eta; }
  // Get results
  double GetA() { return m_A; }
  double GetB() { return m_B; }
  double GetC() { return m_C; }
  double GetD() { return m_D; }
  double GetE() { return m_E; }
  double GetF() { return m_F; }
  double GetY() { return m_Y; }
  double GetZ() { return m_Z; }
  EvtComplex GetC1() { return m_C1; }
  EvtComplex GetC2() { return m_C2; }
  EvtComplex GetC3() { return m_C3; }
  EvtComplex GetC4() { return m_C4; }
  EvtComplex GetC5() { return m_C5; }
  EvtComplex GetC6() { return m_C6; }
  EvtComplex GetC7() { return m_C7; }
  EvtComplex GetC8() { return m_C8; }
  EvtComplex GetC9() { return m_C9; }
  EvtComplex GetC10() { return m_C10; }
  EvtComplex GetC7eff0() { return m_C7eff0; }
  EvtComplex GetC8eff0() { return m_C8eff0; }
  EvtComplex GetC9tilda() { return m_C9tilda; }
  EvtComplex GetC10tilda() { return m_C10tilda; }
  EvtComplex GetP0() { return m_P0; }
  double GetPE() { return m_PE; }

private:

  int m_n_f,m_ksi;
  double m_Lambda,m_alphaMZ,m_mu,m_M_Z,m_M_t,m_M_W,m_alphaS,m_eta,m_sin2W,m_ialpha;
  EvtComplex m_C1,m_C2,m_C3,m_C4,m_C5,m_C6,m_C7,m_C7eff0,m_C8,m_C8eff0,m_C9,m_C9tilda,m_C10,m_C10tilda,m_P0;
  double m_A,m_B,m_C,m_D,m_E,m_F,m_Y,m_Z,m_PE;

  double k[6][8],a[8],h[8],p[8],r[2][8],s[8],q[8],g[8];

};

#endif
