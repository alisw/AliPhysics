
//////////////////////////////////////////////////////////////////////
//
// Module: EvtVubBLNPHybrid.hh
//
//
// tripleDiff from BLNPHybrid's notebook, based on hep-ph/0504071
// 
// Description:
// Class to generate inclusive B to X_u l nu decays.
// This class is based on EvtVubBLNP by Sheila Mclachlin with an update to
// generate the inclusive decays in such a way that the right
// mix of inclusive and exclusive decays is obtained.
// 
//NOTE:
// - A set of weights (for bins in the kinematic variables mX, q2, El)
//   is read from DECAY.DEC. This set of weights must be consistent
//   with the other parameters specified (excl. BF, non-res BF, mb, a).
// - If no binning/weights are specified in DECAY.DEC the hybrid
//   reweighting is not activated
//
//
//////////////////////////////////////////////////////////////////


#ifndef EVTVUBBLNPHYBRID_HH
#define EVTVUBBLNPHYBRID_HH

#include <vector>
#include "EvtGenBase/EvtDecayIncoherent.hh"

class EvtParticle;

class EvtVubBLNPHybrid:public  EvtDecayIncoherent  {

public:
  
  EvtVubBLNPHybrid(); 
  virtual ~EvtVubBLNPHybrid();
  
  std::string getName();
  
  EvtDecayBase* clone();

  void initProbMax();

  void init();

  void decay(EvtParticle *Bmeson); 

  void readWeights(int startArg=0);
  
  double getWeight(double mX, double q2, double El);

private:
  //Input for hybrid modell
  enum { nParameters = 10, nVariables = 3 };
  bool _noHybrid;
  bool _storeWhat;
  int    _nbins_mX;
  int    _nbins_q2;
  int    _nbins_El;
  int    _nbins;
  double _masscut;
  double * _bins_mX;
  double * _bins_q2;
  double * _bins_El;
  double * _weights;
  
  // Input parameters
  double mBB;
  double lambda2;

  // Shape function parameters
  double b;
  double Lambda;
  double Ecut;
  double wzero;

  // SF and SSF modes
  int itype;
  double dtype;
  int isubl;

  // flags
  int flag1;
  int flag2;
  int flag3;

  // Quark mass
  double mb;

  // Matching scales
  double muh;
  double mui;
  double mubar;

  // Perturbative quantities
  double CF;
  double CA;

  double beta0;
  double beta1;
  double beta2;

  double zeta3;

  double Gamma0;
  double Gamma1;
  double Gamma2;

  double gp0;
  double gp1;

  double Lbar;
  double mupisq;
  double moment2;

  int flagpower;
  int flag2loop;

  int maxLoop;
  double precision;

  std::vector<double> gvars;
  
  double rate3(double Pp, double Pl, double Pm);
  double F1(double Pp, double Pm, double muh, double mui, double mubar, double doneJS, double done1);
  double F2(double Pp, double Pm, double muh, double mui, double mubar, double done3);
  double F3(double Pp, double Pm, double muh, double mui, double mubar, double done2);
  double DoneJS(double Pp, double Pm, double mui);
  double Done1(double Pp, double Pm, double mui);
  double Done2(double Pp, double Pm, double mui);
  double Done3(double Pp, double Pm, double mui);
  static double IntJS(double what, const std::vector<double> &vars);
  static double Int1(double what, const std::vector<double> &vars);
  static double Int2(double what, const std::vector<double> &vars);
  static double Int3(double what, const std::vector<double> &vars);
  static double g1(double w, const std::vector<double> &vars);
  static double g2(double w, const std::vector<double> &vars);
  static double g3(double w, const std::vector<double> &vars);
  static double Shat(double w, const std::vector<double> &vars);
  static double Mzero(double muf, double mu, double mupisq, const std::vector<double> &vars);
  double wS(double w);
  double t(double w);
  double u(double w);
  double v(double w);
  double myfunction(double w, double Lbar, double mom2);
  double myfunctionBIK(double w, double Lbar, double mom2);
  double dU1nlo(double muh, double mui);
  double U1lo(double muh, double mui);
  double Sfun(double mu1, double mu2, double epsilon);
  double S0(double a1, double r);
  double S1(double a1, double r);
  double S2(double a1, double r);
  double aGamma(double mu1, double mu2, double epsilon);
  double agp(double mu1, double mu2, double epsilon);
  double alo(double muh, double mui);
  double anlo(double muh, double mui);   // d/depsilon of aGamma
  static double alphas(double mu, const std::vector<double> &vars);
  double PolyLog(double v, double z);
  static double Gamma(double z);
  static double Gamma(double a, double x);
  static double gamser(double a, double x, double LogGamma);
  static double gammcf(double a, double x, double LogGamma);
  double findBLNPWhat();
  std::vector<double> _pf;
};

#endif

