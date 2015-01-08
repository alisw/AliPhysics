//////////////////////////////////////////////////////////////////////
//
// Module: EvtVubAC.hh
//
//////////////////////////////////////////////////////////////////

#ifndef EVTVUBAC_HH
#define EVTVUBAC_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"
#include <vector>

class EvtParticle;

class EvtVubAC:public  EvtDecayIncoherent  {

public:
  
  EvtVubAC() {}
  virtual ~EvtVubAC();

  std::string getName();

  EvtDecayBase* clone();

  void initProbMax();

  void init();

  void decay(EvtParticle *Bmeson); 

private:
  // Input parameters
  double mB;
  double lambda2; 
  
  double alphaSmZ;
  double alphaSmB;
  double c;
  double q;
  double k;
  
  double CF;
  double CA;

  double beta0;
  
  std::vector<double> gvars;

  double rate(double u, double w, double xb);
  double wreg(double w);
  double alphaS(double Q);
  double PolyLog(double v, double z);
  double ureg(double u);
  double ularge(double u);
  double Coeff(double u, double w, double xb);
  double Coeff1(double w, double xb);
  double Coeff0(double w, double xb);
  double Sigma(double x1, double x2);
  double max(double ub, double lb);
  double d1(double u, double w, double xb);
  double d(double u, double w, double xb);
  double f(double w);
  double Lambda2(double x, double alphaSmZ);
  int Bisect(double x1, double x2,double precision,double& root,const double alphaSmZ);
  double FindRoot(const double alphaSmZ);

};

#endif


