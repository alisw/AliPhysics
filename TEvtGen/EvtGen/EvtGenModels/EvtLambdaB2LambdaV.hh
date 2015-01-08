//------------------------------------------------------------------------
//
// Module: EvtGen/EvtLambdaB2LambdaV.hh
//
// Description:
//   Class to generate LambdaB -> Lambda(p pi) V(Vp Vm) decays
//   with V a vector meson such as J/psi (mu+mu-)
//                                 Rho (pi+pi-)
//                                 Omega (pi+pi-)
//                                 Rho-omega mixing (pi+pi-)
//
// Author : Eric Conte (LPC Clermont-Ferrand)
//          econte@clermont.in2p3.fr / ziad@clermont.in2p3.fr
//
// Modification history:
//
//    E. Conte        April 13, 2006         Module created
//    E. Conte        February 5, 2006       First draft
//
//------------------------------------------------------------------------

#ifndef EVTLAMBDAB2LAMBDAV_HH
#define EVTLAMBDAB2LAMBDAV_HH

#include <stdlib.h>
#include <string>
#include "EvtGenBase/EvtDecayProb.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"

namespace VID
{
  enum VectorMesonType{JPSI, OMEGA, RHO, RHO_OMEGA_MIXING};
}

//*******************************************************************
//*                                                                 *
//*                   Class EvtLambdaB2LambdaV                      *
//*                                                                 *
//*******************************************************************
//
// DECAY : LambdaB -> Lambda + vector meson
//
// d(Sigma)
// -------- = 1 + A*B*cos(theta) + 2*A*Re(C*exp(i*phi))*sin(theta)
// d(Omega)
//
// with A (real)    : lambdaB  helicity asymmetry parameter
//      B (real)    : lambdaB polarisation
//      C (complex) : lambdaB density matrix element rho+-
//
// cf : O. Leitner, Z.J Ajaltouni, E. Conte, 
//      PCCF RI 0601, ECT-05-15, LPNHE/2006-01, hep-ph/0602043
 
class EvtLambdaB2LambdaV:public  EvtDecayProb
{

public:

  EvtLambdaB2LambdaV();
  virtual ~EvtLambdaB2LambdaV();
  EvtDecayBase* clone();

  virtual std::string getName();
  void init();
  void initProbMax();
  void decay(EvtParticle *lambdab);
 
private:

  //class name for report method
  std::string fname;

  //meson vector identity 
  VID::VectorMesonType Vtype;

  //decay dynamics parameters 
  double A;
  double B;
  EvtComplex C;

  //V mass generator method
  double getVMass(double MASS_LAMBDAB, double MASS_LAMBDA);
  
  //PDF generator method
  double BreitWignerRelPDF(double m,double _m0, double _g0);
  double RhoOmegaMixingPDF(double m, double _mr, double _gr, double _mo, double _go);
};





//*******************************************************************
//*                                                                 *
//*             Class EvtLambda2PPiForLambdaB2LambdaV               *
//*                                                                 *
//*******************************************************************
//
// DECAY : Lambda -> p + pi-
//
// d(Sigma)
// -------- = 1 + A*B*cos(theta) + 2*A*Re(D*exp(i*phi))*sin(theta)
// d(Omega)
//
// with A (real)    : lambda asymmetry parameter
//      B (real)    : lambda polarisation
//      C (real)    : lambdaB polarisation
//      D (complex) : lambda density matrix element rho+-
//
// cf : O. Leitner, Z.J Ajaltouni, E. Conte
//      PCCF RI 0601, ECT-05-15, LPNHE/2006-01, hep-ph/0602043

class EvtLambda2PPiForLambdaB2LambdaV:public  EvtDecayProb
{

public:

  EvtLambda2PPiForLambdaB2LambdaV();
  virtual ~EvtLambda2PPiForLambdaB2LambdaV();
  EvtDecayBase* clone();

  virtual std::string getName();
  void init();
  void initProbMax();
  void decay(EvtParticle *lambda);  

private :

  //class name for report method
  std::string fname;

  //meson vector identity 
  VID::VectorMesonType Vtype;

  //decay dynamics parameters
  double A;
  double B;
  double C;
  EvtComplex D;
};





//*******************************************************************
//*                                                                 *
//*               Class EvtV2VpVmForLambdaB2LambdaV                 *
//*                                                                 *
//*******************************************************************
//
// DECAY : vector meson V -> Vp + Vm
//
// d(Sigma)
// -------- = (1-3A)*cos(theta)^2 + (1+A)   //leptonic decays
// d(Omega)
//
// d(Sigma)
// -------- = (3A-1)*cos(theta)^2 + (1-A)   //hadronic decays
// d(Omega)
//
// with A (real)    : V density matrix element indicating the
//                    probability to be longitudinally polarized
//
// cf : O. Leitner, Z.J Ajaltouni, E. Conte
//      PCCF RI 0601, ECT-05-15, LPNHE/2006-01, hep-ph/0602043

class EvtV2VpVmForLambdaB2LambdaV:public  EvtDecayProb
{

public:

  EvtV2VpVmForLambdaB2LambdaV();
  virtual ~EvtV2VpVmForLambdaB2LambdaV();
  EvtDecayBase* clone();

  virtual std::string getName();
  void init();
  void initProbMax();
  void decay(EvtParticle *V);
  
private:

  //class name for report method
  std::string fname;

  //meson vector identity 
  VID::VectorMesonType Vtype;
  //decay dynamics parameters 
  double A;
};


#endif
