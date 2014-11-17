//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 2001      Brunel University, University of Wisconsin
//
// Module: EvtGen/EvtBtoXsgammaKagan.hh
//
// Description:
//       Implimentation of the Kagan-Neubert model for non-resonant 
//       B->Xs,gamma decays.
//
// Modification history:
//
//       Jane Tinslay, Francesca Di Lodovico     March 27, 2001  Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOXSGAMMAKAGAN_HH
#define EVTBTOXSGAMMAKAGAN_HH

#include <vector>
#include "EvtGenModels/EvtBtoXsgammaAbsModel.hh"

class EvtBtoXsgammaKagan : public EvtBtoXsgammaAbsModel {

public:
  
  EvtBtoXsgammaKagan() {}

  virtual ~EvtBtoXsgammaKagan();

  void init(int, double*);

  void computeHadronicMass(int, double*);
  
  void getDefaultHadronicMass();

  double GetMass(int code);

  double CalcAlphaS(double);

  void CalcWilsonCoeffs();
  void CalcDelta();
  double Fz(double);

private:

  //Input parameters
  double _mb;
  double _mB;
  double _delta;
  double _nIntervalS;
  double _nIntervalmH;
  double _lambdabar;
  double _lam1;
  double _mHmin;
  double _mHmax;  
  //Other parameters
  double _r7;
  double _gam77;
  double _gam27;
  double _gam87;
  double _beta0;
  double _beta1;
  double _alphasmZ;
  double _mZ;
  double _z;
  double _fz;
  double _lam2;
  double _kappabar;
  double _rer2;
  double _rer8;
  double _kSLemmu;
  double _mW;
  double _mt;
  double _ms;
  double _mu;

  double _c2mu;
  double _c70mu;
  double _c80mu;
  double _c71mu;
  double _c7emmu;

  double _cDeltatot;

  double _alpha;
  double _alphasmW;
  double _alphasmt;
  double _alphasmu;
  double _alphasmubar;
  double _etamu;

  std::vector<double> _mHVect;

  static double ReG(double);
  static double ImG(double);
  static double s77(double);
  static double s88(double, double, double);
  static double s78(double);
  static double s22Func(double var, const std::vector<double> &coeffs);
  static double s27Func(double var, const std::vector<double> &coeffs);

  static double Delta(double, double);
  static double DeltaFermiFunc(double, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2, const std::vector<double> &coeffs3);
  static double s77FermiFunc(double, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2);
  static double s88FermiFunc(double, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2, const std::vector<double> &coeffs3);
  static double s78FermiFunc(double, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2);
  static double s22FermiFunc(double, std::vector<double> &coeffs);
  static double s27FermiFunc(double, std::vector<double> &coeffs);
  static double s28FermiFunc(double, std::vector<double> &coeffs);
  static double GetArrayVal(double, double, double, double, std::vector<double>);
  static double sFermiFunc(double, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2, 
			   const std::vector<double> &coeffs3, const std::vector<double> &coeffs4);
  static double FermiFunc(double, const std::vector<double> &coeffs);
  static double diLogFunc(double);
  static double diLogMathematica(double);
  double *massHad; double *brHad;
  static double intervalMH;
  static bool bbprod;
};

#endif



