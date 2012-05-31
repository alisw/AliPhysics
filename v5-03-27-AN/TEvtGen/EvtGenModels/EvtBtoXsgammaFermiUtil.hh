//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2001      Brunel University
//
// Module: EvtGen/EvtBtoXsgammaFermiUtil.hh
//
// Description:
//       Class to hold various fermi functions and their helper functions. The 
//       fermi functions are used in EvtBtoXsgammaKagan.
//
// Modification history:
//
//       Jane Tinslay     March 21, 2001         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOXSGAMMAFERMIUTIL_HH
#define EVTBTOXSGAMMAFERMIUTIL_HH

#include <vector>

class EvtBtoXsgammaFermiUtil {

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  EvtBtoXsgammaFermiUtil() { };
  virtual ~EvtBtoXsgammaFermiUtil() { };

  //Exponential function
  static double FermiExpFunc(double var, const std::vector<double> &coeffs);

  //Gaussian function and its helper functions
  static double FermiGaussFunc(double, std::vector<double> const &coeffs);
  static double FermiGaussFuncRoot(double, double, double, std::vector<double> &coeffs);
  static double FermiGaussRootFcnA(double, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2);
  static double FermiGaussRootFcnB(double, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2);
  static double Gamma(double, const std::vector<double> &coeffs);

  //Roman function and its helper functions
  static double BesselI1(double);
  static double BesselK1(double);
  static double FermiRomanFuncRoot(double, double);
  static double FermiRomanRootFcnA(double);
  static double FermiRomanFunc(double, std::vector<double> const &coeffs);

};

#endif
