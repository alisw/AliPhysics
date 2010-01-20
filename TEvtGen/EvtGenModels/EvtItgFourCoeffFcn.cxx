//--------------------------------------------------------------------------
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtItgFourCoeffFcn.cc
//
// Description:
//      Class describing a function with four vectors of coefficients. 
//
// Modification history:
//
//    Jane Tinslay                March 21, 2001       Module created
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenModels/EvtItgFourCoeffFcn.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

EvtItgFourCoeffFcn::EvtItgFourCoeffFcn( double (*theFunction)(double, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, const std::vector<double> &), double lowerRange, double upperRange, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2, const std::vector<double> &coeffs3, const std::vector<double> &coeffs4):
  EvtItgAbsFunction(lowerRange, upperRange),
  _myFunction(theFunction),
  _coeffs1(coeffs1),
  _coeffs2(coeffs2), 
  _coeffs3(coeffs3),
  _coeffs4(coeffs4)
{}

EvtItgFourCoeffFcn::~EvtItgFourCoeffFcn( )
{}


double
EvtItgFourCoeffFcn::myFunction(double x) const{
  return _myFunction(x, _coeffs1, _coeffs2, _coeffs3, _coeffs4);
}

void
EvtItgFourCoeffFcn::setCoeff(int vect, int which, double value)
{
  if (vect == 1) _coeffs1[which] = value;
  else if (vect == 2) _coeffs2[which] = value;
  else if (vect == 3) _coeffs3[which] = value;
  else if (vect == 4) _coeffs4[which] = value;
}

double
EvtItgFourCoeffFcn::getCoeff(int vect, int which)
{
  if (vect == 1) return _coeffs1[which];
  else if (vect == 2) return _coeffs2[which];
  else if (vect == 3) return _coeffs3[which];
  else if (vect == 4) return _coeffs4[which];
  else {return 0;}
}
