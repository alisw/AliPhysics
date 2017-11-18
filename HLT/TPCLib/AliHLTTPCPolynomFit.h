//-*- Mode: C++ -*-
//*************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************


#ifndef AliHLTTPCPolynomFit_H
#define AliHLTTPCPolynomFit_H


/**
 * The class AliHLTTPCPolynomFit allows one to fit polynomial coefficients c0..cN
 *
 * with given measurements m_i, i=0, ... :
 *
 * m_i = c0*f0_i + c1*f1_i + .. + cN*fN_i + measurement_error
 *
 */

class AliHLTTPCPolynomFit
{
public:

  AliHLTTPCPolynomFit(): fN(0), fA(0), fB(0)
  { }
  
  AliHLTTPCPolynomFit( int nCoefficients ): fN(0), fA(0), fB(0)
  {
    Reset(nCoefficients);
  }
  
  ~AliHLTTPCPolynomFit(){
    delete[] fA;
    delete[] fB;
  }
  
  static int invS( long double M[], int N );

  void Reset( int nCoefficients = -1 );

  void AddMeasurement( float f[], float m);
  
  int Fit( float Coefficients[] );    

private:

  int fN;
  long double *fA;
  long double *fB;
};

#endif
