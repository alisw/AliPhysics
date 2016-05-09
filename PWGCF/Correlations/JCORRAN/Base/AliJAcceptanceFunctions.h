/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#ifndef ALIJACCEPTANCEFUNCTIONS_H
#define ALIJACCEPTANCEFUNCTIONS_H

/*
 * Necessary function definitions for calculating 3D near side acceptance correction
 */
class AliJAcceptanceFunctions{
public:
  AliJAcceptanceFunctions();       // constructor
  ~AliJAcceptanceFunctions(){;}    // destructor
  double AcceptanceCorrection3DNearSide(double *x, double *par);
  double SideBoundaryDistance(double *x, double *par);
};

#endif
