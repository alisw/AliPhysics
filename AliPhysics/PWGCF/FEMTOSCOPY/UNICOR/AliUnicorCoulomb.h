#ifndef ALIUNICORCOULOMB_H
#define ALIUNICORCOULOMB_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2010

//=============================================================================
// Coulomb correlation function
//=============================================================================

#include <TComplex.h>
#include <TGraph.h>

//=============================================================================
class AliUnicorCoulomb : public TGraph {

 public: 
  AliUnicorCoulomb(TRootIOCtor *) : TGraph() {}             // default constructor
  AliUnicorCoulomb(int sign, double mass, double R);        // constructor
  virtual ~AliUnicorCoulomb() {}                            // destructor
  double Cf(double qinv) const {return Eval(qinv);} // value of the correlation function
  static double Gamow(int zz, double m, double k);  // poin source case - Gamow function
  static void Makehist(int zz, double m, const char *outfil); // make TH2(R,Q)

 protected:
  static double   WaveFunction2(int zz, double mass, double k, double x, double y, double z);
  static TComplex WaveFunction(int zz, double mass, double k, double x, double y, double z);
  static TComplex F1(TComplex alpha, TComplex gamma, TComplex z);

  ClassDef(AliUnicorCoulomb,1)
};
//=============================================================================
#endif
