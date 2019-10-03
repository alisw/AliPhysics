/**
 * \file AliTLorentzVector.h
 * \brief Declaration of class AliTLorentzVector
 *
 * In this header file the class AliTLorentzVector is declared.
 * It is a simple reimplementation of TLorentzVector with an added method Phi_0_TwoPi() which returns Phi in 0-2Pi range.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Jan 29, 2016
 */
#ifndef ALITLORENTZVECTOR_H
#define ALITLORENTZVECTOR_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TLorentzVector.h>

class AliTLorentzVector : public TLorentzVector {
public:
  AliTLorentzVector();
  AliTLorentzVector(Double_t x, Double_t y, Double_t z, Double_t t);
  AliTLorentzVector(const Double_t *carray);
  AliTLorentzVector(const Float_t *carray);
  AliTLorentzVector(const TVector3 &vector3, Double_t t);
  AliTLorentzVector(const TLorentzVector &lorentzvector);
  virtual ~AliTLorentzVector();

  Double_t Phi_0_2pi() const;

  ClassDef(AliTLorentzVector, 1) // AliTLorentzVector class
};

#endif /* ALITLORENTZVECTOR_H */
