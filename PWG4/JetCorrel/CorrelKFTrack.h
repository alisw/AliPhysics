#ifndef CORRELKFTRACK_H
#define CORRELKFTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//__________________________________________________
// Container class for a track with AliKF parameters
//-- Author: Paul Constantin

#include "CorrelParticle.h"

class CorrelKFTrack_t : public CorrelParticle_t {
 public:
  
  CorrelKFTrack_t();
  CorrelKFTrack_t(Float_t pt, Float_t p, Float_t e, Float_t m, cPartType_t i, 
		  Double_t* par, Double_t* cov);
  CorrelKFTrack_t(const CorrelKFTrack_t &p);
  virtual ~CorrelKFTrack_t() {;}
  CorrelKFTrack_t& operator=(const CorrelKFTrack_t& rhs);
  virtual CorrelKFTrack_t* Copy();
  
  void SetParam(const Double_t* v) {fParam=(Double_t*)v;}
  void SetCovar(const Double_t* v) {fCovar=(Double_t*)v;}
  Double_t* Param() const {return fParam;}
  Double_t* Covar() const {return fCovar;}
  
  virtual void Show() const;
  
 private:
  Double_t* fParam; // Param[6] = {X, Y, Z, Px, Py, Pz} - position and momentum
  Double_t* fCovar; // Covar[21] = lower-triangular part of the covariance matrix
};

#endif
