#ifndef ALIANALYSISTASKHYPERONCORR_H
#define ALIANALYSISTASKHYPERONCORR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to prepeare lists of photons in calorimeters,
// converted photon etc. and fill few QA histograms
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskSigma0Spectra.h"

class AliAnalysisTaskSigma0Corr : public AliAnalysisTaskSigma0Spectra {
	
 public:
  AliAnalysisTaskSigma0Corr();
  AliAnalysisTaskSigma0Corr(const char* name);
  virtual ~AliAnalysisTaskSigma0Corr() ;// virtual destructor
		
  // Implementation of interface methods
  virtual void Init();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
		
 protected:
//   void ProcessMC();

  void FillCorrelation(AliCaloParticle * trig, const char * type="PHOS-Track") ;

  void FillDeltaPhi(AliCaloParticle * trig, const Int_t itype = 100) ;

 private:
  AliAnalysisTaskSigma0Corr(const AliAnalysisTaskSigma0Corr&); // Not implemented
  AliAnalysisTaskSigma0Corr& operator=(const AliAnalysisTaskSigma0Corr&); // Not implemented
		
 protected:
   		
  ClassDef(AliAnalysisTaskSigma0Corr, 1); // Analysis task for conversion + calorimeters
};

#endif //ALIANALYSISTASKHYPERON_H
