/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef AliPriorsTask_H
#define AliPriorsTask_H

#include "AliAnalysisTaskSE.h"
#include "AliPID.h"
class TH1I;
class TH2D;

class AliPriorsTask : public AliAnalysisTaskSE {
  public:


  AliPriorsTask();
  AliPriorsTask(const Char_t* name);
  AliPriorsTask& operator= (const AliPriorsTask& c);
  AliPriorsTask(const AliPriorsTask& c);
  virtual ~AliPriorsTask();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  
  void SetPriors(Double_t conc[AliPID::kSPECIES])        {for(Int_t i=0; i<AliPID::kSPECIES; i++) fPriors[i]=conc[i];}
  void GetPriors(Double_t conc[AliPID::kSPECIES]) const  {for(Int_t i=0; i<AliPID::kSPECIES; i++) conc[i]=fPriors[i];}

  void SetNiterations(Int_t nIter) {fNiterMax = nIter;}

  Bool_t NotifyBinChange(); // method called at the end of the event loop

 protected:
 
  // Histograms

  TH1I  *fHistEventsProcessed; // simple histo for monitoring the number of events processed
  TH2D  *fCalcPriors;          // histo monitoring priors during iterations

  Double_t fPriors[AliPID::kSPECIES];  // Priors
  Double_t fRecId[AliPID::kSPECIES];   // Reconstructed Id
  Double_t fMCId[AliPID::kSPECIES];    // MC Id

  Int_t fNiterations;          //counter  
  Int_t fNiterMax;             //max number of iterations
  ClassDef(AliPriorsTask,1);
};

#endif
