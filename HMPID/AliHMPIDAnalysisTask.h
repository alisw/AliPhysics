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


#ifndef AliHMPIDAnalysisTASK_H
#define AliHMPIDAnalysisTASK_H

#include "AliAnalysisTaskSE.h"

class TH1I;
class TParticle ;
class TFile ;
class AliStack ;
class AliESDtrack;

class AliHMPIDAnalysisTask : public AliAnalysisTaskSE {
  public:

  enum {kChamber = 7};

  AliHMPIDAnalysisTask();
  AliHMPIDAnalysisTask(const Char_t* name);
  AliHMPIDAnalysisTask& operator= (const AliHMPIDAnalysisTask& c);
  AliHMPIDAnalysisTask(const AliHMPIDAnalysisTask& c);
  virtual ~AliHMPIDAnalysisTask();

  // ANALYSIS FRAMEWORK STUFF
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);

  // setter
  void   SetTrigger(Int_t trigger) {fTrigger = trigger;}

 protected:
  TList          *fHistList ;  // list of histograms
  
  //Number of events
  TH1I  *fHistEventsProcessed ;  // simple histo for monitoring number of processed events
  Int_t          fNevts       ;  //event numbering
  Int_t          fTrigNevts   ;  //event numbering with the requested trigger
  Int_t          fTrigger     ;  //requested trigger
  
  ClassDef(AliHMPIDAnalysisTask,1);
};

#endif
