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
class AliESDEvent;

class AliHMPIDAnalysisTask : public AliAnalysisTaskSE {
  public:

  enum {kChamber = 7};

  AliHMPIDAnalysisTask();
  AliHMPIDAnalysisTask(const Char_t* name);
  virtual ~AliHMPIDAnalysisTask();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  
 protected:
     
 private:     
 
  void   SetTrigger(Int_t trigger) {fTrigger = trigger;}
  AliESDEvent *fESD;              //! ESD object
  TList          *fHmpHistList ;  // list of histograms 
  Int_t          fNevts       ;  //event numbering
  Int_t          fTrigNevts   ;  //event numbering with the requested trigger
  Int_t          fTrigger     ;  //requested trigger
  TH2F          *fHmpInner;
  TH2F          *fHmpPesdPhmp;
  TH2F          *fHmpCkovPesd;
  TH2F          *fHmpCkovPhmp;
  
  TH1F          *fHmpMipTrkDist;
  TH1F          *fHmpMipTrkDistX;
  TH1F          *fHmpMipTrkDistY;
  TH1F          *fHmpMipCharge3cm;
  TH1F          *fHmpMipCharge1cm;
  TH1F          *fHmpNumPhots;
  
  TH1F          *fHmpTrkFlags;
  
  ClassDef(AliHMPIDAnalysisTask,2);
};

#endif
