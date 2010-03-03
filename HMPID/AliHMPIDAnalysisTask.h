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

//==============================================================================
// AliHMPIDAnalysysTask - Class representing a basic analysis tool of HMPID data at  
// level of ESD.
// A set of histograms is created.
//==============================================================================

#ifndef ALIHMPIDANALYSISTASK_H
#define ALIHMPIDANALYSISTASK_H

#include "AliAnalysisTaskSE.h"

class TH1;
class TH2;
class TParticle ;
class AliStack ;
class AliESDtrack;
class AliESDEvent;

class AliHMPIDAnalysisTask : public AliAnalysisTaskSE {
  public:

  enum {kChamber = 7};

  AliHMPIDAnalysisTask();
  AliHMPIDAnalysisTask(const Char_t* name);
  AliHMPIDAnalysisTask& operator= (const AliHMPIDAnalysisTask& c);
  AliHMPIDAnalysisTask(const AliHMPIDAnalysisTask& c);
  virtual ~AliHMPIDAnalysisTask();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
          Bool_t Equal(Double_t x, Double_t y, Double_t tolerance);
  
 protected:
     
 private:     
 
  void   SetTrigger(Int_t trigger) {fTrigger = trigger;}
  AliESDEvent   *fESD;             //! ESD object
  TList         *fHmpHistList ;    // list of histograms 
  Int_t          fNevts       ;    //event numbering
  Int_t          fTrigNevts   ;    //event numbering with the requested trigger
  Int_t          fTrigger     ;    //requested trigger
  TH2F          *fHmpPesdPhmp;     // HMP momentum vs ESD momentum
  TH2F          *fHmpCkovPesd;     // Ckov angle vs ESD momentum
  TH2F          *fHmpCkovPhmp;     // Ckov angle vs HMP momenutm
  TH1F          *fHmpMipTrkDist;   // Track-Mip distance distribution
  TH1F          *fHmpMipTrkDistX;  // Xtrk - Xmip
  TH1F          *fHmpMipTrkDistY;  // Ytrk - Ymip
  TH1F          *fHmpMipCharge3cm; // Mip charge with 3 cm distance cut 
  TH1F          *fHmpMipCharge1cm; // Mip charge with 1 cm distance cut
  TH1F          *fHmpNumPhots;     // Number of reconstructed photo-electrons
  TH1F          *fHmpTrkFlags;     // track flags
  
  ClassDef(AliHMPIDAnalysisTask,2);
};

#endif
