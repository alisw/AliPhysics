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
#include "AliStack.h"

class TH1;
class TParticle;
class TFile;
class AliESDtrack;
class AliESDtrackCuts;
class AliAnalysisFilter;
class AliESDEvent;
class AliVEvent;
class AliPIDResposne;

class AliHMPIDAnalysisTask : public AliAnalysisTaskSE {
  public:

  enum {kChamber = 7};

  AliHMPIDAnalysisTask();
  AliHMPIDAnalysisTask(const Char_t* name);
  AliHMPIDAnalysisTask& operator= (const AliHMPIDAnalysisTask& c);
  AliHMPIDAnalysisTask(const AliHMPIDAnalysisTask& c);
  virtual ~AliHMPIDAnalysisTask();
  
  virtual void   ConnectInputData(Option_t *);
 // virtual void AliHMPIDAnalysisTask::UserCreateObject(Option_t *)
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

          void   SetUseMC(Bool_t useMC) { fUseMC = useMC; }
          Bool_t Equal(Double_t x, Double_t y, Double_t tolerance);

 protected:
     
 private:     
 
  AliESDEvent        *fESD;             //! ESD object
  AliMCEvent         *fMC;              //! MC event

  Bool_t             fUseMC;            // decide whether use or not the MC information

  TList              *fHmpHistList ;    // list of histograms

  TH1F               *fHmpNevents;
  TH1F               *fZvertex;

  AliPIDResponse     *fPIDResponse;
  AliESDtrackCuts    *fTrackCuts;
  AliAnalysisFilter  *fTrackFilter;
  
  TTree              *fTree;            // tree with useful data for subsequent analysis
  Float_t            fVar[69];          // array of data to fill the tree

  ClassDef(AliHMPIDAnalysisTask,4);
};

#endif
