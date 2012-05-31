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
//
// HMPID PID efficiency using the VO.
// 
// Francesco Barile
//
//==============================================================================

#ifndef ALIHMPIDPIDEfficiencyV0ANALYSISTASK_H
#define ALIHMPIDPIDEfficiencyV0ANALYSISTASK_H

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"

class TH1;
class TParticle;
class TFile;
class AliESDtrack;
class AliESDEvent;

class AliHMPIDPIDEfficiencyV0AnalysisTask : public AliAnalysisTaskSE {
  public:

  enum {kChamber = 7};

  AliHMPIDPIDEfficiencyV0AnalysisTask();
  AliHMPIDPIDEfficiencyV0AnalysisTask(const Char_t* name);
  AliHMPIDPIDEfficiencyV0AnalysisTask& operator= (const AliHMPIDPIDEfficiencyV0AnalysisTask& c);
  AliHMPIDPIDEfficiencyV0AnalysisTask(const AliHMPIDPIDEfficiencyV0AnalysisTask& c);
  virtual ~AliHMPIDPIDEfficiencyV0AnalysisTask();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

          void   SetUseMC(Bool_t useMC) { fUseMC = useMC; }
          Bool_t Equal(Double_t x, Double_t y, Double_t tolerance);

 protected:
     
 private:     
 
  AliESDEvent *fESD;               //! ESD object
  AliMCEvent  *fMC;                //! MC event

  Bool_t       fUseMC;             // decide whether use or not the MC information

  TList         *fHmpHistList ;    // list of histograms

  TH1F          *fHmpNevents;
  TH2F          *fThetavsMom;
  TH1D          *fmass;
  TH1D          *fpdg;
  TH2F          *farme;
  TH2F          *farmeb;
  TH1D          *hangle;
  TH1D          *massap;
  TH1D          *massan;
  TH1D          *fmassaHMPID;  
  TH1D          *hdiff;
  TH1D          *hdiffn;
  TTree         *fTree;            // tree with useful data for subsequent analysis
  Float_t        fVar[90];         // array of data to fill the tree

  ClassDef(AliHMPIDPIDEfficiencyV0AnalysisTask,4);
};

#endif
