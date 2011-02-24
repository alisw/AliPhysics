#ifndef ALIANALYSISTASKSEMONITNORM_cxx
#define ALIANALYSISTASKSEMONITNORM_cxx

/* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// 
// Class for monitoring of information used for open charm normalization
// (triggers, candles, ...)
//
// Origin: davide.caffarri@pd.infn.it
//
//*************************************************************************


class TH1F;
class TH2F;
class AliESDEvent;
class AliCounterCollection;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSEMonitNorm : public AliAnalysisTaskSE 
{
 public:

  AliAnalysisTaskSEMonitNorm(const char *name = "AliAnalysisTaskSEMonitNorm");
  virtual ~AliAnalysisTaskSEMonitNorm(); 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  

 protected:
  
  AliESDEvent *fESD;            // ESD object
  TList       *fOutput;         //! list send on output slot 0
  AliCounterCollection *fCounterTrigg; //! counter for the differents triggered events. 
  AliCounterCollection *fCounterNotTrigg; //! counter for the differents not triggered events.
  AliCounterCollection *fCounterCandleTrig; //! counter for candles in the triggered events (esd)
  AliCounterCollection *fCounterCandleNotTrig; //! counter for candles in the triggered events (esd)  

 private:    

  AliAnalysisTaskSEMonitNorm(const AliAnalysisTaskSEMonitNorm&); // not implemented
  AliAnalysisTaskSEMonitNorm& operator=(const AliAnalysisTaskSEMonitNorm&); // not implemented

  
  ClassDef(AliAnalysisTaskSEMonitNorm,1); // class for monitoring of normalization information
};

#endif
