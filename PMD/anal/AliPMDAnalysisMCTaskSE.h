#ifndef AliPMDAnalysisMCTaskSE_cxx
#define AliPMDAnalysisMCTaskSE_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**************************************************************************
 * Analysis Class Implimentation for the MC Truth and responce Matrix 
 *       Auther:   Satyajit Jena, IIT Bombay |  sjena@cern.ch
 * 
 *                     Mon Nov 22 19:54:27 CET 2010
 *
 **************************************************************************/

class TList;
class TFile;
class TH1F;
class TH1D;
class TH2F;
class AliESDEvent;
class AliMCEvent;
class AliESDPmdTrack;

#include "AliAnalysisTaskSE.h"

class AliPMDAnalysisMCTaskSE : public AliAnalysisTaskSE {
 public:
    
  AliPMDAnalysisMCTaskSE(const char *name = "AliPMDAnalysisMCTaskSE");
  virtual ~AliPMDAnalysisMCTaskSE() {}
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
 
  void EventByEvent(AliESDEvent* esd, AliMCEvent* mcEvent);
  
 private:

  TList* fPhysList;
      
  // TH1F* fhCounter; //!Event Counter Book Keeping
  // TH1F* fhVtx;     //!Vertex Cut X
  // TH1F* fhVty;
  // TH1F* fhVtz;

  // Int_t fCntr;     //!Event Counter

  TH2F *fhResponseAll;  //!2D Responce Matrix
  TH1F *fhTrueAll;      //!1D True Multiplicity
  TH1F *fhMeasuredAll;

  TH2F *fhResponse[10];
  TH1F *fhTrue[10];
  TH1F *fhMeasured[10];
  


   //___________________________________________________
  AliPMDAnalysisMCTaskSE(const AliPMDAnalysisMCTaskSE&); 
  AliPMDAnalysisMCTaskSE& operator=(const AliPMDAnalysisMCTaskSE&); 
  
  ClassDef(AliPMDAnalysisMCTaskSE, 1); 
};

#endif

