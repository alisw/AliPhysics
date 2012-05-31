#ifndef ALIZDCCALIBTASK_cxx
#define ALIZDCCALIBTASK_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#define MAXENVALUE  10000 // Max. en. value for histos

#include "AliAnalysisTaskSE.h"

class TTree;
class TFile;
class TList;
class TH1F;
class AliESDEvent;

class AliZDCCalibTask : public AliAnalysisTaskSE {

 public:
  
  AliZDCCalibTask(const char *name = "ZDCCalibTask");
  virtual ~AliZDCCalibTask();
  AliZDCCalibTask(const AliZDCCalibTask&); 
  AliZDCCalibTask& operator=(const AliZDCCalibTask&); 
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *);
  virtual void  Terminate(Option_t *);
 
 private:

  void    BookHistos();
  void    DrawHistos();

  AliESDEvent* fESD ;	      //! Declaration of leave types
  TList*       fListOfHistos; //! Output list of histograms
  TH1F*        fZNCHisto;     //! ZNC histograms  
  TH1F*        fZPCHisto;     //! ZPC histograms  
  TH1F*        fZNAHisto;     //! ZNA histograms  
  TH1F*        fZPAHisto;     //! ZPA histograms  
  
  ClassDef(AliZDCCalibTask, 1); 

};

#endif
