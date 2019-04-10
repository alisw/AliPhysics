#ifndef ALIANALYSISTASKCOUNTEVENTS
#define ALIANALYSISTASKCOUNTEVENTS

/* Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskCountEvents
// AliAnalysisTaskSE to count events per trigger mask after various selections
// 
//
// Author:
//          F. Prino, prino@to.infn.it
//          
//*************************************************************************

class TList;
class TH1F;
class TH2F;
class TH3F;
class TString;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCountEvents : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskCountEvents();
  virtual ~AliAnalysisTaskCountEvents();

  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);

 private:
  void ConfigureXaxis(TH1* histo);
  
  AliAnalysisTaskCountEvents(const AliAnalysisTaskCountEvents &source);
  AliAnalysisTaskCountEvents& operator=(const AliAnalysisTaskCountEvents &source);
  
  TList*  fOutput;                    //!<!  list of output histos
  TH1F* fHistNEventsPhysSel;          //!<!  histo with N of events  
  TH1F* fHistNEventsSPDVert;          //!<!  histo with N of events  
  TH1F* fHistNEventsTrackVert;        //!<!  histo with N of events  
  TH1F* fHistNEventsZvert10cm;        //!<!  histo with N of events  
  TH2F* fHistNEventsPhysSelVsCent;    //!<!  histo with N of events vs. centr.  
  TH2F* fHistNEventsSPDVertVsCent;    //!<!  histo with N of events vs. centr.  
  TH2F* fHistNEventsTrackVertVsCent;  //!<!  histo with N of events vs. centr.  
  TH2F* fHistNEventsZvert10cmVsCent;  //!<!  histo with N of events vs. centr.  

  ClassDef(AliAnalysisTaskCountEvents,2);
};


#endif
