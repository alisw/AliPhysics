#ifndef ALIJCORRANTASK_H
#define ALIJCORRANTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

//#include <TTree.h>
//#include <TList.h>
#include <TVectorT.h>
#include "AliAnalysisTaskSE.h"
//#include "AliAnalysisFilter.h"
//#include "AliMCEvent.h"
#include "AliJRunHeader.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
//#include "AliESDFMD.h"
#include "AliESDZDC.h"
#include "AliJConst.h"
#include "AliESDpid.h"
#include "AliPHOSGeoUtils.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliJFilter.h"
#include "AliJEfficiencyScanner.h"

//==============================================================

using namespace std;

class TH1D;
class TH2D;
class TNtuple;
class TList;
class TTree;
class TFormula;

class AliMCEvent; 
class AliESDEvent; 
class AliAODEvent; 
class AliAODTrack; 
class AliESDtrackCuts;
class AliESDVZERO;
class AliESDCentrality;

class AliJRunHeader;
class AliMCEvent;
class AliAnalysisFilter;
class AliJTrack;
class AliJEventHeader;



class AliJCORRANTask : public AliAnalysisTaskSE {

 public:
  AliJCORRANTask();
  AliJCORRANTask(const char *name,  TString inputformat);
  AliJCORRANTask(const AliJCORRANTask& ap);   
  AliJCORRANTask& operator = (const AliJCORRANTask& ap);
  virtual ~AliJCORRANTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t* );

  AliJFilter *GetFilter() { return fFilter; }

  AliJRunHeader * GetJRunHeader(){ return fAliJRunHeader; }
  void SetJRunHeader( AliJRunHeader * hdr ){ fAliJRunHeader = hdr ; }

 private:
  
  AliJFilter *fFilter; // filter object
  AliJRunHeader * fAliJRunHeader;

  ClassDef(AliJCORRANTask, 3); 
};
#endif // AliJCORRANTask_H
