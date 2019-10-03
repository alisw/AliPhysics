#ifndef AliAnalysisTaskVtXY_h
#define AliAnalysisTaskVtXY_h
/* Copyright (c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * see cxx source for full Copyright notice         */
//-------------------------------------------------------
//
// ESD based analysis of the main vertex resolution in XY in order
// to estimate the beam interaction spot location and size 
//
//------------------------------------------------------- 
class TStyle;
class TH2F;
class TProfile;
class AliESDEvent;
class AliVertex;
class AliESDVertex;
class AliVertexerTracks;
class AliESDVertexer;
#include "AliAnalysisTask.h"

class AliAnalysisTaskVtXY : public AliAnalysisTask {
 public:
  AliAnalysisTaskVtXY(const char *name = "AliAnalysisTaskVtXY");
  virtual ~AliAnalysisTaskVtXY() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDEvent *fESD;        //ESD object
  TList *fList;             //TList output object
  TProfile       *fHistVtx; //Vtx spectrum
  TProfile       *fHistVty; //Vty spectrum
  AliAnalysisTaskVtXY(const AliAnalysisTaskVtXY&); //not implemented
  AliAnalysisTaskVtXY& operator=(const AliAnalysisTaskVtXY&); //not implemented
  
  ClassDef(AliAnalysisTaskVtXY, 1); //example of analysis
};

#endif
