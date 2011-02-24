#ifndef AliAnalysisTaskMEVertexingHF_cxx
#define AliAnalysisTaskMEVertexingHF_cxx

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *  * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskMEVertexingHF
// AliAnalysisTaskME for event mixing, building the background for 
// heavy-flavour decay candidates
// Author: R.Romita, r.romita@gsi.de
//*************************************************************************


class TH1F;
class TList;
class AliESDEvent;

#include "AliAnalysisTaskME.h"
#include "AliAnalysisVertexingHF.h"
#include "AliMixedEvent.h"
#include <TClonesArray.h>


class AliAnalysisTaskMEVertexingHF : public AliAnalysisTaskME {
 public:
  AliAnalysisTaskMEVertexingHF(const char *name = "AliAnalysisTaskMEVertexingHF");
  virtual ~AliAnalysisTaskMEVertexingHF() {}
  
  virtual void   Init();
  virtual void   LocalInit() {Init();}
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliAnalysisVertexingHF *fvHF;      // Vertexer heavy flavour
  AliMixedEvent *fMixedEvent;        // Mixed event
  TClonesArray  *fVerticesHFTClArr;  //Array of heavy-flavour vertices
  TClonesArray  *fD0toKpiTClArr;     //Array of D0->Kpi
  TClonesArray  *fJPSItoEleTClArr;   //Array of Jpsi->ee
  TClonesArray  *fCharm3ProngTClArr; //Array of D+,Ds,Lc
  TClonesArray  *fCharm4ProngTClArr; // Array of D0->Kpipipi
  TClonesArray  *fDstarTClArr;          // Array of D*->D0pi
  TClonesArray  *fCascadesTClArr;       // Array of Cascades : V0 + track (lambda_c)
  TClonesArray  *fLikeSign2ProngTClArr; // Array of LikeSign2Prong
  TClonesArray  *fLikeSign3ProngTClArr; // Array of LikeSign3Prong


  AliAnalysisTaskMEVertexingHF(const AliAnalysisTaskMEVertexingHF&); // not implemented
  AliAnalysisTaskMEVertexingHF& operator=(const AliAnalysisTaskMEVertexingHF&); // not implemented
  
  ClassDef(AliAnalysisTaskMEVertexingHF,1); // event mixing for heavy-flavour vertexing
};

#endif
