#ifndef ALIANALYSISTASKSEVERTEXINGHF_H
#define ALIANALYSISTASKSEVERTEXINGHF_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSEVertexingHF
// AliAnalysisTaskSE for the reconstruction of heavy-flavour decay candidates
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//*************************************************************************

#include <TClonesArray.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODHFUtil.h"

class TList;
class TString;

class AliAnalysisTaskSEVertexingHF : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEVertexingHF();
  AliAnalysisTaskSEVertexingHF(const char *name);
  virtual ~AliAnalysisTaskSEVertexingHF();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void SetDeltaAODFileName(const char* name) {fDeltaAODFileName=name;}
  const char* GetDeltaAODFileName() const {return fDeltaAODFileName.Data();}
  AliAnalysisVertexingHF *GetVertexingHF() const {return fVHF;}
  
 private:

  AliAnalysisTaskSEVertexingHF(const AliAnalysisTaskSEVertexingHF &source);
  AliAnalysisTaskSEVertexingHF& operator=(const AliAnalysisTaskSEVertexingHF& source); 

  AliAnalysisVertexingHF *fVHF;        // Vertexer heavy flavour
  TList        *fListOfCuts;           //! List of analysis cuts
  TString       fDeltaAODFileName;     // Name of output file
  TClonesArray *fVerticesHFTClArr;     // Array of heavy-flavour vertices
  TClonesArray *fD0toKpiTClArr;        // Array of D0->Kpi
  TClonesArray *fJPSItoEleTClArr;      // Array of Jpsi->ee
  TClonesArray *fCharm3ProngTClArr;    // Array of D+,Ds,Lc
  TClonesArray *fCharm4ProngTClArr;    // Array of D0->Kpipipi
  TClonesArray *fDstarTClArr;          // Array of D*->D0pi
  TClonesArray *fCascadesTClArr;       // Array of Cascades : V0 + track (lambda_c)
  TClonesArray *fLikeSign2ProngTClArr; // Array of LikeSign2Prong
  TClonesArray *fLikeSign3ProngTClArr; // Array of LikeSign3Prong 
  AliAODHFUtil *fHFUtilInfo;              // VZERO branch (to be removed)

  ClassDef(AliAnalysisTaskSEVertexingHF,6); // AliAnalysisTaskSE for the reconstruction of heavy-flavour decay candidates
};

#endif

