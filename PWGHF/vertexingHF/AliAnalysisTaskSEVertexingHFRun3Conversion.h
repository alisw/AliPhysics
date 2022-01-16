#ifndef ALIANALYSISTASKSEVERTEXINGHFRUN3CONVERSION_H
#define ALIANALYSISTASKSEVERTEXINGHFRUN3CONVERSION_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
/// \class Class AliAnalysisTaskSEVertexingHFRun3Conversion
/// \brief AliAnalysisTaskSE for the reconstruction of heavy-flavour decay candidates for conversion to AO2D
/// \author Author: F.Prino, prino@to.infn.it
//*************************************************************************

#include <TClonesArray.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class TList;
class TString;
class TTree;

class AliAnalysisTaskSEVertexingHFRun3Conversion : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEVertexingHFRun3Conversion();
  AliAnalysisTaskSEVertexingHFRun3Conversion(const char *name);
  virtual ~AliAnalysisTaskSEVertexingHFRun3Conversion();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void SetMakeReducedCandidates(Bool_t opt){fMakeReducedCandidates=opt;}
  void ResetTreeAtEachEvent(Bool_t opt){fResetTreeAtEachEv=opt;}
  AliAnalysisVertexingHF *GetVertexingHF() const {return fVHF;}
  
 private:

  AliAnalysisTaskSEVertexingHFRun3Conversion(const AliAnalysisTaskSEVertexingHFRun3Conversion &source);
  AliAnalysisTaskSEVertexingHFRun3Conversion& operator=(const AliAnalysisTaskSEVertexingHFRun3Conversion& source); 

  AliAnalysisVertexingHF *fVHF;        /// Vertexer heavy flavour
  Bool_t       fMakeReducedCandidates; /// flag to use reduced size candidates
  Bool_t       fResetTreeAtEachEv;     /// flag to reset trees at each event
  TTree        *fD0CandidateTree;      //!<! output tree
  TTree        *f3ProngCandidateTree;  //!<! output tree
  TTree        *fDstarCandidateTree;   //!<! output tree
  TTree        *fCascadeCandidateTree; //!<! output tree
  Int_t        fEventIndex;            //!<! tree variable
  Int_t        fD0track0;              //!<! tree variable
  Int_t        fD0track1;              //!<! tree variable
  Int_t        f3ptrack0;              //!<! tree variable
  Int_t        f3ptrack1;              //!<! tree variable
  Int_t        f3ptrack2;              //!<! tree variable
  Int_t        fDstD0tr0;              //!<! tree variable
  Int_t        fDstD0tr1;              //!<! tree variable
  Int_t        fDstSofPi;              //!<! tree variable
  Int_t        fCasV0ind;              //!<! tree variable
  Int_t        fCasV0tr0;              //!<! tree variable
  Int_t        fCasV0tr1;              //!<! tree variable
  Int_t        fCasBachl;              //!<! tree variable
  TList        *fListOfCuts;           //!<! List of analysis cuts
  TClonesArray *fVerticesHFTClArr;     /// Array of heavy-flavour vertices
  TClonesArray *fD0toKpiTClArr;        /// Array of D0->Kpi
  TClonesArray *fJPSItoEleTClArr;      /// Array of Jpsi->ee
  TClonesArray *fCharm3ProngTClArr;    /// Array of D+,Ds,Lc
  TClonesArray *fCharm4ProngTClArr;    /// Array of D0->Kpipipi
  TClonesArray *fDstarTClArr;          /// Array of D*->D0pi
  TClonesArray *fCascadesTClArr;       /// Array of Cascades : V0 + track (lambda_c)
  TClonesArray *fLikeSign2ProngTClArr; /// Array of LikeSign2Prong
  TClonesArray *fLikeSign3ProngTClArr; /// Array of LikeSign3Prong

  /// \cond CLASSIMP     
  ClassDef(AliAnalysisTaskSEVertexingHFRun3Conversion,3); /// AliAnalysisTaskSE for the reconstruction of heavy-flavour decay candidates
  /// \endcond
};

#endif

