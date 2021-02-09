#ifndef ALIANALYSISTASKFINDABLEHYPERTRITON3_H
#define ALIANALYSISTASKFINDABLEHYPERTRITON3_H

/**************************************************************************
 *                                                                        *
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

///////////////////////////////////////////////////////////////////////////
//
// New class generating the tree of findable hypertritons in the 3 body
// channel. This class reimplements the hypertriton 3 case of the
// PWGLF/STRANGENESS/Cascades/Run2/AliAnalysisTaskStrEffStudy class
// with some improvements for this specific analysis.
//
// Author:
// P. Fecchio,  pfecchio@cern.ch
///////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVertexerTracks.h"

#include <TMath.h>

class AliESDtrackCuts;
class AliESDtrack;
class TTree;
class TH1D;
class TH3D;
class TList;
class TTree;

class AliAnalysisTaskFindableHypertriton3 : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskFindableHypertriton3(TString taskname = "FindableHypertriton3Task");
  virtual ~AliAnalysisTaskFindableHypertriton3();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  static  bool  HasTOF(AliESDtrack *t);

  AliEventCuts fEventCuts; /// event cuts class

private:
  AliAnalysisTaskFindableHypertriton3(const AliAnalysisTaskFindableHypertriton3 &);            // not implemented
  AliAnalysisTaskFindableHypertriton3 &operator=(const AliAnalysisTaskFindableHypertriton3 &); // not implemented

  AliPIDResponse *fPIDResponse;   //! PID response object
  AliESDtrackCuts *fESDtrackCuts; //! ESD track cuts
  AliESDVertex *fPrimaryVertex;   //! Primary vertex of the current event

  // Settings
  float fCosPoiningAngleLimit;

  // output object
  TList *fOutputList;   //! Output list
  TTree *fFindableTree; //!

  // Findable Tree
  AliESDtrack *fTreeHyp3BodyVarTracks[3]; //!
  Int_t fTreeHyp3BodyVarPDGcodes[3];
  Float_t fTreeHyp3BodyVarNsigmaTPC[3];
  Float_t fTreeHyp3BodyVarNsigmaTOF[3];

  ULong64_t fTreeHyp3BodyVarEventId;
  Int_t fTreeHyp3BodyVarMotherId;
  Bool_t fTreeHyp3BodyVarIsFakeCand;

  Float_t fTreeHyp3BodyVarTruePx;
  Float_t fTreeHyp3BodyVarTruePy;
  Float_t fTreeHyp3BodyVarTruePz;

  Float_t fTreeHyp3BodyVarDecayVx;
  Float_t fTreeHyp3BodyVarDecayVy;
  Float_t fTreeHyp3BodyVarDecayVz;
  Float_t fTreeHyp3BodyVarDecayT;

  Float_t fTreeHyp3BodyVarPVx;
  Float_t fTreeHyp3BodyVarPVy;
  Float_t fTreeHyp3BodyVarPVz;
  Float_t fTreeHyp3BodyVarPVt;

  Float_t fTreeHyp3BodyVarMagneticField;
  Float_t fTreeHyp3BodyVarCentrality;

  TH1D *fHistEventCounter; //!
  TH1D *fHistCentrality;   //!

  TH3D *fHistGeneratedPtVsCtVsCentralityHypTrit3;     //!
  TH3D *fHistGeneratedPtVsCtVsCentralityAntiHypTrit3; //!

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskFindableHypertriton3, 1); // analysisclass
  /// \endcond
};

#endif