#ifndef ALIANALYSISTASKHYPERTRITON3NEW_H
#define ALIANALYSISTASKHYPERTRITON3NEW_H

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
#include "AliVertexerTracks.h"

#include <TMath.h>

class AliPIDResponse;
class AliESDtrackCuts;
class AliESDtrack;
class TTree;
class TH1D;
class TH3D;
class TList;
class TTree;

/// Define bit flags for datacompression
const unsigned char g = 0x1; // on if is a good candidate
const unsigned char r = 0x2; // on if is reflection candidate
/// if is not good and not reflection is a fake

class AliAnalysisTaskHypertriton3New : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskHypertriton3New(TString taskname = "Hypertriton3TaskNew");
  virtual ~AliAnalysisTaskHypertriton3New();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  AliEventCuts fEventCuts; /// event cuts class

private:
  AliAnalysisTaskHypertriton3New(const AliAnalysisTaskHypertriton3New &);            // not implemented
  AliAnalysisTaskHypertriton3New &operator=(const AliAnalysisTaskHypertriton3New &); // not implemented

  AliPIDResponse *fPIDResponse;   // PID response object
  AliESDtrackCuts *fESDtrackCuts; // ESD track cuts
  AliESDVertex *fPrimaryVertex;   // Primary vertex of the current event

  // Settings
  float fCosPoiningAngleLimit;

  // output object
  TList *fOutputList;   //! Output list
  TTree *fTree; //!

  // Findable Tree
  AliESDtrack *fTreeHyp3BodyVarTracks[3];
  Int_t fTreeHyp3BodyVarPDGcodes[3];

  ULong64_t fTreeHyp3BodyVarEventId;
  Int_t fTreeHyp3BodyVarMotherId;
  UChar_t fTreeHyp3BodyVarCandStat;

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

  TH1D *fHistEventCounter; //!
  TH1D *fHistCentrality;   //!

  TH3D *fHistGeneratedPtVsYVsCentralityHypTrit;     //!
  TH3D *fHistGeneratedPtVsYVsCentralityAntiHypTrit; //!

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskHypertriton3New, 1); // analysisclass
  /// \endcond
};

#endif
