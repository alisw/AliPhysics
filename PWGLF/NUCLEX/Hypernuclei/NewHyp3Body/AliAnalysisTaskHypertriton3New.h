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
// New class for the analysis of the hypertriton 3 prongs decay.
//
// Author:
// P. Fecchio,  pfecchio@cern.ch
///////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPID.h"
#include "AliVertexerTracks.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <TMath.h>

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> FourVector_t;

class AliPIDResponse;
class AliESDEvent;
class AliESDtrack;
class TTree;

// class HyperTritonCandidate {

// private:
//   short mass;
// };

class AliAnalysisTaskHypertriton3New : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskHypertriton3New();
  AliAnalysisTaskHypertriton3New(bool isMC, TString taskname);
  virtual ~AliAnalysisTaskHypertriton3New();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  AliEventCuts fEventCuts;     /// event cuts class
  AliESDtrackCuts *fTrackCuts; /// ESD track cuts

private:
  AliAnalysisTaskHypertriton3New(const AliAnalysisTaskHypertriton3New &);            // not implemented
  AliAnalysisTaskHypertriton3New &operator=(const AliAnalysisTaskHypertriton3New &); // not implemented

  bool PassPIDSelection(AliESDtrack *track, AliPID::EParticleType specie, float nSigmaCut);

  // settings
  bool fIsMC;                  ///<  Switch between MC and data
  float fNSigma[3];            ///<
  float fCosPoiningAngleLimit; ///<

  // support objects
  AliESDEvent *fESDevent;             ///< ESD event
  const AliESDVertex *fPrimaryVertex; //! Primary vertex of the current event
  AliPIDResponse *fPIDResponse;       //! PID response object

  std::vector<AliESDtrack *> fParticles[3]; ///<

  // output list object
  TList *fOutputList; //! Output list
  TTree *fHypTree;    //! Output tree

  // objects in the tree
  float fMagneticField;

  TH1D *fMInvBackgroundStd;  //!<!
  TH1D *fMInvBackgroundCuts; //!<!
  TH1D *fHistPtStd[3];       //!<!
  TH1D *fHistPtCuts[3];      //!<!
  TH1D *fHistDCAXY[3];       //!<!
  TH1D *fHistDCAZ[3];        //!<!
  TH1D *fHistCosPAngle;      //!<!

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskHypertriton3New, 1); // analysisclass
  /// \endcond
};

#endif
