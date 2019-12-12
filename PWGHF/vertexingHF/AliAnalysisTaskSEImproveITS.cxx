/*************************************************************************
* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
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

#include <TObjArray.h>
#include <TClonesArray.h>
#include <TGraph.h>
#include <TFile.h>
#include <TList.h>
#include <TNtuple.h>

#include "AliVertex.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliVertexerTracks.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisVertexingHF.h"
#include "AliNeutralTrackParam.h"
#include "AliAnalysisTaskSEImproveITS.h"
#include "AliVertexingHFUtils.h"
//
// Implementation of the "hybrid-approach" for ITS upgrade studies.
// The tastk smears the track parameters according to estimations
// from single-track upgrade studies. Afterwards it recalculates
// the parameters of the reconstructed decays.
//
// WARNING: This will affect all tasks in a train after this one
// (which is typically desired, though).
//

AliAnalysisTaskSEImproveITS::AliAnalysisTaskSEImproveITS() 
  :AliAnalysisTaskSE(),
   fD0ZResPCur  (0),
   fD0ZResKCur  (0),
   fD0ZResPiCur (0),
   fD0ZResECur (0),
   fD0RPResPCur (0),
   fD0RPResKCur (0),
   fD0RPResPiCur(0),
   fD0RPResECur(0),
   fD0RPSigmaPullRatioP(0),
   fD0RPSigmaPullRatioK(0),
   fD0RPSigmaPullRatioPi(0),
   fD0RPSigmaPullRatioE(0),
   fPt1ResPCur  (0),
   fPt1ResKCur  (0),
   fPt1ResPiCur (0),
   fPt1ResECur (0),
   fD0ZResPUpg  (0),
   fD0ZResKUpg  (0),
   fD0ZResPiUpg (0),
   fD0ZResEUpg (0),
   fD0RPResPUpg (0),
   fD0RPResKUpg (0),
   fD0RPResPiUpg(0),
   fD0RPResEUpg(0),
   fPt1ResPUpg  (0),
   fPt1ResKUpg  (0),
   fPt1ResPiUpg (0),
   fPt1ResEUpg (0),
   fD0ZResPCurSA  (0),
   fD0ZResKCurSA  (0),
   fD0ZResPiCurSA (0),
   fD0ZResECurSA (0),
   fD0RPResPCurSA (0),
   fD0RPResKCurSA (0),
   fD0RPResPiCurSA(0),
   fD0RPResECurSA(0),
   fPt1ResPCurSA  (0),
   fPt1ResKCurSA  (0),
   fPt1ResPiCurSA (0),
   fPt1ResECurSA (0),
   fD0ZResPUpgSA  (0),
   fD0ZResKUpgSA  (0),
   fD0ZResPiUpgSA (0),
   fD0ZResEUpgSA (0),
   fD0RPResPUpgSA (0),
   fD0RPResKUpgSA (0),
   fD0RPResPiUpgSA(0),
   fD0RPResEUpgSA(0),
   fPt1ResPUpgSA  (0),
   fPt1ResKUpgSA  (0),
   fPt1ResPiUpgSA (0),
   fPt1ResEUpgSA (0),
   // specific stuff for PbPb 2018
   fIsPbPb2018(kFALSE)
   // kFirst
  ,fD0ZResPCur_PbPb2018_kFirst(0)
  ,fD0ZResKCur_PbPb2018_kFirst(0) 
  ,fD0ZResPiCur_PbPb2018_kFirst(0) 
  ,fD0ZResECur_PbPb2018_kFirst(0) 
  ,fD0RPResPCur_PbPb2018_kFirst(0) 
  ,fD0RPResKCur_PbPb2018_kFirst(0) 
  ,fD0RPResPiCur_PbPb2018_kFirst(0) 
  ,fD0RPResECur_PbPb2018_kFirst(0)
  ,fD0RPSigmaPullRatioP_PbPb2018_kFirst(0) 
  ,fD0RPSigmaPullRatioK_PbPb2018_kFirst(0) 
  ,fD0RPSigmaPullRatioPi_PbPb2018_kFirst(0) 
  ,fD0RPSigmaPullRatioE_PbPb2018_kFirst(0)
  ,fPt1ResPCur_PbPb2018_kFirst(0)
  ,fPt1ResKCur_PbPb2018_kFirst(0)
  ,fPt1ResPiCur_PbPb2018_kFirst(0) 
  ,fPt1ResECur_PbPb2018_kFirst(0) 
  ,fD0ZResPUpg_PbPb2018_kFirst(0) 
  ,fD0ZResKUpg_PbPb2018_kFirst(0) 
  ,fD0ZResPiUpg_PbPb2018_kFirst(0) 
  ,fD0ZResEUpg_PbPb2018_kFirst(0)
  ,fD0RPResPUpg_PbPb2018_kFirst(0)
  ,fD0RPResKUpg_PbPb2018_kFirst(0) 
  ,fD0RPResPiUpg_PbPb2018_kFirst(0) 
  ,fD0RPResEUpg_PbPb2018_kFirst(0) 
  ,fPt1ResPUpg_PbPb2018_kFirst(0)
  ,fPt1ResKUpg_PbPb2018_kFirst(0) 
  ,fPt1ResPiUpg_PbPb2018_kFirst(0)
  ,fPt1ResEUpg_PbPb2018_kFirst(0) 
  ,fD0ZResPCurSA_PbPb2018_kFirst(0) 
  ,fD0ZResKCurSA_PbPb2018_kFirst(0)
  ,fD0ZResPiCurSA_PbPb2018_kFirst(0) 
  ,fD0ZResECurSA_PbPb2018_kFirst(0) 
  ,fD0RPResPCurSA_PbPb2018_kFirst(0) 
  ,fD0RPResKCurSA_PbPb2018_kFirst(0) 
  ,fD0RPResPiCurSA_PbPb2018_kFirst(0) 
  ,fD0RPResECurSA_PbPb2018_kFirst(0) 
  ,fPt1ResPCurSA_PbPb2018_kFirst(0) 
  ,fPt1ResKCurSA_PbPb2018_kFirst(0) 
  ,fPt1ResPiCurSA_PbPb2018_kFirst(0) 
  ,fPt1ResECurSA_PbPb2018_kFirst(0)
  ,fD0ZResPUpgSA_PbPb2018_kFirst(0) 
  ,fD0ZResKUpgSA_PbPb2018_kFirst(0) 
  ,fD0ZResPiUpgSA_PbPb2018_kFirst(0) 
  ,fD0ZResEUpgSA_PbPb2018_kFirst(0) 
  ,fD0RPResPUpgSA_PbPb2018_kFirst(0) 
  ,fD0RPResKUpgSA_PbPb2018_kFirst(0) 
  ,fD0RPResPiUpgSA_PbPb2018_kFirst(0) 
  ,fD0RPResEUpgSA_PbPb2018_kFirst(0) 
  ,fPt1ResPUpgSA_PbPb2018_kFirst(0)
  ,fPt1ResKUpgSA_PbPb2018_kFirst(0) 
  ,fPt1ResPiUpgSA_PbPb2018_kFirst(0)
  ,fPt1ResEUpgSA_PbPb2018_kFirst(0) 
  // kOnlySecond
  ,fD0ZResPCur_PbPb2018_kOnlySecond(0) 
  ,fD0ZResKCur_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPiCur_PbPb2018_kOnlySecond(0) 
  ,fD0ZResECur_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPCur_PbPb2018_kOnlySecond(0) 
  ,fD0RPResKCur_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPiCur_PbPb2018_kOnlySecond(0) 
  ,fD0RPResECur_PbPb2018_kOnlySecond(0) 
  ,fD0RPSigmaPullRatioP_PbPb2018_kOnlySecond(0) 
  ,fD0RPSigmaPullRatioK_PbPb2018_kOnlySecond(0) 
  ,fD0RPSigmaPullRatioPi_PbPb2018_kOnlySecond(0) 
  ,fD0RPSigmaPullRatioE_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPCur_PbPb2018_kOnlySecond(0) 
  ,fPt1ResKCur_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPiCur_PbPb2018_kOnlySecond(0)
  ,fPt1ResECur_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPUpg_PbPb2018_kOnlySecond(0) 
  ,fD0ZResKUpg_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPiUpg_PbPb2018_kOnlySecond(0) 
  ,fD0ZResEUpg_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPUpg_PbPb2018_kOnlySecond(0) 
  ,fD0RPResKUpg_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPiUpg_PbPb2018_kOnlySecond(0) 
  ,fD0RPResEUpg_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPUpg_PbPb2018_kOnlySecond(0) 
  ,fPt1ResKUpg_PbPb2018_kOnlySecond(0)  
  ,fPt1ResPiUpg_PbPb2018_kOnlySecond(0) 
  ,fPt1ResEUpg_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResKCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPiCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResECurSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResKCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPiCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResECurSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPCurSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResKCurSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPiCurSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResECurSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResKUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPiUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResEUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResKUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPiUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResEUpgSA_PbPb2018_kOnlySecond(0)  
  ,fPt1ResPUpgSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResKUpgSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPiUpgSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResEUpgSA_PbPb2018_kOnlySecond(0),
   fRunInVertexing(kFALSE),
   fImproveTracks(kTRUE),
   fUpdateSecVertCovMat(kTRUE),
   fUpdateSTCovMatrix(kTRUE),
   fUpdatePulls(kTRUE),
   fMimicData(kFALSE),
   fIsAOD        (kTRUE),
   fSmearOnlySignal(kFALSE),
   fMCs         (0),
   fDebugOutput (0),
   fDebugNtuple (0),
   fDebugVars   (0), 
   fNDebug      (0)
{
  //
  // Default constructor.
  for(Int_t jh=0; jh<2; jh++){
    // templates of mean (original)
    for(Int_t ih=0; ih<4; ih++){
      fD0RPMeanPCur[jh][ih]=0x0;
      fD0RPMeanKCur[jh][ih]=0x0;
      fD0RPMeanPiCur[jh][ih]=0x0;
      fD0RPMeanECur[jh][ih]=0x0;
      fD0RPMeanPUpg[jh][ih]=0x0;
      fD0RPMeanKUpg[jh][ih]=0x0;
      fD0RPMeanPiUpg[jh][ih]=0x0;
      fD0RPMeanEUpg[jh][ih]=0x0;
    }
    // templates of mean for Pb-Pb 2018 
    for(UInt_t ih = 0; ih < 24; ih++)
    {
      fD0RPMeanPCur_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanKCur_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanPiCur_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanECur_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanPUpg_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanKUpg_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanPiUpg_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanEUpg_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanPCur_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanKCur_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanPiCur_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanECur_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanPUpg_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanKUpg_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanPiUpg_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanEUpg_PbPb2018_kOnlySecond[jh][ih]=0x0;
    }
    
  }
  //
}

AliAnalysisTaskSEImproveITS::AliAnalysisTaskSEImproveITS(const char *name,
                           const char *period,
                           const char *systematic,
                           Bool_t isRunInVertexing,
			   Int_t ndebug)
  :AliAnalysisTaskSE(name),
   fD0ZResPCur  (0),
   fD0ZResKCur  (0),
   fD0ZResPiCur (0),
   fD0ZResECur (0),
   fD0RPResPCur (0),
   fD0RPResKCur (0),
   fD0RPResPiCur(0),
   fD0RPResECur(0),
   fD0RPSigmaPullRatioP(0),
   fD0RPSigmaPullRatioK(0),
   fD0RPSigmaPullRatioPi(0),
   fD0RPSigmaPullRatioE(0),
   fPt1ResPCur  (0),
   fPt1ResKCur  (0),
   fPt1ResPiCur (0),
   fPt1ResECur (0),
   fD0ZResPUpg  (0),
   fD0ZResKUpg  (0),
   fD0ZResPiUpg (0),
   fD0ZResEUpg (0),
   fD0RPResPUpg (0),
   fD0RPResKUpg (0),
   fD0RPResPiUpg(0),
   fD0RPResEUpg(0),
   fPt1ResPUpg  (0),
   fPt1ResKUpg  (0),
   fPt1ResPiUpg (0),
   fPt1ResEUpg (0),
   fD0ZResPCurSA  (0),
   fD0ZResKCurSA  (0),
   fD0ZResPiCurSA (0),
   fD0ZResECurSA (0),
   fD0RPResPCurSA (0),
   fD0RPResKCurSA (0),
   fD0RPResPiCurSA(0),
   fD0RPResECurSA(0),
   fPt1ResPCurSA  (0),
   fPt1ResKCurSA  (0),
   fPt1ResPiCurSA (0),
   fPt1ResECurSA (0),
   fD0ZResPUpgSA  (0),
   fD0ZResKUpgSA  (0),
   fD0ZResPiUpgSA (0),
   fD0ZResEUpgSA (0),
   fD0RPResPUpgSA (0),
   fD0RPResKUpgSA (0),
   fD0RPResPiUpgSA(0),
   fD0RPResEUpgSA(0),
   fPt1ResPUpgSA  (0),
   fPt1ResKUpgSA  (0),
   fPt1ResPiUpgSA (0),
   fPt1ResEUpgSA (0),
   // specific stuff for PbPb 2018
   fIsPbPb2018(kFALSE)
   // kFirst
  ,fD0ZResPCur_PbPb2018_kFirst(0)
  ,fD0ZResKCur_PbPb2018_kFirst(0) 
  ,fD0ZResPiCur_PbPb2018_kFirst(0) 
  ,fD0ZResECur_PbPb2018_kFirst(0) 
  ,fD0RPResPCur_PbPb2018_kFirst(0) 
  ,fD0RPResKCur_PbPb2018_kFirst(0) 
  ,fD0RPResPiCur_PbPb2018_kFirst(0) 
  ,fD0RPResECur_PbPb2018_kFirst(0)
  ,fD0RPSigmaPullRatioP_PbPb2018_kFirst(0) 
  ,fD0RPSigmaPullRatioK_PbPb2018_kFirst(0) 
  ,fD0RPSigmaPullRatioPi_PbPb2018_kFirst(0) 
  ,fD0RPSigmaPullRatioE_PbPb2018_kFirst(0)
  ,fPt1ResPCur_PbPb2018_kFirst(0)
  ,fPt1ResKCur_PbPb2018_kFirst(0)
  ,fPt1ResPiCur_PbPb2018_kFirst(0) 
  ,fPt1ResECur_PbPb2018_kFirst(0) 
  ,fD0ZResPUpg_PbPb2018_kFirst(0) 
  ,fD0ZResKUpg_PbPb2018_kFirst(0) 
  ,fD0ZResPiUpg_PbPb2018_kFirst(0) 
  ,fD0ZResEUpg_PbPb2018_kFirst(0)
  ,fD0RPResPUpg_PbPb2018_kFirst(0)
  ,fD0RPResKUpg_PbPb2018_kFirst(0) 
  ,fD0RPResPiUpg_PbPb2018_kFirst(0) 
  ,fD0RPResEUpg_PbPb2018_kFirst(0) 
  ,fPt1ResPUpg_PbPb2018_kFirst(0)
  ,fPt1ResKUpg_PbPb2018_kFirst(0) 
  ,fPt1ResPiUpg_PbPb2018_kFirst(0)
  ,fPt1ResEUpg_PbPb2018_kFirst(0) 
  ,fD0ZResPCurSA_PbPb2018_kFirst(0) 
  ,fD0ZResKCurSA_PbPb2018_kFirst(0)
  ,fD0ZResPiCurSA_PbPb2018_kFirst(0) 
  ,fD0ZResECurSA_PbPb2018_kFirst(0) 
  ,fD0RPResPCurSA_PbPb2018_kFirst(0) 
  ,fD0RPResKCurSA_PbPb2018_kFirst(0) 
  ,fD0RPResPiCurSA_PbPb2018_kFirst(0) 
  ,fD0RPResECurSA_PbPb2018_kFirst(0) 
  ,fPt1ResPCurSA_PbPb2018_kFirst(0) 
  ,fPt1ResKCurSA_PbPb2018_kFirst(0) 
  ,fPt1ResPiCurSA_PbPb2018_kFirst(0) 
  ,fPt1ResECurSA_PbPb2018_kFirst(0)
  ,fD0ZResPUpgSA_PbPb2018_kFirst(0) 
  ,fD0ZResKUpgSA_PbPb2018_kFirst(0) 
  ,fD0ZResPiUpgSA_PbPb2018_kFirst(0) 
  ,fD0ZResEUpgSA_PbPb2018_kFirst(0) 
  ,fD0RPResPUpgSA_PbPb2018_kFirst(0) 
  ,fD0RPResKUpgSA_PbPb2018_kFirst(0) 
  ,fD0RPResPiUpgSA_PbPb2018_kFirst(0) 
  ,fD0RPResEUpgSA_PbPb2018_kFirst(0) 
  ,fPt1ResPUpgSA_PbPb2018_kFirst(0)
  ,fPt1ResKUpgSA_PbPb2018_kFirst(0) 
  ,fPt1ResPiUpgSA_PbPb2018_kFirst(0)
  ,fPt1ResEUpgSA_PbPb2018_kFirst(0) 
  // kOnlySecond
  ,fD0ZResPCur_PbPb2018_kOnlySecond(0) 
  ,fD0ZResKCur_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPiCur_PbPb2018_kOnlySecond(0) 
  ,fD0ZResECur_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPCur_PbPb2018_kOnlySecond(0) 
  ,fD0RPResKCur_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPiCur_PbPb2018_kOnlySecond(0) 
  ,fD0RPResECur_PbPb2018_kOnlySecond(0) 
  ,fD0RPSigmaPullRatioP_PbPb2018_kOnlySecond(0) 
  ,fD0RPSigmaPullRatioK_PbPb2018_kOnlySecond(0) 
  ,fD0RPSigmaPullRatioPi_PbPb2018_kOnlySecond(0) 
  ,fD0RPSigmaPullRatioE_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPCur_PbPb2018_kOnlySecond(0) 
  ,fPt1ResKCur_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPiCur_PbPb2018_kOnlySecond(0)
  ,fPt1ResECur_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPUpg_PbPb2018_kOnlySecond(0) 
  ,fD0ZResKUpg_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPiUpg_PbPb2018_kOnlySecond(0) 
  ,fD0ZResEUpg_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPUpg_PbPb2018_kOnlySecond(0) 
  ,fD0RPResKUpg_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPiUpg_PbPb2018_kOnlySecond(0) 
  ,fD0RPResEUpg_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPUpg_PbPb2018_kOnlySecond(0) 
  ,fPt1ResKUpg_PbPb2018_kOnlySecond(0)  
  ,fPt1ResPiUpg_PbPb2018_kOnlySecond(0) 
  ,fPt1ResEUpg_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResKCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPiCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResECurSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResKCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPiCurSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResECurSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPCurSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResKCurSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPiCurSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResECurSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResKUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResPiUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0ZResEUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResKUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResPiUpgSA_PbPb2018_kOnlySecond(0) 
  ,fD0RPResEUpgSA_PbPb2018_kOnlySecond(0)  
  ,fPt1ResPUpgSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResKUpgSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResPiUpgSA_PbPb2018_kOnlySecond(0) 
  ,fPt1ResEUpgSA_PbPb2018_kOnlySecond(0),
   fRunInVertexing(isRunInVertexing),
   fImproveTracks(kTRUE),
   fUpdateSecVertCovMat(kTRUE),
   fUpdateSTCovMatrix(kTRUE),
   fUpdatePulls(kTRUE),
   fMimicData(kFALSE),
   fIsAOD        (kTRUE),
   fSmearOnlySignal(kFALSE),
   fMCs         (0),
   fDebugOutput (0),
   fDebugNtuple (0),
   fDebugVars   (0),
   fNDebug      (ndebug)
{
  //
  // Constructor to be used to create the task.
  // The the URIs specify the resolution files to be used. 
  // They are expected to contain TGraphs with the resolutions
  // for the current and the upgraded ITS (see code for details).
  // One may also specify for how many tracks debug information
  // is written to the output.
  //
  for(Int_t jh=0; jh<2; jh++){
    // templates of mean (original)
    for(Int_t ih=0; ih<4; ih++){
      fD0RPMeanPCur[jh][ih]=0x0;
      fD0RPMeanKCur[jh][ih]=0x0;
      fD0RPMeanPiCur[jh][ih]=0x0;
      fD0RPMeanECur[jh][ih]=0x0;
      fD0RPMeanPUpg[jh][ih]=0x0;
      fD0RPMeanKUpg[jh][ih]=0x0;
      fD0RPMeanPiUpg[jh][ih]=0x0;
      fD0RPMeanEUpg[jh][ih]=0x0;
    }
    // templates of mean for Pb-Pb 2018 
    for(UInt_t ih = 0; ih < 24; ih++)
    {
      fD0RPMeanPCur_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanKCur_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanPiCur_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanECur_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanPUpg_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanKUpg_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanPiUpg_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanEUpg_PbPb2018_kFirst[jh][ih]=0x0;
      fD0RPMeanPCur_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanKCur_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanPiCur_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanECur_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanPUpg_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanKUpg_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanPiUpg_PbPb2018_kOnlySecond[jh][ih]=0x0;
      fD0RPMeanEUpg_PbPb2018_kOnlySecond[jh][ih]=0x0;
    }
  }
  
  TString resfileCurURI = Form("alien:///alice/cern.ch/user/p/pwg_hf/common/Improver/%s/%s/ITSgraphs_Current.root",period,systematic);
  TString resfileUpgURI = Form("alien:///alice/cern.ch/user/p/pwg_hf/common/Improver/%s/%s/ITSgraphs_NewAll-X0.3-Res4um.root",period,systematic);

  printf("\n### reading file %s ...\n",resfileCurURI.Data());
  TFile *resfileCur=TFile::Open(resfileCurURI.Data());
  if(resfileCur)  printf("... READ ###\n");
  if( (resfileCurURI.Contains("LHC18r") || resfileCurURI.Contains("LHC18q")) && 
      (resfileUpgURI.Contains("LHC18r") || resfileUpgURI.Contains("LHC18q")))     fIsPbPb2018=kTRUE; 
  printf("\n\n===\n=== fIsPbPb2018: %s\n===\n\n",fIsPbPb2018?"kTRUE":"kFALSE");
  if(resfileCur) {
    if(!fIsPbPb2018){    
      if(resfileCur->Get("D0RPResP" )) {
        fD0RPResPCur =(TGraph*)(resfileCur->Get("D0RPResP" )->Clone("D0RPResPCur" ));
      }
      if(resfileCur->Get("D0RPResK" )) {
        fD0RPResKCur =(TGraph*)(resfileCur->Get("D0RPResK" )->Clone("D0RPResKCur" ));
      }
      if(resfileCur->Get("D0RPResPi")) {
        fD0RPResPiCur=(TGraph*)(resfileCur->Get("D0RPResPi")->Clone("D0RPResPiCur"));
      }
      if(resfileCur->Get("D0RPResE")) {
        fD0RPResECur=(TGraph*)(resfileCur->Get("D0RPResE")->Clone("D0RPResECur"));
      }
      if(resfileCur->Get("D0RPSigmaPullRatioP" )) {
        fD0RPSigmaPullRatioP =(TGraph*)(resfileCur->Get("D0RPSigmaPullRatioP" ));
      }
      if(resfileCur->Get("D0RPSigmaPullRatioK" )) {
        fD0RPSigmaPullRatioK =(TGraph*)(resfileCur->Get("D0RPSigmaPullRatioK" ));
      }
      if(resfileCur->Get("D0RPSigmaPullRatioPi")) {
        fD0RPSigmaPullRatioPi=(TGraph*)(resfileCur->Get("D0RPSigmaPullRatioPi"));
      }
      if(resfileCur->Get("D0RPSigmaPullRatioE")) {
        fD0RPSigmaPullRatioE=(TGraph*)(resfileCur->Get("D0RPSigmaPullRatioE"));
      }
      for(Int_t j=0; j<2; j++){
        for(Int_t i=0; i<4; i++){
	        if(resfileCur->Get(Form("D0RPMeanP_B%d_phi%d",j,i))) {
	          fD0RPMeanPCur[j][i]=(TGraph*)(resfileCur->Get(Form("D0RPMeanP_B%d_phi%d",j,i))->Clone(Form("D0RPMeanPCur_B%d_phi%d",j,i)));
	        }
	        if(resfileCur->Get(Form("D0RPMeanK_B%d_phi%d",j,i))) {
	          fD0RPMeanKCur[j][i]=(TGraph*)(resfileCur->Get(Form("D0RPMeanK_B%d_phi%d",j,i))->Clone(Form("D0RPMeanKCur_B%d_phi%d",j,i)));
	        }
	        if(resfileCur->Get(Form("D0RPMeanPi_B%d_phi%d",j,i))) {
	          fD0RPMeanPiCur[j][i]=(TGraph*)(resfileCur->Get(Form("D0RPMeanPi_B%d_phi%d",j,i))->Clone(Form("D0RPMeanPiCur_B%d_phi%d",j,i)));
	        }
	        if(resfileCur->Get(Form("D0RPMeanE_B%d_phi%d",j,i))) {
	          fD0RPMeanECur[j][i]=(TGraph*)(resfileCur->Get(Form("D0RPMeanE_B%d_phi%d",j,i))->Clone(Form("D0RPMeanECur_B%d_phi%d",j,i)));
	        }
        }
      }
      if(resfileCur->Get("D0ZResP"  )) {
        fD0ZResPCur  =(TGraph*)(resfileCur->Get("D0ZResP"  )->Clone("D0ZResPCur"  ));
      }
      if(resfileCur->Get("D0ZResK"  )) {
        fD0ZResKCur  =(TGraph*)(resfileCur->Get("D0ZResK"  )->Clone("D0ZResKCur"  ));
      }
      if(resfileCur->Get("D0ZResPi" )) {
        fD0ZResPiCur =(TGraph*)(resfileCur->Get("D0ZResPi" )->Clone("D0ZResPiCur" ));
      }
      if(resfileCur->Get("D0ZResE" )) {
        fD0ZResECur =(TGraph*)(resfileCur->Get("D0ZResE" )->Clone("D0ZResECur" ));
      }
      if(resfileCur->Get("Pt1ResP"  )) {
        fPt1ResPCur  =(TGraph*)(resfileCur->Get("Pt1ResP"  )->Clone("Pt1ResPCur"  ));
      }
      if(resfileCur->Get("Pt1ResK"  )) {
        fPt1ResKCur  =(TGraph*)(resfileCur->Get("Pt1ResK"  )->Clone("Pt1ResKCur"  ));
      }
      if(resfileCur->Get("Pt1ResPi" )) {
        fPt1ResPiCur =(TGraph*)(resfileCur->Get("Pt1ResPi" )->Clone("Pt1ResPiCur" ));
      }
      if(resfileCur->Get("Pt1ResE" )) {
        fPt1ResECur =(TGraph*)(resfileCur->Get("Pt1ResE" )->Clone("Pt1ResECur" ));
      }
      if(resfileCur->Get("D0RPResPSA" )) {
        fD0RPResPCurSA =(TGraph*)(resfileCur->Get("D0RPResPSA" )->Clone("D0RPResPCurSA" ));
      }
      if(resfileCur->Get("D0RPResKSA" )) {
        fD0RPResKCurSA =(TGraph*)(resfileCur->Get("D0RPResKSA" )->Clone("D0RPResKCurSA" ));
      }
      if(resfileCur->Get("D0RPResPiSA")) {
        fD0RPResPiCurSA=(TGraph*)(resfileCur->Get("D0RPResPiSA")->Clone("D0RPResPiCurSA"));
      }
      if(resfileCur->Get("D0RPResESA")) {
        fD0RPResECurSA=(TGraph*)(resfileCur->Get("D0RPResESA")->Clone("D0RPResECurSA"));
      }
      if(resfileCur->Get("D0ZResPSA"  )) {
        fD0ZResPCurSA  =(TGraph*)(resfileCur->Get("D0ZResPSA"  )->Clone("D0ZResPCurSA"  ));
      }
      if(resfileCur->Get("D0ZResKSA"  )) {
        fD0ZResKCurSA  =(TGraph*)(resfileCur->Get("D0ZResKSA"  )->Clone("D0ZResKCurSA"  ));
      }
      if(resfileCur->Get("D0ZResPiSA" )) {
        fD0ZResPiCurSA =(TGraph*)(resfileCur->Get("D0ZResPiSA" )->Clone("D0ZResPiCurSA" ));
      }
      if(resfileCur->Get("D0ZResESA" )) {
        fD0ZResECurSA =(TGraph*)(resfileCur->Get("D0ZResESA" )->Clone("D0ZResECurSA" ));
      }
      if(resfileCur->Get("Pt1ResPSA"  )) {
        fPt1ResPCurSA  =(TGraph*)(resfileCur->Get("Pt1ResPSA"  )->Clone("Pt1ResPCurSA"  ));
      }
      if(resfileCur->Get("Pt1ResKSA"  )) {
        fPt1ResKCurSA  =(TGraph*)(resfileCur->Get("Pt1ResKSA"  )->Clone("Pt1ResKCurSA"  ));
      }
      if(resfileCur->Get("Pt1ResPiSA" )) {
        fPt1ResPiCurSA =(TGraph*)(resfileCur->Get("Pt1ResPiSA" )->Clone("Pt1ResPiCurSA" ));
      }
      if(resfileCur->Get("Pt1ResESA" )) {
        fPt1ResECurSA =(TGraph*)(resfileCur->Get("Pt1ResESA" )->Clone("Pt1ResECurSA" ));
      }
      delete resfileCur;
    }
    else  // analysing PbPb 2018 periods
    {
      if(resfileCur->Get("kFirst_D0RPResP") && resfileCur->Get("kOnlySecond_D0RPResP")){
        fD0RPResPCur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPResP")->Clone("kFirst_D0RPResPCur"));
        fD0RPResPCur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPResP")->Clone("kOnlySecond_D0RPResPCur"));
      }
      if(resfileCur->Get("kFirst_D0RPResK") && resfileCur->Get("kOnlySecond_D0RPResK")){
        fD0RPResKCur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPResK")->Clone("kFirst_D0RPResKCur"));
        fD0RPResKCur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPResK")->Clone("kOnlySecond_D0RPResKCur"));
      }
      if(resfileCur->Get("kFirst_D0RPResPi") && resfileCur->Get("kOnlySecond_D0RPResPi")){
        fD0RPResPiCur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPResPi")->Clone("kFirst_D0RPResPiCur"));
        fD0RPResPiCur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPResPi")->Clone("kOnlySecond_D0RPResPiCur"));
      }
      if(resfileCur->Get("kFirst_D0RPResE") && resfileCur->Get("kOnlySecond_D0RPResE")){
        fD0RPResECur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPResE")->Clone("kFirst_D0RPResECur"));
        fD0RPResECur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPResE")->Clone("kOnlySecond_D0RPResECur"));
      }
      if(resfileCur->Get("kFirst_D0RPSigmaPullRatioP") && resfileCur->Get("kOnlySecond_D0RPSigmaPullRatioP")){
        fD0RPSigmaPullRatioP_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPSigmaPullRatioP"));
        fD0RPSigmaPullRatioP_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPSigmaPullRatioP"));
      }
      if(resfileCur->Get("kFirst_D0RPSigmaPullRatioK") && resfileCur->Get("kOnlySecond_D0RPSigmaPullRatioK")){
        fD0RPSigmaPullRatioK_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPSigmaPullRatioK"));
        fD0RPSigmaPullRatioK_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPSigmaPullRatioK"));
      }
      if(resfileCur->Get("kFirst_D0RPSigmaPullRatioPi") && resfileCur->Get("kOnlySecond_D0RPSigmaPullRatioPi")){
        fD0RPSigmaPullRatioPi_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPSigmaPullRatioPi"));
        fD0RPSigmaPullRatioPi_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPSigmaPullRatioPi"));
      }
      if(resfileCur->Get("kFirst_D0RPSigmaPullRatioE") && resfileCur->Get("kOnlySecond_D0RPSigmaPullRatioE")){
        fD0RPSigmaPullRatioE_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPSigmaPullRatioE"));
        fD0RPSigmaPullRatioE_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPSigmaPullRatioE"));
      }
      for(UInt_t j = 0; j < 2; j++)
      {
        for(Int_t i=0; i<24; i++){
          if(resfileCur->Get(Form("kFirst_D0RPMeanP_B%d_phi%d",j,i)) && resfileCur->Get(Form("kOnlySecond_D0RPMeanP_B%d_phi%d",j,i))){
            fD0RPMeanPCur_PbPb2018_kFirst[j][i]=(TGraph*)resfileCur->Get(Form("kFirst_D0RPMeanP_B%d_phi%d",j,i))->Clone(Form("kFirst_D0RPMeanPCur_B%d_phi%d",j,i));
            fD0RPMeanPCur_PbPb2018_kOnlySecond[j][i]=(TGraph*)resfileCur->Get(Form("kOnlySecond_D0RPMeanP_B%d_phi%d",j,i))->Clone(Form("kOnlySecond_D0RPMeanPCur_B%d_phi%d",j,i));
          }
          if(resfileCur->Get(Form("kFirst_D0RPMeanK_B%d_phi%d",j,i)) && resfileCur->Get(Form("kOnlySecond_D0RPMeanK_B%d_phi%d",j,i))){
            fD0RPMeanKCur_PbPb2018_kFirst[j][i]=(TGraph*)resfileCur->Get(Form("kFirst_D0RPMeanK_B%d_phi%d",j,i))->Clone(Form("kFirst_D0RPMeanKCur_B%d_phi%d",j,i));
            fD0RPMeanKCur_PbPb2018_kOnlySecond[j][i]=(TGraph*)resfileCur->Get(Form("kOnlySecond_D0RPMeanK_B%d_phi%d",j,i))->Clone(Form("kOnlySecond_D0RPMeanKCur_B%d_phi%d",j,i));
          }
          if(resfileCur->Get(Form("kFirst_D0RPMeanPi_B%d_phi%d",j,i)) && resfileCur->Get(Form("kOnlySecond_D0RPMeanPi_B%d_phi%d",j,i))){
            fD0RPMeanPiCur_PbPb2018_kFirst[j][i]=(TGraph*)resfileCur->Get(Form("kFirst_D0RPMeanPi_B%d_phi%d",j,i))->Clone(Form("kFirst_D0RPMeanPiCur_B%d_phi%d",j,i));
            fD0RPMeanPiCur_PbPb2018_kOnlySecond[j][i]=(TGraph*)resfileCur->Get(Form("kOnlySecond_D0RPMeanPi_B%d_phi%d",j,i))->Clone(Form("kOnlySecond_D0RPMeanPiCur_B%d_phi%d",j,i));
          }      
          if(resfileCur->Get(Form("kFirst_D0RPMeanE_B%d_phi%d",j,i)) && resfileCur->Get(Form("kOnlySecond_D0RPMeanE_B%d_phi%d",j,i))){
            fD0RPMeanECur_PbPb2018_kFirst[j][i]=(TGraph*)resfileCur->Get(Form("kFirst_D0RPMeanE_B%d_phi%d",j,i))->Clone(Form("kFirst_D0RPMeanECur_B%d_phi%d",j,i));
            fD0RPMeanECur_PbPb2018_kOnlySecond[j][i]=(TGraph*)resfileCur->Get(Form("kOnlySecond_D0RPMeanE_B%d_phi%d",j,i))->Clone(Form("kOnlySecond_D0RPMeanECur_B%d_phi%d",j,i));
          }    
        }
      }
      if(resfileCur->Get("kFirst_D0ZResP") && resfileCur->Get("kOnlySecond_D0ZResP")){
        fD0ZResPCur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0ZResP")->Clone("kFirst_D0ZResPCur"));
        fD0ZResPCur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0ZResP")->Clone("kOnlySecond_D0ZResPCur"));
      }
      if(resfileCur->Get("kFirst_D0ZResK") && resfileCur->Get("kOnlySecond_D0ZResK")){
        fD0ZResKCur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0ZResK")->Clone("kFirst_D0ZResKCur"));
        fD0ZResKCur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0ZResK")->Clone("kOnlySecond_D0ZResKCur"));
      }
      if(resfileCur->Get("kFirst_D0ZResPi") && resfileCur->Get("kOnlySecond_D0ZResPi")){
        fD0ZResPiCur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0ZResPi")->Clone("kFirst_D0ZResPiCur"));
        fD0ZResPiCur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0ZResPi")->Clone("kOnlySecond_D0ZResPiCur"));
      }
      if(resfileCur->Get("kFirst_D0ZResE") && resfileCur->Get("kOnlySecond_D0ZResE")){
        fD0ZResECur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0ZResE")->Clone("kFirst_D0ZResECur"));
        fD0ZResECur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0ZResE")->Clone("kOnlySecond_D0ZResECur"));
      }
      if(resfileCur->Get("kFirst_Pt1ResP") && resfileCur->Get("kOnlySecond_Pt1ResP")){
        fPt1ResPCur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_Pt1ResP")->Clone("kFirst_Pt1ResPCur"));
        fPt1ResPCur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_Pt1ResP")->Clone("kOnlySecond_Pt1ResPCur"));
      }
      if(resfileCur->Get("kFirst_Pt1ResK") && resfileCur->Get("kOnlySecond_Pt1ResK")){
        fPt1ResKCur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_Pt1ResK")->Clone("kFirst_Pt1ResKCur"));
        fPt1ResKCur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_Pt1ResK")->Clone("kOnlySecond_Pt1ResKCur"));
      }
      if(resfileCur->Get("kFirst_Pt1ResPi") && resfileCur->Get("kOnlySecond_Pt1ResPi")){
        fPt1ResPiCur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_Pt1ResPi")->Clone("kFirst_Pt1ResPiCur"));
        fPt1ResPiCur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_Pt1ResPi")->Clone("kOnlySecond_Pt1ResPiCur"));
      }
      if(resfileCur->Get("kFirst_Pt1ResE") && resfileCur->Get("kOnlySecond_Pt1ResE")){
        fPt1ResECur_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_Pt1ResE")->Clone("kFirst_Pt1ResECur"));
        fPt1ResECur_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_Pt1ResE")->Clone("kOnlySecond_Pt1ResECur"));
      }
      if(resfileCur->Get("kFirst_D0RPResPSA") && resfileCur->Get("kOnlySecond_D0RPResPSA")){
        fD0RPResPCurSA_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPResPSA")->Clone("kFirst_D0RPResPCurSA"));
        fD0RPResPCurSA_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPResPSA")->Clone("kOnlySecond_D0RPResPCurSA"));
      }
      if(resfileCur->Get("kFirst_D0RPResKSA") && resfileCur->Get("kOnlySecond_D0RKResPSA")){
        fD0RPResKCurSA_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPResKSA")->Clone("kFirst_D0RPResKCurSA"));
        fD0RPResKCurSA_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPResKSA")->Clone("kOnlySecond_D0RPResKCurSA"));
      }
      if(resfileCur->Get("kFirst_D0RPResPiSA") && resfileCur->Get("kOnlySecond_D0RPResPiSA")){
        fD0RPResPiCurSA_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPResPiSA")->Clone("kFirst_D0RPResPiCurSA"));
        fD0RPResPiCurSA_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPResPiSA")->Clone("kOnlySecond_D0RPResPiCurSA"));
      }
      if(resfileCur->Get("kFirst_D0RPResESA") && resfileCur->Get("kOnlySecond_D0RPResESA")){
        fD0RPResECurSA_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0RPResESA")->Clone("kFirst_D0RPResECurSA"));
        fD0RPResECurSA_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0RPResESA")->Clone("kOnlySecond_D0RPResECurSA"));
      }
      if(resfileCur->Get("kFirst_D0ZResPSA") && resfileCur->Get("kOnlySecond_D0ZResPSA")){
        fD0ZResPCurSA_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0ZResPSA")->Clone("kFirst_D0ZResPCurSA"));
        fD0ZResPCurSA_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0ZResPSA")->Clone("kOnlySecond_D0ZResPCurSA"));
      }
      if(resfileCur->Get("kFirst_D0ZResKSA") && resfileCur->Get("kOnlySecond_D0ZResKSA")){
        fD0ZResKCurSA_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0ZResKSA")->Clone("kFirst_D0ZResKCurSA"));
        fD0ZResKCurSA_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0ZResKSA")->Clone("kOnlySecond_D0ZResKCurSA"));
      }
      if(resfileCur->Get("kFirst_D0ZResPiSA") && resfileCur->Get("kOnlySecond_D0ZResPiSA")){
        fD0ZResPiCurSA_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0ZResPiSA")->Clone("kFirst_D0ZResPiCurSA"));
        fD0ZResPiCurSA_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0ZResPiSA")->Clone("kOnlySecond_D0ZResPiCurSA"));
      }
      if(resfileCur->Get("kFirst_D0ZResESA") && resfileCur->Get("kOnlySecond_D0ZResESA")){
        fD0ZResECurSA_PbPb2018_kFirst=(TGraph*)(resfileCur->Get("kFirst_D0ZResESA")->Clone("kFirst_D0ZResECurSA"));
        fD0ZResECurSA_PbPb2018_kOnlySecond=(TGraph*)(resfileCur->Get("kOnlySecond_D0ZResESA")->Clone("kOnlySecond_D0ZResECurSA"));
      }
      if(resfileCur->Get("kFirst_Pt1ResPSA") && resfileCur->Get("kOnlySecond_Pt1ResPSA")){
        fPt1ResPCurSA_PbPb2018_kFirst=(TGraph*) (resfileCur->Get("kFirst_Pt1ResPSA")->Clone("kFirst_Pt1ResPCurSA"));
        fPt1ResPCurSA_PbPb2018_kOnlySecond=(TGraph*) (resfileCur->Get("kOnlySecond_Pt1ResPSA")->Clone("kOnlySecond_Pt1ResPCurSA"));
      }
      if(resfileCur->Get("kFirst_Pt1ResKSA") && resfileCur->Get("kOnlySecond_Pt1ResKSA")){
        fPt1ResKCurSA_PbPb2018_kFirst=(TGraph*) (resfileCur->Get("kFirst_Pt1ResKSA")->Clone("kFirst_Pt1ResKCurSA"));
        fPt1ResKCurSA_PbPb2018_kOnlySecond=(TGraph*) (resfileCur->Get("kOnlySecond_Pt1ResKSA")->Clone("kOnlySecond_Pt1ResKCurSA"));
      }
      if(resfileCur->Get("kFirst_Pt1ResPiSA") && resfileCur->Get("kOnlySecond_Pt1ResPiSA")){
        fPt1ResPiCurSA_PbPb2018_kFirst=(TGraph*) (resfileCur->Get("kFirst_Pt1ResPiSA")->Clone("kFirst_Pt1ResPiCurSA"));
        fPt1ResPiCurSA_PbPb2018_kOnlySecond=(TGraph*) (resfileCur->Get("kOnlySecond_Pt1ResPiSA")->Clone("kOnlySecond_Pt1ResPiCurSA"));
      }
      if(resfileCur->Get("kFirst_Pt1ResESA") && resfileCur->Get("kOnlySecond_Pt1ResESA")){
        fPt1ResECurSA_PbPb2018_kFirst=(TGraph*) (resfileCur->Get("kFirst_Pt1ResESA")->Clone("kFirst_Pt1ResECurSA"));
        fPt1ResECurSA_PbPb2018_kOnlySecond=(TGraph*) (resfileCur->Get("kOnlySecond_Pt1ResESA")->Clone("kOnlySecond_Pt1ResECurSA"));
      }
      delete resfileCur;
    }    
  }
    
  //TString resfileUpgURI = Form("alien:///alice/cern.ch/user/p/pwg_hf/common/Improver/%s/%s/ITSgraphs_NewAll-X0.3-Res4um.root",period,systematic);
  printf("\n### reading file %s ...\n",resfileUpgURI.Data());
  TFile *resfileUpg=TFile::Open(resfileUpgURI.Data());
  if(resfileUpg)  printf("... READ ###\n");
  if(resfileUpg) {
    if(!fIsPbPb2018){ 
      if(resfileUpg->Get("D0RPResP" )) {
        fD0RPResPUpg =(TGraph*)(resfileUpg->Get("D0RPResP" )->Clone("D0RPResPUpg" ));
      }
      if(resfileUpg->Get("D0RPResK" )) {
        fD0RPResKUpg =(TGraph*)(resfileUpg->Get("D0RPResK" )->Clone("D0RPResKUpg" ));
      }
      if(resfileUpg->Get("D0RPResPi")) {
        fD0RPResPiUpg=(TGraph*)(resfileUpg->Get("D0RPResPi")->Clone("D0RPResPiUpg"));
      }
      if(resfileUpg->Get("D0RPResE")) {
        fD0RPResEUpg=(TGraph*)(resfileUpg->Get("D0RPResE")->Clone("D0RPResEUpg"));
      }
      for(Int_t j=0; j<2; j++){
        for(Int_t i=0; i<4; i++){
	        if(resfileUpg->Get(Form("D0RPMeanP_B%d_phi%d",j,i))) {
	          fD0RPMeanPUpg[j][i]=(TGraph*)(resfileUpg->Get(Form("D0RPMeanP_B%d_phi%d",j,i))->Clone(Form("D0RPMeanPUpg_B%d_phi%d",j,i)));
	        }
	        if(resfileUpg->Get(Form("D0RPMeanK_B%d_phi%d",j,i))) {
	          fD0RPMeanKUpg[j][i]=(TGraph*)(resfileUpg->Get(Form("D0RPMeanK_B%d_phi%d",j,i))->Clone(Form("D0RPMeanKUpg_B%d_phi%d",j,i)));
	        }
	        if(resfileUpg->Get(Form("D0RPMeanPi_B%d_phi%d",j,i))) {
	          fD0RPMeanPiUpg[j][i]=(TGraph*)(resfileUpg->Get(Form("D0RPMeanPi_B%d_phi%d",j,i))->Clone(Form("D0RPMeanPiUpg_B%d_phi%d",j,i)));
	        }
	        if(resfileUpg->Get(Form("D0RPMeanE_B%d_phi%d",j,i))) {
	          fD0RPMeanEUpg[j][i]=(TGraph*)(resfileUpg->Get(Form("D0RPMeanE_B%d_phi%d",j,i))->Clone(Form("D0RPMeanEUpg_B%d_phi%d",j,i)));
	        }
        }
      }
      if(resfileUpg->Get("D0ZResP"  )) {
        fD0ZResPUpg  =(TGraph*)(resfileUpg->Get("D0ZResP"  )->Clone("D0ZResPUpg"  ));
      }
      if(resfileUpg->Get("D0ZResK"  )) {
        fD0ZResKUpg  =(TGraph*)(resfileUpg->Get("D0ZResK"  )->Clone("D0ZResKUpg"  ));
      }
      if(resfileUpg->Get("D0ZResPi" )) {
        fD0ZResPiUpg =(TGraph*)(resfileUpg->Get("D0ZResPi" )->Clone("D0ZResPiUpg" ));
      }
      if(resfileUpg->Get("D0ZResE" )) {
        fD0ZResEUpg =(TGraph*)(resfileUpg->Get("D0ZResE" )->Clone("D0ZResEUpg" ));
      }
      if(resfileUpg->Get("Pt1ResP"  )) {
        fPt1ResPUpg  =(TGraph*)(resfileUpg->Get("Pt1ResP"  )->Clone("Pt1ResPUpg"  ));
      }
      if(resfileUpg->Get("Pt1ResK"  )) {
        fPt1ResKUpg  =(TGraph*)(resfileUpg->Get("Pt1ResK"  )->Clone("Pt1ResKUpg"  ));
      }
      if(resfileUpg->Get("Pt1ResPi" )) {
        fPt1ResPiUpg =(TGraph*)(resfileUpg->Get("Pt1ResPi" )->Clone("Pt1ResPiUpg" ));
      }
      if(resfileUpg->Get("Pt1ResE" )) {
        fPt1ResEUpg =(TGraph*)(resfileUpg->Get("Pt1ResE" )->Clone("Pt1ResEUpg" ));
      }
      if(resfileUpg->Get("D0RPResPSA" )) {
        fD0RPResPUpgSA =(TGraph*)(resfileUpg->Get("D0RPResPSA" )->Clone("D0RPResPUpgSA" ));
      }
      if(resfileUpg->Get("D0RPResKSA" )) {
        fD0RPResKUpgSA =(TGraph*)(resfileUpg->Get("D0RPResKSA" )->Clone("D0RPResKUpgSA" ));
      }
      if(resfileUpg->Get("D0RPResPiSA")) {
        fD0RPResPiUpgSA=(TGraph*)(resfileUpg->Get("D0RPResPiSA")->Clone("D0RPResPiUpgSA"));
      }
      if(resfileUpg->Get("D0RPResESA")) {
        fD0RPResEUpgSA=(TGraph*)(resfileUpg->Get("D0RPResESA")->Clone("D0RPResEUpgSA"));
      }
      if(resfileUpg->Get("D0ZResPSA"  )) {
        fD0ZResPUpgSA  =(TGraph*)(resfileUpg->Get("D0ZResPSA"  )->Clone("D0ZResPUpgSA"  ));
      }
      if(resfileUpg->Get("D0ZResKSA"  )) {
        fD0ZResKUpgSA  =(TGraph*)(resfileUpg->Get("D0ZResKSA"  )->Clone("D0ZResKUpgSA"  ));
      }
      if(resfileUpg->Get("D0ZResPiSA" )) {
        fD0ZResPiUpgSA =(TGraph*)(resfileUpg->Get("D0ZResPiSA" )->Clone("D0ZResPiUpgSA" ));
      }
      if(resfileUpg->Get("D0ZResESA" )) {
        fD0ZResEUpgSA =(TGraph*)(resfileUpg->Get("D0ZResESA" )->Clone("D0ZResEUpgSA" ));
      }
      if(resfileUpg->Get("Pt1ResPSA"  )) {
        fPt1ResPUpgSA  =(TGraph*)(resfileUpg->Get("Pt1ResPSA"  )->Clone("Pt1ResPUpgSA"  ));
      }
      if(resfileUpg->Get("Pt1ResKSA"  )) {
        fPt1ResKUpgSA  =(TGraph*)(resfileUpg->Get("Pt1ResKSA"  )->Clone("Pt1ResKUpgSA"  ));
      }
      if(resfileUpg->Get("Pt1ResPiSA" )) {
        fPt1ResPiUpgSA =(TGraph*)(resfileUpg->Get("Pt1ResPiSA" )->Clone("Pt1ResPiUpgSA" ));
      }
      if(resfileUpg->Get("Pt1ResESA" )) {
        fPt1ResEUpgSA =(TGraph*)(resfileUpg->Get("Pt1ResESA" )->Clone("Pt1ResEUpgSA" ));
      }
      delete resfileUpg;
    }
    else  // analysing PbPb 2018 periods
    {
      if(resfileUpg->Get("kFirst_D0RPResP") && resfileUpg->Get("kOnlySecond_D0RPResP")){
        fD0RPResPUpg_PbPb2018_kFirst=(TGraph*)(resfileUpg->Get("kFirst_D0RPResP")->Clone("kFirst_D0RPResPUpg"));
        fD0RPResPUpg_PbPb2018_kOnlySecond=(TGraph*)(resfileUpg->Get("kOnlySecond_D0RPResP")->Clone("kOnlySecond_D0RPResPUpg"));
      }
      if(resfileUpg->Get("kFirst_D0RPResK") && resfileUpg->Get("kOnlySecond_D0RPResK")){
        fD0RPResKUpg_PbPb2018_kFirst=(TGraph*)(resfileUpg->Get("kFirst_D0RPResK")->Clone("kFirst_D0RPResKUpg"));
        fD0RPResKUpg_PbPb2018_kOnlySecond=(TGraph*)(resfileUpg->Get("kOnlySecond_D0RPResK")->Clone("kOnlySecond_D0RPResKUpg"));
      }
      if(resfileUpg->Get("kFirst_D0RPResPi") && resfileUpg->Get("kOnlySecond_D0RPResPi")){
        fD0RPResPiUpg_PbPb2018_kFirst=(TGraph*)(resfileUpg->Get("kFirst_D0RPResPi")->Clone("kFirst_D0RPResPiUpg"));
        fD0RPResPiUpg_PbPb2018_kOnlySecond=(TGraph*)(resfileUpg->Get("kOnlySecond_D0RPResPi")->Clone("kOnlySecond_D0RPResPiUpg"));
      }
      if(resfileUpg->Get("kFirst_D0RPResE") && resfileUpg->Get("kOnlySecond_D0RPResE")){
        fD0RPResEUpg_PbPb2018_kFirst=(TGraph*)(resfileUpg->Get("kFirst_D0RPResE")->Clone("kFirst_D0RPResEUpg"));
        fD0RPResEUpg_PbPb2018_kOnlySecond=(TGraph*)(resfileUpg->Get("kOnlySecond_D0RPResE")->Clone("kOnlySecond_D0RPResEUpg"));
      }
      for(UInt_t j = 0; j < 2; j++){
        for(UInt_t i = 0; i < 24; i++){
          if(resfileUpg->Get(Form("kFirst_D0RPMeanP_B%d_phi%d",j,i)) && resfileUpg->Get(Form("kOnlySecond_D0RPMeanP_B%d_phi%d",j,i))){
            fD0RPMeanPUpg_PbPb2018_kFirst[j][i]=(TGraph*)(resfileUpg->Get(Form("kFirst_D0RPMeanP_B%d_phi%d",j,i))->Clone(Form("kFirst_D0RPMeanP_B%d_phi%d",j,i)));
            fD0RPMeanPUpg_PbPb2018_kOnlySecond[j][i]=(TGraph*)(resfileUpg->Get(Form("kOnlySecond_D0RPMeanP_B%d_phi%d",j,i))->Clone(Form("kOnlySecond_D0RPMeanP_B%d_phi%d",j,i)));          
          }
          if(resfileUpg->Get(Form("kFirst_D0RPMeanK_B%d_phi%d",j,i)) && resfileUpg->Get(Form("kOnlySecond_D0RPMeanK_B%d_phi%d",j,i))){
            fD0RPMeanKUpg_PbPb2018_kFirst[j][i]=(TGraph*)(resfileUpg->Get(Form("kFirst_D0RPMeanK_B%d_phi%d",j,i))->Clone(Form("kFirst_D0RPMeanK_B%d_phi%d",j,i)));
            fD0RPMeanKUpg_PbPb2018_kOnlySecond[j][i]=(TGraph*)(resfileUpg->Get(Form("kOnlySecond_D0RPMeanK_B%d_phi%d",j,i))->Clone(Form("kOnlySecond_D0RPMeanK_B%d_phi%d",j,i)));          
          }
          if(resfileUpg->Get(Form("kFirst_D0RPMeanPi_B%d_phi%d",j,i)) && resfileUpg->Get(Form("kOnlySecond_D0RPMeanPi_B%d_phi%d",j,i))){
            fD0RPMeanPiUpg_PbPb2018_kFirst[j][i]=(TGraph*)(resfileUpg->Get(Form("kFirst_D0RPMeanPi_B%d_phi%d",j,i))->Clone(Form("kFirst_D0RPMeanPi_B%d_phi%d",j,i)));
            fD0RPMeanPiUpg_PbPb2018_kOnlySecond[j][i]=(TGraph*)(resfileUpg->Get(Form("kOnlySecond_D0RPMeanPi_B%d_phi%d",j,i))->Clone(Form("kOnlySecond_D0RPMeanPi_B%d_phi%d",j,i)));          
          }
          if(resfileUpg->Get(Form("kFirst_D0RPMeanE_B%d_phi%d",j,i)) && resfileUpg->Get(Form("kOnlySecond_D0RPMeanE_B%d_phi%d",j,i))){
            fD0RPMeanEUpg_PbPb2018_kFirst[j][i]=(TGraph*)(resfileUpg->Get(Form("kFirst_D0RPMeanE_B%d_phi%d",j,i))->Clone(Form("kFirst_D0RPMeanE_B%d_phi%d",j,i)));
            fD0RPMeanEUpg_PbPb2018_kOnlySecond[j][i]=(TGraph*)(resfileUpg->Get(Form("kOnlySecond_D0RPMeanE_B%d_phi%d",j,i))->Clone(Form("kOnlySecond_D0RPMeanE_B%d_phi%d",j,i)));          
          }
        } 
      }
      if(resfileUpg->Get("kFirst_D0ZResP") && resfileUpg->Get("kOnlySecond_D0ZResP")){
        fD0ZResPUpg_PbPb2018_kFirst=(TGraph*)(resfileUpg->Get("kFirst_D0ZResP")->Clone("kFirst_D0ZResPUpg"));
        fD0ZResPUpg_PbPb2018_kOnlySecond=(TGraph*)(resfileUpg->Get("kOnlySecond_D0ZResP")->Clone("kOnlySecond_D0ZResPUpg"));
      }
      if(resfileUpg->Get("kFirst_D0ZResK") && resfileUpg->Get("kOnlySecond_D0ZResK")){
        fD0ZResKUpg_PbPb2018_kFirst=(TGraph*)(resfileUpg->Get("kFirst_D0ZResK")->Clone("kFirst_D0ZResKUpg"));
        fD0ZResKUpg_PbPb2018_kOnlySecond=(TGraph*)(resfileUpg->Get("kOnlySecond_D0ZResK")->Clone("kOnlySecond_D0ZResKUpg"));
      }
      if(resfileUpg->Get("kFirst_D0ZResPi") && resfileUpg->Get("kOnlySecond_D0ZResPi")){
        fD0ZResPiUpg_PbPb2018_kFirst=(TGraph*)(resfileUpg->Get("kFirst_D0ZResPi")->Clone("kFirst_D0ZResPiUpg"));
        fD0ZResPiUpg_PbPb2018_kOnlySecond=(TGraph*)(resfileUpg->Get("kOnlySecond_D0ZResPi")->Clone("kOnlySecond_D0ZResPiUpg"));
      }
      if(resfileUpg->Get("kFirst_D0ZResE") && resfileUpg->Get("kOnlySecond_D0ZResE")){
        fD0ZResEUpg_PbPb2018_kFirst=(TGraph*)(resfileUpg->Get("kFirst_D0ZResE")->Clone("kFirst_D0ZResEUpg"));
        fD0ZResEUpg_PbPb2018_kOnlySecond=(TGraph*)(resfileUpg->Get("kOnlySecond_D0ZResE")->Clone("kOnlySecond_D0ZResEUpg"));
      }
      if(resfileUpg->Get("kFirst_Pt1ResP") && resfileUpg->Get("kOnlySecond_Pt1ResP")){
        fPt1ResPUpg_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_Pt1ResP")->Clone("kFirst_Pt1ResPUpg"));
        fPt1ResPUpg_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_Pt1ResP")->Clone("kOnlySecond_Pt1ResPUpg"));
      }
      if(resfileUpg->Get("kFirst_Pt1ResK") && resfileUpg->Get("kOnlySecond_Pt1ResK")){
        fPt1ResKUpg_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_Pt1ResK")->Clone("kFirst_Pt1ResKUpg"));
        fPt1ResKUpg_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_Pt1ResK")->Clone("kOnlySecond_Pt1ResKUpg"));
      }
      if(resfileUpg->Get("kFirst_Pt1ResPi") && resfileUpg->Get("kOnlySecond_Pt1ResPi")){
        fPt1ResPiUpg_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_Pt1ResPi")->Clone("kFirst_Pt1ResPiUpg"));
        fPt1ResPiUpg_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_Pt1ResPi")->Clone("kOnlySecond_Pt1ResPiUpg"));
      }
      if(resfileUpg->Get("kFirst_Pt1ResE") && resfileUpg->Get("kOnlySecond_Pt1ResE")){
        fPt1ResEUpg_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_Pt1ResE")->Clone("kFirst_Pt1ResEUpg"));
        fPt1ResEUpg_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_Pt1ResE")->Clone("kOnlySecond_Pt1ResEUpg"));
      }
      if(resfileUpg->Get("kFirst_D0RPResPSA") && resfileUpg->Get("kOnlySecond_D0RPResPSA")){
        fD0RPResPUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_D0RPResPSA")->Clone("kFirst_D0RPResPUpgSA"));
        fD0RPResPUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_D0RPResPSA")->Clone("kOnlySecond_D0RPResPUpgSA"));
      }
      if(resfileUpg->Get("kFirst_D0RPResKSA") && resfileUpg->Get("kOnlySecond_D0RPResKSA")){
        fD0RPResKUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_D0RPResKSA")->Clone("kFirst_D0RPResKUpgSA"));
        fD0RPResKUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_D0RPResKSA")->Clone("kOnlySecond_D0RPResKUpgSA"));
      }
      if(resfileUpg->Get("kFirst_D0RPResPiSA") && resfileUpg->Get("kOnlySecond_D0RPResPiSA")){
        fD0RPResPiUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_D0RPResPiSA")->Clone("kFirst_D0RPResPiUpgSA"));
        fD0RPResPiUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_D0RPResPiSA")->Clone("kOnlySecond_D0RPResPiUpgSA"));
      }
      if(resfileUpg->Get("kFirst_D0RPResESA") && resfileUpg->Get("kOnlySecond_D0RPResESA")){
        fD0RPResEUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_D0RPResESA")->Clone("kFirst_D0RPResEUpgSA"));
        fD0RPResEUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_D0RPResESA")->Clone("kOnlySecond_D0RPResEUpgSA"));
      }
      if(resfileUpg->Get("kFirst_D0ZResPSA") && resfileUpg->Get("kOnlyFirst_D0ZResPSA")){
        fD0ZResPUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_D0ZResPSA")->Clone("kFirst_D0ZResPUpgSA"));
        fD0ZResPUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_D0ZResPSA")->Clone("kOnlySecond_D0ZResPUpgSA"));
      }
      if(resfileUpg->Get("kFirst_D0ZResKSA") && resfileUpg->Get("kOnlyFirst_D0ZResKSA")){
        fD0ZResKUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_D0ZResKSA")->Clone("kFirst_D0ZResKUpgSA"));
        fD0ZResKUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_D0ZResKSA")->Clone("kOnlySecond_D0ZResKUpgSA"));
      }
      if(resfileUpg->Get("kFirst_D0ZResPiSA") && resfileUpg->Get("kOnlyFirst_D0ZResPiSA")){
        fD0ZResPiUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_D0ZResPiSA")->Clone("kFirst_D0ZResPiUpgSA"));
        fD0ZResPiUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_D0ZResPiSA")->Clone("kOnlySecond_D0ZResPiUpgSA"));
      }
      if(resfileUpg->Get("kFirst_D0ZResESA") && resfileUpg->Get("kOnlyFirst_D0ZResESA")){
        fD0ZResEUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_D0ZResESA")->Clone("kFirst_D0ZResEUpgSA"));
        fD0ZResEUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_D0ZResESA")->Clone("kOnlySecond_D0ZResEUpgSA"));
      }
      if(resfileUpg->Get("kFirst_Pt1ResPSA") && resfileUpg->Get("kOnlySecond_Pt1ResPSA")){
        fPt1ResPUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_Pt1ResPSA")->Clone("kFirst_Pt1ResPUpgSA"));
        fPt1ResPUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_Pt1ResPSA")->Clone("kOnlySecond_Pt1ResPUpgSA"));
      }
      if(resfileUpg->Get("kFirst_Pt1ResKSA") && resfileUpg->Get("kOnlySecond_Pt1ResKSA")){
        fPt1ResKUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_Pt1ResKSA")->Clone("kFirst_Pt1ResKUpgSA"));
        fPt1ResKUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_Pt1ResKSA")->Clone("kOnlySecond_Pt1ResKUpgSA"));
      }
      if(resfileUpg->Get("kFirst_Pt1ResPiSA") && resfileUpg->Get("kOnlySecond_Pt1ResPiSA")){
        fPt1ResPiUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_Pt1ResPiSA")->Clone("kFirst_Pt1ResPiUpgSA"));
        fPt1ResPiUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_Pt1ResPiSA")->Clone("kOnlySecond_Pt1ResPiUpgSA"));
      }
      if(resfileUpg->Get("kFirst_Pt1ResESA") && resfileUpg->Get("kOnlySecond_Pt1ResESA")){
        fPt1ResEUpgSA_PbPb2018_kFirst=(TGraph*) (resfileUpg->Get("kFirst_Pt1ResESA")->Clone("kFirst_Pt1ResEUpgSA"));
        fPt1ResEUpgSA_PbPb2018_kOnlySecond=(TGraph*) (resfileUpg->Get("kOnlySecond_Pt1ResESA")->Clone("kOnlySecond_Pt1ResEUpgSA"));
      }
      delete resfileUpg;
    }
  
  }
    
  DefineOutput(1,TList::Class());
}

AliAnalysisTaskSEImproveITS::~AliAnalysisTaskSEImproveITS() {
  //
  // Destructor.
  //
  if (fDebugOutput) delete fDebugOutput;
}

void AliAnalysisTaskSEImproveITS::UserCreateOutputObjects() {
  //
  // Creation of user output objects.
  //
  fDebugOutput=new TList();
  fDebugOutput->SetOwner();
  fDebugOutput->SetName("debug");
  fDebugNtuple=new TNtuple("fDebugNtuple","Smearing","pdg:ptmc:d0rpo:d0zo:pt1o:sd0rpo:sd0zo:spt1o:d0rpn:d0zn:pt1n:sd0rpn:sd0zn:spt1n:d0rpmc:d0zmc:pt1mc:pullcorr:d0zoinsigma:d0zninsigma:d0rpoinsigma:d0rpninsigma");
  fDebugVars=new Float_t[fDebugNtuple->GetNvar()];
  
  fDebugOutput->Add(fDebugNtuple );

  if(fD0RPResPCur) fDebugOutput->Add(fD0RPResPCur );
  if(fD0RPResKCur) fDebugOutput->Add(fD0RPResKCur );
  if(fD0RPResPiCur) fDebugOutput->Add(fD0RPResPiCur);
  if(fD0RPResECur) fDebugOutput->Add(fD0RPResECur);
  if(fD0RPSigmaPullRatioP) fDebugOutput->Add(fD0RPSigmaPullRatioP);
  if(fD0RPSigmaPullRatioK) fDebugOutput->Add(fD0RPSigmaPullRatioK);
  if(fD0RPSigmaPullRatioPi) fDebugOutput->Add(fD0RPSigmaPullRatioPi);
  if(fD0RPSigmaPullRatioE) fDebugOutput->Add(fD0RPSigmaPullRatioE);
  for(Int_t j=0; j<2; j++){
    for(Int_t i=0; i<4; i++){
      if(fD0RPMeanPCur[j][i]) fDebugOutput->Add(fD0RPMeanPCur[j][i]);
      if(fD0RPMeanKCur[j][i]) fDebugOutput->Add(fD0RPMeanKCur[j][i]);
      if(fD0RPMeanPiCur[j][i]) fDebugOutput->Add(fD0RPMeanPiCur[j][i]);
      if(fD0RPMeanECur[j][i]) fDebugOutput->Add(fD0RPMeanECur[j][i]);
      if(fD0RPMeanPUpg[j][i]) fDebugOutput->Add(fD0RPMeanPUpg[j][i]);
      if(fD0RPMeanKUpg[j][i]) fDebugOutput->Add(fD0RPMeanKUpg[j][i]);
      if(fD0RPMeanPiUpg[j][i]) fDebugOutput->Add(fD0RPMeanPiUpg[j][i]);
      if(fD0RPMeanEUpg[j][i]) fDebugOutput->Add(fD0RPMeanEUpg[j][i]);
    }
  }
  if(fD0ZResPCur) fDebugOutput->Add(fD0ZResPCur  ); 
  if(fD0ZResKCur) fDebugOutput->Add(fD0ZResKCur  );
  if(fD0ZResPiCur) fDebugOutput->Add(fD0ZResPiCur );
  if(fD0ZResECur) fDebugOutput->Add(fD0ZResECur );
  if(fPt1ResPCur) fDebugOutput->Add(fPt1ResPCur  );
  if(fPt1ResKCur) fDebugOutput->Add(fPt1ResKCur  );
  if(fPt1ResPiCur) fDebugOutput->Add(fPt1ResPiCur );
  if(fPt1ResECur) fDebugOutput->Add(fPt1ResECur );
  if(fD0RPResPUpg) fDebugOutput->Add(fD0RPResPUpg );
  if(fD0RPResKUpg) fDebugOutput->Add(fD0RPResKUpg );
  if(fD0RPResPiUpg) fDebugOutput->Add(fD0RPResPiUpg);
  if(fD0RPResEUpg) fDebugOutput->Add(fD0RPResEUpg);
  if(fD0ZResPUpg) fDebugOutput->Add(fD0ZResPUpg  );
  if(fD0ZResKUpg) fDebugOutput->Add(fD0ZResKUpg  );
  if(fD0ZResPiUpg) fDebugOutput->Add(fD0ZResPiUpg );
  if(fD0ZResEUpg) fDebugOutput->Add(fD0ZResEUpg );
  if(fPt1ResPUpg) fDebugOutput->Add(fPt1ResPUpg  );
  if(fPt1ResKUpg) fDebugOutput->Add(fPt1ResKUpg  );
  if(fPt1ResPiUpg) fDebugOutput->Add(fPt1ResPiUpg );
  if(fPt1ResEUpg) fDebugOutput->Add(fPt1ResEUpg );

  // stuff for PbPb 2018 periods
  // Cur
  if(fD0RPResPCur_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResPCur_PbPb2018_kFirst);
  if(fD0RPResPCur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResPCur_PbPb2018_kOnlySecond);
  if(fD0RPResKCur_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResKCur_PbPb2018_kFirst);
  if(fD0RPResKCur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResKCur_PbPb2018_kOnlySecond);
  if(fD0RPResPiCur_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResPiCur_PbPb2018_kFirst);
  if(fD0RPResPiCur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResPiCur_PbPb2018_kOnlySecond);
  if(fD0RPResECur_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResECur_PbPb2018_kFirst);
  if(fD0RPResECur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResECur_PbPb2018_kOnlySecond);
  if(fD0RPSigmaPullRatioP_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPSigmaPullRatioP_PbPb2018_kFirst);
  if(fD0RPSigmaPullRatioP_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPSigmaPullRatioP_PbPb2018_kOnlySecond);
  if(fD0RPSigmaPullRatioK_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPSigmaPullRatioK_PbPb2018_kFirst);
  if(fD0RPSigmaPullRatioK_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPSigmaPullRatioK_PbPb2018_kOnlySecond);
  if(fD0RPSigmaPullRatioPi_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPSigmaPullRatioPi_PbPb2018_kFirst);
  if(fD0RPSigmaPullRatioPi_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPSigmaPullRatioPi_PbPb2018_kOnlySecond);
  if(fD0RPSigmaPullRatioE_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPSigmaPullRatioE_PbPb2018_kFirst);
  if(fD0RPSigmaPullRatioE_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPSigmaPullRatioE_PbPb2018_kOnlySecond);
  for(UInt_t j = 0; j < 2; j++){
    for(Int_t i=0; i<24; i++){
      if(fD0RPMeanPCur_PbPb2018_kFirst[j][i])  fDebugOutput->Add(fD0RPMeanPCur_PbPb2018_kFirst[j][i]);
      if(fD0RPMeanPCur_PbPb2018_kOnlySecond[j][i])  fDebugOutput->Add(fD0RPMeanPCur_PbPb2018_kOnlySecond[j][i]);
      if(fD0RPMeanKCur_PbPb2018_kFirst[j][i])  fDebugOutput->Add(fD0RPMeanKCur_PbPb2018_kFirst[j][i]);
      if(fD0RPMeanKCur_PbPb2018_kOnlySecond[j][i])  fDebugOutput->Add(fD0RPMeanKCur_PbPb2018_kOnlySecond[j][i]);
      if(fD0RPMeanPiCur_PbPb2018_kFirst[j][i])  fDebugOutput->Add(fD0RPMeanPiCur_PbPb2018_kFirst[j][i]);
      if(fD0RPMeanPiCur_PbPb2018_kOnlySecond[j][i])  fDebugOutput->Add(fD0RPMeanPiCur_PbPb2018_kOnlySecond[j][i]);
      if(fD0RPMeanECur_PbPb2018_kFirst[j][i])  fDebugOutput->Add(fD0RPMeanECur_PbPb2018_kFirst[j][i]);
      if(fD0RPMeanECur_PbPb2018_kOnlySecond[j][i])  fDebugOutput->Add(fD0RPMeanECur_PbPb2018_kOnlySecond[j][i]);
    }
  }
  if(fD0ZResPCur_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResPCur_PbPb2018_kFirst);
  if(fD0ZResPCur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResPCur_PbPb2018_kOnlySecond);
  if(fD0ZResKCur_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResKCur_PbPb2018_kFirst);
  if(fD0ZResKCur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResKCur_PbPb2018_kOnlySecond);
  if(fD0ZResPiCur_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResPiCur_PbPb2018_kFirst);
  if(fD0ZResPiCur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResPiCur_PbPb2018_kOnlySecond);
  if(fD0ZResECur_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResECur_PbPb2018_kFirst);
  if(fD0ZResECur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResECur_PbPb2018_kOnlySecond);
  if(fPt1ResPCur_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResPCur_PbPb2018_kFirst);
  if(fPt1ResPCur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResPCur_PbPb2018_kOnlySecond);
  if(fPt1ResKCur_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResKCur_PbPb2018_kFirst);
  if(fPt1ResKCur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResKCur_PbPb2018_kOnlySecond);
  if(fPt1ResPiCur_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResPiCur_PbPb2018_kFirst);
  if(fPt1ResPiCur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResPiCur_PbPb2018_kOnlySecond);
  if(fPt1ResECur_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResECur_PbPb2018_kFirst);
  if(fPt1ResECur_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResECur_PbPb2018_kOnlySecond);
  if(fD0RPResPCurSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResPCurSA_PbPb2018_kFirst);
  if(fD0RPResPCurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResPCurSA_PbPb2018_kOnlySecond);
  if(fD0RPResKCurSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResKCurSA_PbPb2018_kFirst);
  if(fD0RPResKCurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResKCurSA_PbPb2018_kOnlySecond);
  if(fD0RPResPiCurSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResPiCurSA_PbPb2018_kFirst);
  if(fD0RPResPiCurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResPiCurSA_PbPb2018_kOnlySecond);
  if(fD0RPResECurSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResECurSA_PbPb2018_kFirst);
  if(fD0RPResECurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResECurSA_PbPb2018_kOnlySecond);
  if(fD0ZResPCurSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResPCurSA_PbPb2018_kFirst);
  if(fD0ZResPCurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResPCurSA_PbPb2018_kOnlySecond);
  if(fD0ZResKCurSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResKCurSA_PbPb2018_kFirst);
  if(fD0ZResKCurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResKCurSA_PbPb2018_kOnlySecond);
  if(fD0ZResPiCurSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResPiCurSA_PbPb2018_kFirst);
  if(fD0ZResPiCurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResPiCurSA_PbPb2018_kOnlySecond);
  if(fD0ZResECurSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResECurSA_PbPb2018_kFirst);
  if(fD0ZResECurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResECurSA_PbPb2018_kOnlySecond);
  if(fPt1ResPCurSA_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResPCurSA_PbPb2018_kFirst);
  if(fPt1ResPCurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResPCurSA_PbPb2018_kOnlySecond);
  if(fPt1ResKCurSA_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResKCurSA_PbPb2018_kFirst);
  if(fPt1ResKCurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResKCurSA_PbPb2018_kOnlySecond);
  if(fPt1ResPiCurSA_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResPiCurSA_PbPb2018_kFirst);
  if(fPt1ResPiCurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResPiCurSA_PbPb2018_kOnlySecond);
  if(fPt1ResECurSA_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResECurSA_PbPb2018_kFirst);
  if(fPt1ResECurSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResECurSA_PbPb2018_kOnlySecond);
  // New
  if(fD0RPResPUpg_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResPUpg_PbPb2018_kFirst);
  if(fD0RPResPUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResPUpg_PbPb2018_kOnlySecond);
  if(fD0RPResKUpg_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResKUpg_PbPb2018_kFirst);
  if(fD0RPResKUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResKUpg_PbPb2018_kOnlySecond);
  if(fD0RPResPiUpg_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResPiUpg_PbPb2018_kFirst);
  if(fD0RPResPiUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResPiUpg_PbPb2018_kOnlySecond);
  if(fD0RPResEUpg_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResEUpg_PbPb2018_kFirst);
  if(fD0RPResEUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResEUpg_PbPb2018_kOnlySecond);
  for(UInt_t j = 0; j < 2; j++){
    for(Int_t i=0; i<24; i++){
      if(fD0RPMeanPUpg_PbPb2018_kFirst[j][i])  fDebugOutput->Add(fD0RPMeanPUpg_PbPb2018_kFirst[j][i]);
      if(fD0RPMeanPUpg_PbPb2018_kOnlySecond[j][i])  fDebugOutput->Add(fD0RPMeanPUpg_PbPb2018_kOnlySecond[j][i]);
      if(fD0RPMeanKUpg_PbPb2018_kFirst[j][i])  fDebugOutput->Add(fD0RPMeanKUpg_PbPb2018_kFirst[j][i]);
      if(fD0RPMeanKUpg_PbPb2018_kOnlySecond[j][i])  fDebugOutput->Add(fD0RPMeanKUpg_PbPb2018_kOnlySecond[j][i]);
      if(fD0RPMeanPiUpg_PbPb2018_kFirst[j][i])  fDebugOutput->Add(fD0RPMeanPiUpg_PbPb2018_kFirst[j][i]);
      if(fD0RPMeanPiUpg_PbPb2018_kOnlySecond[j][i])  fDebugOutput->Add(fD0RPMeanPiUpg_PbPb2018_kOnlySecond[j][i]);
      if(fD0RPMeanEUpg_PbPb2018_kFirst[j][i])  fDebugOutput->Add(fD0RPMeanEUpg_PbPb2018_kFirst[j][i]);
      if(fD0RPMeanEUpg_PbPb2018_kOnlySecond[j][i])  fDebugOutput->Add(fD0RPMeanEUpg_PbPb2018_kOnlySecond[j][i]);
    }
  }
  if(fD0ZResPUpg_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResPUpg_PbPb2018_kFirst);
  if(fD0ZResPUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResPUpg_PbPb2018_kOnlySecond);
  if(fD0ZResKUpg_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResKUpg_PbPb2018_kFirst);
  if(fD0ZResKUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResKUpg_PbPb2018_kOnlySecond);
  if(fD0ZResPiUpg_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResPiUpg_PbPb2018_kFirst);
  if(fD0ZResPiUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResPiUpg_PbPb2018_kOnlySecond);
  if(fD0ZResEUpg_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResEUpg_PbPb2018_kFirst);
  if(fD0ZResEUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResEUpg_PbPb2018_kOnlySecond);
  if(fPt1ResPUpg_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResPUpg_PbPb2018_kFirst);
  if(fPt1ResPUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResPUpg_PbPb2018_kOnlySecond);
  if(fPt1ResKUpg_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResKUpg_PbPb2018_kFirst);
  if(fPt1ResKUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResKUpg_PbPb2018_kOnlySecond);
  if(fPt1ResPiUpg_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResPiUpg_PbPb2018_kFirst);
  if(fPt1ResPiUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResPiUpg_PbPb2018_kOnlySecond);
  if(fPt1ResEUpg_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResEUpg_PbPb2018_kFirst);
  if(fPt1ResEUpg_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResEUpg_PbPb2018_kOnlySecond);
  if(fD0RPResPUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResPUpgSA_PbPb2018_kFirst);
  if(fD0RPResPUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResPUpgSA_PbPb2018_kOnlySecond);
  if(fD0RPResKUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResKUpgSA_PbPb2018_kFirst);
  if(fD0RPResKUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResKUpgSA_PbPb2018_kOnlySecond);
  if(fD0RPResPiUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResPiUpgSA_PbPb2018_kFirst);
  if(fD0RPResPiUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResPiUpgSA_PbPb2018_kOnlySecond);
  if(fD0RPResEUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0RPResEUpgSA_PbPb2018_kFirst);
  if(fD0RPResEUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0RPResEUpgSA_PbPb2018_kOnlySecond);
  if(fD0ZResPUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResPUpgSA_PbPb2018_kFirst);
  if(fD0ZResPUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResPUpgSA_PbPb2018_kOnlySecond);
  if(fD0ZResKUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResKUpgSA_PbPb2018_kFirst);
  if(fD0ZResKUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResKUpgSA_PbPb2018_kOnlySecond);
  if(fD0ZResPiUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResPiUpgSA_PbPb2018_kFirst);
  if(fD0ZResPiUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResPiUpgSA_PbPb2018_kOnlySecond);
  if(fD0ZResEUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fD0ZResEUpgSA_PbPb2018_kFirst);
  if(fD0ZResEUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fD0ZResEUpgSA_PbPb2018_kOnlySecond);
  if(fPt1ResPUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResPUpgSA_PbPb2018_kFirst);
  if(fPt1ResPUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResPUpgSA_PbPb2018_kOnlySecond);
  if(fPt1ResKUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResKUpgSA_PbPb2018_kFirst);
  if(fPt1ResKUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResKUpgSA_PbPb2018_kOnlySecond);
  if(fPt1ResPiUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResPiUpgSA_PbPb2018_kFirst);
  if(fPt1ResPiUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResPiUpgSA_PbPb2018_kOnlySecond);
  if(fPt1ResEUpgSA_PbPb2018_kFirst)  fDebugOutput->Add(fPt1ResEUpgSA_PbPb2018_kFirst);
  if(fPt1ResEUpgSA_PbPb2018_kOnlySecond)  fDebugOutput->Add(fPt1ResEUpgSA_PbPb2018_kOnlySecond);


  PostData(1,fDebugOutput);
}

void AliAnalysisTaskSEImproveITS::UserExec(Option_t*) {
  //
  // The event loop
  //
  AliAODEvent *ev=0x0;
  AliESDEvent *evesd=0x0;
  Double_t bz=0.;

  if(fIsAOD) {
    if(!fRunInVertexing) {
      ev=dynamic_cast<AliAODEvent*>(InputEvent());
    } else {
      if(AODEvent() && IsStandardAOD()) ev = dynamic_cast<AliAODEvent*> (AODEvent());
    }  
    if(!ev) return;
    bz=ev->GetMagneticField();
  }
  else {
    evesd = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!evesd) {
      AliError("event not found. Nothing done!");
      return;
    }
    bz=evesd->GetMagneticField();
  }

  if(fIsAOD) {
    
    fMCs=static_cast<TClonesArray*>(ev->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
    AliAODMCHeader *mcHeader = (AliAODMCHeader*)ev->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!fMCs || !mcHeader) return;

    // first loop on candidates to fill them in case of reduced AODs
    // this is done to have the same behaviour of the improver with full (pp, p-Pb) and recuced (Pb-Pb) candidates
    AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
    
    // D0->Kpi
    TClonesArray *array2Prong=static_cast<TClonesArray*>(ev->GetList()->FindObject("D0toKpi"));
    if (array2Prong) {
      for (Int_t icand=0;icand<array2Prong->GetEntriesFast();++icand) {
	AliAODRecoDecayHF2Prong *decay=static_cast<AliAODRecoDecayHF2Prong*>(array2Prong->At(icand));
	vHF->GetProng(ev,decay,0); // needed to fill fAODMap in AliAnalysisVertexingHF
	if(fSmearOnlySignal && AliVertexingHFUtils::IsCandidateInjected(decay,ev,mcHeader,fMCs)==kFALSE) continue;
	vHF->FillRecoCand(ev,(AliAODRecoDecayHF2Prong*)decay);
      }
    }
    // Dstar->Kpipi
    TClonesArray *arrayCascade=static_cast<TClonesArray*>(ev->GetList()->FindObject("Dstar"));
    if (arrayCascade) {
      for (Int_t icand=0;icand<arrayCascade->GetEntriesFast();++icand) {
	AliAODRecoCascadeHF *decayDstar=static_cast<AliAODRecoCascadeHF*>(arrayCascade->At(icand));
	vHF->GetProng(ev,decayDstar,0); // needed to fill fAODMap in AliAnalysisVertexingHF
	if(fSmearOnlySignal && AliVertexingHFUtils::IsCandidateInjected(decayDstar,ev,mcHeader,fMCs)==kFALSE) continue;
	vHF->FillRecoCasc(ev,((AliAODRecoCascadeHF*)decayDstar),kTRUE);
      }
    }
    // Three prong
    TClonesArray *array3Prong=static_cast<TClonesArray*>(ev->GetList()->FindObject("Charm3Prong"));
    if (array3Prong) {
      for (Int_t icand=0;icand<array3Prong->GetEntriesFast();++icand) {
	AliAODRecoDecayHF3Prong *decay=static_cast<AliAODRecoDecayHF3Prong*>(array3Prong->At(icand));
	vHF->GetProng(ev,decay,0); // needed to fill fAODMap in AliAnalysisVertexingHF
	if(fSmearOnlySignal && AliVertexingHFUtils::IsCandidateInjected(decay,ev,mcHeader,fMCs)==kFALSE) continue;
	vHF->FillRecoCand(ev,(AliAODRecoDecayHF3Prong*)decay);
      }
    }
  
    
    // Smear all tracks
    if (fImproveTracks) {
      for(Int_t itrack=0;itrack<ev->GetNumberOfTracks();++itrack) {
	AliAODTrack * trk = static_cast<AliAODTrack*>(ev->GetTrack(itrack));
	if(!trk) AliFatal("Not a standard AOD");
	if(fSmearOnlySignal && AliVertexingHFUtils::IsTrackInjected(trk,mcHeader,fMCs)==kFALSE) continue;
	SmearTrack(trk,bz);
      }
    }

    // TODO: recalculated primary vertex
    AliVVertex *primaryVertex=ev->GetPrimaryVertex();
    
    // Recalculate all candidates
    // D0->Kpi
    if (array2Prong) {
      for (Int_t icand=0;icand<array2Prong->GetEntriesFast();++icand) {
	AliAODRecoDecayHF2Prong *decay=static_cast<AliAODRecoDecayHF2Prong*>(array2Prong->At(icand));
	
	if(fSmearOnlySignal && AliVertexingHFUtils::IsCandidateInjected(decay,ev,mcHeader,fMCs)==kFALSE) continue;
	if(!vHF->FillRecoCand(ev,(AliAODRecoDecayHF2Prong*)decay))continue;
	
	// recalculate vertices
	AliVVertex *oldSecondaryVertex=decay->GetSecondaryVtx();
		
	AliExternalTrackParam et1; et1.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(0)));
	AliExternalTrackParam et2; et2.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(1)));
	
	TObjArray ta12;
	
	ta12.Add(&et1); ta12.Add(&et2); 
	AliESDVertex *v12 =RecalculateVertex(oldSecondaryVertex,&ta12 ,bz);
	
	
	// update secondary vertex
	Double_t pos[3];
	Double_t covpos[6];
	v12->GetXYZ(pos);
	v12->GetCovMatrix(covpos);
	decay->GetSecondaryVtx()->SetPosition(pos[0],pos[1],pos[2]);
	if(fUpdateSecVertCovMat) decay->GetSecondaryVtx()->SetCovMatrix(covpos);
	decay->GetSecondaryVtx()->SetChi2perNDF(v12->GetChi2toNDF()); 
	
	// update d0 
	Double_t d0z0[2],covd0z0[3];
	Double_t d0[2],d0err[2];
	et1.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
	d0[0]=d0z0[0];
	d0err[0] = TMath::Sqrt(covd0z0[0]);
	et2.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
	d0[1]=d0z0[0];
	d0err[1] = TMath::Sqrt(covd0z0[0]);   
	decay->Setd0Prongs(2,d0);
	decay->Setd0errProngs(2,d0err);
	// 
	
	
	Double_t xdummy=0.,ydummy=0.;
	Double_t dca;
	dca=et1.GetDCA(&et2,bz,xdummy,ydummy);
	decay->SetDCA(dca);
	
	
	
	Double_t px[2],py[2],pz[2];
	for (Int_t i=0;i<2;++i) {
	  AliExternalTrackParam et;
	  et.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(i)));
	  et.PropagateToDCA(v12,bz,100.,d0z0,covd0z0);
	  px[i]=et.Px();
	  py[i]=et.Py();
	  pz[i]=et.Pz();
	}
	decay->SetPxPyPzProngs(2,px,py,pz);
	delete v12;
      }
    }
    
    
    // Dstar->Kpipi
    if (arrayCascade) {
      for (Int_t icand=0;icand<arrayCascade->GetEntriesFast();++icand) {
	AliAODRecoCascadeHF *decayDstar=static_cast<AliAODRecoCascadeHF*>(arrayCascade->At(icand));
	if(fSmearOnlySignal && AliVertexingHFUtils::IsCandidateInjected(decayDstar,ev,mcHeader,fMCs)==kFALSE) continue;
	if(!vHF->FillRecoCasc(ev,((AliAODRecoCascadeHF*)decayDstar),kTRUE))continue;
	//Get D0 from D*
	AliAODRecoDecayHF2Prong* decay=(AliAODRecoDecayHF2Prong*)decayDstar->Get2Prong();
	
	// recalculate vertices
	//AliVVertex *oldSecondaryVertex=decay->GetSecondaryVtx();
	
	//soft pion
	AliExternalTrackParam et3; et3.CopyFromVTrack(static_cast<AliAODTrack*>(decayDstar->GetBachelor()));
	
	//track D0
	AliNeutralTrackParam *trackD0 = new AliNeutralTrackParam(decay);
	
	//!!!!TODO: covariance matrix
	
	// update d0 
	Double_t d0z0[2],covd0z0[3];
	Double_t d01[2],d01err[2];
	
	//the D*
	et3.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
	d01[0]=d0z0[0];
	d01err[0] = TMath::Sqrt(covd0z0[0]); 
	trackD0->PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
	d01[1]=d0z0[0];
	d01err[1] = TMath::Sqrt(covd0z0[0]);  
	decayDstar->Setd0Prongs(2,d01);
	decayDstar->Setd0errProngs(2,d01err);
        
	// delete v12;
	delete trackD0; trackD0=NULL;
	
	// a run for D*
	Double_t px1[2],py1[2],pz1[2];
	for (Int_t i=0;i<2;++i) {
	  const AliAODTrack *t1=static_cast<AliAODTrack*>(decayDstar->GetDaughter(i));
	  px1[i]=t1->Px();
	  py1[i]=t1->Py();
	  pz1[i]=t1->Pz();
	}
	decayDstar->SetPxPyPzProngs(2,px1,py1,pz1);
	
      }
    }
    
    
    // Three prong
    if (array3Prong) {
      for (Int_t icand=0;icand<array3Prong->GetEntriesFast();++icand) {
	AliAODRecoDecayHF3Prong *decay=static_cast<AliAODRecoDecayHF3Prong*>(array3Prong->At(icand));
	if(fSmearOnlySignal && AliVertexingHFUtils::IsCandidateInjected(decay,ev,mcHeader,fMCs)==kFALSE) continue;
	if(!vHF->FillRecoCand(ev,(AliAODRecoDecayHF3Prong*)decay))continue;
	
	// recalculate vertices
	AliVVertex *oldSecondaryVertex=decay->GetSecondaryVtx();
	AliExternalTrackParam et1; et1.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(0)));
	AliExternalTrackParam et2; et2.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(1)));
	AliExternalTrackParam et3; et3.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(2)));
	TObjArray ta123,ta12,ta23;
	ta123.Add(&et1);ta123.Add(&et2);ta123.Add(&et3);
	ta12. Add(&et1);ta12 .Add(&et2);
	ta23 .Add(&et2);ta23 .Add(&et3);
	AliESDVertex *v123=RecalculateVertex(oldSecondaryVertex,&ta123,bz);
	AliESDVertex *v12 =RecalculateVertex(oldSecondaryVertex,&ta12 ,bz);
	AliESDVertex *v23 =RecalculateVertex(oldSecondaryVertex,&ta23 ,bz);
	
	// update secondary vertex
	Double_t pos[3];
	Double_t covpos[6];
	v123->GetXYZ(pos);
	v123->GetCovMatrix(covpos);
	decay->GetSecondaryVtx()->SetPosition(pos[0],pos[1],pos[2]);
	if(fUpdateSecVertCovMat) decay->GetSecondaryVtx()->SetCovMatrix(covpos);
	decay->GetSecondaryVtx()->SetChi2perNDF(v123->GetChi2toNDF()); 
	
	// update d0 for all progs
	Double_t d0z0[2],covd0z0[3];
	Double_t d0[3],d0err[3];
	et1.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
	d0[0]=d0z0[0];
	d0err[0] = TMath::Sqrt(covd0z0[0]);
	et2.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
	d0[1]=d0z0[0];
	d0err[1] = TMath::Sqrt(covd0z0[0]);
	et3.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
	d0[2]=d0z0[0];
	d0err[2] = TMath::Sqrt(covd0z0[0]);
	decay->Setd0Prongs   (3,d0   );
	decay->Setd0errProngs(3,d0err);
	// TODO: setter missing
	
	// update dca for prong combinations
	Double_t xdummy=0.,ydummy=0.;
	Double_t dca[3];
	dca[0]=et1.GetDCA(&et2,bz,xdummy,ydummy);
	dca[1]=et3.GetDCA(&et2,bz,xdummy,ydummy);
	dca[2]=et1.GetDCA(&et3,bz,xdummy,ydummy);
	decay->SetDCAs(3,dca);
	
	// update sigmavertex = dispersion
	Float_t sigmaV=v123->GetDispersion(); 
	decay->SetSigmaVert(sigmaV); 
	// update dist12 and dist23
	primaryVertex->GetXYZ(pos);
	decay->SetDist12toPrim(TMath::Sqrt((v12->GetX()-pos[0])*(v12->GetX()-pos[0])
					   +(v12->GetY()-pos[1])*(v12->GetY()-pos[1])
					   +(v12->GetZ()-pos[2])*(v12->GetZ()-pos[2])));
	decay->SetDist23toPrim(TMath::Sqrt((v23->GetX()-pos[0])*(v23->GetX()-pos[0])
					   +(v23->GetY()-pos[1])*(v23->GetY()-pos[1])
					   +(v23->GetZ()-pos[2])*(v23->GetZ()-pos[2])));
	
	
	Double_t px[3],py[3],pz[3];
	for (Int_t i=0;i<3;++i) {
	  AliExternalTrackParam et;
	  et.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(i)));
	  et.PropagateToDCA(v123,bz,100.,d0z0,covd0z0);
	  px[i]=et.Px();
	  py[i]=et.Py();
	  pz[i]=et.Pz();
	}
	decay->SetPxPyPzProngs(3,px,py,pz);
	
	delete v123;delete v12;delete v23;
      }
    }
    delete vHF;
    
  } // end AOD
  else {
    //
    // In case of ESD: only smear all tracks
    //
    if (!fMCEvent) return;
    if (fImproveTracks) {
      for(Int_t itrack=0;itrack<evesd->GetNumberOfTracks();++itrack) {
	AliESDtrack * trk = static_cast<AliESDtrack*>(evesd->GetTrack(itrack));
	if(!trk) AliFatal("No a standard ESD");
	SmearTrack(trk,bz);
      }
    }    
  }// end ESD
  
}

void AliAnalysisTaskSEImproveITS::SmearTrack(AliVTrack *track,Double_t bz) {

  // flags for PbPb 2018 
  Bool_t is_kFirst_Trk      = kFALSE;
  Bool_t is_kOnlySecond_Trk = kFALSE;

  // Early exit, if this track has nothing in common with the ITS
  if(!fIsPbPb2018){
    if (!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)))  return;
  }
  else{ // PbPb 2018 analysed
    // kFirst tracks
    if(track->HasPointOnITSLayer(0))  is_kFirst_Trk = kTRUE;
    // kOnlySecond tracks  
    else if( !(track->HasPointOnITSLayer(0)) && track->HasPointOnITSLayer(1) )  is_kOnlySecond_Trk = kTRUE; 
    if( !is_kFirst_Trk && !is_kOnlySecond_Trk ) return;
  }
  

  // Check if the track was already "improved" (this is done with a trick using layer 7 (ie the 8th))
  if (TESTBIT(track->GetITSClusterMap(),7)) return;
  //

  // Get reconstructed track parameters
  AliExternalTrackParam et; et.CopyFromVTrack(track);
  Double_t *param=const_cast<Double_t*>(et.GetParameter());


  Double_t *covar=const_cast<Double_t*>(et.GetCovariance());

  // Get MC info
  Int_t imc=track->GetLabel();
  if (imc<=0) return;
  Double_t mcx[3];
  Double_t mcp[3];
  Double_t mccv[36]={0.};
  Short_t  mcc;
  const AliVParticle *mc= 0x0;
  if(fIsAOD) {
    mc = static_cast<AliVParticle*>(fMCs->At(imc));
  }
  else {
    mc = static_cast<AliVParticle*>(fMCEvent->GetTrack(imc));
  }
  if(!mc) return;
  mc->XvYvZv(mcx);
  mc->PxPyPz(mcp);
  mcc=mc->Charge();
  AliExternalTrackParam mct(mcx,mcp,mccv,mcc);
  const Double_t *parammc=mct.GetParameter();
//TODO:  const Double_t *covermc=mct.GetCovariance();
  AliVertex vtx(mcx,1.,1);

  // Correct reference points and frames according to MC
  // TODO: B-Field correct?
  // TODO: failing propagation....
  et.PropagateToDCA(&vtx,bz,10.);
  et.Rotate(mct.GetAlpha());

  // Select appropriate smearing functions
  Double_t ptmc=TMath::Abs(mc->Pt());
  Double_t phimc=mc->Phi();
  Int_t phiBin=PhiBin(phimc);
  Int_t magfield=0;
  if(bz<0.) magfield=0;
  else if(bz>0.)magfield=1;
  Double_t sd0rpn=0.;
  Double_t sd0mrpn=0.;
  Double_t sd0zn =0.;
  Double_t spt1n =0.;
  Double_t sd0rpo=0.;
  Double_t sd0mrpo=0.;
  Double_t sd0zo =0.;
  Double_t spt1o =0.;
  Double_t pullcorr=1.;
  Int_t pdgcode = 0;
  if(fIsAOD){
    const AliAODMCParticle *amcpart = static_cast<const AliAODMCParticle *>(mc);
    if(!amcpart) return;
    pdgcode = amcpart->GetPdgCode();
  } else {
    const AliMCParticle *emcpart = static_cast<const AliMCParticle *>(mc);
    if(!emcpart) return;
    pdgcode = emcpart->PdgCode();
  }
    
  if(!fIsPbPb2018){  
    switch (pdgcode) {
    case 2212: case -2212:
      sd0rpo=EvalGraph(ptmc,fD0RPResPCur,fD0RPResPCurSA);
      sd0zo =EvalGraph(ptmc,fD0ZResPCur,fD0ZResPCurSA);
      spt1o =EvalGraph(ptmc,fPt1ResPCur,fPt1ResPCurSA);
      sd0rpn=EvalGraph(ptmc,fD0RPResPUpg,fD0RPResPUpgSA);
      sd0zn =EvalGraph(ptmc,fD0ZResPUpg,fD0ZResPUpgSA);
      spt1n =EvalGraph(ptmc,fPt1ResPUpg,fPt1ResPUpgSA);
      sd0mrpo=EvalGraph(ptmc,fD0RPMeanPCur[magfield][phiBin],fD0RPMeanPCur[magfield][phiBin]);
      sd0mrpn=EvalGraph(ptmc,fD0RPMeanPUpg[magfield][phiBin],fD0RPMeanPUpg[magfield][phiBin]);
      if(fD0RPSigmaPullRatioP) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioP,fD0RPSigmaPullRatioP);
      break;
    case 321: case -321:
      sd0rpo=EvalGraph(ptmc,fD0RPResKCur,fD0RPResKCurSA);
      sd0zo =EvalGraph(ptmc,fD0ZResKCur,fD0ZResKCurSA);
      spt1o =EvalGraph(ptmc,fPt1ResKCur,fPt1ResKCurSA);
      sd0rpn=EvalGraph(ptmc,fD0RPResKUpg,fD0RPResKUpgSA);
      sd0zn =EvalGraph(ptmc,fD0ZResKUpg,fD0ZResKUpgSA);
      spt1n =EvalGraph(ptmc,fPt1ResKUpg,fPt1ResKUpgSA);
      sd0mrpo=EvalGraph(ptmc,fD0RPMeanKCur[magfield][phiBin],fD0RPMeanKCur[magfield][phiBin]);
      sd0mrpn=EvalGraph(ptmc,fD0RPMeanKUpg[magfield][phiBin],fD0RPMeanKUpg[magfield][phiBin]);
      if(fD0RPSigmaPullRatioK) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioK,fD0RPSigmaPullRatioK);
      break;
    case 211: case -211:
      sd0rpo=EvalGraph(ptmc,fD0RPResPiCur,fD0RPResPiCurSA);
      sd0zo =EvalGraph(ptmc,fD0ZResPiCur,fD0ZResPiCurSA);
      spt1o =EvalGraph(ptmc,fPt1ResPiCur,fPt1ResPiCurSA);
      sd0rpn=EvalGraph(ptmc,fD0RPResPiUpg,fD0RPResPiUpgSA);
      sd0zn =EvalGraph(ptmc,fD0ZResPiUpg,fD0ZResPiUpgSA);
      spt1n =EvalGraph(ptmc,fPt1ResPiUpg,fPt1ResPiUpgSA);
      sd0mrpo=EvalGraph(ptmc,fD0RPMeanPiCur[magfield][phiBin],fD0RPMeanPiCur[magfield][phiBin]);
      sd0mrpn=EvalGraph(ptmc,fD0RPMeanPiUpg[magfield][phiBin],fD0RPMeanPiUpg[magfield][phiBin]);
      if(fD0RPSigmaPullRatioPi) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioPi,fD0RPSigmaPullRatioPi);
      break;
    case 11: case -11:
      sd0rpo=EvalGraph(ptmc,fD0RPResECur,fD0RPResECurSA);
      sd0zo =EvalGraph(ptmc,fD0ZResECur,fD0ZResECurSA);
      spt1o =EvalGraph(ptmc,fPt1ResECur,fPt1ResECurSA);
      sd0rpn=EvalGraph(ptmc,fD0RPResEUpg,fD0RPResEUpgSA);
      sd0zn =EvalGraph(ptmc,fD0ZResEUpg,fD0ZResEUpgSA);
      spt1n =EvalGraph(ptmc,fPt1ResEUpg,fPt1ResEUpgSA);
      sd0mrpo=EvalGraph(ptmc,fD0RPMeanECur[magfield][phiBin],fD0RPMeanECur[magfield][phiBin]);
      sd0mrpn=EvalGraph(ptmc,fD0RPMeanEUpg[magfield][phiBin],fD0RPMeanEUpg[magfield][phiBin]);
      if(fD0RPSigmaPullRatioE) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioE,fD0RPSigmaPullRatioE);
      break;
    default:
      return;
    }
  }
  else{ // PbPb 2018 periods templates
    // kFirst tracks
    if(is_kFirst_Trk){
      switch (pdgcode) {
      case 2212: case -2212:
        sd0rpo=EvalGraph(ptmc,fD0RPResPCur_PbPb2018_kFirst,fD0RPResPCurSA_PbPb2018_kFirst);
        sd0zo =EvalGraph(ptmc,fD0ZResPCur_PbPb2018_kFirst,fD0ZResPCurSA_PbPb2018_kFirst);
        spt1o =EvalGraph(ptmc,fPt1ResPCur_PbPb2018_kFirst,fPt1ResPCurSA_PbPb2018_kFirst);
        sd0rpn=EvalGraph(ptmc,fD0RPResPUpg_PbPb2018_kFirst,fD0RPResPUpgSA_PbPb2018_kFirst);
        sd0zn =EvalGraph(ptmc,fD0ZResPUpg_PbPb2018_kFirst,fD0ZResPUpgSA_PbPb2018_kFirst);
        spt1n =EvalGraph(ptmc,fPt1ResPUpg_PbPb2018_kFirst,fPt1ResPUpgSA_PbPb2018_kFirst);
        sd0mrpo=EvalGraph(ptmc,fD0RPMeanKCur_PbPb2018_kFirst[magfield][phiBin],fD0RPMeanKCur_PbPb2018_kFirst[magfield][phiBin]);
        sd0mrpn=EvalGraph(ptmc,fD0RPMeanKUpg_PbPb2018_kFirst[magfield][phiBin],fD0RPMeanKUpg_PbPb2018_kFirst[magfield][phiBin]);
        if(fD0RPSigmaPullRatioP_PbPb2018_kFirst) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioP_PbPb2018_kFirst,fD0RPSigmaPullRatioP_PbPb2018_kFirst);
        break;
      case 321: case -321:
        sd0rpo=EvalGraph(ptmc,fD0RPResKCur_PbPb2018_kFirst,fD0RPResKCurSA_PbPb2018_kFirst);
        sd0zo =EvalGraph(ptmc,fD0ZResKCur_PbPb2018_kFirst,fD0ZResKCurSA_PbPb2018_kFirst);
        spt1o =EvalGraph(ptmc,fPt1ResKCur_PbPb2018_kFirst,fPt1ResKCurSA_PbPb2018_kFirst);
        sd0rpn=EvalGraph(ptmc,fD0RPResKUpg_PbPb2018_kFirst,fD0RPResKUpgSA_PbPb2018_kFirst);
        sd0zn =EvalGraph(ptmc,fD0ZResKUpg_PbPb2018_kFirst,fD0ZResKUpgSA_PbPb2018_kFirst);
        spt1n =EvalGraph(ptmc,fPt1ResKUpg_PbPb2018_kFirst,fPt1ResKUpgSA_PbPb2018_kFirst);
        sd0mrpo=EvalGraph(ptmc,fD0RPMeanKCur_PbPb2018_kFirst[magfield][phiBin],fD0RPMeanKCur_PbPb2018_kFirst[magfield][phiBin]);
        sd0mrpn=EvalGraph(ptmc,fD0RPMeanKUpg_PbPb2018_kFirst[magfield][phiBin],fD0RPMeanKUpg_PbPb2018_kFirst[magfield][phiBin]);
        if(fD0RPSigmaPullRatioK_PbPb2018_kFirst) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioK_PbPb2018_kFirst,fD0RPSigmaPullRatioK_PbPb2018_kFirst);
        break;
      case 211: case -211:
        sd0rpo=EvalGraph(ptmc,fD0RPResPiCur_PbPb2018_kFirst,fD0RPResPiCurSA_PbPb2018_kFirst);
        sd0zo =EvalGraph(ptmc,fD0ZResPiCur_PbPb2018_kFirst,fD0ZResPiCurSA_PbPb2018_kFirst);
        spt1o =EvalGraph(ptmc,fPt1ResPiCur_PbPb2018_kFirst,fPt1ResPiCurSA_PbPb2018_kFirst);
        sd0rpn=EvalGraph(ptmc,fD0RPResPiUpg_PbPb2018_kFirst,fD0RPResPiUpgSA_PbPb2018_kFirst);
        sd0zn =EvalGraph(ptmc,fD0ZResPiUpg_PbPb2018_kFirst,fD0ZResPiUpgSA_PbPb2018_kFirst);
        spt1n =EvalGraph(ptmc,fPt1ResPiUpg_PbPb2018_kFirst,fPt1ResPiUpgSA_PbPb2018_kFirst); 
        sd0mrpo=EvalGraph(ptmc,fD0RPMeanPiCur_PbPb2018_kFirst[magfield][phiBin],fD0RPMeanPiCur_PbPb2018_kFirst[magfield][phiBin]);
        sd0mrpn=EvalGraph(ptmc,fD0RPMeanPiUpg_PbPb2018_kFirst[magfield][phiBin],fD0RPMeanPiUpg_PbPb2018_kFirst[magfield][phiBin]);
        if(fD0RPSigmaPullRatioPi_PbPb2018_kFirst) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioPi_PbPb2018_kFirst,fD0RPSigmaPullRatioPi_PbPb2018_kFirst);
        break;
      case 11: case -11:
        sd0rpo=EvalGraph(ptmc,fD0RPResECur_PbPb2018_kFirst,fD0RPResECurSA_PbPb2018_kFirst);
        sd0zo =EvalGraph(ptmc,fD0ZResECur_PbPb2018_kFirst,fD0ZResECurSA_PbPb2018_kFirst);
        spt1o =EvalGraph(ptmc,fPt1ResECur_PbPb2018_kFirst,fPt1ResECurSA_PbPb2018_kFirst);
        sd0rpn=EvalGraph(ptmc,fD0RPResEUpg_PbPb2018_kFirst,fD0RPResEUpgSA_PbPb2018_kFirst);
        sd0zn =EvalGraph(ptmc,fD0ZResEUpg_PbPb2018_kFirst,fD0ZResEUpgSA_PbPb2018_kFirst);
        spt1n =EvalGraph(ptmc,fPt1ResEUpg_PbPb2018_kFirst,fPt1ResEUpgSA_PbPb2018_kFirst);
        sd0mrpo=EvalGraph(ptmc,fD0RPMeanECur_PbPb2018_kFirst[magfield][phiBin],fD0RPMeanECur_PbPb2018_kFirst[magfield][phiBin]);
        sd0mrpn=EvalGraph(ptmc,fD0RPMeanEUpg_PbPb2018_kFirst[magfield][phiBin],fD0RPMeanEUpg_PbPb2018_kFirst[magfield][phiBin]);
        if(fD0RPSigmaPullRatioE_PbPb2018_kFirst) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioE_PbPb2018_kFirst,fD0RPSigmaPullRatioE_PbPb2018_kFirst);
        break;
      default:
        return;
      }
    }
    // kOnlySecond tracks
    else if(is_kOnlySecond_Trk){
      switch (pdgcode) {
      case 2212: case -2212:
        sd0rpo=EvalGraph(ptmc,fD0RPResPCur_PbPb2018_kOnlySecond,fD0RPResPCurSA_PbPb2018_kOnlySecond);
        sd0zo =EvalGraph(ptmc,fD0ZResPCur_PbPb2018_kOnlySecond,fD0ZResPCurSA_PbPb2018_kOnlySecond);
        spt1o =EvalGraph(ptmc,fPt1ResPCur_PbPb2018_kOnlySecond,fPt1ResPCurSA_PbPb2018_kOnlySecond);
        sd0rpn=EvalGraph(ptmc,fD0RPResPUpg_PbPb2018_kOnlySecond,fD0RPResPUpgSA_PbPb2018_kOnlySecond);
        sd0zn =EvalGraph(ptmc,fD0ZResPUpg_PbPb2018_kOnlySecond,fD0ZResPUpgSA_PbPb2018_kOnlySecond);
        spt1n =EvalGraph(ptmc,fPt1ResPUpg_PbPb2018_kOnlySecond,fPt1ResPUpgSA_PbPb2018_kOnlySecond);
        sd0mrpo=EvalGraph(ptmc,fD0RPMeanKCur_PbPb2018_kOnlySecond[magfield][phiBin],fD0RPMeanKCur_PbPb2018_kOnlySecond[magfield][phiBin]);
        sd0mrpn=EvalGraph(ptmc,fD0RPMeanKUpg_PbPb2018_kOnlySecond[magfield][phiBin],fD0RPMeanKUpg_PbPb2018_kOnlySecond[magfield][phiBin]);
        if(fD0RPSigmaPullRatioP_PbPb2018_kOnlySecond) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioP_PbPb2018_kOnlySecond,fD0RPSigmaPullRatioP_PbPb2018_kOnlySecond);
        break;
      case 321: case -321:
        sd0rpo=EvalGraph(ptmc,fD0RPResKCur_PbPb2018_kOnlySecond,fD0RPResKCurSA_PbPb2018_kOnlySecond);
        sd0zo =EvalGraph(ptmc,fD0ZResKCur_PbPb2018_kOnlySecond,fD0ZResKCurSA_PbPb2018_kOnlySecond);
        spt1o =EvalGraph(ptmc,fPt1ResKCur_PbPb2018_kOnlySecond,fPt1ResKCurSA_PbPb2018_kOnlySecond);
        sd0rpn=EvalGraph(ptmc,fD0RPResKUpg_PbPb2018_kOnlySecond,fD0RPResKUpgSA_PbPb2018_kOnlySecond);
        sd0zn =EvalGraph(ptmc,fD0ZResKUpg_PbPb2018_kOnlySecond,fD0ZResKUpgSA_PbPb2018_kOnlySecond);
        spt1n =EvalGraph(ptmc,fPt1ResKUpg_PbPb2018_kOnlySecond,fPt1ResKUpgSA_PbPb2018_kOnlySecond);
        sd0mrpo=EvalGraph(ptmc,fD0RPMeanKCur_PbPb2018_kOnlySecond[magfield][phiBin],fD0RPMeanKCur_PbPb2018_kOnlySecond[magfield][phiBin]);
        sd0mrpn=EvalGraph(ptmc,fD0RPMeanKUpg_PbPb2018_kOnlySecond[magfield][phiBin],fD0RPMeanKUpg_PbPb2018_kOnlySecond[magfield][phiBin]);
        if(fD0RPSigmaPullRatioK_PbPb2018_kOnlySecond) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioK_PbPb2018_kOnlySecond,fD0RPSigmaPullRatioK_PbPb2018_kOnlySecond);
        break;
      case 211: case -211:
        sd0rpo=EvalGraph(ptmc,fD0RPResPiCur_PbPb2018_kOnlySecond,fD0RPResPiCurSA_PbPb2018_kOnlySecond);
        sd0zo =EvalGraph(ptmc,fD0ZResPiCur_PbPb2018_kOnlySecond,fD0ZResPiCurSA_PbPb2018_kOnlySecond);
        spt1o =EvalGraph(ptmc,fPt1ResPiCur_PbPb2018_kOnlySecond,fPt1ResPiCurSA_PbPb2018_kOnlySecond);
        sd0rpn=EvalGraph(ptmc,fD0RPResPiUpg_PbPb2018_kOnlySecond,fD0RPResPiUpgSA_PbPb2018_kOnlySecond);
        sd0zn =EvalGraph(ptmc,fD0ZResPiUpg_PbPb2018_kOnlySecond,fD0ZResPiUpgSA_PbPb2018_kOnlySecond);
        spt1n =EvalGraph(ptmc,fPt1ResPiUpg_PbPb2018_kOnlySecond,fPt1ResPiUpgSA_PbPb2018_kOnlySecond); 
        sd0mrpo=EvalGraph(ptmc,fD0RPMeanPiCur_PbPb2018_kOnlySecond[magfield][phiBin],fD0RPMeanPiCur_PbPb2018_kOnlySecond[magfield][phiBin]);
        sd0mrpn=EvalGraph(ptmc,fD0RPMeanPiUpg_PbPb2018_kOnlySecond[magfield][phiBin],fD0RPMeanPiUpg_PbPb2018_kOnlySecond[magfield][phiBin]);
        if(fD0RPSigmaPullRatioPi_PbPb2018_kOnlySecond) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioPi_PbPb2018_kOnlySecond,fD0RPSigmaPullRatioPi_PbPb2018_kOnlySecond);
        break;
      case 11: case -11:
        sd0rpo=EvalGraph(ptmc,fD0RPResECur_PbPb2018_kOnlySecond,fD0RPResECurSA_PbPb2018_kOnlySecond);
        sd0zo =EvalGraph(ptmc,fD0ZResECur_PbPb2018_kOnlySecond,fD0ZResECurSA_PbPb2018_kOnlySecond);
        spt1o =EvalGraph(ptmc,fPt1ResECur_PbPb2018_kOnlySecond,fPt1ResECurSA_PbPb2018_kOnlySecond);
        sd0rpn=EvalGraph(ptmc,fD0RPResEUpg_PbPb2018_kOnlySecond,fD0RPResEUpgSA_PbPb2018_kOnlySecond);
        sd0zn =EvalGraph(ptmc,fD0ZResEUpg_PbPb2018_kOnlySecond,fD0ZResEUpgSA_PbPb2018_kOnlySecond);
        spt1n =EvalGraph(ptmc,fPt1ResEUpg_PbPb2018_kOnlySecond,fPt1ResEUpgSA_PbPb2018_kOnlySecond);
        sd0mrpo=EvalGraph(ptmc,fD0RPMeanECur_PbPb2018_kOnlySecond[magfield][phiBin],fD0RPMeanECur_PbPb2018_kOnlySecond[magfield][phiBin]);
        sd0mrpn=EvalGraph(ptmc,fD0RPMeanEUpg_PbPb2018_kOnlySecond[magfield][phiBin],fD0RPMeanEUpg_PbPb2018_kOnlySecond[magfield][phiBin]);
        if(fD0RPSigmaPullRatioE_PbPb2018_kOnlySecond) pullcorr=EvalGraph(ptmc,fD0RPSigmaPullRatioE_PbPb2018_kOnlySecond,fD0RPSigmaPullRatioE_PbPb2018_kOnlySecond);
        break;
      default:
        return;
      }
    }
  }
  

  // Use the same units (i.e. cm and GeV/c)! TODO: pt!
  sd0rpo*=1.e-4;
  sd0zo *=1.e-4;
  sd0rpn*=1.e-4;
  sd0zn *=1.e-4;
  sd0mrpo*=1.e-4;
  sd0mrpn*=1.e-4;

  // Apply the smearing
  Double_t d0zo  =param  [1];
  Double_t d0zmc =parammc[1];
  Double_t d0rpo =param  [0];
  Double_t d0rpmc=parammc[0];
  Double_t pt1o  =param  [4];
  Double_t pt1mc =parammc[4];
  Double_t dd0zo =d0zo-d0zmc;
  Double_t dd0zn =dd0zo *(sd0zo >0. ? (sd0zn /sd0zo ) : 1.);
  Double_t d0zn  =d0zmc+dd0zn;
  Double_t dd0rpo=d0rpo-d0rpmc;
  Double_t dd0rpn=dd0rpo*(sd0rpo>0. ? (sd0rpn/sd0rpo) : 1.);
  Double_t dd0mrpn=TMath::Abs(sd0mrpn)-TMath::Abs(sd0mrpo);
  Double_t d0rpn =d0rpmc+dd0rpn-dd0mrpn;
  Double_t d0zoinsigma = 0.;
  if(covar[0] > 0.) d0zoinsigma = d0zo/TMath::Sqrt(covar[2]);
  Double_t d0rpoinsigma = 0.;
  if(covar[2] > 0.) d0rpoinsigma = d0rpo/TMath::Sqrt(covar[0]);
  
    if(fMimicData){
       dd0mrpn=sd0mrpn-sd0mrpo;
       d0rpn =d0rpmc+dd0rpn+dd0mrpn;
    }
    
  Double_t dpt1o =pt1o-pt1mc;
  Double_t dpt1n =dpt1o *(spt1o >0. ? (spt1n /spt1o ) : 1.);
  Double_t pt1n  =pt1mc+dpt1n;
  param[0]=d0rpn;
  param[1]=d0zn ;
  param[4]=pt1n ;

   //cov matrix update
   if(fUpdateSTCovMatrix){
    if(sd0rpo>0.)            covar[0]*=(sd0rpn/sd0rpo)*(sd0rpn/sd0rpo);//yy
    if(sd0zo>0. && sd0rpo>0.)covar[1]*=(sd0rpn/sd0rpo)*(sd0zn/sd0zo);//yz
    if(sd0zo>0.)             covar[2]*=(sd0zn/sd0zo)*(sd0zn/sd0zo);//zz
    if(sd0rpo>0.)            covar[3]*=(sd0rpn/sd0rpo);//yl
    if(sd0zo>0.)             covar[4]*=(sd0zn/sd0zo);//zl
    if(sd0rpo>0.)            covar[6]*=(sd0rpn/sd0rpo);//ysenT
    if(sd0zo>0.)             covar[7]*=(sd0zn/sd0zo);//zsenT
    if(sd0rpo>0. && spt1o>0.)covar[10]*=(sd0rpn/sd0rpo)*(spt1n/spt1o);//ypt
    if(sd0zo>0. && spt1o>0.) covar[11]*=(sd0zn/sd0zo)*(spt1n/spt1o);//zpt
    if(spt1o>0.)             covar[12]*=(spt1n/spt1o);//sinPhipt
    if(spt1o>0.)             covar[13]*=(spt1n/spt1o);//tanTpt
    if(spt1o>0.)             covar[14]*=(spt1n/spt1o)*(spt1n/spt1o);//ptpt
  }
  if(fUpdatePulls){
  
      covar[0]*=pullcorr*pullcorr;//yy
      covar[1]*=pullcorr;//yz
      covar[3]*=pullcorr;//yl
      covar[6]*=pullcorr;//ysenT
      covar[10]*=pullcorr;//ypt
      
  }

  Double_t d0zninsigma = 0.;
  if(covar[0] > 0.) d0zninsigma = d0zn/TMath::Sqrt(covar[2]);
  Double_t d0rpninsigma = 0.;
  if(covar[2] > 0.) d0rpninsigma = d0rpn/TMath::Sqrt(covar[0]);

  // Copy the smeared parameters to the AOD track
  Double_t x[3];
  Double_t p[3];
  et.GetXYZ(x);
  et.GetPxPyPz(p);
  Double_t cv[21];
  et.GetCovarianceXYZPxPyPz(cv);

  //if(fIsPbPb2018){
  //  printf("--- is_kFirst_Trk: %d\n", is_kFirst_Trk);
  //  printf("--- is_kOnlySecond_Trk: %d\n", is_kOnlySecond_Trk);
  //  for(UInt_t ibin=0; ibin<21; ibin++) printf("%f\n",cv[ibin]);
  //  printf("\n");
  //}

  if(fIsAOD) {
    AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
    aodtrack->SetPosition(x,kFALSE);
    aodtrack->SetP(p,kTRUE);
    aodtrack->SetCovMatrix(cv);
    // Mark the track as "improved" with a trick (this is done with a trick using layer 7 (ie the 8th))
    UChar_t itsClusterMap = aodtrack->GetITSClusterMap();
    SETBIT(itsClusterMap,7);
    aodtrack->SetITSClusterMap(itsClusterMap);
    //
    
  } else {
    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
    Short_t sign = esdtrack->Charge();
    esdtrack->Set(x,p,cv,sign);
    esdtrack->RelateToVVertex(InputEvent()->GetPrimaryVertex(), bz,100.);
    // Mark the track as "improved" with a trick (this is done with a trick using layer 7 (ie the 8th))
    UChar_t itsClusterMap = esdtrack->GetITSClusterMap();
    SETBIT(itsClusterMap,7);
    esdtrack->SetITSClusterMap(itsClusterMap);
  }


  // write out debug infos
  if (fDebugNtuple->GetEntriesFast()<fNDebug) {
    Int_t idbg=0;
    fDebugVars[idbg++]=pdgcode;
    fDebugVars[idbg++]=ptmc  ;
    fDebugVars[idbg++]=d0rpo ;
    fDebugVars[idbg++]=d0zo  ;
    fDebugVars[idbg++]=pt1o  ;
    fDebugVars[idbg++]=sd0rpo;
    fDebugVars[idbg++]=sd0zo ;
    fDebugVars[idbg++]=spt1o ;
    fDebugVars[idbg++]=d0rpn ;
    fDebugVars[idbg++]=d0zn  ;
    fDebugVars[idbg++]=pt1n  ;
    fDebugVars[idbg++]=sd0rpn;
    fDebugVars[idbg++]=sd0zn ;
    fDebugVars[idbg++]=spt1n ;
    fDebugVars[idbg++]=d0rpmc;
    fDebugVars[idbg++]=d0zmc ;
    fDebugVars[idbg++]=pt1mc ;
    fDebugVars[idbg++]=pullcorr ;
    fDebugVars[idbg++]=d0zoinsigma ;
    fDebugVars[idbg++]=d0zninsigma ;
    fDebugVars[idbg++]=d0rpoinsigma ;
    fDebugVars[idbg++]=d0rpninsigma ;
    fDebugNtuple->Fill(fDebugVars);
    PostData(1,fDebugOutput);
  }
}

AliESDVertex* AliAnalysisTaskSEImproveITS::RecalculateVertex(const AliVVertex *old,TObjArray *tracks,Double_t bField) {
  //
  // Helper function to recalculate a vertex.
  //

  static UShort_t ids[]={1,2,3}; //TODO: unsave...
  AliVertexerTracks vertexer(bField);
  vertexer.SetVtxStart(old->GetX(),old->GetY(),old->GetZ());
  AliESDVertex *vertex=vertexer.VertexForSelectedTracks(tracks,ids);
  return vertex;
}

Double_t AliAnalysisTaskSEImproveITS::EvalGraph(Double_t x,const TGraph *graph,const TGraph *graphSA) const {
  //
  // Evaluates a TGraph without linear extrapolation. Instead the last
  // valid point of the graph is used when out of range.
  // The function assumes an ascending order of X.
  //

  if(!graph){
    printf("\tEvalGraph fails !\n");
    return 0.;
  } 

  // TODO: find a pretty solution for this:
  Int_t    n   =graph->GetN();
  Double_t xmin=graph->GetX()[0  ];
  Double_t xmax=graph->GetX()[n-1];
  if (x<xmin) {
    if(!graphSA) return graph->Eval(xmin);
    Double_t xminSA=graphSA->GetX()[0];
    if(x<xminSA) return graphSA->Eval(xminSA);
    return graphSA->Eval(x);
  }
  if (x>xmax) return graph->Eval(xmax);
  return graph->Eval(x);
}

//________________________________________________________________________

Int_t AliAnalysisTaskSEImproveITS::PhiBin(Double_t phi) const { 
  Double_t pi=TMath::Pi();
  if(phi>2.*pi || phi<0.) return -1;

  if(!fIsPbPb2018){
    if((phi<=(pi/4.)) || (phi>7.*(pi/4.))) return 0;
    if((phi>(pi/4.)) && (phi<=3.*(pi/4.))) return 1;
    if((phi>3.*(pi/4.)) && (phi<=5.*(pi/4.))) return 2;
    if((phi>(5.*pi/4.)) && (phi<=7.*(pi/4.))) return 3;
  }
  else{ // correction performed in 24 φ bins for PbPb 2018 periods
    Double_t width = 2.*pi/24;
    Int_t jBin=TMath::Floor(phi/width);
    return jBin; // by construction is 0<=jBin<23 since phi is in 0-2pi
  }
  
  return -1;
}

ClassImp(AliAnalysisTaskSEImproveITS);

