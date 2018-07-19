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
#include "AliAnalysisVertexingHF.h"
#include "AliNeutralTrackParam.h"
#include "AliAnalysisTaskSEImproveITS.h"

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
   fRunInVertexing(kFALSE),
   fImproveTracks(kTRUE),
   fUpdateSecVertCovMat(kTRUE),
   fUpdateSTCovMatrix(kTRUE),
   fUpdatePulls(kTRUE),
   fMimicData(kFALSE),
   fIsAOD        (kTRUE),
   fMCs         (0),
   fDebugOutput (0),
   fDebugNtuple (0),
   fDebugVars   (0), 
   fNDebug      (0)
{
  //
  // Default constructor.
  for(Int_t jh=0; jh<2; jh++){
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
   fRunInVertexing(isRunInVertexing),
   fImproveTracks(kTRUE),
   fUpdateSecVertCovMat(kTRUE),
   fUpdateSTCovMatrix(kTRUE),
   fUpdatePulls(kTRUE),
   fMimicData(kFALSE),
   fIsAOD        (kTRUE),
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
  }
  
  TString resfileCurURI = Form("alien:///alice/cern.ch/user/p/pwg_hf/common/Improver/%s/%s/ITSgraphs_Current.root",period,systematic);
  TFile *resfileCur=TFile::Open(resfileCurURI.Data());
  if(resfileCur) {    
    if(resfileCur->Get("D0RPResP" )) {
      fD0RPResPCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResP" )));
      fD0RPResPCur ->SetName("D0RPResPCur" );
    }
    if(resfileCur->Get("D0RPResK" )) {
      fD0RPResKCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResK" )));
      fD0RPResKCur ->SetName("D0RPResKCur" );
    }
    if(resfileCur->Get("D0RPResPi")) {
      fD0RPResPiCur=new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResPi")));
      fD0RPResPiCur->SetName("D0RPResPiCur");
    }
    if(resfileCur->Get("D0RPResE")) {
      fD0RPResECur=new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResE")));
      fD0RPResECur->SetName("D0RPResECur");
    }
    if(resfileCur->Get("D0RPSigmaPullRatioP" )) {
      fD0RPSigmaPullRatioP =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPSigmaPullRatioP" )));
    }
    if(resfileCur->Get("D0RPSigmaPullRatioK" )) {
      fD0RPSigmaPullRatioK =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPSigmaPullRatioK" )));
    }
    if(resfileCur->Get("D0RPSigmaPullRatioPi")) {
      fD0RPSigmaPullRatioPi=new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPSigmaPullRatioPi")));
    }
    if(resfileCur->Get("D0RPSigmaPullRatioE")) {
      fD0RPSigmaPullRatioE=new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPSigmaPullRatioE")));
    }
    for(Int_t j=0; j<2; j++){
      for(Int_t i=0; i<4; i++){
	if(resfileCur->Get(Form("D0RPMeanP_B%d_phi%d",j,i))) {
	  fD0RPMeanPCur[j][i]=new TGraph(*static_cast<TGraph*>(resfileCur->Get(Form("D0RPMeanP_B%d_phi%d",j,i))));
	  fD0RPMeanPCur[j][i]->SetName(Form("D0RPMeanPCur_B%d_phi%d",j,i));
	}
	if(resfileCur->Get(Form("D0RPMeanK_B%d_phi%d",j,i))) {
	  fD0RPMeanKCur[j][i]=new TGraph(*static_cast<TGraph*>(resfileCur->Get(Form("D0RPMeanK_B%d_phi%d",j,i))));
	  fD0RPMeanKCur[j][i]->SetName(Form("D0RPMeanKCur_B%d_phi%d",j,i));
	}
	if(resfileCur->Get(Form("D0RPMeanPi_B%d_phi%d",j,i))) {
	  fD0RPMeanPiCur[j][i]=new TGraph(*static_cast<TGraph*>(resfileCur->Get(Form("D0RPMeanPi_B%d_phi%d",j,i))));
	  fD0RPMeanPiCur[j][i]->SetName(Form("D0RPMeanPiCur_B%d_phi%d",j,i));
	}
	if(resfileCur->Get(Form("D0RPMeanE_B%d_phi%d",j,i))) {
	  fD0RPMeanECur[j][i]=new TGraph(*static_cast<TGraph*>(resfileCur->Get(Form("D0RPMeanE_B%d_phi%d",j,i))));
	  fD0RPMeanECur[j][i]->SetName(Form("D0RPMeanECur_B%d_phi%d",j,i));
	}
      }
    }
    if(resfileCur->Get("D0ZResP"  )) {
      fD0ZResPCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResP"  )));
      fD0ZResPCur  ->SetName("D0ZResPCur"  );
    }
    if(resfileCur->Get("D0ZResK"  )) {
      fD0ZResKCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResK"  )));
      fD0ZResKCur  ->SetName("D0ZResKCur"  );
    }
    if(resfileCur->Get("D0ZResPi" )) {
      fD0ZResPiCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResPi" )));
      fD0ZResPiCur ->SetName("D0ZResPiCur" );
    }
    if(resfileCur->Get("D0ZResE" )) {
      fD0ZResECur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResE" )));
      fD0ZResECur ->SetName("D0ZResECur" );
    }
    if(resfileCur->Get("Pt1ResP"  )) {
      fPt1ResPCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResP"  )));
      fPt1ResPCur  ->SetName("Pt1ResPCur"  );
    }
    if(resfileCur->Get("Pt1ResK"  )) {
      fPt1ResKCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResK"  )));
      fPt1ResKCur  ->SetName("Pt1ResKCur"  );
    }
    if(resfileCur->Get("Pt1ResPi" )) {
      fPt1ResPiCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResPi" )));
      fPt1ResPiCur ->SetName("Pt1ResPiCur" );
    }
    if(resfileCur->Get("Pt1ResE" )) {
      fPt1ResECur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResE" )));
      fPt1ResECur ->SetName("Pt1ResECur" );
    }
    if(resfileCur->Get("D0RPResPSA" )) {
      fD0RPResPCurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResPSA" )));
      fD0RPResPCurSA ->SetName("D0RPResPCurSA" );
    }
    if(resfileCur->Get("D0RPResKSA" )) {
      fD0RPResKCurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResKSA" )));
      fD0RPResKCurSA ->SetName("D0RPResKCurSA" );
    }
    if(resfileCur->Get("D0RPResPiSA")) {
      fD0RPResPiCurSA=new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResPiSA")));
      fD0RPResPiCurSA->SetName("D0RPResPiCurSA");
    }
    if(resfileCur->Get("D0RPResESA")) {
      fD0RPResECurSA=new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResESA")));
      fD0RPResECurSA->SetName("D0RPResECurSA");
    }
    if(resfileCur->Get("D0ZResPSA"  )) {
      fD0ZResPCurSA  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResPSA"  )));
      fD0ZResPCurSA  ->SetName("D0ZResPCurSA"  );
    }
    if(resfileCur->Get("D0ZResKSA"  )) {
      fD0ZResKCurSA  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResKSA"  )));
      fD0ZResKCurSA  ->SetName("D0ZResKCurSA"  );
    }
    if(resfileCur->Get("D0ZResPiSA" )) {
      fD0ZResPiCurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResPiSA" )));
      fD0ZResPiCurSA ->SetName("D0ZResPiCurSA" );
    }
    if(resfileCur->Get("D0ZResESA" )) {
      fD0ZResECurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResESA" )));
      fD0ZResECurSA ->SetName("D0ZResECurSA" );
    }
    if(resfileCur->Get("Pt1ResPSA"  )) {
      fPt1ResPCurSA  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResPSA"  )));
      fPt1ResPCurSA  ->SetName("Pt1ResPCurSA"  );
    }
    if(resfileCur->Get("Pt1ResKSA"  )) {
      fPt1ResKCurSA  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResKSA"  )));
      fPt1ResKCurSA  ->SetName("Pt1ResKCurSA"  );
    }
    if(resfileCur->Get("Pt1ResPiSA" )) {
      fPt1ResPiCurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResPiSA" )));
      fPt1ResPiCurSA ->SetName("Pt1ResPiCurSA" );
    }
    if(resfileCur->Get("Pt1ResESA" )) {
      fPt1ResECurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResESA" )));
      fPt1ResECurSA ->SetName("Pt1ResECurSA" );
    }
    delete resfileCur;
  }
    
  TString resfileUpgURI = Form("alien:///alice/cern.ch/user/p/pwg_hf/common/Improver/%s/%s/ITSgraphs_NewAll-X0.3-Res4um.root",period,systematic);
  TFile *resfileUpg=TFile::Open(resfileUpgURI.Data());
  if(resfileUpg) {
    if(resfileUpg->Get("D0RPResP" )) {
      fD0RPResPUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResP" )));
      fD0RPResPUpg ->SetName("D0RPResPUpg" );
    }
    if(resfileUpg->Get("D0RPResK" )) {
      fD0RPResKUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResK" )));
      fD0RPResKUpg ->SetName("D0RPResKUpg" );
    }
    if(resfileUpg->Get("D0RPResPi")) {
      fD0RPResPiUpg=new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResPi")));
      fD0RPResPiUpg->SetName("D0RPResPiUpg");
    }
    if(resfileUpg->Get("D0RPResE")) {
      fD0RPResEUpg=new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResE")));
      fD0RPResEUpg->SetName("D0RPResEUpg");
    }
    for(Int_t j=0; j<2; j++){
      for(Int_t i=0; i<4; i++){
	if(resfileUpg->Get(Form("D0RPMeanP_B%d_phi%d",j,i))) {
	  fD0RPMeanPUpg[j][i]=new TGraph(*static_cast<TGraph*>(resfileUpg->Get(Form("D0RPMeanP_B%d_phi%d",j,i))));
	  fD0RPMeanPUpg[j][i]->SetName(Form("D0RPMeanPUpg_B%d_phi%d",j,i));
	}
	if(resfileUpg->Get(Form("D0RPMeanK_B%d_phi%d",j,i))) {
	  fD0RPMeanKUpg[j][i]=new TGraph(*static_cast<TGraph*>(resfileUpg->Get(Form("D0RPMeanK_B%d_phi%d",j,i))));
	  fD0RPMeanKUpg[j][i]->SetName(Form("D0RPMeanKUpg_B%d_phi%d",j,i));
	}
	if(resfileUpg->Get(Form("D0RPMeanPi_B%d_phi%d",j,i))) {
	  fD0RPMeanPiUpg[j][i]=new TGraph(*static_cast<TGraph*>(resfileUpg->Get(Form("D0RPMeanPi_B%d_phi%d",j,i))));
	  fD0RPMeanPiUpg[j][i]->SetName(Form("D0RPMeanPiUpg_B%d_phi%d",j,i));
	}
	if(resfileUpg->Get(Form("D0RPMeanE_B%d_phi%d",j,i))) {
	  fD0RPMeanEUpg[j][i]=new TGraph(*static_cast<TGraph*>(resfileUpg->Get(Form("D0RPMeanE_B%d_phi%d",j,i))));
	  fD0RPMeanEUpg[j][i]->SetName(Form("D0RPMeanEUpg_B%d_phi%d",j,i));
	}
      }
    }
    if(resfileUpg->Get("D0ZResP"  )) {
      fD0ZResPUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResP"  )));
      fD0ZResPUpg  ->SetName("D0ZResPUpg"  );
    }
    if(resfileUpg->Get("D0ZResK"  )) {
      fD0ZResKUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResK"  )));
      fD0ZResKUpg  ->SetName("D0ZResKUpg"  );
    }
    if(resfileUpg->Get("D0ZResPi" )) {
      fD0ZResPiUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResPi" )));
      fD0ZResPiUpg ->SetName("D0ZResPiUpg" );
    }
    if(resfileUpg->Get("D0ZResE" )) {
      fD0ZResEUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResE" )));
      fD0ZResEUpg ->SetName("D0ZResEUpg" );
    }
    if(resfileUpg->Get("Pt1ResP"  )) {
      fPt1ResPUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResP"  )));
      fPt1ResPUpg  ->SetName("Pt1ResPUpg"  );
    }
    if(resfileUpg->Get("Pt1ResK"  )) {
      fPt1ResKUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResK"  )));
      fPt1ResKUpg  ->SetName("Pt1ResKUpg"  );
    }
    if(resfileUpg->Get("Pt1ResPi" )) {
      fPt1ResPiUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResPi" )));
      fPt1ResPiUpg ->SetName("Pt1ResPiUpg" );
    }
    if(resfileUpg->Get("Pt1ResE" )) {
      fPt1ResEUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResE" )));
      fPt1ResEUpg ->SetName("Pt1ResEUpg" );
    }
    if(resfileUpg->Get("D0RPResPSA" )) {
      fD0RPResPUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResPSA" )));
      fD0RPResPUpgSA ->SetName("D0RPResPUpgSA" );
    }
    if(resfileUpg->Get("D0RPResKSA" )) {
      fD0RPResKUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResKSA" )));
      fD0RPResKUpgSA ->SetName("D0RPResKUpgSA" );
    }
    if(resfileUpg->Get("D0RPResPiSA")) {
      fD0RPResPiUpgSA=new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResPiSA")));
      fD0RPResPiUpgSA->SetName("D0RPResPiUpgSA");
    }
    if(resfileUpg->Get("D0RPResESA")) {
      fD0RPResEUpgSA=new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResESA")));
      fD0RPResEUpgSA->SetName("D0RPResEUpgSA");
    }
    if(resfileUpg->Get("D0ZResPSA"  )) {
      fD0ZResPUpgSA  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResPSA"  )));
      fD0ZResPUpgSA  ->SetName("D0ZResPUpgSA"  );
    }
    if(resfileUpg->Get("D0ZResKSA"  )) {
      fD0ZResKUpgSA  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResKSA"  )));
      fD0ZResKUpgSA  ->SetName("D0ZResKUpgSA"  );
    }
    if(resfileUpg->Get("D0ZResPiSA" )) {
      fD0ZResPiUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResPiSA" )));
      fD0ZResPiUpgSA ->SetName("D0ZResPiUpgSA" );
    }
    if(resfileUpg->Get("D0ZResESA" )) {
      fD0ZResEUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResESA" )));
      fD0ZResEUpgSA ->SetName("D0ZResEUpgSA" );
    }
    if(resfileUpg->Get("Pt1ResPSA"  )) {
      fPt1ResPUpgSA  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResPSA"  )));
      fPt1ResPUpgSA  ->SetName("Pt1ResPUpgSA"  );
    }
    if(resfileUpg->Get("Pt1ResKSA"  )) {
      fPt1ResKUpgSA  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResKSA"  )));
      fPt1ResKUpgSA  ->SetName("Pt1ResKUpgSA"  );
    }
    if(resfileUpg->Get("Pt1ResPiSA" )) {
      fPt1ResPiUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResPiSA" )));
      fPt1ResPiUpgSA ->SetName("Pt1ResPiUpgSA" );
    }
    if(resfileUpg->Get("Pt1ResESA" )) {
      fPt1ResEUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResESA" )));
      fPt1ResEUpgSA ->SetName("Pt1ResEUpgSA" );
    }
    delete resfileUpg;
  }
    
  DefineOutput(1,TNtuple::Class());
}

AliAnalysisTaskSEImproveITS::~AliAnalysisTaskSEImproveITS() {
  //
  // Destructor.
  //
  delete fDebugOutput;
}

void AliAnalysisTaskSEImproveITS::UserCreateOutputObjects() {
  //
  // Creation of user output objects.
  //
  fDebugOutput=new TList();
  fDebugOutput->SetOwner();
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
    
    // first loop on candidates to fill them in case of reduced AODs
    // this is done to have the same behaviour of the improver with full (pp, p-Pb) and recuced (Pb-Pb) candidates
    AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
    
    // D0->Kpi
    TClonesArray *array2Prong=static_cast<TClonesArray*>(ev->GetList()->FindObject("D0toKpi"));
    if (array2Prong) {
      for (Int_t icand=0;icand<array2Prong->GetEntries();++icand) {
	AliAODRecoDecayHF2Prong *decay=static_cast<AliAODRecoDecayHF2Prong*>(array2Prong->At(icand));
	vHF->FillRecoCand(ev,(AliAODRecoDecayHF2Prong*)decay);
      }
    }
    // Dstar->Kpipi
    TClonesArray *arrayCascade=static_cast<TClonesArray*>(ev->GetList()->FindObject("Dstar"));
    if (arrayCascade) {
      for (Int_t icand=0;icand<arrayCascade->GetEntries();++icand) {
	AliAODRecoCascadeHF *decayDstar=static_cast<AliAODRecoCascadeHF*>(arrayCascade->At(icand));
	vHF->FillRecoCasc(ev,((AliAODRecoCascadeHF*)decayDstar),kTRUE);
      }
    }
    // Three prong
    TClonesArray *array3Prong=static_cast<TClonesArray*>(ev->GetList()->FindObject("Charm3Prong"));
    if (array3Prong) {
      for (Int_t icand=0;icand<array3Prong->GetEntries();++icand) {
	AliAODRecoDecayHF3Prong *decay=static_cast<AliAODRecoDecayHF3Prong*>(array3Prong->At(icand));
	vHF->FillRecoCand(ev,(AliAODRecoDecayHF3Prong*)decay);
      }
    }
    
    
    // Smear all tracks
    fMCs=static_cast<TClonesArray*>(ev->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
    if (!fMCs) return;
    if (fImproveTracks) {
      for(Int_t itrack=0;itrack<ev->GetNumberOfTracks();++itrack) {
	AliAODTrack * trk = dynamic_cast<AliAODTrack*>(ev->GetTrack(itrack));
	if(!trk) AliFatal("Not a standard AOD");
	SmearTrack(trk,bz);
      }
    }
    
    // TODO: recalculated primary vertex
    AliVVertex *primaryVertex=ev->GetPrimaryVertex();
    
    
    // Recalculate all candidates
    // D0->Kpi
    if (array2Prong) {
      for (Int_t icand=0;icand<array2Prong->GetEntries();++icand) {
	AliAODRecoDecayHF2Prong *decay=static_cast<AliAODRecoDecayHF2Prong*>(array2Prong->At(icand));
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
      for (Int_t icand=0;icand<arrayCascade->GetEntries();++icand) {
	AliAODRecoCascadeHF *decayDstar=static_cast<AliAODRecoCascadeHF*>(arrayCascade->At(icand));
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
      for (Int_t icand=0;icand<array3Prong->GetEntries();++icand) {
	AliAODRecoDecayHF3Prong *decay=static_cast<AliAODRecoDecayHF3Prong*>(array3Prong->At(icand));
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
	AliESDtrack * trk = dynamic_cast<AliESDtrack*>(evesd->GetTrack(itrack));
	if(!trk) AliFatal("No a standard ESD");
	SmearTrack(trk,bz);
      }
    }    
  }// end ESD
  
}

void AliAnalysisTaskSEImproveITS::SmearTrack(AliVTrack *track,Double_t bz) {
  // Early exit, if this track has nothing in common with the ITS

  if (!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)))
    return;

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
    mc = dynamic_cast<AliVParticle*>(fMCs->At(imc));
  }
  else {
    mc = dynamic_cast<AliVParticle*>(fMCEvent->GetTrack(imc));
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
  et.PropagateToDCA(&vtx,track->GetBz(),10.);
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
  if (fDebugNtuple->GetEntries()<fNDebug) {
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

  if(!graph) return 0.;

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
  if((phi<=(pi/4.)) || (phi>7.*(pi/4.))) return 0;
  if((phi>(pi/4.)) && (phi<=3.*(pi/4.))) return 1;
  if((phi>3.*(pi/4.)) && (phi<=5.*(pi/4.))) return 2;
  if((phi>(5.*pi/4.)) && (phi<=7.*(pi/4.))) return 3;
  return -1;
}

ClassImp(AliAnalysisTaskSEImproveITS);

