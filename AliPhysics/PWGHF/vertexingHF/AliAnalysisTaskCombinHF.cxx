/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: $ */

//*************************************************************************
// Class AliAnalysisTaskCombinHF
// AliAnalysisTaskSE to build D meson candidates by combining tracks
//  background is computed LS and track rotations is
// Authors: F. Prino, A. Rossi
/////////////////////////////////////////////////////////////

#include <TList.h>
#include <TH1F.h>
#include <TDatabasePDG.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAnalysisTaskCombinHF.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCombinHF);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskCombinHF::AliAnalysisTaskCombinHF():
  AliAnalysisTaskSE(),
  fOutput(0x0),
  fListCuts(0x0),
  fHistNEvents(0x0),
  fHistEventMultCent(0x0),
  fHistEventMultCentEvSel(0x0),
  fHistEventMultZv(0x0),
  fHistEventMultZvEvSel(0x0),
  fHistXsecVsPtHard(0x0),
  fHistTrackStatus(0x0),
  fHistTrackEtaMultZv(0x0),
  fHistCheckOrigin(0x0),
  fHistCheckOriginRecoD(0x0),
  fHistCheckOriginRecoVsGen(0x0),
  fHistCheckDecChan(0x0),
  fHistCheckDecChanAcc(0x0),
  fPtVsYVsMultGenPrompt(0x0),
  fPtVsYVsMultGenLargeAccPrompt(0x0),
  fPtVsYVsMultGenLimAccPrompt(0x0),
  fPtVsYVsMultGenAccPrompt(0x0),
  fPtVsYVsMultGenAccEvSelPrompt(0x0),
  fPtVsYVsMultRecoPrompt(0x0),
  fPtVsYVsMultGenFeeddw(0x0),
  fPtVsYVsMultGenLargeAccFeeddw(0x0),
  fPtVsYVsMultGenLimAccFeeddw(0x0),
  fPtVsYVsMultGenAccFeeddw(0x0),
  fPtVsYVsMultGenAccEvSelFeeddw(0x0),
  fPtVsYVsMultRecoFeeddw(0x0),
  fPtVsYVsPtBGenFeeddw(0x0),
  fPtVsYVsPtBGenLargeAccFeeddw(0x0),
  fPtVsYVsPtBGenLimAccFeeddw(0x0),
  fPtVsYVsPtBGenAccFeeddw(0x0),
  fPtVsYVsPtBGenAccEvSelFeeddw(0x0),
  fPtVsYVsPtBRecoFeeddw(0x0),
  fMassVsPtVsY(0x0),
  fMassVsPtVsYRot(0x0),
  fMassVsPtVsYLSpp(0x0),
  fMassVsPtVsYLSmm(0x0),
  fMassVsPtVsYSig(0x0),
  fMassVsPtVsYRefl(0x0),
  fMassVsPtVsYBkg(0x0),
  fBMohterPtGen(0x0),
  fNSelected(0x0),
  fNormRotated(0x0),
  fDeltaMass(0x0),
  fDeltaMassFullAnalysis(0x0),
  fMassVsPtVsYME(0x0),
  fMassVsPtVsYMELSpp(0x0),
  fMassVsPtVsYMELSmm(0x0),
  fEventsPerPool(0x0),
  fMixingsPerPool(0x0),
  fMassVsPtVsCosthSt(0x0),
  fMassVsPtVsCosthStRot(0x0),
  fMassVsPtVsCosthStLSpp(0x0),
  fMassVsPtVsCosthStLSmm(0x0),
  fMassVsPtVsCosthStSig(0x0),
  fMassVsPtVsCosthStRefl(0x0),
  fMassVsPtVsCosthStBkg(0x0),
  fMassVsPtVsCosthStME(0x0),
  fMassVsPtVsCosthStMELSpp(0x0),
  fMassVsPtVsCosthStMELSmm(0x0),
  fHistonSigmaTPCPion(0x0),
  fHistonSigmaTPCPionGoodTOF(0x0),
  fHistonSigmaTOFPion(0x0),
  fHistonSigmaTPCKaon(0x0),
  fHistonSigmaTPCKaonGoodTOF(0x0),
  fHistonSigmaTOFKaon(0x0),
  fHistoPtKPtPiPtD(0x0),
  fHistoPtKPtPiPtDSig(0x0),
  fFilterMask(BIT(4)),
  fTrackCutsAll(0x0),
  fTrackCutsPion(0x0),
  fTrackCutsKaon(0x0),
  fFillHistosVsCosThetaStar(kFALSE),
  fApplyCutCosThetaStar(kFALSE),
  fCutCosThetaStar(999.),
  fPhiMassCut(99999.),
  fCutCos3PiKPhiRFrame(-1.1),
  fCutCosPiDsLabFrame(1.1),
  fPidHF(0x0),
  fAnalysisCuts(0x0),
  fMinMass(1.720),
  fMaxMass(2.150),
  fMaxPt(10.),
  fPtBinWidth(0.5),
  fEtaAccCut(0.9),
  fPtAccCut(0.1),
  fNRotations(9),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fNRotations3(9),
  fMinAngleForRot3(2*TMath::Pi()/6),
  fMaxAngleForRot3(4*TMath::Pi()/6),
  fCounter(0x0),
  fMeson(kDzero),
  fReadMC(kFALSE),
  fEnforceMBTrigMaskInMC(kTRUE),
  fGoUpToQuark(kTRUE),
  fFullAnalysis(0),
  fSignalOnlyMC(kFALSE),
  fSelectPtHardRange(kFALSE),
  fMinPtHard(0.),
  fMaxPtHard(999999.),
  fPIDstrategy(knSigma),
  fmaxPforIDPion(0.8),
  fmaxPforIDKaon(2.),
  fKeepNegID(kFALSE),
  fPIDselCaseZero(0),
  fBayesThresKaon(0.4),
  fBayesThresPion(0.4),
  fDoEventMixing(1),
  fNumberOfEventsForMixing(20),
  fMaxzVertDistForMix(5.),
  fMaxMultDiffForMix(5.),
  fNzVertPools(1),
  fNzVertPoolsLimSize(2),
  fzVertPoolLims(0x0),
  fNMultPools(1),
  fNMultPoolsLimSize(2),
  fMultPoolLims(0x0),
  fNOfPools(1),
  fEventBuffer(0x0),
  fEventInfo(new TObjString("")),
  fVtxZ(0),
  fMultiplicity(0),
  fNumOfMultBins(200),
  fMinMultiplicity(-0.5),
  fMaxMultiplicity(199.5),
  fKaonTracks(0x0),
  fPionTracks(0x0)
{
  /// default constructor
}

//________________________________________________________________________
AliAnalysisTaskCombinHF::AliAnalysisTaskCombinHF(Int_t meson, AliRDHFCuts* analysiscuts):
  AliAnalysisTaskSE("DmesonCombin"),
  fOutput(0x0),
  fListCuts(0x0),
  fHistNEvents(0x0),
  fHistEventMultCent(0x0),
  fHistEventMultCentEvSel(0x0),
  fHistEventMultZv(0x0),
  fHistEventMultZvEvSel(0x0),
  fHistXsecVsPtHard(0x0),
  fHistTrackStatus(0x0),
  fHistTrackEtaMultZv(0x0),
  fHistCheckOrigin(0x0),
  fHistCheckOriginRecoD(0x0),
  fHistCheckOriginRecoVsGen(0x0),
  fHistCheckDecChan(0x0),
  fHistCheckDecChanAcc(0x0),
  fPtVsYVsMultGenPrompt(0x0),
  fPtVsYVsMultGenLargeAccPrompt(0x0),
  fPtVsYVsMultGenLimAccPrompt(0x0),
  fPtVsYVsMultGenAccPrompt(0x0),
  fPtVsYVsMultGenAccEvSelPrompt(0x0),
  fPtVsYVsMultRecoPrompt(0x0),
  fPtVsYVsMultGenFeeddw(0x0),
  fPtVsYVsMultGenLargeAccFeeddw(0x0),
  fPtVsYVsMultGenLimAccFeeddw(0x0),
  fPtVsYVsMultGenAccFeeddw(0x0),
  fPtVsYVsMultGenAccEvSelFeeddw(0x0),
  fPtVsYVsMultRecoFeeddw(0x0),
  fPtVsYVsPtBGenFeeddw(0x0),
  fPtVsYVsPtBGenLargeAccFeeddw(0x0),
  fPtVsYVsPtBGenLimAccFeeddw(0x0),
  fPtVsYVsPtBGenAccFeeddw(0x0),
  fPtVsYVsPtBGenAccEvSelFeeddw(0x0),
  fPtVsYVsPtBRecoFeeddw(0x0),
  fMassVsPtVsY(0x0),
  fMassVsPtVsYRot(0x0),
  fMassVsPtVsYLSpp(0x0),
  fMassVsPtVsYLSmm(0x0),
  fMassVsPtVsYSig(0x0),
  fMassVsPtVsYRefl(0x0),
  fMassVsPtVsYBkg(0x0),
  fBMohterPtGen(0x0),
  fNSelected(0x0),
  fNormRotated(0x0),
  fDeltaMass(0x0),
  fDeltaMassFullAnalysis(0x0),
  fMassVsPtVsYME(0x0),
  fMassVsPtVsYMELSpp(0x0),
  fMassVsPtVsYMELSmm(0x0),
  fEventsPerPool(0x0),
  fMixingsPerPool(0x0),
  fMassVsPtVsCosthSt(0x0),
  fMassVsPtVsCosthStRot(0x0),
  fMassVsPtVsCosthStLSpp(0x0),
  fMassVsPtVsCosthStLSmm(0x0),
  fMassVsPtVsCosthStSig(0x0),
  fMassVsPtVsCosthStRefl(0x0),
  fMassVsPtVsCosthStBkg(0x0),
  fMassVsPtVsCosthStME(0x0),
  fMassVsPtVsCosthStMELSpp(0x0),
  fMassVsPtVsCosthStMELSmm(0x0),
  fHistonSigmaTPCPion(0x0),
  fHistonSigmaTPCPionGoodTOF(0x0),
  fHistonSigmaTOFPion(0x0),
  fHistonSigmaTPCKaon(0x0),
  fHistonSigmaTPCKaonGoodTOF(0x0),
  fHistonSigmaTOFKaon(0x0),
  fHistoPtKPtPiPtD(0x0),
  fHistoPtKPtPiPtDSig(0x0),
  fFilterMask(BIT(4)),
  fTrackCutsAll(0x0),
  fTrackCutsPion(0x0),
  fTrackCutsKaon(0x0),
  fFillHistosVsCosThetaStar(kFALSE),
  fApplyCutCosThetaStar(kFALSE),
  fCutCosThetaStar(999.),
  fPhiMassCut(99999.),
  fCutCos3PiKPhiRFrame(-1),
  fCutCosPiDsLabFrame(1.1),
  fPidHF(0x0),
  fAnalysisCuts(analysiscuts),
  fMinMass(1.720),
  fMaxMass(2.150),
  fMaxPt(10.),
  fPtBinWidth(0.5),
  fEtaAccCut(0.9),
  fPtAccCut(0.1),
  fNRotations(9),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fNRotations3(9),
  fMinAngleForRot3(2*TMath::Pi()/6),
  fMaxAngleForRot3(4*TMath::Pi()/6),
  fCounter(0x0),
  fMeson(meson),
  fReadMC(kFALSE),
  fEnforceMBTrigMaskInMC(kTRUE),
  fGoUpToQuark(kTRUE),
  fFullAnalysis(0),
  fSignalOnlyMC(kFALSE),
  fSelectPtHardRange(kFALSE),
  fMinPtHard(0.),
  fMaxPtHard(999999.),
  fPIDstrategy(knSigma),
  fmaxPforIDPion(0.8),
  fmaxPforIDKaon(2.),
  fKeepNegID(kFALSE),
  fPIDselCaseZero(0),
  fBayesThresKaon(0.4),
  fBayesThresPion(0.4),
  fDoEventMixing(1),
  fNumberOfEventsForMixing(20),
  fMaxzVertDistForMix(5.),
  fMaxMultDiffForMix(5.),
  fNzVertPools(1),
  fNzVertPoolsLimSize(2),
  fzVertPoolLims(0x0),
  fNMultPools(1),
  fNMultPoolsLimSize(2),
  fMultPoolLims(0x0),
  fNOfPools(1),
  fEventBuffer(0x0),
  fEventInfo(new TObjString("")),
  fVtxZ(0),
  fMultiplicity(0),
  fNumOfMultBins(200),
  fMinMultiplicity(-0.5),
  fMaxMultiplicity(199.5),
  fKaonTracks(0x0),
  fPionTracks(0x0)
{
  /// standard constructor
  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,AliNormalizationCounter::Class());
  DefineOutput(3,TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskCombinHF::~AliAnalysisTaskCombinHF()
{
  //
  /// Destructor
  //
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistEventMultCent;
    delete fHistEventMultCentEvSel;
    delete fHistEventMultZv;
    delete fHistEventMultZvEvSel;
    delete fHistXsecVsPtHard;
    delete fHistTrackStatus;
    delete fHistTrackEtaMultZv;
    delete fHistCheckOrigin;
    delete fHistCheckOriginRecoD;
    delete fHistCheckOriginRecoVsGen;
    delete fHistCheckDecChan;
    delete fHistCheckDecChanAcc;
    delete fPtVsYVsMultGenPrompt;
    delete fPtVsYVsMultGenLargeAccPrompt;
    delete fPtVsYVsMultGenLimAccPrompt;
    delete fPtVsYVsMultGenAccPrompt;
    delete fPtVsYVsMultGenAccEvSelPrompt;
    delete fPtVsYVsMultRecoPrompt;
    delete fPtVsYVsMultGenFeeddw;
    delete fPtVsYVsMultGenLargeAccFeeddw;
    delete fPtVsYVsMultGenLimAccFeeddw;
    delete fPtVsYVsMultGenAccFeeddw;
    delete fPtVsYVsMultGenAccEvSelFeeddw;
    delete fPtVsYVsMultRecoFeeddw;
    delete fPtVsYVsPtBGenFeeddw;
    delete fPtVsYVsPtBGenLargeAccFeeddw;
    delete fPtVsYVsPtBGenLimAccFeeddw;
    delete fPtVsYVsPtBGenAccFeeddw;
    delete fPtVsYVsPtBGenAccEvSelFeeddw;
    delete fPtVsYVsPtBRecoFeeddw;
    delete fMassVsPtVsY;
    delete fMassVsPtVsYLSpp;
    delete fMassVsPtVsYLSmm;
    delete fMassVsPtVsYRot;
    delete fMassVsPtVsYSig;
    delete fMassVsPtVsYRefl;
    delete fMassVsPtVsYBkg;
    delete fBMohterPtGen;
    delete fNSelected;
    delete fNormRotated;
    delete fDeltaMass;
    delete fDeltaMassFullAnalysis;
    delete fMassVsPtVsYME;
    delete fMassVsPtVsYMELSpp;
    delete fMassVsPtVsYMELSmm;
    delete fMassVsPtVsCosthSt;
    delete fMassVsPtVsCosthStRot;
    delete fMassVsPtVsCosthStLSpp;
    delete fMassVsPtVsCosthStLSmm;
    delete fMassVsPtVsCosthStSig;
    delete fMassVsPtVsCosthStRefl;
    delete fMassVsPtVsCosthStBkg;
    delete fMassVsPtVsCosthStME;
    delete fMassVsPtVsCosthStMELSpp;
    delete fMassVsPtVsCosthStMELSmm;
    delete fHistonSigmaTPCPion;
    delete fHistonSigmaTPCPionGoodTOF;
    delete fHistonSigmaTOFPion;
    delete fHistonSigmaTPCKaon;
    delete fHistonSigmaTPCKaonGoodTOF;
    delete fHistonSigmaTOFKaon;
    delete fHistoPtKPtPiPtD;
    delete fHistoPtKPtPiPtDSig;
  }

  delete fOutput;
  if (fListCuts) delete fListCuts;
  delete fCounter;
  delete fTrackCutsAll;
  delete fTrackCutsPion;
  delete fTrackCutsKaon;
  delete fAnalysisCuts;
  if(fKaonTracks) fKaonTracks->Delete();
  if(fPionTracks) fPionTracks->Delete();
  delete fKaonTracks;
  delete fPionTracks;

  if(fEventBuffer){
    for(Int_t i=0; i<fNOfPools; i++) delete fEventBuffer[i];
    delete fEventBuffer;
  }
  delete fEventInfo;
  delete [] fzVertPoolLims;
  delete [] fMultPoolLims;
}

//________________________________________________________________________
void AliAnalysisTaskCombinHF::ConfigureZVertPools(Int_t nPools, Double_t*  zVertLimits)
{
  /// sets the pools for event mizing in zvertex
  if(fzVertPoolLims) delete [] fzVertPoolLims;
  fNzVertPools=nPools;
  fNzVertPoolsLimSize=nPools+1;
  fzVertPoolLims = new Double_t[fNzVertPoolsLimSize];
  for(Int_t ib=0; ib<fNzVertPoolsLimSize; ib++) fzVertPoolLims[ib]=zVertLimits[ib];
  return;  
}
//________________________________________________________________________
void AliAnalysisTaskCombinHF::ConfigureMultiplicityPools(Int_t nPools, Double_t*  multLimits)
{
  // sets the pools for event mizing in zvertex
  if(fMultPoolLims) delete [] fMultPoolLims;
  fNMultPools=nPools;
  fNMultPoolsLimSize=nPools+1;
  fMultPoolLims = new Double_t[fNMultPoolsLimSize];
  for(Int_t ib=0; ib<nPools+1; ib++) fMultPoolLims[ib]=multLimits[ib];
  return;  
}
//________________________________________________________________________
void AliAnalysisTaskCombinHF::UserCreateOutputObjects()
{
  /// Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskCombinHF::UserCreateOutputObjects() \n");
  
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");
  
  fHistNEvents = new TH1F("hNEvents", "number of events ",10,-0.5,9.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"n. passing IsEvSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"n. rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"n. rejected due to phys sel");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"n. rejected due to not reco vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"n. rejected for contr vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"n. rejected for zSPD-zTrack");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"n. rejected for vertex out of accept");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"n. rejected for pileup");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"no. of out centrality events");
  
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fHistEventMultCent = new TH2F("hEventMultCent"," ; Centrality (V0M) ; N_{tracklets} (|#eta|<1)",100,0.,100.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
  fOutput->Add(fHistEventMultCent);

  fHistEventMultCentEvSel = new TH2F("hEventMultCentEvSel"," ; Centrality (V0M) ; N_{tracklets} (|#eta|<1)",100,0.,100.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
  fOutput->Add(fHistEventMultCentEvSel);

  fHistEventMultZv = new TH2F("hEventMultZv"," ; z_{vertex} (cm) ; N_{tracklets} (|#eta|<1)",30,-15.,15.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
  fOutput->Add(fHistEventMultZv);

  fHistEventMultZvEvSel = new TH2F("hEventMultZvEvSel"," ; z_{vertex} (cm) ; N_{tracklets} (|#eta|<1)",30,-15.,15.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
  fOutput->Add(fHistEventMultZvEvSel);

  fHistXsecVsPtHard = new TH1F("hXsecVsPtHard", " ; pthard (GeV/c) ; Xsec", 200,0.,100.);
  fOutput->Add(fHistXsecVsPtHard);
  
  fHistTrackStatus  = new TH1F("hTrackStatus", "",8,-0.5,7.5);
  fHistTrackStatus->GetXaxis()->SetBinLabel(1,"Not OK");
  fHistTrackStatus->GetXaxis()->SetBinLabel(2,"Track OK");
  fHistTrackStatus->GetXaxis()->SetBinLabel(3,"Kaon, Not OK");
  fHistTrackStatus->GetXaxis()->SetBinLabel(4,"Kaon OK");
  fHistTrackStatus->GetXaxis()->SetBinLabel(5,"Pion, Not OK");
  fHistTrackStatus->GetXaxis()->SetBinLabel(6,"Pion OK");
  fHistTrackStatus->GetXaxis()->SetBinLabel(7,"Kaon||Pion, Not OK");
  fHistTrackStatus->GetXaxis()->SetBinLabel(8,"Kaon||Pion OK");
  fHistTrackStatus->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistTrackStatus->SetMinimum(0);
  fOutput->Add(fHistTrackStatus);
  
  fHistTrackEtaMultZv = new TH3F("hTrackEtaMultZv","",40,-1.,1.,30,-15.,15.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
  fOutput->Add(fHistTrackEtaMultZv);

  Int_t nPtBins = (Int_t)(fMaxPt/fPtBinWidth+0.001);
  Double_t maxPt=fPtBinWidth*nPtBins;

  if(fReadMC){
    
    fHistCheckOrigin=new TH2F("hCheckOrigin"," ; origin ; generator",7,-1.5,5.5,2,-0.5,1.5);
    fHistCheckOrigin->GetYaxis()->SetBinLabel(1,"Hijing");
    fHistCheckOrigin->GetYaxis()->SetBinLabel(2,"Injected");
    fHistCheckOriginRecoD=new TH2F("hCheckOriginRecoD"," ; origin ; generator",7,-1.5,5.5,2,-0.5,1.5);
    fHistCheckOriginRecoD->GetYaxis()->SetBinLabel(1,"Hijing");
    fHistCheckOriginRecoD->GetYaxis()->SetBinLabel(2,"Injected");
    fHistCheckOriginRecoVsGen=new TH2F("hCheckOriginRecoVsGen"," ; Origin (reco D) ; Origin (gen D)",7,-1.5,5.5,7,-1.5,5.5);
    fOutput->Add(fHistCheckOrigin);
    fOutput->Add(fHistCheckOriginRecoD);
    fOutput->Add(fHistCheckOriginRecoVsGen);
  
    fHistCheckDecChan=new TH1F("hCheckDecChan","",7,-2.5,4.5);
    fHistCheckDecChan->SetMinimum(0);
    fOutput->Add(fHistCheckDecChan);
    
    fHistCheckDecChanAcc=new TH1F("hCheckDecChanAcc","",7,-2.5,4.5);
    fHistCheckDecChanAcc->SetMinimum(0);
    fOutput->Add(fHistCheckDecChanAcc);
    
    fPtVsYVsMultGenPrompt = new TH3F("hPtVsYVsMultGenPrompt","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenPrompt);
    
    fPtVsYVsMultGenLargeAccPrompt = new TH3F("hPtVsYVsMultGenLargeAccPrompt","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenLargeAccPrompt);
    
    fPtVsYVsMultGenLimAccPrompt = new TH3F("hPtVsYVsMultGenLimAccPrompt","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenLimAccPrompt);
    
    fPtVsYVsMultGenAccPrompt = new TH3F("hPtVsYVsMultGenAccPrompt","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenAccPrompt);
    
    fPtVsYVsMultGenAccEvSelPrompt = new TH3F("hPtVsYVsMultGenAccEvSelPrompt","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenAccEvSelPrompt);
 
    fPtVsYVsMultRecoPrompt = new TH3F("hPtVsYVsMultRecoPrompt","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultRecoPrompt);

    fPtVsYVsMultGenFeeddw = new TH3F("hPtVsYVsMultGenFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenFeeddw);
    
    fPtVsYVsMultGenLargeAccFeeddw = new TH3F("hPtVsYVsMultGenLargeAccFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenLargeAccFeeddw);
    
    fPtVsYVsMultGenLimAccFeeddw = new TH3F("hPtVsYVsMultGenLimAccFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenLimAccFeeddw);
    
    fPtVsYVsMultGenAccFeeddw = new TH3F("hPtVsYVsMultGenAccFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenAccFeeddw);
    
    fPtVsYVsMultGenAccEvSelFeeddw = new TH3F("hPtVsYVsMultGenAccEvSelFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultGenAccEvSelFeeddw);
 
    fPtVsYVsMultRecoFeeddw = new TH3F("hPtVsYVsMultRecoFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,fNumOfMultBins,fMinMultiplicity,fMaxMultiplicity);
    fOutput->Add(fPtVsYVsMultRecoFeeddw);
 
    fPtVsYVsPtBGenFeeddw = new TH3F("hPtVsYVsPtBGenFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,100,0.,50.);
    fOutput->Add(fPtVsYVsPtBGenFeeddw);
    
    fPtVsYVsPtBGenLargeAccFeeddw = new TH3F("hPtVsYVsPtBGenLargeAccFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,100,0.,50.);
    fOutput->Add(fPtVsYVsPtBGenLargeAccFeeddw);
    
    fPtVsYVsPtBGenLimAccFeeddw = new TH3F("hPtVsYVsPtBGenLimAccFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,100,0.,50.);
    fOutput->Add(fPtVsYVsPtBGenLimAccFeeddw);
    
    fPtVsYVsPtBGenAccFeeddw = new TH3F("hPtVsYVsPtBGenAccFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,100,0.,50.);
    fOutput->Add(fPtVsYVsPtBGenAccFeeddw);
    
    fPtVsYVsPtBGenAccEvSelFeeddw = new TH3F("hPtVsYVsPtBGenAccEvSelFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,100,0.,50.);
    fOutput->Add(fPtVsYVsPtBGenAccEvSelFeeddw);
 
    fPtVsYVsPtBRecoFeeddw = new TH3F("hPtVsYVsPtBRecoFeeddw","",nPtBins,0.,maxPt,20,-1.,1.,100,0.,50.);
    fOutput->Add(fPtVsYVsPtBRecoFeeddw);
 }
  
  
  Int_t nMassBins=static_cast<Int_t>(fMaxMass*1000.-fMinMass*1000.);
  Double_t maxm=fMinMass+nMassBins*0.001;
  fMassVsPtVsY=new TH3F("hMassVsPtVsY","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsY);
  
  fMassVsPtVsYRot=new TH3F("hMassVsPtVsYRot","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsYRot);
  
  fMassVsPtVsYLSpp=new TH3F("hMassVsPtVsYLSpp","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsYLSpp);
  fMassVsPtVsYLSmm=new TH3F("hMassVsPtVsYLSmm","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsYLSmm);
  
  fMassVsPtVsYSig=new TH3F("hMassVsPtVsYSig","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsYSig);
  
  fMassVsPtVsYRefl=new TH3F("hMassVsPtVsYRefl","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsYRefl);
  
  fMassVsPtVsYBkg=new TH3F("hMassVsPtVsYBkg","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsYBkg);

  fBMohterPtGen=new TH1F("hBMohterPtGen","",100,0.,50.);
  fOutput->Add(fBMohterPtGen);
  
  fNSelected=new TH1F("hNSelected","",100,-0.5,99.5);
  fOutput->Add(fNSelected);
  
  fNormRotated=new TH1F("hNormRotated","",11,-0.5,10.5);
  fOutput->Add(fNormRotated);
  
  fDeltaMass=new TH1F("hDeltaMass","",100,-0.4,0.4);
  fOutput->Add(fDeltaMass);
  
  Int_t binSparseDMassRot[5]={nMassBins,100,24,40,20};
  Double_t edgeLowSparseDMassRot[5]={fMinMass,-0.4,0.,-4.,0};
  Double_t edgeHighSparseDMassRot[5]={maxm,0.4,12.,4.,3.14};
  fDeltaMassFullAnalysis=new THnSparseF("fDeltaMassFullAnalysis","fDeltaMassFullAnalysis;inv mass (GeV/c);#Delta inv mass (GeV/c) ; p_{T}^{D} (GeV/c); #Delta p_{T} (GeV/c); daughter angle (2prongs) (rad);",5,binSparseDMassRot,edgeLowSparseDMassRot,edgeHighSparseDMassRot);
  fOutput->Add(fDeltaMassFullAnalysis);
  
  fMassVsPtVsYME=new TH3F("hMassVsPtVsYME","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsYME);

  fMassVsPtVsYMELSpp=new TH3F("hMassVsPtVsYMELSpp","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsYMELSpp);

  fMassVsPtVsYMELSmm=new TH3F("hMassVsPtVsYMELSmm","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,20,-1.,1.);
  fOutput->Add(fMassVsPtVsYMELSmm);

  fNOfPools=fNzVertPools*fNMultPools;
  if(!fzVertPoolLims || !fMultPoolLims) fNOfPools=1;
  if(fDoEventMixing==2) fNOfPools=1;
  if(fNOfPools>1 && fzVertPoolLims && fMultPoolLims){
    fEventsPerPool=new TH2F("hEventsPerPool","hEventsPerPool",fNzVertPools,fzVertPoolLims,fNMultPools,fMultPoolLims);
    fMixingsPerPool=new TH2F("hMixingsPerPool","hMixingsPerPool",fNzVertPools,fzVertPoolLims,fNMultPools,fMultPoolLims);
  }else{
    fEventsPerPool=new TH2F("hEventsPerPool","hEventsPerPool",1,-10.,10.,1,-0.5,2000.5);
    fMixingsPerPool=new TH2F("hMixingsPerPool","hMixingsPerPool",1,-10.,10.,1,-0.5,2000.5);
  }
  fOutput->Add(fEventsPerPool);
  fOutput->Add(fMixingsPerPool);

  fMassVsPtVsCosthSt=new TH3F("hMassVsPtVsCosthSt","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fMassVsPtVsCosthStRot=new TH3F("hMassVsPtVsCosthStRot","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fMassVsPtVsCosthStLSpp=new TH3F("hMassVsPtVsCosthStLSpp","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fMassVsPtVsCosthStLSmm=new TH3F("hMassVsPtVsCosthStLSmm","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fMassVsPtVsCosthStSig=new TH3F("hMassVsPtVsCosthStSig","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fMassVsPtVsCosthStRefl=new TH3F("hMassVsPtVsCosthStRefl","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fMassVsPtVsCosthStBkg=new TH3F("hMassVsPtVsCosthStBkg","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fMassVsPtVsCosthStME=new TH3F("hMassVsPtVsCosthStME","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fMassVsPtVsCosthStMELSpp=new TH3F("hMassVsPtVsCosthStMELSpp","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fMassVsPtVsCosthStMELSmm=new TH3F("hMassVsPtVsCosthStMELSmm","",nMassBins,fMinMass,maxm,nPtBins,0.,maxPt,6,0.4,1.);
  fOutput->Add(fMassVsPtVsCosthSt);
  fOutput->Add(fMassVsPtVsCosthStRot);
  fOutput->Add(fMassVsPtVsCosthStLSpp);
  fOutput->Add(fMassVsPtVsCosthStLSmm);
  fOutput->Add(fMassVsPtVsCosthStSig);
  fOutput->Add(fMassVsPtVsCosthStRefl);
  fOutput->Add(fMassVsPtVsCosthStBkg);
  fOutput->Add(fMassVsPtVsCosthStME);
  fOutput->Add(fMassVsPtVsCosthStMELSpp);
  fOutput->Add(fMassVsPtVsCosthStMELSmm);

  fHistonSigmaTPCPion=new TH2F("hnSigmaTPCPion"," ; p (GeV/c) ; n#sigma^{#pi}_{TPC}",20,0.,10.,100,-5.,5.);
  fHistonSigmaTPCPionGoodTOF=new TH2F("hnSigmaTPCPionGoodTOF"," ; p (GeV/c) ; n#sigma^{#pi}_{TPC}",20,0.,10.,100,-5.,5.);
  fHistonSigmaTOFPion=new TH2F("hnSigmaTOFPion"," ; p (GeV/c) ; n#sigma^{#pi}_{TOF}",20,0.,10.,100,-5.,5.);
  fHistonSigmaTPCKaon=new TH2F("hnSigmaTPCKaon"," ; p (GeV/c) ; n#sigma^{K}_{TPC}",20,0.,10.,100,-5.,5.);
  fHistonSigmaTPCKaonGoodTOF=new TH2F("hnSigmaTPCKaonGoodTOF"," ; p (GeV/c) ; n#sigma^{K}_{TPC}",20,0.,10.,100,-5.,5.);
  fHistonSigmaTOFKaon=new TH2F("hnSigmaTOFKaon"," ; p (GeV/c) ; n#sigma^{K}_{TOF}",20,0.,10.,100,-5.,5.);
  fOutput->Add(fHistonSigmaTPCPion);
  fOutput->Add(fHistonSigmaTPCPionGoodTOF);
  fOutput->Add(fHistonSigmaTOFPion);
  fOutput->Add(fHistonSigmaTPCKaon);
  fOutput->Add(fHistonSigmaTPCKaonGoodTOF);
  fOutput->Add(fHistonSigmaTOFKaon);

  fHistoPtKPtPiPtD = new TH3F("hPtKPtPiPtD"," ; p_{T}(D) ; p_{T}(K) ; p_{T}(#pi)",32,0.,16.,60,0.,15.,60,0.,15.);
  fHistoPtKPtPiPtDSig = new TH3F("hPtKPtPiPtDSig"," ; p_{T}(D) ; p_{T}(K) ; p_{T}(#pi)",32,0.,16.,60,0.,15.,60,0.,15.);
  fOutput->Add(fHistoPtKPtPiPtD);
  fOutput->Add(fHistoPtKPtPiPtDSig);
  
  //Counter for Normalization
  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();

  fListCuts = new TList();
  fListCuts->SetOwner();
  if(fTrackCutsAll){
    AliESDtrackCuts* tatosave=new AliESDtrackCuts(*fTrackCutsAll);
    fListCuts->Add(tatosave);
  }
  if(fTrackCutsPion){
    AliESDtrackCuts* tptosave=new AliESDtrackCuts(*fTrackCutsPion);
    tptosave->SetName(Form("%sForPions",fTrackCutsPion->GetName()));
    fListCuts->Add(tptosave);
  }
  if(fTrackCutsKaon){
    AliESDtrackCuts* tktosave=new AliESDtrackCuts(*fTrackCutsKaon);
    tktosave->SetName(Form("%sForKaons",fTrackCutsKaon->GetName()));
    fListCuts->Add(tktosave);
  }
  
  if(fAnalysisCuts->GetPidHF()){
    AliAODPidHF* pidtosave=new AliAODPidHF(*(fAnalysisCuts->GetPidHF()));
    fListCuts->Add(pidtosave);
  }
  TH1F* hCutValues = new TH1F("hCutValues","",7,0.5,7.5);
  hCutValues->SetBinContent(1,fFilterMask);
  hCutValues->GetXaxis()->SetBinLabel(1,"Filter bit");
  hCutValues->SetBinContent(2,(Float_t)fApplyCutCosThetaStar);
  hCutValues->GetXaxis()->SetBinLabel(2,"Use costhetastar (D0)");
  hCutValues->SetBinContent(3,fCutCosThetaStar);
  hCutValues->GetXaxis()->SetBinLabel(3,"costhetastar (D0)");
  hCutValues->SetBinContent(4,fPhiMassCut);
  hCutValues->GetXaxis()->SetBinLabel(4,"phi mass (Ds)");
  hCutValues->SetBinContent(5,fCutCos3PiKPhiRFrame);
  hCutValues->GetXaxis()->SetBinLabel(5,"cos3piK (Ds)");
  hCutValues->SetBinContent(6,fCutCosPiDsLabFrame);
  hCutValues->GetXaxis()->SetBinLabel(6,"cospiDs (Ds)");
  hCutValues->SetBinContent(7,fAnalysisCuts->GetUseTimeRangeCutForPbPb2018());
  hCutValues->GetXaxis()->SetBinLabel(7,"TimeRangeCut");
  fListCuts->Add(hCutValues);
  PostData(3, fListCuts);

  
  fKaonTracks = new TObjArray();
  fPionTracks=new TObjArray();
  fKaonTracks->SetOwner();
  fPionTracks->SetOwner();

  fEventBuffer = new TTree*[fNOfPools];
  for(Int_t i=0; i<fNOfPools; i++){
    fEventBuffer[i]=new TTree(Form("EventBuffer_%d",i), "Temporary buffer for event mixing");
    fEventBuffer[i]->Branch("zVertex", &fVtxZ);
    fEventBuffer[i]->Branch("multiplicity", &fMultiplicity);
    fEventBuffer[i]->Branch("eventInfo", "TObjString",&fEventInfo);
    fEventBuffer[i]->Branch("karray", "TObjArray", &fKaonTracks);
    fEventBuffer[i]->Branch("parray", "TObjArray", &fPionTracks);
  }

  PostData(1,fOutput);
  PostData(2,fCounter);
}

//________________________________________________________________________
void AliAnalysisTaskCombinHF::UserExec(Option_t */*option*/){
  /// Build the 3-track combinatorics (+-+ and -+-) for D+->Kpipi decays
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
  }
  if(!aod){
    printf("AliAnalysisTaskCombinHF::UserExec: AOD not found!\n");
    return;
  }
  
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  if(fReadMC){
    // Reject events with trigger mask 0 of the LHC13d3 production
    // For these events the ITS layers are skipped in the trakcing
    // and the vertex reconstruction efficiency from tracks is biased
    Int_t runnumber = aod->GetRunNumber();
    if(aod->GetTriggerMask()==0 &&
       (runnumber>=195344 && runnumber<=195677)){
      return;
    }
    // Set the trigger mask for physics selection to kMB in the MC
    if(fEnforceMBTrigMaskInMC){
      //      printf("Enforce trigger mask to kMB, previous mask = %d\n",fAnalysisCuts->GetTriggerMask());
      fAnalysisCuts->SetTriggerMask(AliVEvent::kMB);
    }
    for(Int_t j=0; j<200000; j++) fOrigContainer[j]=-1;
  }

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  
  
  fHistNEvents->Fill(0); // count event
  // Post the data already here
  PostData(1,fOutput);
  
  fCounter->StoreEvent(aod,fAnalysisCuts,fReadMC);
  
  Bool_t isEvSel=fAnalysisCuts->IsEventSelected(aod);
  if(fAnalysisCuts->IsEventRejectedDueToTrigger() || fAnalysisCuts->IsEventRejectedDuePhysicsSelection()){
    if(fAnalysisCuts->IsEventRejectedDueToTrigger()) fHistNEvents->Fill(2);
    if(fAnalysisCuts->IsEventRejectedDuePhysicsSelection()) fHistNEvents->Fill(3);
  }else{
    if(fAnalysisCuts->IsEventRejectedDueToCentrality()){
      fHistNEvents->Fill(9);
    }else{
      if(fAnalysisCuts->IsEventRejectedDueToBadPrimaryVertex()){
	if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex())fHistNEvents->Fill(4);
	if(fAnalysisCuts->IsEventRejectedDueToVertexContributors())fHistNEvents->Fill(5);
	if(fAnalysisCuts->IsEventRejectedDueToBadTrackVertex())fHistNEvents->Fill(6);
      }else{
	if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())fHistNEvents->Fill(7);
	else if(fAnalysisCuts->IsEventRejectedDueToPileup())fHistNEvents->Fill(8);
      }
    }
  }

  // PID object should be taken AFTER IsEventSelected to have proper call to SetupPid !!
  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
  fPidHF = fAnalysisCuts->GetPidHF();
  fPidHF->SetPidResponse(pidResp);
  //
  
  Int_t ntracks=aod->GetNumberOfTracks();
  fVtxZ = aod->GetPrimaryVertex()->GetZ();
  fMultiplicity = AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.); 
  Float_t evCentr=fAnalysisCuts->GetCentrality(aod);
  if(!fAnalysisCuts->IsEventRejectedDueToTrigger() && !fAnalysisCuts->IsEventRejectedDuePhysicsSelection() &&
     !fAnalysisCuts->IsEventRejectedDueToBadPrimaryVertex() && !fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()){
    fHistEventMultCent->Fill(evCentr,fMultiplicity);
  }

  if(fAnalysisCuts->GetUseCentrality()>0 && fAnalysisCuts->IsEventSelectedInCentrality(aod)!=0) return;
  // events not passing the centrality selection can be removed immediately. For the others we must count the generated D mesons


  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;
  if(fReadMC){
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskCombinHF::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskCombinHF::UserExec: MC header branch not found!\n");
      return;
    }
    // selection on pt hard bins in Pb-Pb
    if(fSelectPtHardRange){
      TList *lh=mcHeader->GetCocktailHeaders();
      if(lh){
	Int_t nh=lh->GetEntries();
	for(Int_t i=0;i<nh;i++){
	  AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
	  TString genname=gh->GetName();
	  if(genname.Contains("ythia") || genname.Contains("YTHIA")){
	    AliGenPythiaEventHeader* pyth=(AliGenPythiaEventHeader*)lh->At(i);
	    Double_t ptha=pyth->GetPtHard();
	    Double_t xsec=pyth->GetXsection();
	    if(ptha<fMinPtHard || ptha>fMaxPtHard) return;
	    fHistXsecVsPtHard->SetBinContent(fHistXsecVsPtHard->GetXaxis()->FindBin(ptha),xsec);
	  }
	}
      }
    }
    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) < fAnalysisCuts->GetMaxVtxZ()){ // only cut on zVertex applied to count the signal
      FillGenHistos(arrayMC,mcHeader,isEvSel);
    }
    fHistEventMultZv->Fill(zMCVertex,fMultiplicity);
    if(isEvSel) fHistEventMultZvEvSel->Fill(zMCVertex,fMultiplicity);
    // switch off event mixing in case of signal only MC
    if(fSignalOnlyMC) fDoEventMixing=0;
  }else{
    fHistEventMultZv->Fill(fVtxZ,fMultiplicity);
    if(isEvSel) fHistEventMultZvEvSel->Fill(fVtxZ,fMultiplicity);
  }


  if(!isEvSel)return;
  
  fHistNEvents->Fill(1);
  fHistEventMultCentEvSel->Fill(evCentr,fMultiplicity);


  Int_t pdgOfD=421;
  if(fMeson==kDplus) pdgOfD=411;
  else if(fMeson==kDs) pdgOfD=431;

  // select and flag tracks
  UChar_t* status = new UChar_t[ntracks];
  for(Int_t iTr=0; iTr<ntracks; iTr++){
    status[iTr]=0;
    AliAODTrack* track=dynamic_cast<AliAODTrack*>(aod->GetTrack(iTr));
    if(!track){
      AliWarning("Error in casting track to AOD track. Not a standard AOD?");
      continue;
    }
    if(fReadMC && fSignalOnlyMC && arrayMC){
      // for fast MC analysis we skip tracks not coming from charm hadrons
      Bool_t isCharm=AliVertexingHFUtils::IsTrackFromHadronDecay(pdgOfD,track,arrayMC);
      if(!isCharm) continue;
    }
    if(IsTrackSelected(track)) status[iTr]+=1;
    
    // PID
    if (fPIDstrategy == knSigma) {
      // nsigma PID
      Double_t trmom=track->P();
      Bool_t okTOF=fPidHF->CheckTOFPIDStatus(track);
      if(IsKaon(track)){
	Double_t nstpc,nstof;
	fPidHF->GetnSigmaTPC(track,AliPID::kKaon,nstpc);
	fPidHF->GetnSigmaTOF(track,AliPID::kKaon,nstof);
	fHistonSigmaTPCKaon->Fill(trmom,nstpc);
	if(okTOF) fHistonSigmaTPCKaonGoodTOF->Fill(trmom,nstpc);
	fHistonSigmaTOFKaon->Fill(trmom,nstof);
	status[iTr]+=2;
      }
      if(IsPion(track)){
	Double_t nstpc,nstof;
	fPidHF->GetnSigmaTPC(track,AliPID::kPion,nstpc);
	fPidHF->GetnSigmaTOF(track,AliPID::kPion,nstof);
	fHistonSigmaTPCPion->Fill(trmom,nstpc);
	if(okTOF) fHistonSigmaTPCPionGoodTOF->Fill(trmom,nstpc);
	fHistonSigmaTOFPion->Fill(trmom,nstof);
	status[iTr]+=4;
      }
    }
    else if (fPIDstrategy == kBayesianMaxProb || fPIDstrategy == kBayesianThres) {
      // Bayesian PID
      Double_t *weights = new Double_t[AliPID::kSPECIES];
      fPidHF->GetPidCombined()->ComputeProbabilities(track, fPidHF->GetPidResponse(), weights);
      if (fPIDstrategy == kBayesianMaxProb) {
        if (TMath::MaxElement(AliPID::kSPECIES, weights) == weights[AliPID::kKaon]) status[iTr] += 2;
        if (TMath::MaxElement(AliPID::kSPECIES, weights) == weights[AliPID::kPion]) status[iTr] += 4;
      }
      if (fPIDstrategy == kBayesianThres) {
        if (weights[AliPID::kKaon] > fBayesThresKaon) status[iTr] += 2;
        if (weights[AliPID::kPion] > fBayesThresPion) status[iTr] += 4;
      }
      delete[] weights;
    }
    
    fHistTrackStatus->Fill(status[iTr]);
    fHistTrackEtaMultZv->Fill(track->Eta(),fVtxZ,fMultiplicity);
  }
  
  // build the combinatorics
  Int_t nSelected=0;
  Int_t nFiltered=0;
  Double_t dummypos[3]={0.,0.,0.};
  AliAODVertex* v2=new AliAODVertex(dummypos,999.,-1,2);
  AliAODVertex* v3=new AliAODVertex(dummypos,999.,-1,3);
  // dummy values of track impact parameter, needed to build an AliAODRecoDecay object
  Double_t d02[2]={0.,0.};
  Double_t d03[3]={0.,0.,0.};
  AliAODRecoDecay* tmpRD2 = new AliAODRecoDecay(0x0,2,0,d02);
  AliAODRecoDecay* tmpRD3 = new AliAODRecoDecay(0x0,3,1,d03);
  UInt_t pdg0[2]={321,211};
  UInt_t pdgp[3]={321,211,211};
  UInt_t pdgs[3]={321,211,321};
  Double_t tmpp[3];
  Double_t px[3],py[3],pz[3];
  Int_t dgLabels[3];
  fKaonTracks->Delete();
  fPionTracks->Delete();
 
  for(Int_t iTr1=0; iTr1<ntracks; iTr1++){
    AliAODTrack* trK=dynamic_cast<AliAODTrack*>(aod->GetTrack(iTr1));
    if(!trK){
      AliWarning("Error in casting track to AOD track. Not a standard AOD?");
      continue;
    }
    if((status[iTr1] & 1)==0) continue;
    if(fDoEventMixing>0){
      if(status[iTr1] & 2) fKaonTracks->AddLast(new TLorentzVector(trK->Px(),trK->Py(),trK->Pz(),trK->Charge()));
      if(status[iTr1] & 4) fPionTracks->AddLast(new TLorentzVector(trK->Px(),trK->Py(),trK->Pz(),trK->Charge()));
    }
    if((status[iTr1] & 2)==0) continue;
    Int_t chargeK=trK->Charge();
    trK->GetPxPyPz(tmpp);
    px[0] = tmpp[0];
    py[0] = tmpp[1];
    pz[0] = tmpp[2];
    dgLabels[0]=trK->GetLabel();
    for(Int_t iTr2=0; iTr2<ntracks; iTr2++){
      if((status[iTr2] & 1)==0) continue;
      if((status[iTr2] & 4)==0) continue;
      if(iTr1==iTr2) continue;
      AliAODTrack* trPi1=dynamic_cast<AliAODTrack*>(aod->GetTrack(iTr2));
      if(!trPi1){
	AliWarning("Error in casting track to AOD track. Not a standard AOD?");
	continue;
      }
      Int_t chargePi1=trPi1->Charge();
      trPi1->GetPxPyPz(tmpp);
      px[1] = tmpp[0];
      py[1] = tmpp[1];
      pz[1] = tmpp[2];
      dgLabels[1]=trPi1->GetLabel();
      if(fMeson==kDzero){
        if(chargePi1==chargeK){
	  // LS candidate
	  FillLSHistos(421,2,tmpRD2,px,py,pz,pdg0,chargePi1);
	}else{
	  // OS candidate
	  nFiltered++;
	  v2->AddDaughter(trK);
	  v2->AddDaughter(trPi1);
	  tmpRD2->SetSecondaryVtx(v2);
	  Bool_t ok=FillHistos(421,2,tmpRD2,px,py,pz,pdg0,arrayMC,mcHeader,dgLabels);
	  v2->RemoveDaughters();
	  if(ok) nSelected++;
	}
      }else{
        for(Int_t iTr3=iTr2+1; iTr3<ntracks; iTr3++){
          if((status[iTr3] & 1)==0) continue;
          if(fMeson==kDplus && (status[iTr3] & 4)==0) continue;
          if(fMeson==kDs && (status[iTr3] & 2)==0) continue;
          if(iTr1==iTr3) continue;
          AliAODTrack* trPi2=dynamic_cast<AliAODTrack*>(aod->GetTrack(iTr3));
	  if(!trPi2){
	    AliWarning("Error in casting track to AOD track. Not a standard AOD?");
	    continue;
	  }
          Int_t chargePi2=trPi2->Charge();
          trPi2->GetPxPyPz(tmpp);
          px[2] = tmpp[0];
          py[2] = tmpp[1];
          pz[2] = tmpp[2];
          dgLabels[2]=trPi2->GetLabel();
	  if(fMeson==kDs){
	    Double_t massKK=ComputeInvMassKK(trK,trPi2);
	    Double_t deltaMass=massKK-TDatabasePDG::Instance()->GetParticle(333)->Mass();
	    if(TMath::Abs(deltaMass)>fPhiMassCut) continue;
	    tmpRD3->SetPxPyPzProngs(3,px,py,pz);
	    Double_t cos1=((AliAODRecoDecayHF3Prong*)tmpRD3)->CosPiKPhiRFrameKpiK();
	    Double_t kincutPiKPhi=TMath::Abs(cos1*cos1*cos1);
	    if(kincutPiKPhi<fCutCos3PiKPhiRFrame) continue;
	    Double_t cosPiDsLabFrame=((AliAODRecoDecayHF3Prong*)tmpRD3)->CosPiDsLabFrameKpiK();
	    if(cosPiDsLabFrame>fCutCosPiDsLabFrame) continue;
	  }
	  Bool_t isThreeLS=kFALSE;
	  if(chargePi1==chargeK && chargePi2==chargeK){
	    isThreeLS=kTRUE;
	    if(fMeson==kDplus)FillLSHistos(411,3,tmpRD3,px,py,pz,pdgp,chargePi1);
	    else if(fMeson==kDs)FillLSHistos(431,3,tmpRD3,px,py,pz,pdgs,chargePi1);
	  }
	  Bool_t acceptOS=kFALSE;	  
	  if(fMeson==kDplus){
	    if(chargePi1!=chargeK && chargePi2!=chargeK)acceptOS=kTRUE;
	  }else if(fMeson==kDs){
	    if(chargePi2!=chargeK && !isThreeLS) acceptOS=kTRUE;
	  }
	  if(acceptOS){
	    nFiltered++;
  	    v3->AddDaughter(trK);
	    v3->AddDaughter(trPi1);
	    v3->AddDaughter(trPi2);
	    tmpRD3->SetSecondaryVtx(v3);
	    Bool_t ok=kFALSE;
	    if(fMeson==kDplus) ok=FillHistos(411,3,tmpRD3,px,py,pz,pdgp,arrayMC,mcHeader,dgLabels);
	    else if(fMeson==kDs) ok=FillHistos(431,3,tmpRD3,px,py,pz,pdgs,arrayMC,mcHeader,dgLabels);
	    v3->RemoveDaughters();
	    if(ok) nSelected++;
	  }
	}
      }
    }
  }
  
  delete [] status;
  delete v2;
  delete v3;
  delete tmpRD2;
  delete tmpRD3;
  
  fNSelected->Fill(nSelected);
  
  fCounter->StoreCandidates(aod,nFiltered,kTRUE);
  fCounter->StoreCandidates(aod,nSelected,kFALSE);
  fEventInfo->SetString(Form("Ev%d_esd%d_Pi%d_K%d",mgr->GetNcalls(),((AliAODHeader*)aod->GetHeader())->GetEventNumberESDFile(),fPionTracks->GetEntries(),fKaonTracks->GetEntries()));
  if(fDoEventMixing==1){
    Int_t ind=GetPoolIndex(fVtxZ,fMultiplicity);
    if(ind>=0 && ind<fNOfPools){
      fEventsPerPool->Fill(fVtxZ,fMultiplicity);
      fEventBuffer[ind]->Fill();
      if(fEventBuffer[ind]->GetEntries() >= fNumberOfEventsForMixing){
	fMixingsPerPool->Fill(fVtxZ,fMultiplicity);
	  DoMixingWithPools(ind);
	  ResetPool(ind);
      }
    }
  }else if(fDoEventMixing==2){ // mix with cuts, no pools
      fEventBuffer[0]->Fill();
  }
  PostData(1,fOutput);
  PostData(2,fCounter);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskCombinHF::FillLSHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, Int_t charge){
  /// Fill histos for LS candidates

  if(fReadMC && fSignalOnlyMC) return;
  tmpRD->SetPxPyPzProngs(nProngs,px,py,pz);
  Double_t pt = tmpRD->Pt();
  Double_t minv2 = tmpRD->InvMass2(nProngs,pdgdau);
  if(minv2>fMinMass*fMinMass && minv2<fMaxMass*fMaxMass){
    Double_t rapid = tmpRD->Y(pdgD);
    if(fAnalysisCuts->IsInFiducialAcceptance(pt,rapid)){
      Bool_t fillLS=kTRUE;
      Double_t costhst=0;
      Double_t absCosThSt=0;
      if(TMath::Abs(pdgD)==421 && (fApplyCutCosThetaStar || fFillHistosVsCosThetaStar)){
	costhst=tmpRD->CosThetaStar(0,421,321,211); // kaon is the first daughter
	absCosThSt=TMath::Abs(costhst);
	if(fApplyCutCosThetaStar && absCosThSt>fCutCosThetaStar) fillLS=kFALSE;
      }
      if(fillLS){
	Double_t invMass=TMath::Sqrt(minv2);
	if(charge>0) fMassVsPtVsYLSpp->Fill(invMass,pt,rapid);
	else fMassVsPtVsYLSmm->Fill(invMass,pt,rapid);
	if(fFillHistosVsCosThetaStar){
	  if(charge>0) fMassVsPtVsCosthStLSpp->Fill(invMass,pt,absCosThSt);
	  else fMassVsPtVsCosthStLSmm->Fill(invMass,pt,absCosThSt);
	}
      }
    }
  }
  return;
}

//________________________________________________________________________
void AliAnalysisTaskCombinHF::FillGenHistos(TClonesArray* arrayMC, AliAODMCHeader *mcHeader, Bool_t isEvSel){
  /// Fill histos with generated quantities
  Int_t totPart=arrayMC->GetEntriesFast();
  Int_t thePDG=411;
  Int_t nProng=3;
  if(fMeson==kDzero){
    thePDG=421;
    nProng=2;
  }else if(fMeson==kDs){
    thePDG=431;
    nProng=3;
   }
  for(Int_t ip=0; ip<totPart; ip++){
    AliAODMCParticle *part = (AliAODMCParticle*)arrayMC->At(ip);
    if(TMath::Abs(part->GetPdgCode())==thePDG){
      Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,part,fGoUpToQuark);
      if(ip<200000) fOrigContainer[ip]=orig;
      Bool_t isInj=AliVertexingHFUtils::IsTrackInjected(ip,mcHeader,arrayMC);
      Int_t deca=0;
      Bool_t isGoodDecay=kFALSE;
      Int_t labDau[4]={-1,-1,-1,-1};
      if(fMeson==kDzero){
        deca=AliVertexingHFUtils::CheckD0Decay(arrayMC,part,labDau);
        if(part->GetNDaughters()!=2) continue;
        if(deca==1) isGoodDecay=kTRUE;
      }else if(fMeson==kDplus){
        deca=AliVertexingHFUtils::CheckDplusDecay(arrayMC,part,labDau);
        if(deca>0) isGoodDecay=kTRUE;
      }else if(fMeson==kDs){
        deca=AliVertexingHFUtils::CheckDsDecay(arrayMC,part,labDau);
        if(deca==1) isGoodDecay=kTRUE;
      }
      if(labDau[0]==-1){
        //	printf(Form("Meson %d Label of daughters not filled correctly -- %d\n",fMeson,isGoodDecay));
        continue; //protection against unfilled array of labels
      }
      fHistCheckDecChan->Fill(deca);
      Bool_t isInAcc=CheckAcceptance(arrayMC,nProng,labDau);
      if(isInAcc) fHistCheckDecChanAcc->Fill(deca);
      if(isGoodDecay){
        Double_t ptgen=part->Pt();
        Double_t ygen=part->Y();
	Double_t ptbmoth=0.;
	if(orig==5) ptbmoth=AliVertexingHFUtils::GetBeautyMotherPt(arrayMC,part);
	if(fAnalysisCuts->IsInFiducialAcceptance(ptgen,ygen)){
	  fHistCheckOrigin->Fill(orig,isInj);
	  if(orig==4){
	    fPtVsYVsMultGenPrompt->Fill(ptgen,ygen,fMultiplicity);
	    if(TMath::Abs(ygen)<0.5) fPtVsYVsMultGenLimAccPrompt->Fill(ptgen,ygen,fMultiplicity);
	    if(isInAcc) fPtVsYVsMultGenAccPrompt->Fill(ptgen,ygen,fMultiplicity);
	    if(isEvSel && isInAcc) fPtVsYVsMultGenAccEvSelPrompt->Fill(ptgen,ygen,fMultiplicity);
	  }else if(orig==5){
	    fPtVsYVsMultGenFeeddw->Fill(ptgen,ygen,fMultiplicity);
	    fPtVsYVsPtBGenFeeddw->Fill(ptgen,ygen,ptbmoth);
	    if(TMath::Abs(ygen)<0.5){
	      fPtVsYVsMultGenLimAccFeeddw->Fill(ptgen,ygen,fMultiplicity);
	      fPtVsYVsPtBGenLimAccFeeddw->Fill(ptgen,ygen,ptbmoth);
	      fBMohterPtGen->Fill(ptbmoth);
	    }
	    if(isInAcc){
	      fPtVsYVsMultGenAccFeeddw->Fill(ptgen,ygen,fMultiplicity);
	      fPtVsYVsPtBGenAccFeeddw->Fill(ptgen,ygen,ptbmoth);
	    }
	    if(isEvSel && isInAcc){
	      fPtVsYVsMultGenAccEvSelFeeddw->Fill(ptgen,ygen,fMultiplicity);
	      fPtVsYVsMultGenAccEvSelFeeddw->Fill(ptgen,ygen,ptbmoth);
	    }
	  }
	}
        if(TMath::Abs(ygen)<0.9){
	  if(orig==4) fPtVsYVsMultGenLargeAccPrompt->Fill(ptgen,ygen,fMultiplicity);
	  else if(orig==5){
	    fPtVsYVsMultGenLargeAccFeeddw->Fill(ptgen,ygen,fMultiplicity);
	    fPtVsYVsPtBGenLargeAccFeeddw->Fill(ptgen,ygen,ptbmoth);
	  }
	}
      }
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::FillHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, TClonesArray *arrayMC, AliAODMCHeader *mcHeader, Int_t* dgLabels){
  /// Fill histos for candidates with proper charge sign
  
  Bool_t accept=kFALSE;
  
  tmpRD->SetPxPyPzProngs(nProngs,px,py,pz);
  Double_t pt = tmpRD->Pt();
  Double_t minv2 = tmpRD->InvMass2(nProngs,pdgdau);
  Double_t mass=TMath::Sqrt(minv2);
  
  if(minv2>fMinMass*fMinMass && minv2<fMaxMass*fMaxMass){
    Double_t rapid = tmpRD->Y(pdgD);
    if(fAnalysisCuts->IsInFiducialAcceptance(pt,rapid)){
      accept=kTRUE;
      Double_t costhst=0;
      Double_t absCosThSt=0;
      Double_t ptK=0;
      Double_t ptPi=0;
      if(TMath::Abs(pdgD)==421 && (fApplyCutCosThetaStar || fFillHistosVsCosThetaStar)){
	costhst=tmpRD->CosThetaStar(0,421,321,211); // kaon is the first daughter
	absCosThSt=TMath::Abs(costhst);
	if(fApplyCutCosThetaStar && absCosThSt>fCutCosThetaStar) accept=kFALSE;
      }
      if(accept){
	fMassVsPtVsY->Fill(mass,pt,rapid);
	if(fFillHistosVsCosThetaStar) fMassVsPtVsCosthSt->Fill(mass,pt,absCosThSt);
	if(pdgD==421){
	  ptK=TMath::Sqrt(px[0]*px[0]+py[0]*py[0]);
	  ptPi=TMath::Sqrt(px[1]*px[1]+py[1]*py[1]);
	  if(TMath::Abs(mass-1.865)<0.025) fHistoPtKPtPiPtD->Fill(pt,ptK,ptPi);
	}
	if(fReadMC){
	  Int_t signPdg[3]={0,0,0};
	  for(Int_t iii=0; iii<nProngs; iii++) signPdg[iii]=pdgdau[iii];
	  Int_t labD = tmpRD->MatchToMC(pdgD,arrayMC,nProngs,signPdg);
	  if(labD>=0){
	    AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(TMath::Abs(dgLabels[0])));
	    if(part){
	      Int_t pdgCode = TMath::Abs( part->GetPdgCode() );
	      if(pdgCode==321){ // if the first daughter is a Kaon, this is signal with correct mass assignment
		fMassVsPtVsYSig->Fill(mass,pt,rapid);
		if(fFillHistosVsCosThetaStar) fMassVsPtVsCosthStSig->Fill(mass,pt,absCosThSt);
		if(pdgD==421) fHistoPtKPtPiPtDSig->Fill(pt,ptK,ptPi);
		AliAODMCParticle* dmes =  dynamic_cast<AliAODMCParticle*>(arrayMC->At(labD));
		if(dmes){
		  Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,dmes,fGoUpToQuark);
		  Bool_t isInj=AliVertexingHFUtils::IsTrackInjected(labD,mcHeader,arrayMC);
		  fHistCheckOriginRecoD->Fill(orig,isInj);
		  if(labD<200000) fHistCheckOriginRecoVsGen->Fill(fOrigContainer[labD],orig);
		  if(orig==4) fPtVsYVsMultRecoPrompt->Fill(dmes->Pt(),dmes->Y(),fMultiplicity);
		  else if(orig==5){
		    Double_t ptbmoth=AliVertexingHFUtils::GetBeautyMotherPt(arrayMC,dmes);
		    fPtVsYVsMultRecoFeeddw->Fill(dmes->Pt(),dmes->Y(),fMultiplicity);
		    fPtVsYVsPtBRecoFeeddw->Fill(dmes->Pt(),dmes->Y(),ptbmoth);
		  }
		}
	      }else{ // if the first daughter is not a kaon, it is a reflection
		fMassVsPtVsYRefl->Fill(mass,pt,rapid);
		if(fFillHistosVsCosThetaStar) fMassVsPtVsCosthStRefl->Fill(mass,pt,absCosThSt);
	      }
	    }
	  }else{
	    if(fSignalOnlyMC) accept=kFALSE;
	    else fMassVsPtVsYBkg->Fill(mass,pt,rapid);
	    if(fFillHistosVsCosThetaStar) fMassVsPtVsCosthStBkg->Fill(mass,pt,absCosThSt);
	  }
	}
      }
    }
  }
  // skip track rotations in case of signal only MC
  if(fReadMC && fSignalOnlyMC) return accept;

  // Track rotations to estimate the background
  Int_t nRotated=0;
  Double_t massRot=0;// calculated later only if candidate is acceptable
  Double_t angleProngXY;
  if(TMath::Abs(pdgD)==421)angleProngXY=TMath::ACos((px[0]*px[1]+py[0]*py[1])/TMath::Sqrt((px[0]*px[0]+py[0]*py[0])*(px[1]*px[1]+py[1]*py[1])));
  else if(TMath::Abs(pdgD)==431) {
    Double_t px_phi = px[0]+px[2];
    Double_t py_phi = py[0]+py[2];
    Double_t pz_phi = pz[0]+pz[2];
    angleProngXY=TMath::ACos((pz_phi*px[1]+py_phi*py[1])/TMath::Sqrt((px_phi*px_phi+py_phi*py_phi)*(px[1]*px[1]+py[1]*py[1])));
  }
  else {//angle between pion and phi meson
    angleProngXY=TMath::ACos(((px[0]+px[1])*px[2]+(py[0]+py[1])*py[2])/TMath::Sqrt(((px[0]+px[1])*(px[0]+px[1])+(py[0]+py[1])*(py[0]+py[1]))*(px[2]*px[2]+py[2]*py[2])));
  }
  Double_t ptOrig=pt;
  
  
  Double_t rotStep=0.;
  if(fNRotations>1) rotStep=(fMaxAngleForRot-fMinAngleForRot)/(fNRotations-1); // -1 is to ensure that the last rotation is done with angle=fMaxAngleForRot
  if(TMath::Abs(pdgD)==421 || TMath::Abs(pdgD)==431) fNRotations3=1;
  Double_t rotStep3=0.;
  if(fNRotations3>1) rotStep3=(fMaxAngleForRot3-fMinAngleForRot3)/(fNRotations3-1); // -1 is to ensure that the last rotation is done with angle=fMaxAngleForRot

  for(Int_t irot=0; irot<fNRotations; irot++){
    Double_t phirot=fMinAngleForRot+rotStep*irot;
    Double_t tmpx=px[0];
    Double_t tmpy=py[0];
    Double_t tmpx2=px[2];
    Double_t tmpy2=py[2];
    if(pdgD==431) {
      //rotate pion w.r.t. phi meson
      tmpx=px[1];
      tmpy=py[1];
      px[1]=tmpx*TMath::Cos(phirot)-tmpy*TMath::Sin(phirot);
      py[1]=tmpx*TMath::Sin(phirot)+tmpy*TMath::Cos(phirot);
    }
    else {
      px[0]=tmpx*TMath::Cos(phirot)-tmpy*TMath::Sin(phirot);
      py[0]=tmpx*TMath::Sin(phirot)+tmpy*TMath::Cos(phirot);
    }
    for(Int_t irot3=0; irot3<fNRotations3; irot3++){
      if(pdgD==411){
	Double_t phirot2=fMaxAngleForRot3-rotStep3*irot;
	px[2]=tmpx*TMath::Cos(phirot2)-tmpy*TMath::Sin(phirot2);
	py[2]=tmpx*TMath::Sin(phirot2)+tmpy*TMath::Cos(phirot2);
      }
      tmpRD->SetPxPyPzProngs(nProngs,px,py,pz);
      pt = tmpRD->Pt();
      minv2 = tmpRD->InvMass2(nProngs,pdgdau);
      if(minv2>fMinMass*fMinMass && minv2<fMaxMass*fMaxMass){
	Double_t rapid = tmpRD->Y(pdgD);
	if(fAnalysisCuts->IsInFiducialAcceptance(pt,rapid)){
	  Bool_t fillRotCase=kTRUE;
	  Double_t costhst=0;
	  Double_t absCosThSt=0;
	  if(TMath::Abs(pdgD)==421 && (fApplyCutCosThetaStar || fFillHistosVsCosThetaStar)){
	    costhst=tmpRD->CosThetaStar(0,421,321,211); // kaon is the first daughter
	    absCosThSt=TMath::Abs(costhst);
	    if(fApplyCutCosThetaStar && absCosThSt>fCutCosThetaStar) fillRotCase=kFALSE;
	  }
	  if(fillRotCase){
	    massRot=TMath::Sqrt(minv2);
	    fMassVsPtVsYRot->Fill(massRot,pt,rapid);
	    if(fFillHistosVsCosThetaStar) fMassVsPtVsCosthStRot->Fill(massRot,pt,absCosThSt);
	    nRotated++;
	    fDeltaMass->Fill(massRot-mass);
	    if(fFullAnalysis){
	      Double_t pointRot[5]={mass,massRot-mass,ptOrig,pt-ptOrig,angleProngXY};
	      fDeltaMassFullAnalysis->Fill(pointRot);
	    }
	  }
	}
      }
    }
    if(pdgD==431) {
      px[1]=tmpx;
      py[1]=tmpy;
    }
    else {
      px[0]=tmpx;
      py[0]=tmpy;
      if(pdgD==411){
        px[2]=tmpx2;
        py[2]=tmpy2;
      }
    }
  }
  fNormRotated->Fill(nRotated);
  
  return accept;
  
}
//________________________________________________________________________
void AliAnalysisTaskCombinHF::FillMEHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau){
  /// Fill histos for candidates in MixedEvents
    
  tmpRD->SetPxPyPzProngs(nProngs,px,py,pz);
  Double_t pt = tmpRD->Pt();
  Double_t minv2 = tmpRD->InvMass2(nProngs,pdgdau);
  Double_t mass=TMath::Sqrt(minv2);

  if(minv2>fMinMass*fMinMass && minv2<fMaxMass*fMaxMass){
    Double_t rapid = tmpRD->Y(pdgD);
    if(fAnalysisCuts->IsInFiducialAcceptance(pt,rapid)){
      Bool_t fillME=kTRUE;
      Double_t costhst=0;
      Double_t absCosThSt=0;
      if(TMath::Abs(pdgD)==421 && (fApplyCutCosThetaStar || fFillHistosVsCosThetaStar)){
	costhst=tmpRD->CosThetaStar(0,421,321,211); // kaon is the first daughter
	absCosThSt=TMath::Abs(costhst);
	if(fApplyCutCosThetaStar && absCosThSt>fCutCosThetaStar) fillME=kFALSE;
      }
      if(fillME){
	fMassVsPtVsYME->Fill(mass,pt,rapid);
	if(fFillHistosVsCosThetaStar) fMassVsPtVsCosthStME->Fill(mass,pt,absCosThSt);
      }
    }
  }
  return;
}
//________________________________________________________________________
void AliAnalysisTaskCombinHF::FillMEHistosLS(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, Int_t charge){
  /// Fill histos for candidates in MixedEvents
    
  tmpRD->SetPxPyPzProngs(nProngs,px,py,pz);
  Double_t pt = tmpRD->Pt();
  Double_t minv2 = tmpRD->InvMass2(nProngs,pdgdau);
  Double_t mass=TMath::Sqrt(minv2);

  if(minv2>fMinMass*fMinMass && minv2<fMaxMass*fMaxMass){
    Double_t rapid = tmpRD->Y(pdgD);
    if(fAnalysisCuts->IsInFiducialAcceptance(pt,rapid)){
      Bool_t fillME=kTRUE;
      Double_t costhst=0;
      Double_t absCosThSt=0;
      if(TMath::Abs(pdgD)==421 && (fApplyCutCosThetaStar || fFillHistosVsCosThetaStar)){
	costhst=tmpRD->CosThetaStar(0,421,321,211); // kaon is the first daughter
	absCosThSt=TMath::Abs(costhst);
	if(fApplyCutCosThetaStar && absCosThSt>fCutCosThetaStar) fillME=kFALSE;
      }
      if(fillME){
	if(charge>0) fMassVsPtVsYMELSpp->Fill(mass,pt,rapid);
	else if(charge<0) fMassVsPtVsYMELSmm->Fill(mass,pt,rapid);
	if(fFillHistosVsCosThetaStar){
	  if(charge>0) fMassVsPtVsCosthStMELSpp->Fill(mass,pt,absCosThSt);
	  else if(charge<0) fMassVsPtVsCosthStMELSmm->Fill(mass,pt,absCosThSt);
	}
      }
    }
  }
  return;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::IsTrackSelected(AliAODTrack* track){
  /// track selection cuts
  
  if(track->Charge()==0) return kFALSE;
  if(track->GetID()<0&&!fKeepNegID)return kFALSE;
  if(fFilterMask>=0){
    if(!(track->TestFilterMask(fFilterMask))) return kFALSE;
  }
  if(!SelectAODTrack(track,fTrackCutsAll)) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::IsKaon(AliAODTrack* track){
  /// kaon selection cuts
  
  if(!fPidHF) return kTRUE;
  Int_t isKaon=fPidHF->MakeRawPid(track,AliPID::kKaon);
  Double_t mom=track->P();
  if(SelectAODTrack(track,fTrackCutsKaon)) {
    if(isKaon>=1)    return kTRUE;
    if(isKaon<=-1) return kFALSE;
    switch(fPIDselCaseZero){// isKaon=0
      case 0:
      {
        return kTRUE;// always accept
      }
        break;
      case 1:
      {
        if(isKaon>=0 && track->P()>fmaxPforIDKaon)return kTRUE;// accept only if in a compatibility band starting from p=fmaxPforIDKaon
      }
        break;
      case 2:
      {
        if(track->P()>fmaxPforIDKaon)return kTRUE;
        AliPIDResponse *pidResp=fPidHF->GetPidResponse();// more elaborated strategy: asymmetric cuts, with fix momenta and nsigma ranges for the moment
        Double_t nsigma=pidResp->NumberOfSigmasTPC(track,AliPID::kKaon);
        if(nsigma>-2.&& nsigma<3. && mom<0.6)isKaon=1;
        else if(nsigma>-1.&& nsigma<3.&& mom<0.8)isKaon=1;
        if(isKaon==1)return kTRUE;
      }
        break;
      default:
      {
        AliWarning(Form("WRONG CASE OF PID STRATEGY SELECTED: %d (can range from 0 to 2)",fPIDselCaseZero));
        return kFALSE;// actually case 0 could be set as the default and return kTRUE
      }
    }
  }
  
  return kFALSE;
}
//_______________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::IsPion(AliAODTrack* track){
  /// pion selection cuts
  
  if(!fPidHF) return kTRUE;
  Int_t isPion=fPidHF->MakeRawPid(track,AliPID::kPion);
  Double_t mom=track->P();
  if(SelectAODTrack(track,fTrackCutsPion)) {
    if(isPion>=1)    return kTRUE;
    if(isPion<=-1) return kFALSE;
    switch(fPIDselCaseZero){// isPion=0
      case 0:
      {
        return kTRUE;// always accept
      }
        break;
      case 1:
      {
        if(track->P()>fmaxPforIDPion)return kTRUE;// accept only if in a compatibility band starting from p=fmaxPforIDPion
      }
        break;
      case 2:
      {
        // more elaborated strategy: asymmetric cuts, with fix momenta and nsigma ranges for the moment
        if(track->P()>fmaxPforIDPion)return kTRUE;
        AliPIDResponse *pidResp=fPidHF->GetPidResponse();
        Double_t nsigma=pidResp->NumberOfSigmasTPC(track,AliPID::kPion);
        if(nsigma<2.&& nsigma>-3. && mom<0.6)isPion=1;
        else if(nsigma<1. && nsigma> -3. && mom<0.8)isPion=1;
        if(isPion==1)return kTRUE;
      }
        break;
      default:
      {
        AliWarning(Form("WRONG CASE OF PID STRATEGY SELECTED: %d (can range from 0 to 2)",fPIDselCaseZero));
        return kFALSE;// actually case 0 could be set as the default and return kTRUE
      }
    }
  }
  
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::SelectAODTrack(AliAODTrack *track, AliESDtrackCuts *cuts){
  /// AOD track selection
  
  if(!cuts) return kTRUE;
  
  AliESDtrack esdTrack(track);
  // set the TPC cluster info
  esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
  esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
  esdTrack.SetTPCPointsF(track->GetTPCNclsF());
  if(!cuts->IsSelected(&esdTrack)) return kFALSE;
  
  return kTRUE;
}

//_________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::CheckAcceptance(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau){
  /// check if the decay products are in the good eta and pt range
  for (Int_t iProng = 0; iProng<nProng; iProng++){
    AliAODMCParticle* mcPartDaughter=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDau[iProng]));
    if(!mcPartDaughter) return kFALSE;
    Double_t eta = mcPartDaughter->Eta();
    Double_t pt = mcPartDaughter->Pt();
    if (TMath::Abs(eta) > fEtaAccCut || pt < fPtAccCut) return kFALSE;
  }
  return kTRUE;
}
//_________________________________________________________________
Int_t AliAnalysisTaskCombinHF::GetPoolIndex(Double_t zvert, Double_t mult){
  /// check in which of the pools the current event falls
  if(!fzVertPoolLims || !fMultPoolLims) return 0;
  Int_t theBinZ=TMath::BinarySearch(fNzVertPoolsLimSize,fzVertPoolLims,zvert);
  if(theBinZ<0 || theBinZ>=fNzVertPoolsLimSize) return -1;
  Int_t theBinM=TMath::BinarySearch(fNMultPoolsLimSize,fMultPoolLims,mult);
  if(theBinM<0 || theBinM>=fNMultPoolsLimSize) return -1;
  return fNMultPools*theBinZ+theBinM;
}
//_________________________________________________________________
void AliAnalysisTaskCombinHF::ResetPool(Int_t poolIndex){
  /// delete the contets of the pool
  if(poolIndex<0 || poolIndex>=fNOfPools) return;
  delete fEventBuffer[poolIndex];
  fEventBuffer[poolIndex]=new TTree(Form("EventBuffer_%d",poolIndex), "Temporary buffer for event mixing");
  fEventBuffer[poolIndex]->Branch("zVertex", &fVtxZ);
  fEventBuffer[poolIndex]->Branch("multiplicity", &fMultiplicity);
  fEventBuffer[poolIndex]->Branch("eventInfo", "TObjString",&fEventInfo);
  fEventBuffer[poolIndex]->Branch("karray", "TObjArray", &fKaonTracks);
  fEventBuffer[poolIndex]->Branch("parray", "TObjArray", &fPionTracks);
  return;
}
//_________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::CanBeMixed(Double_t zv1, Double_t zv2, Double_t mult1, Double_t mult2){
  /// check mixing
  if(TMath::Abs(zv2-zv1)>fMaxzVertDistForMix) return kFALSE;
  if(TMath::Abs(mult2-mult1)>fMaxMultDiffForMix) return kFALSE;
  return kTRUE;
}
//_________________________________________________________________
void AliAnalysisTaskCombinHF::DoMixingWithCuts(){
  /// perform mixed event analysis

  if(fDoEventMixing==0) return;
  Int_t nEvents=fEventBuffer[0]->GetEntries();
  if(fDebug > 1) printf("AnalysisTaskCombinHF::DoMixingWithCuts Start Event Mixing of %d events\n",nEvents);

  TObjArray* karray=0x0;
  TObjArray* parray=0x0;
  Double_t zVertex,mult;
  TObjString* eventInfo=0x0;
  fEventBuffer[0]->SetBranchAddress("karray", &karray);
  fEventBuffer[0]->SetBranchAddress("parray", &parray);
  fEventBuffer[0]->SetBranchAddress("eventInfo",&eventInfo);
  fEventBuffer[0]->SetBranchAddress("zVertex", &zVertex);
  fEventBuffer[0]->SetBranchAddress("multiplicity", &mult);
  Double_t d02[2]={0.,0.};
  Double_t d03[3]={0.,0.,0.};
  AliAODRecoDecay* tmpRD2 = new AliAODRecoDecay(0x0,2,0,d02);
  AliAODRecoDecay* tmpRD3 = new AliAODRecoDecay(0x0,3,1,d03);
  UInt_t pdg0[2]={321,211};
  //  UInt_t pdgp[3]={321,211,211};
  Double_t px[3],py[3],pz[3];
  Int_t evId1,esdId1,nk1,np1;
  Int_t evId2,esdId2,nk2,np2;

  for(Int_t iEv1=0; iEv1<nEvents; iEv1++){
    fEventBuffer[0]->GetEvent(iEv1);
    TObjArray* karray1=(TObjArray*)karray->Clone();
    Double_t zVertex1=zVertex;
    Double_t mult1=mult;
    Int_t nKaons=karray1->GetEntries();
    Int_t nPionsForCheck=parray->GetEntries();
    sscanf((eventInfo->String()).Data(),"Ev%d_esd%d_Pi%d_K%d",&evId1,&esdId1,&np1,&nk1);
    if(nk1!=nKaons || np1!=nPionsForCheck){ 
      printf("AnalysisTaskCombinHF::DoMixingWithCuts ERROR: read event does not match to the stored one\n");
      delete karray1;
      continue;
    }
    for(Int_t iEv2=0; iEv2<fNumberOfEventsForMixing; iEv2++){
      Int_t iToMix=iEv1+iEv2+1;
      if(iEv1>=(nEvents-fNumberOfEventsForMixing)) iToMix=iEv1-iEv2-1;
      if(iToMix<0) continue;
      if(iToMix==iEv1) continue;
      if(iToMix<iEv1) continue;
      fEventBuffer[0]->GetEvent(iToMix);
      Double_t zVertex2=zVertex;
      Double_t mult2=mult;
      if(TMath::Abs(zVertex2-zVertex1)<0.0001 && TMath::Abs(mult2-mult1)<0.001){
	printf("AnalysisTaskCombinHF::DoMixingWithCuts ERROR: same event in mixing??? %d %d   %f %f  %f %f\n",iEv1,iEv2,zVertex1,zVertex2,mult1,mult2);
	continue;
      }
      TObjArray* parray2=(TObjArray*)parray->Clone();
      Int_t nPions=parray2->GetEntries();
      Int_t nKaonsForCheck=karray->GetEntries();
      sscanf((eventInfo->String()).Data(),"Ev%d_esd%d_Pi%d_K%d",&evId2,&esdId2,&np2,&nk2);
      if(nk2!=nKaonsForCheck || np2!=nPions){ 
	printf("AnalysisTaskCombinHF::DoMixingWithCuts ERROR: read event does not match to the stored one\n");
	delete parray2;
	continue;
      }
      if(evId2==evId1 && esdId2==esdId1){
 	printf("AnalysisTaskCombinHF::DoMixingWithCuts ERROR: same event in mixing??? %d %d   nK=%d %d  nPi=%d %d\n",evId1,evId2,nKaons,nKaonsForCheck,nPionsForCheck,nPions);
	delete parray2;
	continue;
      }
      if(CanBeMixed(zVertex1,zVertex2,mult1,mult2)){
	for(Int_t iTr1=0; iTr1<nKaons; iTr1++){
	  TLorentzVector* trK=(TLorentzVector*)karray1->At(iTr1);
	  Double_t chargeK=trK->T();
	  px[0] = trK->Px();
	  py[0] = trK->Py();
	  pz[0] = trK->Pz();
	  for(Int_t iTr2=0; iTr2<nPions; iTr2++){
	    TLorentzVector* trPi1=(TLorentzVector*)parray2->At(iTr2);
	    Double_t chargePi1=trPi1->T();
	    px[1] = trPi1->Px();
	    py[1] = trPi1->Py();
	    pz[1] = trPi1->Pz();
	    if(fMeson==kDzero && chargePi1*chargeK<0){
	      FillMEHistos(421,2,tmpRD2,px,py,pz,pdg0);
	    }	  
	    if(fMeson==kDzero && chargePi1*chargeK>0){
	      FillMEHistosLS(421,2,tmpRD2,px,py,pz,pdg0,(Int_t)chargePi1);
	    }
	  }
	}
      }
      delete parray2;
    }
    delete karray1;
  }
  delete tmpRD2;
  delete tmpRD3;
}
//_________________________________________________________________
void AliAnalysisTaskCombinHF::DoMixingWithPools(Int_t poolIndex){
  /// perform mixed event analysis

  if(fDoEventMixing==0) return;
  if(poolIndex<0 || poolIndex>fNzVertPools*fNMultPools) return;

  Int_t nEvents=fEventBuffer[poolIndex]->GetEntries();
  if(fDebug > 1) printf("AliAnalysisTaskCombinHF::DoMixingWithPools Start Event Mixing of %d events\n",nEvents);
  TObjArray* karray=0x0;
  TObjArray* parray=0x0;
  Double_t zVertex,mult;
  TObjString* eventInfo=0x0;
  fEventBuffer[poolIndex]->SetBranchAddress("karray", &karray);
  fEventBuffer[poolIndex]->SetBranchAddress("parray", &parray);
  fEventBuffer[poolIndex]->SetBranchAddress("eventInfo",&eventInfo);
  fEventBuffer[poolIndex]->SetBranchAddress("zVertex", &zVertex);
  fEventBuffer[poolIndex]->SetBranchAddress("multiplicity", &mult);

  // dummy values of track impact parameter, needed to build an AliAODRecoDecay object
  Double_t d02[2]={0.,0.};
  Double_t d03[3]={0.,0.,0.};
  AliAODRecoDecay* tmpRD2 = new AliAODRecoDecay(0x0,2,0,d02);
  AliAODRecoDecay* tmpRD3 = new AliAODRecoDecay(0x0,3,1,d03);
  UInt_t pdg0[2]={321,211};
  UInt_t pdgp[3]={321,211,211};
  UInt_t pdgs[3]={321,211,321};
  Double_t px[3],py[3],pz[3];
  Int_t evId1,esdId1,nk1,np1;
  Int_t evId2,esdId2,nk2,np2;

  for(Int_t iEv1=0; iEv1<nEvents; iEv1++){
    fEventBuffer[poolIndex]->GetEvent(iEv1);
    TObjArray* karray1=(TObjArray*)karray->Clone();
    Double_t zVertex1=zVertex;
    Double_t mult1=mult;
    Int_t nKaons=karray1->GetEntries();
    Int_t nPionsForCheck=parray->GetEntries();
    sscanf((eventInfo->String()).Data(),"Ev%d_esd%d_Pi%d_K%d",&evId1,&esdId1,&np1,&nk1);
    if(nk1!=nKaons || np1!=nPionsForCheck){ 
      printf("AliAnalysisTaskCombinHF::DoMixingWithPools ERROR: read event does not match to the stored one\n");
      delete karray1;
      continue;
    }
    for(Int_t iEv2=0; iEv2<nEvents; iEv2++){
      if(iEv2==iEv1) continue;
      fEventBuffer[poolIndex]->GetEvent(iEv2);
      Double_t zVertex2=zVertex;
      Double_t mult2=mult;
      if(TMath::Abs(zVertex2-zVertex1)<0.0001 && TMath::Abs(mult2-mult1)<0.001){
	printf("AliAnalysisTaskCombinHF::DoMixingWithPools ERROR: same event in mixing??? %d %d   %f %f  %f %f\n",iEv1,iEv2,zVertex1,zVertex2,mult1,mult2);
	continue;
      }
      TObjArray* parray2=(TObjArray*)parray->Clone();
      Int_t nPions=parray2->GetEntries();
      Int_t nKaonsForCheck=karray->GetEntries();
      sscanf((eventInfo->String()).Data(),"Ev%d_esd%d_Pi%d_K%d",&evId2,&esdId2,&np2,&nk2);
      if(nk2!=nKaonsForCheck || np2!=nPions){ 
	printf("AliAnalysisTaskCombinHF::DoMixingWithPools ERROR: read event does not match to the stored one\n");
	delete parray2;
	continue;
      }
      if(evId2==evId1 && esdId2==esdId1){
 	printf("AliAnalysisTaskCombinHF::DoMixingWithPools ERROR: same event in mixing??? %d %d   nK=%d %d  nPi=%d %d\n",evId1,evId2,nKaons,nKaonsForCheck,nPionsForCheck,nPions);
	delete parray2;
	continue;
      }     
      TObjArray* parray3=0x0;
      Int_t nPions3=0;
      if(fMeson!=kDzero && fMeson!=kDs){
	Int_t iEv3=iEv2+1;
	if(iEv3==iEv1) iEv3=iEv2+2;
	if(iEv3>=nEvents) iEv3=iEv2-3;
	if(nEvents==2) iEv3=iEv1;
	if(iEv3<0) iEv3=iEv2-1;
	fEventBuffer[poolIndex]->GetEvent(iEv3);
	parray3=(TObjArray*)parray->Clone();
	nPions3=parray3->GetEntries();
      }
      for(Int_t iTr1=0; iTr1<nKaons; iTr1++){
	TLorentzVector* trK=(TLorentzVector*)karray1->At(iTr1);
	Double_t chargeK=trK->T();
	px[0] = trK->Px();
	py[0] = trK->Py();
	pz[0] = trK->Pz();
	for(Int_t iTr2=0; iTr2<nPions; iTr2++){
	  TLorentzVector* trPi1=(TLorentzVector*)parray2->At(iTr2);
	  Double_t chargePi1=trPi1->T();
	  px[1] = trPi1->Px();
	  py[1] = trPi1->Py();
	  pz[1] = trPi1->Pz();
	  if(chargePi1*chargeK<0){
	    if(fMeson==kDzero){
	      FillMEHistos(421,2,tmpRD2,px,py,pz,pdg0);
      }else if(fMeson==kDs) {
        for(Int_t iTr3=iTr1+1; iTr3<nKaons; iTr3++){
          TLorentzVector* trK2=(TLorentzVector*)karray1->At(iTr3);
          Double_t chargeK2=trK2->T();
          px[2] = trK2->Px();
          py[2] = trK2->Py();
          pz[2] = trK2->Pz();
          Double_t massKK=ComputeInvMassKK(trK,trK2);
          Double_t deltaMass=massKK-TDatabasePDG::Instance()->GetParticle(333)->Mass();
          Double_t cos1=CosPiKPhiRFrame(trK,trK2,trPi1);
          Double_t kincutPiKPhi=TMath::Abs(cos1*cos1*cos1);
          Double_t cosPiDsLabFrame=CosPiDsLabFrame(trK,trK2,trPi1);
          if(chargeK2*chargeK<0 && TMath::Abs(deltaMass)<fPhiMassCut && kincutPiKPhi>fCutCos3PiKPhiRFrame && cosPiDsLabFrame<fCutCosPiDsLabFrame){
            FillMEHistos(431,3,tmpRD3,px,py,pz,pdgs);
          }
        }
      } else{
	      if(parray3){
		for(Int_t iTr3=iTr2+1; iTr3<nPions3; iTr3++){
		  TLorentzVector* trPi2=(TLorentzVector*)parray3->At(iTr3);
		  Double_t chargePi2=trPi2->T();
		  px[2] = trPi2->Px();
		  py[2] = trPi2->Py();
		  pz[2] = trPi2->Pz();
		  if(chargePi2*chargeK<0){
		    if(fMeson==kDplus) FillMEHistos(411,3,tmpRD3,px,py,pz,pdgp);
		  }
		}
	      }
	    }
	  }else if(chargePi1*chargeK>0){
	    if(fMeson==kDzero) FillMEHistosLS(421,2,tmpRD2,px,py,pz,pdg0,(Int_t)chargePi1);
	  }
	}
      }
      delete parray3;
      delete parray2;
    }
    delete karray1;
  }
  delete tmpRD2;
  delete tmpRD3;
}
//_________________________________________________________________
void AliAnalysisTaskCombinHF::FinishTaskOutput()
{
  /// perform mixed event analysis
  if(fDoEventMixing==0) return;
  printf("AliAnalysisTaskCombinHF: FinishTaskOutput\n");

  if(fDoEventMixing==1){
    for(Int_t i=0; i<fNOfPools; i++){
      Int_t nEvents=fEventBuffer[i]->GetEntries();
      if(nEvents>1) DoMixingWithPools(i);	  
    }
  }else if(fDoEventMixing==2){
    DoMixingWithCuts();
  }
}
//_________________________________________________________________
Double_t AliAnalysisTaskCombinHF::ComputeInvMassKK(AliAODTrack* tr1, AliAODTrack* tr2) const{
  /// inv mass of KK
  Double_t massK=TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t p1=tr1->P();
  Double_t p2=tr2->P();
  Double_t pxtot=tr1->Px()+tr2->Px();
  Double_t pytot=tr1->Py()+tr2->Py();
  Double_t pztot=tr1->Pz()+tr2->Pz();
  Double_t e1=TMath::Sqrt(massK*massK+p1*p1);
  Double_t e2=TMath::Sqrt(massK*massK+p2*p2);
  Double_t etot=e1+e2;
  Double_t m2=etot*etot-(pxtot*pxtot+pytot*pytot+pztot*pztot);
  return TMath::Sqrt(m2);
}
//_________________________________________________________________
Double_t AliAnalysisTaskCombinHF::ComputeInvMassKK(TLorentzVector* tr1, TLorentzVector* tr2) const{
  /// inv mass of KK
  Double_t massK=TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t p1=tr1->P();
  Double_t p2=tr2->P();
  Double_t pxtot=tr1->Px()+tr2->Px();
  Double_t pytot=tr1->Py()+tr2->Py();
  Double_t pztot=tr1->Pz()+tr2->Pz();
  Double_t e1=TMath::Sqrt(massK*massK+p1*p1);
  Double_t e2=TMath::Sqrt(massK*massK+p2*p2);
  Double_t etot=e1+e2;
  Double_t m2=etot*etot-(pxtot*pxtot+pytot*pytot+pztot*pztot);
  return TMath::Sqrt(m2);
}

//----------------------------------------------------------------------
Double_t AliAnalysisTaskCombinHF::CosPiKPhiRFrame(TLorentzVector* dauK1, TLorentzVector* dauK2, TLorentzVector* daupi) const {
  /// computes cosine of angle between pi and K in the phi rest frame
  
  Double_t massK=TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t masspi=TDatabasePDG::Instance()->GetParticle(211)->Mass();

  Double_t eK1 = TMath::Sqrt(massK*massK+dauK1->P()*dauK1->P());
  Double_t eK2 = TMath::Sqrt(massK*massK+dauK2->P()*dauK2->P());
  Double_t epi = TMath::Sqrt(masspi*masspi+daupi->P()*daupi->P());

  Double_t ePhi=eK1+eK2;
  Double_t pxPhi=dauK1->Px()+dauK2->Px();
  Double_t pyPhi=dauK1->Py()+dauK2->Py();
  Double_t pzPhi=dauK1->Pz()+dauK2->Pz();
  Double_t bxPhi=pxPhi/ePhi;
  Double_t byPhi=pyPhi/ePhi;
  Double_t bzPhi=pzPhi/ePhi;
  
  TVector3 vecK1Phiframe;
  TLorentzVector vecK1(dauK1->Px(),dauK1->Py(),dauK1->Pz(),eK1);
  vecK1.Boost(-bxPhi,-byPhi,-bzPhi);
  vecK1.Boost(vecK1Phiframe);
  vecK1Phiframe=vecK1.BoostVector();
  
  TVector3 vecPiPhiframe;
  TLorentzVector vecPi(daupi->Px(),daupi->Py(),daupi->Pz(),epi);
  vecPi.Boost(-bxPhi,-byPhi,-bzPhi);
  vecPi.Boost(vecPiPhiframe);
  vecPiPhiframe=vecPi.BoostVector();
  
  Double_t innera=vecPiPhiframe.Dot(vecK1Phiframe);
  Double_t norm1a=TMath::Sqrt(vecPiPhiframe.Dot(vecPiPhiframe));
  Double_t norm2a=TMath::Sqrt(vecK1Phiframe.Dot(vecK1Phiframe));
  Double_t cosK1PhiFrame=innera/(norm1a*norm2a);
  
  return cosK1PhiFrame;
}

//----------------------------------------------------------------------
Double_t AliAnalysisTaskCombinHF::CosPiDsLabFrame(TLorentzVector* dauK1, TLorentzVector* dauK2, TLorentzVector* daupi) const {
  /// computes cosine of angle between pi and Ds in the Ds rest frame
  
  Double_t massDs=TDatabasePDG::Instance()->GetParticle(431)->Mass();
  Double_t masspi=TDatabasePDG::Instance()->GetParticle(211)->Mass();

  Double_t pxD=dauK1->Px()+dauK2->Px()+daupi->Px();
  Double_t pyD=dauK1->Px()+dauK2->Px()+daupi->Px();
  Double_t pzD=dauK1->Px()+dauK2->Px()+daupi->Px();
  Double_t pD=TMath::Sqrt(pxD*pxD+pyD*pyD+pzD*pzD);
  Double_t eD=TMath::Sqrt(massDs*massDs+pD*pD);
  Double_t epi = TMath::Sqrt(masspi*masspi+daupi->P()*daupi->P());
  Double_t bxD=pxD/eD;
  Double_t byD=pyD/eD;
  Double_t bzD=pzD/eD;
  
  TVector3 piDsframe;
  TLorentzVector vecPi(daupi->Px(),daupi->Py(),daupi->Pz(),epi);
  vecPi.Boost(-bxD,-byD,-bzD);
  vecPi.Boost(piDsframe);
  piDsframe=vecPi.BoostVector();
  
  TVector3 vecDs(pxD,pyD,pzD);
  
  Double_t inner=vecDs.Dot(piDsframe);
  Double_t norm1=TMath::Sqrt(vecDs.Dot(vecDs));
  Double_t norm2=TMath::Sqrt(piDsframe.Dot(piDsframe));
  Double_t cosPiDsFrame=inner/(norm1*norm2);
  
  return cosPiDsFrame;
}

//_________________________________________________________________
void AliAnalysisTaskCombinHF::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  //
  if(fDebug > 1) printf("AliAnalysisTaskCombinHF: Terminate() \n");
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  if(fHistNEvents){
    printf("Number of analyzed events = %d\n",(Int_t)fHistNEvents->GetBinContent(2));
  }else{
    printf("ERROR: fHistNEvents not available\n");
    return;
  }
  return;
}

