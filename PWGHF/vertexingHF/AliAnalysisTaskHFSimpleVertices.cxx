#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliVertexerTracks.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAnalysisTaskHFSimpleVertices.h"
#include "AliMultSelection.h"
#include <TH1F.h>
#include <TSystem.h>
#include <TChain.h>
#include <TDatabasePDG.h>


/**************************************************************************
 * Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Implementation of class AliAnalysisTaskHFSimpleVertices
// AliAnalysisTaskSE to extract D meson candidates from ESDs
//          
//*************************************************************************

ClassImp(AliAnalysisTaskHFSimpleVertices)
//______________________________________________________________________________
AliAnalysisTaskHFSimpleVertices::AliAnalysisTaskHFSimpleVertices() : 
  AliAnalysisTaskSE("HFSimpleVertices"), 
  fOutput{nullptr},
  fHistNEvents{nullptr},
  fHistPtAllTracks{nullptr},
  fHistPtSelTracks{nullptr},
  fHistTglAllTracks{nullptr},
  fHistTglSelTracks{nullptr},
  fHistImpParAllTracks{nullptr},
  fHistImpParSelTracks{nullptr},
  fHistITSmapAllTracks{nullptr},
  fHistITSmapSelTracks{nullptr},
  fHistPrimVertX{nullptr},
  fHistPrimVertY{nullptr},
  fHistPrimVertZ{nullptr},
  fHist2ProngVertX{nullptr},
  fHist2ProngVertY{nullptr},
  fHist2ProngVertZ{nullptr},
  fHistDplusVertX{nullptr},
  fHistDplusVertY{nullptr},
  fHistDplusVertZ{nullptr},
  fHistInvMassD0{nullptr},
  fHistPtD0{nullptr},
  fHistPtD0Dau0{nullptr},
  fHistPtD0Dau1{nullptr},
  fHistImpParD0Dau0{nullptr},
  fHistImpParD0Dau1{nullptr},
  fHistd0Timesd0{nullptr},
  fHistDecLenD0{nullptr},
  fHistDecLenXYD0{nullptr},
  fHistInvMassDplus{nullptr},
  fUsePhysSel(kTRUE),
  fTriggerMask(AliVEvent::kAny),
  fSelectOnCentrality(kFALSE),
  fMinCentrality(-1.),
  fMaxCentrality(110.),
  fCentrEstimator("V0M"),
  fDo3Prong(kFALSE),
  fMaxDecVertRadius2(8),
  fMassDzero(0.),
  fMassDplus(0.),
  fMassDs(0.),
  fMassLambdaC(0.),
  fTrackCuts2pr{nullptr},
  fTrackCuts3pr{nullptr},
  fMaxTracksToProcess(9999999),
  fNPtBins(25),
  fMinPtDzero(0.),
  fMaxPtDzero(9999.),
  fSelectD0(1),
  fSelectD0bar(1)
{
  //
  

  InitDefault();
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskHFSimpleVertices::~AliAnalysisTaskHFSimpleVertices(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;

  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistPtAllTracks;
    delete fHistPtSelTracks;
    delete fHistTglAllTracks;
    delete fHistTglSelTracks;
    delete fHistImpParAllTracks;
    delete fHistImpParSelTracks;
    delete fHistITSmapAllTracks;
    delete fHistITSmapSelTracks;
    delete fHistPrimVertX;
    delete fHistPrimVertY;
    delete fHistPrimVertZ;
    delete fHist2ProngVertX;
    delete fHist2ProngVertY;
    delete fHist2ProngVertZ;
    delete fHistDplusVertX;
    delete fHistDplusVertY;
    delete fHistDplusVertZ;
    delete fHistInvMassD0;
    delete fHistPtD0;
    delete fHistPtD0Dau0;
    delete fHistPtD0Dau1;
    delete fHistImpParD0Dau0;
    delete fHistImpParD0Dau1;
    delete fHistd0Timesd0;
    delete fHistDecLenD0;
    delete fHistDecLenXYD0;
    delete fHistInvMassDplus;
  }
  delete fOutput;
  delete fTrackCuts2pr;
  delete fTrackCuts3pr;
}
 
//___________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::InitDefault(){
  /// initialization with default values
  
  fMassDzero = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  fMassDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  fMassDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  fMassLambdaC = TDatabasePDG::Instance()->GetParticle(4122)->Mass();

  fTrackCuts2pr = new AliESDtrackCuts("AliESDtrackCuts", "default");
  fTrackCuts2pr->SetPtRange(0., 1.e10);
  // fTrackCuts->SetEtaRange(-0.8, +0.8);
  fTrackCuts2pr->SetMinNClustersTPC(50);
  fTrackCuts2pr->SetRequireITSRefit(kTRUE);
  fTrackCuts2pr->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					  AliESDtrackCuts::kAny);
  // fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  // fTrackCuts->SetMaxDCAToVertexZ(3.2);
  // fTrackCuts->SetMaxDCAToVertexXY(2.4);
  // fTrackCuts->SetDCAToVertex2D(kTRUE);
  
  fTrackCuts3pr = new AliESDtrackCuts("AliESDtrackCuts", "default3p");
  fTrackCuts3pr->SetPtRange(0., 1.e10);
  // fTrackCuts->SetEtaRange(-0.8, +0.8);
  fTrackCuts3pr->SetMinNClustersTPC(50);
  fTrackCuts3pr->SetRequireITSRefit(kTRUE);
  fTrackCuts3pr->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					  AliESDtrackCuts::kAny);
  
  fNPtBins=25;
  Double_t defaultPtBins[26] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0, 16.0, 20.0, 24.0, 36.0, 50.0, 100.0};
  for(Int_t ib=0; ib<fNPtBins+1; ib++) fPtBinLims[ib]=defaultPtBins[ib];
  
  Double_t defaultD0Cuts[25][kNCutVarsDzero] =
    {{0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   /* pt<0.5*/
     {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   /* 0.5<pt<1*/
     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  /* 1<pt<1.5 */
     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  /* 1.5<pt<2 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  /* 2<pt<2.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  /* 2.5<pt<3 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  /* 3<pt<3.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  /* 3.5<pt<4 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 4<pt<4.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 4.5<pt<5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 5<pt<5.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 5.5<pt<6 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 6<pt<6.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 6.5<pt<7 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   /* 7<pt<7.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   /* 7.5<pt<8 */
     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   /* 8<pt<9 */
     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   /* 9<pt<10 */
     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   /* 10<pt<12 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 10000. * 1E-8, 0.85, 0., 0.},   /* 12<pt<16 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 16<pt<20 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 20<pt<24 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 24<pt<36 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 36<pt<50 */
     {0.400, 300. * 1E-4, 1.0, 0.6, 0.6, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.80, 0., 0.}}; /* pt>50 */
  for(Int_t ib=0; ib<fNPtBins; ib++){
   for(Int_t jc=0; jc<kNCutVarsDzero; jc++){
     fDzeroCuts[ib][jc]=defaultD0Cuts[ib][jc];
   }
  }
}

//___________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::InitFromJson(TString filename){
  /// read configuration from json file
  if (filename != "" && gSystem->Exec(Form("ls %s > /dev/null", filename.Data())) == 0) {
    printf("------Read configuration from JSON file------\n");
    Double_t ptmintrack2 = GetJsonFloat(filename.Data(), "ptmintrack_2prong");
    printf("Min pt track (2 prong)= %f\n", ptmintrack2);
    if(ptmintrack2>0) fTrackCuts2pr->SetPtRange(ptmintrack2, 1.e10);
    Double_t ptmintrack3 = GetJsonFloat(filename.Data(), "ptmintrack_3prong");
    printf("Min pt track (3 prong)= %f\n", ptmintrack3);
    if(ptmintrack3>0) fTrackCuts3pr->SetPtRange(ptmintrack3, 1.e10);
    Int_t do3Prongs = GetJsonInteger(filename.Data(), "do3prong");
    printf("do3prong     = %d\n", do3Prongs);
    if(do3Prongs>0) fDo3Prong=kTRUE;
    Int_t selectD0 = GetJsonInteger(filename.Data(), "d_selectionFlagD0");
    printf("d_selectionFlagD0 = %d\n",selectD0);
    if(selectD0>=0) fSelectD0=selectD0;
    Int_t selectD0bar = GetJsonInteger(filename.Data(), "d_selectionFlagD0bar");
    printf("d_selectionFlagD0bar = %d\n",selectD0bar);
    if(selectD0>=0) fSelectD0bar=selectD0bar;
    Int_t minncluTPC = GetJsonInteger(filename.Data(), "d_tpcnclsfound");
    if(minncluTPC>0) printf("minncluTPC   = %d\n", minncluTPC);
    fTrackCuts2pr->SetMinNClustersTPC(minncluTPC);
    fTrackCuts3pr->SetMinNClustersTPC(minncluTPC);
    Double_t dcatoprimxymin2 = GetJsonFloat(filename.Data(), "dcatoprimxymin_2prong");
    printf("dcatoprimxymin  (2 prong) = %f\n", dcatoprimxymin2);
    if(dcatoprimxymin2>0) fTrackCuts2pr->SetMinDCAToVertexXY(dcatoprimxymin2);
    Double_t dcatoprimxymin3 = GetJsonFloat(filename.Data(), "dcatoprimxymin_3prong");
    printf("dcatoprimxymin  (3 prong) = %f\n", dcatoprimxymin3);
    if(dcatoprimxymin3>0) fTrackCuts3pr->SetMinDCAToVertexXY(dcatoprimxymin3);
    Double_t etamax2 = GetJsonFloat(filename.Data(), "etamax_2prong");
    printf("Max eta  (2 prong) = %f\n", etamax2);
    if(etamax2>0) fTrackCuts2pr->SetEtaRange(-etamax2, +etamax2);
    Double_t etamax3 = GetJsonFloat(filename.Data(), "etamax_3prong");
    printf("Max eta  (3 prong) = %f\n", etamax3);
    if(etamax3>0) fTrackCuts3pr->SetEtaRange(-etamax3, +etamax3);
    
    Double_t d_maxr = GetJsonFloat(filename.Data(), "d_maxr");
    if(d_maxr>0) fMaxDecVertRadius2=d_maxr*d_maxr;
    Double_t ptMinCand = GetJsonFloat(filename.Data(), "d_pTCandMin");
    printf("Min pt Dzero cand = %f\n", ptMinCand);
    if(ptMinCand>=0.) fMinPtDzero=ptMinCand;
    Double_t ptMaxCand = GetJsonFloat(filename.Data(), "d_pTCandMax");
    printf("Max pt Dzero cand = %f\n", ptMaxCand);
    if(ptMaxCand>=0. && ptMaxCand>fMinPtDzero) fMaxPtDzero=ptMaxCand;
    printf("---------------------------------------------\n");
  }else{
    AliError(Form("Json configuration file %s not found\n",filename.Data()));
  }
}

//___________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::UserCreateOutputObjects() {
  // create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",15,-0.5,14.5);
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"PhysSel");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"InCentralityClass");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Good vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Pass zSPD-zTrk vert sel");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"|zvert|<10");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Pileup cut");
  fOutput->Add(fHistNEvents);

  // single track histos
  fHistPtAllTracks = new TH1F("hPtAllTracks", " All tracks ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtSelTracks = new TH1F("hPtSelTracks", " Selected tracks ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistTglAllTracks = new TH1F("hTglAllTracks", " All tracks ; tg#lambda", 100, -5., 5.);
  fHistTglSelTracks = new TH1F("hTglSelTracks", " Selected tracks ; tg#lambda", 100, -5., 5.);
  fHistImpParAllTracks = new TH1F("hImpParAllTracks", " All tracks ; d_{0}^{xy} (cm)", 100, -1, 1.);
  fHistImpParSelTracks = new TH1F("hImpParSelTracks", " Selected tracks ; d_{0}^{xy} (cm)", 100, -1, 1.);
  fHistITSmapAllTracks = new TH1F("hITSmapAllTracks", " All tracks ; ITS cluster map", 64, -0.5, 63.5);
  fHistITSmapSelTracks = new TH1F("hITSmapSelTracks", " Selected tracks ; ITS cluster map", 64, -0.5, 63.5);
  fOutput->Add(fHistPtAllTracks);
  fOutput->Add(fHistPtSelTracks);
  fOutput->Add(fHistTglAllTracks);
  fOutput->Add(fHistTglSelTracks);
  fOutput->Add(fHistImpParAllTracks);
  fOutput->Add(fHistImpParSelTracks);
  fOutput->Add(fHistITSmapAllTracks);
  fOutput->Add(fHistITSmapSelTracks);

  // vertex histos
  fHistPrimVertX = new TH1F("hPrimVertX"," Primary Vertex ; x (cm)",100, -0.5, 0.5);
  fHistPrimVertY = new TH1F("hPrimVertY"," Primary Vertex ; y (cm)",100, -0.5, 0.5);
  fHistPrimVertZ = new TH1F("hPrimVertZ"," Primary Vertex ; z (cm)",100, -20.0, 20.0);
  fHist2ProngVertX = new TH1F("h2ProngVertX"," Secondary Vertex ; x (cm)",1000, -2., 2.);
  fHist2ProngVertY = new TH1F("h2ProngVertY"," Secondary Vertex ; y (cm)",1000, -2., 2.);
  fHist2ProngVertZ = new TH1F("h2ProngVertZ"," Secondary Vertex ; z (cm)",1000, -20.0, 20.0);
  fHistDplusVertX = new TH1F("hDplusVertX"," Secondary Vertex ; x (cm)",100, -1., 1.);
  fHistDplusVertY = new TH1F("hDplusVertY"," Secondary Vertex ; y (cm)",100, -1., 1.);
  fHistDplusVertZ = new TH1F("hDplusVertZ"," Secondary Vertex ; z (cm)",100, -20.0, 20.0);
  fOutput->Add(fHistPrimVertX);
  fOutput->Add(fHistPrimVertY);
  fOutput->Add(fHistPrimVertZ);
  fOutput->Add(fHist2ProngVertX);
  fOutput->Add(fHist2ProngVertY);
  fOutput->Add(fHist2ProngVertZ);
  fOutput->Add(fHistDplusVertX);
  fOutput->Add(fHistDplusVertY);
  fOutput->Add(fHistDplusVertZ);
  
 
  // D meson candidate histos
  fHistInvMassD0 = new TH1F("hInvMassD0" , " ; M_{K#pi} (GeV/c^{2})",500, 0, 5.0);
  fHistPtD0  = new TH1F("hPtD0" , " ; D^{0} p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtD0Dau0 = new TH1F("hPtD0Dau0" , " D^{0} prong0 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtD0Dau1 = new TH1F("hPtD0Dau1" , " D^{0} prong1 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistImpParD0Dau0 = new TH1F("hImpParD0Dau0" , " D^{0} prong0 ; d_{0}^{xy} (cm)", 100, -1.0, 1.0);
  fHistImpParD0Dau1 = new TH1F("hImpParD0Dau1" , " D^{0} prong1 ; d_{0}^{xy} (cm)", 100, -1.0, 1.0);
  fHistd0Timesd0 = new TH1F("hd0Timesd0" , " d_{0}^{xy}x d_{0}^{xy} (cm^{2})", 500, -1.0, 1.0);
  fHistDecLenD0 = new TH1F("hDecLenD0" , " ; Decay Length (cm)",200, 0., 2.0);
  fHistDecLenXYD0 = new TH1F("hDecLenXYD0" , " ; Decay Length xy (cm)",200, 0., 2.0);
  fHistInvMassDplus = new TH1F("hInvMassDplus" , " ; M_{K#pi#pi} (GeV/c^{2})",500, 1.6, 2.1);
  fOutput->Add(fHistPtD0);
  fOutput->Add(fHistPtD0Dau0);
  fOutput->Add(fHistPtD0Dau1);
  fOutput->Add(fHistImpParD0Dau0);
  fOutput->Add(fHistImpParD0Dau1);
  fOutput->Add(fHistd0Timesd0);
  fOutput->Add(fHistInvMassD0);
  fOutput->Add(fHistDecLenD0);
  fOutput->Add(fHistDecLenXYD0);
  fOutput->Add(fHistInvMassDplus);
  
  PostData(1,fOutput);

}
//______________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::UserExec(Option_t *)
{
  //

  AliESDEvent *esd = (AliESDEvent*) (InputEvent());
  if(!esd) {
    printf("AliAnalysisTaskHFSimpleVertices::UserExec(): bad ESD\n");
    return;
  } 

  // AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  // AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
  
  // AliMCEvent* mcEvent = nullptr;

  // if(fReadMC){
  //   AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  //   if (!eventHandler) {
  //     Printf("ERROR: Could not retrieve MC event handler");
  //     return;
  //   }
  //   mcEvent = eventHandler->MCEvent();
  //   if (!mcEvent) {
  //     Printf("ERROR: Could not retrieve MC event");
  //     return;
  //   }
  // }


  fHistNEvents->Fill(0);
  if(fUsePhysSel){
    Bool_t isPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
    if(!isPhysSel) return;
  }
  fHistNEvents->Fill(1);

  if(fSelectOnCentrality){
    AliMultSelection* mulSel = (AliMultSelection*)esd->FindListObject("MultSelection");
    if(mulSel){
      Double_t centr=mulSel->GetMultiplicityPercentile(fCentrEstimator.Data());
      if(centr<fMinCentrality || centr>fMaxCentrality) return;
    }
  }
  fHistNEvents->Fill(2);
  
  AliESDVertex* primVtxTrk = (AliESDVertex*)esd->GetPrimaryVertex();
  AliESDVertex* primVtxSPD = (AliESDVertex*)esd->GetPrimaryVertexSPD();
  TString titTrc=primVtxTrk->GetTitle();
  if(titTrc.IsNull())return;
  if (primVtxTrk->IsFromVertexer3D() || primVtxTrk->IsFromVertexerZ()) return;
  if (primVtxTrk->GetNContributors() < 2) return;
  if (primVtxSPD->GetNContributors()<1) return;
  fHistNEvents->Fill(3);

  double covTrc[6],covSPD[6];
  primVtxTrk->GetCovarianceMatrix(covTrc);
  primVtxSPD->GetCovarianceMatrix(covSPD);
  double dz = primVtxTrk->GetZ()-primVtxSPD->GetZ();
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
  if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return; // bad vertexing
  fHistNEvents->Fill(4);

  Float_t zvert=primVtxTrk->GetZ();
  if(TMath::Abs(zvert)>10) return;
  fHistNEvents->Fill(5);
  
  fHistPrimVertX->Fill(primVtxTrk->GetX());
  fHistPrimVertY->Fill(primVtxTrk->GetY());
  fHistPrimVertZ->Fill(primVtxTrk->GetZ());

  AliAODVertex *vertexAODp = ConvertToAODVertex(primVtxTrk);
  Double_t bzkG = (Double_t)esd->GetMagneticField();
  Int_t totTracks = TMath::Min(fMaxTracksToProcess, esd->GetNumberOfTracks());
  Double_t d0track[2],covd0track[2];
  
  // Apply single track cuts and flag them
  UChar_t* status = new UChar_t[totTracks];
  for (Int_t iTrack = 0; iTrack < totTracks; iTrack++) {
    status[iTrack] = 0;
    AliESDtrack* track = esd->GetTrack(iTrack);
    track->PropagateToDCA(primVtxTrk, bzkG, 100., d0track, covd0track);
    fHistPtAllTracks->Fill(track->Pt());
    fHistTglAllTracks->Fill(track->GetTgl());
    fHistImpParAllTracks->Fill(d0track[0]);
    fHistITSmapAllTracks->Fill(track->GetITSClusterMap());
    status[iTrack] = SingleTrkCuts(track,primVtxTrk,bzkG);
  }

  Double_t d03[3] = {0., 0., 0.};
  AliAODRecoDecay* rd4massCalc3 = new AliAODRecoDecay(0x0, 3, 1, d03);
  TObjArray* twoTrackArray = new TObjArray(2);
  TObjArray* threeTrackArray = new TObjArray(3);
  AliVertexerTracks* vt = new AliVertexerTracks(bzkG);
  Double_t mom0[3], mom1[3], mom2[3];
  for (Int_t iPosTrack_0 = 0; iPosTrack_0 < totTracks; iPosTrack_0++) {
    AliESDtrack* track_p0 = esd->GetTrack(iPosTrack_0);
    track_p0->GetPxPyPz(mom0);
    if (status[iPosTrack_0] == 0) continue;
    track_p0->PropagateToDCA(primVtxTrk, bzkG, 100., d0track, covd0track);
    fHistPtSelTracks->Fill(track_p0->Pt());
    fHistTglSelTracks->Fill(track_p0->GetTgl());
    fHistImpParSelTracks->Fill(d0track[0]);
    fHistITSmapSelTracks->Fill(track_p0->GetITSClusterMap());
    if (track_p0->Charge() < 0) continue;
    for (Int_t iNegTrack_0 = 0; iNegTrack_0 < totTracks; iNegTrack_0++) {
      AliESDtrack* track_n0 = esd->GetTrack(iNegTrack_0);
      track_n0->GetPxPyPz(mom1);
      if (track_n0->Charge() > 0) continue;
      if (status[iNegTrack_0] == 0) continue;
      twoTrackArray->AddAt(track_p0, 0);
      twoTrackArray->AddAt(track_n0, 1);
      AliESDVertex* trkv = ReconstructSecondaryVertex(vt, twoTrackArray, primVtxTrk);
      if (trkv == 0x0) {
	twoTrackArray->Clear();
	continue;
      }
      fHist2ProngVertX->Fill(trkv->GetX());
      fHist2ProngVertY->Fill(trkv->GetY());
      fHist2ProngVertZ->Fill(trkv->GetZ());

      double deltax = trkv->GetX() - primVtxTrk->GetX();
      double deltay = trkv->GetY() - primVtxTrk->GetY();
      double deltaz = trkv->GetZ() - primVtxTrk->GetZ();
      double decaylength = TMath::Sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
      double decaylengthxy = TMath::Sqrt(deltax * deltax + deltay * deltay);

      AliAODVertex* vertexAOD = ConvertToAODVertex(trkv);
      delete trkv;
      AliAODRecoDecayHF2Prong* the2Prong = Make2Prong(twoTrackArray, vertexAOD, bzkG);
      the2Prong->SetOwnPrimaryVtx(vertexAODp);
      Int_t deroSel = 3;
      if(fSelectD0 + fSelectD0bar > 0){
	deroSel = DzeroSelectionCuts(the2Prong);
      }
      if(deroSel>0){
	Double_t m0 = the2Prong->InvMassD0();
	Double_t m0b = the2Prong->InvMassD0bar();
	Double_t ptD = the2Prong->Pt();
	Double_t ptDau0 = the2Prong->PtProng(0);
	Double_t ptDau1 = the2Prong->PtProng(1);
	Double_t ipDau0 = the2Prong->Getd0Prong(0);
	Double_t ipDau1 = the2Prong->Getd0Prong(1);
	Double_t d0xd0 = the2Prong->Prodd0d0();
	if (fSelectD0 == 0 || deroSel == 1 || deroSel == 3) fHistInvMassD0->Fill(m0);
	if (fSelectD0bar == 0 || deroSel == 2 || deroSel == 3) fHistInvMassD0->Fill(m0b);
	fHistPtD0->Fill(ptD);
	fHistPtD0Dau0->Fill(ptDau0);
	fHistPtD0Dau1->Fill(ptDau1);
	fHistImpParD0Dau0->Fill(ipDau0);
	fHistImpParD0Dau1->Fill(ipDau1);
	fHistd0Timesd0->Fill(d0xd0);
	fHistDecLenD0->Fill(decaylength);
	fHistDecLenXYD0->Fill(decaylengthxy);
      }
      delete the2Prong;
      delete vertexAOD;
      
      if (fDo3Prong) {
	if(status[iPosTrack_0]<=1) continue;
	if(status[iNegTrack_0]<=1) continue;
	for (Int_t iPosTrack_1 = iPosTrack_0 + 1; iPosTrack_1 < totTracks; iPosTrack_1++) {
	  AliESDtrack* track_p1 = esd->GetTrack(iPosTrack_1);
	  if (!track_p1) continue;
	  if (track_p1->Charge() < 0) continue;
	  track_p1->GetPxPyPz(mom2);
	  if (status[iPosTrack_1] <= 1) continue;
	  // order tracks according to charge: +-+
	  threeTrackArray->AddAt(track_p0, 0);
	  threeTrackArray->AddAt(track_n0, 1);
	  threeTrackArray->AddAt(track_p1, 2);
	  Int_t massSel = SelectInvMassAndPt3prong(threeTrackArray, rd4massCalc3);
	  if (massSel == 0) {
	    threeTrackArray->Clear();
	    continue;
	  }
	  AliESDVertex* trkv3 = ReconstructSecondaryVertex(vt, threeTrackArray, primVtxTrk);
	  if (trkv3 == 0x0) {
	    threeTrackArray->Clear();
	    continue;
	  }
	  AliAODVertex* vertexAOD3 = ConvertToAODVertex(trkv3);
	  AliAODRecoDecayHF3Prong* the3Prong = Make3Prong(threeTrackArray, vertexAOD3, bzkG);
	  //  the3Prong->SetOwnPrimaryVtx(vertexAODp);
	  if (massSel & (1 << kbitDplus)) {
	    Double_t mp = the3Prong->InvMassDplus();
	    fHistDplusVertX->Fill(trkv3->GetX());
	    fHistDplusVertY->Fill(trkv3->GetY());
	    fHistDplusVertZ->Fill(trkv3->GetZ());
	    fHistInvMassDplus->Fill(mp);
	  }
	  delete trkv3;
	  delete the3Prong;
	  delete vertexAOD3;
	  threeTrackArray->Clear();
	}
	for (Int_t iNegTrack_1 = iNegTrack_0 + 1; iNegTrack_1 < totTracks; iNegTrack_1++) {
	  AliESDtrack* track_n1 = esd->GetTrack(iNegTrack_1);
	  if (!track_n1) continue;
	  if (track_n1->Charge() > 0) continue;
	  if (status[iNegTrack_1] <= 1) continue;
	  track_n1->GetPxPyPz(mom2);
	  // order tracks according to charge: -+-
	  threeTrackArray->AddAt(track_n0, 0);
	  threeTrackArray->AddAt(track_p0, 1);
	  threeTrackArray->AddAt(track_n1, 2);
	  Int_t massSel = SelectInvMassAndPt3prong(threeTrackArray, rd4massCalc3);
	  if (massSel == 0) {
	    threeTrackArray->Clear();
	    continue;
	  }
	  AliESDVertex* trkv3 = ReconstructSecondaryVertex(vt, threeTrackArray, primVtxTrk);
	  if (trkv3 == 0x0) {
	    threeTrackArray->Clear();
	    continue;
	  }
	  AliAODVertex* vertexAOD3 = ConvertToAODVertex(trkv3);
	  AliAODRecoDecayHF3Prong* the3Prong = Make3Prong(threeTrackArray, vertexAOD3, bzkG);
	  //  the3Prong->SetOwnPrimaryVtx(vertexAODp);
	  if (massSel & (1 << kbitDplus)) {
	    Double_t mp = the3Prong->InvMassDplus();
	    fHistDplusVertX->Fill(trkv3->GetX());
	    fHistDplusVertY->Fill(trkv3->GetY());
	    fHistDplusVertZ->Fill(trkv3->GetZ());
	    fHistInvMassDplus->Fill(mp);
	  }
	  delete trkv3;
	  delete the3Prong;
	  delete vertexAOD3;
	  threeTrackArray->Clear();
	}
      }
      twoTrackArray->Clear();
    }
    //  delete vertexAODp;
  }
  delete[] status;
  delete vt;
  delete twoTrackArray;
  delete threeTrackArray;
  delete rd4massCalc3;
  delete vertexAODp;
  
  PostData(1,fOutput);
  
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskHFSimpleVertices::GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3], float bzkG)
{
  /// fast calculation (no covariance matrix treatment) of track momentum at secondary vertex

  Double_t alpha = tr->GetAlpha();
  Double_t sn = TMath::Sin(alpha), cs = TMath::Cos(alpha);
  Double_t x = tr->GetX(), y = tr->GetParameter()[0], snp = tr->GetParameter()[2];
  Double_t xv = secVert->GetX() * cs + secVert->GetY() * sn;
  Double_t yv = -secVert->GetX() * sn + secVert->GetY() * cs;
  x -= xv;
  y -= yv;
  Double_t crv = tr->GetC(bzkG);
  if (TMath::Abs(bzkG) < 0.000001)
    crv = 0.;
  double csp = TMath::Sqrt((1. - snp) * (1. + snp));

  Double_t tgfv = -(crv * x - snp) / (crv * y + csp);
  cs = 1. / TMath::Sqrt(1 + tgfv * tgfv);
  sn = cs < 1. ? tgfv * cs : 0.;

  x = xv * cs + yv * sn;
  Double_t alpNew = alpha + TMath::ASin(sn);
  Double_t ca = TMath::Cos(alpNew - alpha), sa = TMath::Sin(alpNew - alpha);
  Double_t p2 = tr->GetSnp();
  Double_t xNew = tr->GetX() * ca + tr->GetY() * sa;
  Double_t p2New = p2 * ca - TMath::Sqrt((1. - p2) * (1. + p2)) * sa;
  momentum[0] = tr->GetSigned1Pt();
  momentum[1] = p2New * (x - xNew) * tr->GetC(bzkG);
  momentum[2] = tr->GetTgl();
  Bool_t retCode = tr->Local2GlobalMomentum(momentum, alpNew);
  return retCode;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::SingleTrkCuts(AliESDtrack* trk, AliESDVertex* primVert, Double_t bzkG)
{
  if (!trk->PropagateToDCA(primVert, bzkG, kVeryBig)) return kFALSE;
  trk->RelateToVertex(primVert, bzkG, kVeryBig);
  Int_t retCode=0;
  if(fTrackCuts2pr->AcceptTrack(trk)) retCode+=1;
  if(fTrackCuts3pr->AcceptTrack(trk)) retCode+=2;
  return retCode;
}
//______________________________________________________________________________
AliESDVertex* AliAnalysisTaskHFSimpleVertices::ReconstructSecondaryVertex(AliVertexerTracks* vt, TObjArray* trkArray, AliESDVertex* primvtx)
{

  vt->SetVtxStart(primvtx);

  AliESDVertex* trkv = (AliESDVertex*)vt->VertexForSelectedESDTracks(trkArray);
  if (trkv->GetNContributors() != trkArray->GetEntriesFast())
    return 0x0;
  Double_t vertRadius2 = trkv->GetX() * trkv->GetX() + trkv->GetY() * trkv->GetY();
  if(vertRadius2>fMaxDecVertRadius2) return 0x0;
  return trkv;
}
//______________________________________________________________________________
AliAODVertex* AliAnalysisTaskHFSimpleVertices::ConvertToAODVertex(AliESDVertex* trkv)
{
  Double_t pos[3], cov[6], chi2perNDF;
  trkv->GetXYZ(pos);       // position
  trkv->GetCovMatrix(cov); // covariance matrix
  chi2perNDF = trkv->GetChi2toNDF();
  //  double dispersion = trkv->GetDispersion();
  AliAODVertex* vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, 2);
  return vertexAOD;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::DzeroSelectionCuts(AliAODRecoDecayHF2Prong* cand)
{
  bool isD0 = true;
  bool isD0bar = true;
  Double_t ptCand = cand->Pt();
  if (ptCand < fMinPtDzero || ptCand > fMaxPtDzero) return 0;
  Int_t jPtBin = GetPtBin(ptCand);
  if (jPtBin==-1) return 0;
  if (cand->Prodd0d0() > fDzeroCuts[jPtBin][7]) return 0;
  if (cand->CosPointingAngle() < fDzeroCuts[jPtBin][8]) return 0;
  if (cand->CosPointingAngleXY() < fDzeroCuts[jPtBin][9]) return 0;
  if (cand->NormalizedDecayLengthXY() < fDzeroCuts[jPtBin][10]) return 0;
  Double_t decayLengthCut = TMath::Min((cand->P() * 0.0066) + 0.01, 0.06);
  if (TMath::Abs(cand->Normalizedd0Prong(0)) < 0.5 || TMath::Abs(cand->Normalizedd0Prong(1)) < 0.5) return 0;
  if (cand->DecayLength() * cand->DecayLength() < decayLengthCut * decayLengthCut) return 0;
  // if (cand->NormalizedDecayLength() * cand->NormalizedDecayLength() < 1.0) return 0;
  if (TMath::Abs(cand->InvMassD0()-fMassDzero) > fDzeroCuts[jPtBin][0] ) isD0=false;
  if (TMath::Abs(cand->InvMassD0bar()-fMassDzero) > fDzeroCuts[jPtBin][0] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;

  if (cand->Pt2Prong(0) < fDzeroCuts[jPtBin][4]*fDzeroCuts[jPtBin][4] || cand->Pt2Prong(1) < fDzeroCuts[jPtBin][3]*fDzeroCuts[jPtBin][3] ) isD0=false;
  if (cand->Pt2Prong(0) < fDzeroCuts[jPtBin][3]*fDzeroCuts[jPtBin][3] || cand->Pt2Prong(1) < fDzeroCuts[jPtBin][4]*fDzeroCuts[jPtBin][4] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;

  if (TMath::Abs(cand->Getd0Prong(0)) > fDzeroCuts[jPtBin][6] || TMath::Abs(cand->Getd0Prong(1)) > fDzeroCuts[jPtBin][5] ) isD0=false;
  if (TMath::Abs(cand->Getd0Prong(0)) > fDzeroCuts[jPtBin][5] || TMath::Abs(cand->Getd0Prong(1)) > fDzeroCuts[jPtBin][6] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;

  Double_t cosThetaStarD0,cosThetaStarD0bar;
  cand->CosThetaStarD0(cosThetaStarD0,cosThetaStarD0bar);
  if (TMath::Abs(cosThetaStarD0) > fDzeroCuts[jPtBin][2] ) isD0=false;
  if (TMath::Abs(cosThetaStarD0bar) > fDzeroCuts[jPtBin][2] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;


  Int_t returnValue=0;
  if(isD0) returnValue+=1;
  if(isD0bar) returnValue+=2;
  return returnValue;

}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::GetPtBin(Double_t ptCand)
{
  for (Int_t i = 0; i < fNPtBins; i++) {
    if (ptCand>=fPtBinLims[i] && ptCand<fPtBinLims[i+1]){
      return i;
    }
  }
  return -1;
}

//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::SelectInvMassAndPt3prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc3)
{

  Int_t retval = 0;
  Double_t momentum[3];
  Double_t px[3], py[3], pz[3];
  for (Int_t iTrack = 0; iTrack < 3; iTrack++) {
    AliESDtrack* track = (AliESDtrack*)trkArray->UncheckedAt(iTrack);
    track->GetPxPyPz(momentum);
    px[iTrack] = momentum[0];
    py[iTrack] = momentum[1];
    pz[iTrack] = momentum[2];
  }
  UInt_t pdg3[3];
  Int_t nprongs = 3;
  rd4massCalc3->SetPxPyPzProngs(nprongs, px, py, pz);
  Double_t minv2, mrange;
  Double_t lolim, hilim;
  mrange = 0.1;
  lolim = fMassDplus - mrange;
  hilim = fMassDplus + mrange;
  pdg3[0] = 211;
  pdg3[1] = 321;
  pdg3[2] = 211;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if (minv2 > lolim * lolim && minv2 < hilim * hilim)
    retval += (1 << kbitDplus);
  lolim = fMassDs - mrange;
  hilim = fMassDs + mrange;
  for (Int_t ih = 0; ih < 2; ih++) {
    Int_t k = ih * 2;
    pdg3[k] = 321;
    pdg3[1] = 321;
    pdg3[2 - k] = 211;
    minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
    if (minv2 > lolim * lolim && minv2 < hilim * hilim && (retval & (1 << kbitDs)) == 0)
      retval += (1 << kbitDs);
  }
  lolim = fMassLambdaC - mrange;
  hilim = fMassLambdaC + mrange;
  pdg3[0] = 2212;
  pdg3[1] = 321;
  pdg3[2] = 211;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if (minv2 > lolim * lolim && minv2 < hilim * hilim && (retval & (1 << kbitLc)) == 0)
    retval += (1 << kbitLc);
  pdg3[0] = 211;
  pdg3[1] = 321;
  pdg3[2] = 2212;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if (minv2 > lolim * lolim && minv2 < hilim * hilim && (retval & (1 << kbitLc)) == 0)
    retval += (1 << kbitLc);

  return retval;
}
//______________________________________________________________________________
AliAODRecoDecayHF2Prong* AliAnalysisTaskHFSimpleVertices::Make2Prong(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t bzkG)
{

  AliESDtrack* track_0 = (AliESDtrack*)twoTrackArray->UncheckedAt(0);
  AliESDtrack* track_1 = (AliESDtrack*)twoTrackArray->UncheckedAt(1);

  Double_t px[2], py[2], pz[2], d0[2], d0err[2];
  Double_t momentum[3];
  GetTrackMomentumAtSecVert(track_0, secVert, momentum, bzkG);
  px[0] = momentum[0];
  py[0] = momentum[1];
  pz[0] = momentum[2];
  GetTrackMomentumAtSecVert(track_1, secVert, momentum, bzkG);
  px[1] = momentum[0];
  py[1] = momentum[1];
  pz[1] = momentum[2];

  Float_t d0z0f[2], covd0z0f[3];
  track_0->GetImpactParameters(d0z0f, covd0z0f);
  d0[0] = d0z0f[0];
  d0err[0] = TMath::Sqrt(covd0z0f[0]);
  track_1->GetImpactParameters(d0z0f, covd0z0f);
  d0[1] = d0z0f[0];
  d0err[1] = TMath::Sqrt(covd0z0f[0]);

  Double_t xdummy, ydummy;
  float dcap1n1 = track_0->GetDCA(track_1, bzkG, xdummy, ydummy);

  AliAODRecoDecayHF2Prong* the2Prong = new AliAODRecoDecayHF2Prong(0x0, px, py, pz, d0, d0err, dcap1n1);
  AliAODVertex* ownsecv=secVert->CloneWithoutRefs();
  the2Prong->SetOwnSecondaryVtx(ownsecv);
  return the2Prong;
}
//______________________________________________________________________________
AliAODRecoDecayHF3Prong* AliAnalysisTaskHFSimpleVertices::Make3Prong(TObjArray* threeTrackArray, AliAODVertex* secVert, Double_t bzkG)
{

  AliESDtrack* track_0 = (AliESDtrack*)threeTrackArray->UncheckedAt(0);
  AliESDtrack* track_1 = (AliESDtrack*)threeTrackArray->UncheckedAt(1);
  AliESDtrack* track_2 = (AliESDtrack*)threeTrackArray->UncheckedAt(2);

  Double_t px[3], py[3], pz[3], d0[3], d0err[3];
  Double_t momentum[3];
  GetTrackMomentumAtSecVert(track_0, secVert, momentum, bzkG);
  px[0] = momentum[0];
  py[0] = momentum[1];
  pz[0] = momentum[2];
  GetTrackMomentumAtSecVert(track_1, secVert, momentum, bzkG);
  px[1] = momentum[0];
  py[1] = momentum[1];
  pz[1] = momentum[2];
  GetTrackMomentumAtSecVert(track_2, secVert, momentum, bzkG);
  px[2] = momentum[0];
  py[2] = momentum[1];
  pz[2] = momentum[2];
  Float_t d0z0f[2], covd0z0f[3];
  track_0->GetImpactParameters(d0z0f, covd0z0f);
  d0[0] = d0z0f[0];
  d0err[0] = TMath::Sqrt(covd0z0f[0]);
  track_1->GetImpactParameters(d0z0f, covd0z0f);
  d0[1] = d0z0f[0];
  d0err[1] = TMath::Sqrt(covd0z0f[0]);
  track_2->GetImpactParameters(d0z0f, covd0z0f);
  d0[2] = d0z0f[0];
  d0err[2] = TMath::Sqrt(covd0z0f[0]);

  Double_t xdummy, ydummy;
  float dcap1n1 = track_0->GetDCA(track_1, bzkG, xdummy, ydummy);
  float dcap2n1 = track_2->GetDCA(track_1, bzkG, xdummy, ydummy);
  float dcap1p2 = track_0->GetDCA(track_2, bzkG, xdummy, ydummy);
  Double_t dca[3] = {dcap1n1, dcap2n1, dcap1p2};
  Double_t dispersion = 0;
  Double_t dist12 = 0.;
  Double_t dist23 = 0.;
  Short_t charge = (Short_t)(track_0->Charge() + track_1->Charge() + track_2->Charge());

  // construct the candidate passing a NULL pointer for the secondary vertex to avoid creation of TRef
  AliAODRecoDecayHF3Prong* the3Prong = new AliAODRecoDecayHF3Prong(0x0, px, py, pz, d0, d0err, dca, dispersion, dist12, dist23, charge);
  // add a pointer to the secondary vertex via SetOwnSecondaryVtx (no TRef created)
  // AliAODVertex* ownsecv=secVert->CloneWithoutRefs();
  // the3Prong->SetOwnSecondaryVtx(ownsecv);
  return the3Prong;
}
//______________________________________________________________________________
char* AliAnalysisTaskHFSimpleVertices::GetJsonString(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  char* value=0x0;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      value=strtok(line, ":");
      value=strtok(NULL, ":");
      break;
    }
  }
  fclose(fj);
  return value;
}
//______________________________________________________________________________
bool AliAnalysisTaskHFSimpleVertices::GetJsonBool(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  bool value=false;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      char* token=strtok(line, ":");
      token=strtok(NULL, ":");
      TString temp=token;
      temp.ReplaceAll("\"","");
      temp.ReplaceAll(",","");
      if(temp.Contains("true")) value=true;
      break;
    }
  }
  fclose(fj);
  return value;
}

//______________________________________________________________________________
int AliAnalysisTaskHFSimpleVertices::GetJsonInteger(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  int value=-999;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      char* token=strtok(line, ":");
      token=strtok(NULL, ":");
      TString temp=token;
      temp.ReplaceAll("\"","");
      temp.ReplaceAll(",","");
      value=temp.Atoi();
      break;
    }
  }
  fclose(fj);
  return value;
}
//______________________________________________________________________________
float AliAnalysisTaskHFSimpleVertices::GetJsonFloat(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  float value=-999.;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      char* token=strtok(line, ":");
      token=strtok(NULL, ":");
      TString temp=token;
      temp.ReplaceAll("\"","");
      temp.ReplaceAll(",","");
      value=temp.Atof();
      break;
    }
  }
  fclose(fj);
  return value;
}

//______________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  printf("AliAnalysisTaskHFSimpleVertices::Terminate --- Number of events: read = %.0f\n",fHistNEvents->GetBinContent(1));
  return;
}





