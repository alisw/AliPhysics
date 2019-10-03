/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//         AliAnalysisTaskV2AllChAOD class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskV2AllChAOD.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODEventCuts.h"
#include "AliPIDCombined.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include <TMCProcess.h>
#include <TRandom.h>

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskV2AllChAOD)

//________________________________________________________________________
AliAnalysisTaskV2AllChAOD::AliAnalysisTaskV2AllChAOD(const char *name) : AliAnalysisTaskSE(name),  
fAOD(0x0),
fTrackCuts(0x0),
fEventCuts(0x0),
fIsMC(0),
fCharge(0),
fVZEROside(0),
fOutput(0x0),
fOutput_lq(0x0),
fOutput_sq(0x0),
fnCentBins(20),
fnQvecBins(100),
fQvecUpperLim(100),
fCutLargeQperc(90.),
fCutSmallQperc(10.),
fEtaGapMin(-0.5),
fEtaGapMax(0.5),
fTrkBit(128),
fEtaCut(0.8),
fMinPt(0.2),
fMaxPt(20.0),
fMinTPCNcls(70),
fFillTHn(kFALSE),
fCentrality(0),
fQvector(0),
fQvector_lq(0),
fQvector_sq(0),
fResSP(0),
fResSP_vs_Cent(0),
fEta_vs_Phi_bef(0),
fEta_vs_PhiA(0),
fEta_vs_PhiB(0),
fResSP_lq(0),
fResSP_vs_Cent_lq(0),
fResSP_sq(0),
fResSP_vs_Cent_sq(0),
fResSP_inclusive(0),
fv2SPGap1A_inclusive_mb(0),
fv2SPGap1B_inclusive_mb(0),
fv2SPGap1A_inclusive_lq(0),
fv2SPGap1B_inclusive_lq(0),
fv2SPGap1A_inclusive_sq(0),
fv2SPGap1B_inclusive_sq(0),
fResSPmc_inclusive(0),
fv2SPGap1Amc_inclusive_mb(0),
fv2SPGap1Bmc_inclusive_mb(0),
fv2SPGap1Amc_inclusive_lq(0),
fv2SPGap1Bmc_inclusive_lq(0),
fv2SPGap1Amc_inclusive_sq(0),
fv2SPGap1Bmc_inclusive_sq(0),
fResGap1w(0),
fV2IntGap1w(0),
fResSP_qbin(0),
fIsRecoEff(0),
fRecoEffList(0),
fQvecGen(0),
fQgenType(0),
fnNchBins(400),
fDoCentrSystCentrality(0)
{

  for (Int_t i = 0; i< 9; i++){

    fv2SPGap1A[i] = 0;
    fv2SPGap1B[i] = 0;

    fSinGap1Aq[i] = 0;
    fCosGap1Aq[i] = 0;
    fSinGap1Bq[i] = 0;
    fCosGap1Bq[i] = 0;

    fSinGap1A[i] = 0;
    fCosGap1A[i] = 0;
    fSinGap1B[i] = 0;
    fCosGap1B[i] = 0;

    //large q
    fv2SPGap1A_lq[i] = 0;
    fv2SPGap1B_lq[i] = 0;

    fSinGap1Aq_lq[i] = 0;
    fCosGap1Aq_lq[i] = 0;
    fSinGap1Bq_lq[i] = 0;
    fCosGap1Bq_lq[i] = 0;

    fSinGap1A_lq[i] = 0;
    fCosGap1A_lq[i] = 0;
    fSinGap1B_lq[i] = 0;
    fCosGap1B_lq[i] = 0;

    //small q
    fv2SPGap1A_sq[i] = 0;
    fv2SPGap1B_sq[i] = 0;

    fSinGap1Aq_sq[i] = 0;
    fCosGap1Aq_sq[i] = 0;
    fSinGap1Bq_sq[i] = 0;
    fCosGap1Bq_sq[i] = 0;

    fSinGap1A_sq[i] = 0;
    fCosGap1A_sq[i] = 0;
    fSinGap1B_sq[i] = 0;
    fCosGap1B_sq[i] = 0;

    fResSP_vs_Qvec[i] = 0;
    fV2IntGap1wq[i] = 0;  

  }
    
  for(Int_t j=0; j<10; j++){
    fv2SPGap1A_qbin[j]=0x0;
    fv2SPGap1B_qbin[j]=0x0;
  }

  fRecoEffList=new TList();
  fRecoEffList->SetOwner();
  fRecoEffList->SetName("fRecoEffList");

  // Default constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliSpectraAODEventCuts::Class());
  DefineOutput(3, AliSpectraAODTrackCuts::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskV2AllChAOD::UserCreateOutputObjects()
{
  // create output objects
  fOutput=new TList();
  fOutput->SetOwner();
  fOutput->SetName("fOutput");

  fOutput_lq=new TList();
  fOutput_lq->SetOwner();
  fOutput_lq->SetName("fOutput_lq");

  fOutput_sq=new TList();
  fOutput_sq->SetOwner();
  fOutput_sq->SetName("fOutput_sq");

  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");

  if( fFillTHn ){ 
    //dimensions of THnSparse for Q vector checks
    const Int_t nvarev=6;
    //                                             cent         q-rec_perc        qvec_v0a       q-rec_v0c        qvec-gen_tpc          Nch
    Int_t    binsHistRealEv[nvarev] = {     fnCentBins,              100,        fnQvecBins,     fnQvecBins,       fnQvecBins,     fnNchBins};
    Double_t xminHistRealEv[nvarev] = {             0.,               0.,                0.,             0.,               0.,            0.};
    Double_t xmaxHistRealEv[nvarev] = {           100.,             100.,     fQvecUpperLim,  fQvecUpperLim,    fQvecUpperLim,         2000.};

    THnSparseF* NSparseHistEv = new THnSparseF("NSparseHistEv","NSparseHistEv",nvarev,binsHistRealEv,xminHistRealEv,xmaxHistRealEv);
    NSparseHistEv->GetAxis(0)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
    NSparseHistEv->GetAxis(0)->SetName(Form("%s_cent",fEventCuts->GetCentralityMethod().Data()));

    NSparseHistEv->GetAxis(1)->SetTitle("q-vec rec percentile");
    NSparseHistEv->GetAxis(1)->SetName("Qrec_perc");

    NSparseHistEv->GetAxis(2)->SetTitle("q-vec V0A");
    NSparseHistEv->GetAxis(2)->SetName("Qrec_V0A");

    NSparseHistEv->GetAxis(3)->SetTitle("q-vec V0C");
    NSparseHistEv->GetAxis(3)->SetName("Qrec_V0C");

    NSparseHistEv->GetAxis(4)->SetTitle("q-vec TPC");
    NSparseHistEv->GetAxis(4)->SetName("Qgen_TPC");

    NSparseHistEv->GetAxis(5)->SetTitle("Ncharged");
    NSparseHistEv->GetAxis(5)->SetName("Nch");
    fOutput->Add(NSparseHistEv);
  }

  fCentrality = new TH1D("fCentrality", "centrality distribution; centrality", 200, 0., 100);
  fOutput->Add(fCentrality);

  fQvector = new TH1D("fQvector", "q-vector distribution; q-vector", fnQvecBins, 0., fQvecUpperLim);
  fOutput->Add(fQvector);

  fQvector_lq = new TH1D("fQvector_lq", "q-vector distribution; q-vector", fnQvecBins, 0., fQvecUpperLim);
  fOutput_lq->Add(fQvector_lq);

  fQvector_sq = new TH1D("fQvector_sq", "q-vector distribution; q-vector", fnQvecBins, 0., fQvecUpperLim);
  fOutput_sq->Add(fQvector_sq);

  // binning common to all the THn
  //change it according to your needs + move it to global variables -> setter/getter
  //   Double_t ptBins[] = {0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0};
  //   const Int_t nptBins = 31;
  Double_t ptBins[] = {0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0};
  const Int_t nptBins = 33;

  fResSP = new TProfile("fResSP", "Resolution; centrality; Resolution", 9, -0.5, 8.5);
  fOutput->Add(fResSP);

  fResSP_vs_Cent = new TProfile("fResSP_vs_Cent", "Resolution; centrality; Resolution", 20., 0., 100.);
  fOutput->Add(fResSP_vs_Cent);

  fEta_vs_Phi_bef = new TH2D("fEta_vs_Phi_bef","eta vs phi distribution before eta gap;#eta;#phi",200.,-1.,1.,175.,0.,7.);
  fOutput->Add(fEta_vs_Phi_bef);

  fEta_vs_PhiA = new TH2D("fEta_vs_PhiA","eta vs phi distribution;#eta;#phi",200.,-1.,1.,175.,0.,7.);
  fOutput->Add(fEta_vs_PhiA);

  fEta_vs_PhiB = new TH2D("fEta_vs_PhiB","eta vs phi distribution;#eta;#phi",200.,-1.,1.,175.,0.,7.);
  fOutput->Add(fEta_vs_PhiB);

  // MC closure test
  fResSP_inclusive = new TProfile("fResSP_inclusive", "Resolution; ese; Resolution", 3, 0., 3.);
  fOutput->Add(fResSP_inclusive);

  fv2SPGap1A_inclusive_mb = new TProfile("fv2SPGap1A_inclusive_mb", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
  fOutput->Add(fv2SPGap1A_inclusive_mb);

  fv2SPGap1B_inclusive_mb = new TProfile("fv2SPGap1B_inclusive_mb", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
  fOutput->Add(fv2SPGap1B_inclusive_mb);


  //large q
  fResSP_lq = new TProfile("fResSP_lq", "Resolution; centrality; Resolution", 9, -0.5, 8.5);
  fOutput_lq->Add(fResSP_lq);

  fResSP_vs_Cent_lq = new TProfile("fResSP_vs_Cent_lq", "Resolution; centrality; Resolution", 20., 0., 100.);
  fOutput_lq->Add(fResSP_vs_Cent_lq);

  // MC closure test
  fv2SPGap1A_inclusive_lq = new TProfile("fv2SPGap1A_inclusive_lq", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
  fOutput_lq->Add(fv2SPGap1A_inclusive_lq);

  fv2SPGap1B_inclusive_lq = new TProfile("fv2SPGap1B_inclusive_lq", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
  fOutput_lq->Add(fv2SPGap1B_inclusive_lq);

  //small q resolution
  fResSP_sq = new TProfile("fResSP_sq", "Resolution; centrality; Resolution", 9, -0.5, 8.5);
  fOutput_sq->Add(fResSP_sq);

  fResSP_vs_Cent_sq = new TProfile("fResSP_vs_Cent_sq", "Resolution; centrality; Resolution", 20., 0., 100.);
  fOutput_sq->Add(fResSP_vs_Cent_sq);

  // MC closure test
  fv2SPGap1A_inclusive_sq = new TProfile("fv2SPGap1A_inclusive_sq", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
  fOutput_sq->Add(fv2SPGap1A_inclusive_sq);

  fv2SPGap1B_inclusive_sq = new TProfile("fv2SPGap1B_inclusive_sq", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
  fOutput_sq->Add(fv2SPGap1B_inclusive_sq);

  for (Int_t iC = 0; iC < 9; iC++){

    fv2SPGap1A[iC] = new TProfile(Form("fv2SPGap1A_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1A[iC]);

    fv2SPGap1B[iC] = new TProfile(Form("fv2SPGap1B_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1B[iC]);

    fSinGap1Aq[iC] = new TProfile(Form("fSinGap1Aq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fSinGap1Aq[iC]);

    fCosGap1Aq[iC] = new TProfile(Form("fCosGap1Aq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fCosGap1Aq[iC]);

    fSinGap1Bq[iC] = new TProfile(Form("fSinGap1Bq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fSinGap1Bq[iC]);

    fCosGap1Bq[iC] = new TProfile(Form("fCosGap1Bq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fCosGap1Bq[iC]);

    fSinGap1A[iC] = new TProfile(Form("fSinGap1A_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fSinGap1A[iC]);

    fCosGap1A[iC] = new TProfile(Form("fCosGap1A_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fCosGap1A[iC]);

    fSinGap1B[iC] = new TProfile(Form("fSinGap1B_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fSinGap1B[iC]);

    fCosGap1B[iC] = new TProfile(Form("fCosGap1B_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput->Add(fCosGap1B[iC]);

    //large q
    fv2SPGap1A_lq[iC] = new TProfile(Form("fv2SPGap1A_lq_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_lq->Add(fv2SPGap1A_lq[iC]);

    fv2SPGap1B_lq[iC] = new TProfile(Form("fv2SPGap1B_lq_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_lq->Add(fv2SPGap1B_lq[iC]);

    fSinGap1Aq_lq[iC] = new TProfile(Form("fSinGap1Aq_lq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fSinGap1Aq_lq[iC]);

    fCosGap1Aq_lq[iC] = new TProfile(Form("fCosGap1Aq_lq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fCosGap1Aq_lq[iC]);

    fSinGap1Bq_lq[iC] = new TProfile(Form("fSinGap1Bq_lq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fSinGap1Bq_lq[iC]);

    fCosGap1Bq_lq[iC] = new TProfile(Form("fCosGap1Bq_lq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fCosGap1Bq_lq[iC]);

    fSinGap1A_lq[iC] = new TProfile(Form("fSinGap1A_lq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fSinGap1A_lq[iC]);

    fCosGap1A_lq[iC] = new TProfile(Form("fCosGap1A_lq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fCosGap1A_lq[iC]);

    fSinGap1B_lq[iC] = new TProfile(Form("fSinGap1B_lq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fSinGap1B_lq[iC]);

    fCosGap1B_lq[iC] = new TProfile(Form("fCosGap1B_lq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_lq->Add(fCosGap1B_lq[iC]);

    //small q
    fv2SPGap1A_sq[iC] = new TProfile(Form("fv2SPGap1A_sq_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_sq->Add(fv2SPGap1A_sq[iC]);

    fv2SPGap1B_sq[iC] = new TProfile(Form("fv2SPGap1B_sq_%d", iC), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_sq->Add(fv2SPGap1B_sq[iC]);

    fSinGap1Aq_sq[iC] = new TProfile(Form("fSinGap1Aq_sq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fSinGap1Aq_sq[iC]);

    fCosGap1Aq_sq[iC] = new TProfile(Form("fCosGap1Aq_sq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fCosGap1Aq_sq[iC]);

    fSinGap1Bq_sq[iC] = new TProfile(Form("fSinGap1Bq_sq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fSinGap1Bq_sq[iC]);

    fCosGap1Bq_sq[iC] = new TProfile(Form("fCosGap1Bq_sq_%d", iC), "p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fCosGap1Bq_sq[iC]);

    fSinGap1A_sq[iC] = new TProfile(Form("fSinGap1A_sq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fSinGap1A_sq[iC]);

    fCosGap1A_sq[iC] = new TProfile(Form("fCosGap1A_sq_%d", iC), ";p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fCosGap1A_sq[iC]);

    fSinGap1B_sq[iC] = new TProfile(Form("fSinGap1B_sq_%d", iC), ";p_{T} (GeV/c);#LT sin(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fSinGap1B_sq[iC]);

    fCosGap1B_sq[iC] = new TProfile(Form("fCosGap1B_sq_%d", iC), "p_{T} (GeV/c);#LT cos(2*#phi) #GT", nptBins, ptBins);
    fOutput_sq->Add(fCosGap1B_sq[iC]);
  }

  // MC closure test
  
  if(fIsMC){
    fResSPmc_inclusive = new TProfile("fResSPmc_inclusive", "Resolution; ese; Resolution", 3, 0., 3.);
    fOutput->Add(fResSPmc_inclusive);

    fv2SPGap1Amc_inclusive_mb = new TProfile("fv2SPGap1Amc_inclusive_mb", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1Amc_inclusive_mb);

    fv2SPGap1Bmc_inclusive_mb = new TProfile("fv2SPGap1Bmc_inclusive_mb", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1Bmc_inclusive_mb);

    //large-q
    fv2SPGap1Amc_inclusive_lq = new TProfile("fv2SPGap1Amc_inclusive_lq", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_lq->Add(fv2SPGap1Amc_inclusive_lq);

    fv2SPGap1Bmc_inclusive_lq = new TProfile("fv2SPGap1Bmc_inclusive_lq", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_lq->Add(fv2SPGap1Bmc_inclusive_lq);

    //small-q
    fv2SPGap1Amc_inclusive_sq = new TProfile("fv2SPGap1Amc_inclusive_sq", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_sq->Add(fv2SPGap1Amc_inclusive_sq);

    fv2SPGap1Bmc_inclusive_sq = new TProfile("fv2SPGap1Bmc_inclusive_sq", "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput_sq->Add(fv2SPGap1Bmc_inclusive_sq);
  }
  
  //v2 vs qvec...
  
  fResGap1w = new TProfile("fResGap1w", "Resolution; centrality; Resolution", 9, -0.5, 8.5);
  fResGap1w->Sumw2();
  fOutput->Add(fResGap1w);

  fV2IntGap1w = new TProfile("fV2IntGap1w", "; centrality; v_{2}", 9, -0.5, 8.5);
  fV2IntGap1w->Sumw2();
  fOutput->Add(fV2IntGap1w);
  
  for (Int_t iC = 0; iC < 9; iC++){
    
    fResSP_vs_Qvec[iC] = new TProfile(Form("fResSP_vs_Qvec_%d", iC), "Resolution; Qvec (V0A); Resolution", 20., 0., 100.);
    fResSP_vs_Qvec[iC]->Sumw2();
    fOutput->Add(fResSP_vs_Qvec[iC]);

    fV2IntGap1wq[iC] = new TProfile(Form("fV2IntGap1wq_%d", iC), "integrated v_{2} vs q-vector; Qvec (V0A); v_{2}", 20., 0., 100.);
    fV2IntGap1wq[iC]->Sumw2();
    fOutput->Add(fV2IntGap1wq[iC]);
    
  }
  
  // v2 vs pt in q-vec bins

  fResSP_qbin = new TProfile("fResSP_qbin", "Resolution; centrality; Resolution", 10, -0.5, 9.5);
  fOutput->Add(fResSP_qbin);
  
  for (Int_t iQ = 0; iQ < 10; iQ++){
    
    fv2SPGap1A_qbin[iQ] = new TProfile(Form("fv2SPGap1A_qbin_%d", iQ), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1A_qbin[iQ]);

    fv2SPGap1B_qbin[iQ] = new TProfile(Form("fv2SPGap1B_qbin_%d", iQ), "v_{2}{2} vs p_{T}; p_{T} (GeV/c); v_{2}{2}", nptBins, ptBins);
    fOutput->Add(fv2SPGap1B_qbin[iQ]);

  }
  
  
  PostData(1, fOutput  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fOutput_lq  );
  PostData(5, fOutput_sq  );
}

//________________________________________________________________________

void AliAnalysisTaskV2AllChAOD::UserExec(Option_t *)
{
  //Printf("An event");
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (!fAOD) {
    AliWarning("ERROR: AliAODEvent not available \n");
    return;
  }

  if (strcmp(fAOD->ClassName(), "AliAODEvent"))
  {
    AliFatal("Not processing AODs");
  }

  if(!fEventCuts->IsSelected(fAOD,fTrackCuts))return;//event selection

  //Get q-vector percentile.
  Double_t Qvec=0.;
  if(fIsMC && fQvecGen) Qvec = fEventCuts->GetQvecPercentileMC(fVZEROside, fQgenType);
  else Qvec = fEventCuts->GetQvecPercentile(fVZEROside);

  fQvector->Fill(Qvec);
  if (Qvec > fCutLargeQperc && Qvec < 100.) fQvector_lq->Fill(Qvec);
  if (Qvec > 0. && Qvec < fCutSmallQperc) fQvector_sq->Fill(Qvec);
  
    Int_t qvecClass = -1;
  if ((Qvec > 0) && (Qvec <= 10.0))
    qvecClass = 0;
  else if ((Qvec > 10.0) && (Qvec <= 20.0))
    qvecClass = 1;
  else if ((Qvec > 20.0) && (Qvec <= 30.0))
    qvecClass = 2;
  else if ((Qvec > 30.0) && (Qvec <= 40.0))
    qvecClass = 3;
  else if ((Qvec > 40.0) && (Qvec <= 50.0))
    qvecClass = 4;
  else if ((Qvec > 50.0) && (Qvec <= 60.0))
    qvecClass = 5;
  else if ((Qvec > 60.0) && (Qvec <= 70.0))
    qvecClass = 6;
  else if ((Qvec > 70.0) && (Qvec <= 80.0))
    qvecClass = 7;
  else if ((Qvec > 80.0) && (Qvec <= 90.0))
    qvecClass = 8;
  else if ((Qvec > 90.0) && (Qvec <= 100.0))
    qvecClass = 9;

  if(qvecClass==-1)return; // FIXME if the qvec is not defined... return!!!
  

  Double_t Cent=(fDoCentrSystCentrality)?1.01*fEventCuts->GetCent():fEventCuts->GetCent();
  fCentrality->Fill(Cent);

  Int_t centV0 = -1;
  if ((Cent > 0) && (Cent <= 5.0))
    centV0 = 0;
  else if ((Cent > 5.0) && (Cent <= 10.0))
    centV0 = 1;
  else if ((Cent > 10.0) && (Cent <= 20.0))
    centV0 = 2;
  else if ((Cent > 20.0) && (Cent <= 30.0))
    centV0 = 3;
  else if ((Cent > 30.0) && (Cent <= 40.0))
    centV0 = 4;
  else if ((Cent > 40.0) && (Cent <= 50.0))
    centV0 = 5;
  else if ((Cent > 50.0) && (Cent <= 60.0))
    centV0 = 6;
  else if ((Cent > 60.0) && (Cent <= 70.0))
    centV0 = 7;
  else if ((Cent > 70.0) && (Cent <= 80.0))
    centV0 = 8;

  if(centV0==-1)return; // FIXME if the centrality is not defined or >80%... return!!!

  if(fIsMC) MCclosure(Qvec); // fill mc histograms for montecarlo closure

  Double_t QxGap1A = 0., QyGap1A = 0.;
  Double_t QxGap1B = 0., QyGap1B = 0.;
  Double_t multGap1A = 0, multGap1B = 0;

  for (Int_t loop = 0; loop < 2; loop++){

    //main loop on tracks
    for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {

      AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));

      if(!track) AliFatal("Not a standard AOD");

      if(fCharge != 0 && track->Charge() != fCharge) continue;//if fCharge != 0 only select fCharge 

      if (!fTrackCuts->IsSelected(track,kTRUE)) continue; //track selection (rapidity selection NOT in the standard cuts)

      if ( track->Pt() < fMinPt || track->Pt() > fMaxPt ) continue;
      
      fEta_vs_Phi_bef->Fill( track->Eta(), track->Phi() );

      if (loop == 0) {

        if (track->Eta() > fEtaGapMax){

          if(fIsRecoEff){
            Double_t receff = GetRecoEff(track->Pt(), centV0);
            QxGap1A += (TMath::Cos(2.*track->Phi()))/receff;
            QyGap1A += (TMath::Sin(2.*track->Phi()))/receff;
            multGap1A+=1./receff;
          } else  {
            QxGap1A += TMath::Cos(2.*track->Phi());
            QyGap1A += TMath::Sin(2.*track->Phi());
            multGap1A++;
          }

          fSinGap1Aq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
          fCosGap1Aq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));

          fEta_vs_PhiA->Fill( track->Eta(), track->Phi() );

          if (Qvec > fCutLargeQperc && Qvec < 100.){
            fSinGap1Aq_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
            fCosGap1Aq_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
          }

          if (Qvec > 0. && Qvec < fCutSmallQperc){
            fSinGap1Aq_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
            fCosGap1Aq_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
          }

        }

        if (track->Eta() < fEtaGapMin){

          if(fIsRecoEff){
            Double_t receff = GetRecoEff(track->Pt(), centV0);
            QxGap1B += (TMath::Cos(2.*track->Phi()))/receff;
            QyGap1B += (TMath::Sin(2.*track->Phi()))/receff;
            multGap1B+=1./receff;
          } else {
            QxGap1B += TMath::Cos(2.*track->Phi());
            QyGap1B += TMath::Sin(2.*track->Phi());
            multGap1B++;
          }

          fCosGap1Bq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
          fSinGap1Bq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));

          fEta_vs_PhiB->Fill( track->Eta(), track->Phi() );

          if (Qvec > fCutLargeQperc && Qvec < 100.){
            fSinGap1Bq_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
            fCosGap1Bq_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
          }

          if (Qvec > 0. && Qvec < fCutSmallQperc){
            fSinGap1Bq_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
            fCosGap1Bq_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
          }

        }

      } else {

        //eval v2 scalar product
        if (track->Eta() < fEtaGapMin && multGap1A > 0){

          Double_t v2SPGap1A = 0.;
          if(fIsRecoEff){
            Double_t receff = GetRecoEff(track->Pt(), centV0);
            Double_t uxGap1A = (TMath::Cos(2.*track->Phi()))/receff;
            Double_t uyGap1A = (TMath::Sin(2.*track->Phi()))/receff;
            //	    multGap1A = multGap1A/receff;
            v2SPGap1A = (uxGap1A*QxGap1A + uyGap1A*QyGap1A)/(Double_t)multGap1A;
          } else v2SPGap1A = (TMath::Cos(2.*track->Phi())*QxGap1A + TMath::Sin(2.*track->Phi())*QyGap1A)/(Double_t)multGap1A;

          //           Double_t v2SPGap1A = (TMath::Cos(2.*track->Phi())*QxGap1A + TMath::Sin(2.*track->Phi())*QyGap1A)/(Double_t)multGap1A;

          fv2SPGap1A[centV0]->Fill(track->Pt(), v2SPGap1A);
          fv2SPGap1A_qbin[qvecClass]->Fill(track->Pt(), v2SPGap1A);

          fv2SPGap1A_inclusive_mb->Fill(track->Pt(), v2SPGap1A); //mb v2 for mc closure

          fSinGap1A[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
          fCosGap1A[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));

          if (Qvec > fCutLargeQperc && Qvec < 100.){
            fv2SPGap1A_lq[centV0]->Fill(track->Pt(), v2SPGap1A);
            fSinGap1A_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
            fCosGap1A_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));

            fv2SPGap1A_inclusive_lq->Fill(track->Pt(), v2SPGap1A); //lq v2 for mc closure
          }

          if (Qvec > 0. && Qvec < fCutSmallQperc){
            fv2SPGap1A_sq[centV0]->Fill(track->Pt(), v2SPGap1A);
            fSinGap1A_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
            fCosGap1A_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));

            fv2SPGap1A_inclusive_sq->Fill(track->Pt(), v2SPGap1A); //sq v2 for mc closure
          }

        }

        if (track->Eta() > fEtaGapMax && multGap1B > 0){

          Double_t v2SPGap1B = 0;
          if(fIsRecoEff){
            Double_t receff = GetRecoEff(track->Pt(), centV0);
            Double_t uxGap1B = (TMath::Cos(2.*track->Phi()))/receff;
            Double_t uyGap1B = (TMath::Sin(2.*track->Phi()))/receff;
//            multGap1B = multGap1B/receff;
            v2SPGap1B = (uxGap1B*QxGap1B + uyGap1B*QyGap1B)/(Double_t)multGap1B;
          }
          else v2SPGap1B = (TMath::Cos(2.*track->Phi())*QxGap1B + TMath::Sin(2.*track->Phi())*QyGap1B)/(Double_t)multGap1B;

          // 	  Double_t v2SPGap1B = (TMath::Cos(2.*track->Phi())*QxGap1B + TMath::Sin(2.*track->Phi())*QyGap1B)/(Double_t)multGap1B;

          fv2SPGap1B[centV0]->Fill(track->Pt(), v2SPGap1B);
          fv2SPGap1B_qbin[qvecClass]->Fill(track->Pt(), v2SPGap1B);

          fv2SPGap1B_inclusive_mb->Fill(track->Pt(), v2SPGap1B); //mb v2 for mc closure

          fCosGap1B[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));
          fSinGap1B[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));

          if (Qvec > fCutLargeQperc && Qvec < 100.){
            fv2SPGap1B_lq[centV0]->Fill(track->Pt(), v2SPGap1B);
            fSinGap1B_lq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
            fCosGap1B_lq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));

            fv2SPGap1B_inclusive_lq->Fill(track->Pt(), v2SPGap1B); //lq v2 for mc closure
          }

          if (Qvec > 0. && Qvec < fCutSmallQperc){
            fv2SPGap1B_sq[centV0]->Fill(track->Pt(), v2SPGap1B);
            fSinGap1B_sq[centV0]->Fill(track->Pt(), TMath::Sin(2.*track->Phi()));
            fCosGap1B_sq[centV0]->Fill(track->Pt(), TMath::Cos(2.*track->Phi()));

            fv2SPGap1B_inclusive_sq->Fill(track->Pt(), v2SPGap1B); //sq v2 for mc closure
          }

        }
      }// end else 
    } // end loop on tracks
  } // end loop


  if (multGap1A > 0 && multGap1B > 0){
    Double_t res = (QxGap1A*QxGap1B + QyGap1A*QyGap1B)/(Double_t)multGap1A/(Double_t)multGap1B;
    fResSP->Fill((Double_t)centV0, res);
    fResSP_qbin->Fill((Double_t)qvecClass, res);

    fResSP_inclusive->Fill(0., res); //mb v2 for mc closure
    fResSP_vs_Cent->Fill(Cent, res);

    if (Qvec > fCutLargeQperc && Qvec < 100.){
      fResSP_lq->Fill((Double_t)centV0, res);
      fResSP_vs_Cent_lq->Fill(Cent, res);

      fResSP_inclusive->Fill(1., res); //lq v2 for mc closure
    }

    if (Qvec > 0. && Qvec < fCutSmallQperc){
      fResSP_sq->Fill((Double_t)centV0, res);
      fResSP_vs_Cent_sq->Fill(Cent, res);

      fResSP_inclusive->Fill(2., res); //sq v2 for mc closure
    }
    
    //_________________________________________________________________
    //v2 vs qvec
    fResGap1w->Fill(Double_t(centV0), (QxGap1A*QxGap1B + QyGap1A*QyGap1B)/multGap1A/multGap1B, (Double_t)(multGap1A*multGap1B));

    Double_t nGap1 = multGap1A*multGap1B;
    Double_t qqGap1 = (QxGap1A*QxGap1B + QyGap1A*QyGap1B)/nGap1;

    fV2IntGap1w->Fill(Double_t(centV0), qqGap1, nGap1);

    fResSP_vs_Qvec[centV0]->Fill(Double_t(Qvec), (QxGap1A*QxGap1B + QyGap1A*QyGap1B)/multGap1A/multGap1B, (Double_t)(multGap1A*multGap1B));
    fV2IntGap1wq[centV0]->Fill(Double_t(Qvec), qqGap1, nGap1);
    
  }// end multiplicity if


  if( fFillTHn ){ 

    Double_t varEv[6];
    varEv[0]=Cent;

    varEv[1]=(Double_t)Qvec; // qvec_rec_perc

    varEv[2]=(Double_t)fEventCuts->GetqV0A();

    varEv[3]=(Double_t)fEventCuts->GetqV0C();

    varEv[4]=(Double_t)fEventCuts->GetqTPC();

    varEv[5]=(Double_t)fEventCuts->GetNch(); // Nch

    ((THnSparseF*)fOutput->FindObject("NSparseHistEv"))->Fill(varEv);//event loop

  }


  PostData(1, fOutput  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fOutput_lq  );
  PostData(5, fOutput_sq  );
}

//_________________________________________________________________
Bool_t  AliAnalysisTaskV2AllChAOD::GetDCA(const AliAODTrack* trk, Double_t * p){

  //AliAODTrack::DCA(): for newest AOD fTrack->DCA() always gives -999. This should fix.
  //FIXME should update EventCuts?
  //FIXME add track->GetXYZ(p) method

  double xyz[3],cov[3];

  if (!trk->GetXYZ(xyz)) { // dca is not stored
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(trk);
    AliVEvent* ev = (AliVEvent*)trk->GetEvent();
    if (!ev) {/*printf("Event is not connected to the track\n");*/ return kFALSE;}
    if (!etp.PropagateToDCA(ev->GetPrimaryVertex(), ev->GetMagneticField(),999,xyz,cov)) return kFALSE; // failed, track is too far from vertex
  }
  p[0] = xyz[0];
  p[1] = xyz[1];
  return kTRUE;

}

//_________________________________________________________________
void  AliAnalysisTaskV2AllChAOD::MCclosure(Double_t qvec){
  // First do MC to fill up the MC particle array

  TClonesArray *arrayMC = 0;
  if (fIsMC)
  {
    arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if (!arrayMC) {
      AliFatal("Error: MC particles branch not found!\n");
    }

    Double_t QxGap1Amc = 0., QyGap1Amc = 0.;
    Double_t QxGap1Bmc = 0., QyGap1Bmc = 0.;
    Int_t multGap1Amc = 0, multGap1Bmc = 0;

    for (Int_t loop = 0; loop < 2; loop++){

      Int_t nMC = arrayMC->GetEntries();

      for (Int_t iMC = 0; iMC < nMC; iMC++)
      {
        AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
        if(!partMC->Charge()) continue;//Skip neutrals
        if(fCharge != 0 && partMC->Charge()*fCharge < 0.) continue;//if fCharge != 0 only select fCharge

        if (!(partMC->IsPhysicalPrimary()))
          continue;

        if(partMC->Eta()<fTrackCuts->GetEtaMin() || partMC->Eta()>fTrackCuts->GetEtaMax()) continue;

        //Printf("a particle");


        if (loop == 0) {

          if (partMC->Eta() > fEtaGapMax){
            QxGap1Amc += TMath::Cos(2.*partMC->Phi());
            QyGap1Amc += TMath::Sin(2.*partMC->Phi());
            multGap1Amc++;
          }

          if (partMC->Eta() < fEtaGapMin){
            QxGap1Bmc += TMath::Cos(2.*partMC->Phi());
            QyGap1Bmc += TMath::Sin(2.*partMC->Phi());
            multGap1Bmc++;
          }

        } else {

          //eval v2 scalar product
          if (partMC->Eta() < fEtaGapMin && multGap1Amc > 0){
            Double_t v2SPGap1Amc = (TMath::Cos(2.*partMC->Phi())*QxGap1Amc + TMath::Sin(2.*partMC->Phi())*QyGap1Amc)/(Double_t)multGap1Amc;
            fv2SPGap1Amc_inclusive_mb->Fill(partMC->Pt(), v2SPGap1Amc);

            if (qvec > fCutLargeQperc && qvec < 100.){
              fv2SPGap1Amc_inclusive_lq->Fill(partMC->Pt(), v2SPGap1Amc);
            }

            if (qvec > 0. && qvec < fCutSmallQperc){
              fv2SPGap1Amc_inclusive_sq->Fill(partMC->Pt(), v2SPGap1Amc);
            }

          }

          if (partMC->Eta() > fEtaGapMax && multGap1Bmc > 0){
            Double_t v2SPGap1Bmc = (TMath::Cos(2.*partMC->Phi())*QxGap1Bmc + TMath::Sin(2.*partMC->Phi())*QyGap1Bmc)/(Double_t)multGap1Bmc;
            fv2SPGap1Bmc_inclusive_mb->Fill(partMC->Pt(), v2SPGap1Bmc);

            if (qvec > fCutLargeQperc && qvec < 100.){
              fv2SPGap1Bmc_inclusive_lq->Fill(partMC->Pt(), v2SPGap1Bmc);
            }

            if (qvec > 0. && qvec < fCutSmallQperc){
              fv2SPGap1Bmc_inclusive_sq->Fill(partMC->Pt(), v2SPGap1Bmc);
            }

          }


        }// end else
      } // end loop on partMCs
    } // end loop

    if (multGap1Amc > 0 && multGap1Bmc > 0){
      Double_t resmc = (QxGap1Amc*QxGap1Bmc + QyGap1Amc*QyGap1Bmc)/(Double_t)multGap1Amc/(Double_t)multGap1Bmc;
      fResSPmc_inclusive->Fill(0.,resmc);

      if (qvec > fCutLargeQperc && qvec < 100.){
        fResSPmc_inclusive->Fill(1.,resmc);
      }

      if (qvec > 0. && qvec < fCutSmallQperc){
        fResSPmc_inclusive->Fill(2.,resmc);
      }

    }

  }// end if MC
}

//_________________________________________________________________
Double_t AliAnalysisTaskV2AllChAOD::GetRecoEff(Double_t pt, Int_t iC){

  if(iC>8) return 1.;

  if(pt<0.2 || pt>100.) return 1.;

  // // //   //spectra ese binning
  // // //   //const Double_t ptBins[] = {0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,5.0,6.0,7.0,8.0,9.0,10.,12.,15.,20.,25.,30.,35.,40.,50.,75.,100.};
  // // //   //const Int_t nptBins=34;

  const Double_t fEpsilon=0.000001;
  //   cout<<"list: "<<endl;
  //   fRecoEffList->ls();
  TH1D *h = (TH1D*)fRecoEffList->At(iC);

  Int_t bin = h->FindBin(pt);

  Double_t lowlim = h->GetBinLowEdge(bin);
  Double_t uplim = h->GetBinLowEdge(bin) + h->GetBinWidth(bin);

  Double_t eff = 1.;

  if( pt>lowlim && pt<uplim ) eff = h->GetBinContent(bin);
  if( pt == lowlim ) eff = h->GetBinContent( h->FindBin(pt+fEpsilon) );
  if( pt == uplim ) eff = h->GetBinContent( h->FindBin(pt-fEpsilon) );

  return eff;

}
//_________________________________________________________________
void   AliAnalysisTaskV2AllChAOD::Terminate(Option_t *)
{
  // Terminate
}
