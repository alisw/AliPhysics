/*************************************************************************
* Copyright(c) 1998-2008,ALICE Experiment at CERN,All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use,copy,modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee,provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/


/**************************************************************************************************
* AliAnalysisTaskSEpPbCorrelationsJetV2:
* This task is developed from AliAnalysisTaskSEpPbCorrelationsYS.cxx, aimed to calculate the TPC pairs 
* - FMD correlation, which can obtain the V2 of jet particles.
**************************************************************************************************/


#include "AliAnalysisManager.h"
#include "TGrid.h"
#include "AliLog.h"
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TH3.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMap.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TTree.h>
#include <TProfile.h>
#include <TObjectTable.h>
#include "AliAnalysisUtils.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliCFContainer.h"
#include "AliGenEventHeader.h"
#include "AliTHn.h"

#include "AliAODEvent.h"
#include "AliESDAD.h"
//#include "AliESDEvent.h"
#include "AliVAD.h"
#include "AliAODForwardMult.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"
#include "AliESDVertex.h"
//#include "AliAODPid.h"
#include "AliAODHandler.h"
//#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliAODcascade.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliCentrality.h"
#include "AliEventPoolManager.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMultiplicity.h"
//#include "AliPID.h"
//#include "AliPIDResponse.h"
#include "AliVParticle.h"
#include "Riostream.h"

//#include "AliFlowEventSimple.h"
///#include "AliFlowVector.h"/
//#include "AliFlowTrackSimple.h"
//#include "AliAnalysisTaskMSP.h"
//#include "AliFlowAnalysisWithMSP.h"
#include "AliMultSelection.h"
#include "AliFMDCorrSecondaryMap.h"
#include "AliForwardCorrectionManager.h"
#include "AliTrigAssoPairST.h"
//#include "AliForwardUtil.h"

#include "AliAnalysisTaskSEpPbCorrelationsJetV2.h"
#include "AliAnalysisTaskSEpPbCorrelationsYS.h"
using namespace std;

ClassImp(AliAnalysisTaskSEpPbCorrelationsJetV2)
ClassImp(AliAssociatedTrackYS)
ClassImp(AliTrigAssoPairST)
ClassImp(AliMixTrackYS)
ClassImp(AliAssociatedVZEROYS)

AliAnalysisTaskSEpPbCorrelationsJetV2::AliAnalysisTaskSEpPbCorrelationsJetV2()
    : AliAnalysisTaskSE(),
      fcollisiontype("pPb"),
      fDataType(kTRUE),
      frun2(kTRUE),
      fQA(kTRUE),
      fFMDcut(kTRUE),
      fFMDcutmode(1),
      fReduceDphi(1.5707),
      fAnaMode("V0AV0C"),
      fasso("Phi"),
      fSymmetricFMD(kFALSE),
      fLikeSign(kTRUE),
      fCentType("ZNA"),
      fprimFMD(kTRUE),
      fprimTPC(kTRUE),
      fcentcalib(kTRUE),
      fNEntries(0),
      lCentrality(0),
      bSign(0),
      fZVertex(10.),
      fOutputList(0),
      fOutputList1(0),
      fOutputList2(0),
      //fPIDResponse(0),
      ffilterbit(5),
      fnoClusters(70),
      fCutChargedDCAzMax(5),
      fCutChargedDCAxyMax(15),
      fPtMin(0.2),
      fPtMax(3.0),
      fAsscoptCut(0.5),
      fCenMin(0.),
      fCenMax(10.),
      fEtaMax(0.8),
      fEtaMaxExtra(0.),
      fEtaMinExtra(-0.8),
      fEventCuts(0),
      fUtils(),
      fEvent(0),
      mcEvent(0),
      lPrimaryBestVtx(0),
      fPrimaryZVtx(0),
      fvzero(0),
      fPoolMgr(0),
      fPoolMaxNEvents(2000),
      fPoolMinNTracks(50000),
      fMinEventsToMix(5),
      fNzVtxBins(0),
      fNCentBins(0),
      fh2_pt_trig_asso(0),
      fHistzvertex(0),
      fHistCentrality(0),
      fHistCentrality_beforecut(0),
      fHistCentrality_beforeFMDMulcut(0),
      fHistCentzvertex(0),
      mixedDist(0),
      mixedDist2(0),
      //fHistLeadQA(0),      
      frefetaa(0),
      frefetac(0),
      frefvz(0),
      fh2_FMD_eta_phi(0),
      fh2_FMD_eta_phi_afterCut(0),
      fHist_NeventRun(0),
      fHist_V0AMultRun(0),
      fHist_V0CMultRun(0),
      fHist_FMDAMultRun(0),
      fHist_FMDCMultRun(0),
      fhistfmdphiacc(0),  
      fhFMDmultchannel(0),
      fhFMDmultchannel_actual(0),
      fhistfmd(0),
      fhistits(0), 
      fhSecFMD(0),
      fFMDV0(0),
      fFMDV0_post(0),
      fFMDV0A(0),
      fFMDV0A_post(0),
      fFMDV0C(0),
      fFMDV0C_post(0),
      fHist_Stat(0),
      fHistTriggerTrack(0),
      fHistReconstTrack(0),
      fHistTriggerTrackMix(0),
      fHistReconstTrackMix(0){
      for (Int_t iBin = 0; iBin < 100; iBin++) {
    fZvtxBins[iBin] = 0.;
    fCentBins[iBin] = 0.;
  }
  for(Int_t i=0;i<4;i++){
    fhrefetaFMD[i]=0;
    fhrefphiFMD[i]=0;
  }
  for(Int_t i=0;i<10;i++){
    fhcorr[i]=0;
  }
    
}
AliAnalysisTaskSEpPbCorrelationsJetV2::AliAnalysisTaskSEpPbCorrelationsJetV2(const char *name)
    : AliAnalysisTaskSE(name),
      fcollisiontype("pPb"),
      fDataType(kTRUE),
      frun2(kTRUE),
      fQA(kTRUE),
      fFMDcut(kTRUE),
      fFMDcutmode(1),
      fReduceDphi(1.5707),
      fAnaMode("V0AV0C"),
      fasso("Phi"),
      fSymmetricFMD(kFALSE), 
      fLikeSign(kTRUE), 
      fCentType("ZNA"),
      fprimFMD(kTRUE),
      fprimTPC(kTRUE),
      fcentcalib(kTRUE),
      fNEntries(0),
      lCentrality(0),
      bSign(0),
      fZVertex(10.),
      fOutputList(0),
      fOutputList1(0),
      fOutputList2(0),
      //fPIDResponse(0),
      ffilterbit(5),
      fnoClusters(70),
      fCutChargedDCAzMax(5),
      fCutChargedDCAxyMax(15),
      fPtMin(0.2),
      fPtMax(3.0),
      fAsscoptCut(0.5),
      fCenMin(0.),
      fCenMax(10.),
      fEtaMax(0.8),
      fEtaMaxExtra(0.),
      fEtaMinExtra(-0.8),
      fEventCuts(0),
      fUtils(),
      fEvent(0),
      mcEvent(0),
      lPrimaryBestVtx(0),
      fPrimaryZVtx(0),
      fvzero(0),
      fPoolMgr(0),
      fPoolMaxNEvents(2000),
      fPoolMinNTracks(50000),
      fMinEventsToMix(5),
      fNzVtxBins(0),
      fNCentBins(0),
      fh2_pt_trig_asso(0),
      fHistzvertex(0),
      fHistCentrality(0),
      fHistCentrality_beforecut(0),
      fHistCentrality_beforeFMDMulcut(0),
      fHistCentzvertex(0),
      mixedDist(0),
      mixedDist2(0),
      //fHistLeadQA(0),      
      frefetaa(0),
      frefetac(0),
      frefvz(0),
      fh2_FMD_eta_phi(0),
      fh2_FMD_eta_phi_afterCut(0),
      fHist_NeventRun(0),
      fHist_V0AMultRun(0),
      fHist_V0CMultRun(0),
      fHist_FMDAMultRun(0),
      fHist_FMDCMultRun(0),
      fhistfmdphiacc(0),  
      fhFMDmultchannel(0),
      fhFMDmultchannel_actual(0),
      fhistfmd(0),
      fhistits(0), 
      fhSecFMD(0),
      fFMDV0(0),
      fFMDV0_post(0),
      fFMDV0A(0),
      fFMDV0A_post(0),
      fFMDV0C(0),
      fFMDV0C_post(0),
      fHist_Stat(0),
      fHistTriggerTrack(0),
      fHistReconstTrack(0),
      fHistTriggerTrackMix(0),
      fHistReconstTrackMix(0){
      for (Int_t iBin = 0; iBin < 100; iBin++) {
          fZvtxBins[iBin] = 0.;
          fCentBins[iBin] = 0.;
        }
        for(Int_t i=0;i<4;i++){
          fhrefetaFMD[i]=0;
          fhrefphiFMD[i]=0;
        }
        for(Int_t i=0;i<10;i++){
          fhcorr[i]=0;
        }
        DefineOutput(1, TList::Class());
        DefineOutput(2, TList::Class());
        DefineOutput(3, TList::Class());
      }

AliAnalysisTaskSEpPbCorrelationsJetV2::~AliAnalysisTaskSEpPbCorrelationsJetV2()
{
  if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList;
    fOutputList = 0x0;
  }
  
  if (fOutputList1 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList1;
    fOutputList1 = 0x0;
  }
  
  if (fOutputList2 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList2;
    fOutputList2 = 0x0;
  }
}
void AliAnalysisTaskSEpPbCorrelationsJetV2::UserCreateOutputObjects() {
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  fOutputList->SetName("global");
  DefineGeneralOutput();
  PostData(1, fOutputList);

  fOutputList1 = new TList();
  fOutputList1->SetOwner(kTRUE);
  fOutputList1->SetName("anahistos");
  DefineCorrOutput();

  //DefineVZEROOutput();
  PostData(2, fOutputList1);

  fOutputList2 = new TList();
  fOutputList2->SetOwner(kTRUE);
  fOutputList2->SetName("QA");
  DefinedQAHistos();
  PostData(3, fOutputList2);

  fEventCuts.AddQAplotsToList(fOutputList);

  if(fcentcalib){

   TGrid::Connect("alien://");
   //TFile*file=TFile::Open("alien:///alice/cern.ch/user/y/ysekiguc/fcalibration_centrality_AMPT_modireco.root");
   TFile*file=TFile::Open("alien:///alice/cern.ch/user/s/sitang/AMPT_Centrality_Calibration/fcalibration_centrality_AMPT_modireco.root");

   if(!file) AliError("No correction factor");
   fhcorr[0]=(TH1D*)file->Get("hcent");
   fOutputList->Add(fhcorr[0]);
  }

  frefetac=new TH1F("frefetac","frefetac",30,-3.4,-1.9);
  fOutputList2->Add(frefetac);
  frefetaa=new TH1F("frefetaa","frefetaa",62,1.9,5.0);
  fOutputList2->Add(frefetaa);
  
  frefvz=new TH1F("frefvz","z-vertex",10,-10,10);
  fOutputList2->Add(frefvz);

  fOutputList2->Add(new TH2D("hPt_Cen_mid","", 100, 0., 20., 100, 0., 100.)); // For Rcp analysis
  fOutputList2->Add(new TH2D("hPt_Cen_mid_prim","", 100, 0., 20., 100, 0., 100.)); // For Rcp analysis
 
   fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins, fNzVtxBins, fZvtxBins);
   if (!fPoolMgr) return;
   fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);
 }
 void AliAnalysisTaskSEpPbCorrelationsJetV2::DefineGeneralOutput() {

   fHist_Stat = new TH1F("fHist_Stat", "Stat Histogram", 14, -0.5, 13.5);
   fHist_Stat->GetXaxis()->SetBinLabel(1, "All Events");
   fHist_Stat->GetXaxis()->SetBinLabel(2, "Analyzed Events");
   fHist_Stat->GetXaxis()->SetBinLabel(3, "MultSelection OK");
   fHist_Stat->GetXaxis()->SetBinLabel(4, "Vertex OK");
   fHist_Stat->GetXaxis()->SetBinLabel(5, "Centrality OK");
   fHist_Stat->GetXaxis()->SetBinLabel(6, "FMD OK");
   fHist_Stat->GetXaxis()->SetBinLabel(7, "FMD multi cut");
   fHist_Stat->GetXaxis()->SetBinLabel(8, "--");
   fHist_Stat->GetXaxis()->SetBinLabel(9, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(10, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(11, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(12, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(13, "-");
   fHist_Stat->GetXaxis()->SetBinLabel(14, "-");
   fOutputList->Add(fHist_Stat);
   
   fHistzvertex = new TH1F("fHistzvertex", ";VZ;count", 60, -15, 15);
   fOutputList->Add(fHistzvertex);

   Double_t fmaxcent;
   if(fCentType=="Manual"){
     fmaxcent=200;
   }else{
     if (fcollisiontype=="HMPP") fmaxcent=1.;
     else fmaxcent=100.;
   }
   fHistCentrality = new TH1F("fHistCentrality", ";centrality;count", 100, 0, fmaxcent);
   fOutputList->Add(fHistCentrality);

   fHistCentrality_beforecut = new TH1F("fHistCentrality_beforecut", ";centrality;count", 100, 0, fmaxcent);  // Before FMD multiplicity > 0 cut
   fOutputList->Add(fHistCentrality_beforecut);
 
   fHistCentrality_beforeFMDMulcut = new TH1F("fHistCentrality_beforeFMDMulcut", ";centrality;count", 100, 0, fmaxcent); // Before FMD/V0 cut
   fOutputList->Add(fHistCentrality_beforeFMDMulcut);
 
   fHistCentzvertex = new TH2F("fHistCentzvertex", "Cent;VZ;count", 100,0, fmaxcent, 60, -15, 15);
   fOutputList->Add(fHistCentzvertex);

  }

 void AliAnalysisTaskSEpPbCorrelationsJetV2::DefinedQAHistos() {
   Int_t ncentmax;
   if(fCentType=="Manual") ncentmax=200;
   else ncentmax=100;
      
   mixedDist=new TH2F("mixedDist", ";centrality;tracks;events", 101, 0, ncentmax, 200, 0, fPoolMinNTracks*1.5 );
   mixedDist2=new TH2F("mixedDist2", ";centrality;events;events", 101, 0,ncentmax, 100, 0, 1000) ;
   fOutputList2->Add(mixedDist);
   fOutputList2->Add(mixedDist2);

	 fHist_NeventRun=new TH1F("fHist_NeventRun","fHist_NeventRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_NeventRun);

	 fHist_V0AMultRun=new TH1F("fHist_V0AMultRun","fHist_V0AMultRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_V0AMultRun);
	 fHist_V0CMultRun=new TH1F("fHist_V0CMultRun","fHist_V0CMultRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_V0CMultRun);
	 fHist_FMDAMultRun=new TH1F("fHist_FMDAMultRun","fHist_FMDAMultRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_FMDAMultRun);
	 fHist_FMDCMultRun=new TH1F("fHist_FMDCMultRun","fHist_FMDCMultRun",200,-0.5,199.5);
	 fOutputList2->Add(fHist_FMDCMultRun);


	 Float_t nmutFoward = 1000.;
	 
	 fFMDV0 = new TH2F("FMDV0", "FMD vs V0 pre cut;FMD;V0;",2000, 0, 2*nmutFoward, 2000, 0, 2*nmutFoward);
	 fOutputList2->Add(fFMDV0);
	 fFMDV0_post=new TH2F("FMDV0_post", "FMD vs V0 post cut;FMD;V0;",2000, 0, 2*nmutFoward, 2000, 0, 2*nmutFoward);
	 fOutputList2->Add(fFMDV0_post);
	 fFMDV0A = new TH2F("FMDV0A", "FMD vs V0A;FMD;V0A;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
	 fOutputList2->Add(fFMDV0A);
	 fFMDV0A_post = new TH2F("FMDV0A_post", "FMD vs V0A post cut;FMD;V0A;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
	 fOutputList2->Add(fFMDV0A_post);
	 fFMDV0C = new TH2F("FMDV0C", "FMD vs V0C;FMD;V0C;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
	 fOutputList2->Add(fFMDV0C);
	 fFMDV0C_post = new TH2F("FMDV0C_post", "FMD vs V0C post cut;FMD;V0C;",1000, 0, nmutFoward, 1000, 0, nmutFoward);
	 fOutputList2->Add(fFMDV0C_post);
	 
	 fh2_FMD_eta_phi=new TH2D("fh2_FMD_eta_phi","fh2_FMD_eta_phi",200,-4,6,20,0,2*TMath::Pi());
	 fOutputList2->Add(fh2_FMD_eta_phi);

         fh2_FMD_eta_phi_afterCut=new TH2D("fh2_FMD_eta_phi_After_FMDV0_cut","fh2_FMD_eta_phi_After_FMDV0_cut",200,-4,6,20,0,2*TMath::Pi());
         fOutputList2->Add(fh2_FMD_eta_phi_afterCut);
	 
	 fhistfmdphiacc=new TH2D("fhistfmdphiacc","fhistfmdphiacc",200,-4,6,20,0,ncentmax);
	 fOutputList2->Add(fhistfmdphiacc);
	 
	 fhFMDmultchannel=new TH2F("fhFMDmultchannel","fhFMDmultchannel",200,-4,6,100,0,100);
	 fOutputList2->Add(fhFMDmultchannel);
           
         fhFMDmultchannel_actual=new TH2F("fhFMDmultchannel_actual","fhFMDmultchannel_actual",200,-4,6,100,0,100);
         fOutputList2->Add(fhFMDmultchannel_actual);       
	 
	 	 
	 const Int_t ifmdbin[4]={200,20,20,20};
	 fhistfmd=new AliTHn("fhistfmd","fhistfmd",1,4,ifmdbin);
	 fhistfmd->SetBinLimits(0,-4.,6.);
	 fhistfmd->SetBinLimits(1,0.,2*TMath::Pi());
	 if(fCentType=="Manual")fhistfmd->SetBinLimits(2,0.,200.);
	 else fhistfmd->SetBinLimits(2,0.,100.);
	 fhistfmd->SetBinLimits(3,-10.,10.);
	 fhistfmd->SetVarTitle(0,"eta");
	 fhistfmd->SetVarTitle(1,"phi");
	 fhistfmd->SetVarTitle(2,"centrality");
	 fhistfmd->SetVarTitle(3,"vzy");
// old	 fOutputList2->Add(fhistfmd);
	 
	 const Int_t iitsbin[4]={200,20,20,40};
	 const	 Double_t MinITS[4]={-4,0.,0.,-10.};
	 const	 Double_t MaxITS[4]={6,2*TMath::Pi(),100.,10.};
	 fhistits=new THnSparseF("fhistits","fhistits",4,iitsbin,MinITS,MaxITS);
	 // old  fOutputList2->Add(fhistits);
	 
	 const Double_t binning_etafmd[51]={
					    -3.4,-3.3,-3.2,-3.1,-3.0,
					    -2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,
					    -1.9,-1.8,-1.7,
					    1.7,1.8,1.9,
					    2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
					    3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
					    4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9};
	 
	  }

 void AliAnalysisTaskSEpPbCorrelationsJetV2::DefineCorrOutput() {

//========================= For Jet v2
   Double_t binning_pt_assoc[] = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 200.};
   Int_t nbins_pt_assoc = sizeof(binning_pt_assoc) / sizeof(Double_t) - 1; // For TPC-TPC

   Double_t binning_pt_lead[] = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0};
   Int_t nbins_pt_lead  = sizeof(binning_pt_lead) / sizeof(Double_t) - 1;
//==========================

//=========================== For Inc v2
   Double_t binning_pt_leadTPC[] = {0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 5.0, 6.0, 7.0, 10.0};
   Int_t nbins_pt_leadTPC  = sizeof(binning_pt_leadTPC) / sizeof(Double_t) - 1;

   Double_t binning_dphi[73] = {
       -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464,
       -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865,
       -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266,
       0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332,
       0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931,
       1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530,
       1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129,
       2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727,
       2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326,
       3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925,
       3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524,
       4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123,
       4.712389};

//=========================

    Int_t nCFStepstrig=1;

    if(fAnaMode=="TPCTPCFMDA" || fAnaMode=="TPCTPCFMDC")
    {
     Double_t binning_deta_tpctpc_trig[] = {-1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1., -0.9, -0.8, -0.7 ,-0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
     Int_t ndetatpctpc_trig = sizeof(binning_deta_tpctpc_trig)/sizeof(Double_t)-1;

     const Int_t iEvtBinFMD[] = {nbins_pt_lead, ndetatpctpc_trig, 18, 10};
     const Int_t nEvtVarsFMD = sizeof(iEvtBinFMD) / sizeof(Int_t); 
     
     fHistTriggerTrack = new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVarsFMD, iEvtBinFMD);
     
     fHistTriggerTrack->SetBinLimits(0, binning_pt_lead);
     //fHistTriggerTrack->SetBinLimits(1, binning_pt_assoc);
     fHistTriggerTrack->SetBinLimits(1, binning_deta_tpctpc_trig);
     fHistTriggerTrack->SetBinLimits(2, -0.5*TMath::Pi(),1.5*TMath::Pi());
     fHistTriggerTrack->SetBinLimits(3, -10.,10.);

     fHistTriggerTrack->SetVarTitle(0, "leading p_{T} GeV/c");
     //fHistTriggerTrack->SetVarTitle(1, "associate p_{T} GeV/c");
     fHistTriggerTrack->SetVarTitle(1, "#Delta#eta");
     fHistTriggerTrack->SetVarTitle(2, "#Delta#phi");
     fHistTriggerTrack->SetVarTitle(3, "zvertex");
    }     

//======================= AliTHn to Store Trigger particle in same event
   if(fAnaMode=="TPCFMDA" || fAnaMode=="TPCFMDC")
   {
    const Int_t iEvtBinTPCFMD[] = {nbins_pt_leadTPC, 10};
    const Int_t nEvtVarsTPCFMD = sizeof(iEvtBinTPCFMD) / sizeof(Int_t);
    fHistTriggerTrack = new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVarsTPCFMD, iEvtBinTPCFMD);
    fHistTriggerTrack->SetBinLimits(0, binning_pt_leadTPC);
    fHistTriggerTrack->SetBinLimits(1, -10.,10.);
   
    fHistTriggerTrack->SetVarTitle(0, "leading p_{T} GeV/c");
    fHistTriggerTrack->SetVarTitle(1, "zvertex");
   }

   if(fAnaMode=="FMDAFMDC")
   {
    Double_t binning_cent_fmdfmd[] = {0.,10., 60.,100.};
    Int_t ncentbin = sizeof(binning_cent_fmdfmd)/sizeof(Double_t) - 1;

    const Int_t iEvtBinFMDAFMDC[] = {10};
    const Int_t nEvtVarsFMDAFMDC = sizeof(iEvtBinFMDAFMDC) / sizeof(Int_t);
    fHistTriggerTrack = new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVarsFMDAFMDC, iEvtBinFMDAFMDC);
    //fHistTriggerTrack->SetBinLimits(0, binning_cent_fmdfmd);
    fHistTriggerTrack->SetBinLimits(0, -10.,10.);
    
    //fHistTriggerTrack->SetVarTitle(0, "Centrality");
    fHistTriggerTrack->SetVarTitle(0, "zvertex");
   }
 

   if(fAnaMode=="TPCTPC")
   {
    const Int_t nEvtVarsTPC = 3;
    const Int_t iEvtBinTPC[] = {nbins_pt_lead, 10, 12};
    
    fHistTriggerTrack = new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFStepstrig, nEvtVarsTPC, iEvtBinTPC);

    fHistTriggerTrack->SetBinLimits(0, binning_pt_lead);
    fHistTriggerTrack->SetBinLimits(1, -10.,10.);
    fHistTriggerTrack->SetBinLimits(2, 0.,12.);

    fHistTriggerTrack->SetVarTitle(0, "leading p_{T} GeV/c");
    fHistTriggerTrack->SetVarTitle(1, "zvertex");
    fHistTriggerTrack->SetVarTitle(2, "Random Number");
   }
//=======================

//======================= AliTHn to Store Trigger particle in mix event

   const Int_t nEvtVars = 2;
   const Int_t iEvtBin[2] = {nbins_pt_lead, 12};
    
   fHistTriggerTrackMix = new AliTHn("fHistTriggerTrackMix", "fHistTriggerTrackMix", nCFStepstrig, nEvtVars, iEvtBin);
   fHistTriggerTrackMix->SetBinLimits(0, binning_pt_lead);
   fHistTriggerTrackMix->SetBinLimits(1, 0.,12.);
   fHistTriggerTrackMix->SetVarTitle(0, "leading p_{T} GeV/c");
   fHistTriggerTrackMix->SetVarTitle(1, "Random Number");


   fOutputList1->Add(fHistTriggerTrack);
   fOutputList1->Add(fHistTriggerTrackMix);
//==========================   
   //////////////////////////////////////////
   //Containers two particle correlation
   //////////////////////////////////////////

   Int_t nCFSteps = 1;

   if(fAnaMode=="TPCTPC") {
     Double_t binning_deta_tpctpc[33] = {-1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5,
					 -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5, 
                                         0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5, 1.6};
     const Int_t iTrackBin_TPCTPC[5] = {32,nbins_pt_assoc,nbins_pt_lead, 72, 10};
     fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, 5, iTrackBin_TPCTPC);
     fHistReconstTrack->SetBinLimits(0, binning_deta_tpctpc);
     fHistReconstTrack->SetBinLimits(1, binning_pt_assoc);
     fHistReconstTrack->SetBinLimits(2, binning_pt_lead);
     fHistReconstTrack->SetBinLimits(3, binning_dphi);
     fHistReconstTrack->SetBinLimits(4, -10,10);
     fHistReconstTrack->SetVarTitle(0, "#Delta#eta");
     fHistReconstTrack->SetVarTitle(1, "p_{T} GeV/c");
     fHistReconstTrack->SetVarTitle(2, "leading p_{T} GeV/c");
     fHistReconstTrack->SetVarTitle(3, "#Delta#phi");
     fHistReconstTrack->SetVarTitle(4, "Vz");
     fHistReconstTrackMix =  new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, 5, iTrackBin_TPCTPC);
     fHistReconstTrackMix->SetBinLimits(0, binning_deta_tpctpc);
     fHistReconstTrackMix->SetBinLimits(1, binning_pt_assoc);
     fHistReconstTrackMix->SetBinLimits(2, binning_pt_lead);
     fHistReconstTrackMix->SetBinLimits(3, binning_dphi);
     fHistReconstTrackMix->SetBinLimits(4, -10,10);
     fHistReconstTrackMix->SetVarTitle(0, "#Delta#eta");
     fHistReconstTrackMix->SetVarTitle(1, "p_{T} GeV/c");
     fHistReconstTrackMix->SetVarTitle(2, "leading p_{T} GeV/c");
     fHistReconstTrackMix->SetVarTitle(3, "#Delta#phi");
     fHistReconstTrackMix->SetVarTitle(4, "Vz");
   }else if (fAnaMode=="TPCTPCFMDA" || fAnaMode=="TPCTPCFMDC"){

    Double_t binning_deta_tpctpc[] = {-1.6, -1.3, -1.1, -0.9, -0.7 ,-0.5, -0.3, -0.1,  0.1,  0.3,  0.5, 0.7, 0.9, 1.1, 1.3, 1.6};
    Int_t ndetatpctpc = sizeof(binning_deta_tpctpc)/sizeof(Double_t)-1;            
     
   
   Double_t binning_dphi_reduce[] = {
       -1.570796,  -1.221730,
       -0.872665, -0.523599,
       -0.174533,  0.174533,  0.523599, 
       0.872665,  1.221730,   1.570796,  
       1.919862,  2.268928,   2.617994,  
       2.967060,  3.316126,  3.665191,  
       4.014257,  4.363323,   4.712389};
    Int_t nbinning_dphi_reduce = sizeof(binning_dphi_reduce)/sizeof(Double_t) - 1;
    	
  
     Double_t binning_detaFMDTPC[]={-6.,-5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4., -3.8, -3.6, -3.4, -3.2, -3., -2.8, -2.6, -2.4, -2.2, -2., -1.8, -1.6, -1.4, -1.2, -1., -0.8};
     Double_t binning_detaFMDCTPC[]={ 1., 1.2, 1.4, 1.6, 1.8, 2. , 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4.};

     Int_t ndetatpcfmd;
     if(fAnaMode=="TPCTPCFMDA") {
       ndetatpcfmd= sizeof(binning_detaFMDTPC)/sizeof(Double_t) - 1;
     }else{
       ndetatpcfmd = sizeof(binning_detaFMDCTPC)/sizeof(Double_t) - 1;
     }

     Int_t nCFStepstpcfmd=1;
     //Int_t iTrackBin_tpcfmd[7]={ndetatpctpc, nbins_pt_assoc, nbins_pt_lead, 10, ndetatpcfmd, nbinning_dphi_reduce, 24};
     Int_t iTrackBin_tpcfmd[6]={ndetatpctpc, nbins_pt_lead, 10, ndetatpcfmd, nbinning_dphi_reduce, 24};

     fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFStepstpcfmd, 6, iTrackBin_tpcfmd);    
     fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFStepstpcfmd, 6, iTrackBin_tpcfmd);

 //only for test

     if(fAnaMode=="TPCTPCFMDA") 
     {
      fHistReconstTrack->SetBinLimits(3,binning_detaFMDTPC);
      fHistReconstTrackMix->SetBinLimits(3,binning_detaFMDTPC);
     }
     else if(fAnaMode=="TPCTPCFMDC") 
     { 
      fHistReconstTrack->SetBinLimits(3,binning_detaFMDCTPC);
      fHistReconstTrackMix->SetBinLimits(3,binning_detaFMDCTPC);
     }
     
     fHistReconstTrack->SetBinLimits(0, binning_deta_tpctpc);
     fHistReconstTrack->SetBinLimits(1, binning_pt_lead);
     fHistReconstTrack->SetBinLimits(2,-10.,10.);
     fHistReconstTrack->SetBinLimits(4,binning_dphi_reduce);
     fHistReconstTrack->SetBinLimits(5, -0.5*TMath::Pi(),1.5*TMath::Pi());
     
     fHistReconstTrack->SetVarTitle(0,"#Delta#eta_{TPC-TPC}");
     fHistReconstTrack->SetVarTitle(1,"p_{T}_{TPC-trig} GeV/c");
     fHistReconstTrack->SetVarTitle(2,"z vertex");
     fHistReconstTrack->SetVarTitle(3,"#Delta#eta_{TPC-FMD}");
     fHistReconstTrack->SetVarTitle(4,"#Delta#phi_{TPC-FMD}");
     fHistReconstTrack->SetVarTitle(5,"#Delta#phi_{TPC-TPC}");

     fHistReconstTrackMix->SetBinLimits(0,binning_deta_tpctpc);
     fHistReconstTrackMix->SetBinLimits(1,binning_pt_lead);
     fHistReconstTrackMix->SetBinLimits(2,-10.,10.);
     fHistReconstTrackMix->SetBinLimits(4,binning_dphi_reduce);
     fHistReconstTrackMix->SetBinLimits(5, -0.5*TMath::Pi(),1.5*TMath::Pi());

     fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta_{TPC-TPC}");
     fHistReconstTrackMix->SetVarTitle(1,"p_{T}_{TPC-trig} GeV/c");
     fHistReconstTrackMix->SetVarTitle(2,"z vertex");
     fHistReconstTrackMix->SetVarTitle(3,"#Delta#eta_{TPC-FMD}");
     fHistReconstTrackMix->SetVarTitle(4,"#Delta#phi_{TPC-FMD}");
     fHistReconstTrackMix->SetVarTitle(5,"#Delta#phi_{TPC-TPC}");
   }
   else if(fAnaMode=="TPCFMDA" || fAnaMode=="TPCFMDC")
   {
    Double_t binning_detaTPCFMDA[]={-6.,-5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4., -3.8, -3.6, -3.4, -3.2, -3., -2.8, -2.6, -2.4, -2.2, -2., -1.8, -1.6, -1.4, -1.2, -1., -0.8};
    Double_t binning_detaTPCFMDC[]={ 1., 1.2, 1.4, 1.6, 1.8, 2. , 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4.};

    Int_t nCFStepstpcfmd=1;

    Int_t iTrackBin_tpcfmdA[] = {nbins_pt_leadTPC, 10, 26, 72};
    Int_t iTrackBin_tpcfmdC[] = {nbins_pt_leadTPC, 10, 15, 72};
    Int_t nTrackBin_tpcfmd = sizeof(iTrackBin_tpcfmdA) / sizeof(Int_t);
    
    if(fAnaMode=="TPCFMDA")
    {
     fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFStepstpcfmd, nTrackBin_tpcfmd, iTrackBin_tpcfmdA);
     fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFStepstpcfmd, nTrackBin_tpcfmd, iTrackBin_tpcfmdA);
     fHistReconstTrack->SetBinLimits(2, -6., -0.8);
     fHistReconstTrackMix->SetBinLimits(2, -6., -0.8);
    }
    else
    {
     fHistReconstTrack = new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFStepstpcfmd, nTrackBin_tpcfmd, iTrackBin_tpcfmdC);
     fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFStepstpcfmd, nTrackBin_tpcfmd, iTrackBin_tpcfmdC);
     fHistReconstTrack->SetBinLimits(2, 1., 4.);
     fHistReconstTrackMix->SetBinLimits(2, 1., 4.);
    }
      
    fHistReconstTrack->SetBinLimits(0, binning_pt_leadTPC);
    fHistReconstTrack->SetBinLimits(1, -10, 10);
    fHistReconstTrack->SetBinLimits(3, binning_dphi);
  
    fHistReconstTrack->SetVarTitle(0,"p_{T}_{TPC-trig} GeV/c");
    fHistReconstTrack->SetVarTitle(1,"z vertex");
    fHistReconstTrack->SetVarTitle(2,"#Delta#eta_{TPC-FMD}");
    fHistReconstTrack->SetVarTitle(3,"#Delta#phi_{TPC-FMD}");
    
    fHistReconstTrackMix->SetBinLimits(0, binning_pt_leadTPC);
    fHistReconstTrackMix->SetBinLimits(1, -10, 10);
    fHistReconstTrackMix->SetBinLimits(3, binning_dphi);
    
    fHistReconstTrackMix->SetVarTitle(0,"p_{T}_{TPC-trig} GeV/c");
    fHistReconstTrackMix->SetVarTitle(1,"z vertex");
    fHistReconstTrackMix->SetVarTitle(2,"#Delta#eta_{TPC-FMD}");
    fHistReconstTrackMix->SetVarTitle(3,"#Delta#phi_{TPC-FMD}");
   }
   else if(fAnaMode == "FMDAFMDC")
   {
    Double_t binning_cent_fmdfmd[] = {0.,10., 60.,100.};
    Int_t ncentbin=sizeof(binning_cent_fmdfmd)/sizeof(Double_t) - 1;

    Int_t nCFStepstpcfmd=1;

    Int_t iTrackBin_fmdAfmdC[] = {10, 24, 20};
    Int_t nTrackBin_fmdAfmdC = sizeof(iTrackBin_fmdAfmdC) / sizeof(Int_t);

    fHistReconstTrack = new AliTHn("fHistReconstTrack","fHistReconstTrack",nCFStepstpcfmd,nTrackBin_fmdAfmdC,iTrackBin_fmdAfmdC);
    fHistReconstTrackMix = new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFStepstpcfmd, nTrackBin_fmdAfmdC, iTrackBin_fmdAfmdC);

    fHistReconstTrack->SetBinLimits(0, -10, 10);
    fHistReconstTrack->SetBinLimits(1, 3.4,8.2);
    fHistReconstTrack->SetBinLimits(2, -0.55*TMath::Pi(),1.45*TMath::Pi());

    fHistReconstTrack->SetVarTitle(0,"z vertex");
    fHistReconstTrack->SetVarTitle(1,"#Delta#eta_{FMDA-FMDC}");
    fHistReconstTrack->SetVarTitle(2,"#Delta#phi_{FMDA-FMDC}");   

    fHistReconstTrackMix->SetBinLimits(0, -10, 10);
    fHistReconstTrackMix->SetBinLimits(1, 3.4,8.2);
    fHistReconstTrackMix->SetBinLimits(2, -0.5*TMath::Pi(),1.5*TMath::Pi());
   
    fHistReconstTrackMix->SetVarTitle(0,"z vertex");
    fHistReconstTrackMix->SetVarTitle(1,"#Delta#eta_{FMDA-FMDC}");
    fHistReconstTrackMix->SetVarTitle(2,"#Delta#phi_{FMDA-FMDC}");
   }
  

   fOutputList1->Add(fHistReconstTrack);
   fOutputList1->Add(fHistReconstTrackMix);
 }

 void AliAnalysisTaskSEpPbCorrelationsJetV2::UserExec(Option_t *) {

   DumpTObjTable("Start analysis");
   AliAnalysisManager *mgr        = AliAnalysisManager::GetAnalysisManager();

// Bottommost Selection for Events
   AliInputEventHandler *inEvMain = (AliInputEventHandler *)(mgr->GetInputEventHandler());
   if (!inEvMain)    return;

   fEvent = dynamic_cast<AliAODEvent *>(inEvMain->GetEvent());
   if (!fEvent) 
   {
    AliWarning("ERROR: fEvent not available \n");
    return;
   }


   fHist_Stat->Fill(0);

// Default Events cuts for Run-2
   if(fcollisiontype=="pPb")
   {
    if(!fEventCuts.AcceptEvent(fEvent))   // Events Cut 
    {
     PostData(1, fOutputList);
     PostData(2, fOutputList1);
     PostData(3, fOutputList2);
     return;
    }
   }
   
   fHist_Stat->Fill(1);      

// Multiplicity Selection for Run-2
   AliMultSelection *multSelection = (AliMultSelection *)fEvent->FindListObject("MultSelection");
   if(!multSelection) return;
   fHist_Stat->Fill(2);
 
   

// The Vz cut for events
   lPrimaryBestVtx = fEvent->GetPrimaryVertex();
   fPrimaryZVtx = lPrimaryBestVtx->GetZ();
/* 
   if(lPrimaryBestVtx->IsFromVertexer3D()) cout << "It is 3D reconstruction" << endl;
   if(lPrimaryBestVtx->IsFromVertexerZ())  cout << "It is Vz reconstruction" << endl;
   if(!lPrimaryBestVtx->IsFromVertexer3D() && !lPrimaryBestVtx->IsFromVertexerZ()) cout << "It is TPC reconstruction" << endl; 
*/
   //lPrimaryBestVtx->Print();
   if((TMath::Abs(fPrimaryZVtx)) >= fZVertex)   
   {
    PostData(1, fOutputList);
    PostData(2, fOutputList1);
    PostData(3, fOutputList2);
    return;
   }
   fHist_Stat->Fill(3);

// Magnetic Field of the event
   bSign = 0.;
   bSign = (InputEvent()->GetMagneticField() > 0) ? 1 : -1;

   // Multiplicity Calibration for MC
   fvzero = fEvent->GetVZEROData();
   if(fcentcalib)
   { 
    Float_t sum = 0., max = 0.;
    for(Int_t i = 32; i < 64; ++i)
    {      
     sum +=fvzero->GetMultiplicity(i);
     if (fvzero->GetMultiplicity(i) > max) max = fvzero->GetMultiplicity(i);
    }
    sum -= max;
     //     fV0Amultmodi->Fill(sum);
    
    Int_t nbinmult= fhcorr[0]->GetXaxis()->FindBin(sum);
    lCentrality=fhcorr[0]->GetBinContent(nbinmult);
   }
   else
   {
    lCentrality = multSelection->GetMultiplicityPercentile(fCentType);
    Int_t qual = multSelection->GetEvSelCode();
    if (qual == 199)  lCentrality = -999;
    if (lCentrality < 0. || lCentrality > 100. - 0.0000001)   return;
   }

   fHist_Stat->Fill(4);  // Basic centrality

   //===================== Particle Level Information
   if(!fDataType && fprimFMD)
   {
    AliAODMCHeader* aodMCheader=(AliAODMCHeader*)fEvent->FindListObject(AliAODMCHeader::StdBranchName());
    TClonesArray *mcArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
    if(!mcArray)
    {
     Printf("No MC particle branch found");
     return;
    }
    
    Int_t nMCAllTracks = mcArray->GetEntriesFast();
    for (Int_t i = 0; i < nMCAllTracks; i++)
    {
     AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(i);
     if (!mcTrack) 
     {
      Error("ReadEventAODMC", "Could not receive particle %d", i);
      continue;
     }
     Bool_t TrIsPrim=mcTrack->IsPhysicalPrimary();
     Float_t mcTrackEta = mcTrack->Eta();
     Float_t mcTrackPt  = mcTrack->Pt();
     Bool_t TrCharge=mcTrack->Charge()!=0;
     if(!TrCharge)        continue;
     if(!TrIsPrim)        continue;
     if(mcTrackEta<0.8 && mcTrackEta>-0.8) (static_cast<TH2D*>(fOutputList2 ->FindObject("hPt_Cen_mid_prim")))->Fill(mcTrackPt,lCentrality);
    }
   }
  //===================== For Rcp
  Int_t nTracks = fEvent->GetNumberOfTracks();
  for(Int_t i = 0; i < nTracks; i++)
  {
   AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(fEvent->GetTrack(i));
   if (!aodTrack)      continue;
   if (!IsAcceptedTrack(aodTrack))      continue;
   if (aodTrack->Charge() == 0)      continue; 
   Double_t aodTrack_pt = aodTrack->Pt();
   (static_cast<TH2D*>(fOutputList2 ->FindObject("hPt_Cen_mid")))->Fill(aodTrack_pt,lCentrality);
  }
  if(fAnaMode=="Rcp")
  {
   PostData(1, fOutputList);
   PostData(2, fOutputList1);
   PostData(3, fOutputList2);
   return;
  }
  //=====================

   fHistCentrality_beforecut->Fill(lCentrality);
   
   //   fHist_Stat->Fill(5);
   DumpTObjTable("After event selection");

   MakeAna();
      
   PostData(1, fOutputList);
   PostData(2, fOutputList1);
   PostData(3, fOutputList2);
 }

 void AliAnalysisTaskSEpPbCorrelationsJetV2::Terminate(Option_t *) {
   //  AliInfo(Form("Number of Correlation
   DumpTObjTable("End of the analysis");
   Printf("Entries======================%d",fNEntries);
   if (fPoolMgr)    delete fPoolMgr;   // PoolMgr->ClearPools();
 }

void AliAnalysisTaskSEpPbCorrelationsJetV2::MakeAna() {

   DumpTObjTable("start correlation analysis");

   TObjArray *selectedTracksLeading = new TObjArray;
   selectedTracksLeading->SetOwner(kTRUE);

   TObjArray *selectedTracksAssociated = new TObjArray;
   selectedTracksAssociated->SetOwner(kTRUE);

   TObjArray* selectedTracksMC1=new TObjArray;
   selectedTracksMC1->SetOwner(kTRUE);
   TObjArray* selectedTracksMC2=new TObjArray;
   selectedTracksMC2->SetOwner(kTRUE);


   Float_t nFMD_fwd_hits=0;
   Float_t nFMD_bwd_hits=0;
   Float_t nFMD_fwdV0acc_hits=0;
   Float_t nFMD_bwdV0acc_hits=0;
   AliAODForwardMult*aodForward=static_cast<AliAODForwardMult*>(fEvent->FindListObject("Forward"));

   // Shape of d2Ndetadphi: 200, -4, 6, 20, 0, 2pi
   Int_t ivzbin=frefvz->GetXaxis()->FindBin(fPrimaryZVtx);
   const TH2D& d2Ndetadphi = aodForward->GetHistogram();
   TH1*hphiacceptance=aodForward->GetPhiAcceptance();

   Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
   Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
   //cout << "n eta "<< nEta << "n phi "<<nPhi << endl;
   Double_t pt = 0;

   for (Int_t iEta = 1; iEta <= nEta; iEta++) 
   {
     Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
     if (!valid) continue;       
     Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
     Float_t phiacc=hphiacceptance->GetBinContent(iEta);
     fhistfmdphiacc->Fill(eta,lCentrality,phiacc);
     
     for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) 
     {
      // Bin content is most likely number of particles!
      Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);	 
      Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
	 
	 if (mostProbableN > 0) {
	   if(eta>0){
	     nFMD_fwd_hits+=mostProbableN;
	     if(2.8<eta && eta<5.03) nFMD_fwdV0acc_hits+=mostProbableN;
	   }else{
	     nFMD_bwd_hits+=mostProbableN;
	     if(-3.4<eta && eta<-2.01) nFMD_bwdV0acc_hits+=mostProbableN;
	   }
	 }

	 if (mostProbableN > 0) {
	   if(eta < 4.8 && eta > 1.8){
             if (!fSymmetricFMD || (eta < 3.2))
             {                     
	      Int_t nfmdetabin1=frefetaa->FindBin(eta);
	      if(fAnaMode=="TPCTPCFMDA" || fAnaMode=="TPCFMDA") selectedTracksAssociated->Add(new AliAssociatedTrackYS(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));	
              if(fAnaMode=="FMDAFMDC") selectedTracksLeading->Add(new AliAssociatedTrackYS(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));
              fhFMDmultchannel_actual->Fill(eta,mostProbableN);
             }
	   }else if( eta < -1.8 && eta > -3.2){
	      Int_t nfmdetabin=frefetac->FindBin(eta);
	      if(fAnaMode=="TPCTPCFMDC" || fAnaMode=="ITSFMDC" ||fAnaMode=="FMDAFMDC" || fAnaMode == "TPCFMDC") selectedTracksAssociated->Add(new AliAssociatedTrackYS(-999,eta,phi,-999,-999,-999,-999,-999,mostProbableN));
              fhFMDmultchannel_actual->Fill(eta,mostProbableN);
           }
	   fhFMDmultchannel->Fill(eta,mostProbableN);
	   Double_t cont[4]={eta,phi,lCentrality,fPrimaryZVtx};
	   fhistfmd->Fill(cont,0,mostProbableN);
	   fh2_FMD_eta_phi->Fill(eta,phi,mostProbableN);
	 }
     }
    }
     delete hphiacceptance;
     
     //FMD events cuts: FMD hits > 0
    if(fFMDcut){ 
     if(nFMD_fwd_hits==0 || nFMD_bwd_hits==0){
       selectedTracksLeading->Clear();
       delete selectedTracksLeading;
       selectedTracksAssociated->Clear();
       delete selectedTracksAssociated;
       selectedTracksMC1->Clear();
       delete selectedTracksMC1;
       selectedTracksMC2->Clear();
       delete selectedTracksMC2;
       PostData(1, fOutputList);
       PostData(2, fOutputList1);
       PostData(3, fOutputList2);
       return;
     } //events cuts
    }  

     fHist_Stat->Fill(5);
     fHistCentrality_beforeFMDMulcut->Fill(lCentrality);  
    
     DumpTObjTable("End of fill fmd tracks");
     
     Float_t nV0A_hits = fvzero->GetMTotV0A();
     Float_t nV0C_hits = fvzero->GetMTotV0C();
     
     fFMDV0->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
     fFMDV0A->Fill(nFMD_fwd_hits, nV0A_hits);
     fFMDV0C->Fill(nFMD_bwd_hits, nV0C_hits);
     
     if(fFMDcut){
       Double_t FMDcutapar0=0.;
       Double_t FMDcutapar1=0.;
       Double_t FMDcutcpar0=0.;
       Double_t FMDcutcpar1=0.;
       switch(fFMDcutmode){
       case 1:
	 FMDcutapar0=1.3;
	 FMDcutapar1=200;
	 FMDcutcpar0=2.;
	 FMDcutcpar1=200;
       	 break;
       case 2:
	 FMDcutapar0=1.3;
	 FMDcutapar1=600;
	 FMDcutcpar0=2.;
	 FMDcutcpar1=600;
	 break;
       case 3:
	 FMDcutapar0=1.5;
	 FMDcutapar1=100;
	 FMDcutcpar0=2.3;
	 FMDcutcpar1=100;
	 break;
       case 4:
	 FMDcutapar0=1.5;
	 FMDcutapar1=300;
	 FMDcutcpar0=2.3;
	 FMDcutcpar1=300;
	 break;
       case 5:
	 FMDcutapar0=1.3;
	 FMDcutapar1=400;
	 FMDcutcpar0=2.;
	 FMDcutcpar1=400;
	 break;
       case 6:
         FMDcutapar0=1.75;
         FMDcutapar1=150;
         FMDcutcpar0=1.4;
         FMDcutcpar1=120;
         break;
       case 7:
         FMDcutapar0=1.64755;
         FMDcutapar1=119.602;
         FMDcutcpar0=2.73426;
         FMDcutcpar1=150.31;
         break;
       case 8:
         FMDcutapar0=1.64755;
         FMDcutapar1=159.47;
         FMDcutcpar0=2.73426;
         FMDcutcpar1=200.413;
         break;
       case 9:
         FMDcutapar0=1.2031;
         FMDcutapar1=73.123;
         FMDcutcpar0=2.25453;
         FMDcutcpar1=104.941;
         break;
       case 10:
         FMDcutapar0=1.2031;
         FMDcutapar1=97.4973;
         FMDcutcpar0=2.25453;
         FMDcutcpar1=139.921;
         break;
       case 11://pp 2 sigma cut
         FMDcutapar0=1.2031;
         FMDcutapar1=48.7486;
         FMDcutcpar0=2.25453;
         FMDcutcpar1=69.9606;
         break;
       case 12://Pbp 2 sigma cut old
         FMDcutapar0=1.64755;
         FMDcutapar1=79.7346;
         FMDcutcpar0=2.73426;
         break;
       case 13://pPb 1 sigma cut
         FMDcutapar0=1.64755;
         FMDcutapar1=40.907;
         FMDcutcpar0=2.73426;
         FMDcutcpar1=36.5323;
         break;
       case 14://pPb 2 sigma cut
         FMDcutapar0=1.64755;
         FMDcutapar1=80.7744;
         FMDcutcpar0=2.73426;
         FMDcutcpar1=86.6355;
         break;
       case 15://pPb 2 sigma cut
         FMDcutapar0=1.64755;
         FMDcutapar1=80.7744;
         FMDcutcpar0=2.73426;
         FMDcutcpar1=86.6355;
         break;
       case 16://PbPb 3sigma cut
         FMDcutapar0=1.26713;
         FMDcutapar1=211.719;
         FMDcutcpar0=2.47928;
         FMDcutcpar1=250.645;
         break;
       case 17://PbPb 4sigma cut
         FMDcutapar0=1.26713;
         FMDcutapar1=274.066;
         FMDcutcpar0=2.47928;
         FMDcutcpar1=338.14;
         break;
     default: break;
     }

// only for pPb
     if((nV0A_hits<(FMDcutapar0*nFMD_fwd_hits-FMDcutapar1)) || (nV0C_hits<(FMDcutcpar0*nFMD_bwd_hits-FMDcutcpar1))){
	   selectedTracksLeading->Clear();
	   delete selectedTracksLeading;
	   selectedTracksAssociated->Clear();
	   delete selectedTracksAssociated;
           selectedTracksMC1->Clear();
           delete selectedTracksMC1;
           selectedTracksMC2->Clear();
           delete selectedTracksMC2;
	   PostData(1, fOutputList);
	   PostData(2, fOutputList1);
	   PostData(3, fOutputList2);
	   return;
       }
     
     fFMDV0_post->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
     fFMDV0A_post->Fill(nFMD_fwd_hits, nV0A_hits);
     fFMDV0C_post->Fill(nFMD_bwd_hits, nV0C_hits);

//==================== Fill the 2D map after the FMD/V0 cut 

     for (Int_t iEta = 1; iEta <= nEta; iEta++)
     {
      Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
      if (!valid) continue;
      Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
      for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++)
      {
       Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
       Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
       if (mostProbableN > 0) fh2_FMD_eta_phi_afterCut->Fill(eta,phi,mostProbableN);
      }
     }
  }
   fHist_Stat->Fill(6);

  
   fHistCentrality->Fill(lCentrality);
   //Centrality Selection  
   if(lCentrality>fCenMax || lCentrality<fCenMin) 
   {
    selectedTracksLeading->Clear();
    delete selectedTracksLeading;
    selectedTracksAssociated->Clear();
    delete selectedTracksAssociated;
    selectedTracksMC1->Clear();
    delete selectedTracksMC1;
    selectedTracksMC2->Clear();
    delete selectedTracksMC2;
    PostData(1, fOutputList);
    PostData(2, fOutputList1);
    PostData(3, fOutputList2);
    return;
   }

   //=============================================================================== primary FMD
  if(!fDataType && fprimFMD)
  {
   Int_t pdgcode=0;
   Bool_t TrIsPrim=kFALSE;
   Bool_t TrIsSecondMate=kFALSE;
   Bool_t TrIsSecondWeak=kFALSE;
   Bool_t TrIsOthers=kFALSE;
   Double_t mcTrackEta=-999.;
   Double_t mcTrackPt=-999.;
   Double_t mcTrackPhi=-999.;
   Bool_t TrCharge=kFALSE;

   AliAODMCHeader* aodMCheader=(AliAODMCHeader*)fEvent->FindListObject(AliAODMCHeader::StdBranchName());
   TClonesArray *mcArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
   if(!mcArray){
    Printf("No MC particle branch found");
    return;
   }

   Int_t nMCAllTracks = mcArray->GetEntriesFast();
   Int_t nMCtrackssamecut=0;
   Int_t ntrackv0aall=0;
   Int_t ntrackv0cprimary=0;
   
   for (Int_t i = 0; i < nMCAllTracks; i++)
   {
     AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(i);
     if (!mcTrack) {
      Error("ReadEventAODMC", "Could not receive particle %d", i);
      continue;
     }

     TrIsPrim=mcTrack->IsPhysicalPrimary();
     TrIsSecondMate=mcTrack->IsSecondaryFromMaterial();
     TrIsSecondWeak=mcTrack->IsSecondaryFromWeakDecay();
     TrIsOthers=!TrIsPrim && !TrIsSecondMate && !TrIsSecondWeak;

     mcTrackEta = mcTrack->Eta();
     mcTrackPt  = mcTrack->Pt();
     mcTrackPhi = mcTrack->Phi();
     TrCharge=mcTrack->Charge()!=0;
     pdgcode=TMath::Abs(mcTrack->PdgCode());

     if(!TrCharge)        continue;
     if(mcTrackPt>fPtMax) continue;

     if(TrIsPrim)
     {
      if(fAnaMode=="TPCTPC") 
      {
       if(mcTrackEta<-0.8 || mcTrackEta>0.8) continue;
       if(mcTrackPt<fPtMin || mcTrackPt>fPtMax) continue;
       selectedTracksMC1->Add(new AliAssociatedTrackYS(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
       selectedTracksMC2->Add(new AliAssociatedTrackYS(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
      }
      else if(fAnaMode=="TPCTPCFMDA"|| fAnaMode=="TPCFMDA" || fAnaMode=="TPCFMDC" || fAnaMode=="TPCTPCFMDC")
      {
       if(mcTrackEta>-0.8 && mcTrackEta<0.8)
       {
        if(mcTrackPt<fPtMin || mcTrackPt>fPtMax) continue;
        selectedTracksMC1->Add(new AliAssociatedTrackYS(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
       }
       else
       {
        if(fAnaMode=="TPCTPCFMDA" || fAnaMode=="TPCFMDA")
        {
         if(mcTrackEta>1.8 && mcTrackEta<4.8)  selectedTracksMC2->Add(new AliAssociatedTrackYS(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
        }
        else if(fAnaMode=="TPCTPCFMDC" || fAnaMode=="TPCFMDC" ||fAnaMode=="FMDAFMDC")
        {
          if(mcTrackEta>-3.2  && mcTrackEta<-1.8) selectedTracksMC2->Add(new AliAssociatedTrackYS(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
        }
       }
      }
      else if(fAnaMode == "FMDAFMDC")
      {
       if(mcTrackEta> 1.8 && mcTrackEta<4.8 ) selectedTracksMC1->Add(new AliAssociatedTrackYS(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
       if(mcTrackEta>-3.2 && mcTrackEta<-1.8) selectedTracksMC2->Add(new AliAssociatedTrackYS(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
      }
     }
   }
  }
//===============================================================================
   fHistzvertex->Fill(fPrimaryZVtx);
   fHistCentzvertex->Fill(lCentrality, fPrimaryZVtx);
   
   DumpTObjTable("End of FMD vs V0 cuts");


// Associate Particle   
if(fAnaMode=="TPCTPC"){
    selectedTracksAssociated=GetAcceptedTracksLeading(fEvent,kFALSE,selectedTracksAssociated);
 }
// Leading Particle
 if(fAnaMode=="TPCTPCFMDA" || fAnaMode=="TPCTPC" || fAnaMode=="TPCTPCFMDC" || fAnaMode =="TPCFMDC" || fAnaMode == "TPCFMDA"){
   selectedTracksLeading=GetAcceptedTracksLeading(fEvent,kTRUE,selectedTracksLeading);
 }
 
 DumpTObjTable("End of TPC/ITS track fill");


//================== For Jet v2
 TObjArray *selectedTracksAssociated_TPC = new TObjArray; selectedTracksAssociated_TPC->SetOwner(kTRUE);
 selectedTracksAssociated_TPC = GetAcceptedTracksLeading(fEvent,kFALSE,selectedTracksAssociated_TPC);
 TObjArray *selected_TPC_Pairs = new TObjArray; selected_TPC_Pairs->SetOwner(kTRUE);
//================= 

 if(!fprimFMD){
  FillCorrelationTracks(lCentrality,selectedTracksLeading,selectedTracksAssociated,selectedTracksAssociated_TPC,fHistTriggerTrack,fHistReconstTrack,bSign,0, selected_TPC_Pairs);
  FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),selectedTracksLeading,selectedTracksAssociated,selectedTracksAssociated_TPC,fHistTriggerTrackMix,fHistReconstTrackMix,bSign,0,selected_TPC_Pairs);
  DumpTObjTable("End of fill  Correlation");
 }
 else if(fprimTPC)
 {
  FillCorrelationTracks(lCentrality,selectedTracksMC1,selectedTracksMC2,selectedTracksMC1,fHistTriggerTrack,fHistReconstTrack,bSign,0, selected_TPC_Pairs);
  FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),selectedTracksMC1,selectedTracksMC2,selectedTracksMC1,fHistTriggerTrackMix,fHistReconstTrackMix,bSign,0,selected_TPC_Pairs);
  DumpTObjTable("End of fill  Correlation");
 }
 else
 {
  FillCorrelationTracks(lCentrality,selectedTracksLeading,selectedTracksMC2,selectedTracksAssociated_TPC,fHistTriggerTrack,fHistReconstTrack,bSign,0, selected_TPC_Pairs);
  FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),selectedTracksLeading,selectedTracksMC2,selectedTracksAssociated_TPC,fHistTriggerTrackMix,fHistReconstTrackMix,bSign,0,selected_TPC_Pairs);
  DumpTObjTable("End of fill  Correlation");
 }

 selected_TPC_Pairs->Clear();
 delete selected_TPC_Pairs;
 selectedTracksLeading->Clear();
 delete selectedTracksLeading;
 selectedTracksAssociated->Clear();
 delete selectedTracksAssociated;
 selectedTracksAssociated_TPC->Clear();
 delete selectedTracksAssociated_TPC;
 selectedTracksMC1->Clear();
 delete selectedTracksMC1;
 selectedTracksMC2->Clear();
 delete selectedTracksMC2; 
 //DumpTObjTable("after delete TObjects");
 fNEntries++;
}

TObjArray* AliAnalysisTaskSEpPbCorrelationsJetV2::GetFMDhitsYS(Bool_t Aside){
  TObjArray *tracks1 = new TObjArray;
    tracks1->SetOwner(kTRUE);
    AliAODForwardMult* aodForward =static_cast<AliAODForwardMult*>(fEvent->FindListObject("Forward"));
    // Shape of d2Ndetadphi: 200, -4, 6, q20, 0, 2pi
    const TH2D& d2Ndetadphi = aodForward->GetHistogram();
    Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
    Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
    //AliAnalysisTaskValidation::Tracks ret_vector;
    // FMD has no pt resolution!
    Float_t pt = 0;
    for (Int_t iEta = 1; iEta <= nEta; iEta++) {
      Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
      if (!valid) {
	// No data expected for this eta
	continue;
        }
      Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
      for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
	// Bin content is most likely number of particles!
	Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
	if (mostProbableN > 0) {
	  Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
	  //ret_vector.push_back(AliAnalysisTaskValidation::Track(eta, phi, pt, mostProbableN));
	  if(Aside){
	    if(eta<0) continue;
	  } else{
	    if(eta>0) continue;
	  }
	  tracks1->Add(new AliAssociatedVZEROYS(mostProbableN,eta,phi,0,0,0));
	  Double_t cont[3]={eta,phi,lCentrality};
	  fhistfmd->Fill(cont,0,mostProbableN);
	  fh2_FMD_eta_phi->Fill(eta,phi,mostProbableN);
	}
      }
    }    
   return tracks1;
  }



TObjArray *AliAnalysisTaskSEpPbCorrelationsJetV2::GetAcceptedTracksLeading(AliAODEvent *fAOD,Bool_t leading,TObjArray*tracks) {
  //TObjArray *tracks = new TObjArray;
  //tracks->SetOwner(kTRUE);
  Int_t nTracks = fAOD->GetNumberOfTracks();
  Double_t pidqa[4];
  for (Int_t i = 0; i < nTracks; i++) {
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i));
    if (!aodTrack)      continue;
    if (!IsAcceptedTrack(aodTrack))      continue;
    if (aodTrack->Charge() == 0)      continue;
    if(leading){
     pidqa[0]=aodTrack->Pt();
     pidqa[1]=aodTrack->Eta();
     pidqa[2]=RangePhi(aodTrack->Phi());
     pidqa[3]=lCentrality;
     //fHistLeadQA->Fill(pidqa,0);
    }
    Int_t SpAsso=0;
    tracks->Add(new AliAssociatedTrackYS(aodTrack->Charge(), aodTrack->Eta(), aodTrack->Phi(), aodTrack->Pt(), aodTrack->GetID(), -999, -999, SpAsso, 1));
  }
  return tracks;
}

Bool_t AliAnalysisTaskSEpPbCorrelationsJetV2::IsAcceptedTrack(const AliAODTrack *aodTrack) {
  if (!aodTrack)
    return kFALSE;
  //  tight DCA cut
  if(fcollisiontype=="PbPb"){
    if (!aodTrack->TestFilterBit(ffilterbit))  return kFALSE; 
    Float_t nCrossedRowsTPC =aodTrack->GetTPCClusterInfo(2,1);
    if (nCrossedRowsTPC < fnoClusters) return kFALSE;
    if(fCutChargedDCAzMax > 0. ||fCutChargedDCAxyMax > 0.){
      Double_t dTrackXYZ[3]  = {0.,0.,0.};
      Double_t dVertexXYZ[3] = {0.,0.,0.};
      Double_t dDCAXYZ[3]    = {0.,0.,0.};
      aodTrack->GetXYZ(dTrackXYZ);
      lPrimaryBestVtx->GetXYZ(dVertexXYZ);
      for(Short_t i(0); i < 3; i++) { dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i]; }
      
      if(fCutChargedDCAzMax > 0. && TMath::Abs(dDCAXYZ[2]) >fCutChargedDCAzMax) return kFALSE;
      if(fCutChargedDCAxyMax > 0. && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]) > fCutChargedDCAxyMax) return kFALSE;
    }
  }else{
    //if (!aodTrack->TestFilterMask(BIT(ffilterbit)))  return kFALSE; 
    if (!aodTrack->TestFilterBit(ffilterbit))  return kFALSE; 
  }
  if (aodTrack->Pt() < fPtMin) return kFALSE;
  if (aodTrack->Pt() > fPtMax) return kFALSE;  
  if (TMath::Abs(aodTrack->Eta()) > fEtaMax) return kFALSE;

  return kTRUE;
}

void AliAnalysisTaskSEpPbCorrelationsJetV2::FillCorrelationTracks( Double_t centrality, TObjArray *triggerArray, TObjArray *selectedTrackArray, TObjArray *selectedTrackArray_TPC, AliTHn *triggerHist, AliTHn *associateHist, Float_t bSign, Int_t step, TObjArray *selected_TPC_Pairs)
{
 bSign=0;// default
 step=1;//default
 
 if (!triggerHist || !associateHist)    return;

//========= For TPC-TPC
 Double_t binscontTrigTPCTPC[3];
 Double_t binscontTPCTPC[5] = {0.};

//========== For TPC-TPC-FMD
 //Double_t binscontTrig[6];
 Double_t binscontTrig[5];
 //Double_t binscont[7] = {0.}; 
 Double_t binscont[6] = {0.}; 

//========== For TPC-FMDA or TPC-FMDC
 Double_t binscontTrigTPCFMD[2];
 Double_t binscontTPCFMD[4] = {0.};
 
//========= For FMDA-FMDC
 Double_t binscontTrigFMDA[1];
 Double_t binscontFMDAFMDC[3] = {0.};

 if(fAnaMode == "FMDAFMDC")
 {
  for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
      AliAssociatedTrackYS* trigger = (AliAssociatedTrackYS*) triggerArray->At(i);
      if(!trigger)continue;
      Double_t  triggerMultiplicity= trigger->Multiplicity();
      Double_t  triggerEta  = trigger->Eta();
      Double_t  triggerPhi  = trigger->Phi();
      binscontTrigFMDA[0]=fPrimaryZVtx;
      triggerHist->Fill(binscontTrigFMDA,0,(Double_t)triggerMultiplicity);
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
        AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) selectedTrackArray->At(j);
        if(!associate)continue;
        Double_t assophi=associate->Phi();
        Double_t assoeta=associate->Eta();
        Double_t mult=associate->Multiplicity()*triggerMultiplicity;
        binscontFMDAFMDC[0] = fPrimaryZVtx;
        binscontFMDAFMDC[1] = triggerEta-associate->Eta();
        binscontFMDAFMDC[2] = RangePhi_FMD(triggerPhi-associate->Phi());
        //cout << binscontFMDAFMDC[2] << endl;
        if(triggerPhi==assophi && triggerEta==assoeta) continue;
        associateHist->Fill(binscontFMDAFMDC,0,(Double_t)mult);
      }
    }
 }
 else
 { 
 for(Int_t i = 0; i < triggerArray->GetEntriesFast(); i++)
 {  
  AliAssociatedTrackYS *trigger = (AliAssociatedTrackYS *)triggerArray->At(i);
  if (!trigger)    continue;
  Int_t trigID = trigger->GetID();
  Double_t triggerPt = trigger->Pt();
  Double_t triggerEta = trigger->Eta();
  Double_t triggerPhi = trigger->Phi();

  binscontTrigTPCTPC[0] = triggerPt;
  binscontTrigTPCTPC[1] = fPrimaryZVtx;
  binscontTrigTPCTPC[2] = rand()%12 + 0.5;

  binscontTrigTPCFMD[0] = triggerPt;
  binscontTrigTPCFMD[1] = fPrimaryZVtx;

  if(fAnaMode=="TPCTPC") triggerHist->Fill(binscontTrigTPCTPC, 0);
  if(fAnaMode=="TPCFMDA"||fAnaMode=="TPCFMDC") triggerHist->Fill(binscontTrigTPCFMD, 0);
  
  if(fAnaMode=="TPCTPCFMDA" || fAnaMode=="TPCTPCFMDC"||fAnaMode=="TPCTPC")
  {
   for (Int_t j = 0; j < selectedTrackArray_TPC->GetEntriesFast(); j++) 
   {
    AliAssociatedTrackYS *associate_TPC =   (AliAssociatedTrackYS*)selectedTrackArray_TPC->At(j);
    if (!associate_TPC)        continue;
    if (associate_TPC->Pt()<0.5) continue;
    if (triggerPt < associate_TPC->Pt())          continue;
    if (trigID == associate_TPC->GetID())          continue;
    if (fLikeSign)
    {
     if((trigger->Charge())*(associate_TPC->Charge())<0) continue;
    }
    Double_t dTPC_Pairs_Eta = triggerEta-associate_TPC->Eta();
    Double_t dTPC_Pairs_phi = RangePhi(triggerPhi-associate_TPC->Phi());

    binscontTrig[0] = triggerPt;
    //binscontTrig[1] = associate_TPC->Pt();
    binscontTrig[1] = dTPC_Pairs_Eta;
    binscontTrig[2] = dTPC_Pairs_phi;
    binscontTrig[3] = fPrimaryZVtx;

    if(fAnaMode=="TPCTPC")
    {
     binscontTPCTPC[0] = dTPC_Pairs_Eta;
     binscontTPCTPC[1] = associate_TPC->Pt();
     binscontTPCTPC[2] = triggerPt;
     binscontTPCTPC[3] = dTPC_Pairs_phi;
     binscontTPCTPC[4] = fPrimaryZVtx;

     associateHist->Fill(binscontTPCTPC, 0);
    }
	
 //   cout << Form("f2pc_double_%d_%d",fh2_pt_trig_asso->GetXaxis()->FindBin(triggerPt)-1,fh2_pt_trig_asso->GetYaxis()->FindBin(associate_TPC->Pt())-1) << endl;
 
    if(fAnaMode=="TPCTPCFMDA" || fAnaMode=="TPCTPCFMDC")
    {
     if(fReduceDphi>0.) 
     {
      if((dTPC_Pairs_phi<-1*fReduceDphi)||(dTPC_Pairs_phi>fReduceDphi)) continue;
     }
     triggerHist->Fill(binscontTrig, 0);
    }

    selected_TPC_Pairs->Add(new AliTrigAssoPairST(trigger->Charge(), triggerEta, triggerPhi, triggerPt, associate_TPC->Pt(), trigger->GetID(), -999, -999, 0, 1, dTPC_Pairs_Eta,dTPC_Pairs_phi));
   }
  }
  if(fAnaMode == "TPCFMDA" || fAnaMode == "TPCFMDC")
  {
   for(Int_t k = 0; k < selectedTrackArray->GetEntriesFast(); k++)
   {
    AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) selectedTrackArray->At(k);
    binscontTPCFMD[0] = triggerPt; // TPC pt  
    binscontTPCFMD[1] = fPrimaryZVtx; // Vz  
    binscontTPCFMD[2] = triggerEta - associate->Eta(); // TPC-FMD eta1-eta2   
    binscontTPCFMD[3] = RangePhi(triggerPhi - associate->Phi());; // TPC-FMD phi1-phi2   
    associateHist->Fill(binscontTPCFMD, 0, (Double_t)associate->Multiplicity());
   }
  }
 }

 if(fAnaMode=="TPCTPCFMDA" || fAnaMode=="TPCTPCFMDC")
 {
  for(Int_t j2 = 0;j2 < selected_TPC_Pairs->GetEntriesFast(); j2++)
  {
   AliTrigAssoPairST *associate_TPC_Pair = (AliTrigAssoPairST*)selected_TPC_Pairs->At(j2);
   if(!associate_TPC_Pair) continue;
   if(associate_TPC_Pair->Pt_Asso()< fAsscoptCut) continue;
   for(Int_t k = 0; k < selectedTrackArray->GetEntriesFast(); k++)
   {
    binscont[0] = associate_TPC_Pair->Getdeta_pairs(); // TPC eta1-eta2  
    //binscont[1] = associate_TPC_Pair->Pt_Asso(); // associate TPC pt  
    binscont[1] = associate_TPC_Pair->Pt();; // trigger TPC pt
    binscont[2] = fPrimaryZVtx; //Vz
    AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) selectedTrackArray->At(k);
    if(!associate) continue;
    Double_t dFMD_Eta = associate->Eta();
    binscont[3] = associate_TPC_Pair->Eta() - dFMD_Eta; // eta1-eta3
    binscont[4] = RangePhi(associate_TPC_Pair->Phi()-associate->Phi()); // phi1-phi3
    binscont[5] = associate_TPC_Pair->Getdphi_pairs(); // TPC phi1-phi2

    associateHist->Fill(binscont, 0, (Double_t)associate->Multiplicity()); 
   }
  }
  }
 }
}

void AliAnalysisTaskSEpPbCorrelationsJetV2::FillCorrelationTracksMixing(Double_t centrality, Double_t pvxMix, TObjArray *triggerArray, TObjArray *selectedTrackArray, TObjArray *selectedTrackArray_TPC, AliTHn *triggerHist, AliTHn *associateHist, Float_t bSign, Int_t step, TObjArray *selected_TPC_Pairs)
{
  Float_t bSign1=bSign;
  Int_t step1=step;

  AliEventPool *pool = fPoolMgr->GetEventPool(centrality, pvxMix);
  if (!pool){
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality,
                  pvxMix));
  }
  //cout<< "event"<<pool->GetCurrentNEvents() << endl;
  if (pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
      Double_t binscontTrig[2];
      mixedDist ->Fill(centrality, pool->NTracksInPool());
      mixedDist2->Fill(centrality, pool->GetCurrentNEvents());      
      Int_t nMix = pool->GetCurrentNEvents();

      if(fAnaMode=="TPCTPCFMDA" || fAnaMode=="TPCTPCFMDC")
      {
       //Double_t binscont[7];
       Double_t binscont[6];
       for(Int_t j2 = 0; j2 < selected_TPC_Pairs->GetEntriesFast(); j2++)
       {
        AliTrigAssoPairST *associate_TPC_Pair = (AliTrigAssoPairST*)selected_TPC_Pairs->At(j2);
        if (!associate_TPC_Pair) continue;
        if (associate_TPC_Pair->Pt_Asso()< fAsscoptCut) continue;

        binscontTrig[0] = associate_TPC_Pair->Pt();
        binscontTrig[1] = rand()%12 + 0.5; 
        binscont[0] = associate_TPC_Pair->Getdeta_pairs(); // TPC eta1-eta2  
//        binscont[1] = associate_TPC_Pair->Pt_Asso(); // associate TPC pt  
        binscont[1] = associate_TPC_Pair->Pt(); // trigger TPC pt
        binscont[2] = pvxMix; //Vz

        for (Int_t jMix = 0; jMix < nMix; jMix++) {
         TObjArray *mixEvents = pool->GetEvent(jMix);
         for(Int_t k = 0; k < mixEvents->GetEntriesFast(); k++)
         {
          AliAssociatedTrackYS* associate=(AliAssociatedTrackYS*)  mixEvents->At(k);
          if(!associate) continue;
          Double_t dFMD_Eta = associate->Eta();
          binscont[3] = associate_TPC_Pair->Eta() - dFMD_Eta; // eta1-eta3
          binscont[4] = RangePhi(associate_TPC_Pair->Phi()-associate->Phi()); // phi1-phi3
          binscont[5] = associate_TPC_Pair->Getdphi_pairs(); // TPC phi1-phi2
          associateHist->Fill(binscont, 0,(Double_t)associate->Multiplicity()/(Double_t)nMix);
         }
        }
       }
      }

     if(fAnaMode=="TPCFMDA" || fAnaMode=="TPCFMDC")
     {
      Double_t binscont[4];
      for(Int_t j2 = 0; j2 < triggerArray->GetEntriesFast(); ++j2)
      {
       AliAssociatedTrackYS *trig = (AliAssociatedTrackYS *)triggerArray->At(j2);
       if (!trig)          continue;
       Double_t triggerPhi = trig->Phi();
       Double_t triggerEta = trig->Eta();
       Double_t triggerPt = trig->Pt();
      
       binscontTrig[0] = triggerPhi;
       binscontTrig[1] = pvxMix;

        for (Int_t jMix = 0; jMix < nMix; jMix++) {
         TObjArray *mixEvents = pool->GetEvent(jMix);
         for(Int_t k = 0; k < mixEvents->GetEntriesFast(); k++)
         {
          AliAssociatedTrackYS* associate=(AliAssociatedTrackYS*)  mixEvents->At(k);
          if(!associate) continue;
          binscont[0] = triggerPt; // Trigger pT
          binscont[1] = pvxMix;    // Vz
          binscont[2] = triggerEta - associate->Eta(); // TPC-FMD eta1-eta2
          binscont[3] = RangePhi(triggerPhi - associate->Phi()); // TPC-FMD phi1-phi2
          associateHist->Fill(binscont, 0,(Double_t)associate->Multiplicity()/(Double_t)nMix);
         }
        }
       }
      }
 
      if(fAnaMode=="TPCTPC")
      {
       Double_t binscont[5];
       for(Int_t j2 = 0; j2 < triggerArray->GetEntriesFast(); j2++)
       {
        AliAssociatedTrackYS *trig = (AliAssociatedTrackYS *)triggerArray->At(j2);
        if (!trig)          continue;
        Double_t triggerPhi = trig->Phi();
        Double_t triggerEta = trig->Eta();
        Double_t triggerPt = trig->Pt();
        for(Int_t jMix = 0; jMix < nMix; jMix++)
        {
         TObjArray *mixEvents = pool->GetEvent(jMix);
         for(Int_t k = 0; k < mixEvents->GetEntriesFast(); k++)
         {
          AliAssociatedTrackYS* associate=(AliAssociatedTrackYS*)  mixEvents->At(k);
          if(!associate) continue;
          if(fLikeSign)
          {
           if((trig->Charge())*(associate->Charge())<0) continue;
          }
          binscont[0] = triggerEta - associate->Eta();              // TPC-TPC dPhi
          binscont[1] = associate->Pt();
          binscont[2] = triggerPt;                                  
          binscont[3] = RangePhi(triggerPhi - associate->Phi());    // TPC-TPC dPhi
          binscont[4] = pvxMix;                                     // Vz
          associateHist->Fill(binscont, 0,1./(Double_t)nMix);
         }
        }
       } 
      }
     if(fAnaMode == "FMDAFMDC")
     {
      Double_t binscont[3]; 
      for(Int_t j2=0;j2<triggerArray->GetEntriesFast();j2++)
      {
        AliAssociatedTrackYS* trigger = (AliAssociatedTrackYS*) triggerArray->At(j2);
        if(!trigger)continue;
        Float_t triggerMultiplicity= trigger->Multiplicity();
        Float_t  triggerEta  = trigger->Eta();
        Float_t  triggerPhi  = trigger->Phi();
        for (Int_t jMix = 0; jMix < nMix; jMix++){
          TObjArray *mixEvents = pool->GetEvent(jMix);
          for(Int_t k = 0; k < mixEvents->GetEntriesFast(); k++)
          {     
           AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) mixEvents->At(k);
           if(!associate)continue;
           //          Double_t associatemultiplicity=associate->Multiplicity();
           Float_t assophi=associate->Phi();
           Float_t assoeta=associate->Eta();
           Float_t mult=triggerMultiplicity*associate->Multiplicity();
           binscont[0] = pvxMix;
           binscont[1] = triggerEta-associate->Eta();
           binscont[2] = RangePhi_FMD(triggerPhi-associate->Phi());
           associateHist->Fill(binscont,0,mult/(Double_t)nMix); 
          }
        } 
      }  
     }
  }
  TObjArray* tracksClone=CloneTrack(selectedTrackArray);
  pool->UpdatePool(tracksClone);
}


TObjArray* AliAnalysisTaskSEpPbCorrelationsJetV2::CloneTrack(TObjArray*selectedTrackArray){
  TObjArray *tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i = 0; i < selectedTrackArray->GetEntriesFast(); i++) {
    AliAssociatedTrackYS *particle =  (AliAssociatedTrackYS *)selectedTrackArray->At(i);
    tracksClone->Add(new AliAssociatedTrackYS(particle->Charge(), particle->Eta(), particle->Phi(), particle->Pt(),
					      particle->GetID(), particle->GetIDFirstDaughter(),
					      particle->GetIDSecondDaughter(), particle->WhichCandidate(),
					      particle->Multiplicity()));
  }
  
  return tracksClone;
}

Double_t AliAnalysisTaskSEpPbCorrelationsJetV2::RangePhi(Double_t DPhi) {
  if (DPhi < -TMath::Pi() / 2)   DPhi += 2 * TMath::Pi();
  if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
  return DPhi;
}

Double_t AliAnalysisTaskSEpPbCorrelationsJetV2::RangePhi_FMD(Double_t DPhi) {
  //if (DPhi < (-TMath::Pi() / 2 -0.0001))   DPhi += 2 * TMath::Pi();
  //if (DPhi > (3 * TMath::Pi() / 2-0.0001)) DPhi -= 2*TMath::Pi();
  DPhi = TMath::ATan2(TMath::Sin(DPhi), TMath::Cos(DPhi));
  if (DPhi < (-0.5*TMath::Pi()-0.0001))    DPhi += 2 * TMath::Pi();
  //if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
  return DPhi;
}

Double_t AliAnalysisTaskSEpPbCorrelationsJetV2::RangePhi2(Double_t DPhi) {
  DPhi = TMath::ATan2(TMath::Sin(DPhi), TMath::Cos(DPhi));
  if (DPhi < -1.178097)    DPhi += 2 * TMath::Pi();
  return DPhi;
}



void AliAnalysisTaskSEpPbCorrelationsJetV2::DumpTObjTable(const char* note)
{
  if(note) {
    //    printf("TObjectTable::%s",note);
  }
  //  gObjectTable->Print();
}

 

