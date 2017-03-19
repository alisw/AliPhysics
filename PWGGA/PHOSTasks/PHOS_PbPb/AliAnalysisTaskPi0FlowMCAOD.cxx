/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// Extension to Pi0FLOw, mimicing AliPHOSHijingEfficiency
// by Dmitri Peressounko, 05.02.2013
// Authors: Henrik Qvigstad, Dmitri Peressounko
// Date   : 05.04.2013
// Adopted for AOD analysis by Paul Batzing and Boris Polishchuk (10.03.2014).
/* $Id$ */

#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "THashList.h"
#include "TArray.h"
#include "TArrayD.h"


#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliPHOSHijingEfficiency.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSCalibData.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliCDBManager.h"
#include "AliCentrality.h" 
#include "AliEventplane.h"
#include "TProfile.h"
#include <TPDGCode.h>
#include "AliOADBContainer.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisTaskPi0FlowMCAOD.h"
ClassImp(AliAnalysisTaskPi0FlowMCAOD);

const Double_t AliAnalysisTaskPi0FlowMCAOD::kRCut = 1.;

AliAnalysisTaskPi0FlowMCAOD::AliAnalysisTaskPi0FlowMCAOD(const char* name, AliAnalysisTaskPi0Flow::Period period)
: AliAnalysisTaskPi0Flow(name, period),fMcArray(0x0),kOffVertexCutSet(kTRUE)
{
}

AliAnalysisTaskPi0FlowMCAOD::~AliAnalysisTaskPi0FlowMCAOD()
{
}

void AliAnalysisTaskPi0FlowMCAOD::UserCreateOutputObjects()
{
  // Do standard Pi0Flow CreateOuputObjects
  AliAnalysisTaskPi0Flow::UserCreateOutputObjects();
  
  // MC Generated histograms
  char key[55];
  for(Int_t cent=0; cent < fCentEdges.GetSize()-1; cent++){
    snprintf(key,55,"hMC_rap_gamma_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity pi0",200,-1.,1.)) ;
    snprintf(key,55,"hMC_rap_pi0_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity pi0",200,-1.,1.)) ;
    snprintf(key,55,"hMC_rap_eta_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity eta",200,-1.,1.)) ;
    snprintf(key,55,"hMC_phi_gamma_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Phi pi0",200,0.,TMath::TwoPi())) ;
    snprintf(key,55,"hMC_phi_pi0_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Phi pi0",200,0.,TMath::TwoPi())) ;
    snprintf(key,55,"hMC_phi_eta_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Phi eta",200,0.,TMath::TwoPi())) ;
    snprintf(key,55,"hMC_all_gamma_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity photon",250,0.,25.)) ;
    snprintf(key,55,"hMC_all_pi0_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity pi0",250,0.,25.)) ;
    snprintf(key,55,"hMC_all_eta_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity eta",250,0.,25.)) ;
    snprintf(key,55,"hMC_unitEta_gamma_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Pt photon",250,0.,25.)) ;
    snprintf(key,55,"hMC_unitEta_pi0_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity eta",250,0.,25.)) ;
    snprintf(key,55,"hMC_unitEta_eta_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity eta",250,0.,25.)) ;

    snprintf(key,55,"hMC_all_K0S_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity photon",250,0.,25.)) ;

    snprintf(key,55,"hMC_unitEta_K0S_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity eta",250,0.,25.)) ;

  }
  fOutputContainer->Add(new TH2F("hMC_gamma_vertex","Creation vertex",25,0.,25.,1000,0.,500.)) ;
  fOutputContainer->Add(new TH2F("hMC_pi0_vertex","Creation vertex",25,0.,25.,1000,0.,500.)) ;
  fOutputContainer->Add(new TH2F("hMC_eta_vertex","Creation vertex",25,0.,25.,1000,0.,500.)) ;
 
 
  Int_t nPt      = 200;
  Double_t ptMin = 0;
  Double_t ptMax = 20; 
  fOutputContainer->Add(new TH2F("Vertex","Pi0 creation vertex",nPt,ptMin,ptMax,5000,0.,500.));
  fOutputContainer->Add(new TH3F("hSecondPi0RphiZ","Secondary pi0 vertex",450,0.,450.,100,0.,TMath::TwoPi(),200,-100.,100.));
  fOutputContainer->Add(new TH2F("hSecondPi0RE","Secondary pi0 vertex",450,0.,450.,200,0.,20.));
  fOutputContainer->Add(new TH3F("hMass_R","Mass vs radius any parent",50,0.,0.25,100,0.,10.,300,0.,600.));
  fOutputContainer->Add(new TH3F("Real_pi_R","All clusters",50,0.,0.25,100,0.,10.,250,0.,500.));
  fOutputContainer->Add(new TH3F("Real_pi_Z","All clusters",50,0.,0.25,100,0.,10.,100,-100.,100.));
//  fOutputContainer->Add(new TH2F(Form("Real_npnp_RZ"),"All clusters",250,0.,500.,100,-100.,100.));
//  fOutputContainer->Add(new TH3F(Form("Real_mass_R"),"All clusters",50,0.,0.25,100,0.,10.,300,0.,600.));

  const Int_t nM       = 500;
  const Double_t mMin  = 0.0;
  const Double_t mMax  = 1.0;

  for(Int_t cen=0; cen < fCentEdges.GetSize()-1; cen++){
    fOutputContainer->Add(new TH1F(Form("hPrimPhot_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimEl_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimPi0_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimEta_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimPipm_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimP_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimPbar_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimN_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimNbar_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimK0S_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimK0L_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimKpm_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimOther_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));

    //pairs from common parents
    fOutputContainer->Add(new TH2F(Form("hParentAll_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentK0s_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentGamma_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentEl_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentOther_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentDirPi0_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
   
    //common parent - pi0
    fOutputContainer->Add(new TH2F(Form("hParentPi0NoPrim_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Eta_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Omega_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Pipm_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Kpm_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Ks_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Kl_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0pn_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0antipn_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    
  }
  
  
  //Photon contaminations
  fOutputContainer->Add(new TH2F("hPipmGammaConvR","Conversion radius" ,200,0.,20.,1000,0.,500.));
  fOutputContainer->Add(new TH2F("hPipmElConvR","Conversion radius" ,200,0.,20.,1000,0.,500.));
  fOutputContainer->Add(new TH2F("hPipmNConvR","Conversion radius" ,200,0.,20.,1000,0.,500.));
  fOutputContainer->Add(new TH2F("hPipmOtherConvR","Conversion radius" ,200,0.,20.,1000,0.,500.));
  fOutputContainer->Add(new TH2F("hPipmGammaConvRZ","Conversion radius" ,400,-200.,200.,1000,0.,500.)); 
   
   const Int_t nTypes=24 ;
   char partTypes[nTypes][55] ;
   snprintf(partTypes[0],55,"hGammaNoPrim") ; //
   snprintf(partTypes[1],55,"hGammaPhot") ; //
   snprintf(partTypes[2],55,"hGammaEl") ; //
   snprintf(partTypes[3],55,"hGammaPi0") ; //
   snprintf(partTypes[4],55,"hGammaEta") ; //
   snprintf(partTypes[5],55,"hhGammaOmega") ; //
   snprintf(partTypes[6],55,"hGammaPipm") ; //
   snprintf(partTypes[7],55,"hGammaP") ; //
   snprintf(partTypes[8],55,"hGammaPbar") ; //
   snprintf(partTypes[9],55,"hGammaN") ; //
   snprintf(partTypes[10],55,"hGammaNbar") ; //
   snprintf(partTypes[11],55,"hGammaK0S") ; //
   snprintf(partTypes[12],55,"hGammaK0L") ; //
   snprintf(partTypes[13],55,"hGammaKpm") ; //
   snprintf(partTypes[14],55,"hGammaKstar") ; //
   snprintf(partTypes[15],55,"hGammaDelta") ; //
   snprintf(partTypes[16],55,"hGammaOtherCharged") ; //
   snprintf(partTypes[17],55,"hGammaOtherNeutral") ; //
   snprintf(partTypes[18],55,"hGammaPipmGamma") ; //
   snprintf(partTypes[19],55,"hGammaPipmEl") ; //
   snprintf(partTypes[20],55,"hGammaPipmOther") ; //
   snprintf(partTypes[21],55,"hGammaPipmDirect") ; //
   snprintf(partTypes[22],55,"hGammaPipmp") ; //
   snprintf(partTypes[23],55,"hGammaPipmn") ; //
 
   const Int_t nPID=12 ;
   char cPID[12][25] ;
   snprintf(cPID[0],25,"All") ;
   snprintf(cPID[1],25,"Allcore") ;
   snprintf(cPID[2],25,"CPV") ;
   snprintf(cPID[3],25,"CPVcore") ;
   snprintf(cPID[4],25,"CPV2") ;
   snprintf(cPID[5],25,"CPV2core") ;
   snprintf(cPID[6],25,"Disp") ;
   snprintf(cPID[7],25,"Dispcore") ;
   snprintf(cPID[8],25,"Disp2") ;
   snprintf(cPID[9],25,"Disp2core") ;
   snprintf(cPID[10],25,"Both") ;
   snprintf(cPID[11],25,"Bothcore") ;
 
   for(Int_t itype=0; itype<nTypes; itype++){
     for(Int_t iPID=0; iPID<nPID; iPID++){
       for(Int_t cen=0; cen < fCentEdges.GetSize()-1; cen++){
         fOutputContainer->Add(new TH1F(Form("%s_%s_cen%d",partTypes[itype],cPID[iPID],cen),"Cluster parents",nPt,ptMin,ptMax));
       }
     }
   }

  PostData(1, fOutputContainer);
}

void AliAnalysisTaskPi0FlowMCAOD::UserExec(Option_t* option)
{
  fMcArray = GetMCArray();
  if(!fMcArray) { PostData(1, fOutputContainer); return; }
 
  AliAnalysisTaskPi0Flow::UserExec(option);
}


void AliAnalysisTaskPi0FlowMCAOD::SelectPhotonClusters()
{
  AliAnalysisTaskPi0Flow::SelectPhotonClusters();
  
  for (Int_t i1=0; i1<fCaloPhotonsPHOS->GetEntriesFast(); i1++) {
    AliCaloPhoton * photon = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1);
    AliVCluster* cluster = photon->GetCluster();
    Bool_t sure=  kTRUE;
    Int_t primary=FindPrimary(cluster,sure) ;
    photon->SetPrimary(primary);
    photon->SetWeight(PrimaryWeight(primary)) ;
  }

    if(kOffVertexCutSet) {
        for (Int_t i1=0; i1<fCaloPhotonsPHOS->GetEntriesFast(); i1++) {
            AliCaloPhoton * photon = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1);
            Int_t primary = photon->GetPrimary();
            AliAODMCParticle* p = GetParticle(primary);
            if(R(p) >kRCut) {
                if(p->PdgCode()==11 || p->PdgCode()==-11) continue;
                else { fCaloPhotonsPHOS->Remove(photon); fCaloPhotonsPHOS->Compress(); }
            }
        }
    }
    
}

void AliAnalysisTaskPi0FlowMCAOD::FillSelectedClusterHistograms()
{
  for (Int_t i1=0; i1<fCaloPhotonsPHOS->GetEntriesFast(); i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;

    Double_t dphiA=ph1->Phi()-fRPV0A ;
    while(dphiA<0)dphiA+=TMath::Pi() ;
    while(dphiA>TMath::Pi())dphiA-=TMath::Pi() ;

    Double_t dphiC=ph1->Phi()-fRPV0C ;
    while(dphiC<0)dphiC+=TMath::Pi() ;
    while(dphiC>TMath::Pi())dphiC-=TMath::Pi() ;

    Double_t dphiT=ph1->Phi()-fRP ;
    while(dphiT<0)dphiT+=TMath::Pi() ;
    while(dphiT>TMath::Pi())dphiT-=TMath::Pi() ;
    
    Double_t pt = ph1->Pt() ;
    Double_t ptcore = ph1->GetMomV2()->Pt() ;
    Double_t w = ph1->GetWeight();

    FillHistogram(Form("hPhotPhiV0AAll_cen%d",fCentBin),pt,dphiA, w) ;
    FillHistogram(Form("hPhotPhiV0CAll_cen%d",fCentBin),pt,dphiC, w) ;
    if(fHaveTPCRP)
      FillHistogram(Form("hPhotPhiTPCAll_cen%d",fCentBin),pt,dphiT, w) ;
    FillHistogram(Form("hPhotPhiV0AAllcore_cen%d",fCentBin),ptcore,dphiA, w) ;
    FillHistogram(Form("hPhotPhiV0CAllcore_cen%d",fCentBin),ptcore,dphiC, w) ;
    if(fHaveTPCRP)
      FillHistogram(Form("hPhotPhiTPCAllcore_cen%d",fCentBin),ptcore,dphiT, w) ;

    FillHistogram(Form("hPhotAll_cen%d",fCentBin),pt, w) ;
    FillHistogram(Form("hPhotAllcore_cen%d",fCentBin),ptcore, w) ;
    if(ph1->IsntUnfolded()){
      FillHistogram(Form("hPhotAllwou_cen%d",fCentBin),pt, w) ;
      FillHistogram(Form("hPhotPhiV0AAllwou_cen%d",fCentBin),pt,dphiA, w) ;
      FillHistogram(Form("hPhotPhiV0CAllwou_cen%d",fCentBin),pt,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCAllwou_cen%d",fCentBin),pt,dphiT, w) ;
    }
    if(ph1->IsCPVOK()){
      FillHistogram(Form("hPhotPhiV0ACPV_cen%d",fCentBin),pt,dphiA, w) ;
      FillHistogram(Form("hPhotPhiV0CCPV_cen%d",fCentBin),pt,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCCPV_cen%d",fCentBin),pt,dphiT, w) ;

      FillHistogram(Form("hPhotPhiV0ACPVcore_cen%d",fCentBin),ptcore,dphiA, w) ;
      FillHistogram(Form("hPhotPhiV0CCPVcore_cen%d",fCentBin),ptcore,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCCPVcore_cen%d",fCentBin),ptcore,dphiT, w) ;

      FillHistogram(Form("hPhotCPV_cen%d",fCentBin),pt, w) ;
      FillHistogram(Form("hPhotCPVcore_cen%d",fCentBin),ptcore, w) ;
    }
    if(ph1->IsCPV2OK()){
      FillHistogram(Form("hPhotPhiV0ACPV2_cen%d",fCentBin),pt,dphiA, w) ;
      FillHistogram(Form("hPhotPhiV0CCPV2_cen%d",fCentBin),pt,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCCPV2_cen%d",fCentBin),pt,dphiT, w) ;

      FillHistogram(Form("hPhotPhiV0ACPV2core_cen%d",fCentBin),ptcore,dphiA, w) ;
      FillHistogram(Form("hPhotPhiV0CCPV2core_cen%d",fCentBin),ptcore,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCCPV2core_cen%d",fCentBin),ptcore,dphiT, w) ;
      FillHistogram(Form("hPhotCPV2_cen%d",fCentBin),pt, w) ;
      FillHistogram(Form("hPhotCPV2core_cen%d",fCentBin),ptcore, w) ;
    }
    if(ph1->IsDispOK()){
      FillHistogram(Form("hPhotPhiV0ADisp_cen%d",fCentBin),pt,dphiA, w) ;
      FillHistogram(Form("hPhotPhiV0CDisp_cen%d",fCentBin),pt,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCDisp_cen%d",fCentBin),pt,dphiT, w) ;

      FillHistogram(Form("hPhotPhiV0ADispcore_cen%d",fCentBin),ptcore,dphiA, w) ;
      FillHistogram(Form("hPhotPhiV0CDispcore_cen%d",fCentBin),ptcore,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCDispcore_cen%d",fCentBin),ptcore,dphiT, w) ;

      if(ph1->IsntUnfolded()){
        FillHistogram(Form("hPhotPhiV0ADispwou_cen%d",fCentBin),pt,dphiA, w) ;
        FillHistogram(Form("hPhotPhiV0CDispwou_cen%d",fCentBin),pt,dphiC, w) ;
        if(fHaveTPCRP)
          FillHistogram(Form("hPhotPhiTPCDispwou_cen%d",fCentBin),pt,dphiT, w) ;

      }
      FillHistogram(Form("hPhotDisp_cen%d",fCentBin),pt, w) ;
      FillHistogram(Form("hPhotDispcore_cen%d",fCentBin),ptcore, w) ;
      if(ph1->IsntUnfolded()){
        FillHistogram(Form("hPhotDispwou_cen%d",fCentBin),pt, w) ;
      }
      if(ph1->IsCPVOK()){
	FillHistogram(Form("hPhotPhiV0ABoth_cen%d",fCentBin),pt,dphiA, w) ;
	FillHistogram(Form("hPhotPhiV0CBoth_cen%d",fCentBin),pt,dphiC, w) ;
        if(fHaveTPCRP)
  	  FillHistogram(Form("hPhotPhiTPCBoth_cen%d",fCentBin),pt,dphiT, w) ;

	FillHistogram(Form("hPhotPhiV0ABothcore_cen%d",fCentBin),ptcore,dphiA, w) ;
	FillHistogram(Form("hPhotPhiV0CBothcore_cen%d",fCentBin),ptcore,dphiC, w) ;
        if(fHaveTPCRP)
	  FillHistogram(Form("hPhotPhiTPCBothcore_cen%d",fCentBin),ptcore,dphiT, w) ;

	FillHistogram(Form("hPhotBoth_cen%d",fCentBin),pt, w) ;
	FillHistogram(Form("hPhotBothcore_cen%d",fCentBin),ptcore, w) ;
      }
    }
    if(ph1->IsDisp2OK()){
      FillHistogram(Form("hPhotPhiV0ADisp2_cen%d",fCentBin),pt,dphiA, w) ;
      FillHistogram(Form("hPhotPhiV0CDisp2_cen%d",fCentBin),pt,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCDisp2_cen%d",fCentBin),pt,dphiT, w) ;
      FillHistogram(Form("hPhotPhiV0ADisp2core_cen%d",fCentBin),ptcore,dphiA, w) ;
      FillHistogram(Form("hPhotPhiV0CDisp2core_cen%d",fCentBin),ptcore,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCDisp2core_cen%d",fCentBin),ptcore,dphiT, w) ;

      FillHistogram(Form("hPhotDisp2_cen%d",fCentBin),pt, w) ;
      FillHistogram(Form("hPhotDisp2core_cen%d",fCentBin),ptcore, w) ;
      if(ph1->IsCPVOK()){
	FillHistogram(Form("hPhotPhiV0ABoth2_cen%d",fCentBin),pt,dphiA, w) ;
	FillHistogram(Form("hPhotPhiV0CBoth2_cen%d",fCentBin),pt,dphiC, w) ;
        if(fHaveTPCRP)
  	  FillHistogram(Form("hPhotPhiTPCBoth2_cen%d",fCentBin),pt,dphiT, w) ;

	FillHistogram(Form("hPhotPhiV0ABoth2core_cen%d",fCentBin),ptcore,dphiA, w) ;
	FillHistogram(Form("hPhotPhiV0CBoth2core_cen%d",fCentBin),ptcore,dphiC, w) ;
        if(fHaveTPCRP)
	  FillHistogram(Form("hPhotPhiTPCBoth2core_cen%d",fCentBin),ptcore,dphiT, w) ;

	FillHistogram(Form("hPhotBoth2_cen%d",fCentBin),pt, w) ;
	FillHistogram(Form("hPhotBoth2core_cen%d",fCentBin),ptcore, w) ;
      }
    }
  }
}

void AliAnalysisTaskPi0FlowMCAOD::ConsiderPi0s()
{
  char key[55];
  for (Int_t i1=0; i1 < fCaloPhotonsPHOS->GetEntriesFast()-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;
    const Double_t w1 = ph1->GetWeight();
    for (Int_t i2=i1+1; i2<fCaloPhotonsPHOS->GetEntriesFast(); i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i2) ;
      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      
      const Double_t w2 = ph2->GetWeight();
      Double_t w = TMath::Sqrt(w1*w2);
      
      FillHistogram("hPHOSphi",fCentrality,p12.Pt(),p12.Phi(), w) ;
      Double_t dphiA=p12.Phi()-fRPV0A ;
      while(dphiA<0)dphiA+=TMath::Pi() ;
      while(dphiA>TMath::Pi())dphiA-=TMath::Pi() ;

      Double_t dphiC=p12.Phi()-fRPV0C ;
      while(dphiC<0)dphiC+=TMath::Pi() ;
      while(dphiC>TMath::Pi())dphiC-=TMath::Pi() ;

      Double_t dphiT=p12.Phi()-fRP ;
      while(dphiT<0)dphiT+=TMath::Pi() ;
      while(dphiT>TMath::Pi())dphiT-=TMath::Pi() ;

      Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;
      Double_t m=p12.M() ;
      Double_t mcore=pv12.M() ;
      Double_t pt=p12.Pt() ;
      Double_t ptcore=pv12.Pt() ;
      Double_t pt1=ph1->Pt() ;
      Double_t pt2=ph2->Pt() ;
      Double_t ptcore1=ph1->GetMomV2()->Pt() ;
      Double_t ptcore2=ph2->GetMomV2()->Pt() ;

      FillHistogram(Form("hMassPtV0AAll_cen%d",fCentBin),m,pt,dphiA, w) ;
      FillHistogram(Form("hMassPtV0CAll_cen%d",fCentBin),m,pt,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hMassPtTPCAll_cen%d",fCentBin),m,pt,dphiT, w) ;

      FillHistogram(Form("hMassPtV0AAllcore_cen%d",fCentBin),mcore,ptcore,dphiA, w) ;
      FillHistogram(Form("hMassPtV0CAllcore_cen%d",fCentBin),mcore,ptcore,dphiC, w) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hMassPtTPCAllcore_cen%d",fCentBin),mcore,ptcore,dphiT, w) ;


      FillHistogram(Form("hPi0All_cen%d",fCentBin),m,pt, w) ;
      FillHistogram(Form("hPi0Allcore_cen%d",fCentBin),mcore,ptcore, w) ;
      if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
        FillHistogram(Form("hPi0Allwou_cen%d",fCentBin),m,pt, w) ;
        FillHistogram(Form("hMassPtV0AAllwou_cen%d",fCentBin),m,pt,dphiA, w) ;
        FillHistogram(Form("hMassPtV0CAllwou_cen%d",fCentBin),m,pt,dphiC, w) ;
        if(fHaveTPCRP)
          FillHistogram(Form("hMassPtTPCAllwou_cen%d",fCentBin),m,pt,dphiT, w) ;
      }

      FillHistogram(Form("hSingleAll_cen%d",fCentBin),m,pt1, w) ;
      FillHistogram(Form("hSingleAll_cen%d",fCentBin),m,pt2, w) ;
      FillHistogram(Form("hSingleAllcore_cen%d",fCentBin),mcore,ptcore1, w) ;
      FillHistogram(Form("hSingleAllcore_cen%d",fCentBin),mcore,ptcore2, w) ;
      if(ph1->IsntUnfolded())
        FillHistogram(Form("hSingleAllwou_cen%d",fCentBin),m,pt1, w) ;
      if(ph2->IsntUnfolded())
        FillHistogram(Form("hSingleAllwou_cen%d",fCentBin),m,pt2, w) ;
      if(ph1->IsCPVOK()){
        FillHistogram(Form("hSingleCPV_cen%d",fCentBin),m,pt1, w) ;
        FillHistogram(Form("hSingleCPVcore_cen%d",fCentBin),mcore,ptcore1, w) ;
      }
      if(ph2->IsCPVOK()){
        FillHistogram(Form("hSingleCPV_cen%d",fCentBin),m,pt2, w) ;
        FillHistogram(Form("hSingleCPVcore_cen%d",fCentBin),mcore,ptcore2, w) ;
      }
      if(ph1->IsCPV2OK()){
        FillHistogram(Form("hSingleCPV2_cen%d",fCentBin),m,pt1, w) ;
        FillHistogram(Form("hSingleCPV2core_cen%d",fCentBin),mcore,ptcore2, w) ;
      }
      if(ph2->IsCPV2OK()){
        FillHistogram(Form("hSingleCPV2_cen%d",fCentBin),m,pt2, w) ;
        FillHistogram(Form("hSingleCPV2core_cen%d",fCentBin),mcore,ptcore2, w) ;
      }
      if(ph1->IsDispOK()){
        FillHistogram(Form("hSingleDisp_cen%d",fCentBin),m,pt1, w) ;
        if(ph1->IsntUnfolded()){
          FillHistogram(Form("hSingleDispwou_cen%d",fCentBin),m,pt1, w) ;
	}
        FillHistogram(Form("hSingleDispcore_cen%d",fCentBin),mcore,ptcore1, w) ;
      }
      if(ph2->IsDispOK()){
        FillHistogram(Form("hSingleDisp_cen%d",fCentBin),m,pt2, w) ;
        if(ph1->IsntUnfolded()){
          FillHistogram(Form("hSingleDispwou_cen%d",fCentBin),m,pt2, w) ;
	}
        FillHistogram(Form("hSingleDispcore_cen%d",fCentBin),mcore,ptcore2, w) ;
      }
      if(ph1->IsDisp2OK()){
        FillHistogram(Form("hSingleDisp2_cen%d",fCentBin),m,pt1, w) ;
        FillHistogram(Form("hSingleDisp2core_cen%d",fCentBin),mcore,ptcore1, w) ;
      }
      if(ph2->IsDisp2OK()){
        FillHistogram(Form("hSingleDisp2_cen%d",fCentBin),m,pt2, w) ;
        FillHistogram(Form("hSingleDisp2core_cen%d",fCentBin),mcore,ptcore1, w) ;
      }
      if(ph1->IsDispOK() && ph1->IsCPVOK()){
        FillHistogram(Form("hSingleBoth_cen%d",fCentBin),m,pt1, w) ;
        FillHistogram(Form("hSingleBothcore_cen%d",fCentBin),mcore,ptcore1, w) ;
      }
      if(ph2->IsDispOK() && ph2->IsCPVOK()){
        FillHistogram(Form("hSingleBoth_cen%d",fCentBin),m,pt2, w) ;
        FillHistogram(Form("hSingleBothcore_cen%d",fCentBin),mcore,ptcore2, w) ;
      }
      if(ph1->IsDisp2OK() && ph1->IsCPVOK()){
        FillHistogram(Form("hSingleBoth2_cen%d",fCentBin),m,pt1, w) ;
        FillHistogram(Form("hSingleBoth2core_cen%d",fCentBin),mcore,ptcore1, w) ;
      }
      if(ph2->IsDisp2OK() && ph2->IsCPVOK()){
        FillHistogram(Form("hSingleBoth2_cen%d",fCentBin),m,pt2, w) ;
        FillHistogram(Form("hSingleBoth2core_cen%d",fCentBin),mcore,ptcore2, w) ;
      }


      if(a<kAlphaCut){
        FillHistogram(Form("hPi0All_a07_cen%d",fCentBin),m,pt, w) ;
      }

      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hMassPtCPV_cen%d",fCentBin) ;
	FillHistogram(Form("hMassPtV0ACPV_cen%d",fCentBin),m,pt,dphiA, w) ;
	FillHistogram(Form("hMassPtV0CCPV_cen%d",fCentBin),m,pt,dphiC, w) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCCPV_cen%d",fCentBin),m,pt,dphiT, w) ;

	FillHistogram(Form("hMassPtV0ACPVcore_cen%d",fCentBin),mcore,ptcore,dphiA, w) ;
	FillHistogram(Form("hMassPtV0CCPVcore_cen%d",fCentBin),mcore,ptcore,dphiC, w) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCCPVcore_cen%d",fCentBin),mcore,ptcore,dphiT, w) ;

	FillHistogram(Form("hPi0CPV_cen%d",fCentBin),m,pt, w) ;
	FillHistogram(Form("hPi0CPVcore_cen%d",fCentBin),mcore, ptcore, w) ;

        if(a<kAlphaCut){
          FillHistogram(Form("hPi0CPV_a07_cen%d",fCentBin),m,pt, w) ;
        }
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	FillHistogram(Form("hMassPtV0ACPV2_cen%d",fCentBin),m,pt,dphiA, w) ;
	FillHistogram(Form("hMassPtV0CCPV2_cen%d",fCentBin),m,pt,dphiC, w) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCCPV2_cen%d",fCentBin),m,pt,dphiT, w) ;
	FillHistogram(Form("hMassPtV0ACPV2core_cen%d",fCentBin),mcore,ptcore,dphiA, w) ;
	FillHistogram(Form("hMassPtV0CCPV2core_cen%d",fCentBin),mcore,ptcore,dphiC, w) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCCPV2core_cen%d",fCentBin),mcore,ptcore,dphiT, w) ;
	
	FillHistogram(Form("hPi0CPV2_cen%d",fCentBin),m,pt, w) ;
	FillHistogram(Form("hPi0CPV2core_cen%d",fCentBin),mcore, ptcore, w) ;
        if(a<kAlphaCut){
          FillHistogram(Form("hPi0CPV2_a07_cen%d",fCentBin),m,pt, w) ;
        }
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hMassPtDisp_cen%d",fCentBin) ;
	FillHistogram(Form("hMassPtV0ADisp_cen%d",fCentBin),m,pt,dphiA, w) ;
	FillHistogram(Form("hMassPtV0CDisp_cen%d",fCentBin),m,pt,dphiC, w) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCDisp_cen%d",fCentBin),m,pt,dphiT, w) ;
	
	FillHistogram(Form("hMassPtV0ADispcore_cen%d",fCentBin),mcore, ptcore,dphiA, w) ;
	FillHistogram(Form("hMassPtV0CDispcore_cen%d",fCentBin),mcore, ptcore,dphiC, w) ;
	if(fHaveTPCRP)
	  FillHistogram(Form("hMassPtTPCDispcore_cen%d",fCentBin),mcore, ptcore,dphiT, w) ;

	FillHistogram(Form("hPi0Disp_cen%d",fCentBin),m,pt, w) ;
	FillHistogram(Form("hPi0Dispcore_cen%d",fCentBin),mcore, ptcore, w) ;
	
	if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
	  FillHistogram(Form("hPi0Dispwou_cen%d",fCentBin),m,pt, w) ;

	  FillHistogram(Form("hMassPtV0ADispwou_cen%d",fCentBin),m,pt,dphiA, w) ;
 	  FillHistogram(Form("hMassPtV0CDispwou_cen%d",fCentBin),m,pt,dphiC, w) ;
	  if(fHaveTPCRP)
  	    FillHistogram(Form("hMassPtTPCDispwou_cen%d",fCentBin),m,pt,dphiT, w) ;
	}

        if(a<kAlphaCut){
          FillHistogram(Form("hPi0Disp_a07_cen%d",fCentBin),m,pt, w) ;
        }
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMassPtV0ABoth_cen%d",fCentBin),m,pt,dphiA, w) ;
	  FillHistogram(Form("hMassPtV0CBoth_cen%d",fCentBin),m,pt,dphiC, w) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMassPtTPCBoth_cen%d",fCentBin),m,pt,dphiT, w) ;

	  FillHistogram(Form("hMassPtV0ABothcore_cen%d",fCentBin),mcore,ptcore,dphiA, w) ;
	  FillHistogram(Form("hMassPtV0CBothcore_cen%d",fCentBin),mcore,ptcore,dphiC, w) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMassPtTPCBothcore_cen%d",fCentBin),mcore,ptcore,dphiT, w) ;

	  FillHistogram(Form("hPi0Both_cen%d",fCentBin),m,pt, w) ;
	  FillHistogram(Form("hPi0Bothcore_cen%d",fCentBin),mcore,ptcore, w) ;

          if(a<kAlphaCut){
            snprintf(key,55,"hPi0Both_a07_cen%d",fCentBin) ;
            FillHistogram(Form("hPi0Both_a07_cen%d",fCentBin),m,pt, w) ;
          }
          if(ph1->Module()==1 && ph2->Module()==1)
	    FillHistogram("hPi0M11", m, pt, w) ;
          else if(ph1->Module()==2 && ph2->Module()==2)
	    FillHistogram("hPi0M22", m, pt, w) ;
          else if(ph1->Module()==3 && ph2->Module()==3)
	    FillHistogram("hPi0M33", m, pt, w) ;
          else if(ph1->Module()==1 && ph2->Module()==2)
	    FillHistogram("hPi0M12", m, pt, w) ;
          else if(ph1->Module()==1 && ph2->Module()==3)
	    FillHistogram("hPi0M13", m, pt, w) ;
          else if(ph1->Module()==2 && ph2->Module()==3)
	    FillHistogram("hPi0M23", m, pt, w) ;

        }
	
      }
      
      
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	FillHistogram(Form("hPi0Disp2_cen%d",fCentBin),m,pt,w) ;
  	FillHistogram(Form("hPi0Disp2core_cen%d",fCentBin),mcore, ptcore, w) ;	

	FillHistogram(Form("hMassPtV0ADisp2_cen%d",fCentBin),m,pt,dphiA, w) ;
	FillHistogram(Form("hMassPtV0CDisp2_cen%d",fCentBin),m,pt,dphiC, w) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCDisp2_cen%d",fCentBin),m,pt,dphiT, w) ;

	FillHistogram(Form("hMassPtV0ADisp2core_cen%d",fCentBin),mcore, ptcore,dphiA, w) ;
	FillHistogram(Form("hMassPtV0CDisp2core_cen%d",fCentBin),mcore, ptcore,dphiC, w) ;
	if(fHaveTPCRP)
	  FillHistogram(Form("hMassPtTPCDisp2core_cen%d",fCentBin),mcore, ptcore,dphiT, w) ;
	  
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMassPtV0ABoth2_cen%d",fCentBin),m,pt,dphiA, w) ;
	  FillHistogram(Form("hMassPtV0CBoth2_cen%d",fCentBin),m,pt,dphiC, w) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMassPtTPCBoth2_cen%d",fCentBin),m,pt,dphiT, w) ;

	  FillHistogram(Form("hMassPtV0ABoth2core_cen%d",fCentBin),mcore,ptcore,dphiA, w) ;
	  FillHistogram(Form("hMassPtV0CBoth2core_cen%d",fCentBin),mcore,ptcore,dphiC, w) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMassPtTPCBoth2core_cen%d",fCentBin),mcore,ptcore,dphiT, w) ;

	  FillHistogram(Form("hPi0Both2_cen%d",fCentBin),m,pt, w) ;
	  FillHistogram(Form("hPi0Both2core_cen%d",fCentBin),mcore,ptcore, w) ;
	}

      }
    } // end of loop i2
  } // end of loop i1
}

//________________________________________________________________________
void AliAnalysisTaskPi0FlowMCAOD::ConsiderPi0sMix()
{
  char key[55];

  TList * arrayList = GetCaloPhotonsPHOSList(fVtxBin, fCentBin, fEMRPBin);

  for (Int_t i1=0; i1<fCaloPhotonsPHOS->GetEntriesFast(); i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;
    const Double_t w1 = ph1->GetWeight();
    for(Int_t evi=0; evi<arrayList->GetEntries();evi++){
      TObjArray * mixPHOS = static_cast<TObjArray*>(arrayList->At(evi));
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
	TLorentzVector p12  = *ph1  + *ph2;
	TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
	
	const Double_t w2 = ph2->GetWeight();
	Double_t w = TMath::Sqrt(w1*w2);

	Double_t dphiA=p12.Phi()-fRPV0A ;
	while(dphiA<0)dphiA+=TMath::Pi() ;
	while(dphiA>TMath::Pi())dphiA-=TMath::Pi() ;

	Double_t dphiC=p12.Phi()-fRPV0C ;
	while(dphiC<0)dphiC+=TMath::Pi() ;
	while(dphiC>TMath::Pi())dphiC-=TMath::Pi() ;

	Double_t dphiT=p12.Phi()-fRP ;
	while(dphiT<0)dphiT+=TMath::Pi() ;
	while(dphiT>TMath::Pi())dphiT-=TMath::Pi() ;


        Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;
        Double_t m=p12.M() ;
        Double_t mcore=pv12.M() ;
        Double_t pt=p12.Pt() ;
        Double_t ptcore=pv12.Pt() ;
        Double_t pt1=ph1->Pt() ;
        Double_t pt2=ph2->Pt() ;
        Double_t ptcore1=ph1->GetMomV2()->Pt() ;
        Double_t ptcore2=ph2->GetMomV2()->Pt() ;


	snprintf(key,55,"hMiMassPtAll_cen%d",fCentBin) ;
	FillHistogram(Form("hMiMassPtV0AAll_cen%d",fCentBin),m,pt,dphiA, w) ;
	FillHistogram(Form("hMiMassPtV0CAll_cen%d",fCentBin),m,pt,dphiC, w) ;
	if(fHaveTPCRP)
 	  FillHistogram(Form("hMiMassPtTPCAll_cen%d",fCentBin),m,pt,dphiT, w) ;

	FillHistogram(Form("hMiMassPtV0AAllcore_cen%d",fCentBin),mcore, ptcore, dphiA, w) ;
	FillHistogram(Form("hMiMassPtV0CAllcore_cen%d",fCentBin),mcore, ptcore, dphiC, w) ;
        if(fHaveTPCRP)
	  FillHistogram(Form("hMiMassPtTPCAllcore_cen%d",fCentBin),mcore, ptcore, dphiT, w) ;

	FillHistogram(Form("hMiPi0All_cen%d",fCentBin),m,pt, w) ;
	FillHistogram(Form("hMiPi0Allcore_cen%d",fCentBin),mcore,ptcore, w) ;
	if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
	  FillHistogram(Form("hMiPi0Allwou_cen%d",fCentBin),m,pt, w) ;
          FillHistogram(Form("hMiMassPtV0AAllwou_cen%d",fCentBin),m,pt,dphiA, w) ;
          FillHistogram(Form("hMiMassPtV0CAllwou_cen%d",fCentBin),m,pt,dphiC, w) ;
          if(fHaveTPCRP)
            FillHistogram(Form("hMiMassPtTPCAllwou_cen%d",fCentBin),m,pt,dphiT, w) ;
	}

        FillHistogram(Form("hMiSingleAll_cen%d",fCentBin),m,pt1, w) ;
        FillHistogram(Form("hMiSingleAll_cen%d",fCentBin),m,pt2, w) ;
        FillHistogram(Form("hMiSingleAllcore_cen%d",fCentBin),mcore,ptcore1, w) ;
        FillHistogram(Form("hMiSingleAllcore_cen%d",fCentBin),mcore,ptcore2, w) ;
        if(ph1->IsntUnfolded())
          FillHistogram(Form("hMiSingleAllwou_cen%d",fCentBin),m,pt1, w) ;
        if(ph2->IsntUnfolded())
          FillHistogram(Form("hMiSingleAllwou_cen%d",fCentBin),m,pt2, w) ;
        if(ph1->IsCPVOK()){
          FillHistogram(Form("hMiSingleCPV_cen%d",fCentBin),m,pt1, w) ;
          FillHistogram(Form("hMiSingleCPVcore_cen%d",fCentBin),mcore,ptcore1, w) ;
        }
        if(ph2->IsCPVOK()){
          FillHistogram(Form("hMiSingleCPV_cen%d",fCentBin),m,pt2, w) ;
          FillHistogram(Form("hMiSingleCPVcore_cen%d",fCentBin),mcore,ptcore2, w) ;
        }
        if(ph1->IsCPV2OK()){
          FillHistogram(Form("hMiSingleCPV2_cen%d",fCentBin),m,pt1, w) ;
          FillHistogram(Form("hMiSingleCPV2core_cen%d",fCentBin),mcore,ptcore1, w) ;
        }
        if(ph2->IsCPV2OK()){
          FillHistogram(Form("hMiSingleCPV2_cen%d",fCentBin),m,pt2, w) ;
          FillHistogram(Form("hMiSingleCPV2core_cen%d",fCentBin),mcore,ptcore2, w) ;
        }
        if(ph1->IsDispOK()){
          FillHistogram(Form("hMiSingleDisp_cen%d",fCentBin),m,pt1, w) ;
          if(ph1->IsntUnfolded()){
            FillHistogram(Form("hMiSingleDispwou_cen%d",fCentBin),m,pt1, w) ;
	  }
          FillHistogram(Form("hMiSingleDispcore_cen%d",fCentBin),mcore,ptcore1, w) ;
        }
        if(ph2->IsDispOK()){
          FillHistogram(Form("hMiSingleDisp_cen%d",fCentBin),m,pt2, w) ;
          if(ph1->IsntUnfolded()){
            FillHistogram(Form("hMiSingleDispwou_cen%d",fCentBin),m,pt2, w) ;
	  }
          FillHistogram(Form("hMiSingleDispcore_cen%d",fCentBin),mcore,ptcore2, w) ;
        }
        if(ph1->IsDisp2OK()){
          FillHistogram(Form("hMiSingleDisp2_cen%d",fCentBin),m,pt1, w) ;
          FillHistogram(Form("hMiSingleDisp2core_cen%d",fCentBin),mcore,ptcore1, w) ;
        }
        if(ph2->IsDisp2OK()){
          FillHistogram(Form("hMiSingleDisp2_cen%d",fCentBin),m,pt2, w) ;
          FillHistogram(Form("hMiSingleDisp2core_cen%d",fCentBin),mcore,ptcore2, w) ;
        }
        if(ph1->IsDispOK() && ph1->IsCPVOK()){
          snprintf(key,55,"hMiSingleBoth_cen%d",fCentBin) ;
          FillHistogram(key,m,pt1, w) ;
          snprintf(key,55,"hMiSingleBothcore_cen%d",fCentBin);
          FillHistogram(key,mcore,ptcore1, w) ;
        }
        if(ph2->IsDispOK() && ph2->IsCPVOK()){
          snprintf(key,55,"hMiSingleBoth_cen%d",fCentBin);
          FillHistogram(key,m,pt2, w) ;
          snprintf(key,55,"hMiSingleBothcore_cen%d",fCentBin);
          FillHistogram(key,mcore,ptcore2, w) ;
        }
        if(ph1->IsDisp2OK() && ph1->IsCPVOK()){
          FillHistogram(Form("hMiSingleBoth2_cen%d",fCentBin),m,pt1, w) ;
          FillHistogram(Form("hMiSingleBoth2core_cen%d",fCentBin),mcore,ptcore1, w) ;
        }
        if(ph2->IsDisp2OK() && ph2->IsCPVOK()){
          FillHistogram(Form("hMiSingleBoth2_cen%d",fCentBin),m,pt2, w) ;
          FillHistogram(Form("hMiSingleBoth2core_cen%d",fCentBin),mcore,ptcore2, w) ;
        }



        if(a<kAlphaCut){
          FillHistogram(Form("hMiPi0All_a07_cen%d",fCentBin),m,pt, w) ;
        }
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMiMassPtV0ACPV_cen%d",fCentBin),m,pt,dphiA, w) ;
	  FillHistogram(Form("hMiMassPtV0CCPV_cen%d",fCentBin),m,pt,dphiC, w) ;
	  if(fHaveTPCRP)
 	    FillHistogram(Form("hMiMassPtTPCCPV_cen%d",fCentBin),m,pt,dphiT, w) ;

	  FillHistogram(Form("hMiMassPtV0ACPVcore_cen%d",fCentBin),mcore, ptcore,dphiA, w) ;
	  FillHistogram(Form("hMiMassPtV0CCPVcore_cen%d",fCentBin),mcore, ptcore,dphiC, w) ;
	  if(fHaveTPCRP)
 	    FillHistogram(Form("hMiMassPtTPCCPVcore_cen%d",fCentBin),mcore, ptcore,dphiT, w) ;

	  FillHistogram(Form("hMiPi0CPV_cen%d",fCentBin),m,pt, w) ;
	  FillHistogram(Form("hMiPi0CPVcore_cen%d",fCentBin),mcore, ptcore, w) ;

	  if(a<kAlphaCut){
            FillHistogram(Form("hMiPi0CPV_a07_cen%d",fCentBin),m,pt, w) ;
          }
	}
	if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	  FillHistogram(Form("hMiPi0CPV2_cen%d",fCentBin),m,pt, w) ;
	  FillHistogram(Form("hMiPi0CPV2core_cen%d",fCentBin),mcore, ptcore, w) ;

	  FillHistogram(Form("hMiMassPtV0ACPV2_cen%d",fCentBin),m,pt,dphiA, w) ;
	  FillHistogram(Form("hMiMassPtV0CCPV2_cen%d",fCentBin),m,pt,dphiC, w) ;
	  if(fHaveTPCRP)
 	    FillHistogram(Form("hMiMassPtTPCCPV2_cen%d",fCentBin),m,pt,dphiT, w) ;
	  FillHistogram(Form("hMiMassPtV0ACPV2core_cen%d",fCentBin),mcore,ptcore,dphiA, w) ;
	  FillHistogram(Form("hMiMassPtV0CCPV2core_cen%d",fCentBin),mcore,ptcore,dphiC, w) ;
	  if(fHaveTPCRP)
 	    FillHistogram(Form("hMiMassPtTPCCPV2core_cen%d",fCentBin),mcore,ptcore,dphiT, w) ;

	  if(a<kAlphaCut){
            FillHistogram(Form("hMiPi0CPV2_a07_cen%d",fCentBin),m,pt, w) ;
          }
	}
	if(ph1->IsDispOK() && ph2->IsDispOK()){
	  FillHistogram(Form("hMiMassPtV0ADisp_cen%d",fCentBin),m,pt,dphiA, w) ;
	  FillHistogram(Form("hMiMassPtV0CDisp_cen%d",fCentBin),m,pt,dphiC, w) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMiMassPtTPCDisp_cen%d",fCentBin),m,pt,dphiT, w) ;

	  FillHistogram(Form("hMiMassPtV0ADispcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiA, w) ;
	  FillHistogram(Form("hMiMassPtV0CDispcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiC, w) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMiMassPtTPCDispcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiT, w) ;


	  FillHistogram(Form("hMiPi0Disp_cen%d",fCentBin),m,pt, w) ;
	  FillHistogram(Form("hMiPi0Dispcore_cen%d",fCentBin),pv12.M(),pv12.Pt(), w) ;
          if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
	    FillHistogram(Form("hMiPi0Dispwou_cen%d",fCentBin),m,pt, w) ;
	    FillHistogram(Form("hMiMassPtV0ADispwou_cen%d",fCentBin),m,pt,dphiA, w) ;
	    FillHistogram(Form("hMiMassPtV0CDispwou_cen%d",fCentBin),m,pt,dphiC, w) ;
            if(fHaveTPCRP)
	      FillHistogram(Form("hMiMassPtTPCDispwou_cen%d",fCentBin),m,pt,dphiT, w) ;
	  }

	  if(a<kAlphaCut){
            FillHistogram(Form("hMiPi0Disp_a07_cen%d",fCentBin),m,pt, w) ;
          }
	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	    FillHistogram(Form("hMiMassPtV0ABoth_cen%d",fCentBin),m,pt,dphiA, w) ;
	    FillHistogram(Form("hMiMassPtV0CBoth_cen%d",fCentBin),m,pt,dphiC, w) ;
	    if(fHaveTPCRP)
  	      FillHistogram(Form("hMiMassPtTPCBoth_cen%d",fCentBin),m,pt,dphiT, w) ;

	    FillHistogram(Form("hMiMassPtV0ABothcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiA, w) ;
	    FillHistogram(Form("hMiMassPtV0CBothcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiC, w) ;
	    if(fHaveTPCRP)
  	      FillHistogram(Form("hMiMassPtTPCBothcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiT, w) ;

	    FillHistogram(Form("hMiPi0Both_cen%d",fCentBin),m,pt, w) ;
	    FillHistogram(Form("hMiPi0Bothcore_cen%d",fCentBin),pv12.M(),pv12.Pt(), w) ;

	    if(a<kAlphaCut){
              FillHistogram(Form("hMiPi0Both_a07_cen%d",fCentBin),m,pt, w) ;
            }
	  }
	}
	
  	if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	  FillHistogram(Form("hMiMassPtV0ADisp2_cen%d",fCentBin),m,pt,dphiA, w) ;
	  FillHistogram(Form("hMiMassPtV0CDisp2_cen%d",fCentBin),m,pt,dphiC, w) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMiMassPtTPCDisp2_cen%d",fCentBin),m,pt,dphiT, w) ;

	  FillHistogram(Form("hMiMassPtV0ADisp2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiA, w) ;
	  FillHistogram(Form("hMiMassPtV0CDisp2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiC, w) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMiMassPtTPCDisp2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiT, w) ;


	  FillHistogram(Form("hMiPi0Disp2_cen%d",fCentBin),m,pt, w) ;
	  FillHistogram(Form("hMiPi0Disp2core_cen%d",fCentBin),pv12.M(),pv12.Pt(), w) ;

	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	    FillHistogram(Form("hMiMassPtV0ABoth2_cen%d",fCentBin),m,pt,dphiA, w) ;
	    FillHistogram(Form("hMiMassPtV0CBoth2_cen%d",fCentBin),m,pt,dphiC, w) ;
	    if(fHaveTPCRP)
  	      FillHistogram(Form("hMiMassPtTPCBoth2_cen%d",fCentBin),m,pt,dphiT, w) ;

	    FillHistogram(Form("hMiMassPtV0ABoth2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiA, w) ;
	    FillHistogram(Form("hMiMassPtV0CBoth2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiC, w) ;
	    if(fHaveTPCRP)
  	      FillHistogram(Form("hMiMassPtTPCBoth2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiT, w) ;

	    FillHistogram(Form("hMiPi0Both2_cen%d",fCentBin),m,pt, w) ;
	    FillHistogram(Form("hMiPi0Both2core_cen%d",fCentBin),pv12.M(),pv12.Pt(), w) ;

	  }
	}
      } // end of loop i2
    }
  } // end of loop i1
}


void AliAnalysisTaskPi0FlowMCAOD::ProcessMC()
{
  FillMCHist();
  FillSecondaries() ;
}



//___________________________________________________________________________
void AliAnalysisTaskPi0FlowMCAOD::FillMCHist(){
  //fill histograms for efficiensy etc. calculation
  
  //---------First pi0/eta-----------------------------
  char partName[10] ; char hkey[55] ;
  Int_t ntrack = fMcArray->GetEntriesFast(); 
  
  for(Int_t i=0;i<ntrack;i++){
    AliAODMCParticle* particle = GetParticle(i);
      
    if(particle->PdgCode() == kPi0)
      snprintf(partName,10,"pi0") ;
    else
      if(particle->PdgCode() == kEta)
	snprintf(partName,10,"eta") ;
      else
	if(particle->PdgCode() == kGamma)
	  snprintf(partName,10,"gamma") ;
	else
	  if(particle->PdgCode() == 310)
	    snprintf(partName,10,"K0S") ;
	  else
	    continue ;
    //Primary particle
    Double_t r  = R(particle) ;
    Double_t pt = particle->Pt() ;
    //Distribution over vertex
    FillHistogram(Form("hMC_%s_vertex",partName),pt,r) ;
    
    if(r >kRCut)
      continue ;

    Double_t phi=particle->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    
    Double_t phig = 180./TMath::Pi()*phi; // phi in deg
    
    //Total number of pi0 with creation radius <1 cm
    Double_t weight = PrimaryParticleWeight(particle) ;  
    snprintf(hkey,55,"hMC_all_%s_cen%d",partName,fCentBin) ;
    
    FillHistogram(hkey,pt,weight) ;
    
    if(TMath::Abs(particle->Y())<0.135 && phig>240. && phig<320.){
      snprintf(hkey,55,"hMC_unitEta_%s_cen%d",partName,fCentBin) ;
      FillHistogram(hkey,pt,weight) ;
      
      snprintf(hkey,55,"hMC_rap_%s_cen%d",partName,fCentBin) ;
      FillHistogram(hkey,particle->Y(),weight) ;
      
      snprintf(hkey,55,"hMC_phi_%s_cen%d",partName,fCentBin) ;
      FillHistogram(hkey,phi,weight) ;
    }  
  }   
  
}
 


//________________________________________________________________________
void AliAnalysisTaskPi0FlowMCAOD::FillSecondaries()
{
  //Sort secondaires
  
  //Fill spectra of primary particles 
  //with proper weight
  if( fDebug )
    AliInfo("start");
  
  Int_t ntrack = fMcArray->GetEntriesFast(); 
  
  for(Int_t i=0; i<ntrack; i++){
    AliAODMCParticle* p = GetParticle(i);
    
    if(R(p)>kRCut)
      continue ;
    if( TMath::Abs(p->Pt())<1.e-6 )
      continue;
    if(TMath::Abs(p->Y())>0.5)
      continue ;
    Double_t w = PrimaryParticleWeight(p) ;  
    Int_t primPdgCode=p->PdgCode() ;
    switch(primPdgCode){
    case  kGamma: FillHistogram(Form("hPrimPhot_cen%d",fCentBin),p->Pt(),w); 
      break ;
    case  kElectron: 
    case -kElectron: 
      FillHistogram(Form("hPrimEl_cen%d",fCentBin),p->Pt(),w); 
      break ;
    case  kPi0: 
      FillHistogram(Form("hPrimPi0_cen%d",fCentBin),p->Pt(),w); 
      break ;
    case  kEta: 
      FillHistogram(Form("hPrimEta_cen%d",fCentBin),p->Pt(),w); 
      break ;
    case  kPiPlus: 
    case  kPiMinus: 
      FillHistogram(Form("hPrimPipm_cen%d",fCentBin),p->Pt(),w); 
      break ;		  
    case  kProton:  //p 
      FillHistogram(Form("hPrimP_cen%d",fCentBin),p->Pt(),w); 
      break ;		  
    case kProtonBar:  //pbar
      FillHistogram(Form("hPrimPbar_cen%d",fCentBin),p->Pt(),w); 
      break ;		  
    case  kNeutron:  //n 
      FillHistogram(Form("hPrimN_cen%d",fCentBin),p->Pt(),w); 
      break ;		  
    case  kNeutronBar:  //nbar
      FillHistogram(Form("hPrimNbar_cen%d",fCentBin),p->Pt(),w); 
      break ;
    case  310:  //nbar
      FillHistogram(Form("hPrimK0S_cen%d",fCentBin),p->Pt(),w); 
      break ;
    case  130:  //nbar
      FillHistogram(Form("hPrimK0L_cen%d",fCentBin),p->Pt(),w); 
      break ;
    case  321:  //K+
    case -321:  //K-
      FillHistogram(Form("hPrimKpm_cen%d",fCentBin),p->Pt(),w); 
      break ;
    default:	   //other
      FillHistogram(Form("hPrimOther_cen%d",fCentBin),p->Pt(),w);    
    }
  }
  if(fDebug)
    AliInfo("Origins of secondary pi0s");
  //Origins of secondary pi0s
  for(Int_t i=0; i<ntrack; i++){
    AliAODMCParticle* p = GetParticle(i);
    if(p->PdgCode()!=111)
      continue ;
    FillHistogram("Vertex",p->Pt(),R(p));
    if(R(p)<kRCut)
      continue ;
    Double_t phi=p->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    FillHistogram("hSecondPi0RphiZ",R(p),phi,p->Zv()) ;   
    Double_t w = PrimaryParticleWeight(p) ;  
    FillHistogram("hSecondPi0RE",R(p),p->Pt(),w) ;   
  }

  TLorentzVector p1;

  Int_t inPHOS=fCaloPhotonsPHOS->GetEntries() ;
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;
    Double_t w1=ph1->GetWeight() ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i2) ;
      TLorentzVector p12  = *ph1  + *ph2;
      Double_t w2=ph2->GetWeight() ;
      Double_t w = TMath::Sqrt(w1*w2) ;
      FillHistogram(Form("hParentAll_cen%d",fCentBin),p12.M(),p12.Pt(),w) ;
      Int_t prim=FindCommonParent(ph1->GetPrimary(),ph2->GetPrimary()) ;
      if(prim>-1){
        AliAODMCParticle * particle = GetParticle(prim);
        FillHistogram("hMass_R",p12.M(),p12.Pt(),R(particle)) ;
		
	
        Int_t pdgCode=particle->PdgCode() ;
        if(pdgCode!=111){ //common parent - not pi0
          if(pdgCode==22)  
            FillHistogram(Form("hParentGamma_cen%d",fCentBin),p12.M(),p12.Pt(),w) ;
	  else{		    
            if(pdgCode==11 || pdgCode==-11){   
              FillHistogram(Form("hParentEl_cen%d",fCentBin),p12.M(),p12.Pt(),w) ;
	    }
  	    else{
              if(InPi0mass(p12.M() ,p12.Pt())){
	        if(fDebug >1) AliInfo(Form("Common parent: %d",pdgCode));
	      }
              FillHistogram(Form("hParentOther_cen%d",fCentBin),p12.M(),p12.Pt(),w) ;
	    }
	  }//Not photons
        }//Common parent not pi0
        else{ //common parent - pi0
          FillHistogram(Form("hParentPi0_cen%d",fCentBin),p12.M(),p12.Pt(),w) ;	
          FillHistogram(Form("Real_pi_R"),p12.M(),p12.Pt(),R(particle),w) ;	
          FillHistogram(Form("Real_pi_Z"),p12.M(),p12.Pt(),particle->Zv(),w) ;	
	  if(R(particle)<kRCut && TMath::Abs(particle->Zv())<fMaxAbsVertexZ){
            FillHistogram(Form("hParentDirPi0_cen%d",fCentBin),p12.M(),p12.Pt(),w) ;
	    continue ;
	  }
	  //Common particle pi0, created off-vertex
  	  Int_t primPi0 = particle->GetMother();
	  if(primPi0==-1){
            FillHistogram(Form("hParentPi0NoPrim_cen%d",fCentBin),p12.M(),p12.Pt(),w) ;
	  }
	  else{
    	    Int_t primPdgCode=((AliAODMCParticle *)fMcArray->At(primPi0))->PdgCode();
            switch(primPdgCode){
            case 221: FillHistogram(Form("hParentPi0Eta_cen%d",fCentBin),p12.M(),p12.Pt(),w) ; //eta
	      break ;
            case 223: FillHistogram(Form("hParentPi0Omega_cen%d",fCentBin),p12.M(),p12.Pt(),w) ; //omega
	      break ;
	    case  211:  //pi+-
	    case -211: FillHistogram(Form("hParentPi0Pipm_cen%d",fCentBin),p12.M(),p12.Pt(),w) ; //
	      break ;
	    case  321:  //K+-
	    case -321: FillHistogram(Form("hParentPi0Kpm_cen%d",fCentBin),p12.M(),p12.Pt(),w) ; //
	      break ;
	    case 310: FillHistogram(Form("hParentPi0Ks_cen%d",fCentBin),p12.M(),p12.Pt(),w) ; // K0S
	      break ;
	    case 130: FillHistogram(Form("hParentPi0Kl_cen%d",fCentBin),p12.M(),p12.Pt(),w) ; // K0L
	      break ;
	    case  2212:  //p 
	    case  2112:  //n 
	      FillHistogram(Form("hParentPi0pn_cen%d",fCentBin),p12.M(),p12.Pt(),w) ; // pn
	      break ;
	    case -2212:  //pbar
	    case -2112:  //nbar
	      FillHistogram(Form("hParentPi0antipn_cen%d",fCentBin),p12.M(),p12.Pt(),w) ; // pn
	      break ;
	    default:	   //other
	      FillHistogram(Form("hParentPi0Other_cen%d",fCentBin),p12.M(),p12.Pt(),w) ; //
	    }//switch	  
          }//pi0 with primary
        }//common parent - pi0
      }//there is common primary 
    }//seond photon loop
  }//first photon loop
  
  
  //Now look at photon contaiminations
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
  
    AliCaloPhoton * ph1=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;
    Int_t iprim = ph1->GetPrimary() ;
    if(iprim<0)
      FillAllHistograms(Form("hGammaNoPrim_cen%d",fCentBin),ph1) ; //
    else{
      //Find primary at vertex
      AliAODMCParticle * primPHOS = GetParticle(iprim);
      Int_t iprimV = primPHOS->GetMother();
      AliAODMCParticle * primVtx = primPHOS ;
      while((iprimV>-1) && R(primVtx)>kRCut){
	primVtx = GetParticle(iprimV);
        iprimV = primVtx->GetMother();
      }
    
      //photon
      Int_t primPdgCode=primVtx->PdgCode() ;
      switch(primPdgCode){
      case  22: FillAllHistograms("hGammaPhot",ph1); 
	break ;
      case  11: 
      case -11: 
	FillAllHistograms("hGammaEl",ph1); 
	break ;
      case  111: 
	FillAllHistograms("hGammaPi0",ph1); 
	break ;
      case  221: 
	FillAllHistograms("hGammaEta",ph1); 
	break ;
      case 223: FillAllHistograms("hGammaOmega",ph1) ; //omega
	break ;
      case  211: 
      case -211: 
	FillAllHistograms("hGammaPipm",ph1); 
	//Find particle entered PHOS
	if(primVtx == primPHOS)
	  FillAllHistograms("hGammaPipmDirect",ph1); 
	else{
	  Int_t primPdgPHOS=primPHOS->PdgCode() ;
	  if(primPdgPHOS==22){
	    FillAllHistograms("hGammaPipmGamma",ph1); 
	    FillHistogram("hPipmGammaConvR",ph1->Pt(),R(primPHOS));
	    FillHistogram("hPipmGammaConvRZ",primPHOS->Zv(),R(primPHOS));
	    break ;		  
	  }
	  if(TMath::Abs(primPdgPHOS)==11){
	    FillAllHistograms("hGammaPipmEl",ph1); 
	    FillHistogram("hPipmElConvR",ph1->Pt(),R(primPHOS));
	    break ;		  
	  }
	  if(TMath::Abs(primPdgPHOS)==2212){
	    FillAllHistograms("hGammaPipmp",ph1); 
	    FillHistogram("hPipmNConvR",ph1->Pt(),R(primPHOS));
	    break ;		  
	  }
	  if(TMath::Abs(primPdgPHOS)==2112){
	    FillAllHistograms("hGammaPipmn",ph1); 
	    FillHistogram("hPipmNConvR",ph1->Pt(),R(primPHOS));
	    break ;		  
	  }
	  FillAllHistograms("hGammaPipmOther",ph1); 
	  FillHistogram("hPipmOtherConvR",ph1->Pt(),R(primPHOS));		    
	}
	break ;		  
      case  2212:  //p 
	FillAllHistograms("hGammaP",ph1); 
	break ;		  
      case -2212:  //pbar
	FillAllHistograms("hGammaPbar",ph1); 
	break ;		  
      case  2112:  //n 
	FillAllHistograms("hGammaN",ph1); 
	break ;		  
      case -2112:  //nbar
	FillAllHistograms("hGammaNbar",ph1) ; // pn
	break ;
      case  310:  //nbar
	FillAllHistograms("hGammaK0S",ph1) ; // pn
	break ;
      case  130:  //nbar
	FillAllHistograms("hGammaK0L",ph1) ; // pn
	break ;
      case  321:  //K+
      case -321:  //K-
	FillAllHistograms("hGammaKpm",ph1) ; // pn
	break ;
      case -323: 
      case  323: 
      case -313: 
      case  313: FillAllHistograms("hGammaKstar",ph1) ; // K*(892)
	break ;
		  
      case -2224 : //Deltas
      case  2224 : //Deltas
      case -2214 : //Deltas
      case  2214 : //Deltas
      case -2114 : //Deltas
      case  2114 : //Deltas
      case -1114 : //Deltas
      case  1114 : //Deltas
	FillAllHistograms("hGammaDelta",ph1) ; // pn
	break ;		  
      default:	   //other
	if(primVtx->Charge())
	  FillAllHistograms("hGammaOtherCharged",ph1) ; //
	else
	  FillAllHistograms("hGammaOtherNeutral",ph1) ; //
      }
    }
  
  }//single photons
  
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0FlowMCAOD::FillAllHistograms(const char * particleName,AliCaloPhoton * ph)
{
  //Fill All PID histograms
        
  Double_t w=ph->GetWeight() ;
  Double_t pt = ph->Pt() ;
  Double_t ptC=ph->GetMomV2()->Pt() ;
  FillHistogram(Form("%s_All_cen%d",particleName,fCentBin),pt,w) ;
  FillHistogram(Form("%s_Allcore_cen%d",particleName,fCentBin),ptC,w) ;
  if(ph->IsCPVOK()){
    FillHistogram(Form("%s_CPV_cen%d",particleName,fCentBin),pt,w) ;
    FillHistogram(Form("%s_CPVcore_cen%d",particleName,fCentBin),ptC,w) ;
  }
  if(ph->IsCPV2OK()){
    FillHistogram(Form("%s_CPV2_cen%d",particleName,fCentBin),pt,w) ;
    FillHistogram(Form("%s_CPV2core_cen%d",particleName,fCentBin),ptC,w) ;
  }
  if(ph->IsDispOK()){     
    FillHistogram(Form("%s_Disp_cen%d",particleName,fCentBin),pt,w) ;
    FillHistogram(Form("%s_Dispcore_cen%d",particleName,fCentBin),ptC,w) ;
    if(ph->IsDisp2OK()){
      FillHistogram(Form("%s_Disp2_cen%d",particleName,fCentBin),pt,w) ;
      FillHistogram(Form("%s_Disp2core_cen%d",particleName,fCentBin),ptC,w) ;
    }
    if(ph->IsCPVOK()){
      FillHistogram(Form("%s_Both_cen%d",particleName,fCentBin),pt,w) ;
      FillHistogram(Form("%s_Bothcore_cen%d",particleName,fCentBin),ptC,w) ;
    }
  }  
}


//___________________________________________________________________________
Double_t AliAnalysisTaskPi0FlowMCAOD::PrimaryWeight(Int_t /* primary */){
  return 1.;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPi0FlowMCAOD::PrimaryParticleWeight(AliAODMCParticle * /* particle */){
  return 1.;
}

//________________________________________________________________________
TClonesArray* AliAnalysisTaskPi0FlowMCAOD::GetMCArray()
{
  fMcArray = 0;
  AliAODInputHandler* aodHandler=dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  
  if (aodHandler){
    AliAODEvent *aod=aodHandler->GetEvent();
    if (aod) {
      fMcArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fMcArray) AliError("Could not retrieve MC array!");
    }
    else AliError("Could not retrieve AOD event!");
  }
  
  return fMcArray;
}

AliAODMCParticle* AliAnalysisTaskPi0FlowMCAOD::GetParticle(Int_t particlepos) {
  //Returns particle at given position.
  
  if(!fMcArray){
    AliError("MC array is not initialized, run GetMCArray() first!");
    return 0;
  }

  return (AliAODMCParticle *) fMcArray->At(particlepos);
}


Int_t AliAnalysisTaskPi0FlowMCAOD::FindPrimary(AliVCluster*clu,  Bool_t&sure){
  //Finds primary and estimates if it unique one?
  //First check can it be photon/electron

  const Double_t emFraction=0.9; //part of energy of cluster to be assigned to EM particle
  Int_t n = clu->GetNLabels();

  for(Int_t i=0;  i<n;  i++){
    Int_t label = clu->GetLabelAt(i);
    AliAODMCParticle*  p=  GetParticle(label) ;
    Int_t pdg = p->PdgCode() ;
    if(pdg==22  ||  pdg==11 || pdg == -11){
      if(p->E()>emFraction*clu->E()){
	sure=kTRUE ;
	return label;
      }
    }
  }

  Double_t*  Ekin=  new  Double_t[n] ;

  for(Int_t i=0;  i<n;  i++){
    AliAODMCParticle*  p = GetParticle(clu->GetLabelAt(i)) ;
    Ekin[i]=p->P() ;  // estimate of kinetic energy
    if(p->PdgCode()==-2212  ||  p->PdgCode()==-2112){
      Ekin[i]+=1.8  ;  //due to annihilation
    }
  }
  Int_t iMax=0;
  Double_t eMax=0.,eSubMax=0. ;
  for(Int_t i=0;  i<n;  i++){
    if(Ekin[i]>eMax){
      eSubMax=eMax;
      eMax=Ekin[i];
      iMax=i;
    }
  }
  if(eSubMax>0.8*eMax)//not obvious primary
    sure=kFALSE;
  else
    sure=kTRUE;
  delete[]  Ekin;
  return clu->GetLabelAt(iMax);
}

//________________________________________________________________________
Int_t AliAnalysisTaskPi0FlowMCAOD::FindCommonParent(Int_t iPart, Int_t jPart){
  //check if there is a common parent for particles i and j
  // -1: no common parent or wrong iPart/jPart
  
  Int_t ntrack = fMcArray->GetEntriesFast();
  if(iPart==-1 || iPart>=ntrack || jPart==-1 || jPart>=ntrack) return -1;
  
  Int_t iprim1 = iPart;
  
  while(iprim1>-1){  
    Int_t iprim2=jPart;
    
    while(iprim2>-1){
      if(iprim1==iprim2)
	return iprim1 ;
      iprim2 = GetParticle(iprim2)->GetMother();
    }
    
    iprim1 = GetParticle(iprim1)->GetMother();
  }
  
  return -1;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPi0FlowMCAOD::HaveParent(Int_t iPart, Int_t pdgParent){
  //check if there is a common parent for particles i and j
  // -1: no common parent or wrong iPart/jPart
  
  Int_t ntrack = fMcArray->GetEntriesFast();
  if(iPart==-1 || iPart>=ntrack) return -1;
  
  Int_t iprim1=iPart;
  while(iprim1>-1){  
    AliAODMCParticle * tmp = GetParticle(iprim1) ;
    if(tmp->PdgCode()==pdgParent)
      return kTRUE ;
    iprim1 = tmp->GetMother();
  }
  return kFALSE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPi0FlowMCAOD::InPi0mass(Double_t m, Double_t /*pt*/){

 return TMath::Abs(m-0.135)<0.007*2.5 ;
}

//________________________________________________________________________
Double32_t AliAnalysisTaskPi0FlowMCAOD::R(AliAODMCParticle* p) {
  //Radius of vertex.

  Double32_t x = p->Xv();
  Double32_t y = p->Yv();
  Double32_t z = p->Zv();

  return x*x + y*y + z*z;
}
