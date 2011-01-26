#include "AliAnalysisTaskSE.h"
#include "AliTrackPointArray.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliStack.h"
#include "AliPID.h"
#include "AliITSPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisTaskITSsaTracks.h"

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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
// Implementation of class AliAnalysisTaskITSsaTracks
// AliAnalysisTaskSE to extract QA and performance histos for ITS standalone tracks
// 
//
// Authors: L. Milano, milano@to.infn.it
//          F. Prino, prino@to.infn.it
//          
//*************************************************************************

ClassImp(AliAnalysisTaskITSsaTracks)
//______________________________________________________________________________
AliAnalysisTaskITSsaTracks::AliAnalysisTaskITSsaTracks() : AliAnalysisTaskSE("ITSsa resolution"), 
  fOutput(0),
  fHistNEvents(0),
  fHistPtTPCITSAll(0),
  fHistPtTPCITSGood(0),
  fHistPtTPCITSFake(0),
  fHistPtITSsaAll(0),
  fHistPtITSsaGood(0),
  fHistPtITSsaFake(0),
  fHistPtITSpureSAAll(0),
  fHistPtITSpureSAGood(0),
  fHistPtITSpureSAFake(0),
  fHistdedxvsPtITSpureSA3cls(0),
  fHistdedxvsPITSpureSA3cls(0),
  fHistdedxvsPtITSpureSA4cls(0),
  fHistdedxvsPITSpureSA4cls(0),
  fNPtBins(kMaxPtBins),
  fMinITSpts(4),
  fMinSPDpts(1),
  fMinITSptsForPid(3),
  fMinITSptsForMatch(6),
  fMinTPCpts(50),
  fReadMC(kFALSE),
  fUseMCId(kFALSE)
{
  //
  for(Int_t ilay=0; ilay<6;ilay++) fRequirePoint[ilay]=kFALSE;
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  const Int_t nbins = 29;
  Double_t xbins[nbins+1]={0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,
			   0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,
			   0.85,0.90,0.95,1.00,1.20,1.40,1.60,1.80,1.90,2.00};
  SetPtBins(nbins,xbins);
}


//___________________________________________________________________________
AliAnalysisTaskITSsaTracks::~AliAnalysisTaskITSsaTracks(){
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}
   
//________________________________________________________________________
void AliAnalysisTaskITSsaTracks::SetPtBins(Int_t n, Double_t* lim){
  // define pt bins for analysis
  if(n>kMaxPtBins){
    printf("Max. number of Pt bins = %d\n",kMaxPtBins);
    return;
  }else{
    fNPtBins=n;
    for(Int_t i=0; i<fNPtBins+1; i++) fPtLimits[i]=lim[i];
    for(Int_t i=fNPtBins+1; i<kMaxPtBins+1; i++) fPtLimits[i]=99999999.;
  }
}
//___________________________________________________________________________
void AliAnalysisTaskITSsaTracks::UserCreateOutputObjects() {
  // create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",3,-0.5,2.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  //binning for the dedx histogram
  const Int_t hnbinsdedx=400;
  Double_t hxmindedx = 0.01;
  Double_t hxmaxdedx = 10;
  Double_t hlogxmindedx = TMath::Log10(hxmindedx);
  Double_t hlogxmaxdedx = TMath::Log10(hxmaxdedx);
  Double_t hbinwidthdedx = (hlogxmaxdedx-hlogxmindedx)/hnbinsdedx;
  Double_t hxbinsdedx[hnbinsdedx+1];
  hxbinsdedx[0] = 0.01; 
  for (Int_t i=1;i<=hnbinsdedx;i++) {
    hxbinsdedx[i] = hxmindedx + TMath::Power(10,hlogxmindedx+i*hbinwidthdedx);
  }
  
  fHistdedxvsPtITSpureSA3cls = new TH2F("hdedxvsPtITSpureSA3cls","",hnbinsdedx,hxbinsdedx,900,0,1000);
  fHistdedxvsPtITSpureSA3cls->Sumw2();
  fOutput->Add(fHistdedxvsPtITSpureSA3cls);
  
  fHistdedxvsPITSpureSA3cls = new TH2F("hdedxvsPITSpureSA3cls","",hnbinsdedx,hxbinsdedx,900,0,1000);
  fHistdedxvsPITSpureSA3cls->Sumw2();
  fOutput->Add(fHistdedxvsPITSpureSA3cls);
  
  fHistdedxvsPtITSpureSA4cls = new TH2F("hdedxvsPtITSpureSA4cls","",hnbinsdedx,hxbinsdedx,900,0,1000);
  fHistdedxvsPtITSpureSA4cls->Sumw2();
  fOutput->Add(fHistdedxvsPtITSpureSA4cls);
  
  fHistdedxvsPITSpureSA4cls = new TH2F("hdedxvsPITSpureSA4cls","",hnbinsdedx,hxbinsdedx,900,0,1000);
  fHistdedxvsPITSpureSA4cls->Sumw2();
  fOutput->Add(fHistdedxvsPITSpureSA4cls);
  
  TString spname[3]={"Pion","Kaon","Proton"};
  TString hisname;
  const Int_t nbins = fNPtBins;
  Double_t xbins[nbins+1];
  for(Int_t ibin=0; ibin<=nbins; ibin++) xbins[ibin]=fPtLimits[ibin];


  fHistPtTPCITSAll = new TH1F("hPtTPCITSAll","",100,0.,2.);
  fHistPtTPCITSAll->Sumw2();
  fOutput->Add(fHistPtTPCITSAll);
  fHistPtTPCITSGood = new TH1F("hPtTPCITSGood","",100,0.,2.);
  fHistPtTPCITSGood->Sumw2();
  fOutput->Add(fHistPtTPCITSGood);
  fHistPtTPCITSFake = new TH1F("hPtTPCITSFake","",100,0.,2.);
  fHistPtTPCITSFake->Sumw2();
  fOutput->Add(fHistPtTPCITSFake);

  fHistPtITSsaAll  = new TH1F("hPtITSsaAll","",100,0.,2.);
  fHistPtITSsaAll->Sumw2();
  fOutput->Add(fHistPtITSsaAll);
  fHistPtITSsaGood  = new TH1F("hPtITSsaGood","",100,0.,2.);
  fHistPtITSsaGood->Sumw2();
  fOutput->Add(fHistPtITSsaGood);
  fHistPtITSsaFake  = new TH1F("hPtITSsaFake","",100,0.,2.);
  fHistPtITSsaFake->Sumw2();
  fOutput->Add(fHistPtITSsaFake);

  fHistPtITSpureSAAll = new TH1F("hPtITSpureSAAll","",100,0.,2.);
  fHistPtITSpureSAAll->Sumw2();
  fOutput->Add(fHistPtITSpureSAAll);
  fHistPtITSpureSAGood = new TH1F("hPtITSpureSAGood","",100,0.,2.);
  fHistPtITSpureSAGood->Sumw2();
  fOutput->Add(fHistPtITSpureSAGood);
  fHistPtITSpureSAFake = new TH1F("hPtITSpureSAFake","",100,0.,2.);
  fHistPtITSpureSAFake->Sumw2();
  fOutput->Add(fHistPtITSpureSAFake);

  for(Int_t iSpec=0; iSpec<kNspecies; iSpec++){

    hisname.Form("hPtTPCITS%s",spname[iSpec].Data());
    fHistPtTPCITS[iSpec] = new TH1F(hisname.Data(),"",100,0.,2.);
    fHistPtTPCITS[iSpec]->Sumw2();
    fOutput->Add(fHistPtTPCITS[iSpec]);

    hisname.Form("hPtITSsa%s",spname[iSpec].Data());
    fHistPtITSsa[iSpec] = new TH1F(hisname.Data(),"",100,0.,2.);
    fHistPtITSsa[iSpec]->Sumw2();
    fOutput->Add(fHistPtITSsa[iSpec]);

    hisname.Form("hPtITSpureSA%s",spname[iSpec].Data());
    fHistPtITSpureSA[iSpec] = new TH1F(hisname.Data(),"",100,0.,2.);
    fHistPtITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistPtITSpureSA[iSpec]);

    //---

    hisname.Form("hEtaPhiTPCITS%s",spname[iSpec].Data());
    fHistEtaPhiTPCITS[iSpec] = new TH2F(hisname.Data(),"",50,-1,1.,50,0.,2.*TMath::Pi());
    fHistEtaPhiTPCITS[iSpec]->Sumw2();
    fOutput->Add(fHistEtaPhiTPCITS[iSpec]);

    hisname.Form("hEtaPhiITSsa%s",spname[iSpec].Data());
    fHistEtaPhiITSsa[iSpec] = new TH2F(hisname.Data(),"",50,-1,1.,50,0.,2.*TMath::Pi());
    fHistEtaPhiITSsa[iSpec]->Sumw2();
    fOutput->Add(fHistEtaPhiITSsa[iSpec]);

    hisname.Form("hEtaPhiITSpureSA%s",spname[iSpec].Data());
    fHistEtaPhiITSpureSA[iSpec] = new TH2F(hisname.Data(),"",50,-1,1.,50,0.,2.*TMath::Pi());
    fHistEtaPhiITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistEtaPhiITSpureSA[iSpec]);

    //---

    hisname.Form("hNcluTPCITS%s",spname[iSpec].Data());
    fHistNcluTPCITS[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-0.5,6.5);
    fHistNcluTPCITS[iSpec]->Sumw2();
    fOutput->Add(fHistNcluTPCITS[iSpec]);

    hisname.Form("hNcluITSsa%s",spname[iSpec].Data());
    fHistNcluITSsa[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-0.5,6.5);
    fHistNcluITSsa[iSpec]->Sumw2();
    fOutput->Add(fHistNcluITSsa[iSpec]);

    hisname.Form("hNcluITSpureSA%s",spname[iSpec].Data());
    fHistNcluITSpureSA[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-0.5,6.5);
    fHistNcluITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistNcluITSpureSA[iSpec]);

    //---

    hisname.Form("hd0rphiITSpureSA%s",spname[iSpec].Data());
    fHistd0rphiITSpureSA[iSpec] = new TH2F(hisname.Data(),"",nbins,xbins,2000,-1,1);
    fHistd0rphiITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistd0rphiITSpureSA[iSpec]);

    hisname.Form("hd0zITSpureSA%s",spname[iSpec].Data());
    fHistd0zITSpureSA[iSpec] = new TH2F(hisname.Data(),"",nbins,xbins,2000,-1,1);
    fHistd0zITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistd0zITSpureSA[iSpec]);

    //---

    hisname.Form("hCluInLayTPCITS%s",spname[iSpec].Data());
    fHistCluInLayTPCITS[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-1.5,5.5);
    fHistCluInLayTPCITS[iSpec]->Sumw2();
    fOutput->Add(fHistCluInLayTPCITS[iSpec]);
    
    hisname.Form("hCluInLayITSsa%s",spname[iSpec].Data());
    fHistCluInLayITSsa[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-1.5,5.5);
    fHistCluInLayITSsa[iSpec]->Sumw2();
    fOutput->Add(fHistCluInLayITSsa[iSpec]);

    hisname.Form("hCluInLayITSpureSA%s",spname[iSpec].Data());
    fHistCluInLayITSpureSA[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-1.5,5.5);
    fHistCluInLayITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistCluInLayITSpureSA[iSpec]);
    
    hisname.Form("hOuterLayITSpureSA%s",spname[iSpec].Data());
    fHistOuterLayITSpureSA[iSpec] = new TH2F(hisname.Data(),"",100,0.,2.,7,-1.5,5.5);
    fHistOuterLayITSpureSA[iSpec]->Sumw2();
    fOutput->Add(fHistOuterLayITSpureSA[iSpec]);
    
    //---

    hisname.Form("hPtResid%s",spname[iSpec].Data());
    fHistPtResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-1.,1.);
    fHistPtResid[iSpec]->Sumw2();
    fOutput->Add(fHistPtResid[iSpec]);
    hisname.Form("hPtRelResid%s",spname[iSpec].Data());
    fHistPtRelResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-0.5,0.5);
    fHistPtRelResid[iSpec]->Sumw2();
    fOutput->Add(fHistPtRelResid[iSpec]);
    
    hisname.Form("hInvPtResid%s",spname[iSpec].Data());
    fHistInvPtResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-1.,1.);
    fHistInvPtResid[iSpec]->Sumw2();
    fOutput->Add(fHistInvPtResid[iSpec]);
    hisname.Form("hInvPtRelResid%s",spname[iSpec].Data());
    fHistInvPtRelResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-0.5,0.5);
    fHistInvPtRelResid[iSpec]->Sumw2();
    fOutput->Add(fHistInvPtRelResid[iSpec]);
    
    hisname.Form("hMCPtResid%s",spname[iSpec].Data());
    fHistMCPtResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-1.,1.);
    fHistMCPtResid[iSpec]->Sumw2();
    fOutput->Add(fHistMCPtResid[iSpec]);
    hisname.Form("hMCPtRelResid%s",spname[iSpec].Data());
    fHistMCPtRelResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-0.5,0.5);
    fHistMCPtRelResid[iSpec]->Sumw2();
    fOutput->Add(fHistMCPtRelResid[iSpec]);
    
    hisname.Form("hMCInvPtResid%s",spname[iSpec].Data());
    fHistMCInvPtResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-1.,1.);
    fHistMCInvPtResid[iSpec]->Sumw2();
    fOutput->Add(fHistMCInvPtResid[iSpec]);
    hisname.Form("hMCInvPtRelResid%s",spname[iSpec].Data());
    fHistMCInvPtRelResid[iSpec]=new TH2F(hisname.Data(),hisname.Data(),nbins,xbins,100,-0.5,0.5);
    fHistMCInvPtRelResid[iSpec]->Sumw2();
    fOutput->Add(fHistMCInvPtRelResid[iSpec]);
    
  }

}
//______________________________________________________________________________
void AliAnalysisTaskITSsaTracks::UserExec(Option_t *)
{
  //

  AliESDEvent *esd = (AliESDEvent*) (InputEvent());


  if(!esd) {
    printf("AliAnalysisTaskSDDRP::Exec(): bad ESD\n");
    return;
  } 


  if(!ESDfriend()) {
    printf("AliAnalysisTaskSDDRP::Exec(): bad ESDfriend\n");
    return;
  }
  
  AliStack* stack=0;

  if(fReadMC){
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    stack = mcEvent->Stack();
    if (!stack) {
      Printf("ERROR: stack not available");
      return;
    }
  }

  PostData(1, fOutput);

  fHistNEvents->Fill(0);
  const AliESDVertex *spdv=esd->GetPrimaryVertexSPD();
  if(spdv->GetNContributors()<=0) return;
  fHistNEvents->Fill(1);

  Int_t ntracks = esd->GetNumberOfTracks();
  for (Int_t iTrack=0; iTrack < ntracks; iTrack++) {
    AliESDtrack * track = esd->GetTrack(iTrack);
    if (!track) continue;
    Int_t status=track->GetStatus();
    if(!(status&AliESDtrack::kITSrefit)) continue;
    Bool_t isSA=kTRUE;
    if(status&AliESDtrack::kTPCin) isSA=kFALSE;
    Int_t nTPCclus=track->GetNcls(1);
    Int_t nITSclus=track->GetNcls(0);
    if(nITSclus<fMinITSpts)continue;

    UChar_t clumap=track->GetITSClusterMap();
    Int_t nSPDPoints=0;
    for(Int_t i=0; i<2; i++){
      if(clumap&(1<<i)) ++nSPDPoints;
    }
    if(nSPDPoints<fMinSPDpts) continue;

    Int_t nPointsForPid=0;
    for(Int_t i=2; i<6; i++){
      if(clumap&(1<<i)) ++nPointsForPid;
    }
    
    Float_t pttrack=track->Pt();
    Float_t ptrack=track->P();
    Float_t dedx=track->GetITSsignal();
    Int_t hadronSpecie=-1;
    if(status&AliESDtrack::kTPCin){ 
      fHistPtTPCITSAll->Fill(pttrack);
    }else{
      if(status&AliESDtrack::kITSpureSA){
	fHistPtITSpureSAAll->Fill(pttrack);
      }else{
	fHistPtITSsaAll->Fill(pttrack);
      }
    }

    if(fReadMC && fUseMCId){
      Int_t trlabel=track->GetLabel();
      if(trlabel>=0){
	if(status&AliESDtrack::kTPCin){ 
	  fHistPtTPCITSGood->Fill(pttrack);
	}else{
	  if(status&AliESDtrack::kITSpureSA){
	    fHistPtITSpureSAGood->Fill(pttrack);
	  }else{
	    fHistPtITSsaGood->Fill(pttrack);
	  }
	}
      }else{
	if(status&AliESDtrack::kTPCin){ 
	  fHistPtTPCITSFake->Fill(pttrack);
	}else{
	  if(status&AliESDtrack::kITSpureSA){
	    fHistPtITSpureSAFake->Fill(pttrack);
	  }else{
	    fHistPtITSsaFake->Fill(pttrack);
	  }
	}
	continue; // fake track
      }
      TParticle* part = stack->Particle(trlabel);
      Int_t pdg=TMath::Abs(part->GetPdgCode());
      if(pdg==211) hadronSpecie=kPion;
      else if(pdg==321) hadronSpecie=kKaon;
      else if(pdg==2212) hadronSpecie=kProton;
    }else{
      if(nPointsForPid<fMinITSptsForPid)continue;
      AliITSPIDResponse pidits(fReadMC);
      Float_t nSigmaPion=pidits.GetNumberOfSigmas(ptrack,dedx,AliPID::kPion,nPointsForPid,isSA);
      if(nSigmaPion>-2. && nSigmaPion<2.){
	hadronSpecie=kPion;
      }else{
	Float_t nSigmaKaon=pidits.GetNumberOfSigmas(ptrack,dedx,AliPID::kKaon,nPointsForPid,isSA);
	if(nSigmaKaon>-2. && nSigmaKaon<2.){
	  hadronSpecie=kKaon;
	}else{
	  Float_t nSigmaProton=pidits.GetNumberOfSigmas(ptrack,dedx,AliPID::kProton,nPointsForPid,isSA);
	  if(nSigmaProton>-2. && nSigmaProton<2.){
	    hadronSpecie=kProton;
	  }
	}
      } 
    }
    if(hadronSpecie<0) continue;
    if(status&AliESDtrack::kTPCin){ // TPC+ITS tracks
      fHistPtTPCITS[hadronSpecie]->Fill(pttrack);
      fHistEtaPhiTPCITS[hadronSpecie]->Fill(track->Eta(),track->Phi());
      fHistNcluTPCITS[hadronSpecie]->Fill(pttrack,nITSclus);
      fHistCluInLayTPCITS[hadronSpecie]->Fill(-1.);
      for(Int_t iBit=0; iBit<6; iBit++){
	if(clumap&(1<<iBit)) fHistCluInLayTPCITS[hadronSpecie]->Fill(pttrack,iBit);
      }
    }else{ // ITS standalone and pureSA tracks
      if(status&AliESDtrack::kITSpureSA){
	Float_t impactXY=-999, impactZ=-999;
	track->GetImpactParameters(impactXY, impactZ);
	fHistPtITSpureSA[hadronSpecie]->Fill(pttrack);
	fHistEtaPhiITSpureSA[hadronSpecie]->Fill(track->Eta(),track->Phi());
	fHistNcluITSpureSA[hadronSpecie]->Fill(pttrack,nITSclus);
	fHistd0rphiITSpureSA[hadronSpecie]->Fill(pttrack,impactXY);
	fHistd0zITSpureSA[hadronSpecie]->Fill(pttrack,impactZ);
	if(nPointsForPid==3){
	  fHistdedxvsPtITSpureSA3cls->Fill(pttrack,dedx);
	  fHistdedxvsPITSpureSA3cls->Fill(ptrack,dedx);
	}
	if(nPointsForPid==4){
	  fHistdedxvsPtITSpureSA4cls->Fill(pttrack,dedx);
	  fHistdedxvsPITSpureSA4cls->Fill(ptrack,dedx);
	}
	fHistCluInLayITSpureSA[hadronSpecie]->Fill(-1.);
	Int_t outerLay=-1;
	for(Int_t iBit=0; iBit<6; iBit++){
	  if(clumap&(1<<iBit)){
	    fHistCluInLayITSpureSA[hadronSpecie]->Fill(pttrack,iBit);
	    if(iBit>outerLay) outerLay=iBit;
	  }
	}
	fHistOuterLayITSpureSA[hadronSpecie]->Fill(pttrack,outerLay);	
	
	
	if(fReadMC){  
	  Int_t trlabel=track->GetLabel();
	  if(trlabel<0)continue; // fake track
	  TParticle* part = stack->Particle(trlabel);
	  Float_t ptgen=part->Pt();
	  Float_t invpttrack=track->OneOverPt();
	  Float_t invptgen=0.;
	  if(ptgen>0.) invptgen=1./ptgen;
	  fHistMCPtResid[hadronSpecie]->Fill(pttrack,pttrack-ptgen);
	  fHistMCPtRelResid[hadronSpecie]->Fill(pttrack,(pttrack-ptgen)/ptgen);
	  fHistMCInvPtResid[hadronSpecie]->Fill(pttrack,invpttrack-invptgen);
	  fHistMCInvPtRelResid[hadronSpecie]->Fill(pttrack,(invpttrack-invptgen)/invptgen);	  
	}
      }else{
	fHistPtITSsa[hadronSpecie]->Fill(pttrack);
	fHistEtaPhiITSsa[hadronSpecie]->Fill(track->Eta(),track->Phi());
	fHistNcluITSsa[hadronSpecie]->Fill(pttrack,nITSclus);
	fHistCluInLayITSsa[hadronSpecie]->Fill(-1.);
	for(Int_t iBit=0; iBit<6; iBit++){
	  if(clumap&(1<<iBit)) fHistCluInLayITSsa[hadronSpecie]->Fill(pttrack,iBit);
	}
      }
    }

    if(nITSclus<fMinITSptsForMatch || nTPCclus<fMinTPCpts) continue;      
    Int_t idxMI[12],idxSA[12];
    for(Int_t icl=0; icl<12; icl++){ 
      idxMI[icl]=-1; 
      idxSA[icl]=-1;
    }
    Int_t ncls=track->GetClusters(0,idxMI);
    if(fMinITSpts<6){
      Bool_t accept=kTRUE;
      for(Int_t ilay=0; ilay<6; ilay++){
	if(fRequirePoint[ilay]){
	  Int_t mask = 1<<ilay;
	  if(!(clumap & mask)){ 
	    accept=kFALSE;
	    break;
	  }
	}
      }
      if(!accept) continue;
    }
    // Sort
    for(Int_t i=0;i<12;i++){
      for(Int_t j=i+1;j<12;j++){
	if(idxMI[j]>idxMI[i]){
	  Int_t tmp=idxMI[j];
	  idxMI[j]=idxMI[i];
	  idxMI[i]=tmp;
	}
      }
    }
    //    for(Int_t i=0; i<12; i++) printf("%d ",idxMI[i]);
    //    printf("\n");
    if(idxMI[0]<0 && idxMI[0]==idxMI[1]) continue;
    Bool_t matched=kFALSE;
    for (Int_t iTrack2 = 0; iTrack2 < ntracks; iTrack2++) {
      if(matched) break;
      if(iTrack2==iTrack) continue;
      AliESDtrack* track2 = esd->GetTrack(iTrack2);
      Int_t status2=track2->GetStatus();
      if(!(status2&AliESDtrack::kITSrefit)) continue;
      if(!(status2&AliESDtrack::kITSpureSA)) continue;      
      Int_t clumap2=track2->GetITSClusterMap();
      Int_t nITSclus2=track2->GetNcls(0);
      Int_t nTPCclus2=track2->GetNcls(1);
      if(nITSclus2<fMinITSpts || nTPCclus2>0) continue; 
      Int_t ncls2=track2->GetClusters(0,idxSA);
      if(ncls2!=ncls) continue;
      if(fMinITSpts<6){
	Bool_t accept=kTRUE;
	for(Int_t ilay=0; ilay<6; ilay++){
	  if(fRequirePoint[ilay]){
	    Int_t mask = 1<<ilay;
	    if(!(clumap2 & mask)){ 
	      accept=kFALSE;
	      break;
	    }
	  }
	}
	if(!accept) continue;
      }
      // Sort
      for(Int_t i=0;i<12;i++){
	for(Int_t j=i+1;j<12;j++){
	  if(idxSA[j]>idxSA[i]){
	    Int_t tmp=idxSA[j];
	    idxSA[j]=idxSA[i];
	    idxSA[i]=tmp;
	  }
	}
      }
      Int_t match=0;
      for(Int_t icl=0; icl<ncls; icl++){
	if(idxSA[icl]!=idxMI[icl]){
	  match=0; 
	  break;
	}
	else match++;
      }
      if(match==ncls && match>0){
	matched=kTRUE;
	Float_t pt1=track->Pt();
	Float_t pt2=track2->Pt();
	Float_t ptm1=track->OneOverPt();
	Float_t ptm2=track2->OneOverPt();
	fHistPtResid[hadronSpecie]->Fill(pt1,pt2-pt1);
	fHistPtRelResid[hadronSpecie]->Fill(pt1,(pt2-pt1)/pt1);
	fHistInvPtResid[hadronSpecie]->Fill(pt1,ptm2-ptm1);
	fHistInvPtRelResid[hadronSpecie]->Fill(pt1,(ptm2-ptm1)/ptm1);	
      }
    }
  }

  PostData(1,fOutput);
  
}
//______________________________________________________________________________
void AliAnalysisTaskITSsaTracks::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  fHistdedxvsPtITSpureSA3cls= dynamic_cast<TH2F*>(fOutput->FindObject("hdedxvsPtITSpureSA3cls"));
  fHistdedxvsPITSpureSA3cls= dynamic_cast<TH2F*>(fOutput->FindObject("hdedxvsPITSpureSA3cls"));
  fHistdedxvsPtITSpureSA4cls= dynamic_cast<TH2F*>(fOutput->FindObject("hdedxvsPtITSpureSA4cls"));
  fHistdedxvsPITSpureSA4cls= dynamic_cast<TH2F*>(fOutput->FindObject("hdedxvsPITSpureSA4cls"));
    
  TString spname[3]={"Pion","Kaon","Proton"};

  fHistPtTPCITSAll =dynamic_cast<TH1F*>(fOutput->FindObject("hPtTPCITSAll"));
  fHistPtTPCITSGood =dynamic_cast<TH1F*>(fOutput->FindObject("hPtTPCITSGood"));
  fHistPtTPCITSFake =dynamic_cast<TH1F*>(fOutput->FindObject("hPtTPCITSFake"));

  fHistPtITSsaAll =dynamic_cast<TH1F*>(fOutput->FindObject("hPtITSsaAll"));
  fHistPtITSsaGood =dynamic_cast<TH1F*>(fOutput->FindObject("hPtITSsaGood"));
  fHistPtITSsaFake =dynamic_cast<TH1F*>(fOutput->FindObject("hPtITSsaFake"));

  fHistPtITSpureSAAll =dynamic_cast<TH1F*>(fOutput->FindObject("hPtITSpureSAAll"));
  fHistPtITSpureSAGood =dynamic_cast<TH1F*>(fOutput->FindObject("hPtITSpureSAGood"));
  fHistPtITSpureSAFake =dynamic_cast<TH1F*>(fOutput->FindObject("hPtITSpureSAFake"));

  for(Int_t iSpec=0; iSpec<kNspecies; iSpec++){
    fHistPtTPCITS[iSpec]= dynamic_cast<TH1F*>(fOutput->FindObject(Form("hPtTPCITS%s",spname[iSpec].Data())));
    fHistPtITSsa[iSpec]= dynamic_cast<TH1F*>(fOutput->FindObject(Form("hPtITSsa%s",spname[iSpec].Data())));
    fHistPtITSpureSA[iSpec]= dynamic_cast<TH1F*>(fOutput->FindObject(Form("hPtITSpureSA%s",spname[iSpec].Data())));

    fHistEtaPhiTPCITS[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hEtaPhiTPCITS%s",spname[iSpec].Data())));
    fHistEtaPhiITSsa[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hEtaPhiITSsa%s",spname[iSpec].Data())));
    fHistEtaPhiITSpureSA[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hEtaPhiITSpureSA%s",spname[iSpec].Data())));
    
    fHistNcluTPCITS[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hNcluTPCITS%s",spname[iSpec].Data())));
    fHistNcluITSsa[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hNcluITSsa%s",spname[iSpec].Data())));
    fHistNcluITSpureSA[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hNcluITSpureSA%s",spname[iSpec].Data())));
    
    fHistd0rphiITSpureSA[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hd0rphiITSpureSA%s",spname[iSpec].Data())));
    fHistd0zITSpureSA[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hd0zITSpureSA%s",spname[iSpec].Data())));
    
    fHistCluInLayTPCITS[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hCluInLayTPCITS%s",spname[iSpec].Data())));
    fHistCluInLayITSsa[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hCluInLayITSsa%s",spname[iSpec].Data())));
    fHistCluInLayITSpureSA[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hCluInLayITSpureSA%s",spname[iSpec].Data())));
    
    fHistOuterLayITSpureSA[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hOuterLayITSpureSA%s",spname[iSpec].Data())));
  
    fHistPtResid[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hPtResid%s",spname[iSpec].Data())));
    fHistPtRelResid[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hPtRelResid%s",spname[iSpec].Data())));
    fHistInvPtResid[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hInvPtResid%s",spname[iSpec].Data())));
    fHistInvPtRelResid[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hInvPtRelResid%s",spname[iSpec].Data())));
    fHistMCPtResid[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hMCPtResid%s",spname[iSpec].Data())));
    fHistMCPtRelResid[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hMCPtRelResid%s",spname[iSpec].Data())));
    fHistMCInvPtResid[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hMCInvPtResid%s",spname[iSpec].Data())));
    fHistMCInvPtRelResid[iSpec]= dynamic_cast<TH2F*>(fOutput->FindObject(Form("hMCInvPtRelResid%s",spname[iSpec].Data())));
  
  }
  return;
}






