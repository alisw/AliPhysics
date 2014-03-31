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
#include "AliVertexingHFUtils.h"
#include "AliAnalysisTaskCombinHF.h"

ClassImp(AliAnalysisTaskCombinHF)


//________________________________________________________________________
AliAnalysisTaskCombinHF::AliAnalysisTaskCombinHF():
  AliAnalysisTaskSE(),
  fOutput(0x0), 
  fHistNEvents(0x0),
  fHistTrackStatus(0x0),
  fHistCheckOrigin(0x0),
  fHistCheckOriginSel(0x0),
  fHistCheckDecChan(0x0),
  fHistCheckDecChanAcc(0x0),
  fPtVsYGen(0x0),
  fPtVsYGenLimAcc(0x0),
  fPtVsYGenAcc(0x0),
  fPtVsYReco(0x0),
  fMassVsPtVsY(0x0),
  fMassVsPtVsYRot(0x0),
  fMassVsPtVsYLSpp(0x0),
  fMassVsPtVsYLSmm(0x0),
  fMassVsPtVsYSig(0x0),
  fMassVsPtVsYRefl(0x0),
  fMassVsPtVsYBkg(0x0),
  fNSelected(0x0),
  fNormRotated(0x0),
  fDeltaMass(0x0),
  fDeltaMassFullAnalysis(0x0),
  fFilterMask(BIT(4)),
  fTrackCutsAll(0x0),
  fTrackCutsPion(0x0),
  fTrackCutsKaon(0x0),
  fPidHF(new AliAODPidHF()),
  fAnalysisCuts(0x0),
  fMinMass(1.750),
  fMaxMass(2.150),
  fEtaAccCut(0.9),
  fPtAccCut(0.1),
  fNRotations(9),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fMinAngleForRot3(2*TMath::Pi()/6),
  fMaxAngleForRot3(4*TMath::Pi()/6),
  fCounter(0x0),
  fMeson(kDzero),
  fReadMC(kFALSE),
  fPromptFeeddown(kPrompt),
  fGoUpToQuark(kTRUE),
  fFullAnalysis(0)  
{
  // default constructor
}

//________________________________________________________________________
AliAnalysisTaskCombinHF::AliAnalysisTaskCombinHF(Int_t meson, AliRDHFCuts* analysiscuts):
  AliAnalysisTaskSE("DmesonCombin"),
  fOutput(0x0), 
  fHistNEvents(0x0),
  fHistTrackStatus(0x0),
  fHistCheckOrigin(0x0),
  fHistCheckOriginSel(0x0),
  fHistCheckDecChan(0x0),
  fHistCheckDecChanAcc(0x0),
  fPtVsYGen(0x0),
  fPtVsYGenLimAcc(0x0),
  fPtVsYGenAcc(0x0),
  fPtVsYReco(0x0),
  fMassVsPtVsY(0x0),
  fMassVsPtVsYRot(0x0),
  fMassVsPtVsYLSpp(0x0),
  fMassVsPtVsYLSmm(0x0),
  fMassVsPtVsYSig(0x0),
  fMassVsPtVsYRefl(0x0),
  fMassVsPtVsYBkg(0x0),
  fNSelected(0x0),
  fNormRotated(0x0),
  fDeltaMass(0x0),
  fDeltaMassFullAnalysis(0x0),
  fFilterMask(BIT(4)),
  fTrackCutsAll(0x0),
  fTrackCutsPion(0x0),
  fTrackCutsKaon(0x0),
  fPidHF(new AliAODPidHF()),
  fAnalysisCuts(analysiscuts),
  fMinMass(1.750),
  fMaxMass(2.150),
  fEtaAccCut(0.9),
  fPtAccCut(0.1),
  fNRotations(9),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fMinAngleForRot3(2*TMath::Pi()/6),
  fMaxAngleForRot3(4*TMath::Pi()/6),
  fCounter(0x0),
  fMeson(meson),
  fReadMC(kFALSE),
  fPromptFeeddown(1),
  fGoUpToQuark(kTRUE),
  fFullAnalysis(0)

{
  // standard constructor

  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,AliNormalizationCounter::Class());
 }

//________________________________________________________________________
AliAnalysisTaskCombinHF::~AliAnalysisTaskCombinHF()
{
  //
  // Destructor
  //
  delete fOutput;
  delete fHistNEvents;
  delete fHistTrackStatus;
  delete fHistCheckOrigin;
  delete fHistCheckOriginSel;
  delete fHistCheckDecChan;
  delete fHistCheckDecChanAcc;
  delete fPtVsYGen;
  delete fPtVsYGenLimAcc;
  delete fPtVsYGenAcc;
  delete fPtVsYReco;
  delete fMassVsPtVsY; 
  delete fMassVsPtVsYLSpp;
  delete fMassVsPtVsYLSmm;
  delete fMassVsPtVsYRot;
  delete fMassVsPtVsYSig;
  delete fMassVsPtVsYRefl;
  delete fMassVsPtVsYBkg;
  delete fNSelected;
  delete fNormRotated;
  delete fDeltaMass;
  delete fCounter;
  delete fTrackCutsAll;
  delete fTrackCutsPion;
  delete fTrackCutsKaon;
  delete fPidHF;
  delete fAnalysisCuts;
}

//________________________________________________________________________
void AliAnalysisTaskCombinHF::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskCombinHF::UserCreateOutputObjects() \n");

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "number of events ",8,-0.5,7.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"n. passing IsEvSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"n. rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"n. rejected due to not reco vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"n. rejected for contr vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"n. rejected for vertex out of accept");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"n. rejected for pileup events");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"no. of out centrality events");

  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

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
  fHistTrackStatus->Sumw2();
  fHistTrackStatus->SetMinimum(0);
  fOutput->Add(fHistTrackStatus);

  if(fReadMC){

    fHistCheckOrigin=new TH1F("hCheckOrigin","",7,-1.5,5.5);
    fHistCheckOrigin->Sumw2();
    fHistCheckOrigin->SetMinimum(0);
    fOutput->Add(fHistCheckOrigin);

    fHistCheckOriginSel=new TH1F("hCheckOriginSel","",7,-1.5,5.5);
    fHistCheckOriginSel->Sumw2();
    fHistCheckOriginSel->SetMinimum(0);
    fOutput->Add(fHistCheckOriginSel);

    fHistCheckDecChan=new TH1F("hCheckDecChan","",7,-2.5,4.5);
    fHistCheckDecChan->Sumw2();
    fHistCheckDecChan->SetMinimum(0);
    fOutput->Add(fHistCheckDecChan);

    fHistCheckDecChanAcc=new TH1F("hCheckDecChanAcc","",7,-2.5,4.5);
    fHistCheckDecChanAcc->Sumw2();
    fHistCheckDecChanAcc->SetMinimum(0);
    fOutput->Add(fHistCheckDecChanAcc);

    fPtVsYGen= new TH2F("hPtVsYGen","",20,0.,10.,20,-1.,1.);
    fPtVsYGen->Sumw2();
    fPtVsYGen->SetMinimum(0);
    fOutput->Add(fPtVsYGen);

    fPtVsYGenLimAcc= new TH2F("hPtVsYGenLimAcc","",20,0.,10.,20,-1.,1.);
    fPtVsYGenLimAcc->Sumw2();
    fPtVsYGenLimAcc->SetMinimum(0);
    fOutput->Add(fPtVsYGenLimAcc);

    fPtVsYGenAcc= new TH2F("hPtVsYGenAcc","",20,0.,10.,20,-1.,1.);
    fPtVsYGenAcc->Sumw2();
    fPtVsYGenAcc->SetMinimum(0);
    fOutput->Add(fPtVsYGenAcc);

    fPtVsYReco= new TH2F("hPtVsYReco","",20,0.,10.,20,-1.,1.);
    fPtVsYReco->Sumw2();
    fPtVsYReco->SetMinimum(0);
    fOutput->Add(fPtVsYReco);
  }


 Int_t nMassBins=fMaxMass*1000.-fMinMass*1000.;
  Double_t maxm=fMinMass+nMassBins*0.001;
  fMassVsPtVsY=new TH3F("hMassVsPtVsY","",nMassBins,fMinMass,maxm,20,0.,10.,20,-1.,1.);
  fMassVsPtVsY->Sumw2();
  fMassVsPtVsY->SetMinimum(0);
  fOutput->Add(fMassVsPtVsY);

  fMassVsPtVsYRot=new TH3F("hMassVsPtVsYRot","",nMassBins,fMinMass,maxm,20,0.,10.,20,-1.,1.);
  fMassVsPtVsYRot->Sumw2();
  fMassVsPtVsYRot->SetMinimum(0);
  fOutput->Add(fMassVsPtVsYRot);

  if(fMeson==kDzero){
    fMassVsPtVsYLSpp=new TH3F("hMassVsPtVsYLSpp","",nMassBins,fMinMass,maxm,20,0.,10.,20,-1.,1.);
    fMassVsPtVsYLSpp->Sumw2();
    fMassVsPtVsYLSpp->SetMinimum(0);
    fOutput->Add(fMassVsPtVsYLSpp);
    fMassVsPtVsYLSmm=new TH3F("hMassVsPtVsYLSmm","",nMassBins,fMinMass,maxm,20,0.,10.,20,-1.,1.);
    fMassVsPtVsYLSmm->Sumw2();
    fMassVsPtVsYLSmm->SetMinimum(0);
    fOutput->Add(fMassVsPtVsYLSmm);
  }

  fMassVsPtVsYSig=new TH3F("hMassVsPtVsYSig","",nMassBins,fMinMass,maxm,20,0.,10.,20,-1.,1.);
  fMassVsPtVsYSig->Sumw2();
  fMassVsPtVsYSig->SetMinimum(0);
  fOutput->Add(fMassVsPtVsYSig);

  fMassVsPtVsYRefl=new TH3F("hMassVsPtVsYRefl","",nMassBins,fMinMass,maxm,20,0.,10.,20,-1.,1.);
  fMassVsPtVsYRefl->Sumw2();
  fMassVsPtVsYRefl->SetMinimum(0);
  fOutput->Add(fMassVsPtVsYRefl);

  fMassVsPtVsYBkg=new TH3F("hMassVsPtVsYBkg","",nMassBins,fMinMass,maxm,20,0.,10.,20,-1.,1.);
  fMassVsPtVsYBkg->Sumw2();
  fMassVsPtVsYBkg->SetMinimum(0);
  fOutput->Add(fMassVsPtVsYBkg);

  fNSelected=new TH1F("hNSelected","",100,-0.5,99.5);
  fNSelected->Sumw2();
  fNSelected->SetMinimum(0);
  fOutput->Add(fNSelected);

  fNormRotated=new TH1F("hNormRotated","",11,-0.5,10.5);
  fNormRotated->Sumw2();
  fNormRotated->SetMinimum(0);
  fOutput->Add(fNormRotated);

  fDeltaMass=new TH1F("hDeltaMass","",100,-0.4,0.4);
  fDeltaMass->Sumw2();
  fDeltaMass->SetMinimum(0);
  fOutput->Add(fDeltaMass);

  Int_t binSparseDMassRot[5]={nMassBins,100,24,40,20};
  Double_t edgeLowSparseDMassRot[5]={fMinMass,-0.4,0.,-4.,0};
  Double_t edgeHighSparseDMassRot[5]={maxm,0.4,12.,4.,3.14};
  fDeltaMassFullAnalysis=new THnSparseF("fDeltaMassFullAnalysis","fDeltaMassFullAnalysis;inv mass (GeV/c);#Delta inv mass (GeV/c) ; p_{T}^{D} (GeV/c); #Delta p_{T} (GeV/c); daughter angle (2prongs) (rad);",5,binSparseDMassRot,edgeLowSparseDMassRot,edgeHighSparseDMassRot);
  fOutput->Add(fDeltaMassFullAnalysis);

  //Counter for Normalization
  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();

  PostData(1,fOutput); 
  PostData(2,fCounter);   
}

//________________________________________________________________________
void AliAnalysisTaskCombinHF::UserExec(Option_t */*option*/){
  //Build the 3-track combinatorics (+-+ and -+-) for D+->Kpipi decays

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
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
  fPidHF->SetPidResponse(pidResp);
  
 
  fHistNEvents->Fill(0); // count event
  // Post the data already here
  PostData(1,fOutput);

  fCounter->StoreEvent(aod,fAnalysisCuts,fReadMC);
 
  Bool_t isEvSel=fAnalysisCuts->IsEventSelected(aod);
  if(fAnalysisCuts->IsEventRejectedDueToTrigger())fHistNEvents->Fill(2);
  if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex())fHistNEvents->Fill(3);
  if(fAnalysisCuts->IsEventRejectedDueToVertexContributors())fHistNEvents->Fill(4);
  if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())fHistNEvents->Fill(5);
  //  if(fAnalysisCuts->IsEventRejectedDueToPileup())fHistNEvents->Fill(6);
  if(fAnalysisCuts->IsEventRejectedDueToCentrality()) fHistNEvents->Fill(7); 

  if(!isEvSel)return;
 
  fHistNEvents->Fill(1);
  
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
    FillGenHistos(arrayMC);
  }


  Int_t ntracks=aod->GetNTracks();

  // select and flag tracks
  UChar_t* status = new UChar_t[ntracks];
  for(Int_t iTr=0; iTr<ntracks; iTr++){
    status[iTr]=0;
    AliAODTrack* track=aod->GetTrack(iTr);
    if(IsTrackSelected(track)) status[iTr]+=1;
    if(IsKaon(track)) status[iTr]+=2;
    if(IsPion(track)) status[iTr]+=4;
    fHistTrackStatus->Fill(status[iTr]);
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
  //  UInt_t pdgs[3]={321,321,211};
  Double_t tmpp[3];
  Double_t px[3],py[3],pz[3];
  Int_t dgLabels[3];

  for(Int_t iTr1=0; iTr1<ntracks; iTr1++){
    AliAODTrack* trK=aod->GetTrack(iTr1);
    if((status[iTr1] & 1)==0) continue;
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
      AliAODTrack* trPi1=aod->GetTrack(iTr2);
      Int_t chargePi1=trPi1->Charge();
      trPi1->GetPxPyPz(tmpp);
      px[1] = tmpp[0]; 
      py[1] = tmpp[1]; 
      pz[1] = tmpp[2]; 
      dgLabels[1]=trPi1->GetLabel();
      if(chargePi1==chargeK){
	if(fMeson==kDzero) FillLSHistos(421,2,tmpRD2,px,py,pz,pdg0,chargePi1);
	continue;
      }
      if(fMeson==kDzero){
	nFiltered++;
	v2->AddDaughter(trK);
	v2->AddDaughter(trPi1);
	tmpRD2->SetSecondaryVtx(v2);
	Bool_t ok=FillHistos(421,2,tmpRD2,px,py,pz,pdg0,arrayMC,dgLabels);
	v2->RemoveDaughters();
	if(ok) nSelected++;
      }else{
	for(Int_t iTr3=iTr2+1; iTr3<ntracks; iTr3++){
	  if((status[iTr3] & 1)==0) continue;
	  if((status[iTr3] & 4)==0) continue;
	  if(iTr1==iTr3) continue;
	  AliAODTrack* trPi2=aod->GetTrack(iTr3);
	  Int_t chargePi2=trPi2->Charge();
	  if(chargePi2==chargeK) continue;
	  nFiltered++;
	  trPi2->GetPxPyPz(tmpp);
	  px[2] = tmpp[0]; 
	  py[2] = tmpp[1]; 
	  pz[2] = tmpp[2]; 
	  dgLabels[2]=trPi2->GetLabel();
	  v3->AddDaughter(trK);
	  v3->AddDaughter(trPi1);
	  v3->AddDaughter(trPi2);
	  tmpRD3->SetSecondaryVtx(v3);
	  Bool_t ok=FillHistos(411,3,tmpRD3,px,py,pz,pdgp,arrayMC,dgLabels);
	  v3->RemoveDaughters();
	  if(ok) nSelected++;
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

  PostData(1,fOutput); 
  PostData(2,fCounter);    

  return;
}

//________________________________________________________________________
void AliAnalysisTaskCombinHF::FillLSHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, Int_t charge){
  // Fill histos for LS candidates
    
  tmpRD->SetPxPyPzProngs(nProngs,px,py,pz);
  Double_t pt = tmpRD->Pt();
  Double_t minv2 = tmpRD->InvMass2(nProngs,pdgdau);
  if(minv2>fMinMass*fMinMass && minv2<fMaxMass*fMaxMass){
    Double_t rapid = tmpRD->Y(pdgD);
    if(fAnalysisCuts->IsInFiducialAcceptance(pt,rapid)){
      if(charge>0) fMassVsPtVsYLSpp->Fill(TMath::Sqrt(minv2),pt,rapid);
      else fMassVsPtVsYLSmm->Fill(TMath::Sqrt(minv2),pt,rapid);
    }
  }
  return;
}

//________________________________________________________________________
void AliAnalysisTaskCombinHF::FillGenHistos(TClonesArray* arrayMC){
  // Fill histos with generated quantities
  Int_t totPart=arrayMC->GetEntriesFast();
  Int_t thePDG=411;
  Int_t nProng=3;
  if(fMeson==kDzero){
    thePDG=421;
    nProng=2;
  }
  for(Int_t ip=0; ip<totPart; ip++){
    AliAODMCParticle *part = (AliAODMCParticle*)arrayMC->At(ip);
    if(TMath::Abs(part->GetPdgCode())==thePDG){
      Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,part,fGoUpToQuark);
      fHistCheckOrigin->Fill(orig);
      if(fPromptFeeddown==kFeeddown && orig!=5) continue;
      else if(fPromptFeeddown==kPrompt && orig!=4) continue;
      else if(fPromptFeeddown==kBoth && orig<4) continue;
      fHistCheckOriginSel->Fill(orig);
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
      }
      fHistCheckDecChan->Fill(deca);
      if(labDau[0]==-1){
	//	printf(Form("Meson %d Label of daughters not filled correctly -- %d\n",fMeson,isGoodDecay));
	continue; //protection against unfilled array of labels
      }
      Bool_t isInAcc=CheckAcceptance(arrayMC,nProng,labDau);
      if(isInAcc) fHistCheckDecChanAcc->Fill(deca);
      if(isGoodDecay){
	Double_t ptgen=part->Pt();
	Double_t ygen=part->Y();
	if(fAnalysisCuts->IsInFiducialAcceptance(ptgen,ygen)){
	  fPtVsYGen->Fill(ptgen,ygen);
	  if(TMath::Abs(ygen)<0.5) fPtVsYGenLimAcc->Fill(ptgen,ygen);
	  if(isInAcc) fPtVsYGenAcc->Fill(ptgen,ygen);
	}
      }
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::FillHistos(Int_t pdgD,Int_t nProngs, AliAODRecoDecay* tmpRD, Double_t* px, Double_t* py, Double_t* pz, UInt_t *pdgdau, TClonesArray *arrayMC, Int_t* dgLabels){
  // Fill histos for candidates with proper charge sign

  Bool_t accept=kFALSE;

  tmpRD->SetPxPyPzProngs(nProngs,px,py,pz);
  Double_t pt = tmpRD->Pt();
  Double_t minv2 = tmpRD->InvMass2(nProngs,pdgdau);
  Double_t mass=TMath::Sqrt(minv2);

  if(minv2>fMinMass*fMinMass && minv2<fMaxMass*fMaxMass){
    Double_t rapid = tmpRD->Y(pdgD);
    if(fAnalysisCuts->IsInFiducialAcceptance(pt,rapid)){
      fMassVsPtVsY->Fill(mass,pt,rapid);
      accept=kTRUE;
      if(fReadMC){
	Int_t signPdg[3]={0,0,0};
	for(Int_t iii=0; iii<nProngs; iii++) signPdg[iii]=pdgdau[iii];
	Int_t labD = tmpRD->MatchToMC(pdgD,arrayMC,nProngs,signPdg);
	if(labD>=0){
	  AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(TMath::Abs(dgLabels[0])));
	  if(part){
	    Int_t pdgCode = TMath::Abs( part->GetPdgCode() );
	    if(pdgCode==321){
	      AliAODMCParticle* dmes =  dynamic_cast<AliAODMCParticle*>(arrayMC->At(labD));
	      if(dmes){
		Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,dmes,fGoUpToQuark);
		if((fPromptFeeddown==kFeeddown && orig==5)|| (fPromptFeeddown==kPrompt && orig==4) || (fPromptFeeddown==kBoth && orig>=4)) {
		  fPtVsYReco->Fill(dmes->Pt(),dmes->Y());
		}
	      }
	      fMassVsPtVsYSig->Fill(mass,pt,rapid);
	    }else{
	      fMassVsPtVsYRefl->Fill(mass,pt,rapid);
	    }
	  }
	}else{
	  fMassVsPtVsYBkg->Fill(mass,pt,rapid);
	}
      }
    }
  }

  Int_t nRotated=0;
  Double_t massRot=0;// calculated later only if candidate is acceptable
  Double_t angleProngXY;
  if(TMath::Abs(pdgD)==421)angleProngXY=TMath::ACos((px[0]*px[1]+py[0]*py[1])/TMath::Sqrt((px[0]*px[0]+py[0]*py[0])*(px[1]*px[1]+py[1]*py[1])));
  else {
    angleProngXY=TMath::ACos(((px[0]+px[1])*px[2]+(py[0]+py[1])*py[2])/TMath::Sqrt(((px[0]+px[1])*(px[0]+px[1])+(py[0]+py[1])*(py[0]+py[1]))*(px[2]*px[2]+py[2]*py[2])));
  }
  Double_t ptOrig=pt;


  Double_t rotStep=(fMaxAngleForRot-fMinAngleForRot)/(fNRotations-1); // -1 is to ensure that the last rotation is done with angle=fMaxAngleForRot
  Double_t rotStep3=(fMaxAngleForRot3-fMinAngleForRot3)/(fNRotations-1); // -1 is to ensure that the last rotation is done with angle=fMaxAngleForRot

  for(Int_t irot=0; irot<fNRotations; irot++){
    Double_t phirot=fMinAngleForRot+rotStep*irot;
    Double_t tmpx=px[0];
    Double_t tmpy=py[0];
    Double_t tmpx2=px[2];
    Double_t tmpy2=py[2];
    px[0]=tmpx*TMath::Cos(phirot)-tmpy*TMath::Sin(phirot);
    py[0]=tmpx*TMath::Sin(phirot)+tmpy*TMath::Cos(phirot);
    if(pdgD==411 || pdgD==431){
      Double_t phirot2=fMinAngleForRot3+rotStep3*irot;
      px[2]=tmpx*TMath::Cos(phirot2)-tmpy*TMath::Sin(phirot2);
      py[2]=tmpx*TMath::Sin(phirot2)+tmpy*TMath::Cos(phirot2);
    }
    tmpRD->SetPxPyPzProngs(nProngs,px,py,pz);
    pt = tmpRD->Pt();
    minv2 = tmpRD->InvMass2(nProngs,pdgdau);
    if(minv2>fMinMass*fMinMass && minv2<fMaxMass*fMaxMass){
      Double_t rapid = tmpRD->Y(pdgD);
      if(fAnalysisCuts->IsInFiducialAcceptance(pt,rapid)){
	massRot=TMath::Sqrt(minv2);
	fMassVsPtVsYRot->Fill(massRot,pt,rapid);
	nRotated++;
	fDeltaMass->Fill(massRot-mass);
	if(fFullAnalysis){
	  Double_t pointRot[5]={mass,massRot-mass,ptOrig,pt-ptOrig,angleProngXY};
	  fDeltaMassFullAnalysis->Fill(pointRot);
	}
      }
    }
    px[0]=tmpx;
    py[0]=tmpy;
    if(pdgD==411 || pdgD==431){
      px[2]=tmpx2;
      py[2]=tmpy2;
    }
  }
  fNormRotated->Fill(nRotated);

  return accept;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::IsTrackSelected(AliAODTrack* track){
  // track selection cuts

  if(track->Charge()==0) return kFALSE;
  if(!(track->TestFilterMask(fFilterMask))) return kFALSE;
  if(!SelectAODTrack(track,fTrackCutsAll)) return kFALSE;    
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::IsKaon(AliAODTrack* track){
  // kaon selection cuts

  if(!fPidHF) return kTRUE;
  Int_t isKaon=fPidHF->MakeRawPid(track,AliPID::kKaon);  
  if(isKaon>=0 && SelectAODTrack(track,fTrackCutsKaon)) return kTRUE;
  return kFALSE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::IsPion(AliAODTrack* track){
  // pion selection cuts

  if(!fPidHF) return kTRUE;
  Int_t isPion=fPidHF->MakeRawPid(track,AliPID::kPion);
  if(isPion>=0&& SelectAODTrack(track,fTrackCutsPion)) return kTRUE;
  return kFALSE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskCombinHF::SelectAODTrack(AliAODTrack *track, AliESDtrackCuts *cuts){
  // AOD track selection

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
  // check if the decay products are in the good eta and pt range
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
void AliAnalysisTaskCombinHF::Terminate(Option_t */*option*/)
{
  // Terminate analysis
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

