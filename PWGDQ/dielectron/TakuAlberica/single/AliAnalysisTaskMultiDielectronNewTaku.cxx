/*************************************************************************
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                        Basic Analysis Task                            //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include "TVector2.h"


#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliTriggerAnalysis.h>


#include <TGeoGlobalMagField.h>
#include "TGeoManager.h"
#include "AliGeomManager.h"
#include <AliMagF.h>


#include "AliDielectronTaku.h"
#include "AliDielectronHistosTaku.h"
#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliDielectronEventCuts.h"

#include "AliAnalysisTaskMultiDielectronNewTaku.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliVVertex.h"
#include "AliESDVZERO.h"
#include "AliEMCALTrack.h"
//#include "AliFlowEventSimple.h"
//#include "AliFlowLYZEventPlane.h"


ClassImp(AliAnalysisTaskMultiDielectronNewTaku)

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectronNewTaku::AliAnalysisTaskMultiDielectronNewTaku() :
  AliAnalysisTaskSE(),
  fListDielectron(),
  fListHistos(),
  fListTree(NULL),
  fListCF(),
  fTree(NULL),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kMB),
  fTriggerOnV0AND(kFALSE),
  fRejectPileup(kFALSE),
  fTriggerAnalysis(0x0),
  fEventFilter(0x0),
  fCutsEvent(0x0),
  fCutsMother(0x0),
  fEventStat(0x0),
  fEvent(0x0),
  fdEdXvsPt(0x0),
  fdEdXnSigmaElecvsPt(0x0),
  fTOFbetavsPt(0x0),
  fTOFnSigmaElecvsPt(0x0),
  fTPCcrossedRowsvsPt(0x0),
  fTPCchi2vsPt(0x0),
  fNEvent(0),
  fkTriggerMask(0),
  fkTriggerCent(0),
  fkNCut(0),
  fkRunNumber(0),
  fkCentrality(0),
  fkXvPrim(0),
  fkYvPrim(0),
  fkZvPrim(0),
  fkXRes(0),
  fkYRes(0),
  fkZRes(0),
  fkNTrk(0),
  fkTracks(0),
  fkNacc(0),
  fkNaccTrcklts(0),
  fkNch(0),
  fkZDCN1E(0),
  fkZDCP1E(0),
  fkZDCN2E(0),
  fkZDCP2E(0),
  fkV0A(0),
  fkV0C(0),
  fkNPar(0),
  //  fFlowEvent(0),
  //  fLyzEp(0),
  fQsum(NULL),
  fQ2sum(0),
  fkTriggerInfo(0),
  fMag(0)
{
  //
  // Constructor
  //


}

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectronNewTaku::AliAnalysisTaskMultiDielectronNewTaku(const char *name,
								     AliDielectronEventCuts* cutsEvent
								     ) :
  AliAnalysisTaskSE(name),
  fListDielectron(),
  fListHistos(),
  fListTree(NULL),
  fListCF(),
  fTree(NULL),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kMB),
  fTriggerOnV0AND(kFALSE),
  fRejectPileup(kFALSE),
  fTriggerAnalysis(0x0),
  fEventFilter(0x0),
  fCutsEvent(cutsEvent),
  fCutsMother(0x0),
  fEventStat(0x0),
  fEvent(0x0),
  fdEdXvsPt(0x0),
  fdEdXnSigmaElecvsPt(0x0),
  fTOFbetavsPt(0x0),
  fTOFnSigmaElecvsPt(0x0),
  fTPCcrossedRowsvsPt(0x0),
  fTPCchi2vsPt(0x0),
  fNEvent(0),
  fkTriggerMask(0),
  fkTriggerCent(0),
  fkNCut(0),
  fkRunNumber(0),
  fkCentrality(0),
  fkXvPrim(0),
  fkYvPrim(0),
  fkZvPrim(0),
  fkXRes(0),
  fkYRes(0),
  fkZRes(0),
  fkNTrk(0),
  fkTracks(0),
  fkNacc(0),
  fkNaccTrcklts(0),
  fkNch(0),
  fkZDCN1E(0),
  fkZDCP1E(0),
  fkZDCN2E(0),
  fkZDCP2E(0),
  fkV0A(0),
  fkV0C(0),
  fkNPar(0),
  //  fFlowEvent(0),
  //  fLyzEp(0),
  fQsum(NULL),
  fQ2sum(0),
  fkTriggerInfo(0),
  fMag(0)
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
  //DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TH1D::Class());
  DefineOutput(4, TTree::Class());

  fName  = name;
  // Constructor.
  fQsum = new TVector2();        // flow vector sum
  fQ2sum = 0;


  cout<<" ************** AliAnalysisTaskMultiDielectron::AliAnalysisTaskMultiDielectron  ***********"<<endl;
}


//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronNewTaku::UserCreateOutputObjects()
{
  cout<<" ************** AliAnalysisTaskMultiDielectron::Start ***********"<<endl;
  fListHistos.SetName("Dielectron_Histos_Multi");
  fListTree.SetName("Dielectron_Trees_Multi");
  fListCF.SetName("Dielectron_CF_Multi");
  fListDielectron.SetOwner();
  fListHistos.SetOwner();
  fListTree.SetOwner();
  fListCF.SetOwner();

  //
  // Add all histogram manager histogram lists to the output TList
  //

  if (!fListHistos.IsEmpty()||!fListCF.IsEmpty()) return; //already initialised

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
//   Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  TIter nextDie(&fListDielectron);
  AliDielectronTaku *die=0;
  while ( (die=static_cast<AliDielectronTaku*>(nextDie())) ){
    die->Init();
    if (die->GetHistogramList()) fListHistos.Add(const_cast<THashList*>(die->GetHistogramList()));
    if (die->GetCFManagerPair()) fListCF.Add(const_cast<AliCFContainer*>(die->GetCFManagerPair()->GetContainer()));
    //if (die->GetTreeList()) fListTree.Add(const_cast<THashList*>(die->GetTreeList()));
    ///if (die->GetTreeList()) fTree = (TTree*)(die->GetTreeManager()->GetSingleTree());
  }

  Int_t cuts=fListDielectron.GetEntries();
  Int_t nbins=kNbinsEvent+2*cuts;
  if (!fEventStat){
    
    fEventStat=new TH1D(Form("hEventStat_%s",fName.Data()),"Event statistics",nbins,0,nbins);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");

    //default names
    fEventStat->GetXaxis()->SetBinLabel(3,"Bin3 not used");
    fEventStat->GetXaxis()->SetBinLabel(4,"Bin4 not used");
    fEventStat->GetXaxis()->SetBinLabel(5,"Bin5 not used");
    
    if (fTriggerOnV0AND&&isESD) fEventStat->GetXaxis()->SetBinLabel(3,"V0and triggers");
    if (fEventFilter) fEventStat->GetXaxis()->SetBinLabel(4,"After Event Filter");
    if (fRejectPileup) fEventStat->GetXaxis()->SetBinLabel(5,"After Pileup rejection");
    
    for (Int_t i=0; i<cuts; ++i){
      fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+1)+2*i,Form("#splitline{1 candidate}{%s}",fListDielectron.At(i)->GetName()));
      fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+2)+2*i,Form("#splitline{With >1 candidate}{%s}",fListDielectron.At(i)->GetName()));
    }
  }

  if (!fTriggerAnalysis) fTriggerAnalysis=new AliTriggerAnalysis;
  fTriggerAnalysis->EnableHistograms();
  fTriggerAnalysis->SetAnalyzeMC(AliDielectronMC::Instance()->HasMC());


  //// this is just test
  //// I would like to see dE/dx, pT, and so on...

  TList *tQAElectron = new TList();
  tQAElectron->SetName("QAElectron");
  tQAElectron->SetOwner();


  int nbinx=400;
  float max_x=20;
  float min_x=0.2;
  float binw = (TMath::Log(max_x)-TMath::Log(min_x))/nbinx;
  double xbin[401];
  for(int ii=0;ii<nbinx+1;ii++){
    xbin[ii] = TMath::Exp(TMath::Log(min_x) + 0.5*binw+binw*ii);
  }


  //// my histo 
  fEvent = new TH1D("Event","Number of Events",   60,0,60);
  tQAElectron->Add(fEvent);
  //fdEdXvsPt = new TH2D("dEdXvsPt","dE/dX vs. PT of TPC", 400,0, 20, 2000,0,200);
  fdEdXvsPt = new TH2D("dEdXvsPt","dE/dX vs. PT of TPC", nbinx, xbin, 2000,0,200);
  tQAElectron->Add(fdEdXvsPt);
  fdEdXnSigmaElecvsPt = new TH2D("fdEdXnSigmaElecvsPt"," dE/dX normalized to electron vs. pT of TPC",
				 //				 400, 0, 20, 2000, -10, 10);
				 nbinx, xbin, 2000, -10, 10);
  tQAElectron->Add(fdEdXnSigmaElecvsPt);
  fTOFbetavsPt = new TH2D("fTOFbetavsPt","TOF beta vs. p", 400, 0, 20, 1200, 0, 1.2);
  tQAElectron->Add(fTOFbetavsPt);
  fTOFnSigmaElecvsPt = new TH2D("fTOFnSigmaElecvsPt","TOF nsigma for electron", 400, 0, 20, 2000, -10, 10);
  tQAElectron->Add(fTOFnSigmaElecvsPt);

  fTPCcrossedRowsvsPt = new TH2D("fTTPCcrossedRowsvsPt","TPC crossed rows", 400, 0, 20, 160,0,160);
  tQAElectron->Add(fTPCcrossedRowsvsPt);
  fTPCchi2vsPt = new TH2D("fTTPCchi2RowsvsPt","TPC chi2", 400, 0, 20, 1000,0,200);
  tQAElectron->Add(fTPCchi2vsPt);


  fListHistos.Add(tQAElectron);
  
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
  PostData(3, fEventStat);
  //PostData(4, &fListTree);

  //PostData(1,&fListTree);
  fTree = new TTree(Form("tree_%s",fName.Data()),"single");
  SetBranches(fTree);

  PostData(4,fTree);


  ////// basic cuts 
  fCutsMother = new AliESDtrackCuts;
  fCutsMother->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  fCutsMother->SetRequireTPCRefit(kTRUE);
  fCutsMother->SetRequireITSRefit(kTRUE);
  fCutsMother->SetMaxDCAToVertexZ(3.0);
  fCutsMother->SetMaxDCAToVertexXY(1.0);
  fCutsMother->SetEtaRange( -0.9 , 0.9 );
  fCutsMother->SetAcceptKinkDaughters(kFALSE);
  fCutsMother->SetPtRange(0.2,10);
  fCutsMother->SetMinNClustersTPC(70);
  fCutsMother->SetMaxChi2PerClusterTPC(4);
  

  //lee yang zeros event plane
  //fLyzEp = new AliFlowLYZEventPlane() ;
  //fLyzEp-> Init();

  fkTriggerInfo = new TObjArray(100);

  fNEvent = 0;
  fkNPar=0;
  cout<<" ************** AliAnalysisTaskMultiDielectron::UserCreateOutputObjects End ***********"<<endl;
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronNewTaku::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //
  if(fNEvent%100==0){
    cout<<"Processing event "<<fNEvent<<endl;
  }
  fNEvent++;

  //cout<<"Start UserExec"<<endl;

  if (fListHistos.IsEmpty()&&fListCF.IsEmpty()) return;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if (!inputHandler) return;

  if ( inputHandler->GetPIDResponse() ){
    AliDielectronVarManager::SetPIDResponse( inputHandler->GetPIDResponse() );
  } else {
    //load esd pid bethe bloch parameters depending on the existance of the MC handler
    // yes: MC parameters
    // no:  data parameters
    //ESD case
    if (isESD){
      if (!AliDielectronVarManager::GetESDpid()){
        if (AliDielectronMC::Instance()->HasMC()) {
          AliDielectronVarManager::InitESDpid();
        } else {
          AliDielectronVarManager::InitESDpid(1);
	  /*
	  Double_t fAlephParam[5]={2.11543/50,
				   20.3394,
				   5.0411e-11,
				   2.15543,
				   2.88663};
	  */
	  Double_t fAlephParam[5]={2.11543/122,
				   42.3394,
				   2.0411e-22,
				   2.25543,
				   6.89
	  };
                                   
                                   
                                   


	  AliESDpid *fESDpid = new AliESDpid();
	  fESDpid->GetTPCResponse().SetBetheBlochParameters(fAlephParam[0],
							    fAlephParam[1],
							    fAlephParam[2],
							    fAlephParam[3],
							    fAlephParam[4]);
	  AliDielectronVarManager::SetESDpid(fESDpid);
        }
      }
    }
    //AOD case
    if (isAOD){
      if (!AliDielectronVarManager::GetAODpidUtil()){
        if (AliDielectronMC::Instance()->HasMC()) {
          AliDielectronVarManager::InitAODpidUtil();
        } else {
          AliDielectronVarManager::InitAODpidUtil(1);
        }
      }
    }
  } 



  /*
  AliESDInputHandler *esdHandler=0x0;
  if ( (esdHandler=dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler())) && esdHandler->GetESDpid() ){
    AliDielectronVarManager::SetESDpid(esdHandler->GetESDpid());
  } else {
    //load esd pid bethe bloch parameters depending on the existance of the MC handler
    // yes: MC parameters
    // no:  data parameters
    if (!AliDielectronVarManager::GetESDpid()){
      if (AliDielectronMC::Instance()->HasMC()) {
        AliDielectronVarManager::InitESDpid();
      } else {
        AliDielectronVarManager::InitESDpid(1);
      }
    }
  } 
  */

  //////////////////////////////////////////
  /////////////////////////////////////////
  // just copy from hor
  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());

  AliESDRun *fRun = (AliESDRun*)fESD->GetESDRun();
  fMag = fRun->GetMagneticField();
  if(fRun){
    for(int itr=0; itr<AliESDRun::kNTriggerClasses;itr++){
      fkTrigName = new TObjString();
      fkTrigName->SetString(fRun->GetTriggerClass(itr));
      //cout<<itr<<" "<<triggername[itr]<<endl;
      fkTriggerInfo->AddAt(fkTrigName, itr);
      fkTrigName->Clear();
    }
  }



  //just dump the fired trigger class
  if(fESD){
    AliESDHeader* head = (AliESDHeader*)fESD->GetHeader();
    //cout<<" trigger "<<head->GetFiredTriggerInputs()<<" : "<<head->GetTriggerMask()<<" : "<<fESD->GetFiredTriggerClasses()<<" "<<endl;
    fkTriggerMask = head->GetTriggerMask();
    
    for(int itr=0; itr<AliESDRun::kNTriggerClasses;itr++){
      if( (head->GetTriggerMask() >> itr) & 0x1 == 1 ){
	fEvent->Fill(10+itr);
      }
    }
  }

  if(fESD) {
    fEvent->Fill(0);
    Bool_t  isEvT = kFALSE;
    /*
    isEvT = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() 
	     & AliVEvent::kMB);
    */
    isEvT = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() 
	     //& AliVEvent::kSemiCentral);
	     & AliVEvent::kCentral);

    //AliTriggerAnalysis *fTrigAna = new AliTriggerAnalysis();
    if (isEvT){
      fEvent->Fill(1);
      if(fCutsEvent->IsSelected(fESD)){
        fEvent->Fill(2);
        const AliESDVertex *vertex = fESD->GetPrimaryVertex(); 
        if((vertex->GetNContributors()>2)&&(TMath::Abs(vertex->GetZ())<7)) {
          fEvent->Fill(3);
	}
      }
    }

    //////////////////////////////////////////
    //// fill the number of fired triggers for each trigger types
    /*
    isEvT = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
    //AliTriggerAnalysis *fTrigAna = new AliTriggerAnalysis();
    if (isEvT){
      if(isEvT & AliVEvent::kMB) fEvent->Fill(10);
      if(isEvT & AliVEvent::kINT7) fEvent->Fill(11);
      if(isEvT & AliVEvent::kMUON) fEvent->Fill(12);
      if(isEvT & AliVEvent::kHighMult) fEvent->Fill(13);
      if(isEvT & AliVEvent::kEMC1) fEvent->Fill(14);
      if(isEvT & AliVEvent::kCINT5) fEvent->Fill(15);
      if(isEvT & AliVEvent::kMUSPB) fEvent->Fill(16);
      if(isEvT & AliVEvent::kMUSHPB) fEvent->Fill(17);
      if(isEvT & AliVEvent::kMuonLikePB) fEvent->Fill(18);
      if(isEvT & AliVEvent::kMuonUnlikePB) fEvent->Fill(19);
      if(isEvT & AliVEvent::kEMC7) fEvent->Fill(20);
      if(isEvT & AliVEvent::kMUS7) fEvent->Fill(21);
      if(isEvT & AliVEvent::kPHI1) fEvent->Fill(22);
      if(isEvT & AliVEvent::kPHOSPb) fEvent->Fill(23);
      if(isEvT & AliVEvent::kEMCEJE) fEvent->Fill(24);
      if(isEvT & AliVEvent::kEMCEGA) fEvent->Fill(25);
      if(isEvT & AliVEvent::kCentral) fEvent->Fill(26);
      if(isEvT & AliVEvent::kSemiCentral) fEvent->Fill(27);
      if(isEvT & AliVEvent::kDG5) fEvent->Fill(28);
      if(isEvT & AliVEvent::kZED) fEvent->Fill(29);
      if(isEvT & AliVEvent::kAny) fEvent->Fill(30);
    }
    */
  }
  //////////////////////////////////////////
  //////////////////////////////////////////


  // Was event selected ?
  UInt_t isSelected = AliVEvent::kAny;
  if( fSelectPhysics && inputHandler && inputHandler->GetEventSelection() ) {
    isSelected = inputHandler->IsEventSelected();
    isSelected&=fTriggerMask;
  }
  
  //Before physics selection
  fEventStat->Fill(kAllEvents);
  if (isSelected==0) {
    PostData(3,fEventStat);
    return;
  }
  //after physics selection
  fEventStat->Fill(kSelectedEvents);

  //V0and
  if (fTriggerOnV0AND&&isESD){
    //if (!fTriggerAnalysis->IsOfflineTriggerFired(static_cast<AliESDEvent*>(InputEvent()), AliTriggerAnalysis::kV0AND)) return;
  }
  fEventStat->Fill(kV0andEvents);
  
  //event filter
  if (fEventFilter) {
    //if (!fEventFilter->IsSelected(InputEvent())) return;
  }
  fEventStat->Fill(kFilteredEvents);
  
  //pileup
  if (fRejectPileup){
    if (InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.)) return;
  }
  fEventStat->Fill(kPileupEvents);


  Bool_t isEvT1 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() 
		   & AliVEvent::kCentral);

  Bool_t isEvT2 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() 
		   & AliVEvent::kSemiCentral);


  Bool_t isEvT3 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() 
		   & AliVEvent::kMB);
  
  fkTriggerCent = 0;
  if(!isEvT1 && !isEvT2 && !isEvT3){
    return ;
  }else{
    if(isEvT1){
      fkTriggerCent += 100;
    }
    if(isEvT2){
      fkTriggerCent += 10;
    }
    if(isEvT3){
      fkTriggerCent += 1;
    }
  }

  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );

  AliDielectronPID::SetCorrVal((Double_t)InputEvent()->GetRunNumber());

  FillEvent(InputEvent());
  fkNPar = 0;

  //Process event in all AliDielectron instances
  TIter nextDie(&fListDielectron);
  AliDielectronTaku *die=0;
  Int_t idie=0;
  while ( (die=static_cast<AliDielectronTaku*>(nextDie())) ){
    die->Process(InputEvent());
    if (die->HasCandidates()){
      Int_t ncandidates=die->GetPairArray(1)->GetEntriesFast();
      if (ncandidates==1) fEventStat->Fill((kNbinsEvent)+2*idie);
      else if (ncandidates>1) fEventStat->Fill((kNbinsEvent+1)+2*idie);
    }
    AliDielectronVarManager::Fill(InputEvent(), fgValues);
    //Fill track information, separately for the track array candidates                                                         
    for (Int_t ii=0; ii<2; ++ii){ 
      TObjArray *obj = (TObjArray*)die->GetTrackArray(ii);
      Int_t ntracks=obj->GetEntriesFast();
      for (Int_t itrack=0; itrack<ntracks; ++itrack){
	////////////////////////////////////////////////////////////////////
	AliDielectronVarManager::Fill(obj->UncheckedAt(itrack), fgValues);
	////////////////////////////////////////////////////////////////////
	AliVTrack *trk = static_cast<AliVTrack*>(obj->UncheckedAt(itrack));
	double par[3] ;
	MomentumEnergyMatch(trk, par);
	fgValues[AliDielectronVarManager::kE] = par[0];
	fgValues[AliDielectronVarManager::kDeltaEta] = par[1];
	fgValues[AliDielectronVarManager::kDeltaPhi] = par[2];
	/*
	if(par[0]<0){
	  continue;
	}
	*/
	
	AliESDtrack *esdtrack = static_cast<AliESDtrack*>(obj->UncheckedAt(itrack));
	Double_t dca[2]={-999.,-999.};
	Double_t cov[3]={-999.,-999.,-999.};
	Double_t kBeampiperadius=3.;
	esdtrack->PropagateToDCA(InputEvent()->GetPrimaryVertex(), 
				 InputEvent()->GetMagneticField(), kBeampiperadius, dca, cov);
	fgValues[AliDielectronVarManager::kLegDist] = dca[1];
	fgValues[AliDielectronVarManager::kLegDistXY] = dca[0];
	fgValues[AliDielectronVarManager::kNclsTPC] = esdtrack->GetTPCCrossedRows();
	
	/*
	///// KF analysis 
	AliKFParticle kfTrack = AliKFParticle(*trk, 11); //assuming electron //charge is from trk
	/// i would like to store
	/// X, Y, Z, Px, Py, Pz, S(decay length/mom), Chi2, NDF
	if(kfTrack.GetNDF()!=0) fgValues[AliDielectronVarManager::kChi2NDF] = kfTrack.GetChi2()/kfTrack.GetNDF();
	fgValues[AliDielectronVarManager::kDecayLength] = kfTrack.GetS();
	fgValues[AliDielectronVarManager::kR] = kfTrack.GetR();
	fgValues[AliDielectronVarManager::kThetaHE] = kfTrack.GetX();
	fgValues[AliDielectronVarManager::kPhiHE] = kfTrack.GetY();
	fgValues[AliDielectronVarManager::kThetaCS] = kfTrack.GetZ();
	fgValues[AliDielectronVarManager::kPhiCS] = kfTrack.GetPhi();
	fgValues[AliDielectronVarManager::kITSsignalSSD1] = kfTrack.GetPx();
	fgValues[AliDielectronVarManager::kITSsignalSSD2] = kfTrack.GetPy();
	fgValues[AliDielectronVarManager::kITSsignalSDD1] = kfTrack.GetPz();
	fgValues[AliDielectronVarManager::kITSsignalSDD2] = kfTrack.GetP();
	*/
	////////////////////////////////////////////////////////////////////
	for(int ich=0;ich<AliDielectronVarManager::kNMaxValues; ich++){
	  fgData[ich][fkNPar] = fgValues[ich];
	}
	fkNPar++;
      }
    }
    fkNCut = idie;
    ++idie;
    fTree->Fill();
  }

  //////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //// fillback to single particle without trackcuts
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    if(!fCutsMother->AcceptTrack(track)) continue;
    fdEdXvsPt->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
    fdEdXnSigmaElecvsPt->Fill(track->GetTPCmomentum(), 
			      AliDielectronVarManager::GetESDpid()->NumberOfSigmasTPC(track,
										      AliPID::kElectron)
			      -AliDielectronPID::GetCorrVal());

    fTPCcrossedRowsvsPt->Fill(track->GetTPCmomentum(), track->GetTPCCrossedRows());
    fTPCchi2vsPt->Fill(track->GetTPCmomentum(), track->GetTPCchi2());

    /// for beta caliculaton 
    Double_t l = track->GetIntegratedLength();  // cm
    Double_t t = track->GetTOFsignal();
    Double_t t0 = AliDielectronVarManager::GetESDpid()->GetTOFResponse().GetTimeZero(); // ps
    Double_t beta = 0;
    if( (l < 360. || l > 800.) || (t <= 0.) || (t0 >999990.0) ) {
      beta;
    }
    else {
      t -= t0; // subtract the T0
      l *= 0.01;  // cm ->m
      t *= 1e-12; //ps -> s
    
      Double_t v = l / t;
      beta = v / TMath::C();
    }

    fTOFbetavsPt->Fill(track->GetTPCmomentum(), beta);
    //// electron cuts //////
    if(fabs(AliDielectronVarManager::GetESDpid()->NumberOfSigmasTPC(track,AliPID::kElectron)
	    -AliDielectronPID::GetCorrVal())<3){
      fTOFnSigmaElecvsPt->Fill(track->GetTPCmomentum(), 
			       AliDielectronVarManager::GetESDpid()->NumberOfSigmasTOF(track,
										       AliPID::kElectron));
    }
  }
  


  //cout<<idie<<endl;
  
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
  PostData(3,fEventStat);
  //PostData(4, &fListTree);
  //PostData(1, &fListTree);
  PostData(4, fTree);
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronNewTaku::FinishTaskOutput()
{
  //
  // Write debug tree
  //
  cout<<" ******* AliAnalysisTaskMultiDielectronNewTaku::FinishTaskOutput() *****"<<endl;
  TIter nextDie(&fListDielectron);
  AliDielectronTaku *die=0;
  while ( (die=static_cast<AliDielectronTaku*>(nextDie())) ){
    die->SaveDebugTree();
  }
}

//_________________________________________________________________________________              
void AliAnalysisTaskMultiDielectronNewTaku::MomentumEnergyMatch(const AliVParticle *track, double *par){

  Float_t  clsPos[3];                                                                                          
  Double_t trkPos[3];  
  Double_t matchclsE = -9999.9;               
  
  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);                                      
  AliESDEvent *evt = (AliESDEvent*)esdtrack->GetESDEvent();              
  Double_t  magF = evt->GetMagneticField();         
  Double_t magSign = 1.0;                                                                                    
  if(magF<0)magSign = -1.0;           
  if (!TGeoGlobalMagField::Instance()->GetField()) {   
    printf("Loading field map...\n");                                                                    
    //AliMagF* field = new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG);                                 
    AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG); // for 10d              
    TGeoGlobalMagField::Instance()->SetField(field);                                                      
  }                                                   
                                                                                                            
  AliEMCALTrack *emctrack = new AliEMCALTrack(*esdtrack);                                                 
  Double_t fieldB[3];                                                                                        
  emctrack->GetBxByBz(fieldB);                                                                               
  //printf("%g %g %g \n", fieldB[0], fieldB[1], fieldB[2]);                                                    
  double min_r=99999.0;
  double min_dphi=-9999.0;
  double min_deta=-9999.0;
  
  for(Int_t icl=0; icl<evt->GetNumberOfCaloClusters(); icl++){
    AliVCluster *cluster = (AliVCluster*) evt->GetCaloCluster(icl);                                          
    if(!cluster->IsEMCAL()) continue;          
    cluster->GetPosition(clsPos);                                                                            
    if(!emctrack->PropagateToGlobal(clsPos[0],clsPos[1],clsPos[2],0.,0.) )  continue;                        
    emctrack->GetXYZ(trkPos);          
    TVector3 clsPosVec(clsPos[0],clsPos[1],clsPos[2]);                                                       
    TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);   
    Double_t delEmcphi = clsPosVec.Phi()-trkPosVec.Phi();  // track cluster matching                         
    Double_t delEmceta = clsPosVec.Eta()-trkPosVec.Eta();  // track cluster matching      
    double rmatch = sqrt(pow(delEmcphi,2)+pow(delEmceta,2));                                                 
    /*
    if(rmatch<min_r){
      min_r = rmatch;
      min_dphi = delEmcphi;     
      min_deta = delEmceta;
      matchclsE = cluster->E();
    }
    */
    if(rmatch<0.02 && rmatch<min_r){
      min_r = rmatch;
      min_dphi = delEmcphi; 
      min_deta = delEmceta; 
      matchclsE = cluster->E();
    }
  }
  delete emctrack;        

  par[0] = matchclsE;
  par[1] = min_dphi;
  par[2] = min_deta;

}


//_________________________________________________________________________________              
void AliAnalysisTaskMultiDielectronNewTaku::FillEvent(AliVEvent * const ev){

  //AliKFVertex *fgKFVertex = new fgKFVertex();
  //if (ev && ev->GetPrimaryVertex()) fgKFVertex=new AliKFVertex(*ev->GetPrimaryVertex());

  if(ev){
    fkRunNumber = ev->GetRunNumber();
    fkXvPrim = ev->GetPrimaryVertex()->GetX();
    fkYvPrim = ev->GetPrimaryVertex()->GetY();
    fkZvPrim = ev->GetPrimaryVertex()->GetZ();
    fkNTrk = ev->GetNumberOfTracks();
    fkNacc = AliDielectronHelper::GetNacc(ev);
    fkNaccTrcklts = AliDielectronHelper::GetNaccTrcklts(ev);
    
    fkZDCN1E = ev->GetZDCN1Energy();
    fkZDCP1E = ev->GetZDCP1Energy();
    fkZDCN2E = ev->GetZDCN2Energy();
    fkZDCP2E = ev->GetZDCP2Energy();
  }
  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(ev);

  fkCentrality=-1;
  AliCentrality *esdCentrality = const_cast<AliESDEvent*>(fESD)->GetCentrality();
  if (esdCentrality) fkCentrality = esdCentrality->GetCentralityPercentile("V0M");

  if(fESD->GetPrimaryVertex()){
    fkXRes = fESD->GetPrimaryVertex()->GetXRes();
    fkYRes = fESD->GetPrimaryVertex()->GetYRes();
    fkZRes = fESD->GetPrimaryVertex()->GetZRes();
  }
  
  if(fESD->GetVZEROData()){
    fkV0A=0;
    fkV0C=0;
    for(int ich=0;ich<32;ich++){
      fkV0A += fESD->GetVZEROData()->GetMultiplicityV0A(ich);
      fkV0C += fESD->GetVZEROData()->GetMultiplicityV0C(ich);
    }
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  /////////////// get event plane ///////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  AliEventplane *fEventPlane = (AliEventplane*)fESD->GetEventplane();
  if (fEventPlane)
    fkRP = fEventPlane->GetEventplane("V0", fESD,2);

  /*
  fkRP=-999; fkRPQx=-999; fkRPQy=-999;
  fkRPsub1=-999; fkRPsub1Qx=-999; fkRPsub1Qy=-999;
  fkRPsub2=-999; fkRPsub2Qx=-999; fkRPsub2Qy=-999;

  if(fEventPlane){

    fkRP = fEventPlane->GetEventplane("Q"); 
    if(fEventPlane->GetQVector()){
      fkRPQx = fEventPlane->GetQVector()->X();
      fkRPQy = fEventPlane->GetQVector()->Y();
    }
    if(fEventPlane->GetQsub1()){
      fkRPsub1 = fEventPlane->GetQsub1()->Phi()/2;
      fkRPsub1Qx = fEventPlane->GetQsub1()->X();
      fkRPsub1Qy = fEventPlane->GetQsub1()->Y();
    }
    if(fEventPlane->GetQsub2()){
      fkRPsub2 = fEventPlane->GetQsub2()->Phi()/2;
      fkRPsub2Qx = fEventPlane->GetQsub2()->X();
      fkRPsub2Qy = fEventPlane->GetQsub2()->Y();
    }

    fkQsubRes = fEventPlane->GetQsubRes();
     cout<<fkRP<<" "<<fkRPQx<<" "<<fkRPQy<<" "<<fkRPsub1<<" "<<fkRPsub1Qx<<" "<<fkRPsub1Qy<<" "
	<<fkRPsub2<<" "<<fkRPsub2Qx<<" "<<fkRPsub2Qy<<endl;
 
    cout<<fEventPlane->GetEventplane("Q")<<endl;
    TVector2 *vec2 = (TVector2*)fEventPlane->GetQVector();
    if(!vec2){ cout<<" no EP information "<<endl;}
    else{
      Double_t dRPEP = vec2->Phi()/2; 
      cout<<dRPEP<<endl;
    }
 

  }
 
  fFlowEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));
  
  if (fFlowEvent) {
    //get the Q vector from the FlowEvent
    AliFlowVector vQ = fFlowEvent->GetQ(); 
    //if (vQ.X()== 0. && vQ.Y()== 0. ) { cout<<"Q vector is NULL!"<<endl; } //coding violation
    //Weight with the multiplicity
    Double_t dQX = 0.;
    Double_t dQY = 0.;
    if (TMath::AreEqualAbs(vQ.GetMult(),0.0,1e-10)) {
      dQX = vQ.X()/vQ.GetMult();
      dQY = vQ.Y()/vQ.GetMult();
    } else {cerr<<"vQ.GetMult() is zero!"<<endl; }
    vQ.Set(dQX,dQY);
    *fQsum += vQ;
    fQ2sum += vQ.Mod2();
    fLyzEp->CalculateRPandW(vQ);
    Double_t dWR = fLyzEp->GetWR();     
    
    Double_t dRP = fLyzEp->GetPsi();
    //plot difference between event plane from EP-method and LYZ-method
    Double_t dRPEP = vQ.Phi()/2;                              //gives distribution from (0 to pi)
    //Double_t dRPEP = 0.5*TMath::ATan2(vQ.Y(),vQ.X());       //gives distribution from (-pi/2 to pi/2)

    Double_t dDeltaPhi = dRPEP - dRP;
    if (dDeltaPhi < 0.) { dDeltaPhi += TMath::Pi(); }        //to shift distribution from (-pi/2 to pi/2) to (0 to pi)

    cout<<dRPEP<<" "<<dRP<<" "<<endl;

  }else{
    cout<<" no fFlowEvent "<<endl;
  }
  */





}


//_________________________________________________________________________________              
void AliAnalysisTaskMultiDielectronNewTaku::SetBranches(TTree *t){

  t->Branch("kNEvent",&fNEvent,"kNEvent/I"); 
  t->Branch("kMag",&fMag,"kMag/D"); 
  t->Branch("fkTriggerInfo","TObjArray",&fkTriggerInfo);
  t->Branch("kTriggerMask",&fkTriggerMask,"kTriggerMask/D"); 
  t->Branch("kTriggerCent",&fkTriggerCent,"kTriggerCent/I"); 
  t->Branch("fkNCut", &fkNCut, "fkNCut/D");
  t->Branch("fkRunNumber", &fkRunNumber, "fkRunNumber/D");
  t->Branch("fkCentrality", &fkCentrality, "fkCentrality/D");
  t->Branch("fkXvPrim", &fkXvPrim, "fkXvPrim/D");
  t->Branch("fkYvPrim", &fkYvPrim, "fkYvPrim/D");
  t->Branch("fkZvPrim", &fkZvPrim, "fkZvPrim/D");
  t->Branch("fkXRes", &fkXRes, "fkXRes/D");
  t->Branch("fkYRes", &fkYRes, "fkYRes/D");
  t->Branch("fkZRes", &fkZRes, "fkZRes/D");
  t->Branch("fkNTrk", &fkNTrk, "fkNTrk/D");
  t->Branch("fkTracks", &fkTracks, "fkTracks/D");
  t->Branch("fkNacc", &fkNacc, "fkNacc/D");
  t->Branch("fkNaccTrcklts", &fkNaccTrcklts, "fkNaccTrcklts/D");
  t->Branch("fkNch", &fkNch, "fkNch/D");
  t->Branch("fkZDCN1E", &fkZDCN1E, "fkZDCN1E/D");
  t->Branch("fkZDCP1E", &fkZDCP1E, "fkZDCP1E/D");
  t->Branch("fkZDCN2E", &fkZDCN2E, "fkZDCN2E/D");
  t->Branch("fkZDCP2E", &fkZDCP2E, "fkZDCP2E/D");
  t->Branch("fkV0A", &fkV0A, "fkV0A/D");
  t->Branch("fkV0C", &fkV0C, "fkV0C/D");
  
  t->Branch("fkRP",&fkRP,"fkRP/D");
  t->Branch("fkRPQx",&fkRPQx,"fkRPQx/D");
  t->Branch("fkRPQy",&fkRPQy,"fkRPQy/D");
  t->Branch("fkRPsub1",&fkRPsub1,"fkRPsub1/D");
  t->Branch("fkRPsub1Qx",&fkRPsub1Qx,"fkRPsub1Qx/D");
  t->Branch("fkRPsub1Qy",&fkRPsub1Qy,"fkRPsub1Qy/D");
  t->Branch("fkRPsub2",&fkRPsub2,"fkRPsub2/D");
  t->Branch("fkRPsub2Qx",&fkRPsub2Qx,"fkRPsub2Qx/D");
  t->Branch("fkRPsub2Qy",&fkRPsub2Qy,"fkRPsub2Qy/D");
  t->Branch("fkQsubRes",&fkQsubRes,"fkQsubRes/D");


  t->Branch("fkNPar", &fkNPar, "fkNPar/I");
  t->Branch("kPx",fgData[AliDielectronVarManager::kPx],"kPx[fkNPar]/D"); 
  t->Branch("kPy",fgData[AliDielectronVarManager::kPy],"kPy[fkNPar]/D"); 
  t->Branch("kPz",fgData[AliDielectronVarManager::kPz],"kPz[fkNPar]/D"); 
  t->Branch("kPt",fgData[AliDielectronVarManager::kPt],"kPt[fkNPar]/D"); 
  t->Branch("kP",fgData[AliDielectronVarManager::kP],"kP[fkNPar]/D"); 
  t->Branch("kXv",fgData[AliDielectronVarManager::kXv],"kXv[fkNPar]/D"); 
  t->Branch("kYv",fgData[AliDielectronVarManager::kYv],"kYv[fkNPar]/D"); 
  t->Branch("kZv",fgData[AliDielectronVarManager::kZv],"kZv[fkNPar]/D"); 
  t->Branch("kOneOverPt",fgData[AliDielectronVarManager::kOneOverPt],"kOneOverPt[fkNPar]/D"); 
  t->Branch("kPhi",fgData[AliDielectronVarManager::kPhi],"kPhi[fkNPar]/D"); 
  t->Branch("kTheta",fgData[AliDielectronVarManager::kTheta],"kTheta[fkNPar]/D"); 
  t->Branch("kEta",fgData[AliDielectronVarManager::kEta],"kEta[fkNPar]/D"); 
  t->Branch("kY",fgData[AliDielectronVarManager::kY],"kY[fkNPar]/D"); 
  t->Branch("kE",fgData[AliDielectronVarManager::kE],"kE[fkNPar]/D"); 
  t->Branch("kM",fgData[AliDielectronVarManager::kM],"kM[fkNPar]/D"); 
  t->Branch("kCharge",fgData[AliDielectronVarManager::kCharge],"kCharge[fkNPar]/D"); 
  t->Branch("kNclsITS",fgData[AliDielectronVarManager::kNclsITS],"kNclsITS[fkNPar]/D"); 
  t->Branch("kNclsTPC",fgData[AliDielectronVarManager::kNclsTPC],"kNclsTPC[fkNPar]/D"); 
  t->Branch("kNclsTPCiter1",fgData[AliDielectronVarManager::kNclsTPCiter1],"kNclsTPCiter1[fkNPar]/D"); 
  t->Branch("kNFclsTPC",fgData[AliDielectronVarManager::kNFclsTPC],"kNFclsTPC[fkNPar]/D"); 
  t->Branch("kNFclsTPCr",fgData[AliDielectronVarManager::kNFclsTPCr],"kNFclsTPCr[fkNPar]/D"); 
  t->Branch("kNFclsTPCrFrac",fgData[AliDielectronVarManager::kNFclsTPCrFrac],"kNFclsTPCrFrac[fkNPar]/D"); 
  t->Branch("kTPCsignalN",fgData[AliDielectronVarManager::kTPCsignalN],"kTPCsignalN[fkNPar]/D"); 
  t->Branch("kTPCsignalNfrac",fgData[AliDielectronVarManager::kTPCsignalNfrac],"kTPCsignalNfrac[fkNPar]/D"); 
  t->Branch("kTPCchi2Cl",fgData[AliDielectronVarManager::kTPCchi2Cl],"kTPCchi2Cl[fkNPar]/D"); 
  t->Branch("kTrackStatus",fgData[AliDielectronVarManager::kTrackStatus],"kTrackStatus[fkNPar]/D"); 
  t->Branch("kNclsTRD",fgData[AliDielectronVarManager::kNclsTRD],"kNclsTRD[fkNPar]/D"); 
  t->Branch("kTRDntracklets",fgData[AliDielectronVarManager::kTRDntracklets],"kTRDntracklets[fkNPar]/D"); 
  t->Branch("kTRDpidQuality",fgData[AliDielectronVarManager::kTRDpidQuality],"kTRDpidQuality[fkNPar]/D"); 
  t->Branch("kTRDprobEle",fgData[AliDielectronVarManager::kTRDprobEle],"kTRDprobEle[fkNPar]/D"); 
  t->Branch("kTRDprobPio",fgData[AliDielectronVarManager::kTRDprobPio],"kTRDprobPio[fkNPar]/D"); 
  t->Branch("kImpactParXY",fgData[AliDielectronVarManager::kImpactParXY],"kImpactParXY[fkNPar]/D"); 
  t->Branch("kImpactParZ",fgData[AliDielectronVarManager::kImpactParZ],"kImpactParZ[fkNPar]/D"); 
  t->Branch("kTrackLength",fgData[AliDielectronVarManager::kTrackLength],"kTrackLength[fkNPar]/D"); 
  t->Branch("kPdgCode",fgData[AliDielectronVarManager::kPdgCode],"kPdgCode[fkNPar]/D"); 
  t->Branch("kPdgCodeMother",fgData[AliDielectronVarManager::kPdgCodeMother],"kPdgCodeMother[fkNPar]/D"); 
  t->Branch("kPdgCodeGrandMother",fgData[AliDielectronVarManager::kPdgCodeGrandMother],"kPdgCodeGrandMother[fkNPar]/D"); 
  t->Branch("kNumberOfDaughters",fgData[AliDielectronVarManager::kNumberOfDaughters],"kNumberOfDaughters[fkNPar]/D"); 
  t->Branch("kHaveSameMother",fgData[AliDielectronVarManager::kHaveSameMother],"kHaveSameMother[fkNPar]/D"); 
  t->Branch("kIsJpsiPrimary",fgData[AliDielectronVarManager::kIsJpsiPrimary],"kIsJpsiPrimary[fkNPar]/D"); 
  t->Branch("kITSsignal",fgData[AliDielectronVarManager::kITSsignal],"kITSsignal[fkNPar]/D"); 
  t->Branch("kITSsignalSSD1",fgData[AliDielectronVarManager::kITSsignalSSD1],"kITSsignalSSD1[fkNPar]/D"); 
  t->Branch("kITSsignalSSD2",fgData[AliDielectronVarManager::kITSsignalSSD2],"kITSsignalSSD2[fkNPar]/D"); 
  t->Branch("kITSsignalSDD1",fgData[AliDielectronVarManager::kITSsignalSDD1],"kITSsignalSDD1[fkNPar]/D"); 
  t->Branch("kITSsignalSDD2",fgData[AliDielectronVarManager::kITSsignalSDD2],"kITSsignalSDD2[fkNPar]/D"); 
  t->Branch("kITSclusterMap",fgData[AliDielectronVarManager::kITSclusterMap],"kITSclusterMap[fkNPar]/D"); 
  t->Branch("kITSnSigmaEle",fgData[AliDielectronVarManager::kITSnSigmaEle],"kITSnSigmaEle[fkNPar]/D"); 
  t->Branch("kITSnSigmaPio",fgData[AliDielectronVarManager::kITSnSigmaPio],"kITSnSigmaPio[fkNPar]/D"); 
  t->Branch("kITSnSigmaMuo",fgData[AliDielectronVarManager::kITSnSigmaMuo],"kITSnSigmaMuo[fkNPar]/D"); 
  t->Branch("kITSnSigmaKao",fgData[AliDielectronVarManager::kITSnSigmaKao],"kITSnSigmaKao[fkNPar]/D"); 
  t->Branch("kITSnSigmaPro",fgData[AliDielectronVarManager::kITSnSigmaPro],"kITSnSigmaPro[fkNPar]/D"); 
  t->Branch("kPIn",fgData[AliDielectronVarManager::kPIn],"kPIn[fkNPar]/D"); 
  t->Branch("kTPCsignal",fgData[AliDielectronVarManager::kTPCsignal],"kTPCsignal[fkNPar]/D"); 
  t->Branch("kTOFsignal",fgData[AliDielectronVarManager::kTOFsignal],"kTOFsignal[fkNPar]/D"); 
  t->Branch("kTOFbeta",fgData[AliDielectronVarManager::kTOFbeta],"kTOFbeta[fkNPar]/D"); 
  t->Branch("kTPCnSigmaEle",fgData[AliDielectronVarManager::kTPCnSigmaEle],"kTPCnSigmaEle[fkNPar]/D"); 
  t->Branch("kTPCnSigmaPio",fgData[AliDielectronVarManager::kTPCnSigmaPio],"kTPCnSigmaPio[fkNPar]/D"); 
  t->Branch("kTPCnSigmaMuo",fgData[AliDielectronVarManager::kTPCnSigmaMuo],"kTPCnSigmaMuo[fkNPar]/D"); 
  t->Branch("kTPCnSigmaKao",fgData[AliDielectronVarManager::kTPCnSigmaKao],"kTPCnSigmaKao[fkNPar]/D"); 
  t->Branch("kTPCnSigmaPro",fgData[AliDielectronVarManager::kTPCnSigmaPro],"kTPCnSigmaPro[fkNPar]/D"); 
  t->Branch("kTOFnSigmaEle",fgData[AliDielectronVarManager::kTOFnSigmaEle],"kTOFnSigmaEle[fkNPar]/D"); 
  t->Branch("kTOFnSigmaPio",fgData[AliDielectronVarManager::kTOFnSigmaPio],"kTOFnSigmaPio[fkNPar]/D"); 
  t->Branch("kTOFnSigmaMuo",fgData[AliDielectronVarManager::kTOFnSigmaMuo],"kTOFnSigmaMuo[fkNPar]/D"); 
  t->Branch("kTOFnSigmaKao",fgData[AliDielectronVarManager::kTOFnSigmaKao],"kTOFnSigmaKao[fkNPar]/D"); 
  t->Branch("kTOFnSigmaPro",fgData[AliDielectronVarManager::kTOFnSigmaPro],"kTOFnSigmaPro[fkNPar]/D"); 
  t->Branch("kKinkIndex0",fgData[AliDielectronVarManager::kKinkIndex0],"kKinkIndex0[fkNPar]/D"); 
  t->Branch("kChi2NDF",fgData[AliDielectronVarManager::kChi2NDF],"kChi2NDF[fkNPar]/D"); 
  t->Branch("kDecayLength",fgData[AliDielectronVarManager::kDecayLength],"kDecayLength[fkNPar]/D"); 
  t->Branch("kR",fgData[AliDielectronVarManager::kR],"kR[fkNPar]/D"); 
  t->Branch("kOpeningAngle",fgData[AliDielectronVarManager::kOpeningAngle],"kOpeningAngle[fkNPar]/D"); 
  t->Branch("kThetaHE",fgData[AliDielectronVarManager::kThetaHE],"kThetaHE[fkNPar]/D"); 
  t->Branch("kPhiHE",fgData[AliDielectronVarManager::kPhiHE],"kPhiHE[fkNPar]/D"); 
  t->Branch("kThetaCS",fgData[AliDielectronVarManager::kThetaCS],"kThetaCS[fkNPar]/D"); 
  t->Branch("kPhiCS",fgData[AliDielectronVarManager::kPhiCS],"kPhiCS[fkNPar]/D"); 
  t->Branch("kLegDist",fgData[AliDielectronVarManager::kLegDist],"kLegDist[fkNPar]/D"); 
  t->Branch("kLegDistXY",fgData[AliDielectronVarManager::kLegDistXY],"kLegDistXY[fkNPar]/D"); 
  t->Branch("kDeltaEta",fgData[AliDielectronVarManager::kDeltaEta],"kDeltaEta[fkNPar]/D"); 
  t->Branch("kDeltaPhi",fgData[AliDielectronVarManager::kDeltaPhi],"kDeltaPhi[fkNPar]/D"); 
  t->Branch("kMerr",fgData[AliDielectronVarManager::kMerr],"kMerr[fkNPar]/D"); 
  t->Branch("kDCA",fgData[AliDielectronVarManager::kDCA],"kDCA[fkNPar]/D"); 
  t->Branch("kPairType",fgData[AliDielectronVarManager::kPairType],"kPairType[fkNPar]/D"); 
  t->Branch("kPseudoProperTime",fgData[AliDielectronVarManager::kPseudoProperTime],"kPseudoProperTime[fkNPar]/D"); 
  t->Branch("kXvPrim",fgData[AliDielectronVarManager::kXvPrim],"kXvPrim=kPairMax[fkNPar]/D"); 
  t->Branch("kYvPrim",fgData[AliDielectronVarManager::kYvPrim],"kYvPrim[fkNPar]/D"); 
  t->Branch("kZvPrim",fgData[AliDielectronVarManager::kZvPrim],"kZvPrim[fkNPar]/D"); 
  t->Branch("kXRes",fgData[AliDielectronVarManager::kXRes],"kXRes[fkNPar]/D"); 
  t->Branch("kYRes",fgData[AliDielectronVarManager::kYRes],"kYRes[fkNPar]/D"); 
  t->Branch("kZRes",fgData[AliDielectronVarManager::kZRes],"kZRes[fkNPar]/D"); 
  t->Branch("kNTrk",fgData[AliDielectronVarManager::kNTrk],"kNTrk[fkNPar]/D"); 
  t->Branch("kTracks",fgData[AliDielectronVarManager::kTracks],"kTracks[fkNPar]/D"); 
  t->Branch("kNacc",fgData[AliDielectronVarManager::kNacc],"kNacc[fkNPar]/D"); 
  t->Branch("kNaccTrcklts",fgData[AliDielectronVarManager::kNaccTrcklts],"kNaccTrcklts[fkNPar]/D"); 
  t->Branch("kNch",fgData[AliDielectronVarManager::kNch],"kNch[fkNPar]/D"); 
  t->Branch("kCentrality",fgData[AliDielectronVarManager::kCentrality],"kCentrality[fkNPar]/D"); 
  t->Branch("kNevents",fgData[AliDielectronVarManager::kNevents],"kNevents[fkNPar]/D"); 


}

