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

#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliTriggerAnalysis.h>
#include <AliPIDResponse.h>
#include <AliTPCPIDResponse.h>

#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliDielectronMixingHandler.h"
#include "AliAnalysisTaskMultiDielectronTG.h"

#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliLog.h"

#include <vector>
#include <deque>
#include <cstdlib> 
#include <TRandom.h>

const char* kPairClassNamesTG[7] = {
  "unlike",
  "like_pp",
  "like_ee",
  "mixunlike_pe",
  "mixunlike_ep",
  "mixlike_pp",
  "mixlike_ee"
};


ClassImp(AliDielectronSingleTG)
ClassImp(AliAnalysisTaskMultiDielectronTG)

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectronTG::AliAnalysisTaskMultiDielectronTG() :
  AliAnalysisTaskSE(),
  fListDielectron(),
  fListHistos(),
  fListCF(),
  fQAElectron(),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kMB),
  fExcludeTriggerMask(0),
  fTriggerOnV0AND(kFALSE),
  fRejectPileup(kFALSE),
  fTriggerLogic(kAny),
  fTriggerAnalysis(0x0),
  fEventFilter(0x0),
  fCutsMother(0x0),   
  fEventStat(0x0),
  fEvent(0x0),
  fdEdXvsPt(0x0),
  fdEdXnSigmaElecvsPt(0x0),
  fdEdXvsPtTOF(0x0),
  fdEdXnSigmaElecvsPtTOF(0x0),
  fTOFbetavsPt(0x0),
  fTOFnSigmaElecvsPt(0x0),
  fNCrossedRowsTPC(0x0),
  fChi2ClusTPC(0x0),
  fRatioCrossClusTPC(0x0),
  fVem(0x0),
  fVep(0x0),
  fVemtmp(0x0),
  fVeptmp(0x0),
  fdconvphiv(acos(-1.0)),
  fdconvMee(100),
  fdop(0),
  fbz(0),
  fdv0mixing(kTRUE),
  fBGRejUnlike(kFALSE),
  fBGRejLike(kTRUE)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectronTG::AliAnalysisTaskMultiDielectronTG(const char *name) :
  AliAnalysisTaskSE(name),
  fListDielectron(),
  fListHistos(),
  fListCF(),
  fQAElectron(),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kMB),
  fExcludeTriggerMask(0),
  fTriggerOnV0AND(kFALSE),
  fRejectPileup(kFALSE),
  fTriggerLogic(kAny),
  fTriggerAnalysis(0x0),
  fEventFilter(0x0),
  fCutsMother(0x0),   
  fEventStat(0x0),
  fEvent(0x0),
  fdEdXvsPt(0x0),
  fdEdXnSigmaElecvsPt(0x0),
  fdEdXvsPtTOF(0x0),
  fdEdXnSigmaElecvsPtTOF(0x0),
  fTOFbetavsPt(0x0),
  fTOFnSigmaElecvsPt(0x0),
  fNCrossedRowsTPC(0x0),
  fChi2ClusTPC(0x0),
  fRatioCrossClusTPC(0x0),
  fVem(0x0),
  fVep(0x0),
  fVemtmp(0x0),
  fVeptmp(0x0),
  fdconvphiv(acos(-1.0)),
  fdconvMee(100),
  fdop(0),
  fbz(0),
  fdv0mixing(kTRUE),
  fBGRejUnlike(kFALSE),
  fBGRejLike(kTRUE)
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TH1D::Class());
  fListHistos.SetName("Dielectron_Histos_Multi");
  fListCF.SetName("Dielectron_CF_Multi");
  fListDielectron.SetOwner();
  fListHistos.SetOwner();
  fListCF.SetOwner();

  ///////////////
  for(int i=0;i<fgkNDIE; i++){
    for(int j=0;j<fgkNZBIN;j++){
      for(int k=0;k<fgkNCENT;k++){
	for(int l=0; l<fgkNRPBIN; l++){
	  fibuf[i][j][k][l] = 0;
	  fpoolp[i][j][k][l].clear();
	  fpoolm[i][j][k][l].clear();
	  for(int ib=0;ib<fgkNBUF; ib++){
	    fvep[ib][i][j][k][l].clear();
	    fvem[ib][i][j][k][l].clear();
	  }
	}
      }
    }
  }



}

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectronTG::~AliAnalysisTaskMultiDielectronTG()
{
  //
  // Destructor
  //

  //histograms and CF are owned by the dielectron framework.
  //however they are streamed to file, so in the first place the
  //lists need to be owner...
  fListHistos.SetOwner(kFALSE);
  fListCF.SetOwner(kFALSE);
  
  for(int i=0;i<fgkNDIE; i++){
    for(int j=0;j<fgkNZBIN;j++){
      for(int k=0;k<fgkNCENT;k++){
	for(int l=0; l<fgkNRPBIN; l++){
	  fibuf[i][j][k][l] = 0;
	  fpoolp[i][j][k][l].clear();
	  fpoolm[i][j][k][l].clear();
	  for(int ib=0;ib<fgkNBUF; ib++){
	    fvep[ib][i][j][k][l].clear();
	    fvem[ib][i][j][k][l].clear();
	  }
	}
      }
    }
  }
}
//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //

  if (!fListHistos.IsEmpty()||!fListCF.IsEmpty()) return; //already initialised

//   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
//   Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
//   Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    die->Init();
    if (die->GetHistogramList()) fListHistos.Add(const_cast<THashList*>(die->GetHistogramList()));
    if (die->GetCFManagerPair()) fListCF.Add(const_cast<AliCFContainer*>(die->GetCFManagerPair()->GetContainer()));
  }

  Int_t cuts=fListDielectron.GetEntries();
  Int_t nbins=kNbinsEvent+2*cuts;
  if (!fEventStat){
    fEventStat=new TH1D("hEventStat","Event statistics",nbins,0,nbins);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");

    //default names
    fEventStat->GetXaxis()->SetBinLabel(3,"Bin3 not used");
    fEventStat->GetXaxis()->SetBinLabel(4,"Bin4 not used");
    fEventStat->GetXaxis()->SetBinLabel(5,"Bin5 not used");
    
    if(fTriggerOnV0AND) fEventStat->GetXaxis()->SetBinLabel(3,"V0and triggers");
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
  

  
  int nbinx=400;
  float maxx=20;
  float minx=0.2;
  float binw = (TMath::Log(maxx)-TMath::Log(minx))/nbinx;
  double xbin[401];
  for(int ii=0;ii<nbinx+1;ii++){
    xbin[ii] = TMath::Exp(TMath::Log(minx) + 0.5*binw+binw*ii);
  }

  
  fQAElectron = new TList();
  fQAElectron->SetName("QAElectron");
  fQAElectron->SetOwner();


  fEvent = new TH1D("Event","centrality",   100,0,100);
  fQAElectron->Add(fEvent);
  fdEdXvsPt = new TH2D("dEdXvsPt","dE/dX vs. PT of TPC", nbinx, xbin, 2000,0,200);
  fQAElectron->Add(fdEdXvsPt);
  fdEdXnSigmaElecvsPt = new TH2D("fdEdXnSigmaElecvsPt"," dE/dX normalized to electron vs. pT of TPC",
                                 nbinx, xbin, 2000, -10, 10);
  fQAElectron->Add(fdEdXnSigmaElecvsPt);

  fdEdXvsPtTOF = new TH2D("dEdXvsPtTOF","dE/dX vs. PT of TPC", nbinx, xbin, 2000,0,200);
  fQAElectron->Add(fdEdXvsPtTOF);
  fdEdXnSigmaElecvsPtTOF = new TH2D("fdEdXnSigmaElecvsPtTOF"," dE/dX normalized to electron vs. pT of TPC",
                                 nbinx, xbin, 2000, -10, 10);
  fQAElectron->Add(fdEdXnSigmaElecvsPtTOF);



  fTOFbetavsPt = new TH2D("fTOFbetavsPt","TOF beta vs. p", 400, 0, 20, 1200, 0, 1.2);
  fQAElectron->Add(fTOFbetavsPt);
  fTOFnSigmaElecvsPt = new TH2D("fTOFnSigmaElecvsPt","TOF nsigma for electron", 400, 0, 20, 2000, -10, 10);
  fQAElectron->Add(fTOFnSigmaElecvsPt);

  fNCrossedRowsTPC = new TH2F("fNCrossedRowsTPC", "TPC nCrossed Rows vs. pT", 200, 0, 20, 200, 0, 200);
  fQAElectron->Add(fNCrossedRowsTPC);
  fChi2ClusTPC = new TH2F("fChi2ClusTPC", "hChi2ClusTPC vs. pT", 200, 0, 20, 200, 0, 10);
  fQAElectron->Add(fChi2ClusTPC);
  
  fRatioCrossClusTPC = new TH2F("fRatioCrossClusTPC", "hRatioCrossClusTPC vs. pT", 200, 0, 20, 200, 0, 10);     
  fQAElectron->Add(fRatioCrossClusTPC);

  fListHistos.Add(fQAElectron);

  fListHistos.SetOwner();  
  
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
  PostData(3, fEventStat);

  fCutsMother = new AliESDtrackCuts;
  fCutsMother->SetDCAToVertex2D(kTRUE);
  fCutsMother->SetMaxDCAToVertexZ(3.0);
  fCutsMother->SetMaxDCAToVertexXY(1.0);
  fCutsMother->SetPtRange(  0.05 , 200.0);
  fCutsMother->SetEtaRange( -0.84 , 0.84 );
  fCutsMother->SetAcceptKinkDaughters(kFALSE);
  fCutsMother->SetRequireITSRefit(kTRUE);
  fCutsMother->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  fCutsMother->SetMinNClustersITS(3);
  fCutsMother->SetRequireTPCRefit(kTRUE);


  AliInfo(Form("PairCutType = %d %d %d %d %d",
	       fRejectPairFlag[0],
	       fRejectPairFlag[1],
	       fRejectPairFlag[2],
	       fRejectPairFlag[3],
	       fRejectPairFlag[4]));

}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //

  if (fListHistos.IsEmpty()&&fListCF.IsEmpty()) return;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if (!inputHandler) return;
  
//   AliPIDResponse *pidRes=inputHandler->GetPIDResponse();
  if ( inputHandler->GetPIDResponse() ){
    // for the 2.76 pass2 MC private train. Together with a sigma shift of -0.169
//    pidRes->GetTPCResponse().SetSigma(4.637e-3,2.41332105409873257e+04);
    AliDielectronVarManager::SetPIDResponse( inputHandler->GetPIDResponse() );
  } else {
    AliFatal("This task needs the PID response attached to the input event handler!");
  }
  
  // Was event selected ?
  ULong64_t isSelected = AliVEvent::kAny;
  Bool_t isRejected = kFALSE;
  if( fSelectPhysics && inputHandler){
    if((isESD && inputHandler->GetEventSelection()) || isAOD){
      isSelected = inputHandler->IsEventSelected();
      if (fExcludeTriggerMask && (isSelected&fExcludeTriggerMask)) isRejected=kTRUE;
      if (fTriggerLogic==kAny) isSelected&=fTriggerMask;
      else if (fTriggerLogic==kExact) isSelected=((isSelected&fTriggerMask)==fTriggerMask);
    }
   }
 
 
  //Before physics selection
  fEventStat->Fill(kAllEvents);
  if (isSelected==0||isRejected) {
    PostData(3,fEventStat);
    return;
  }
  //after physics selection
  fEventStat->Fill(kSelectedEvents);

  //V0and
  if(fTriggerOnV0AND){
    if(isESD){if (!fTriggerAnalysis->IsOfflineTriggerFired(static_cast<AliESDEvent*>(InputEvent()), AliTriggerAnalysis::kV0AND))
            return;}
    if(isAOD){if(!((static_cast<AliAODEvent*>(InputEvent()))->GetVZEROData()->GetV0ADecision() == AliVVZERO::kV0BB &&
		   (static_cast<AliAODEvent*>(InputEvent()))->GetVZEROData()->GetV0CDecision() == AliVVZERO::kV0BB) )
	return;}
  }
  

  fEventStat->Fill(kV0andEvents);

  //Fill Event histograms before the event filter
  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  Double_t values[AliDielectronVarManager::kNMaxValues]={0};
  Double_t valuesMC[AliDielectronVarManager::kNMaxValues]={0};
  AliDielectronVarManager::SetEvent(InputEvent());
  AliDielectronVarManager::Fill(InputEvent(),values);
  AliDielectronVarManager::Fill(InputEvent(),valuesMC);
  Bool_t hasMC=AliDielectronMC::Instance()->HasMC();
  if (hasMC) {
    if (AliDielectronMC::Instance()->ConnectMCEvent())
      AliDielectronVarManager::Fill(AliDielectronMC::Instance()->GetMCEvent(),valuesMC);
  }

  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    AliDielectronHistos *h=die->GetHistoManager();
    if (h){
      if (h->GetHistogramList()->FindObject("Event_noCuts"))
        h->FillClass("Event_noCuts",AliDielectronVarManager::kNMaxValues,values);
      if (hasMC && h->GetHistogramList()->FindObject("MCEvent_noCuts"))
        h->FillClass("Event_noCuts",AliDielectronVarManager::kNMaxValues,valuesMC);
    }
  }
  nextDie.Reset();
  
  //event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(InputEvent())) return;
  }
  fEventStat->Fill(kFilteredEvents);
  
  //pileup
  if (fRejectPileup){
    if (InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.)) return;
  }
  fEventStat->Fill(kPileupEvents);
  
  //bz for AliKF
  fbz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( fbz );
  AliDielectronPID::SetCorrVal((Double_t)InputEvent()->GetRunNumber());
  
  //Process event in all AliDielectron instances
  //   TIter nextDie(&fListDielectron);
  //   AliDielectron *die=0;
  Int_t idie=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    //AliInfo(Form(" **************** die->Process(InputEvent()) : %d **************************", idie));
    die->SetDontClearArrays(kTRUE);
    die->Process(InputEvent());
    if (die->HasCandidates()){
      Int_t ncandidates=die->GetPairArray(1)->GetEntriesFast();
      if (ncandidates==1) fEventStat->Fill((kNbinsEvent)+2*idie);
      else if (ncandidates>1) fEventStat->Fill((kNbinsEvent+1)+2*idie);
    }


    AliDielectronVarManager::Fill(InputEvent(), fgValues);
    for (Int_t ii=0; ii<2; ++ii){ 
      TObjArray *obj = (TObjArray*)die->GetTrackArray(ii);
      Int_t ntracks=obj->GetEntriesFast();
      //AliInfo(Form(" ************** # of tracks = %d", ntracks));
      for (Int_t itrack=0; itrack<ntracks; ++itrack){
	
	////////////////////////////////////////////////////////////////////
	AliDielectronVarManager::Fill(obj->UncheckedAt(itrack), fgValues);
        ////////////////////////////////////////////////////////////////////

	if(fgValues[AliDielectronVarManager::kCharge]>0){
	  fVeptmp.push_back(new  AliDielectronSingleTG(1, 
						       fgValues[AliDielectronVarManager::kCentrality],
						       fgValues[AliDielectronVarManager::kXv],
						       fgValues[AliDielectronVarManager::kYv],
						       fgValues[AliDielectronVarManager::kZv],
						       fgValues[AliDielectronVarManager::kPx],
						       fgValues[AliDielectronVarManager::kPy],
						       fgValues[AliDielectronVarManager::kPz],
						       fgValues[AliDielectronVarManager::kPt],
						       fgValues[AliDielectronVarManager::kEta],
						       fgValues[AliDielectronVarManager::kPhi],
						       fgValues[AliDielectronVarManager::kTheta],
						       1, 1, static_cast<AliVTrack*>(obj->UncheckedAt(itrack)))
			    );
	}else if(fgValues[AliDielectronVarManager::kCharge]<0){
	  fVemtmp.push_back(new  AliDielectronSingleTG(-1, 
						       fgValues[AliDielectronVarManager::kCentrality],
						       fgValues[AliDielectronVarManager::kXv],
						       fgValues[AliDielectronVarManager::kYv],
						       fgValues[AliDielectronVarManager::kZv],
						       fgValues[AliDielectronVarManager::kPx],
						       fgValues[AliDielectronVarManager::kPy],
						       fgValues[AliDielectronVarManager::kPz],
						       fgValues[AliDielectronVarManager::kPt],
						       fgValues[AliDielectronVarManager::kEta],
						       fgValues[AliDielectronVarManager::kPhi],
						       fgValues[AliDielectronVarManager::kTheta],
						       1, 1, static_cast<AliVTrack*>(obj->UncheckedAt(itrack)))
			    );
	}
      }
    }
    //AliInfo(Form("size of e and p = %d %d", (int)fVeptmp.size(), (int)fVemtmp.size()));
    

    CheckGhostPairs(fVeptmp);
    CheckGhostPairs(fVemtmp);
    if(fRejectPairFlag[idie]==1 || fRejectPairFlag[idie]==2){
      RejectPairs(fVeptmp, fVemtmp, idie);
    }
    RandomizePool(fVeptmp, fVemtmp);    
    CalcPair(fVep, fVem, die, idie);

    //    AliInfo(Form("size of e and p (after) = %d %d", (int)fVep.size(), (int)fVem.size()));

    double dwcent = 100.0/fgkNCENT;
    double dwiz = 20.0/fgkNZBIN;
    double dwrp = acos(-1.0)/fgkNRPBIN;

    int icent = (int)(fgValues[AliDielectronVarManager::kCentrality]/dwcent);
    int izbin = (int)((fgValues[AliDielectronVarManager::kZvPrim]+10)/dwiz);
    int irp = (int)((fgValues[AliDielectronVarManager::kV0ACrpH2])/dwrp);
    
    if(icent<0) icent=0;
    if(icent>=fgkNCENT) icent=fgkNCENT-1;
    if(izbin<0) izbin=0;
    if(izbin>=fgkNZBIN) izbin=fgkNZBIN-1;
    if(irp<0) irp=0;
    if(irp>=fgkNRPBIN) irp=fgkNRPBIN-1;
    
    fvep[fibuf[idie][izbin][icent][irp]][idie][izbin][icent][irp].clear();
    for(int iep = 0; iep<(int)fVep.size();iep++) {
      fvep[fibuf[idie][izbin][icent][irp]][idie][izbin][icent][irp].push_back(fVep[iep]);
      fpoolp[idie][izbin][icent][irp].push_back(fVep[iep]);
      if(fpoolp[idie][izbin][icent][irp].size()>fgkMAXPOOL) {
	fpoolp[idie][izbin][icent][irp].pop_front();
      }
    }
    fvem[fibuf[idie][izbin][icent][irp]][idie][izbin][icent][irp].clear();
    for(int iem = 0; iem<(int)fVem.size();iem++) {
      fvem[fibuf[idie][izbin][icent][irp]][idie][izbin][icent][irp].push_back(fVem[iem]);
      fpoolm[idie][izbin][icent][irp].push_back(fVem[iem]);
      if(fpoolm[idie][izbin][icent][irp].size()>fgkMAXPOOL) {
	fpoolm[idie][izbin][icent][irp].pop_front();
      }
    }


    fibuf[idie][izbin][icent][irp]++;
    if(fibuf[idie][izbin][icent][irp]>= fgkNBUF) fibuf[idie][izbin][icent][irp]=0; 


    fVeptmp.clear();
    fVemtmp.clear();
    fVep.clear();
    fVem.clear();


    ++idie;
  }


  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fEvent->Fill(values[AliDielectronVarManager::kCentrality]);
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    if(!fCutsMother->AcceptTrack(track)) continue;
    fdEdXvsPt->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
    fdEdXnSigmaElecvsPt->Fill(track->GetTPCmomentum(), 
                              AliDielectronVarManager::GetPIDResponse()->NumberOfSigmasTPC(track,
                                                                                      AliPID::kElectron)
                              -AliDielectronPID::GetCorrVal());
    /// for beta caliculaton 
    Double_t l = track->GetIntegratedLength();  // cm
    Double_t t = track->GetTOFsignal();
    Double_t t0 = AliDielectronVarManager::GetPIDResponse()->GetTOFResponse().GetTimeZero(); // ps
    Double_t beta = 0;
    if( (l < 360. || l > 800.) || (t <= 0.) || (t0 >999990.0) ) {
      beta=-9999;
    }
    else {
      t -= t0; // subtract the T0
      l *= 0.01;  // cm ->m
      t *= 1e-12; //ps -> s
    
      Double_t v = l / t;
      beta = v / TMath::C();
    }

    fTOFbetavsPt->Fill(track->GetTPCmomentum(), beta);
    fTOFnSigmaElecvsPt->Fill(track->GetTPCmomentum(), 
                             AliDielectronVarManager::GetPIDResponse()->NumberOfSigmasTOF(track,
                                                                                     AliPID::kElectron));
    ////rough PID is required 
    if( fabs(AliDielectronVarManager::GetPIDResponse()->NumberOfSigmasTOF(track, AliPID::kElectron))<3){
      fdEdXvsPtTOF->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
      fdEdXnSigmaElecvsPtTOF->Fill(track->GetTPCmomentum(), 
				   AliDielectronVarManager::GetPIDResponse()->NumberOfSigmasTPC(track,
												AliPID::kElectron)
				   -AliDielectronPID::GetCorrVal());

      
      if(track->GetTPCsignal()>70 && track->GetTPCsignal()<90){
	fNCrossedRowsTPC->Fill(track->GetTPCmomentum(), track->GetTPCCrossedRows());
	//fChi2ClusTPC->Fill(track->GetTPCmomentum(), track->GetTPCchi2()/track->GetTPCNcls());
	fChi2ClusTPC->Fill(track->GetTPCmomentum(), track->GetTPCchi2());
	fRatioCrossClusTPC->Fill(track->GetTPCmomentum(), track->GetTPCCrossedRows()/track->GetTPCNclsF());
      }
    }
  }
  
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
  PostData(3,fEventStat);
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::FinishTaskOutput()
{
  //
  // Write debug tree
  //
  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    die->SaveDebugTree();
    AliDielectronMixingHandler *mix=die->GetMixingHandler();
//    printf("\n\n\n===============\ncall mix in Terminate: %p (%p)\n=================\n\n",mix,die);
    if(!mix) continue;
    for (Int_t ipool=0; ipool<mix->GetNumberOfBins(); ++ipool){
      mix->MixRemaining(die, ipool);
    }
  }
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::CheckGhostPairs(vector<AliDielectronSingleTG*> e1)
{
  ////// Check whether the pairs are in ghost pairs 
  /////  ghost means that one of the pairs is highly associated but reconstructed as a true particle


  bool reject = false;
  if(e1.size()>1){
    for(int i1=0; i1<(int)e1.size(); i1++){
      reject = false;
      for(int i2=i1+1; i2<(int)e1.size(); i2++){
        if( fabs(e1[i1]->Phi() - e1[i2]->Phi())<0.01 && 
	    fabs(e1[i1]->Eta() - e1[i2]->Eta())<0.005
	    ){
          reject = true;
          e1[i2]->SetGstFlag(0);
        }
      }
      if(reject==true)e1[i1]->SetGstFlag(0);
    }
  }
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskMultiDielectronTG::CheckGhost(vector<AliDielectronSingleTG*> e1, vector<AliDielectronSingleTG*> e2)
{
  ////// To be sure whether there are no ghost pairs in h event mixing 

  if(e1.size()>0 && e2.size()>0){
    for(int i1=0; i1<(int)e1.size(); i1++){
      for(int i2=0; i2<(int)e2.size(); i2++){
        if( fabs(e1[i1]->Phi() - e2[i2]->Phi())<0.01 && 
	    fabs(e1[i1]->Eta() - e2[i2]->Eta())<0.005
	    ){
	  return true;
        }
      }
    }
  }
  return false;
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::RandomizePool(vector<AliDielectronSingleTG*> e1, vector<AliDielectronSingleTG*> e2)
{
  
  ///// Randomize the pool constent to cancel the filling scheme of the single tracks
  //// 


  int size1 = e1.size();
  int usedindex[1000];
  for(int i=0;i<1000;i++){
    usedindex[i] = -1;
  }
  for(int i=0;i<size1;i++){
    usedindex[i] = 0;
  }

  for(int i=0;i<size1;i++){
    int j = (int)(gRandom->Uniform(0,size1));
    while(usedindex[j]==1){
      j = (int)(gRandom->Uniform(0,size1));
    }
    if( (e1[j]->GetGstFlag()==1) &&
	(e1[j]->GetConvFlag()==1)
        ){
      fVep.push_back(e1[j]);
    }
    usedindex[j] = 1;
  }
  

  int size2 = e2.size();
  for(int i=0;i<1000;i++){
    usedindex[i] = -1;
  }
  for(int i=0;i<size2;i++){
    usedindex[i] = 0;
  }

  for(int i=0;i<size2;i++){
    int j = (int)(gRandom->Uniform(0,size2));
    while(usedindex[j]==1){
      j = (int)(gRandom->Uniform(0,size2));
    }
    if( (e2[j]->GetGstFlag()==1) &&
	(e2[j]->GetConvFlag()==1)
        ){
      fVem.push_back(e2[j]);
    }
    usedindex[j] = 1;
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::CalcPair(vector<AliDielectronSingleTG*> ve1, 
						vector<AliDielectronSingleTG*> ve2, AliDielectron *die, Int_t idie)
{

  //
  // main routine for the pair procesing 
  //


  for(int i1=0; i1<(int)ve1.size(); i1++){
    for(int i2=0; i2<(int)ve2.size(); i2++){
      FillPair(ve1[i1], ve2[i2], 0, die, idie);  
    }    
  }

  for(int i1=0; i1<(int)ve1.size(); i1++){
    for(int i2=i1+1; i2<(int)ve1.size(); i2++){
      FillPair(ve1[i1], ve1[i2], 1, die, idie);  
    }    
  }

  for(int i1=0; i1<(int)ve2.size(); i1++){
    for(int i2=i1+1; i2<(int)ve2.size(); i2++){
      FillPair(ve2[i1], ve2[i2], 2, die, idie );  
    }    
  }


  double dwcent = 100.0/fgkNCENT;
  double dwiz = 20.0/fgkNZBIN;
  double dwrp = acos(-1.0)/fgkNRPBIN;

  int icent = (int)(fgValues[AliDielectronVarManager::kCentrality]/dwcent);
  int izbin = (int)((fgValues[AliDielectronVarManager::kZvPrim]+10)/dwiz);
  int irp = (int)((fgValues[AliDielectronVarManager::kV0ACrpH2])/dwrp);

  if(icent<0) icent=0;
  if(icent>=fgkNCENT) icent=fgkNCENT-1;
  if(izbin<0) izbin=0;
  if(izbin>=fgkNZBIN) izbin=fgkNZBIN-1;
  if(irp<0) irp=0;
  if(irp>=fgkNRPBIN) irp=fgkNRPBIN-1;

  int nmixed;
  if(ve1.size()>0) {
    //
    // Now mixed event for +- pairs
    //
    nmixed = 0;
    for(int ibuf=0;(nmixed<fgkNMix);ibuf++) {
      int ntry = 0;
      while((fBGRejUnlike && CheckGhost(ve1, fvem[ibuf][idie][izbin][icent][irp])) &&  ntry<fgkMAXTRY) {
        ReshuffleBuffer(fvem[ibuf][idie][izbin][icent][irp],fpoolm[idie][izbin][icent][irp]);
        ntry++;
      }
      for(int i1=0; i1<(int)ve1.size(); i1++){
	for(int i2=0; i2<(int)fvem[ibuf][idie][izbin][icent][irp].size(); i2++){
          FillPair(ve1[i1],fvem[ibuf][idie][izbin][icent][irp][i2], 3, die, idie);
        }
      }
      ++nmixed;
    }//for(ibuf)
  }
  if(ve2.size()>0) {
    //
    // Now mixed event for -+ pairs
    //
    nmixed = 0;
    for(int ibuf=0;(nmixed<fgkNMix);ibuf++) {
      int ntry = 0;
      while((fBGRejUnlike && CheckGhost(ve2, fvep[ibuf][idie][izbin][icent][irp])) &&  ntry<fgkMAXTRY) {
        ReshuffleBuffer(fvep[ibuf][idie][izbin][icent][irp],fpoolp[idie][izbin][icent][irp]);
        ntry++;
      }
      for(int i1=0; i1<(int)ve2.size(); i1++){
        for(int i2=0; i2<(int)fvep[ibuf][idie][izbin][icent][irp].size(); i2++){
	  FillPair(fvep[ibuf][idie][izbin][icent][irp][i2],ve2[i1],4, die, idie);
        }
      }
      ++nmixed;
    }//for(ibuf)
  }

  if(ve1.size()>0) {
    //
    // Now mixed event for ++ pairs
    //
    nmixed = 0;
    for(int ibuf=0;(nmixed<fgkNMix);ibuf++) {
      int ntry = 0;
      while((fBGRejLike && CheckGhost(ve1, fvep[ibuf][idie][izbin][icent][irp])) &&  ntry<fgkMAXTRY) {
        ReshuffleBuffer(fvep[ibuf][idie][izbin][icent][irp],fpoolp[idie][izbin][icent][irp]);
        ntry++;
      }
      for(int i1=0; i1<(int)ve1.size(); i1++){
        for(int i2=0;i2<(int)fvep[ibuf][idie][izbin][icent][irp].size(); i2++){
          FillPair(ve1[i1],fvep[ibuf][idie][izbin][icent][irp][i2], 5, die, idie);
        }
      }
      ++nmixed;
    }//for(ibuf)
  }

  if(ve2.size()>0) {
    //
    // Now mixed event for +- pairs
    //
    nmixed = 0;
    for(int ibuf=0;(nmixed<fgkNMix);ibuf++) {
      int ntry = 0;
      while((fBGRejLike && CheckGhost(ve2, fvem[ibuf][idie][izbin][icent][irp])) &&  ntry<fgkMAXTRY) {
        ReshuffleBuffer(fvem[ibuf][idie][izbin][icent][irp],fpoolm[idie][izbin][icent][irp]);
        ntry++;
      }
      for(int i1=0; i1<(int)ve2.size(); i1++){
        for(int i2=0; i2<(int)fvem[ibuf][idie][izbin][icent][irp].size(); i2++){
          FillPair(ve2[i1],fvem[ibuf][idie][izbin][icent][irp][i2],6, die, idie);
        }
      }
      ++nmixed;
    }//for(ibuf)
  }

}


//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::FillPair(AliDielectronSingleTG *iep, 
						AliDielectronSingleTG *iem, int type, AliDielectron *die, Int_t idie)
{                      

  //
  // main routine for filling kinematics of pairs
  //


  
  double dmass, dphiv, dpxpair, dpypair, dpzpair;
  double dptpair, depair, dphipair, detapair, dcos, dpsi;

  CalcVars(iep, iem, dmass, dphiv, dpxpair, dpypair, dpzpair, 
            dptpair, depair, dphipair, detapair, dcos, dpsi);


  double dopeningangle =  -9999;
  double dv0mass =  -9999;
  double dv0pxpair = -9999;
  double dv0pypair = -9999;
  double dv0pzpair = -9999;
  double dv0ptpair = -9999;
  double dv0epair = -9999;
  double dv0xvpair = -9999;
  double dv0yvpair = -9999;
  double dv0zvpair = -9999;
  double dv0phipair = -9999;
  double dv0etapair = -9999;
  double dv0rpair = -9999;
  double dpsipair =  -9999;
  double dphivpair =  -9999;

  ////////////////////////////
  ///// calculate v0 ////////
  ///////////////////////////
  Bool_t v0OFF=kFALSE;
  /// for the moment, this doesn't work for MC
  if(fdv0mixing == kFALSE && (type==3 || type==4 || type==5 || type==6)){
    v0OFF = kTRUE;
  }
  if(die->GetHasMC()==kTRUE && (type==3 || type==4 || type==5 || type==6)){
    v0OFF = kTRUE;
  }
  if(type==0 || type==1 || type==2){
    v0OFF = kFALSE;
  }
  

  if(v0OFF==kFALSE){
    AliVTrack* trackob1= iep->GetTrack();    
    AliVTrack* trackob2= iem->GetTrack();    

    if(!trackob1 || !trackob2){
      return; 
    }

    AliDielectronPair candidate;
    candidate.SetTracks(trackob1, (int)(11*iep->Charge()), 
			trackob2, (int)(11*iem->Charge()));
    
    if(trackob1==trackob2){
      AliInfo("Same Objects!!");
      return; 
    }

    //const AliKFParticle &kfPairLeg1 = candidate.GetKFFirstDaughter();
    //const AliKFParticle &kfPairLeg2 = candidate.GetKFSecondDaughter();
    
    const AliKFParticle &kfPair = candidate.GetKFParticle();

    /*
    dv0mass = candidate.M();
    dv0pxpair = candidate.Px();
    dv0pypair = candidate.Py();
    dv0pzpair = candidate.Pz();
    dv0ptpair = candidate.Pt();
    dv0epair = candidate.E();
    dv0xvpair = candidate.Xv();
    dv0yvpair = candidate.Yv();
    dv0zvpair = candidate.Zv();
    dv0phipair = candidate.Phi();
    dv0etapair = candidate.Eta();
    */
    dv0mass = kfPair.GetMass();
    dv0pxpair = kfPair.GetPx();
    dv0pypair = kfPair.GetPy();
    dv0pzpair = kfPair.GetPz();
    dv0ptpair = kfPair.GetPt();
    dv0epair = kfPair.GetE();
    dv0xvpair = kfPair.GetX();
    dv0yvpair = kfPair.GetY();
    dv0zvpair = kfPair.GetZ();
    dv0phipair = kfPair.GetPhi();
    dv0etapair = kfPair.GetEta();
    dv0rpair =  kfPair.GetR();

    dopeningangle = candidate.OpeningAngle();    
    dpsipair = candidate.PsiPair(fbz);
    dphivpair = candidate.PhivPair(fbz);


  }

  Double_t values[AliDielectronVarManager::kNMaxValues];
  TString  className1;
  TString  className2;
  className1.Form("MyPair_%s",kPairClassNamesTG[type]);
  className2.Form("MyPairV0_%s",kPairClassNamesTG[type]);
  
  AliDielectronHistos *fHistos = die->GetHistoManager();
  Bool_t pairClass1=fHistos->GetHistogramList()->FindObject(className1.Data())!=0x0;
  Bool_t pairClass2=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;

  if (pairClass1 && PairTrackcut(dphiv, dcos, dmass, idie)==true){
    ///import pair variables to values!!
    values[AliDielectronVarManager::kPx] = dpxpair;
    values[AliDielectronVarManager::kPy] = dpypair;
    values[AliDielectronVarManager::kPz] = dpzpair;
    values[AliDielectronVarManager::kPt] = dptpair;
    values[AliDielectronVarManager::kXv] = dv0xvpair;
    values[AliDielectronVarManager::kYv] = dv0yvpair;
    values[AliDielectronVarManager::kZv] = dv0zvpair;
    values[AliDielectronVarManager::kR] = dv0rpair;
    values[AliDielectronVarManager::kE] = depair;
    values[AliDielectronVarManager::kEta] = detapair;
    values[AliDielectronVarManager::kM] = dmass;
    values[AliDielectronVarManager::kPsiPair] = dphiv;
    values[AliDielectronVarManager::kPhivPair] = dphiv;
    values[AliDielectronVarManager::kPhi]  = dphipair;
    values[AliDielectronVarManager::kOpeningAngle]  = dcos;
    values[AliDielectronVarManager::kCosPointingAngle]  = TMath::Abs(TMath::ATan2(TMath::Sin(iep->Phi()-iem->Phi()),TMath::Cos(iep->Phi()-iem->Phi())));
    fHistos->FillClass(className1, AliDielectronVarManager::kNMaxValues, values);
  }


  if (pairClass2 && PairTrackcut(dphiv, dopeningangle, dv0mass, idie)==true){
    values[AliDielectronVarManager::kPx] = dv0pxpair;
    values[AliDielectronVarManager::kPy] = dv0pypair;
    values[AliDielectronVarManager::kPz] = dv0pzpair;
    values[AliDielectronVarManager::kPt] = dv0ptpair;
    values[AliDielectronVarManager::kXv] = dv0xvpair;
    values[AliDielectronVarManager::kYv] = dv0yvpair;
    values[AliDielectronVarManager::kZv] = dv0zvpair;
    values[AliDielectronVarManager::kR] = dv0rpair;
    values[AliDielectronVarManager::kE] = dv0epair;
    values[AliDielectronVarManager::kEta] = dv0etapair;
    values[AliDielectronVarManager::kM] = dv0mass;
    values[AliDielectronVarManager::kPsiPair] = dpsipair;
    values[AliDielectronVarManager::kPhivPair] = dphivpair;
    values[AliDielectronVarManager::kPhi]  = dv0phipair;
    values[AliDielectronVarManager::kOpeningAngle]  = dopeningangle;
    values[AliDielectronVarManager::kCosPointingAngle]  = TMath::Abs(TMath::ATan2(TMath::Sin(iep->Phi()-iem->Phi()),TMath::Cos(iep->Phi()-iem->Phi())));
    fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
  }


  
}

//_________________________________________________________________________________
bool AliAnalysisTaskMultiDielectronTG::PairTrackcut(double phiv, double op, double mass, int idie)
{

  //
  // pair-by-pair cuts
  //
  

  bool pairOK = true;

  if(fRejectPairFlag[idie] == 1 || fRejectPairFlag[idie] == 2 ||
     fRejectPairFlag[idie] == 3 || fRejectPairFlag[idie] == 4 ){
    if(fRejectPairFlag[idie] == 2 || fRejectPairFlag[idie] == 4 ){
      if(fbz>0 && (phiv>fdconvphiv && mass < fdconvMee) ){
	pairOK = false;
      }else if(fbz<0 && phiv<acos(-1.0)-fdconvphiv  && mass < fdconvMee){
	pairOK = false;
      }
    }else if(fRejectPairFlag[idie] == 1 || fRejectPairFlag[idie] == 3){
      if(op<fdop){
	pairOK = false;
      }
    }
  }

  return pairOK;

}


//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::CalcVars(AliDielectronSingleTG *iep, AliDielectronSingleTG *iem, 
						double &mass, double &phiv, double &px, double &py, double&pz,
						double &pt, double &e, double &phi, 
						double &eta, double &cos, double &psi)
{


  //
  // standalone calculator for the pair variables
  //

  px = iep->Px()+iem->Px();
  py = iep->Py()+iem->Py();
  pz = iep->Pz()+iem->Pz();
  pt = sqrt(px*px+py*py);
  double dppair = sqrt(pt*pt+pz*pz);
  static const double me=0.0005109989;
  e = sqrt(me*me+iep->Px()*iep->Px()+iep->Py()*iep->Py()+iep->Pz()*iep->Pz())
    + sqrt(me*me+iem->Px()*iem->Px()+iem->Py()*iem->Py()+iem->Pz()*iem->Pz());
  
  mass =  e*e-px*px-py*py-pz*pz;
  if(mass>=0){
    mass = sqrt(mass);
  }
   
  
  phi = atan2(py, px);
  eta = -0.5*TMath::Log((dppair+pz)/(dppair-pz));
  double p1 = sqrt(pow(iep->Px(),2)+pow(iep->Py(),2)+pow(iep->Pz(),2));
  double p2 = sqrt(pow(iem->Px(),2)+pow(iem->Py(),2)+pow(iem->Pz(),2));
  cos = acos((iep->Px()*iem->Px()+iep->Py()*iem->Py()+iep->Pz()*iem->Pz())/(p1*p2));

  double dtheta = iep->Theta()-iem->Theta();
  psi = asin(dtheta/cos);


  //unit vector of (pep+pem) 
  double pl = dppair;
  double ux = px/pl;
  double uy = py/pl;
  double uz = pz/pl;
  double ax = uy/sqrt(ux*ux+uy*uy);
  double ay = -ux/sqrt(ux*ux+uy*uy); 
  
  //momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by 
  //definition. 
  //double ptep = iep->Px()*ax + iep->Py()*ay; 
  //double ptem = iem->Px()*ax + iem->Py()*ay; 
  
  double pxep = iep->Px();
  double pyep = iep->Py();
  double pzep = iep->Pz();
  double pxem = iem->Px();
  double pyem = iem->Py();
  double pzem = iem->Pz();
  
  
  //vector product of pep X pem 
  double vpx = pyep*pzem - pzep*pyem; 
  double vpy = pzep*pxem - pxep*pzem; 
  double vpz = pxep*pyem - pyep*pxem; 
  double vp = sqrt(vpx*vpx+vpy*vpy+vpz*vpz); 
  //double thev = acos(vpz/vp); 
  
  //unit vector of pep X pem 
  double vx = vpx/vp; 
  double vy = vpy/vp; 
  double vz = vpz/vp; 

  //The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz) 
  double wx = uy*vz - uz*vy; 
  double wy = uz*vx - ux*vz; 
  double wz = ux*vy - uy*vx; 
  double wl = sqrt(wx*wx+wy*wy+wz*wz); 
  // by construction, (wx,wy,wz) must be a unit vector. 
  if(fabs(wl - 1.0) > 0.00001) std::cout << "Calculation error in W vector"<<std::endl; 
  // measure angle between (wx,wy,wz) and (ax,ay,0). The angle between them 
  // should be small if the pair is conversion 
  //
  double cosPhiV = wx*ax + wy*ay; 
  phiv = acos(cosPhiV); 
  
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::ReshuffleBuffer(vector<AliDielectronSingleTG*> ve, deque<AliDielectronSingleTG*> pool)
{
  //If there is not enough electron in the pool, give up
  //
  // ReshuffleBuffer for th event mixing 
  //

  unsigned int ne = ve.size();
  unsigned int poolsize = pool.size();
  int used[fgkMAXPOOL];
  for(int i=0;i<(int)fgkMAXPOOL;i++){
    used[i]=0;
  }

  if(poolsize < ne) {
    std::cout <<" pool size="<<poolsize<<" ne"<<ne<<std::endl;
    return;
  }
  for(unsigned int ie=0; ie < ne; ie++) {
    int j = rand()%poolsize;
    while(used[j]==1){
      j = rand()%poolsize;    
    }
    ve[ie] = pool[j];
    used[j]=1;
  }

}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::RejectPairs(vector<AliDielectronSingleTG*> e1, vector<AliDielectronSingleTG*> e2, Int_t idie){

  ////////////////////////////////////
  ///// to reject pairs from track arrays 
  ///// conversions, ghost ..
  ///////////////////////////////////
  if(e1.size()>0 && e2.size()>0){
    for(int i1=0; i1<(int)e1.size(); i1++){
      for(int i2=0; i2<(int)e2.size(); i2++){
	if(fRejectPairFlag[idie]==1){
	  Double_t cosine = GetOpeningAngle(e1[i1], e2[i2]);
	  if(cosine<fdop){
	    e1[i1]->SetConvFlag(0);
	    e2[i2]->SetConvFlag(0);
	  }
	}else if(fRejectPairFlag[idie]==2){
	  Double_t phiv = GetPhiv(e1[i1], e2[i2]);
	  Double_t mee = GetMass(e1[i1], e2[i2]);
	  if(fbz>0 && ( phiv>fdconvphiv && mee < fdconvMee) ){
	    e1[i1]->SetConvFlag(0);
	    e2[i2]->SetConvFlag(0);
	  }else if(fbz<0 && phiv<acos(-1.0)-fdconvphiv && mee < fdconvMee){
	    e1[i1]->SetConvFlag(0);
	    e2[i2]->SetConvFlag(0);
	  }
	}
      }
    }
  }
  if(e1.size()>0){
    for(int i1=0; i1<(int)e1.size(); i1++){
      for(int i2=i1+1; i2<(int)e1.size(); i2++){
	if(fRejectPairFlag[idie]==1){
	  Double_t cosine = GetOpeningAngle(e1[i1], e1[i2]);
	  if(cosine<fdop){
	    e1[i1]->SetConvFlag(0);
	    e1[i2]->SetConvFlag(0);
	  }
	}else if(fRejectPairFlag[idie]==2){
	  Double_t phiv = GetPhiv(e1[i1], e1[i2]);
	  Double_t mee = GetMass(e1[i1], e1[i2]);
	  if(fbz>0 && phiv>fdconvphiv && mee < fdconvMee){
	    e1[i1]->SetConvFlag(0);
	    e1[i2]->SetConvFlag(0);
	  }else if(fbz<0 && phiv<acos(-1.0)-fdconvphiv && mee < fdconvMee){
	    e1[i1]->SetConvFlag(0);
	    e1[i2]->SetConvFlag(0);
	  }
	}
      }
    }
  }

  if(e2.size()>0){
    for(int i1=0; i1<(int)e2.size(); i1++){
      for(int i2=i1+1; i2<(int)e2.size(); i2++){
	if(fRejectPairFlag[idie]==1){
	  Double_t cosine = GetOpeningAngle(e2[i1], e2[i2]);
	  if(cosine<fdop){
	    e2[i1]->SetConvFlag(0);
	    e2[i2]->SetConvFlag(0);
	  }
	}else if(fRejectPairFlag[idie]==2){
	  Double_t phiv = GetPhiv(e2[i1], e2[i2]);
	  Double_t mee = GetMass(e2[i1], e2[i2]);
	  if(fbz>0 && phiv>fdconvphiv && mee < fdconvMee){
	    e2[i1]->SetConvFlag(0);
	    e2[i2]->SetConvFlag(0);
	  }else if(fbz<0 && phiv<acos(-1.0)-fdconvphiv && mee < fdconvMee){
	    e2[i1]->SetConvFlag(0);
	    e2[i2]->SetConvFlag(0);
	  }
	}
      }
    }
  }
}


//_________________________________________________________________________________
Double_t AliAnalysisTaskMultiDielectronTG::GetOpeningAngle(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2){

  //////////////////////
  //////// calculate pairs and get opening angle 
  //////////////////////

  double dmass, dphiv, dpxpair, dpypair, dpzpair;
  double dptpair, depair, dphipair, detapair, dcos, dpsi;

  CalcVars(e1, e2, dmass, dphiv, dpxpair, dpypair, dpzpair, 
            dptpair, depair, dphipair, detapair, dcos, dpsi);

  return dcos;
}

//_________________________________________________________________________________
Double_t AliAnalysisTaskMultiDielectronTG::GetPhiv(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2){

  //////////////////////
  //////// calculate pairs and get phiv
  //////////////////////

  double dmass, dphiv, dpxpair, dpypair, dpzpair;
  double dptpair, depair, dphipair, detapair, dcos, dpsi;

  CalcVars(e1, e2, dmass, dphiv, dpxpair, dpypair, dpzpair, 
            dptpair, depair, dphipair, detapair, dcos, dpsi);

  return dphiv;
}

//_________________________________________________________________________________
Double_t AliAnalysisTaskMultiDielectronTG::GetMass(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2){

  //////////////////////
  //////// calculate pairs and get mass
  //////////////////////

  double dmass, dphiv, dpxpair, dpypair, dpzpair;
  double dptpair, depair, dphipair, detapair, dcos, dpsi;
  
  CalcVars(e1, e2, dmass, dphiv, dpxpair, dpypair, dpzpair, 
            dptpair, depair, dphipair, detapair, dcos, dpsi);

  return dmass;
}
