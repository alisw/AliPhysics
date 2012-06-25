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

const char* fgkPairClassNamesTG[7] = {
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
  tQAElectron(),
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
  hNCrossedRowsTPC(0x0),
  hChi2ClusTPC(0x0),
  hRatioCrossClusTPC(0x0),
  vem(0x0),
  vep(0x0),
  vem_tmp(0x0),
  vep_tmp(0x0),
  d_conv_phiv(acos(-1.0)),
  bz(0),
  d_v0_mixing(kTRUE)
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
  tQAElectron(),
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
  hNCrossedRowsTPC(0x0),
  hChi2ClusTPC(0x0),
  hRatioCrossClusTPC(0x0),
  vem(0x0),
  vep(0x0),
  vem_tmp(0x0),
  vep_tmp(0x0),
  d_conv_phiv(acos(-1.0)),
  bz(0),
  d_v0_mixing(kTRUE)
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
  for(int i=0;i<NDIE; i++){
    for(int j=0;j<NZBIN;j++){
      for(int k=0;k<NCENT;k++){
	for(int l=0; l<NRPBIN; l++){
	  d_ibuf[i][j][k][l] = 0;
	  d_poolp[i][j][k][l].clear();
	  d_poolm[i][j][k][l].clear();
	  for(int ib=0;ib<NBUF; ib++){
	    d_vep[ib][i][j][k][l].clear();
	    d_vem[ib][i][j][k][l].clear();
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
  
  for(int i=0;i<NDIE; i++){
    for(int j=0;j<NZBIN;j++){
      for(int k=0;k<NCENT;k++){
	for(int l=0; l<NRPBIN; l++){
	  d_ibuf[i][j][k][l] = 0;
	  d_poolp[i][j][k][l].clear();
	  d_poolm[i][j][k][l].clear();
	  for(int ib=0;ib<NBUF; ib++){
	    d_vep[ib][i][j][k][l].clear();
	    d_vem[ib][i][j][k][l].clear();
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
  float max_x=20;
  float min_x=0.2;
  float binw = (TMath::Log(max_x)-TMath::Log(min_x))/nbinx;
  double xbin[401];
  for(int ii=0;ii<nbinx+1;ii++){
    xbin[ii] = TMath::Exp(TMath::Log(min_x) + 0.5*binw+binw*ii);
  }

  
  tQAElectron = new TList();
  tQAElectron->SetName("QAElectron");
  tQAElectron->SetOwner();


  fEvent = new TH1D("Event","centrality",   100,0,100);
  tQAElectron->Add(fEvent);
  fdEdXvsPt = new TH2D("dEdXvsPt","dE/dX vs. PT of TPC", nbinx, xbin, 2000,0,200);
  tQAElectron->Add(fdEdXvsPt);
  fdEdXnSigmaElecvsPt = new TH2D("fdEdXnSigmaElecvsPt"," dE/dX normalized to electron vs. pT of TPC",
                                 nbinx, xbin, 2000, -10, 10);
  tQAElectron->Add(fdEdXnSigmaElecvsPt);

  fdEdXvsPtTOF = new TH2D("dEdXvsPtTOF","dE/dX vs. PT of TPC", nbinx, xbin, 2000,0,200);
  tQAElectron->Add(fdEdXvsPtTOF);
  fdEdXnSigmaElecvsPtTOF = new TH2D("fdEdXnSigmaElecvsPtTOF"," dE/dX normalized to electron vs. pT of TPC",
                                 nbinx, xbin, 2000, -10, 10);
  tQAElectron->Add(fdEdXnSigmaElecvsPtTOF);



  fTOFbetavsPt = new TH2D("fTOFbetavsPt","TOF beta vs. p", 400, 0, 20, 1200, 0, 1.2);
  tQAElectron->Add(fTOFbetavsPt);
  fTOFnSigmaElecvsPt = new TH2D("fTOFnSigmaElecvsPt","TOF nsigma for electron", 400, 0, 20, 2000, -10, 10);
  tQAElectron->Add(fTOFnSigmaElecvsPt);

  hNCrossedRowsTPC = new TH2F("hNCrossedRowsTPC", "TPC nCrossed Rows vs. pT", 200, 0, 20, 200, 0, 200);
  tQAElectron->Add(hNCrossedRowsTPC);
  hChi2ClusTPC = new TH2F("hChi2ClusTPC", "hChi2ClusTPC vs. pT", 200, 0, 20, 200, 0, 10);
  tQAElectron->Add(hChi2ClusTPC);
  
  hRatioCrossClusTPC = new TH2F("hRatioCrossClusTPC", "hRatioCrossClusTPC vs. pT", 200, 0, 20, 200, 0, 10);     
  tQAElectron->Add(hRatioCrossClusTPC);

  fListHistos.Add(tQAElectron);

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
  bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );

  AliDielectronPID::SetCorrVal((Double_t)InputEvent()->GetRunNumber());
  
  //Process event in all AliDielectron instances
  //   TIter nextDie(&fListDielectron);
  //   AliDielectron *die=0;
  Int_t idie=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    //AliInfo(" **************** die->Process(InputEvent()) **************************");
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
	  vep_tmp.push_back(new  AliDielectronSingleTG(1, 
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
	  vem_tmp.push_back(new  AliDielectronSingleTG(-1, 
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
    //AliInfo(Form("size of e and p = %d %d", (int)vep_tmp.size(), (int)vem_tmp.size()));


    check_ghost_pairs(vep_tmp);
    check_ghost_pairs(vem_tmp);
    randomize_pool(vep_tmp, vem_tmp);    
    calc_pair(vep, vem, die, idie);

    //    AliInfo(Form("size of e and p (after) = %d %d", (int)vep.size(), (int)vem.size()));

    double dw_cent = 100.0/NCENT;
    double dw_iz = 20.0/NZBIN;
    double dw_rp = acos(-1.0)/NRPBIN;

    int icent = (int)(fgValues[AliDielectronVarManager::kCentrality]/dw_cent);
    int izbin = (int)((fgValues[AliDielectronVarManager::kZvPrim]+10)/dw_iz);
    int irp = (int)((fgValues[AliDielectronVarManager::kV0ACrpH2])/dw_rp);
    
    if(icent<0) icent=0;
    if(icent>=NCENT) icent=NCENT-1;
    if(izbin<0) izbin=0;
    if(izbin>=NZBIN) izbin=NZBIN-1;
    if(irp<0) irp=0;
    if(irp>=NRPBIN) irp=NRPBIN-1;
    
    d_vep[d_ibuf[idie][izbin][icent][irp]][idie][izbin][icent][irp].clear();
    for(int iep = 0; iep<(int)vep.size();iep++) {
      d_vep[d_ibuf[idie][izbin][icent][irp]][idie][izbin][icent][irp].push_back(vep[iep]);
      d_poolp[idie][izbin][icent][irp].push_back(vep[iep]);
      if(d_poolp[idie][izbin][icent][irp].size()>MAXPOOL) {
	d_poolp[idie][izbin][icent][irp].pop_front();
      }
    }
    d_vem[d_ibuf[idie][izbin][icent][irp]][idie][izbin][icent][irp].clear();
    for(int iem = 0; iem<(int)vem.size();iem++) {
      d_vem[d_ibuf[idie][izbin][icent][irp]][idie][izbin][icent][irp].push_back(vem[iem]);
      d_poolm[idie][izbin][icent][irp].push_back(vem[iem]);
      if(d_poolm[idie][izbin][icent][irp].size()>MAXPOOL) {
	d_poolm[idie][izbin][icent][irp].pop_front();
      }
    }


    d_ibuf[idie][izbin][icent][irp]++;
    if(d_ibuf[idie][izbin][icent][irp]>= NBUF) d_ibuf[idie][izbin][icent][irp]=0; 


    vep_tmp.clear();
    vem_tmp.clear();
    vep.clear();
    vem.clear();


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
	hNCrossedRowsTPC->Fill(track->GetTPCmomentum(), track->GetTPCCrossedRows());
	hChi2ClusTPC->Fill(track->GetTPCmomentum(), track->GetTPCchi2()/track->GetTPCNcls());
	hRatioCrossClusTPC->Fill(track->GetTPCmomentum(), track->GetTPCCrossedRows()/track->GetTPCNclsF());
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
    if (mix) mix->MixRemaining(die);
  }
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::check_ghost_pairs(vector<AliDielectronSingleTG*> e1){
  bool reject = false;
  if(e1.size()>1){
    for(int i1=0; i1<(int)e1.size(); i1++){
      reject = false;
      for(int i2=i1+1; i2<(int)e1.size(); i2++){
        if( fabs(e1[i1]->Phi() - e1[i2]->Phi())<0.01 ){
          reject = true;
          e1[i2]->SetGstFlag(0);
        }
      }
      if(reject==true)e1[i1]->SetGstFlag(0);
    }
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::randomize_pool(vector<AliDielectronSingleTG*> e1, vector<AliDielectronSingleTG*> e2){
  
  int size1 = e1.size();
  int used_index[1000];
  for(int i=0;i<1000;i++){
    used_index[i] = -1;
  }
  for(int i=0;i<size1;i++){
    used_index[i] = 0;
  }

  for(int i=0;i<size1;i++){
    int j = (int)(gRandom->Uniform(0,size1));
    while(used_index[j]==1){
      j = (int)(gRandom->Uniform(0,size1));
    }
    if( (e1[j]->GetGstFlag()==1) &&
	(e1[j]->GetConvFlag()==1)
        ){
      vep.push_back(e1[j]);
    }
    used_index[j] = 1;
  }
  

  int size2 = e2.size();
  for(int i=0;i<1000;i++){
    used_index[i] = -1;
  }
  for(int i=0;i<size2;i++){
    used_index[i] = 0;
  }

  for(int i=0;i<size2;i++){
    int j = (int)(gRandom->Uniform(0,size2));
    while(used_index[j]==1){
      j = (int)(gRandom->Uniform(0,size2));
    }
    if( (e2[j]->GetGstFlag()==1) &&
	(e2[j]->GetConvFlag()==1)
        ){
      vem.push_back(e2[j]);
    }
    used_index[j] = 1;
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::calc_pair(vector<AliDielectronSingleTG*> ve1, 
						vector<AliDielectronSingleTG*> ve2, AliDielectron *die, Int_t idie){

  for(int i1=0; i1<(int)ve1.size(); i1++){
    for(int i2=0; i2<(int)ve2.size(); i2++){
      fill_pair(ve1[i1], ve2[i2], 0, die);  
    }    
  }

  for(int i1=0; i1<(int)ve1.size(); i1++){
    for(int i2=i1+1; i2<(int)ve1.size(); i2++){
      fill_pair(ve1[i1], ve1[i2], 1, die);  
    }    
  }

  for(int i1=0; i1<(int)ve2.size(); i1++){
    for(int i2=i1+1; i2<(int)ve2.size(); i2++){
      fill_pair(ve2[i1], ve2[i2], 2, die);  
    }    
  }


  double dw_cent = 100.0/NCENT;
  double dw_iz = 20.0/NZBIN;
  double dw_rp = acos(-1.0)/NRPBIN;

  int icent = (int)(fgValues[AliDielectronVarManager::kCentrality]/dw_cent);
  int izbin = (int)((fgValues[AliDielectronVarManager::kZvPrim]+10)/dw_iz);
  int irp = (int)((fgValues[AliDielectronVarManager::kV0ACrpH2])/dw_rp);

  if(icent<0) icent=0;
  if(icent>=NCENT) icent=NCENT-1;
  if(izbin<0) izbin=0;
  if(izbin>=NZBIN) izbin=NZBIN-1;
  if(irp<0) irp=0;
  if(irp>=NRPBIN) irp=NRPBIN-1;

  int nmixed;
  if(ve1.size()>0) {
    //
    // Now mixed event for +- pairs
    //
    nmixed = 0;
    for(int ibuf=0;(nmixed<NMix);ibuf++) {
      int ntry = 0;
      while(ntry<MAX_TRY) {
        reshuffle_buffer(d_vem[ibuf][idie][izbin][icent][irp],d_poolm[idie][izbin][icent][irp]);
        ntry++;
      }
      for(int i1=0; i1<(int)ve1.size(); i1++){
	for(int i2=0; i2<(int)d_vem[ibuf][idie][izbin][icent][irp].size(); i2++){
          fill_pair(ve1[i1],d_vem[ibuf][idie][izbin][icent][irp][i2], 3, die);
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
    for(int ibuf=0;(nmixed<NMix);ibuf++) {
      int ntry = 0;
      while(ntry<MAX_TRY) {
        reshuffle_buffer(d_vep[ibuf][idie][izbin][icent][irp],d_poolp[idie][izbin][icent][irp]);
        ntry++;
      }
      for(int i1=0; i1<(int)ve2.size(); i1++){
        for(int i2=0; i2<(int)d_vep[ibuf][idie][izbin][icent][irp].size(); i2++){
          fill_pair(ve2[i1],d_vep[ibuf][idie][izbin][icent][irp][i2],4, die);
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
    for(int ibuf=0;(nmixed<NMix);ibuf++) {
      int ntry = 0;
      while(ntry<MAX_TRY) {
        reshuffle_buffer(d_vep[ibuf][idie][izbin][icent][irp],d_poolp[idie][izbin][icent][irp]);
        ntry++;
      }
      for(int i1=0; i1<(int)ve1.size(); i1++){
        for(int i2=0;i2<(int)d_vep[ibuf][idie][izbin][icent][irp].size(); i2++){
          fill_pair(ve1[i1],d_vep[ibuf][idie][izbin][icent][irp][i2], 5, die);
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
    for(int ibuf=0;(nmixed<NMix);ibuf++) {
      int ntry = 0;
      while(ntry<MAX_TRY) {
        reshuffle_buffer(d_vem[ibuf][idie][izbin][icent][irp],d_poolm[idie][izbin][icent][irp]);
        ntry++;
      }
      for(int i1=0; i1<(int)ve2.size(); i1++){
        for(int i2=0; i2<(int)d_vem[ibuf][idie][izbin][icent][irp].size(); i2++){
          fill_pair(ve2[i1],d_vem[ibuf][idie][izbin][icent][irp][i2],6, die);
        }
      }
      ++nmixed;
    }//for(ibuf)
  }

}


//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::fill_pair(AliDielectronSingleTG *iep, 
						AliDielectronSingleTG *iem, int type, AliDielectron *die){                      
  
  double d_mass, d_phiv, d_pxpair, d_pypair, d_pzpair;
  double d_ptpair, d_epair, d_phipair, d_etapair, d_cos, d_psi;

  calc_vars(iep, iem, d_mass, d_phiv, d_pxpair, d_pypair, d_pzpair, 
            d_ptpair, d_epair, d_phipair, d_etapair, d_cos, d_psi);


  double d_openingangle =  -9999;
  double d_v0_mass =  -9999;
  double d_v0_pxpair = -9999;
  double d_v0_pypair = -9999;
  double d_v0_pzpair = -9999;
  double d_v0_ptpair = -9999;
  double d_v0_epair = -9999;
  double d_v0_xv_pair = -9999;
  double d_v0_yv_pair = -9999;
  double d_v0_zv_pair = -9999;
  double d_v0_phipair = -9999;
  double d_v0_etapair = -9999;
  double d_v0_r_pair = -9999;
  double d_psi_pair =  -9999;

  ////////////////////////////
  ///// calculate v0 ////////
  ///////////////////////////
  Bool_t V0OFF=kFALSE;
  /// for the moment, this doesn't work for MC
  if(d_v0_mixing == kFALSE && (type==3 || type==4 || type==5 || type==6)){
    V0OFF = kTRUE;
  }
  if(die->GetHasMC()==kTRUE && (type==3 || type==4 || type==5 || type==6)){
    V0OFF = kTRUE;
  }
  if(type==0 || type==1 || type==2){
    V0OFF = kFALSE;
  }
  

  if(V0OFF==kFALSE){
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
    const AliKFParticle &kfPair = candidate.GetKFParticle();
    d_openingangle = candidate.OpeningAngle();
    d_v0_mass = candidate.M();
    d_v0_pxpair = candidate.Px();
    d_v0_pypair = candidate.Py();
    d_v0_pzpair = candidate.Pz();
    d_v0_ptpair = candidate.Pt();
    d_v0_epair = candidate.E();
    d_v0_xv_pair = candidate.Xv();
    d_v0_yv_pair = candidate.Yv();
    d_v0_zv_pair = candidate.Zv();
    d_v0_phipair = candidate.Phi();
    // d_v0_theta_pair = candidate.Theta();
    d_v0_etapair = candidate.Eta();
    d_v0_r_pair =  kfPair.GetR();
    
    d_psi_pair = candidate.PsiPair(bz);
  }

  Double_t values[AliDielectronVarManager::kNMaxValues];
  TString  className1;
  TString  className2;
  className1.Form("MyPair_%s",fgkPairClassNamesTG[type]);
  className2.Form("MyPairV0_%s",fgkPairClassNamesTG[type]);
  
  AliDielectronHistos *fHistos = die->GetHistoManager();
  Bool_t pairClass1=fHistos->GetHistogramList()->FindObject(className1.Data())!=0x0;
  Bool_t pairClass2=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;

  if (pairClass1 && PairTrackcut(d_phiv)==true){
    ///import pair variables to values!!
    values[AliDielectronVarManager::kPx] = d_pxpair;
    values[AliDielectronVarManager::kPy] = d_pypair;
    values[AliDielectronVarManager::kPz] = d_pzpair;
    values[AliDielectronVarManager::kPt] = d_ptpair;
    values[AliDielectronVarManager::kXv] = d_v0_xv_pair;
    values[AliDielectronVarManager::kYv] = d_v0_yv_pair;
    values[AliDielectronVarManager::kZv] = d_v0_zv_pair;
    values[AliDielectronVarManager::kR] = d_v0_r_pair;
    values[AliDielectronVarManager::kE] = d_epair;
    values[AliDielectronVarManager::kEta] = d_etapair;
    values[AliDielectronVarManager::kM] = d_mass;
    values[AliDielectronVarManager::kPsiPair] = d_phiv;
    values[AliDielectronVarManager::kPhi]  = d_phipair;
    values[AliDielectronVarManager::kOpeningAngle]  = d_cos;
    fHistos->FillClass(className1, AliDielectronVarManager::kNMaxValues, values);
  }


  if (pairClass2 && PairTrackcut(d_phiv)==true){
    values[AliDielectronVarManager::kPx] = d_v0_pxpair;
    values[AliDielectronVarManager::kPy] = d_v0_pypair;
    values[AliDielectronVarManager::kPz] = d_v0_pzpair;
    values[AliDielectronVarManager::kPt] = d_v0_ptpair;
    values[AliDielectronVarManager::kXv] = d_v0_xv_pair;
    values[AliDielectronVarManager::kYv] = d_v0_yv_pair;
    values[AliDielectronVarManager::kZv] = d_v0_zv_pair;
    values[AliDielectronVarManager::kR] = d_v0_r_pair;
    values[AliDielectronVarManager::kE] = d_v0_epair;
    values[AliDielectronVarManager::kEta] = d_v0_etapair;
    values[AliDielectronVarManager::kM] = d_v0_mass;
    values[AliDielectronVarManager::kPsiPair] = d_psi_pair;
    values[AliDielectronVarManager::kPhi]  = d_v0_phipair;
    values[AliDielectronVarManager::kOpeningAngle]  = d_openingangle;
    fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
  }


  
}

//_________________________________________________________________________________
bool AliAnalysisTaskMultiDielectronTG::PairTrackcut(double phiv){
  

  bool pairOK = true;

  //var is phiv for the moment 
  if(bz>0 && phiv>d_conv_phiv){
    pairOK = false;
  }else if(bz<0 && phiv<acos(-1.0)-d_conv_phiv){
    pairOK = false;
  }

  return pairOK;

}


//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectronTG::calc_vars(AliDielectronSingleTG *iep, AliDielectronSingleTG *iem, 
						double &mass, double &phiv, double &px, double &py, double&pz,
						double &pt, double &e, double &phi, 
						double &eta, double &cos, double &psi){

  px = iep->Px()+iem->Px();
  py = iep->Py()+iem->Py();
  pz = iep->Pz()+iem->Pz();
  pt = sqrt(px*px+py*py);
  double d_ppair = sqrt(pt*pt+pz*pz);
  static const double me=0.0005109989;
  e = sqrt(me*me+iep->Px()*iep->Px()+iep->Py()*iep->Py()+iep->Pz()*iep->Pz())
    + sqrt(me*me+iem->Px()*iem->Px()+iem->Py()*iem->Py()+iem->Pz()*iem->Pz());
  
  mass =  e*e-px*px-py*py-pz*pz;
  if(mass<0){
    mass = mass;
  }else{
    mass = sqrt(mass);
  }
   
  
  phi = atan2(py, px);
  eta = -0.5*TMath::Log((d_ppair+pz)/(d_ppair-pz));
  double p1 = sqrt(pow(iep->Px(),2)+pow(iep->Py(),2)+pow(iep->Pz(),2));
  double p2 = sqrt(pow(iem->Px(),2)+pow(iem->Py(),2)+pow(iem->Pz(),2));
  cos = acos((iep->Px()*iem->Px()+iep->Py()*iem->Py()+iep->Pz()*iem->Pz())/(p1*p2));

  double dtheta = iep->Theta()-iem->Theta();
  psi = asin(dtheta/cos);


  //unit vector of (pep+pem) 
  double pl = d_ppair;
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
void AliAnalysisTaskMultiDielectronTG::reshuffle_buffer(vector<AliDielectronSingleTG*> ve, deque<AliDielectronSingleTG*> pool){
  //If there is not enough electron in the pool, give up

  unsigned int ne = ve.size();
  unsigned int poolsize = pool.size();
  int used[MAXPOOL];
  for(int i=0;i<(int)MAXPOOL;i++){
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
