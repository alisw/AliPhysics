#include "AliAnalysisNuclMult.h"

// ROOT includes
#include <TMath.h>
#include "TChain.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TFile.h"
#include "TList.h"
#include "TH3F.h"

// AliRoot includes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisFilter.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliCentrality.h"
#include "AliPPVsMultUtils.h"

ClassImp(AliAnalysisNuclMult)

//_____________________________________________________________________________
AliAnalysisNuclMult::AliAnalysisNuclMult():
  AliAnalysisTaskSE(),
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fTrackFilter(0x0),
  fPPVsMultUtils(0),
  fPIDResponse(NULL),
  fList(new TList()),
  iMultEstimator(0),
  multiplicityMin(-999.),
  multiplicityMax(999.),
  htriggerMask(NULL),
  htriggerMask_noMB(NULL),
  hNspdVertex(NULL),
  hzvertex(NULL),
  hpileUp(NULL),
  hmult(NULL),
  hNtrack(NULL)
{
  for(Int_t i=0;i<9;i++) stdFlagPid[i]=0;
  fList->SetName("results");
  //fList->SetOwner();
}
//______________________________________________________________________________
AliAnalysisNuclMult::AliAnalysisNuclMult(const char *name):
  AliAnalysisTaskSE(name),
  fAOD(NULL), 
  fESD(NULL),
  fEvent(NULL),
  fTrackFilter(0x0),
  fPPVsMultUtils(0),
  fPIDResponse(NULL),
  fList(new TList()),
  iMultEstimator(0),
  multiplicityMin(-999.),
  multiplicityMax(999.),
  htriggerMask(NULL),
  htriggerMask_noMB(NULL),
  hNspdVertex(NULL),
  hzvertex(NULL),
  hpileUp(NULL),
  hmult(NULL),
  hNtrack(NULL) 
{
  for(Int_t i=0;i<9;i++) stdFlagPid[i]=0;
  DefineOutput(1, TList::Class());
  fList->SetName("results");
}
//_____________________________________________________________________________
AliAnalysisNuclMult::~AliAnalysisNuclMult()
{
  if(fList) delete fList;
}
//______________________________________________________________________________
void AliAnalysisNuclMult::UserCreateOutputObjects()
{

  //Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,mu,pi,K,p,d,t,He3,He4
  stdFlagPid[0]=1; stdFlagPid[1]=2; stdFlagPid[2]=4; stdFlagPid[3]=8; stdFlagPid[4]=16; stdFlagPid[5]=32; stdFlagPid[6]=64; stdFlagPid[7]=128; stdFlagPid[8]=256; 
  
  Char_t nameSpec[18][30];
  snprintf(nameSpec[0],20,"e^{+}"); snprintf(nameSpec[1],20,"#mu^{+}"); snprintf(nameSpec[2],20,"#pi^{+}"); snprintf(nameSpec[3],20,"K^{+}"); snprintf(nameSpec[4],20,"p"); snprintf(nameSpec[5],20,"d"); snprintf(nameSpec[6],20,"t"); snprintf(nameSpec[7],20,"^{3}He"); snprintf(nameSpec[8],20,"^{4}He");
  snprintf(nameSpec[9],20,"e^{-}"); snprintf(nameSpec[10],20,"#mu^{-}"); snprintf(nameSpec[11],20,"#pi^{-}"); snprintf(nameSpec[12],20,"K^{-}"); snprintf(nameSpec[13],20,"#bar{p}"); snprintf(nameSpec[14],20,"#bar{d}"); snprintf(nameSpec[15],20,"#bar{t}"); snprintf(nameSpec[16],20,"^{3}#bar{He}"); snprintf(nameSpec[17],20,"^{4}#bar{He}");

  htriggerMask = new TH1I("htriggerMask","Trigger mask. Attention: before to cut on multiplicity",34,0,34);
  const Char_t *xaxisTitle[34]={"kMB","kINT7","kMUON","kHighMult","kEMC1","kCINT5","kCMUS5","kMUSH7","kMUL7","kMUU7","kEMC7","kEMC8","kMUS7","kPHI1","kPHI7","kEMCEJE","kEMCEGA","kCentral","kSemiCentral","kDG5","kZED","kSPI7","kSPI","kINT8","kMuonSingleLowPt8","kMuonSingleHighPt8","kMuonLikeLowPt8","kMuonUnlikeLowPt8","kMuonUnlikeLowPt0","kTRD","kFastOnly","kUserDefined","kAny","kAnyINT"};
  for(Int_t i=0;i<34;i++) {
    htriggerMask->Fill(xaxisTitle[i],0);
  }
  
  htriggerMask_noMB = new TH1I("htriggerMask_noMB","Trigger mask excluding MB events. Attention: before to cut on multiplicity",34,0,34);
  for(Int_t i=0;i<34;i++) {
    htriggerMask_noMB->Fill(xaxisTitle[i],0);
  }
    
  hNspdVertex = new TH1I("hNspdVertex","Number of vertices determined in the SPD. Attention: before to cut on multiplicity;N_{vtx}^{SPD}",220,-20,200);

  hzvertex = new TH1F("hzvertex","z-vertex distribution. Attention: before to cut on multiplicity;z_{vtx} (cm)",1000,-50,50);

  hpileUp = new TH1I("hpileUp","If the event is tagged as pileup in the SPD. Attention: before to cut on multiplicity",2,0,2);
  const Char_t *xaxisTitle2[2]={"kFALSE","kTRUE"};
  for(Int_t i=0;i<2;i++) {
    hpileUp->Fill(xaxisTitle2[i],0);
  }
  
  hmult = new TH1I("hmult","Multiplicity distribution (after cuts on event)",30000,-100,200);
  if(iMultEstimator==0) hmult->GetXaxis()->SetTitle("V0M Multiplicity Percentile");
  else if(iMultEstimator==1) hmult->GetXaxis()->SetTitle("kTrackletsITSTPC");

  hNtrack = new TH1I("hNtrack","Number of tracks per event (after cuts on event);N_{tracks}",2100,-100,2000);

  fdEdxVSp[0] = new TH2F("fdEdxVSp_pos","TPC dE/dx (positive charge); p_{TPC}/|z| (GeV/c); dE/dx_{TPC} (a.u.)",200,0,10,1000,0,2000);//100,500//500,2000
  fdEdxVSp[1] = new TH2F("fdEdxVSp_neg","TPC dE/dx (negative charge); p_{TPC}/|z| (GeV/c); dE/dx_{TPC} (a.u.)",200,0,10,1000,0,2000);

  Char_t name_hDeDxExp[9][200];
  Char_t title_hDeDxExp[9][200];
  for(Int_t iS=0;iS<9;iS++) {
    snprintf(name_hDeDxExp[iS],200,"hDeDxExp_%s",nameSpec[iS]);
    snprintf(title_hDeDxExp[iS],200,"Expected TPC dE/dx of %s;p_{TPC}/|z| (GeV/c);dE/dx_{TPC} (a.u.)",nameSpec[iS]);
    hDeDxExp[iS] = new TProfile(name_hDeDxExp[iS],title_hDeDxExp[iS],200,0,10,0,1000,"");//,500,0,5,0,1000,"");
  }

  Char_t name_fNsigmaTPC[18][200];
  Char_t title_fNsigmaTPC[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fNsigmaTPC[iS],200,"NsigmaTPC_%s",nameSpec[iS]);
    snprintf(title_fNsigmaTPC[iS],200,"n#sigma_{TPC} (%s);p_{T}/|z| (GeV/c);n_{#sigma_{TPC}}^{%s}",nameSpec[iS],nameSpec[iS]);
    if(iS==5 || iS==5+9) fNsigmaTPC[iS] = new TH2F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],200,0,10,200,-10,10);
    else fNsigmaTPC[iS] = new TH2F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],1,0,10,1,-10,10);
  }
  
  Char_t name_fDca[2][18][200];
  Char_t title_fDca[2][18][200];
  for(Int_t i=0;i<2;i++) {
    for(Int_t iS=0;iS<18;iS++) {
      if(i==0){
	snprintf(name_fDca[i][iS],200,"fDca_%s",nameSpec[iS]);
	snprintf(title_fDca[i][iS],200,"DCA (%s) (before DCA_{z} cut);DCA_{xy} (cm);DCA_{z} (cm);p_{T}/|z| (GeV/c)",nameSpec[iS]);
      }
      else {
	snprintf(name_fDca[i][iS],200,"fDca_withTOF_%s",nameSpec[iS]);
	snprintf(title_fDca[i][iS],200,"DCA (%s) with TOF matching (before DCA_{z} cut);DCA_{xy} (cm);DCA_{z} (cm);p_{T}/|z| (GeV/c)",nameSpec[iS]);
      }
      if(iS==5 || iS==5+9) fDca[i][iS] = new TH3F(name_fDca[i][iS],title_fDca[i][iS],250,-1.0,1.0,68,-3.4,3.4,30,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
      else fDca[i][iS] = new TH3F(name_fDca[i][iS],title_fDca[i][iS],1,-1.0,1.0,1,-3.4,3.4,1,0,3);
    }
  }
  
  fBetaTOFvspt[0] = new TH2F("fBetaTOFvspt_pos","#beta_{TOF} (positive charge);p_{T}/|z| (GeV/c);#beta_{TOF}",200,0,10,260,0.4,1.05);//500,520//200,260
  fBetaTOFvspt[1] = new TH2F("fBetaTOFvspt_neg","#beta_{TOF} (negative charge);p_{T}/|z| (GeV/c);#beta_{TOF}",200,0,10,260,0.4,1.05);

  Char_t name_hBetaExp[9][200];
  Char_t title_hBetaExp[9][200];
  for(Int_t iS=0;iS<9;iS++) {
    snprintf(name_hBetaExp[iS],200,"hBetaTOFVSpt_Exp_%s",nameSpec[iS]);
    snprintf(title_hBetaExp[iS],200,"Expected #beta_{TOF} (%s);p_{T}/|z| (GeV/c); #beta_{TOF}",nameSpec[iS]);
    hBetaExp[iS] = new TProfile(name_hBetaExp[iS],title_hBetaExp[iS],200,0,10,0.4,1.05,"");
  }

  fM2tof[0] = new TH2F("fM2tof_pos","m^{2}_{TOF} (positive charge);p_{T}/|z| (GeV/c); m^{2}_{TOF} (GeV^{2}/c^{4})",120,0,6,300,0,6);
  fM2tof[1] = new TH2F("fM2tof_neg","m^{2}_{TOF} (negative charge);p_{T}/|z| (GeV/c); m^{2}_{TOF} (GeV^{2}/c^{4})",120,0,6,300,0,6);
  
  Char_t name_fM2vspt[18][200];
  Char_t title_fM2vspt[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fM2vspt[iS],200,"fM2vspt_%s",nameSpec[iS]);
    snprintf(title_fM2vspt[iS],200,"m^{2}_{TOF} (3#sigma TPC dE/dx cut on %s);p_{T}/|z| (GeV/c); m^{2}_{TOF} (GeV^{2}/c^{4})",nameSpec[iS]);
    if(iS==5 || iS==5+9) fM2vspt[iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],120,0,6,300,0,6);
    else fM2vspt[iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],1,0,6,1,0,6);
  }

  Char_t name_fpTcorr[200];
  snprintf(name_fpTcorr,200,"fpTcorr_%s",nameSpec[5]);
  fpTcorr[0] = new TF1(name_fpTcorr,"[0]-[1]*TMath::Exp(-[2]*x)",0,10);
  fpTcorr[0]->FixParameter(0,-3.39633e-03);
  fpTcorr[0]->FixParameter(1,4.38176e-01);
  fpTcorr[0]->FixParameter(2,3.04490e+00);
  fpTcorr[0]->GetXaxis()->SetTitle("p_{T}^{reco} (GeV/c)");
  fpTcorr[0]->GetYaxis()->SetTitle("p_{T}^{reco} - p_{T}^{true} (GeV/c)");

  snprintf(name_fpTcorr,200,"fpTcorr_%s",nameSpec[5+9]);
  fpTcorr[1] = new TF1(name_fpTcorr,"[0]-[1]*TMath::Exp(-[2]*x)",0,10);
  fpTcorr[1]->FixParameter(0,-2.73165e-03);
  fpTcorr[1]->FixParameter(1,4.66999e-01);
  fpTcorr[1]->FixParameter(2,3.09068e+00);
  fpTcorr[1]->GetXaxis()->SetTitle("p_{T}^{reco} (GeV/c)");
  fpTcorr[1]->GetYaxis()->SetTitle("p_{T}^{reco} - p_{T}^{true} (GeV/c)");
  
  fList->Add(htriggerMask);
  fList->Add(htriggerMask_noMB);
  fList->Add(hNspdVertex);
  fList->Add(hpileUp);
  fList->Add(hzvertex);
  fList->Add(hmult);
  fList->Add(hNtrack);
  for(Int_t i=0;i<2;i++) fList->Add(fdEdxVSp[i]);
  for(Int_t i=0;i<9;i++) fList->Add(hDeDxExp[i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fBetaTOFvspt[i]);
  for(Int_t i=0;i<9;i++) fList->Add(hBetaExp[i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[1][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2tof[i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[5+9*i]);

  // Post output data.
  PostData(1, fList);

  fPPVsMultUtils = new AliPPVsMultUtils();

}
//______________________________________________________________________________
void AliAnalysisNuclMult::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fAOD && !fESD){
    Printf("%s:%d AODEvent and ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    return;
  }
  if(fESD) fEvent = fESD;
  else fEvent = fAOD;
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler=(AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse();
  
  //----------------- Event selection --------------:
  
  //physics selection:
  Bool_t isPhysSelected = kFALSE;
  //Trigger mask filled:
  for(Int_t i=0;i<32;i++) {
    unsigned bit=(1<<i);
    if(inputHandler->IsEventSelected() & bit) htriggerMask->Fill(i);
  }
  if(inputHandler->IsEventSelected() & AliVEvent::kAny) htriggerMask->Fill(33-1);
  if(inputHandler->IsEventSelected() & AliVEvent::kAnyINT) htriggerMask->Fill(34-1);
  
  if(!(inputHandler->IsEventSelected() & AliVEvent::kMB)) {
    for(Int_t i=0;i<32;i++) {
      unsigned bit=(1<<i);
      if(inputHandler->IsEventSelected() & bit) htriggerMask_noMB->Fill(i);
    }
    if(inputHandler->IsEventSelected() & AliVEvent::kAny) htriggerMask_noMB->Fill(33-1);
    if(inputHandler->IsEventSelected() & AliVEvent::kAnyINT) htriggerMask_noMB->Fill(34-1);
  }
  
  if(multiplicityMin!=-999) isPhysSelected = ((inputHandler->IsEventSelected() & AliVEvent::kMB) || (inputHandler->IsEventSelected() & AliVEvent::kHighMult));
  else isPhysSelected = (inputHandler->IsEventSelected() & AliVEvent::kMB);
  if(!isPhysSelected) return;

  const AliVVertex* vtxEVENT = fEvent->GetPrimaryVertex();
  
  //event must have at least an SPD-determined primary vertex
  Int_t NevSpd=-2;
  NevSpd=vtxEVENT->GetNContributors();
  hNspdVertex->Fill(NevSpd);
  if(NevSpd<1) return;

  //event must not be tagged as pileup
  Bool_t isPileUpSpd=kFALSE;
  if(fESD) {
    isPileUpSpd=((AliESDEvent *)fEvent)->IsPileupFromSPD();
  }
  else if(fAOD) { 
    isPileUpSpd=((AliAODEvent *)fEvent)->IsPileupFromSPD();
  }
  hpileUp->Fill(isPileUpSpd);
  if(isPileUpSpd) return;
  
   //event must have a primary vertex located within |z| < 10 cm
  Float_t zvtx=999.9;
  zvtx = vtxEVENT->GetZ();
  hzvertex->Fill(zvtx);
  if(TMath::Abs(zvtx)>10) return;

  //--------------------Multiplicity determination (Mid-pseudo-rapidity estimator):
  Float_t mult=-99;
  if(iMultEstimator==0) {
    mult=fPPVsMultUtils->GetMultiplicityPercentile(fEvent,"V0M");
  }
  else if(iMultEstimator==1) {
    if(fESD) {
      mult=AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent*)fESD,AliESDtrackCuts::kTrackletsITSTPC,0.8);
    }
    else if(fAOD) {
      mult=((AliAODHeader * )fEvent->GetHeader())->GetRefMultiplicityComb08();
    }
  }
  
  if(multiplicityMin!=-999) {
    if(mult<multiplicityMin || mult>multiplicityMax) return;
  }
  
  hmult->Fill(mult);

  Int_t nTracks = fEvent->GetNumberOfTracks();
  hNtrack->Fill(nTracks);

  //----------------------loop on the TRACKS-----------------------------
  for(Int_t iT = 0; iT < nTracks; iT++) { 
    
    AliVTrack* track = (AliVTrack *) fEvent->GetTrack(iT);
    
    if (!track){
      continue;
    }
    
    //------------------------- Track cuts:

    //track selection:
    if(fESD) {
      if(!fTrackFilter->IsSelected(track)) continue;
    }
    else if(fAOD) {
      if(!(((AliAODTrack *) track)->TestFilterBit(16))) continue;
      //TestFilterBit(16) -- Standard Cuts with very loose DCA: GetStandardITSTPCTrackCuts2010(kFALSE) && SetMaxDCAToVertexXY(2.4) && SetMaxDCAToVertexZ(3.2) && SetDCaToVertex2D(kTRUE)
    }
    
     //For the geometrical cuts
    Double_t eta = track->Eta();
    Double_t rapidity = track->Y();
    
    if(TMath::Abs(eta)>0.8) continue;
    if(TMath::Abs(rapidity)>0.5) continue;
    
    //track extrapolation to the primary vertex
    Double_t b[2] = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
    if (!track->PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 100., b, bCov))
      continue;
    
    Double_t DCAxy = b[0];
    Double_t DCAz = b[1];
    
    //Cut on the DCAxy
    if(TMath::Abs(DCAxy)>1.0) continue;
    this->FillDca(DCAxy,DCAz,track);

    //Cut on the DCAz
    if(TMath::Abs(DCAz)>1.0) continue;
        
    //------------------------- Track info:

    Double_t phi= track->Phi();
    Short_t charge = track->Charge();
    Double_t p = track->P();
    Double_t pt = track->Pt();
    Double_t dedx = track->GetTPCsignal();
    Double_t pTPC = track->GetTPCmomentum();
    
    //------------------------- TPC info:
    
    if(charge>0) fdEdxVSp[0]->Fill(pTPC,dedx);
    else if(charge<0) fdEdxVSp[1]->Fill(pTPC,dedx);

    Double_t nsigmaTPC[9];
    Double_t expdedx[9];
    
    //Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,mu,pi,K,p,d,t,He3,He4
    Int_t FlagPid = 0;
    
    for(Int_t iS=0;iS<9;iS++){
      expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
      hDeDxExp[iS]->Fill(pTPC,expdedx[iS]);
      nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
      if(TMath::Abs(nsigmaTPC[iS])<3) {
	FlagPid += ((Int_t)TMath::Power(2,iS));
      }
    }
    
    this->PtCorrection(pt, FlagPid, charge);

    for(Int_t iS=0;iS<9;iS++){
      if(charge>0) fNsigmaTPC[iS]->Fill(pt,nsigmaTPC[iS]);
      else if(charge<0) fNsigmaTPC[iS+9]->Fill(pt,nsigmaTPC[iS]);
    }

    //-------------------------- TOF info:

    Bool_t kTOF=kFALSE;
    kTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    
    if(!kTOF) continue;

    Double_t tof  = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
    Double_t beta = 0.0;
    Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};

    Double_t exptimes[9];
    track->GetIntegratedTimes(exptimes);
    //Integrated times of nuclei:
    for(Int_t iN=5;iN<9;iN++) {
      exptimes[iN] = exptimes[4]*exptimes[4]*(massOverZ[iN]*massOverZ[iN]/p/p+1)/(massOverZ[4]*massOverZ[4]/p/p+1);
      exptimes[iN] = TMath::Sqrt(exptimes[iN]);
    }  
    
    for(Int_t iS=0;iS<9;iS++){
      if(exptimes[iS]<1e-12) continue;
      hBetaExp[iS]->Fill(pt,exptimes[0]/exptimes[iS]);
    }
    
    beta=exptimes[0];
    if(tof<1e-12) continue;
    beta=beta/tof;//beta = L/tof/c = t_e/tof
    
    if(charge>0) fBetaTOFvspt[0]->Fill(pt,beta);
    else if(charge<0) fBetaTOFvspt[1]->Fill(pt,beta);
    
    Double_t gamma2=1.-(beta*beta);
    if(gamma2<1e-12) continue;
    gamma2=1./gamma2;
    
    Double_t mass2=0.;
    mass2=beta*beta*gamma2;
    if(mass2<1e-12) continue;
    mass2=(p*p)/mass2;
    
    if(charge>0) fM2tof[0]->Fill(pt,mass2);
    else if(charge<0) fM2tof[1]->Fill(pt,mass2);
    
    for(Int_t iS=0;iS<9;iS++){
      if(FlagPid & stdFlagPid[iS]) {
	if(charge>0) {
	  fM2vspt[iS]->Fill(pt,mass2);
	}	  
	else if(charge<0) {
	  fM2vspt[iS+9]->Fill(pt,mass2);
	}
      }
    }
    
    
  }//end track loop
}//end loop on the events
//_____________________________________________________________________________
void AliAnalysisNuclMult::FillDca(Double_t DCAxy, Double_t DCAz, AliVTrack *track) {
  
  Double_t nsigmaTPC[9];

  //Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,mu,pi,K,p,d,t,He3,He4
  Int_t FlagPid = 0;
  
  Bool_t kTOF=kFALSE;
  kTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);

  Short_t charge = track->Charge();
  Double_t pt = track->Pt();

  for(Int_t iS=0;iS<9;iS++){
    nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
    if(TMath::Abs(nsigmaTPC[iS])<3) {
      FlagPid += ((Int_t)TMath::Power(2,iS));
    }
  }
  
  for(Int_t iS=0;iS<9;iS++){
    if(FlagPid & stdFlagPid[iS]) {
      if(charge>0) {
	fDca[0][iS]->Fill(DCAxy,DCAz,pt);
	if(kTOF) fDca[1][iS]->Fill(DCAxy,DCAz,pt);
      }      
      else if(charge<0) {
	fDca[0][iS+9]->Fill(DCAxy,DCAz,pt);
	if(kTOF) fDca[1][iS+9]->Fill(DCAxy,DCAz,pt);
      }    
    }
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisNuclMult::PtCorrection(Double_t &pt, Int_t FlagPid, Short_t charge) {

  if(FlagPid & stdFlagPid[5]) {
    if(charge>0) {
      pt=pt-fpTcorr[0]->Eval(pt);
    }
    else if(charge<0) {
      pt=pt-fpTcorr[1]->Eval(pt);
    }
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisNuclMult::Terminate(Option_t *)
{ 
  printf("Terminate()\n");
}
