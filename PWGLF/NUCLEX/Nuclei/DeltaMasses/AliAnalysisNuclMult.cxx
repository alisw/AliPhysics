#include "AliAnalysisNuclMult.h"

// ROOT includes
#include <TMath.h>
#include "TChain.h"
#include "TH2F.h"
#include "TH2I.h"
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
#include "TTree.h"

// AliRoot includes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
//#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliESDtrackCuts.h"
//#include "AliAODTrack.h"
//#include "AliAODPid.h"
#include "AliCentrality.h"
#include "AliPPVsMultUtils.h"
// for Monte Carlo:
//#include "AliMCEventHandler.h"
//#include "AliMCEvent.h"
//#include "AliStack.h"
//#include "AliVParticle.h"

ClassImp(AliAnalysisNuclMult)

//_____________________________________________________________________________
AliAnalysisNuclMult::AliAnalysisNuclMult():
  AliAnalysisTaskSE(),
  fESD(NULL),
  fEvent(NULL),
  fESDtrackCuts(NULL),
  fPPVsMultUtils(NULL),
  fPIDResponse(NULL),
  fList(new TList()),
  multMin(0.),
  multMax(100.),
  htriggerMask(NULL),
  htriggerMask_noMB(NULL),
  hzvertex(NULL),
  hNevent(NULL),
  hV0mult(NULL),
  hTrackMult(NULL),
  hCheckTrackSel(NULL)
{
  fList->SetName("results");
  //fList->SetOwner();
}
//______________________________________________________________________________
AliAnalysisNuclMult::AliAnalysisNuclMult(const char *name):
  AliAnalysisTaskSE(name),
  fESD(NULL),
  fEvent(NULL),
  fESDtrackCuts(NULL),
  fPPVsMultUtils(NULL),
  fPIDResponse(NULL),
  fList(new TList()),
  multMin(0.),
  multMax(100.),
  htriggerMask(NULL),
  htriggerMask_noMB(NULL),
  hzvertex(NULL),
  hNevent(NULL),
  hV0mult(NULL),
  hTrackMult(NULL),
  hCheckTrackSel(NULL)
{
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

  Char_t nameSpec[18][30];
  snprintf(nameSpec[0],20,"e^{+}"); snprintf(nameSpec[1],20,"#mu^{+}"); snprintf(nameSpec[2],20,"#pi^{+}"); snprintf(nameSpec[3],20,"K^{+}"); snprintf(nameSpec[4],20,"p"); snprintf(nameSpec[5],20,"d"); snprintf(nameSpec[6],20,"t"); snprintf(nameSpec[7],20,"^{3}He"); snprintf(nameSpec[8],20,"^{4}He");
  snprintf(nameSpec[9],20,"e^{-}"); snprintf(nameSpec[10],20,"#mu^{-}"); snprintf(nameSpec[11],20,"#pi^{-}"); snprintf(nameSpec[12],20,"K^{-}"); snprintf(nameSpec[13],20,"#bar{p}"); snprintf(nameSpec[14],20,"#bar{d}"); snprintf(nameSpec[15],20,"#bar{t}"); snprintf(nameSpec[16],20,"^{3}#bar{He}"); snprintf(nameSpec[17],20,"^{4}#bar{He}");
  
  htriggerMask = new TH1I("htriggerMask","Trigger mask. Attention: before to cut on multiplicity",32,0,32);
  const Char_t *xaxisTitle[32]={"kMB","kINT7","kMUON","kHighMult","kEMC1","kCINT5","kCMUS5","kMUSH7","kMUL7","kMUU7","kEMC7","kEMC8","kMUS7","kPHI1","kPHI7","kEMCEJE","kEMCEGA","kCentral","kSemiCentral","kDG5","kZED","kSPI7","kSPI","kINT8","kMuonSingleLowPt8","kMuonSingleHighPt8","kMuonLikeLowPt8","kMuonUnlikeLowPt8","kMuonUnlikeLowPt0","kTRD","kFastOnly","kUserDefined"};
  for(Int_t i=0;i<32;i++) {
    htriggerMask->Fill(xaxisTitle[i],0);
  }
  
  htriggerMask_noMB = new TH1I("htriggerMask_noMB","Trigger mask excluding MB events. Attention: before to cut on multiplicity",32,0,32);
  for(Int_t i=0;i<32;i++) {
    htriggerMask_noMB->Fill(xaxisTitle[i],0);
  }

  hzvertex = new TH1F("hzvertex","z-vertex distribution. After (only) the trigger selection and the selection of INEL. collisions. Therefore it's filled before the other event cuts;z_{vtx}",1000,-50,50);
 
  hNevent = new TH1I("hNevent","Event counter. To check the event selection compare the last bin with the number of hV0mult entries",7,0,7);
  const Char_t *xaxisTitle2[7]={"kMB || kHighMult","IsINELgtZERO","IsAcceptedVertexPosition","HasNoInconsistentSPDandTrackVertices","IsNotPileupSPDInMultBins","IsEventSelected","InsideMultiplicityBin"};
for(Int_t i=0;i<7;i++) {
    hNevent->Fill(xaxisTitle2[i],0);
  }
  
  hV0mult = new TH1I("hV0mult","Multiplicity distribution (after cuts on event);V0M Multiplicity Percentile",1000,0,100);

  hTrackMult = new TH1I("hTrackMult","Mid-pseudo-rapidity estimator of multiplicity;kTrackletsITSTPC (|#eta|<0.8)",1000,0,1000);
  
  hCheckTrackSel = new TH1I("hCheckTrackSel","Number of tracks per event after the track selection",12,0,12);
  const Char_t *xaxisTitle3[12]={"|y|<0.5","nTPCclusters>=70","chi2perTPCcluster<=4","isTPCrefit","isITSrefit","nSPD>0","NoKinkDaughters","chi2perITScluster<=36","isPropagatedToDca","|DCAxy|<1","|DCAz|<1","|#eta|<0.8"};
  for(Int_t i=0;i<12;i++) {
    hCheckTrackSel->Fill(xaxisTitle3[i],0);
  }

  hnTPCclusters[0] = new TH1I("hnTPCclusters_0","Number of TPC clusters (before track cuts)",200,0,200);
  hchi2TPC[0] = new TH1D("hchi2TPC_0","#chi^{2} per TPC cluster (before track cuts)",1000,0,100);
  hisTPCrefit[0] = new TH1I("hisTPCrefit_0","kTPCrefit (before track cuts)",2,0,2);
  hisITSrefit[0] = new TH1I("hisITSrefit_0","kITSrefit (before track cuts)",2,0,2);
  hnSPD[0] = new TH1I("hnSPD_0","Number of SPD rec. points (before track cuts)",3,0,3);
  hnKinkDaughters[0] = new TH1I("hnKinkDaughters_0","Number of Kink Daughters (before track cuts)",40,0,40);
  hsigmaToVtx[0] = new TH1D("hsigmaToVtx_0","Number of sigma to the vertex (before track cuts)",400,0,40);
  hchi2ITS[0] = new TH1D("hchi2ITS_0","#chi^{2} per ITS cluster (before track cuts)",1000,0,100);
  
  heta[0] = new TH1D("heta_0","#eta (before track cuts);#eta",200,-1,1);
  hisPropagatedToDca[0] = new TH1I("hisPropagatedToDca_0","kPropagatedToDca (before track cuts)",2,0,2);
  
  hnTPCclusters[1] = new TH1I("hnTPCclusters_1","Number of TPC clusters (after track cuts)",200,0,200);
  hchi2TPC[1] = new TH1D("hchi2TPC_1","#chi^{2} per TPC cluster (after track cuts)",1000,0,100);
  hisTPCrefit[1] = new TH1I("hisTPCrefit_1","kTPCrefit (after track cuts)",2,0,2);
  hisITSrefit[1] = new TH1I("hisITSrefit_1","kITSrefit (after track cuts)",2,0,2);
  hnSPD[1] = new TH1I("hnSPD_1","Number of SPD rec. points (after track cuts)",3,0,3);
  hnKinkDaughters[1] = new TH1I("hnKinkDaughters_1","Number of Kink Daughters (after track cuts)",40,0,40);
  hsigmaToVtx[1] = new TH1D("hsigmaToVtx_1","Number of sigma to the vertex (after track cuts)",400,0,40);
  hchi2ITS[1] = new TH1D("hchi2ITS_1","#chi^{2} per ITS cluster (after track cuts)",1000,0,100);
    
  heta[1] = new TH1D("heta_1","#eta (after track cuts);#eta",200,-1,1);
  hisPropagatedToDca[1] = new TH1I("hisPropagatedToDca_1","kPropagatedToDca (after track cuts)",2,0,2);
  
  fptCorr[0] = new TF1("fptCorr_d","[0]-[1]*TMath::Exp(-[2]*x)",0,10);
  fptCorr[0]->FixParameter(0,-3.97081e-03);
  fptCorr[0]->FixParameter(1,4.94023e-01);
  fptCorr[0]->FixParameter(2,3.21308e+00);
  fptCorr[0]->GetXaxis()->SetTitle("p_{T}^{reco} (GeV/c)");
  fptCorr[0]->GetYaxis()->SetTitle("p_{T}^{reco} - p_{T}^{true} (GeV/c)");

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
    snprintf(title_fNsigmaTPC[iS],200,"n#sigma_{TPC} (%s);p_{T} (GeV/c);n_{#sigma_{TPC}}^{%s}",nameSpec[iS],nameSpec[iS]);
    if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9) fNsigmaTPC[iS] = new TH2F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],200,0,10,1000,-50,50);
    else fNsigmaTPC[iS] = new TH2F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],1,0,10,1,-10,10);
  }

  Char_t name_fDca[18][200];
  Char_t title_fDca[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of %s;DCA_{xy} (cm);DCA_{z} (cm);p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fDca[iS]=new TH3F(name_fDca[iS],title_fDca[iS],200,-1.0,1.0,100,-1.0,1.0,100,0,5);
    else fDca[iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-2.0,2.0,1,0,5);
  }
    
  Char_t name_fNsigmaTOF[18][200];
  Char_t title_fNsigmaTOF[18][200];    
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fNsigmaTOF[iS],200,"NsigmaTOF_%s",nameSpec[iS]);
    snprintf(title_fNsigmaTOF[iS],200,"n#sigma_{TOF} (%s);p_{T} (GeV/c);n_{#sigma_{TOF}}^{%s}",nameSpec[iS],nameSpec[iS]);
    if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9) fNsigmaTOF[iS] = new TH2F(name_fNsigmaTOF[iS],title_fNsigmaTOF[iS],200,0,10,1000,-50,50);
    else fNsigmaTOF[iS] = new TH2F(name_fNsigmaTOF[iS],title_fNsigmaTOF[iS],1,0,10,1,-10,10);
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

  fM2tof[0] = new TH2F("fM2tof_pos","m^{2}_{TOF} (positive charge);p_{T} (GeV/c); m^{2}_{TOF} (GeV^{2}/c^{4})",200,0,10,500,0,10);
  fM2tof[1] = new TH2F("fM2tof_neg","m^{2}_{TOF} (negative charge);p_{T} (GeV/c); m^{2}_{TOF} (GeV^{2}/c^{4})",200,0,10,500,0,10);
  
  Char_t name_fM2vspt[18][200];
  Char_t title_fM2vspt[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fM2vspt[iS],200,"fM2vspt_%s",nameSpec[iS]);
    snprintf(title_fM2vspt[iS],200,"m^{2}_{TOF} (3#sigma TPC dE/dx cut on %s);p_{T} (GeV/c); m^{2}_{TOF} (GeV^{2}/c^{4})",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fM2vspt[iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],200,0,10,500,0,10);
    else fM2vspt[iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],1,0,10,1,0,10);
  }
    
  if(multMin<1e-18 && multMax>99) {//plots added only for integrated multiplicity
    fList->Add(htriggerMask);
    fList->Add(htriggerMask_noMB);
    fList->Add(hzvertex);
    fList->Add(hNevent);
  }

  fList->Add(hV0mult);
  fList->Add(hTrackMult);
  
  fList->Add(hCheckTrackSel);
  
  for(Int_t i=0;i<2;i++) fList->Add(hnTPCclusters[i]);
  for(Int_t i=0;i<2;i++) fList->Add(hchi2TPC[i]);
  for(Int_t i=0;i<2;i++) fList->Add(hisTPCrefit[i]);
  for(Int_t i=0;i<2;i++) fList->Add(hisITSrefit[i]);
  for(Int_t i=0;i<2;i++) fList->Add(hnSPD[i]);
  for(Int_t i=0;i<2;i++) fList->Add(hnKinkDaughters[i]);
  for(Int_t i=0;i<2;i++) fList->Add(hsigmaToVtx[i]);
  for(Int_t i=0;i<2;i++) fList->Add(hchi2ITS[i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(heta[i]);
  for(Int_t i=0;i<2;i++) fList->Add(hisPropagatedToDca[i]);

  //fList->Add(fptCorr[0]);

  //temp:
  for(Int_t i=0;i<2;i++) fList->Add(fdEdxVSp[i]);
  for(Int_t i=0;i<9;i++) fList->Add(hDeDxExp[i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[2+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[3+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fDca[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[5+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTOF[2+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTOF[3+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTOF[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTOF[5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fBetaTOFvspt[i]);
  for(Int_t i=0;i<9;i++) fList->Add(hBetaExp[i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2tof[i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[5+9*i]);
  
  // Post output data.
  PostData(1, fList);

}
//______________________________________________________________________________
void AliAnalysisNuclMult::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  //fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fESD){
    Printf("%s:%d AODEvent and ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    return;
  }
  fEvent = fESD;
      
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler=(AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse();
  
  //Trigger mask filled:
  for(Int_t i=0;i<32;i++) {
    unsigned bit=(1<<i);//shift of 1 of i positions
    if(inputHandler->IsEventSelected() & bit) htriggerMask->Fill(i);
  }
  //for no MB events:
  if(!(inputHandler->IsEventSelected() & AliVEvent::kMB)) {
    for(Int_t i=0;i<32;i++) {
      unsigned bit=(1<<i);
      if(inputHandler->IsEventSelected() & bit) htriggerMask_noMB->Fill(i);
    }
  }
  
  //------------------------- Event selection:
  this->EventSelectionMonitor();
  
  if(!fPPVsMultUtils->IsEventSelected(fEvent, AliVEvent::kMB) && !fPPVsMultUtils->IsEventSelected(fEvent, AliVEvent::kHighMult)) return;
  
  Float_t mult = fPPVsMultUtils->GetMultiplicityPercentile(fEvent, "V0M", kFALSE);//kFALSE because I made the event selection before
  
  if(!this->IsInsideMultiplicityBin(mult)) return;
  hV0mult->Fill(mult);
  //------------------------- Event selection (end)

  Int_t Ntracklets=fESDtrackCuts->GetReferenceMultiplicity((AliESDEvent*)fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);
  hTrackMult->Fill(Ntracklets);
  
  Int_t nTracks = fEvent->GetNumberOfTracks();
  
  //-------------------------------- Loop on reconstructed TRACKS
  for(Int_t iT=0;iT<nTracks;iT++) { 

    AliVTrack* track = (AliVTrack *) fEvent->GetTrack(iT);
    if (!track){
      continue;
    }
    
    Double_t rapidity = track->Y();
    
    //------------------------- Track cuts:
   
    //rapidity cut:
    if(TMath::Abs(rapidity)>0.5) continue;
    
    Double_t DCAxy, DCAz;
    if(!this->AcceptTrack(track, DCAxy, DCAz)) continue;
    
    //------------------------- Track cuts (end)
        
    Short_t charge = track->Charge();
    Double_t pt = track->Pt();

    Double_t nsigmaTPC[9];
    this->GetNsigmaTPC(track, nsigmaTPC);

    //pT correction applied:
    this->PtCorr(pt, nsigmaTPC);
    
    for(Int_t iS=0;iS<9;iS++){
      if(charge>0) fNsigmaTPC[iS]->Fill(pt,nsigmaTPC[iS]);
      else if(charge<0) fNsigmaTPC[iS+9]->Fill(pt,nsigmaTPC[iS]);
    }

    //DCA filled:
    for(Int_t iS=0;iS<9;iS++){
      if(TMath::Abs(nsigmaTPC[iS])<3) {
	if(charge>0) fDca[iS]->Fill(DCAxy,DCAz,pt);
	else if(charge<0) fDca[iS+9]->Fill(DCAxy,DCAz,pt);
      }
    }
  
    //TOF matching required:
    if(!this->IsTOFmatching(track)) continue;
    
    Double_t beta = this->GetBeta(track, pt);
    
    Double_t p = track->P();
    
    //mass determination:
    Double_t m2 = this->GetM2(p, beta);
    
    if(charge>0) fM2tof[0]->Fill(pt,m2);
    else if(charge<0) fM2tof[1]->Fill(pt,m2);
    
    for(Int_t iS=0;iS<9;iS++){
      if(TMath::Abs(nsigmaTPC[iS])<3) {
	if(charge>0) fM2vspt[iS]->Fill(pt,m2);
	else if(charge<0) fM2vspt[iS+9]->Fill(pt,m2);
      }
    }

  }//----------------------loop on reconstructed TRACKS (end)
  
}//end loop on the events
//_____________________________________________________________________________
void AliAnalysisNuclMult::Terminate(Option_t *) { 
  printf("Terminate()\n");
}
//_____________________________________________________________________________
void AliAnalysisNuclMult::EventSelectionMonitor() {
  
  if(fPPVsMultUtils->IsEventSelected(fEvent, AliVEvent::kMB) || fPPVsMultUtils->IsEventSelected(fEvent, AliVEvent::kHighMult)) {
    hNevent->Fill(5);
  }
  
  if(fPPVsMultUtils->IsSelectedTrigger(fEvent, AliVEvent::kMB) || fPPVsMultUtils->IsSelectedTrigger(fEvent, AliVEvent::kHighMult)) {
    hNevent->Fill(0);
    
    if(fPPVsMultUtils->IsINELgtZERO(fEvent)) {
      hNevent->Fill(1);
      
      hzvertex->Fill(fEvent->GetPrimaryVertex()->GetZ());
      if(fPPVsMultUtils->IsAcceptedVertexPosition(fEvent)) {
	hNevent->Fill(2);

	if(fPPVsMultUtils->HasNoInconsistentSPDandTrackVertices(fEvent)) {
	  hNevent->Fill(3);
	  	  
	  if(fPPVsMultUtils->IsNotPileupSPDInMultBins(fEvent)) {
	    hNevent->Fill(4);
	  }
	  
	  if(this->IsInsideMultiplicityBin(fPPVsMultUtils->GetMultiplicityPercentile(fEvent,"V0M"))) {
	    hNevent->Fill(6);
	  }
	  
	}
      }
    }
  }
  
  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisNuclMult::IsInsideMultiplicityBin(Float_t multiplicity) {

  if(multiplicity < multMin || multiplicity > multMax+1e-18) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisNuclMult::AcceptTrack(AliVTrack *track, Double_t &DCAxy, Double_t &DCAz) {
 
  //-----------------where the std. cuts act:
  Int_t nTPCclusters=((AliESDtrack *)track)->GetTPCclusters(0);
  Double_t chi2TPC=-1;
  if(nTPCclusters!=0) chi2TPC=track->GetTPCchi2()/Double_t(nTPCclusters);
  Bool_t isTPCrefit=(track->GetStatus() & AliVTrack::kTPCrefit);
  Bool_t isITSrefit=(track->GetStatus() & AliVTrack::kITSrefit);
  Int_t nSPD=0;
  for(int i=0;i<2;i++) {
    if(track->HasPointOnITSLayer(i)) nSPD++;
  }
  Int_t nKinkDaughters=track->GetKinkIndex(0);
  Double_t sigmaToVtx=fESDtrackCuts->GetSigmaToVertex((AliESDtrack *)track);//we don't cut on this variable
  Int_t nITSclusters=track->GetITSclusters(0);
  Double_t chi2ITS=-1;
  if(nITSclusters!=0) chi2ITS=track->GetITSchi2()/Double_t(nITSclusters);
  
  //----------------- other cuts:
  Double_t eta = track->Eta();
  Double_t dca[2], cov[3];
  Bool_t isPropagatedToDca=track->PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 100., dca, cov);
  DCAxy = dca[0], DCAz = dca[1];

  //before to cut:
  hnTPCclusters[0]->Fill(nTPCclusters);
  hchi2TPC[0]->Fill(chi2TPC);
  hisTPCrefit[0]->Fill(isTPCrefit);
  hisITSrefit[0]->Fill(isITSrefit);
  hnSPD[0]->Fill(nSPD);
  hnKinkDaughters[0]->Fill(nKinkDaughters);
  hsigmaToVtx[0]->Fill(sigmaToVtx);
  hchi2ITS[0]->Fill(chi2ITS);
  heta[0]->Fill(eta);
  hisPropagatedToDca[0]->Fill(isPropagatedToDca);

  this->TrackSelectionMonitor(nTPCclusters, chi2TPC, isTPCrefit, isITSrefit, nSPD, nKinkDaughters, chi2ITS, isPropagatedToDca, DCAxy, DCAz, eta);
  
  if(!fESDtrackCuts->AcceptTrack((AliESDtrack *)track)) return kFALSE;
  
  if(!isPropagatedToDca) return kFALSE;

  if(TMath::Abs(DCAxy)>1.0) return kFALSE;

  if(TMath::Abs(DCAz)>1.0) return kFALSE;

  if(TMath::Abs(eta)>0.8) return kFALSE;

  //after the cuts:
  hnTPCclusters[1]->Fill(nTPCclusters);
  hchi2TPC[1]->Fill(chi2TPC);
  hisTPCrefit[1]->Fill(isTPCrefit);
  hisITSrefit[1]->Fill(isITSrefit);
  hnSPD[1]->Fill(nSPD);
  hnKinkDaughters[1]->Fill(nKinkDaughters);
  hsigmaToVtx[1]->Fill(sigmaToVtx);
  hchi2ITS[1]->Fill(chi2ITS);
  heta[1]->Fill(eta);
  hisPropagatedToDca[1]->Fill(isPropagatedToDca);
  
  return kTRUE; 
}
//_____________________________________________________________________________
void AliAnalysisNuclMult::TrackSelectionMonitor(Int_t nTPCclusters, Double_t chi2TPC, Bool_t isTPCrefit, Bool_t isITSrefit, Int_t nSPD, Int_t nKinkDaughters, Double_t chi2ITS, Bool_t isPropagatedToDca, Double_t DCAxy, Double_t DCAz, Double_t eta) {
  
  hCheckTrackSel->Fill(0);
  if(nTPCclusters>=70) {
    hCheckTrackSel->Fill(1);
    
    if(chi2TPC<=4) {
      hCheckTrackSel->Fill(2);
      
      if(isTPCrefit) {
	hCheckTrackSel->Fill(3);
        
	if(isITSrefit) {
	  hCheckTrackSel->Fill(4);
	  
	  if(nSPD>0) {
	    hCheckTrackSel->Fill(5);
	    
	    if(nKinkDaughters<1) {
	      hCheckTrackSel->Fill(6);
	      
	      if(chi2ITS<=36) {
		hCheckTrackSel->Fill(7);
		
		if(isPropagatedToDca) {
		  hCheckTrackSel->Fill(8);
		  
		  if(TMath::Abs(DCAxy)<1.0) {
		    hCheckTrackSel->Fill(9);
		    
		    if(TMath::Abs(DCAz)<1.0) {//temporary set to 2.0
		      hCheckTrackSel->Fill(10);
		      
		      if(TMath::Abs(eta)<0.8) {
			hCheckTrackSel->Fill(11);
			
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisNuclMult::PtCorr(Double_t &pt, Double_t nsigmaTPC[9]) {

  if(TMath::Abs(nsigmaTPC[5])<3) pt=pt-fptCorr[0]->Eval(pt);
   
  return;
}
//_____________________________________________________________________________
void AliAnalysisNuclMult::GetNsigmaTPC(AliVTrack *track, Double_t nsigmaTPC[9]) {

  Short_t charge = track->Charge();
  Double_t pTPC = track->GetTPCmomentum();

  Double_t dedx = track->GetTPCsignal();

  Double_t expdedx[9];
  for(Int_t iS=0;iS<9;iS++){
    expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
  }
    
  if(charge>0) fdEdxVSp[0]->Fill(pTPC,dedx);
  else if(charge<0) fdEdxVSp[1]->Fill(pTPC,dedx);
  
  for(Int_t iS=0;iS<9;iS++){
    hDeDxExp[iS]->Fill(pTPC,expdedx[iS]);
  }
  
  for(Int_t iS=0;iS<9;iS++){
    nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
  }
  
  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisNuclMult::IsTOFmatching(AliVTrack *track) {

  Bool_t kTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
  
  if(!kTOF) return kFALSE;
  
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisNuclMult::GetExpTOFtime(AliVTrack *track, Double_t p, Double_t exptimes[9]) {

  Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};

  track->GetIntegratedTimes(exptimes);
  //Integrated times of nuclei:
  for(Int_t iN=5;iN<9;iN++) {
    if(p>1e-18) exptimes[iN] = exptimes[4]*exptimes[4]*(massOverZ[iN]*massOverZ[iN]/p/p+1)/(massOverZ[4]*massOverZ[4]/p/p+1);
    exptimes[iN] = TMath::Sqrt(exptimes[iN]);
  }  
  
  return;
}
//_____________________________________________________________________________
Double_t AliAnalysisNuclMult::GetBeta(AliVTrack *track, Double_t pt) {
  
  Short_t charge = track->Charge();
  Double_t p = track->P();
  Double_t tof = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
  Double_t len = track->GetIntegratedLength();

  Double_t exptimes[9];
  this->GetExpTOFtime(track, p, exptimes);
    
  Double_t tofres[9];
  for(Int_t iS=0;iS<9;iS++) tofres[iS] = fPIDResponse->GetTOFResponse().GetExpectedSigma(p, exptimes[iS], (AliPID::EParticleType) iS);

  Double_t nsigmaTOF[9];
  for(Int_t iS=0;iS<9;iS++) {
    nsigmaTOF[iS] = -99999.9;
    if(tofres[iS]>1e-18) nsigmaTOF[iS] = (tof-exptimes[iS])/tofres[iS];
  }

  for(Int_t iS=0;iS<9;iS++) {
    if(charge>0) fNsigmaTOF[iS]->Fill(pt,nsigmaTOF[iS]);
    else if(charge<0) fNsigmaTOF[iS+9]->Fill(pt,nsigmaTOF[iS]);
  }
  
  Double_t beta = -99.9;
  if(tof>1e-18) beta = len/(tof*2.99792457999999984e-02);//beta=L/(tc)
  
  if(charge>0) fBetaTOFvspt[0]->Fill(pt,beta);
  else if(charge<0) fBetaTOFvspt[1]->Fill(pt,beta);

  Double_t betaexp[9];
  for(Int_t iS=0;iS<9;iS++) {
    if(exptimes[iS]>1e-18) betaexp[iS] = len/(exptimes[iS]*2.99792457999999984e-02);
    else betaexp[iS] = -99.9;
  }

  for(Int_t iS=0;iS<9;iS++){
    hBetaExp[iS]->Fill(pt,betaexp[iS]);
  }
  
  return beta;
}
//_____________________________________________________________________________
Double_t AliAnalysisNuclMult::GetM2(Double_t p, Double_t beta) {

  if(beta<1e-18) return -1;
  
  Double_t m2 = (p*p)*(1-beta*beta)/(beta*beta);
  
  return m2;
}
