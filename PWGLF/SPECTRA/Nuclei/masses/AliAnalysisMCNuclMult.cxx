#include "AliAnalysisMCNuclMult.h"

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
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisFilter.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliCentrality.h"
#include "AliPPVsMultUtils.h"
// for Monte Carlo:
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliVParticle.h"

ClassImp(AliAnalysisMCNuclMult)

//_____________________________________________________________________________
AliAnalysisMCNuclMult::AliAnalysisMCNuclMult():
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
  fNtrVsMult(NULL),
  prNtrVsMult(NULL),
  hmult_tot(NULL),
  hmult(NULL),
  hNtrack(NULL),
  hpdg(NULL)
{
  for(Int_t i=0;i<9;i++) stdPdg[i] = 0;
  fList->SetName("results");
  //fList->SetOwner();
}
//______________________________________________________________________________
AliAnalysisMCNuclMult::AliAnalysisMCNuclMult(const char *name):
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
  fNtrVsMult(NULL),
  prNtrVsMult(NULL),
  hmult_tot(NULL),
  hmult(NULL),
  hNtrack(NULL),
  hpdg(NULL)
{
  for(Int_t i=0;i<9;i++) stdPdg[i] = 0;
  DefineOutput(1, TList::Class());
  fList->SetName("results");
}
//_____________________________________________________________________________
AliAnalysisMCNuclMult::~AliAnalysisMCNuclMult()
{
  if(fList) delete fList;
}
//______________________________________________________________________________
void AliAnalysisMCNuclMult::UserCreateOutputObjects()
{

  //e,mu,pi,K,p,d,t,He3,He4
  stdPdg[0] = 11; stdPdg[1] = 13; stdPdg[2] = 211; stdPdg[3] = 321; stdPdg[4] = 2212; stdPdg[5] = 10020; stdPdg[6] = 10030; stdPdg[7] = 20030; stdPdg[8] = 20040;
  
  Char_t nameSpec[18][30];
  snprintf(nameSpec[0],20,"e^{+}"); snprintf(nameSpec[1],20,"#mu^{+}"); snprintf(nameSpec[2],20,"#pi^{+}"); snprintf(nameSpec[3],20,"K^{+}"); snprintf(nameSpec[4],20,"p"); snprintf(nameSpec[5],20,"d"); snprintf(nameSpec[6],20,"t"); snprintf(nameSpec[7],20,"^{3}He"); snprintf(nameSpec[8],20,"^{4}He");
  snprintf(nameSpec[9],20,"e^{-}"); snprintf(nameSpec[10],20,"#mu^{-}"); snprintf(nameSpec[11],20,"#pi^{-}"); snprintf(nameSpec[12],20,"K^{-}"); snprintf(nameSpec[13],20,"#bar{p}"); snprintf(nameSpec[14],20,"#bar{d}"); snprintf(nameSpec[15],20,"#bar{t}"); snprintf(nameSpec[16],20,"^{3}#bar{He}"); snprintf(nameSpec[17],20,"^{4}#bar{He}");
  
  //-----Event characterization---------
  
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

  hzvertex = new TH1F("hzvertex","z-vertex distribution. Attention: before to cut on multiplicity;z_{vtx}",1000,-50,50);

  hpileUp = new TH1I("hpileUp","If the event is tagged as pileup in the SPD. Attention: before to cut on multiplicity",2,0,2);
  const Char_t *xaxisTitle2[2]={"kFALSE","kTRUE"};
  for(Int_t i=0;i<2;i++) {
    hpileUp->Fill(xaxisTitle2[i],0);
  }
  
  fNtrVsMult = new TH2I("fNtrVsMult","N_{tracks} vs. multiplicity (after cuts on event);multiplicity;N_{tracks}",1200,-200,1000,1000,-100,900);
  prNtrVsMult = new TProfile("prNtrVsMult","N_{tracks} vs. multiplicity (after cuts on event);multiplicity;N_{tracks}",1200,-200,1000,-100,900,"");

  hmult_tot = new TH1I("hmult_tot","Multiplicity distribution (after cuts on event)",30000,-100,200);
  if(iMultEstimator==0) hmult_tot->GetXaxis()->SetTitle("V0M Multiplicity Percentile");
  else if(iMultEstimator==1) hmult_tot->GetXaxis()->SetTitle("kTrackletsITSTPC");

  hmult = new TH1I("hmult","Multiplicity distribution (after cuts on event)",30000,-100,200);
  if(iMultEstimator==0) hmult->GetXaxis()->SetTitle("V0M Multiplicity Percentile");
  else if(iMultEstimator==1) hmult->GetXaxis()->SetTitle("kTrackletsITSTPC");

  hNtrack = new TH1I("hNtrack","Number of tracks per event (after cuts on event);N_{tracks}",2100,-100,2000);

  hpdg = new TH1I("hpdg","Pdg label of generated particles;Pdg label",50082,-25041,25041);

  Char_t name_hpt[200];
  Char_t title_hpt[200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_gen_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of generated %s (|y|<0.5). Attention: before to cut on multiplicity;p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[0][iS] = new TH1F(name_hpt,title_hpt,100,0,5);
    else hpt[0][iS] = new TH1F(name_hpt,title_hpt,1,0,5);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_trigger_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of %s after the trigger selection. Attention: before to cut on multiplicity;p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[1][iS] = new TH1F(name_hpt,title_hpt,100,0,5);
    else hpt[1][iS] = new TH1F(name_hpt,title_hpt,1,0,5);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_vtx_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of %s for events within |z_{vtx}|<10cm. Attention: before to cut on multiplicity;p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[2][iS] = new TH1F(name_hpt,title_hpt,100,0,5);
    else hpt[2][iS] = new TH1F(name_hpt,title_hpt,1,0,5);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_mult_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of %s for events after multiplicity selection;p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[3][iS] = new TH1F(name_hpt,title_hpt,100,0,5);
    else hpt[3][iS] = new TH1F(name_hpt,title_hpt,1,0,5);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_acc_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of %s for events after |#eta|<0.8 cut;p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[4][iS] = new TH1F(name_hpt,title_hpt,100,0,5);
    else hpt[4][iS] = new TH1F(name_hpt,title_hpt,1,0,5);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_reco_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reconstructed %s;p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[5][iS] = new TH1F(name_hpt,title_hpt,100,0,5);
    else hpt[5][iS] = new TH1F(name_hpt,title_hpt,1,0,5);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_recoTof_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reconstructed %s and matched to the TOF;p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[6][iS] = new TH1F(name_hpt,title_hpt,100,0,5);
    else hpt[6][iS] = new TH1F(name_hpt,title_hpt,1,0,5);
  }

  Char_t name_fptRecoVsTrue[200];
  Char_t title_fptRecoVsTrue[200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fptRecoVsTrue,200,"fptRecoVsTrue_%s",nameSpec[iS]);
    snprintf(title_fptRecoVsTrue,200,"p_{T}^{reco.} vs p_{T}^{true} of %s;p_{T}^{reco.} (GeV/c);p_{T}^{reco.} - p_{T}^{true} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fptRecoVsTrue[0][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,100,0,5,300,-0.6,0.3);//320,-0.6,0.2
    else fptRecoVsTrue[0][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,1,0,5,1,-0.6,0.3);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fptRecoVsTrue,200,"fptRecoVsTrue_prim_%s",nameSpec[iS]);
    snprintf(title_fptRecoVsTrue,200,"p_{T}^{reco.} vs p_{T}^{true} of %s (primaries);p_{T}^{reco.} (GeV/c);p_{T}^{reco.} - p_{T}^{true} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fptRecoVsTrue[1][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,100,0,5,300,-0.6,0.3);//320,-0.6,0.2
    else fptRecoVsTrue[1][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,1,0,5,1,-0.6,0.3);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fptRecoVsTrue,200,"fptRecoVsTrue_sec_%s",nameSpec[iS]);
    snprintf(title_fptRecoVsTrue,200,"p_{T}^{reco.} vs p_{T}^{true} of %s (secondaries);p_{T}^{reco.} (GeV/c);p_{T}^{reco.} - p_{T}^{true} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fptRecoVsTrue[2][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,100,0,5,300,-0.6,0.3);//320,-0.6,0.2
    else fptRecoVsTrue[2][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,1,0,5,1,-0.6,0.3);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fptRecoVsTrue,200,"fptRecoVsTrue_done_%s",nameSpec[iS]);
    snprintf(title_fptRecoVsTrue,200,"p_{T}^{reco.} vs p_{T}^{true} of %s (correction already applied);p_{T}^{reco.} (GeV/c);p_{T}^{reco.} - p_{T}^{true} (GeV/c)",nameSpec[iS]);
    if(iS==5 || iS==5+9) fptRecoVsTrue[3][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,100,0,5,300,-0.6,0.3);//320,-0.6,0.2
    else fptRecoVsTrue[3][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,1,0,5,1,-0.6,0.3);
  }

  Char_t name_prptRecoVsTrue[200];
  Char_t title_prptRecoVsTrue[200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_prptRecoVsTrue,200,"prptRecoVsTrue_%s",nameSpec[iS]);
    snprintf(title_prptRecoVsTrue,200,"p_{T}^{reco.} vs p_{T}^{true} of %s;p_{T}^{reco.} (GeV/c);#LT p_{T}^{reco.} - p_{T}^{true} #GT (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) prptRecoVsTrue[0][iS] = new TProfile(name_prptRecoVsTrue,title_prptRecoVsTrue,100,0,5,-0.3,0.3,"");
    else prptRecoVsTrue[0][iS] = new TProfile(name_prptRecoVsTrue,title_prptRecoVsTrue,1,0,5,-0.3,0.3,"");
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_prptRecoVsTrue,200,"prptRecoVsTrue_prim_%s",nameSpec[iS]);
    snprintf(title_prptRecoVsTrue,200,"p_{T}^{reco.} vs p_{T}^{true} of %s (primaries);p_{T}^{reco.} (GeV/c);#LT p_{T}^{reco.} - p_{T}^{true} #GT (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) prptRecoVsTrue[1][iS] = new TProfile(name_prptRecoVsTrue,title_prptRecoVsTrue,100,0,5,-0.3,0.3,"");
    else prptRecoVsTrue[1][iS] = new TProfile(name_prptRecoVsTrue,title_prptRecoVsTrue,1,0,5,-0.3,0.3,"");
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_prptRecoVsTrue,200,"prptRecoVsTrue_sec_%s",nameSpec[iS]);
    snprintf(title_prptRecoVsTrue,200,"p_{T}^{reco.} vs p_{T}^{true} of %s (secondaries);p_{T}^{reco.} (GeV/c);#LT p_{T}^{reco.} - p_{T}^{true} #GT (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) prptRecoVsTrue[2][iS] = new TProfile(name_prptRecoVsTrue,title_prptRecoVsTrue,100,0,5,-0.3,0.3,"");
    else prptRecoVsTrue[2][iS] = new TProfile(name_prptRecoVsTrue,title_prptRecoVsTrue,1,0,5,-0.3,0.3,"");
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_prptRecoVsTrue,200,"prptRecoVsTrue_done_%s",nameSpec[iS]);
    snprintf(title_prptRecoVsTrue,200,"p_{T}^{reco.} vs p_{T}^{true} of %s (correction already applied);p_{T}^{reco.} (GeV/c);#LT p_{T}^{reco.} - p_{T}^{true} #GT (GeV/c)",nameSpec[iS]);
    if(iS==5 || iS==5+9) prptRecoVsTrue[3][iS] = new TProfile(name_prptRecoVsTrue,title_prptRecoVsTrue,100,0,5,-0.3,0.3,"");
    else prptRecoVsTrue[3][iS] = new TProfile(name_prptRecoVsTrue,title_prptRecoVsTrue,1,0,5,-0.3,0.3,"");
  }
  
  Char_t name_fDca[18][200];
  Char_t title_fDca[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_prim_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of primary %s;DCA_{xy} (cm);DCA_{z} (cm);p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==5 || iS==5+9) fDca[0][0][iS]=new TH3F(name_fDca[iS],title_fDca[iS],250,-1.0,1.0,68,-3.4,3.4,30,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
    else fDca[0][0][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-3.4,3.4,1,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_secMat_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of secondary (from material) %s;DCA_{xy} (cm);DCA_{z} (cm);p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==5 || iS==5+9) fDca[0][1][iS]=new TH3F(name_fDca[iS],title_fDca[iS],250,-1.0,1.0,68,-3.4,3.4,30,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
    else fDca[0][1][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-3.4,3.4,1,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_secWeak_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of secondary (from weak decay) %s;DCA_{xy} (cm);DCA_{z} (cm);p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==5 || iS==5+9) fDca[0][2][iS]=new TH3F(name_fDca[iS],title_fDca[iS],250,-1.0,1.0,68,-3.4,3.4,30,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
    else fDca[0][2][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-3.4,3.4,1,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
  }
  
  //TOF matching required:
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_prim_tof_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of primary %s (matched to the TOF);DCA_{xy} (cm);DCA_{z} (cm);p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==5 || iS==5+9) fDca[1][0][iS]=new TH3F(name_fDca[iS],title_fDca[iS],250,-1.0,1.0,68,-3.4,3.4,30,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
    else fDca[1][0][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-3.4,3.4,1,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_secMat_tof_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of secondary (from material) %s (matched to the TOF);DCA_{xy} (cm);DCA_{z} (cm);p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==5 || iS==5+9) fDca[1][1][iS]=new TH3F(name_fDca[iS],title_fDca[iS],250,-1.0,1.0,68,-3.4,3.4,30,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
    else fDca[1][1][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-3.4,3.4,1,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_secWeak_tof_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of secondary (from weak decay) %s (matched to the TOF);DCA_{xy} (cm);DCA_{z} (cm);p_{T}/|z| (GeV/c)",nameSpec[iS]);
    if(iS==5 || iS==5+9) fDca[1][2][iS]=new TH3F(name_fDca[iS],title_fDca[iS],250,-1.0,1.0,68,-3.4,3.4,30,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
    else fDca[1][2][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-3.4,3.4,1,0,3);//2000,-3.4,3.4,64,-3.4,3.4,30,0,3);
  }
  
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
    snprintf(title_fNsigmaTPC[iS],200,"n#sigma_{TPC} (%s) filled with all particles;p_{T}/|z| (GeV/c);n_{#sigma_{TPC}}^{%s}",nameSpec[iS],nameSpec[iS]);
    if(iS==5 || iS==5+9) fNsigmaTPC[0][iS] = new TH2F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],200,0,10,200,-10,10);
    else fNsigmaTPC[0][iS] = new TH2F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],1,0,10,1,-10,10);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fNsigmaTPC[iS],200,"NsigmaTPC_id_%s",nameSpec[iS]);
    snprintf(title_fNsigmaTPC[iS],200,"n#sigma_{TPC} (%s) filled only with %s (primaries and secondaries);p_{T}/|z| (GeV/c);n_{#sigma_{TPC}}^{%s}",nameSpec[iS],nameSpec[iS],nameSpec[iS]);
    if(iS==5 || iS==5+9) fNsigmaTPC[1][iS] = new TH2F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],200,0,10,200,-10,10);
    else fNsigmaTPC[1][iS] = new TH2F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],1,0,10,1,-10,10);
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
    if(iS==5 || iS==5+9) fM2vspt[0][iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],120,0,6,300,0,6);
    else fM2vspt[0][iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],1,0,6,1,0,6);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fM2vspt[iS],200,"fM2vspt_id_%s",nameSpec[iS]);
    snprintf(title_fM2vspt[iS],200,"m^{2}_{TOF} of %s (3#sigma TPC dE/dx cut on %s);p_{T}/|z| (GeV/c); m^{2}_{TOF} (GeV^{2}/c^{4})",nameSpec[iS],nameSpec[iS]);
    if(iS==5 || iS==5+9) fM2vspt[1][iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],120,0,6,300,0,6);
    else fM2vspt[1][iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],1,0,6,1,0,6);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fM2vspt[iS],200,"fM2vspt_id_noTPC_%s",nameSpec[iS]);
    snprintf(title_fM2vspt[iS],200,"m^{2}_{TOF} of %s;p_{T}/|z| (GeV/c); m^{2}_{TOF} (GeV^{2}/c^{4})",nameSpec[iS]);
    if(iS==5 || iS==5+9) fM2vspt[2][iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],120,0,6,300,0,6);
    else fM2vspt[2][iS] = new TH2F(name_fM2vspt[iS],title_fM2vspt[iS],1,0,6,1,0,6);
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
  
  //event characterization:
  fList->Add(htriggerMask);
  fList->Add(htriggerMask_noMB);
  fList->Add(hpdg);
  fList->Add(hNspdVertex);
  fList->Add(hpileUp);
  fList->Add(hzvertex);
  fList->Add(fNtrVsMult);
  fList->Add(prNtrVsMult);
  fList->Add(hmult_tot);
  fList->Add(hmult);
  fList->Add(hNtrack);

  //generated particles:
  for(Int_t i=0;i<2;i++) fList->Add(hpt[0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(hpt[0][5+9*i]);
  //for triggered events:
  for(Int_t i=0;i<2;i++) fList->Add(hpt[1][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(hpt[1][5+9*i]);
  //after v-vertex cut (10cm):
  for(Int_t i=0;i<2;i++) fList->Add(hpt[2][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(hpt[2][5+9*i]);
  //after multiplicity cut:
  for(Int_t i=0;i<2;i++) fList->Add(hpt[3][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(hpt[3][5+9*i]);
  //after eta cut (0.8):
  for(Int_t i=0;i<2;i++) fList->Add(hpt[4][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(hpt[4][5+9*i]);
  //reconstructed tracks:
  for(Int_t i=0;i<2;i++) fList->Add(hpt[5][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(hpt[5][5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fDca[0][0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[0][1][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[0][2][4+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fDca[0][0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[0][1][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[0][2][5+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(prptRecoVsTrue[0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(prptRecoVsTrue[0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[1][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[1][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(prptRecoVsTrue[1][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(prptRecoVsTrue[1][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[2][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[2][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(prptRecoVsTrue[2][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(prptRecoVsTrue[2][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[3][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(prptRecoVsTrue[3][5+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fdEdxVSp[i]);
  for(Int_t i=0;i<9;i++) fList->Add(hDeDxExp[i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[1][5+9*i]);
   
  //TOF matched:
  for(Int_t i=0;i<2;i++) fList->Add(hpt[6][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(hpt[6][5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fDca[1][0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[1][1][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[1][2][4+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fDca[1][0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[1][1][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[1][2][5+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fBetaTOFvspt[i]);
  for(Int_t i=0;i<9;i++) fList->Add(hBetaExp[i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2tof[i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[1][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[2][5+9*i]);
   
  // Post output data.
  PostData(1, fList);

  fPPVsMultUtils = new AliPPVsMultUtils();

}
//______________________________________________________________________________
void AliAnalysisMCNuclMult::UserExec(Option_t *) 
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
  
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if(!eventHandler) return;
  AliMCEvent* mcEvent = eventHandler->MCEvent();
  AliStack* fStack = mcEvent->Stack();
  
  //---------------------- EVENT analysis
  Bool_t isPhysSelected = kFALSE;
  /*for(Int_t i=0;i<32;i++) {
    Int_t bit=(1<<i);
    if(inputHandler->IsEventSelected() & bit) htriggerMask->Fill(i);
    }*/
  
  if(multiplicityMin!=-999) isPhysSelected = ((inputHandler->IsEventSelected() & AliVEvent::kMB) || (inputHandler->IsEventSelected() & AliVEvent::kHighMult));
  else isPhysSelected = (inputHandler->IsEventSelected() & AliVEvent::kMB);
  //if(!isPhysSelected) return;

  const AliVVertex* vtxEVENT = fEvent->GetPrimaryVertex();
  //event must have at least an SPD-determined primary vertex
  Int_t NevSpd=-1;
  NevSpd=vtxEVENT->GetNContributors();
  //hNspdVertex->Fill(NevSpd);
  //if(NevSpd<1) return;

  //event must not be tagged as pileup
  Bool_t isPileUpSpd=kFALSE;
  isPileUpSpd=((AliESDEvent *)fEvent)->IsPileupFromSPD();
  //hpileUp->Fill(isPileUpSpd);
  //if(isPileUpSpd) return;

  //event must have a primary vertex located within |z| < 10 cm
  Float_t zvtx=999.9;
  zvtx = vtxEVENT->GetZ();
  //hzvertex->Fill(zvtx);
  //if(TMath::Abs(zvtx)>10) return;
  
  /*//--------------------Multiplicity determination (Mid-pseudo-rapidity estimator):
  Int_t mult=AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent*)fESD,AliESDtrackCuts::kTrackletsITSTPC,0.8);
  if(multiplicityMin!=-999) {
    if(mult<multiplicityMin || mult>multiplicityMax) return;
  }
  hmult->Fill(mult);*/

  Int_t nMCpart = mcEvent->GetNumberOfTracks();
  //Int_t nMCprimary =  mcEvent->GetNumberOfPrimaries();
  
  //Int_t stdPdg[9] = {11,13,211,321,2212,10020,10030,20030,20040};//e,mu,pi,K,p,d,t,He3,He4
  
  //------------------------------- Loop on MC particles (I)

  for(Int_t iM=0;iM<nMCpart;iM++){
    
    Int_t Pdg=-999999;
    Short_t charge=0;
    Int_t label;
    Bool_t isPrimary = kFALSE;
    Int_t kSpec=-1;
    Double_t rapidity;
    Double_t eta;
    Double_t pt;
    
    AliVParticle *mcpart = (AliVParticle *) mcEvent->GetTrack(iM);
    Pdg = mcpart->PdgCode();
    //for nuclei: e.g deuteron: 1000010020-1000000000=10020
    if(Pdg>1000000000)Pdg -= 1000000000;
    else if (Pdg<-1000000000) Pdg += 1000000000;
    
    for(Int_t iS=0;iS<9;iS++) {
      if(TMath::Abs(Pdg)!=stdPdg[iS]) continue;
      kSpec=iS;
    }
    if(kSpec<0) continue;
    
    charge = mcpart->Charge();
    label = mcpart->GetLabel();

    isPrimary = fStack->IsPhysicalPrimary(label);
    //isPrimary = (iM < nMCprimary);
    //cut on secondary particles:
    if(!isPrimary) continue;
    rapidity = mcpart->Y();
    eta = mcpart->Eta();
    //rapidity cut:
    if(TMath::Abs(rapidity)>0.5) continue;
    pt=mcpart->Pt();

    hpdg->Fill(Pdg);

    if(charge>0) {
      hpt[0][kSpec]->Fill(pt);
    }
    else if(charge<0) {
      hpt[0][kSpec+9]->Fill(pt);
    }
      
    //----------Trigger condition:
    if(!isPhysSelected) continue;
    if(NevSpd<1) continue;
    if(isPileUpSpd) continue;

    if(charge>0) {
      hpt[1][kSpec]->Fill(pt);
    }
    else if(charge<0) {
      hpt[1][kSpec+9]->Fill(pt);
    }
    
    //----------beam-induced background cut:
    if(TMath::Abs(zvtx)>10) continue;

    if(charge>0) {
      hpt[2][kSpec]->Fill(pt);
    }
    else if(charge<0) {
      hpt[2][kSpec+9]->Fill(pt);
    }
  
  }//-------------------------------------- Loop on MC particles (I) (end)
  
  //---------------------------- Cut on events -----------

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
  
  if(!isPhysSelected) return;
  
  hNspdVertex->Fill(NevSpd);
  if(NevSpd<1) return;

  hpileUp->Fill(isPileUpSpd);
  if(isPileUpSpd) return;

  hzvertex->Fill(zvtx);
  if(TMath::Abs(zvtx)>10) return;

  //---------------------------- Cut on events (end) -----------
  
  //--------------------Multiplicity determination (Mid-pseudo-rapidity estimator):
 
  Float_t mult=-99;
  if(iMultEstimator==0) {
    //fPPVsMultUtils = new AliPPVsMultUtils();
    mult=fPPVsMultUtils->GetMultiplicityPercentile(fEvent,"V0M");
  }
  else if(iMultEstimator==1) mult=(Float_t)AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent*)fESD,AliESDtrackCuts::kTrackletsITSTPC,0.8);
  
  Int_t nTracks = fEvent->GetNumberOfTracks();
  fNtrVsMult->Fill(mult,nTracks);
  prNtrVsMult->Fill(mult,nTracks);
  
  hmult_tot->Fill(mult);  
  if(multiplicityMin!=-999) {
    if(mult<multiplicityMin || mult>multiplicityMax) return;
  }
  hmult->Fill(mult);

  //------------------------------- Loop on MC particles (II)

  for(Int_t iM=0;iM<nMCpart;iM++){
    
    Int_t Pdg=-999999;
    Short_t charge=0;
    Int_t label;
    Bool_t isPrimary = kFALSE;
    Int_t kSpec=-1;
    Double_t rapidity;
    Double_t eta;
    Double_t pt;
    
    AliVParticle *mcpart = (AliVParticle *) mcEvent->GetTrack(iM);
    Pdg = mcpart->PdgCode();
    //for nuclei: e.g deuteron: 1000010020-1000000000=10020
    if(Pdg>1000000000)Pdg -= 1000000000;
    else if (Pdg<-1000000000) Pdg += 1000000000;
    
    for(Int_t iS=0;iS<9;iS++) {
      if(TMath::Abs(Pdg)!=stdPdg[iS]) continue;
      kSpec=iS;
    }
    if(kSpec<0) continue;
    
    charge = mcpart->Charge();
    label = mcpart->GetLabel();

    isPrimary = fStack->IsPhysicalPrimary(label);
    //isPrimary = (iM < nMCprimary);
    //cut on secondary particles:
    if(!isPrimary) continue;
    rapidity = mcpart->Y();
    eta = mcpart->Eta();
    //rapidity cut:
    if(TMath::Abs(rapidity)>0.5) continue;
    pt=mcpart->Pt();

    //----------multiplicity cut:
    if(charge>0) {
      hpt[3][kSpec]->Fill(pt);
    }
    else if(charge<0) {
      hpt[3][kSpec+9]->Fill(pt);
    }
    //----------pseudorapidity cut (acceptance):
    if(TMath::Abs(eta)>0.8) continue;

    if(charge>0) {
      hpt[4][kSpec]->Fill(pt);
    }
    else if(charge<0) {
      hpt[4][kSpec+9]->Fill(pt);
    }
      
  }//-------------------------------------- Loop on MC particles (II) (end)
  
  //Int_t nTracks = fEvent->GetNumberOfTracks();
  hNtrack->Fill(nTracks);
  
  //----------------------loop on reconstructed TRACKS-----------------------------
  for(Int_t iT = 0; iT < nTracks; iT++) { 
    
    AliVTrack* track = (AliVTrack *) fEvent->GetTrack(iT);
    
    if (!track){
      continue;
    }
    
    //---------- info on MC true:
    Int_t label = track->GetLabel();
    //Int_t labelT=label;
    label=TMath::Abs(label);
    AliVParticle* mcpart = (AliVParticle *) mcEvent->GetTrack(label);
    
    Int_t Pdg=-999999;
    Short_t t_charge;
    Int_t t_label;
    Bool_t isPrimary = kFALSE;
    Bool_t isSecMat = kFALSE;
    Bool_t isSecWeak = kFALSE;
    Int_t kSpec=-1;
    Double_t t_rapidity=-9999.; 
    Double_t t_eta=-9999.;
    Double_t t_pt=-9999.;
    Pdg = mcpart->PdgCode();
    //for nuclei: e.g deuteron: 01000010020-1000000000=10020
    if(Pdg>1000000000)Pdg -= 1000000000;
    else if (Pdg<-1000000000) Pdg += 1000000000;
    for(Int_t iS=0;iS<9;iS++) {
      if(TMath::Abs(Pdg)!=stdPdg[iS]) continue;
      kSpec=iS;
    }
    if(kSpec<0) continue;
    
    t_charge = mcpart->Charge();
    t_label = mcpart->GetLabel();
    isPrimary = fStack->IsPhysicalPrimary(t_label);
    isSecMat = fStack->IsSecondaryFromMaterial(t_label);
    isSecWeak = fStack->IsSecondaryFromWeakDecay(t_label);
    //isPrimary = (label < nMCprimary);
    t_rapidity = mcpart->Y();
    t_eta = mcpart->Eta();
    t_pt=mcpart->Pt();
    //---------- info on MC true (end)
    
    //rapidity cut:
    if(TMath::Abs(t_rapidity)>0.5) continue;
    
    //----------pseudorapidity cut (acceptance):
    if(TMath::Abs(t_eta)>0.8) continue;
    
    //------------------------- Track cuts ----------------------------
    if(!fTrackFilter->IsSelected(track)) continue;
   
    //track extrapolation to the primary vertex
    Double_t b[2] = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
    if (!track->PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 100., b, bCov))
      continue;
    Double_t DCAxy = b[0];
    Double_t DCAz = b[1];
   
    //---------- Track info:
    Double_t rapidity = track->Y();
    Double_t eta = track->Eta();
    Double_t phi= track->Phi();
    Short_t charge = track->Charge();
    Double_t p = track->P();
    Double_t pt = track->Pt();
    //tpc:
    Double_t dedx = fPIDResponse->GetTPCsignalTunedOnData(track);//track->GetTPCsignal();
    Double_t pTPC = track->GetTPCmomentum();
    //tof:
    Bool_t kTOF=kFALSE;
    kTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    Double_t tof  = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
    Double_t beta = 0.0;
    Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};

    //Cut on the DCAxy
    if(TMath::Abs(DCAxy)>1.0) continue;

    this->FillDca(DCAxy, DCAz, Pdg, t_charge, isPrimary, isSecMat, isSecWeak, t_pt, kTOF);

    //Cut on the DCAz
    if(TMath::Abs(DCAz)>1.0) continue;

    //------------------------- Track cuts (end) -------------------------
    
    this->ForPtCorr(pt, t_pt, Pdg, t_charge, isPrimary);
    this->PtCorrection(pt, Pdg, t_charge);
    this->CheckPtCorr(pt, t_pt, Pdg, t_charge);
    
    //-------------------------- TPC info ---------------------------------
    
    //reconstructed PRIMARY particles:
    if(isPrimary) {

      if(t_charge>0) {
	hpt[5][kSpec]->Fill(t_pt);
      }
      else if(t_charge<0) {
	hpt[5][kSpec+9]->Fill(t_pt);
      }
            
    }// reconstructed PRIMARY particles (end)

    if(t_charge>0) fdEdxVSp[0]->Fill(pTPC,dedx);
    else if(t_charge<0) fdEdxVSp[1]->Fill(pTPC,dedx);
    
    Double_t nsigmaTPC[9];
    Double_t expdedx[9];
    
    Int_t stdFlagPid[9] = {1,2,4,8,16,32,64,128,256};//e,mu,pi,K,p,d,t,He3,He4
    Int_t FlagPid = 0;
    
    for(Int_t iS=0;iS<9;iS++){
      expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
      hDeDxExp[iS]->Fill(pTPC,expdedx[iS]);
      nsigmaTPC[iS] = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) iS);
      
      if(t_charge>0) fNsigmaTPC[0][iS]->Fill(t_pt,nsigmaTPC[iS]);
      else if(t_charge<0) fNsigmaTPC[0][iS+9]->Fill(t_pt,nsigmaTPC[iS]);

      if(TMath::Abs(nsigmaTPC[iS])<3) {
	FlagPid += ((Int_t)TMath::Power(2,iS));
      }
      
    }

    //using Pdg label:
    if(t_charge>0) fNsigmaTPC[1][kSpec]->Fill(t_pt,nsigmaTPC[kSpec]);
    else if(t_charge<0) fNsigmaTPC[1][kSpec+9]->Fill(t_pt,nsigmaTPC[kSpec]);

    //-------------------------- TOF matching required---------------------
    if(!kTOF) continue;
      
    //reconstructed PRIMARY particles-II:
    if(isPrimary) {
      if(t_charge>0) {
	hpt[6][kSpec]->Fill(t_pt);
      }
      else if(t_charge<0) {
	hpt[6][kSpec+9]->Fill(t_pt);
      }
    }// reconstructed PRIMARY particles-II(end)

    Double_t exptimes[9];
    track->GetIntegratedTimes(exptimes);
    //Integrated times of nuclei:
    for(Int_t iN=5;iN<9;iN++) {
      exptimes[iN] = exptimes[4]*exptimes[4]*(massOverZ[iN]*massOverZ[iN]/p/p+1)/(massOverZ[4]*massOverZ[4]/p/p+1);
      exptimes[iN] = TMath::Sqrt(exptimes[iN]);
    }  
    
    for(Int_t iS=0;iS<9;iS++){
      if(exptimes[iS]<1e-12) continue;
      hBetaExp[iS]->Fill(t_pt,exptimes[0]/exptimes[iS]);
    }
    
    beta=exptimes[0];
    if(tof<1e-12) continue;
    beta=beta/tof;//beta = L/tof/c = t_e/tof
    
    if(charge>0) fBetaTOFvspt[0]->Fill(t_pt,beta);
    else fBetaTOFvspt[1]->Fill(t_pt,beta);
    
    Double_t gamma2=1.-(beta*beta);
    if(gamma2<1e-12) continue;
    gamma2=1./gamma2;
    
    Double_t mass2=0.;
    mass2=beta*beta*gamma2;
    if(mass2<1e-12) continue;
    mass2=(p*p)/mass2;
    
    if(charge>0) fM2tof[0]->Fill(t_pt,mass2);
    else fM2tof[1]->Fill(t_pt,mass2);
    
    for(Int_t iS=0;iS<9;iS++){
      if(FlagPid & stdFlagPid[iS]) {
	if(t_charge>0) {
	  fM2vspt[0][iS]->Fill(t_pt,mass2);
	}	  
	else if(t_charge<0) {
	  fM2vspt[0][iS+9]->Fill(t_pt,mass2);
	}
      }
    }
    
    //TPC info and Pdg label used:
    if(FlagPid & stdFlagPid[kSpec]) {
      if(t_charge>0) {
	fM2vspt[1][kSpec]->Fill(t_pt,mass2);
      }	  
      else if(t_charge<0) {
	fM2vspt[1][kSpec+9]->Fill(t_pt,mass2);
      }
    }
    
    //using Pdg label:
    if(t_charge>0) {
      fM2vspt[2][kSpec]->Fill(t_pt,mass2);
    }	  
    else if(t_charge<0) {
      fM2vspt[2][kSpec+9]->Fill(t_pt,mass2);
    }
    
  }//loop on reconstructed TRACKS (end)
}//end loop on the events
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::ForPtCorr(Double_t pt, Double_t t_pt, Int_t Pdg, Short_t t_charge, Bool_t isPrimary) { 

  Int_t kSpec=-1;
  for(Int_t iS=0;iS<9;iS++) {
    if(TMath::Abs(Pdg)!=stdPdg[iS]) continue;
    kSpec=iS;
  }
  if(kSpec<0) return;

  if(kSpec>6) {//He3, He4
    pt*=2;
    t_pt*=2;
  }
  
  if(t_charge>0) {
    fptRecoVsTrue[0][kSpec]->Fill(pt,pt-t_pt);
    prptRecoVsTrue[0][kSpec]->Fill(pt,pt-t_pt);
    if(isPrimary) {
      fptRecoVsTrue[1][kSpec]->Fill(pt,pt-t_pt);
      prptRecoVsTrue[1][kSpec]->Fill(pt,pt-t_pt);
    }
    else {
      fptRecoVsTrue[2][kSpec]->Fill(pt,pt-t_pt);
      prptRecoVsTrue[2][kSpec]->Fill(pt,pt-t_pt);
    }
  }
  else if(t_charge<0) {
    fptRecoVsTrue[0][kSpec+9]->Fill(pt,pt-t_pt);
    prptRecoVsTrue[0][kSpec+9]->Fill(pt,pt-t_pt);
    if(isPrimary) {
      fptRecoVsTrue[1][kSpec+9]->Fill(pt,pt-t_pt);
      prptRecoVsTrue[1][kSpec+9]->Fill(pt,pt-t_pt);
    }
    else {
      fptRecoVsTrue[2][kSpec+9]->Fill(pt,pt-t_pt);
      prptRecoVsTrue[2][kSpec+9]->Fill(pt,pt-t_pt);
    }
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::PtCorrection(Double_t &pt, Int_t Pdg, Short_t t_charge) {

  Int_t kSpec=-1;
  for(Int_t iS=0;iS<9;iS++) {
    if(TMath::Abs(Pdg)!=stdPdg[iS]) continue;
    kSpec=iS;
  }
  if(kSpec<0) return;

  if(kSpec==5) {
    if(t_charge>0) {
      pt=pt-fpTcorr[0]->Eval(pt);
    }
    else if(t_charge<0) {
      pt=pt-fpTcorr[1]->Eval(pt);
    }
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::CheckPtCorr(Double_t pt, Double_t t_pt, Int_t Pdg, Short_t t_charge) {

  Int_t kSpec=-1;
  for(Int_t iS=0;iS<9;iS++) {
    if(TMath::Abs(Pdg)!=stdPdg[iS]) continue;
    kSpec=iS;
  }
  if(kSpec<0) return;
  
  if(kSpec==5) {
    if(t_charge>0) {
      fptRecoVsTrue[3][kSpec]->Fill(pt,pt-t_pt);
      prptRecoVsTrue[3][kSpec]->Fill(pt,pt-t_pt);
    }
    else if(t_charge<0) {
      fptRecoVsTrue[3][kSpec+9]->Fill(pt,pt-t_pt);
      prptRecoVsTrue[3][kSpec+9]->Fill(pt,pt-t_pt);
    }
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::FillDca(Double_t DCAxy, Double_t DCAz, Int_t Pdg, Short_t t_charge, Bool_t isPrimary, Bool_t isSecMat, Bool_t isSecWeak, Double_t t_pt, Bool_t kTOF) {
  
  //Int_t stdPdg[9] = {11,13,211,321,2212,10020,10030,20030,20040};//e,mu,pi,K,p,d,t,He3,He4
  Int_t kSpec=-1;
  for(Int_t iS=0;iS<9;iS++) {
    if(TMath::Abs(Pdg)!=stdPdg[iS]) continue;
    kSpec=iS;
  }
  if(kSpec<0) return;
  
  if(isPrimary) {
    if(t_charge>0) {
      fDca[0][0][kSpec]->Fill(DCAxy,DCAz,t_pt);
    }     
    else if(t_charge<0) {
      fDca[0][0][kSpec+9]->Fill(DCAxy,DCAz,t_pt);
    }    
  }
  else if(isSecMat){
    if(t_charge>0) {
      fDca[0][1][kSpec]->Fill(DCAxy,DCAz,t_pt);
    }     
    else if(t_charge<0) {
      fDca[0][1][kSpec+9]->Fill(DCAxy,DCAz,t_pt);
    }    
  }
  else if(isSecWeak){
    if(t_charge>0) {
      fDca[0][2][kSpec]->Fill(DCAxy,DCAz,t_pt);
    }     
    else if(t_charge<0) {
      fDca[0][2][kSpec+9]->Fill(DCAxy,DCAz,t_pt);
    }    
  }

  if(!kTOF) return; //TOF matching required
  if(isPrimary) {
    if(t_charge>0) {
      fDca[1][0][kSpec]->Fill(DCAxy,DCAz,t_pt);
    }     
    else if(t_charge<0) {
      fDca[1][0][kSpec+9]->Fill(DCAxy,DCAz,t_pt);
    }    
  }
  else if(isSecMat){
    if(t_charge>0) {
      fDca[1][1][kSpec]->Fill(DCAxy,DCAz,t_pt);
    }     
    else if(t_charge<0) {
      fDca[1][1][kSpec+9]->Fill(DCAxy,DCAz,t_pt);
    }    
  }
  else if(isSecWeak){
    if(t_charge>0) {
      fDca[1][2][kSpec]->Fill(DCAxy,DCAz,t_pt);
    }     
    else if(t_charge<0) {
      fDca[1][2][kSpec+9]->Fill(DCAxy,DCAz,t_pt);
    }    
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::Terminate(Option_t *) { 
  printf("Terminate()\n");
}
