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
  hpdg(NULL),
  htemp(NULL),
  hCheckTrackSel(NULL)
{
  for(Int_t i=0;i<18;i++) stdPdg[i] = 0;
  fList->SetName("results");
  //fList->SetOwner();
}
//______________________________________________________________________________
AliAnalysisMCNuclMult::AliAnalysisMCNuclMult(const char *name):
  AliAnalysisTaskSE(name),
  fAOD(NULL), 
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
  hpdg(NULL),
  htemp(NULL),
  hCheckTrackSel(NULL)
{
  for(Int_t i=0;i<18;i++) stdPdg[i] = 0;
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
  stdPdg[0+9] = -11; stdPdg[1+9] = -13; stdPdg[2+9] = -211; stdPdg[3+9] = -321; stdPdg[4+9] = -2212; stdPdg[5+9] = -10020; stdPdg[6+9] = -10030; stdPdg[7+9] = -20030; stdPdg[8+9] = -20040;

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
  const Char_t *xaxisTitle2[7]={"kMB || kHighMult","IsINELgtZERO","IsAcceptedVertexPosition","IsNotPileupSPDInMultBins","HasNoInconsistentSPDandTrackVertices","IsEventSelected","InsideMultiplicityBin"};
for(Int_t i=0;i<7;i++) {
    hNevent->Fill(xaxisTitle2[i],0);
  }
  
  hV0mult = new TH1I("hV0mult","Multiplicity distribution (after cuts on event);V0M Multiplicity Percentile",1000,0,100);

  hTrackMult = new TH1I("hTrackMult","Mid-pseudo-rapidity estimator of multiplicity;kTrackletsITSTPC (|#eta|<0.8)",1000,0,1000);
  
  hpdg = new TH2I("hpdg","Pdg label of generated particles (after the event selection);Pdg label;isPrimary           isSecMat           isSecWeak   ",50082,-25041,25041,3,0,3);
  hpdg->GetYaxis()->SetNdivisions(105);
  
  htemp = new TH1I("htemp","Number of tracks per event (after cuts on event);N_{tracks}",2100,-100,2000);
    
  hCheckTrackSel = new TH1I("hCheckTrackSel","Number of tracks per event after the track selection",12,0,12);
  const Char_t *xaxisTitle3[12]={"Ntracks","nTPCclusters>=70","chi2perTPCcluster<=4","isTPCrefit","isITSrefit","nSPD>0","NoKinkDaughters","chi2perITScluster<=36","isPropagatedToDca","|DCAxy|<1","|DCAz|<2","|eta|<0.8"};
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
  
  hnTPCclusters[1] = new TH1I("hnTPCclusters_1","Number of TPC clusters (after track cuts, DCAz cut not yet applied)",200,0,200);
  hchi2TPC[1] = new TH1D("hchi2TPC_1","#chi^{2} per TPC cluster (after track cuts, DCAz cut not yet applied)",1000,0,100);
  hisTPCrefit[1] = new TH1I("hisTPCrefit_1","kTPCrefit (after track cuts, DCAz cut not yet applied)",2,0,2);
  hisITSrefit[1] = new TH1I("hisITSrefit_1","kITSrefit (after track cuts, DCAz cut not yet applied)",2,0,2);
  hnSPD[1] = new TH1I("hnSPD_1","Number of SPD rec. points (after track cuts, DCAz cut not yet applied)",3,0,3);
  hnKinkDaughters[1] = new TH1I("hnKinkDaughters_1","Number of Kink Daughters (after track cuts, DCAz cut not yet applied)",40,0,40);
  hsigmaToVtx[1] = new TH1D("hsigmaToVtx_1","Number of sigma to the vertex (after track cuts, DCAz cut not yet applied)",400,0,40);
  hchi2ITS[1] = new TH1D("hchi2ITS_1","#chi^{2} per ITS cluster (after track cuts, DCAz cut not yet applied)",1000,0,100);
    
  heta[1] = new TH1D("heta_1","#eta (after track cuts, DCAz cut not yet applied);#eta",200,-1,1);
  hisPropagatedToDca[1] = new TH1I("hisPropagatedToDca_1","kPropagatedToDca (after track cuts, DCAz cut not yet applied)",2,0,2);

  fdEdxVSp[0] = new TH2F("fdEdxVSp_pos","TPC dE/dx (positive charge); p_{TPC}/|z| (GeV/c); dE/dx_{TPC} (a.u.)",200,0,10,1000,0,2000);//100,500//500,2000
  fdEdxVSp[1] = new TH2F("fdEdxVSp_neg","TPC dE/dx (negative charge); p_{TPC}/|z| (GeV/c); dE/dx_{TPC} (a.u.)",200,0,10,1000,0,2000);

  Char_t name_hDeDxExp[9][200];
  Char_t title_hDeDxExp[9][200];
  for(Int_t iS=0;iS<9;iS++) {
    snprintf(name_hDeDxExp[iS],200,"hDeDxExp_%s",nameSpec[iS]);
    snprintf(title_hDeDxExp[iS],200,"Expected TPC dE/dx of %s;p_{TPC}/|z| (GeV/c);dE/dx_{TPC} (a.u.)",nameSpec[iS]);
    hDeDxExp[iS] = new TProfile(name_hDeDxExp[iS],title_hDeDxExp[iS],200,0,10,0,1000,"");//,500,0,5,0,1000,"");
  }

  //temporary plot:
  fBetaTOFvspt[0] = new TH2F("fBetaTOFvspt_pos","#beta_{TOF} (positive charge);p_{T}^{reco.}/|z| (GeV/c);#beta_{TOF}",200,0,10,260,0.4,1.05);//500,520//200,260
  fBetaTOFvspt[1] = new TH2F("fBetaTOFvspt_neg","#beta_{TOF} (negative charge);p_{T}^{reco.}/|z| (GeV/c);#beta_{TOF}",200,0,10,260,0.4,1.05);

  Char_t name_hBetaExp[9][200];
  Char_t title_hBetaExp[9][200];
  for(Int_t iS=0;iS<9;iS++) {
    snprintf(name_hBetaExp[iS],200,"hBetaTOFVSpt_Exp_%s",nameSpec[iS]);
    snprintf(title_hBetaExp[iS],200,"Expected #beta_{TOF} (%s);p_{T}^{reco.}/|z| (GeV/c); #beta_{TOF}",nameSpec[iS]);
    hBetaExp[iS] = new TProfile(name_hBetaExp[iS],title_hBetaExp[iS],200,0,10,0.4,1.05,"");
  } 
  //-------------

  Char_t name_hpt[200];
  Char_t title_hpt[200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_gen_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of generated primary %s for |y|<0.5.;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[0][0][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[0][0][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_gen_secMat_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of generated secondary %s from mat. for |y|<0.5.;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[0][1][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[0][1][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_gen_secWeak_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of generated secondary %s from w.d. for |y|<0.5.;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[0][2][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[0][2][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }

  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_reco_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reco. primary %s;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[1][0][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[1][0][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_reco_secMat_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from mat.;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[1][1][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[1][1][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_reco_secWeak_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from w.d.;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[1][2][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[1][2][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }

  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_matchedTof_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reco. primary %s matched at TOF;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[2][0][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[2][0][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_matchedTof_secMat_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from mat. matched at TOF;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[2][1][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[2][1][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_matchedTof_secWeak_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from w.d. matched at TOF;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[2][2][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[2][2][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }

  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_goodmatchTof_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reco. primary %s matched correctly at TOF;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[3][0][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[3][0][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_goodmatchTof_secMat_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from mat. matched correctly at TOF;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[3][1][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[3][1][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_hpt,200,"hp_goodmatchTof_secWeak_%s",nameSpec[iS]);
    snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from w.d. matched correctly at TOF;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) hpt[3][2][iS] = new TH1F(name_hpt,title_hpt,200,0,10);
    else hpt[3][2][iS] = new TH1F(name_hpt,title_hpt,1,0,10);
  }
    
  Char_t name_fptRecoVsTrue[200];
  Char_t title_fptRecoVsTrue[200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fptRecoVsTrue,200,"fptRecoVsTrue_%s",nameSpec[iS]);
    snprintf(title_fptRecoVsTrue,200,"%s;p_{T}^{reco.}/z (GeV/c);p_{T}^{reco.} - p_{T}^{true} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fptRecoVsTrue[0][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,100,0,5,300,-0.6,0.3);
    else fptRecoVsTrue[0][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,1,0,5,1,-0.6,0.3);
  }
  
  Char_t name_fDca[18][200];
  Char_t title_fDca[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_prim_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of primary %s;DCA_{xy} (cm);DCA_{z} (cm);p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fDca[0][iS]=new TH3F(name_fDca[iS],title_fDca[iS],200,-1.0,1.0,200,-2.0,2.0,100,0,5);
    else fDca[0][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-2.0,2.0,1,0,5);
  }
  
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_secMat_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of secondary (from material) %s;DCA_{xy} (cm);DCA_{z} (cm);p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fDca[1][iS]=new TH3F(name_fDca[iS],title_fDca[iS],200,-1.0,1.0,200,-2.0,2.0,100,0,5);
    else fDca[1][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-2.0,2.0,1,0,5);
  }
  
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDca[iS],200,"fDca_secWeak_%s",nameSpec[iS]);
    snprintf(title_fDca[iS],200,"DCA of secondary (from weak decay) %s;DCA_{xy} (cm);DCA_{z} (cm);p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fDca[2][iS]=new TH3F(name_fDca[iS],title_fDca[iS],200,-1.0,1.0,200,-2.0,2.0,100,0,5);
    else fDca[2][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-1.0,1.0,1,-2.0,2.0,1,0,5);
  }
  
  if(multMin<1e-12 && multMax>99) {//plots added only for integrated multiplicity
    fList->Add(htriggerMask);
    fList->Add(htriggerMask_noMB);
    fList->Add(hzvertex);
    fList->Add(hNevent);
  }

  fList->Add(hV0mult);
  fList->Add(hTrackMult);
  
  fList->Add(hpdg);
  fList->Add(htemp);
  
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
  
  for(Int_t i=0;i<2;i++) fList->Add(fdEdxVSp[i]);
  for(Int_t i=0;i<9;i++) fList->Add(hDeDxExp[i]);
  for(Int_t i=0;i<2;i++) fList->Add(fBetaTOFvspt[i]);
  for(Int_t i=0;i<9;i++) fList->Add(hBetaExp[i]);
  
  for(Int_t j=0;j<3;j++) {
    for(Int_t i=0;i<2;i++) fList->Add(hpt[0][j][4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(hpt[0][j][5+9*i]);
    
    for(Int_t i=0;i<2;i++) fList->Add(hpt[1][j][4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(hpt[1][j][5+9*i]);
    
    for(Int_t i=0;i<2;i++) fList->Add(hpt[2][j][4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(hpt[2][j][5+9*i]);
    
    for(Int_t i=0;i<2;i++) fList->Add(hpt[3][j][4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(hpt[3][j][5+9*i]);
  }
  
  for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[0][5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fDca[0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[1][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[2][4+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fDca[0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[1][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDca[2][5+9*i]);
  
  // Post output data.
  PostData(1, fList);

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
  
  Int_t nMCpart = mcEvent->GetNumberOfTracks();
  //------------------------------- Loop on MC particles
  for(Int_t iM=0;iM<nMCpart;iM++){
    
    AliVParticle *mcpart = (AliVParticle *) mcEvent->GetTrack(iM);
    
    Int_t Pdg = mcpart->PdgCode();
    //for nuclei: e.g. deuteron: 1000010020-1000000000=10020
    if(Pdg>1e+9) Pdg -= 1e+9;
    else if (Pdg<-1e+9) Pdg += 1e+9;
    Int_t kSpec=-1;
    for(Int_t iS=0;iS<18;iS++) {
      if(Pdg!=stdPdg[iS]) continue;
      kSpec=iS;
    }
    
    Int_t label = mcpart->GetLabel();
    Bool_t isPrimary = fStack->IsPhysicalPrimary(label);
    Bool_t isSecMat = fStack->IsSecondaryFromMaterial(label);
    Bool_t isSecWeak = fStack->IsSecondaryFromWeakDecay(label);
    Double_t rapidity = mcpart->Y();
    Double_t pt = mcpart->Pt();
    
    if(isPrimary) hpdg->Fill(Pdg,0);
    else if(isSecMat) hpdg->Fill(Pdg,1);
    else if(isSecWeak) hpdg->Fill(Pdg,2);
    
    if(TMath::Abs(rapidity)>0.5) continue;
    
    if(kSpec<0) continue;
    
    if(isPrimary) {
      hpt[0][0][kSpec]->Fill(pt);
    }
    else if(isSecMat) {
      hpt[0][1][kSpec]->Fill(pt); 
    }
    else if(isSecWeak) {
      hpt[0][2][kSpec]->Fill(pt); 
    }
    
  }//------------------------------- Loop on MC particles (end)
   
  Int_t nTracks = fEvent->GetNumberOfTracks();
  htemp->Fill(nTracks);
  //-------------------------------- Loop on reconstructed TRACKS
  for(Int_t iT=0;iT<nTracks;iT++) { 

    AliVTrack* track = (AliVTrack *) fEvent->GetTrack(iT);
    if (!track){
      continue;
    }
    
    Int_t label = track->GetLabel();
    label=TMath::Abs(label);
    
    //---------- info on MC true:
    AliVParticle *mcpart = (AliVParticle *) mcEvent->GetTrack(label);
    
    Int_t Pdg = mcpart->PdgCode();
    //for nuclei: e.g. deuteron: 1000010020-1000000000=10020
    if(Pdg>1e+9) Pdg -= 1e+9;
    else if (Pdg<-1e+9) Pdg += 1e+9;
    Int_t kSpec=-1;
    for(Int_t iS=0;iS<18;iS++) {
      if(Pdg!=stdPdg[iS]) continue;
      kSpec=iS;
    }

    Int_t t_label = mcpart->GetLabel();
    Bool_t isPrimary = fStack->IsPhysicalPrimary(label);
    Bool_t isSecMat = fStack->IsSecondaryFromMaterial(t_label);
    Bool_t isSecWeak = fStack->IsSecondaryFromWeakDecay(t_label);
    Double_t t_rapidity = mcpart->Y();
    Double_t t_pt = mcpart->Pt();
    
    //---------- info on MC true (end)

     //rapidity cut:
    if(TMath::Abs(t_rapidity)>0.5) continue;
    
    //------------------------- Track cuts:
    Double_t DCAxy, DCAz;
    if(!this->AcceptTrack(track, DCAxy, DCAz)) continue;
    
    this->FillDca(DCAxy, DCAz, t_pt, kSpec, isPrimary, isSecMat, isSecWeak);

    if(TMath::Abs(DCAz)>1.0) continue;
    //------------------------- Track cuts (end)

    Double_t pt = track->Pt(); //use only on ForPtCorr !
    this->ForPtCorr(pt, t_pt, kSpec);

    this->CheckTPCsignal(track);
    
    //reconstructed particles:
    if(kSpec>-1){
      if(isPrimary) {
	hpt[1][0][kSpec]->Fill(t_pt);
      }
      else if(isSecMat) {
	hpt[1][1][kSpec]->Fill(t_pt); 
      }
      else if(isSecWeak) {
	hpt[1][2][kSpec]->Fill(t_pt); 
      }
    }// reconstructed particles (end)

    Bool_t kTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    //-------------------------- TOF matching required---------------------
    if(!kTOF) continue;
    
    Int_t *toflabel = new Int_t[3];
    ((AliESDtrack *)track)->GetTOFLabel(toflabel);
    Bool_t isTOFgoodMatch = kFALSE;
    if(toflabel[0]==label) isTOFgoodMatch = kTRUE;
    
    this->CheckTOFsignal(track);
    
    if(kSpec>-1) {
      
      if(isPrimary) {
	hpt[2][0][kSpec]->Fill(t_pt);
	//good match at TOF
	if(isTOFgoodMatch) {
	  hpt[3][0][kSpec]->Fill(t_pt);
	}
      }
      else if(isSecMat) {
	hpt[2][1][kSpec]->Fill(t_pt);
	//good match at TOF
	if(isTOFgoodMatch) {
	  hpt[3][1][kSpec]->Fill(t_pt);
	}
      }
      else if(isSecWeak) {
	hpt[2][2][kSpec]->Fill(t_pt);
	//good match at TOF
	if(isTOFgoodMatch) {
	  hpt[3][2][kSpec]->Fill(t_pt);
	}
      }
    
    }
        
  }//----------------------loop on reconstructed TRACKS (end)
  
}//end loop on the events
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::Terminate(Option_t *) { 
  printf("Terminate()\n");
}
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::EventSelectionMonitor() {
  
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

	if(fPPVsMultUtils->IsNotPileupSPDInMultBins(fEvent)) {
	  hNevent->Fill(3);

	  if(fPPVsMultUtils->HasNoInconsistentSPDandTrackVertices(fEvent)) {
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
Bool_t AliAnalysisMCNuclMult::IsInsideMultiplicityBin(Float_t multiplicity) {

  if(multiplicity < multMin || multiplicity > multMax+1e-12) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisMCNuclMult::AcceptTrack(AliVTrack *track, Double_t &DCAxy, Double_t &DCAz) {
 
  //for the moment it checks for std. tracks only

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

  if(TMath::Abs(DCAz)>2.0) return kFALSE;    // temporarily 1cm cut moved in exec !

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
void AliAnalysisMCNuclMult::TrackSelectionMonitor(Int_t nTPCclusters, Double_t chi2TPC, Bool_t isTPCrefit, Bool_t isITSrefit, Int_t nSPD, Int_t nKinkDaughters, Double_t chi2ITS, Bool_t isPropagatedToDca, Double_t DCAxy, Double_t DCAz, Double_t eta) {

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
		 
		    if(TMath::Abs(DCAz)<2.0) {//temporary set to 2.0
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
void AliAnalysisMCNuclMult::FillDca(Double_t DCAxy, Double_t DCAz, Double_t t_pt, Int_t kSpec, Bool_t isPrimary, Bool_t isSecMat, Bool_t isSecWeak) {
  
  if(kSpec<0) return;
  
  if(isPrimary) {
    fDca[0][kSpec]->Fill(DCAxy,DCAz,t_pt);
  }     
  else if(isSecMat){
    fDca[1][kSpec]->Fill(DCAxy,DCAz,t_pt);
  }
  else if(isSecWeak){
    fDca[2][kSpec]->Fill(DCAxy,DCAz,t_pt);
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::ForPtCorr(Double_t pt, Double_t t_pt, Int_t kSpec) { 
  
  if(kSpec<0) return;
  
  if(kSpec==7 || kSpec==8 || kSpec==7+9 || kSpec==8+9) t_pt=t_pt/2;
  
  fptRecoVsTrue[0][kSpec]->Fill(pt,pt-t_pt);
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::CheckTPCsignal(AliVTrack *track) {

  Short_t charge = track->Charge();
  Double_t pTPC = track->GetTPCmomentum();
  Double_t dedx = fPIDResponse->GetTPCsignalTunedOnData(track);//track->GetTPCsignal();
  
  if(charge>0) fdEdxVSp[0]->Fill(pTPC,dedx);
  else if(charge<0) fdEdxVSp[1]->Fill(pTPC,dedx);
  Double_t expdedx[9];
  for(Int_t iS=0;iS<9;iS++){
    expdedx[iS] = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, (AliPID::EParticleType) iS, AliTPCPIDResponse::kdEdxDefault, kTRUE);
    hDeDxExp[iS]->Fill(pTPC,expdedx[iS]);
  }
  
  return;
}
//_____________________________________________________________________________
void AliAnalysisMCNuclMult::CheckTOFsignal(AliVTrack *track) {
  //-------------------------- TOF plots:
  //To do: p->t_p pt->t_pt, remember that the true variables are not divided by the charge Z!
  //For the moment I use only the reconstructed variables (are only "QA" plots)!
    
  Double_t pt = track->Pt();
  Double_t p = track->P();
  Short_t charge = track->Charge();

  Double_t tof  = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
  Double_t beta = 0.;
  Double_t massOverZ[9] = {0.000511,0.105658,0.139570,0.493677,0.938272,1.875612859,2.808921005,1.404195741,1.863689620};
  Double_t exptimes[9];
  track->GetIntegratedTimes(exptimes);
  //Integrated times of nuclei:
  for(Int_t iN=5;iN<9;iN++) {
    if(p<1e-12) continue;
    exptimes[iN] = exptimes[4]*exptimes[4]*(massOverZ[iN]*massOverZ[iN]/p/p+1)/(massOverZ[4]*massOverZ[4]/p/p+1);
    exptimes[iN] = TMath::Sqrt(exptimes[iN]);
  } 
  for(Int_t iS=0;iS<9;iS++){
    if(exptimes[iS]<1e-12) continue;
    hBetaExp[iS]->Fill(pt,exptimes[0]/exptimes[iS]);
  }
  if(tof>1e-12) {
    beta=exptimes[0];
    beta=beta/tof;//beta = L/tof/c = t_e/tof
  }
  if(charge>0) fBetaTOFvspt[0]->Fill(pt,beta);
  else if(charge<0) fBetaTOFvspt[1]->Fill(pt,beta);
  
  return;
}
