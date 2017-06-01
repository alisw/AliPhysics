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
#include "THnSparse.h"

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
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
//#include "AliVParticle.h"
#include "AliGenEventHeader.h"

ClassImp(AliAnalysisNuclMult)

//Pdg code of e,mu,pi,K,p,d,t,He3,He4:
const Int_t PdgStd[18]={ 11, 13, 211, 321, 2212, 10020, 10030, 20030, 20040,
			-11,-13,-211,-321,-2212,-10020,-10030,-20030,-20040};

const Float_t multMin[7]={  0, 0,  5, 10, 20, 40,  70};
const Float_t multMax[7]={100, 5, 10, 20, 40, 70, 100};

const Float_t trackletsMin[14]={   0, 0, 1, 4,  7, 10, 15, 20, 25, 30, 40, 50, 60,   70};
const Float_t trackletsMax[14]={1000, 1, 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 1000};

//_____________________________________________________________________________
AliAnalysisNuclMult::AliAnalysisNuclMult():
  AliAnalysisTaskSE(),
  isMC(kFALSE),
  fESD(NULL),
  fEvent(NULL),
  eventHandler(NULL),
  mcEvent(NULL),
  fStack(NULL),
  fESDtrackCuts(NULL),
  fPPVsMultUtils(NULL),
  fPIDResponse(NULL),
  fList(new TList()),
  nTPCclustersMin(70),
  nsigmaTPCMax(3.),
  DCAxyMax(0.5),
  DCAzMax(1.),
  hzvertex(NULL),
  hNevent(NULL),
  hV0mult(NULL),
  hTracklets(NULL),
  hTrackletsVsV0mult(NULL),
  hCheckTrackSel(NULL)
{
  fList->SetName("results");
}
//______________________________________________________________________________
AliAnalysisNuclMult::AliAnalysisNuclMult(const char *name):
  AliAnalysisTaskSE(name),
  isMC(kFALSE),
  fESD(NULL),
  fEvent(NULL),
  eventHandler(NULL),
  mcEvent(NULL),
  fStack(NULL),
  fESDtrackCuts(NULL),
  fPPVsMultUtils(NULL),
  fPIDResponse(NULL),
  fList(new TList()),
  nTPCclustersMin(70),
  nsigmaTPCMax(3.),
  DCAxyMax(0.5),
  DCAzMax(1.),
  hzvertex(NULL),
  hNevent(NULL),
  hV0mult(NULL),
  hTracklets(NULL),
  hTrackletsVsV0mult(NULL),
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

  fList->SetOwner(kTRUE);
  
  AliESDtrackCuts* esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE, 0);
  esdTrackCuts->SetMinNClustersTPC(nTPCclustersMin);
  esdTrackCuts->SetMaxDCAToVertexZ(1e+9);
  //DCA cuts set in the task
  fESDtrackCuts = esdTrackCuts;

  const Char_t nameSpec[18][30]={"e^{+}","#mu^{+}","#pi^{+}","K^{+}",     "p",      "d",      "t",      "^{3}He",      "^{4}He",
				 "e^{-}","#mu^{-}","#pi^{-}","K^{-}","#bar{p}","#bar{d}","#bar{t}","^{3}#bar{He}","^{4}#bar{He}"};
    
  const Int_t Nmultbins=20;
  const Float_t multV0Min=0, multV0Max=100;

  const Int_t Ntrackletsbins=13;
  const Float_t trackletsbins[Ntrackletsbins+1]={0, 1, 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 1000};

  const Int_t Nptbins=100;
  const Float_t ptMin=0, ptMax=10;

  const Int_t Nsigmabins=500;
  const Float_t sigmaMin=-50, sigmaMax=50;

  const Int_t NDCAzbins=100;
  const Int_t NDCAxybins=100;
  
  htriggerMask[0] = new TH1F("htriggerMask","Trigger mask (before the event selection)",32,0,32);
  const Char_t *xaxisTitle[32]={"kMB","kINT7","kMUON","kHighMult","kEMC1","kCINT5","kCMUS5","kMUSH7","kMUL7","kMUU7","kEMC7","kEMC8","kMUS7","kPHI1","kPHI7","kEMCEJE","kEMCEGA","kCentral","kSemiCentral","kDG5","kZED","kSPI7","kSPI","kINT8","kMuonSingleLowPt8","kMuonSingleHighPt8","kMuonLikeLowPt8","kMuonUnlikeLowPt8","kMuonUnlikeLowPt0","kTRD","kFastOnly","kUserDefined"};
  for(Int_t i=0;i<32;i++) {
    htriggerMask[0]->Fill(xaxisTitle[i],0);
  }
  
  htriggerMask[1] = new TH1F("htriggerMask_noMB","Trigger mask excluding MB events (before the event selection)",32,0,32);
  for(Int_t i=0;i<32;i++) {
    htriggerMask[1]->Fill(xaxisTitle[i],0);
  }

  hzvertex = new TH1F("hzvertex","z-vertex distribution. After (only) the trigger selection and the selection of INEL. collisions. Therefore it's filled before the other event cuts;z_{vtx} (cm)",200,-20,20);
 
  hNevent = new TH1F("hNevent","Event counter. To check the event selection compare the last bin with the number of hV0mult entries",7,0,7);
  const Char_t *xaxisTitle2[7]={"kMB","IsINELgtZERO","IsAcceptedVertexPosition","HasNoInconsistentSPDandTrackVertices","IsNotPileupSPDInMultBins","IsEventSelected","InsideMultiplicityBin"};
for(Int_t i=0;i<7;i++) {
    hNevent->Fill(xaxisTitle2[i],0);
  }
  
  hV0mult = new TH1F("hV0mult","Multiplicity distribution (after the event selection);V0M Multiplicity Percentile",100,0,100);

  hTracklets = new TH1F("hTracklets","Mid-pseudo-rapidity estimator of multiplicity;Number of tracklets",1000,0,1000);
  
  hTrackletsVsV0mult = new TH2F("hTrackletsVsV0mult",";V0M Multiplicity Percentile;Number of tracklets",100,0,100,1000,0,1000);

  hrapidity[0] = new TH1F("hrapidity_0","Rapidity (assuming deuteron mass) before |y|<0.5 cut;rapidity",1000,-5,5);
  hrapidity[1] = new TH1F("hrapidity_1","Rapidity (assuming deuteron mass) after |y|<0.5 cut;rapidity",1000,-5,5);
  if(isMC) {
    hrapidity[0]->SetTitle("Rapidity before |y|<0.5 cut");
    hrapidity[1]->SetTitle("Rapidity after |y|<0.5 cut");
  }
    
  hCheckTrackSel = new TH1F("hCheckTrackSel","Number of tracks per event after the track selection",12,0,12);
  const Char_t *xaxisTitle3[12]={"|y|<0.5",Form("nTPCclusters>=%d",nTPCclustersMin),"chi2perTPCcluster<=4","isTPCrefit","isITSrefit","nSPD>0","NoKinkDaughters","chi2perITScluster<=36","isPropagatedToDca",Form("|DCAxy|<%.1f",DCAxyMax),Form("|DCAz|<%.1f",DCAzMax),"|#eta|<0.8"};
  for(Int_t i=0;i<12;i++) {
    hCheckTrackSel->Fill(xaxisTitle3[i],0);
  }
  
  hnTPCclusters[0] = new TH1F("hnTPCclusters_0","Number of TPC clusters (before track cuts)",200,0,200);
  hchi2TPC[0] = new TH1F("hchi2TPC_0","#chi^{2} per TPC cluster (before track cuts)",1000,0,100);
  hisTPCrefit[0] = new TH1F("hisTPCrefit_0","kTPCrefit (before track cuts)",2,0,2);
  hisITSrefit[0] = new TH1F("hisITSrefit_0","kITSrefit (before track cuts)",2,0,2);
  hnSPD[0] = new TH1F("hnSPD_0","Number of SPD rec. points (before track cuts)",3,0,3);
  hnKinkDaughters[0] = new TH1F("hnKinkDaughters_0","Number of Kink Daughters (before track cuts)",40,0,40);
  hsigmaToVtx[0] = new TH1F("hsigmaToVtx_0","Number of sigma to the vertex (before track cuts)",400,0,40);
  hchi2ITS[0] = new TH1F("hchi2ITS_0","#chi^{2} per ITS cluster (before track cuts)",1000,0,100);
  
  heta[0] = new TH1F("heta_0","#eta (before track cuts);#eta",200,-1,1);
  hisPropagatedToDca[0] = new TH1F("hisPropagatedToDca_0","kPropagatedToDca (before track cuts)",2,0,2);
  
  hnTPCclusters[1] = new TH1F("hnTPCclusters_1","Number of TPC clusters (after track cuts)",200,0,200);
  hchi2TPC[1] = new TH1F("hchi2TPC_1","#chi^{2} per TPC cluster (after track cuts)",1000,0,100);
  hisTPCrefit[1] = new TH1F("hisTPCrefit_1","kTPCrefit (after track cuts)",2,0,2);
  hisITSrefit[1] = new TH1F("hisITSrefit_1","kITSrefit (after track cuts)",2,0,2);
  hnSPD[1] = new TH1F("hnSPD_1","Number of SPD rec. points (after track cuts)",3,0,3);
  hnKinkDaughters[1] = new TH1F("hnKinkDaughters_1","Number of Kink Daughters (after track cuts)",40,0,40);
  hsigmaToVtx[1] = new TH1F("hsigmaToVtx_1","Number of sigma to the vertex (after track cuts)",400,0,40);
  hchi2ITS[1] = new TH1F("hchi2ITS_1","#chi^{2} per ITS cluster (after track cuts)",1000,0,100);
    
  heta[1] = new TH1F("heta_1","#eta (after track cuts);#eta",200,-1,1);
  hisPropagatedToDca[1] = new TH1F("hisPropagatedToDca_1","kPropagatedToDca (after track cuts)",2,0,2);
  
  fptCorr[0] = new TF1("fptCorr_d","[0]-[1]*TMath::Exp(-[2]*x)",0,10);
  fptCorr[0]->FixParameter(0,-3.97081e-03);
  fptCorr[0]->FixParameter(1,4.94023e-01);
  fptCorr[0]->FixParameter(2,3.21308e+00);
  fptCorr[0]->GetXaxis()->SetTitle("p_{T}^{reco} (GeV/c)");
  fptCorr[0]->GetYaxis()->SetTitle("p_{T}^{reco} - p_{T}^{true} (GeV/c)");

  fdEdxVSp[0] = new TH2F("fdEdxVSp_pos","TPC dE/dx (positive charges); p_{TPC}/|z| (GeV/c); dE/dx_{TPC} (a.u.)",200,0,10,1000,0,2000);
  fdEdxVSp[1] = new TH2F("fdEdxVSp_neg","TPC dE/dx (negative charges); p_{TPC}/|z| (GeV/c); dE/dx_{TPC} (a.u.)",200,0,10,1000,0,2000);

  Char_t name_hDeDxExp[9][200];
  Char_t title_hDeDxExp[9][200];
  for(Int_t iS=0;iS<9;iS++) {
    snprintf(name_hDeDxExp[iS],200,"hDeDxExp_%s",nameSpec[iS]);
    snprintf(title_hDeDxExp[iS],200,"Expected TPC dE/dx of %s;p_{TPC}/|z| (GeV/c);dE/dx_{TPC} (a.u.)",nameSpec[iS]);
    hDeDxExp[iS] = new TProfile(name_hDeDxExp[iS],title_hDeDxExp[iS],200,0,10,0,1000,"");//,500,0,5,0,1000,"");
  }

  if(1) {
    const Int_t Ndim=4;
    Int_t bins[Ndim] = {Nsigmabins, Nptbins, Nmultbins, Ntrackletsbins};
    Double_t xmin[Ndim] = {sigmaMin, ptMin, multV0Min,    0};
    Double_t xmax[Ndim] = {sigmaMax, ptMax, multV0Max, 1000};
    Char_t name_fSparse[200];
    Char_t title_fSparse[200];
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fSparse,200,"fSparseNsigmaTPC_%s",nameSpec[iS]);
      snprintf(title_fSparse,200,"%s",nameSpec[iS]);
      TString axis[Ndim] = {Form("n_{#sigma_{TPC}}^{%s}",nameSpec[iS]),"p_{T} (GeV/c)","V0M Multiplicity Percentile","Number of tracklets"};
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) {
	fSparseNsigmaTPC[iS]=new THnSparseF(name_fSparse, title_fSparse, Ndim, bins, xmin, xmax);
	fSparseNsigmaTPC[iS]->GetAxis(3)->Set(Ntrackletsbins, trackletsbins);
      }
      else {
	Int_t Ebins[Ndim] = {1, 1, 1, 1};
	fSparseNsigmaTPC[iS]=new THnSparseF(name_fSparse, title_fSparse, Ndim, Ebins, xmin, xmax);
      }
      for(Int_t j=0;j<Ndim;j++) {
	fSparseNsigmaTPC[iS]->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
      }
    }
  }

  //fNsigmaTPC[1] with only one bin per axis (temporarily)
  Char_t name_fNsigmaTPC[18][200];
  Char_t title_fNsigmaTPC[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fNsigmaTPC[iS],200,"NsigmaTPC_0_%s",nameSpec[iS]);
    snprintf(title_fNsigmaTPC[iS],200,"n#sigma_{TPC} (%s);p_{T} (GeV/c);V0M Multiplicity Percentile;n_{#sigma_{TPC}}^{%s}",nameSpec[iS],nameSpec[iS]);
    if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9) fNsigmaTPC[0][iS] = new TH3F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],Nptbins,ptMin,ptMax,Nmultbins,multV0Min,multV0Max,Nsigmabins,sigmaMin,sigmaMax);
    else fNsigmaTPC[0][iS] = new TH3F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],1,ptMin,ptMax,1,multV0Min,multV0Max,1,sigmaMin,sigmaMax);
  
    snprintf(name_fNsigmaTPC[iS],200,"NsigmaTPC_1_%s",nameSpec[iS]);
    snprintf(title_fNsigmaTPC[iS],200,"n#sigma_{TPC} (%s);p_{T} (GeV/c);Number of tracklets;n_{#sigma_{TPC}}^{%s}",nameSpec[iS],nameSpec[iS]);
    if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9) {
      fNsigmaTPC[1][iS] = new TH3F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],1,ptMin,ptMax,1,0,1000,1,sigmaMin,sigmaMax);//Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nsigmabins,sigmaMin,sigmaMax
      //fNsigmaTPC[1][iS]->GetYaxis()->Set(Ntrackletsbins, trackletsbins);
    }
    else fNsigmaTPC[1][iS] = new TH3F(name_fNsigmaTPC[iS],title_fNsigmaTPC[iS],1,ptMin,ptMax,1,0,1000,1,sigmaMin,sigmaMax);
  }
  
  Char_t name_fNsigmaTPCwTOF[18][200];
  Char_t title_fNsigmaTPCwTOF[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fNsigmaTPCwTOF[iS],200,"NsigmaTPCwTOF_%s",nameSpec[iS]);
    snprintf(title_fNsigmaTPCwTOF[iS],200,"n#sigma_{TPC} (%s) (TOF matching required);p_{T} (GeV/c);n_{#sigma_{TPC}}^{%s}",nameSpec[iS],nameSpec[iS]);
    if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9) fNsigmaTPCwTOF[iS] = new TH2F(name_fNsigmaTPCwTOF[iS],title_fNsigmaTPCwTOF[iS],Nptbins,ptMin,ptMax,Nsigmabins,sigmaMin,sigmaMax);
    else fNsigmaTPCwTOF[iS] = new TH2F(name_fNsigmaTPCwTOF[iS],title_fNsigmaTPCwTOF[iS],1,ptMin,ptMax,1,sigmaMin,sigmaMax);
  }
    
  if(1) {
    const Int_t Ndim=4;
    Int_t bins[Ndim] = {NDCAxybins, 50, Nmultbins, Ntrackletsbins};
    Double_t xmin[Ndim] = {-0.5, 0, multV0Min,    0};
    Double_t xmax[Ndim] = { 0.5, 5, multV0Max, 1000};
    Char_t name_fSparse[200];
    Char_t title_fSparse[200];
    TString axis[Ndim] = {"DCA_{xy} (cm)","p_{T} (GeV/c)","V0M Multiplicity Percentile","Number of tracklets"};
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fSparse,200,"fSparseDcaxy_%s",nameSpec[iS]);
      snprintf(title_fSparse,200,"%s",nameSpec[iS]);
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) {
	fSparseDcaxy[iS]=new THnSparseF(name_fSparse, title_fSparse, Ndim, bins, xmin, xmax);
	fSparseDcaxy[iS]->GetAxis(3)->Set(Ntrackletsbins, trackletsbins);
      }
      else {
	Int_t Ebins[Ndim] = {1, 1, 1, 1};
	fSparseDcaxy[iS]=new THnSparseF(name_fSparse, title_fSparse, Ndim, Ebins, xmin, xmax);
      }
      for(Int_t j=0;j<Ndim;j++) {
	fSparseDcaxy[iS]->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
      }
    }
  }
  
  Char_t name_fDcaxy[18][200];
  Char_t title_fDcaxy[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDcaxy[iS],200,"fDca_xy_0_%s",nameSpec[iS]);
    snprintf(title_fDcaxy[iS],200,"%s;DCA_{xy} (cm);V0M Multiplicity Percentile;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9)
      fDcaxy[0][iS]=new TH3F(name_fDcaxy[iS],title_fDcaxy[iS],NDCAxybins,-0.5,0.5,Nmultbins,multV0Min,multV0Max,50,0,5);
    else fDcaxy[0][iS]=new TH3F(name_fDcaxy[iS],title_fDcaxy[iS],1,-0.5,0.5,1,multV0Min,multV0Max,1,0,5);
    
    snprintf(name_fDcaxy[iS],200,"fDca_xy_1_%s",nameSpec[iS]);
    snprintf(title_fDcaxy[iS],200,"%s;DCA_{xy} (cm);Number of tracklets;p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) {
      fDcaxy[1][iS]=new TH3F(name_fDcaxy[iS],title_fDcaxy[iS],NDCAxybins,-0.5,0.5,Ntrackletsbins,0,1000,50,0,5);
      fDcaxy[1][iS]->GetYaxis()->Set(Ntrackletsbins, trackletsbins);
    }
    else fDcaxy[1][iS]=new TH3F(name_fDcaxy[iS],title_fDcaxy[iS],1,-0.5,0.5,1,0,1000,1,0,5);
  }
  
  //fDcawTOF with only one bin per axis (temporarily)
  Char_t name_fDca[18][200];
  Char_t title_fDca[18][200];
  for(Int_t iM=0;iM<7;iM++) {
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fDca[iS],200,"fDca_xy_wTOF_0_multMin=%.0f_multMax=%.0f_%s",multMin[iM],multMax[iM],nameSpec[iS]);
      snprintf(title_fDca[iS],200,"%s (wTOF) V0M %.0f-%.0f%%;DCA_{xy} (cm);m^{2}_{TOF} (GeV^{2}/c^{4});p_{T} (GeV/c)",
	       nameSpec[iS],multMin[iM],multMax[iM]);
      if(iS==5 || iS==5+9) fDcawTOF_0[iM][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-0.5,0.5,1,0,1,50,0,5);//NDCAxybins,-0.5,0.5,100,0,10,50,0,5
      else fDcawTOF_0[iM][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-0.5,0.5,1,ptMin,ptMax,1,0,5);
    }
  }
  for(Int_t iM=0;iM<14;iM++) {
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fDca[iS],200,"fDca_xy_wTOF_1_multMin=%.0f_multMax=%.0f_%s",trackletsMin[iM],trackletsMax[iM],nameSpec[iS]);
      snprintf(title_fDca[iS],200,"%s (wTOF) Number of tracklets [%.0f-%.0f[;DCA_{xy} (cm);m^{2}_{TOF} (GeV^{2}/c^{4});p_{T} (GeV/c)",
	       nameSpec[iS],trackletsMin[iM],trackletsMax[iM]);
      if(iS==5 || iS==5+9) fDcawTOF_1[iM][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-0.5,0.5,1,0,10,1,0,5);//NDCAxybins,-0.5,0.5,100,0,10,50,0,5
      else fDcawTOF_1[iM][iS]=new TH3F(name_fDca[iS],title_fDca[iS],1,-0.5,0.5,1,ptMin,ptMax,1,0,5);
    }
  }
  
  Char_t name_fDcaz[18][200];
  Char_t title_fDcaz[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fDcaz[iS],200,"fDca_z_%s",nameSpec[iS]);
    snprintf(title_fDcaz[iS],200,"%s;DCA_{z} (cm);p_{T} (GeV/c)",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fDcaz[iS]=new TH2F(name_fDcaz[iS],title_fDcaz[iS],NDCAzbins,-DCAzMax,DCAzMax,50,0,5);
    else fDcaz[iS]=new TH2F(name_fDcaz[iS],title_fDcaz[iS],1,-DCAzMax,DCAzMax,1,0,5);
  }

  Char_t name_fNsigmaTOF[18][200];
  Char_t title_fNsigmaTOF[18][200];    
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fNsigmaTOF[iS],200,"NsigmaTOF_%s",nameSpec[iS]);
    snprintf(title_fNsigmaTOF[iS],200,"n#sigma_{TOF} (%s);p_{T} (GeV/c);n_{#sigma_{TOF}}^{%s}",nameSpec[iS],nameSpec[iS]);
    if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9) fNsigmaTOF[iS] = new TH2F(name_fNsigmaTOF[iS],title_fNsigmaTOF[iS],Nptbins,ptMin,ptMax,Nsigmabins,sigmaMin,sigmaMax);
    else fNsigmaTOF[iS] = new TH2F(name_fNsigmaTOF[iS],title_fNsigmaTOF[iS],1,ptMin,ptMax,1,sigmaMin,sigmaMax);
  }

  Char_t name_fNsigmas[18][200];
  Char_t title_fNsigmas[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fNsigmas[iS],200,"NsigmaTPCvsTOF_%s",nameSpec[iS]);
    snprintf(title_fNsigmas[iS],200,"%s;p_{T} (GeV/c);n_{#sigma_{TOF}}^{%s};n_{#sigma_{TPC}}^{%s}",nameSpec[iS],nameSpec[iS],nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fNsigmas[iS] = new TH3F(name_fNsigmas[iS],title_fNsigmas[iS],Nptbins,ptMin,ptMax,100,-10,10,100,-10,10);
    else fNsigmas[iS] = new TH3F(name_fNsigmas[iS],title_fNsigmas[iS],1,ptMin,ptMax,1,-10,10,1,-10,10);
  }

  fBetaTOFvspt[0] = new TH2F("fBetaTOFvspt_pos","#beta_{TOF} (positive charges);p_{T}/|z| (GeV/c);#beta_{TOF}",200,0,10,260,0.4,1.05);//500,520//200,260
  fBetaTOFvspt[1] = new TH2F("fBetaTOFvspt_neg","#beta_{TOF} (negative charges);p_{T}/|z| (GeV/c);#beta_{TOF}",200,0,10,260,0.4,1.05);

  Char_t name_hBetaExp[9][200];
  Char_t title_hBetaExp[9][200];
  for(Int_t iS=0;iS<9;iS++) {
    snprintf(name_hBetaExp[iS],200,"hBetaTOFVSpt_Exp_%s",nameSpec[iS]);
    snprintf(title_hBetaExp[iS],200,"Expected #beta_{TOF} (%s);p_{T}/|z| (GeV/c); #beta_{TOF}",nameSpec[iS]);
    hBetaExp[iS] = new TProfile(name_hBetaExp[iS],title_hBetaExp[iS],200,0,10,0.4,1.05,"");
  }

  fM2tof[0][0] = new TH3F("fM2tof_pos_0","m^{2}_{TOF} (positive charges);p_{T} (GeV/c);V0M Multiplicity Percentile;m^{2}_{TOF} (GeV^{2}/c^{4})",Nptbins,ptMin,ptMax,Nmultbins,multV0Min,multV0Max,500,0,10);
  fM2tof[0][1] = new TH3F("fM2tof_neg_0","m^{2}_{TOF} (negative charges);p_{T} (GeV/c);V0M Multiplicity Percentile;m^{2}_{TOF} (GeV^{2}/c^{4})",Nptbins,ptMin,ptMax,Nmultbins,multV0Min,multV0Max,500,0,10);
  
  fM2tof[1][0] = new TH3F("fM2tof_pos_1","m^{2}_{TOF} (positive charges);p_{T} (GeV/c);Number of tracklets;m^{2}_{TOF} (GeV^{2}/c^{4})",Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,500,0,10);
  fM2tof[1][1] = new TH3F("fM2tof_neg_1","m^{2}_{TOF} (negative charges);p_{T} (GeV/c);Number of tracklets;m^{2}_{TOF} (GeV^{2}/c^{4})",Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,500,0,10);
  for(Int_t i=0;i<2;i++) {
    fM2tof[1][i]->GetYaxis()->Set(Ntrackletsbins, trackletsbins);
  }

  if(1) {
    const Int_t Ndim=5;
    Int_t bins[Ndim] = {500, 50, NDCAxybins, Nmultbins, Ntrackletsbins};
    Double_t xmin[Ndim] = { 0, 0, -0.5, multV0Min,    0};
    Double_t xmax[Ndim] = {10, 5,  0.5, multV0Max, 1000};
    Char_t name_fSparse[200];
    Char_t title_fSparse[200];
    TString axis[Ndim] = {"m^{2}_{TOF} (GeV^{2}/c^{4})", "p_{T} (GeV/c)", "DCA_{xy} (cm)", "V0M Multiplicity Percentile","Number of tracklets"};
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fSparse,200,"fSparseM2vspt_%s",nameSpec[iS]);
      snprintf(title_fSparse,200,"%s",nameSpec[iS]);
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) {
	fSparseM2vspt[iS]=new THnSparseF(name_fSparse, title_fSparse, Ndim, bins, xmin, xmax);
	fSparseM2vspt[iS]->GetAxis(4)->Set(Ntrackletsbins, trackletsbins);
      }
      else {
	Int_t Ebins[Ndim] = {1, 1, 1, 1, 1};
	fSparseM2vspt[iS]=new THnSparseF(name_fSparse, title_fSparse, Ndim, Ebins, xmin, xmax);
      }
      for(Int_t j=0;j<Ndim;j++) {
	fSparseM2vspt[iS]->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
      }
    }
  }

  //fM2vspt with only one bins per axis (temporarily)
  Char_t name_fM2vspt[18][200];
  Char_t title_fM2vspt[18][200];
  for(Int_t iS=0;iS<18;iS++) {
    snprintf(name_fM2vspt[iS],200,"fM2vspt_0_%s",nameSpec[iS]);
    snprintf(title_fM2vspt[iS],200,"m^{2}_{TOF} (3#sigma TPC dE/dx cut on %s);p_{T} (GeV/c);V0M Multiplicity Percentile;m^{2}_{TOF} (GeV^{2}/c^{4})",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fM2vspt[0][iS] = new TH3F(name_fM2vspt[iS],title_fM2vspt[iS],1,ptMin,ptMax,1,multV0Min,multV0Max,1,0,10);//Nptbins,ptMin,ptMax,Nmultbins,multV0Min,multV0Max,500,0,10)
    else fM2vspt[0][iS] = new TH3F(name_fM2vspt[iS],title_fM2vspt[iS],1,ptMin,ptMax,1,multV0Min,multV0Max,1,0,10);
        
    snprintf(name_fM2vspt[iS],200,"fM2vspt_1_%s",nameSpec[iS]);
    snprintf(title_fM2vspt[iS],200,"m^{2}_{TOF} (3#sigma TPC dE/dx cut on %s);p_{T} (GeV/c);Number of tracklets;m^{2}_{TOF} (GeV^{2}/c^{4})",nameSpec[iS]);
    if(iS==4 || iS==4+9 || iS==5 || iS==5+9) {
      fM2vspt[1][iS] = new TH3F(name_fM2vspt[iS],title_fM2vspt[iS],1,ptMin,ptMax,1,0,1,500,0,10);//Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,500,0,10
      //fM2vspt[1][iS]->GetYaxis()->Set(Ntrackletsbins, trackletsbins);
    }
    else fM2vspt[1][iS] = new TH3F(name_fM2vspt[iS],title_fM2vspt[iS],1,ptMin,ptMax,1,0,1000,1,0,10);
  }

  htemp[0]=new TH1F("htemp0","Ntracklets for identified deuterons with TPC (pt<5)", Ntrackletsbins, trackletsbins);
  
  //Only for MC:
  if(isMC) {
    
    hpdg[0] = new TH2F("hpdg","Pdg label of generated particles (after the event selection);Pdg label;isPrimary           isSecMat           isSecWeak   ",50082,-25041,25041,3,0,3);
    hpdg[0]->GetYaxis()->SetNdivisions(105);

    //bes=before event selection
    const Int_t Ntrackletsbins_bes=14;
    const Float_t trackletsbins_bes[Ntrackletsbins_bes+1]={-1, 0, 1, 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 1000};

    const Int_t Nmultbins_bes=22;
    const Float_t multbins_bes[Nmultbins_bes+1]={-1, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 101};

    if(1) {
      const Int_t Ndim=4;
      Int_t bins[Ndim] = {Nptbins, Ntrackletsbins_bes, Nmultbins_bes, 200};
      Double_t xmin[Ndim] = {ptMin,   -1,  -1, -20};
      Double_t xmax[Ndim] = {ptMax, 1000, 101,  20};
      Char_t name_fSparse[200];
      Char_t title_fSparse[200];
      TString axis[Ndim] = {"p_{T} (GeV/c)", "Number of tracklets", "V0M Multiplicity Percentile", "z_{vtx} (cm)"};
      for(Int_t iS=0;iS<18;iS++) {
	snprintf(name_fSparse,200,"fSparsehpt_%s",nameSpec[iS]);
	snprintf(title_fSparse,200,"p_{T} distribution of generated primary %s (at least one charged particle in |#eta|<1)",nameSpec[iS]);
	if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9) {
	  fSparsehpt[iS]=new THnSparseF(name_fSparse, title_fSparse, Ndim, bins, xmin, xmax);
	  fSparsehpt[iS]->GetAxis(1)->Set(Ntrackletsbins_bes, trackletsbins_bes);
	  fSparsehpt[iS]->GetAxis(2)->Set(Nmultbins_bes, multbins_bes);
	}
	else {
	  Int_t Ebins[Ndim] = {1, 1, 1, 1};
	  fSparsehpt[iS]=new THnSparseF(name_fSparse, title_fSparse, Ndim, Ebins, xmin, xmax);
	}
	for(Int_t j=0;j<Ndim;j++) {
	  fSparsehpt[iS]->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
	}
      }
    }

    hMCTrackletsVsV0mult[0] = new TH3F("hMCTrackletsVsV0mult","Event counter for events with at least one charged particle in |#eta|<1;z_{vtx} (cm);V0M Multiplicity Percentile;Number of tracklets",200,-20,20,Nmultbins_bes,-1,101,Ntrackletsbins_bes,-1,1000);
    hMCTrackletsVsV0mult[0]->GetYaxis()->Set(Nmultbins_bes, multbins_bes);
    hMCTrackletsVsV0mult[0]->GetZaxis()->Set(Ntrackletsbins_bes, trackletsbins_bes);

    htest[0] = new TH1F("htest","htest;pT",Nptbins,ptMin,ptMax);

    //variable binning on Ntracklets set below
    Char_t name_hpt[200];
    Char_t title_hpt[200];
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_gen_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of generated primary %s for |y|<0.5.;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[0][0][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[0][0][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_gen_secMat_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of generated secondary %s from mat. for |y|<0.5.;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[0][1][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[0][1][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_gen_secWeak_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of generated secondary %s from w.d. for |y|<0.5.;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[0][2][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[0][2][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_reco_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of reco. primary %s;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[1][0][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[1][0][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_reco_secMat_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from mat.;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[1][1][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[1][1][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_reco_secWeak_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from w.d.;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[1][2][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[1][2][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_matchedTof_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of reco. primary %s matched at TOF;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[2][0][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[2][0][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_matchedTof_secMat_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from mat. matched at TOF;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[2][1][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[2][1][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_matchedTof_secWeak_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from w.d. matched at TOF;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[2][2][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[2][2][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_goodmatchTof_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of reco. primary %s matched correctly at TOF;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[3][0][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[3][0][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_goodmatchTof_secMat_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from mat. matched correctly at TOF;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[3][1][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[3][1][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_hpt,200,"hpt_goodmatchTof_secWeak_%s",nameSpec[iS]);
      snprintf(title_hpt,200,"p_{T} distribution of reco. secondary %s from w.d. matched correctly at TOF;p_{T} (GeV/c);Number of tracklets;V0M Multiplicity Percentile",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	hpt[3][2][iS] = new TH3F(name_hpt,title_hpt,Nptbins,ptMin,ptMax,Ntrackletsbins,0,1000,Nmultbins,multV0Min,multV0Max);
      else hpt[3][2][iS] = new TH3F(name_hpt,title_hpt,1,ptMin,ptMax,1,0,1000,1,multV0Min,multV0Max);
    }
    
    for(Int_t i=0;i<4;i++) {
      for(Int_t j=0;j<3;j++) {
	for(Int_t iS=0;iS<18;iS++) {
	  if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	    hpt[i][j][iS]->GetYaxis()->Set(Ntrackletsbins, trackletsbins);
	}
      }
    }
    
    Char_t name_fptRecoVsTrue[200];
    Char_t title_fptRecoVsTrue[200];
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fptRecoVsTrue,200,"fptRecoVsTrue_%s",nameSpec[iS]);
      snprintf(title_fptRecoVsTrue,200,"%s;p_{T}^{reco.}/z (GeV/c);p_{T}^{reco.} - p_{T}^{true} (GeV/c)",nameSpec[iS]);
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fptRecoVsTrue[0][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,100,0,5,300,-0.6,0.3);
      else fptRecoVsTrue[0][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,1,0,5,1,-0.6,0.3);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fptRecoVsTrue,200,"fptRecoVsTrue_2_%s",nameSpec[iS]);
      snprintf(title_fptRecoVsTrue,200,"%s, after pT correction;p_{T}^{reco.}/z (GeV/c);p_{T}^{reco.} - p_{T}^{true} (GeV/c)",nameSpec[iS]);
      if(iS==5 || iS==5+9) fptRecoVsTrue[1][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,100,0,5,300,-0.6,0.3);
      else fptRecoVsTrue[1][iS] = new TH2F(name_fptRecoVsTrue,title_fptRecoVsTrue,1,0,5,1,-0.6,0.3);
    }
    
    Char_t name_fmcDca[18][200];
    Char_t title_fmcDca[18][200];
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fmcDca[iS],200,"fmcDca_xy_prim_%s",nameSpec[iS]);
      snprintf(title_fmcDca[iS],200,"primary %s;DCA_{xy};p_{T} (GeV/c)",nameSpec[iS]);
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fmcDca[0][0][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],NDCAxybins,-0.5,0.5,50,0,5);
      else fmcDca[0][0][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],1,-0.5,0.5,1,0,5);
      
      snprintf(name_fmcDca[iS],200,"fmcDca_z_prim_%s",nameSpec[iS]);
      snprintf(title_fmcDca[iS],200,"primary %s;DCA_{z};p_{T} (GeV/c)",nameSpec[iS]);
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fmcDca[0][1][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],NDCAzbins,-DCAzMax,DCAzMax,50,0,5);
      else fmcDca[0][1][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],1,-DCAzMax,DCAzMax,1,0,5);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fmcDca[iS],200,"fmcDca_xy_secMat_%s",nameSpec[iS]);
      snprintf(title_fmcDca[iS],200,"secondary %s from mat.;DCA_{xy};p_{T} (GeV/c)",nameSpec[iS]);
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fmcDca[1][0][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],NDCAxybins,-0.5,0.5,50,0,5);
      else fmcDca[1][0][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],1,-0.5,0.5,1,0,5);
      
      snprintf(name_fmcDca[iS],200,"fmcDca_z_secMat_%s",nameSpec[iS]);
      snprintf(title_fmcDca[iS],200,"secondary %s from mat.;DCA_{z};p_{T} (GeV/c)",nameSpec[iS]);
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fmcDca[1][1][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],NDCAzbins,-DCAzMax,DCAzMax,50,0,5);
      else fmcDca[1][1][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],1,-DCAzMax,DCAzMax,1,0,5);
    }
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fmcDca[iS],200,"fmcDca_xy_secWeak_%s",nameSpec[iS]);
      snprintf(title_fmcDca[iS],200,"secondary %s from w.d.;DCA_{xy};p_{T} (GeV/c)",nameSpec[iS]);
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fmcDca[2][0][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],NDCAxybins,-0.5,0.5,50,0,5);
      else fmcDca[2][0][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],1,-0.5,0.5,1,0,5);
      
      snprintf(name_fmcDca[iS],200,"fmcDca_z_secWeak_%s",nameSpec[iS]);
      snprintf(title_fmcDca[iS],200,"secondary %s from w.d.;DCA_{z};p_{T} (GeV/c)",nameSpec[iS]);
      if(iS==4 || iS==4+9 || iS==5 || iS==5+9) fmcDca[2][1][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],NDCAzbins,-DCAzMax,DCAzMax,50,0,5);
      else fmcDca[2][1][iS]=new TH2F(name_fmcDca[iS],title_fmcDca[iS],1,-DCAzMax,DCAzMax,1,0,5);
    }
    
    Char_t name_fmcNsigmaTOF[18][200];
    Char_t title_fmcNsigmaTOF[18][200];
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fmcNsigmaTOF[iS],200,"fmcNsigmaTOF_prim_equalLabel_%s",nameSpec[iS]);
      snprintf(title_fmcNsigmaTOF[iS],200,"primary %s, TOFlabel == label;p_{T} (GeV/c);n_{#sigma}^{TOF}",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	fmcNsigmaTOF[0][0][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],Nptbins,ptMin,ptMax,Nsigmabins,sigmaMin,sigmaMax);
      else fmcNsigmaTOF[0][0][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],1,ptMin,ptMax,1,sigmaMin,sigmaMax);
      
      snprintf(name_fmcNsigmaTOF[iS],200,"fmcNsigmaTOF_secMat_equalLabel_%s",nameSpec[iS]);
      snprintf(title_fmcNsigmaTOF[iS],200,"secondary %s from mat., TOFlabel == label;p_{T} (GeV/c);n_{#sigma}^{TOF}",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	fmcNsigmaTOF[1][0][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],Nptbins,ptMin,ptMax,Nsigmabins,sigmaMin,sigmaMax);
      else fmcNsigmaTOF[1][0][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],1,ptMin,ptMax,1,sigmaMin,sigmaMax);
      
      snprintf(name_fmcNsigmaTOF[iS],200,"fmcNsigmaTOF_secWeak_equalLabel_%s",nameSpec[iS]);
      snprintf(title_fmcNsigmaTOF[iS],200,"secondary %s from w.d., TOFlabel == label;p_{T} (GeV/c);n_{#sigma}^{TOF}",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	fmcNsigmaTOF[2][0][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],Nptbins,ptMin,ptMax,Nsigmabins,sigmaMin,sigmaMax);
      else fmcNsigmaTOF[2][0][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],1,ptMin,ptMax,1,sigmaMin,sigmaMax);
    }
    
    for(Int_t iS=0;iS<18;iS++) {
      snprintf(name_fmcNsigmaTOF[iS],200,"fmcNsigmaTOF_prim_differentLabel_%s",nameSpec[iS]);
      snprintf(title_fmcNsigmaTOF[iS],200,"primary %s, TOFlabel != label;p_{T} (GeV/c);n_{#sigma}^{TOF}",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	fmcNsigmaTOF[0][1][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],Nptbins,ptMin,ptMax,Nsigmabins,sigmaMin,sigmaMax);
      else fmcNsigmaTOF[0][1][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],1,ptMin,ptMax,1,sigmaMin,sigmaMax);
      
      snprintf(name_fmcNsigmaTOF[iS],200,"fmcNsigmaTOF_secMat_differentLabel_%s",nameSpec[iS]);
      snprintf(title_fmcNsigmaTOF[iS],200,"secondary %s from mat., TOFlabel != label;p_{T} (GeV/c);n_{#sigma}^{TOF}",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	fmcNsigmaTOF[1][1][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],Nptbins,ptMin,ptMax,Nsigmabins,sigmaMin,sigmaMax);
      else fmcNsigmaTOF[1][1][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],1,ptMin,ptMax,1,sigmaMin,sigmaMax);
      
      snprintf(name_fmcNsigmaTOF[iS],200,"fmcNsigmaTOF_secWeak_differentLabel_%s",nameSpec[iS]);
      snprintf(title_fmcNsigmaTOF[iS],200,"secondary %s from w.d., TOFlabel != label;p_{T} (GeV/c);n_{#sigma}^{TOF}",nameSpec[iS]);
      if(iS==2 || iS==2+9 || iS==3 || iS==3+9 || iS==4 || iS==4+9 || iS==5 || iS==5+9)
	fmcNsigmaTOF[2][1][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],Nptbins,ptMin,ptMax,Nsigmabins,sigmaMin,sigmaMax);
      else fmcNsigmaTOF[2][1][iS] = new TH2F(name_fmcNsigmaTOF[iS],title_fmcNsigmaTOF[iS],1,ptMin,ptMax,1,sigmaMin,sigmaMax);
    }
    
  }//Only for MC (end)

  fList->Add(htriggerMask[0]);
  fList->Add(htriggerMask[1]);
  fList->Add(hzvertex);
  fList->Add(hNevent);
  
  fList->Add(hV0mult);
  fList->Add(hTracklets);
  fList->Add(hTrackletsVsV0mult);
  
  for(Int_t i=0;i<2;i++) fList->Add(hrapidity[i]);
  
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
  
  for(Int_t i=0;i<2;i++) fList->Add(fSparseNsigmaTPC[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fSparseNsigmaTPC[5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[0][2+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[0][3+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[1][2+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[1][3+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[1][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPC[1][5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fSparseDcaxy[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fSparseDcaxy[5+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fDcaxy[0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDcaxy[0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDcaxy[1][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDcaxy[1][5+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fDcaz[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fDcaz[5+9*i]);
  
  for(Int_t i=0;i<7;i++) {
    fList->Add(fDcawTOF_0[i][5]);
    fList->Add(fDcawTOF_0[i][5+9]);
  }
  for(Int_t i=0;i<14;i++) {
    fList->Add(fDcawTOF_1[i][5]);
    fList->Add(fDcawTOF_1[i][5+9]);
  }

  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTOF[2+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTOF[3+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTOF[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTOF[5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fNsigmas[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmas[5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fBetaTOFvspt[i]);
  for(Int_t i=0;i<9;i++) fList->Add(hBetaExp[i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2tof[0][i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2tof[1][i]);

  for(Int_t i=0;i<2;i++) fList->Add(fSparseM2vspt[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fSparseM2vspt[5+9*i]);

  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[0][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[0][5+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[1][4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fM2vspt[1][5+9*i]);
  
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPCwTOF[2+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPCwTOF[3+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPCwTOF[4+9*i]);
  for(Int_t i=0;i<2;i++) fList->Add(fNsigmaTPCwTOF[5+9*i]);
  
  fList->Add(htemp[0]);
  
  //Only for MC:
  if(isMC) {
    fList->Add(hpdg[0]);
    
    fList->Add(hMCTrackletsVsV0mult[0]);

    fList->Add(htest[0]);

    for(Int_t i=0;i<2;i++) fList->Add(fSparsehpt[2+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fSparsehpt[3+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fSparsehpt[4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fSparsehpt[5+9*i]);

    for(Int_t j=0;j<3;j++) {
      for(Int_t i=0;i<2;i++) fList->Add(hpt[0][j][2+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[0][j][3+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[0][j][4+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[0][j][5+9*i]);
      
      for(Int_t i=0;i<2;i++) fList->Add(hpt[1][j][2+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[1][j][3+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[1][j][4+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[1][j][5+9*i]);
      
      for(Int_t i=0;i<2;i++) fList->Add(hpt[2][j][2+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[2][j][3+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[2][j][4+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[2][j][5+9*i]);
      
      for(Int_t i=0;i<2;i++) fList->Add(hpt[3][j][2+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[3][j][3+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[3][j][4+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(hpt[3][j][5+9*i]);
    }
  
    for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[0][4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[0][5+9*i]);
    
    for(Int_t i=0;i<2;i++) fList->Add(fptRecoVsTrue[1][5+9*i]);
    
    //DCAxy:
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[0][0][4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[1][0][4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[2][0][4+9*i]);
    
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[0][0][5+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[1][0][5+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[2][0][5+9*i]);
    
    //DCAz:
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[0][1][4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[1][1][4+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[2][1][4+9*i]);
    
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[0][1][5+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[1][1][5+9*i]);
    for(Int_t i=0;i<2;i++) fList->Add(fmcDca[2][1][5+9*i]);
    
    for(Int_t j=0;j<3;j++) {
      for(Int_t i=0;i<2;i++) fList->Add(fmcNsigmaTOF[j][0][2+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(fmcNsigmaTOF[j][1][2+9*i]);
      
      for(Int_t i=0;i<2;i++) fList->Add(fmcNsigmaTOF[j][0][3+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(fmcNsigmaTOF[j][1][3+9*i]);
      
      for(Int_t i=0;i<2;i++) fList->Add(fmcNsigmaTOF[j][0][4+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(fmcNsigmaTOF[j][1][4+9*i]);
      
      for(Int_t i=0;i<2;i++) fList->Add(fmcNsigmaTOF[j][0][5+9*i]);
      for(Int_t i=0;i<2;i++) fList->Add(fmcNsigmaTOF[j][1][5+9*i]);
    }
    
  }//Only for MC (end):
  
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

  if(isMC) {
    eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!eventHandler) return;
    mcEvent = eventHandler->MCEvent();
    fStack = mcEvent->Stack();
  }
  
  //Trigger mask filled:
  for(Int_t i=0;i<32;i++) {
    unsigned bit=(1<<i);//shift of 1 of i positions
    if(inputHandler->IsEventSelected() & bit) htriggerMask[0]->Fill(i);
  }
  //for no MB events:
  if(!(inputHandler->IsEventSelected() & AliVEvent::kMB)) {
    for(Int_t i=0;i<32;i++) {
      unsigned bit=(1<<i);
      if(inputHandler->IsEventSelected() & bit) htriggerMask[1]->Fill(i);
    }
  }
  
  // Fill of MC particles pT spectrum for the true INEL>0:
  // at least one primary charged particle in |eta|<1 and |zvtx_mc|<10 cm
  // See AliAnalysisTaskPPVsMultCrossCheckMC.cxx
  if(isMC) {
       
    Bool_t isTrueINEL=kFALSE;

    for(Int_t iM=0;iM<mcEvent->GetNumberOfTracks();iM++){
      
      AliVParticle *mcpart = (AliVParticle *) mcEvent->GetTrack(iM);
      Int_t label = mcpart->GetLabel();
      Bool_t isPrimary = fStack->IsPhysicalPrimary(label);
      Double_t eta =  mcpart->Eta();
      Double_t charge = ((Double_t) mcpart->Charge())/3;
      
      if(TMath::Abs(charge)<0.001) continue;
      if(!isPrimary) continue;
      if(TMath::Abs(eta)>1) continue;
      
      isTrueINEL=kTRUE;  
    }
    
    //z_vtx (MC)
    Float_t zvtx_mc=mcEvent->GetPrimaryVertex()->GetZ();
    
    Float_t Ntracklets_bes = (Float_t) fPPVsMultUtils->GetStandardReferenceMultiplicity(fEvent, kFALSE);
    if(Ntracklets_bes<0) Ntracklets_bes=-1;

    Float_t mult_bes = fPPVsMultUtils->GetMultiplicityPercentile(fEvent, "V0M", kFALSE);
    if(mult_bes<0) mult_bes=-1;
    if(mult_bes>100) mult_bes=100;
    
    for(Int_t iM=0;iM<mcEvent->GetNumberOfTracks();iM++){
      
      if(!isTrueINEL) break;

      AliVParticle *mcpart = (AliVParticle *) mcEvent->GetTrack(iM);
      
      Int_t kSpec = this->GetMCSpec(mcpart);
      Int_t label = mcpart->GetLabel();
      Bool_t isPrimary = fStack->IsPhysicalPrimary(label);
      Double_t rapidity = mcpart->Y();
      Double_t pt = mcpart->Pt();
      
      if(TMath::Abs(rapidity)>0.5) continue;
      
      if(kSpec<0) continue;
      
      if(!isPrimary) continue;
      
      Double_t value[4]={pt, Ntracklets_bes, mult_bes, zvtx_mc};
      fSparsehpt[kSpec]->Fill(value);
      
      Int_t Pdg=this->GetPdgCode(mcpart);
      if(TMath::Abs(zvtx_mc)<10 && TMath::Abs(Pdg)==211) htest[0]->Fill(pt);
      
    }//loop on MC particles closed
      
    if(isTrueINEL) hMCTrackletsVsV0mult[0]->Fill(zvtx_mc, mult_bes, Ntracklets_bes);

  }//stop MC
  
  //------------------------- Event selection:
  this->EventSelectionMonitor();
  
  if(!fPPVsMultUtils->IsEventSelected(fEvent, AliVEvent::kMB)) return;
  
  Float_t mult = fPPVsMultUtils->GetMultiplicityPercentile(fEvent, "V0M", kFALSE);//kFALSE because I made the event selection before
  
  if(!this->IsInsideFullMultRange(mult)) return;
  hV0mult->Fill(mult);
  //------------------------- Event selection (end)

  Float_t Ntracklets = (Float_t) fPPVsMultUtils->GetStandardReferenceMultiplicity(fEvent, kFALSE);//kFALSE because I made the event selection before
  hTracklets->Fill(Ntracklets);
  hTrackletsVsV0mult->Fill(mult,Ntracklets);
  
  //It is used only to fill fDCAwTOF
  Int_t iMultBin=this->GetMultiplicityBin(mult);
  Int_t iTrackletsBin=this->GetTrackletsBin(Ntracklets);
  
  //------------------------------- Loop on MC particles
  if(isMC) {

    for(Int_t iM=0;iM<mcEvent->GetNumberOfTracks();iM++){
      
      AliVParticle *mcpart = (AliVParticle *) mcEvent->GetTrack(iM);
      
      Int_t Pdg = this->GetPdgCode(mcpart);
      Int_t kSpec = this->GetMCSpec(mcpart);
      Int_t label = mcpart->GetLabel();
      Bool_t isPrimary = fStack->IsPhysicalPrimary(label);
      Bool_t isSecMat = fStack->IsSecondaryFromMaterial(label);
      Bool_t isSecWeak = fStack->IsSecondaryFromWeakDecay(label);
      Double_t rapidity = mcpart->Y();
      Double_t pt = mcpart->Pt();
      
      if(isPrimary) hpdg[0]->Fill(Pdg,0);
      else if(isSecMat) hpdg[0]->Fill(Pdg,1);
      else if(isSecWeak) hpdg[0]->Fill(Pdg,2);
      
      if(TMath::Abs(rapidity)>0.5) continue;
      
      if(kSpec<0) continue;
      
      if(isPrimary) {
	hpt[0][0][kSpec]->Fill(pt,Ntracklets,mult);
      }
      else if(isSecMat) {
	hpt[0][1][kSpec]->Fill(pt,Ntracklets,mult); 
      }
      else if(isSecWeak) {
	hpt[0][2][kSpec]->Fill(pt,Ntracklets,mult); 
      }
    }
  }//------------------------------- Loop on MC particles (end)

  Int_t nTracks = fEvent->GetNumberOfTracks();
  
  //-------------------------------- Loop on reconstructed TRACKS
  for(Int_t iT=0;iT<nTracks;iT++) { 

    AliVTrack* track = (AliVTrack *) fEvent->GetTrack(iT);
    if (!track){
      continue;
    }
    
    //---------- info on MC true:
    AliVParticle *mcpart;
    Int_t label=-1, Pdg=-1, kSpec=-1, t_label=-1;
    Bool_t isPrimary=kFALSE, isSecMat=kFALSE, isSecWeak=kFALSE;
    Double_t t_pt=-1;
    if(isMC) {
      
      label = TMath::Abs(track->GetLabel());
      
      mcpart = (AliVParticle *) mcEvent->GetTrack(label);
      
      Pdg = this->GetPdgCode(mcpart);
      kSpec = this->GetMCSpec(mcpart);
      t_label = mcpart->GetLabel();
      isPrimary = fStack->IsPhysicalPrimary(t_label);
      isSecMat = fStack->IsSecondaryFromMaterial(t_label);
      isSecWeak = fStack->IsSecondaryFromWeakDecay(t_label);
      //Double_t t_rapidity = mcpart->Y();
      t_pt = mcpart->Pt();
            
    }      
    //---------- info on MC true (end)
    
    Double_t rapidity;
    if(isMC) rapidity = mcpart->Y();
    else rapidity = this->GetRapidity(track);
    
    //------------------------- Track cuts:
   
    //rapidity cut:
    hrapidity[0]->Fill(rapidity);
    if(TMath::Abs(rapidity)>0.5) continue;
    hrapidity[1]->Fill(rapidity);

    Double_t DCAxy, DCAz;
    if(!this->AcceptTrack(track, DCAxy, DCAz)) continue;
        
    //------------------------- Track cuts (end)
        
    //reconstructed MC particles:
    if(isMC) {
      this->ForPtCorr(track->Pt(), t_pt, kSpec);
      
      if(kSpec>-1){
	if(isPrimary) {
	  hpt[1][0][kSpec]->Fill(t_pt,Ntracklets,mult);
	  
	  fmcDca[0][0][kSpec]->Fill(DCAxy,t_pt);
	  fmcDca[0][1][kSpec]->Fill(DCAz,t_pt);
	}
	else if(isSecMat) {
	  hpt[1][1][kSpec]->Fill(t_pt,Ntracklets,mult);

	  fmcDca[1][0][kSpec]->Fill(DCAxy,t_pt);
	  fmcDca[1][1][kSpec]->Fill(DCAz,t_pt);
	}
	else if(isSecWeak) {
	  hpt[1][2][kSpec]->Fill(t_pt,Ntracklets,mult);

	  fmcDca[2][0][kSpec]->Fill(DCAxy,t_pt);
	  fmcDca[2][1][kSpec]->Fill(DCAz,t_pt);
	}
      }
    }// reconstructed MC particles (end)
    
    Short_t charge = track->Charge();
    Double_t pt = track->Pt();

    Double_t nsigmaTPC[9];
    this->GetNsigmaTPC(track, nsigmaTPC);

    //pT correction applied:
    this->PtCorr(pt, nsigmaTPC);
    
    for(Int_t iS=0;iS<9;iS++){
      if(charge>0) {
	fNsigmaTPC[0][iS]->Fill(pt,mult,nsigmaTPC[iS]);
	fNsigmaTPC[1][iS]->Fill(pt,Ntracklets,nsigmaTPC[iS]);
	Double_t value[4]={nsigmaTPC[iS], pt, mult, Ntracklets};
	fSparseNsigmaTPC[iS]->Fill(value);
      }
      else if(charge<0) {
	fNsigmaTPC[0][iS+9]->Fill(pt,mult,nsigmaTPC[iS]);
	fNsigmaTPC[1][iS+9]->Fill(pt,Ntracklets,nsigmaTPC[iS]);
	Double_t value[4]={nsigmaTPC[iS], pt, mult, Ntracklets};
	fSparseNsigmaTPC[iS+9]->Fill(value);
      }
    }

    //DCA filled:
    for(Int_t iS=0;iS<9;iS++){
      if(TMath::Abs(nsigmaTPC[iS])<nsigmaTPCMax) {
	if(charge>0) {
	  fDcaxy[0][iS]->Fill(DCAxy,mult,pt);
	  fDcaxy[1][iS]->Fill(DCAxy,Ntracklets,pt);
	  
	  Double_t value[4]={DCAxy, pt, mult, Ntracklets};
	  fSparseDcaxy[iS]->Fill(value);
	  
	  fDcaz[iS]->Fill(DCAz,pt);

	  if(iS==5 && pt<5) htemp[0]->Fill(Ntracklets);
	}
	else if(charge<0) {
	  fDcaxy[0][iS+9]->Fill(DCAxy,mult,pt);
	  fDcaxy[1][iS+9]->Fill(DCAxy,Ntracklets,pt);

	  Double_t value[4]={DCAxy, pt, mult, Ntracklets};
	  fSparseDcaxy[iS+9]->Fill(value);
	  
	  fDcaz[iS+9]->Fill(DCAz,pt);
	}
      }
    }
    
    //TOF matching required:
    if(!this->IsTOFmatching(track)) continue;
    
    for(Int_t iS=0;iS<9;iS++){
      if(charge>0) fNsigmaTPCwTOF[iS]->Fill(pt,nsigmaTPC[iS]);
      else if(charge<0) fNsigmaTPCwTOF[iS+9]->Fill(pt,nsigmaTPC[iS]);
    }
    
    Double_t nsigmaTOF[9];
    Double_t beta = this->GetBeta(track, pt, nsigmaTOF);
    
    for(Int_t iS=0;iS<9;iS++){
      if(charge>0) fNsigmas[iS]->Fill(pt,nsigmaTOF[iS],nsigmaTPC[iS]);
      else if(charge<0) fNsigmas[iS+9]->Fill(pt,nsigmaTOF[iS],nsigmaTPC[iS]);
    }
    
    Double_t p = track->P();
    
    //mass determination:
    Double_t m2 = this->GetM2(p, beta);
    
    if(charge>0) {
      fM2tof[0][0]->Fill(pt,mult,m2);
      fM2tof[1][0]->Fill(pt,Ntracklets,m2);
    }
    else if(charge<0) {
      fM2tof[0][1]->Fill(pt,mult,m2);
      fM2tof[1][1]->Fill(pt,Ntracklets,m2);
    }
    
    for(Int_t iS=0;iS<9;iS++){
      if(TMath::Abs(nsigmaTPC[iS])<nsigmaTPCMax) {
	if(charge>0) {
	  fM2vspt[0][iS]->Fill(pt,mult,m2);
	  fM2vspt[1][iS]->Fill(pt,Ntracklets,m2);

	  Double_t value[5]={m2, pt, DCAxy, mult, Ntracklets};
	  fSparseM2vspt[iS]->Fill(value);

	  fDcawTOF_0[0][iS]->Fill(DCAxy,m2,pt);
	  if(iMultBin>0) fDcawTOF_0[iMultBin][iS]->Fill(DCAxy,m2,pt);
	  
	  fDcawTOF_1[0][iS]->Fill(DCAxy,m2,pt);
	  if(iTrackletsBin>0) fDcawTOF_1[iTrackletsBin][iS]->Fill(DCAxy,m2,pt);
	}
	else if(charge<0) {
	  fM2vspt[0][iS+9]->Fill(pt,mult,m2);
	  fM2vspt[1][iS+9]->Fill(pt,Ntracklets,m2);
	  
	  Double_t value[5]={m2, pt, DCAxy, mult, Ntracklets};
	  fSparseM2vspt[iS+9]->Fill(value);

	  fDcawTOF_0[0][iS+9]->Fill(DCAxy,m2,pt);
	  if(iMultBin>0) fDcawTOF_0[iMultBin][iS+9]->Fill(DCAxy,m2,pt);

	  fDcawTOF_1[0][iS+9]->Fill(DCAxy,m2,pt);
	  if(iTrackletsBin>0) fDcawTOF_1[iTrackletsBin][iS+9]->Fill(DCAxy,m2,pt);
	}
      }
    }
    
    //good TOF match reconstructed MC particles 
    if(isMC) {
      Bool_t isTOFgoodMatch = this->IsTOFgoodmatching(track, label, nsigmaTOF, kSpec, t_pt, isPrimary, isSecMat, isSecWeak);
      
      if(kSpec>-1) {
	
	if(isPrimary) {
	  hpt[2][0][kSpec]->Fill(t_pt,Ntracklets,mult);
	  //good match at TOF
	  if(isTOFgoodMatch) {
	    hpt[3][0][kSpec]->Fill(t_pt,Ntracklets,mult);
	  }
	}
	else if(isSecMat) {
	  hpt[2][1][kSpec]->Fill(t_pt,Ntracklets,mult);
	  //good match at TOF
	  if(isTOFgoodMatch) {
	    hpt[3][1][kSpec]->Fill(t_pt,Ntracklets,mult);
	  }
	}
	else if(isSecWeak) {
	  hpt[2][2][kSpec]->Fill(t_pt,Ntracklets,mult);
	  //good match at TOF
	  if(isTOFgoodMatch) {
	    hpt[3][2][kSpec]->Fill(t_pt,Ntracklets,mult);
	  }
	}
	
      }
      
    }//good TOF match reconstructed MC particles (end)
    
  }//----------------------loop on reconstructed TRACKS (end)
  
}//end loop on the events
//_____________________________________________________________________________
void AliAnalysisNuclMult::Terminate(Option_t *) { 
  printf("Terminate()\n");
}
//_____________________________________________________________________________
void AliAnalysisNuclMult::EventSelectionMonitor() {
  
  if(fPPVsMultUtils->IsEventSelected(fEvent, AliVEvent::kMB)) {
    hNevent->Fill(5);
  }
  
  if(fPPVsMultUtils->IsSelectedTrigger(fEvent, AliVEvent::kMB)) {
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
	  
	  if(this->IsInsideFullMultRange(fPPVsMultUtils->GetMultiplicityPercentile(fEvent,"V0M"))) {
	    hNevent->Fill(6);
	  }
	  
	}
      }
    }
  }
  
  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisNuclMult::IsInsideFullMultRange(Float_t multiplicity) {

  if(multiplicity < 0 || multiplicity > 100) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________________
Int_t AliAnalysisNuclMult::GetMultiplicityBin(Float_t multiplicity) {

  Int_t iMultBinT=-1;
  for(Int_t iM=1;iM<7;iM++) {
    if(multiplicity > multMin[iM] && multiplicity < multMax[iM]) iMultBinT=iM;
  }
  Int_t iMultBin=iMultBinT;
  
  return iMultBin;
}
//_____________________________________________________________________________
Int_t AliAnalysisNuclMult::GetTrackletsBin(Int_t Ntracklets) {

  Int_t iTrackletsBinT=-1;
  for(Int_t iM=1;iM<14;iM++) {
    if(Ntracklets >= trackletsMin[iM] && Ntracklets < trackletsMax[iM]) iTrackletsBinT=iM;
  }
  Int_t iTrackletsBin=iTrackletsBinT;
  
  return iTrackletsBin;
}
//_____________________________________________________________________________
Double_t AliAnalysisNuclMult::GetRapidity(AliVTrack *track) {

  Double_t p = track->P();
  Double_t pz = track->Pz();
  Double_t massOverZ = AliPID::ParticleMassZ(AliPID::kDeuteron);
  
  Double_t E = TMath::Sqrt(p*p+massOverZ*massOverZ);
  Double_t rapidity;
  
  if((E-pz)>1e-18) rapidity = 0.5*TMath::Log((E+pz)/(E-pz));
  else rapidity = -99;
    
  return rapidity;
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

  if(TMath::Abs(DCAxy)>DCAxyMax) return kFALSE;

  if(TMath::Abs(DCAz)>DCAzMax) return kFALSE;

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
  if(nTPCclusters>=nTPCclustersMin) {
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
		  
		  if(TMath::Abs(DCAxy)<0.5) {
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

  track->GetIntegratedTimes(exptimes);
  //Integrated times of nuclei:
  Double_t m_proton = AliPID::ParticleMass(4);
  for(Int_t iN=5;iN<9;iN++) {
    Double_t massOverZ = AliPID::ParticleMassZ(iN);
    if(p>1e-18) exptimes[iN] = exptimes[4]*exptimes[4]*(massOverZ*massOverZ/p/p+1)/(m_proton*m_proton/p/p+1);
    exptimes[iN] = TMath::Sqrt(exptimes[iN]);
  }  
  
  return;
}
//_____________________________________________________________________________
Double_t AliAnalysisNuclMult::GetBeta(AliVTrack *track, Double_t pt, Double_t nsigmaTOF[9]) {
  
  Short_t charge = track->Charge();
  Double_t p = track->P();
  Double_t tof = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
  Double_t len = track->GetIntegratedLength();

  Double_t exptimes[9];
  this->GetExpTOFtime(track, p, exptimes);
    
  Double_t tofres[9];
  for(Int_t iS=0;iS<9;iS++) tofres[iS] = fPIDResponse->GetTOFResponse().GetExpectedSigma(p, exptimes[iS], (AliPID::EParticleType) iS);

  //Double_t nsigmaTOF[9];
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
Int_t AliAnalysisNuclMult::GetPdgCode(AliVParticle *mcpart) {
  
  Int_t Pdg = mcpart->PdgCode();
  //for nuclei: e.g. deuteron: 1000010020-1000000000=10020
  if(Pdg>1e+9) Pdg -= 1e+9;
  else if (Pdg<-1e+9) Pdg += 1e+9;
  
  return Pdg;
}
//_____________________________________________________________________________
Int_t AliAnalysisNuclMult::GetMCSpec(AliVParticle *mcpart) {
  
  Int_t Pdg = this->GetPdgCode(mcpart);

  Int_t kSpec=-1;
  for(Int_t iS=0;iS<18;iS++) {
    if(Pdg!=PdgStd[iS]) continue;
    kSpec=iS;
  }

  return kSpec;
}
//_____________________________________________________________________________
Double_t AliAnalysisNuclMult::GetM2(Double_t p, Double_t beta) {

  if(beta<1e-18) return -1;
  
  Double_t m2 = (p*p)*(1-beta*beta)/(beta*beta);
  
  return m2;
}
//_____________________________________________________________________________
void AliAnalysisNuclMult::ForPtCorr(Double_t pt, Double_t t_pt, Int_t kSpec) { 
  
  if(kSpec<0) return;
  
  //The true pt is not divided by the charge Z:
  if(kSpec==7 || kSpec==8 || kSpec==7+9 || kSpec==8+9) t_pt=t_pt/2;
  
  fptRecoVsTrue[0][kSpec]->Fill(pt,pt-t_pt);
  
  //To check the pt correction:
  fptRecoVsTrue[1][kSpec]->Fill(pt,(pt-fptCorr[0]->Eval(pt))-t_pt);

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisNuclMult::IsTOFgoodmatching(AliVTrack *track, Int_t label, Double_t nsigmaTOF[9], Int_t kSpec, Double_t t_pt, Bool_t isPrimary, Bool_t isSecMat, Bool_t isSecWeak) {
 
  /*
    return ( TOFlabel == label || TMath::Abs(nsigmaTOF)<3 )
   */

  Int_t *toflabel = new Int_t[3];
  ((AliESDtrack *)track)->GetTOFLabel(toflabel);
  Bool_t isTOFgoodMatch = kFALSE;
  if(toflabel[0]==label) isTOFgoodMatch = kTRUE;

  if(kSpec<0) return isTOFgoodMatch;
  
  Int_t kSpec_2=kSpec;
  if(kSpec_2 > 8) {//e.g. kSpec={5,5+9}->kSpec_2={5}
    kSpec_2-=9;
  }  
 
  if(isPrimary) {
    if(isTOFgoodMatch) fmcNsigmaTOF[0][0][kSpec]->Fill(t_pt,nsigmaTOF[kSpec_2]);
    else fmcNsigmaTOF[0][1][kSpec]->Fill(t_pt,nsigmaTOF[kSpec_2]);
  }
  else if(isSecMat) {
    if(isTOFgoodMatch) fmcNsigmaTOF[1][0][kSpec]->Fill(t_pt,nsigmaTOF[kSpec_2]);
    else fmcNsigmaTOF[1][1][kSpec]->Fill(t_pt,nsigmaTOF[kSpec_2]);
  }
  else if(isSecWeak) {
   if(isTOFgoodMatch) fmcNsigmaTOF[2][0][kSpec]->Fill(t_pt,nsigmaTOF[kSpec_2]);
   else fmcNsigmaTOF[2][1][kSpec]->Fill(t_pt,nsigmaTOF[kSpec_2]);
  }
  
  Bool_t isTOFgoodMatch_2 = isTOFgoodMatch;

  if(TMath::Abs(nsigmaTOF[kSpec_2])<3) isTOFgoodMatch_2 = kTRUE;

  return isTOFgoodMatch_2;
}
