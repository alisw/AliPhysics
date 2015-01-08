#include <iostream>
#include <math.h>
#include "TChain.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjString.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THn.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TObjectTable.h"
#include <vector>

//#include "TStopwatch.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"


#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
//#include "AliAnalysisUtils.h"
#include "AliHelperPID.h"
#include "AliEventPoolManager.h"

#include "AliAnalysisTaskFemtoESE.h"

//#include "AliSpectraAODEventCuts.h"
//#include "AliSpectraAODTrackCuts.h"
//#include "/opt/alice/aliroot/master/src/PWGLF/SPECTRA/PiKaPr/TestAOD/AliSpectraAODEventCuts.h"
//#include "/opt/alice/aliroot/master/src/PWGLF/SPECTRA/PiKaPr/TestAOD/AliSpectraAODTrackCuts.h"
//#include "AliSpectraAODEventCuts.h"
//#include "AliSpectraAODTrackCuts.h"

ClassImp(AliAnalysisTaskFemtoESE)


//________________________________________________________________________
// Default constructor
AliAnalysisTaskFemtoESE::AliAnalysisTaskFemtoESE() : 
AliAnalysisTaskSE(),
    fAOD(0x0), 
    fOutputList(0x0),
    fHelperPID(0x0),
    fPoolMgr(0x0),
    fEventCuts(0x0),
    fTrackCuts(0x0),
    fFilterBit(128),
    fSelectBit(AliVEvent::kMB),
    bIsLHC10h(1),
    fEventCounter(0),
    fMixingTracks(10000),
    fBfield(0.),
    fMinSepPairEta(0.),
    fMinSepPairPhi(0.),
    fQinvMin(-1.),
    fMaxDcaXY(1000.),
    fMaxDcaZ(1000.),
    fPtMin(0.14),
    fPtMax(1.5),
    fEtaMax(0.8),
    fShareQuality(0.5),
    fShareFraction(0.05),
    nCountSamePairs(0),
    nCountMixedPairs(0),
    nCountTracks(0),
    fMinQPerc(-1000),
    fMaxQPerc(1000),
    qlimit(0.2),
    qbins(40),
    fQPercDet(0),
    fEPDet(0),
    nKtBins(0),
    nKtBins1(1),
    ktBins(0),
    nEPBins(0),
    nEPBins1(1),
    epBins(0),
    nEPBinsMix(0),
    nEPBinsMix1(1),
    epBinsMix(0),
    nCentBins(0),
    nCentBins1(1),
    centBins(0),
    nVzBins(0),
    nVzBins1(1),
    vzBins(0),
    hq(0x0),
    hqmix(0x0),
    hqinv(0x0),
    nqPercBinsLHC11h(1), 
    qPercBinsLHC11h(0),
    hpx(0x0),
    hpy(0x0),
    hpz(0x0),
    hpt(0x0),
    hE(0x0),
    hphieta(0x0),
    hphieta_pid(0x0),
    hpt_pid(0x0),
    hvzcent(0x0),
    hcent(0x0),
    hcentUnweighted(0x0),
    hcentn(0x0),
    hphistaretapair10(0x0),
    hphistaretapair16(0x0),
    hphistaretapair10a(0x0),
    hphistaretapair16a(0x0),
    hphistaretapair10b(0x0),
    hphistaretapair16b(0x0),
    hphietapair(0x0),
    hphietapair2(0x0),
    hpidid(0x0),
    hkt(0x0),
    hktcheck(0x0),
    hkt3(0x0),
    hdcaxy(0x0),
    hdcaz(0x0),
    hsharequal(0x0),
    hsharefrac(0x0),
    hsharequalmix(0x0),
    hsharefracmix(0x0),
    hPsiTPC(0x0),
    hPsiV0A(0x0),
    hPsiV0C(0x0),
    hShiftTPC(0x0),
    hShiftV0A(0x0),
    hShiftV0C(0x0),
    hPsiMix(0x0),
    hCheckEPA(0x0),
    hCheckEPC(0x0),
    hCheckEPmix(0x0),
    hAvDphi(0x0),
    hNpairs(0x0),
    hPairDphi(0x0),
    hPairDphiMix(0x0),
    hcentq(0x0),
    hMixedDistTracks(0x0),
    hMixedDistEvents(0x0),
    hQvecV0A(0x0),
    hQvecV0C(0x0),
    hresV0ATPC(0x0),
    hresV0CTPC(0x0),
    hresV0AV0C(0x0),
    hqinvcheck(0x0),
    hktbins(0x0),
    hcentbins(0x0),
    hepbins(0x0)
{
}

//________________________________________________________________________
// Constructor
AliAnalysisTaskFemtoESE::AliAnalysisTaskFemtoESE(const char* name) : 
  AliAnalysisTaskSE(name),
  fAOD(0x0), 
  fOutputList(0x0),
  fHelperPID(0x0),
  fPoolMgr(0x0),
  fEventCuts(0x0),
  fTrackCuts(0x0),
  fFilterBit(128),
  fSelectBit(AliVEvent::kMB),
  bIsLHC10h(1),
  fEventCounter(0),
  fMixingTracks(10000),
  fBfield(0.),
  fMinSepPairEta(0.),
  fMinSepPairPhi(0.),
  fQinvMin(-1.),
  fMaxDcaXY(1000.),
  fMaxDcaZ(1000.),
  fPtMin(0.14),
  fPtMax(1.5),
  fEtaMax(0.8),
  fShareQuality(0.5),
  fShareFraction(0.05),
  nCountSamePairs(0),
  nCountMixedPairs(0),
  nCountTracks(0),
  fMinQPerc(-1000),
  fMaxQPerc(1000),
  qlimit(0.2),
  qbins(40),
  fQPercDet(0),
  fEPDet(0),
  nKtBins(0),
  nKtBins1(1),
  ktBins(0),
  nEPBins(0),
  nEPBins1(1),
  epBins(0),
  nEPBinsMix(0),
  nEPBinsMix1(1),
  epBinsMix(0),
  nCentBins(0),
  nCentBins1(1),
  centBins(0),
  nVzBins(0),
  nVzBins1(1),
  vzBins(0),
  hq(0x0),
  hqmix(0x0),
  hqinv(0x0),
  nqPercBinsLHC11h(1), 
  qPercBinsLHC11h(0),
  hpx(0x0),
  hpy(0x0),
  hpz(0x0),
  hpt(0x0),
  hE(0x0),
  hphieta(0x0),
  hphieta_pid(0x0),
  hpt_pid(0x0),
  hvzcent(0x0),
  hcent(0x0),
  hcentUnweighted(0x0),
  hcentn(0x0),
  hphistaretapair10(0x0),
  hphistaretapair16(0x0),
  hphistaretapair10a(0x0),
  hphistaretapair16a(0x0),
  hphistaretapair10b(0x0),
  hphistaretapair16b(0x0),
  hphietapair(0x0),
  hphietapair2(0x0),
  hpidid(0x0),
  hkt(0x0),
  hktcheck(0x0),
  hkt3(0x0),
  hdcaxy(0x0),
  hdcaz(0x0),
  hsharequal(0x0),
  hsharefrac(0x0),
  hsharequalmix(0x0),
  hsharefracmix(0x0),
  hPsiTPC(0x0),
  hPsiV0A(0x0),
  hPsiV0C(0x0),
  hShiftTPC(0x0),
  hShiftV0A(0x0),
  hShiftV0C(0x0),
  hPsiMix(0x0),
  hCheckEPA(0x0),
  hCheckEPC(0x0),
  hCheckEPmix(0x0),
  hAvDphi(0x0),
  hNpairs(0x0),
  hPairDphi(0x0),
  hPairDphiMix(0x0),
  hcentq(0x0),
  hMixedDistTracks(0x0),
  hMixedDistEvents(0x0),
  hQvecV0A(0x0),
  hQvecV0C(0x0),
  hresV0ATPC(0x0),
  hresV0CTPC(0x0),
  hresV0AV0C(0x0),
  hqinvcheck(0x0),
  hktbins(0x0),
  hcentbins(0x0),
  hepbins(0x0)
{
  
  Printf("*******************************************");
  Printf("AliAnalysisTaskFemtoESE named %s",name);
  Printf("*******************************************");

  // default binning
  //SetEPBins(12,-TMath::Pi()/12.,2*TMath::Pi()-TMath::Pi()/12.);
  SetEPBins(12);
  Double_t ktBinsTemp[5] = {0.2,0.3,0.4,0.5,0.7};
  SetKtBins(4,ktBinsTemp);
  Double_t centBinsTemp[7] = {0,5,10,20,30,40,50};
  SetCentBins(6,centBinsTemp);
  Double_t vzBinsTemp[9] = {-8,-6,-4,-2,0,2,4,6,8};
  SetVzBins(8,vzBinsTemp);

  nEPBinsMix = 6;
  nEPBinsMix1 = 7;
  epBinsMix = new Double_t[nEPBinsMix1];
  for(Int_t ee = 0; ee < nEPBinsMix1; ee++)
    {epBinsMix[ee] = (TMath::Pi()/(Double_t)nEPBinsMix)*((Double_t)ee)-TMath::Pi()/2.;}

  vertex[0] = vertex[1] = vertex[2] = 0.;

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliHelperPID::Class());
  DefineOutput(3, AliSpectraAODEventCuts::Class());
  DefineOutput(4, AliSpectraAODTrackCuts::Class());

}

//________________________________________________________________________
// destructor
AliAnalysisTaskFemtoESE::~AliAnalysisTaskFemtoESE()
{
  /* do nothing yet */  
}
//________________________________________________________________________
// copy constructor
AliAnalysisTaskFemtoESE::AliAnalysisTaskFemtoESE(const AliAnalysisTaskFemtoESE &/*obj*/) :
  AliAnalysisTaskSE(),
  fAOD(0x0), 
  fOutputList(0x0),
  fHelperPID(0x0),
  fPoolMgr(0x0),
  fEventCuts(0x0),
  fTrackCuts(0x0),
  fFilterBit(128),
  fSelectBit(AliVEvent::kMB),
  bIsLHC10h(1),
  fEventCounter(0),
  fMixingTracks(10000),
  fBfield(0.),
  fMinSepPairEta(0.),
  fMinSepPairPhi(0.),
  fQinvMin(-1.),
  fMaxDcaXY(1000.),
  fMaxDcaZ(1000.),
  fPtMin(0.14),
  fPtMax(1.5),
  fEtaMax(0.8),
  fShareQuality(0.5),
  fShareFraction(0.05),
  nCountSamePairs(0),
  nCountMixedPairs(0),
  nCountTracks(0),
  fMinQPerc(-1000),
  fMaxQPerc(1000),
  qlimit(0.2),
  qbins(40),
  fQPercDet(0),
  fEPDet(0),
  nKtBins(0),
  nKtBins1(1),
  ktBins(0),
  nEPBins(0),
  nEPBins1(1),
  epBins(0),
  nEPBinsMix(0),
  nEPBinsMix1(1),
  epBinsMix(0),
  nCentBins(0),
  nCentBins1(1),
  centBins(0),
  nVzBins(0),
  nVzBins1(1),
  vzBins(0),
  hq(0x0),
  hqmix(0x0),
  hqinv(0x0),
  nqPercBinsLHC11h(1), 
  qPercBinsLHC11h(0),
  hpx(0x0),
  hpy(0x0),
  hpz(0x0),
  hpt(0x0),
  hE(0x0),
  hphieta(0x0),
  hphieta_pid(0x0),
  hpt_pid(0x0),
  hvzcent(0x0),
  hcent(0x0),
  hcentUnweighted(0x0),
  hcentn(0x0),
  hphistaretapair10(0x0),
  hphistaretapair16(0x0),
  hphistaretapair10a(0x0),
  hphistaretapair16a(0x0),
  hphistaretapair10b(0x0),
  hphistaretapair16b(0x0),
  hphietapair(0x0),
  hphietapair2(0x0),
  hpidid(0x0),
  hkt(0x0),
  hktcheck(0x0),
  hkt3(0x0),
  hdcaxy(0x0),
  hdcaz(0x0),
  hsharequal(0x0),
  hsharefrac(0x0),
  hsharequalmix(0x0),
  hsharefracmix(0x0),
  hPsiTPC(0x0),
  hPsiV0A(0x0),
  hPsiV0C(0x0),
  hShiftTPC(0x0),
  hShiftV0A(0x0),
  hShiftV0C(0x0),
  hPsiMix(0x0),
  hCheckEPA(0x0),
  hCheckEPC(0x0),
  hCheckEPmix(0x0),
  hAvDphi(0x0),
  hNpairs(0x0),
  hPairDphi(0x0),
  hPairDphiMix(0x0),
  hcentq(0x0),
  hMixedDistTracks(0x0),
  hMixedDistEvents(0x0),
  hQvecV0A(0x0),
  hQvecV0C(0x0),
  hresV0ATPC(0x0),
  hresV0CTPC(0x0),
  hresV0AV0C(0x0),
  hqinvcheck(0x0),
  hktbins(0x0),
  hcentbins(0x0),
  hepbins(0x0)
{
  /* do nothing yet */  
} 
//________________________________________________________________________
// assignment operator
AliAnalysisTaskFemtoESE& AliAnalysisTaskFemtoESE::operator=(const AliAnalysisTaskFemtoESE &/*obj*/)
{
  /* do nothing yet */  
  return *this;
}

//________________________________________________________________________
void AliAnalysisTaskFemtoESE::UserCreateOutputObjects()
{
  // Create histograms

  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");
  if (!fHelperPID)  AliFatal("HelperPID should be set in the steering macro");

  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->SetName("fOutputList");
  
  hpx = new TH1D("hpx","px",200,-2,2);
  hpx->GetXaxis()->SetTitle("p_{x}");
  fOutputList->Add(hpx);
  hpy = new TH1D("hpy","py",200,-2,2);
  hpy->GetXaxis()->SetTitle("p_{y}");
  fOutputList->Add(hpy);
  hpz = new TH1D("hpz","pz",200,-2,2);
  hpz->GetXaxis()->SetTitle("p_{z}");
  fOutputList->Add(hpz);
  hpt = new TH1D("hpt","pt",100,0,2);
  hpt->GetXaxis()->SetTitle("p_{t}");
  fOutputList->Add(hpt);
  hE = new TH1D("hE","E",100,0,2);
  hE->GetXaxis()->SetTitle("E");
  fOutputList->Add(hE);
  hphieta = new TH2D("hphieta","track #varphi vs #eta",100,0,2*TMath::Pi(),80,-0.8,0.8);
  hphieta->GetXaxis()->SetTitle("#varphi");
  hphieta->GetYaxis()->SetTitle("#eta");
  fOutputList->Add(hphieta);
  hphieta_pid = new TH2D("hphieta_pid","PID check -- #Delta#varphi vs #Delta#eta",100,-0.3,0.3,100,-0.3,0.3);
  hphieta_pid->GetXaxis()->SetTitle("#Delta#varphi");
  hphieta_pid->GetYaxis()->SetTitle("#Delta#eta");
  fOutputList->Add(hphieta_pid);
  hpt_pid = new TH1D("hpt_pid","PID check -- #Delta p_{t}",100,-0.5,0.5);
  hpt_pid->GetXaxis()->SetTitle("#Delta p_{t}");
  fOutputList->Add(hpt_pid);
  hvzcent = new TH2D("hvzcent","vz vs cent",20,-10,10,nCentBins,centBins);
  hvzcent->GetXaxis()->SetTitle("v_{z}");
  hvzcent->GetYaxis()->SetTitle("centrality");
  fOutputList->Add(hvzcent);
  hcent = new TH1D("hcent","cent",200,0,50);
  hcent->GetXaxis()->SetTitle("centrality");
  fOutputList->Add(hcent);
  hcentUnweighted = new TH1D("hcentUnweighted","cent - unweighted",200,0,50);
  hcentUnweighted->GetXaxis()->SetTitle("centrality");
  fOutputList->Add(hcentUnweighted);
  hcentn = new TH2D("hcentn","cent vs npions",50,0,50,100,0,2000);
  hcentn->GetXaxis()->SetTitle("Centrality");
  hcentn->GetYaxis()->SetTitle("Number of pions");
  fOutputList->Add(hcentn);
  hphistaretapair10 = new TH3D("hphistaretapair10","pair #Delta#varphi* vs #Delta#eta at r=1.0m",100,-0.15,0.15,100,-0.1,0.1,10,0,1);
  hphistaretapair10->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair10->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair10->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair10);
  hphistaretapair16 = new TH3D("hphistaretapair16","pair #Delta#varphi* vs #Delta#eta at r=1.6m",100,-0.15,0.15,100,-0.1,0.1,10,0,1);
  hphistaretapair16->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair16->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair16->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair16);
  hphistaretapair10a = new TH3D("hphistaretapair10a","pair #Delta#varphi* vs #Delta#eta at r=1.0m",100,-0.15,0.15,100,-0.1,0.1,10,0,1);
  hphistaretapair10a->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair10a->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair10a->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair10a);
  hphistaretapair16a = new TH3D("hphistaretapair16a","pair #Delta#varphi* vs #Delta#eta at r=1.6m",100,-0.15,0.15,100,-0.1,0.1,10,0,1);
  hphistaretapair16a->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair16a->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair16a->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair16a);
  hphistaretapair10b = new TH3D("hphistaretapair10b","pair #Delta#varphi* vs #Delta#eta at r=1.0m",100,-TMath::Pi(),TMath::Pi(),100,-1.6,1.6,10,0,1);
  hphistaretapair10b->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair10b->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair10b->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair10b);
  hphistaretapair16b = new TH3D("hphistaretapair16b","pair #Delta#varphi* vs #Delta#eta at r=1.6m",100,-TMath::Pi(),TMath::Pi(),100,-1.6,1.6,10,0,1);
  hphistaretapair16b->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair16b->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair16b->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair16b);
  hphietapair = new TH3D("hphietapair","pair #Delta#varphi vs #Delta#eta",100,-0.1,0.1,100,-0.1,0.1,10,0,1);
  hphietapair->GetXaxis()->SetTitle("#Delta#varphi");
  hphietapair->GetYaxis()->SetTitle("#Delta#eta");
  hphietapair->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphietapair);
  hphietapair2 = new TH3D("hphietapair2","pair #varphi vs #eta",100,-TMath::Pi(),TMath::Pi(),100,-1.6,1.6,10,0,1);
  hphietapair2->GetXaxis()->SetTitle("#Delta#varphi");
  hphietapair2->GetYaxis()->SetTitle("#eta");
  hphietapair2->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphietapair2);
  hpidid = new TH1D("hpidid","pid id",9,-4.5,4.5);
  hpidid->GetXaxis()->SetTitle("track PID ID");
  fOutputList->Add(hpidid);
  hkt = new TH2D("hkt","k_{T}",nCentBins,centBins,100,0,2);
  hkt->GetXaxis()->SetTitle("centrality");
  hkt->GetYaxis()->SetTitle("k_{T}");
  fOutputList->Add(hkt);
  hktcheck = new TH1D("hktcheck","k_{T} check",50,0,1);
  hktcheck->GetXaxis()->SetTitle("k_{T}");
  fOutputList->Add(hktcheck);
  hkt3 = new TH3D("hkt3","kt vs pt",50,0,1,50,0,5,50,0,5);
  hkt3->GetXaxis()->SetTitle("k_{T}");
  hkt3->GetYaxis()->SetTitle("p_{T,1}");
  hkt3->GetZaxis()->SetTitle("p_{T,2}");
  fOutputList->Add(hkt3);
  hdcaxy = new TH2D("hdcaxy","DCA xy",100,-5,5,100,-5,5);
  hdcaxy->GetXaxis()->SetTitle("DCA x");
  hdcaxy->GetYaxis()->SetTitle("DCA y");
  fOutputList->Add(hdcaxy);
  hdcaz = new TH1D("hdcaz","DCA z",100,-5,5);
  hdcaz->GetXaxis()->SetTitle("DCA z");
  fOutputList->Add(hdcaz);
  hsharequal = new TH1D("hsharequal","Share Quality",102,-1.02,1.02);
  hsharequal->GetXaxis()->SetTitle("Share Quality");
  fOutputList->Add(hsharequal);
  hsharefrac = new TH1D("hsharefrac","Share Fraction",100,0,1);
  hsharefrac->GetXaxis()->SetTitle("Share Fraction");
  fOutputList->Add(hsharefrac);
  hsharequalmix = new TH1D("hsharequalmix","Share Quality -- mixed events",102,-1.02,1.02);
  hsharequalmix->GetXaxis()->SetTitle("Share Quality");
  fOutputList->Add(hsharequalmix);
  hsharefracmix = new TH1D("hsharefracmix","Share Fraction -- mixed events",100,0,1);
  hsharefracmix->GetXaxis()->SetTitle("Share Fraction");
  fOutputList->Add(hsharefracmix);
  hPsiTPC = new TH2D("hPsiTPC","TPC EP",nCentBins,centBins,100,-1*TMath::Pi(),TMath::Pi());
  hPsiTPC->GetXaxis()->SetTitle("Centrality");
  hPsiTPC->GetYaxis()->SetTitle("#Psi{TPC}");
  fOutputList->Add(hPsiTPC);
  hPsiV0A = new TH2D("hPsiV0A","V0A EP",nCentBins,centBins,100,-1*TMath::Pi(),TMath::Pi());
  hPsiV0A->GetXaxis()->SetTitle("Centrality");
  hPsiV0A->GetYaxis()->SetTitle("#Psi{V0A}");
  fOutputList->Add(hPsiV0A);
  hPsiV0C = new TH2D("hPsiV0C","V0C EP",nCentBins,centBins,100,-1*TMath::Pi(),TMath::Pi());
  hPsiV0C->GetXaxis()->SetTitle("Centrality");
  hPsiV0C->GetYaxis()->SetTitle("#Psi{V0C}");
  fOutputList->Add(hPsiV0C);
  hShiftTPC = new TH2D("hShiftTPC","TPC shifting values",nCentBins,centBins,6,-0.5,5.5);
  hShiftTPC->GetXaxis()->SetTitle("Centrality");
  fOutputList->Add(hShiftTPC);
  hShiftV0A = new TH2D("hShiftV0A","V0A shifting values",nCentBins,centBins,6,-0.5,5.5);
  hShiftV0A->GetXaxis()->SetTitle("Centrality");
  fOutputList->Add(hShiftV0A);
  hShiftV0C = new TH2D("hShiftV0C","V0C shifting values",nCentBins,centBins,6,-0.5,5.5);
  hShiftV0C->GetXaxis()->SetTitle("Centrality");
  fOutputList->Add(hShiftV0C);
  hPsiMix = new TH2D("hPsiMix","Mix EP",nCentBins,centBins,100,-1*TMath::Pi(),TMath::Pi());
  hPsiMix->GetXaxis()->SetTitle("Centrality");
  hPsiMix->GetYaxis()->SetTitle("#Psi{Mix}");
  fOutputList->Add(hPsiMix);
  hCheckEPA = new TH2D("hCheckEPA","Check EP V0A",nCentBins,centBins,100,-1*TMath::Pi(),TMath::Pi());
  hCheckEPA->GetXaxis()->SetTitle("Centrality");
  hCheckEPA->GetYaxis()->SetTitle("PsiV0A - PsiTPC");
  fOutputList->Add(hCheckEPA);
  hCheckEPC = new TH2D("hCheckEPC","Check EP V0C",nCentBins,centBins,100,-1*TMath::Pi(),TMath::Pi());
  hCheckEPC->GetXaxis()->SetTitle("Centrality");
  hCheckEPC->GetYaxis()->SetTitle("PsiV0C - PsiTPC");
  fOutputList->Add(hCheckEPC);
  hCheckEPmix = new TH2D("hCheckEPmix","Check EP mixed events",100,-1*TMath::Pi(),TMath::Pi(),100,-1*TMath::Pi(),TMath::Pi());
  hCheckEPmix->GetXaxis()->SetTitle("Psi1 - Psi_mix");
  hCheckEPmix->GetYaxis()->SetTitle("Psi1 - Psi2");
  fOutputList->Add(hCheckEPmix);
  hAvDphi = new TH3D("hAvDphi","average event plane angle",nEPBins,epBins,nCentBins,centBins,nKtBins,ktBins);
  hAvDphi->GetXaxis()->SetTitle("EP bins");
  hAvDphi->GetYaxis()->SetTitle("cent bins");
  hAvDphi->GetZaxis()->SetTitle("kt bins");
  fOutputList->Add(hAvDphi);
  hNpairs = new TH3D("hNpairs","number of pairs",nEPBins,epBins,nCentBins,centBins,nKtBins,ktBins);
  hNpairs->GetXaxis()->SetTitle("EP bins");
  hNpairs->GetYaxis()->SetTitle("cent bins");
  hNpairs->GetZaxis()->SetTitle("kt bins");
  fOutputList->Add(hNpairs);
  hPairDphi = new TH1D("hPairDphi","pairs wrt EP",100,0,2*TMath::Pi());
  hPairDphi->GetXaxis()->SetTitle("#Delta#phi");
  fOutputList->Add(hPairDphi);
  hPairDphiMix = new TH1D("hPairDphiMix","mixed pairs wrt EP",100,0,2*TMath::Pi());
  hPairDphiMix->GetXaxis()->SetTitle("#Delta#phi");
  fOutputList->Add(hPairDphiMix);
  hcentq = new TH2D("hcentq","qvec vs cent",100,0,100,50,0,50);
  hcentq->GetXaxis()->SetTitle("q_{2} percentile");
  hcentq->GetYaxis()->SetTitle("centrality");
  fOutputList->Add(hcentq);
  hMixedDistTracks = new TH2D("hMixedDistTracks", ";centrality;tracks", nCentBins, centBins, 200, 0, fMixingTracks * 2.);
  fOutputList->Add(hMixedDistTracks);
  hMixedDistEvents = new TH2D("hMixedDistEvents", ";centrality;events", nCentBins, centBins, 41, -0.5, 40.5);
  fOutputList->Add(hMixedDistEvents);
  hQvecV0A = new TH2D("hQvecV0A","Qvector in V0A",50,0,50,200,0,5);
  hQvecV0A->GetXaxis()->SetTitle("Centrality");
  hQvecV0A->GetYaxis()->SetTitle("normalized Qvector");
  fOutputList->Add(hQvecV0A);
  hQvecV0C = new TH2D("hQvecV0C","Qvector in V0C",50,0,50,200,0,5);
  hQvecV0C->GetXaxis()->SetTitle("Centrality");
  hQvecV0C->GetYaxis()->SetTitle("normalized Qvector");
  fOutputList->Add(hQvecV0C);

  // resolution histograms
  hresV0ATPC = new TH1D("hresV0ATPC","cent vs cos(2*(V0A-TPC))",nCentBins,centBins);
  hresV0ATPC->GetXaxis()->SetTitle("centrality");
  fOutputList->Add(hresV0ATPC);
  hresV0CTPC = new TH1D("hresV0CTPC","cent vs cos(2*(V0C-TPC))",nCentBins,centBins);
  hresV0CTPC->GetXaxis()->SetTitle("centrality");
  fOutputList->Add(hresV0CTPC);
  hresV0AV0C = new TH1D("hresV0AV0C","cent vs cos(2*(V0A-V0C))",nCentBins,centBins);
  hresV0AV0C->GetXaxis()->SetTitle("centrality");
  fOutputList->Add(hresV0AV0C);

  hqinvcheck = new TH3F("hqinvcheck","Qinv vs kt vs cent",100,0,1,50,0,1,10,0,50);
  hqinvcheck->GetXaxis()->SetTitle("qinv");
  hqinvcheck->GetYaxis()->SetTitle("kt");
  hqinvcheck->GetZaxis()->SetTitle("centrality");
  fOutputList->Add(hqinvcheck);

  Double_t limit1 = 8.71488;
  Double_t limit2 = 12.0126;
  hq = new TH3F****[2];
  hqmix = new TH3F****[2];
  hqinv = new TH3F****[2];
  for(Int_t z = 0; z < 2; z++) // charge combinations
    {
      hq[z] = new TH3F***[nKtBins];
      hqmix[z] = new TH3F***[nKtBins];
      hqinv[z] = new TH3F***[nKtBins];
      for(Int_t k = 0; k < nKtBins; k++)
	{
	  hq[z][k] = new TH3F**[nEPBins];
	  hqmix[z][k] = new TH3F**[nEPBins];
	  hqinv[z][k] = new TH3F**[nEPBins];
	  for(Int_t e = 0; e < nEPBins; e++)
	    {
	      hq[z][k][e] = new TH3F*[nCentBins];
	      hqmix[z][k][e] = new TH3F*[nCentBins];
	      hqinv[z][k][e] = new TH3F*[nCentBins];
	      for(Int_t c = 0; c < nCentBins; c++)
		{
		  //hq[z][k][e][c] = new TH3F(Form("hq%i_k%i_e%i_c%i",z,k,e,c),Form("hq%i_k%i_e%i_c%i",z,k,e,c),60,-0.21,0.21,60,-0.21,0.21,60,-0.21,0.21);
		  hq[z][k][e][c] = new TH3F(Form("hq%i_k%i_e%i_c%i",z,k,e,c),Form("hq%i_k%i_e%i_c%i",z,k,e,c),qbins,-1.*qlimit,qlimit,qbins,-1.*qlimit,qlimit,qbins,-1.*qlimit,qlimit);
		  hq[z][k][e][c]->GetXaxis()->SetTitle("qout"); hq[z][k][e][c]->GetYaxis()->SetTitle("qside"); hq[z][k][e][c]->GetZaxis()->SetTitle("qlong");
		  if(!(centBins[c]>limit2 || centBins[c+1]<limit1))
		    hq[z][k][e][c]->Sumw2();
		  // set sumw2 only for the correlation histograms which are filled with centrality weights (around cent=10)
		  fOutputList->Add(hq[z][k][e][c]);
		  //hqmix[z][k][e][c] = new TH3F(Form("hqmix%i_k%i_e%i_c%i",z,k,e,c),Form("hqmix%i_k%i_e%i_c%i",z,k,e,c),60,-0.21,0.21,60,-0.21,0.21,60,-0.21,0.21);
		  hqmix[z][k][e][c] = new TH3F(Form("hqmix%i_k%i_e%i_c%i",z,k,e,c),Form("hqmix%i_k%i_e%i_c%i",z,k,e,c),qbins,-1.*qlimit,qlimit,qbins,-1.*qlimit,qlimit,qbins,-1.*qlimit,qlimit);
		  if(!(centBins[c]>limit2 || centBins[c+1]<limit1))
		    hqmix[z][k][e][c]->Sumw2();
		  fOutputList->Add(hqmix[z][k][e][c]);
		  //hqinv[z][k][e][c] = new TH3F(Form("hqinv%i_k%i_e%i_c%i",z,k,e,c),Form("hqinv%i_k%i_e%i_c%i",z,k,e,c),60,-0.21,0.21,60,-0.21,0.21,60,-0.21,0.21);
		  hqinv[z][k][e][c] = new TH3F(Form("hqinv%i_k%i_e%i_c%i",z,k,e,c),Form("hqinv%i_k%i_e%i_c%i",z,k,e,c),qbins,-1.*qlimit,qlimit,qbins,-1.*qlimit,qlimit,qbins,-1.*qlimit,qlimit);
		  //hqinv[z][k][e][c]->Sumw2();
		  fOutputList->Add(hqinv[z][k][e][c]);
		}
	    }
	}
    }

  // create dummy histograms which just hold the values of the kt, cent, ep bin edges
  hktbins = new TH1F("hktbins","kt bins",nKtBins,ktBins);
  fOutputList->Add(hktbins);
  hcentbins = new TH1F("hcentbins","cent bins",nCentBins,centBins);
  fOutputList->Add(hcentbins);
  hepbins = new TH1F("hepbins","ep bins",nEPBins,epBins);
  fOutputList->Add(hepbins);

  Printf("************************");
  Printf("using the %s detector for event plane determination!",fEPDet ? "V0C" : "V0A");
  Printf("using the %s detector for q-vector determination!",fQPercDet ? "V0C" : "V0A");
  Printf("************************");

  vertex[0] = vertex[1] = vertex[2] = 0.;

  // event mixing pool
  fPoolMgr = new AliEventPoolManager*[2];
  Int_t poolsize = 1000;
  fPoolMgr[0] = new AliEventPoolManager(poolsize, fMixingTracks, nCentBins, centBins, nVzBins, vzBins,nEPBinsMix,epBinsMix);
  fPoolMgr[0]->SetTargetValues(fMixingTracks, 0.01, 5); // check these values
  fPoolMgr[1] = new AliEventPoolManager(poolsize, fMixingTracks, nCentBins, centBins, nVzBins, vzBins,nEPBinsMix,epBinsMix);
  fPoolMgr[1]->SetTargetValues(fMixingTracks, 0.01, 5); // check these values

  nCountSamePairs = 0;
  nCountMixedPairs = 0;
  nCountTracks = 0;

  //stopwatch = new TStopwatch();
  //stopwatch->Start();

  PostData(1, fOutputList);
  PostData(2, fHelperPID);
  PostData(3, fEventCuts);
  PostData(4, fTrackCuts);

}

//________________________________________________________________________
void AliAnalysisTaskFemtoESE::UserExec(Option_t *) 
{ 
  // Main loop
  // Called for each event

  //if(!fAODcase) {cout<<"ESDs not supported"<<endl; return;}

  fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!fAOD) {Printf("ERROR: fAOD not available"); return;}

  //if(!EventCut(fAOD)) return;  
  if(!EventCut()) return; 

  fEventCounter++;
  if(fEventCounter%1000==0) Printf("===========  Event # %i  ===========",fEventCounter);

  AliCentrality *centrality;// for AODs and ESDs
  const AliAODVertex *primaryVertexAOD;
 
  //AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  //AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  //fPIDResponse = inputHandler->GetPIDResponse();

  fBfield = fAOD->GetMagneticField();

  /////////////////////////////////////////////////
  
  Float_t centralityPercentile=0;
  
  centrality = fAOD->GetCentrality();
  centralityPercentile = centrality->GetCentralityPercentile("V0M");
  if(centralityPercentile <= centBins[0]) return;
  if(centralityPercentile > centBins[nCentBins]) return;
  Double_t centWeight = GetCentralityWeight(centralityPercentile);
  
  primaryVertexAOD = fAOD->GetPrimaryVertex();
  vertex[0]=primaryVertexAOD->GetX(); vertex[1]=primaryVertexAOD->GetY(); vertex[2]=primaryVertexAOD->GetZ();
  Double_t zvtx = vertex[2];
  if(zvtx < vzBins[0] || zvtx > vzBins[nVzBins]) return; // Z-Vertex Cut 
  //cout<<"Centrality % = " << centralityPercentile << "  z-vertex = " << zvtx << endl;

  //stopwatch->Stop();
  //Printf("%lf   %lf",stopwatch->RealTime(),stopwatch->CpuTime());
  //stopwatch->Start();

  // get event plane from V0's
  if(!fEventCuts->IsSelected(fAOD,fTrackCuts)) {Printf("Error! Event not accepted by AliAODSpectraEventCuts!"); return;}
  Double_t psiV0A = fEventCuts->GetPsiV0A();
  Double_t psiV0C = fEventCuts->GetPsiV0C();
  Double_t qperc = -999;
  qperc = fEventCuts->GetQvecPercentile(fQPercDet);//0: VZERO-A 1: VZERO-C // now works for both LHC10h and LHC11h (for 0-50% centrality only!)
  if(!bIsLHC10h)
    {
      //Printf("%lf   %lf",qperc,GetQPercLHC11h(fEventCuts->GetqV0A()));
      hQvecV0A->Fill(centralityPercentile,fEventCuts->GetqV0A(),centWeight);
      hQvecV0C->Fill(centralityPercentile,fEventCuts->GetqV0C(),centWeight);
      //Printf("q vector = %lf  percentile = %lf",(fQPercDet==0 ? fEventCuts->GetqV0A() : fEventCuts->GetqV0C()),qperc);
    }

  if(psiV0A == -999) return;
  if(psiV0C == -999) return;
  if(qperc < fMinQPerc || qperc > fMaxQPerc) return;

  //ProcInfo_t procInfo;
  //gSystem->GetProcInfo(&procInfo);
  //Printf("beginning of event: ResMem %ld VMem %ld", procInfo.fMemResident, procInfo.fMemVirtual);

  Double_t psiEP = psiV0A;
  if(fEPDet==1) psiEP = psiV0C;

  Double_t sin2phi = 0, cos2phi = 0;

  TObjArray* tracksPos = new TObjArray();
  tracksPos->SetOwner(kTRUE);
  TObjArray* tracksNeg = new TObjArray();
  tracksNeg->SetOwner(kTRUE);

  // Track loop -- select pions
  for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
    AliAODTrack* aodtrack = (AliAODTrack*)fAOD->GetTrack(i);
    if (!aodtrack) continue;
    if(!TrackCut(aodtrack)) continue;

    // check event plane angle using tracks in the TPC
    if(aodtrack->Pt() < 2 && aodtrack->Pt() > 0.2)
      {
	sin2phi += (aodtrack->Pt())*sin(2*aodtrack->Phi());
	cos2phi += (aodtrack->Pt())*cos(2*aodtrack->Phi());
      }

    // filter bit 7 PID method...
    Int_t trackPID=999;
    for(Int_t m = 0; m < fAOD->GetNumberOfTracks(); m++) {
      AliAODTrack* aodtrack2 = (AliAODTrack*)fAOD->GetTrack(m);
      if (!aodtrack2) continue;
      if(aodtrack->GetID() != (-aodtrack2->GetID() - 1)) continue;// (-aodTrack2->GetID() - 1)
      trackPID=fHelperPID->GetParticleSpecies((AliVTrack*)aodtrack2,kTRUE);
      hphieta_pid->Fill(aodtrack->Phi()-aodtrack2->Phi(),aodtrack->Eta()-aodtrack2->Eta());
      hpt_pid->Fill(aodtrack->Pt()-aodtrack2->Pt());
      //cout << aodtrack->Phi() << "   " << aodtrack->Eta() << "   " << aodtrack->Pt() << endl;
      //cout << aodtrack2->Phi() << "   " << aodtrack2->Eta() << "   " << aodtrack2->Pt() << "   " << dataID2 << endl;
    }

    hpidid->Fill((trackPID+1)*aodtrack->Charge());

    // select pions
    if(trackPID==0)
      {
	AliFemtoESEBasicParticle* particle = new AliFemtoESEBasicParticle(sqrt(pow(aodtrack->P(),2)+pow(0.13957, 2)),aodtrack->Px(),aodtrack->Py(),aodtrack->Pz(),aodtrack->Charge(),aodtrack->Phi(),aodtrack->Eta());
	particle->SetPsiEP(psiEP);
	particle->SetTPCClusterMap(aodtrack->GetTPCClusterMap());
	particle->SetTPCSharedMap(aodtrack->GetTPCSharedMap());

	if(particle->Charge()>0)
	  tracksPos->Add(particle);
	if(particle->Charge()<0)
	  tracksNeg->Add(particle);

	// track qa plots
	hpx->Fill(particle->Px());
	hpy->Fill(particle->Py());
	hpz->Fill(particle->Pz());
	hpt->Fill(particle->Pt());
	hE->Fill(particle->E());
	hphieta->Fill(particle->Phi(),particle->Eta());
      }

  }
  // end track loop

  Int_t ntracks = tracksPos->GetEntriesFast()+tracksNeg->GetEntriesFast();

  // get EP from TPC, just to check
  Double_t psiTPC = 0.;
  if(sin2phi != 0 && cos2phi != 0)
    psiTPC = 0.5*atan2(sin2phi,cos2phi);
  else return;

  hPsiTPC->Fill(centralityPercentile,psiTPC,centWeight);
  hPsiV0A->Fill(centralityPercentile,psiV0A,centWeight);
  hPsiV0C->Fill(centralityPercentile,psiV0C,centWeight);
  Double_t dphiEP = psiTPC-psiV0A;
  if(dphiEP>TMath::Pi()) dphiEP-=2*TMath::Pi();
  if(dphiEP<-TMath::Pi()) dphiEP+=2*TMath::Pi();
  hCheckEPA->Fill(centralityPercentile,dphiEP,centWeight);
  dphiEP = psiTPC-psiV0C;
  if(dphiEP>TMath::Pi()) dphiEP-=2*TMath::Pi();
  if(dphiEP<-TMath::Pi()) dphiEP+=2*TMath::Pi();
  hCheckEPC->Fill(centralityPercentile,dphiEP,centWeight);

  // values for EP shifting method  
  hShiftTPC->Fill(centralityPercentile,0.,cos(2*psiTPC)*centWeight);
  hShiftTPC->Fill(centralityPercentile,1.,sin(2*psiTPC)*centWeight);
  hShiftTPC->Fill(centralityPercentile,2.,cos(4*psiTPC)*centWeight);
  hShiftTPC->Fill(centralityPercentile,3.,sin(4*psiTPC)*centWeight);
  hShiftTPC->Fill(centralityPercentile,4.,cos(6*psiTPC)*centWeight);
  hShiftTPC->Fill(centralityPercentile,5.,sin(6*psiTPC)*centWeight);

  hShiftV0A->Fill(centralityPercentile,0.,cos(2*psiV0A)*centWeight);
  hShiftV0A->Fill(centralityPercentile,1.,sin(2*psiV0A)*centWeight);
  hShiftV0A->Fill(centralityPercentile,2.,cos(4*psiV0A)*centWeight);
  hShiftV0A->Fill(centralityPercentile,3.,sin(4*psiV0A)*centWeight);
  hShiftV0A->Fill(centralityPercentile,4.,cos(6*psiV0A)*centWeight);
  hShiftV0A->Fill(centralityPercentile,5.,sin(6*psiV0A)*centWeight);

  hShiftV0C->Fill(centralityPercentile,0.,cos(2*psiV0C)*centWeight);
  hShiftV0C->Fill(centralityPercentile,1.,sin(2*psiV0C)*centWeight);
  hShiftV0C->Fill(centralityPercentile,2.,cos(4*psiV0C)*centWeight);
  hShiftV0C->Fill(centralityPercentile,3.,sin(4*psiV0C)*centWeight);
  hShiftV0C->Fill(centralityPercentile,4.,cos(6*psiV0C)*centWeight);
  hShiftV0C->Fill(centralityPercentile,5.,sin(6*psiV0C)*centWeight);


  hcentq->Fill(qperc,centralityPercentile,centWeight);

  nCountTracks += ntracks;
  //cout << "Found " << ntracks << " pion tracks..." << endl;

  hvzcent->Fill(zvtx,centralityPercentile,centWeight);
  hcentUnweighted->Fill(centralityPercentile);
  hcent->Fill(centralityPercentile,centWeight);
  hcentn->Fill(centralityPercentile,ntracks,centWeight);

  // resolution histograms
  hresV0ATPC->Fill(centralityPercentile,cos(2*(psiV0A-psiTPC))*centWeight);
  hresV0CTPC->Fill(centralityPercentile,cos(2*(psiV0C-psiTPC))*centWeight);
  hresV0AV0C->Fill(centralityPercentile,cos(2*(psiV0A-psiV0C))*centWeight);

  AliEventPool* poolPos = fPoolMgr[0]->GetEventPool(centralityPercentile,zvtx,psiEP);
  if (!poolPos) AliFatal(Form("No pool found for centrality = %f, vz = %f, ep = %f", centralityPercentile, zvtx,psiEP));
  //Printf("Positive pool found for centrality = %f, vz = %f, ep = %f with %i events in it", centralityPercentile, zvtx,psiEP,poolPos->GetCurrentNEvents());
  AliEventPool* poolNeg = fPoolMgr[1]->GetEventPool(centralityPercentile,zvtx,psiEP);
  if (!poolNeg) AliFatal(Form("No pool found for centrality = %f, vz = %f, ep = %f", centralityPercentile, zvtx,psiEP));
  //Printf("Negative pool found for centrality = %f, vz = %f, ep = %f with %i events in it", centralityPercentile, zvtx,psiEP,poolNeg->GetCurrentNEvents());

  TrackLoop(tracksPos,poolPos,0,psiEP,centralityPercentile);
  TrackLoop(tracksNeg,poolNeg,1,psiEP,centralityPercentile);
  //TObjArray* clonedtracks = CloneAndReduceTrackList(tracks,psiEP);
  poolPos->UpdatePool(tracksPos);
  poolNeg->UpdatePool(tracksNeg);
  //cout << "pool contains " << pool->GetCurrentNEvents() << " events and " << pool->NTracksInPool() << " tracks." << endl;
  //tracks->Clear();
  
  //delete tracks;

  //gSystem->GetProcInfo(&procInfo);
  //Printf("end of event: ResMem %ld VMem %ld", procInfo.fMemResident, procInfo.fMemVirtual);

  // Post output data.
  PostData(1, fOutputList);
  PostData(2, fHelperPID);
  PostData(3, fEventCuts);
  PostData(4, fTrackCuts);
}


void AliAnalysisTaskFemtoESE::TrackLoop(TObjArray *tracks, AliEventPool *pool, Int_t z, Double_t psiEP, Float_t centralityPercentile)
{
  // z=0 for positive charges, z=1 for negative charges
  Double_t kt = 0;
  Double_t qout=0, qside=0, qlong=0;
  Double_t pVect1[4] = {0,0,0,0};
  Double_t pVect2[4] = {0,0,0,0};
  Int_t k, e, c; //bin indices for histograms

  Double_t centWeight = GetCentralityWeight(centralityPercentile);

  Int_t ntracks = tracks->GetEntriesFast();

  Int_t countMixedEvents = 0, countMixedTracks = 0;

  // same event loop
  for(Int_t j = 0; j < ntracks; j++)
    {
      //cout << endl << j << "   ";
      AliFemtoESEBasicParticle* track1 = (AliFemtoESEBasicParticle*)tracks->At(j);
      pVect1[0]=track1->E();
      pVect1[1]=track1->Px();
      pVect1[2]=track1->Py();
      pVect1[3]=track1->Pz();
      //cout << pVect1[0] << "   " << pVect1[1] << "   " <<  pVect1[2] << "   " << pVect1[3] << endl;

      for(Int_t i = j+1; i < ntracks; i++)
	{
	  AliFemtoESEBasicParticle* track2 = (AliFemtoESEBasicParticle*)tracks->At(i);

	  Double_t deltaphistar10 = DeltaPhiStar(track1,track2,1.0);
	  Double_t deltaphistar16 = DeltaPhiStar(track1,track2,1.6);

	  hphistaretapair10->Fill(deltaphistar10,track1->Eta()-track2->Eta(),kt);
	  hphistaretapair16->Fill(deltaphistar16,track1->Eta()-track2->Eta(),kt);

	  if(!PairCut(track1,track2,kFALSE)) continue;

	  hphistaretapair10a->Fill(deltaphistar10,track1->Eta()-track2->Eta(),kt);
	  hphistaretapair16a->Fill(deltaphistar16,track1->Eta()-track2->Eta(),kt);
	  hphistaretapair10b->Fill(deltaphistar10,track1->Eta()-track2->Eta(),kt);
	  hphistaretapair16b->Fill(deltaphistar16,track1->Eta()-track2->Eta(),kt);

	  pVect2[0]=track2->E();
	  pVect2[1]=track2->Px();
	  pVect2[2]=track2->Py();
	  pVect2[3]=track2->Pz();

	  //qinv = GetQinv(pVect1, pVect2); // = qinv**2 = (P1x-P2x)**2 + (P1y-P2y)**2 + (P1z-P2z)**2 - (P1t-P2t)**2 
	  GetQosl(pVect1, pVect2, qout, qside, qlong); // qout, qside, qlong = components of Q=P1-P2 in the P=P1+P2 frame
	  if(!(fabs(qout)<qlimit && fabs(qside)<qlimit && fabs(qlong)<qlimit)) continue;
	  nCountSamePairs++;
	  kt = sqrt(pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.; // = Kt = |pT1+pT2|/2
	  hkt->Fill(centralityPercentile,kt,centWeight);
	  hkt3->Fill(kt,track1->Pt(),track2->Pt());
	  Double_t deltaphi = GetDeltaPhiEP(pVect1[1],pVect1[2],pVect2[1],pVect2[2],psiEP); // angle to event plane in correct range
	  if(!FindBin(kt,deltaphi,centralityPercentile,k,e,c)) continue;
	  hktcheck->Fill(kt);
	  if(deltaphi > TMath::Pi())
	    {
	      hq[z][k][e][c]->Fill(qout,-1.*qside,-1.*qlong,centWeight);
	      hqinv[z][k][e][c]->Fill(qout,-1.*qside,-1.*qlong,GetQinv(pVect1, pVect2)*centWeight);
	    }
	  else
	    {
	      hq[z][k][e][c]->Fill(qout,qside,qlong,centWeight);
	      hqinv[z][k][e][c]->Fill(qout,qside,qlong,GetQinv(pVect1, pVect2)*centWeight);
	    }
	  hNpairs->Fill(deltaphi,centralityPercentile,kt);
	  hAvDphi->Fill(deltaphi,centralityPercentile,kt,deltaphi);
	  hPairDphi->Fill(deltaphi);
	  hqinvcheck->Fill(GetQinv(pVect1, pVect2),kt,centralityPercentile,centWeight);
	  Double_t dphi = track1->Phi()-track2->Phi();
	  if(dphi<-TMath::Pi()) dphi += 2*TMath::Pi();
	  if(dphi>TMath::Pi()) dphi -= 2*TMath::Pi();
	  hphietapair->Fill(dphi,track1->Eta()-track2->Eta(),kt);
	  hphietapair2->Fill(dphi,track1->Eta()-track2->Eta(),kt);
	  //cout << k << "  ";
	}
    }

  // mixed event loop
  countMixedTracks = 0;
  countMixedEvents = 0;
  if (pool->IsReady()) 
    {
      Int_t nMix = pool->GetCurrentNEvents();
      
      for (Int_t jMix=0; jMix<nMix; jMix++) // loop over events in pool
	{
	  TObjArray* bgTracks = pool->GetEvent(jMix);
	  Int_t ntracksmix = bgTracks->GetEntriesFast();
	  if(ntracksmix <= 0) continue;
	  
	  // compare event planes of two events
	  AliFemtoESEBasicParticle* tracktest = (AliFemtoESEBasicParticle*)bgTracks->UncheckedAt(0);
	  Double_t psiEPmixtemp = tracktest->GetPsiEP();
	  /*Double_t dphiEPtest = fPsiEPmixtemp-psiEP;
	  while(dphiEPtest>2*TMath::Pi()) dphiEPtest-=2*TMath::Pi();
	  while(dphiEPtest<0) dphiEPtest+=2*TMath::Pi();
	  if(dphiEPtest>TMath::Pi()) dphiEPtest-=TMath::Pi();
	  if(dphiEPtest>TMath::Pi()/2.) dphiEPtest = TMath::Pi()-dphiEPtest;
	  if(dphiEPtest > TMath::Pi()/6.) continue; // event planes must be within pi/6
	  */
	  
	  Double_t psiEPmix = 0.5*atan2(sin(2*psiEP)+sin(2*psiEPmixtemp),cos(2*psiEP)+cos(2*psiEPmixtemp));
	  hPsiMix->Fill(centralityPercentile,psiEPmix,centWeight);

	  Double_t dphimix = psiEP-psiEPmix;
	  if(dphimix < -TMath::Pi()) dphimix += 2*TMath::Pi();
	  if(dphimix > TMath::Pi()) dphimix -= 2*TMath::Pi();
	  Double_t dphi12 = psiEP-psiEPmixtemp;
	  if(dphi12 < -TMath::Pi()) dphi12 += 2*TMath::Pi();
	  if(dphi12 > TMath::Pi()) dphi12 -= 2*TMath::Pi();
	  hCheckEPmix->Fill(dphimix,dphi12);
	    
	  countMixedTracks += ntracksmix;
	  countMixedEvents += 1;

	  for(Int_t j = 0; j < ntracks; j++)
	    {
	      //cout << endl << j << "   ";
	      AliFemtoESEBasicParticle* track1 = (AliFemtoESEBasicParticle*)tracks->At(j);
	      pVect1[0]=track1->E();
	      pVect1[1]=track1->Px();
	      pVect1[2]=track1->Py();
	      pVect1[3]=track1->Pz();
	      
	      for(Int_t i = 0; i < ntracksmix; i++)
		{
		  AliFemtoESEBasicParticle* track2 = (AliFemtoESEBasicParticle*)bgTracks->UncheckedAt(i);

		  if(!PairCut(track1,track2,kTRUE)) continue;

		  pVect2[0]=track2->E();
		  pVect2[1]=track2->Px();
		  pVect2[2]=track2->Py();
		  pVect2[3]=track2->Pz();

		  if(psiEPmixtemp != track2->GetPsiEP()) AliFatal("Error! Event plane angles are wrong in mixing!!");

		  //qinv = GetQinv(pVect1, pVect2); // qinv**2 = (P1x-P2x)**2 + (P1y-P2y)**2 + (P1z-P2z)**2 - (P1t-P2t)**2 
		  GetQosl(pVect1, pVect2, qout, qside, qlong); // qout, qside, qlong = components of Q=P1-P2 in the P=P1+P2 frame
		  if(!(fabs(qout)<qlimit && fabs(qside)<qlimit && fabs(qlong)<qlimit)) continue;
		  nCountMixedPairs++;
		  kt = sqrt(pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.; // = Kt = |pT1+pT2|/2
		  Double_t deltaphi = GetDeltaPhiEP(pVect1[1],pVect1[2],pVect2[1],pVect2[2],psiEPmix); // angle to event plane in correct range
		  //Double_t weight = 1./(Double_t)nMix;
		  if(!FindBin(kt,deltaphi,centralityPercentile,k,e,c)) continue;
		  if(deltaphi > TMath::Pi())
		    hqmix[z][k][e][c]->Fill(qout,-1.*qside,-1.*qlong,centWeight);
		  else
		    hqmix[z][k][e][c]->Fill(qout,qside,qlong,centWeight);

		  hPairDphiMix->Fill(deltaphi);

		}
	    }
	}
    }

  hMixedDistTracks->Fill(centralityPercentile, countMixedTracks);
  hMixedDistEvents->Fill(centralityPercentile, countMixedEvents);
  //Printf("mixed tracks: %i   mixed events: %i",countMixedTracks,countMixedEvents);
}





//________________________________________________________________________
void AliAnalysisTaskFemtoESE::Terminate(Option_t *) 
{

  if(ktBins) delete [] ktBins;
  if(epBins) delete [] epBins;
  if(centBins) delete [] centBins;
  if(vzBins) delete [] vzBins;
  if(qPercBinsLHC11h) delete [] qPercBinsLHC11h;

  // Called once at the end of the query
 
  Printf("Done");

}


/*

//________________________________________________________________________
Bool_t AliAnalysisTaskFemtoESE::AcceptPair(AliChaoticityTrackStruct *first, AliChaoticityTrackStruct *second)
{
  
if(fabs(first->fEta-second->fEta) > fMinSepPairEta) return kTRUE;
  
// propagate through B field to r=1m
Float_t phi1 = first->fPhi - asin(first->fCharge*(0.1*fBfield)*0.15/first->fPt);// 0.15 for D=1m
if(phi1 > 2*PI) phi1 -= 2*PI;
if(phi1 < 0) phi1 += 2*PI;
Float_t phi2 = second->fPhi - asin(second->fCharge*(0.1*fBfield)*0.15/second->fPt);// 0.15 for D=1m 
if(phi2 > 2*PI) phi2 -= 2*PI;
if(phi2 < 0) phi2 += 2*PI;
  
Float_t deltaphi = phi1 - phi2;
if(deltaphi > PI) deltaphi -= 2*PI;
if(deltaphi < -PI) deltaphi += 2*PI;
deltaphi = fabs(deltaphi);

if(deltaphi < fMinSepPairPhi) return kFALSE;// Min Separation
    
  
// propagate through B field to r=1.6m
phi1 = first->fPhi - asin(first->fCharge*(0.1*fBfield)*0.24/first->fPt);// mine. 0.24 for D=1.6m
if(phi1 > 2*PI) phi1 -= 2*PI;
if(phi1 < 0) phi1 += 2*PI;
phi2 = second->fPhi - asin(second->fCharge*(0.1*fBfield)*0.24/second->fPt);// mine. 0.24 for D=1.6m 
if(phi2 > 2*PI) phi2 -= 2*PI;
if(phi2 < 0) phi2 += 2*PI;
  
deltaphi = phi1 - phi2;
if(deltaphi > PI) deltaphi -= 2*PI;
if(deltaphi < -PI) deltaphi += 2*PI;
deltaphi = fabs(deltaphi);

if(deltaphi < fMinSepPairPhi) return kFALSE;// Min Separation
  
  
   
//
  
Int_t ncl1 = first->fClusterMap.GetNbits();
Int_t ncl2 = second->fClusterMap.GetNbits();
Int_t sumCls = 0; Int_t sumSha = 0; Int_t sumQ = 0;
Double_t shfrac = 0; Double_t qfactor = 0;
for(Int_t imap = 0; imap < ncl1 && imap < ncl2; imap++) {
if (first->fClusterMap.TestBitNumber(imap) && second->fClusterMap.TestBitNumber(imap)) {// Both clusters
if (first->fSharedMap.TestBitNumber(imap) && second->fSharedMap.TestBitNumber(imap)) { // Shared
sumQ++;
sumCls+=2;
sumSha+=2;}
else {sumQ--; sumCls+=2;}
}
else if (first->fClusterMap.TestBitNumber(imap) || second->fClusterMap.TestBitNumber(imap)) {// Non shared
sumQ++;
sumCls++;}
}
if (sumCls>0) {
qfactor = sumQ*1.0/sumCls;
shfrac = sumSha*1.0/sumCls;
}
  
if(qfactor > fShareQuality || shfrac > fShareFraction) return kFALSE;
  
  
return kTRUE;
  

}
//________________________________________________________________________
Float_t AliAnalysisTaskFemtoESE::GamovFactor(Int_t chargeBin1, Int_t chargeBin2, Float_t qinv)
{
Float_t arg = G_Coeff/qinv;
  
if(chargeBin1==chargeBin2) return (exp(arg)-1)/(arg);
else {return (exp(-arg)-1)/(-arg);}
  
}
//________________________________________________________________________
void AliAnalysisTaskFemtoESE::Shuffle(Int_t *iarr, Int_t i1, Int_t i2)
{
Int_t j, k;
Int_t a = i2 - i1;
for (Int_t i = i1; i < i2+1; i++) {
j = (Int_t) (gRandom->Rndm() * a);
k = iarr[j];
iarr[j] = iarr[i];
iarr[i] = k;
}
}
*/
//________________________________________________________________________
Double_t AliAnalysisTaskFemtoESE::GetQinv(Double_t track1[], Double_t track2[]){
  
  Double_t qinv2=1.0;
  qinv2 = pow(track1[1]-track2[1],2.) + pow(track1[2]-track2[2],2.) + pow(track1[3]-track2[3],2.) - pow(track1[0]-track2[0],2.);
  if(qinv2 >= 0.) return sqrt(qinv2);
  else return -1.*sqrt(-1.*qinv2);
}

//________________________________________________________________________
void AliAnalysisTaskFemtoESE::GetQosl(Double_t track1[], Double_t track2[], Double_t& qout, Double_t& qside, Double_t& qlong){
 
  Double_t p0 = track1[0] + track2[0];
  Double_t px = track1[1] + track2[1];
  Double_t py = track1[2] + track2[2];
  Double_t pz = track1[3] + track2[3];
  
  Double_t mt = sqrt(p0*p0 - pz*pz);
  Double_t pt = sqrt(px*px + py*py);
  
  Double_t v0 = track1[0] - track2[0];
  Double_t vx = track1[1] - track2[1];
  Double_t vy = track1[2] - track2[2];
  Double_t vz = track1[3] - track2[3];
   
  if(gRandom->Rndm()<0.5)
    {
      v0 = -v0;
      vx = -vx;
      vy = -vy;
      vz = -vz;
    }

  //cout << p0 << "  " << px << "  " << py << "  " << pz << "  " << v0 << "  " << vx << "  " << vy << "  " << vz << "  " << mt << "  " << pt << endl;
  
  qout = (px*vx + py*vy)/pt;
  qside = (px*vy - py*vx)/pt;
  qlong = (p0*vz - pz*v0)/mt;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskFemtoESE::EventCut(/*AliAODEvent* fevent*/){

  // Trigger Cut

  Bool_t isSelected1 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  Bool_t isSelected2 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
  Bool_t isSelected3 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
  if(!isSelected1 && !isSelected2 && !isSelected3) {return kFALSE;}

  /*  
  // Pile-up rejection
  AliAnalysisUtils *AnaUtil=new AliAnalysisUtils();
  if(!fPbPbcase) AnaUtil->SetUseMVPlpSelection(kTRUE);// use Multi-Vertex tool for pp and pPb
  else AnaUtil->SetUseMVPlpSelection(kFALSE);
  Bool_t pileUpCase=AnaUtil->IsPileUpEvent(fAOD); 
  if(pileUpCase) return;

  // Vertexing
  primaryVertexAOD = fAOD->GetPrimaryVertex();
  vertex[0]=primaryVertexAOD->GetX(); vertex[1]=primaryVertexAOD->GetY(); vertex[2]=primaryVertexAOD->GetZ();
    
  if(fabs(vertex[2]) > 10) {cout<<"Zvertex Out of Range. Skip Event"<<endl; return;} // Z-Vertex Cut 
  ((TH3D*)fOutputList->FindObject("fVertexDist"))->Fill(vertex[0], vertex[1], vertex[2]);
    
  for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
  AliAODTrack* aodtrack = fAOD->GetTrack(i);
  if (!aodtrack) continue;
  AODTracks++;
  if(!aodtrack->TestFilterBit(BIT(fFilterBit))) continue;// AOD filterBit cut
  FBTracks++;
  }
  ((TH1D*)fOutputList->FindObject("fMultDist2"))->Fill(FBTracks);

  //if(fAOD->IsPileupFromSPD()) {cout<<"PileUpEvent. Skip Event"<<endl; return;} // Old Pile-up cut
  if(primaryVertexAOD->GetNContributors() < 1) {cout<<"Bad Vertex. Skip Event"<<endl; return;}
   
  ((TH1D*)fOutputList->FindObject("fMultDist3"))->Fill(FBTracks);
 
  fBfield = fAOD->GetMagneticField();
    
  for(Int_t i=0; i<fZvertexBins; i++){
  if( (vertex[2] >= zstart+i*zstep) && (vertex[2] < zstart+(i+1)*zstep) ){
  zbin=i;
  break;
  }
  }
  */

  return kTRUE;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskFemtoESE::TrackCut(AliAODTrack* ftrack){

  if (!ftrack->TestFilterBit(fFilterBit)) return kFALSE;

  if(ftrack->Pt() < fPtMin) return kFALSE;
  if(ftrack->Pt() > fPtMax) return kFALSE;
  if(fabs(ftrack->Eta()) > fEtaMax) return kFALSE;
 
  if(ftrack->GetTPCNcls() < 80) return kFALSE;// TPC nCluster cut

  //Double_t trackdca[3] = {ftrack->XAtDCA(),ftrack->YAtDCA(),ftrack->ZAtDCA()};
  //Double_t dcaxy = sqrt( pow(trackdca[0] - vertex[0],2) + pow(trackdca[1] - vertex[1],2));
  //Double_t dcaz = sqrt( pow(trackdca[2] - vertex[2],2));
  Double_t dcaxy = sqrt( pow(ftrack->XAtDCA(),2) + pow(ftrack->YAtDCA(),2));
  Double_t dcaz = fabs(ftrack->ZAtDCA());
  if(dcaxy > fMaxDcaXY) return kFALSE;
  if(dcaz > fMaxDcaZ) return kFALSE;
  hdcaxy->Fill(ftrack->XAtDCA(),ftrack->YAtDCA());
  hdcaz->Fill(ftrack->ZAtDCA());


  //// FilterBit Overlap Check
  //if(fFilterBit != 7){
  //	Bool_t goodTrackOtherFB = kFALSE;
  //	if(fMCcase && fAOD->GetRunNumber()<=126437) goodTrackOtherFB=kTRUE;// FB7 to FB5 mapping in 10f6a MC does not work
  //	
  //	for (Int_t j = 0; j < fAOD->GetNumberOfTracks(); j++) {
  //	  AliAODTrack* aodtrack2 = fAOD->GetTrack(randomIndex[j]);
  //	  if(!aodtrack2) continue;
  //	  if(!aodtrack2->TestFilterBit(BIT(fFilterBit))) continue;
  //	  
  //	  if(-(aodtrack->GetID()+1)==aodtrack2->GetID()) {goodTrackOtherFB=kTRUE; break;}
  //	  
  //	}
  //	if(!goodTrackOtherFB) continue;
  //    }

  return kTRUE;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskFemtoESE::PairCut(AliFemtoESEBasicParticle* ftrack1, AliFemtoESEBasicParticle* ftrack2, Bool_t mix){

  // check for same charge
  if(ftrack2->Charge() != ftrack1->Charge()) return kFALSE;
  
  // qinv cut
  if(fQinvMin > 0)
    {
      Double_t trackvec1[4] = {ftrack1->E(),ftrack1->Px(),ftrack1->Py(),ftrack1->Pz()};
      Double_t trackvec2[4] = {ftrack2->E(),ftrack2->Px(),ftrack2->Py(),ftrack2->Pz()};
      Double_t qinv = GetQinv(trackvec1,trackvec2);
      if(qinv < fQinvMin) return kFALSE; // qinv < 0.005
    }
  // deltaEta x deltaPhi* cut
  if(fabs(ftrack1->Eta()-ftrack2->Eta()) < fMinSepPairEta)
    {
      //Double_t deltaphistar = DeltaPhiStar(ftrack1,ftrack2,1.0); // angular separation at r=1m
      //deltaphistar = fabs(deltaphistar);
      //if(deltaphistar < fMinSepPairPhi) return kFALSE;
      Double_t deltaphistar = DeltaPhiStar(ftrack1,ftrack2,1.6); // angular separation at r=1.6m
      deltaphistar = fabs(deltaphistar);
      if(deltaphistar < fMinSepPairPhi) return kFALSE;
    }

  // share fraction & share quality cut
  TBits clusterMap1 = (TBits)(ftrack1->GetTPCClusterMap());
  TBits sharedMap1 = (TBits)(ftrack1->GetTPCSharedMap());
  TBits clusterMap2 = (TBits)(ftrack2->GetTPCClusterMap());
  TBits sharedMap2 = (TBits)(ftrack2->GetTPCSharedMap());

  Int_t ncl1 = clusterMap1.GetNbits();
  Int_t ncl2 = clusterMap2.GetNbits();
  Int_t sumCls = 0; Int_t sumSha = 0; Int_t sumQ = 0;
  Double_t shfrac = 0; Double_t qfactor = 0;
  for(Int_t imap = 0; imap < ncl1 && imap < ncl2; imap++)
    {
      if (clusterMap1.TestBitNumber(imap) && clusterMap2.TestBitNumber(imap)) // Both clusters
	{
	  if (sharedMap1.TestBitNumber(imap) && sharedMap2.TestBitNumber(imap)) // Shared
	    {
	      sumQ++;
	      sumCls+=2;
	      sumSha+=2;}
	  else {sumQ--; sumCls+=2;}
	}
      else if (clusterMap1.TestBitNumber(imap) || clusterMap2.TestBitNumber(imap)) // Non shared
	{
	  sumQ++;
	  sumCls++;
	}
    }
  if (sumCls>0)
    {
      qfactor = sumQ*1.0/sumCls;
      shfrac = sumSha*1.0/sumCls;
    }
  
  // sumCls -- number of clusters in track 1 + number of clusters in track 2 (clusters in both tracks counted twice)
  // sumSha -- number of shared clusters (counted twice)
  // sumQ -- ?

  if(!mix)
    {
      hsharequal->Fill(qfactor);
      hsharefrac->Fill(shfrac);
    }
  else
    {
      hsharequalmix->Fill(qfactor);
      hsharefracmix->Fill(shfrac);
    }

  
  if(qfactor > fShareQuality || shfrac > fShareFraction) return kFALSE;
  
  return kTRUE;
}

Double_t AliAnalysisTaskFemtoESE::DeltaPhiStar(AliAODTrack* ftrack1, AliAODTrack* ftrack2, Double_t r)
{
  Double_t phi1 = ftrack1->Phi() - asin(0.3*ftrack1->Charge()*(0.1*fBfield)*r/(2.*ftrack1->Pt())); // magnetic field converted from kilogauss to tesla
  Double_t phi2 = ftrack2->Phi() - asin(0.3*ftrack2->Charge()*(0.1*fBfield)*r/(2.*ftrack2->Pt())); 
  
  Double_t deltaphi = phi1 - phi2;
  while(deltaphi > TMath::Pi()) deltaphi -= 2*TMath::Pi();
  while(deltaphi < -TMath::Pi()) deltaphi += 2*TMath::Pi();

  return deltaphi;
}

/*TObjArray* AliAnalysisTaskFemtoESE::CloneAndReduceTrackList(TObjArray* tracks, Double_t psi)
  {
  // clones a track list by using AliDPhiBasicParticle which uses much less memory (used for event mixing)

  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
  AliAODTrack* particle = (AliAODTrack*) tracks->UncheckedAt(i);
  AliFemtoESEBasicParticle* copy = new AliFemtoESEBasicParticle(sqrt(pow(particle->P(),2)+pow(0.13957, 2)),particle->Px(),particle->Py(),particle->Pz(),particle->Charge(),particle->Phi(),particle->Eta());
  copy->SetPsiEP(psi);
  copy->SetTPCClusterMap(particle->GetTPCClusterMap());
  copy->SetTPCSharedMap(particle->GetTPCSharedMap());

  tracksClone->Add(copy);
  }
  
  return tracksClone;
  }*/

Double_t AliAnalysisTaskFemtoESE::GetDeltaPhiEP(Double_t px1, Double_t py1, Double_t px2, Double_t py2, Double_t psi)
{
  // angle of pair
  Double_t px = px1+px2;
  Double_t py = py1+py2;
  Double_t phi = atan2(py,px);

  Double_t dphi = phi-psi;
  while(dphi < 0) dphi += 2*TMath::Pi();
  while(dphi > 2*TMath::Pi()) dphi -= 2*TMath::Pi();

  return dphi;
}

Bool_t AliAnalysisTaskFemtoESE::FindBin(Double_t kt, Double_t phi, Double_t cent, Int_t& a, Int_t& b, Int_t&c)
{
  a = b = c = -1;
  for(Int_t i = 0; i < nKtBins; i++)
    {
      if(kt >= ktBins[i] && kt < ktBins[i+1])
	{
	  a = i;
	  break;
	}
    }
  if(a==-1) return kFALSE;

  if(phi > TMath::Pi()) phi = 2*TMath::Pi()-phi;
  for(Int_t i = 0; i < nEPBins; i++)
    {
      if(phi >= epBins[i] && phi < epBins[i+1])
	{
	  b = i;
	  break;
	}
    }
  if(b==-1) return kFALSE;

  for(Int_t i = 0; i < nCentBins; i++)
    {
      if(cent >= centBins[i] && cent < centBins[i+1])
	{
	  c = i;
	  break;
	}
    }
  if(c==-1) return kFALSE;

  return kTRUE;
}



void AliAnalysisTaskFemtoESE::SetKtBins(Int_t n, Double_t* bins)
{
  if(ktBins) delete [] ktBins;
  nKtBins = n;
  nKtBins1 = n+1;
  ktBins = new Double_t[nKtBins+1];
  for(Int_t i = 0; i < nKtBins+1; i++)
    ktBins[i]=bins[i];
  Printf("Setting %i kt bins: ",nKtBins);
  for(Int_t i = 0; i < nKtBins+1; i++) Printf("%lf",ktBins[i]);
}
void AliAnalysisTaskFemtoESE::SetCentBins(Int_t n, Double_t* bins)
{
  if(centBins) delete [] centBins;
  nCentBins = n;
  nCentBins1 = n+1;
  centBins = new Double_t[nCentBins+1];
  for(Int_t i = 0; i < nCentBins+1; i++)
    centBins[i]=bins[i];
  Printf("Setting %i centrality bins: ",nCentBins);
  for(Int_t i = 0; i < nCentBins+1; i++) Printf("%lf",centBins[i]);
}
void AliAnalysisTaskFemtoESE::SetVzBins(Int_t n, Double_t* bins)
{
  if(vzBins) delete [] vzBins;
  nVzBins = n;
  nVzBins1 = n+1;
  vzBins = new Double_t[nVzBins+1];
  for(Int_t i = 0; i < nVzBins+1; i++)
    vzBins[i]=bins[i];
  Printf("Setting %i vz bins: ",nVzBins);
  for(Int_t i = 0; i < nVzBins+1; i++) Printf("%lf",vzBins[i]);
}
void AliAnalysisTaskFemtoESE::SetEPBins(Int_t n)
{
  if(epBins) delete [] epBins;
  nEPBins = n;
  nEPBins1 = n+1;
  epBins = new Double_t[nEPBins+1];
  for(Int_t y = 0; y < nEPBins+1; y++)
    epBins[y] = (TMath::Pi()/(Double_t)n)*((Double_t)y);
  Printf("Setting %i EP bins: ",nEPBins);
  for(Int_t i = 0; i < nEPBins+1; i++) Printf("%lf",epBins[i]);
}

/*Double_t AliAnalysisTaskFemtoESE::GetQPercLHC11h(Double_t qvec)
{
  // preliminary attemp at calculating qvector percentile in LHC11h -- still very approximate and only works in 5% bins
  if(!qPercBinsLHC11h)
    {
      Double_t tempArray[21] = {0.0, 0.22995, 0.33047, 0.410831, 0.480728, 0.545566, 0.606841, 0.66634, 0.725193, 0.783813, 0.843311, 0.904185, 0.96796, 1.03522, 1.10768, 1.18774, 1.27808, 1.3857, 1.52438, 1.73633, 4.95};
      nqPercBinsLHC11h = 21;
      qPercBinsLHC11h = new Double_t[nqPercBinsLHC11h];
      for(Int_t n = 0; n < nqPercBinsLHC11h; n++) qPercBinsLHC11h[n] = tempArray[n];
    }

  for(Int_t t = 0; t < nqPercBinsLHC11h-1; t++)
    if(qvec > qPercBinsLHC11h[t] && qvec < qPercBinsLHC11h[t+1]) return 50.*(2*t+1)/(Double_t)(nqPercBinsLHC11h-1);

  if(qvec < qPercBinsLHC11h[0]) return 0.0;
  if(qvec > qPercBinsLHC11h[nqPercBinsLHC11h-1]) return 100.0;

  return 0.0;

  }*/

/*Double_t AliAnalysisTaskFemtoESE::GetCentralityWeight(Double_t cent)
{

  // use makeCentWeighting.C to fit centrality distribution to obtain this parameterization
  Double_t par1 = 2.60629;
  Double_t par2 = 0.579333;
  Double_t limit1 = 8.71488;
  Double_t limit2 = 12.0126;

  if(cent < limit1) return 1./par1;
  if(cent > limit2) return 1./par2;

  Double_t slope = (par2-par1)/(limit2-limit1);
  Double_t b = par1-slope*limit1;

  return 1./(slope*cent+b);
  }*/

Double_t AliAnalysisTaskFemtoESE::GetCentralityWeight(Double_t cent)
{

  // use makeCentWeighting.C to fit centrality distribution to obtain this parameterization
  Double_t par1 = 2.60629;
  Double_t par2 = 0.579333;
  Double_t limit1 = 8.71488;
  Double_t limit2 = 12.0126;

  if(cent < limit1) return 1.;
  if(cent > limit2) return 1.;

  Double_t slope = (par2-par1)/(limit2-limit1);
  Double_t b = par1-slope*limit1;

  if(cent < 10.)
    return par1/(slope*cent+b);
  else
    return par2/(slope*cent+b);
}
