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
  fEventCounter(0),
  fMixingTracks(10000),
  fBfield(0.),
  fMinSepPairEta(0.),
  fMinSepPairPhi(0.),
  fShareQuality(0.5),
  fShareFraction(0.05),
  nCountSamePairs(0),
  nCountMixedPairs(0),
  nCountTracks(0),
  fMinQPerc(-1000),
  fMaxQPerc(1000),
  fQPercDet(0),
  fEPDet(0),
  fPsiEPmix(0),
  fPsiEPmixtemp(0),
  nKtBins(0),
  nKtBins1(1),
  ktBins(0),
  nEPBins(0),
  nEPBins1(1),
  epBins(0),
  nCentBins(0),
  nCentBins1(1),
  centBins(0),
  nVzBins(0),
  nVzBins1(1),
  vzBins(0),
  hq(0x0),
  hqmix(0x0)
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
  fEventCounter(0),
  fMixingTracks(10000),
  fBfield(0.),
  fMinSepPairEta(0.),
  fMinSepPairPhi(0.),
  fShareQuality(0.5),
  fShareFraction(0.05),
  nCountSamePairs(0),
  nCountMixedPairs(0),
  nCountTracks(0),
  fMinQPerc(-1000),
  fMaxQPerc(1000),
  fQPercDet(0),
  fEPDet(0),
  fPsiEPmix(0),
  fPsiEPmixtemp(0),
  nKtBins(0),
  nKtBins1(1),
  ktBins(0),
  nEPBins(0),
  nEPBins1(1),
  epBins(0),
  nCentBins(0),
  nCentBins1(1),
  centBins(0),
  nVzBins(0),
  nVzBins1(1),
  vzBins(0),
  hq(0x0),
  hqmix(0x0)
{
  
  Printf("*******************************************");
  Printf("AliAnalysisTaskFemtoESE named %s",name);
  Printf("*******************************************");

  // default binning
  SetEPBins(12,-TMath::Pi()/12.,2*TMath::Pi()-TMath::Pi()/12.);
  Double_t ktBinsTemp[5] = {0.2,0.3,0.4,0.5,0.7};
  SetKtBins(4,ktBinsTemp);
  Double_t centBinsTemp[6] = {0,10,20,30,40,50};
  SetCentBins(5,centBinsTemp);
  Double_t vzBinsTemp[11] = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
  SetVzBins(10,vzBinsTemp);

  vertex[0] = vertex[1] = vertex[2] = 0.;

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliHelperPID::Class());

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
  fEventCounter(0),
  fMixingTracks(10000),
  fBfield(0.),
  fMinSepPairEta(0.),
  fMinSepPairPhi(0.),
  fShareQuality(0.5),
  fShareFraction(0.05),
  nCountSamePairs(0),
  nCountMixedPairs(0),
  nCountTracks(0),
  fMinQPerc(-1000),
  fMaxQPerc(1000),
  fQPercDet(0),
  fEPDet(0),
  fPsiEPmix(0),
  fPsiEPmixtemp(0),
  nKtBins(0),
  nKtBins1(1),
  ktBins(0),
  nEPBins(0),
  nEPBins1(1),
  epBins(0),
  nCentBins(0),
  nCentBins1(1),
  centBins(0),
  nVzBins(0),
  nVzBins1(1),
  vzBins(0),
  hq(0x0),
  hqmix(0x0)
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
  
  TH1D *hpx = new TH1D("hpx","px",200,-2,2);
  hpx->GetXaxis()->SetTitle("p_{x}");
  fOutputList->Add(hpx);
  TH1D *hpy = new TH1D("hpy","py",200,-2,2);
  hpy->GetXaxis()->SetTitle("p_{y}");
  fOutputList->Add(hpy);
  TH1D *hpz = new TH1D("hpz","pz",200,-2,2);
  hpz->GetXaxis()->SetTitle("p_{z}");
  fOutputList->Add(hpz);
  TH1D *hpt = new TH1D("hpt","pt",100,0,2);
  hpt->GetXaxis()->SetTitle("p_{t}");
  fOutputList->Add(hpt);
  TH1D *hE = new TH1D("hE","E",100,0,2);
  hE->GetXaxis()->SetTitle("E");
  fOutputList->Add(hE);
  TH2D *hphieta = new TH2D("hphieta","track #varphi vs #eta",100,0,2*TMath::Pi(),80,-0.8,0.8);
  hphieta->GetXaxis()->SetTitle("#varphi");
  hphieta->GetYaxis()->SetTitle("#eta");
  fOutputList->Add(hphieta);
  TH2D *hphieta_pid = new TH2D("hphieta_pid","PID check -- #Delta#varphi vs #Delta#eta",100,-0.3,0.3,100,-0.3,0.3);
  hphieta_pid->GetXaxis()->SetTitle("#Delta#varphi");
  hphieta_pid->GetYaxis()->SetTitle("#Delta#eta");
  fOutputList->Add(hphieta_pid);
  TH1D *hpt_pid = new TH1D("hpt_pid","PID check -- #Delta p_{t}",100,-0.5,0.5);
  hpt_pid->GetXaxis()->SetTitle("#Delta p_{t}");
  fOutputList->Add(hpt_pid);
  TH2D *hvzcent = new TH2D("hvzcent","vz vs cent",nVzBins,vzBins,nCentBins,centBins);
  hvzcent->GetXaxis()->SetTitle("v_{z}");
  hvzcent->GetYaxis()->SetTitle("centrality");
  fOutputList->Add(hvzcent);
  TH1D *hcent = new TH1D("hcent","cent",10,0,100);
  hcent->GetXaxis()->SetTitle("centrality");
  fOutputList->Add(hcent);
  TH2D *hcentn = new TH2D("hcentn","cent vs npions",50,0,50,100,0,2000);
  hcentn->GetXaxis()->SetTitle("Centrality");
  hcentn->GetYaxis()->SetTitle("Number of pions");
  fOutputList->Add(hcentn);
  TH3D *hphistaretapair10 = new TH3D("hphistaretapair10","pair #Delta#varphi* vs #Delta#eta at r=1.0m",100,-0.15,0.15,100,-0.1,0.1,10,0,1);
  hphistaretapair10->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair10->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair10->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair10);
  TH3D *hphistaretapair16 = new TH3D("hphistaretapair16","pair #Delta#varphi* vs #Delta#eta at r=1.6m",100,-0.15,0.15,100,-0.1,0.1,10,0,1);
  hphistaretapair16->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair16->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair16->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair16);
  TH3D *hphistaretapair10a = new TH3D("hphistaretapair10a","pair #Delta#varphi* vs #Delta#eta at r=1.0m",100,-0.15,0.15,100,-0.1,0.1,10,0,1);
  hphistaretapair10a->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair10a->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair10a->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair10a);
  TH3D *hphistaretapair16a = new TH3D("hphistaretapair16a","pair #Delta#varphi* vs #Delta#eta at r=1.6m",100,-0.15,0.15,100,-0.1,0.1,10,0,1);
  hphistaretapair16a->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair16a->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair16a->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair16a);
  TH3D *hphistaretapair10b = new TH3D("hphistaretapair10b","pair #Delta#varphi* vs #Delta#eta at r=1.0m",100,-TMath::Pi(),TMath::Pi(),100,-1.6,1.6,10,0,1);
  hphistaretapair10b->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair10b->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair10b->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair10b);
  TH3D *hphistaretapair16b = new TH3D("hphistaretapair16b","pair #Delta#varphi* vs #Delta#eta at r=1.6m",100,-TMath::Pi(),TMath::Pi(),100,-1.6,1.6,10,0,1);
  hphistaretapair16b->GetXaxis()->SetTitle("#Delta#varphi*");
  hphistaretapair16b->GetYaxis()->SetTitle("#Delta#eta");
  hphistaretapair16b->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphistaretapair16b);
  TH3D *hphietapair = new TH3D("hphietapair","pair #Delta#varphi vs #Delta#eta",100,-0.1,0.1,100,-0.1,0.1,10,0,1);
  hphietapair->GetXaxis()->SetTitle("#Delta#varphi");
  hphietapair->GetYaxis()->SetTitle("#Delta#eta");
  hphietapair->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphietapair);
  TH3D *hphietapair2 = new TH3D("hphietapair2","pair #varphi vs #eta",100,-TMath::Pi(),TMath::Pi(),100,-1.6,1.6,10,0,1);
  hphietapair2->GetXaxis()->SetTitle("#Delta#varphi");
  hphietapair2->GetYaxis()->SetTitle("#eta");
  hphietapair2->GetZaxis()->SetTitle("k_{T}");
  fOutputList->Add(hphietapair2);
  TH1D *hpidid = new TH1D("hpidid","pid id",9,-4.5,4.5);
  hpidid->GetXaxis()->SetTitle("track PID ID");
  fOutputList->Add(hpidid);
  TH1D *hkt = new TH1D("hkt","k_{T}",100,0,2);
  hkt->GetXaxis()->SetTitle("k_{T}");
  fOutputList->Add(hkt);
  TH1D *hktcheck = new TH1D("hktcheck","k_{T} check",100,0,2);
  hktcheck->GetXaxis()->SetTitle("k_{T}");
  fOutputList->Add(hktcheck);
  TH3D *hkt3 = new TH3D("hkt3","kt vs pt",50,0,1,50,0,5,50,0,5);
  hkt3->GetXaxis()->SetTitle("k_{T}");
  hkt3->GetYaxis()->SetTitle("p_{T,1}");
  hkt3->GetZaxis()->SetTitle("p_{T,2}");
  fOutputList->Add(hkt3);
  TH2D *hdcaxy = new TH2D("hdcaxy","DCA xy",100,-5,5,100,-5,5);
  hdcaxy->GetXaxis()->SetTitle("DCA x");
  hdcaxy->GetYaxis()->SetTitle("DCA y");
  fOutputList->Add(hdcaxy);
  TH1D *hdcaz = new TH1D("hdcaz","DCA z",100,-5,5);
  hdcaz->GetXaxis()->SetTitle("DCA z");
  fOutputList->Add(hdcaz);
  TH1D *hsharequal = new TH1D("hsharequal","Share Quality",102,-1.02,1.02);
  hsharequal->GetXaxis()->SetTitle("Share Quality");
  fOutputList->Add(hsharequal);
  TH1D *hsharefrac = new TH1D("hsharefrac","Share Fraction",100,0,1);
  hsharefrac->GetXaxis()->SetTitle("Share Fraction");
  fOutputList->Add(hsharefrac);
  TH1D *hsharequalmix = new TH1D("hsharequalmix","Share Quality -- mixed events",102,-1.02,1.02);
  hsharequalmix->GetXaxis()->SetTitle("Share Quality");
  fOutputList->Add(hsharequalmix);
  TH1D *hsharefracmix = new TH1D("hsharefracmix","Share Fraction -- mixed events",100,0,1);
  hsharefracmix->GetXaxis()->SetTitle("Share Fraction");
  fOutputList->Add(hsharefracmix);
  TH1D *hPsiTPC = new TH1D("hPsiTPC","TPC EP",100,-1*TMath::Pi(),TMath::Pi());
  hPsiTPC->GetXaxis()->SetTitle("#Psi{TPC}");
  fOutputList->Add(hPsiTPC);
  TH1D *hPsiV0A = new TH1D("hPsiV0A","V0A EP",100,-1*TMath::Pi(),TMath::Pi());
  hPsiV0A->GetXaxis()->SetTitle("#Psi{V0A}");
  fOutputList->Add(hPsiV0A);
  TH1D *hPsiV0C = new TH1D("hPsiV0C","V0C EP",100,-1*TMath::Pi(),TMath::Pi());
  hPsiV0C->GetXaxis()->SetTitle("#Psi{V0C}");
  fOutputList->Add(hPsiV0C);
  TH1D *hCheckEPA = new TH1D("hCheckEPA","Check EP V0A",100,-1*TMath::Pi(),TMath::Pi());
  hCheckEPA->GetXaxis()->SetTitle("PsiV0A - PsiTPC");
  fOutputList->Add(hCheckEPA);
  TH1D *hCheckEPC = new TH1D("hCheckEPC","Check EP V0C",100,-1*TMath::Pi(),TMath::Pi());
  hCheckEPC->GetXaxis()->SetTitle("PsiV0C - PsiTPC");
  fOutputList->Add(hCheckEPC);
  TH2D* hCheckEPmix = new TH2D("hCheckEPmix","Check EP mixed events",100,-1*TMath::Pi(),TMath::Pi(),100,-1*TMath::Pi(),TMath::Pi());
  hCheckEPmix->GetXaxis()->SetTitle("Psi1 - Psi_mix");
  hCheckEPmix->GetYaxis()->SetTitle("Psi1 - Psi2");
  fOutputList->Add(hCheckEPmix);
  TH2D *hcentq = new TH2D("hcentq","qvec vs cent",100,0,100,5,0,50);
  hcentq->GetXaxis()->SetTitle("q_{2} percentile");
  hcentq->GetYaxis()->SetTitle("centrality");
  fOutputList->Add(hcentq);
  TH2D* hMixedDist = new TH2D("hMixedDist", ";centrality;tracks;events", 101, 0, 101, 200, 0, fMixingTracks * 1.5);
  fOutputList->Add(hMixedDist);

  // resolution histograms (same binning as hvzcent)
  TH2D *hresV0ATPC = new TH2D("hresV0ATPC","vz vs cent vs cos(2*(V0A-TPC))",nVzBins,vzBins,nCentBins,centBins);
  hresV0ATPC->GetXaxis()->SetTitle("v_{z}");
  hresV0ATPC->GetYaxis()->SetTitle("centrality");
  fOutputList->Add(hresV0ATPC);
  TH2D *hresV0CTPC = new TH2D("hresV0CTPC","vz vs cent vs cos(2*(V0C-TPC))",nVzBins,vzBins,nCentBins,centBins);
  hresV0CTPC->GetXaxis()->SetTitle("v_{z}");
  hresV0CTPC->GetYaxis()->SetTitle("centrality");
  fOutputList->Add(hresV0CTPC);
  TH2D *hresV0AV0C = new TH2D("hresV0AV0C","vz vs cent vs cos(2*(V0A-V0C))",nVzBins,vzBins,nCentBins,centBins);
  hresV0AV0C->GetXaxis()->SetTitle("v_{z}");
  hresV0AV0C->GetYaxis()->SetTitle("centrality");
  fOutputList->Add(hresV0AV0C);

  hq = new TH3F****[nKtBins];
  hqmix = new TH3F****[nKtBins];
  for(Int_t k = 0; k < nKtBins; k++)
    {
      hq[k] = new TH3F***[nEPBins];
      hqmix[k] = new TH3F***[nEPBins];
      for(Int_t e = 0; e < nEPBins; e++)
	{
	  hq[k][e] = new TH3F**[nCentBins];
	  hqmix[k][e] = new TH3F**[nCentBins];
	  for(Int_t c = 0; c < nCentBins; c++)
	    {
	      hq[k][e][c] = new TH3F*[nVzBins];
	      hqmix[k][e][c] = new TH3F*[nVzBins];
	      for(Int_t v = 0; v < nVzBins; v++)
		{
		  hq[k][e][c][v] = new TH3F(Form("hq_%i_%i_%i_%i",k,e,c,v),Form("hq_%i_%i_%i_%i",k,e,c,v),20,-0.2,0.2,20,-0.2,0.2,20,-0.2,0.2);
		  fOutputList->Add(hq[k][e][c][v]);
		  hqmix[k][e][c][v] = new TH3F(Form("hqmix_%i_%i_%i_%i",k,e,c,v),Form("hqmix_%i_%i_%i_%i",k,e,c,v),20,-0.2,0.2,20,-0.2,0.2,20,-0.2,0.2);
		  fOutputList->Add(hqmix[k][e][c][v]);
		}
	    }
	}
    }

  // create dummy histograms which just hold the values of the kt, cent, vz, ep bin edges
  TH1F* hktbins = new TH1F("hktbins","kt bins",nKtBins,ktBins);
  fOutputList->Add(hktbins);
  TH1F* hcentbins = new TH1F("hcentbins","cent bins",nCentBins,centBins);
  fOutputList->Add(hcentbins);
  TH1F* hepbins = new TH1F("hepbins","ep bins",nEPBins,epBins);
  fOutputList->Add(hepbins);
  TH1F* hvzbins = new TH1F("hvzbins","vz bins",nVzBins,vzBins);
  fOutputList->Add(hvzbins);

  Printf("************************");
  Printf("using the %s detector for event plane determination",fEPDet ? "V0C" : "V0A");
  Printf("using the %s detector for q-vector determination",fQPercDet ? "V0C" : "V0A");
  Printf("************************");

  vertex[0] = vertex[1] = vertex[2] = 0.;

  // event mixing pool
  Int_t poolsize = 1000;
  fPoolMgr = new AliEventPoolManager(poolsize, fMixingTracks, nCentBins, centBins, nVzBins, vzBins);
  fPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5); // check these values


  nCountSamePairs = 0;
  nCountMixedPairs = 0;
  nCountTracks = 0;

  PostData(1, fOutputList);
  PostData(2, fHelperPID);
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
  
  primaryVertexAOD = fAOD->GetPrimaryVertex();
  vertex[0]=primaryVertexAOD->GetX(); vertex[1]=primaryVertexAOD->GetY(); vertex[2]=primaryVertexAOD->GetZ();
  Double_t zvtx = vertex[2];
  if(zvtx < vzBins[0] || zvtx > vzBins[nVzBins]) return; // Z-Vertex Cut 
  //cout<<"Centrality % = " << centralityPercentile << "  z-vertex = " << zvtx << endl;

  //ProcInfo_t procInfo;
  //cout << "beginning of event" << endl;
  //gSystem->GetProcInfo(&procInfo);
  //printf("ResMem %ld VMem %ld\n", procInfo.fMemResident, procInfo.fMemVirtual);

 // get event plane from V0's
  if(!fEventCuts->IsSelected(fAOD,fTrackCuts)) Printf("Error! Event not accepted by AliAODSpectraEventCuts!");
  //TVector2* qvecA = fEventCuts->GetqV0A();
  //TVector2* qvecC = fEventCuts->GetqV0C();
  //Double_t psiV0 = qvecA->Phi()/2.;
  //fEventCuts->CalculateQVectorLHC10h();
  Double_t psiV0A = fEventCuts->GetPsiV0A();
  Double_t psiV0C = fEventCuts->GetPsiV0C();
  //Double_t qV0a = fEventCuts->GetqV0A();
  Double_t qperc = fEventCuts->GetQvecPercentile(fQPercDet);//0: VZERO-A 1: VZERO-C
  //cout << "Psi = " << psiV0 << "   qV0a = " << qV0a << "   percentile = " << qperc << endl;

  if(psiV0A == -999) return;
  if(psiV0C == -999) return;
  if(qperc < fMinQPerc || qperc > fMaxQPerc) return;


  TH1D* hpx = (TH1D*)fOutputList->FindObject("hpx");
  TH1D* hpy = (TH1D*)fOutputList->FindObject("hpy");
  TH1D* hpz = (TH1D*)fOutputList->FindObject("hpz");
  TH1D* hpt = (TH1D*)fOutputList->FindObject("hpt");
  TH1D* hE = (TH1D*)fOutputList->FindObject("hE");
  TH2D* hphieta = (TH2D*)fOutputList->FindObject("hphieta");
  TH2D* hphieta_pid = (TH2D*)fOutputList->FindObject("hphieta_pid");
  TH1D* hpt_pid = (TH1D*)fOutputList->FindObject("hpt_pid");
  TH2D* hvzcent = (TH2D*)fOutputList->FindObject("hvzcent");
  TH1D* hcent = (TH1D*)fOutputList->FindObject("hcent");
  TH2D* hcentn = (TH2D*)fOutputList->FindObject("hcentn");
  TH3D* hphistaretapair10 = (TH3D*)fOutputList->FindObject("hphistaretapair10");
  TH3D* hphistaretapair16 = (TH3D*)fOutputList->FindObject("hphistaretapair16");
  TH3D* hphistaretapair10a = (TH3D*)fOutputList->FindObject("hphistaretapair10a");
  TH3D* hphistaretapair16a = (TH3D*)fOutputList->FindObject("hphistaretapair16a");
  TH3D* hphistaretapair10b = (TH3D*)fOutputList->FindObject("hphistaretapair10b");
  TH3D* hphistaretapair16b = (TH3D*)fOutputList->FindObject("hphistaretapair16b");
  TH3D* hphietapair = (TH3D*)fOutputList->FindObject("hphietapair");
  TH3D* hphietapair2 = (TH3D*)fOutputList->FindObject("hphietapair2");
  TH1D* hpidid = (TH1D*)fOutputList->FindObject("hpidid");
  TH1D* hkt = (TH1D*)fOutputList->FindObject("hkt");
  TH1D* hktcheck = (TH1D*)fOutputList->FindObject("hktcheck");
  TH3D* hkt3 = (TH3D*)fOutputList->FindObject("hkt3");
  TH1D* hPsiTPC = (TH1D*)fOutputList->FindObject("hPsiTPC");
  TH1D* hPsiV0A = (TH1D*)fOutputList->FindObject("hPsiV0A");
  TH1D* hPsiV0C = (TH1D*)fOutputList->FindObject("hPsiV0C");
  TH1D* hCheckEPA = (TH1D*)fOutputList->FindObject("hCheckEPA");
  TH1D* hCheckEPC = (TH1D*)fOutputList->FindObject("hCheckEPC");
  TH2D* hCheckEPmix = (TH2D*)fOutputList->FindObject("hCheckEPmix");
  TH2D* hcentq = (TH2D*)fOutputList->FindObject("hcentq");
  TH2D* hMixedDist = (TH2D*)fOutputList->FindObject("hMixedDist");
  TH2D* hresV0ATPC = (TH2D*)fOutputList->FindObject("hresV0ATPC");
  TH2D* hresV0CTPC = (TH2D*)fOutputList->FindObject("hresV0CTPC");
  TH2D* hresV0AV0C = (TH2D*)fOutputList->FindObject("hresV0AV0C");

  Double_t sin2phi = 0, cos2phi = 0;

  TObjArray* tracks = new TObjArray();
  //tracks->SetOwner(kTRUE);

  // Track loop -- select pions
  for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
    AliAODTrack* aodtrack = fAOD->GetTrack(i);
    if (!aodtrack) continue;
    if(!TrackCut(aodtrack)) continue;

    // filter bit 7 PID method...
    Int_t trackPID=999;
    for(Int_t m = 0; m < fAOD->GetNumberOfTracks(); m++) {
      AliAODTrack* aodtrack2 = fAOD->GetTrack(m);
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
      tracks->Add(aodtrack);

      // check event plane angle using tracks in the TPC
      if(aodtrack->Pt() < 2 && aodtrack->Pt() > 0.2)
	{
	  sin2phi += (aodtrack->Pt())*sin(2*aodtrack->Phi());
	  cos2phi += (aodtrack->Pt())*cos(2*aodtrack->Phi());
	}

  }
  // end track loop

  Int_t ntracks = tracks->GetEntriesFast();

  // get EP from TPC, just to check
  Double_t psiTPC = 0.;
  if(ntracks > 0)
    psiTPC = 0.5*atan2(sin2phi,cos2phi);
  else return;

  Double_t psiEP = psiV0A;
  if(fEPDet==1) psiEP = psiV0C;

  hPsiTPC->Fill(psiTPC);
  hPsiV0A->Fill(psiV0A);
  hPsiV0C->Fill(psiV0C);
  Double_t dphiEP = psiTPC-psiV0A;
  if(dphiEP>TMath::Pi()) dphiEP-=2*TMath::Pi();
  if(dphiEP<-TMath::Pi()) dphiEP+=2*TMath::Pi();
  hCheckEPA->Fill(dphiEP);
  dphiEP = psiTPC-psiV0C;
  if(dphiEP>TMath::Pi()) dphiEP-=2*TMath::Pi();
  if(dphiEP<-TMath::Pi()) dphiEP+=2*TMath::Pi();
  hCheckEPC->Fill(dphiEP);


  hcentq->Fill(qperc,centralityPercentile);


  Double_t kt = 0;
  Double_t qout=0, qside=0, qlong=0;
  Double_t pVect1[4] = {0,0,0,0};
  Double_t pVect2[4] = {0,0,0,0};
  Int_t k, e, c, v; //bin indices for histograms

  nCountTracks += ntracks;
  //cout << "Found " << ntracks << " pion tracks..." << endl;

  hvzcent->Fill(zvtx,centralityPercentile);
  hcent->Fill(centralityPercentile);
  hcentn->Fill(centralityPercentile,ntracks);

  // resolution histograms
  hresV0ATPC->Fill(zvtx,centralityPercentile,cos(2*(psiV0A-psiTPC)));
  hresV0CTPC->Fill(zvtx,centralityPercentile,cos(2*(psiV0C-psiTPC)));
  hresV0AV0C->Fill(zvtx,centralityPercentile,cos(2*(psiV0A-psiV0C)));

  AliEventPool* pool = fPoolMgr->GetEventPool(centralityPercentile,zvtx);
  if (!pool) AliFatal(Form("No pool found for centrality = %f, vz = %f", centralityPercentile, zvtx));
  //if (!pool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f, Psi_EP = %f", centralityPercentile, zvtx, psiEP));
  //if(pool->IsReady()) hMixedDist->Fill(centralityPercentile, pool->NTracksInPool());


  for(Int_t j = 0; j < ntracks; j++)
    {
      //cout << endl << j << "   ";
      AliAODTrack* track1 = (AliAODTrack*)tracks->At(j);
      pVect1[0]=sqrt(pow(track1->P(),2)+pow(0.13957, 2));
      pVect1[1]=track1->Px();
      pVect1[2]=track1->Py();
      pVect1[3]=track1->Pz();
      //cout << pVect1[0] << "   " << pVect1[1] << "   " <<  pVect1[2] << "   " << pVect1[3] << endl;

      // track qa plots
      hpx->Fill(pVect1[1]);
      hpy->Fill(pVect1[2]);
      hpz->Fill(pVect1[3]);
      hpt->Fill(track1->Pt());
      hE->Fill(pVect1[0]);
      hphieta->Fill(track1->Phi(),track1->Eta());

      // same event
      for(Int_t i = j+1; i < ntracks; i++)
	{
	  AliAODTrack* track2 = (AliAODTrack*)tracks->At(i);

	  hphistaretapair10->Fill(DeltaPhiStar(track1,track2,1.0),track1->Eta()-track2->Eta(),kt);
	  hphistaretapair16->Fill(DeltaPhiStar(track1,track2,1.6),track1->Eta()-track2->Eta(),kt);

	  if(!PairCut(track1,track2,kFALSE)) continue;

	  hphistaretapair10a->Fill(DeltaPhiStar(track1,track2,1.0),track1->Eta()-track2->Eta(),kt);
	  hphistaretapair16a->Fill(DeltaPhiStar(track1,track2,1.6),track1->Eta()-track2->Eta(),kt);
	  hphistaretapair10b->Fill(DeltaPhiStar(track1,track2,1.0),track1->Eta()-track2->Eta(),kt);
	  hphistaretapair16b->Fill(DeltaPhiStar(track1,track2,1.6),track1->Eta()-track2->Eta(),kt);

	  pVect2[0]=sqrt(pow(track2->P(),2)+pow(0.13957, 2));
	  pVect2[1]=track2->Px();
	  pVect2[2]=track2->Py();
	  pVect2[3]=track2->Pz();

	  //qinv = GetQinv(pVect1, pVect2); // = qinv**2 = (P1x-P2x)**2 + (P1y-P2y)**2 + (P1z-P2z)**2 - (P1t-P2t)**2 
	  GetQosl(pVect1, pVect2, qout, qside, qlong); // qout, qside, qlong = components of Q=P1-P2 in the P=P1+P2 frame
	  kt = sqrt(pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.; // = Kt = |pT1+pT2|/2
	  hkt->Fill(kt);
	  hkt3->Fill(kt,track1->Pt(),track2->Pt());
	  Double_t deltaphi = GetDeltaPhiEP(pVect1[1],pVect1[2],pVect2[1],pVect2[2],psiEP); // angle to event plane in correct range
	  if(fabs(qout)<0.2 && fabs(qside)<0.2 && fabs(qlong)<0.2) nCountSamePairs++;
	  if(kt < ktBins[0] || kt > ktBins[nKtBins]) continue;
	  if(!FindBin(kt,deltaphi,centralityPercentile,zvtx,k,e,c,v)) continue;
	  hktcheck->Fill(kt);
	  hq[k][e][c][v]->Fill(qout,qside,qlong);
	  Double_t dphi = track1->Phi()-track2->Phi();
	  if(dphi<-TMath::Pi()) dphi += 2*TMath::Pi();
	  if(dphi>TMath::Pi()) dphi -= 2*TMath::Pi();
	  hphietapair->Fill(dphi,track1->Eta()-track2->Eta(),kt);
	  hphietapair2->Fill(dphi,track1->Eta()-track2->Eta(),kt);
	  //cout << k << "  ";
	}

      // mixed event
      if (pool->IsReady()) 
	{
	  Int_t nMix = pool->GetCurrentNEvents();
	  Int_t countmix = 0;

	  for (Int_t jMix=0; jMix<nMix; jMix++) 
	    {
	      TObjArray* bgTracks = pool->GetEvent(jMix);
	      Int_t ntracksmix = bgTracks->GetEntriesFast();

	      if(ntracksmix > 0)
		{
		  AliFemtoESEBasicParticle* tracktest = (AliFemtoESEBasicParticle*)bgTracks->UncheckedAt(0);
		  fPsiEPmixtemp = tracktest->GetPsiEP();
		  Double_t dphiEPtest = fPsiEPmixtemp-psiEP;
		  while(dphiEPtest>2*TMath::Pi()) dphiEPtest-=2*TMath::Pi();
		  while(dphiEPtest<0) dphiEPtest+=2*TMath::Pi();
		  if(dphiEPtest>TMath::Pi()) dphiEPtest-=TMath::Pi();
		  if(dphiEPtest>TMath::Pi()/2.) dphiEPtest = TMath::Pi()-dphiEPtest;
		  if(dphiEPtest > TMath::Pi()/6.) continue;
		  countmix += ntracksmix;

		  fPsiEPmix = 0.5*atan2(sin(2*psiEP)+sin(2*fPsiEPmixtemp),cos(2*psiEP)+cos(2*fPsiEPmixtemp));
		  Double_t dphimix = psiEP-fPsiEPmix;
		  if(dphimix < -TMath::Pi()) dphimix += 2*TMath::Pi();
		  if(dphimix > TMath::Pi()) dphimix -= 2*TMath::Pi();
		  Double_t dphi12 = psiEP-fPsiEPmixtemp;
		  if(dphi12 < -TMath::Pi()) dphi12 += 2*TMath::Pi();
		  if(dphi12 > TMath::Pi()) dphi12 -= 2*TMath::Pi();
		  hCheckEPmix->Fill(dphimix,dphi12);
		}


	      //cout << "mixing with " << ntracksmix << " tracks" << endl;
	      for(Int_t i = 0; i < ntracksmix; i++)
		{
		  AliFemtoESEBasicParticle* track2 = (AliFemtoESEBasicParticle*)bgTracks->UncheckedAt(i);

		  if(!PairCut(track1,track2,kTRUE)) continue;

		  pVect2[0]=track2->E();
		  pVect2[1]=track2->Px();
		  pVect2[2]=track2->Py();
		  pVect2[3]=track2->Pz();

		  if(fPsiEPmixtemp != track2->GetPsiEP()) AliFatal("Error! Event plane angles are wrong in mixing!!");

		  //qinv = GetQinv(pVect1, pVect2); // qinv**2 = (P1x-P2x)**2 + (P1y-P2y)**2 + (P1z-P2z)**2 - (P1t-P2t)**2 
		  GetQosl(pVect1, pVect2, qout, qside, qlong); // qout, qside, qlong = components of Q=P1-P2 in the P=P1+P2 frame
		  kt = sqrt(pow(pVect1[1]+pVect2[1],2) + pow(pVect1[2]+pVect2[2],2))/2.; // = Kt = |pT1+pT2|/2
		  Double_t deltaphi = GetDeltaPhiEP(pVect1[1],pVect1[2],pVect2[1],pVect2[2],fPsiEPmix); // angle to event plane in correct range

		  //Double_t weight = 1./(Double_t)nMix;
		  if(fabs(qout)<0.2 && fabs(qside)<0.2 && fabs(qlong)<0.2) nCountMixedPairs++;
		  if(kt < ktBins[0] || kt > ktBins[nKtBins]) continue;
		  if(!FindBin(kt,deltaphi,centralityPercentile,zvtx,k,e,c,v)) continue;
		  hqmix[k][e][c][v]->Fill(qout,qside,qlong);

		}
	    }

	  hMixedDist->Fill(centralityPercentile, countmix);

	}
    }

  TObjArray* clonedtracks = CloneAndReduceTrackList(tracks,psiEP);
  pool->UpdatePool(clonedtracks);
  //cout << "pool contains " << pool->GetCurrentNEvents() << " events and " << pool->NTracksInPool() << " tracks." << endl;
  //tracks->Clear();
  
  delete tracks;

  // Post output data.
  PostData(1, fOutputList);
  PostData(2, fHelperPID);
}
//________________________________________________________________________
void AliAnalysisTaskFemtoESE::Terminate(Option_t *) 
{

  if(ktBins) delete [] ktBins;
  if(epBins) delete [] epBins;
  if(centBins) delete [] centBins;
  if(vzBins) delete [] vzBins;

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
  
  Double_t qinv=1.0;
  qinv = sqrt( pow(track1[1]-track2[1],2) + pow(track1[2]-track2[2],2) + pow(track1[3]-track2[3],2) - pow(track1[0]-track2[0],2));
  return qinv;
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

  if(ftrack->Pt() < 0.14) return kFALSE;
  if(ftrack->Pt() > 1.5) return kFALSE;
  if(fabs(ftrack->Eta()) > 0.8) return kFALSE;
 
  if(ftrack->GetTPCNcls() < 80) return kFALSE;// TPC nCluster cut

  Double_t trackdca[3] = {ftrack->XAtDCA(),ftrack->YAtDCA(),ftrack->ZAtDCA()};
  //ftrack->XYZAtDCA(trackdca);
  //Double_t dcaxy = sqrt( pow(trackpos[0] - vertex[0],2) + pow(trackpos[1] - vertex[1],2));
  //Double_t dcaz = sqrt( pow(trackpos[2] - vertex[2],2));
  ((TH2D*)fOutputList->FindObject("hdcaxy"))->Fill(trackdca[0],trackdca[1]);
  ((TH1D*)fOutputList->FindObject("hdcaz"))->Fill(trackdca[2]);
  //if(dcaxy > 0.2) return kFALSE;
  //if(dcaz > 0.15) return kFALSE;


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
Bool_t AliAnalysisTaskFemtoESE::PairCut(AliAODTrack* ftrack1, AliAODTrack* ftrack2, Bool_t mix){

  // check for same charge
  if(ftrack2->Charge() != ftrack1->Charge()) return kFALSE;
  
  // qinv cut
  Double_t trackvec1[4] = {sqrt(pow(ftrack1->P(),2)+pow(0.13957, 2)),ftrack1->Px(),ftrack1->Py(),ftrack1->Pz()};
  Double_t trackvec2[4] = {sqrt(pow(ftrack2->P(),2)+pow(0.13957, 2)),ftrack2->Px(),ftrack2->Py(),ftrack2->Pz()};
  Double_t qinv = GetQinv(trackvec1,trackvec2);
  if(qinv < 0.005) return kFALSE;

  // deltaEta x deltaPhi* cut
  if(fabs(ftrack1->Eta()-ftrack2->Eta()) < fMinSepPairEta)
    {
      Double_t deltaphistar = DeltaPhiStar(ftrack1,ftrack2,1.0); // angular separation at r=1m
      deltaphistar = fabs(deltaphistar);
      if(deltaphistar < fMinSepPairPhi) return kFALSE;
      deltaphistar = DeltaPhiStar(ftrack1,ftrack2,1.6); // angular separation at r=1.6m
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
  
  TH1D* hsharequal;
  TH1D* hsharefrac;

  if(!mix)
    {
      hsharequal = (TH1D*)fOutputList->FindObject("hsharequal");
      hsharefrac = (TH1D*)fOutputList->FindObject("hsharefrac");
    }
  else
    {
     hsharequal = (TH1D*)fOutputList->FindObject("hsharequalmix");
     hsharefrac = (TH1D*)fOutputList->FindObject("hsharefracmix");
    }
  hsharequal->Fill(qfactor);
  hsharefrac->Fill(shfrac);
  
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

TObjArray* AliAnalysisTaskFemtoESE::CloneAndReduceTrackList(TObjArray* tracks, Double_t psi)
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
}

Double_t AliAnalysisTaskFemtoESE::GetDeltaPhiEP(Double_t px1, Double_t py1, Double_t px2, Double_t py2, Double_t psi)
{
  // angle of pair
  Double_t px = px1+px2;
  Double_t py = py1+py2;
  Double_t phi = atan2(py,px);

  Double_t dphi = phi-psi;
  while(dphi < epBins[0]) dphi += 2*TMath::Pi();
  while(dphi > epBins[nEPBins]) dphi -= 2*TMath::Pi();

  return dphi;
}

Bool_t AliAnalysisTaskFemtoESE::FindBin(Double_t kt, Double_t phi, Double_t cent, Double_t vz, Int_t& a, Int_t& b, Int_t&c, Int_t& d)
{
  a = b = c = d = -1;
  for(Int_t i = 0; i < nKtBins; i++)
    {
      if(kt >= ktBins[i] && kt < ktBins[i+1])
	{
	  a = i;
	  break;
	}
    }
  if(a==-1) return kFALSE;

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

  for(Int_t i = 0; i < nVzBins; i++)
    {
      if(vz >= vzBins[i] && vz < vzBins[i+1])
	{
	  d = i;
	  break;
	}
    }
  if(d==-1) return kFALSE;

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
void AliAnalysisTaskFemtoESE::SetEPBins(Int_t n, Double_t min, Double_t max)
{
  if(epBins) delete [] epBins;
  nEPBins = n;
  nEPBins1 = n+1;
  epBins = new Double_t[nEPBins+1];
  for(Int_t y = 0; y < nEPBins+1; y++)
    epBins[y] = min+((max-min)/(Double_t)n)*((Double_t)y);
  Printf("Setting %i EP bins: ",nEPBins);
  for(Int_t i = 0; i < nEPBins+1; i++) Printf("%lf",epBins[i]);
}

