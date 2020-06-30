/*   This macro produces:  pT spectra in different multiplicity and Delta phi bins
     Aditya Nath Mishra Wigner RCP, Budapest, Hungary
     Please report bugs to: amishra@cern.ch / aditya.nath.mishra@wigner.hu
     last update: 16/06/2020

*/

#include "AliAnalysisTaskUeSpectraDphi.h"

// ROOT includes
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2D.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TFile.h>


// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliEventCuts.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include <TTreeStream.h>
#include <AliHeader.h>
#include <AliAnalysisUtils.h>
#include <AliMultiplicity.h>
#include <AliMultSelection.h>
#include <AliMultSelectionTask.h>
#include <AliAODInputHandler.h>
#include <AliAODHandler.h>
#include <AliAODVertex.h>
#include <AliAODTrack.h>
#include <AliAODPid.h>
#include <AliDataFile.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>

#include <iostream>
using namespace std;

const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
const Int_t ndPhiBins =18;
Double_t dPhiBinsOA[ndPhiBins] ={0.174533,0.349066,0.523599,0.698132,0.872665,1.0472,1.22173,1.39626,1.5708,1.74533,1.91986,2.0944,2.26893,2.44346,2.61799,2.79253,2.96706,3.15};
Double_t dPhiBinsSA[ndPhiBins+1] ={0,0.174533,0.349066,0.523599,0.698132,0.872665,1.0472,1.22173,1.39626,1.5708,1.74533,1.91986,2.0944,2.26893,2.44346,2.61799,2.79253,2.96706,3.15};

ClassImp(AliAnalysisTaskUeSpectraDphi)

//_____________________________________________________________________________
AliAnalysisTaskUeSpectraDphi::AliAnalysisTaskUeSpectraDphi():
AliAnalysisTaskSE(),
  fESD(0x0),
  fEventCuts(0x0),
  fMCEvent(0x0),
  fMCStack(0x0),
  fAnalysisMC(kFALSE),
  fTrackFilter(0x0),
  ftrigBit(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),
  fEtaCut(0.8),
  fNcl(70),
  fTriggeredEventMB(-999),
  fListOfObjects(0x0),
  fisPS(kFALSE),
  fZvtx_SPD(0x0),
  fMCVzCut(0x0),
  fisTracklet(kFALSE),
  fVtxBeforeCuts(0x0),
  fVtxAfterCuts(0x0),
  hSelEv(0x0),
  hINEL0(0x0),
  hPS(0x0),
  hVtxPS(0x0),
  hpT(0x0),
  hEta(0x0),
  hPhi(0x0),
  hPtL(0x0),
  hPhiL(0x0),
  hEtaL(0x0),
  hDphi(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  hpTDphiBinsSAWLP(0x0),
  hpTvsDphiOA(0x0),
  hpTvsDphiSA(0x0),
  hMultvsDphiOA(0x0),
  hMultvsDphiSA(0x0),
  hMultvspTvsDphi(0x0),
  isINEL0True(0x0),
  isINEL0Rec(0x0),
  hINEL0MCTrig(0x0),
  hINEL0MCTrue(0x0),
  hPS_MC(0x0),
  hVtxPS_MC(0x0),
  hpTMCTrue(0x0),
  hEtaMCTrue(0x0),
  hPhiMCTrue(0x0),
  hPtLMCTrue(0x0),
  hPhiLMCTrue(0x0),
  hEtaLMCTrue(0x0),
  hDphiMCTrue(0x0),
  hpTDphiBinsSAWLPMCTrue(0x0),
  hpTvsDphiOAMCTrue(0x0),
  hpTvsDphiSAMCTrue(0x0),
  hMultvsDphiOAMCTrue(0x0),
  hMultvsDphiSAMCTrue(0x0),
  hMultvspTvsDphiMCTrue(0x0),
  fMultSelection(0x0),
  ftrackmult08(-999),
  fv0mpercentile(-999)
{
  for(Int_t i=0; i<18; i++){
    hpTDphiBinsOA[i] = 0;
    hpTDphiBinsSA[i] = 0;
    hMultDphiBinsOA[i] = 0;
    hMultDphiBinsSA[i] = 0;

    hpTDphiBinsOAMCTrue[i] = 0;
    hpTDphiBinsSAMCTrue[i] = 0;
    hMultDphiBinsOAMCTrue[i] = 0;
    hMultDphiBinsSAMCTrue[i] = 0;
  }
  // Default constructor (should not be used)
}

//______________________________________________________________________________
AliAnalysisTaskUeSpectraDphi::AliAnalysisTaskUeSpectraDphi(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fEventCuts(0x0),
  fMCEvent(0x0),
  fMCStack(0x0),
  fAnalysisMC(kFALSE),
  fTrackFilter(0x0),
  ftrigBit(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),
  fEtaCut(0.8),
  fNcl(70),
  fTriggeredEventMB(-999),
  fListOfObjects(0x0),
  fisPS(kFALSE),
  fZvtx_SPD(0x0),
  fMCVzCut(0x0),
  fisTracklet(kFALSE),
  fVtxBeforeCuts(0x0),
  fVtxAfterCuts(0x0),
  hSelEv(0x0),
  hINEL0(0x0),
  hpT(0x0),
  hEta(0x0),
  hPhi(0x0),
  hPtL(0x0),
  hPhiL(0x0),
  hEtaL(0x0),
  hDphi(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  hpTDphiBinsSAWLP(0x0),
  hpTvsDphiOA(0x0),
  hpTvsDphiSA(0x0),
  hMultvsDphiOA(0x0),
  hMultvsDphiSA(0x0),
  hMultvspTvsDphi(0x0),
  isINEL0True(0x0),
  isINEL0Rec(0x0),
  hINEL0MCTrig(0x0),
  hINEL0MCTrue(0x0),
  hPS_MC(0x0),
  hVtxPS_MC(0x0),
  hpTMCTrue(0x0),
  hEtaMCTrue(0x0),
  hPhiMCTrue(0x0),
  hPtLMCTrue(0x0),
  hPhiLMCTrue(0x0),
  hEtaLMCTrue(0x0),
  hDphiMCTrue(0x0),
  hpTDphiBinsSAWLPMCTrue(0x0),
  hpTvsDphiOAMCTrue(0x0),
  hpTvsDphiSAMCTrue(0x0),
  hMultvsDphiOAMCTrue(0x0),
  hMultvsDphiSAMCTrue(0x0),
  hMultvspTvsDphiMCTrue(0x0),
  fMultSelection(0x0),
  ftrackmult08(-999),
  fv0mpercentile(-999)
{
  for(Int_t i=0; i<18; i++){
    hpTDphiBinsOA[i] = 0;
    hpTDphiBinsSA[i] = 0;
    hMultDphiBinsOA[i] = 0;
    hMultDphiBinsSA[i] = 0;

    hpTDphiBinsOAMCTrue[i] = 0;
    hpTDphiBinsSAMCTrue[i] = 0;
    hMultDphiBinsOAMCTrue[i] = 0;
    hMultDphiBinsSAMCTrue[i] = 0;
  }
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskUeSpectraDphi::Exit(const char *msg) {

  Printf("%s", msg);
  return;
}


//_____________________________________________________________________________
AliAnalysisTaskUeSpectraDphi::~AliAnalysisTaskUeSpectraDphi()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fListOfObjects;
    fListOfObjects = 0x0;
  }

}

//______________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::UserCreateOutputObjects()
{

  const Int_t nPtBins      = 70;
  Double_t PtBins[nPtBins+1] = {
    0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
    0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,
    1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
    4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20,22.0,
    24.0,26.0,28.0,32.0,36.0,42.0,50.0,60.0,80.0,100.0,130.0,
    160.0,200.0};

   const Int_t nMultBins = 200;
  const Double_t MultBins[nMultBins+1] = {
    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,
    40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,
    77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,
    111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,
    140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,
    169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,
    198,199,200};

  const Int_t nDphiBins = 180;
  const Double_t DphiBins[nDphiBins+1] = {
    0,0.01745,0.0349,0.05235,0.0698,0.08725,0.1047,0.12215,0.1396,0.15705,0.1745,0.19195,0.2094,0.22685,0.2443,0.26175,
    0.2792,0.29665,0.3141,0.33155,0.349,0.36645,0.3839,0.40135,0.4188,0.43625,0.4537,0.47115,0.4886,0.50605,0.5235,0.54095,
    0.5584,0.57585,0.5933,0.61075,0.6282,0.64565,0.6631,0.68055,0.698,0.71545,0.7329,0.75035,0.7678,0.78525,0.8027,0.82015,
    0.8376,0.85505,0.8725,0.88995,0.9074,0.92485,0.9423,0.95975,0.9772,0.99465,1.0121,1.02955,1.047,1.06445,1.0819,1.09935,
    1.1168,1.13425,1.1517,1.16915,1.1866,1.20405,1.2215,1.23895,1.2564,1.27385,1.2913,1.30875,1.3262,1.34365,1.3611,1.37855,
    1.396,1.41345,1.4309,1.44835,1.4658,1.48325,1.5007,1.51815,1.5356,1.55305,1.5705,1.58795,1.6054,1.62285,1.6403,1.65775,
    1.6752,1.69265,1.7101,1.72755,1.745,1.76245,1.7799,1.79735,1.8148,1.83225,1.8497,1.86715,1.8846,1.90205,1.9195,1.93695,
    1.9544,1.97185,1.9893,2.00675,2.0242,2.04165,2.0591,2.07655,2.094,2.11145,2.1289,2.14635,2.1638,2.18125,2.1987,2.21615,
    2.2336,2.25105,2.2685,2.28595,2.3034,2.32085,2.3383,2.35575,2.3732,2.39065,2.4081,2.42555,2.443,2.46045,2.4779,2.49535,
    2.5128,2.53025,2.5477,2.56515,2.5826,2.60005,2.6175,2.63495,2.6524,2.66985,2.6873,2.70475,2.7222,2.73965,2.7571,2.77455,
    2.792,2.80945,2.8269,2.84435,2.8618,2.87925,2.8967,2.91415,2.9316,2.94905,2.9665,2.98395,3.0014,3.01885,3.0363,3.05375,
    3.0712,3.08865,3.1061,3.12355,3.15
   };

  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested

  // Definition of trackcuts
  if(!fTrackFilter){
    fTrackFilter = new AliAnalysisFilter("trackFilter2015");
    SetTrackCuts(fTrackFilter);
  }

  OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  //
  // Histograms
  //

  hSelEv = 0;
  hSelEv = new TH1D("hSelEv","Number of events; Events; Counts", 10, 0, 10);
  fListOfObjects->Add(hSelEv);

  hINEL0 = 0;
  hINEL0 = new TH1D("hINEL0","Number of events; Events; Counts", 10, 0, 10);
  fListOfObjects->Add(hINEL0);

  hPS = new TH1D("hPS","",5,0,5);
  fListOfObjects->Add(hPS);
  
  hVtxPS = new TH1D("hVtxPS","",5,0,5);
  fListOfObjects->Add(hVtxPS);

  fVtxBeforeCuts = 0;
  fVtxBeforeCuts = new TH1D("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 300, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);

  fVtxAfterCuts = 0;
  fVtxAfterCuts = new TH1D("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 300, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);

  hRefMult08 = 0;
  hRefMult08 = new TH1D("hRefMult08","Multiplicity (-0.8 < #eta < 0.8);N_{ch};count",nMultBins,MultBins);
  fListOfObjects->Add(hRefMult08);
  
  hV0Mmult = 0;
  hV0Mmult = new TH1D("hV0Mmult","V0M ;V0M percentile;count",110,0,110);
  fListOfObjects->Add(hV0Mmult);

  if (fAnalysisMC){
    hINEL0MCTrig = 0;
    hINEL0MCTrig = new TH1D("hINEL0MCTrig","Number of events; Events; Counts", 10, 0, 10);
    fListOfObjects->Add(hINEL0MCTrig);
    
    hINEL0MCTrue = 0;
    hINEL0MCTrue = new TH1D("hINEL0MCTrue","Number of events; Events; Counts", 10, 0, 10);
    fListOfObjects->Add(hINEL0MCTrue);

    hPS_MC = new TH1D("hPS_MC","",5,0,5);
    fListOfObjects->Add(hPS_MC);
    
    hVtxPS_MC = new TH1D("hVtxPS_MC","",5,0,5);
    fListOfObjects->Add(hVtxPS_MC);

    hpTMCTrue = new TH1D("hpTMCTrue","",nPtBins,PtBins);
    fListOfObjects->Add(hpTMCTrue);

    hEtaMCTrue = new TH1D("hEtaMCTrue","; #eta^{leading};counts",20,-1,1);
    fListOfObjects->Add(hEtaMCTrue);

    hPhiMCTrue = new TH1D("hPhiMCTrue", ";#phi (rad); count", 64,0,2.0*TMath::Pi());
    fListOfObjects->Add(hPhiMCTrue);

    hPtLMCTrue = 0;
    hPtLMCTrue = new TH1D("hPtLMCTrue",";#it{p}_{T}^{leading} (GeV/#it{c});counts",nPtBins,PtBins);
    fListOfObjects->Add(hPtLMCTrue);

    hEtaLMCTrue = 0;
    hEtaLMCTrue = new TH1D("hEtaLMCTrue","; #eta^{leading};counts",20,-1,1);
    fListOfObjects->Add(hEtaLMCTrue);

    hPhiLMCTrue = 0;
    hPhiLMCTrue = new TH1D("hPhiLMCTrue","; #phi^{leading} (rad);counts",64,0,2.0*TMath::Pi());
    fListOfObjects->Add(hPhiLMCTrue);

    hDphiMCTrue = 0;
    hDphiMCTrue = new TH1D("hDphiMCTrue","",64,-2*TMath::Pi(),2*TMath::Pi());
    fListOfObjects->Add(hDphiMCTrue);

    hpTDphiBinsSAWLPMCTrue = 0;
    hpTDphiBinsSAWLPMCTrue = new TH1D("hpTDphiBinsSAWLPMCTrue","Charged particles WLP (sliding #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBins,PtBins);
    fListOfObjects->Add(hpTDphiBinsSAWLPMCTrue);

    hpTvsDphiOAMCTrue = 0;
    hpTvsDphiOAMCTrue= new TH2D("hpTvsDphiOAMCTrue","p_{T} vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",360,-3.15,3.15,nPtBins,PtBins);
    fListOfObjects->Add(hpTvsDphiOAMCTrue);

    hpTvsDphiSAMCTrue = 0;
    hpTvsDphiSAMCTrue = new TH2D("hpTvsDphiSAMCTrue","p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",nDphiBins,DphiBins,nPtBins,PtBins);
    fListOfObjects->Add(hpTvsDphiSAMCTrue);

    hMultvsDphiOAMCTrue = 0;
    hMultvsDphiOAMCTrue= new TH2D("hMultvsDphiOAMCTrue","Multiplicity vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);Multiplicity",360,-3.15,3.15,nMultBins,MultBins);
    fListOfObjects->Add(hMultvsDphiOAMCTrue);
    
    hMultvsDphiSAMCTrue = 0;
    hMultvsDphiSAMCTrue = new TH2D("hMultvsDphiSAMCTrue","Multiplicity vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity",nDphiBins,DphiBins,nMultBins,MultBins);
    fListOfObjects->Add(hMultvsDphiSAMCTrue);

    hMultvspTvsDphiMCTrue = 0;
    hMultvspTvsDphiMCTrue = new TH3D("hMultvspTvsDphiMCTrue","Multiplicity vs p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity;p_{T} (GeV/c)",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
    fListOfObjects->Add(hMultvspTvsDphiMCTrue);

    for(Int_t i=0;i<18;i++){
      //True---
      hpTDphiBinsOAMCTrue[i] = 0;
      hpTDphiBinsOAMCTrue[i] = new TH1D(Form("hpTdPhiBinsOAMCTrue%d",i),"Charged particles (opening #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBins,PtBins);
      fListOfObjects->Add(hpTDphiBinsOAMCTrue[i]);

      hpTDphiBinsSAMCTrue[i] = 0;
      hpTDphiBinsSAMCTrue[i] = new TH1D(Form("hpTdPhiBinsSAMCTrue%d",i),"Charged particles (sliding #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBins,PtBins);
      fListOfObjects->Add(hpTDphiBinsSAMCTrue[i]);

      hMultDphiBinsOAMCTrue[i] = 0;
      hMultDphiBinsOAMCTrue[i] = new TH1D(Form("hMultdPhiBinsOAMCTrue%d",i),"Charged particles (opening #Delta#phi);Multiplicity;count",nMultBins,MultBins);
      fListOfObjects->Add(hMultDphiBinsOAMCTrue[i]);

      hMultDphiBinsSAMCTrue[i] = 0;
      hMultDphiBinsSAMCTrue[i] = new TH1D(Form("hMultdPhiBinsSAMCTrue%d",i),"Charged particles (sliding #Delta#phi);Multiplicity;count",nMultBins,MultBins);
      fListOfObjects->Add(hMultDphiBinsSAMCTrue[i]);
    }

  }

      hpT = new TH1D("hpT","",nPtBins,PtBins);
      fListOfObjects->Add(hpT);

      hEta = new TH1D("hEta","; #eta^{leading};counts",20,-1,1);
      fListOfObjects->Add(hEta);

      hPhi = new TH1D("hPhi", ";#phi (rad); count", 64,0,2.0*TMath::Pi());
      fListOfObjects->Add(hPhi);

      hPtL = 0;
      hPtL = new TH1D("hPtL",";#it{p}_{T}^{leading} (GeV/#it{c});counts",nPtBins,PtBins);
      fListOfObjects->Add(hPtL);

      hEtaL = 0;
      hEtaL = new TH1D("hEtaL","; #eta^{leading};counts",20,-1,1);
      fListOfObjects->Add(hEtaL);

      hPhiL = 0;
      hPhiL = new TH1D("hPhiL","; #phi^{leading} (rad);counts",64,0,2.0*TMath::Pi());
      fListOfObjects->Add(hPhiL);

      hDphi = 0;
      hDphi = new TH1D("hDphi","",64,-2*TMath::Pi(),2*TMath::Pi());
      fListOfObjects->Add(hDphi);

      hpTvsDphiOA = 0;
      hpTvsDphiOA = new TH2D("hpTvsDphiOA","p_{T} vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",360,-3.15,3.15,nPtBins,PtBins);
      fListOfObjects->Add(hpTvsDphiOA);

      hpTvsDphiSA = 0;
      hpTvsDphiSA = new TH2D("hpTvsDphiSA","p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",nDphiBins,DphiBins,nPtBins,PtBins);
      fListOfObjects->Add(hpTvsDphiSA);

      hMultvsDphiOA = 0;
      hMultvsDphiOA= new TH2D("hMultvsDphiOA","Multiplicity vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);Multiplicity",360,-3.15,3.15,nMultBins,MultBins);
      fListOfObjects->Add(hMultvsDphiOA);
      
      hpTDphiBinsSAWLP = 0;
      hpTDphiBinsSAWLP = new TH1D("hpTDphiBinsSAWLP","Charged particles WLP (sliding #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBins,PtBins);
      fListOfObjects->Add(hpTDphiBinsSAWLP);
      
      hMultvsDphiSA = 0;
      hMultvsDphiSA = new TH2D("hMultvsDphiSA","Multiplicity vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity",nDphiBins,DphiBins,nMultBins,MultBins);
      fListOfObjects->Add(hMultvsDphiSA);

      hMultvspTvsDphi = 0;
      hMultvspTvsDphi = new TH3D("hMultvspTvsDphi","Multiplicity vs p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);Multiplicity;p_{T} (GeV/c)",nDphiBins,DphiBins,nMultBins,MultBins,nPtBins,PtBins);
      fListOfObjects->Add(hMultvspTvsDphi);

      for(Int_t i=0;i<18;i++){
	hpTDphiBinsOA[i] = 0;
	hpTDphiBinsOA[i] = new TH1D(Form("hpTdPhiBinsOA%d",i),"Charged particles (opening #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBins,PtBins);
	fListOfObjects->Add(hpTDphiBinsOA[i]);

	hpTDphiBinsSA[i] = 0;
	hpTDphiBinsSA[i] = new TH1D(Form("hpTdPhiBinsSA%d",i),"Charged particles (sliding #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBins,PtBins);
	fListOfObjects->Add(hpTDphiBinsSA[i]);

	hMultDphiBinsOA[i] = 0;
	hMultDphiBinsOA[i] = new TH1D(Form("hMultdPhiBinsOA%d",i),"Charged particles (opening #Delta#phi);Multiplicity;count",nMultBins,MultBins);
	fListOfObjects->Add(hMultDphiBinsOA[i]);

	hMultDphiBinsSA[i] = 0;
	hMultDphiBinsSA[i] = new TH1D(Form("hMultdPhiBinsSA%d",i),"Charged particles (sliding #Delta#phi);Multiplicity;count",nMultBins,MultBins);
	fListOfObjects->Add(hMultDphiBinsSA[i]);
      }

  fEventCuts.AddQAplotsToList(fListOfObjects);
  PostData(1, fListOfObjects);

}

//______________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::UserExec(Option_t *)
{

  // -----------------------------------------------------
  //			 InputEvent
  // -----------------------------------------------------

  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  // -----------------------------------------------------
  //			 E S D
  // -----------------------------------------------------
  fESD = dynamic_cast<AliESDEvent*>(event);
  
  if(!fESD){
    Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }
  
  // -----------------------------------------------------
  //			 MC
  // -----------------------------------------------------
  
  if (fAnalysisMC){
    fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!fMCEvent)  {
      Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }
    
    fMCStack = fMCEvent->Stack();
    if(!fMCStack){
      Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }

       
    AliHeader *header = fMCEvent->Header();
    if(!header) {AliDebug( AliLog::kError , "Header not avaible" ); return; }
  }
  
  
  /************ BEGINNING OF EVENT SELECTION *******************/
  // Get trigger decision
  fTriggeredEventMB = 0; //init
  if (!fAnalysisMC){  // for data check event selection as well
    if((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & ftrigBit )  fTriggeredEventMB = 1;  //event triggered as minimum bias
  }
  else {
    if (ftrigBit)  fTriggeredEventMB = 1;
  }

  Bool_t SPDvsClustersBG = kFALSE;

  AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
  if (!AnalysisUtils)
    {
      cout<<"------- No AnalysisUtils Object Found --------"<<AnalysisUtils<<endl;
      return;
    }
  else SPDvsClustersBG = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD); // We want NO BG

  Bool_t isNotPileUp = !fESD->IsPileupFromSPD(5,0.8);
  Bool_t IncompleteDAQ = fESD->IsIncompleteDAQ(); // we want is not incomplete DAQ

  // vertex
  const AliESDVertex * vertex    =    fESD->GetPrimaryVertex(); // tracks vertex, if not -> spd vertex, if not TPC vertex
  Bool_t isVtxGood = vertex->GetStatus() && selectVertex2015pp( fESD ,kTRUE,kFALSE,kTRUE); // requires Tracks and spd vertex, and Zconsistency of 5mm
  Double_t vertex_z = vertex->GetZ();
  Bool_t isVtxInZCut = (TMath::Abs(vertex_z) <= fVtxCut); // Zvtx in +- 10

  // Implement INEL>0
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  Bool_t isINEL0 = kFALSE;
  for (Int_t i = 0; i < mult->GetNumberOfTracklets(); ++i){ if (TMath::Abs(mult->GetEta(i)) < 1.) isINEL0 = kTRUE;}

  /********** IS PHYSICS SELECTION FLAG ****************************/
  fisPS = fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp;

  // recontructed INEL > 0 is PS + vtx + Zvtx inside +-10 ------
  isINEL0Rec = kFALSE;
  if ( isINEL0 && fisPS && isVtxGood && isVtxInZCut) isINEL0Rec = kTRUE;

  if (fisPS) fVtxBeforeCuts->Fill(vertex_z);              // VZ hack
  if (isINEL0Rec) fVtxAfterCuts->Fill(vertex_z);
  hSelEv->Fill(0); // all events
  if(fTriggeredEventMB) hSelEv->Fill(1); // triggered events
  if (fTriggeredEventMB && !IncompleteDAQ ) hSelEv->Fill(2); // trigger + IsIncompleteDAQ
  if (fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG) hSelEv->Fill(3); // trigger + IsIncompleteDAQ + BG rejection
  if(fPileUpRej)
    {
      if(fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp) //
	hSelEv->Fill(4); // trigger + IsIncompleteDAQ + BG rejection + PileUp
    }
  if (fisPS)hSelEv->Fill(5); //PS: trigger + IsIncompleteDAQ + BG rejection + PileUp + 1 Tracklet in eta +-1
  if (fisPS && isVtxGood) hSelEv->Fill(6); //PS + GetPrimaryVertex
  if (isINEL0Rec) hSelEv->Fill(7); //PS + GetPrimaryVertex + isVtxInZCut

  // -------------------------------------- multiplcity estimators section ------------------------------------------ //
  ftrackmult08 = -999;
  fv0mpercentile = -999;

  //ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);     //reference
  ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);     //tracklets
  //ftrackmult08 = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD); //Combined estimator

  hRefMult08->Fill(ftrackmult08);

  fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection"); // Esto es para 13 TeV
  if (!fMultSelection)
    cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
  else
    fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  hV0Mmult->Fill(fv0mpercentile);

  // cout<<"------- V0M mult ==  "<<fv0mpercentile<<"--------"<<endl;

  // ------------------------------------------ end of mult estimators -------------------------------------------------//
   if (fAnalysisMC){ // analysis for generated MC
    const AliVVertex *vertexMC = (AliVVertex*) fMCEvent->GetPrimaryVertex();
    fMCVzCut     = (TMath::Abs(vertexMC->GetZ()) <= fVtxCut);   // ZMCvtx in +- 10
    isINEL0True = isMCEventTrueINEL0(fMCEvent);

    // for trigger efficiency PS / INEL > 0 true
    if ( fisPS ){ 
      hINEL0MCTrig->Fill(0); // for trigger efficiency Trig/True
      hPS_MC->Fill(0);       // for missing vtx correction hPS_MC/hVtxPS_MC
    }
    if (fisPS && isVtxGood) hVtxPS_MC->Fill(0);
    if (isINEL0True) {
      hINEL0MCTrue->Fill(0);
      AnalyzeMC(); // analysis for MC
    }
   }
   else
     {
       if (isINEL0Rec) hINEL0->Fill(0);
       // Two histos for missing vtx correction
       if ( fisPS ) hPS->Fill(0);   
       if ( fisPS && isVtxGood ) hVtxPS->Fill(0);
     }

  if (isINEL0Rec) {AnalyzeESD(fESD);}
  //  cout<<"hello!!!"<<endl;

  // Post output data.
  PostData(1, fListOfObjects);

}
//________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::AnalyzeMC(){
  
  // selection on leading particle
  Double_t pt_leading    = 0;
  Double_t p_leading    = 0;
  Double_t eta_leading    = 0;
  Double_t phi_leading    = 0;
  Int_t    i_leading = 0;

   Int_t mult = 0;
  Int_t multOA[18];
  Int_t multSA[18];

  for ( int iT = 0 ; iT < fMCStack->GetNtrack(); iT++ ){ // loop over TRUE MC
    
    TParticle *mcParticle = fMCStack->Particle(iT);
    
    if (!mcParticle){
      cout<<"no mcParticle"<<endl;
      continue;
    }
    
    if(!(fMCStack->IsPhysicalPrimary(iT))) continue;
    
    Double_t eta = mcParticle->Eta();
    Double_t pt = mcParticle->Pt();
    Double_t phi = mcParticle->Phi();
    
    int partPDG = TMath::Abs(mcParticle->GetPdgCode());
    if ( TMath::Abs(eta) > fEtaCut ) continue;
    if ( pt < 0.15 ) continue;
    mult++;  
    
    if(pt>pt_leading){
      pt_leading      = pt;
      eta_leading     = eta;
      phi_leading     = phi;
      i_leading = iT;
    }
  
    hpTMCTrue->Fill(pt);
    hEtaMCTrue->Fill(eta);
    hPhiMCTrue->Fill(phi);

  }// end loop over tracks
  
  hPtLMCTrue->Fill(pt_leading);
  hEtaLMCTrue->Fill(eta_leading);
  hPhiLMCTrue->Fill(phi_leading);
  
  for(Int_t j=0;j<ndPhiBins;j++){
    multOA[j] =0;
    multSA[j]=0;
  }

  for ( int iT = 0 ; iT < fMCStack->GetNtrack(); iT++ ){ // loop over TRUE MC
    
    TParticle *mcParticle = fMCStack->Particle(iT);
    
    if (!mcParticle){
      cout<<"no mcParticle"<<endl;
      continue;
    }
    
    if(!fMCStack->IsPhysicalPrimary(iT)) continue;
    
    Double_t eta = mcParticle->Eta();
    Double_t pt = mcParticle->Pt();
    Double_t phi = mcParticle->Phi();
    
    int partPDG = TMath::Abs(mcParticle->GetPdgCode());
    if ( TMath::Abs(eta) > fEtaCut ) continue;
    if ( pt < 0.15 ) continue;
  
     Double_t DPhiOA = DeltaPhi( phi, phi_leading );
     Double_t DPhiSA = TMath::Abs(DPhiOA);
     
     hDphiMCTrue->Fill(DPhiOA);
     hpTvsDphiOAMCTrue->Fill(DPhiOA,pt);
     hpTvsDphiSAMCTrue->Fill(DPhiSA,pt);
     hMultvsDphiOAMCTrue->Fill(DPhiOA,mult);
     hMultvsDphiSAMCTrue->Fill(DPhiSA,mult);
     hMultvspTvsDphiMCTrue->Fill(DPhiSA,mult,pt);
     
     for(Int_t j=0;j<ndPhiBins;j++){
       if(TMath::Abs(DPhiOA)<= dPhiBinsOA[j]){
	 hpTDphiBinsOAMCTrue[j]->Fill(pt);
	 multOA[j]++;
       }
       if(DPhiSA > dPhiBinsSA[j] && DPhiSA <= dPhiBinsSA[j+1]){
	 hpTDphiBinsSAMCTrue[j]->Fill(pt);
	 multSA[j]++;
       }
     }
     if(DPhiSA >= dPhiBinsSA[0] && DPhiSA <= dPhiBinsSA[1]) hpTDphiBinsSAWLPMCTrue->Fill(pt);
  }// end loop over tracks
  for(Int_t k=0;k<ndPhiBins;k++){
    hMultDphiBinsOAMCTrue[k]->Fill(multOA[k]);
    hMultDphiBinsSAMCTrue[k]->Fill(multSA[k]);
  }
  
}
//________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::AnalyzeESD(AliESDEvent* fESD){
  
   // selection on leading particle
  Double_t pt_leading    = 0;
  Double_t p_leading    = 0;
  Double_t eta_leading    = 0;
  Double_t phi_leading    = 0;

  Int_t    i_leading = 0;

  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

    AliESDtrack* esdTrack = fESD->GetTrack(i);

    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t momentum = esdTrack->P();
    Double_t pt       = esdTrack->Pt();

    if(TMath::Abs(eta) > fEtaCut) continue;
    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack)) continue;
    if(pt<0.15) continue;
    Short_t ncl = esdTrack->GetTPCsignalN();
    if(ncl<fNcl) continue;

    if(pt>pt_leading){
      pt_leading      = pt;
      p_leading       = momentum;
      eta_leading     = eta;
      phi_leading     = phi;
      i_leading = i;
    }

    hpT->Fill(pt);
    hEta->Fill(eta);
    hPhi->Fill(phi);

  }// end loop over tracks

  if(pt_leading<0.15) return;

  hPtL->Fill(pt_leading);
  hEtaL->Fill(eta_leading);
  hPhiL->Fill(phi_leading);

  Int_t mult = 0;
  Int_t multOA[18];
  Int_t multSA[18];

  for(Int_t j=0;j<ndPhiBins;j++){
    multOA[j] =0;
    multSA[j]=0;
  }

  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

    AliESDtrack* esdTrack = fESD->GetTrack(i);
    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t pt       = esdTrack->Pt();

    if(TMath::Abs(eta) > fEtaCut) continue;
    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack)) continue;
    if(pt<0.15) continue;
    Short_t ncl = esdTrack->GetTPCsignalN();
    if(ncl<fNcl) continue;
    mult++;

    Double_t DPhiOA = DeltaPhi( phi, phi_leading );
    Double_t DPhiSA = TMath::Abs(DPhiOA);

    hDphi->Fill(DPhiOA);
    hpTvsDphiOA->Fill(DPhiOA,pt);
    hpTvsDphiSA->Fill(DPhiSA,pt);
    hMultvsDphiOA->Fill(DPhiOA,mult);
    hMultvsDphiSA->Fill(DPhiSA,mult);
    hMultvspTvsDphi->Fill(DPhiSA,mult,pt);

    for(Int_t j=0;j<ndPhiBins;j++){
      if(TMath::Abs(DPhiOA)<= dPhiBinsOA[j]){
	hpTDphiBinsOA[j]->Fill(pt);
	multOA[j]++;
      }
      if(DPhiSA > dPhiBinsSA[j] && DPhiSA <= dPhiBinsSA[j+1]){
	hpTDphiBinsSA[j]->Fill(pt);
	multSA[j]++;
      }
    }
    if(DPhiSA >= dPhiBinsSA[0] && DPhiSA <= dPhiBinsSA[1]) hpTDphiBinsSAWLP->Fill(pt);
  }// end loop over tracks
  for(Int_t k=0;k<ndPhiBins;k++){
    hMultDphiBinsOA[k]->Fill(multOA[k]);
    hMultDphiBinsSA[k]->Fill(multSA[k]);
  }
}

//______________________________________________________________________
Bool_t AliAnalysisTaskUeSpectraDphi::selectVertex2015pp(AliESDEvent *esd,
							Bool_t checkSPDres, //enable check on vtx resolution
							Bool_t requireSPDandTrk, //ask for both trk and SPD vertex
							Bool_t checkProximity) //apply cut on relative position of spd and trk verteces
{

  if (!esd) return kFALSE;

  const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
  const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
  Bool_t hasSPD = spdVertex->GetStatus();
  Bool_t hasTrk = trkVertex->GetStatus();

  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;

  //reject events if none between the SPD or track verteces are available
  //if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (!hasSPD) return kFALSE;
    //on demand check the spd vertex resolution and reject if not satisfied
    if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
  } else {
    if (hasSPD) {
      //if enabled check the spd vertex resolution and reject if not satisfied
      //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
      if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
      if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE;
    }
  }
  return kTRUE;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskUeSpectraDphi::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}

//____________________________________________________________
void AliAnalysisTaskUeSpectraDphi::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  // TPC
  esdTrackCuts->SetMinNCrossedRowsTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);

  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  // 7*(0.0015+0.0050/pt^1.1)
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");

  esdTrackCuts->SetMaxDCAToVertexZ(2);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetMaxChi2PerClusterITS(36);

  fTrackFilter->AddCuts(esdTrackCuts);

}

Bool_t AliAnalysisTaskUeSpectraDphi::isMCEventTrueINEL0(AliMCEvent* fMCEvent)
{
  Bool_t isINEL0 = kFALSE;
  for ( int iT = 0 ; iT < fMCStack->GetNtrack(); iT++ ) // loop over TRUE MC
  {
    TParticle *mcParticle = fMCStack->Particle(iT);

    if (!mcParticle)
    {
      cout<<"no mcParticle"<<endl;
      continue;
    }

    if(!fMCStack->IsPhysicalPrimary(iT))
      continue;

    if(!(mcParticle->Pt()>0.0))
      continue;

    Double_t eta = mcParticle->Eta();
    if ( TMath::Abs(eta) > 1.0 )
      continue;

    if ( !( TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3 ) )
      continue;

    isINEL0 = kTRUE;
    break;

  }

  return isINEL0;

}

Double_t AliAnalysisTaskUeSpectraDphi::DeltaPhi(Double_t phia, Double_t phib,
						Double_t rangeMin, Double_t rangeMax)
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();

  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi > pi)        dphi -= 2*pi;
  else if (dphi < -pi)  dphi += 2*pi;

  return dphi;
}
