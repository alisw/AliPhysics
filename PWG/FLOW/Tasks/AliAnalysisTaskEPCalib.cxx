#include <iostream>
#include <cstdlib>
#include <sys/time.h>
// ROOT classes
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TExMap.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TComplex.h"
#include "TSpline.h"
#include "TGrid.h"
// Alice analysis base class
#include "AliAnalysisTaskSE.h"
// Alice analysis additional classes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
// Alice AOD classes
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
// Alice MC classes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
// Alice "V" classes
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"
// Alice PID classes
#include "AliAODPid.h"
#include "AliAODpidUtil.h"
#include "AliPID.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
#include "AliAnalysisTaskEPCalib.h"
#include "AliCentrality.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEPCalib);

//---------------------------------------------------
AliAnalysisTaskEPCalib::AliAnalysisTaskEPCalib() : 
  AliAnalysisTaskSE(),
  fTPCEstOn(false),
  fTPCNUAWeight(false),
  fFillTPCQMean(false),
  fFillTPCShift(false),
  fTPCCalib (false),
  fFillVZEROQMean(false),
  fVZEROEstOn (false),
  fVZEROGainEq(false),
  fVZEROCalib(false),
  fQAV0(false),
  fDebug(0),
  fHarmonic(2),
  fTrigger(1),
  fFltbit(1),
  fNclsCut(70),
  fChi2Hg(4.0),
  fChi2Lo(0.1),
  fDcaCutz(3.2),
  fDcaCutxy(2.4),
  fPtMin(0.2),
  fPtMax(5.0),
  fCbinHg(8),
  fCbinLo(0),
  fPeriod("LHC10h"),
  fMultComp("pileupByEDSTPC128"),
  fCentCut(7.5),
  fRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCentBin(-999),
  fCent(-999),
  fCentSPD(-999),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fZvtxCut(10.0),
  fListNUA1(NULL),
  fListNUA2(NULL),
  fListNUA3(NULL),
  hNUAweightPlus(NULL),
  hNUAweightMinus(NULL),
  hCorrectNUAPos(NULL),
  hCorrectNUANeg(NULL),
  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),
  fMultCutPU(NULL),
  fOutputList(NULL),
  hEvtCount(NULL),
  hRunNumBin(NULL),
  hPt(NULL),
  hPDedx(NULL)
{
  // QA
  for (int i = 0; i < NRUNNUM; ++i) fRunNumList[i]="0";
  for (int i = 0; i < 2; ++i) hCent[i]=NULL;
  for (int i = 0; i < 2; ++i) hVz[i]=NULL;
  for (int i = 0; i < 8; ++i) hCentQA[i]=NULL;
  for (int i = 0; i < 2; ++i) hMultCentQA[i]=NULL;
  for (int i = 0; i < 6; ++i) hMultMultQA[i]=NULL;
  for (int i = 0; i < 2; ++i) hEta[i]=NULL;
  for (int i = 0; i < 2; ++i) hPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hEtaPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaXy[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaZ[i]=NULL;
  for (int i = 0; i < 2; ++i) hNhits[i]=NULL;

  // V0
  for (int iCent = 0; iCent < NCENTBINS; ++iCent){
    for (int i = 0; i < 3; ++i){
      hPsiVZERODirectGet[iCent][i]=NULL;
      hPsiV0Cor[iCent][i]=NULL;
    }     
  };
  for (int iRun = 0; iRun < NRUNNUM; ++iRun){
    hMultV0Fill[iRun]=NULL;
    hMultV0Read[iRun]=NULL;
    for (int i = 0; i < 3; ++i){
      pV0XMeanFill[iRun][i]=NULL;
      pV0YMeanFill[iRun][i]=NULL;
      hQx3mV0[iRun][i]=NULL;
      hQy3mV0[iRun][i]=NULL;
    }
  };
  for (int i = 0; i < 3; ++i){
    hQnCentCor[i]=NULL;
    hQxCentCor[i]=NULL;
    hQyCentCor[i]=NULL;
    hQxVtxCor[i]=NULL;
    hQyVtxCor[i]=NULL;
  }

  //TPC
  for (int iCent = 0; iCent < NCENTBINS; ++iCent){
    for (int i = 0; i < 12; ++i)hPsiTPC[iCent][i]=NULL;
  };
  for (int i = 0; i < 9; ++i){
    hQxTPCCent[i]=NULL;
    hQyTPCCent[i]=NULL;
    hQnTPCCent[i]=NULL;
  };
  for (int iRun = 0; iRun < NRUNNUM; ++iRun){
    for (int i = 0; i < 4; ++i){
      pTPCCosMeanFill[iRun][i]=NULL;
      pTPCSinMeanFill[iRun][i]=NULL;
      hTPCCosMeanRead[iRun][i]=NULL;
      hTPCSinMeanRead[iRun][i]=NULL;
    }
  };
  for (int i = 0; i < 3; ++i){
    pTPCShiftFillCoeffCos[i]=NULL;
    pTPCShiftFillCoeffSin[i]=NULL;
    hTPCShiftReadCoeffCos[i]=NULL;
    hTPCShiftReadCoeffSin[i]=NULL;
  }
}

//---------------------------------------------------
AliAnalysisTaskEPCalib::AliAnalysisTaskEPCalib(const char *name) : 
  AliAnalysisTaskSE(name),
  fTPCEstOn(false),
  fTPCNUAWeight(false),
  fFillTPCQMean(false),
  fFillTPCShift(false),
  fTPCCalib (false),
  fFillVZEROQMean(false),
  fVZEROEstOn (false),
  fVZEROGainEq(false),
  fVZEROCalib(false),
  fQAV0(false),
  fDebug(0),
  fHarmonic(2),
  fTrigger(1),
  fFltbit(1),
  fNclsCut(70),
  fChi2Hg(4.0),
  fChi2Lo(0.1),
  fDcaCutz(3.2),
  fDcaCutxy(2.4),
  fPtMin(0.2),
  fPtMax(5.0),
  fCbinHg(8),
  fCbinLo(0),
  fPeriod("LHC10h"),
  fMultComp("pileupByEDSTPC128"),
  fCentCut(7.5),
  fRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCentBin(-999),
  fCent(-999),
  fCentSPD(-999),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fZvtxCut(10.0),
  fListNUA1(NULL),
  fListNUA2(NULL),
  fListNUA3(NULL),
  hNUAweightPlus(NULL),
  hNUAweightMinus(NULL),
  hCorrectNUAPos(NULL),
  hCorrectNUANeg(NULL),
  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),
  fMultCutPU(NULL),
  fOutputList(NULL),
  hEvtCount(NULL),
  hRunNumBin(NULL),
  hPt(NULL),
  hPDedx(NULL)
{
  // QA
  for (int i = 0; i < NRUNNUM; ++i) fRunNumList[i]="0";
  for (int i = 0; i < 2; ++i) hCent[i]=NULL;
  for (int i = 0; i < 2; ++i) hVz[i]=NULL;
  for (int i = 0; i < 8; ++i) hCentQA[i]=NULL;
  for (int i = 0; i < 2; ++i) hMultCentQA[i]=NULL;
  for (int i = 0; i < 6; ++i) hMultMultQA[i]=NULL;
  for (int i = 0; i < 2; ++i) hEta[i]=NULL;
  for (int i = 0; i < 2; ++i) hPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hEtaPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaXy[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaZ[i]=NULL;
  for (int i = 0; i < 2; ++i) hNhits[i]=NULL;

  // V0
  for (int iCent = 0; iCent < NCENTBINS; ++iCent){
    for (int i = 0; i < 3; ++i){
      hPsiVZERODirectGet[iCent][i]=NULL;
      hPsiV0Cor[iCent][i]=NULL;
    }     
  };
  for (int iRun = 0; iRun < NRUNNUM; ++iRun){
    hMultV0Fill[iRun]=NULL;
    hMultV0Read[iRun]=NULL;
    for (int i = 0; i < 3; ++i){
      pV0XMeanFill[iRun][i]=NULL;
      pV0YMeanFill[iRun][i]=NULL;
      hQx3mV0[iRun][i]=NULL;
      hQy3mV0[iRun][i]=NULL;
    }
  };
  for (int i = 0; i < 3; ++i){
    hQnCentCor[i]=NULL;
    hQxCentCor[i]=NULL;
    hQyCentCor[i]=NULL;
    hQxVtxCor[i]=NULL;
    hQyVtxCor[i]=NULL;
  }

  //TPC
  for (int iCent = 0; iCent < NCENTBINS; ++iCent){
    for (int i = 0; i < 12; ++i)hPsiTPC[iCent][i]=NULL;
  };
  for (int i = 0; i < 9; ++i){
    hQxTPCCent[i]=NULL;
    hQyTPCCent[i]=NULL;
    hQnTPCCent[i]=NULL;
  };
  for (int iRun = 0; iRun < NRUNNUM; ++iRun){
    for (int i = 0; i < 4; ++i){
      pTPCCosMeanFill[iRun][i]=NULL;
      pTPCSinMeanFill[iRun][i]=NULL;
      hTPCCosMeanRead[iRun][i]=NULL;
      hTPCSinMeanRead[iRun][i]=NULL;
    }
  };
  for (int i = 0; i < 3; ++i){
    pTPCShiftFillCoeffCos[i]=NULL;
    pTPCShiftFillCoeffSin[i]=NULL;
    hTPCShiftReadCoeffCos[i]=NULL;
    hTPCShiftReadCoeffSin[i]=NULL;
  }

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//---------------------------------------------------
AliAnalysisTaskEPCalib::~AliAnalysisTaskEPCalib()
{
}

//---------------------------------------------------
void AliAnalysisTaskEPCalib::Terminate(Option_t *) 
{
}

//---------------------------------------------------
void AliAnalysisTaskEPCalib::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetName(GetName());
  fOutputList->SetOwner(kTRUE);

  // Run Info 10h
  if (fPeriod.EqualTo("LHC10h")){
    TString RunList[NRUNNUM]={
      "139510","139507","139505","139503","139465","139438","139437","139360","139329","139328","139314","139310",
      "139309","139173","139107","139105","139038","139037","139036","139029","139028","138872","138871","138870",
      "138837","138732","138730","138666","138662","138653","138652","138638","138624","138621","138583","138582",
      "138579","138578","138534","138469","138442","138439","138438","138396","138364","138275","138225","138201",
      "138197","138192","138190","137848","137844","137752","137751","137724","137722","137718","137704","137693",
      "137692","137691","137686","137685","137639","137638","137608","137595","137549","137546","137544","137541",
      "137539","137531","137530","137443","137441","137440","137439","137434","137432","137431","137430","137243",
      "137236","137235","137232","137231","137230","137162","137161"};
      hRunNumBin = new TH1I("runNumBin","",100,0,100);
      for (int i=0; i<NRUNNUM; ++i) {    
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,RunList[i].Data());
        fRunNumList[i] = RunList[i];
      }
      fOutputList->Add(hRunNumBin);
  };
  // Run Info 11h
  if (fPeriod.EqualTo("LHC11h")){
    TString RunList[NRUNNUM]={"170387","170040","170268","170228","170207","169838","170159","170204","170311","170084",
      "169835","170088","170593","170203","170270","169846","170163","170388","170155","170083",
      "170572","169837","169855","170306","170269","170089","170309","170091","170081","170230",
      "170085","170315","170027","170193","170312","170313","170308","169858","169859"};
      hRunNumBin = new TH1I("runNumBin","",100,0,100);
      for (int i=0; i<NRUNNUM; ++i) {    
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,RunList[i].Data());
        fRunNumList[i] = RunList[i];
      }
      fOutputList->Add(hRunNumBin);
  };
  // Run Info 15o
  if (fPeriod.EqualTo("LHC15o") ){
    TString RunList[NRUNNUM]={
         "246994","246991","246989","246984","246982","246948","246945","246928","246871","246870","246867","246865", 
         "246864","246859","246858","246851","246847","246846","246845","246844","246810","246809","246808","246807", 
         "246805","246804","246766","246765","246763","246760","246759","246758","246757","246751","246750","246434", 
         "246431","246424","246392","246391","246276","246275","246272","246271","246225","246222","246217","246185", 
         "246182","246181","246180","246178","246153","246152","246151","246148","246115","246113","246089","246087", 
         "246053","246052","246049","246048","246042","246037","246036","246012","246003","246001","245963","245954", 
         "245952","245949","245923","245833","245831","245829","245793","245785","245775","245766","245759","245752", 
         "245731","245729","245705","245702","245692","245683","245554","245545","245544","245543","245542","245540", 
         "245535","245507","245505","245504","245501","245497","245496","245454","245453","245450","245446","245441", 
         "245411","245410","245409","245407","245401","245397","245396","245353","245349","245347","245346","245345", 
         "245343","245259","245233","245232","245231","245152","245151","245146","245145","245068","245066","245064", 
         "244983","244982","244980","244975","244918","244917"};
      hRunNumBin = new TH1I("runNumBin","",150,0,150);
      for (int i=0; i<NRUNNUM; ++i){    
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,RunList[i].Data());
        fRunNumList[i] = RunList[i];
      }
      fOutputList->Add(hRunNumBin);
  };

  //------------------
  // QA
  //------------------
  // event-wise
  hEvtCount = new TH1D("hEvtCount", "Event Count", 25, 1, 26);
  hEvtCount->GetXaxis()->SetBinLabel(1,"all");
  hEvtCount->GetXaxis()->SetBinLabel(2,"readin");
  hEvtCount->GetXaxis()->SetBinLabel(3,"evt");
  hEvtCount->GetXaxis()->SetBinLabel(4,"runNum");
  hEvtCount->GetXaxis()->SetBinLabel(5,"vertex");
  hEvtCount->GetXaxis()->SetBinLabel(6,"centrality");
  hEvtCount->GetXaxis()->SetBinLabel(7,"pileup");
  hEvtCount->GetXaxis()->SetBinLabel(8,"V0FillGE");
  hEvtCount->GetXaxis()->SetBinLabel(9,"V0FillRC");
  hEvtCount->GetXaxis()->SetBinLabel(10,"V0Calib");
  hEvtCount->GetXaxis()->SetBinLabel(11,"TPCFillMean");
  hEvtCount->GetXaxis()->SetBinLabel(12,"TPCFillShift");
  hEvtCount->GetXaxis()->SetBinLabel(13,"TPCCalib");
  hEvtCount->GetXaxis()->SetBinLabel(20,"manager");
  hEvtCount->GetXaxis()->SetBinLabel(21,"handler");
  hEvtCount->GetXaxis()->SetBinLabel(22,"fAOD");
  hEvtCount->GetXaxis()->SetBinLabel(23,"fPID");
  hEvtCount->GetXaxis()->SetBinLabel(24,"fUtils");
  hEvtCount->GetXaxis()->SetBinLabel(25,"fMultSel");
  fOutputList->Add(hEvtCount);

  hCent[0] = new TH1D("hCentBefCut", "", 100, 0., 100.);
  hCent[1] = new TH1D("hCentAftCut", "", 100, 0., 100.);
  fOutputList->Add(hCent[0]);
  fOutputList->Add(hCent[1]);

  hVz[0] = new TH1D("hVzBeforeCut", "", 200, -50., 50.);
  hVz[1] = new TH1D("hVzAfterCut",  "", 200, -50., 50.);
  fOutputList->Add(hVz[0]);
  fOutputList->Add(hVz[1]);

  hCentQA[0] = new TH2D("hCentQAV0MSPD1BefCut", ";centV0M;centSPD1", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[0]);
  hCentQA[1] = new TH2D("hCentQAV0MSPD1AftCut", ";centV0M;centSPD1", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[1]);

  hCentQA[2] = new TH2D("hCentQAV0MTRKBefCut", ";centV0M;centTRK", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[2]);
  hCentQA[3] = new TH2D("hCentQAV0MTRKAftCut", ";centV0M;centTRK", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[3]);

  hCentQA[4] = new TH2D("hCentQAV0MSPD0BefCut", ";centV0M;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[4]);
  hCentQA[5] = new TH2D("hCentQAV0MSPD0AftCut", ";centV0M;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[5]);

  hCentQA[6] = new TH2D("hCentQASPD1SPD0BefCut", ";centSPD1;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[6]);
  hCentQA[7] = new TH2D("hCentQASPD1SPD0AftCut", ";centSPD1;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[7]);

  hMultCentQA[0] = new TH2D("hMultCentQABefCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  hMultCentQA[1] = new TH2D("hMultCentQAAftCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  fOutputList->Add(hMultCentQA[0]);
  fOutputList->Add(hMultCentQA[1]);

  if (fMultComp.EqualTo("pileupByEDSTPC128") || fMultComp.EqualTo("pileupByGlobalTPC1")){
    hMultMultQA[0] = new TH2D("hMultMultQAmTPCmESDPileupBefCut", "befCut;multTPC;multESD", 50, 0, 5000, 160, 0, 16000);
    hMultMultQA[1] = new TH2D("hMultMultQAmClobalmTPCEFPileupBefCut", "befCut;multGlobal;multTPCFE", 50, 0, 5000, 50, 0, 5000);
    hMultMultQA[2] = new TH2D("hMultMultQAmFB32mTrkTOFBefCut", "befCut;multTrkTOF;nTrk", 201, 0, 20000, 201, 0, 20000);
    hMultMultQA[3] = new TH2D("hMultMultQAmTPCmESDPileupAftCut", "aftCut;multTPC;multESD", 50, 0, 5000, 160, 0, 16000);
    hMultMultQA[4] = new TH2D("hMultMultQAmClobalmTPCEFPileupAftCut", "aftCut;multGlobal;multTPCFE", 50, 0, 5000, 50, 0, 5000);
    hMultMultQA[5] = new TH2D("hMultMultQAmFB32mTrkTOFAftCut", "aftCut;multTrkTOF;nTrk", 201, 0, 20000, 201, 0, 20000);
    fOutputList->Add(hMultMultQA[0]);
    fOutputList->Add(hMultMultQA[1]);
    fOutputList->Add(hMultMultQA[2]);
    fOutputList->Add(hMultMultQA[3]);
    fOutputList->Add(hMultMultQA[4]);
    fOutputList->Add(hMultMultQA[5]);
  }

  // track-wise
  hPt = new TH1D("hPt", "", 200, 0., 20.);
  fOutputList->Add(hPt);
  hEta[0] = new TH1D("hEtaBeforeCut", "", 200, -10., 10.);
  hEta[1] = new TH1D("hEtaAfterCut",  "", 200, -10., 10.);
  fOutputList->Add(hEta[0]);
  fOutputList->Add(hEta[1]);
  hPhi[0] = new TH1D("hPhiBeforeCor", "", 50, 0, 6.283185);
  hPhi[1] = new TH1D("hPhiAfterCor", "", 50, 0, 6.283185);
  fOutputList->Add(hPhi[0]);
  fOutputList->Add(hPhi[1]);
  hEtaPhi[0] = new TH2D("hEtaPhiBeforeCor", "", 50, 0, 6.283185, 16,-0.8,0.8);
  hEtaPhi[1] = new TH2D("hEtaPhiAfterCor", "", 50, 0, 6.283185, 16,-0.8,0.8);
  fOutputList->Add(hEtaPhi[0]);
  fOutputList->Add(hEtaPhi[1]);
  hNhits[0] = new TH1D("hNhitsBeforeCut", "", 200, 0., 200.);
  hNhits[1] = new TH1D("hNhitsAfterCut",  "", 200, 0., 200.);
  fOutputList->Add(hNhits[0]);
  fOutputList->Add(hNhits[1]);
  hPDedx = new TH2D("hPDedx", "", 400, -10., 10., 400, 0, 1000);
  fOutputList->Add(hPDedx);
  if (fPeriod.EqualTo("LHC11h") || fPeriod.EqualTo("LHC10h")){
    hDcaXy[0] = new TH1D("hDcaXyBeforeCut", "", 100, 0., 10.);
    hDcaXy[1] = new TH1D("hDcaXyAfterCut",  "", 100, 0., 10.);
    fOutputList->Add(hDcaXy[0]);
    fOutputList->Add(hDcaXy[1]);
    hDcaZ[0] = new TH1D("hDcaZBeforeCut", "", 100, 0., 10.);
    hDcaZ[1] = new TH1D("hDcaZAfterCut",  "", 100, 0., 10.);
    fOutputList->Add(hDcaZ[0]);
    fOutputList->Add(hDcaZ[1]);
  }

  // V0 Calib
  if (fVZEROEstOn) {

    // Prepared for function of "GetEventPlane"
    for (int iCent=0; iCent<NCENTBINS; ++iCent) {
      for (int i=0; i<3; ++i) {
        hPsiVZERODirectGet[iCent][i] = new TH1D(Form("psiVZERO_GetEP_cent%i_%i",iCent,i), "",360,0,2*3.142);
        if (fVZEROCalib) fOutputList->Add(hPsiVZERODirectGet[iCent][i]);
      }
    }
      
    // Gain Eq
    if (fVZEROGainEq){
      AliOADBContainer* contMult = new AliOADBContainer("hMultV0BefCorPfpx");
      for (int iRun = 0; iRun < NRUNNUM; ++iRun ){
        hMultV0Fill[iRun] = new TH1D(Form("hMultV0Run%i", (int)fRunNumList[iRun].Atoi()), "", 64, 0, 64);
        contMult->AppendObject(hMultV0Fill[iRun], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());      
      } 
      contMult->WriteToFile("calibV015oP2.root");
    };

    // Fill Mean for recenter
    if (fFillVZEROQMean) {
      AliOADBContainer* contQx3mm = new AliOADBContainer("fqxm3m");
      AliOADBContainer* contQy3mm = new AliOADBContainer("fqym3m");
      AliOADBContainer* contQx3cm = new AliOADBContainer("fqxc3m");
      AliOADBContainer* contQy3cm = new AliOADBContainer("fqyc3m");
      AliOADBContainer* contQx3am = new AliOADBContainer("fqxa3m");
      AliOADBContainer* contQy3am = new AliOADBContainer("fqya3m");
      for (int iRun = 0; iRun < NRUNNUM; ++iRun ){
        pV0XMeanFill[iRun][0] = new TProfile(Form("hV0QxMeanMRun%i", (int)fRunNumList[iRun].Atoi()),"", 90,0,90);
        pV0YMeanFill[iRun][0] = new TProfile(Form("hV0QyMeanMRun%i", (int)fRunNumList[iRun].Atoi()),"", 90,0,90);
        pV0XMeanFill[iRun][1] = new TProfile(Form("hV0QxMeanCRun%i", (int)fRunNumList[iRun].Atoi()),"", 90,0,90);
        pV0YMeanFill[iRun][1] = new TProfile(Form("hV0QyMeanCRun%i", (int)fRunNumList[iRun].Atoi()),"", 90,0,90);
        pV0XMeanFill[iRun][2] = new TProfile(Form("hV0QxMeanARun%i", (int)fRunNumList[iRun].Atoi()),"", 90,0,90);
        pV0YMeanFill[iRun][2] = new TProfile(Form("hV0QyMeanARun%i", (int)fRunNumList[iRun].Atoi()),"", 90,0,90);
        contQx3mm->AppendObject(pV0XMeanFill[iRun][0], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());  
        contQy3mm->AppendObject(pV0YMeanFill[iRun][0], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());  
        contQx3cm->AppendObject(pV0XMeanFill[iRun][1], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());  
        contQy3cm->AppendObject(pV0YMeanFill[iRun][1], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());  
        contQx3am->AppendObject(pV0XMeanFill[iRun][2], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());  
        contQy3am->AppendObject(pV0YMeanFill[iRun][2], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());  
      }
      contQx3mm->WriteToFile("calibV015oP2.root");
      contQx3mm->WriteToFile("calibV015oP2.root");
      contQx3cm->WriteToFile("calibV015oP2.root");
      contQy3cm->WriteToFile("calibV015oP2.root");
      contQy3am->WriteToFile("calibV015oP2.root");
      contQy3am->WriteToFile("calibV015oP2.root");
    };

    // Read Mean & Calib
    if (fFillVZEROQMean || fVZEROCalib){
      if (!gGrid) {
          TGrid::Connect("alien://");
      }
      
      TFile* fileV0Calib15o = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/calibV015oP2.root");
      if(!fileV0Calib15o){
          printf("OADB V0 calibration file cannot be opened\n");
          return;
      }

      // TFile* fileV0Calib15o = TFile::Open("~/alice/LEGOTrainTest/refData/refLHC15o/calibV015oP2.root");
      // if(!fileV0Calib15o){
      //     printf("OADB V0 calibration file cannot be opened\n");
      //     return;
      // }

      // Mult
      AliOADBContainer* contMult = (AliOADBContainer*) fileV0Calib15o->Get("hMultV0BefCorPfpx");
      if(!contMult){
          printf("OADB object hMultV0BefCorr is not available in the file\n");
          return;
      }
      for (int iRun = 0; iRun < NRUNNUM; ++iRun){
        if(!(contMult->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hMultV0Read[iRun] = ((TH1D*) contMult->GetObject(fRunNumList[iRun].Atoi()));
      }
  
      if (fVZEROCalib){
        // V0C Qx Mean 
       AliOADBContainer* contQx3cm = (AliOADBContainer*) fileV0Calib15o->Get("fqxc3m");
       if(!contQx3cm){
           printf("OADB object fqxc2m is not available in the file\n");
           return;
       }
       
       // V0C Qy Mean 
       AliOADBContainer* contQy3cm = (AliOADBContainer*) fileV0Calib15o->Get("fqyc3m");
       if(!contQy3cm){
           printf("OADB object fqyc2m is not available in the file\n");
           return;
       }
       
       // V0A Qx Mean 
       AliOADBContainer* contQx3am = (AliOADBContainer*) fileV0Calib15o->Get("fqxa3m");
       if(!contQx3am){
           printf("OADB object fqxa2m is not available in the file\n");
           return;
       }
       
       // V0A Qy Mean 
       AliOADBContainer* contQy3am = (AliOADBContainer*) fileV0Calib15o->Get("fqya3m");
       if(!contQy3am){
           printf("OADB object fqya2m is not available in the file\n");
           return;
       }
  
       for (int iRun = 0; iRun < NRUNNUM; ++iRun){
         // V0C Qx Mean 
         if(!(contQx3cm->GetObject(fRunNumList[iRun].Atoi() ) ) ){
             printf("OADB object fqxc3m is not available for run %i\n", fRunNumList[iRun].Atoi() );
             return;
         }
         hQx3mV0[iRun][1] = ((TH1D*) contQx3cm->GetObject(fRunNumList[iRun].Atoi() ) );
         
         // V0C Qy Mean 
         if(!(contQy3cm->GetObject(fRunNumList[iRun].Atoi() ) ) ){
             printf("OADB object fqyc2m is not available for run %i\n", fRunNumList[iRun].Atoi() );
             return;
         }
         hQy3mV0[iRun][1] = ((TH1D*) contQy3cm->GetObject(fRunNumList[iRun].Atoi() ) );
         
         // V0A Qx Mean
         if(!(contQx3am->GetObject(fRunNumList[iRun].Atoi() ) ) ){
             printf("OADB object fqxa2m is not available for run %i\n", fRunNumList[iRun].Atoi() );
             return;
         }
         hQx3mV0[iRun][2] = ((TH1D*) contQx3am->GetObject(fRunNumList[iRun].Atoi() ) );
         
         // V0A Qy Mean
         if(!(contQy3am->GetObject(fRunNumList[iRun].Atoi() ) ) ){
             printf("OADB object fqya2m is not available for run %i\n", fRunNumList[iRun].Atoi() );
             return;
         }
         hQy3mV0[iRun][2] = ((TH1D*) contQy3am->GetObject(fRunNumList[iRun].Atoi() ) );
       }
      }
    };

    for (int iCent = 0; iCent < 10; ++iCent){
      hPsiV0Cor[iCent][0] = new TH1D(Form("hPsiCorV0MCent%i", iCent),"",180, 0, 3.142);
      hPsiV0Cor[iCent][1] = new TH1D(Form("hPsiCorV0CCent%i", iCent),"",180, 0, 3.142);
      hPsiV0Cor[iCent][2] = new TH1D(Form("hPsiCorV0ACent%i", iCent),"",180, 0, 3.142);
      if (fVZEROCalib){
        fOutputList->Add(hPsiV0Cor[iCent][0]);
        fOutputList->Add(hPsiV0Cor[iCent][1]);
        fOutputList->Add(hPsiV0Cor[iCent][2]);
      }
    }

    hQnCentCor[0] = new TH2D("hqnV0M","",100, 0, 100, 72, 0, 12); 
    hQnCentCor[1] = new TH2D("hqnV0C","",100, 0, 100, 72, 0, 12); 
    hQnCentCor[2] = new TH2D("hqnV0A","",100, 0, 100, 72, 0, 12); 
    if (fVZEROCalib){
      fOutputList->Add(hQnCentCor[0]);
      fOutputList->Add(hQnCentCor[1]);
      fOutputList->Add(hQnCentCor[2]);    
    }
  
    if (fQAV0){
      hQxCentCor[0] = new TH2D("hQxCorV0MvsCentSPD","",100, 0, 100, 32, -600, 600);
      hQxCentCor[1] = new TH2D("hQxCorV0CvsCentSPD","",100, 0, 100, 32, -600, 600);
      hQxCentCor[2] = new TH2D("hQxCorV0AvsCentSPD","",100, 0, 100, 32, -600, 600);
      hQyCentCor[0] = new TH2D("hQyCorV0MvsCentSPD","",100, 0, 100, 32, -600, 600);
      hQyCentCor[1] = new TH2D("hQyCorV0CvsCentSPD","",100, 0, 100, 32, -600, 600);
      hQyCentCor[2] = new TH2D("hQyCorV0AvsCentSPD","",100, 0, 100, 32, -600, 600);
      hQxVtxCor[0] = new TH2D("hQxCorV0MvsVz","",20, -10, 10, 32, -600, 600);
      hQxVtxCor[1] = new TH2D("hQxCorV0CvsVz","",20, -10, 10, 32, -600, 600);
      hQxVtxCor[2] = new TH2D("hQxCorV0AvsVz","",20, -10, 10, 32, -600, 600);
      hQyVtxCor[0] = new TH2D("hQyCorV0MvsVz","",20, -10, 10, 32, -600, 600);
      hQyVtxCor[1] = new TH2D("hQyCorV0CvsVz","",20, -10, 10, 32, -600, 600);
      hQyVtxCor[2] = new TH2D("hQyCorV0AvsVz","",20, -10, 10, 32, -600, 600);
      if (fVZEROCalib){
        for (int i = 0; i < 3; ++i) fOutputList->Add(hQxCentCor[i]); 
        for (int i = 0; i < 3; ++i) fOutputList->Add(hQyCentCor[i]); 
        for (int i = 0; i < 3; ++i) fOutputList->Add(hQxVtxCor[i]); 
        for (int i = 0; i < 3; ++i) fOutputList->Add(hQyVtxCor[i]); 
      }
    }
  };

  // TPC Calib
  if (fTPCEstOn) {

    if (fTPCNUAWeight) {
      if (fPeriod.EqualTo("LHC10h")){
        TFile* fileNUACor = TFile::Open("","READ");
        if (!fileNUACor) {
          AliError(Form("%s: Wrong TPCNUAEntry file.", GetName()));
          return;
        }

        fListNUA1 = (TList*)(fileNUACor->Get("listNUA_139510to138653"));
        fListNUA2 = (TList*)(fileNUACor->Get("listNUA_138652to137693"));
        fListNUA2 = (TList*)(fileNUACor->Get("listNUA_138652to137693"));
      }

      if (fPeriod.EqualTo("LHC15o")){
        TFile* fileNUACor = TFile::Open("","READ");
        if (!fileNUACor) {
          AliError(Form("%s: Wrong TPCNUAEntry file.", GetName()));
          return;
        }
        TDirectoryFile* directoryfileNUA = (TDirectoryFile*)fileNUACor->Get("ZDCgains");
        fListNUA1 = (TList*)(directoryfileNUA->Get("fNUA_ChPosChNeg"));
      }

    };

    if (fFillTPCQMean) {
      AliOADBContainer* contQxEtaPVzP= new AliOADBContainer("fqxetapvzp");
      AliOADBContainer* contQyEtaPVzP= new AliOADBContainer("fqyetapvzp");
      AliOADBContainer* contQxEtaPVzN= new AliOADBContainer("fqxetapvzn");
      AliOADBContainer* contQyEtaPVzN= new AliOADBContainer("fqyetapvzn");
      AliOADBContainer* contQxEtaNVzP= new AliOADBContainer("fqxetanvzp");
      AliOADBContainer* contQyEtaNVzP= new AliOADBContainer("fqyetanvzp");
      AliOADBContainer* contQxEtaNVzN= new AliOADBContainer("fqxetanvzn");
      AliOADBContainer* contQyEtaNVzN= new AliOADBContainer("fqyetanvzn");
      for (int iRun = 0; iRun <NRUNNUM; ++iRun){
        pTPCCosMeanFill[iRun][0] = new TProfile(Form("hTPCQxMeanEtaPVzPRun%i", (int)fRunNumList[iRun].Atoi()),"",10,0,10);
        pTPCSinMeanFill[iRun][0] = new TProfile(Form("hTPCQyMeanEtaPVzPRun%i", (int)fRunNumList[iRun].Atoi()),"",10,0,10);
        pTPCCosMeanFill[iRun][1] = new TProfile(Form("hTPCQxMeanEtaPVzNRun%i", (int)fRunNumList[iRun].Atoi()),"",10,0,10);
        pTPCSinMeanFill[iRun][1] = new TProfile(Form("hTPCQyMeanEtaPVzNRun%i", (int)fRunNumList[iRun].Atoi()),"",10,0,10);
        pTPCCosMeanFill[iRun][2] = new TProfile(Form("hTPCQxMeanEtaNVzPRun%i", (int)fRunNumList[iRun].Atoi()),"",10,0,10);
        pTPCSinMeanFill[iRun][2] = new TProfile(Form("hTPCQyMeanEtaNVzPRun%i", (int)fRunNumList[iRun].Atoi()),"",10,0,10);
        pTPCCosMeanFill[iRun][3] = new TProfile(Form("hTPCQxMeanEtaNVzNRun%i", (int)fRunNumList[iRun].Atoi()),"",10,0,10);
        pTPCSinMeanFill[iRun][3] = new TProfile(Form("hTPCQyMeanEtaNVzNRun%i", (int)fRunNumList[iRun].Atoi()),"",10,0,10);
        contQxEtaPVzP->AppendObject(pTPCCosMeanFill[iRun][0], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());  
        contQyEtaPVzP->AppendObject(pTPCSinMeanFill[iRun][0], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());
        contQxEtaPVzN->AppendObject(pTPCCosMeanFill[iRun][1], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());
        contQyEtaPVzN->AppendObject(pTPCSinMeanFill[iRun][1], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());
        contQxEtaNVzP->AppendObject(pTPCCosMeanFill[iRun][2], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());
        contQyEtaNVzP->AppendObject(pTPCSinMeanFill[iRun][2], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());
        contQxEtaNVzN->AppendObject(pTPCCosMeanFill[iRun][3], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());
        contQyEtaNVzN->AppendObject(pTPCSinMeanFill[iRun][3], fRunNumList[iRun].Atoi(), fRunNumList[iRun].Atoi());
        contQxEtaPVzP->WriteToFile("calibTPC15oP2.root");
        contQyEtaPVzP->WriteToFile("calibTPC15oP2.root");
        contQxEtaPVzN->WriteToFile("calibTPC15oP2.root");
        contQyEtaPVzN->WriteToFile("calibTPC15oP2.root");
        contQxEtaNVzP->WriteToFile("calibTPC15oP2.root");
        contQyEtaNVzP->WriteToFile("calibTPC15oP2.root");
        contQxEtaNVzN->WriteToFile("calibTPC15oP2.root");
        contQyEtaNVzN->WriteToFile("calibTPC15oP2.root");
      }
    }

    if (fFillTPCShift || fTPCCalib) {
      TFile* fileTPCCalib15o = TFile::Open("fileTPCCalib15o.root","READ");
      if (!fileTPCCalib15o) {
          printf("OADB TPC calibration file cannot be opened\n");
          return;
      }

      // Read Qx Qy Mean
      // Eta>0 Vz>0
      AliOADBContainer* contQxEtaPVzP = (AliOADBContainer*) fileTPCCalib15o->Get("fqxetapvzp");
      if(!contQxEtaPVzP){
          printf("OADB object fqxetapvzp is not available in the file\n");
          return;
      }
      for (int iRun = 0; iRun < NRUNNUM; ++iRun){
        if(!(contQxEtaPVzP->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object fqxetapvzp is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hTPCCosMeanRead[iRun][0] = ((TH1D*) contQxEtaPVzP->GetObject(fRunNumList[iRun].Atoi()));
      }

      AliOADBContainer* contQyEtaPVzP = (AliOADBContainer*) fileTPCCalib15o->Get("fqyetapvzp");
      if(!contQyEtaPVzP){
          printf("OADB object fqyetapvzp is not available in the file\n");
          return;
      }
      for (int iRun = 0; iRun < NRUNNUM; ++iRun){
        if(!(contQyEtaPVzP->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object fqyetapvzp is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hTPCSinMeanRead[iRun][0] = ((TH1D*) contQyEtaPVzP->GetObject(fRunNumList[iRun].Atoi()));
      }

      // Eta>0 Vz<0
      AliOADBContainer* contQxEtaPVzN = (AliOADBContainer*) fileTPCCalib15o->Get("fqxetapvzn");
      if(!contQxEtaPVzN){
          printf("OADB object fqxetapvzn is not available in the file\n");
          return;
      }
      for (int iRun = 0; iRun < NRUNNUM; ++iRun){
        if(!(contQxEtaPVzN->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object fqxetapvzn is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hTPCCosMeanRead[iRun][1] = ((TH1D*) contQxEtaPVzN->GetObject(fRunNumList[iRun].Atoi()));
      }
      
      AliOADBContainer* contQyEtaPVzN = (AliOADBContainer*) fileTPCCalib15o->Get("fqyetapvzn");
      if(!contQyEtaPVzN){
          printf("OADB object fqyetapvzn is not available in the file\n");
          return;
      }
      for (int iRun = 0; iRun < NRUNNUM; ++iRun){
        if(!(contQyEtaPVzN->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object fqyetapvzn is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hTPCSinMeanRead[iRun][1] = ((TH1D*) contQyEtaPVzN->GetObject(fRunNumList[iRun].Atoi()));
      }

      // Eta<0 Vz>0
      AliOADBContainer* contQxEtaNVzP = (AliOADBContainer*) fileTPCCalib15o->Get("fqxetanvzp");
      if(!contQxEtaNVzP){
          printf("OADB object fqxetanvzp is not available in the file\n");
          return;
      }
      for (int iRun = 0; iRun < NRUNNUM; ++iRun){
        if(!(contQxEtaNVzP->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object fqxetanvzp is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hTPCCosMeanRead[iRun][2] = ((TH1D*) contQxEtaNVzP->GetObject(fRunNumList[iRun].Atoi()));
      }
      
      AliOADBContainer* contQyEtaNVzP = (AliOADBContainer*) fileTPCCalib15o->Get("fqyetanvzp");
      if(!contQyEtaNVzP){
          printf("OADB object fqyetanvzp is not available in the file\n");
          return;
      }
      for (int iRun = 0; iRun < NRUNNUM; ++iRun){
        if(!(contQyEtaNVzP->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object fqyetanvzp is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hTPCSinMeanRead[iRun][2] = ((TH1D*) contQyEtaNVzP->GetObject(fRunNumList[iRun].Atoi()));
      }

      // Eta<0 Vz<0
      AliOADBContainer* contQxEtaNVzN = (AliOADBContainer*) fileTPCCalib15o->Get("fqxetanvzn");
      if(!contQxEtaNVzN){
          printf("OADB object fqxetanvzn is not available in the file\n");
          return;
      }
      for (int iRun = 0; iRun < NRUNNUM; ++iRun){
        if(!(contQxEtaNVzN->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object fqxetanvzn is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hTPCCosMeanRead[iRun][3] = ((TH1D*) contQxEtaNVzN->GetObject(fRunNumList[iRun].Atoi()));
      }
      
      AliOADBContainer* contQyEtaNVzN = (AliOADBContainer*) fileTPCCalib15o->Get("fqyetanvzn");
      if(!contQyEtaNVzN){
          printf("OADB object fqyetanvzn is not available in the file\n");
          return;
      }
      for (int iRun = 0; iRun < NRUNNUM; ++iRun){
        if(!(contQyEtaNVzN->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object fqyetanvzn is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hTPCSinMeanRead[iRun][3] = ((TH1D*) contQyEtaNVzN->GetObject(fRunNumList[iRun].Atoi()));
      }

      if (fFillTPCShift){
        pTPCShiftFillCoeffCos[0]  = new TProfile2D("TPCShiftCoeffCos" ,"",10,0,10,20,0,20);
        pTPCShiftFillCoeffSin[0]   = new TProfile2D("TPCShiftCoeffSin" ,"",10,0,10,20,0,20);
        pTPCShiftFillCoeffCos[1]  = new TProfile2D("TPCShiftCoeffCosC","",10,0,10,20,0,20);
        pTPCShiftFillCoeffSin[1]   = new TProfile2D("TPCShiftCoeffSinC","",10,0,10,20,0,20);
        pTPCShiftFillCoeffCos[2]  = new TProfile2D("TPCShiftCoeffCosA","",10,0,10,20,0,20);
        pTPCShiftFillCoeffSin[2]   = new TProfile2D("TPCShiftCoeffSinA","",10,0,10,20,0,20);
        fOutputList->Add(pTPCShiftFillCoeffCos[0]);
        fOutputList->Add(pTPCShiftFillCoeffSin[0]);
        fOutputList->Add(pTPCShiftFillCoeffCos[1]);
        fOutputList->Add(pTPCShiftFillCoeffSin[1]);
        fOutputList->Add(pTPCShiftFillCoeffCos[2]);
        fOutputList->Add(pTPCShiftFillCoeffSin[2]);
      }
    }

    if (fTPCCalib) {
      TFile* fileTPCShift = TFile::Open("shiftTPC15oP2.root","READ");
      if (!fileTPCShift) {
          printf("OADB TPC shift file cannot be opened\n");
          return;
      }

      TList* listTPCCalib =  (TList*)fileTPCShift->Get("output");
      hTPCShiftReadCoeffCos[0] = (TH2D*)listTPCCalib->FindObject("TPCShiftCoeffCos");
      hTPCShiftReadCoeffCos[1] = (TH2D*)listTPCCalib->FindObject("TPCShiftCoeffCosC");
      hTPCShiftReadCoeffCos[2] = (TH2D*)listTPCCalib->FindObject("TPCShiftCoeffCosA");
      hTPCShiftReadCoeffSin[0] =  (TH2D*)listTPCCalib->FindObject("TPCShiftCoeffSin");
      hTPCShiftReadCoeffSin[1] =  (TH2D*)listTPCCalib->FindObject("TPCShiftCoeffSinC");
      hTPCShiftReadCoeffSin[1] =  (TH2D*)listTPCCalib->FindObject("TPCShiftCoeffSinA");

      //Read shift run-by-run
      // AliOADBContainer* contQxShiftM = (AliOADBContainer*) fileTPCShift->Get("");
      // if(!contQxShiftM){
      //     printf("OADB object hMultV0BefCorr is not available in the file\n");
      //     return;
      // }
      // for (int iRun = 0; iRun < NRUNNUM; ++iRun){
      //   if(!(contQxShiftM->GetObject(fRunNumList[iRun].Atoi() ) ) ){
      //       printf("OADB object  is not available for run %i\n", fRunNumList[iRun].Atoi() );
      //       return;
      //   }
      //   pTPCShiftCoeffCos[0] = ((TH1D*) contQxShiftM->GetObject(fRunNumList[iRun].Atoi()));
      // }

      // AliOADBContainer* contQyShiftM = (AliOADBContainer*) fileTPCShift->Get("");
      // if(!contQyShiftM){
      //     printf("OADB object hMultV0BefCorr is not available in the file\n");
      //     return;
      // }
      // for (int iRun = 0; iRun < NRUNNUM; ++iRun){
      //   if(!(contQyShiftM->GetObject(fRunNumList[iRun].Atoi() ) ) ){
      //       printf("OADB object  is not available for run %i\n", fRunNumList[iRun].Atoi() );
      //       return;
      //   }
      //   pTPCShiftCoeffSin[0] = ((TH1D*) contQyShiftM->GetObject(fRunNumList[iRun].Atoi()));
      // }

      // AliOADBContainer* contQxShiftC = (AliOADBContainer*) fileTPCShift->Get("");
      // if(!contQxShiftC){
      //     printf("OADB object hMultV0BefCorr is not available in the file\n");
      //     return;
      // }
      // for (int iRun = 0; iRun < NRUNNUM; ++iRun){
      //   if(!(contQxShiftC->GetObject(fRunNumList[iRun].Atoi() ) ) ){
      //       printf("OADB object  is not available for run %i\n", fRunNumList[iRun].Atoi() );
      //       return;
      //   }
      //   pTPCShiftCoeffCos[1] = ((TH1D*) contQxShiftC->GetObject(fRunNumList[iRun].Atoi()));
      // }

      // AliOADBContainer* contQyShiftC = (AliOADBContainer*) fileTPCShift->Get("");
      // if(!contQyShiftC){
      //     printf("OADB object hMultV0BefCorr is not available in the file\n");
      //     return;
      // }
      // for (int iRun = 0; iRun < NRUNNUM; ++iRun){
      //   if(!(contQyShiftC->GetObject(fRunNumList[iRun].Atoi() ) ) ){
      //       printf("OADB object  is not available for run %i\n", fRunNumList[iRun].Atoi() );
      //       return;
      //   }
      //   pTPCShiftCoeffSin[1] = ((TH1D*) contQyShiftC->GetObject(fRunNumList[iRun].Atoi()));
      // }

      // AliOADBContainer* contQxShiftA = (AliOADBContainer*) fileTPCShift->Get("");
      // if(!contQxShiftA){
      //     printf("OADB object hMultV0BefCorr is not available in the file\n");
      //     return;
      // }
      // for (int iRun = 0; iRun < NRUNNUM; ++iRun){
      //   if(!(contQxShiftA->GetObject(fRunNumList[iRun].Atoi() ) ) ){
      //       printf("OADB object  is not available for run %i\n", fRunNumList[iRun].Atoi() );
      //       return;
      //   }
      //   pTPCShiftCoeffCos[1] = ((TH1D*) contQxShiftA->GetObject(fRunNumList[iRun].Atoi()));
      // }

      // AliOADBContainer* contQyShiftA = (AliOADBContainer*) fileTPCShift->Get("");
      // if(!contQyShiftA){
      //     printf("OADB object hMultV0BefCorr is not available in the file\n");
      //     return;
      // }
      // for (int iRun = 0; iRun < NRUNNUM; ++iRun){
      //   if(!(contQyShiftA->GetObject(fRunNumList[iRun].Atoi() ) ) ){
      //       printf("OADB object  is not available for run %i\n", fRunNumList[iRun].Atoi() );
      //       return;
      //   }
      //   pTPCShiftCoeffSin[1] = ((TH1D*) contQyShiftA->GetObject(fRunNumList[iRun].Atoi()));
      // }
    }

    for (int iCent=0; iCent<10; ++iCent) {
      for (int i=0; i<12; ++i) {
        hPsiTPC[iCent][i] = new TH1D(Form("psiTPCCent%i_%i",iCent,i), "",360,0,2*3.142);
        fOutputList->Add(hPsiTPC[iCent][i]);
      }
    }

    for (int i=0; i<9; ++i) {
      hQxTPCCent[i] = new TH2D(Form("qxTPCVsCent_%i",i), "",100, 0, 100, 32, -600, 600);
      hQyTPCCent[i] = new TH2D(Form("qyTPCVsCent_%i",i), "",100, 0, 100, 32, -600, 600);
      hQnTPCCent[i] = new TH2D(Form("qnTPCVsCent_%i",i), "",100, 0, 100, 24, 0, 12);
      fOutputList->Add(hQxTPCCent[i]);
      fOutputList->Add(hQyTPCCent[i]);
      fOutputList->Add(hQnTPCCent[i]);
    }
  };

  // Dobrin 15o pass2 V0 Calib
  fSPDCutPU = new TF1("fSPDCutPU", "450. + 3.9*x", 0, 50000);
    
  Double_t parV0[8] = {33.4237, 0.953516, 0.0712137, 227.923, 8.9239, -0.00319679, 0.000306314, -7.6627e-07};
  fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);
    
  Double_t parV0CL0[6] = {0.0193587, 0.975914, 0.675714, 0.0292263, -0.000549509, 5.86421e-06};
  fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);
    
  Double_t parFB32[9] = {-812.822, 6.41796, 5421.83, -0.382601, 0.0299686, -26.6249, 321.388, -0.82615, 0.0167828};
  fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
  fMultCutPU->SetParameters(parFB32);


  // Rihan Pile-up function 18?
  // fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  // Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  // fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  // fV0CutPU->SetParameters(parV0);
  
  // Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  // fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  // fCenCutLowPU->SetParameters(parV0CL0);
  // fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  // fCenCutHighPU->SetParameters(parV0CL0);

  // Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  // fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  // fMultCutPU->SetParameters(parFB32);

  PostData(1,fOutputList);
}

//---------------------------------------------------
void AliAnalysisTaskEPCalib::UserExec(Option_t *)
{
  hEvtCount->Fill(1);

  //----------------------------
  // Info
  //----------------------------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else hEvtCount->Fill(20);
  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else hEvtCount->Fill(21);
  AliAODEvent* fAOD = (AliAODEvent*)InputEvent();
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else hEvtCount->Fill(22);
  AliPIDResponse* fPID = handler->GetPIDResponse();
  if (!fPID) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else hEvtCount->Fill(23);
  AliAnalysisUtils* fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else hEvtCount->Fill(24);
  if (fPeriod.EqualTo("LHC15o")){
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    if (!fMultSel) {
      AliError(Form("%s: Could not get AliMultSelection", GetName()));
    } else hEvtCount->Fill(25);
    if (!manager || !handler || !fAOD || !fPID || !fUtils || !fMultSel) return;
  }
  else if (!manager || !handler || !fAOD || !fPID || !fUtils) return;
  hEvtCount->Fill(2);
  if (fDebug) Printf("Info done!");

  //----------------------------
  // Trigger
  //----------------------------
  ULong64_t mask = handler->IsEventSelected();
  Bool_t isTrigselected = false;
        if (fTrigger == 0) { // Run2
      isTrigselected = mask&AliVEvent::kINT7;
        } else if (fTrigger == 1) {
      isTrigselected = mask&AliVEvent::kMB; //10h
        } else if (fTrigger == 2) {
      isTrigselected = mask&AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral; //11h
        }
  if(isTrigselected == false) return;
  hEvtCount->Fill(3);
  if (fDebug) Printf("trigger done!");

  //------------------
  // event
  //------------------
  // run number
  fRunNum = fAOD->GetRunNumber();
  fRunNumBin = GetRunNumBin(fRunNum);
  if (fRunNumBin<0) return;
  hRunNumBin->Fill(fRunNumBin);
  hEvtCount->Fill(4);
  if (fDebug) Printf("run nummbr done!");

  // vertex
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  double vx    = fVtx->GetX();
  double vy    = fVtx->GetY();
  double vz    = fVtx->GetZ();
  double dz    = vz - fAOD->GetPrimaryVertexSPD()->GetZ();
  if (fabs(vz)>(double)fZvtxCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors()<1) return;
  hVz[0]->Fill(vz);
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
  // fEventCuts.SetCentralityEstimators("V0M","CL1");
  // if (!fEventCuts->AcceptEvent(fAOD) ) return;
  if (fPeriod.EqualTo("LHC10h")  || fPeriod.EqualTo("LHC11h") ) {if (fabs(dz)>0.5) return;}
  else if (fPeriod.EqualTo("LHC15o") ){
      double covTrc[6],covSPD[6];
      fVtx->GetCovarianceMatrix(covTrc);
      fAOD->GetPrimaryVertexSPD()->GetCovarianceMatrix(covSPD);
      double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
      double errTrc = TMath::Sqrt(covTrc[5]);
      double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
      if (fabs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return;
  }
  hVz[1]->Fill(vz); 
  for (int i = 0; i < 20; ++i) {
      if (vz > -10+i*1 && vz < -10+(i+1)*1) {fVzBin = i; break;}
  }
  if (fVzBin<-990) return;
  hEvtCount->Fill(5);
  if (fDebug) Printf("vertex done!");

  // centrality
  // centrality
  double centV0M = -1, centTRK = -1, centSPD0 = -1, centSPD1 = -1, centV0A = -1;
  if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC11h")){
    centV0M    =  fAOD->GetCentrality()->GetCentralityPercentile("V0M");
    centTRK    =  fAOD->GetCentrality()->GetCentralityPercentile("TRK");
    centSPD0  =  fAOD->GetCentrality()->GetCentralityPercentile("CL0");
    centSPD1  =  fAOD->GetCentrality()->GetCentralityPercentile("CL1");
    centV0A     =  fAOD->GetCentrality()->GetCentralityPercentile("V0A");
  }
  else if (fPeriod.EqualTo("LHC15o") ){
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    centV0M    =  fMultSel->GetMultiplicityPercentile("V0M");
    centTRK    =  fMultSel->GetMultiplicityPercentile("TRK");
    centSPD0 =  fMultSel->GetMultiplicityPercentile("CL0");
    centSPD1 =  fMultSel->GetMultiplicityPercentile("CL1");
  }
  fCent  = centV0M;
  fCentSPD = centSPD1;
  hCentQA[0]->Fill(centV0M,centSPD1);
  hCentQA[2]->Fill(centV0M,centTRK);
  hCentQA[4]->Fill(centV0M,centSPD0);
  hCentQA[6]->Fill(centSPD1,centSPD0);
  if (fabs(fCent-centSPD1)>fCentCut) return;
  hCentQA[1]->Fill(centV0M,centSPD1);
  hCentQA[3]->Fill(centV0M,centTRK);
  hCentQA[5]->Fill(centV0M,centSPD0);
  hCentQA[7]->Fill(centSPD1,centSPD0);
  if (fCent<0) return;
  // cent bin
  fCentBin = (int)fCent/10;
  if (fCentBin<fCbinLo || fCentBin>=fCbinHg) return;
  hCent[0]->Fill(fCent);
  hEvtCount->Fill(6);
  if (fDebug) Printf("centrality done!");

  // pile up
  if (fPeriod.EqualTo("LHC10h")) {if(!RemovalForRun1(fAOD, fUtils) ) return;}
  if (fPeriod.EqualTo("LHC11h")) {if(!RemovalForRun1(fAOD, fUtils) ) return;}
  if (fPeriod.EqualTo("LHC15o")) {if (!RejectEvtTFFit (fAOD)) return;} // 15o_pass2
  hCent[1]->Fill(fCent);
  hEvtCount->Fill(7);
  if (fDebug) Printf("pile-up done!");

  if (fTPCEstOn) TPCPlane(fAOD);
  if (fVZEROEstOn) V0Plane(fAOD);

  PostData(1,fOutputList);
}

//---------------------------------------------------
int AliAnalysisTaskEPCalib::GetRunNumBin(int runNum)
{
  int runNBin=-1;
  if (fPeriod.EqualTo("LHC10h") ){
      int fRunNumList[90]={
        139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310,
        139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870,
        138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582,
        138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201,
        138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693,
        137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541,
        137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243,
        137236, 137235, 137232, 137231, 137162, 137161};
      for (int i = 0; i < 90; ++i) {
        if (runNum==fRunNumList[i]) {runNBin=i; break;}
        else continue;
      }
  } else if (fPeriod.EqualTo("LHC11h") ){
      int fRunNumList[39]={
        170387,170040,170268,170228,170207,169838,170159,170204,170311,170084,
        169835,170088,170593,170203,170270,169846,170163,170388,170155,170083,
        170572,169837,169855,170306,170269,170089,170309,170091,170081,170230,
        170085,170315,170027,170193,170312,170313,170308,169858,169859};
        for (int i = 0; i < 39; ++i) {
          if (runNum==fRunNumList[i]) {runNBin=i; break;}
          else continue;
        }
  } else if (fPeriod.EqualTo("LHC15o") ){
      int fRunNumList[138]={
        246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 
        246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 
        246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246434, 
        246431, 246424, 246392, 246391, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 
        246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 
        246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 
        245952, 245949, 245923, 245833, 245831, 245829, 245793, 245785, 245775, 245766, 245759, 245752, 
        245731, 245729, 245705, 245702, 245692, 245683, 245554, 245545, 245544, 245543, 245542, 245540, 
        245535, 245507, 245505, 245504, 245501, 245497, 245496, 245454, 245453, 245450, 245446, 245441, 
        245411, 245410, 245409, 245407, 245401, 245397, 245396, 245353, 245349, 245347, 245346, 245345, 
        245343, 245259, 245233, 245232, 245231, 245152, 245151, 245146, 245145, 245068, 245066, 245064, 
        244983, 244982, 244980, 244975, 244918, 244917};
      for (int i = 0; i < 138; ++i) {
        if (runNum==fRunNumList[i]) {runNBin=i; break;}
        else continue;
      }
  }
  return runNBin;
}

//---------------------------------------------------

bool AliAnalysisTaskEPCalib::RemovalForRun1(AliAODEvent* fAOD, AliAnalysisUtils* fUtils)
{
  // pileup
  fUtils->SetUseOutOfBunchPileUp(true);
  fUtils->SetUseMVPlpSelection(true);
  // fUtils->SetMinPlpContribMV(5);
  bool isPileup = fUtils->IsPileUpEvent(fAOD);
  // bool isPileup = fUtils->IsPileUpMV(fAOD); // pp, p-Pb
  if (isPileup) return false;
  return true;   
}

//---------------------------------------------------

bool AliAnalysisTaskEPCalib::RejectEvtTFFit(AliAODEvent* fAOD)
{
  Float_t centV0M=-1.;
  Float_t centCL1=-1.;
  Float_t centCL0=-1.;

  AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
 
  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }
  centV0M = (Float_t) fCent;
  centCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  Int_t nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
    
  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTrk = 0;
  for (Int_t it = 0; it < nTracks; it++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
      if (!aodTrk){
          delete aodTrk;
          continue;
      }
        
      if (aodTrk->TestFilterBit(32)){
        if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
      }
  }

  hMultCentQA[0]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)

  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  // pile-up cuts
  if (centCL0 < fCenCutLowPU->Eval(centV0M)) return false;        
  if (centCL0 > fCenCutHighPU->Eval(centV0M)) return false;
           
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
           
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;

  if (Float_t(multTrk) < fMultCutPU->Eval(centV0M)) return false;

  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;

  if (fAOD->IsIncompleteDAQ()) return false;

  hMultCentQA[1]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}

//---------------------------------------------------

void AliAnalysisTaskEPCalib::TPCPlane(AliAODEvent* fAOD)
{
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  double mag    = fAOD->GetMagneticField();
  // sum q vector for psi
  double sumCos[3][3]  = {0}; // [M,C,A][Raw, NUA, Shift]
  double sumSin[3][3]  = {0};
  double mult[3][2] = {0};

  // loop tracks
  int nTrk = fAOD->GetNumberOfTracks();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrk);
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
    if (!track->TestFilterBit(fFltbit)) continue; 
    if (!AcceptAODTrack(fAOD, track, fVtx)) continue;

    double phi    = track->Phi();
    double weight=1;
    double vz    = fVtx->GetZ();
    double pt     = track->Pt();
    double eta    = track->Eta();
    int        charge = track->Charge();
    hPhi[0]->Fill(phi);

    double cosnphi = cos(fHarmonic*phi);
    double sinnphi = sin(fHarmonic*phi);
    sumCos[0][0] += cosnphi;
    sumSin[0][0] += sinnphi;
    mult[0][0]++;
    if (eta>0) {
      sumCos[1][0] += cosnphi;
      sumSin[1][0] += sinnphi;
      mult[1][0]++;
    } else {
      sumCos[2][0] += cosnphi;
      sumSin[2][0] += sinnphi;
      mult[2][0]++;
    }

    if (fTPCNUAWeight) {
      double wAcc = GetNUACor(charge, phi, eta, vz);
      if (wAcc<0) continue;
      hPhi[1]->Fill(phi, wAcc);      
      cosnphi = wAcc*cos(fHarmonic*phi);
      sinnphi = wAcc*sin(fHarmonic*phi);
      sumCos[0][1] += cosnphi;
      sumSin[0][1] += sinnphi;
      mult[0][1] += wAcc; 
      if (eta>0) {
        sumCos[1][1] += cosnphi;
        sumSin[1][1] += sinnphi;
        mult[1][1] += wAcc; 
      } else {
        sumCos[2][1] += cosnphi;
        sumSin[2][1] += sinnphi;
        mult[2][1] += wAcc;
      }
    }

    if (fFillTPCQMean) {
      if (eta>0 && vz>0) {
        pTPCCosMeanFill[fRunNumBin][0]->Fill(fCentBin,cosnphi);
        pTPCSinMeanFill[fRunNumBin][0]->Fill(fCentBin,sinnphi);
      } else if (eta>0 && vz<0) {
        pTPCCosMeanFill[fRunNumBin][1]->Fill(fCentBin,cosnphi);
        pTPCSinMeanFill[fRunNumBin][1]->Fill(fCentBin,sinnphi);
      } else if(eta<0 && vz>0) {
        pTPCCosMeanFill[fRunNumBin][2]->Fill(fCentBin,cosnphi);
        pTPCSinMeanFill[fRunNumBin][2]->Fill(fCentBin,sinnphi);
      } else if(eta<0 && vz<0) {
        pTPCCosMeanFill[fRunNumBin][3]->Fill(fCentBin,cosnphi);
        pTPCSinMeanFill[fRunNumBin][3]->Fill(fCentBin,sinnphi);
      }
    }

    if (fFillTPCShift || fTPCCalib) {
      double cosMean = 0, sinMean = 0;
      if (eta>0 && vz>0) {
        cosMean = hTPCCosMeanRead[fRunNumBin][0]->GetBinContent(fCentBin+1);
        sinMean = hTPCSinMeanRead[fRunNumBin][0]->GetBinContent(fCentBin+1);
      } else if (eta>0 && vz<0) {
        cosMean = hTPCCosMeanRead[fRunNumBin][1]->GetBinContent(fCentBin+1);
        sinMean = hTPCSinMeanRead[fRunNumBin][1]->GetBinContent(fCentBin+1);
      } else if(eta<0 && vz>0) {
        cosMean = hTPCCosMeanRead[fRunNumBin][2]->GetBinContent(fCentBin+1);
        sinMean = hTPCSinMeanRead[fRunNumBin][2]->GetBinContent(fCentBin+1);
      } else if(eta<0 && vz<0) {
        cosMean = hTPCCosMeanRead[fRunNumBin][3]->GetBinContent(fCentBin+1);
        sinMean = hTPCSinMeanRead[fRunNumBin][3]->GetBinContent(fCentBin+1);
      }
      sumCos[0][2] += (cosnphi-cosMean); 
      sumSin[0][2] += (sinnphi-sinMean);
      if (eta>0) {
        sumCos[1][2] += (cosnphi-cosMean);
        sumSin[1][2] += (sinnphi-sinMean);
      } else {
        sumCos[2][2] += (cosnphi-cosMean);
        sumSin[2][2] += (sinnphi-sinMean);
      }
    }
  }; // loop track end

  if (fFillTPCQMean) hEvtCount->Fill(11);
  if (fFillTPCShift) hEvtCount->Fill(12);

  double psi[3] = {0}, qn[3] = {0};
  for (int i = 0; i < 3; ++i){ // Raw : M, C, A
    psi[i] = GetEventPlane(sumCos[i][0], sumSin[i][0]);
    hQxTPCCent[i]->Fill(fCent, sumCos[i][0]);
    hQyTPCCent[i]->Fill(fCent, sumSin[i][0]);
    qn[i] = sqrt(sumCos[i][0]*sumCos[i][0] + sumSin[i][0]*sumSin[i][0]) / sqrt(mult[i][0]);
    hQnTPCCent[i]->Fill(fCent, qn[i]);
    hPsiTPC[fCentBin][i]->Fill(psi[i]);
  }

  if (fTPCNUAWeight) {
    for (int i = 0; i < 3; ++i){ // NUA : M, C, A
      psi[i] = GetEventPlane(sumCos[i][1], sumSin[i][1]);
      hQxTPCCent[i+3]->Fill(fCent, sumCos[i][1]);
      hQyTPCCent[i+3]->Fill(fCent, sumSin[i][1]);
      qn[i] = sqrt(sumCos[i][1]*sumCos[i][1] + sumSin[i][1]*sumSin[i][1]) / sqrt(mult[i][1]);
      hQnTPCCent[i+3]->Fill(fCent, qn[i]);
      hPsiTPC[fCentBin][i+3]->Fill(psi[i]);
    }
  }

  if (fFillTPCShift) {
    for (int i = 0; i < 3; ++i){ // Recenter : M, C, A
      psi[i] = GetEventPlane(sumCos[i][2], sumSin[i][2]);
      hQxTPCCent[i+6]->Fill(fCent, sumCos[i][2]);
      hQyTPCCent[i+6]->Fill(fCent, sumSin[i][2]);
      qn[i] = sqrt(sumCos[i][2]*sumCos[i][2] + sumSin[i][2]*sumSin[i][2]) / sqrt(mult[i][2]);
      hQnTPCCent[i+6]->Fill(fCent, qn[i]);
      hPsiTPC[fCentBin][i+6]->Fill(psi[i]);
      for (int j = 1; j <= 20; ++j) {
        pTPCShiftFillCoeffCos[i]->Fill(fCentBin, j-1, cos(j*fHarmonic*psi[i]));
        pTPCShiftFillCoeffSin[i]->Fill(fCentBin, j-1, sin(j*fHarmonic*psi[i]));
      }
    }
  }

  if (fTPCCalib) {
    for (int i = 0; i < 3; ++i){ // Calib : M, C, A
      psi[i] = GetEventPlane(sumCos[i][2], sumSin[i][2]);
      hQxTPCCent[i+6]->Fill(fCent, sumCos[i][2]);
      hQyTPCCent[i+6]->Fill(fCent, sumSin[i][2]);
      qn[i] = sqrt(sumCos[i][2]*sumCos[i][2] + sumSin[i][2]*sumSin[i][2]) / sqrt(mult[i][2]);
      hQnTPCCent[i+6]->Fill(fCent, qn[i]);
      hPsiTPC[fCentBin][i+6]->Fill(psi[i]);
      double psiShift = 0;
      for (int j = 1; j <= 20; ++j) {
        int binRead = hTPCShiftReadCoeffCos[i]->FindBin(fCentBin, j-1);
        double shiftCos   = hTPCShiftReadCoeffCos[i]->GetBinContent(binRead);
        double shiftSin    = hTPCShiftReadCoeffSin[i]->GetBinContent(binRead);
        psiShift += (2/i/fHarmonic)*(-shiftSin*cos(i*psi[i])+shiftCos*sin(i*psi[i]));
      }
      hPsiTPC[fCentBin][i+9]->Fill(psiShift);
    }
    hEvtCount->Fill(13);
  }

  return;
}

//---------------------------------------------------
void AliAnalysisTaskEPCalib::V0Plane(AliAODEvent* fAOD)
{
  AliEventplane* V0EP = fAOD->GetEventplane();
  double psiV0=-999, psiV0A=-999, psiV0C=-999;
  if (fVZEROCalib){
    psiV0  = TVector2::Phi_0_2pi(V0EP->GetEventplane("V0",  fAOD)); if (psiV0>TMath::Pi())  psiV0 -=TMath::Pi();
    psiV0A = TVector2::Phi_0_2pi(V0EP->GetEventplane("V0A", fAOD)); if (psiV0A>TMath::Pi()) psiV0A-=TMath::Pi();
    psiV0C = TVector2::Phi_0_2pi(V0EP->GetEventplane("V0C", fAOD)); if (psiV0C>TMath::Pi()) psiV0C-=TMath::Pi();
    hPsiVZERODirectGet[fCentBin][0]->Fill(psiV0);
    hPsiVZERODirectGet[fCentBin][1]->Fill(psiV0A);
    hPsiVZERODirectGet[fCentBin][2]->Fill(psiV0C);
  }
  if (fDebug) Printf("Directly Get Psi done!");

  Int_t iCentSPD = (Int_t)fCentSPD;
  if (iCentSPD >= 90) return;

  double qx[3] = {0}, qy[3] = {0}, qxGE[3] = {0}, qyGE[3] = {0}, qxRecenter[3] = {0}, qyRecenter[3] = {0};
  double multRing[3] = {0}, multRingGE[3] = {0};
  // [0]: M; [1]: C; [2]: A;
  for(int iCh = 0; iCh < 64; ++iCh) {
    double phi = TMath::Pi()/8. + TMath::Pi()/4.*(iCh%8);
    double mult = 0.;

    if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC11h")) mult= fAOD->GetVZEROEqMultiplicity(iCh);
    else if (fPeriod.EqualTo("LHC15o")) {
      AliAODVZERO* aodV0 = fAOD->GetVZEROData();
      mult = aodV0->GetMultiplicity(iCh);
    }

    qx[0] += mult*TMath::Cos(fHarmonic*phi);
    qy[0] += mult*TMath::Sin(fHarmonic*phi);
    multRing[0] += mult;
    if (fVZEROGainEq) hMultV0Fill[fRunNumBin]->Fill(iCh, mult);

    if (iCh<32) { // C
      qx[1] += mult*TMath::Cos(fHarmonic*phi);
      qy[1] += mult*TMath::Sin(fHarmonic*phi);
      multRing[1] += mult;

      if (fFillVZEROQMean || fVZEROCalib){
        double multCorGEC = -1;

        if (iCh < 8)
            multCorGEC = mult/hMultV0Read[fRunNumBin]->GetBinContent(iCh+1)*hMultV0Read[fRunNumBin]->GetBinContent(1);
        else if (iCh >= 8 && iCh < 16)
            multCorGEC = mult/hMultV0Read[fRunNumBin]->GetBinContent(iCh+1)*hMultV0Read[fRunNumBin]->GetBinContent(9);
        else if (iCh >= 16 && iCh < 24)
           multCorGEC = mult/hMultV0Read[fRunNumBin]->GetBinContent(iCh+1)*hMultV0Read[fRunNumBin]->GetBinContent(17);
        else if (iCh >= 24 && iCh < 32)
           multCorGEC = mult/hMultV0Read[fRunNumBin]->GetBinContent(iCh+1)*hMultV0Read[fRunNumBin]->GetBinContent(25);
        if (multCorGEC<0) continue;

        qxGE[1] += multCorGEC*TMath::Cos(fHarmonic*phi);
        qyGE[1] += multCorGEC*TMath::Sin(fHarmonic*phi);   
        multRingGE[1] += multCorGEC;
        qxGE[0] += multCorGEC*TMath::Cos(fHarmonic*phi);
        qyGE[0] += multCorGEC*TMath::Sin(fHarmonic*phi); 
        multRingGE[0] += multCorGEC;
      }    
    } else if (iCh>=32 && iCh<64) { // A
      qx[2] += mult*TMath::Cos(fHarmonic*phi);
      qy[2] += mult*TMath::Sin(fHarmonic*phi);
      multRing[2] += mult;

      if (fFillVZEROQMean || fVZEROCalib){
        double multCorGEA = -1;

        if (iCh >= 32 && iCh < 40)
            multCorGEA = mult/hMultV0Read[fRunNumBin]->GetBinContent(iCh+1)*hMultV0Read[fRunNumBin]->GetBinContent(33);
        else if (iCh >= 40 && iCh < 48)
            multCorGEA = mult/hMultV0Read[fRunNumBin]->GetBinContent(iCh+1)*hMultV0Read[fRunNumBin]->GetBinContent(41);
        else if (iCh >= 48 && iCh < 56)
            multCorGEA = mult/hMultV0Read[fRunNumBin]->GetBinContent(iCh+1)*hMultV0Read[fRunNumBin]->GetBinContent(49);
        else if (iCh >= 56 && iCh < 64)
            multCorGEA = mult/hMultV0Read[fRunNumBin]->GetBinContent(iCh+1)*hMultV0Read[fRunNumBin]->GetBinContent(57);
        if (multCorGEA<0) continue;

        qxGE[2] += multCorGEA*TMath::Cos(fHarmonic*phi);
        qyGE[2] += multCorGEA*TMath::Sin(fHarmonic*phi);
        multRingGE[2] += multCorGEA;  
        qxGE[0] += multCorGEA*TMath::Cos(fHarmonic*phi);
        qyGE[0] += multCorGEA*TMath::Sin(fHarmonic*phi); 
        multRingGE[0] += multCorGEA;
      }
    }
  };

  if (fDebug) Printf("mult GE done!");

  if (fVZEROGainEq) hEvtCount->Fill(8);
  if (multRing[0] < 1e-6 || multRing[1] < 1e-6 || multRing[2] < 1e-6) return;

  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  double vz = fVtx->GetZ();

  for (int i = 1; i < 3; ++i) {

    // Fill Qx Qy Mean
    if (fFillVZEROQMean){
      pV0XMeanFill[fRunNumBin][i]->Fill(iCentSPD,qxGE[i]);
      pV0XMeanFill[fRunNumBin][i]->Fill(iCentSPD,qyGE[i]);      
    }

    // Calib
    if (fVZEROCalib) {
      // Qn distribution
      double qxMean = hQx3mV0[fRunNumBin][i]->GetBinContent(iCentSPD+1);
      double qyMean = hQy3mV0[fRunNumBin][i]->GetBinContent(iCentSPD+1);

      qxRecenter[i] = qxGE[i] - qxMean;
      qyRecenter[i] = qyGE[i] - qyMean;

      double Qn_thisEvt = sqrt(qxRecenter[i]*qxRecenter[i] + qyRecenter[i]*qyRecenter[i])/ sqrt(multRingGE[i]);
      hQnCentCor[i]->Fill(fCent, Qn_thisEvt); 

      // psi
      double psiRecenter = GetEventPlane(qxRecenter[i], qyRecenter[i]);
      hPsiV0Cor[fCentBin][i]->Fill(psiRecenter); 

      if (fQAV0){
        hQxCentCor[i]->Fill(iCentSPD, qxRecenter[i]);
        hQyCentCor[i]->Fill(iCentSPD, qyRecenter[i]);
        hQxVtxCor[i]->Fill(vz, qxRecenter[i]);
        hQyVtxCor[i]->Fill(vz, qyRecenter[i]);
      }
    }
  };
  if (fDebug) Printf("Calib done!");
  if (fFillVZEROQMean) hEvtCount->Fill(9);
  if (fVZEROCalib) hEvtCount->Fill(10);

  return;
}

//---------------------------------------------------

double AliAnalysisTaskEPCalib::GetEventPlane(double qx, double qy) 
{
  double psi = TMath::ATan2(qy, qx)/fHarmonic;

  if (psi < 0) return psi+TMath::Pi();
  else return psi;
}

//---------------------------------------------------

bool AliAnalysisTaskEPCalib::AcceptAODTrack(AliAODEvent* fAOD, AliAODTrack *track, AliAODVertex* fVtx)
{
    //------------------
    // track cut  
    //------------------
    double pt     = track->Pt();
    double eta    = track->Eta();
    int    nhits  = track->GetTPCNcls();
    double dedx   = track->GetTPCsignal();
    double chi2   = track->Chi2perNDF();
    int    charge = track->Charge();
    hEta[0]->Fill(eta);
    hNhits[0]->Fill(nhits);
    if(pt < fPtMin || pt > fPtMax || fabs(eta)>fEtaCut || fabs(nhits)<fNclsCut || chi2<fChi2Lo || chi2>fChi2Hg || dedx<fDedxCut) return false;
    hPt->Fill(pt);
    hEta[1]->Fill(eta);
    hNhits[1]->Fill(nhits);
    hPDedx->Fill(track->P()*charge, dedx);
    if (fPeriod.EqualTo("LHC11h") || fPeriod.EqualTo("LHC10h")){
      //------------------
      // dca cut
      //------------------
      double mag = fAOD->GetMagneticField(); 
      double dcaxy  = 999.;
      double dcaz   = 999.;
      double r[3];
      double dca[2];
      double cov[3];
      double vx    = fVtx->GetX();
      double vy    = fVtx->GetY();
      double vz    = fVtx->GetZ();
      bool proptodca = track->PropagateToDCA(fVtx, mag, 100., dca, cov);
      if (track->GetXYZ(r)) {
        dcaxy = r[0];
        dcaz  = r[1];
      } else {
        double dcax = r[0] - vx;
        double dcay = r[1] - vy;
        dcaz  = r[2] - vz;
        dcaxy = sqrt(dcax*dcax + dcay*dcay);
        // dcaxy = dca[0];
      }
      hDcaXy[0]->Fill(dcaxy);
      if (fabs(dcaxy)>fDcaCutxy) return false;
      hDcaXy[1]->Fill(dcaxy);
      hDcaZ[0]->Fill(dcaz);
      if (fabs(dcaz)>fDcaCutz) return false;
      hDcaZ[1]->Fill(dcaz);     
    }

    return true;
}

//---------------------------------------------------

double AliAnalysisTaskEPCalib::GetNUACor(int charge, double phi, double eta, double vz)
{
  double weightNUA = 1;
  if (fVzBin<0 || fCentBin<0 || fRunNum<0) return -1;
  if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC11h")){
    if (fRunNumBin<30){
      hNUAweightPlus = (TH2D*)fListNUA1->FindObject(Form("weightdPhidEta_run%i_cent%i_vz%i_plus",fRunNum,fCentBin,fVzBin));
      hNUAweightMinus = (TH2D*)fListNUA1->FindObject(Form("weightdPhidEta_run%i_cent%i_vz%i_minus",fRunNum,fCentBin,fVzBin));
    }
    else if (fRunNumBin>=30 && fRunNumBin<60){
      hNUAweightPlus = (TH2D*)fListNUA2->FindObject(Form("weightdPhidEta_run%i_cent%i_vz%i_plus",fRunNum,fCentBin,fVzBin));
      hNUAweightMinus = (TH2D*)fListNUA2->FindObject(Form("weightdPhidEta_run%i_cent%i_vz%i_minus",fRunNum,fCentBin,fVzBin));
    }
    else if (fRunNumBin>60){
      hNUAweightPlus = (TH2D*)fListNUA3->FindObject(Form("weightdPhidEta_run%i_cent%i_vz%i_plus",fRunNum,fCentBin,fVzBin));
      hNUAweightMinus = (TH2D*)fListNUA3->FindObject(Form("weightdPhidEta_run%i_cent%i_vz%i_minus",fRunNum,fCentBin,fVzBin));
    }
    if(!hNUAweightPlus || !hNUAweightMinus) return -1;
    if (charge>0){ 
      int phiBin = hNUAweightPlus->GetXaxis()->FindBin(phi);
      int etaBin = hNUAweightPlus->GetYaxis()->FindBin(eta);  
      if (hNUAweightPlus->GetBinContent(phiBin, etaBin)>0) weightNUA = hNUAweightPlus->GetBinContent(phiBin, etaBin);
      return weightNUA;
    } else if (charge<0){
      int phiBin = hNUAweightMinus->GetXaxis()->FindBin(phi);
      int etaBin = hNUAweightMinus->GetYaxis()->FindBin(eta);  
      if (hNUAweightMinus->GetBinContent(phiBin, etaBin)>0) weightNUA = hNUAweightMinus->GetBinContent(phiBin, etaBin);
      return weightNUA;
    } 
  } else if (fPeriod.EqualTo("LHC15o")){ // Rihan and Protty 's NUA Results
    if (charge>0){ 
      hCorrectNUAPos = (TH3F*) fListNUA1->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Pos_Cent0_Run%d",fRunNum));
      if (!hCorrectNUAPos) return -1;
      int iBinNUA = hCorrectNUAPos->FindBin(vz,phi,eta); 
      if (hCorrectNUAPos->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUAPos->GetBinContent(iBinNUA);
      return  weightNUA;
    } else if (charge<0){
      hCorrectNUANeg = (TH3F*) fListNUA1->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Neg_Cent0_Run%d",fRunNum));
      if (!hCorrectNUANeg) return -1;
      int iBinNUA = hCorrectNUANeg->FindBin(vz,phi,eta); 
      if (hCorrectNUANeg->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUANeg->GetBinContent(iBinNUA);
      return weightNUA;
    } 
    // In Rihan and Protty 's NUA results, the phi distribution is independent on centrality and particle charge
  }
  return weightNUA;
}
