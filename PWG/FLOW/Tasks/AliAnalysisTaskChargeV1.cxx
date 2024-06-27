/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <random> // std::default_random_engine
#include <chrono> // std::chrono::system_clock
// ROOT classes
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1I.h"
#include "TH3I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TExMap.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TComplex.h"
#include "AliAnalysisTask.h"
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
#include "AliAODZDC.h"
#include "TGrid.h"

#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"

// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
// #include "AliEventCuts.h"
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
#include "AliOADBContainer.h"
#include "AliPID.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliAnalysisTaskChargeV1.h"

class AliAnalysisTaskChargeV1; // your analysis class

using namespace std; // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskChargeV1) // classimp: necessary for root

    AliAnalysisTaskChargeV1::AliAnalysisTaskChargeV1() : AliAnalysisTaskSE(),
                                                         fAOD(nullptr),
                                                         fOutputList(nullptr),
                                                         fQAList(nullptr),
                                                         fPID(nullptr),
                                                         mHarmonic(1.),
                                                         fFilterBit(768),
                                                         fPtMin(0.2),
                                                         fPtMax(5.),
                                                         fEtaMax(0.8),
                                                         fNhitsMin(70),
                                                         fChi2Max(4.0),
                                                         fDeDxMin(10),
                                                         fNSigmaTPCCut(3),
                                                         fNSigmaTOFCut(3),
                                                         fTrigger("kINT7"),
                                                         TriggerIsOn(false),
                                                         fVzCut(10.0),
                                                         TPCcos_t(nullptr),
                                                         TPCcos_p(nullptr),
                                                         px_P(nullptr),
                                                         px_T(nullptr),
                                                         v1_t(nullptr),
                                                         v1_p(nullptr),
                                                         fHist2Psi1ZNCCent(nullptr),
                                                         fHist2Psi1ZNACent(nullptr),
                                                         ZDCpx_P(nullptr),
                                                         ZDCpx_T(nullptr),
                                                         ZDCv1_t(nullptr),
                                                         ZDCv1_p(nullptr),
                                                         ZDCv1_t_15o(nullptr),
                                                         ZDCv1_p_15o(nullptr),
                                                         ZDCResQ(nullptr),
                                                         ZDCcos_t(nullptr),
                                                         ZDCcos_p(nullptr),
                                                         Proton_EtaPhi(nullptr),
                                                         Kion_EtaPhi(nullptr),
                                                         Pion_EtaPhi(nullptr),
                                                         PosHadron_EtaPhi(nullptr),
                                                         PosHadron_PhiPsi_p(nullptr),
                                                         PosHadron_PhiPsi_t(nullptr),
                                                         NegHadron_EtaPhi(nullptr),
                                                         NegHadron_PhiPsi_p(nullptr),
                                                         NegHadron_PhiPsi_t(nullptr)
{
  runNum = -999;
  oldRunNum = -999;
  runNumBin = -999;
  for (int i = 0; i < 3; ++i)
    vtx[i] = -999;
  vzBin = -999;
  cent = -999;
  centSPD1 = -999;
  centBin = -999;
  hEvtCount = nullptr;
  hRunNumBin = nullptr;
  hCent = nullptr;
  for (int i = 0; i < 3; ++i)
    hCentCorr[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    hVxy[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    hVz[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    fHist2DMultCentQA[i] = nullptr;
  // Track-wise
  for (int i = 0; i < 2; ++i)
    hPt[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    hEta[i] = nullptr;
  for (int i = 0; i < 8; ++i)
    hBeforePhi[i] = nullptr;
  for (int i = 0; i < 8; ++i)
    hAfterPhi[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    hNhits[i] = nullptr;
  hPDedx = nullptr;
  // pile up
  fSPDCutPU = nullptr;
  fV0CutPU = nullptr;
  fCenCutLowPU = nullptr;
  fCenCutHighPU = nullptr;
  fMultCutPU = nullptr;
  // NUE
  IsDoNUE = false;
  fListNUE = nullptr;
  hNUEweightPlus = nullptr;
  hNUEweightMinus = nullptr;
  // NUA
  IsDoNUA = false;
  fListNUA = nullptr;
  hCorrectNUAPos = nullptr;
  hCorrectNUANeg = nullptr;
  // TPC Plane
  pos1Plane = nullptr;
  neg1Plane = nullptr;
  Res1Square = nullptr;

  nCentrality = 7;
  ptEta = nullptr;
  ResQ = nullptr;
  Psi_P = nullptr;
  Psi_T = nullptr;
  Psi_PT = nullptr;
  // ZDC
  fListZDCCalib = nullptr;
  fHZDCCparameters = nullptr;
  fHZDCAparameters = nullptr;
  fProfileZDCPsi1Correlation = nullptr;
  fProfileZDCPsi2Correlation = nullptr;

  // ZDC QA
  for (int i = 0; i < 2; i++)
    fProfileZNCTowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNCQxCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNCQyCent[i] = nullptr;
  for (int i = 0; i < 3; i++)
    fHist2CalibPsi1ZNCCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNATowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNAQxCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNAQyCent[i] = nullptr;
  for (int i = 0; i < 3; i++)
    fHist2CalibPsi1ZNACent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZDCQxAQxCCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZDCQxAQyCCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZDCQyAQxCCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZDCQyAQyCCent[i] = nullptr;
  fPsi1ZNC = -999;
  fPsi1ZNA = -999;
  // PID QA
  fHistPIDPt = nullptr;
  fHistPIDEta = nullptr;
  fHistPIDPhi = nullptr;
  fHist2ProtonSigTPC = nullptr;
  fHist2ProtonSigTOF = nullptr;
  fHist2PionSigTPC = nullptr;
  fHist2PionSigTOF = nullptr;
  fHist2KionSigTPC = nullptr;
  fHist2KionSigTOF = nullptr;

  // ZDC v1
  Qtx = -999;
  Qty = -999;
  Qpx = -999;
  Qpy = -999;

  fPeriod = nullptr;
  fZDCGainAlpha = 0.395;
  fUseBadTowerCalib = false;
  fBadTowerCalibList = nullptr;
  fUseZDCSpectraCorr = false;
  fZDCSpectraCorrList = nullptr;
  fhZNSpectra = nullptr;
  fhZNSpectraCor = nullptr;
  fhZNSpectraPow = nullptr;
  fhZNBCCorr = nullptr; 
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] =  NULL;
    }
  }
  for(Int_t c=0; c<100; c++) {
    fBadTowerCalibHist[c] = NULL;
  }
  for(Int_t i=0; i<8; i++) {
    SpecCorMu1[i] = NULL;
    SpecCorMu2[i] = NULL;
    SpecCorSi[i] = NULL;
    SpecCorAv[i] = NULL;
  }
}
//_____________________________________________________________________________
AliAnalysisTaskChargeV1::AliAnalysisTaskChargeV1(const char *name) : AliAnalysisTaskSE(name),
                                                                     fAOD(nullptr),
                                                                     fOutputList(nullptr),
                                                                     fQAList(nullptr),
                                                                     fPID(nullptr),
                                                                     mHarmonic(1.),
                                                                     fFilterBit(768),
                                                                     fPtMin(0.2),
                                                                     fPtMax(5.),
                                                                     fEtaMax(0.8),
                                                                     fNhitsMin(70),
                                                                     fChi2Max(4.),
                                                                     fDeDxMin(10),
                                                                     fNSigmaTPCCut(3),
                                                                     fNSigmaTOFCut(3),
                                                                     fTrigger("kINT7"),
                                                                     TriggerIsOn(false),
                                                                     fVzCut(10.0), 
                                                                     TPCcos_t(nullptr),
                                                                     TPCcos_p(nullptr),
                                                                     px_P(nullptr),
                                                                     px_T(nullptr),
                                                                     v1_t(nullptr),
                                                                     v1_p(nullptr),
                                                                     fHist2Psi1ZNCCent(nullptr),
                                                                     fHist2Psi1ZNACent(nullptr),
                                                                     ZDCpx_P(nullptr),
                                                                     ZDCpx_T(nullptr),
                                                                     ZDCv1_t(nullptr),
                                                                     ZDCv1_p(nullptr),
                                                                     ZDCv1_t_15o(nullptr),
                                                                     ZDCv1_p_15o(nullptr),
                                                                     ZDCResQ(nullptr),
                                                                     ZDCcos_t(nullptr),
                                                                     ZDCcos_p(nullptr),
                                                                     Proton_EtaPhi(nullptr),
                                                                     Kion_EtaPhi(nullptr),
                                                                     Pion_EtaPhi(nullptr),
                                                                     PosHadron_EtaPhi(nullptr),
                                                                     PosHadron_PhiPsi_p(nullptr),
                                                                     PosHadron_PhiPsi_t(nullptr),
                                                                     NegHadron_EtaPhi(nullptr),
                                                                     NegHadron_PhiPsi_p(nullptr),
                                                                     NegHadron_PhiPsi_t(nullptr)
{
  runNum = -999;
  oldRunNum = -999;
  runNumBin = -999;
  for (int i = 0; i < 3; ++i)
    vtx[i] = -999;
  vzBin = -999;
  cent = -999;
  centSPD1 = -999;
  centBin = -999;
  hEvtCount = nullptr;
  hRunNumBin = nullptr;
  hCent = nullptr;
  for (int i = 0; i < 3; ++i)
    hCentCorr[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    hVxy[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    hVz[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    fHist2DMultCentQA[i] = nullptr;
  // Track-wise
  for (int i = 0; i < 2; ++i)
    hPt[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    hEta[i] = nullptr;
  for (int i = 0; i < 8; ++i)
    hBeforePhi[i] = nullptr;
  for (int i = 0; i < 8; ++i)
    hAfterPhi[i] = nullptr;
  for (int i = 0; i < 2; ++i)
    hNhits[i] = nullptr;
  hPDedx = nullptr;
  // pile up
  fSPDCutPU = nullptr;
  fV0CutPU = nullptr;
  fCenCutLowPU = nullptr;
  fCenCutHighPU = nullptr;
  fMultCutPU = nullptr;
  // NUE
  IsDoNUE = true;
  fListNUE = nullptr;
  hNUEweightPlus = nullptr;
  hNUEweightMinus = nullptr;
  // NUA
  IsDoNUA = true;
  fListNUA = nullptr;
  hCorrectNUAPos = nullptr;
  hCorrectNUANeg = nullptr;
  // TPC Plane
  pos1Plane = nullptr;
  neg1Plane = nullptr;
  Res1Square = nullptr;

  nCentrality = 7;
  ptEta = nullptr;
  ResQ = nullptr;
  Psi_P = nullptr;
  Psi_T = nullptr;
  Psi_PT = nullptr;

  // ZDC
  fListZDCCalib = nullptr;
  fHZDCCparameters = nullptr;
  fHZDCAparameters = nullptr;
  fProfileZDCPsi1Correlation = nullptr;
  fProfileZDCPsi2Correlation = nullptr;
  // ZDC QA
  for (int i = 0; i < 2; i++)
    fProfileZNCTowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNCQxCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNCQyCent[i] = nullptr;
  for (int i = 0; i < 3; i++)
    fHist2CalibPsi1ZNCCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNATowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNAQxCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZNAQyCent[i] = nullptr;
  for (int i = 0; i < 3; i++)
    fHist2CalibPsi1ZNACent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZDCQxAQxCCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZDCQxAQyCCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZDCQyAQxCCent[i] = nullptr;
  for (int i = 0; i < 2; i++)
    fProfileZDCQyAQyCCent[i] = nullptr;
  fPsi1ZNC = -999;
  fPsi1ZNA = -999;
  // PID QA
  fHistPIDPt = nullptr;
  fHistPIDEta = nullptr;
  fHistPIDPhi = nullptr;
  fHist2ProtonSigTPC = nullptr;
  fHist2ProtonSigTOF = nullptr;
  fHist2PionSigTPC = nullptr;
  fHist2PionSigTOF = nullptr;
  fHist2KionSigTPC = nullptr;
  fHist2KionSigTOF = nullptr;

  // ZDC v1
  Qtx = -999;
  Qty = -999;
  Qpx = -999;
  Qpy = -999;

  fPeriod = nullptr;
  fZDCGainAlpha = 0.395;
  fUseBadTowerCalib = false;
  fBadTowerCalibList = nullptr;
  fUseZDCSpectraCorr = false;
  fZDCSpectraCorrList = nullptr;
  fhZNSpectra = nullptr;
  fhZNSpectraCor = nullptr;
  fhZNSpectraPow = nullptr;
  fhZNBCCorr = nullptr; 
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] =  NULL;
    }
  }
  for(Int_t c=0; c<100; c++) {
    fBadTowerCalibHist[c] = NULL;
  }
  for(Int_t i=0; i<8; i++) {
    SpecCorMu1[i] = NULL;
    SpecCorMu2[i] = NULL;
    SpecCorSi[i] = NULL;
    SpecCorAv[i] = NULL;
  }  
  // constructor
  DefineInput(0, TChain::Class()); // define the input of the analysis: in this case we take a 'chain' of events
                                   // this chain is created by the analysis manager, so no need to worry about it,
                                   // it does its work automatically
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this case it's a list of histograms
  DefineOutput(2, TList::Class()); // you can add more output objects by calling DefineOutput(2, classname::Class())
                                   // if you add more output objects, make sure to call PostData for all of them, and to
                                   // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskChargeV1::~AliAnalysisTaskChargeV1()
{
  // destructor
  if (fListZDCCalib)  delete fListZDCCalib;
  if (fBadTowerCalibList) delete fBadTowerCalibList;
  if (fQAList)
    delete fQAList;
  if (fOutputList)
  {
    delete fOutputList; // at the end of your task, it is deleted from memory by calling this function
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskChargeV1::UserCreateOutputObjects(){
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output file
  //
  fOutputList = new TList();       // this is a list which will contain all of your histograms
  fOutputList->SetName(GetName()); // at the end of the analysis, the contents of this list are written
                                   // to the output file
  fOutputList->SetOwner(kTRUE);    // memory stuff: the list is owner of all objects it contains and will delete them
                                   // if requested (dont worry about this now)
  TH1::SetDefaultSumw2(kTRUE);
  // event-wise
  hEvtCount = new TH1I("evtCount", "", 20, 1, 21);
  hEvtCount->GetXaxis()->SetBinLabel(1, "All");
  hEvtCount->GetXaxis()->SetBinLabel(2, "Info");
  hEvtCount->GetXaxis()->SetBinLabel(3, "Run number");
  hEvtCount->GetXaxis()->SetBinLabel(4, "Vertex");
  hEvtCount->GetXaxis()->SetBinLabel(5, "Cent");
  hEvtCount->GetXaxis()->SetBinLabel(6, "Pile up");
  hEvtCount->GetXaxis()->SetBinLabel(7,"Trigger");
  hEvtCount->GetXaxis()->SetBinLabel(8, "Get VZERO Plane");
  hEvtCount->GetXaxis()->SetBinLabel(10, "Manager");
  hEvtCount->GetXaxis()->SetBinLabel(11, "Handler");
  hEvtCount->GetXaxis()->SetBinLabel(12, "AOD");
  hEvtCount->GetXaxis()->SetBinLabel(13, "PID");
  hEvtCount->GetXaxis()->SetBinLabel(14, "Utils");
  hEvtCount->GetXaxis()->SetBinLabel(17, "TPC plane");
  hEvtCount->GetXaxis()->SetBinLabel(18, "VZERO plane");
  hEvtCount->GetXaxis()->SetBinLabel(19, "ZDC plane");
  hEvtCount->GetXaxis()->SetBinLabel(20, "loops end");
  fOutputList->Add(hEvtCount);

  if(!fPeriod)
  std::cout << "!!!!!!!!!!!!!!!  fPeriod NOT been set !!!!!!!!!!!!!!!" << std::endl;

  if (fPeriod.EqualTo("LHC15o")) {
  TString runNumList[77] = {"246994","246991","246989","246984","246982","246980","246948","246945","246928","246851",
                            "246847","246846","246845","246844","246810","246809","246808","246807","246805","246804",
                            "246766","246765","246763","246760","246759","246758","246757","246751","246750","246495",
                            "246493","246488","246487","246434","246431","246428","246424","246276","246275","246272",
                            "246271","246225","246222","246217","246185","246182","246181","246180","246178","246153",
                            "246152","246151","246115","246113","246089","246053","246052","246049","246048","246042",
                            "246037","246036","246012","246003","246001","245954","245952","245949","245923","245833",
                            "245831","245829","245705","245702","245700","245692","245683"};
  //  Int_t dRun15opidfix[] = {245145, 245146, 245151, 245152, 245231, 245232, 245259, 245343, 245345, 245346, 245347, 245349, 245353, 245396, 245397, 245401, 245407, 245409, 245441, 245446, 245450, 245454, 245496, 245497, 245501, 245504, 245505, 245507, 245535, 245540, 245542, 245543, 245544, 245545, 245554};

  hRunNumBin = new TH1I("runNumBin", "", 77, 0, 77);
  for (int i = 0; i < 77; ++i)
  {
    hRunNumBin->GetXaxis()->SetBinLabel(i + 1, runNumList[i].Data());
  }
  fOutputList->Add(hRunNumBin);
  }

  if (fPeriod.EqualTo("LHC18q")||fPeriod.EqualTo("LHC18r")) {
  TString runNumList[214] = {"296623", "296622", "296621", "296619", "296618", "296616", "296615", "296594", "296553", "296552",
                             "296551", "296550", "296548", "296547", "296516", "296512", "296511", "296510", "296509", "296472",
                             "296433", "296424", "296423", "296420", "296419", "296415", "296414", "296383", "296381", "296380",
                             "296379", "296378", "296377", "296376", "296375", "296312", "296309", "296304", "296303", "296280",
                             "296279", "296273", "296270", "296269", "296247", "296246", "296244", "296243", "296242", "296241",
                             "296240", "296198", "296197", "296196", "296195", "296194", "296192", "296191", "296143", "296142",
                             "296135", "296134", "296133", "296132", "296123", "296074", "296066", "296065", "296063", "296062",
                             "296060", "296016", "295942", "295941", "295937", "295936", "295913", "295910", "295909", "295861",
                             "295860", "295859", "295856", "295855", "295854", "295853", "295831", "295829", "295826", "295825",
                             "295822", "295819", "295818", "295816", "295791", "295788", "295786", "295763", "295762", "295759",
                             "295758", "295755", "295754", "295725", "295723", "295721", "295719", "295718", "295717", "295714",
                             "295712", "295676", "295675", "295673", "295668", "295667", "295666", "295615", "295612", "295611",
                             "295610", "295589", "295588", "295586", "295585", "297595", "297590", "297588", "297558", "297544", "297542", "297541", "297540", "297537", "297512",
                             "297483", "297479", "297452", "297451", "297450", "297446", "297442", "297441", "297415", "297414",
                             "297413", "297406", "297405", "297380", "297379", "297372", "297367", "297366", "297363", "297336",
                             "297335", "297333", "297332", "297317", "297311", "297310", "297278", "297222", "297221", "297218",
                             "297196", "297195", "297193", "297133", "297132", "297129", "297128", "297124", "297123", "297119",
                             "297118", "297117", "297085", "297035", "297031", "296966", "296941", "296938", "296935", "296934",
                             "296932", "296931", "296930", "296903", "296900", "296899", "296894", "296852", "296851", "296850",
                             "296848", "296839", "296838", "296836", "296835", "296799", "296794", "296793", "296790", "296787",
                             "296786", "296785", "296784", "296781", "296752", "296694", "296693", "296691", "296690"};
  hRunNumBin = new TH1I("runNumBin", "", 214, 0, 214);
  for (int i = 0; i < 214; ++i)
  {
    hRunNumBin->GetXaxis()->SetBinLabel(i + 1, runNumList[i].Data());
  }
  fOutputList->Add(hRunNumBin);
  }

  hCent = new TH1D("centrality", "", 100, 0, 100);
  fOutputList->Add(hCent);
  hCentCorr[0] = new TH2D("centcorr0", "", 100, 0, 100, 100, 0, 100);
  hCentCorr[1] = new TH2D("centcorr1", "", 100, 0, 100, 100, 0, 100);
  hCentCorr[2] = new TH2D("centcorr2", "", 100, 0, 100, 100, 0, 100);
  for (int i = 0; i < 3; ++i)
    fOutputList->Add(hCentCorr[i]);

  hVxy[0] = new TH2D("vxy0", "", 100, -0.5, 0.5, 100, -0.5, 0.5);
  hVxy[1] = new TH2D("vxy1", "", 100, -0.5, 0.5, 100, -0.5, 0.5);
  hVz[0] = new TH1D("vz0", "", 200, -50, 50);
  hVz[1] = new TH1D("vz1", "", 200, -50, 50);
  for (int i = 0; i < 2; ++i)
    fOutputList->Add(hVxy[i]);
  for (int i = 0; i < 2; ++i)
    fOutputList->Add(hVz[i]);

  fHist2DMultCentQA[0] = new TH2D("fHist2DMultCentQA_BfCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  fHist2DMultCentQA[1] = new TH2D("fHist2DMultCentQA_AfCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  fOutputList->Add(fHist2DMultCentQA[0]);
  fOutputList->Add(fHist2DMultCentQA[1]);

  // track-wise
  hPt[0] = new TH1D("hPtBeforeCut", "", 200, 0., 20.);
  hPt[1] = new TH1D("hPtAfterCut", "", 200, 0., 20.);
  for (int i = 0; i < 2; ++i)
    fOutputList->Add(hPt[i]);
  hEta[0] = new TH1D("hEtaBeforeCut", "", 200, -10., 10.);
  hEta[1] = new TH1D("hEtaAfterCut", "", 200, -10., 10.);
  for (int i = 0; i < 2; ++i)
    fOutputList->Add(hEta[i]);
  for (int i = 0; i < 8; ++i)
  {
    hBeforePhi[i] = new TH2D(Form("hPhiBeforeCut_cent%i", i), "", 20, 0, 2 * TMath::Pi(), 3, 0, 3);
    hAfterPhi[i] = new TH2D(Form("hPhiAfterCut_cent%i", i), "", 20, 0, 2 * TMath::Pi(), 3, 0, 3);
    fOutputList->Add(hBeforePhi[i]);
    fOutputList->Add(hAfterPhi[i]);
  }
  hNhits[0] = new TH1D("hNhitsBeforeCut", "", 200, 0., 200.);
  hNhits[1] = new TH1D("hNhitsAfterCut", "", 200, 0., 200.);
  for (int i = 0; i < 2; ++i)
    fOutputList->Add(hNhits[i]);
  hPDedx = new TH2D("hPDedx", "", 400, -10., 10., 400, 0, 1000);
  fOutputList->Add(hPDedx);

  // pile up
  fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);
  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);

  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);

  ////////////////////////
  // NUE
  ////////////////////////
  if (IsDoNUE)
  {
    if (!fListNUE)
    {
      std::cout << ("NUE list not found") << std::endl;
      return;
    }

    hNUEweightPlus = (TH1D *)fListNUE->FindObject("trkEfficiencyChrgPos");
    hNUEweightMinus = (TH1D *)fListNUE->FindObject("trkEfficiencyChrgNeg");
  }

  ////////////////////////
  // NUA
  ////////////////////////
  if (IsDoNUA)
  {
    if (!fListNUA)
    {
      std::cout << ("NUA list not found") << std::endl;
      return;
    }
    hCorrectNUAPos = new TH3F();
    hCorrectNUANeg = new TH3F();
  }

  // TPC Plane
  pos1Plane = new TH1D("pos1Plane", "", 100, 0, 2 * TMath::Pi());
  fOutputList->Add(pos1Plane);
  neg1Plane = new TH1D("neg1Plane", "", 100, 0, 2 * TMath::Pi());
  fOutputList->Add(neg1Plane);
  Res1Square = new TProfile("Res1Square_cent", "", 1, 0, 1);
  fOutputList->Add(Res1Square);
  TPCcos_t = new TProfile2D *[nCentrality];
  TPCcos_p = new TProfile2D *[nCentrality];
  for (int i = 0; i < nCentrality; ++i)
  {
    TPCcos_t[i] = new TProfile2D(Form("TPCcos_t%i", i), "", 8, 0, 8, 5, -0.8, 0.8);
    TPCcos_p[i] = new TProfile2D(Form("TPCcos_p%i", i), "", 8, 0, 8, 5, -0.8, 0.8);
    fOutputList->Add(TPCcos_t[i]);
    fOutputList->Add(TPCcos_p[i]);
  }

  // scalar product method TPC
  px_P = new TProfile2D *[nCentrality];
  px_T = new TProfile2D *[nCentrality];
  v1_t = new TProfile2D *[nCentrality];
  v1_p = new TProfile2D *[nCentrality];
  for (int i = 0; i < nCentrality; ++i)
  {
    px_P[i] = new TProfile2D(Form("px_P%i", i), "", 8, 0, 8, 5, -0.8, 0.8);
    px_T[i] = new TProfile2D(Form("px_T%i", i), "", 8, 0, 8, 5, -0.8, 0.8);
    v1_t[i] = new TProfile2D(Form("v1_t%i", i), "", 8, 0, 8, 5, -0.8, 0.8);
    v1_p[i] = new TProfile2D(Form("v1_p%i", i), "", 8, 0, 8, 5, -0.8, 0.8);
    fOutputList->Add(px_P[i]);
    fOutputList->Add(px_T[i]);
    fOutputList->Add(v1_t[i]);
    fOutputList->Add(v1_p[i]);
  }
  ResQ = new TProfile("ResQ", "", 10, 0, 10);
  ptEta = new TProfile("ptEta", "", 5, -0.8, 0.8);
  fOutputList->Add(ResQ);
  fOutputList->Add(ptEta);

  Psi_P = new TH1D("Psi_P", "", 100, -2 * TMath::Pi(), 2 * TMath::Pi());
  Psi_T = new TH1D("Psi_T", "", 100, -2 * TMath::Pi(), 2 * TMath::Pi());
  Psi_PT = new TH1D("Psi_PT", "", 100, -4 * TMath::Pi(), 4 * TMath::Pi());
  fOutputList->Add(Psi_P);
  fOutputList->Add(Psi_T);
  fOutputList->Add(Psi_PT);

  ////////////////////////
  // ZDC
  ////////////////////////
  fQAList = new TList();
  fQAList->SetName("fQAList");
  fQAList->SetOwner(kTRUE);

  fHZDCCparameters = new TH1D();
  fHZDCAparameters = new TH1D();
  if (!fListZDCCalib)
  {
    std::cout << ("ZDC calibration list not found") << std::endl;
    return;
  }

  fProfileZDCPsi1Correlation = new TProfile("fProfileZDCPsi1Correlation", "fProfileZDCPsi1Correlation;centrality;Res", 10, 0., 10);
  fProfileZDCPsi2Correlation = new TProfile("fProfileZDCPsi2Correlation", "fProfileZDCPsi2Correlation;centrality;Res", 10, 0., 10);
  fOutputList->Add(fProfileZDCPsi1Correlation);
  fOutputList->Add(fProfileZDCPsi2Correlation);
  fHist2Psi1ZNCCent = new TH1D *[nCentrality];
  fHist2Psi1ZNACent = new TH1D *[nCentrality];
  for (int i = 0; i < nCentrality; ++i)
  {
    fHist2Psi1ZNCCent[i] = new TH1D(Form("fHist2Psi1ZNCCent%i", i), "", 100, -TMath::Pi(), TMath::Pi());
    fHist2Psi1ZNACent[i] = new TH1D(Form("fHist2Psi1ZNACent%i", i), "", 100, -TMath::Pi(), TMath::Pi());
    fOutputList->Add(fHist2Psi1ZNCCent[i]);
    fOutputList->Add(fHist2Psi1ZNACent[i]);
  }

  // QA
  std::string charCalibStep;
  for (int i = 0; i < 2; i++)
  {
    if (i == 0)
      charCalibStep = "GE";
    if (i == 1)
      charCalibStep = "RC";
    fProfileZNCQxCent[i] = new TProfile(Form("fProfileZNCQxCent%s", charCalibStep.data()), "", 80, 0, 80.);
    fProfileZNCQyCent[i] = new TProfile(Form("fProfileZNCQyCent%s", charCalibStep.data()), "", 80, 0, 80.);
    fHist2CalibPsi1ZNCCent[i] = new TH1D(Form("fHist2CalibPsi1ZNCCent%s", charCalibStep.data()), "", 100, -TMath::Pi(), TMath::Pi());
    fQAList->Add(fProfileZNCQxCent[i]);
    fQAList->Add(fProfileZNCQyCent[i]);
    fQAList->Add(fHist2CalibPsi1ZNCCent[i]);

    fProfileZNAQxCent[i] = new TProfile(Form("fProfileZNAQxCent%s", charCalibStep.data()), "", 80, 0, 80.);
    fProfileZNAQyCent[i] = new TProfile(Form("fProfileZNAQyCent%s", charCalibStep.data()), "", 80, 0, 80.);
    fHist2CalibPsi1ZNACent[i] = new TH1D(Form("fHist2CalibPsi1ZNACent%s", charCalibStep.data()), "", 100, -TMath::Pi(), TMath::Pi());
    fQAList->Add(fProfileZNAQxCent[i]);
    fQAList->Add(fProfileZNAQyCent[i]);
    fQAList->Add(fHist2CalibPsi1ZNACent[i]);

    fProfileZDCQxAQxCCent[i] = new TProfile(Form("fProfileZDCQxAQxCCent%s", charCalibStep.data()), "", 80, 0, 80.);
    fProfileZDCQxAQyCCent[i] = new TProfile(Form("fProfileZDCQxAQyCCent%s", charCalibStep.data()), "", 80, 0, 80.);
    fProfileZDCQyAQxCCent[i] = new TProfile(Form("fProfileZDCQyAQxCCent%s", charCalibStep.data()), "", 80, 0, 80.);
    fProfileZDCQyAQyCCent[i] = new TProfile(Form("fProfileZDCQyAQyCCent%s", charCalibStep.data()), "", 80, 0, 80.);
    fQAList->Add(fProfileZDCQxAQxCCent[i]);
    fQAList->Add(fProfileZDCQxAQyCCent[i]);
    fQAList->Add(fProfileZDCQyAQxCCent[i]);
    fQAList->Add(fProfileZDCQyAQyCCent[i]);
  }

  fHist2CalibPsi1ZNCCent[2] = new TH1D("fHist2CalibPsi1ZNCCentSF", "", 100, -TMath::Pi(), TMath::Pi());
  fHist2CalibPsi1ZNACent[2] = new TH1D("fHist2CalibPsi1ZNACentSF", "", 100, -TMath::Pi(), TMath::Pi());
  fQAList->Add(fHist2CalibPsi1ZNCCent[2]);
  fQAList->Add(fHist2CalibPsi1ZNACent[2]);

  fProfileZNCTowerMeanEnegry[0] = new TProfile("fProfileZNCTowerMeanEnegryRW", "", 5, 0, 5);
  fProfileZNCTowerMeanEnegry[1] = new TProfile("fProfileZNCTowerMeanEnegryGE", "", 5, 0, 5);
  fProfileZNATowerMeanEnegry[0] = new TProfile("fProfileZNATowerMeanEnegryRW", "", 5, 0, 5);
  fProfileZNATowerMeanEnegry[1] = new TProfile("fProfileZNATowerMeanEnegryGE", "", 5, 0, 5);
  fQAList->Add(fProfileZNCTowerMeanEnegry[0]);
  fQAList->Add(fProfileZNCTowerMeanEnegry[1]);
  fQAList->Add(fProfileZNATowerMeanEnegry[0]);
  fQAList->Add(fProfileZNATowerMeanEnegry[1]);

  // PID QA
  fHistPIDPt = new TH2D("fHistPIDPt", "fHistPIDPt;p_{T}", 8, 0, 8, 200, 0, 20);
  fHistPIDEta = new TH2D("fHistPIDEta", "fHistPIDEta;#eta", 8, 0, 8, 100, -1, 1);
  fHistPIDPhi = new TH2D("fHistPIDPhi", "fHistPIDPhi;#phi", 8, 0, 8, 100, 0, TMath::TwoPi());
  fHist2ProtonSigTPC = new TH2D("fHist2ProtonSigTPC", "fHist2ProtonSigTPC;SigTPC", 25, 0, 5, 100, 0, 10);
  fHist2ProtonSigTOF = new TH2D("fHist2ProtonSigTOF", "fHist2ProtonSigTOF;SigTOF", 25, 0., 5., 100, 0., 10.);
  fHist2PionSigTPC = new TH2D("fHist2PionSigTPC", "fHist2PionSigTPC;SigTPC", 25, 0, 5, 100, 0, 10);
  fHist2PionSigTOF = new TH2D("fHist2PionSigTOF", "fHist2PionSigTOF;SigTOF", 25, 0., 5., 100, 0., 10.);
  fHist2KionSigTPC = new TH2D("fHist2KionSigTPC", "fHist2KionSigTPC;SigTPC", 25, 0, 5, 100, 0, 10);
  fHist2KionSigTOF = new TH2D("fHist2KionSigTOF", "fHist2KionSigTOF;SigTOF", 25, 0., 5., 100, 0., 10.);
  fQAList->Add(fHistPIDPt);
  fQAList->Add(fHistPIDEta);
  fQAList->Add(fHistPIDPhi);
  fQAList->Add(fHist2ProtonSigTPC);
  fQAList->Add(fHist2ProtonSigTOF);
  fQAList->Add(fHist2PionSigTPC);
  fQAList->Add(fHist2PionSigTOF);
  fQAList->Add(fHist2KionSigTPC);
  fQAList->Add(fHist2KionSigTOF);

  /// ZDC v1pt
  ZDCpx_P = new TProfile2D *[nCentrality];
  ZDCpx_T = new TProfile2D *[nCentrality];
  ZDCv1_t = new TProfile2D *[nCentrality];
  ZDCv1_p = new TProfile2D *[nCentrality];
  ZDCv1_t_15o = new TProfile2D *[nCentrality];
  ZDCv1_p_15o = new TProfile2D *[nCentrality];
  ZDCcos_t = new TProfile2D *[nCentrality];
  ZDCcos_p = new TProfile2D *[nCentrality];
  Proton_EtaPhi = new TH2D *[nCentrality];
  Kion_EtaPhi = new TH2D *[nCentrality];
  Pion_EtaPhi = new TH2D *[nCentrality];
  PosHadron_EtaPhi = new TH2D *[nCentrality];
  PosHadron_PhiPsi_p = new TH2D *[nCentrality];
  PosHadron_PhiPsi_t = new TH2D *[nCentrality];
  NegHadron_EtaPhi = new TH2D *[nCentrality];
  NegHadron_PhiPsi_p = new TH2D *[nCentrality];
  NegHadron_PhiPsi_t = new TH2D *[nCentrality];
  for (int i = 0; i < nCentrality; ++i)
  {
    ZDCpx_P[i] = new TProfile2D(Form("ZDCpx_P%i", i), "", 8, 0, 8, 16, -0.8, 0.8);
    ZDCpx_T[i] = new TProfile2D(Form("ZDCpx_T%i", i), "", 8, 0, 8, 16, -0.8, 0.8);
    ZDCv1_t[i] = new TProfile2D(Form("ZDCv1_t%i", i), "", 8, 0, 8, 16, -0.8, 0.8);
    ZDCv1_p[i] = new TProfile2D(Form("ZDCv1_p%i", i), "", 8, 0, 8, 16, -0.8, 0.8);
    ZDCv1_t_15o[i] = new TProfile2D(Form("ZDCv1_t_15o%i", i), "", 8, 0, 8, 5, -0.8, 0.8);
    ZDCv1_p_15o[i] = new TProfile2D(Form("ZDCv1_p_15o%i", i), "", 8, 0, 8, 5, -0.8, 0.8);
    ZDCcos_t[i] = new TProfile2D(Form("ZDCcos_t%i", i), "", 8, 0, 8, 16, -0.8, 0.8);
    ZDCcos_p[i] = new TProfile2D(Form("ZDCcos_p%i", i), "", 8, 0, 8, 16, -0.8, 0.8);
    Proton_EtaPhi[i] = new TH2D(Form("Proton_EtaPhi%i", i),"", 16, -0.8, 0.8, 100, 0, 2*TMath::Pi());
    Pion_EtaPhi[i] = new TH2D(Form("Pion_EtaPhi%i", i),"", 16, -0.8, 0.8, 100, 0, 2*TMath::Pi());
    Kion_EtaPhi[i] = new TH2D(Form("Kion_EtaPhi%i", i),"", 16, -0.8, 0.8, 100, 0, 2*TMath::Pi());
    PosHadron_EtaPhi[i] = new TH2D(Form("PosHadron_EtaPhi%i", i),"", 16, -0.8, 0.8, 100, 0, 2*TMath::Pi());
    PosHadron_PhiPsi_p[i] = new TH2D(Form("PosHadron_PhiPsi_p%i", i),"", 16, -0.8, 0.8, 100, 0, 2*TMath::Pi());
    PosHadron_PhiPsi_t[i] = new TH2D(Form("PosHadron_PhiPsi_t%i", i),"", 16, -0.8, 0.8, 100, 0, 2*TMath::Pi());

    NegHadron_EtaPhi[i] = new TH2D(Form("NegHadron_EtaPhi%i", i),"", 16, -0.8, 0.8, 100, 0, 2*TMath::Pi());
    NegHadron_PhiPsi_p[i] = new TH2D(Form("NegHadron_PhiPsi_p%i", i),"", 16, -0.8, 0.8, 100, 0, 2*TMath::Pi());
    NegHadron_PhiPsi_t[i] = new TH2D(Form("NegHadron_PhiPsi_t%i", i),"", 16, -0.8, 0.8, 100, 0, 2*TMath::Pi());

    
    fOutputList->Add(ZDCpx_P[i]);
    fOutputList->Add(ZDCpx_T[i]);
    fOutputList->Add(ZDCv1_t[i]);
    fOutputList->Add(ZDCv1_p[i]);    
    fOutputList->Add(ZDCv1_t_15o[i]);
    fOutputList->Add(ZDCv1_p_15o[i]);
    fOutputList->Add(ZDCcos_t[i]);
    fOutputList->Add(ZDCcos_p[i]);
    fOutputList->Add(Proton_EtaPhi[i]);
    fOutputList->Add(Kion_EtaPhi[i]);
    fOutputList->Add(Pion_EtaPhi[i]);
    fOutputList->Add(PosHadron_EtaPhi[i]);
    fOutputList->Add(PosHadron_PhiPsi_p[i]);
    fOutputList->Add(PosHadron_PhiPsi_t[i]);
    fOutputList->Add(NegHadron_EtaPhi[i]);
    fOutputList->Add(NegHadron_PhiPsi_p[i]);
    fOutputList->Add(NegHadron_PhiPsi_t[i]);
  }
  ZDCResQ = new TProfile("ZDCResQ", "", 10, 0, 10);
  fOutputList->Add(ZDCResQ);

  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] = new TH1D();
      // fOutputList->Add(fTowerGainEq[c][i]);
    }
  }
  if(fBadTowerCalibList) {
    for(Int_t c=0; c<100; c++) {
      fBadTowerCalibHist[c] = (TH2D*)fBadTowerCalibList->FindObject(Form("TH2Resp[%d]",c));
      fOutputList->Add(fBadTowerCalibHist[c]);
    }
  }
  if(fZDCSpectraCorrList) {
    for(Int_t i=0; i<8; i++) {
      SpecCorMu1[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorMu1[%d]",i));
      fOutputList->Add(SpecCorMu1[i]);
      SpecCorMu2[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorMu2[%d]",i));
      fOutputList->Add(SpecCorMu2[i]);
      SpecCorAv[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorAv[%d]",i));
      fOutputList->Add(SpecCorAv[i]);
      SpecCorSi[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorSi[%d]",i));
      fOutputList->Add(SpecCorSi[i]);
    }
  fhZNSpectra = new TH3D("fhZNSpectra","fhZNSpectra",100,0.,100.,8,0.,8.,1000,0.,1.E5);
  fOutputList->Add(fhZNSpectra);
  fhZNSpectraCor = new TH3D("fhZNSpectraCor","fhZNSpectraCor",100,0.,100.,8,0.,8.,1000,0.,1.E5);
  fOutputList->Add(fhZNSpectraCor);
  fhZNSpectraPow = new TH3D("fhZNSpectraPow","fhZNSpectraPow",100,0.,100.,8,0.,8.,1000,0.,TMath::Power(1.E5,fZDCGainAlpha));
  fOutputList->Add(fhZNSpectraPow);
  }
  fhZNBCCorr = new TH3D("fhZNBCCorr","fhZNBCCorr",100,0.,100.,500,0.,1.E5,500,0.,1.E5);
  fOutputList->Add(fhZNBCCorr);

  PostData(1, fOutputList); // postdata will notify the analysis manager of changes / updates to the
  PostData(2, fQAList);     // fOutputList object. the manager will in the end take care of writing your output to file
                            // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskChargeV1::UserExec(Option_t *)
{
  hEvtCount->Fill(1);
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager)
  {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  }
  else
    hEvtCount->Fill(10);
  AliAODInputHandler *handler = (AliAODInputHandler *)manager->GetInputEventHandler();
  if (!handler)
  {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  }
  else
    hEvtCount->Fill(11);
  fAOD = dynamic_cast<AliAODEvent *>(InputEvent());
  if (!fAOD)
  {
    AliError(Form("%s: Could not get AOD event", GetName()));
    return;
  }
  else
    hEvtCount->Fill(12);
    fPID = handler->GetPIDResponse();
    if (!fPID)
    {
      AliError(Form("%s: Could not get PIDResponse", GetName()));
    }
   else 
    hEvtCount->Fill(13);
  AliAnalysisUtils *fUtils = new AliAnalysisUtils();
  if (!fUtils)
  {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  }
  else
    hEvtCount->Fill(14);
  // AliMultSelection *fMultSel = (AliMultSelection *)InputEvent()->FindListObject("MultSelection");
  // if (!fMultSel)
  // {
  //   AliError(Form("%s: Could not get AliMultSelection", GetName()));
  // }
  // else
  //   hEvtCount->Fill(15);
  if (!manager || !handler || !fAOD || !fUtils)
    return;
  hEvtCount->Fill(2);


  // Trigger
  //----------------------------
  if(TriggerIsOn){
  UInt_t mask = handler->IsEventSelected();
  bool isTrigselected = false;
  if (fTrigger.EqualTo("kMB"))
  isTrigselected = mask & AliVEvent::kMB;
  else if (fTrigger.EqualTo("kINT7"))
  isTrigselected = mask & AliVEvent::kINT7;
  else if (fTrigger.EqualTo("kINT7+kSemiCentral"))
  isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kSemiCentral);
  else if (fTrigger.EqualTo("kINT7+kCentral+kSemiCentral"))
  isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kCentral + AliVEvent::kSemiCentral);
  if (isTrigselected == false) return;
  hEvtCount->Fill(7);}

  //------------------
  // event-wise
  //------------------
  // runNumber
  runNum = fAOD->GetRunNumber();
  if(fPeriod.EqualTo("LHC15o")){      
    if (runNum != oldRunNum){
      if (!LoadCalibHistForThisRun()) return;
      oldRunNum = runNum;
    }   
  }
  else if(fPeriod.EqualTo("LHC18q")||fPeriod.EqualTo("LHC18r")){  
    if (runNum != oldRunNum){
      if (!LoadCalibHistForThisRun())
        return;
      oldRunNum = runNum;
    }
  }
  runNumBin = GetRunNumBin(runNum);
  if (runNumBin < 0)
    return;
  hRunNumBin->Fill(runNumBin);
  hEvtCount->Fill(3);
    

  // vertex
  const AliVVertex *vtTrc = fAOD->GetPrimaryVertex();
  const AliVVertex *vtSPD = fAOD->GetPrimaryVertexSPD();
  vtx[0] = (double)vtTrc->GetX();
  vtx[1] = (double)vtTrc->GetY();
  vtx[2] = (double)vtTrc->GetZ();
  hVxy[0]->Fill(vtx[0], vtx[1]);
  hVz[0]->Fill(vtx[2]);
  if (fabs(vtx[2]) > fVzCut) return;

  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = vtTrc->GetZ() - vtSPD->GetZ();
  double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  if (TMath::Abs(dz) > 0.2 || nsigTot > 10 || nsigTrc > 20)
    return;
  hVxy[1]->Fill(vtx[0], vtx[1]);
  hVz[1]->Fill(vtx[2]);

  for (int i = 0; i < 20; ++i)
  {
    if (vtx[2] > -10 + i * 1 && vtx[2] < -10 + (i + 1) * 1)
    {
      vzBin = i;
      break;
    }
  }
  if (vzBin < -990)
    return;
  hEvtCount->Fill(4);

  // centrality
  AliMultSelection *fMultSel = (AliMultSelection *)InputEvent()->FindListObject("MultSelection");
  cent = fMultSel->GetMultiplicityPercentile("V0M");
  centSPD1 = fMultSel->GetMultiplicityPercentile("CL1");
  hCentCorr[0]->Fill(cent, centSPD1);
  if (fabs(cent - centSPD1) > 7.5)
    return;
  if (cent < 10 || cent >= 50)
    return;
  hCentCorr[1]->Fill(cent, centSPD1);
  if (cent >= 0 && cent < 5)
    centBin = 0;
  else if (cent >= 5 && cent < 10)
    centBin = 1;
  else
    centBin = (int)cent / 10 + 1; // centbin
  hCent->Fill(cent);
  hEvtCount->Fill(5);

  // pileup 18q
  if (!RejectEvtTFFit())
    return;
  hEvtCount->Fill(6);

  ///---------------------------
  // ZDC Plane
  ///----------------------------
  if (!GetZDCPlaneLsFit())
    return;
  hEvtCount->Fill(19);

  //------------------
  //* loop trk
  //------------------

  vector<double> vecPhi;
  vector<int> vecCharge;
  vector<double> vecEta;
  vector<double> vecWeight;
  vector<double> vecPt;
  vector<double> vecPx;
  vector<double> vecPOI;

  // TPC QxQy
  double sumCosPos = 0.;
  double sumCosNeg = 0.;
  double sumSinPos = 0.;
  double sumSinNeg = 0.;
  int etaPos = 0;
  int etaNeg = 0;
//  double Psi2Qx = 0;
//  double Psi2Qy = 0;
  int nTrk = fAOD->GetNumberOfTracks();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk)
  {
    AliAODTrack *track = static_cast<AliAODTrack *>(fAOD->GetTrack(iTrk));
    if (!track)
    {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }

    double pt = track->Pt();
    double eta = track->Eta();
    double phi = track->Phi();
    int charge = track->Charge();
    int nhits = track->GetTPCNcls();
    double dedx = track->GetTPCsignal();
    double chi2 = track->Chi2perNDF();
    double px = track->Px();

    if (!track->TestFilterBit(fFilterBit))
      continue;

    int etaBin = -1;
    if (eta > 0.)
    {
      etaBin = 0;
    }
    else if (eta < 0.)
    {
      etaBin = 1;
    }
    hBeforePhi[centBin]->Fill(phi, etaBin);
    hBeforePhi[centBin]->Fill(phi, 2);
    hPt[0]->Fill(pt);
    hEta[0]->Fill(eta);
    hNhits[0]->Fill(nhits);

    if (pt < fPtMin || pt > fPtMax)
      continue;
    if (fabs(eta) > fEtaMax)
      continue;
    if (fabs(nhits) < fNhitsMin)
      continue;
    if (chi2 > fChi2Max)
      continue;
    if (dedx < fDeDxMin)
      continue;

    //------------------
    // NUE & NUA
    //------------------
    double weight = 1.;
    if (IsDoNUE)
    {
      double wEffi = GetNUECor(charge, pt);
      if (wEffi < 0)
        continue;
      else
        weight *= wEffi;
    }
    if (IsDoNUA)
    {
      double wAcc = GetNUACor(charge, phi, eta, vtx[2]);
      if (wAcc < 0)
        continue;
      else
        weight *= wAcc;
    }
    hAfterPhi[centBin]->Fill(phi, etaBin, weight);
    hAfterPhi[centBin]->Fill(phi, 2, weight);
    hPt[1]->Fill(pt, weight);
    hEta[1]->Fill(eta);
    hNhits[1]->Fill(nhits);
    hPDedx->Fill(track->P() * charge, dedx);
    // PID
    double fPOIBIN = 0;
    bool isItProttrk = CheckPIDofParticle(track, 3); // 3=proton
    bool isItkiontrk = CheckPIDofParticle(track, 2); // 2=kion
    bool isItpiontrk = CheckPIDofParticle(track, 1); // 1=pion
    if (charge > 0)
    {
      if (isItProttrk && !isItkiontrk && !isItpiontrk)
        fPOIBIN = 0.5;
      if (!isItProttrk && isItkiontrk && !isItpiontrk)
        fPOIBIN = 1.5;
      if (!isItProttrk && !isItkiontrk && isItpiontrk)
        fPOIBIN = 2.5;
    }
    else
    {
      if (isItProttrk && !isItkiontrk && !isItpiontrk)
        fPOIBIN = 3.5;
      if (!isItProttrk && isItkiontrk && !isItpiontrk)
        fPOIBIN = 4.5;
      if (!isItProttrk && !isItkiontrk && isItpiontrk)
        fPOIBIN = 5.5;
    }
    if (fPOIBIN == 0)
      continue;
    // PID QA
    if (abs(fPOIBIN - 0.5) < 1E-6 || abs(fPOIBIN - 3.5) < 1E-6)
    {
      float nSigTPC = fPID->NumberOfSigmasTPC(track, AliPID::kProton);
      float nSigTOF = fPID->NumberOfSigmasTOF(track, AliPID::kProton);
      float nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);
      fHist2ProtonSigTPC->Fill(pt, nSigTPC);
      fHist2ProtonSigTOF->Fill(pt, nSigRMS);
    }
    if (abs(fPOIBIN - 1.5) < 1E-6 || abs(fPOIBIN - 4.5) < 1E-6)
    {
      float nSigTPC = fPID->NumberOfSigmasTPC(track, AliPID::kKaon);
      float nSigTOF = fPID->NumberOfSigmasTOF(track, AliPID::kKaon);
      float nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);
      fHist2KionSigTPC->Fill(pt, nSigTPC);
      fHist2KionSigTOF->Fill(pt, nSigRMS);
    }
    if (abs(fPOIBIN - 2.5) < 1E-6 || abs(fPOIBIN - 5.5) < 1E-6)
    {
      float nSigTPC = fPID->NumberOfSigmasTPC(track, AliPID::kPion);
      float nSigTOF = fPID->NumberOfSigmasTOF(track, AliPID::kPion);
      float nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);
      fHist2PionSigTPC->Fill(pt, nSigTPC);
      fHist2PionSigTOF->Fill(pt, nSigRMS);
    }
    vecPhi.push_back(phi);
    vecCharge.push_back(charge);
    vecEta.push_back(eta);
    vecWeight.push_back(weight);
    vecPt.push_back(pt);
    vecPx.push_back(px);
    vecPOI.push_back(fPOIBIN);

    if (eta > 0.)
    {
      sumCosPos += weight * cos(phi);
      sumSinPos += weight * sin(phi);
      etaPos++;
    }
    else if (eta < 0.)
    {
      sumCosNeg += weight * cos(phi);
      sumSinNeg += weight * sin(phi);
      etaNeg++;
    }
  }
  if (etaNeg == 0 || etaPos == 0)
    return;

  // TPC Plane
  TVector2 Q1TPCPos;
  Q1TPCPos.Set(sumCosPos, sumSinPos);
  double psiPos = Q1TPCPos.Phi();
  pos1Plane->Fill(psiPos);
  TVector2 Q1TPCNeg;
  Q1TPCNeg.Set(sumCosNeg, sumSinNeg);
  double psiNeg = Q1TPCNeg.Phi();
  neg1Plane->Fill(psiNeg);
  hEvtCount->Fill(17);
  TComplex negQ(sumCosNeg, sumSinNeg);
  TComplex negQStar = TComplex::Conjugate(negQ);
  TComplex posQ(sumCosPos, sumSinPos);
  TComplex posQStar = TComplex::Conjugate(posQ);
  Psi_P->Fill(posQ.Theta());
  Psi_T->Fill(negQ.Theta());
  Psi_PT->Fill(posQ.Theta() - negQ.Theta());
  // Fill Resolution
  Res1Square->Fill(0.5, cos(psiPos - psiNeg));
  fHist2Psi1ZNCCent[centBin]->Fill(fPsi1ZNC);
  fHist2Psi1ZNACent[centBin]->Fill(fPsi1ZNA);
  fProfileZDCPsi1Correlation->Fill(centBin + 0.5, TMath::Cos(1 * (fPsi1ZNC - fPsi1ZNA)));
  fProfileZDCPsi2Correlation->Fill(centBin + 0.5, TMath::Cos(2 * (fPsi1ZNC - fPsi1ZNA)));
  TComplex ZDCQt(Qtx, Qty);
  TComplex ZDCQp(Qpx, Qpy);
  ZDCResQ->Fill(centBin + 0.5, (ZDCQt * TComplex::Conjugate(ZDCQp)).Re());
  ResQ->Fill(centBin + 0.5, (negQ * posQStar).Re());
  for (vector<double>::size_type iTrk = 0; iTrk < vecPhi.size(); iTrk++)
  {
    double phi = vecPhi[iTrk];
    int charge = vecCharge[iTrk];
    double eta = vecEta[iTrk];
    double weight = vecWeight[iTrk];
    double pt = vecPt[iTrk];
    double px = vecPx[iTrk];
    double iTrkpoi = vecPOI[iTrk];
    ptEta->Fill(eta, pt);

    // PID QA
    fHistPIDPt->Fill(iTrkpoi, pt, weight);
    fHistPIDEta->Fill(iTrkpoi, eta, weight);
    fHistPIDPhi->Fill(iTrkpoi, phi, weight);

    // v1pt TPC
    TComplex u(cos(phi), sin(phi));
    if (eta > 0.)
    {
      px_P[centBin]->Fill(iTrkpoi, eta, pt * ((u * (TComplex::Conjugate(posQ - u))).Re()), weight);
      px_T[centBin]->Fill(iTrkpoi, eta, pt * ((u * negQStar).Re()), weight);
      v1_p[centBin]->Fill(iTrkpoi, eta, (u * (TComplex::Conjugate(posQ - u))).Re(), weight);
      v1_t[centBin]->Fill(iTrkpoi, eta, (u * negQStar).Re(), weight);
      TPCcos_p[centBin]->Fill(iTrkpoi, eta, cos(phi - ((posQ - u).Theta())), weight);
      TPCcos_t[centBin]->Fill(iTrkpoi, eta, cos(phi - (negQ.Theta())), weight);

      if (charge > 0.)
      {
        px_P[centBin]->Fill(6.5, eta, pt * ((u * (TComplex::Conjugate(posQ - u))).Re()), weight);
        px_T[centBin]->Fill(6.5, eta, pt * ((u * negQStar).Re()), weight);
        v1_p[centBin]->Fill(6.5, eta, (u * (TComplex::Conjugate(posQ - u))).Re(), weight);
        v1_t[centBin]->Fill(6.5, eta, (u * negQStar).Re(), weight);
        TPCcos_p[centBin]->Fill(6.5, eta, cos(phi - ((posQ - u).Theta())), weight);
        TPCcos_t[centBin]->Fill(6.5, eta, cos(phi - (negQ.Theta())), weight);
      }
      else if (charge < 0.)
      {
        px_P[centBin]->Fill(7.5, eta, pt * ((u * (TComplex::Conjugate(posQ - u))).Re()), weight);
        px_T[centBin]->Fill(7.5, eta, pt * ((u * negQStar).Re()), weight);
        v1_p[centBin]->Fill(7.5, eta, (u * (TComplex::Conjugate(posQ - u))).Re(), weight);
        v1_t[centBin]->Fill(7.5, eta, (u * negQStar).Re(), weight);
        TPCcos_p[centBin]->Fill(7.5, eta, cos(phi - ((posQ - u).Theta())), weight);
        TPCcos_t[centBin]->Fill(7.5, eta, cos(phi - (negQ.Theta())), weight);
      }
    }
    else if (eta < 0.)
    {
      px_P[centBin]->Fill(iTrkpoi, eta, pt * ((u * posQStar).Re()), weight);
      px_T[centBin]->Fill(iTrkpoi, eta, pt * ((u * (TComplex::Conjugate(negQ - u))).Re()), weight);
      v1_p[centBin]->Fill(iTrkpoi, eta, (u * posQStar).Re(), weight);
      v1_t[centBin]->Fill(iTrkpoi, eta, (u * (TComplex::Conjugate(negQ - u))).Re(), weight);
      TPCcos_p[centBin]->Fill(iTrkpoi, eta, cos(phi - (posQ.Theta())), weight);
      TPCcos_t[centBin]->Fill(iTrkpoi, eta, cos(phi - ((negQ - u).Theta())), weight);

      if (charge > 0.)
      {
        px_P[centBin]->Fill(6.5, eta, pt * ((u * posQStar).Re()), weight);
        px_T[centBin]->Fill(6.5, eta, pt * ((u * (TComplex::Conjugate(negQ - u))).Re()), weight);
        v1_p[centBin]->Fill(6.5, eta, (u * posQStar).Re(), weight);
        v1_t[centBin]->Fill(6.5, eta, (u * (TComplex::Conjugate(negQ - u))).Re(), weight);
        TPCcos_p[centBin]->Fill(6.5, eta, cos(phi - (posQ.Theta())), weight);
        TPCcos_t[centBin]->Fill(6.5, eta, cos(phi - ((negQ - u).Theta())), weight);
      }
      else if (charge < 0.)
      {
        px_P[centBin]->Fill(7.5, eta, pt * ((u * posQStar).Re()), weight);
        px_T[centBin]->Fill(7.5, eta, pt * ((u * (TComplex::Conjugate(negQ - u))).Re()), weight);
        v1_p[centBin]->Fill(7.5, eta, (u * posQStar).Re(), weight);
        v1_t[centBin]->Fill(7.5, eta, (u * (TComplex::Conjugate(negQ - u))).Re(), weight);
        TPCcos_p[centBin]->Fill(7.5, eta, cos(phi - (posQ.Theta())), weight);
        TPCcos_t[centBin]->Fill(7.5, eta, cos(phi - ((negQ - u).Theta())), weight);
      }
    }

    ZDCpx_P[centBin]->Fill(iTrkpoi, eta, pt * ((u * TComplex::Conjugate(ZDCQp)).Re()), weight);
    ZDCpx_T[centBin]->Fill(iTrkpoi, eta, pt * ((u * TComplex::Conjugate(ZDCQt)).Re()), weight);
    ZDCv1_p[centBin]->Fill(iTrkpoi, eta, (u * TComplex::Conjugate(ZDCQp)).Re(), weight);
    ZDCv1_t[centBin]->Fill(iTrkpoi, eta, (u * TComplex::Conjugate(ZDCQt)).Re(), weight);
    ZDCcos_t[centBin]->Fill(iTrkpoi, eta, cos(phi - fPsi1ZNC), weight);
    ZDCcos_p[centBin]->Fill(iTrkpoi, eta, cos(phi - fPsi1ZNA), weight);
    
    if (abs(iTrkpoi - 0.5) < 1E-6 || abs(iTrkpoi - 3.5) < 1E-6)   Proton_EtaPhi[centBin]->Fill(eta, phi);
    
    if (abs(iTrkpoi - 1.5) < 1E-6 || abs(iTrkpoi - 4.5) < 1E-6)   Kion_EtaPhi[centBin]->Fill(eta, phi);

    if (abs(iTrkpoi - 2.5) < 1E-6 || abs(iTrkpoi - 5.5) < 1E-6)   Pion_EtaPhi[centBin]->Fill(eta, phi);

    if (charge > 0.)  {
      PosHadron_EtaPhi[centBin]->Fill(eta, phi);
      PosHadron_PhiPsi_p[centBin]->Fill(eta, phi-fPsi1ZNA);
      PosHadron_PhiPsi_t[centBin]->Fill(eta, phi-fPsi1ZNC);}
    else if (charge < 0.)  {
      NegHadron_EtaPhi[centBin]->Fill(eta, phi);
      NegHadron_PhiPsi_p[centBin]->Fill(eta, phi-fPsi1ZNA);
      NegHadron_PhiPsi_t[centBin]->Fill(eta, phi-fPsi1ZNC);}    

    if (charge > 0.)
    {
      ZDCpx_P[centBin]->Fill(6.5, eta, pt * ((u * TComplex::Conjugate(ZDCQp)).Re()), weight);
      ZDCpx_T[centBin]->Fill(6.5, eta, pt * ((u * TComplex::Conjugate(ZDCQt)).Re()), weight);
      ZDCv1_p[centBin]->Fill(6.5, eta, (u * TComplex::Conjugate(ZDCQp)).Re(), weight);
      ZDCv1_t[centBin]->Fill(6.5, eta, (u * TComplex::Conjugate(ZDCQt)).Re(), weight);
      ZDCv1_p_15o[centBin]->Fill(6.5, eta, (u * TComplex::Conjugate(ZDCQp)).Re(), weight);
      ZDCv1_t_15o[centBin]->Fill(6.5, eta, (u * TComplex::Conjugate(ZDCQt)).Re(), weight);
      ZDCcos_t[centBin]->Fill(6.5, eta, cos(phi - fPsi1ZNC), weight);
      ZDCcos_p[centBin]->Fill(6.5, eta, cos(phi - fPsi1ZNA), weight);
    }
    else if (charge < 0.)
    {
      ZDCpx_P[centBin]->Fill(7.5, eta, pt * ((u * TComplex::Conjugate(ZDCQp)).Re()), weight);
      ZDCpx_T[centBin]->Fill(7.5, eta, pt * ((u * TComplex::Conjugate(ZDCQt)).Re()), weight);
      ZDCv1_p[centBin]->Fill(7.5, eta, (u * TComplex::Conjugate(ZDCQp)).Re(), weight);
      ZDCv1_t[centBin]->Fill(7.5, eta, (u * TComplex::Conjugate(ZDCQt)).Re(), weight);
      ZDCv1_p_15o[centBin]->Fill(7.5, eta, (u * TComplex::Conjugate(ZDCQp)).Re(), weight);
      ZDCv1_t_15o[centBin]->Fill(7.5, eta, (u * TComplex::Conjugate(ZDCQt)).Re(), weight);
      ZDCcos_t[centBin]->Fill(7.5, eta, cos(phi - fPsi1ZNC), weight);
      ZDCcos_p[centBin]->Fill(7.5, eta, cos(phi - fPsi1ZNA), weight);
    }

  }
  hEvtCount->Fill(20);
  PostData(1, fOutputList); // stream the results the analysis of this event to
  PostData(2, fQAList);     // the output manager which will take care of writing
                            // it to a file
}

//---------------------------------------------------
int AliAnalysisTaskChargeV1::GetRunNumBin(int runNum)
{
  int runNumBin = -1;
  // 18q
  if (fPeriod.EqualTo("LHC15o")) {
      int runNumList[77] = {246994,246991,246989,246984,246982,246980,246948,246945,246928,246851,
                            246847,246846,246845,246844,246810,246809,246808,246807,246805,246804,
                            246766,246765,246763,246760,246759,246758,246757,246751,246750,246495,
                            246493,246488,246487,246434,246431,246428,246424,246276,246275,246272,
                            246271,246225,246222,246217,246185,246182,246181,246180,246178,246153,
                            246152,246151,246115,246113,246089,246053,246052,246049,246048,246042,
                            246037,246036,246012,246003,246001,245954,245952,245949,245923,245833,
                            245831,245829,245705,245702,245700,245692,245683};
    for (int i = 0; i < 77; ++i)
    {
      if (runNum == runNumList[i])
      {
        runNumBin = i;
        break;
      }
      else
        continue;
    }
  }
  else if (fPeriod.EqualTo("LHC18q")||fPeriod.EqualTo("LHC18r")) {
  int runNumList[214] = {296623, 296622, 296621, 296619, 296618, 296616, 296615, 296594, 296553, 296552,
                         296551, 296550, 296548, 296547, 296516, 296512, 296511, 296510, 296509, 296472,
                         296433, 296424, 296423, 296420, 296419, 296415, 296414, 296383, 296381, 296380,
                         296379, 296378, 296377, 296376, 296375, 296312, 296309, 296304, 296303, 296280,
                         296279, 296273, 296270, 296269, 296247, 296246, 296244, 296243, 296242, 296241,
                         296240, 296198, 296197, 296196, 296195, 296194, 296192, 296191, 296143, 296142,
                         296135, 296134, 296133, 296132, 296123, 296074, 296066, 296065, 296063, 296062,
                         296060, 296016, 295942, 295941, 295937, 295936, 295913, 295910, 295909, 295861,
                         295860, 295859, 295856, 295855, 295854, 295853, 295831, 295829, 295826, 295825,
                         295822, 295819, 295818, 295816, 295791, 295788, 295786, 295763, 295762, 295759,
                         295758, 295755, 295754, 295725, 295723, 295721, 295719, 295718, 295717, 295714,
                         295712, 295676, 295675, 295673, 295668, 295667, 295666, 295615, 295612, 295611,
                         295610, 295589, 295588, 295586, 295585, 297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512,
                         297483, 297479, 297452, 297451, 297450, 297446, 297442, 297441, 297415, 297414,
                         297413, 297406, 297405, 297380, 297379, 297372, 297367, 297366, 297363, 297336,
                         297335, 297333, 297332, 297317, 297311, 297310, 297278, 297222, 297221, 297218,
                         297196, 297195, 297193, 297133, 297132, 297129, 297128, 297124, 297123, 297119,
                         297118, 297117, 297085, 297035, 297031, 296966, 296941, 296938, 296935, 296934,
                         296932, 296931, 296930, 296903, 296900, 296899, 296894, 296852, 296851, 296850,
                         296848, 296839, 296838, 296836, 296835, 296799, 296794, 296793, 296790, 296787,
                         296786, 296785, 296784, 296781, 296752, 296694, 296693, 296691, 296690};
    for (int i = 0; i < 214; ++i)
    {
      if (runNum == runNumList[i])
      {
        runNumBin = i;
        break;
      }
      else
        continue;
    }
  }
  return runNumBin;
}

//---------------------------------------------------

bool AliAnalysisTaskChargeV1::RejectEvtTFFit()
{
  Float_t centV0M = -999;
  Float_t centCL1 = -999;
  Float_t centCL0 = -999;

  AliMultSelection *fMultSelection = (AliMultSelection *)InputEvent()->FindListObject("MultSelection");
  if (!fMultSelection)
  {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }
  centV0M = (Float_t)cent;
  centCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  Int_t nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets *aodTrkl = (AliAODTracklets *)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTrk = 0;

  for (Int_t it = 0; it < nTracks; it++)
  {
    AliAODTrack *aodTrk = (AliAODTrack *)fAOD->GetTrack(it);
    if (!aodTrk)
    {
      delete aodTrk;
      continue;
    }
    if (aodTrk->TestFilterBit(32))
    {
      if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
    }
  }

  fHist2DMultCentQA[0]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  AliAODVZERO *aodV0 = fAOD->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  // pile-up cuts
  if (centCL0 < fCenCutLowPU->Eval(centV0M))
    return false;
  if (centCL0 > fCenCutHighPU->Eval(centV0M))
    return false;
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls))
    return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot))
    return false;
  if (Float_t(multTrk) < fMultCutPU->Eval(centV0M))
    return false;
  if (((AliAODHeader *)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0)
    return false;
  if (fAOD->IsIncompleteDAQ())
    return false;

  fHist2DMultCentQA[1]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}

//_____________________________________________________________________________
void AliAnalysisTaskChargeV1::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________

bool AliAnalysisTaskChargeV1::LoadCalibHistForThisRun()
{
  if(fPeriod.EqualTo("LHC15o"))
  {    
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[0][i] = (TH1D*)(fListZDCCalib->FindObject(Form("fZNCTower[%d][%d]",runNum,i)));
      fTowerGainEq[1][i] = (TH1D*)(fListZDCCalib->FindObject(Form("fZNATower[%d][%d]",runNum,i)));
      if (!fTowerGainEq[0][i] || !fTowerGainEq[1][i]) return false;
    }
    std::cout << "\n ===========> Info:: ZDC Channel Weights Found for Run " << runNum << std::endl;
    return true;
  }

     
  // 18q/r NUA
  else if(fPeriod.EqualTo("LHC18q")||fPeriod.EqualTo("LHC18r"))
  {
    if (IsDoNUA)
    {
      hCorrectNUAPos->Reset();
      hCorrectNUANeg->Reset();
      hCorrectNUAPos = (TH3F *)fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d", 0, runNum));
      hCorrectNUANeg = (TH3F *)fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d", 0, runNum));
      if (!hCorrectNUAPos)
        return false;
      if (!hCorrectNUANeg)
        return false;
    }
    // ZDC
    fHZDCCparameters = (TH1D *)(fListZDCCalib->FindObject(Form("Run %d", runNum))->FindObject(Form("fZDCCparameters[%d]", runNum)));
    fHZDCAparameters = (TH1D *)(fListZDCCalib->FindObject(Form("Run %d", runNum))->FindObject(Form("fZDCAparameters[%d]", runNum)));
    if (fHZDCCparameters && fHZDCAparameters)
      std::cout << "\n ===========> Info:: ZDC Channel Weights Found for Run " << runNum << std::endl;

    if (!fHZDCCparameters)
      return false;
    if (!fHZDCAparameters)
      return false;

    return true;
  }
}

//---------------------------------------------------

double AliAnalysisTaskChargeV1::GetNUACor(int charge, double phi, double eta, double vz)
{
  double weightNUA = 1;

  if (vzBin < 0 || centBin < 0 || runNum < 0)
    return -1;

  // Rihan and Protty 's NUA Results
  if (charge > 0)
  {
    if (!hCorrectNUAPos)
      return -1;
    int iBinNUA = hCorrectNUAPos->FindBin(vz, phi, eta);
    if (hCorrectNUAPos->GetBinContent(iBinNUA) > 0)
      weightNUA = (double)hCorrectNUAPos->GetBinContent(iBinNUA);
    return weightNUA;
  }
  else if (charge < 0)
  {
    if (!hCorrectNUANeg)
      return -1;
    int iBinNUA = hCorrectNUANeg->FindBin(vz, phi, eta);
    if (hCorrectNUANeg->GetBinContent(iBinNUA) > 0)
      weightNUA = (double)hCorrectNUANeg->GetBinContent(iBinNUA);
    return weightNUA;
  }
  // In Rihan and Protty 's NUA results, the phi distribution is independent on centrality and particle charge

  return weightNUA;
}

//---------------------------------------------------

double AliAnalysisTaskChargeV1::GetNUECor(int charge, double pt)
{
  double weightNUE = 1;

  if (charge > 0)
  {
    //   hNUEweightPlus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgPos");
    if (!hNUEweightPlus)
      return -1;
    int ptBin = hNUEweightPlus->GetXaxis()->FindBin(pt);
    if (hNUEweightPlus->GetBinContent(ptBin) > 0)
    {
      weightNUE = 1. / hNUEweightPlus->GetBinContent(ptBin);
    }
    else
      return -1;
  }
  if (charge < 0)
  {
    //   hNUEweightMinus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgNeg");
    if (!hNUEweightMinus)
      return -1;
    int ptBin = hNUEweightMinus->GetXaxis()->FindBin(pt);
    if (hNUEweightMinus->GetBinContent(ptBin) > 0)
    {
      weightNUE = 1. / hNUEweightMinus->GetBinContent(ptBin);
    }
    else
      return -1;
  }

  return weightNUE;
}

//-----------------------------------------------------

bool AliAnalysisTaskChargeV1::GetZDCPlaneLsFit(){
  AliAODZDC *fZDC = fAOD->GetZDCData();
  float fCent = cent;

  if(fPeriod.EqualTo("LHC15o"))
  {
    const Double_t * towZNCraw = fZDC->GetZNCTowerEnergy();
    const Double_t * towZNAraw = fZDC->GetZNATowerEnergy();

    Double_t Enucl = (runNum < 209122 ? 1380. : 2511.);
    Double_t xyZNC[2]={0.,0.}, xyZNA[2]={0.,0.};
    Double_t towZNC[5]={0.}, towZNA[5]={0.};


    Double_t ZNCcalib=1., ZNAcalib=1.;

    for(Int_t i=0; i<5; i++) {
      if(fTowerGainEq[0][i]) towZNC[i] = towZNCraw[i]*fTowerGainEq[0][i]->GetBinContent(fTowerGainEq[0][i]->FindBin(fCent));
      if(fTowerGainEq[1][i]) towZNA[i] = towZNAraw[i]*fTowerGainEq[1][i]->GetBinContent(fTowerGainEq[1][i]->FindBin(fCent));
      // if(fResetNegativeZDC) {
      //   if(towZNC[i]<0.) towZNC[i] = 0.;
      //   if(towZNA[i]<0.) towZNA[i] = 0.;
      }
        
    if(runNum>=245829) towZNA[2] = 0.;
    Double_t zncEnergy=0., znaEnergy=0.;
    for(Int_t i=0; i<5; i++){
      zncEnergy += towZNC[i];
      znaEnergy += towZNA[i];
    }
    if(runNum>=245829) znaEnergy *= 8./7.;

    Double_t energyZNC = ((AliVAODHeader*)fAOD->GetHeader())->GetZDCN1Energy();
    Double_t energyZNA = ((AliVAODHeader*)fAOD->GetHeader())->GetZDCN2Energy();

    const Double_t x[4] = {-1.75, 1.75, -1.75, 1.75};
    const Double_t y[4] = {-1.75, -1.75, 1.75, 1.75};
    Double_t numXZNC=0., numYZNC=0., denZNC=0., cZNC, wZNC, EZNC, SumEZNC=0.;
    Double_t numXZNA=0., numYZNA=0., denZNA=0., cZNA, wZNA, EZNA, SumEZNA=0., BadChOr;
    Bool_t fAllChONZNC=kTRUE, fAllChONZNA=kTRUE;

    for(Int_t i=0; i<4; i++)
    {
      // get energy
      EZNC = towZNC[i+1];
      fhZNSpectra->Fill(fCent,i+0.5,EZNC);
  //   fhZNSpectraRbR[RunBin]->Fill(fCent,i+0.5,EZNC);
      if(fUseZDCSpectraCorr && EZNC>0.) 
      {
        Double_t mu1 = SpecCorMu1[i]->Interpolate(fCent);
        Double_t mu2 = SpecCorMu2[i]->Interpolate(fCent);
        Double_t av = SpecCorAv[i]->Interpolate(fCent);
        Double_t cor1 = SpecCorSi[i]->Interpolate(fCent);
        EZNC = exp( (log(EZNC) - mu1 + mu2*cor1)/cor1 ) + av;
        fhZNSpectraCor->Fill(fCent,i+0.5,EZNC);
      }
      if(fUseZDCSpectraCorr && EZNC<=0.) fAllChONZNC=kFALSE;

      SumEZNC += EZNC;

      // build centroid
      wZNC = TMath::Power(EZNC, fZDCGainAlpha);
      numXZNC += x[i]*wZNC;
      numYZNC += y[i]*wZNC;
      denZNC += wZNC;
      fhZNSpectraPow->Fill(fCent,i+0.5,wZNC);

      // get energy
      if(i==1) 
      {
        EZNA = towZNA[0]-towZNA[1]-towZNA[3]-towZNA[4];
        if(fUseBadTowerCalib && fBadTowerCalibHist[int(fCent)]) 
        {
          EZNA = GetBadTowerResp(EZNA, fBadTowerCalibHist[int(fCent)]);
        }
      } 
      else {
        EZNA = towZNA[i+1];}
      fhZNSpectra->Fill(fCent,i+4.5,EZNA);

      if(fUseZDCSpectraCorr && EZNA>0.) 
      {
        Double_t mu1 = SpecCorMu1[i+4]->Interpolate(fCent);
        Double_t mu2 = SpecCorMu2[i+4]->Interpolate(fCent);
        Double_t av = SpecCorAv[i+4]->Interpolate(fCent);
        Double_t cor1 = SpecCorSi[i+4]->Interpolate(fCent);
        EZNA = exp( (log(EZNA) - mu1 + mu2*cor1)/cor1 ) + av;
        fhZNSpectraCor->Fill(fCent,i+4.5,EZNA);
      }

      if(fUseZDCSpectraCorr && EZNA<=0.) fAllChONZNA=kFALSE;
      SumEZNA += EZNA;

      // build centroid
      wZNA = TMath::Power(EZNA, fZDCGainAlpha);
      numXZNA += x[i]*wZNA;
      numYZNA += y[i]*wZNA;
      denZNA += wZNA;
      fhZNSpectraPow->Fill(fCent,i+4.5,wZNA);
      }
    // store distribution for unfolding
    if(runNum<245829) 
    {
      Double_t recoE = towZNA[0]-towZNA[1]-towZNA[3]-towZNA[4];
      Double_t trueE = towZNA[2];
      fhZNBCCorr->Fill(fCent,trueE,recoE);
    }
    if(denZNC>0.)
    {
      Double_t nSpecnC = SumEZNC/Enucl;
      cZNC = 1.89358-0.71262/(nSpecnC+0.71789);
      xyZNC[0] = cZNC*numXZNC/denZNC;
      xyZNC[1] = cZNC*numYZNC/denZNC;
      denZNC *= cZNC;
    }
    else{ xyZNC[0] = xyZNC[1] = 0.;}
    if(denZNA>0.)
    {
      Double_t nSpecnA = SumEZNA/Enucl;
      cZNA = 1.89358-0.71262/(nSpecnA+0.71789);
      xyZNA[0] = cZNA*numXZNA/denZNA;
      xyZNA[1] = cZNA*numYZNA/denZNA;
      denZNA *= cZNA;
    }
    else{ xyZNA[0] = xyZNA[1] = 0.;}
    
    if(!fAllChONZNC) denZNC=-1.;
    if(!fAllChONZNA) denZNA=-1.;

    if(denZNC>0. && pow(xyZNC[0]*xyZNC[0]+xyZNC[1]*xyZNC[1],0.5)>1.E-6) 
    {
      double psiZNCGE = GetEventPlane(xyZNC[0], xyZNC[1], 1);
      fPsi1ZNC = psiZNCGE;
      fHist2CalibPsi1ZNCCent[0]->Fill(psiZNCGE);
    }  
    else {return false;}

    if(denZNA>0. && pow(xyZNA[0]*xyZNA[0]+xyZNA[1]*xyZNA[1],0.5)>1.E-6) 
    {
      double psiZNAGE = GetEventPlane(xyZNA[0], xyZNA[1], 1);
      fPsi1ZNA = psiZNAGE;
      fHist2CalibPsi1ZNACent[0]->Fill(psiZNAGE);
    }
    else {return false;}
    Qtx = xyZNC[0];
    Qty = xyZNC[1];
    Qpx = xyZNA[0];
    Qpy = xyZNA[1];

    fProfileZNCQxCent[0]->Fill(fCent, Qtx);
    fProfileZNCQyCent[0]->Fill(fCent, Qty);
    fProfileZNAQxCent[0]->Fill(fCent, Qpx);
    fProfileZNAQyCent[0]->Fill(fCent, Qpy);

    fProfileZDCQxAQxCCent[0]->Fill(fCent, Qpx * Qtx);
    fProfileZDCQxAQyCCent[0]->Fill(fCent, Qpx * Qty);
    fProfileZDCQyAQxCCent[0]->Fill(fCent, Qpy * Qtx);
    fProfileZDCQyAQyCCent[0]->Fill(fCent, Qpy * Qty);
    return true;
  }
  else if(fPeriod.EqualTo("LHC18r")||fPeriod.EqualTo("LHC18q"))
  {
    if (!fZDC)
      return false;
    const double *fZNATowerRawAOD = fZDC->GetZNATowerEnergy();
    const double *fZNCTowerRawAOD = fZDC->GetZNCTowerEnergy();
    for (int iTower = 0; iTower < 5; iTower++)
    {
      if (fZNATowerRawAOD[iTower] < 1.e-6)
        return false;
      if (fZNCTowerRawAOD[iTower] < 1.e-6)
        return false;
    }

    double towZNCraw1GainEq = 0, towZNCraw2GainEq = 0, towZNCraw3GainEq = 0, towZNCraw4GainEq = 0;
    towZNCraw1GainEq = fZNCTowerRawAOD[1] * fHZDCCparameters->GetBinContent(1);
    towZNCraw2GainEq = fZNCTowerRawAOD[2] * fHZDCCparameters->GetBinContent(2);
    towZNCraw3GainEq = fZNCTowerRawAOD[3] * fHZDCCparameters->GetBinContent(3);
    towZNCraw4GainEq = fZNCTowerRawAOD[4] * fHZDCCparameters->GetBinContent(4);

    double towZNAraw1GainEq = 0, towZNAraw2GainEq = 0, towZNAraw3GainEq = 0, towZNAraw4GainEq = 0;
    towZNAraw1GainEq = fZNATowerRawAOD[1] * fHZDCAparameters->GetBinContent(1);
    towZNAraw2GainEq = fZNATowerRawAOD[2] * fHZDCAparameters->GetBinContent(2);
    towZNAraw3GainEq = fZNATowerRawAOD[3] * fHZDCAparameters->GetBinContent(3);
    towZNAraw4GainEq = fZNATowerRawAOD[4] * fHZDCAparameters->GetBinContent(4);

    const double xZDCC[4] = {-1, 1, -1, 1}; // directional vector
    const double yZDCC[4] = {-1, -1, 1, 1};
    const double xZDCA[4] = {1, -1, 1, -1};
    const double yZDCA[4] = {-1, -1, 1, 1};

    double towZNC[5] = {fZNCTowerRawAOD[0], towZNCraw1GainEq, towZNCraw2GainEq, towZNCraw3GainEq, towZNCraw4GainEq};
    double towZNA[5] = {fZNATowerRawAOD[0], towZNAraw1GainEq, towZNAraw2GainEq, towZNAraw3GainEq, towZNAraw4GainEq};

    // double EZNC = 0;
    double wZNC = 0, denZNC = 0, numXZNC = 0, numYZNC = 0;
    // double EZNA = 0;
    double wZNA = 0, denZNA = 0, numXZNA = 0, numYZNA = 0;

    for (int i = 0; i < 4; i++)
    {
      // ZNC
      // get energy
      // EZNC = towZNC[i+1];
      // build ZDCC centroid
      wZNC = TMath::Max(0., 4.0 + TMath::Log(towZNC[i + 1] / fZNCTowerRawAOD[0]));
      numXZNC += xZDCC[i] * wZNC;
      numYZNC += yZDCC[i] * wZNC;
      denZNC += wZNC;

      // ZNA part
      // get energy
      // EZNA = towZNA[i+1];
      // build ZDCA centroid
      wZNA = TMath::Max(0., 4.0 + TMath::Log(towZNA[i + 1] / fZNATowerRawAOD[0]));
      numXZNA += xZDCA[i] * wZNA;
      numYZNA += yZDCA[i] * wZNA;
      denZNA += wZNA;
    }
    if (fabs(denZNC) < 1.e-6)
        return false;
    if (fabs(denZNA) < 1.e-6)
      return false;

    double ZDCCxPosFromLogWeight = numXZNC / denZNC;
    double ZDCCyPosFromLogWeight = numYZNC / denZNC;
    double ZDCAxPosFromLogWeight = numXZNA / denZNA;
    double ZDCAyPosFromLogWeight = numYZNA / denZNA;

    // QA
    fProfileZNCQxCent[0]->Fill(fCent, ZDCCxPosFromLogWeight);
    fProfileZNCQyCent[0]->Fill(fCent, ZDCCyPosFromLogWeight);
    fProfileZNAQxCent[0]->Fill(fCent, ZDCAxPosFromLogWeight);
    fProfileZNAQyCent[0]->Fill(fCent, ZDCAyPosFromLogWeight);

    fProfileZDCQxAQxCCent[0]->Fill(fCent, ZDCAxPosFromLogWeight * ZDCCxPosFromLogWeight);
    fProfileZDCQxAQyCCent[0]->Fill(fCent, ZDCAxPosFromLogWeight * ZDCCyPosFromLogWeight);
    fProfileZDCQyAQxCCent[0]->Fill(fCent, ZDCAyPosFromLogWeight * ZDCCxPosFromLogWeight);
    fProfileZDCQyAQyCCent[0]->Fill(fCent, ZDCAyPosFromLogWeight * ZDCCyPosFromLogWeight);

    double psiZNCGE = GetEventPlane(ZDCCxPosFromLogWeight, ZDCCyPosFromLogWeight, 1);
    double psiZNAGE = GetEventPlane(ZDCAxPosFromLogWeight, ZDCAyPosFromLogWeight, 1);
    if (TMath::IsNaN(psiZNCGE) || TMath::IsNaN(psiZNAGE))
      return false;

    fHist2CalibPsi1ZNCCent[0]->Fill(psiZNCGE);
    fHist2CalibPsi1ZNACent[0]->Fill(psiZNAGE);

    // Recenter
    double ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(5) * fCent + fHZDCCparameters->GetBinContent(6) * vtx[0] + fHZDCCparameters->GetBinContent(7) * vtx[1] + fHZDCCparameters->GetBinContent(8) * vtx[2] + fHZDCCparameters->GetBinContent(9);
    double ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(10) * fCent + fHZDCCparameters->GetBinContent(11) * vtx[0] + fHZDCCparameters->GetBinContent(12) * vtx[1] + fHZDCCparameters->GetBinContent(13) * vtx[2] + fHZDCCparameters->GetBinContent(14);
    double ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(5) * fCent + fHZDCAparameters->GetBinContent(6) * vtx[0] + fHZDCAparameters->GetBinContent(7) * vtx[1] + fHZDCAparameters->GetBinContent(8) * vtx[2] + fHZDCAparameters->GetBinContent(9);
    double ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(10) * fCent + fHZDCAparameters->GetBinContent(11) * vtx[0] + fHZDCAparameters->GetBinContent(12) * vtx[1] + fHZDCAparameters->GetBinContent(13) * vtx[2] + fHZDCAparameters->GetBinContent(14);

    double QxZNC = ZDCCxPosFromLogWeight - ZDCCAvgxPosFromVtxFit;
    double QyZNC = ZDCCyPosFromLogWeight - ZDCCAvgyPosFromVtxFit;

    double QxZNA = ZDCAxPosFromLogWeight - ZDCAAvgxPosFromVtxFit;
    double QyZNA = ZDCAyPosFromLogWeight - ZDCAAvgyPosFromVtxFit;

      // QA
    fProfileZNCQxCent[1]->Fill(fCent, QxZNC);
    fProfileZNCQyCent[1]->Fill(fCent, QyZNC);
    fProfileZNAQxCent[1]->Fill(fCent, QxZNA);
    fProfileZNAQyCent[1]->Fill(fCent, QyZNA);

    fProfileZDCQxAQxCCent[1]->Fill(fCent, QxZNA * QxZNC);
    fProfileZDCQxAQyCCent[1]->Fill(fCent, QxZNA * QyZNC);
    fProfileZDCQyAQxCCent[1]->Fill(fCent, QyZNA * QxZNC);
    fProfileZDCQyAQyCCent[1]->Fill(fCent, QyZNA * QyZNC);

    double psiZNCRC = GetEventPlane(QxZNC, QyZNC, 1);
    double psiZNARC = GetEventPlane(QxZNA, QyZNA, 1);
    if (TMath::IsNaN(psiZNCRC) || TMath::IsNaN(psiZNARC))
      return false;

    fHist2CalibPsi1ZNCCent[1]->Fill(psiZNCRC);
    fHist2CalibPsi1ZNACent[1]->Fill(psiZNARC);
    Qtx = QxZNC;
    Qty = QyZNC;
    Qpx = QxZNA;
    Qpy = QyZNA;
    fPsi1ZNC = psiZNCRC;
    fPsi1ZNA = psiZNARC;

    return true;
  }
}

//--------------------------------------------
inline double AliAnalysisTaskChargeV1::GetEventPlane(double qx, double qy, double harmonic)
{
  double psi = (1. / harmonic) * TMath::ATan2(qy, qx);
    return psi;
}

//---------------------------------------------

bool AliAnalysisTaskChargeV1::CheckPIDofParticle(AliAODTrack *ftrack, int pidToCheck)
{
  if (pidToCheck == 0)
    return kTRUE; //// Charge Particles do not need PID check

  if (!fPID)
  {
    Printf("\n Could Not access PIDResponse Task, please Add the Task...\n return with kFALSE pid\n");
    return kFALSE;
  }

  /// Rihan todo: To set the low pT cuts for nSigmaTPC from AddTaskMacro!
  /// Although someone barely needs to change it given the purity..

  float nSigTPC = 0, nSigTOF = 0, nSigRMS = 0;
  double trkPtPID = ftrack->Pt();

  /// Pion =>
  if (pidToCheck == 1)
  {
    nSigTPC = fPID->NumberOfSigmasTPC(ftrack, AliPID::kPion); // Some warning show here (***TDatabasePDG::AddParicle: particle with PDGcode = 3124 already defind),I don't understand what happended. --chunzheng
    nSigTOF = fPID->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);

    if (trkPtPID <= 0.5 && TMath::Abs(nSigTPC) <= fNSigmaTPCCut)
      return true;
    if (trkPtPID > 0.5 && trkPtPID < 2.0 && TMath::Abs(nSigRMS) <= fNSigmaTOFCut)
      return true;
    return false;
  }
  /// Kaon =>
  else if (pidToCheck == 2)
  {
    nSigTPC = fPID->NumberOfSigmasTPC(ftrack, AliPID::kKaon);
    nSigTOF = fPID->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
    nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);

    if (trkPtPID <= 0.45 && TMath::Abs(nSigTPC) <= fNSigmaTPCCut)
      return true;
    if (trkPtPID > 0.45 && trkPtPID < 2.0 && TMath::Abs(nSigRMS) <= fNSigmaTOFCut)
      return true;
    return false;
  }
  /// proton =>
  else if (pidToCheck == 3)
  { ///
    nSigTPC = fPID->NumberOfSigmasTPC(ftrack, AliPID::kProton);
    nSigTOF = fPID->NumberOfSigmasTOF(ftrack, AliPID::kProton);
    nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);

    if (trkPtPID <= 0.6 && TMath::Abs(nSigTPC) <= fNSigmaTPCCut)
      return true;
    if (trkPtPID > 0.6 && trkPtPID < 2.0 && TMath::Abs(nSigRMS) <= fNSigmaTOFCut)
      return true;
    return false;
  }
  else
  {
    Printf("\n -Ve number not allowed! Choose among: 0,1,2,3 (Charge Pion, Kaon, Proton)\n return with kFALSE \n");
    return false;
  }
}

//------------------------------
Double_t AliAnalysisTaskChargeV1::GetBadTowerResp(Double_t Et, TH2D* BadTowerCalibHist)
{
  Double_t EtC = BadTowerCalibHist->ProjectionY("",BadTowerCalibHist->GetXaxis()->FindBin(Et),BadTowerCalibHist->GetXaxis()->FindBin(Et))->GetRandom();
  return EtC;
}
