/**************************************************************************
 * Copyright(c) 1998-2023, ALICE Experiment at CERN, All rights reserved. *
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

//--------------------------------------------------------------------------------
// Feeddown study of CVE analysis
// Contributor: Chunzheng Wang, <chunzheng.wang@cern.ch>, Shanghai
//--------------------------------------------------------------------------------

#include <sys/time.h>
#include <cstdlib>
#include <cmath>

// ROOT classes
#include "TChain.h"
#include "TF1.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TSpline.h"
#include "TString.h"
// Alice analysis base class
#include "AliAnalysisTaskSE.h"
// Alice analysis additional classes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
// Alice AOD classes
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODVZERO.h"
#include "AliAODVertex.h"
// Alice classes
#include "AliEventCuts.h"
#include "AliEventplane.h"
// Alice MC classes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
// Alice PID classes
#include "AliAnalysisTaskProtonPtDCA.h"
#include "AliPIDResponse.h"


ClassImp(AliAnalysisTaskProtonPtDCA);



//---------------------------------------------------
AliAnalysisTaskProtonPtDCA::AliAnalysisTaskProtonPtDCA()
    : AliAnalysisTaskSE(),
      isNarrowDcaCuts768(true),
      isProtonCustomizedDCACut(true),
      isTightPileUp(false),

      fTrigger("kINT7+kSemiCentral"),
      fPeriod("LHC18q"),
      fVzCut(10.0),
      fFilterBit(768),
      fNclsCut(70),
      fChi2Max(2.5),
      fChi2Min(0.1),

      fNSigmaTPC(3.0),
      fNSigmaRMS(3.0),

      fAOD(nullptr),
      fPIDResponse(nullptr),
      fMultSel(nullptr),
      fRunNum(-999),
      fOldRunNum(-999),
      fRunNumBin(-999),
      fCent(-999),
      fSPDCutPU(nullptr),
      fV0CutPU(nullptr),
      fCenCutLowPU(nullptr),
      fCenCutHighPU(nullptr),
      fMultCutPU(nullptr),
      fQAList(nullptr),
      h2ProtonPtDcaXY(nullptr),
      h2ProtonPtDcaZ(nullptr),
      h2AntiProtonPtDcaXY(nullptr),
      h2AntiProtonPtDcaZ(nullptr)
       {
  for (auto& v : fVertex) v = 0;
}

//---------------------------------------------------
AliAnalysisTaskProtonPtDCA::AliAnalysisTaskProtonPtDCA(const char* name)
    : AliAnalysisTaskSE(name),
      isNarrowDcaCuts768(true),
      isProtonCustomizedDCACut(true),
      isTightPileUp(false),

      fTrigger("kINT7+kSemiCentral"),
      fPeriod("LHC18q"),
      fVzCut(10.0),
      fFilterBit(768),
      fNclsCut(70),
      fChi2Max(2.5),
      fChi2Min(0.1),

      fNSigmaTPC(3.0),
      fNSigmaRMS(3.0),

      fAOD(nullptr),
      fPIDResponse(nullptr),
      fMultSel(nullptr),
      fRunNum(-999),
      fOldRunNum(-999),
      fRunNumBin(-999),
      fCent(-999),
      fSPDCutPU(nullptr),
      fV0CutPU(nullptr),
      fCenCutLowPU(nullptr),
      fCenCutHighPU(nullptr),
      fMultCutPU(nullptr),
      fQAList(nullptr),
      h2ProtonPtDcaXY(nullptr),
      h2ProtonPtDcaZ(nullptr),
      h2AntiProtonPtDcaXY(nullptr),
      h2AntiProtonPtDcaZ(nullptr)
       {
  for (auto& v : fVertex) v = 0;

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//------------------------------------------------

AliAnalysisTaskProtonPtDCA::~AliAnalysisTaskProtonPtDCA() {
  // Destructor
  // histograms are in the output list and deleted when the output
  if (fQAList) delete fQAList;
}

//---------------------------------------------------

void AliAnalysisTaskProtonPtDCA::Terminate(Option_t*) {
  // Terminate loop
  Printf("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskProtonPtDCA::UserCreateOutputObjects() {

  ////////////////////////
  // Pile up Function
  ////////////////////////
  // Rihan 18q/r Pile-up function
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    double parV0[8] = {43.8011,     0.822574,    8.49794e-02,  1.34217e+02,
                       7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
    double parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
    double parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};

    fSPDCutPU = std::unique_ptr<TF1>(new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000));

    fV0CutPU = std::unique_ptr<TF1>(
        new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000));
    fV0CutPU->SetParameters(parV0);

    fCenCutLowPU = std::unique_ptr<TF1>(
        new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100));
    fCenCutLowPU->SetParameters(parV0CL0);

    fCenCutHighPU = std::unique_ptr<TF1>(
        new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100));
    fCenCutHighPU->SetParameters(parV0CL0);

    fMultCutPU = std::unique_ptr<TF1>(
        new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90));
    fMultCutPU->SetParameters(parFB32);
  }

  //------------------
  // QA
  //------------------
  fQAList = new TList();
  fQAList->SetName("fQAList");
  fQAList->SetOwner(kTRUE);

  h2ProtonPtDcaXY = new TH2D("h_pt_dcaxy_proton","h_pt_dcaxy_proton", 23, 0.4, 5, 300, 0, 2);
  h2ProtonPtDcaZ = new TH2D("h_pt_dcaz_proton","h_pt_dcaz_proton", 23, 0.4, 5, 300, 0, 2);
  h2AntiProtonPtDcaXY = new TH2D("h_pt_dcaxy_antiproton","h_pt_dcaxy_antiproton", 23, 0.4, 5, 300, 0, 2);
  h2AntiProtonPtDcaZ = new TH2D("h_pt_dcaz_antiproton","h_pt_dcaz_antiproton", 23, 0.4, 5, 300, 0, 2);

  fQAList->Add(h2ProtonPtDcaXY);
  fQAList->Add(h2ProtonPtDcaZ);
  fQAList->Add(h2AntiProtonPtDcaXY);
  fQAList->Add(h2AntiProtonPtDcaZ);


  PostData(1, fQAList);
  if (fDebug) Printf("Post fResultsList Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskProtonPtDCA::UserExec(Option_t*) {

  if (fDebug) Printf("===============================We are in UserExec!!!===================================");
  //----------------------------
  // Handle
  //----------------------------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  }

  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  }

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  }

  fPIDResponse = handler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  }

  fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  if (!fMultSel) {
    AliError(Form("%s: Could not get AliMultSelection", GetName()));
  }

  if (!manager || !handler || !fAOD || !fPIDResponse || !fMultSel) return;
  if (fDebug) Printf("Handles done!");

  //----------------------------
  // Trigger
  //----------------------------
  unsigned int mask = handler->IsEventSelected();
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
  if (fDebug) Printf("trigger done!");

  //----------------------------
  // Vertex
  //----------------------------
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  fVtx->GetXYZ(fVertex.data());
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  if (fabs(fVertex[0]) < 1e-6 || fabs(fVertex[1]) < 1e-6 || fabs(fVertex[2]) < 1e-6) return;
  if (fabs(fVertex[2]) > fVzCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors() < 1) return;
  if (fDebug) Printf("vertex done!");

  //----------------------------
  // Centrality
  //----------------------------
  float centV0M = fMultSel->GetMultiplicityPercentile("V0M");
  float centTRK = fMultSel->GetMultiplicityPercentile("TRK");
  float centSPD0 = fMultSel->GetMultiplicityPercentile("CL0");
  float centSPD1 = fMultSel->GetMultiplicityPercentile("CL1");

  // we use centV0M as the default centrality
  fCent = centV0M;
  if (fabs(fCent - centSPD1) > 7.5) return;
  if (fCent < 0 || fCent >= 80) return;

  // PF-Preview comment
  if (fCent > 30 && fCent < 50) {
    if (!(mask & (AliVEvent::kINT7 + AliVEvent::kSemiCentral))) return;
  } else {
    if (!(mask & AliVEvent::kINT7)) return;
  }

  if (fDebug) Printf("centrality done!");

  //----------------------------
  // Pile up
  //----------------------------
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    if (!RejectEvtTFFit(centSPD0)) return;
  }
  if (fDebug) Printf("pile-up done!");
  //----------------------------
  // Loop Tracks / Fill Vectors
  //----------------------------
  // must loop tracks becasue we need the TPC plane
  if (!LoopTracks()) return;
  if (fDebug) Printf("Loop Tracks done!");
  //------------------
  // Post output data.
  //------------------
  PostData(1, fQAList);
  if (fDebug) Printf("analysis done!");
}

//---------------------------------------------------

bool AliAnalysisTaskProtonPtDCA::LoopTracks() {
  int nTrks = fAOD->GetNumberOfTracks();
  if (nTrks < 10) return false;
  for (int iTrk = 0; iTrk < nTrks; ++iTrk) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrk);
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
    if (!track->TestFilterBit(fFilterBit)) continue;
    if (!AcceptAODTrack(track)) continue;
    //------------------
    // NUE & NUA
    //------------------
    float phi = track->Phi();
    float pt = track->Pt();
    float eta = track->Eta();
    int charge = track->Charge();
    int id = track->GetID();
    int nhits = track->GetTPCNcls();
    float dedx = track->GetTPCsignal();

    // DCA Cut
    float dcaxy = -999, dcaz = -999;
    if (!GetDCA(dcaxy, dcaz, track)) continue;
    // if FB = 96 or 768, we don't need special DCA cut

    // but we need to set the dca cut for 768 when we start to choose the paiticle for pair
    if (fabs(dcaxy) > 2.4) continue;
    if (fabs(dcaz) > 3.2) continue;
    if (fFilterBit == 768 && isNarrowDcaCuts768) {
      if (fabs(dcaz) > 2.0) continue;
      if (fabs(dcaxy) > 7.0 * (0.0026 + 0.005 / TMath::Power(pt, 1.01))) continue;
    }

    bool isProton = CheckPIDofParticle(track, 3);
    isProton = isProton && (fabs(dcaz) < 1. && fabs(dcaxy) < (0.0105 + 0.035 / TMath::Power(pt, 1.1)));
    if(!isProton) continue;

    if (charge > 0) {
      h2ProtonPtDcaXY->Fill(pt, dcaxy);
      h2ProtonPtDcaZ->Fill(pt, dcaz);
    } else {
      h2AntiProtonPtDcaXY->Fill(pt, dcaxy);
      h2AntiProtonPtDcaZ->Fill(pt, dcaz);
    }
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskProtonPtDCA::RejectEvtTFFit(float centSPD0) {
  int nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  int nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  int nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();

  int nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const int nTracks = fAOD->GetNumberOfTracks();
  int multTrk = 0;
  for (int it = 0; it < nTracks; it++) {
    AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
    if (!aodTrk) continue;
    if (aodTrk->TestFilterBit(32)) {
      if ((fabs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2)) multTrk++;
    }
  }

  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  if(!aodV0) return false;
  float multV0a = aodV0->GetMTotV0A();
  float multV0c = aodV0->GetMTotV0C();
  float multV0Tot = multV0a + multV0c;
  unsigned short multV0aOn = aodV0->GetTriggerChargeA();
  unsigned short multV0cOn = aodV0->GetTriggerChargeC();
  unsigned short multV0On = multV0aOn + multV0cOn;

  // pile-up cuts
  if (centSPD0 < fCenCutLowPU->Eval(fCent)) return false;
  if (centSPD0 > fCenCutHighPU->Eval(fCent)) return false;
  if (nITSCls > fSPDCutPU->Eval(nITSTrkls)) return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;
  if (multTrk < fMultCutPU->Eval(fCent)) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fAOD->IsIncompleteDAQ()) return false;

  if (isTightPileUp == true) {
    int tpcClsTot = fAOD->GetNumberOfTPCClusters();
    float nclsDif = tpcClsTot - (53182.6 + 113.326 * multV0Tot - 0.000831275 * multV0Tot * multV0Tot);
    if (nclsDif > 200000)  // can be varied to 150000, 200000
      return false;
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskProtonPtDCA::AcceptAODTrack(AliAODTrack* track) {
  //------------------
  // track cut
  //------------------
  float pt = track->Pt();
  if (pt < 0.2 || pt > 5.0) return false;
  float eta = track->Eta();
  if (fabs(eta) > 0.8) return false;
  int nhits = track->GetTPCNcls();
  if (nhits < fNclsCut) return false;
  float dedx = track->GetTPCsignal();
  if (dedx < 10.0) return false;
  float chi2 = track->Chi2perNDF();
  if (chi2 < fChi2Min) return false;
  if (chi2 > fChi2Max) return false;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskProtonPtDCA::CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck) {
  if (pidToCheck == 0) return kTRUE;  //// Charge Particles do not need PID check

  if (!fPIDResponse) {
    Printf("\n Could Not access PIDResponse Task, please Add the Task...\n return with kFALSE pid\n");
    return kFALSE;
  }

  /// Rihan todo: To set the low pT cuts for nSigmaTPC from AddTaskMacro!
  /// Although someone barely needs to change it given the purity..

  float nSigTPC = 0, nSigTOF = 0, nSigRMS = 0;
  float trkPtPID = ftrack->Pt();

  /// Pion =>
  if (pidToCheck == 1) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(
        ftrack, AliPID::kPion);  // Some warning show here (***TDatabasePDG::AddParicle: particle with PDGcode = 3124
                                 // already defind),I don't understand what happended. --chunzheng
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);

    if (trkPtPID <= 0.5 && fabs(nSigTPC) <= fNSigmaTPC) return true;
    if (trkPtPID > 0.5 && fabs(nSigRMS) <= fNSigmaRMS) return true;
    return false;
  }
  /// Kaon =>
  else if (pidToCheck == 2) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kKaon);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
    nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);

    if (trkPtPID <= 0.45 && fabs(nSigTPC) <= fNSigmaTPC) return true;
    if (trkPtPID > 0.45 && fabs(nSigRMS) <= fNSigmaRMS) return true;
    return false;
  }
  /// proton =>
  else if (pidToCheck == 3) {  ///
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kProton);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kProton);
    nSigRMS = TMath::Sqrt(nSigTPC * nSigTPC + nSigTOF * nSigTOF);

    bool isProton = fabs(nSigRMS) < fNSigmaRMS;
    return isProton;
  } else {
    Printf("\n -Ve number not allowed! Choose among: 0,1,2,3 (Charge Pion, Kaon, Proton)\n return with kFALSE \n");
    return false;
  }
}

//---------------------------------------------------

bool AliAnalysisTaskProtonPtDCA::GetDCA(float& dcaxy, float& dcaz, AliAODTrack* track) {
  if (!track) return false;
  double r[3];
  if (track->GetXYZ(r)) {
    dcaxy = r[0];
    dcaz = r[1];
  } else {
    float dcax = r[0] - fVertex[0];
    float dcay = r[1] - fVertex[1];
    dcaz = r[2] - fVertex[2];
    dcaxy = sqrt(dcax * dcax + dcay * dcay);
  }
  return true;
}
