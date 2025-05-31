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
#include "AliLog.h"
#include "AliVEventHandler.h"
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
#include "AliAnalysisTaskCVEUtil.h"
#include <bitset>
#include "AliPIDResponse.h"


ClassImp(AliAnalysisTaskCVEUtil);

static float PTBINMAX = 10.f;
static int PTBINS = 200;
static float PTBINMAX2D = 5;
static int PTBINS2D = 50;
static float DCAXYMAX = 1.0f;
static float DCAZMAX = 2.0f;
static int DCABINS = 500;
static float CENTBINMAX = 70.f;
static int CENTBINS = 7;
static int kRunSplit = 296657;

//---------------------------------------------------
AliAnalysisTaskCVEUtil::AliAnalysisTaskCVEUtil()
    : AliAnalysisTaskSE()
    {
}

//---------------------------------------------------
AliAnalysisTaskCVEUtil::AliAnalysisTaskCVEUtil(const char* name)
    : AliAnalysisTaskSE(name)
    {

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//------------------------------------------------

AliAnalysisTaskCVEUtil::~AliAnalysisTaskCVEUtil() {
  // Destructor
  // histograms are in the output list and deleted when the output
  if (fOutputList) delete fOutputList;
  if (runNumList) delete fOutputList;
  if (fMCHists) delete fMCHists;
  if (fDataHists) delete fDataHists;
}

//---------------------------------------------------

void AliAnalysisTaskCVEUtil::Terminate(Option_t*) {
  // Terminate loop
  Printf("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskCVEUtil::UserCreateOutputObjects() {

  ////////////////////////
  // Pile up Function
  ////////////////////////
  // Rihan 18q/r Pile-up function
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    double parV0[8]    = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
    double parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
    double parFB32[8]  = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};

    fSPDCutPU = std::unique_ptr<TF1>(new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000));

    fV0CutPU = std::unique_ptr<TF1>(new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000));
    fV0CutPU->SetParameters(parV0);

    fCenCutLowPU = std::unique_ptr<TF1>(new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100));
    fCenCutLowPU->SetParameters(parV0CL0);

    fCenCutHighPU = std::unique_ptr<TF1>(new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100));
    fCenCutHighPU->SetParameters(parV0CL0);

    fMultCutPU = std::unique_ptr<TF1>(new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90));
    fMultCutPU->SetParameters(parFB32);
  } else {
      if(!fIsMC) {
          AliFatal("Sorry, Sorry only support LHC18q/r dataset and their anchored MC");
      }
  }

  TString runNumList18q[125] = {
      "296623", "296622", "296621", "296619", "296618", "296616", "296615", "296594", "296553", "296552", "296551",
      "296550", "296548", "296547", "296516", "296512", "296511", "296510", "296509", "296472", "296433", "296424",
      "296423", "296420", "296419", "296415", "296414", "296383", "296381", "296380", "296379", "296378", "296377",
      "296376", "296375", "296312", "296309", "296304", "296303", "296280", "296279", "296273", "296270", "296269",
      "296247", "296246", "296244", "296243", "296242", "296241", "296240", "296198", "296197", "296196", "296195",
      "296194", "296192", "296191", "296143", "296142", "296135", "296134", "296133", "296132", "296123", "296074",
      "296066", "296065", "296063", "296062", "296060", "296016", "295942", "295941", "295937", "295936", "295913",
      "295910", "295909", "295861", "295860", "295859", "295856", "295855", "295854", "295853", "295831", "295829",
      "295826", "295825", "295822", "295819", "295818", "295816", "295791", "295788", "295786", "295763", "295762",
      "295759", "295758", "295755", "295754", "295725", "295723", "295721", "295719", "295718", "295717", "295714",
      "295712", "295676", "295675", "295673", "295668", "295667", "295666", "295615", "295612", "295611", "295610",
      "295589", "295588", "295586", "295585"};
  TString runNumList18r[89] = {"297595", "297590", "297588", "297558", "297544", "297542", "297541", "297540", "297537",
                               "297512", "297483", "297479", "297452", "297451", "297450", "297446", "297442", "297441",
                               "297415", "297414", "297413", "297406", "297405", "297380", "297379", "297372", "297367",
                               "297366", "297363", "297336", "297335", "297333", "297332", "297317", "297311", "297310",
                               "297278", "297222", "297221", "297218", "297196", "297195", "297193", "297133", "297132",
                               "297129", "297128", "297124", "297123", "297119", "297118", "297117", "297085", "297035",
                               "297031", "296966", "296941", "296938", "296935", "296934", "296932", "296931", "296930",
                               "296903", "296900", "296899", "296894", "296852", "296851", "296850", "296848", "296839",
                               "296838", "296836", "296835", "296799", "296794", "296793", "296790", "296787", "296786",
                               "296785", "296784", "296781", "296752", "296694", "296693", "296691", "296690"};
  runNumList = new std::map<int, int>;
  if (fPeriod.EqualTo("LHC18q"))
    for (int i = 0; i < 125; i++) runNumList->insert(std::pair<int, int>(runNumList18q[i].Atoi(), i + 1));
  else if (fPeriod.EqualTo("LHC18r"))
    for (int i = 0; i < 89; i++) runNumList->insert(std::pair<int, int>(runNumList18r[i].Atoi(), i + 1));
  else return;
  fHistRunNumBin = new TH1I("runNumBin", "", (int)runNumList->size(), 1, (int)runNumList->size() + 1);
  std::map<int, int>::iterator iter;
  for (auto runNum : *runNumList) fHistRunNumBin->GetXaxis()->SetBinLabel(runNum.second, Form("%i", runNum.first));

  //------------------
  // Output list
  //------------------
  fOutputList = new TList();
  fOutputList->SetName("fOutputList");
  fOutputList->SetOwner(true);

  CreateAllHistograms();

  fOutputList->Add(fHistRunNumBin);

  PostData(1, fOutputList);
  if (fDebug) Printf("Post fOutputList Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskCVEUtil::UserExec(Option_t*) {

  if (fDebug) Printf("===============================We are in UserExec!!!===================================");
  //----------------------------
  // Handle
  //----------------------------
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
    return;
  }

  if(!fInputHandler) {
    AliError(Form("%s: Could not get InputHandler", GetName()));
    return;
  }

  fPIDResponse = fInputHandler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
    return;
  }

  fMultSel = dynamic_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  if (!fMultSel) {
    AliError(Form("%s: Could not get AliMultSelection", GetName()));
    return;
  }

  fEventCuts = new AliEventCuts();
  if (!fEventCuts) {
    AliError(Form("%s: Could not create AliEventCuts", GetName()));
    return;
  }

  if (fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!fMCEvent) {
      AliError(Form("%s: Could not get MCEvent", GetName()));
      return;
    }
  }

  if (fDebug) Printf("Handles done!");

  //----------------------------
  // Trigger
  //----------------------------
  unsigned int mask = fInputHandler->IsEventSelected();
  bool isTrigselected = false;
  if (fTrigger.EqualTo("kMB")) isTrigselected = mask & AliVEvent::kMB;
  else if (fTrigger.EqualTo("kINT7")) isTrigselected = mask & AliVEvent::kINT7;
  else if (fTrigger.EqualTo("kINT7+kSemiCentral")) isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kSemiCentral);
  else if (fTrigger.EqualTo("kINT7+kCentral+kSemiCentral")) isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kCentral + AliVEvent::kSemiCentral);
  if (isTrigselected == false) return;
  if (fDebug) Printf("trigger done!");

  //----------------------------
  // Vertex
  //----------------------------
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  if (!fVtx) return;
  fVtx->GetXYZ(fVertex.data());
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  if (fabs(fVertex[0]) < 1e-6 || fabs(fVertex[1]) < 1e-6 || fabs(fVertex[2]) < 1e-6) return;
  if (fabs(fVertex[2]) > fVzCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors() < 1) return;
  if (fDebug) Printf("vertex done!");

  //----------------------------
  // Run number
  //----------------------------
  if (!runNumList || runNumList->empty()) {
    AliError("runNumList not initialized!");
    return;
  }
  int runNum = fAOD->GetRunNumber();
  if (fPeriod.EqualTo("LHC18q")) {
    if (runNum > kRunSplit) return;
  }
  else if (fPeriod.EqualTo("LHC18r")) {
    if (runNum < kRunSplit) return;
  }
  auto it = runNumList->find(runNum);
  if (it == runNumList->end()) {
    AliError(Form("Run number %d not in runNumList! No calib files", runNum));
    return;
  }
  int runNumBin = it->second;
  fHistRunNumBin->Fill(runNumBin);

  //----------------------------
  // Centrality
  //----------------------------
  float centV0M = fMultSel->GetMultiplicityPercentile("V0M");
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
  if (fIsMC) {
    if (!fEventCuts->AcceptEvent(fAOD)) return; //TODO: Chun
    // if (!RejectEvtTFFit(centSPD0)) return;
  } else {
    if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      if (!RejectEvtTFFit(centSPD0)) return;
    }
  }
  if (fDebug) Printf("pile-up done!");

  //----------------------------
  // Process MC Particles (only for MC)
  //----------------------------
  if(fIsMC && !ProcessMCParticles()) {
    AliError("MC particles processing failed");
  }

  //----------------------------
  // Loop Tracks / Fill Vectors
  //----------------------------
  if (!LoopTracks()) return;
  if (fDebug) Printf("Loop Tracks done!");

  //----------------------------
  // Loop V0s / Fill Vectors
  //----------------------------
  if (!LoopV0s()) return;
  if (fDebug) Printf("Loop V0s done!");

  //------------------
  // Post output data.
  //------------------
  PostData(1, fOutputList);
  if (fDebug) Printf("analysis done!");
}

//---------------------------------------------------

bool AliAnalysisTaskCVEUtil::LoopTracks() {
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

    float pt = track->Pt();
    int charge = track->Charge();

    // DCA Cut
    float dcaxy = -999, dcaz = -999;
    if (!GetDCA(dcaxy, dcaz, track)) continue;
    if (fabs(dcaxy) > 2.4) continue;
    if (fabs(dcaz) > 3.2) continue;
    if (fFilterBit == 768 && fUseNarrowDCA768) {
      if (fabs(dcaz) > 2.0) continue;
      if (fabs(dcaxy) > 7.0 * (0.0026 + 0.005 / TMath::Power(pt, 1.01))) continue;
    }

    bool isPosHadron = (charge > 0 && CheckPIDofParticle(track, 0));
    bool isNegHadron = (charge < 0 && CheckPIDofParticle(track, 0));
    bool isPosPion   = (charge > 0 && CheckPIDofParticle(track, 1));
    bool isNegPion   = (charge < 0 && CheckPIDofParticle(track, 1));
    bool isProton    = (charge > 0 && CheckPIDofParticle(track, 3));
    bool isAntiProton = (charge < 0 && CheckPIDofParticle(track, 3));
    isProton = isProton && (fabs(dcaz) < 1. && fabs(dcaxy) < (0.0105 + 0.035 / TMath::Power(pt, 1.1)));
    isProton = isProton && pt > 0.4 && pt < 5.0;
    isAntiProton = isAntiProton && (fabs(dcaz) < 1. && fabs(dcaxy) < (0.0105 + 0.035 / TMath::Power(pt, 1.1)));
    isAntiProton = isAntiProton && pt > 0.4 && pt < 5.0;

    std::bitset<static_cast<size_t>(ParticleType::COUNT)> mask_data;
    mask_data.set(0, isPosHadron);
    mask_data.set(1, isNegHadron);
    mask_data.set(2, isPosPion);
    mask_data.set(3, isNegPion);
    mask_data.set(4, isProton);
    mask_data.set(5, isAntiProton);
    mask_data.set(6, false);
    mask_data.set(7, false);

    if (fIsMC) {
      int label = (track->GetLabel());
      if (label < 0) continue;
      auto MCtrack = (AliAODMCParticle *)fMCEvent->GetTrack(label);
      if (!MCtrack) continue;

      int pdg = MCtrack->PdgCode();

      std::bitset<static_cast<size_t>(ParticleType::COUNT)> mask_mc;
      mask_mc.set(0, charge > 0);
      mask_mc.set(1, charge < 0);
      mask_mc.set(2, pdg == 211);
      mask_mc.set(3, pdg == -211);
      mask_mc.set(4, pdg == 2212);
      mask_mc.set(5, pdg == -2212);
      mask_mc.set(6, false);
      mask_mc.set(7, false);

      std::bitset<static_cast<size_t>(ParticleType::COUNT)> mask_rc{mask_data};
      std::bitset<static_cast<size_t>(ParticleType::COUNT)> mask_rc_mc;
      mask_rc_mc = mask_mc & mask_rc;

      bool isFromOrigin = MCtrack->IsPhysicalPrimary();
      bool isFromMaterial = MCtrack->IsSecondaryFromMaterial();
      bool isFromLambda = false;
      bool isFromOthers = false;
      if (MCtrack->IsSecondaryFromWeakDecay()) {
        int motherLabal = MCtrack->GetMother();
        auto mother = (AliAODMCParticle *)fMCEvent->GetTrack(motherLabal);
        int motherPdg = mother ? mother->PdgCode() : 0;
        if (motherPdg == 3122 || motherPdg == -3122) {
          isFromLambda = true;
        } else {
          isFromOthers = true;
        }
      }

      for (auto p : fParticles) {
        std::size_t bit = static_cast<std::size_t>(p);
        if (!mask_rc.test(bit)) continue;
        (*fMCHists)[p].h2_pt_rc->Fill(fCent,pt);
        if (isFromOrigin) {
            (*fMCHists)[p].h3_pt_dcaXY_origin_rc->Fill(fCent,pt,dcaxy);
            (*fMCHists)[p].h3_pt_dcaZ_origin_rc->Fill(fCent,pt,dcaz);
        }
        if (isFromMaterial) {
            (*fMCHists)[p].h3_pt_dcaXY_material_rc->Fill(fCent,pt,dcaxy);
            (*fMCHists)[p].h3_pt_dcaZ_material_rc->Fill(fCent,pt,dcaz);
        }
        if (isFromLambda) {
            (*fMCHists)[p].h3_pt_dcaXY_lambda_rc->Fill(fCent,pt,dcaxy);
            (*fMCHists)[p].h3_pt_dcaZ_lambda_rc->Fill(fCent,pt,dcaz);
        }
        if (isFromOthers) {
            (*fMCHists)[p].h3_pt_dcaXY_other_rc->Fill(fCent,pt,dcaxy);
            (*fMCHists)[p].h3_pt_dcaZ_other_rc->Fill(fCent,pt,dcaz);
        }
      }

      for (auto p : fParticles) {
        std::size_t bit = static_cast<std::size_t>(p);
        if (!mask_rc_mc.test(bit)) continue;
        (*fMCHists)[p].h2_pt_rc_real->Fill(fCent,pt);
        if (isFromOrigin) {
            (*fMCHists)[p].h3_pt_dcaXY_origin_rc_real->Fill(fCent,pt,dcaxy);
            (*fMCHists)[p].h3_pt_dcaZ_origin_rc_real->Fill(fCent,pt,dcaz);
        }
        if (isFromMaterial) {
            (*fMCHists)[p].h3_pt_dcaXY_material_rc_real->Fill(fCent,pt,dcaxy);
            (*fMCHists)[p].h3_pt_dcaZ_material_rc_real->Fill(fCent,pt,dcaz);
        }
        if (isFromLambda) {
            (*fMCHists)[p].h3_pt_dcaXY_lambda_rc_real->Fill(fCent,pt,dcaxy);
            (*fMCHists)[p].h3_pt_dcaZ_lambda_rc_real->Fill(fCent,pt,dcaz);
        }
        if (isFromOthers) {
            (*fMCHists)[p].h3_pt_dcaXY_other_rc_real->Fill(fCent,pt,dcaxy);
            (*fMCHists)[p].h3_pt_dcaZ_other_rc_real->Fill(fCent,pt,dcaz);
        }
      }
    } else {
      for (auto p : fParticles) {
        std::size_t bit = static_cast<std::size_t>(p);
        if (!mask_data.test(bit)) continue;
        (*fDataHists)[p].h2_pt -> Fill(fCent, pt);
        (*fDataHists)[p].h3_pt_dcaXY->Fill(fCent, pt, dcaxy);
        (*fDataHists)[p].h3_pt_dcaZ->Fill(fCent, pt, dcaz);
      }
    }
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEUtil::RejectEvtTFFit(float centSPD0) {
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

  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEUtil::AcceptAODTrack(AliAODTrack* track) {
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

bool AliAnalysisTaskCVEUtil::CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck) {
  if (pidToCheck == 0) return true;  //// Charge Particles do not need PID check

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

bool AliAnalysisTaskCVEUtil::GetDCA(float& dcaxy, float& dcaz, AliAODTrack* track) {
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

bool AliAnalysisTaskCVEUtil::ProcessMCParticles() {
  auto AODMCTrackArray = dynamic_cast<TClonesArray *>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!AODMCTrackArray) return false;

  for (long i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
      AliAODMCParticle *trackMC = static_cast<AliAODMCParticle *>(AODMCTrackArray->At(i));
      if (!trackMC) continue;
      if (!trackMC->IsPhysicalPrimary()) continue; //must be physical primary
      if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMCEvent)) continue;

      float vzMC = trackMC->Zv(); //?
      if (fabs(vzMC) > fVzCut) continue;

      int pdg = trackMC->GetPdgCode();
      float pt = trackMC->Pt();
      float rap = trackMC->Y();
      float eta = trackMC->Eta();
      float charge = trackMC->Charge();

      auto In = [](float x, float a, float b) { return x > a && x < b; };

      bool isPosHadron = charge > 0   && In(pt, 0.2f, 5.0f) && In(eta, -0.8f,0.8f);
      bool isNegHadron = charge < 0   && In(pt, 0.2f, 5.0f) && In(eta, -0.8f,0.8f);
      bool isPosPion   = pdg ==  211  && In(pt, 0.2f, 5.0f) && In(eta, -0.8f,0.8f);
      bool isNegPion   = pdg == -211  && In(pt, 0.2f, 5.0f) && In(eta, -0.8f,0.8f);
      bool isProton    = pdg ==  2212 && In(pt, 0.4f, 5.0f) && In(eta, -0.8f,0.8f);
      bool isAntiProton= pdg == -2212 && In(pt, 0.4f, 5.0f) && In(eta, -0.8f,0.8f);
      bool isLambda    = pdg == 3122  && In(pt, 0.5f, 10.f) && In(rap, -0.5f, 0.5f);
      bool isAntiLambda= pdg == -3122 && In(pt, 0.5f, 10.f) && In(rap, -0.5f, 0.5f);

      std::bitset<static_cast<size_t>(ParticleType::COUNT)> mask_mc;
      mask_mc.set(0, isPosHadron);
      mask_mc.set(1, isNegHadron);
      mask_mc.set(2, isPosPion);
      mask_mc.set(3, isNegPion);
      mask_mc.set(4, isProton);
      mask_mc.set(5, isAntiProton);
      mask_mc.set(6, isLambda);
      mask_mc.set(7, isAntiLambda);

      for (auto p : fParticles) {
          std::size_t bit = static_cast<std::size_t>(p);
          if (!mask_mc.test(bit)) continue;
          (*fMCHists)[p].h2_pt_mc->Fill(fCent, pt);
      }
  }
  return true;
}

//---------------------------------------------------
void AliAnalysisTaskCVEUtil::CreateAllHistograms() {
 if(fIsMC) {
     fMCHists = new std::map<ParticleType, MCParticleHists>();
     for (auto particle : fParticles) {
         (*fMCHists)[particle].h2_pt_mc                = new TH2F(Form("h2_pt_mc_%s", ParticleName(particle)), Form("p_{T} distribution of %s in MC", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS, 0.f, PTBINMAX);
         (*fMCHists)[particle].h2_pt_rc                = new TH2F(Form("h2_pt_rc_%s", ParticleName(particle)), Form("p_{T} distribution of %s in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS, 0.f, PTBINMAX);
         (*fMCHists)[particle].h2_pt_rc_real           = new TH2F(Form("h2_pt_rc_real_%s", ParticleName(particle)), Form("p_{T} distribution of %s in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS, 0.f, PTBINMAX);

         (*fMCHists)[particle].h3_pt_dcaXY_origin_rc   = new TH3F(Form("h3_pt_dcaXY_origin_rc_%s", ParticleName(particle)), Form("p_{T}, dcaXY, origin, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAXYMAX);
         (*fMCHists)[particle].h3_pt_dcaXY_material_rc = new TH3F(Form("h3_pt_dcaXY_material_rc_%s", ParticleName(particle)), Form("p_{T}, dcaXY, material, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAXYMAX);
         (*fMCHists)[particle].h3_pt_dcaXY_lambda_rc   = new TH3F(Form("h3_pt_dcaXY_lambda_rc_%s", ParticleName(particle)), Form("p_{T}, dcaXY, lambda, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAXYMAX);
         (*fMCHists)[particle].h3_pt_dcaXY_other_rc    = new TH3F(Form("h3_pt_dcaXY_other_rc_%s", ParticleName(particle)), Form("p_{T}, dcaXY, other, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAXYMAX);
         (*fMCHists)[particle].h3_pt_dcaZ_origin_rc    = new TH3F(Form("h3_pt_dcaZ_origin_rc_%s", ParticleName(particle)), Form("p_{T}, dcaZ, origin, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);
         (*fMCHists)[particle].h3_pt_dcaZ_material_rc  = new TH3F(Form("h3_pt_dcaZ_material_rc_%s", ParticleName(particle)), Form("p_{T}, dcaZ, material, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);
         (*fMCHists)[particle].h3_pt_dcaZ_lambda_rc    = new TH3F(Form("h3_pt_dcaZ_lambda_rc_%s", ParticleName(particle)), Form("p_{T}, dcaZ, lambda, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);
         (*fMCHists)[particle].h3_pt_dcaZ_other_rc     = new TH3F(Form("h3_pt_dcaZ_other_rc_%s", ParticleName(particle)), Form("p_{T}, dcaZ, other, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);

         (*fMCHists)[particle].h3_pt_dcaXY_origin_rc_real   = new TH3F(Form("h3_pt_dcaXY_origin_rc_real_%s", ParticleName(particle)), Form("p_{T}, dcaXY, origin, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAXYMAX);
         (*fMCHists)[particle].h3_pt_dcaXY_material_rc_real = new TH3F(Form("h3_pt_dcaXY_material_rc_real_%s", ParticleName(particle)), Form("p_{T}, dcaXY, material, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAXYMAX);
         (*fMCHists)[particle].h3_pt_dcaXY_lambda_rc_real   = new TH3F(Form("h3_pt_dcaXY_lambda_rc_real_%s", ParticleName(particle)), Form("p_{T}, dcaXY, lambda, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAXYMAX);
         (*fMCHists)[particle].h3_pt_dcaXY_other_rc_real    = new TH3F(Form("h3_pt_dcaXY_other_rc_real_%s", ParticleName(particle)), Form("p_{T}, dcaXY, other, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAXYMAX);
         (*fMCHists)[particle].h3_pt_dcaZ_origin_rc_real    = new TH3F(Form("h3_pt_dcaZ_origin_rc_real_%s", ParticleName(particle)), Form("p_{T}, dcaZ, origin, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);
         (*fMCHists)[particle].h3_pt_dcaZ_material_rc_real  = new TH3F(Form("h3_pt_dcaZ_material_rc_real_%s", ParticleName(particle)), Form("p_{T}, dcaZ, material, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);
         (*fMCHists)[particle].h3_pt_dcaZ_lambda_rc_real    = new TH3F(Form("h3_pt_dcaZ_lambda_rc_real_%s", ParticleName(particle)), Form("p_{T}, dcaZ, lambda, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);
         (*fMCHists)[particle].h3_pt_dcaZ_other_rc_real     = new TH3F(Form("h3_pt_dcaZ_other_rc_real_%s", ParticleName(particle)), Form("p_{T}, dcaZ, other, %s distribution in real collision", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);
         (*fMCHists)[particle].AddToList(fOutputList);
     }
 } else {
     fDataHists = new std::map<ParticleType, DataParticleHists>();
     for (auto particle : fParticles) {
         (*fDataHists)[particle].h2_pt  = new TH2F(Form("h2_pt_%s", ParticleName(particle)), Form("p_{T}, %s distribution in data", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS, 0.f, PTBINMAX);
         (*fDataHists)[particle].h3_pt_dcaXY = new TH3F(Form("h3_pt_dcaXY_%s", ParticleName(particle)), Form("p_{T}, dcaXY, %s distribution in data", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);
         (*fDataHists)[particle].h3_pt_dcaZ = new TH3F(Form("h3_pt_dcaZ_%s", ParticleName(particle)), Form("p_{T}, dcaZ, %s distribution in data", ParticleName(particle)), CENTBINS, 0.f, CENTBINMAX, PTBINS2D,0.f,PTBINMAX2D, DCABINS, 0.f, DCAZMAX);
         (*fDataHists)[particle].AddToList(fOutputList);
     }
  }
}




//---------------------------------------------------
bool AliAnalysisTaskCVEUtil::IsGoodV0(AliAODv0 *thisV0) {
  // Offline reconstructed V0 only
  if (thisV0->GetOnFlyStatus()) return false;
  // Cosinus of pointing angle < 0.997
  double dCPA = thisV0->CosPointingAngle(fVertex.data());
  if (dCPA < 0.997) return false;
  // DCA of V0 < 1.5 cm
  double dV0Dca = thisV0->DcaV0ToPrimVertex();
  if (TMath::Abs(dV0Dca) > 1.5) return false;
  // V0 path length before decay 3-100 cm
  double dDecayLength = thisV0->DecayLengthV0(fVertex.data());
  if (dDecayLength > 100.) return false;
  if (dDecayLength < 3.) return false;
  // DCA between daughters < 0.5cm
  double dDCA = thisV0->DcaV0Daughters();
  if (dDCA > 0.5) return false;
  // DCA of daughters to PV > 0.05 cm
  float nDcaPV = thisV0->DcaNegToPrimVertex();
  float pDcaPV = thisV0->DcaPosToPrimVertex();
  if (nDcaPV < 0.05 || pDcaPV < 0.05) return false;
  return kTRUE;
}

//---------------------------------------------------
bool AliAnalysisTaskCVEUtil::IsGoodDaughterTrack(
    const AliAODTrack *track) {
  // TPC refit
  if (!track->IsOn(AliVTrack::kTPCrefit)) return false;
  // No kinks
  if (int(track->GetProdVertex()->GetType()) == AliAODVertex::kKink) return false;
  // Maximum value of transverse momentum
  double dPt = track->Pt();
  if (dPt > 20) return false;
  // Maximum value of pseudorapidity
  double dEta = track->Eta();
  if (TMath::Abs(dEta) > 0.8) return false;
  // Minimum number of clusters
  float nCrossedRowsTPC = track->GetTPCClusterInfo(2, 1);
  if (nCrossedRowsTPC < 70) return false;
  // Findable clusters > 0
  int findable = track->GetTPCNclsF();
  if (findable <= 0) return false;
  // [number of crossed rows]>0.8  [number of findable clusters].
  if (nCrossedRowsTPC / findable < 0.8) return false;
  return true;
}

//---------------------------------------------------
int AliAnalysisTaskCVEUtil::GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *nTrack) {
  bool isLambda = kFALSE;
  bool isAntiLambda = kFALSE;
  int code = 0;

  // Λ-->(p+)+(π-)
  float nSigTPCPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)); // TPC p+
  float nSigTPCNegPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)); // TPC π-
  //(Λ-)-->(p-)+(π+)
  float nSigTPCPosPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)); // TPC π+
  float nSigTPCNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)); // TPC p-

  isLambda = (nSigTPCPosProton < 3.) && (nSigTPCNegPion < 3.);
  isAntiLambda = (nSigTPCNegProton < 3.) && (nSigTPCPosPion < 3.);

  if (isLambda) code = 3122;
  if (isAntiLambda) code = -3122;
  if (isLambda && isAntiLambda) code = 0;
  return code;
}


bool AliAnalysisTaskCVEUtil::LoopV0s() {
  int nV0s = fAOD->GetNumberOfV0s();
  for (int iV0 = 0; iV0 < nV0s; iV0++) {
    AliAODv0 *v0 = fAOD->GetV0(iV0);
    // topological cuts
    if (!v0) continue;
    if (!IsGoodV0(v0)) continue;
    AliAODTrack *pTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(0));
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(1));
    if (!pTrack || !nTrack) continue;
    if (pTrack->Charge() * nTrack->Charge() > 0) continue;
    if (!(IsGoodDaughterTrack(nTrack)) || !(IsGoodDaughterTrack(pTrack))) continue;

    float rap = v0->RapLambda();
    if (rap < -0.5 || rap > 0.5) continue;
    float pt = v0->Pt();
    if (pt < 0.5 || pt > 10) continue;

    // reconstruct lambda and anti-lambda
    int code = GetLambdaCode(pTrack, nTrack);
    if (code != 3122 && code != -3122) continue;

    float mass;
    if (code == 3122)       mass = v0->MassLambda();
    else if (code == -3122) mass = v0->MassAntiLambda();
    else continue;
    if (fabs(mass - 1.115683) > 0.02) continue;

    if(fIsMC) {
    // get pdg code from one-to-one MC matching
    if (pTrack->GetLabel() < 0 || nTrack->GetLabel() < 0) continue;
    auto pMC = (AliAODMCParticle *)fMCEvent->GetTrack(pTrack->GetLabel());
    auto nMC = (AliAODMCParticle *)fMCEvent->GetTrack(nTrack->GetLabel());
    if (!pMC || !nMC) continue;
    if (pMC->GetMother() != nMC->GetMother()) continue;
    int v0MotherLabel = pMC->GetMother();
    if (v0MotherLabel < 0) continue;
    auto MCtrack = (AliAODMCParticle *)fMCEvent->GetTrack(v0MotherLabel);
    int pdg = MCtrack->GetPdgCode();

     if (code == 3122) {
         (*fMCHists)[ParticleType::kLambda].h2_pt_rc->Fill(fCent,pt);
         if (pdg == 3122) {
             (*fMCHists)[ParticleType::kLambda].h2_pt_rc_real->Fill(fCent,pt);
         }
     } else if (code == -3122) {
         (*fMCHists)[ParticleType::kAntiLambda].h2_pt_rc->Fill(fCent,pt);
         if (pdg == -3122) {
             (*fMCHists)[ParticleType::kAntiLambda].h2_pt_rc_real->Fill(fCent,pt);
         }
     }
    } else {
       // Fill histograms for data
       if (code == 3122) {
           (*fDataHists)[ParticleType::kLambda].h2_pt->Fill(fCent,pt);
       } else if (code == -3122) {
           (*fDataHists)[ParticleType::kAntiLambda].h2_pt->Fill(fCent,pt);
       }
    }
  }

  return true;
}
