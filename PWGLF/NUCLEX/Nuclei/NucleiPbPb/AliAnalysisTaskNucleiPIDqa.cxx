#include "AliAnalysisTaskNucleiPIDqa.h"
#include "AliAnalysisTaskNucleiYield.h"

#include <algorithm>
#include <cmath>
#include <string>
using std::string;

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH2F.h>
#include <TList.h>

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliVTrack.h"

///\cond CLASSIMP
ClassImp(AliAnalysisTaskNucleiPIDqa);
///\endcond

const string                kNames[4]= {"2H","3H","3He","4He"};
const AliPID::EParticleType kSpecies[4] = {AliPID::kDeuteron, AliPID::kTriton, AliPID::kHe3, AliPID::kAlpha};
const char                  kLetter[2] = {'A','M'};
const string                kPIDmethod[5] = {"ITS","TPC","TOF","ITSTPC","TPCTOF"};

void BinLogAxis(TH1 *h) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis* axis = h->GetXaxis();
  int bins = axis->GetNbins();
  double from = axis->GetXmin();
  double to = axis->GetXmax();
  double *newBins = new double[bins + 1];

  newBins[0] = from;
  double factor = pow(to / from, 1. / bins);

  for (int i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i - 1];
  }
  axis->Set(bins, newBins);
  delete[] newBins;
}

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskNucleiPIDqa::AliAnalysisTaskNucleiPIDqa(TString taskname) :  AliAnalysisTaskSE(taskname.Data()),
  fEventCut{false},
  fFilterBit{BIT(4)},
  fNsigmaITS{3.5f},
  fNsigmaTPC{3.5f},
  fNsigmaTOF{3.5f},
  fITSsignalN{3},
  fTPCsignalN{50},
  fList{nullptr},
  fPID{nullptr},
  fITSsignal{nullptr},
  fTPCsignal{nullptr},
  fTOFsignal{nullptr},
  fITSsignalSelected{{{nullptr}}},
  fTPCsignalSelected{{{nullptr}}},
  fTOFsignalSelected{{{nullptr}}},
  fITSnSigmaSelected{{{nullptr}}},
  fTPCnSigmaSelected{{{nullptr}}},
  fTOFnSigmaSelected{{{nullptr}}},
  fUseCustomBethe{false,false,false,false},
  fCustomBethe{},
  fCustomResolution{}
  {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

/// Standard destructor
///
AliAnalysisTaskNucleiPIDqa::~AliAnalysisTaskNucleiPIDqa(){
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fList) delete fList;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskNucleiPIDqa::UserCreateOutputObjects() {

  fList = new TList();
  fList->SetOwner(kTRUE);

  fTPCperformance = new TH2D("fTPCperformance",";#it{p}/ |#it{z}| (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units);Entries",600,0.1,11.,1400,0,1400);
  fTPCperformanceTwoCharges = new TH2D("fTPCperformanceTwoCharges",";#it{p}/ #it{z} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units);Entries",600,-3,3,1400,0,1400);
  BinLogAxis(fTPCperformance);
  fList->Add(fTPCperformance);
  fList->Add(fTPCperformanceTwoCharges);

  for (int iC = 0; iC < 2; ++iC) {
    fITSsignal[iC] = new TH2F(Form("f%cITSsignal",kLetter[iC]),";#it{p} (GeV/#it{c});ITS d#it{E}/d#it{x} (keV / 300#mum);Entries",490,0.2,10.,700,0,1400);
    fTPCsignal[iC] = new TH2F(Form("f%cTPCsignal",kLetter[iC]),";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (a.u.);Entries",490,0.2,10.,1000,0,2000);
    fTOFsignal[iC] = new TH2F(Form("f%cTOFsignal",kLetter[iC]),";#it{p} (GeV/#it{c});#beta;Entries",490,0.2,10.,550,0.,1.1);
    fList->Add(fITSsignal[iC]);
    fList->Add(fTPCsignal[iC]);
    fList->Add(fTOFsignal[iC]);
    for (int iS = 0; iS < 4; ++iS) {
      for (int iPid = 0; iPid < 5; ++iPid) {
        fITSsignalSelected[iC][iS][iPid] = new TH2F(Form("fITSsignalSelected%c%s_%s",kLetter[iC],kNames[iS].data(),kPIDmethod[iPid].data()),
                                                    ";#it{p} (GeV/#it{c});ITS d#it{E}/d#it{x} (keV / 300#mum);Entries",
                                                    490,0.2,10.,700,0,1400);
        fTPCsignalSelected[iC][iS][iPid] = new TH2F(Form("fTPCsignalSelected%c%s_%s",kLetter[iC],kNames[iS].data(),kPIDmethod[iPid].data()),
                                                    ";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (a.u.);Entries",
                                                    490,0.2,10.,1000,0,2000.);
        fTOFsignalSelected[iC][iS][iPid] = new TH2F(Form("fTOFsignalSelected%c%s_%s",kLetter[iC],kNames[iS].data(),kPIDmethod[iPid].data()),
                                                    ";#it{p} (GeV/#it{c});#beta;Entries",
                                                    490,0.2,10.,550,0.,1.1);
        fITSnSigmaSelected[iC][iS][iPid] = new TH2F(Form("fITSnSigmaSelected%c%s_%s",kLetter[iC],kNames[iS].data(),kPIDmethod[iPid].data()),
                                                    ";#it{p} (GeV/#it{c});n#sigma_{ITS};Entries",
                                                    490,0.2,10.,1000,-5.,5.);
        fTPCnSigmaSelected[iC][iS][iPid] = new TH2F(Form("fTPCnSigmaSelected%c%s_%s",kLetter[iC],kNames[iS].data(),kPIDmethod[iPid].data()),
                                                    ";#it{p} (GeV/#it{c});n#sigma_{TPC};Entries",
                                                    490,0.2,10.,1000,-5.,5.);
        fTOFnSigmaSelected[iC][iS][iPid] = new TH2F(Form("fTOFnSigmaSelected%c%s_%s",kLetter[iC],kNames[iS].data(),kPIDmethod[iPid].data()),
                                                    ";#it{p} (GeV/#it{c});n#sigma_{TOF};Entries",
                                                    490,0.2,10.,1000,-5.,5.);
        fList->Add(fITSsignalSelected[iC][iS][iPid]);
        fList->Add(fTOFsignalSelected[iC][iS][iPid]);
        fList->Add(fTPCsignalSelected[iC][iS][iPid]);
        fList->Add(fITSnSigmaSelected[iC][iS][iPid]);
        fList->Add(fTPCnSigmaSelected[iC][iS][iPid]);
        fList->Add(fTOFnSigmaSelected[iC][iS][iPid]);
      }
    }
  }
  fEventCut.AddQAplotsToList(fList);

  PostData(1, fList);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskNucleiPIDqa::UserExec(Option_t *) {
  AliVEvent *ev = InputEvent();
  if (!fEventCut.AcceptEvent(ev)) {
    PostData(1, fList);
    return;
  }

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPID = handl->GetPIDResponse();

  /// Checking how many deuterons in acceptance are reconstructed well
  for (Int_t iT = 0; iT < (Int_t)ev->GetNumberOfTracks(); ++iT) {
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(ev->GetTrack(iT));

    if (track->GetID() <= 0) continue;
    if (!track->TestFilterBit(fFilterBit)) continue;
    if (track->GetTPCsignalN() < fTPCsignalN) continue;
    const float beta = AliAnalysisTaskNucleiYield::HasTOF(track,fPID);
    const int hasTOF = beta > 1.e-24 ? 1 : 0;
    const int iC = (track->Charge() > 0) ? 1 : 0;

    if (track->GetTPCsignalN() > 69) {
      fTPCperformance->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
      fTPCperformanceTwoCharges->Fill(track->GetTPCmomentum() * (2 * iC - 1), track->GetTPCsignal());
    }

    int nSPD = 0u, nSDD = 0u, nSSD = 0u;
    int nITS = AliAnalysisTaskNucleiYield::GetNumberOfITSclustersPerLayer(track, nSPD, nSDD, nSSD);
    const int hasITS = int(nITS - nSPD >= fITSsignalN);

    if (hasITS) fITSsignal[iC]->Fill(track->P(), track->GetITSsignal());
    fTPCsignal[iC]->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
    if (hasTOF) fTOFsignal[iC]->Fill(track->P(),beta);

    for (int iS = 0; iS < 4; ++iS) {
      const float nSigmaITS = fPID->NumberOfSigmasITS(track,kSpecies[iS]);
      float nSigmaTPC = fPID->NumberOfSigmasTPC(track,kSpecies[iS]);
      if (fUseCustomBethe[iS]) {
        const float betaGamma = track->GetTPCmomentum() / AliPID::ParticleMass(kSpecies[iS]);
        const float* pars = fCustomBethe[iS];
        const float expSignal = AliExternalTrackParam::BetheBlochAleph(betaGamma, pars[0], pars[1], pars[2], pars[3], pars[4]);
        nSigmaTPC  = (track->GetTPCsignal() - expSignal) / (fCustomResolution[iS] * expSignal);
      }
      const float nSigmaTOF = fPID->NumberOfSigmasTOF(track,kSpecies[iS]);

      const bool pidITS = std::abs(nSigmaITS) < fNsigmaITS;
      const bool pidTPC = std::abs(nSigmaTPC) < fNsigmaTPC;
      const bool pidTOF = std::abs(nSigmaTOF) < fNsigmaTOF;
      const bool pidConditions[5] = {pidITS, pidTPC, pidTOF, pidITS&&pidTPC, pidTPC&&pidTOF};

      for (int iPid = 0; iPid < 5; ++iPid) {
        if (pidConditions[iPid]) {
          fITSsignalSelected[iC][iS][iPid]->Fill(track->P(),track->GetITSsignal());
          fTPCsignalSelected[iC][iS][iPid]->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
          fTOFsignalSelected[iC][iS][iPid]->Fill(track->P(),beta);
          fITSnSigmaSelected[iC][iS][iPid]->Fill(track->P(),nSigmaITS);
          fTPCnSigmaSelected[iC][iS][iPid]->Fill(track->GetTPCmomentum(),nSigmaTPC);
          fTOFnSigmaSelected[iC][iS][iPid]->Fill(track->P(),nSigmaTOF);
        }
      }
    }
  } // End AOD track loop

  //  Post output data.
  PostData(1,fList);
}

/// Merge the output. Called once at the end of the query.
///
/// \return void
///
void AliAnalysisTaskNucleiPIDqa::Terminate(Option_t *) {
  return;
}

void AliAnalysisTaskNucleiPIDqa::SetCustomBetheBloch(int iSpecies, float res, const float* bethe) {
  fUseCustomBethe[iSpecies] = true;
  fCustomResolution[iSpecies] = res;
  std::copy(bethe, bethe+5, fCustomBethe[iSpecies]);
}
