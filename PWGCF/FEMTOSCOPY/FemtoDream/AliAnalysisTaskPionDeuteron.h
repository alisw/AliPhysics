#ifndef ALIANALYSISTASKPIONDEUTERON_H
#define ALIANALYSISTASKPIONDEUTERON_H

#include "AliEasyFemto.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "TList.h"
#include "AliPIDResponse.h"

class AliAnalysisTaskPionDeuteron : public AliAnalysisTaskSE, public AliEasyFemto
{
public:
  AliAnalysisTaskPionDeuteron(const char *name = "PioDuteronTask");
  virtual ~AliAnalysisTaskPionDeuteron();

  // General Functions
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  void SetMixingDepth(unsigned int depth);

  void SetCentralityEstimator(int est) { fEstimator = est; }

  // Selections
  AliEventCuts fEventCuts;
  CutContainer fPionCuts;
  CutContainer fProtonCuts;
  CutContainer fDeuteronCuts;

private:
  TList *fOutputList; ///< Output list

  AliPIDResponse *fPID; //!<! PID response class
  int fEstimator;       ///< Choose the centrality estimator from AliEventCut

  TH1F *fNormalisationHist;       //!<! Event selection
  TH1F *hPionTrackSelections;     //!<! Track selection for pions
  TH1F *hProtonTrackSelections;   //!<! Track selection for protons
  TH1F *hDeuteronTrackSelections; //!<! Track selection for deuterons
  TH1F *hDeltaP;                  //!<! delta-P distribution

  TH1F *hPionSpectrum[2];         //!<! proton transverse momentum spectrum {0: positive, 1: negative}
  TH1F *hProtonSpectrum[2];       //!<! proton transverse momentum spectrum {0: positive, 1: negative}
  TH1F *hDeuteronSpectrum[2];     //!<! deuteron transverse momentum spectrum {0: positive, 1: negative}
  TH1F *hFakeDeuteronSpectrum[2]; //!<! fake deuteron (obtained via coalescence) transverse momentum spectrum {0: positive, 1: negative}

  TH2F *hNparticles[2];     //!<! number of pions and deuterons per event {0: positive, 1: negative}
  TH2F *hNparticlesFake[2]; //!<! number of pions and fake deuterons per event {0: positive, 1: negative}

  TH1F *hSameEventPionDeuteronKstarLS[2]; //!<! same-event k* distribution for pion-deuteron (like-sign) {0: positive, 1: negative}
  TH1F *hSameEventPionDeuteronKstarUS[2]; //!<! same-event k* distribution for pion-deuteron (unlike-sign) {0: positive pion, 1: negative pion}

  TH1F *hMixedEventPionDeuteronKstarLS[2]; //!<! mixed-event k* distribution for pion-deuteron (like-sign) {0: positive, 1: negative}
  TH1F *hMixedEventPionDeuteronKstarUS[2]; //!<! mixed-event k* distribution for pion-deuteron (unlike-sign) {0: positive pion, 1: negative pion}

  TH1F *hSameEventPionFakeDeuteronKstarLS[2]; //!<! same-event k* distribution for pion-(fake)deuteron (like-sign) {0: positive, 1: negative}
  TH1F *hSameEventPionFakeDeuteronKstarUS[2]; //!<! same-event k* distribution for pion-(fake)deuteron (unlike-sign) {0: positive pion, 1: negative pion}

  TH1F *hMixedEventPionFakeDeuteronKstarLS[2]; //!<! mixed-event k* distribution for pion-(fake)deuteron (like-sign) {0: positive, 1: negative}
  TH1F *hMixedEventPionFakeDeuteronKstarUS[2]; //!<! mixed-event k* distribution for pion-(fake)deuteron (unlike-sign) {0: positive pion, 1: negative pion}

  TH1F *hSameEventPionProtonKstarLS[2]; //!<! same-event k* distribution for pion-proton (like-sign) {0: positive, 1: negative}
  TH1F *hSameEventPionProtonKstarUS[2]; //!<! same-event k* distribution for pion-proton (unlike-sign) {0: positive pion, 1: negative pion}

  TH1F *hMixedEventPionProtonKstarLS[2]; //!<! mixed-event k* distribution for pion-proton (like-sign) {0: positive, 1: negative}
  TH1F *hMixedEventPionProtonKstarUS[2]; //!<! mixed-event k* distribution for pion-proton (unlike-sign) {0: positive pion, 1: negative pion}

  std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferProtons;     ///< Mixing buffer for protons
  std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferAntiprotons; ///< Mixing buffer for antiprotons

  std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferDeuterons;     ///< Mixing buffer for deuterons
  std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferAntideuterons; ///< Mixing buffer for antideuterons

  std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferFakeDeuterons;     ///< Mixing buffer for fake deuterons
  std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferFakeAntideuterons; ///< Mixing buffer for fake antideuterons

  bool applyDeuteronPID(AliAODTrack *track, AliPIDResponse *fPID, CutContainer &cuts, float &tofBeta); ///< apply PID selection for Deuterons

  ClassDef(AliAnalysisTaskPionDeuteron, 2);
};

#endif
