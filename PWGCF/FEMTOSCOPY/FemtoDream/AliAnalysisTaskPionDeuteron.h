#ifndef ALIANALYSISTASKPIONDEUTERON_H
#define ALIANALYSISTASKPIONDEUTERON_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliPID.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TArrayF.h"
#include "AliAODTrack.h"
#include "TLorentzVector.h"
#include "CustomQueue.h"

#include <vector>

#define LIGHT_SPEED 2.99792457999999984e-02 // in the units that TOF likes
#define EPS 1.e-16

struct CutContainer
{
  int filterBit = BIT(7);
  float minPt = 0.5;
  float maxPt = 4.0;
  int nTPCcls = 80;
  int nCrossedRows = 70;
  float eta = 0.8;
  float dcaZ = 0.3;
  float dcaXY = 0.3;
  float crossedRowsOverFindable = 0.83;
  float nSigmaTPC = 3.;
  float nSigmaTOF = 3.;
  float nSigmaComb = 3.;
  float ptPIDthreshold = 0.75;
};

class AliAnalysisTaskPionDeuteron : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskPionDeuteron(const char *name = "PioDuteronTask");
  virtual ~AliAnalysisTaskPionDeuteron();

  // General Functions
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  // Mixing buffer
  void SetCentralityEstimator(int est) { fEstimator = est; }

  // Utils
  void SetPrimaryPtBins(int nbins, float min, float max);
  void SetPrimaryPtBins(int nbins, float *bins);
  void SetKstarBins(int nbins, float min, float max);
  void SetKstarBins(int nbins, float *bins);
  void SetP0(float p0) { fP0 = p0; }

  void SetZvtxArray(std::vector<float> &vec);
  void SetMultiplicityArray(std::vector<float> &vec);

  int FindBin(float zvtx, float mult);

  void SetMixingDepth(unsigned int depth);

  void SetSimpleCoalescence(bool makeItSimple) { fSimpleCoalescence = makeItSimple; }
  void SetTwoGauss(bool useTwoGauss) { fTwoGauss = useTwoGauss; }

  static float ComputeRadius(float mt);

  // Selections
  AliEventCuts fEventCuts;
  CutContainer fPionCuts;
  CutContainer fProtonCuts;
  CutContainer fDeuteronCuts;

private:
  TList *fOutputList; ///< Output list

  AliPIDResponse *fPID; //!<! PID response class

  int fEstimator;            ///< Choose the centrality estimator from AliEventCut
  float fP0;                 ///< Coalescence momentum p0
  unsigned int fMixingDepth; /// Depth of the mixing buffer

  bool fSimpleCoalescence; ///< If true use simple coalescence, otherwise Wigner approach
  bool fTwoGauss;          ///< If true use two-Gauss deuteron waver function, otherwise simple Gauss

  TArrayF fPrimaryPtBins; ///<  Transverse momentum bins
  TArrayF fKstarBins;     ///<  realtive momentum bins

  TH1F *fNormalisationHist;       //!<! Event selection
  TH1F *hPionSpectrum[2];         //!<! proton transverse momentum spectrum {0: positive, 1: negative}
  TH1F *hProtonSpectrum[2];       //!<! proton transverse momentum spectrum {0: positive, 1: negative}
  TH1F *hDeuteronSpectrum[2];     //!<! deuteron transverse momentum spectrum {0: positive, 1: negative}
  TH1F *hFakeDeuteronSpectrum[2]; //!<! fake deuteron (obtained via coalescence) transverse momentum spectrum {0: positive, 1: negative}
  TH1F *hDeltaP;                  //!<! delta-P distribution

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

  std::vector<float> fZvtxArray;         ///< Arrays with the z of the primary vertex for binning
  std::vector<float> fMultiplicityArray; ///< Arrays with multiplicity for binning
  int fNZvtxBins;                        ///< Number of z vertex bins
  int fNMultiplicityBins;                ///< Number of multiplicity bins

  float GetKstar(TLorentzVector &p1, TLorentzVector &p2);                                                               ///< return relative momentum in the rest frame of the pair
  void FillMixedEvent(std::vector<TLorentzVector> &vec, CustomQueue<std::vector<TLorentzVector>> &buffer, TH1F *histo); ///< fill mixed event distribution

  void DoFakeCoalescence(std::vector<TLorentzVector> &v_proton, std::vector<TLorentzVector> &v_fakedeuteron, std::vector<std::pair<int, int>> &v_fakedeuteronID, TH1F *histo); ///< create deuterons from protons

  bool applyTrackSelection(AliAODTrack *track, CutContainer &cuts, float &dcaXY); ///< apply track selections (no dcaXY, no PID)

  bool applyStandardPID(AliAODTrack *track, CutContainer &cuts, AliPID::EParticleType kType); ///< apply PID standard selection

  bool applyDeuteronPID(AliAODTrack *track, CutContainer &cuts, float &tofBeta); ///< apply PID selection for Deuterons

  float hasTOF(AliAODTrack *track); ///< returns TOF beta after checks

  ClassDef(AliAnalysisTaskPionDeuteron, 1);
};

#endif
