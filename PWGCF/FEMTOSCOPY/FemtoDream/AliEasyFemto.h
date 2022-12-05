#ifndef ALIEASYFEMTO_H
#define ALIEASYFEMTO_H

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

class CutContainer
{
public:
  int filterBit = 128;
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

  void Print()
  {
    printf("filterBit: %d\n", filterBit);
    printf("minPt: %f\n", minPt);
    printf("maxPt: %f\n", maxPt);
    printf("nTPCcls: %d\n", nTPCcls);
    printf("nCrossedRows: %d\n", nCrossedRows);
    printf("eta: %f\n", eta);
    printf("dcaZ: %f\n", dcaZ);
    printf("dcaXY: %f\n", dcaXY);
    printf("crossedRowsOverFindable: %f\n", crossedRowsOverFindable);
    printf("nSigmaTPC: %f\n", nSigmaTPC);
    printf("nSigmaTOF: %f\n", nSigmaTOF);
    printf("nSigmaComb: %f\n", nSigmaComb);
    printf("ptPIDthreshold: %f\n", ptPIDthreshold);
  }

  ClassDefNV(CutContainer, 1);
};

class AliEasyFemto
{
public:
  AliEasyFemto();
  virtual ~AliEasyFemto(){};

  // Utils
  void SetPrimaryPtBins(int nbins, float min, float max);
  void SetPrimaryPtBins(int nbins, float *bins);
  void SetKstarBins(int nbins, float min, float max);
  void SetKstarBins(int nbins, float *bins);
  void SetP0(float p0) { fP0 = p0; }

  void SetZvtxArray(std::vector<float> &vec);
  void SetMultiplicityArray(std::vector<float> &vec);

  // virtual void SetMixingDepth(unsigned int depth);
  int FindBin(float zvtx, float mult);
  static float GetKstar(TLorentzVector &p1, TLorentzVector &p2);                                                        ///< return relative momentum in the rest frame of the pair
  void FillMixedEvent(std::vector<TLorentzVector> &vec, CustomQueue<std::vector<TLorentzVector>> &buffer, TH1F *histo); ///< fill mixed event distribution

  // Selections
  enum kTrackSelection
  {
    kNoSelection = 0,
    kFilterBit,
    kEta,
    kPtMin,
    kPtMax,
    kTPCnCls,
    kTPCnClsS,
    kTPCnCrossedRows,
    kTPCnCrossedOverFindable,
    kDCAz,
    kDCAxy,
    kPID,
    kNselections
  };

  const char *vSelLabels[kTrackSelection::kNselections] = {
      "No selections",
      "Filter bit",
      "#eta",
      "Min #it{p}_{T}",
      "Max #it{p}_{T}",
      "TPC cls",
      "Shared TPC cls",
      "TPC crossed rows",
      "Crossed rows / findable",
      "DCA_{z}",
      "DCA_{xy}",
      "PID"};

  float hasTOF(AliAODTrack *track, AliPIDResponse *fPID);                                                           ///< returns TOF beta after checks
  bool applyTrackSelection(AliAODTrack *track, CutContainer &cuts, float &dcaXY, TH1F *hTrackSelections = nullptr); ///< apply track selections (no dcaXY, no PID)
  bool applyStandardPID(AliAODTrack *track, AliPIDResponse *fPID, CutContainer &cuts, AliPID::EParticleType kType); ///< apply PID standard selection

  // Coalescence utils
  static float ComputeRadius(float mt);
  static float ComputeSimpleGaussProbability(float q, float radius, float d);
  static float ComputeTwoGaussProbability(float q, float radius, float d1, float d2, float Delta);
  static TLorentzVector BaseCoalescence(TLorentzVector &a, TLorentzVector &b, float &deltaP, float &mt); ///< perform coalescence of two nucleons

  // Coalescence
  void SetSimpleCoalescence(bool makeItSimple) { fSimpleCoalescence = makeItSimple; }
  void SetTwoGauss(bool useTwoGauss) { fTwoGauss = useTwoGauss; }

  void DoCoalescence(std::vector<TLorentzVector> &v_proton, std::vector<TLorentzVector> &v_neutron, std::vector<TLorentzVector> &v_deuteron, TH1F *hDeuteron, TH1F *hDeltaP);                     ///< create deuterons from protons and neutrons
  void DoFakeCoalescence(std::vector<TLorentzVector> &v_proton, std::vector<TLorentzVector> &v_fakedeuteron, std::vector<std::pair<int, int>> &v_fakedeuteronID, TH1F *hDeuteron, TH1F *hDeltaP); ///< create fake deuterons from protons

protected:
  float fP0;                 ///< Coalescence momentum p0
  unsigned int fMixingDepth; /// Depth of the mixing buffer

  bool fSimpleCoalescence; ///< If true use simple coalescence, otherwise Wigner approach
  bool fTwoGauss;          ///< If true use two-Gauss deuteron waver function, otherwise simple Gauss

  TArrayF fPrimaryPtBins; ///<  Transverse momentum bins
  TArrayF fKstarBins;     ///<  realtive momentum bins

  std::vector<float> fZvtxArray;         ///< Arrays with the z of the primary vertex for binning
  std::vector<float> fMultiplicityArray; ///< Arrays with multiplicity for binning
  int fNZvtxBins;                        ///< Number of z vertex bins
  int fNMultiplicityBins;                ///< Number of multiplicity bins

  ClassDef(AliEasyFemto, 1);
};

#endif
