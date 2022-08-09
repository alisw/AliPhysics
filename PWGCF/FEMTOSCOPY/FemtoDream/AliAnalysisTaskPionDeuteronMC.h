#ifndef ALIANALYSISTASKPIONDEUTERONMC_H
#define ALIANALYSISTASKPIONDEUTERONMC_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TArrayF.h"
#include "TLorentzVector.h"
#include "CustomQueue.h"

//____________________________________________________________________________________________________________________________________
class AliAnalysisTaskPionDeuteronMC : public AliAnalysisTaskSE
{

public:
    AliAnalysisTaskPionDeuteronMC(const char *name = "PioDuteronMCTask");
    virtual ~AliAnalysisTaskPionDeuteronMC();

    // General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    void SetSourceSizeRadius(TH1F *hSourceR0) { hSourceSize = hSourceR0; }
    void SetCentralityEstimator(int est) { fEstimator = est; }

    // Utils
    void SetPrimaryPtBins(int nbins, float min, float max);
    void SetPrimaryPtBins(int nbins, float *bins);
    void SetKstarBins(int nbins, float min, float max);
    void SetKstarBins(int nbins, float *bins);
    void SetP0(float p0) { fP0 = p0; }

    void SetZvtxArray(std::vector<int> &vec);
    void SetMultiplicityArray(std::vector<int> &vec);

    int FindBin(float zvtx, float mult);

    void SetMixingDepth(unsigned int depth);

    AliEventCuts fEventCuts;

private:
    TList *fOutputList; ///< Output list

    int fEstimator;            ///< Choose the centrality estimator from AliEventCut
    float fP0;                 ///< Coalescence momentum p0
    unsigned int fMixingDepth; /// Depth of the mixing buffer

    TArrayF fPrimaryPtBins; ///<  Transverse momentum bins
    TArrayF fKstarBins;     ///<  realtive momentum bins

    TH1F *fNormalisationHist;   //!<! Event selection
    TH1F *hSourceSize;          //!<! distribution of the source radius
    TH1F *hPionSpectrum[2];     //!<! proton trasnverse momentum spectrum {0: positive, 1: negative}
    TH1F *hProtonSpectrum[2];   //!<! proton trasnverse momentum spectrum {0: positive, 1: negative}
    TH1F *hNeutronSpectrum[2];  //!<! neutron trasnverse momentum spectrum {0: positive, 1: negative}
    TH1F *hDeuteronSpectrum[2]; //!<! deuteron trasnverse momentum spectrum {0: positive, 1: negative}
    TH1F *hDeltaP;              //!<! delta-P distribution

    TH2F *hNparticles[2]; //!<! number of pions and deuterons per event {0: positive, 1: negative}

    TH1F *hSameEventKstarLS[2]; //!<! same-event k* distribution (like-sign) {0: positive, 1: negative}
    TH1F *hSameEventKstarUS[2]; //!<! same-event k* distribution (unlike-sign) {0: positive pion, 1: negative pion}

    TH1F *hMixedEventKstarLS[2]; //!<! same-event k* distribution (like-sign) {0: positive, 1: negative}
    TH1F *hMixedEventKstarUS[2]; //!<! same-event k* distribution (unlike-sign) {0: positive pion, 1: negative pion}

    std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferDeuterons;     ///< Mixing buffer for deuterons
    std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferAntideuterons; ///< Mixing buffer for antideuterons

    std::vector<float> fZvtxArray;         ///< Arrays with the z of the primary vertex for binning
    std::vector<float> fMultiplicityArray; ///< Arrays with multiplicity for binning
    int fNZvtxBins;                        ///< Number of z vertex bins
    int fNMultiplicityBins;                ///< Number of multiplicity bins

    float GetKstar(TLorentzVector &p1, TLorentzVector &p2);                                                               ///< return relative momentum in the rest frame of the pair
    void FillMixedEvent(std::vector<TLorentzVector> &vec, CustomQueue<std::vector<TLorentzVector>> &buffer, TH1F *histo); ///< fill mixed event distribution

    void DoSimpleCoalescence(std::vector<TLorentzVector> &v_proton, std::vector<TLorentzVector> &v_neutron, std::vector<TLorentzVector> &v_deuteron, TH1F *histo); ///< create deuterons from protons and neutrons

    ClassDef(AliAnalysisTaskPionDeuteronMC, 3);
};
//____________________________________________________________________________________________________________________________________

#endif
