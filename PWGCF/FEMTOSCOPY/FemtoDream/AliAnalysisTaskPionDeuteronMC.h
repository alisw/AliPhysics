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

    AliEventCuts fEventCuts;

private:
    TList *fOutputList; ///< Output list

    float fP0;                 ///< Coalescence momentum p0
    unsigned int fMixingDepth; /// Depth of the mixing buffer

    bool fSimpleCoalescence; ///< If true use simple coalescence, otherwise Wigner approach
    bool fTwoGauss;          ///< If true use two-Gauss deuteron waver function, otherwise simple Gauss

    TArrayF fPrimaryPtBins; ///<  Transverse momentum bins
    TArrayF fKstarBins;     ///<  realtive momentum bins

    TH1F *fNormalisationHist;   //!<! Event selection

    TH1F *hPionSpectrum[2];     //!<! proton transverse momentum spectrum {0: positive, 1: negative}
    TH1F *hProtonSpectrum[2];   //!<! proton transverse momentum spectrum {0: positive, 1: negative}
    TH1F *hNeutronSpectrum[2];  //!<! neutron transverse momentum spectrum {0: positive, 1: negative}
    TH1F *hDeuteronSpectrum[2]; //!<! deuteron transverse momentum spectrum {0: positive, 1: negative}
    TH1F *hDeltaP;              //!<! delta-P distribution

    TH2F *hNparticles[2]; //!<! number of pions and deuterons per event {0: positive, 1: negative}

    TH1F *hSameEventPionDeuteronKstarLS[2]; //!<! same-event k* distribution for pion-deuteron (like-sign) {0: positive, 1: negative}
    TH1F *hSameEventPionDeuteronKstarUS[2]; //!<! same-event k* distribution for pion-deuteron (unlike-sign) {0: positive pion, 1: negative pion}

    TH1F *hMixedEventPionDeuteronKstarLS[2]; //!<! mixed-event k* distribution for pion-deuteron (like-sign) {0: positive, 1: negative}
    TH1F *hMixedEventPionDeuteronKstarUS[2]; //!<! mixed-event k* distribution for pion-deuteron (unlike-sign) {0: positive pion, 1: negative pion}

    TH1F *hSameEventPionProtonKstarLS[2]; //!<! same-event k* distribution for pion-proton (like-sign) {0: positive, 1: negative}
    TH1F *hSameEventPionProtonKstarUS[2]; //!<! same-event k* distribution for pion-proton (unlike-sign) {0: positive pion, 1: negative pion}

    TH1F *hMixedEventPionProtonKstarLS[2]; //!<! mixed-event k* distribution for pion-proton (like-sign) {0: positive, 1: negative}
    TH1F *hMixedEventPionProtonKstarUS[2]; //!<! mixed-event k* distribution for pion-proton (unlike-sign) {0: positive pion, 1: negative pion}

    std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferProtons;     ///< Mixing buffer for protons
    std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferAntiprotons; ///< Mixing buffer for antiprotons

    std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferDeuterons;     ///< Mixing buffer for deuterons
    std::vector<CustomQueue<std::vector<TLorentzVector>>> fMixingBufferAntideuterons; ///< Mixing buffer for antideuterons

    std::vector<float> fZvtxArray;         ///< Arrays with the z of the primary vertex for binning
    std::vector<float> fMultiplicityArray; ///< Arrays with multiplicity for binning
    int fNZvtxBins;                        ///< Number of z vertex bins
    int fNMultiplicityBins;                ///< Number of multiplicity bins

    float GetKstar(TLorentzVector &p1, TLorentzVector &p2);                                                               ///< return relative momentum in the rest frame of the pair
    void FillMixedEvent(std::vector<TLorentzVector> &vec, CustomQueue<std::vector<TLorentzVector>> &buffer, TH1F *histo); ///< fill mixed event distribution

    void DoCoalescence(std::vector<TLorentzVector> &v_proton, std::vector<TLorentzVector> &v_neutron, std::vector<TLorentzVector> &v_deuteron, TH1F *histo); ///< create deuterons from protons and neutrons

    ClassDef(AliAnalysisTaskPionDeuteronMC, 4);
};
//____________________________________________________________________________________________________________________________________

#endif
