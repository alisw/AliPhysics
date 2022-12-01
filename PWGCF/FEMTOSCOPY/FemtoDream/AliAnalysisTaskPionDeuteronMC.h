#ifndef ALIANALYSISTASKPIONDEUTERONMC_H
#define ALIANALYSISTASKPIONDEUTERONMC_H


#include "AliEasyFemto.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TArrayF.h"
#include "TLorentzVector.h"
#include "CustomQueue.h"
class AliAnalysisTaskPionDeuteronMC : public AliAnalysisTaskSE, public AliEasyFemto
{

public:
    AliAnalysisTaskPionDeuteronMC(const char *name = "PioDuteronMCTask");
    virtual ~AliAnalysisTaskPionDeuteronMC();

    // General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    void SetMixingDepth(unsigned int depth);

    AliEventCuts fEventCuts;

private:
    TList *fOutputList; ///< Output list

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

    ClassDef(AliAnalysisTaskPionDeuteronMC, 5);
};
//____________________________________________________________________________________________________________________________________

#endif
