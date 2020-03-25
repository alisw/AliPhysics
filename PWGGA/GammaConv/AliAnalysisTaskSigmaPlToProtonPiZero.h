/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskSigmaPlToProtonPiZero_cxx
#define AliAnalysisTaskSigmaPlToProtonPiZero_cxx

class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliAnalysisTaskConvJet.h"
#include "AliAnalysisTaskJetOutlierRemoval.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TRandom3.h"
#include <vector>
#include <map>
#include "THn.h"

class AliAnalysisTaskSigmaPlToProtonPiZero : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskSigmaPlToProtonPiZero();
                                AliAnalysisTaskSigmaPlToProtonPiZero(const char *name);
        virtual                 ~AliAnalysisTaskSigmaPlToProtonPiZero();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        // Function to set correction task setting
    void InitBack();

        // base functions for selecting photon and meson candidates in reconstructed data
    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    void SetIsHeavyIon(Int_t flag){
      fIsHeavyIon = flag;
    }

    // // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    // switches for additional analysis streams or outputs
    void SetLightOutput(Bool_t flag){fDoLightOutput = flag;}

    void SetInOutTimingCluster(Double_t min, Double_t max){
      fDoInOutTimingCluster = kTRUE; fMinTimingCluster = min; fMaxTimingCluster = max;
      return;
    }
      // Setting the cut lists for the conversion photons
    void SetEventCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fEventCutArray = CutArray;
    }
      // Setting the cut lists for the calo photons
    void SetCaloCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fClusterCutArray = CutArray;
    }
    // Setting the cut lists for the meson
    void SetMesonCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fMesonCutArray = CutArray;
    }
    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}
    void SetTrackMatcherRunningMode(Int_t mode){fTrackMatcherRunningMode = mode;}

    Bool_t IsPi0Selected(AliAODConversionMother* pi0cand, TF1* fFitPi0MassDataLowPt, TF1* fFitPi0MassDataHighPt, TF1* fFitWidthData);
    Bool_t IsPi0SelectedMC(AliAODConversionMother* pi0cand, TF1* fFitPi0MassMCLowPt, TF1* fFitPi0MassMCHighPt, TF1* fFitWidthMC);
    Double_t GetPodAlpha(TLorentzVector sigmaVektor, TLorentzVector protonVektor, TLorentzVector rekombinatedPi0);
    Double_t GetQT(TLorentzVector sigmaVektor, TLorentzVector rekombinatedPi0);

    Bool_t IsRealProtonKink(AliESDkink* kink, AliStack* MCstack);
    Bool_t IsRealPhoton(AliAODConversionPhoton PhotonCandidate, AliMCEvent* fMCEvent);

    protected:
        AliVEvent *             fEvent;
        AliESDEvent*            fESD;           //! input event
        AliMCEvent*             fMCEvent;           //! input MC event
        AliStack*               fMCStack;
        AliPIDResponse*         fPIDResponse; //! pid response object
        TList*                  fOutputList;    //! output list
        TList**                 fESDList;
        THnD**                    fHistSigmaPlus;
        THnD**                    fHistSigmaPlusMC;
        TH1F*                   fHistSigmaPlusMCGen;
        TH2F**                   fHistReconstructedMassPi0;
        TH2F**                   fHistReconstructedMassPi0MC;
        TH2F**                   fHistPodolanski;
        TH1F**                   fHistAllTracksPt;
        TH1F**                   fHistProtonPt;
        TH2F**                   fHistThetaPhi;
        TH2F**                   fHistThetaPhiProton;
        TH1F**                   fHistClusterE;
        TH1F**                   fHistClusterEWOCuts;
        TH1F**                   fHistNClusWoCuts;
        TH1F**                   fHistNClusWCuts;
        TH2F**                   fHistDEDx;
        TH2F**                   fHistTPCSignal;
        TF1**                    fFitPi0MassDataLowPt;
        TF1**                   fFitPi0MassDataHighPt;
        TF1**                    fFitPi0MassMCLowPt;
        TF1**                    fFitPi0MassMCHighPt;
        TF1**                    fFitWidthData;
        TF1**                    fFitWidthMC;

        AliV0ReaderV1*        fV0Reader;                                            // basic photon Selection Task
        TString               fV0ReaderName;
        TString               fCorrTaskSetting;
        TList*                fEventCutArray;
        AliConvEventCuts*     fEventCuts;
        TList*                fClusterCutArray;
        AliCaloPhotonCuts*    fCaloPhotonCuts;                                    // List with Cluster Cuts
        TList*                fMesonCutArray;
        AliConversionMesonCuts*   fMesonCuts;                                       // MesonCutObject
        AliConversionPhotonCuts*  fConversionCuts;                                  // ConversionPhotonCutObject                                     // If a jet is near the EMCal in the current event
        TH1F**                fHistNEvents;                                        //! array of histos with event information
        Int_t                 fnCuts;                                               // number of cuts to be analysed in parallel
        Int_t                 fIsHeavyIon;                                          // switch for pp = 0, PbPb = 1, pPb = 2
        Bool_t                fDoLightOutput;                                       // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
        Int_t                 fIsMC;
        Bool_t                fDoInOutTimingCluster;                                // manual timing cut for cluster to combine cluster within timing cut and without
        Double_t              fMinTimingCluster;                                    // corresponding ranges, min
        Double_t              fMaxTimingCluster;                                    // corresponding ranges, max
        Int_t                 fTrackMatcherRunningMode;                             // CaloTrackMatcher running mode

                                              // flag for MC information
                                           // switches ranges of histograms and binnings to pi0 specific analysis
    private:
        AliAnalysisTaskSigmaPlToProtonPiZero(const AliAnalysisTaskSigmaPlToProtonPiZero&); // not implemented
        AliAnalysisTaskSigmaPlToProtonPiZero& operator=(const AliAnalysisTaskSigmaPlToProtonPiZero&); // not implemented

        void FillfHistNEvents(Int_t icut, Float_t in) { if(fHistNEvents[icut]) fHistNEvents[icut]->Fill(in); }

        ClassDef(AliAnalysisTaskSigmaPlToProtonPiZero, 15);
};

#endif
