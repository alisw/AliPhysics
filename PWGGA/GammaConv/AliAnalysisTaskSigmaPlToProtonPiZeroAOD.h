/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskSigmaPlToProtonPiZeroAOD_cxx
#define AliAnalysisTaskSigmaPlToProtonPiZeroAOD_cxx

class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliCaloSigmaCuts.h"
#include "AliAnalysisTaskConvJet.h"
#include "AliAnalysisTaskJetOutlierRemoval.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include <vector>
#include <map>
#include "THn.h"

class AliAnalysisTaskSigmaPlToProtonPiZeroAOD : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskSigmaPlToProtonPiZeroAOD();
                                AliAnalysisTaskSigmaPlToProtonPiZeroAOD(const char *name);
        virtual                 ~AliAnalysisTaskSigmaPlToProtonPiZeroAOD();

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
    // Setting the cut lists for the meson
    void SetSigmaCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fSigmaCutArray = CutArray;
    }
    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}
    void SetTrackMatcherRunningMode(Int_t mode){fTrackMatcherRunningMode = mode;}

    Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);
    Double_t GetPodAlpha(TLorentzVector sigmaVektor, TLorentzVector protonVektor, TLorentzVector rekombinatedPi0);
    Double_t GetQT(TLorentzVector sigmaVektor, TLorentzVector rekombinatedPi0);

    Int_t IsRealProton(AliAODTrack* track, TClonesArray *fAODMCTrackArray, Int_t iCut, Double_t fWeightJetJetMC, Int_t fill);
    Int_t IsRealPhoton(AliAODConversionPhoton *PhotonCandidate, TClonesArray *fAODMCTrackArray, Int_t iCut, Double_t fWeightJetJetMC);
    Int_t IsProtonFromXi0(AliAODTrack* track, TClonesArray *fAODMCTrackArray, Int_t iCut, Double_t fWeightJetJetMC, Int_t fill);
    Int_t IsPhotonFromXi0(AliAODConversionPhoton *PhotonCandidate, TClonesArray *fAODMCTrackArray, Int_t iCut, Double_t fWeightJetJetMC);
    void CalculateBackgroundSwappWGammaGamma(vector < AliVCluster* > photon, vector < AliAODTrack* > proton, Double_t vpos[3], Int_t iCut, Double_t fWeightJetJetMC);
    void CalculateBackgroundSwappWProtonPion(vector < TLorentzVector > pions, vector < AliAODTrack* > proton, Int_t iCut, Double_t fWeightJetJetMC);

    protected:
        AliVEvent *             fEvent;         //!
        AliAODEvent*            fAOD;           //! input event
        AliMCEvent*             fMCEvent;           //! input MC event
        AliPIDResponse*         fPIDResponse; //! pid response object
        TClonesArray*			fAODMCTrackArray;   //!
        TList*                  fOutputList;    //! output list
        TList**                 fAODList; //!
        TH2F**                  fHistSigmaPlus; //!
        TH2F**                  fHistSigmaPlusMC; //! 
        TH2F**                  fHistSigmaPlusMCTrueProtonGamma; //!
        TH2F**                  fHistSigmaPlusMCTrueProton; //!
        TH2F**                  fHistSigmaPlusMCTruePion; //!
        TH2F**                  fHistDoubleCountTrueSigmaInvMassPt; //!
        TH1F**                  fHistGenSigmaPt; //!
        TH1F**                  fHistGenSigmaPerEvent; //!
        TH1F**                  fHistGenProtonPt;//!
        TH1F**                  fHistGenPiZeroPt; //!
        TH2F**                  fHistSigmaPtEta; //!
        TH2F**                  fHistProtonPtEta; //!
        TH2F**                  fHistPi0PtEta; //!
        TH1F**                  fHistGenAngleProtonPiZero; //!
        TH2F**                  fHistReconstructedMassPi0; //!
        TH2F**                  fHistReconstructedMassPi0MC; //!
        TH2F**                  fHistReconstructedMassPi0wCut; //!
        TH2F**                  fHistReconstructedMassPi0MCwCut; //!
        TH2F**                  fHistPodolanski; //!
        TH2F**                  fHistPodolanskiWCut; //!
        TH2F**                  fHistPodolanskiWCutTrue; //!
        TH2F**                  fHistPodolanskiGenTrue; //!
        TH1F**                  fHistAllTracksPt; //!
        TH1F**                  fHistProtonPt; //!
        TH1F**                  fHistTrueProtonPt; //!
        TH2F**                  fHistThetaPhiTrueSigmaPl; //!
        TH2F**                  fHistThetaPhi; //!
        TH2F**                  fHistThetaPhiProton; //!
        TH1F**                  fHistClusterE; //!
        TH1F**                  fHistClusterEWOCuts; //!
        TH1F**                  fHistNClusWoCuts; //!
        TH1F**                  fHistNClusWCuts; //!
        TH1F**                  fHistNProtonsPerEvent; //!
        TH2F**                  fHistDecayangle; //!
        TH2F**                  fHistDecayangleTrue; //!
        TH2F**                  fHistDecayanglewCut; //!
        TH2F**                  fHistDecayangleTruewCut; //!
        TH1F**                  fHistTrackDCAXY; //!
        TH1F**                  fHistTrackDCAZ; //!
        TH1F**                  fHistTrackDCAXYTrue; //!
        TH1F**                  fHistTrackDCAZTrue; //!
        TH1F**                  fHistTrackDCAXYwCuts; //!
        TH1F**                  fHistTrackDCAZwCuts; //!
        TH1F**                  fHistTrackDCAXYTruewCuts; //!
        TH1F**                  fHistTrackDCAZTruewCuts; //!
        TH2F**                  fHistDEDx; //!
        TH2F**                  fHistTOFBeta; //!
        TH2F**                  fHistTPCSignal; //!
        TH1F**                  fHistTPCCluster; //!
        TH1F**                  fHistTPCClusterTrue; //!
        TH1F**                  fHistTPCchi2; //!
        TH1F**                  fHistTPCchi2True; //!
        TH1F**                  fHistITSCluster; //!
        TH1F**                  fHistITSClusterTrue; //!
        TH1F**                  fHistITSchi2; //!
        TH1F**                  fHistITSchi2True; //!
        TH1F**                  fHistTPCClusterwCut; //!
        TH1F**                  fHistTPCClusterTruewCut; //!
        TH1F**                  fHistTPCchi2wCut; //!
        TH1F**                  fHistTPCchi2TruewCut; //!
        TH1F**                  fHistITSClusterwCut; //!
        TH1F**                  fHistITSClusterTruewCut; //!
        TH1F**                  fHistITSchi2wCut; //!
        TH1F**                  fHistITSchi2TruewCut; //!
        TH2F**                  fHistRotationWGammaGamma; //!
        TH2F**                  fHistRotationWProtonPion; //!
        TH2F**                  fHistSigmaMassPtWoPodCut; //!
        TH2F**                  fHistSigmaMassPtWoPodCutMC; //!
        TH1F**                  fHistNLoopsProton; //!
        TH1F**                  fHistNLoopsGamma; //!
        TH2F**                  fHistXi0MC; //! 
        TGenPhaseSpace          fGenPhaseSpace;                                       // For generation of decays into two gammas
        AliV0ReaderV1*          fV0Reader;                                            // basic photon Selection Task
        TString                 fV0ReaderName;
        TString                 fCorrTaskSetting;
        TList*                  fEventCutArray;
        TList*                  fClusterCutArray;
        TList*                  fMesonCutArray;
        TList*                  fSigmaCutArray;
        TH1F**                  fHistNEvents;                                        //! array of histos with event information
        TH1F**                  fHistNEventsWOWeight;                                        //! array of histos with event information
        Int_t                   fnCuts;                                               // number of cuts to be analysed in parallel
        Int_t                   fIsHeavyIon;                                          // switch for pp = 0, PbPb = 1, pPb = 2
        Bool_t                  fDoLightOutput;                                       // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
        Int_t                   fIsMC;
        Bool_t                  fDoInOutTimingCluster;                                // manual timing cut for cluster to combine cluster within timing cut and without
        Double_t                fMinTimingCluster;                                    // corresponding ranges, min
        Double_t                fMaxTimingCluster;                                    // corresponding ranges, max
        Int_t                   fTrackMatcherRunningMode;  
        Double_t                fQTCutUpper;                           // CaloTrackMatcher running mode
        Double_t                fQTCutLower;                           // CaloTrackMatcher running mode
        Double_t                fWeightJetJetMC;                                    //! weight for Jet-Jet MC
        AliAnalysisTaskJetOutlierRemoval*   fOutlierJetReader;      //! JetReader
                                              // flag for MC information
                                           // switches ranges of histograms and binnings to pi0 specific analysis
    private:
        AliAnalysisTaskSigmaPlToProtonPiZeroAOD(const AliAnalysisTaskSigmaPlToProtonPiZeroAOD&); // not implemented
        AliAnalysisTaskSigmaPlToProtonPiZeroAOD& operator=(const AliAnalysisTaskSigmaPlToProtonPiZeroAOD&); // not implemented

        void FillfHistNEvents(Int_t icut, Float_t in, Double_t fWeightJetJetMC) { if(fHistNEvents[icut]) fHistNEvents[icut]->Fill(in, fWeightJetJetMC); }

        ClassDef(AliAnalysisTaskSigmaPlToProtonPiZeroAOD, 27);
};

#endif
