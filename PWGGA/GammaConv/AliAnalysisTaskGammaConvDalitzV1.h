#ifndef ALIANALYSISTASKGAMMACONVDALITZV1_H
#define ALIANALYSISTASKGAMMACONVDALITZV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// Analysis task for pi0->e+e-gamma (Dalitz decay)

#include "AliAnalysisTaskSE.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliDalitzElectronSelector.h"
#include "AliConvEventCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliGammaConversionAODBGHandler.h"
#include "TProfile2D.h"
#include "TH3F.h"
#include <vector>
#include <memory>
using  std::unique_ptr;
#include "AliDalitzAODESD.h"
#include "AliDalitzData.h"
#include "AliDalitzAODESDMC.h"
#include "AliDalitzEventMC.h"

class AliVEvent;
class AliVTrack;
class AliESDInputHandler;
class AliMCEventHandler;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpidCuts;
class AliV0Reader;
class AliGammaConversionHistograms;
class AliTriggerAnalysis;
class AliAODEvent;
class AliAODtrack;
class AliAODtrackCuts;
class AliAODpidCuts;


class AliAnalysisTaskGammaConvDalitzV1: public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskGammaConvDalitzV1();
    AliAnalysisTaskGammaConvDalitzV1( const char* name );
    virtual ~AliAnalysisTaskGammaConvDalitzV1();

    virtual void UserExec(Option_t *);
    virtual void UserCreateOutputObjects();
    virtual Bool_t Notify();
    virtual void Terminate(const Option_t *);

    void SetLogBinningXTH2(TH2* histoRebin);

    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    void SetMoveParticleAccordingToVertex(Bool_t flag){fMoveParticleAccordingToVertex = flag;}

    void SetIsHeavyIon(Int_t flag){
      fIsHeavyIon = flag;
    }

    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void SetEventCutList(Int_t nCuts, TList *CutArray){
      fnCuts= nCuts;
      fCutEventArray = CutArray;
    }
    void SetConversionCutList(Int_t nCuts, TList *CutArray){
      fnCuts= nCuts;
      fCutGammaArray = CutArray;
    }
    void SetElectronCutList(TList *CutArray){
      fCutElectronArray = CutArray;
    }
    void SetMesonCutList(TList *CutArray){
      fCutMesonArray = CutArray;
    }
    void SetDoChicAnalysis(Bool_t flag){ fDoChicAnalysis = flag; }
    void SetDoMesonQA(Int_t flag){ fDoMesonQA = flag; }
    void SetDoTHnSparse(Bool_t flag){fDoTHnSparse = flag;}
    void SetDoHistoDalitzMassLog(Bool_t flag){fDoHistoDalitzMassLog = flag;}
    void SetProductionVertextoVGamma(Bool_t flag) { fSetProductionVertextoVGamma = flag; }
    void SetDoMaterialBudgetWeightingOfGammasForTrueMesons(Bool_t flag) {fDoMaterialBudgetWeightingOfGammasForTrueMesons = flag;}
    void SetDoLightVersion(Bool_t flag) { fDoLightVersion= flag;}

  private:
    void InitBack();
    void ProcessPhotonCandidates();
    void ProcessTruePhotonCandidates(AliAODConversionPhoton*);
    void ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate);
    void ProcessTrueChicCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate, AliAODConversionPhoton *TruejpsiCandidate);
    void RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP);
    void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    void ProcessElectronCandidates();
    void ProcessVirtualGammasCandidates();
    void ProcessMCParticles();
    void CountESDTracks();
    void CalculatePi0DalitzCandidates();
    void CalculateBackground();
    void UpdateEventByEventData();
    void FillElectronQAHistos(AliAODConversionPhoton *Vgamma) const;
    Double_t GetPsiPair(AliDalitzAODESD *trackPos, AliDalitzAODESD *trackNeg ) const;
    Double_t GetPsiPairMC( AliDalitzAODESDMC *fMCPosParticle, AliDalitzAODESDMC *fMCNegParticle) const;
    Double_t GetdeltaPhi(AliDalitzAODESD *trackelectronVgamma, AliDalitzAODESD *trackpositronVgamma ) const;

    Bool_t IsDalitz(AliDalitzAODESDMC *fMCMother) const;
    Bool_t IsPi0DalitzDaughter( Int_t label ) const;
    Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);

    AliV0ReaderV1                     *fV0Reader;
    TString                           fV0ReaderName;
    AliDalitzElectronSelector         *fElecSelector;
    AliGammaConversionAODBGHandler    **fBGHandler;
    AliVEvent                         *fDataEvent;
    AliDalitzData                     *fAODESDEvent;
    AliESDEvent                       *fESDEvent;
    AliMCEvent                        *fMCEvent;
    AliAODEvent                       *fAODEvent;
    AliDalitzEventMC                  *fAODESDEventMC;
    TList                             **fCutFolder;
    TList                             **fESDList;
    TList                             **fBackList;
    TList                             **fMotherList;
    TList                             **fTrueList;
    TList                             **fMCList;
    TList                             **fQAFolder;
    TList                             *fOutputContainer;
    TClonesArray                      *fReaderGammas;
    vector<Int_t>                     fSelectorElectronIndex;
    vector<Int_t>                     fSelectorPositronIndex;
    TList                             *fGoodGammas;
    TList                             *fGoodVirtualGammas;
    TList                             *fGoodElectrons;
    TList                             *fGoodPositrons;
    TList                             *fCutEventArray;
    TList                             *fCutGammaArray;
    TList                             *fCutElectronArray;
    TList                             *fCutMesonArray;
    //TList                           **fGammasPool;
    AliConvEventCuts                  *fEventCuts;
    AliConversionPhotonCuts           *fConversionCuts;
    TH1F                              **hESDConvGammaPt;
    TH1F                              **hESDConvGammaEta;
    //TH2F                            **hESDConvGammaZR;
    THnSparseF                        **sESDConvGammaZR;
    THnSparseF                        **sESDConvGammaXY;
    TH1F                              **hESDDalitzElectronPt;
    TH1F                              **hESDDalitzPositronPt;
    TH1F                              **hESDDalitzElectronPhi;
    TH1F                              **hESDDalitzPositronPhi;
    TH1F                              **hESDDalitzElectronAfterPt;
    TH1F                              **hESDDalitzPositronAfterPt;
    TH1F                              **hESDDalitzElectronAfterEta;
    TH1F                              **hESDDalitzElectronAfterEtaPCut;
    TH1F                              **hESDDalitzPositronAfterEta;
    TH1F                              **hESDDalitzPositronAfterEtaPCut;
    TH1F                              **hESDDalitzElectronAfterPhi;
    TH1F                              **hESDDalitzPositronAfterPhi;
    TH1F                              **hESDDalitzElectronAfterNClsITS;
    TH1F                              **hESDDalitzElectronAfterNClsITSPCut;
    TH1F                              **hESDDalitzPositronAfterNClsITS;
    TH1F                              **hESDDalitzPositronAfterNClsITSPCut;
    TH2F                              **hESDDalitzElectronAfterNFindClsTPC;
    TH1F                              **hESDDalitzElectronAfterNFindClsTPCPCut;
    TH2F                              **hESDDalitzPositronAfterNFindClsTPC;
    TH1F                              **hESDDalitzPositronAfterNFindClsTPCPCut;
    TH2F                              **hESDDalitzElectronAfterNClsTPC;
    TH1F                              **hESDDalitzElectronAfterNClsTPCPCut;
    TH2F                              **hESDDalitzPositronAfterNClsTPC;
    TH1F                              **hESDDalitzPositronAfterNClsTPCPCut;
    TH2F                              **hESDDalitzElectronAfterNCrossedRowsTPC;
    TH1F                              **hESDDalitzElectronAfterNCrossedRowsTPCPCut;
    TH2F                              **hESDDalitzPositronAfterNCrossedRowsTPC;
    TH1F                              **hESDDalitzPositronAfterNCrossedRowsTPCPCut;
    TH2F                              **hESDDalitzPosEleAfterDCAxy;
    TH2F                              **hESDDalitzPosEleAfterDCAz;
    TH2F                              **hESDDalitzElectronAfterTPCdEdxVsP;
    TH2F                              **hESDDalitzPositronAfterTPCdEdxVsP;
    TH2F                              **hESDDalitzElectronAfterTPCdEdxSignalVsP;
    TH2F                              **hESDDalitzPositronAfterTPCdEdxSignalVsP;
    TH2F                              **hESDDalitzElectronAfterTPCdEdxVsEta;
    TH2F                              **hESDDalitzPositronAfterTPCdEdxVsEta;
    TH2F                              **hESDDalitzElectronAfterTPCdEdxVsPhi;
    TH2F                              **hESDDalitzPositronAfterTPCdEdxVsPhi;
    TH2F                              **hESDEposEnegPsiPairDPhi;
    TH2F                              **hESDEposEnegPsiPairEta;
    TH2F                              **hESDEposEnegDPhiEta;
    TH2F                              **hESDEposEnegInvMassPt;
    TH2F                              **hESDTruePi0EposEnegInvMassPi0Pt;
    TH2F                              **hESDEposEnegLikeSignBackInvMassPt;
    TH2F                              **hESDMotherInvMassPt;
    TH2F                              **hESDPi0MotherInvMassPt;
    TH2F                              **hESDPi0MotherDiffInvMassPt;
    TH2F                              **hESDPi0MotherDiffLimInvMassPt;
    TH2F                              **hESDEposEnegInvMassPi0MotherPt;
    TH1F                              **hESDMotherInvMassOpeningAngleGammaElectron;
    THnSparseF                        **sESDMotherInvMassPtZM;
    TH2F                              **hESDMotherBackInvMassPt;
    THnSparseF                        **sESDMotherBackInvMassPtZM;
    TH2F                              **hESDMotherPi0PtY;
    TH2F                              **hESDMotherPi0PtAlpha;
    TH2F                              **hESDMotherPi0PtOpenAngle;
    THnSparseF                        **sESDMotherDalitzPlot;
    TH1F                              **hMCAllGammaPt;
    TH1F                              **hMCAllGammaPi0Pt;
    TH1F                              **hMCConvGammaPt;
    TH2F                              **hMCConvGammaPtR;
    TH1F                              **hMCConvGammaRSPt;
    TH1F                              **hMCConvGammaPi0Pt;
    TH1F                              **hMCAllPositronsPt;
    TH1F                              **hMCDecayPositronPi0Pt;
    TH1F                              **hMCAllElectronsPt;
    TH1F                              **hMCDecayElectronPi0Pt;
    TH1F                              **hMCConvGammaEta;
    TH1F                              **hMCConvGammaR;
    TH1F                              **hMCAllPositronsEta;
    TH1F                              **hMCAllElectronsEta;
    TH1F                              **hMCPi0DalitzGammaPt;
    TH1F                              **hMCPi0DalitzElectronPt;
    TH1F                              **hMCPi0DalitzPositronPt;
    TH1F                              **hMCPi0Pt;
    TH1F                              **hMCPi0WOWeightPt;
    TH1F                              **hMCPi0GGPt;
    TH1F                              **hMCEtaPt;
    TH1F                              **hMCEtaWOWeightPt;
    TH1F                              **hMCEtaGGPt;
    TH1F                              **hMCPi0InAccPt;
    TH1F                              **hMCPi0WOWeightInAccPt;
    TH1F                              **hMCPi0InAccOpeningAngleGammaElectron;
    THnSparseF                        **sMCPi0DalitzPlot;
    TH1F                              **hMCEtaInAccPt;
    TH1F                              **hMCEtaWOWeightInAccPt;
    TH1F                              **hMCChiCPt;
    TH1F                              **hMCChiCInAccPt;
    TH2F                              **hMCPi0EposEnegInvMassPt;
    TH2F                              **hMCEtaEposEnegInvMassPt;
    TH2F                              **hESDEposEnegTruePi0DalitzInvMassPt;
    TH1F                              **hESDEposEnegTruePrimPi0DalitzInvMass;
    TH2F                              **hESDEposEnegTruePi0DalitzPsiPairDPhi;
    TH1F                              **hESDEposEnegTruePi0DalitzPsiPairMC;
    TH2F                              **hESDEposEnegTruePi0DalitzPsiPairEta;
    TH2F                              **hESDEposEnegTruePi0DalitzDPhiEta;
    TH2F                              **hESDEposEnegTrueEtaDalitzInvMassPt;
    TH1F                              **hESDEposEnegTruePrimEtaDalitzInvMass;
    TH2F                              **hESDEposEnegTrueEtaDalitzPsiPairDPhi;
    TH2F                              **hESDEposEnegTruePhotonInvMassPt;
    TH2F                              **hESDEposEnegTrueInvMassPt;
    TH2F                              **hESDEposEnegTrueMotherInvMassPt;
    TH2F                              **hESDEposEnegTruePhotonPsiPairDPhi;
    TH2F                              **hESDEposEnegTruePhotonPsiPairDPhiPtCut;
    TH2F                              **hESDEposEnegTrueJPsiInvMassPt;
    TH2F                              **hESDTrueMotherChiCInvMassPt;
    TH2F                              **hESDTrueMotherChiCDiffInvMassPt;
    TH2F                              **hESDTrueMotherInvMassPt;
    TH2F                              **hESDTrueMotherW0WeightsInvMassPt;
    TH2F                              **hESDTrueMotherDalitzInvMassPt;
    TH2F                              **hESDTrueMotherPi0GGInvMassPt;
    TH2F                              **hESDTrueMotherPi0GGW0WeightsInvMassPt;
    TH2F                              **hESDTruePi0PtY;
    TH2F                              **hESDTruePi0PtAlpha;
    TH2F                              **hESDTruePi0PtOpenAngle;
    THnSparseF                        **sESDTruePi0DalitzPlot;
    TH2F                              **hESDTruePrimaryMotherPi0GGInvMassPt;
    TH2F                              **hESDTrueSecondaryMotherPi0GGInvMassPt;
    TH2F                              **hESDTruePrimaryMotherInvMassMCPt;
    TH2F                              **hESDTruePrimaryMotherInvMassPt;
    TH2F                              **hESDTruePrimaryMotherW0WeightingInvMassPt;
    TH2F                              **hESDTruePrimaryPi0DalitzESDPtMCPt;
    TH2F                              **hESDTrueSecondaryMotherInvMassPt;
    TH2F                              **hESDTrueSecondaryMotherFromK0sInvMassPt;
    TH2F                              **hESDTrueBckGGInvMassPt;
    TH2F                              **hESDTrueBckContInvMassPt;
    TH2F                              **hESDTrueMotherGGInvMassPt;
    TH1F                              **hESDTrueConvGammaPt;
    TH1F                              **hESDTrueConvGammaPtMC;
    TH1F                              **hESDTrueConvGammaR;
    TH1F                              **hESDTrueConvGammaRMC;
    TH1F                              **hESDTruePositronPt;
    TH1F                              **hESDTrueElectronPt;
    TH1F                              **hESDTrueSecConvGammaPt;
    TH1F                              **hESDTrueSecPositronPt;
    TH1F                              **hESDTrueSecElectronPt;
    TH1F                              **hESDTruePi0DalitzConvGammaPt;
    TH1F                              **hESDTruePi0DalitzConvGammaR;
    TH1F                              **hESDTruePi0DalitzPositronPt;
    TH1F                              **hESDTruePi0DalitzPositronPtMB;
    TH1F                              **hESDTruePi0DalitzElectronPt;
    TH1F                              **hESDTruePi0DalitzElectronPtMB;
    TH1F                              **hESDTruePi0DalitzSecConvGammaPt;
    TH1F                              **hESDTruePi0DalitzSecPositronPt;
    TH1F                              **hESDTruePi0DalitzSecElectronPt;
    TH1F                              **hNEvents;
    TH1F                              **hNEventsWOWeight;
    TH1I                              **hNGoodESDTracks;
    TH2F                              **hNGoodESDTracksVsNGoodGammas;
    TH2F                              **hNGoodESDTracksVsNGoodVGammas;
    TH2F                              **fHistoSPDClusterTrackletBackground;        //! array of histos with SPD tracklets vs SPD clusters for background rejection
    TH1I                              **hNV0Tracks;
    //TH3F                              **hESDEposEnegPsiPairpTleptonsDPhi;
    TH3F                              **hESDEposEnegPsiPairpTPionDPhi;
    TH3F                              **hESDEposEnegPsiPairpTEtaDPhi;
    TH3F                              **hESDEposEnegPsiPairpTPhotonDPhi;
    TProfile                          **hEtaShift;
    TH2F                              **fHistoDoubleCountTruePi0InvMassPt;      //! array of histos with double counted pi0s, invMass, pT
    TH2F                              **fHistoDoubleCountTrueEtaInvMassPt;      //! array of histos with double counted etas, invMass, pT
    TH2F                              **fHistoDoubleCountTrueConvGammaRPt;        //! array of histos with double counted photons, R, pT
    TProfile**                        fProfileJetJetXSection;                     //! array of profiles with xsection for jetjet
    TH1F**                            fhJetJetNTrials;                            //! array of histos with ntrials for jetjet
    vector<Int_t>                     fVectorDoubleCountTruePi0s;            //! vector containing labels of validated pi0
    vector<Int_t>                     fVectorDoubleCountTrueEtas;            //! vector containing labels of validated eta
    vector<Int_t>                     fVectorDoubleCountTrueConvGammas;          //! vector containing labels of validated photons

    TRandom3                          fRandom;
    Double_t                          fEventPlaneAngle;
    Double_t                          *fUnsmearedPx;
    Double_t                          *fUnsmearedPy;
    Double_t                          *fUnsmearedPz;
    Double_t                          *fUnsmearedE;
    Double_t                          *fUnsmearedVPx;
    Double_t                          *fUnsmearedVPy;
    Double_t                          *fUnsmearedVPz;
    Double_t                          *fUnsmearedVE;

    Int_t                             fnCuts;
    Int_t                             fiCut;
    Int_t                             fNumberOfESDTracks;
    Int_t                             fNumberOfESDTrackskBoth;
    Int_t                             fNVirtualGammas;
    Bool_t                            fMoveParticleAccordingToVertex;
    Int_t                             fIsHeavyIon;
    Bool_t                            fDoMesonAnalysis;
    Bool_t                            fDoChicAnalysis;
    Int_t                             fDoMesonQA;
    Bool_t                            fSetProductionVertextoVGamma;
    Bool_t                            fIsFromMBHeader;
    Int_t                             fIsMC;
    Bool_t                            fDoTHnSparse;
    Double_t                          fWeightJetJetMC;
    Bool_t                            fDoHistoDalitzMassLog;
    Bool_t                            fDoMaterialBudgetWeightingOfGammasForTrueMesons;
    Bool_t                            fDoLightVersion;

  private:
    AliAnalysisTaskGammaConvDalitzV1( const AliAnalysisTaskGammaConvDalitzV1& ); // Not implemented
    AliAnalysisTaskGammaConvDalitzV1& operator=( const AliAnalysisTaskGammaConvDalitzV1& ); // Not implemented

  ClassDef( AliAnalysisTaskGammaConvDalitzV1, 8 );
};

#endif // ALIANALYSISTASKGAMMACONVDALITZV1_H
