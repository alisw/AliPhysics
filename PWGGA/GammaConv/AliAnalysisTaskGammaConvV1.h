#ifndef ALIANLYSISTASKGAMMACONVV1_cxx
#define ALIANLYSISTASKGAMMACONVV1_cxx

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliConversionMesonCuts.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include <vector>
#include <map>

class AliAnalysisTaskGammaConvV1 : public AliAnalysisTaskSE {

  public:
    AliAnalysisTaskGammaConvV1();
    AliAnalysisTaskGammaConvV1(const char *name);
    virtual ~AliAnalysisTaskGammaConvV1();

    virtual void   UserCreateOutputObjects();
    virtual Bool_t Notify();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);
    void InitBack();

    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    void SetLightOutput(Bool_t flag ){fDoLightOutput = flag;}

    void SetIsHeavyIon(Int_t flag)                                { fIsHeavyIon                 = flag    ;}
    void SetIsMC(Int_t isMC)                                      { fIsMC                       = isMC    ;}
    void SetDoMesonAnalysis(Bool_t flag)                          { fDoMesonAnalysis            = flag    ;}
    void SetDoMesonQA(Int_t flag)                                 { fDoMesonQA                  = flag    ;}
    void SetDoPhotonQA(Int_t flag)                                { fDoPhotonQA                 = flag    ;}
    void SetDoClusterSelectionForTriggerNorm(Bool_t flag)         { fEnableClusterCutsForTrigger= flag    ;}
    void SetDoChargedPrimary(Bool_t flag)                         { fDoChargedPrimary           = flag    ;}
    void SetDoPlotVsCentrality(Bool_t flag)                       { fDoPlotVsCentrality         = flag    ;}
    void SetDoTHnSparse(Bool_t flag)                              { fDoTHnSparse                = flag    ;}
    void SetDoCentFlattening(Int_t flag)                          { fDoCentralityFlat           = flag    ;}
    void ProcessPhotonCandidates();
    void ProcessClusters();
    void CalculatePi0Candidates();
    void CalculateBackground();
    void CalculateBackgroundRP();
    void ProcessMCParticles();
    void ProcessAODMCParticles();
    void RelabelAODPhotonCandidates(Bool_t mode);
    void ProcessTruePhotonCandidates( AliAODConversionPhoton* TruePhotonCandidate);
    void ProcessTruePhotonCandidatesAOD( AliAODConversionPhoton* TruePhotonCandidate);
    void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
    void ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
    void RotateParticle(AliAODConversionPhoton *gamma);
    void RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP);
    void SetEventCutList(Int_t nCuts, TList *CutArray)          { fnCuts                        = nCuts     ;
                                                                  fEventCutArray                = CutArray  ;}
    void SetConversionCutList(Int_t nCuts, TList *CutArray)     { fnCuts                        = nCuts     ;
                                                                  fCutArray                     = CutArray  ;}
                                                                
    void SetMesonCutList(Int_t nCuts, TList *CutArray)          { fnCuts                        = nCuts     ;
                                                                  fMesonCutArray                = CutArray  ;}
    void SetClusterCutList(Int_t nCuts, TList *CutArray)        { fnCuts                        = nCuts     ;
                                                                  fClusterCutArray              = CutArray  ;}

    // BG HandlerSettings
    void SetMoveParticleAccordingToVertex(Bool_t flag)            {fMoveParticleAccordingToVertex = flag;}
    void FillPhotonCombinatorialBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode[], Int_t fDoPhotonQA, Double_t PhiParticle[]);
    void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    void UpdateEventByEventData();
    void SetLogBinningXTH2(TH2* histoRebin);
    Int_t GetSourceClassification(Int_t daughter, Int_t pdgCode);

    // Additional functions
    Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);
    void FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked);
    void FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist);
    
  protected:
    AliV0ReaderV1*                    fV0Reader;                                  //
    TString                           fV0ReaderName;
    Bool_t                            fDoLightOutput;                             // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    AliGammaConversionAODBGHandler**  fBGHandler;                                 //
    AliConversionAODBGHandlerRP**     fBGHandlerRP;                               //
    AliVEvent*                        fInputEvent;                                //
    AliMCEvent*                       fMCEvent;                                   //
    AliStack*                         fMCStack;                                   //
    TList**                           fCutFolder;                                 //
    TList**                           fESDList;                                   //
    TList**                           fBackList;                                  //
    TList**                           fMotherList;                                //
    TList**                           fPhotonDCAList;                             //
    TList**                           fMesonDCAList;                              //
    TList**                           fTrueList;                                  //
    TList**                           fMCList;                                    //
    TList**                           fHeaderNameList;                            //
    TList*                            fOutputContainer;                           //
    TClonesArray*                     fReaderGammas;                              //
    TList*                            fGammaCandidates;                           //
    TList*                            fEventCutArray;                             //
    TList*                            fCutArray;                                  //
    TList*                            fMesonCutArray;                             //
    TList*                            fClusterCutArray;                           //
    TH1F**                            hESDCaloGammaPt;                            //!
    TH1F**                            hESDConvGammaPt;                            //!
    TH1F**                            hESDConvGammaR;                             //!
    TH1F**                            hESDConvGammaEta;                           //!
    TH1F**                            hESDConvGammaPhi;                           //!
    TH1F**                            hESDConvGammaPsiPair;                       //!
    TH2F**                            hESDConvGammaPsiPairPt;                     //!
    TTree**                           tESDConvGammaPtDcazCat;                     //!
    Float_t                           fPtGamma;                                   //!
    Float_t                           fDCAzPhoton;                                //!
    Float_t                           fRConvPhoton;                               //!
    Float_t                           fEtaPhoton;                                 //!
    UChar_t                           iCatPhoton;                                 //!
    UChar_t                           iPhotonMCInfo;                              //!
                                      // 0: garbage,
                                      // 1: background
                                      // 2: secondary photon not from eta or k0s,
                                      // 3: secondary photon from eta, 
                                      // 4: secondary photon from k0s, 
                                      // 5: dalitz
                                      // 6: primary gamma
    TH2F**                            hESDMotherInvMassPt;                        //!
    THnSparseF**                      sESDMotherInvMassPtZM;                      //!
    TH2F**                            hESDMotherBackInvMassPt;                    //!
    THnSparseF**                      sESDMotherBackInvMassPtZM;                  //!
    TH2F**                            hESDMotherInvMassEalpha;                    //!
    TH2F**                            hESDMotherPi0PtY;                           //!
    TH2F**                            hESDMotherEtaPtY;                           //!
    TH2F**                            hESDMotherPi0PtAlpha;                       //!
    TH2F**                            hESDMotherEtaPtAlpha;                       //!
    TH2F**                            hESDMotherPi0PtOpenAngle;                   //!
    TH2F**                            hESDMotherEtaPtOpenAngle;                   //!
    TH2F**                            hESDMotherPi0LowPt;                         //!
    TH2F**                            hESDMotherPi0HighPt;                        //!
    THnSparseF**                      sPtRDeltaROpenAngle;                        //!
    TH1I**                            hMCHeaders;                                 //!
    TH1F**                            hMCAllGammaPt;                              //!
    TH1F**                            hMCDecayGammaPi0Pt;                         //!
    TH1F**                            hMCDecayGammaRhoPt;                         //!
    TH1F**                            hMCDecayGammaEtaPt;                         //!
    TH1F**                            hMCDecayGammaOmegaPt;                       //!
    TH1F**                            hMCDecayGammaEtapPt;                        //!
    TH1F**                            hMCDecayGammaPhiPt;                         //!
    TH1F**                            hMCDecayGammaSigmaPt;                       //!
    TH1F**                            hMCConvGammaPt;                             //!
    TH1F**                            hMCConvGammaR;                              //!
    TH1F**                            hMCConvGammaEta;                            //!
    TH1F**                            hMCPi0Pt;                                   //!
    TH1F**                            hMCPi0WOWeightPt;                           //! array of histos with unweighted pi0, pT
    TH1F**                            hMCPi0WOEvtWeightPt;                        //! array of histos without event weights pi0, pT
    TH1F**                            hMCEtaWOEvtWeightPt;                        //! array of histos without event weights eta, pT
    TH1F**                            hMCEtaPt;                                   //!
    TH1F**                            hMCEtaWOWeightPt;                           //!
    TH1F**                            hMCPi0WOWeightInAccPt;                      //!
    TH1F**                            hMCEtaWOWeightInAccPt;                      //!
    TH1F**                            hMCPi0InAccPt;                              //!
    TH1F**                            hMCEtaInAccPt;                              //!
    TH1F**                            hMCPi0WOEvtWeightInAccPt;                   //!
    TH1F**                            hMCEtaWOEvtWeightInAccPt;                   //!
    TH2F**                            hMCPi0PtY;                                  //!
    TH2F**                            hMCEtaPtY;                                  //!
    TH2F**                            hMCPi0PtAlpha;                              //!
    TH2F**                            hMCEtaPtAlpha;                              //!
    TH1F**                            hMCK0sPt;                                   //!
    TH2F**                            hMCSecPi0PtvsSource;                        //!
    TH2F**                            hMCSecPi0RvsSource;                         //!
    TH1F**                            hMCSecPi0Source;                            //!
    TH2F**                            hMCSecPi0InAccPtvsSource;                   //!
    TH1F**                            hMCSecEtaPt;                                //!
    TH1F**                            hMCSecEtaSource;                            //!
    TH2F**                            hMCPi0PtJetPt;                              //! array of histos with weighted pi0, pT, hardest jet pt
    TH2F**                            hMCEtaPtJetPt;                              //! array of histos with weighted eta, pT, hardest jet pt
    TH1F**                            hMCPhysicalPrimariesPt;                     //!
    TH1F**                            hMCPrimaryPionPlusPt;                       //!
    TH1F**                            hMCPrimaryPionMinusPt;                      //!
    TH1F**                            hMCPrimaryKaonPlusPt;                       //!
    TH1F**                            hMCPrimaryKaonMinusPt;                      //!
    TH1F**                            hMCPrimaryProtonPt;                         //!
    TH1F**                            hMCPrimaryAntiprotonPt;                     //!
    TH1F**                            hMCPrimaryPi0Pt;                            //!
    TH1F**                            hMCPrimaryEtaPt;                            //!
    TH2F**                            hESDTrueMotherInvMassPt;                    //!
    TH2F**                            hESDTruePrimaryMotherInvMassPt;             //!
    TH2F**                            hESDTruePrimaryMotherW0WeightingInvMassPt;  //!
    TProfile2D**                      pESDTruePrimaryMotherWeightsInvMassPt;      //!
    TH2F**                            hESDTruePrimaryPi0MCPtResolPt;              //!
    TH2F**                            hESDTruePrimaryEtaMCPtResolPt;              //!
    TH2F**                            hESDTrueSecondaryMotherInvMassPt;           //!
    TH2F**                            hESDTrueSecondaryMotherFromK0sInvMassPt;    //!
    TH1F**                            hESDTrueK0sWithPi0DaughterMCPt;             //!
    TH2F**                            hESDTrueSecondaryMotherFromEtaInvMassPt;    //!
    TH1F**                            hESDTrueEtaWithPi0DaughterMCPt;             //!
    TH2F**                            hESDTrueSecondaryMotherFromLambdaInvMassPt; //!
    TH1F**                            hESDTrueLambdaWithPi0DaughterMCPt;          //!
    TH2F**                            hESDTrueBckGGInvMassPt;                     //!
    TH2F**                            hESDTrueBckContInvMassPt;                   //!
    TH2F**                            hESDTruePi0PtY;                             //!
    TH2F**                            hESDTrueEtaPtY;                             //!
    TH2F**                            hESDTruePi0PtAlpha;                         //!
    TH2F**                            hESDTrueEtaPtAlpha;                         //!
    TH2F**                            hESDTruePi0PtOpenAngle;                     //!
    TH2F**                            hESDTrueEtaPtOpenAngle;                     //!
    TH2F**                            hESDTruePi0LowPt;                           //!
    TH2F**                            hESDTruePi0HighPt;                          //!
    TH2F**                            hESDTrueMotherDalitzInvMassPt;              //!
    TH1F**                            hESDTrueConvGammaPt;                        //!
    TH1F**                            hESDTrueConvGammaR;                         //!
    TH1F**                            hESDTrueConvGammaPtMC;                      //!
    TH1F**                            hESDTrueConvGammaRMC;                       //!
    TH1F**                            hESDTrueConvGammaEta;                       //!
    TH1F**                            hESDTrueConvGammaPsiPair;                   //!
    TH2F**                            hESDTrueConvGammaPsiPairPt;                 //!
    TH2F**                            hESDCombinatorialPt;                        //!
    TH2F**                            hESDCombinatorialPtDeltaPhi_ek;             //!
    TH2F**                            hESDCombinatorialPtDeltaPhi_ep;             //!
    TH2F**                            hESDCombinatorialPtDeltaPhi_epi;            //!
    TH2F**                            hESDCombinatorialPtDeltaPhi_pik;            //!
    TH2F**                            hESDCombinatorialPtDeltaPhi_pip;            //!
    TH1F**                            hESDTruePrimaryConvGammaPt;                 //!
    TH2F**                            hESDTruePrimaryConvGammaESDPtMCPt;          //!
    TH1F**                            hESDTrueSecondaryConvGammaPt;               //!
    TH1F**                            hESDTrueSecondaryConvGammaFromXFromK0sPt;   //!
    TH1F**                            hESDTrueSecondaryConvGammaFromXFromLambdaPt;//!
    TH2F**                            hESDTrueDalitzPsiPairDeltaPhi;              //!
    TH2F**                            hESDTrueGammaPsiPairDeltaPhi;               //!
    TH2F**                            hDoubleCountTruePi0InvMassPt;               //! array of histos with double counted pi0s, invMass, pT
    TH2F**                            hDoubleCountTrueEtaInvMassPt;               //! array of histos with double counted etas, invMass, pT
    TH2F**                            hDoubleCountTrueConvGammaRPt;               //! array of histos with double counted photons, R, pT
    vector<Int_t>                     vecDoubleCountTruePi0s;                     //! vector containing labels of validated pi0
    vector<Int_t>                     vecDoubleCountTrueEtas;                     //! vector containing labels of validated eta
    vector<Int_t>                     vecDoubleCountTrueConvGammas;               //! vector containing labels of validated photons
    TH1F**                            hMultipleCountTruePi0;                      //! array of histos how often TruePi0s are counted
    TH1F**                            hMultipleCountTrueEta;                      //! array of histos how often TrueEtas are counted
    TH1F**                            hMultipleCountTrueConvGamma;                //! array of histos how often TrueConvGammass are counted
    map<Int_t,Int_t>                  mapMultipleCountTruePi0s;                   //! map containing pi0 labels that are counted at least twice
    map<Int_t,Int_t>                  mapMultipleCountTrueEtas;                   //! map containing eta labels that are counted at least twice
    map<Int_t,Int_t>                  mapMultipleCountTrueConvGammas;             //! map containing photon labels that are counted at least twice
    TH1F**                            hNEvents;                                   //!
    TH1F**                            hNEventsWOWeight;                           //! array of histos with event information without event weights
    TH1F**                            hNGoodESDTracks;                            //!
    TH1F**                            hNEventsWeighted;                           //!
    TH1F**                            hNGoodESDTracksWeighted;                    //!
    TH1F**                            hVertexZ;                                   //!
    TH1F**                            hVertexZWeighted;                           //!
    TH1F**                            hCentrality;                                //!
    TH1F**                            hCentralityFlattened;                       //!
    TH2F**                            hCentralityVsPrimaryTracks;                 //!
    TH1F**                            hNGammaCandidates;                          //!
    TH2F**                            hNGoodESDTracksVsNGammaCandidates;          //!
    TH2F**                            fHistoSPDClusterTrackletBackground;         //! array of histos with SPD tracklets vs SPD clusters for background rejection
    TH1F**                            hNV0Tracks;                                 //!
    TProfile**                        fProfileEtaShift;                           //! array of profiles with eta shift
    TProfile**                        fProfileJetJetXSection;                     //! array of profiles with xsection for jetjet
    TH1F**                            fhJetJetNTrials;                            //! array of histos with ntrials for jetjet
    TProfile**                        hEtaShift;                                  //!
    TTree**                           tESDMesonsInvMassPtDcazMinDcazMaxFlag;      //!
    Float_t                           fInvMass;                                   //!
    Float_t                           fPt;                                        //!
    Float_t                           fDCAzGammaMin;                              //!
    Float_t                           fDCAzGammaMax;                              //!
    UChar_t                           iFlag;                                      //!
    UChar_t                           iMesonMCInfo;                               //!
                                      // 0: garbage,
                                      // 1: background
                                      // 2: secondary meson not from eta or k0s,
                                      // 3: secondary meson from eta, 
                                      // 4: secondary meson from k0s, 
                                      // 5: dalitz
                                      // 6: primary meson gamma-gamma-channel
    Double_t                          fEventPlaneAngle;                           // EventPlaneAngle
    TRandom3                          fRandom;                                    //
    Int_t                             fnGammaCandidates;                          //
    Double_t*                         fUnsmearedPx;                               //[fnGammaCandidates]
    Double_t*                         fUnsmearedPy;                               //[fnGammaCandidates]
    Double_t*                         fUnsmearedPz;                               //[fnGammaCandidates]
    Double_t*                         fUnsmearedE;                                //[fnGammaCandidates]
    Int_t*                            fMCStackPos;                                //[fnGammaCandidates]
    Int_t*                            fMCStackNeg;                                //[fnGammaCandidates]
    Int_t*                            fESDArrayPos;                               //[fnGammaCandidates]
    Int_t*                            fESDArrayNeg;                               //[fnGammaCandidates]
    Int_t                             fnCuts;                                     //
    Int_t                             fiCut;                                      //
    Bool_t                            fMoveParticleAccordingToVertex;             //
    Int_t                             fIsHeavyIon;                                //
    Bool_t                            fDoMesonAnalysis;                           //
    Int_t                             fDoMesonQA;                                 //
    Int_t                             fDoPhotonQA;                                //
    Bool_t                            fDoChargedPrimary;                          //
    Bool_t                            fDoPlotVsCentrality;                        //
    Bool_t                            fIsFromMBHeader;                            //
    Int_t                             fIsMC;                                      //
    Bool_t                            fDoTHnSparse;                               // flag for using THnSparses for background estimation
    Int_t                             fDoCentralityFlat;                          //flag for centrality flattening
    Double_t                          fWeightJetJetMC;                            // weight for Jet-Jet MC
    Double_t*                         fWeightCentrality;                          //[fnCuts], weight for centrality flattening
    Bool_t                            fEnableClusterCutsForTrigger;                //enables ClusterCuts for Trigger
    
  private:

    AliAnalysisTaskGammaConvV1(const AliAnalysisTaskGammaConvV1&); // Prevent copy-construction
    AliAnalysisTaskGammaConvV1 &operator=(const AliAnalysisTaskGammaConvV1&); // Prevent assignment
    ClassDef(AliAnalysisTaskGammaConvV1, 33);
};

#endif
