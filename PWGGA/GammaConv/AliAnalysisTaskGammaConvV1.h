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
#include "AliAnalysisTaskConvJet.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
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
    void SetFileNameBDT(TString filename) { fFileNameBDT = filename.Data() ;}
    void InitializeBDT();
    void ProcessPhotonBDT();
    void ProcessClusters();
    void ProcessJets();
    void ProcessPhotonsHighPtHadronAnalysis();
    void CalculatePi0Candidates();
    void CalculateBackground();
    void CalculateBackgroundSwapp();
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

    void SetDoMaterialBudgetWeightingOfGammasForTrueMesons(Bool_t flag) {fDoMaterialBudgetWeightingOfGammasForTrueMesons = flag;}

    // BG HandlerSettings
    void SetMoveParticleAccordingToVertex(Bool_t flag)            {fMoveParticleAccordingToVertex = flag;}
    void FillPhotonCombinatorialBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode[], Double_t PhiParticle[]);
    void FillPhotonCombinatorialMothersHistESD(TParticle *daughter,TParticle *mother);
    void FillPhotonCombinatorialMothersHistAOD(AliAODMCParticle *daughter, AliAODMCParticle* motherCombPart);
    void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    void UpdateEventByEventData();
    void SetLogBinningXTH2(TH2* histoRebin);
    Int_t GetSourceClassification(Int_t daughter, Int_t pdgCode);

    // Additional functions
    Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);
    void FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked);
    void FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist);
    Double_t GetOriginalInvMass(const AliConversionPhotonBase * photon, AliVEvent * event) const;

  protected:
    AliV0ReaderV1*                    fV0Reader;                                  //
    TString                           fV0ReaderName;
    Bool_t                            fDoLightOutput;                             // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    AliGammaConversionAODBGHandler**  fBGHandler;                                 //
    AliConversionAODBGHandlerRP**     fBGHandlerRP;                               //
    AliVEvent*                        fInputEvent;                                //
    AliMCEvent*                       fMCEvent;                                   //
    TList**                           fCutFolder;                                 //
    TList**                           fESDList;                                   //
    TList**                           fBackList;                                  //
    TList**                           fMotherList;                                //
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
    AliAnalysisTaskConvJet*           fConvJetReader;                             //
    Bool_t                            fDoJetAnalysis;                             //
    Bool_t                            fDoIsolatedAnalysis;                        //
    Bool_t                            fDoHighPtHadronAnalysis;                    //
    Bool_t                            fDoJetQA;                                   //
    TList**                           fJetHistograms;                             //
    TList**                           fTrueJetHistograms;                         //
    Int_t                             fMaxPtNearEMCalPlace;                       //
    Bool_t                            fJetNearEMCal;                              //
    TH1F**                            fHistoCaloGammaPt;                          //!
    TH1F**                            fHistoCaloGammaE;                           //!
    TH1F**                            fHistoConvGammaPt;                          //!
    TH2F**                            fHistoConvGammaPtwithHighPtHadron;          //!
    TH2F**                            fHistoConvGammaPtwithoutHighPtHadron;       //!
    TH1F**                            fHistoConvGammaR;                           //!
    TH1F**                            fHistoConvGammaEta;                         //!
    TH1F**                            fHistoConvGammaPhi;                         //!
    TH2F**                            fHistoConvGammaPsiPairPt;                   //!
    TH1F**                            fHistoConvGammaInvMass;                     //!
    TH1F**                            fHistoConvGammaInvMassReco;                 //!
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
    TH2F**                            fHistoMotherInvMassPt;                        //!
    TH2F**                            fHistoMotherInvMassPtIso;                     //!
    TH2F**                            fHistoMotherInvMassPtNonIso;                  //!
    TH2F**                            fHistoMotherInvMassPtMCRecIsoTrueNonIso;      //!
    TH2F**                            fHistoMotherInvMassPtMCRecNonIsoTrueIso;      //!
    TH2F**                            fHistoMotherEisoPt;                           //!
    TH2F**                            fHistoMotherRisoPt;                           //!
    TH2F**                            fHistoMotherNtracksIsoPt;                     //!
    THnSparseF**                      sESDMotherInvMassPtZM;                        //!
    TH2F**                            fHistoMotherBackInvMassPt;                    //!
    THnSparseF**                      sESDMotherBackInvMassPtZM;                    //!
    TH2F**                            fHistoMotherInvMassEalpha;                    //!
    TH2F**                            fHistoMotherPi0PtY;                           //!
    TH2F**                            fHistoMotherEtaPtY;                           //!
    TH2F**                            fHistoMotherPi0PtAlpha;                       //!
    TH2F**                            fHistoMotherEtaPtAlpha;                       //!
    TH2F**                            fHistoMotherPi0PtOpenAngle;                   //!
    TH2F**                            fHistoMotherEtaPtOpenAngle;                   //!
    THnSparseF**                      sPtRDeltaROpenAngle;                          //!
    TH1I**                            fHistoMCHeaders;                                 //!
    TH1F**                            fHistoMCAllGammaPt;                              //!
    TH2F**                            fHistoMCAllSecondaryGammaPt;                     //!
    TH1F**                            fHistoMCDecayGammaPi0Pt;                         //!
    TH1F**                            fHistoMCDecayGammaRhoPt;                         //!
    TH1F**                            fHistoMCDecayGammaEtaPt;                         //!
    TH1F**                            fHistoMCDecayGammaOmegaPt;                       //!
    TH1F**                            fHistoMCDecayGammaEtapPt;                        //!
    TH1F**                            fHistoMCDecayGammaPhiPt;                         //!
    TH1F**                            fHistoMCDecayGammaSigmaPt;                       //!
    TH1F**                            fHistoMCConvGammaPt;                             //!
    TH2F**                            fHistoMCSecondaryConvGammaPt;                    //!
    TH1F**                            fHistoMCConvGammaR;                              //!
    TH1F**                            fHistoMCConvGammaEta;                            //!
    TH1F**                            fHistoMCPi0Pt;                                   //!
    TH1F**                            fHistoMCPi0WOWeightPt;                           //! array of histos with unweighted pi0, pT
    TH1F**                            fHistoMCPi0WOEvtWeightPt;                        //! array of histos without event weights pi0, pT
    TH1F**                            fHistoMCEtaWOEvtWeightPt;                        //! array of histos without event weights eta, pT
    TH1F**                            fHistoMCEtaPt;                                   //!
    TH1F**                            fHistoMCEtaWOWeightPt;                           //!
    TH1F**                            fHistoMCPi0WOWeightInAccPt;                      //!
    TH1F**                            fHistoMCEtaWOWeightInAccPt;                      //!
    TH1F**                            fHistoMCPi0InAccPt;                              //!
    TH1F**                            fHistoMCEtaInAccPt;                              //!
    TH1F**                            fHistoMCPi0WOEvtWeightInAccPt;                   //!
    TH1F**                            fHistoMCEtaWOEvtWeightInAccPt;                   //!
    TH2F**                            fHistoMCPi0PtY;                                  //!
    TH2F**                            fHistoMCEtaPtY;                                  //!
    TH2F**                            fHistoMCPi0PtAlpha;                              //!
    TH2F**                            fHistoMCEtaPtAlpha;                              //!
    TH2F**                            fHistoMCPrimaryPtvsSource;                       //!
    TH2F**                            fHistoMCSecPi0PtvsSource;                        //!
    TH2F**                            fHistoMCSecPi0RvsSource;                         //!
    TH1F**                            fHistoMCSecPi0Source;                            //!
    TH2F**                            fHistoMCSecPi0InAccPtvsSource;                   //!
    TH1F**                            fHistoMCSecEtaPt;                                //!
    TH1F**                            fHistoMCSecEtaSource;                            //!
    TH2F**                            fHistoMCPi0PtJetPt;                              //! array of histos with weighted pi0, pT, hardest jet pt
    TH2F**                            fHistoMCEtaPtJetPt;                              //! array of histos with weighted eta, pT, hardest jet pt
    TH1F**                            fHistoMCPhysicalPrimariesPt;                  //!
    TH2F**                            fHistoTrueMotherInvMassPt;                    //!
    TH2F**                            fHistoTruePrimaryMotherInvMassPt;             //!
    TH2F**                            fHistoTruePrimaryMotherW0WeightingInvMassPt;  //!
    TProfile2D**                      pESDTruePrimaryMotherWeightsInvMassPt;        //!
    TH2F**                            fHistoTruePrimaryPi0MCPtResolPt;              //!
    TH2F**                            fHistoTruePrimaryEtaMCPtResolPt;              //!
    TH2F**                            fHistoTrueSecondaryMotherInvMassPt;           //!
    TH2F**                            fHistoTrueSecondaryMotherFromK0sInvMassPt;    //!
    TH1F**                            fHistoTrueK0sWithPi0DaughterMCPt;             //!
    TH2F**                            fHistoTrueSecondaryMotherFromK0lInvMassPt;    //!
    TH1F**                            fHistoTrueK0lWithPi0DaughterMCPt;             //!
    TH2F**                            fHistoTrueSecondaryMotherFromEtaInvMassPt;    //!
    TH1F**                            fHistoTrueEtaWithPi0DaughterMCPt;             //!
    TH2F**                            fHistoTrueSecondaryMotherFromLambdaInvMassPt; //!
    TH1F**                            fHistoTrueLambdaWithPi0DaughterMCPt;          //!
    TH2F**                            fHistoTrueBckGGInvMassPt;                     //!
    TH2F**                            fHistoTrueBckContInvMassPt;                   //!
    TH2F**                            fHistoTruePi0PtY;                             //!
    TH2F**                            fHistoTrueEtaPtY;                             //!
    TH2F**                            fHistoTruePi0PtAlpha;                         //!
    TH2F**                            fHistoTrueEtaPtAlpha;                         //!
    TH2F**                            fHistoTruePi0PtOpenAngle;                     //!
    TH2F**                            fHistoTrueEtaPtOpenAngle;                     //!
    TH2F**                            fHistoTrueMotherDalitzInvMassPt;              //!
    TH1F**                            fHistoTrueConvGammaPt;                        //!
    TH1F**                            fHistoTrueConvGammaR;                         //!
    TH1F**                            fHistoTrueConvGammaPtMC;                      //!
    TH1F**                            fHistoTrueConvGammaRMC;                       //!
    TH1F**                            fHistoTrueConvGammaEta;                       //!
    TH2F**                            fHistoTrueConvGammaPsiPairPt;                 //!
    TH1F**                            fHistoTrueConvGammaInvMass;                   //!
    TH1F**                            fHistoTrueConvGammaInvMassReco;               //!
    TH2F**                            fHistoCombinatorialPt;                        //!
    TH3F**                            fHistoCombinatorialMothersPt;                 //!
    TH2F**                            fHistoCombinatorialPtDeltaPhi_ek;             //!
    TH2F**                            fHistoCombinatorialPtDeltaPhi_ep;             //!
    TH2F**                            fHistoCombinatorialPtDeltaPhi_epi;            //!
    TH2F**                            fHistoCombinatorialPtDeltaPhi_pik;            //!
    TH2F**                            fHistoCombinatorialPtDeltaPhi_pip;            //!
    TH1F**                            fHistoTruePrimaryConvGammaPt;                 //!
    TH2F**                            fHistoTrueSecondaryConvGammaPt;                       //!
    TH2F**                            fHistoTrueSecondaryConvGammaMCPt;                     //!
    TH2F**                            fHistoTruePrimaryConvGammaESDPtMCPt;                  //!
    TH2F**                            fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt;    //!
    TH2F**                            fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt;    //!
    TH2F**                            fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt; //!
    TH2F**                            fHistoTrueDalitzPsiPairDeltaPhi;                      //!
    TH2F**                            fHistoTrueGammaPsiPairDeltaPhi;                       //!
    TH2F**                            fHistoDoubleCountTruePi0InvMassPt;               //! array of histos with double counted pi0s, invMass, pT
    TH2F**                            fHistoDoubleCountTrueEtaInvMassPt;               //! array of histos with double counted etas, invMass, pT
    TH2F**                            fHistoDoubleCountTrueConvGammaRPt;               //! array of histos with double counted photons, R, pT
    vector<Int_t>                     vecDoubleCountTruePi0s;                     //! vector containing labels of validated pi0
    vector<Int_t>                     vecDoubleCountTrueEtas;                     //! vector containing labels of validated eta
    vector<Int_t>                     vecDoubleCountTrueConvGammas;               //! vector containing labels of validated photons
    TH1F**                            fHistoMultipleCountTruePi0;                      //! array of histos how often TruePi0s are counted
    TH1F**                            fHistoMultipleCountTrueEta;                      //! array of histos how often TrueEtas are counted
    TH1F**                            fHistoMultipleCountTrueConvGamma;                //! array of histos how often TrueConvGammass are counted
    map<Int_t,Int_t>                  mapMultipleCountTruePi0s;                   //! map containing pi0 labels that are counted at least twice
    map<Int_t,Int_t>                  mapMultipleCountTrueEtas;                   //! map containing eta labels that are counted at least twice
    map<Int_t,Int_t>                  mapMultipleCountTrueConvGammas;             //! map containing photon labels that are counted at least twice
    TH1F**                            fHistoNEvents;                                   //!
    TH1F**                            fHistoNEventsWOWeight;                           //! array of histos with event information without event weights
    TH1F**                            fHistoNGoodESDTracks;                            //!
    TH1F**                            fHistoNEventsWeighted;                           //!
    TH1F**                            fHistoNGoodESDTracksWeighted;                    //!
    TH1F**                            fHistoVertexZ;                                   //!
    TH1F**                            fHistoVertexZWeighted;                           //!
    Int_t                             fDoCentralityFlat;                               //flag for centrality flattening
    TH1F**                            fHistoCentrality;                                //!
    TH1F**                            fHistoCentralityFlattened;                       //!
    TH2F**                            fHistoCentralityVsPrimaryTracks;                 //!
    TH1F**                            fHistoNGammaCandidates;                          //!
    TH2F**                            fHistoNGoodESDTracksVsNGammaCandidates;          //!

    TH1F**                            fHistoPtJet;                                          // Histogram of Jet Pt
    TH1F**                            fHistoJetEta;                                         // Histogram of Jet Eta
    TH1F**                            fHistoJetPhi;                                         // Histogram of Jet Phi
    TH1F**                            fHistoJetArea;                                        // Histogram of Jet Area
    TH1F**                            fHistoNJets;                                          // Histogram of number of jets
    TH1F**                            fHistoEventwJets;                                     // Histogram of number of events with jets > 0
    TH1F**                            fHistoJetPi0PtRatio;                                  // Histogram of PtPi0/PtJet
    TH1F**                            fHistoDoubleCounting;                                 // Histogram if NM candidates are defined within multiple jets
    TH2F**                            fHistoJetMotherInvMassPt;                             // Histogram of NM candidates with a jet in the event
    TH2F**                            fHistoPi0InJetMotherInvMassPt;                        // Histogram of NM candidates that are inside a jet
    TH2F**                            fHistoMotherBackJetInvMassPt;                         // Histogram of Backgrouns candidates that are involved with jets
    TH2F**                            fHistoRJetPi0Cand;                                    // Histogram of RJetPi0Cand vs Pt
    TH2F**                            fHistoEtaPhiJetPi0Cand;                               // Histogram of delta eta and delta phi distr between jet and NM candidates
    TH2F**                            fHistoEtaPhiJetWithPi0Cand;                           // Histogram of delta eta and delta phi distr when pi0 is inside a jet
    TH2F**                            fHistoJetFragmFunc;                                   // Histogram to determine fragmentation function
    TH2F**                            fHistoJetFragmFuncZInvMass;                           // Histogram of Inv Mass distribution with z
    TH2F**                            fHistoTruevsRecJetPt;                                 // Histogram of true jet pt vs reconstructed jet pt
    TH2F**                            fHistoTrueJetMotherInvMassPt;                         // Histogram of true pi0s in an event with a jet
    TH2F**                            fHistoTrueInJetMotherInvMassPt;                       // Histogram of true pi0s in a jet
    TH2F**                            fHistoTruePrimaryJetInvMassPt;                        // Histogram of true primary pi0s in an event with a jet
    TH2F**                            fHistoTruePrimaryinJetInvMassPt;                      // Histogram of true primary pi0s in a jet
    TH2F**                            fHistoTruePrimaryInJetInvMassTruePt;                  // Histogram of true primary pi0s in a jet with their true pt
    TH1F**                            fHistoTrueDoubleCountingJet;                          // Histogram of when a true pi0 is defined to be in multiple jets
    TH2F**                            fHistoTrueJetFragmFunc;                               // Histogram to determine true fragmentation function
    TH2F**                            fHistoTrueJetFragmFuncZInvMass;                       // Histogram to determine true Inv Mass distribution with z
    TH1F**                            fHistoMCPi0JetInAccPt;                                // Histogram with weighted pi0 in a jet event in acceptance, pT
    TH1F**                            fHistoMCPi0inJetInAccPt;                              // Histogram with weighted pi0 in a jet in acceptance, pT
    TH1F**                            fHistoMCEtaJetInAccPt;                                // Histogram with weighted eta in a jet event in acceptance, pT
    TH1F**                            fHistoMCEtainJetInAccPt;                              // Histogram with weighted eta in a jet in acceptance, pT
    TH1F**                            fHistoMCPi0JetEventGenerated;                         // Histogram with mesons in a jet event generated, pT
    TH1F**                            fHistoMCPi0inJetGenerated;                            // Histogram with mesons in a jet generated, pT
    TH1F**                            fHistoMCEtaJetEventGenerated;                         // Histogram with mesons in a jet event generated, pT
    TH1F**                            fHistoMCEtainJetGenerated;                            // Histogram with mesons in a jet generated, pT
    TH2F**                            fHistoTrueSecondaryFromK0sJetInvMassPt;               // Histogram with validated secondary mothers from K0s in an event with a jet, invMass, pt
    TH2F**                            fHistoTrueSecondaryFromK0sinJetInvMassPt;             // Histogram with validated secondary mothers from K0s in a jet, invMass, pt
    TH2F**                            fHistoTrueSecondaryFromLambdaJetInvMassPt;            // Histogram with validated secondary mothers from lambda in an event with a jet, invMass, pt
    TH2F**                            fHistoTrueSecondaryFromLambdainJetInvMassPt;          // Histogram with validated secondary mothers from lambda in a jet, invMass, pt
    TH2F**                            fHistoTrueSecondaryFromK0lJetInvMassPt;               // Histogram with validated secondary mothers from K0l in an event with a jet, invMass, pt
    TH2F**                            fHistoTrueSecondaryFromK0linJetInvMassPt;             // Histogram with validated secondary mothers from K0l in a jet, invMass, pt
    TH2F**                            fHistoTrueSecondaryInvJetMassPt;                      // Histogram with validated secondary mothers in an event with a jet, invMass, pt
    TH2F**                            fHistoTrueSecondaryInvinJetMassPt;                    // Histogram with validated secondary mothers in a jet, invMass, pt
    TH2F**                            fHistoUnfoldingAsData;                                // Histogram to use for jet pi0 unfolding
    TH2F**                            fHistoUnfoldingMissed;                                // Histogram to use for jet pi0 unfolding
    TH2F**                            fHistoUnfoldingReject;                                // Histogram to use for jet pi0 unfolding
    TH2F**                            fHistoUnfoldingAsDataInvMassZ;                        // Histogram to use for jet pi0 unfolding for z inmass
    TH2F**                            fHistoUnfoldingMissedInvMassZ;                        // Histogram to use for jet pi0 unfolding for z inmass
    TH2F**                            fHistoUnfoldingRejectInvMassZ;                        // Histogram to use for jet pi0 unfolding for z inmass

    vector<Double_t>                  fVectorJetPt;                                         // Vector of JetPt
    vector<Double_t>                  fVectorJetPx;                                         // Vector of JetPx
    vector<Double_t>                  fVectorJetPy;                                         // Vector of JetPy
    vector<Double_t>                  fVectorJetPz;                                         // Vector of JetPz
    vector<Double_t>                  fVectorJetEta;                                        // Vector of JetEta
    vector<Double_t>                  fVectorJetPhi;                                        // Vector of JetPhi
    vector<Double_t>                  fVectorJetArea;                                       // Vector of JetArea
    vector<Double_t>                  fTrueVectorJetPt;                                     // Vector of True JetPt
    vector<Double_t>                  fTrueVectorJetPx;                                     // Vector of True JetPx
    vector<Double_t>                  fTrueVectorJetPy;                                     // Vector of True JetPy
    vector<Double_t>                  fTrueVectorJetPz;                                     // Vector of True JetPz
    vector<Double_t>                  fTrueVectorJetEta;                                    // Vector of True JetEta
    vector<Double_t>                  fTrueVectorJetPhi;                                    // Vector of True JetPhi

    TH2F**                            fHistoSPDClusterTrackletBackground;         //! array of histos with SPD tracklets vs SPD clusters for background rejection
    TH2F**                            fHistoV0MultVsNumberTPCoutTracks;           //! correlation V=Mult vs number TPC out Tracks
    TH1F**                            fHistoNV0Tracks;                            //!
    TH1F**                            fHistoBDToutput;                            //!
    TH1F**                            fHistoBDToutputPt;                          //!
    TH1F**                            fHistoBDToutputMCTrue;                      //!
    TProfile**                        fProfileEtaShift;                           //! array of profiles with eta shift
    TProfile**                        fProfileJetJetXSection;                     //! array of profiles with xsection for jetjet
    TH1F**                            fhJetJetNTrials;                            //! array of histos with ntrials for jetjet
    TH1F**                            fHistoEventSphericity;                      //! array of histos with event Sphericity
    TProfile**                        fHistoEtaShift;                             //!
    TTree**                           tESDMesonsInvMassPtDcazMinDcazMaxFlag;      //!
    Float_t                           fInvMass;                                   //!
    Float_t*                          fBDTvariable;                               //!
    TMVA::Reader*                     fBDTreader;                                 //!
    TString                           fFileNameBDT;                               //!
    Bool_t                            fEnableBDT;                                 //!
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
    Int_t*                            fMCEventPos;                                //[fnGammaCandidates]
    Int_t*                            fMCEventNeg;                                //[fnGammaCandidates]
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
    Bool_t                            fIsFromSelectedHeader;                      //
    Int_t                             fIsMC;                                      //
    Bool_t                            fDoTHnSparse;                               // flag for using THnSparses for background estimation
    Double_t                          fWeightJetJetMC;                            // weight for Jet-Jet MC
    Double_t*                         fWeightCentrality;                          //[fnCuts], weight for centrality flattening
    Bool_t                            fEnableClusterCutsForTrigger;               //enables ClusterCuts for Trigger
    Bool_t                            fDoMaterialBudgetWeightingOfGammasForTrueMesons;
    TTree*                            tBrokenFiles;                               // tree for keeping track of broken files
    TObjString*                       fFileNameBroken;                            // string object for broken file name
    Bool_t                            fFileWasAlreadyReported;                    // to store if the current file was already marked broken 
    TClonesArray*                     fAODMCTrackArray;       //! pointer to track array

  private:

    AliAnalysisTaskGammaConvV1(const AliAnalysisTaskGammaConvV1&); // Prevent copy-construction
    AliAnalysisTaskGammaConvV1 &operator=(const AliAnalysisTaskGammaConvV1&); // Prevent assignment
    ClassDef(AliAnalysisTaskGammaConvV1, 50);
};

#endif
