#ifndef ALIANLYSISTASKPI0V2CALO_cxx
#define ALIANLYSISTASKPI0V2CALO_cxx

#include "AliAnalysisTaskSE.h"
#include "AliOADBContainer.h"
#include "AliAODVZERO.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliAnalysisTaskJetOutlierRemoval.h"
#include "AliAnalysisManager.h"
#include "TGrid.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include "TGenPhaseSpace.h"
#include <vector>
#include <map>

class AliAnalysisTaskPi0v2Calo : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskPi0v2Calo();
    AliAnalysisTaskPi0v2Calo(const char *name);
    virtual ~AliAnalysisTaskPi0v2Calo();

    virtual void   UserCreateOutputObjects();
    virtual bool Notify();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);
    void InitBack();

    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    void SetCaloTriggerHelperName(TString name){fCaloTriggerHelperName=name; return;}
    void SetIsHeavyIon(int flag){
      fIsHeavyIon = flag;
    }
    void SetPeriodName(TString name){fPeriod=name; return;}
    void SetInAndOutPlane(bool flag) {fUseInOutPlane=flag; return;}
    void SetFlowQa(bool flag){IsQAVZERO=flag; return;}

    // base functions for selecting photon and meson candidates in reconstructed data
    void ProcessClusters();
    void ProcessConversionCandidates();
    void CalculatePi0Candidates();

    // MC functions
    void SetIsMC(int isMC){fIsMC=isMC;}
    void ProcessMCParticles(int isCurrentEventSelected = 0);
    void ProcessAODMCParticles(int isCurrentEventSelected = 0);
    void ProcessTrueClusterCandidates( AliAODConversionPhoton* TruePhotonCandidate, double clusterM02);
    void ProcessTrueClusterCandidatesAOD( AliAODConversionPhoton* TruePhotonCandidate, double clusterM02);
    void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
    void ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);

    // switches for additional analysis streams or outputs
    void SetLightOutput(int flag){fDoLightOutput = flag;}
    void SetDoMesonAnalysis(bool flag){fDoMesonAnalysis = flag;}
    void SetAllowOverlapHeaders( bool allowOverlapHeader ) {fAllowOverlapHeaders = allowOverlapHeader;}
    void SetDoPi0Only(bool flag){fDoPi0Only = flag;}
    void SetDoPrimaryTrackMatching( bool flag ){ fDoPrimaryTrackMatching = flag;}

    void SetInOutTimingCluster(double min, double max){
      fDoInOutTimingCluster = kTRUE; fMinTimingCluster = min; fMaxTimingCluster = max;
      return;
    }

      // Setting the cut lists for the conversion photons
    void SetEventCutList(int nCuts, TList *CutArray){
      fnCuts = nCuts;
      fEventCutArray = CutArray;
    }

      // Setting the cut lists for the calo photons
    void SetCaloCutList(int nCuts, TList *CutArray){
      fnCuts = nCuts;
      fClusterCutArray = CutArray;
    }

    // Setting the cut lists for the meson
    void SetMesonCutList(int nCuts, TList *CutArray){
      fnCuts = nCuts;
      fMesonCutArray = CutArray;
    }

    // BG HandlerSettings
    void CalculateBackground();
    void CalculateBackgroundSwapp();
    void UpdateEventByEventData();

    // Additional functions for convenience
    int GetSourceClassification(int daughter, int pdgCode);

    bool CheckVectorForDoubleCount(vector<int> &vec, int tobechecked);
    void FillMultipleCountHistoAndClear(map<int,int> &ma, TH1F* hist);

    // Function to enable MC label sorting
    void SetEnableSortingOfMCClusLabels(bool enableSort)  { fEnableSortForClusMC   = enableSort; }

    // Function to enable local debugging mode
    void SetLocalDebugFlag(int iF) {fLocalDebugFlag = iF;}

    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}

    void EventDebugMethod();
    void DebugMethod(AliAODConversionMother *pi0cand, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1);
    void DebugMethodPrint1(AliAODConversionMother *pi0cand, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1);

    void SetTrackMatcherRunningMode(int mode){fTrackMatcherRunningMode = mode;}

    int CheckClustersForMCContribution(int mclabel, bool leading = kFALSE);

    void LoadCalibrationContainers();
    bool LoadVZEROCalibration();
    bool LoadCalibHistForThisRun();
    bool GetVZEROPlane();
    double GetEventPlane(double qx, double qy, double harmonic);
    void FillFlowHistograms(AliAODConversionMother *pi0cand);
    void FillFlowBackHistograms(AliAODConversionMother *pi0cand, double weight);

  protected:
    AliV0ReaderV1*        fV0Reader;                                            // basic photon Selection Task
    TString               fV0ReaderName;
    TString               fCaloTriggerHelperName;
    TString               fCorrTaskSetting;
    TString               fPeriod;                                              // name of the period
    AliGammaConversionAODBGHandler**  fBGHandler;                               // BG handler for Conversion
    AliVEvent*            fInputEvent;                                          // current event
    AliMCEvent*           fMCEvent;                                             // corresponding MC event
    TList**               fCutFolder;                                           // Array of lists for containers belonging to cut
    TList**               fESDList;                                             // Array of lists with histograms with reconstructed properties
    TList**               fQAList;                                              // Array of lists with histograms for QA
    TList**               fBackList;                                            // Array of lists with BG THnSparseF
    TList**               fMotherList;                                          // Array of lists with Signal THnSparseF
    TList**               fTrueList;                                            // Array of lists with histograms with MC validated reconstructed properties
    TList**               fMCList;                                              // Array of lists with histograms with pure MC information
    TList*                fOutputContainer;                                     // Output container
    TClonesArray*         fReaderGammas;                                        // Array with conversion photons selected by V0Reader Cut
    TList*                fGammaCandidates;                                     // current list of photon candidates
    TList*                fClusterCandidates;                                   //! current list of cluster candidates
    TList*                fEventCutArray;                                       // List with Event Cuts
    AliConvEventCuts*     fEventCuts;                                           // EventCutObject
    TList*                fClusterCutArray;                                     // List with Cluster Cuts
    AliCaloPhotonCuts*    fCaloPhotonCuts;                                      // CaloPhotonCutObject
    TList*                fMesonCutArray;                                       // List with Meson Cuts
    AliConversionMesonCuts*   fMesonCuts;                                       // MesonCutObject
    AliAnalysisTaskJetOutlierRemoval*   fOutlierJetReader;                      // JetReader
    AliConversionPhotonCuts*  fConversionCuts;                                  // ConversionPhotonCutObject
    AliCaloTriggerMimicHelper**     fCaloTriggerMimicHelper;                    //!Array wich points to AliCaloTriggerMimicHelper for each Event Cut
    map<TString, bool>  fSetEventCutsOutputlist;                                //! Store, if Output list for Event Cut has already been added

    //histograms for mesons reconstructed quantities
    TH3F**                fHistoMotherInvMassPtPhiV0A;                          //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt, V0A EP angle
    TH3F**                fHistoMotherBackInvMassPtPhiV0A;                      //! array of histogram with BG for mixed event photon pairs, inv Mass, pt, V0A EP angle
    TH3F**                fHistoMotherInvMassPtPhiV0C;                          //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt, V0C EP angle
    TH3F**                fHistoMotherBackInvMassPtPhiV0C;                      //! array of histogram with BG for mixed event photon pairs, inv Mass, pt, V0C EP angle
    TH2F**                fHistoMotherInvMassPtV0CInPlane;                      //! array of histogram with signal + BG for same event photon pairs In Plane V0C, inv Mass, pt
    TH2F**                fHistoMotherInvMassPtV0AInPlane;                      //! array of histogram with signal + BG for same event photon pairs In Plane V0A, inv Mass, pt
    TH2F**                fHistoMotherInvMassPtV0COutPlane;                     //! array of histogram with signal + BG for same event photon pairs Out Plane V0C, inv Mass, pt
    TH2F**                fHistoMotherInvMassPtV0AOutPlane;                     //! array of histogram with signal + BG for same event photon pairs Out Plane V0A, inv Mass, pt
    TH2F**                fHistoMotherBackInvMassPtV0CInPlane;                  //! array of histogram with only BG for same event photon pairs In Plane V0C, inv Mass, pt
    TH2F**                fHistoMotherBackInvMassPtV0AInPlane;                  //! array of histogram with only BG for same event photon pairs In Plane V0A, inv Mass, pt
    TH2F**                fHistoMotherBackInvMassPtV0COutPlane;                 //! array of histogram with only BG for same event photon pairs Out Plane V0C, inv Mass, pt
    TH2F**                fHistoMotherBackInvMassPtV0AOutPlane;                 //! array of histogram with only BG for same event photon pairs Out Plane V0A, inv Mass, pt


    // histograms for rec photon clusters
    TH1F**                fHistoClusGammaPt;                                    //! array of histos with cluster, pt
    TH1F**                fHistoClusGammaE;                                     //! array of histos with cluster, E
    TH1F**                fHistoClusOverlapHeadersGammaPt;                      //! array of histos with cluster, pt overlapping with other headers
    TH1F**                fHistoClusOverlapMBHeaderGammaPt;                     //! array of histos with cluster, pt overlapping with MB header
    TH1F**                fHistoClusAllHeadersGammaPt;                          //! array of histos with cluster, pt all headers
    TH1F**                fHistoClusRejectedHeadersGammaPt;                     //! array of histos with cluster, pt rejected with other headers
    TH2F**                fHistoClusGammaPtM02;                                 //! array of histos with cluster M02 vs. pt
    //histograms for pure MC quantities
    TH1I**                fHistoMCHeaders;                                      //! array of histos for header names
    TH1D**                fHistoMCEventsTrigg;                                  //! array of histos with number of accepted and rejected events selected on MC based trigger (important for mult. dep INEL>0)
    TH1F**                fHistoMCGammaPtNotTriggered;                          //! array of histos with weighted gamm in not triggered events, pT
    TH1F**                fHistoMCGammaPtNoVertex;                              //! array of histos with weighted gamm in not triggered events, pT    
    TH1F**                fHistoMCAllGammaPt;                                   //! array of histos with all gamma, pT
    TH2F**                fHistoMCAllSecondaryGammaPt;                          //! array of histos with all secondary gamma, pT
    TH1F**                fHistoMCDecayGammaPi0Pt;                              //! array of histos with decay gamma from pi0, pT
    TH1F**                fHistoMCDecayGammaRhoPt;                              //! array of histos with decay gamma from rho, pT
    TH1F**                fHistoMCDecayGammaEtaPt;                              //! array of histos with decay gamma from eta, pT
    TH1F**                fHistoMCDecayGammaOmegaPt;                            //! array of histos with decay gamma from omega, pT
    TH1F**                fHistoMCDecayGammaEtapPt;                             //! array of histos with decay gamma from eta', pT
    TH1F**                fHistoMCDecayGammaPhiPt;                              //! array of histos with decay gamma from phi, pT
    TH1F**                fHistoMCDecayGammaSigmaPt;                            //! array of histos with decay gamma from Sigma0, pT
    TH1F**                fHistoMCPi0Pt;                                        //! array of histos with weighted pi0, pT
    TH1F**                fHistoMCPi0PtNotTriggered;                            //! array of histos with weighted pi0 in not triggered events, pT
    TH1F**                fHistoMCPi0PtNoVertex;                                //! array of histos with weighted pi0 in not triggered events, pT
    TH1F**                fHistoMCPi0WOWeightPt;                                //! array of histos with unweighted pi0, pT
    TH1F**                fHistoMCPi0WOEvtWeightPt;                             //! array of histos without event weights pi0, pT
    TH1F**                fHistoMCEtaPt;                                        //! array of histos with weighted eta, pT
    TH1F**                fHistoMCEtaPtNotTriggered;                            //! array of histos with weighted eta in not triggered events, pT
    TH1F**                fHistoMCEtaPtNoVertex;                                //! array of histos with weighted eta in not triggered events, pT    
    TH1F**                fHistoMCEtaWOWeightPt;                                //! array of histos with unweighted eta, pT
    TH1F**                fHistoMCEtaWOEvtWeightPt;                             //! array of histos without event weights eta, pT
    TH1F**                fHistoMCPi0InAccPt;                                   //! array of histos with weighted pi0 in acceptance, pT
    TH1F**                fHistoMCEtaInAccPt;                                   //! array of histos with weighted eta in acceptance, pT
    TH1F**                fHistoMCPi0InAccPtNotTriggered;                       //! array of histos with weighted pi0 in acceptance, pT
    TH1F**                fHistoMCEtaInAccPtNotTriggered;                       //! array of histos with weighted eta in acceptance, pT
    TH1F**                fHistoMCPi0WOEvtWeightInAccPt;                        //! array of histos without evt weight pi0 in acceptance, pT
    TH1F**                fHistoMCEtaWOEvtWeightInAccPt;                        //! array of histos without evt weight eta in acceptance, pT
    TH2F**                fHistoMCPrimaryPtvsSource;                            //! array of histos with weighted primary particles, pT vs source
    TH2F**                fHistoMCSecPi0PtvsSource;                             //! array of histos with secondary pi0, pT, source
    TH1F**                fHistoMCSecPi0Source;                                 //! array of histos with secondary pi0, source
    TH2F**                fHistoMCSecPi0InAccPtvsSource;                        //! array of histos with secondary pi0 in acceptance, pT, source
    TH1F**                fHistoMCSecEtaPt;                                     //! array of histos with secondary eta, pT
    TH1F**                fHistoMCSecEtaSource;                                 //! array of histos with secondary eta, source
    TH2F**                fHistoMCPi0PtJetPt;                                   //! array of histos with weighted pi0, pT, hardest jet pt
    TH2F**                fHistoMCEtaPtJetPt;                                   //! array of histos with weighted eta, pT, hardest jet pt

    // MC validated reconstructed quantities mesons
    TH2F**                fHistoTruePi0InvMassPt;                               //! array of histos with validated mothers, invMass, pt
    TH2F**                fHistoTruePi0InvMassPtAdditional;                     //! array of histos with validated mothers, invMass, pt
    TH2F**                fHistoTrueEtaInvMassPt;                               //! array of histos with validated mothers, invMass, pt
    TH2F**                fHistoTrueEtaInvMassPtAdditional;                     //! array of histos with validated mothers, invMass, pt
    TH2F**                fHistoTruePrimaryPi0InvMassPt;                        //! array of histos with validated weighted primary mothers, invMass, pt
    TH2F**                fHistoTruePrimaryEtaInvMassPt;                        //! array of histos with validated weighted primary mothers, invMass, pt
    TH2F**                fHistoTruePrimaryPi0W0WeightingInvMassPt;             //! array of histos with validated unweighted primary mothers, invMass, pt
    TH2F**                fHistoTruePrimaryEtaW0WeightingInvMassPt;             //! array of histos with validated unweighted primary mothers, invMass, pt
    TProfile2D**          fProfileTruePrimaryPi0WeightsInvMassPt;               //! array of profiles with weights for validated primary mothers, invMass, pt
    TProfile2D**          fProfileTruePrimaryEtaWeightsInvMassPt;               //! array of profiles with weights for validated primary mothers, invMass, pt
    TH2F**                fHistoTrueSecondaryPi0InvMassPt;                      //! array of histos with validated secondary mothers, invMass, pt
    TH2F**                fHistoTrueSecondaryPi0FromK0sInvMassPt;               //! array of histos with validated secondary mothers from K0s, invMass, pt
    TH2F**                fHistoTrueSecondaryPi0FromK0lInvMassPt;               //! array of histos with validated secondary mothers from K0l, invMass, pt
    TH2F**                fHistoTrueSecondaryPi0FromEtaInvMassPt;               //! array of histos with validated secondary mothers from eta, invMass, pt
    TH2F**                fHistoTrueSecondaryPi0FromLambdaInvMassPt;            //! array of histos with validated secondary mothers from Lambda, invMass, pt
    // MC validated reconstructed quantities photons
    TH2F**                fHistoClusPhotonBGPt;                                     //! array of histos with cluster photon BG, pt, source
    TH2F**                fHistoClusPhotonPlusConvBGPt;                             //! array of histos with cluster photon plus conv BG, pt, source
    TH1F**                fHistoTrueClusGammaPt;                                    //! array of histos with validated cluster (electron or photon), pt
    TH2F**                fHistoTrueClusGammaEResE;                                 //! array of histos with photon (validated and conversions) energy resolution ((E_rec-E_true)/E_true) as function of E_rec
    TH2F**                fHistoTrueClusPhotonGammaEResE;                           //! array of histos with validated photon energy resolution ((E_rec-E_true)/E_true) as function of E_rec 
    TH1F**                fHistoTrueNLabelsInClus;                                  //! array of histos with number of labels in cluster
    TH1F**                fHistoTruePrimaryClusGammaPt;                             //! array of histos with validated primary photon cluster, pt
    TH2F**                fHistoTruePrimaryClusGammaESDPtMCPt;                      //! array of histos with validated primary photon cluster, rec Pt, MC pt
    TH1F**                fHistoTruePrimaryClusConvGammaPt;                         //! array of histos with validated primary conv photon cluster, pt
    TH2F**                fHistoTruePrimaryClusConvGammaESDPtMCPt;                  //! array of histos with validated primary conv photon cluster, rec Pt, MC pt
    TH2F**                fHistoTrueSecondaryClusGammaPt;                           //! array of histos with validated secondary photon cluster, pt
    TH2F**                fHistoTrueSecondaryClusConvGammaPt;                       //! array of histos with validated secondary converted photon cluster, pt
    TH2F**                fHistoTrueSecondaryClusGammaMCPt;                         //! array of histos with validated secondary photon cluster, MC pt
    TH2F**                fHistoTrueSecondaryClusConvGammaMCPt;                     //! array of histos with validated secondary converted photon cluster, MC pt
    TH2F**                fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt;        //! array of histos with validated secondary photon cluster from X from K0s, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt;    //! array of histos with validated secondary converted photon cluster from X from K0s, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt;        //! array of histos with validated secondary photon cluster from X from K0l, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt;    //! array of histos with validated secondary converted photon cluster from X from K0l, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt;     //! array of histos with validated secondary photon cluster from X from Lambda, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt; //! array of histos with validated secondary converted photon cluster from X from Lambda, MC pt, rec pt
    TH2F**                fHistoDoubleCountTruePi0InvMassPt;                        //! array of histos with double counted pi0s, invMass, pT
    TH2F**                fHistoDoubleCountTrueEtaInvMassPt;                        //! array of histos with double counted etas, invMass, pT
    TH2F**                fHistoDoubleCountTrueClusterGammaPt;                      //! array of histos with double counted cluster photons
    vector<int>           fVectorDoubleCountTruePi0s;                               //! vector containing labels of validated pi0
    vector<int>           fVectorDoubleCountTrueEtas;                               //! vector containing labels of validated eta
    vector<int>           fVectorDoubleCountTrueClusterGammas;                      //! vector containing labels of validated cluster photons
    TH1F**                fHistoMultipleCountTrueClusterGamma;                      //! array of histos how often TrueClusterGammas are counted
    map<int,int>          fMapMultipleCountTrueClusterGammas;                       //! map containing cluster photon labels that are counted at least twice
    TH2F**                fHistoTruePi0InvMassPtAlpha;                              //! array of histogram with pure pi0 signal inv Mass, energy of cluster
    TH2F**                fHistoTruePi0PureGammaInvMassPtAlpha;                     //! array of histogram with pure pi0 signal (only pure gammas) inv Mass, energy of cluster

    // event histograms
    TH1F**                fHistoNEvents;                                        //! array of histos with event information
    TH1F**                fHistoNEventsWOWeight;                                //! array of histos with event information without event weights
    TH1F**                fHistoNGoodESDTracks;                                 //! array of histos with number of good tracks (2010 Standard track cuts)
    TH1F**                fHistoVertexZ;                                        //! array of histos with vertex z distribution for selected events
    TH1F**                fHistoNGammaCandidates;                               //! array of histos with number of gamma candidates per event
    TH1F**                fHistoNGammaCandidatesBasic;                          //! array of histos with number of gamma candidates per event for basic cluster cut
    TH2F**                fHistoNGoodESDTracksVsNGammaCandidates;               //! array of histos with number of good tracks vs gamma candidates
    TH2F**                fHistoSPDClusterTrackletBackground;                   //! array of histos with SPD tracklets vs SPD clusters for background rejection
    TH1F**                fHistoNV0Tracks;                                      //! array of histos with V0 counts
    TProfile**            fProfileEtaShift;                                     //! array of profiles with eta shift
    TProfile**            fProfileJetJetXSection;                               //! array of profiles with xsection for jetjet
    TH1F**                fHistoJetJetNTrials;                                  //! array of histos with ntrials for jetjet
    TH2F**                fHistoPtHardJJWeight;                                 //! array of histos with ntrials for jetjet
    TH1F**                fHistoPionSpectrum;                                   //! array of histos with charged pion spectrum
    TH1F**                fHistoProtonSpectrum;                                 //! array of histos with proton spectrum
    TH1F**                fHistoKaonSpectrum;                                   //! array of histos with charged kaon spectrum
    TH1F**                fHistoNPionSpectrum;                                  //! array of histos with Neutral pion spectrum
    TH1F**                fHistoEtaSpectrum;                                    //! array of histos with Eta spectrum
    TH1F**                fHistoDMesonSpectrum;                                 //! array of histos with D0 meson spectrum

    // histograms and variables for flow
    TProfile2D** fHistoMotherInvMassPtV0CCos2phi;                               //! array of histogram with signal + BG for same event photon pairs cos(2(phi-Psi)) inv Mass, pt
    TProfile2D** fHistoMotherInvMassPtV0ACos2phi;                               //! array of histogram with signal + BG for same event photon pairs cos(2(phi-Psi)) inv Mass, pt
    TProfile2D** fHistoMotherBackInvMassPtV0CCos2phi;                           //! array of histogram with only BG for same event photon pairs cos(2(phi-Psi)), inv Mass, pt
    TProfile2D** fHistoMotherBackInvMassPtV0ACos2phi;                           //! array of histogram with only BG for same event photon pairs cos(2(phi-Psi)), inv Mass, pt
    int fRunNumber;
    int fOldRunNumber;
    bool IsVZEROCalibOn;                                                        // switch for VZERO qn calib
    bool IsQAVZERO;                                                             // switch for Fill VZERO calib qa
    bool fUseInOutPlane;                                                        // switch to change to in and out of plane method
    TList* fListVZEROCalib;                                                     // read list for V0 Calib
    TFile* fVZEROCalibFile;
    double fPsi2V0C;
    double fPsi2V0A;
    TH2D** fHist2DPsi2V0CCent;
    TH2D** fHist2DPsi2V0ACent;
    TH1D*  hMultV0; // Dobrin
    AliOADBContainer *contMult;
    AliOADBContainer *contQxncm;
    AliOADBContainer *contQyncm;
    AliOADBContainer *contQxnam;
    AliOADBContainer *contQynam;
    // 18q/r
    TH2F* fHCorrectV0ChWeghts;
    // V0C
    TProfile** fProfileV0CQxCentGE;                                           //! array of profiles for QA
    TProfile** fProfileV0CQyCentGE;                                           //! array of profiles for QA
    TProfile** fProfileV0CQxVtxGE;                                            //! array of profiles for QA
    TProfile** fProfileV0CQyVtxGE;                                            //! array of profiles for QA
    TH2D** fHist2CalibPsi2V0CCentGE;                                          //! array of histos for calibration
    TProfile** fProfileV0AQxCentGE;                                           //! array of profiles for QA
    TProfile** fProfileV0AQyCentGE;                                           //! array of profiles for QA
    TProfile** fProfileV0AQxVtxGE;                                            //! array of profiles for QA
    TProfile** fProfileV0AQyVtxGE;                                            //! array of profiles for QA
    TH2D** fHist2CalibPsi2V0ACentGE;                                          //! array of histos for calibration
    TProfile** fProfileV0CQxCentRC;                                           //! array of profiles for QA
    TProfile** fProfileV0CQyCentRC;                                           //! array of profiles for QA
    TProfile** fProfileV0CQxVtxRC;                                            //! array of profiles for QA
    TProfile** fProfileV0CQyVtxRC;                                            //! array of profiles for QA
    TH2D** fHist2CalibPsi2V0CCentRC;                                          //! array of histos for calibration
    TProfile** fProfileV0AQxCentRC;                                           //! array of profiles for QA
    TProfile** fProfileV0AQyCentRC;                                           //! array of profiles for QA
    TProfile** fProfileV0AQxVtxRC;                                            //! array of profiles for QA
    TProfile** fProfileV0AQyVtxRC;                                            //! array of profiles for QA
    TH2D** fHist2CalibPsi2V0ACentRC;                                          //! array of histos for calibration
    TProfile** fHist2V0Res;                                                   //! array of profiles of V0 Resolution

    // additional variables
    float               fCent;                                                // Event centrality
    double              fEventPlaneAngle;                                     // EventPlaneAngle
    TRandom3            fRandom;                                              // random
    int                 fnCuts;                                               // number of cuts to be analysed in parallel
    int                 fiCut;                                                // current cut
    int                 fIsHeavyIon;                                          // switch for pp = 0, PbPb = 1, pPb = 2
    int                 fDoLightOutput;                                       // switch for running light output, 0 -> normal mode, 1 -> light mode, 2 -> minimum
    bool                fDoMesonAnalysis;                                     // flag for meson analysis
    bool                fIsFromDesiredHeader;                                 // flag for MC headers
    bool                fIsOverlappingWithOtherHeader;                        // flag for particles in MC overlapping between headers
    bool                fIsOverlapWithMBHeader;                               // flag for particles overlapping with MB header when MB header is one of the selected headers
    int                 fIsMC;                                                // flag for MC information
    double              fWeightJetJetMC;                                      // weight for Jet-Jet MC
    bool                fDoInOutTimingCluster;                                // manual timing cut for cluster to combine cluster within timing cut and without
    double              fMinTimingCluster;                                    // corresponding ranges, min
    double              fMaxTimingCluster;                                    // corresponding ranges, max
    bool                fEnableSortForClusMC;                                 // switch on sorting for MC labels in cluster
    bool                fDoPrimaryTrackMatching;                              // switch for basic track matching for primaries
    TTree*                tBrokenFiles;                                         // tree for keeping track of broken files
    TObjString*           fFileNameBroken;                                      // string object for broken file name
    TObjString*           fFileNameTrigger;                                     // string object for triggering filename in MB
    TGenPhaseSpace        fGenPhaseSpace;                                       // For generation of decays into two gammas
    TClonesArray*         fAODMCTrackArray;                                     // storage of track array

    int                 fLocalDebugFlag;                                      // debug flag for local running, must be '0' for grid running
    bool                fAllowOverlapHeaders;                                 // enable overlapping headers for cluster selection
    int                 fNCurrentClusterBasic;                                // current number of cluster without minE
    int                 fTrackMatcherRunningMode;                             // CaloTrackMatcher running mode
    bool                fDoPi0Only;                                           // switches ranges of histograms and binnings to pi0 specific analysis

    // for V0 event plane calib
    TH1D* hQx2mV0[2];
    TH1D* hQy2mV0[2];

  private:
    AliAnalysisTaskPi0v2Calo(const AliAnalysisTaskPi0v2Calo&);                  // Prevent copy-construction
    AliAnalysisTaskPi0v2Calo &operator=(const AliAnalysisTaskPi0v2Calo&);       // Prevent assignment

    ClassDef(AliAnalysisTaskPi0v2Calo, 1);
};

#endif
