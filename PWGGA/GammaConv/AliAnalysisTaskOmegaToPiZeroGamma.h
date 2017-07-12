#ifndef ALIANLYSISTASKOMEGATOPIZEROGAMMA_cxx
#define ALIANLYSISTASKOMEGATOPIZEROGAMMA_cxx

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
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include "THnSparse.h"
#include <vector>
#include <map>

class AliAnalysisTaskOmegaToPiZeroGamma : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskOmegaToPiZeroGamma();
    AliAnalysisTaskOmegaToPiZeroGamma(const char *name);
    virtual ~AliAnalysisTaskOmegaToPiZeroGamma();

    virtual void   UserCreateOutputObjects();
    virtual Bool_t Notify();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);
    void InitBack();

    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    void SetIsHeavyIon(Int_t flag){
      fIsHeavyIon = flag;    
    }

    // base functions for selecting photon and meson candidates in reconstructed data
    void ProcessClusters();
    void ProcessPhotonCandidates();
    void CalculateOmegaCandidates();
    
    // MC functions
    void SetIsMC                        ( Int_t isMC)                                       { fIsMC = isMC                              ;}
    void ProcessMCParticles             ();
    void ProcessAODMCParticles          ();
    void ProcessTruePhotonCandidates    ( AliAODConversionPhoton* TruePhotonCandidate);
    void ProcessTruePhotonCandidatesAOD ( AliAODConversionPhoton* TruePhotonCandidate);
    void ProcessTrueClusterCandidates   ( AliAODConversionPhoton* TruePhotonCandidate);
    void ProcessTrueClusterCandidatesAOD( AliAODConversionPhoton* TruePhotonCandidate);
    void RelabelAODPhotonCandidates     ( Bool_t mode);
    void ProcessTrueMesonCandidates     ( AliAODConversionMother *OmegaCandidate,
                                          AliAODConversionPhoton *TrueGammaCandidate0, 
                                          AliAODConversionPhoton *TrueGammaCandidate1, 
                                          AliAODConversionPhoton *TrueGammaCandidate2);
    void ProcessTrueMesonCandidatesAOD  ( AliAODConversionMother *OmegaCandidate,
                                          AliAODConversionPhoton *TrueGammaCandidate0,
                                          AliAODConversionPhoton *TrueGammaCandidate1,
                                          AliAODConversionPhoton *TrueGammaCandidate2);
    
    // switches for additional analysis streams or outputs
    void SetDoMesonQA                   ( Int_t flag )                                      { fDoMesonQA = flag                           ;}
    void SetDoPhotonQA                  ( Int_t flag )                                      { fDoPhotonQA = flag                          ;}
    void SetPlotHistsExtQA              ( Bool_t flag )                                     { fSetPlotHistsExtQA = flag                   ;}

    // Setting the cut lists for the conversion photons
    void SetEventCutList                ( Int_t nCuts, 
                                          TList *CutArray)                                  {
                                                                                              fnCuts = nCuts                              ;
                                                                                              fEventCutArray = CutArray                   ;
                                                                                            }

    // Setting the cut lists for the conversion photons
    void SetConversionCutList           ( Int_t nCuts, 
                                          TList *CutArray)                                  {
                                                                                              fnCuts = nCuts                              ;
                                                                                              fCutArray = CutArray                        ;
                                                                                            }

      // Setting the cut lists for the calo photons
    void SetCaloCutList                 ( Int_t nCuts, 
                                          TList *CutArray)                                  {
                                                                                              fnCuts = nCuts                              ;
                                                                                              fClusterCutArray = CutArray                 ;
                                                                                            }
    
    // Setting cut lists for the neutral pion
    void SetNeutralPionCutList          ( Int_t nCuts,
                                          TList *CutArray)                                  {
                                                                                              fnCuts = nCuts                              ;
                                                                                              fNeutralPionCutArray = CutArray             ;
                                                                                            }

    // Setting cut lists for the (omega) meson
    void SetMesonCutList                ( Int_t nCuts,
                                          TList *CutArray)                                  {
                                                                                              fnCuts = nCuts                              ;
                                                                                              fMesonCutArray = CutArray                   ;
                                                                                            }

    // BG HandlerSettings
    void CalculateBackground            ();
    void SetMoveParticleAccordingToVertex       ( Bool_t flag )                             { fMoveParticleAccordingToVertex = flag       ;}
    void MoveParticleAccordingToVertex          ( AliAODConversionPhoton* particle, 
                                                  const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    void MoveParticleAccordingToVertex          (AliAODConversionMother* particle,
                                                  const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    void UpdateEventByEventData         ();
    
    // Additional functions for convenience
    void SetLogBinningXTH2              ( TH2* histoRebin );
    Int_t GetSourceClassification       (Int_t daughter, 
                                         Int_t pdgCode );
    Bool_t CheckVectorOnly              ( vector<Int_t> &vec, 
                                          Int_t tobechecked );
    Bool_t CheckVectorForDoubleCount    ( vector<Int_t> &vec,
                                          Int_t tobechecked );

    void FillMultipleCountMap           ( map<Int_t,Int_t> &ma,
                                          Int_t tobechecked );
    void FillMultipleCountHistoAndClear ( map<Int_t,
                                          Int_t> &ma,
                                          TH1F* hist );

        // Function to enable MC label sorting
    void SetEnableSortingOfMCClusLabels (Bool_t enableSort) { fEnableSortForClusMC   = enableSort;}

    // set reconstruction method
    void SetReconMethod                 (Int_t ReconMethod) { fReconMethod = ReconMethod;}

    // set flag for doing pi0-gamma angle cut
    void SetDoPiZeroGammaAngleCut         (Bool_t DoPiZeroGammaAngleCut){ fDoPiZeroGammaAngleCut = DoPiZeroGammaAngleCut;}

    // set scaling factors for pi0-gamma angle cut
    void SetlowerFactor                  (Double_t lowerFactor){ flowerFactor = lowerFactor;}
    void SetupperFactor                  (Double_t upperFactor){ fupperFactor = upperFactor;}
    
  protected:
    AliV0ReaderV1*                      fV0Reader;              // basic photon Selection Task
    TString                             fV0ReaderName;
    AliGammaConversionAODBGHandler**    fBGHandler;             // BG handler for Conversion 
    AliGammaConversionAODBGHandler**    fBGClusHandler;         // BG handler for Cluster
    AliGammaConversionAODBGHandler**    fBGPi0Handler;          // BG handler for pi0's
    AliVEvent*                          fInputEvent;            // current event
    AliMCEvent*                         fMCEvent;               // corresponding MC event
    TList**                             fCutFolder;             // Array of lists for containers belonging to cut
    TList**                             fESDList;               // Array of lists with histograms with reconstructed properties
    TList**                             fTrueList;              // Array of lists with histograms with MC validated reconstructed properties
    TList**                             fMCList;                // Array of lists with histograms with pure MC information
    TList**                             fClusterOutputList;     //!Array of lists of output histograms for cluster photons
    TList*                              fOutputContainer;       // Output container
    TClonesArray*                       fReaderGammas;          // Array with conversion photons selected by V0Reader Cut
    TList*                              fGammaCandidates;       // current list of photon candidates
    TList*                              fClusterCandidates;     //! current list of cluster candidates
    TList*                              fPi0Candidates;         //! current list of pi0 candidates
    TList*                              fEventCutArray;         // List with Event Cuts
    AliConvEventCuts*                   fEventCuts;             // EventCutObject
    TList*                              fCutArray;              // List with Conversion Cuts
    AliConversionPhotonCuts*            fConversionCuts;        // ConversionCutObject
    TList*                              fClusterCutArray;       // List with Cluster Cuts
    AliCaloPhotonCuts*                  fCaloPhotonCuts;        // CaloPhotonCutObject
    TList*                              fNeutralPionCutArray;   // List with neutral pion cuts
    TList*                              fMesonCutArray;         // List with meson cuts
    TF1*                                fmaxfit;                // function describing location of max. points in the distribution of pi0-gamma angle vs. pT
    Double_t                            flowerFactor;           // factor maxfit is multiplied by to get lower limit for pi0-gamma angle cut
    Double_t                            fupperFactor;           // factor maxfit is multiplied by to get upper limit for pi0-gamma angle cut
    Double_t                            fMinPi0Pt;              // Min Pi0 Pt cut in GeV
    
    //histograms for Conversions reconstructed quantities
    TH1F**                  fHistoConvGammaPt;                  //! histogram conversion photon pT
    TH1F**                  fHistoConvGammaR;                   //! histogram conversion photon R
    TH1F**                  fHistoConvGammaEta;                 //! histogram conversion photon Eta
                      // 0: garbage,
                      // 1: background
                      // 2: secondary photon not from eta or k0s,
                      // 3: secondary photon from eta, 
                      // 4: secondary photon from k0s, 
                      // 5: dalitz
                      // 6: primary gamma

    //histograms for mesons reconstructed quantities
    TH2F**                  fHistoPhotonPairInvMassPt;              //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
    TH2F**                  fHistoPhotonPairMatchedInvMassPt;       //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
    TH2F**                  fHistoMotherMatchedInvMassPt;           //! array of histograms with invariant mass and Pt of rejected omega candidates due to
                                                                    //! matching between a pcm track and cluster which do not belong to the pair forming pi0
    TH2F**                  fHistoMotherAngleCutRejectedInvMassPt;  //! array of histograms with invariant mass and Pt of omega candidates rejected by pi0-gamma angle cut
    TH2F**                  fHistoPhotonPairYPt;                //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, Y
    TH2F**                  fHistoPhotonPairAlphaPt;            //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, alpha
    TH2F**                  fHistoPhotonPairOpenAnglePt;        //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, openAngle
    TH2F**                  fHistoPhotonPairEtaPhi;             //! array of histograms with Eta, Phi of pi0 candidates
    TH2F**                  fHistoMotherConvPhotonEtaPhi;       //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17 ,eta/phi of conversion photon
    TH2F**                  fHistoMotherInvMassPt;              //! array of histograms for invariant mass of omega candidates
    TH2F**                  fHistoMotherYPt;                    //! array of histograms for Y of omega candidates
    TH2F**                  fHistoMotherAlphaPt;                //! array of histograms with pT, Alpha of omega candidates
    TH2F**                  fHistoMotherEtaPhi;                 //! array of histograms with Eta and Phi of omega candidates
    TH2F**                  fHistoMotherPi0AnglePt;             //! array of histograms with angle between omega candidates and pi0 candidates
    TH2F**                  fHistoMotherGammaAnglePt;           //! array of histograms with angle between omega candidates and combined gammas
    TH2F**                  fHistoPi0GammaAnglePt;              //! array of histograms with angle between pi0 candidates and combined gammas
    TH1F**                  fHistoGammaFromMotherPt;            //! array of histograms with pT of gammas from omega candidates

    // BG histograms
    TH2F**                  fHistoSamePi0DiffGammaBackInvMassPt;//! array of histograms with background generated by combining pi0's from the current event and photons from the handler
    TH2F**                  fHistoDiffPi0SameGammaBackInvMassPt;//! array of histograms with background generated by combining pi0's from the handler and photons from the current event
    
    // histograms for rec photon clusters
    TH1F**                  fHistoClusGammaPt;                  //! array of histos with cluster, pt
    TH1F**                  fHistoClusOverlapHeadersGammaPt;    //! array of histos with cluster, pt overlapping with other headers
                    
    //histograms for pure MC quantities
    TH1F**                  fHistoMCAllGammaPt;                 //! array of histos with all gamma, pT
    TH1F**                  fHistoMCAllGammaEMCALAccPt;         //! array of histos with all gamma in EMCAL acceptance, pT
    TH1F**                  fHistoMCConvGammaPt;                //! array of histos with converted gamma, pT
    TH1F**                  fHistoMCConvGammaR;                 //! array of histos with converted gamma, R
    TH1F**                  fHistoMCConvGammaEta;               //! array of histos with converted gamma, Eta
    TH1F**                  fHistoMCPi0Pt;                      //! array of histos with weighted pi0, pT
    TH1F**                  fHistoMCPi0WOWeightPt;              //! array of histos with unweighted pi0, pT
    TH1F**                  fHistoMCPi0WOEvtWeightPt;           //! array of histos without event weights pi0, pT
    TH1F**                  fHistoMCPi0InAccPt;                 //! array of histos with weighted pi0 in acceptance, pT
    TH1F**                  fHistoMCPi0WOWeightInAccPt;         //! array of histos without weight pi0 in acceptance, pT
    TH2F**                  fHistoMCPi0PtY;                     //! array of histos with weighted pi0, pT, Y
    TH2F**                  fHistoMCPi0PtAlpha;                 //! array of histos with weighted pi0, pT, alpha
    TH2F**                  fHistoMCPi0PtJetPt;                 //! array of histos with weighted pi0, pT, hardest jet pt
    TH1F**                  fHistoMCGammaFromAllOmegaPt;        //! array of histos with pT of photons which are decay products of omegas
    TH1F**                  fHistoMCGammaFromOmegaInAccPt;      //! array of histos with pT of photons from omegas in acceptance
    TH2F**                  fHistoMCPi0FromAllOmegaInvMassPt;   //! array of histos with pT & Inv Mass of pi0s which are daughters of omegas
    TH2F**                  fHistoMCOmegaInAccInvMassPt;        //! array of histos with true omegas in acceptance, Inv Mass, pT
    TH2F**                  fHistoMCOmegaInvMassPt;             //! array of histos with pT of true omegas within rapidity window, for acceptance correction
    TH2F**                  fHistoMCAllOmegaYPt;                //! array of histos with pT, Y of all true omegas
    TH2F**                  fHistoMCOmegaInAccYPt;              //! array of histos with pT, Y of true omegas in acceptance
    TH2F**                  fHistoMCAllOmegaAlphaPt;            //! array of histos with pT, alpha of omegas which decayed into pi0+gamma
    TH2F**                  fHistoMCOmegaInAccAlphaPt;          //! array of histos with pT, alpha of omegas in acceptance which decayed into pi0+gamma
    TH2F**                  fHistoMCPi0FromAllOmegaAlphaPt;     //! array of histos with pT, alpha of pi0s which are daughters of omegas
    TH2F**                  fHistoMCPi0FromOmegaInAccAlphaPt;   //! array of histos with pT, alpha of pi0s which are daughters of omegas in acceptance
    TH2F**                  fHistoMCPi0FromAllOmegaYPt;         //! array of histos with pT, Y of pi0s from all true omegas
    TH2F**                  fHistoMCPi0FromOmegaInAccYPt;       //! array of histos with pT, Y of pi0s from omegas in acceptance
    TH2F**                  fHistoMCPi0FromOmegaInAccInvMassPt; //! array of histos with Inv Mass, pT of pi0s from omegas in acceptance
    TH2F**                  fHistoMCPi0FromAllOmegaEtaPhi;      //! array of histos with eta and phi of pi0s from all true omegas
    TH2F**                  fHistoMCPi0FromOmegaInAccEtaPhi;    //! array of histos with eta and phi of pi0s from omegas in acceptance
    TH2F**                  fHistoMCAllOmegaEtaPhi;             //! array of histos with eta and phi of all true omegas
    TH2F**                  fHistoMCOmegaInAccEtaPhi;           //! array of histos with eta and phi of all true omegas in acceptance
    TH2F**                  fHistoMCAllOmegaPiZeroAnglePt;      //! array of histos with angle between true omega and pi0 daughter and pT of omega
    TH2F**                  fHistoMCAllPiZeroGammaAnglePt;      //! array of histos with angle between pi0 and gamma daughters of true omega and pT of omega
    TH2F**                  fHistoMCAllOmegaGammaAnglePt;       //! array of histos with angle between true omega and gamma daughter and pT of omega
    TH2F**                  fHistoMCInAccOmegaPiZeroAnglePt;    //! array of histos with angle between true omega in acc. and pi0 daughter
    TH2F**                  fHistoMCInAccPiZeroGammaAnglePt;    //! array of histos with angle between pi0 and gamma daughters of true omega in acc.
    TH2F**                  fHistoMCInAccOmegaGammaAnglePt;     //! array of histos with angle between true omega in acc. and gamma daughter
    TH2F**                  fHistoMCAllOmegaInvMassPt;          //! array of histos with all true omegas, invMass, pT
    TH2F**                  fHistoMCAllOmegaPtPi0Pt;            //! array of histos with pT of all omegas and that of their pi0 daughters
    TH2F**                  fHistoMCInAccOmegaPtPi0Pt;          //! array of histos with pT of omegas in acceptance and that of their pi0 daughters
    TH2F**                  fHistoMCAllOmegaPtGammaPt;          //! array of histos with pT of all omegas and those of their gamma daughters
    TH2F**                  fHistoMCInAccOmegaPtGammaPt;        //! array of histos with pT of omegas in acceptance and those of their gamma daughters

    // MC validated reconstructed quantities mesons
    TH2F**                  fHistoTrueOmegaInvMassPt;                           //! array of histos with reconstructed true omegas, invMass, pT
    TH2F**                  fHistoTrueOmegaYPt;                                 //! array of histos with pT, Y of reconstructed true omegas
    TH2F**                  fHistoTrueOmegaAlphaPt;                             //! array of histos with pT, Alpha of reconstructed true omegas
    TH2F**                  fHistoTruePi0FromOmegaYPt;                          //! array of histos with pT, Y of pi0s from reconstructed true omegas
    TH2F**                  fHistoTruePi0FromOmegaInvMassPt;                    //! array of histos with reconstructed true pi0s which are decay products of true omegas, invMass, pT
    TH2F**                  fHistoTruePi0FromOmegaAlphaPt;                      //! array of histos with pT, alpha of pi0s from reconstructed true omegas
    TH2F**                  fHistoTruePi0FromOmegaEtaPhi;                       //! array of histos with eta, phi of pi0s from reconstructed true omegas
    TH2F**                  fHistoTruePi0FromOmegaOpenAnglePt;                  //! array of histos with pT, opening angle of pi0s fro reconstructed true omegas
    TH2F**                  fHistoTrueOmegaPi0AnglePt;                          //! array of histos with pT, angle between reconstructed true omegas and pi0 daughters
    TH2F**                  fHistoTrueOmegaGammaAnglePt;                        //! array of histos with pT, angle between reconstructed true omegas and gamma daughters
    TH2F**                  fHistoTruePi0GammaAnglePt;                          //! array of histos with pT, angle between pi0 and gamma from reconstructed true omegas
    TH2F**                  fHistoTrueOmegaEtaPhi;                              //! array of histos with eta, phi of reconstructed true omegas
    TH2F**                  fHistoTrueOmegaPtPi0Pt;                             //! array of histos with pT of validated omegas against pT of pi0s from validated omegas
    TH2F**                  fHistoTrueOmegaPtGammaPt;                           //! array of histos with pT of validated omegas against pT of gammas from validated omegas
    TH1F**                  fHistoTrueGammaFromOmegaPt;                         //! array of histos with pT of photons from validated omegas

    vector<Int_t>           fVectorRecTruePi0s;                                 //! array of strings containing the stack position of the reconstructed validated pi0
    vector<Int_t>           fVectorDoubleCountTruePi0s;                         //! vector containing labels of validated pi0
    TH1F**                  fHistoMultipleCountTruePi0;                         //! array of histos how often TruePi0s are counted
    map<Int_t,Int_t>        fMapMultipleCountTruePi0s;                          //! map containing pi0 labels that are counted at least twice
    // event histograms
    TH1F**                  fHistoNEvents;                                      //! array of histos with event information
    TH1F**                  fHistoNEventsMinGamma;                              //! array of histos with no. of events containing the minimum number of EMCal/PCM photons for each reconstruction method
    TH1F**                  fHistoMCOmegaDecayChannels;                         //! array of histos with no. of true omegas in each decay channel
    TH1F**                  fHistoNEventsWOWeight;                              //! array of histos with event information without event weights
    TH1F**                  fHistoNGoodESDTracks;                               //! array of histos with number of good tracks (2010 Standard track cuts)
    TH1F**                  fHistoVertexZ;                                      //! array of histos with vertex z distribution for selected events
    TH1F**                  fHistoNGammaCandidates;                             //! array of histos with number of gamma candidates per event
    TH1F**                  fHistoNClusterCandidates;                           //! array of histos with number of cluster candidates per event
    TH2F**                  fHistoNGoodESDTracksVsNGammaCandidates;             //! array of histos with number of good tracks vs gamma candidates
    TH2F**                  fHistoSPDClusterTrackletBackground;                 //! array of histos with SPD tracklets vs SPD clusters for background rejection
    TH1F**                  fHistoNV0Tracks;                                    //! array of histos with V0 counts
    TProfile**              fProfileEtaShift;                                   //! array of profiles with eta shift
    TProfile**              fProfileJetJetXSection;                             //! array of profiles with xsection for jetjet
    TH1F**                  fHistoJetJetNTrials;                                //! array of histos with ntrials for jetjet

    // additional variables
    Double_t                fEventPlaneAngle;                                   // EventPlaneAngle
    Int_t                   fNGammaCandidates;                                  // number of gamma candidates in event
    Double_t*               fUnsmearedPx;                                       //[fNGammaCandidates]
    Double_t*               fUnsmearedPy;                                       //[fNGammaCandidates]
    Double_t*               fUnsmearedPz;                                       //[fNGammaCandidates]
    Double_t*               fUnsmearedE;                                        //[fNGammaCandidates]
    Int_t*                  fMCEventPos;                                        //[fNGammaCandidates]
    Int_t*                  fMCEventNeg;                                        //[fNGammaCandidates]
    Int_t*                  fESDArrayPos;                                       //[fNGammaCandidates]
    Int_t*                  fESDArrayNeg;                                       //[fNGammaCandidates]
    Int_t                   fnCuts;                                             // number of cuts to be analysed in parallel
    Int_t                   fiCut;                                              // current cut  
    Bool_t                  fMoveParticleAccordingToVertex;                     // boolean for BG calculation
    Int_t                   fIsHeavyIon;                                        // switch for pp = 0, PbPb = 1, pPb = 2
    Int_t                   fDoMesonQA;                                         // flag for meson QA
    Int_t                   fDoPhotonQA;                                        // flag for photon QA
    Bool_t                  fIsFromMBHeader;                                    // flag for MC headers
    Bool_t                  fIsOverlappingWithOtherHeader;                      // flag for particles in MC overlapping between headers
    Int_t                   fIsMC;                                              // flag for MC information
    Bool_t                  fSetPlotHistsExtQA;                                 // flag for extended QA hists
    Double_t                fWeightJetJetMC;                                    // weight for Jet-Jet MC 
    Bool_t                  fEnableSortForClusMC;                               // switch on sorting for MC labels in cluster
    Int_t                   fReconMethod;                                       // switch for combining photons: PCM-cal,cal = 0; PCM-cal,PCM = 1; cal-cal,cal = 2;
                                                                                // cal-cal,PCM = 3; PCM-PCM,cal = 4; PCM-PCM,PCM = 5
    Bool_t                  fDoPiZeroGammaAngleCut;                             // flag for pi0-gamma angle cut

  private:
    AliAnalysisTaskOmegaToPiZeroGamma(const AliAnalysisTaskOmegaToPiZeroGamma&); // Prevent copy-construction
    AliAnalysisTaskOmegaToPiZeroGamma &operator=(const AliAnalysisTaskOmegaToPiZeroGamma&); // Prevent assignment

    ClassDef(AliAnalysisTaskOmegaToPiZeroGamma, 11);
};

#endif
