
#ifndef AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson_H
#define AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliAnalysisTaskSE.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliPrimaryPionSelector.h"
#include "AliConversionMesonCuts.h"
#include "AliConvEventCuts.h"
#include "AliCaloPhotonCuts.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliAnalysisTaskJetOutlierRemoval.h"
#include "TProfile2D.h"
#include "TH3F.h"
#include "TF1.h"
#include "TArrayI.h"
#include <vector>

class AliESDInputHandler;
class AliMCEventHandler;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpidCuts;
class AliTriggerAnalysis;

class AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson: public AliAnalysisTaskSE
{
  public:

    AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson();
    AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson( const char* name );
    virtual ~AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson();

    virtual void UserExec(Option_t *);
    virtual void UserCreateOutputObjects();
    virtual Bool_t Notify();
    virtual void Terminate(const Option_t *);

    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}

    void SetMoveParticleAccordingToVertex(Bool_t flag){fMoveParticleAccordingToVertex = flag;}
    void SetIsHeavyIon(Int_t flag){
      if (flag == 1 || flag ==2 ){
        fIsHeavyIon = 1;
      } else {
        fIsHeavyIon = 0;
      }
    }

    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void SetLightOutput(Int_t flag){fDoLightOutput = flag;}
    void SetUnsmearedOutputs(TString unsmearingoutputs){
      fEnableNoCorrOutput           = (unsmearingoutputs.Contains("0") || unsmearingoutputs.Contains("4")) ? kTRUE : kFALSE;
      fEnableSubNDMOutput           = unsmearingoutputs.Contains("1") ? kTRUE : kFALSE;
      fEnableFixedpzOutput          = unsmearingoutputs.Contains("2") ? kTRUE : kFALSE;
      fEnableSubLambdaOutput        = unsmearingoutputs.Contains("3") ? kTRUE : kFALSE;
      fEnablePCMEMCUnsmearing       = unsmearingoutputs.Contains("4") ? kTRUE : kFALSE;
      }
    void SetEventCutList(Int_t nCuts, TList *CutArray){
      fnCuts= nCuts;
      fEventCutArray = CutArray;
    }
    void SetConversionCutList(TList *CutArray){ fGammaCutArray = CutArray;}
    void SetClusterCutList(TList *CutArray){ fClusterCutArray = CutArray;}
    void SetPionCutList(TList *CutArray){ fPionCutArray = CutArray;}
    void SetNeutralPionCutList(TList *CutArray){ fNeutralDecayMesonCutArray = CutArray; }
    void SetMesonCutList(TList *CutArray){ fMesonCutArray = CutArray; }
    void SetDoMesonQA(Int_t flag){ fDoMesonQA = flag; }
    void SetNDMRecoMode(Int_t mode){fNDMRecoMode = mode; }
    void SetTolerance(Double_t tol){fTolerance=tol;}
    void SetSelectedHeavyNeutralMeson(Int_t selectMeson){fSelectedHeavyNeutralMeson=selectMeson;}
    void SetTrackMatcherRunningMode(Int_t mode){fTrackMatcherRunningMode = mode;}
    void SetAllowOverlapHeaders( Bool_t allowOverlapHeader ) {fAllowOverlapHeaders = allowOverlapHeader;}
    void SetPionSelectorName( TString flag ) {fPionSelectorName = flag;}
    TString GetPionSelectorName() {return fPionSelectorName;}
    // Function to enable MC label sorting
    void SetEnableSortingOfMCClusLabels(Bool_t enableSort)  { fEnableSortForClusMC   = enableSort; }
    void SetDoMaterialBudgetWeightingOfGammasForTrueMesons(Bool_t flag) { fDoMaterialBudgetWeightingOfGammasForTrueMesons = flag;}
    void SetDoMaterialBudgetWeightingOfGammasForInvMassHistogram(Bool_t flag) { fUseMatBudWeightsForInvMassHistogram = flag;}
    void SetBckgReductionTree(Int_t cutoff){ 
      fEnableBckgReductionStudy = kTRUE;
      fMLtreeCutOff = cutoff; }

    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}

  private:

    void InitBack();

    // routines for photon selection from conversions
    void ProcessConversionPhotonCandidates();
    void ProcessTrueConversionPhotonCandidates(AliAODConversionPhoton*);
    void ProcessTrueConversionPhotonCandidatesAOD(AliAODConversionPhoton*);

    // routines for photon selection from clusters
    void ProcessCaloPhotonCandidates();
    void ProcessTrueCaloPhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate);
    void ProcessTrueCaloPhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate);

    void ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionMother *TrueNeutralPionCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate, Double_t weightMatBudget);
    void ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionMother *TrueNeutralPionCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate, Double_t weightMatBudget);
    void MoveParticleAccordingToVertex(AliAODConversionMother* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);

    void FixPzToMatchPDGInvMassNDM(AliAODConversionMother* particle);
    void FixPzVecToMatchPDGInvMass(TLorentzVector* track);

    // routines for neutral pion candidates from pure conversion
    void ProcessNeutralDecayMesonCandidatesPureConversions();
    void ProcessTrueNeutralPionCandidatesPureConversions(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
    void ProcessTrueNeutralPionCandidatesPureConversionsAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);

    // routines for neutral pion candidates from pure calo
    void ProcessNeutralPionCandidatesPureCalo();
    void ProcessTrueNeutralPionCandidatesPureCalo(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
    void ProcessTrueNeutralPionCandidatesPureCaloAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);

    // routines for neutral pion candidates from mixed conv + calo
    void ProcessNeutralPionCandidatesMixedConvCalo();
    void ProcessTrueNeutralPionCandidatesMixedConvCalo( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
    void ProcessTrueNeutralPionCandidatesMixedConvCaloAOD( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);

    void ProcessPionCandidates();
    void ProcessPionCandidatesAOD();
    void ProcessMCParticles();
    void ProcessAODMCParticles();
    void CalculateMesonCandidates(AliAODConversionPhoton *vParticle);
    void CalculateBackground(Int_t mode);
    void UpdateEventByEventData();

    Bool_t IsPiPlPiMiPiZeroDecay(AliMCParticle *fMCMother) const;
    Bool_t IsPiPlPiMiPiZeroDecayAOD( TClonesArray* trackArray, AliAODMCParticle *fMCMother) const;
    Bool_t IsPiPlPiMiEtaDecay(AliMCParticle *fMCMother) const;
    Bool_t IsPiPlPiMiEtaDecayAOD( TClonesArray* trackArray, AliAODMCParticle *fMCMother) const;
    Bool_t IsEtaPrimePiPlPiMiEtaDaughter( Int_t label ) const;
    Bool_t IsEtaPrimePiPlPiMiEtaDaughterAOD(TClonesArray* trackArray, Int_t label ) const;
    Bool_t IsEtaPiPlPiMiPiZeroDaughter( Int_t label ) const;
    Bool_t IsEtaPiPlPiMiPiZeroDaughterAOD(TClonesArray* trackArray, Int_t label ) const;
    Bool_t IsOmegaPiPlPiMiPiZeroDaughter( Int_t label ) const;
    Bool_t IsOmegaPiPlPiMiPiZeroDaughterAOD(TClonesArray* trackArray, Int_t label ) const;
    Bool_t IsD0PiPlPiMiPiZeroDaughter( Int_t label ) const;
    Bool_t IsD0PiPlPiMiPiZeroDaughterAOD(TClonesArray* trackArray, Int_t label ) const;
    Bool_t IsRhoDaughterAOD(TClonesArray* trackArray, Int_t label ) const;
    Bool_t IsKaonDaughterAOD(TClonesArray* trackArray, Int_t label ) const;
    AliAODConversionPhoton* ReturnPiPlPiMiOneFromHNMMassClosestToRho(vector<Int_t> lGoodPosPionIndexPrev, vector<Int_t> lGoodNegPionIndexPrev, UInt_t i, UInt_t j);

    Bool_t GammaIsNeutralMesonPiPlPiMiNDMDaughter( Int_t label ) const;
    Bool_t GammaIsNeutralMesonPiPlPiMiNDMDaughterAOD(TClonesArray* trackArray, Int_t label ) const;

    Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);

    Bool_t KinematicCut(AliAODConversionMother *negpion, AliAODConversionMother *pospion, AliAODConversionMother *neutpion, AliAODConversionMother *omega);
    void RelabelAODPhotonCandidates(Bool_t mode);
    AliExternalTrackParam* GetConstrainedParameterAOD(const AliAODTrack* aodTr, const AliAODVertex* vtx, double bz);
    Double32_t CalculateP2(Double_t xyz[3],Double_t pxpypz[3]);
    Double_t GetAlphaLFromLorentz(TLorentzVector Particle1, TLorentzVector Particle2);
    Bool_t MesonIsSelectedByAlphaCut(Double_t alphaValue, Double_t MesonPt, Int_t ParametrizationMode);


    AliV0ReaderV1*                    fV0Reader;                                          //!<! V0Reader for basic conversion photon selection
    TString                           fV0ReaderName;                                      ///< Name of the V0 reader
    TString                           fCorrTaskSetting;                                   ///< Name of the correction framework setting used (if any)
    AliPrimaryPionSelector*           fPionSelector;                                      //!<! primary charged pion selector, basic selection of pi+,pi-
    TString                           fPionSelectorName;                                  ///< Name of the PionSelector
    AliGammaConversionAODBGHandler**  fBGHandlerPiPl;                                     //!<! BG handler Pos Pion
    AliGammaConversionAODBGHandler**  fBGHandlerPiMi;                                     //!<! BG handler Neg Pion
    AliVEvent*                        fInputEvent;                                        //!<! current event
    AliMCEvent*                       fMCEvent;                                           //!<! current MC event
    TList**                           fCutFolder;                                         //!<! list of output folders with main cut name
    TList**                           fESDList;                                           //!<! list with main output histograms for data
    TList**                           fTrueList;                                          //!<! list with validated reconstructed MC histograms
    TList**                           fTrueTreeList;                                      //!<! list containing TTree's for MC True
    TList**                           fMCList;                                            //!<! list with pure MC histograms
    TList*                            fOutputContainer;                                   //!<! output container
    TClonesArray*                     fReaderGammas;                                      //!<! array with photon from fV0Reader
    vector<Int_t>                     fSelectorNegPionIndex;                              //!<! array with pion indices for negative pions from fPionSelector
    vector<Int_t>                     fSelectorPosPionIndex;                              //!<! array with pion indices for positive pions from fPionSelector
    TList*                            fGoodConvGammas;                                    //!<! good conv gammas after selection
    TList*                            fClusterCandidates;                                 //!<! good calo gammas after selection
    TList*                            fNeutralDecayParticleCandidates;                    //!<! good neutral pion candidates
    Bool_t                            fUseMatBudWeightsForInvMassHistogram;               //!<! use material budget weights for Invariant Mass Histograms, needed for rec. Eff
    vector<Float_t>                   fNeutralDecayParticleCandidateMatBudWeights;        //!<! material budget weights for good neutral pion candidates
    TList*                            fNeutralDecayParticleSidebandCandidates;            //!<! good neutral pion candidates from sideband
    TList*                            fNeutralDecayParticleSwappCandidates;               //!<! good neutral pion candidates from Swapp method
    TList*                            fPosPionCandidates;                                 //!<! good positive pion candidates
    TList*                            fNegPionCandidates;                                 //!<! good negative pion candidates
    TList*                            fEventCutArray;                                     ///< array with event cuts
    TList*                            fGammaCutArray;                                     ///< array with Conversion Cuts
    TList*                            fClusterCutArray;                                   ///< array with Cluster Cuts
    TList*                            fPionCutArray;                                      ///< array with charged pion cuts
    TList*                            fNeutralDecayMesonCutArray;                         ///< array with neutral pion cuts
    TList*                            fMesonCutArray;                                     ///< array with neutral meson cuts
    AliConvEventCuts*                 fEventCuts;                                         //!<! current event cuts
    AliConversionPhotonCuts*          fConversionCuts;                                    //!<! current conversion cuts
    AliCaloPhotonCuts*                fClusterCuts;                                       //!<! current cluster cuts
    AliAnalysisTaskJetOutlierRemoval* fOutlierJetReader;                                  // JetReader

    // TTrees                           //!<!
    /**  Tree containing information about an event where eta->pi+pi-pi0 was found if ID isn't covered by current implementations */
    TTree**                           fTreeEventInfoHNM;                                  //!<!
    TTree**                           fTreeTrueNDMFromHNM;                                //!<!
    Short_t                           fCasePiPi;                                          ///< 0:PiPlPiMi 1:PiMiPiZero 1:PiPlPiZero
    Float_t                           fSamePiPiMotherID;                                  ///< ID of mother of two pions
    Float_t                           fSamePiPiMotherInvMass;                             ///< Invariant mass of mother of two pions
    Float_t                           fSamePiPiMotherPt;                                  ///< pT of mother of two pions
    Float_t                           fSamePiPiPiMotherID;                                ///< ID of mother of two pions
    Float_t                           fSamePiPiPiMotherInvMass;                           ///< Invariant mass of mother of two pions
    Float_t                           fSamePiPiPiMotherPt;                                ///< pT of mother of two pions
    Float_t                           fV0MultiplicityHNMEvent;                            ///< V0 multiplicity of an event where a true Eta was found
    Float_t                           fTrackMultiplicityHNMEvent;                         ///< track multiplicity of an event where a true Eta was found
    Float_t                           fZVertexHNMEvent;                                   ///< z position of primary vertex of an event where a true Eta was found
    Float_t                           fPtHNM;                                             ///< pT of a true Eta
    Float_t                           fPtNDM;
    Float_t                           fInvMassNDM;
    Float_t                           fPDGMassNDM;                                        ///< PDG mass of either pi0 or eta
    Float_t                           fNDMMinPtPossible;                                  ///< min pt of NDM measurable for each method
    Float_t                           fPDGMassChargedPion;                                ///< PDG mass of either pi0 or eta
    Int_t                             fPDGCodeNDM;                                        ///< PDG code of either pi0 or eta
    Int_t                             fPDGCodeAnalyzedMeson;                              ///< PDG code of the analyzed heavy netural meson
    Bool_t                            fEnableNoCorrOutput;                                ///< Output the minv-pT histogram without any unsmearing correction
    Bool_t                            fEnableSubNDMOutput;                                ///< Output the minv-pT histogram with the pi0 mass error subtracted for each omega meson
    Bool_t                            fEnableFixedpzOutput;                               ///< Output the minv-pT histogram after correcting the NDM pz to match PDG mass
    Bool_t                            fEnableSubLambdaOutput;                             ///< Output the minv-pT histogram with Lambda*DeltaMpi0 subtracted for each omega meson
    TF1**                             fLambda;                                            ///< Function returning factor which multiplied with DeltaMpi0 will be subracted in the unsmearing
    Bool_t                            fEnablePCMEMCUnsmearing;                            ///< Output the minv-pT histogram with an analytic unsmearing assuming perfect PCM energy resolution
    Bool_t                            fEnableNDMEfficiency;                               ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableNDMInputSpectrum;                            ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableBasicMesonQA;                                ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableBackgroundQA;                                ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnable3DHistoQA;                                   ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableCorrelationTreeQA;                           ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableBackgroundCalculation;                       ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableTreeTrueNDMFromHNM;                          ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableTrueMotherPiPlPiMiNDMAdditionalInvMassPt;    ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableTrueMotherPiPlPiMiNDMInvMassPtBackground;    ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableAsymmetryPlotCombCPionVsNPion;               ///< Turn On or Off if Histograms are created and used
    Bool_t                            fEnableAsymmetryPlot_NotAccepted;                   ///< Turn On or Off if Histograms are created and used
    Bool_t                            enableDalitzAllPt;                                  ///< Turn On or Off if Histograms are created and used
    Bool_t                            enableDalitzLowPt;                                  ///< Turn On or Off if Histograms are created and used
    Bool_t                            enableDalitzMidPt;                                  ///< Turn On or Off if Histograms are created and used
    Bool_t                            enableDalitzHighPt;                                 ///< Turn On or Off if Histograms are created and used
    Double_t                          HistoDalitzPtRangeMin_LowPt;                        ///< Min Range of Dalitz Plots for LowPt
    Double_t                          HistoDalitzPtRangeMax_LowPt;                        ///< Max Range of Dalitz Plots for LowPt
    Double_t                          HistoDalitzPtRangeMin_MidPt;                        ///< Min Range of Dalitz Plots for MidPt
    Double_t                          HistoDalitzPtRangeMax_MidPt;                        ///< Max Range of Dalitz Plots for MidPt
    Double_t                          HistoDalitzPtRangeMin_HighPt;                       ///< Min Range of Dalitz Plots for HighPt
    Double_t                          HistoDalitzPtRangeMax_HighPt;                       ///< Max Range of Dalitz Plots for HighPt
    
    // reconstructed particles
    TH1F**                            fHistoConvGammaPt;                                  //!<! array of histos of conversion photon, pt
    TH1F**                            fHistoConvGammaEta;                                 //!<! array of histos of conversion photon, eta
    TH1F**                            fHistoClusterGammaPt;                               //!<! array of histos of Cluster photon, pt
    TH1F**                            fHistoClusterGammaEta;                              //!<! array of histos of Cluster photon, eta
    TH1F**                            fHistoClusterGammaE;                                //!<! array of histos of Cluster photon, energy
    TH1F**                            fHistoNegPionPt;                                    //!<! array of histos of negative pion, pt
    TH1F**                            fHistoPosPionPt;                                    //!<! array of histos of positive pion, pt
    TH2F**                            fHistoPionPionInvMassPt;                            //!<! array of histos of pion pion, invMass, pT_{pi+pi-}
    TH2F**                            fHistoGammaGammaInvMassPt;                          //!<! array of histos of gamma-gamma, invMass, pT_{gamma gamma}
    TH2F**                            fHistoGammaGammaInvMassPtBeforeCuts;                //!<! array of histos of gamma-gamma, invMass, pT_{gamma gamma}
    TH2F**                            fHistoSwappingGammaGammaInvMassPt;                  //!<! array of histos of gamma-gamma, invMass, pT_{gamma gamma}
    TH2F**                            fHistoMotherInvMassPt;                              //!<! array of histos of pi+pi-pi0 same event, invMass, pT_{pi+pi-pi0}
    TH2F**                            fHistoMotherInvMassPtRejectedKinematic;             //!<! array of histos of rejected pi+pi-pi0 same event, invMass, pT_{pi+pi-pi0}
   
    // general QA histograms 
    TH1I**                            fHistoNumberClusterGamma;                           //!<! array of histos of number of Cluster photons per event
    TH1F**                            fHistoNegPionPhi;                                   //!<! array of histos of negative pion, phi
    TH1F**                            fHistoPosPionPhi;                                   //!<! array of histos of positive pion, phi
    TH1F**                            fHistoNegPionEta;                                   //!<! array of histos of negative pion, eta
    TH1F**                            fHistoPosPionEta;                                   //!<! array of histos of positive pion, eta
    TH2F**                            fHistoNegPionClsTPC;                                //!<! array of histos of negative pion, findable tpc cluster, pT
    TH2F**                            fHistoPosPionClsTPC;                                //!<! array of histos of positive pion, findable tpc cluster, pT
    TH2F**                            fHistoPionDCAxy;                                    //!<! array of histos of pion, dca_xy, pT
    TH2F**                            fHistoPionDCAz;                                     //!<! array of histos of pion, dca_z, pT
    TH2F**                            fHistoPionDCAxyFromOmega;                           //!<! array of histos of pions from omegas, dca vs pT
    TH2F**                            fHistoPionDCAzFromOmega;                            //!<! array of histos of pions from omegas, dca vs pT
    TH2F**                            fHistoPionDCAxyFromRho;                             //!<! array of histos of pions from rho mesons, dca vs pT
    TH2F**                            fHistoPionDCAzFromRho;                              //!<! array of histos of pions from rho mesons, dca vs pT
    TH2F**                            fHistoPionDCAxyFromKaon;                            //!<! array of histos of pions from kaons, dca vs pT
    TH2F**                            fHistoPionDCAzFromKaon;                             //!<! array of histos of pions from kaons, dca vs pT
    TH2F**                            fHistoPionTPCdEdxNSigma;                            //!<! array of histos of pion, p, TPC nSigma dEdx pion
    TH2F**                            fHistoPionTPCdEdx;                                  //!<! array of histos of pion, p, TPC dEdx
    
    // filled only for Dalitz QA flag   
    //Asymmetry Plot
    TH2F**                            fHistoAsymmetryPlotCombCPionVsNPion;                //!<!
    TH2F**                            fHistoAsymmetryPlotCombCPionVsNPion_NotAccepted;    //!<!
    //All Dalitz
    TH2F**                            fHistoDalitzPlotPosFixedPzNDM;                     //!<!
    TH2F**                            fHistoDalitzPlotNegFixedPzNDM;                     //!<!
    TH2F**                            fHistoDalitzPlotPosSubNDM;                         //!<!
    TH2F**                            fHistoDalitzPlotNegSubNDM;                         //!<!
    //Dalitz Low Pt
    TH2F**                            fHistoDalitzPlotPosFixedPzNDM_LowPt;                //!<!
    TH2F**                            fHistoDalitzPlotNegFixedPzNDM_LowPt;                //!<!
    TH2F**                            fHistoDalitzPlotPosSubNDM_LowPt;                    //!<!
    TH2F**                            fHistoDalitzPlotNegSubNDM_LowPt;                    //!<!
    //Dalitz Mid Pt
    TH2F**                            fHistoDalitzPlotPosFixedPzNDM_MidPt;                //!<!
    TH2F**                            fHistoDalitzPlotNegFixedPzNDM_MidPt;                //!<!
    TH2F**                            fHistoDalitzPlotPosSubNDM_MidPt;                    //!<!
    TH2F**                            fHistoDalitzPlotNegSubNDM_MidPt;                    //!<!
    //Dalitz High Pt
    TH2F**                            fHistoDalitzPlotPosFixedPzNDM_HighPt;               //!<!
    TH2F**                            fHistoDalitzPlotNegFixedPzNDM_HighPt;               //!<!
    TH2F**                            fHistoDalitzPlotPosSubNDM_HighPt;                   //!<!
    TH2F**                            fHistoDalitzPlotNegSubNDM_HighPt;                   //!<!

    TH2F**                            fHistoBackInvMassPt      ;                          //!<! Event mixing background group 1 (pi+ and pi- from same event)
    TH2F**                            fHistoMotherLikeSignBackInvMassPt;                  //!<! array of histos of pi+pi+pi0 likesign mixed event, invMass, pT_{pi+pi-pi0}

    // angle distributions
    TH2F**                            fHistoAngleHNMesonPiPlPiMi;                         //!<! angle between combined Pi+ and Pi- and omega
    TH2F**                            fHistoAngleHNMesonNDM;                              //!<! angle between Pi0 and omega
    TH2F**                            fHistoAngleHNMesonPiPl;                             //!<! angle between Pi+ and omega
    TH2F**                            fHistoAngleHNMesonPiMi;                             //!<! angle between Pi- and omega
    TH2F**                            fHistoAnglePiPlPiMi;                                //!<! angle between Pi+ and Pi-
    TH2F**                            fHistoAngleNDMPiMi;                                 //!<! angle between Pi0 and Pi-
    TH2F**                            fHistoAnglePiPlNDM;                                 //!<! angle between Pi+ and Pi0
    TH2F**                            fHistoAngleSum;                                     //!<! angle between omega and Pi0 + angle between Pi+ and Pi- + angle between Pi0 and Pi-
    TH2F**                            fHistoTrueAngleSum;
    TH2F**                            fHistoTrueHNMesonPtvsNDMPt;

    TH2F**                            fHistoMotherInvMassSubNDM;                          //!<! invariant mass of (pi+,pi-,pi0) - invariant mass of pi0
    TH2F**                            fHistoBackInvMassPtSubNDM;                          //!<! background group 1, invMass-invMass(pi0), pT_{pi+pi-pi0} (pi+ and pi- from same event)
    TH2F**                            fHistoMotherLikeSignBackInvMassSubNDMPt;            //!<! array of histos of pi+pi+pi0 likesign mixed event, invMass-invMass(pi0), pT_{pi+pi-pi0}

    TH2F**                            fHistoMotherInvMassFixedPzNDM;                      // invariant mass of (pi+,pi-,pi0) - invariant mass of pi0
    /** background group 1 mixed event, invMass, pT_{pi+pi-pi0}, the Pz of the pi0 was fixed such that its invMass matches the PDG value */
    TH2F**                            fHistoBackInvMassPtFixedPzNDM;                      //!<!
    /** array of histos of pi+pi+pi0 likesign mixed event, invMass, pT_{pi+pi+pi0}, the Pz of the pi0 was fixed such that
     *  its invMass matches the PDG value
     */
    TH2F**                            fHistoMotherLikeSignBackInvMassFixedPzNDMPt;        //!<!
    TH2F**                            fHistoMotherInvMassSubLambda;                       //!<! invariant mass of (pi+,pi-,pi0) - lambda* invariant mass error of pi0
    TH2F**                            fHistoBackInvMassPtSubLambda;                       //!<! background group 1, invMass-invMass(pi0), pT_{pi+pi-pi0} (pi+ and pi- from same event)
    TH2F**                            fHistoMotherLikeSignBackInvMassSubLambdaPt;         //!<! array of histos of pi+pi+pi0 likesign mixed event, invMass-lambda*invMass(pi0), pT_{pi+pi-pi0}
    TH1F**                            fHistoPCMEMCScalingFactor;                          //!<! Scaling factor for cluster photon in PCMEMC unsmearing
    
    // pure MC properties
    TH1F**                            fHistoMCAllGammaPt;                                 //!<! array of histos of all produced gammas in the specified y range
    TH1F**                            fHistoMCConvGammaPt;                                //!<! array of histos of all converted gammas in the specified y range
    TH1F**                            fHistoMCGammaFromNeutralMesonPt;                    //!<! array of histos of all produced gammas from omega or eta via pi+pi-pi0 in the specified y range/
    TH1F**                            fHistoMCPosPionsFromNeutralMesonPt;                 //!<! array of histos of all produced positive pions from omega or eta via pi+pi-pi0 in the specified y range/
    TH1F**                            fHistoMCNegPionsFromNeutralMesonPt;                 //!<! array of histos of all produced negative pions from omega or eta via pi+pi-pi0 in the specified y range/
    TH1F**                            fHistoMCHNMPiPlPiMiNDMPt;                           //!<! array of histos of produced NNM via pi+pi-pi0 in the specified y range
    TH1F**                            fHistoMCHNMPiPlPiMiNDMEta;                          //!<! array of histos of produced HNM via pi+pi-pi0 in the specified y range
    TH1F**                            fHistoMCHNMPiPlPiMiNDMPhi;                          //!<! array of histos of produced HNM via pi+pi-pi0 in the specified y range
    TH1F**                            fHistoMCAllPosPionsPt;                              //!<! array of histos with all produced primary positive pions in the specified y range
    TH1F**                            fHistoMCAllNegPionsPt;                              //!<! array of histos with all produced primary negative pions in the specified y range
    /** array of histos of produced etas via pi+pi-pi0 in the specified y range, with decay products in respective y, eta ranges */
    TH1F**                            fHistoMCHNMPiPlPiMiNDMInAccPt;                      //!<!
    TH1F**                            fHistoMCNDMFromHNMInputPt;                          //!<! array of histos of produced Pi0 coming from
    TH1F**                            fHistoMCNDMFromHNMInputInAccPt;                     //!<!
    TH2F**                            fHistoMCHNMInAccVsNDMPt;                            //!<!
    // reconstructed particles MC validated
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPt;                 //!<! histos with reconstructed validated eta or omega, inv mass, pT
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM;           //!<! histos with reconstructed validated eta or omega, inv mass, pT fixed pi0 mass
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM;       //!<! histos with reconstructed validated eta or omega, inv mass, pT fixed pi0 mass
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPtSubLambda;        //!<! histos with reconstructed validated eta or omega, inv mass, pT fixed pi0 mass
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMAdditionalInvMassPtSubNDM; //!<! histos with reconstructed validated eta or omega, inv mass, pT fixed pi0 mass, only additionally found omegas
    
    TH2F**                          fHistoDoubleCountTruePi0InvMassPt;                    //!<! array of histos with double counted pi0s, invMass, pT
    TH2F**                          fHistoDoubleCountTrueHNMInvMassPt;                    //!<! array of histos with double counted etas, invMass, pT
    TH2F**                          fHistoDoubleCountTrueConvGammaRPt;                    //!<! array of histos with double counted photons, R, pT

    TH2F**                          fHistoTrueMotherGammaGammaInvMassPt;                  //!<! histos with reconstructed validated pi0, inv mass, pT
    TH2F**                          fHistoTrueMotherGammaGammaFromHNMInvMassPt;           //!<! histos with reconstructed validated pi0, inv mass, pT
    TH1F**                          fHistoTrueConvGammaPt;                                //!<! histos with reconstructed validated conv gamma, pT
    TH1F**                          fHistoTrueConvGammaFromNeutralMesonPt;                //!<! histos with reconstructed validated conv gamma from eta or omega via pi0, pT
    TH1F**                          fHistoTrueClusterGammaPt;                             //!<! histos with reconstructed validated cluster gamma, pT
    TH1F**                          fHistoTrueClusterGammaFromNeutralMesonPt;             //!<! histos with reconstructed validated cluster gamma from eta or omega via pi0, pT
    TH1F**                          fHistoTruePosPionPt;                                  //!<! histos with reconstructed validated positive pion, pT
    TH1F**                          fHistoTruePosPionFromNeutralMesonPt;                  //!<! histos with reconstructed validated positive pion from eta or omega, pT
    TH1F**                          fHistoTrueNegPionPt;                                  //!<! histos with reconstructed validated negative pion, pT
    TH1F**                          fHistoTrueNegPionFromNeutralMesonPt;                  //!<! histos with reconstructed validated negative pion from eta or omega, pT

    // reconstructed particles MC validated different mesons
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent;   //!<! histos with all reconstructed validated mesons which are not analyzed (eta or omega), inv mass, pT
    
    TH2F**                          fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt;         //!<! histos with all reconstructed validated pi+pi-NDM that are combinatorical
    TH2F**                          fHistoTruePiPlPiMiNDMContaminationInvMassPt;          //!<! histos with reconstructed pi0 that are not actually pions

    
    // general meson QA (MC)
    TH1F**                            fHistoMCAllMesonPt;                                 //!<! array of histos of pt of all neutral decay mesons in the specified y range
    TH1F**                            fHistoMCAllMesonEta;                                //!<! array of histos of eta of all neutral decay mesons in the specified y range
    TH1F**                            fHistoMCAllMesonPhi;                                //!<! array of histos of phi of all neutral decay mesons in the specified y range
    TH1F**                            fHistoMCMesonFromNeutralMesonPt;                    //!<! array of histos of pt of neutral decay mesons from heavy meson in the specified y range
    TH1F**                            fHistoMCMesonFromNeutralMesonEta;                   //!<! array of histos of eta of neutral decay mesons from heavy meson in the specified y range
    TH1F**                            fHistoMCMesonFromNeutralMesonPhi;                   ///<! array of histos of phi of neutral decay mesons from heavy meson in the specified y range
    TH1F**                            fHistoMCAllPosPionsEta;                             //!<! array of histos with all produced primary positive pions in the specified y range
    TH1F**                            fHistoMCAllPosPionsPhi;                             //!<! array of histos with all produced primary positive pions in the specified y range
    TH1F**                            fHistoMCAllNegPionsEta;                             //!<! array of histos with all produced primary negative pions in the specified y range
    TH1F**                            fHistoMCAllNegPionsPhi;                             //!<! array of histos with all produced primary negative pions in the specified y range
    TH1F**                            fHistoMCPosPionsFromNeutralMesonEta;                //!<! array of histos of all produced positive pions from omega or eta via pi+pi-pi0 in the specified y range/
    TH1F**                            fHistoMCPosPionsFromNeutralMesonPhi;                //!<! array of histos of all produced positive pions from omega or eta via pi+pi-pi0 in the specified y range/
    TH1F**                            fHistoMCNegPionsFromNeutralMesonEta;                //!<! array of histos of all produced negative pions from omega or eta via pi+pi-pi0 in the specified y range/
    TH1F**                            fHistoMCNegPionsFromNeutralMesonPhi;                //!<! array of histos of all produced negative pions from omega or eta via pi+pi-pi0 in the specified y range/
    TH2F**                            fHistoMCHNMPiPlPiMiNDMEtavsPt;                      //!<! array of histos of produced HNM eta vs Pt
     // MC truth properties for heavy meson (and decay products)
    TH1F**                            fHistoMCHeavyAllPt;                                 //!<! array of histos with pt of all heavy mesons
    TH1F**                            fHistoMCHeavyAllEta;                                //!<! array of histos with eta of all heavy mesons
    TH1F**                            fHistoMCHeavyAllPhi;                                //!<! array of histos with phi of all heavy mesons
    TH1F**                            fHistoMCHeavyChannelPt;                             //!<! array of histos with pt of heavy mesons in the decay channel
    TH1F**                            fHistoMCHeavyChannelEta;                            //!<! array of histos with eta of heavy mesons in the decay channel
    TH1F**                            fHistoMCHeavyChannelPhi;                            //!<! array of histos with phi of heavy mesons in the decay channel
    TH1F**                            fHistMCChannelNDMFromHeavyPt;                       //!<! array if histos with pt of the ndm of the heavy meson
    TH1F**                            fHistMCChannelNDMFromHeavyEta;                      //!<! array if histos with eta of the ndm of the heavy meson
    TH1F**                            fHistMCChannelNDMFromHeavyPhi;                      //!<! array if histos with phi of the ndm of the heavy meson
    TH1F**                            fHistMCChannelPiPlusFromHeavyPt;                    //!<! array if histos with pt of the piplus of the heavy meson
    TH1F**                            fHistMCChannelPiPlusFromHeavyEta;                   //!<! array if histos with eta of the piplus of the heavy meson
    TH1F**                            fHistMCChannelPiPlusFromHeavyPhi;                   //!<! array if histos with phi of the piplus of the heavy meson
    TH1F**                            fHistMCChannelPiMinusFromHeavyPt;                   //!<! array if histos with pt of the piminus of the heavy meson
    TH1F**                            fHistMCChannelPiMinusFromHeavyEta;                  //!<! array if histos with eta of the piminus of the heavy meson
    TH1F**                            fHistMCChannelPiPMinusFromHeavyPhi;                 //!<! array if histos with phi of the piminus of the heavy meson
    TH2F**                            fHistMCChannelNDMPtHeavyPt;                         //!<! array of histos of pt correlation ndm - heavy meson
    TH2F**                            fHistMCChannelPiPlusPtHeavyPt;                      //!<! array of histos of pt correlation piplus - heavy meson
    TH2F**                            fHistMCChannelPiMinusPtHeavyPt;                     //!<! array of histos of pt correlation piminus - heavy meson
    TH1F**                            fHistoMCHeavyReconstructiblePt;                     //!<! array of histos with pt of reconstructible heavy mesons
    TH1F**                            fHistoMCHeavyReconstructibleEta;                    //!<! array of histos with eta of reconstructible heavy mesons
    TH1F**                            fHistoMCHeavyReconstructiblePhi;                    //!<! array of histos with phi of reconstructible heavy mesons
    TH1F**                            fHistMCReconstructibleNDMFromHeavyPt;               //!<! array if histos with pt of the ndm of the reconstructible heavy meson
    TH1F**                            fHistMCReconstructibleNDMFromHeavyEta;              //!<! array if histos with eta of the ndm of the reconstructible heavy meson
    TH1F**                            fHistMCReconstructibleNDMFromHeavyPhi;              //!<! array if histos with phi of the ndm of the reconstructible heavy meson
    TH1F**                            fHistMCReconstructiblePiPlusFromHeavyPt;            //!<! array if histos with pt of the piplus of the reconstructible heavy meson
    TH1F**                            fHistMCReconstructiblePiPlusFromHeavyEta;           //!<! array if histos with eta of the piplus of the reconstructible heavy meson
    TH1F**                            fHistMCReconstructiblePiPlusFromHeavyPhi;           //!<! array if histos with phi of the piplus of the reconstructible heavy meson
    TH1F**                            fHistMCReconstructiblePiMinusFromHeavyPt;           //!<! array if histos with pt of the piminus of the reconstructible heavy meson
    TH1F**                            fHistMCReconstructiblePiMinusFromHeavyEta;          //!<! array if histos with eta of the piminus of the reconstructible heavy meson
    TH1F**                            fHistMCReconstructiblePiPMinusFromHeavyPhi;         //!<! array if histos with phi of the piminus of the reconstructible heavy meson
    TH2F**                            fHistMCReconstructibleNDMPtHeavyPt;                 //!<! array of histos of pt correlation ndm - reconstructible heavy meson
    TH2F**                            fHistMCReconstructiblePiPlusPtHeavyPt;              //!<! array of histos of pt correlation piplus - reconstructible heavy meson
    TH2F**                            fHistMCReconstructiblePiMinusPtHeavyPt;             //!<! array of histos of pt correlation piminus - reconstructible heavy meson
    TH1F**                            fHistoTrueMesonFlags;                                 //!<! histo for event counting
    TH2F**                            fHistoTruePionPionInvMassPt;                          //!<! histos with reconstructed validated two pion, invariant mass, pT
    TH2F**                            fHistoTruePionPionFromSameMotherInvMassPt;            //!<! histos with reconstructed validated two pion from same mother, invariant mass, pT
    TH2F**                            fHistoTruePionPionFromHNMInvMassPt;                   //!<! histos with reconstructed validated two pion from eta , invariant mass, pT
    TH2F**                            fHistoTruePionFromHNMInvMassClosestToRhoPt;           //!<! histos filled for each true HNM with the pi+pi- mass closest to rho with one of the pions from the HNM
    TH2F**                            fHistoTruePionFromHNMInvMassPt;                       //!<! histos with reconstructed validated one pion from eta , invariant mass, pT
    // truth for virtual candidate histos
    TH1F**                            fHistoTruevParticleChi2PerNDF;                        //!<!
    TH1F**                            fHistoTruevParticleFromSameMotherChi2PerNDF;          //!<!
    TH1F**                            fHistoTruevParticleFromHNMChi2PerNDF;                 //!<!

    TH1F**                            fHistoTruevParticledS;                                //!<!
    TH1F**                            fHistoTruevParticleFromSameMotherdS;                  //!<!
    TH1F**                            fHistoTruevParticleFromHNMdS;                         //!<!

    // background QA histos (MC)// reconstructed particles MC validated different mesons
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega;    //!<! histos with reconstructed validated eta or omega mesons from  which are not analyzed, inv mass, pT
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho;         //!<! histos with reconstructed validated rho mesons from  which are not analyzed, inv mass, pT
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s;         //!<! histos with reconstructed validated K0s mesons from  which are not analyzed, inv mass, pT
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l;         //!<! histos with reconstructed validated K0l mesons from  which are not analyzed, inv mass, pT
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime;    //!<! histos with reconstructed validated EtaPrime mesons from  which are not analyzed, inv mass, pT
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther;       //!<! histos with reconstructed validated EtaPrime mesons from  which are not analyzed, inv mass, pT
    
    TH2F**                          fHistoTruePionPionArmenteros;                         //!<! histos with reconstructed validated two pion, armenteros (alpha vs qt)
    TH2F**                          fHistoTruePionPionFromRhoArmenteros;                  //!<! histos with reconstructed validated two pions from a rho meson, armenteros (alpha vs qt)
    TH2F**                          fHistoTruePionFromHNMArmenteros;                      //!<! histos with reconstructed validated one pion from HNM, armenteros (alpha vs qt)

    TH2F**                          fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt;         //!<! histos with reconstructed validated pi+ pi-  from omega, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt;       //!<! histos with reconstructed validated pi+ pi-  from eta, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt;         //!<! histos with reconstructed validated pi+ pi-  from rho0, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt;    //!<! histos with reconstructed validated pi+ pi-  from etaprime, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt;         //!<! histos with reconstructed validated pi+ pi-  from K0s, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt;         //!<! histos with reconstructed validated pi+ pi-  from K0s, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt;         //!<! histos with reconstructed validated pi+ pi-  from K0s, invariant mass, pT

    TH2F**                          fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt;       //!<! histos with reconstructed validated pi0 pi-  from omega, invariant mass, pT
    TH2F**                          fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt;     //!<! histos with reconstructed validated pi0 pi-  from eta, invariant mass, pT
    TH2F**                          fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt;       //!<! histos with reconstructed validated pi0 pi-  from rho0, invariant mass, pT
    TH2F**                          fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt;       //!<! histos with reconstructed validated pi0 pi-  from rho0, invariant mass, pT
    TH2F**                          fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt;    //!<! histos with reconstructed validated pi0 pi-  from rho0, invariant mass, pT

    TH2F**                          fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt;       //!<! histos with reconstructed validated pi0 pi+  from omega, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt;     //!<! histos with reconstructed validated pi0 pi+  from eta, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt;       //!<! histos with reconstructed validated pi0 pi+  from rho0, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt;       //!<! histos with reconstructed validated pi0 pi+  from K0l, invariant mass, pT
    TH2F**                          fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt;       //!<! histos with reconstructed validated pi0 pi+  from K0l, invariant mass, pT

    TH2F**                          fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt;     //!<! histos with reconstructed validated pi+pi-NDM that are pure combinatorical (do not share a mother)

    TH2F**                          fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt;     //!<! histos with reconstructed pi+ that are not actually pions
    TH2F**                          fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt;    //!<! histos with reconstructed pi- that are not actually pions
    TH2F**                          fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt;    //!<! histos with reconstructed pi+pi-pi0 that are not actually pions
    TH2F**                          fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt; //!<! histos with reconstructed pi+pi-pi0 that are not actually pions, Crosscheck
    TH2F**                          fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt; //!<! histos with reconstructed pi+pi-pi0 that are not actually pions, Crosscheck
    // --- 

    // Dalitz QA (MC)
    //Asymmetry Plot
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMAsymmetryPlotCombCPionVsNPion;   //!<!
    //Dalitz All Pt
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM;   //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM;   //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM;       //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM;       //!<!
    //Dalitz Low Pt
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt;   //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt;   //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt;       //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt;       //!<!
    //Dalitz Mid Pt
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt;   //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt;   //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt;       //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt;       //!<!
    //Dalitz High Pt
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt;  //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt;  //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt;      //!<!
    TH2F**                          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt;      //!<!

    // 3D histo
    TH3F**                          fHistopi0vsmesonmassshiftangle;                             //!<!

    vector<Int_t>                   fVectorDoubleCountTruePi0s;                           //!<! vector containing labels of validated pi0
    vector<Int_t>                   fVectorDoubleCountTrueHNMs;                           //!<! vector containing labels of validated eta
    vector<Int_t>                   fVectorDoubleCountTrueConvGammas;                     //!<! vector containing labels of validated photons
    // Event properties
    TH1F**                          fHistoNEvents;                                        //!<! histo for event counting
    TH1F**                          fHistoNEventsWOWeight;                                //!<! histo for event counting without weight
    TProfile**                      fProfileJetJetXSection;                               //!<! histo for cross section for jet-jet Monte-Carlo
    TH1F**                          fHistoJetJetNTrials;                                  //!<! histo for number of trials for jet-jet Monte-Carlo
    TH1I**                          fHistoNGoodESDTracks;                                 //!<! histo number of reconstructed primary tracks
    TProfile**                      fProfileEtaShift;                                     //!<! profile for eta shift bookkeeping
    TH2F**                          fHistoSPDClusterTrackletBackground;                   //!<! array of histos with SPD tracklets vs SPD clusters for background rejection

    // not QA histo, but only in not light output
    // virtual candidate histos
    TH1F**                          fHistovParticleChi2PerNDF;                            //!<! 
    TH1F**                          fHistovParticleChi2PerNDFBothConstrained;             //!<!
    TH1F**                          fHistovParticleChi2PerNDFOneConstrained;              //!<!
    TH1F**                          fHistovParticledS;                                    //!<!
    TH1F**                          fHistovParticledSBothConstrained;                     //!<!
    TH1F**                          fHistovParticledSOneConstrained;                      //!<!

    // 2D histo
    TRandom3                        fRandom;                                              ///< random number
    Int_t                           fnCuts;                                               ///< number of cuts to be run in parallel
    Int_t                           fiCut;                                                ///< current cut
    Int_t                           fNumberOfESDTracks;                                   ///< integer with number of primary tracks in this event
    Bool_t                          fMoveParticleAccordingToVertex;                       ///< Flag to move parice to the vertex
    Int_t                           fIsHeavyIon;                                          ///< Flag for collision system 0: pp, 1: PbPb, 2: pPb
    Bool_t                          fDoMesonAnalysis;                                     ///< Flag for switching on meson analysis
    Int_t                           fDoMesonQA;                                           ///< Switching for meson QA 0: no QA, 1: general meson QA, 2: background QA, 3: 3D histogram, 4: Dalitz plots, 5: trees, 23: enable background calculations; 
                                                                                          //    combinations: 6: 1+2, 7: 1+2+3, 8: 1+2+3+5, 9: 2+3, 10: 2+3+5, 11: 1+2+3+4+5
    Bool_t                          fIsFromMBHeader;                                      ///< Flag for particle whether it belongs to accepted header
    Bool_t                          fIsFromDesiredHeader;                                 ///< flag for MC headers
    Bool_t                          fIsOverlappingWithOtherHeader;                        ///< flag for particles in MC overlapping between headers
    Bool_t                          fAllowOverlapHeaders;                                 ///< enable overlapping headers for cluster selection
    Int_t                           fIsMC;                                                ///< Flag for MC
    Int_t                           fSelectedHeavyNeutralMeson;                           ///< Flag for running eta prime
    Int_t                           fDoLightOutput;                                       ///< Flag to turn on light output
    Int_t                           fNDMRecoMode;                                         ///< Flag how neutral pion is reconstructed 0=PCM-PCM, 1=PCM-Calo, 2=Calo-Calo
    Double_t                        fTolerance;                                           ///< tolerance in rad for angle cuts
    Double_t                        fWeightJetJetMC;                                      //!<! Weight for hte jet-jet Monte-Carlo
    Int_t                           fTrackMatcherRunningMode;                             /// CaloTrackMatcher running mode
    Bool_t                          fEnableSortForClusMC;                                 /// switch on sorting for MC labels in cluster
    Bool_t                          fDoMaterialBudgetWeightingOfGammasForTrueMesons;      /// flag to set the usage of material budget weights
    Bool_t                          fDoProfileMaterialBudgetWeights;                      /// Number of MaterialBudgetWeight Bins
    TProfile**                      fProfileMaterialBudgetWeights;                        //!<! histo to track MaterialBudgetWeights
    Int_t                           fNumberOfMaterialBudgetBins;                          /// Number of MaterialBudgetWeight Bins
    Bool_t                          fEnableBckgReductionStudy;                            /// enable the tree for ML based background reduction study
    Int_t                           fMLtreeCutOff;                                        /// the cut-off for reduction in ML background studies
    TH1F*                           fHistoBckReduction;                                   /// the validation histogram for ML background studies

    TArrayI                         fMCEventPos;                                          //!<! Pos. in MC event pos. leg of the photon (for relabelling)
    TArrayI                         fMCEventNeg;                                          //!<! Pos. in MC event neg. leg of the photon (for relabelling)
    TArrayI                         fESDArrayPos;                                         //!<! Pos. in MC AOD array pos. leg of the photon (for relabelling)
    TArrayI                         fESDArrayNeg;                                         //!<! Pos. in MC AOD array pos. leg of the photon (for relabelling)

    /** Tree containing information for the ML based background reduction */
    TTree*                            fTreeBckgReduction;                                 //!<!
    // charged pions branches
    Short_t                           fBuffer_PiPl_px;
    Short_t                           fBuffer_PiPl_py;
    Short_t                           fBuffer_PiPl_pz;
    UShort_t                          fBuffer_PiPl_E;
    Char_t                            fBuffer_PiPl_charge;
    UShort_t                          fBuffer_PiPl_DCAR;
    Short_t                           fBuffer_PiPl_DCAz;
    UShort_t                          fBuffer_PiPl_TPCClus;
    Short_t                           fBuffer_PiPl_dEdxSigma;
    Short_t                           fBuffer_PiPl_TOFdEdxSigma;
    Bool_t                            fBuffer_PiPl_trueID;

    Short_t                           fBuffer_PiMi_px;
    Short_t                           fBuffer_PiMi_py;
    Short_t                           fBuffer_PiMi_pz;
    UShort_t                          fBuffer_PiMi_E;
    Char_t                            fBuffer_PiMi_charge;
    UShort_t                          fBuffer_PiMi_DCAR;
    Short_t                           fBuffer_PiMi_DCAz;
    Short_t                           fBuffer_PiMi_TPCClus;
    Short_t                           fBuffer_PiMi_dEdxSigma;
    Short_t                           fBuffer_PiMi_TOFdEdxSigma;
    Bool_t                            fBuffer_PiMi_trueID;

    Char_t                            fBuffer_PionPair_trueMotherID;        // 0 - not from the same mother, 1 - from the reconstructed NDM, 2 - from NDM but different one <- the plan, for now just the true meson value (1 for neutral pion, 2 for eta, 0 for not true meson)

    // gamma branches
    Short_t                           fBuffer_Gamma1_px;
    Short_t                           fBuffer_Gamma1_py;
    Short_t                           fBuffer_Gamma1_pz;
    UShort_t                          fBuffer_Gamma1_E;
    Short_t                           fBuffer_Gamma1_eta;
    UShort_t                          fBuffer_Gamma1_phi;
    Bool_t                            fBuffer_Gamma1_trueID;

    Short_t                           fBuffer_Gamma2_px;
    Short_t                           fBuffer_Gamma2_py;
    Short_t                           fBuffer_Gamma2_pz;
    UShort_t                          fBuffer_Gamma2_E;
    Short_t                           fBuffer_Gamma2_eta;
    UShort_t                          fBuffer_Gamma2_phi;
    Bool_t                            fBuffer_Gamma2_trueID;

    // ----- PCM standard cuts
    Short_t                           fBuffer_Gamma1_eMomentum;
    Short_t                           fBuffer_Gamma1_eTPCClus;     
    Short_t                           fBuffer_Gamma1_edEdxSigma;  
    Short_t                           fBuffer_Gamma1_epidEdxSigma;
    Short_t                           fBuffer_Gamma1_eTOFPID;     

    Short_t                           fBuffer_Gamma1_pMomentum;
    Short_t                           fBuffer_Gamma1_pTPCClus;     
    Short_t                           fBuffer_Gamma1_pdEdxSigma;  
    Short_t                           fBuffer_Gamma1_ppidEdxSigma;
    Short_t                           fBuffer_Gamma1_pTOFPID;  

    Short_t                           fBuffer_Gamma1_R;
    Short_t                           fBuffer_Gamma1_ArmenterosQt;
    Short_t                           fBuffer_Gamma1_ArmenterosAlpha; 
    Short_t                           fBuffer_Gamma1_chiSquared;
    Short_t                           fBuffer_Gamma1_PsiPair;

    // ***

    Short_t                           fBuffer_Gamma2_eMomentum;
    Short_t                           fBuffer_Gamma2_eTPCClus;     
    Short_t                           fBuffer_Gamma2_edEdxSigma;  
    Short_t                           fBuffer_Gamma2_epidEdxSigma;
    Short_t                           fBuffer_Gamma2_eTOFPID;     

    Short_t                           fBuffer_Gamma2_pMomentum;
    Short_t                           fBuffer_Gamma2_pTPCClus;     
    Short_t                           fBuffer_Gamma2_pdEdxSigma;  
    Short_t                           fBuffer_Gamma2_ppidEdxSigma;
    Short_t                           fBuffer_Gamma2_pTOFPID;  

    Short_t                           fBuffer_Gamma2_R;
    Short_t                           fBuffer_Gamma2_ArmenterosQt;
    Short_t                           fBuffer_Gamma2_ArmenterosAlpha; 
    Short_t                           fBuffer_Gamma2_chiSquared;
    Short_t                           fBuffer_Gamma2_PsiPair;

    // -------- EMC cuts
    Short_t                           fBuffer_Gamma1_M02;
    Short_t                           fBuffer_Gamma2_M02;

    // ---- common 
    Short_t                           fBuffer_GammaPair_OpeningAngle;
    Short_t                           fBuffer_GammaPair_Alpha;
    Short_t                           fBuffer_GammaPair_invMassRec;
    Bool_t                            fBuffer_GammaPair_trueMotherID;         

    // neutral meson 
    Short_t                           fBuffer_NDM_px;
    Short_t                           fBuffer_NDM_py;
    Short_t                           fBuffer_NDM_pz;
    Short_t                           fBuffer_NDM_E;
    Short_t                           fBuffer_NDM_invMassRec;
    Bool_t                            fBuffer_NDM_trueID;




private:
    AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson( const AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson& ); // Not implemented
    AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson& operator=( const AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson& ); // Not implemented

  ClassDef(AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson, 35);
};

#endif // AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson_H
