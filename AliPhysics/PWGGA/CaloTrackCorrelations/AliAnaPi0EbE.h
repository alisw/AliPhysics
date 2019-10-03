#ifndef ALIANAPI0EBE_H
#define ALIANAPI0EBE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaPi0EbE
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Select cluster pairs or single merged clusters with pi0 or eta invariant mass.

/// Class for the analysis of pi0 and eta event by event.
/// Pi0/Eta identified by one of the following:
///  * Invariant mass of 2 cluster in calorimeter
///  * Shower shape analysis in calorimeter
///  * Invariant mass of one cluster in calorimeter and one photon reconstructed in TPC (to be revised)
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaPi0EbE).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///_________________________________________________________________________

// --- ROOT system ---
class TList ;
class TObjString;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaPi0EbE : public AliAnaCaloTrackCorrBaseClass {
  
public:
  
  AliAnaPi0EbE() ; // default ctor
    
  /// Virtual destructor.
  virtual ~AliAnaPi0EbE() { ; }
  
  TObjString *   GetAnalysisCuts();
  
  TList      *   GetCreateOutputObjects();
  
  Int_t          GetMCIndex(Int_t aodTag);
  
  void           Init();
  
  void           InitParameters();
  
  void           MakeAnalysisFillAOD()  ;
  
  void           MakeAnalysisFillHistograms() ;
  
  void           Print(const Option_t * opt) const;
  
  // Main
  
  void           FillEMCALBCHistograms(Float_t energy, Float_t eta, Float_t phi, Float_t time);
  
  void           FillPileUpHistograms(Float_t pt, Float_t time, AliVCluster * c) ;
  
  void           FillRejectedClusterHistograms(Int_t mctag, Int_t nMaxima);
  
  void           FillSelectedClusterHistograms(AliVCluster* cluster, Float_t pt,
                                               Int_t nLocMax,        Int_t tag,
                                               Float_t asy = 0);
  
  void           FillWeightHistograms(AliVCluster *clus);
  
  void           HasPairSameMCMother(Int_t label1 , Int_t label2,
                                     Int_t tag1   , Int_t tag2,
                                     Int_t & label, Int_t & tag);
  
  void           MakeInvMassInCalorimeter() ;
  
  void           MakeInvMassInCalorimeterAndCTS() ;
  
  void           MakeShowerShapeIdentification() ;
  
  // Setters Getters
  
  //
  /// Analysis types
  //
  enum anaTypes
  {
      kIMCalo,        ///<  2 calorimeter clusters invariant mass selection
      kSSCalo,        ///<  1 calorimeter cluster shower shape and split invariatn mass selection
      kIMCaloTracks   ///<  1 calorimeter cluster and 1 photon conversion pair invariant mass selection
  };
    
  anaTypes       GetAnalysisType()                     const { return fAnaType                 ; }
  void           SetAnalysisType(anaTypes ana)               { fAnaType = ana                  ; }
    
  //
  // Used by all analysis types
  //
  TString        GetInputAODGammaConvName()            const { return fInputAODGammaConvName   ; }
  void           SetInputAODGammaConvName(TString name)      { fInputAODGammaConvName = name   ; }

  void           SetNLMCut(Int_t min, Int_t max)             { fNLMCutMin = min;
                                                               fNLMCutMax = max                ; }
  Int_t          GetNLMCutMin()                        const { return fNLMCutMin               ; }
  Int_t          GetNLMCutMax()                        const { return fNLMCutMax               ; }

  void           SwitchOnAllNLMHistoFill()                   { fFillAllNLMHistograms   = kTRUE ; }
  void           SwitchOffAllNLMHistoFill()                  { fFillAllNLMHistograms   = kFALSE; }

  void           SwitchOnTMHistoFill()                       { fFillTMHisto           = kTRUE  ; }
  void           SwitchOffTMHistoFill()                      { fFillTMHisto           = kFALSE ; }
  
  void           SwitchOnSelectedClusterHistoFill()          { fFillSelectClHisto     = kTRUE  ; }
  void           SwitchOffSelectedClusterHistoFill()         { fFillSelectClHisto     = kFALSE ; }

  //
  // Only for invariant mass analysis
  //
  void           SetM02CutForInvMass(Float_t min=0, Float_t max=10)  
                                                             { fM02MinCutForIM         = min   ; 
                                                               fM02MaxCutForIM         = max   ; }

  void           SetR(Float_t r)                             { fR = r                          ; }
  void           SetIsolationCandidateMinPt(Float_t min)     { fIsoCandMinPt = min             ; }
  
  void           SwitchOnSelectIsolatedDecay()               { fSelectIsolatedDecay    = kTRUE ; }
  void           SwitchOffSelectIsolatedDecay()              { fSelectIsolatedDecay    = kFALSE; }
  
  void           SwitchOnSelectPairInIsolationCone()         { fSelectPairInIsoCone    = kTRUE ; }
  void           SwitchOffSelectPairInIsolationCone()        { fSelectPairInIsoCone    = kFALSE; }

  //
  // Only for pi0 SS identification case
  //
  void           SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) 
                                                             { fMinDist  = m1 ; fMinDist2 = m2 ; 
                                                               fMinDist3 = m3                  ; }
  
  
  void           SetNLMMinEnergy(Int_t i, Float_t min)       { if (i < 3 && i >=0 ) fNLMECutMin[i]  = min   ; }
  Float_t        GetNLMMinEnergy(Int_t i) const              { if( i < 3 && i >=0 ) return fNLMECutMin[i]   ;  
                                                               else return 0 ; }
  
  void           SetTimeCut(Double_t min, Double_t max)      { fTimeCutMin = min;
                                                               fTimeCutMax = max               ; }
  Double_t       GetTimeCutMin()                       const { return fTimeCutMin              ; }
  Double_t       GetTimeCutMax()                       const { return fTimeCutMax              ; }
  
  Bool_t         IsTrackMatchRejectionOn()             const { return fRejectTrackMatch        ; }
  void           SwitchOnTrackMatchRejection()               { fRejectTrackMatch      = kTRUE  ; }
  void           SwitchOffTrackMatchRejection()              { fRejectTrackMatch      = kFALSE ; }
  
  void           SwitchOnFillWeightHistograms()              { fFillWeightHistograms  = kTRUE  ; }
  void           SwitchOffFillWeightHistograms()             { fFillWeightHistograms  = kFALSE ; }
  
  void           SwitchOnOnlySimpleSSHistoFill()             { fFillOnlySimpleSSHisto = kTRUE  ; }
  void           SwitchOffOnlySimpleHistoFill()              { fFillOnlySimpleSSHisto = kFALSE ; }
  
  void           SwitchOnFillEMCALBCHistograms()             { fFillEMCALBCHistograms = kTRUE  ; }
  void           SwitchOffFillEMCALBCHistograms()            { fFillEMCALBCHistograms = kFALSE ; }
  
  void           SwitchOnSplitClusterDistToBad()             { fCheckSplitDistToBad   = kTRUE  ; }
  void           SwitchOffSplitClusterDistToBad()            { fCheckSplitDistToBad   = kFALSE ; }
      
  
  /// For MC histograms in arrays, index in the array corresponds to a MC originating particle type.
  enum mcTypes   { kmcPi0      = 0, kmcEta      = 1, kmcPhoton           = 2,
                   kmcPi0Decay = 3, kmcEtaDecay = 4, kmcOtherDecay       = 5,
                   kmcElectron = 6, kmcHadron   = 7                          } ;
  
  /// Total number of MC origin histograms.
  static const Int_t fgkNmcTypes = 8;
  
private:
  
  anaTypes       fAnaType;                                  ///<  Select analysis type.
  
  Int_t          fNLMCutMin  ;                              ///<  Remove clusters/cells with number of local maxima smaller than this value.
  Int_t          fNLMCutMax  ;                              ///<  Remove clusters/cells with number of local maxima larger than this value.
  Bool_t         fFillAllNLMHistograms;                     ///<  Fill all NLM dependent histograms.
  Bool_t         fFillTMHisto;                              ///<  Fill track matching plots.
  Bool_t         fFillSelectClHisto;                        ///<  Fill selected cluster histograms.

  // Invariant mass analysis
  Float_t        fM02MaxCutForIM   ;                        ///<  Study photon clusters with l0 smaller than cut, in inv. mass analysis.
  Float_t        fM02MinCutForIM   ;                        ///<  Study photon clusters with l0 larger than cut, in inv. mass analysis.
  Bool_t         fSelectIsolatedDecay;                      ///<  Select pairs where at least one is declared isolated (run first AliAnaParticleIsolation).
  Bool_t         fSelectPairInIsoCone;                      ///<  Select pair in isolation cone.
  Float_t        fR;                                        ///<  Isolation cone.
  Float_t        fIsoCandMinPt;                             ///<  Isolation candidate minimum pT.
  
  // Only for pi0 SS identification case, kSSCalo
  Float_t        fMinDist ;                                 ///<  Minimal distance to bad channel to accept cluster.
  Float_t        fMinDist2;                                 ///<  Cuts on Minimal distance to study acceptance evaluation.
  Float_t        fMinDist3;                                 ///<  One more cut on distance used for acceptance-efficiency study.
  Float_t        fNLMECutMin[3] ;                           ///<  Minimum energy of the cluster, depending on NLM.
  Double_t       fTimeCutMin  ;                             ///<  Remove clusters/cells with time smaller than this value, in ns.
  Double_t       fTimeCutMax  ;                             ///<  Remove clusters/cells with time larger than this value, in ns.
  Bool_t         fRejectTrackMatch ;                        ///<  Remove clusters which have an associated TPC track.
  Bool_t         fCheckSplitDistToBad;                      ///<  Check the distance to bad channel and to EMCal borders of split clusters.
  Bool_t         fFillWeightHistograms ;                    ///<  Fill weigth histograms.
  Bool_t         fFillOnlySimpleSSHisto;                    ///<  Fill selected cluster histograms, selected SS histograms.
  Bool_t         fFillEMCALBCHistograms;                    ///<  Fill eta-phi BC dependent histograms.
  
  // Only for combination of calorimeter and conversion photons, kIMCaloTracks
  TString        fInputAODGammaConvName;                    ///<  Name of AOD branch with conversion photons.
  
  ///<  cluster/particle kinematic temporal containers
  TLorentzVector fMomentum;                                 //!<! Cluster/pi0 momentum, kinematic temporal containers.
  TLorentzVector fMomentum1;                                //!<! Cluster/photon momentum, kinematic temporal containers.
  TLorentzVector fMomentum2;                                //!<! Cluster/photon momentum, kinematic temporal containers.
  TLorentzVector fMomentum12;                               //!<! Cluster/pi0 momentum, sum 1+2, kinematic temporal containers.
  TLorentzVector fPrimaryMom;                               //!<! Primary momentum, kinematic temporal containers.
  TLorentzVector fGrandMotherMom;                           //!<! Primary momentum, kinematic temporal containers.
  
  // Histograms
  
  TH1F         * fhPt  ;                                    //!<! Number of identified  pi0/eta vs pT
  TH1F         * fhE   ;                                    //!<! Number of identified  pi0/eta vs E
  TH2F         * fhPtEta  ;                                 //!<! Pt vs eta of identified  pi0/eta
  TH2F         * fhPtPhi  ;                                 //!<! Pt vs phi of identified  pi0/eta
  TH2F         * fhEtaPhi  ;                                //!<! eta vs phi of identified  pi0/eta
  TH2F         * fhEtaPhiEMCALBC0  ;                        //!<! Pseudorapidity vs Phi of clusters
  TH2F         * fhEtaPhiEMCALBC1  ;                        //!<! Pseudorapidity vs Phi of clusters
  TH2F         * fhEtaPhiEMCALBCN  ;                        //!<! Pseudorapidity vs Phi of clusters
  
  TH2F         * fhEtaPhiTriggerEMCALBC[11]  ;              //!<! Pseudorapidity vs Phi of pi0 for E > 2
  TH2F         * fhTimeTriggerEMCALBC  [11]  ;              //!<! Time distribution of pi0, when trigger is in a given BC
  TH2F         * fhTimeTriggerEMCALBCPileUpSPD[11] ;        //!<! Time distribution of pi0, when trigger is in a given BC, tagged as pile-up SPD
  TH2F         * fhEtaPhiTriggerEMCALBCUM[11]  ;            //!<! Pseudorapidity vs Phi of pi0 for E > 2, not matched to trigger
  TH2F         * fhTimeTriggerEMCALBCUM[11]  ;              //!<! Time distribution of pi0, when trigger is in a given BC, not matched to trigger
  
  TH2F         * fhTimeTriggerEMCALBC0UMReMatchOpenTime   ; //!<! Time distribution of pi0s in event, when trigger is not found, rematched open time trigger
  TH2F         * fhTimeTriggerEMCALBC0UMReMatchCheckNeigh ; //!<! Time distribution of pi0s in event, when trigger is not found, rematched with neigbour patchs
  TH2F         * fhTimeTriggerEMCALBC0UMReMatchBoth       ; //!<! Time distribution of pi0s in event, when trigger is not found, rematched open both
  
  TH2F         * fhPtCentrality ;                           //!<! Centrality  vs pi0/eta pT
  TH2F         * fhPtEventPlane ;                           //!<! Event plane vs pi0/eta pT
  TH2F         * fhMCPtCentrality[fgkNmcTypes];             //!<! Centrality  vs pi0/eta pT  coming from X
  
  TH1F         * fhPtReject  ;                              //!<! Number of rejected as  pi0/eta vs pT
  TH1F         * fhEReject   ;                              //!<! Number of rejected as  pi0/eta vs E
  TH2F         * fhPtEtaReject  ;                           //!<! pT vs eta of rejected as  pi0/eta
  TH2F         * fhPtPhiReject  ;                           //!<! pT vs phi of rejected as  pi0/eta
  TH2F         * fhEtaPhiReject  ;                          //!<! eta vs phi of rejected as  pi0/eta
  
  TH2F         * fhMass  ;                                  //!<! Pair mass vs E, for all pairs
  TH2F         * fhMassPt  ;                                //!<! Pair mass vs pT, for all pairs
  TH2F         * fhMassPtMaxPair  ;                         //!<! Pair mass vs pT max of the pair, for all pairs
  TH2F         * fhMassPtMinPair  ;                         //!<! Pair mass vs pT min of the pair, for all pairs
  TH2F         * fhMassSplitPt  ;                           //!<! Pair mass vs pT (split), for all pairs
  TH2F         * fhSelectedMass  ;                          //!<! Pair mass vs E, for selected pairs
  TH2F         * fhSelectedMassPt  ;                        //!<! Pair mass vs pT, for selected pairs
  TH2F         * fhSelectedMassSplitPt  ;                   //!<! Pair mass vs pT (split), for selected pairs
  
  TH2F         * fhMassPtIsoRCut  ;                         //!<! Pair mass vs pT, for all pairs when opening angle not larger than iso cone radius
  
  TH2F         * fhMassPtLocMax[3] ;                        //!<! Pair mass vs pT, for all pairs, for each NLM case
  TH2F         * fhSelectedMassPtLocMax[3] ;                //!<! Pair mass vs pT, for selected pairs, for each NLM case
  TH2F         * fhSelectedMassPtLocMaxSM[3][22];           //!<! Pair mass vs pT, for selected pairs, for each NLM case, for each SM
  TH2F         * fhMCSelectedMassPtLocMax[fgkNmcTypes][3] ; //!<! Pair mass vs pT, for selected pairs, vs originating particle
  
  TH2F         * fhSelectedLambda0PtLocMaxSM[3][22];        //!<! Pair mass vs pT, for selected pairs, for each NLM case, for each SM
  
  TH2F         * fhMassNoOverlap  ;                         //!<! Pair mass vs E, for all pairs, no overlap
  TH2F         * fhMassPtNoOverlap  ;                       //!<! Pair mass vs pT, for all pairs, no overlap
  TH2F         * fhMassSplitPtNoOverlap  ;                  //!<! Pair mass vs pT (split), for all pairs, no overlap
  TH2F         * fhSelectedMassNoOverlap  ;                 //!<! Pair mass vs E, for selected pairs, no overlap
  TH2F         * fhSelectedMassPtNoOverlap  ;               //!<! Pair mass vs pT, for selected pairs, no overlap
  TH2F         * fhSelectedMassSplitPtNoOverlap  ;          //!<! Pair mass vs pT (split), for selected pairs, no overlap
    
  TH2F         * fhMCPi0PtRecoPtPrim;                       //!<! pt reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaPtRecoPtPrim;                       //!<! pt reco vs pt prim for eta mother
  TH2F         * fhMCPi0PtRecoPtPrimNoOverlap;              //!<! pt reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaPtRecoPtPrimNoOverlap;              //!<! pt reco vs pt prim for eta mother
  
  TH2F         * fhMCPi0SplitPtRecoPtPrim;                  //!<! pt split reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSplitPtRecoPtPrim;                  //!<! pt split reco vs pt prim for eta mother
  TH2F         * fhMCPi0SplitPtRecoPtPrimNoOverlap;         //!<! pt split reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSplitPtRecoPtPrimNoOverlap;         //!<! pt split reco vs pt prim for eta mother
  
  TH2F         * fhMCPi0SelectedPtRecoPtPrim;               //!<! pt reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSelectedPtRecoPtPrim;               //!<! pt reco vs pt prim for eta mother
  TH2F         * fhMCPi0SelectedPtRecoPtPrimNoOverlap;      //!<! pt reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSelectedPtRecoPtPrimNoOverlap;      //!<! pt reco vs pt prim for eta mother
  
  TH2F         * fhMCPi0SelectedSplitPtRecoPtPrim;          //!<! pt split reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSelectedSplitPtRecoPtPrim;          //!<! pt split reco vs pt prim for eta mother
  TH2F         * fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap; //!<! pt split reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap; //!<! pt split reco vs pt prim for eta mother
  
  TH2F         * fhMCPi0PtRecoPtPrimLocMax[3];              //!<! pt reco vs pt prim for pi0 mother, vs NLM
  TH2F         * fhMCEtaPtRecoPtPrimLocMax[3];              //!<! pt reco vs pt prim for eta mother, vs NLM
  TH2F         * fhMCPi0SplitPtRecoPtPrimLocMax[3];         //!<! pt split reco vs pt prim for pi0 mother, vs NLM
  TH2F         * fhMCEtaSplitPtRecoPtPrimLocMax[3];         //!<! pt split reco vs pt prim for eta mother, vs NLM
  
  TH2F         * fhMCPi0SelectedPtRecoPtPrimLocMax[3];      //!<! pt reco vs pt prim for pi0 mother, vs NLM
  TH2F         * fhMCEtaSelectedPtRecoPtPrimLocMax[3];      //!<! pt reco vs pt prim for eta mother, vs NLM
  TH2F         * fhMCPi0SelectedSplitPtRecoPtPrimLocMax[3]; //!<! pt split reco vs pt prim for pi0 mother, vs NLM
  TH2F         * fhMCEtaSelectedSplitPtRecoPtPrimLocMax[3]; //!<! pt split reco vs pt prim for eta mother, vs NLM
  
  TH2F         * fhAsymmetry ;                              //!<! Cluster pT vs asymmetry of 2 splitted clusters
  TH2F         * fhSelectedAsymmetry  ;                     //!<! Cluster pT vs asymmetry of 2 splitted clusters, for selected pairs
  TH1F         * fhSplitE  ;                                //!<! Split sub-cluster pair energy sum
  TH1F         * fhSplitPt  ;                               //!<! Split sub-cluster pair pT sum
  TH2F         * fhSplitPtEta  ;                            //!<! Split sub-cluster pair pT sum vs eta
  TH2F         * fhSplitPtPhi  ;                            //!<! Split sub-cluster pair pT sum vs phi
  TH2F         * fhNLocMaxSplitPt  ;                        //!<! Split sub-cluster pair pT sum, as a function of n maxima
  
  TH1F         * fhPtDecay  ;                               //!<! Number of identified  pi0/eta decay photons vs pT
  
  TH2F         * fhPtDispersion ;                           //!<! pT vs disp of selected cluster
  TH2F         * fhPtLambda0 ;                              //!<! pT vs lambda0 of selected cluster
  TH2F         * fhPtLambda0NoSplitCut ;                    //!<! pT vs lambda0 of cluster before the split selection.
  TH2F         * fhPtLambda1 ;                              //!<! pT vs lambda1 of selected cluster
  TH2F         * fhPtLambda0NoTRD ;                         //!<! pT vs lambda0 of selected cluster, not behind TRD
  TH2F         * fhPtLambda0FracMaxCellCut ;                //!<! pT vs lambda0 of selected cluster, fraction of cluster energy in max cell cut
  TH2F         * fhPtFracMaxCell ;                          //!<! pT vs frac max cell of selected cluster
  TH2F         * fhPtFracMaxCellNoTRD ;                     //!<! pT vs frac max cell of selected cluster, not behind TRD
  TH2F         * fhPtNCells;                                //!<! pT vs N cells in selected cluster
  TH2F         * fhPtTime;                                  //!<! pT vs Time of selected cluster
  TH2F         * fhEPairDiffTime;                           //!<! E pair vs Pair of clusters time difference vs E
  
  TH2F         * fhPtDispEta ;                              //!<! Shower dispersion in eta direction
  TH2F         * fhPtDispPhi ;                              //!<! Shower dispersion in phi direction
  TH2F         * fhLambda0DispEta[7] ;                      //!<! Shower shape correlation l0 vs disp eta
  TH2F         * fhLambda0DispPhi[7] ;                      //!<! Shower shape correlation l0 vs disp phi
  TH2F         * fhPtSumEta ;                               //!<! Shower dispersion in eta direction
  TH2F         * fhPtSumPhi ;                               //!<! Shower dispersion in phi direction
  TH2F         * fhPtSumEtaPhi ;                            //!<! Shower dispersion in eta and phi direction
  TH2F         * fhPtDispEtaPhiDiff ;                       //!<! Shower dispersion eta - phi
  TH2F         * fhPtSphericity ;                           //!<! Shower sphericity in eta vs phi
  TH2F         * fhDispEtaDispPhi[7] ;                      //!<! Shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F         * fhAsymmetryLambda0[7] ;                    //!<! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhAsymmetryDispEta[7] ;                    //!<! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhAsymmetryDispPhi[7] ;                    //!<! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  
  // MC histograms
  
  TH1F         * fhMCPtDecay            [fgkNmcTypes];      //!<! pT from MC particle
  TH1F         * fhMCPtDecayLostPairPi0;                    //!<! pT for tagged clustres when MC Pi0 Decay, when companion is lost
  TH1F         * fhMCPtDecayLostPairEta;                    //!<! pT for tagged clustres when MC Eta Decay, when companion is lost
  TH2F         * fhMCPtLambda0          [fgkNmcTypes];      //!<! pT vs lambda0 of pi0 pairs but really from MC particle
  TH2F         * fhMCPtLambda1          [fgkNmcTypes];      //!<! pT vs lambda1 of pi0 pairs but really from MC particle
  TH2F         * fhMCPtDispersion       [fgkNmcTypes];      //!<! pT vs dispersion of pi0 pairs but really from MC particle
  TH2F         * fhMCPtLambda0NoTRD     [fgkNmcTypes];      //!<! pT vs lambda0 of pi0 pairs but really from MC particle, not behind TRD
  TH2F         * fhMCPtLambda0FracMaxCellCut[fgkNmcTypes];  //!<! pT vs lambda0 of pi0 pairs but really from MC particle, fraction of cluster energy in max cell cut
  TH2F         * fhMCPtFracMaxCell      [fgkNmcTypes];      //!<! pT vs fraction of max cell
  TH2F         * fhMCPtDispEta          [fgkNmcTypes];      //!<! Shower dispersion in eta direction
  TH2F         * fhMCPtDispPhi          [fgkNmcTypes];      //!<! Shower dispersion in phi direction
  TH2F         * fhMCLambda0DispEta  [7][fgkNmcTypes];      //!<! Shower shape correlation l0 vs disp eta
  TH2F         * fhMCLambda0DispPhi  [7][fgkNmcTypes];      //!<! Shower shape correlation l0 vs disp phi
  TH2F         * fhMCPtSumEtaPhi        [fgkNmcTypes];      //!<! Shower dispersion in eta vs phi direction
  TH2F         * fhMCPtDispEtaPhiDiff   [fgkNmcTypes];      //!<! Shower dispersion in eta -phi direction
  TH2F         * fhMCPtSphericity       [fgkNmcTypes];      //!<! Shower sphericity, eta vs phi
  TH2F         * fhMCDispEtaDispPhi  [7][fgkNmcTypes];      //!<! Shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F         * fhMCPtAsymmetry        [fgkNmcTypes];      //!<! E asymmetry of 2 splitted clusters vs cluster pT
  TH2F         * fhMCAsymmetryLambda0[7][fgkNmcTypes];      //!<! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhMCAsymmetryDispEta[7][fgkNmcTypes];      //!<! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhMCAsymmetryDispPhi[7][fgkNmcTypes];      //!<! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  
  TH1F         * fhMCE                  [fgkNmcTypes];      //!<! Number of identified as pi0 vs E coming from X
  TH1F         * fhMCPt                 [fgkNmcTypes];      //!<! Number of identified as pi0 vs Pt coming from X
  TH2F         * fhMCPtPhi              [fgkNmcTypes];      //!<! pt vs phi of identified as pi0, coming from X
  TH2F         * fhMCPtEta              [fgkNmcTypes];      //!<! pt vs eta of identified as pi0, coming from X
  TH1F         * fhMCEReject            [fgkNmcTypes];      //!<! Number of rejected as pi0 vs E coming from X
  TH1F         * fhMCPtReject           [fgkNmcTypes];      //!<! Number of rejected as pi0 vs Pt coming from X
  
  TH1F         * fhMCSplitE             [fgkNmcTypes];      //!<! Number of identified as pi0 vs sum E  split coming from X
  TH1F         * fhMCSplitPt            [fgkNmcTypes];      //!<! Number of identified as pi0 vs sum Pt split coming from X
  TH2F         * fhMCSplitPtPhi         [fgkNmcTypes];      //!<! pt vs phi of identified as pi0, coming from X
  TH2F         * fhMCSplitPtEta         [fgkNmcTypes];      //!<! pt vs eta of identified as pi0, coming from X
  TH2F         * fhMCNLocMaxSplitPt     [fgkNmcTypes];      //!<! Number of identified as pi0 vs sum Pt split coming from X, for different NLM
  
  TH2F         * fhMCMassPt             [fgkNmcTypes];      //!<! Pair pT vs Mass coming from X
  TH2F         * fhMCMassSplitPt        [fgkNmcTypes];      //!<! Pair pT (split) vs Mass coming from X
  TH2F         * fhMCSelectedMassPt     [fgkNmcTypes];      //!<! selected pair pT vs Mass coming from X
  TH2F         * fhMCSelectedMassSplitPt[fgkNmcTypes];      //!<! selected pair pT (split) vs Mass coming from X
  
  TH2F         * fhMCMassPtNoOverlap     [fgkNmcTypes];     //!<! Pair pT vs Mass coming from X, no random particles overlap
  TH2F         * fhMCMassSplitPtNoOverlap[fgkNmcTypes];     //!<! Pair pT (split) vs Mass coming from X, no random particles overlap
  TH2F         * fhMCSelectedMassPtNoOverlap[fgkNmcTypes];  //!<! selected pair pT vs Mass coming from X, no random particles overlap
  TH2F         * fhMCSelectedMassSplitPtNoOverlap[fgkNmcTypes]; //!<! selected pair pT (split) vs Mass coming from X, no random particles overlap
  
  TH2F         * fhMCPi0PtGenRecoFraction;                  //!<! SS id, clusters id as pi0 (eta), coming from 2 photon, pi0 primary, pt vs E prim pi0 / E reco
  TH2F         * fhMCEtaPtGenRecoFraction;                  //!<! SS id, clusters id as pi0 (eta), coming from 2 photon, eta primary, pt vs E prim eta / E reco
  TH1F         * fhMCPi0DecayPt;                            //!<! SS id, clusters id as pi0 (eta), coming from 1 photon, pi0 decay primary, pt
  TH2F         * fhMCPi0DecayPtFraction;                    //!<! SS id, clusters id as pi0 (eta), coming from 1 photon, pi0 decay primary, pt vs pt decay / pt mother
  TH1F         * fhMCEtaDecayPt;                            //!<! SS id, clusters id as pi0 (eta), coming from 1 photon, eta decay primary, pt
  TH2F         * fhMCEtaDecayPtFraction;                    //!<! SS id, clusters id as pi0 (eta), coming from 1 photon, eta decay primary, pt vs pt decay / pt mother
  TH1F         * fhMCOtherDecayPt;                          //!<! SS id, clusters id as pi0 (eta), coming from 1 photon, other decay primary, pt
  
  TH2F         * fhMassPairMCPi0;                           //!<! Pair mass, origin is same pi0
  TH2F         * fhMassPairMCEta;                           //!<! Pair mass, origin is same eta
  TH2F         * fhAnglePairMCPi0;                          //!<! Pair opening angle, origin is same pi0
  TH2F         * fhAnglePairMCEta;                          //!<! Pair opening angle, origin is same eta
  
  TH2F         * fhMCPi0PtOrigin ;                          //!<! Mass of reoconstructed pi0 pairs  in calorimeter vs mother
  TH2F         * fhMCEtaPtOrigin ;                          //!<! Mass of reoconstructed pi0 pairs  in calorimeter vs mother
  TH2F         * fhMCNotResonancePi0PtOrigin ;              //!<! Mass of reoconstructed pi0 pairs  in calorimeter vs mother
  TH2F         * fhMCPi0PtStatus ;                          //!<! Mass of reoconstructed pi0 pairs  in calorimeter vs mother
  TH2F         * fhMCPi0ProdVertex;                         //!<! Spectrum of selected pi0 vs production vertex
  TH2F         * fhMCEtaProdVertex;                         //!<! Spectrum of selected eta vs production vertex
  
  // Weight studies
  
  TH2F         * fhECellClusterRatio;                       //!<! E cell / e cluster vs e cluster for selected photons
  TH2F         * fhECellClusterLogRatio;                    //!<! Log (e cell / e cluster)  vs e cluster for selected photons
  TH2F         * fhEMaxCellClusterRatio;                    //!<! E max cell / e cluster vs e cluster for selected photons
  TH2F         * fhEMaxCellClusterLogRatio;                 //!<! Log (e max cell / e cluster) vs e cluster for selected photons
  TH2F         * fhLambda0ForW0[14];                        //!<! L0 for 7 defined w0= 3, 3.5 ... 6 for selected photons
  //TH2F         * fhLambda1ForW0[7];                       //!<! L1 for 7 defined w0= 3, 3.5 ... 6 for selected photons
  
  // Track Matching
  TH2F         * fhTrackMatchedDEta     ;                   //!<! Eta distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDPhi     ;                   //!<! Phi distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDEtaDPhi ;                   //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV
  TH2F         * fhTrackMatchedDEtaPos  ;                   //!<! Eta distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDPhiPos  ;                   //!<! Phi distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDEtaDPhiPos ;                //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV
  TH2F         * fhTrackMatchedDEtaNeg  ;                   //!<! Eta distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDPhiNeg  ;                   //!<! Phi distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDEtaDPhiNeg ;                //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV
  
  TH2F         * fhTrackMatchedMCParticlePt;                //!<! Trace origin of matched particle, energy
  TH2F         * fhTrackMatchedMCParticleDEta;              //!<! Trace origin of matched particle, eta residual
  TH2F         * fhTrackMatchedMCParticleDPhi;              //!<! Trace origin of matched particle, phi residual
  TH2F         * fhdEdx  ;                                  //!<! Matched track dEdx vs cluster E
  TH2F         * fhEOverP;                                  //!<! Matched track E cluster over P track vs cluster E
  TH2F         * fhEOverPNoTRD;                             //!<! Matched track E cluster over P track vs cluster E, not behind TRD
  
  // Local maxima
  TH2F         * fhNLocMaxPt;                               //!<! number of maxima in selected clusters
  TH2F         * fhNLocMaxPtSM[22] ;                        //!<! Number of maxima in selected clusters, per super module
  TH2F         * fhMCNLocMaxPt[fgkNmcTypes];                //!<! Number of maxima in selected clusters, vs originating particle
  TH2F         * fhPtLambda0LocMax[3] ;                     //!<! pT vs lambda0 of selected cluster, 1,2,>2 local maxima in cluster
  TH2F         * fhMCPtLambda0LocMax[fgkNmcTypes][3];       //!<! pT vs lambda0 of selected cluster, 1,2,>2 local maxima in cluster, vs originating particle
  TH2F         * fhPtLambda1LocMax[3] ;                     //!<! pT vs lambda1 of selected cluster, 1,2,>2 local maxima in cluster
  TH2F         * fhPtDispersionLocMax[3] ;                  //!<! pT vs lambda1 of selected cluster, 1,2,>2 local maxima in cluster
  TH2F         * fhPtDispEtaLocMax[3] ;                     //!<! pT vs eta dispersion of selected cluster, 1,2,>2 local maxima in cluster
  TH2F         * fhPtDispPhiLocMax[3] ;                     //!<! pT vs phi dispersion of selected cluster, 1,2,>2 local maxima in cluster
  TH2F         * fhPtSumEtaPhiLocMax[3] ;                   //!<! pT vs dispersion in eta and phi direction
  TH2F         * fhPtDispEtaPhiDiffLocMax[3];               //!<! pT vs dispersion eta - phi
  TH2F         * fhPtSphericityLocMax[3] ;                  //!<! pT vs sphericity in eta vs phi
  TH2F         * fhPtAsymmetryLocMax[3] ;                   //!<! E asymmetry of 2 splitted clusters vs cluster E for different NLM
  
  TH2F         * fhMassPairLocMax[8];                       //!<! Pair mass, origin is same pi0, combine clusters depending on number of maxima
  
  TH2F         * fhNLocMaxPtReject;                         //!<! Number of maxima in selected clusters
  TH2F         * fhMCNLocMaxPtReject[fgkNmcTypes];          //!<! Number of maxima in selected clusters
  
  // Pile-up
  TH1F         * fhPtPileUp[7];                             //!<! pT distribution of selected pi0/eta
  TH2F         * fhPtCellTimePileUp[7];                     //!<! pT vs Time inside cluster, before any selection, not max cell
  TH2F         * fhPtTimeDiffPileUp[7];                     //!<! pT vs Time difference inside cluster, before any selection
  TH2F         * fhTimePtNoCut;                             //!<! Time of cluster vs pT, no cut
  TH2F         * fhTimePtSPD;                               //!<! Time of cluster vs pT, IsSPDPileUp
  TH2F         * fhTimePtSPDMulti;                          //!<! Time of cluster vs pT, IsSPDPileUpMulti
  TH2F         * fhTimeNPileUpVertSPD;                      //!<! Time of cluster vs n pile-up vertices from SPD
  TH2F         * fhTimeNPileUpVertTrack;                    //!<! Time of cluster vs n pile-up vertices from Tracks
  TH2F         * fhTimeNPileUpVertContributors;             //!<! Time of cluster vs n pile-up vertex from SPD contributors
  TH2F         * fhTimePileUpMainVertexZDistance;           //!<! Time of cluster vs difference of z main vertex and pile-up vertex
  TH2F         * fhTimePileUpMainVertexZDiamond;            //!<! Time of cluster vs difference of z diamond and pile-up vertex
  
  TH2F         * fhPtNPileUpSPDVtx;                         //!<! Cluster pt vs number of spd pile-up vertices
  TH2F         * fhPtNPileUpTrkVtx;                         //!<! Cluster pt vs number of track pile-up vertices
  TH2F         * fhPtNPileUpSPDVtxTimeCut;                  //!<! Cluster pt vs number of spd pile-up vertices, time cut +-25 ns
  TH2F         * fhPtNPileUpTrkVtxTimeCut;                  //!<! Cluster pt vs number of track pile-up vertices, time cut +- 25 ns
  TH2F         * fhPtNPileUpSPDVtxTimeCut2;                 //!<! Cluster pt vs number of spd pile-up vertices, time cut +-75 ns
  TH2F         * fhPtNPileUpTrkVtxTimeCut2;                 //!<! Cluster pt vs number of track pile-up vertices, time cut +- 75 ns
  
  /// Copy constructor not implemented.
  AliAnaPi0EbE(              const AliAnaPi0EbE & pi0ebe) ;
    
  /// Assignment operator not implemented.
  AliAnaPi0EbE & operator = (const AliAnaPi0EbE & pi0ebe) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaPi0EbE,42) ;
  /// \endcond

} ;

#endif //ALIANAPI0EBE_H




