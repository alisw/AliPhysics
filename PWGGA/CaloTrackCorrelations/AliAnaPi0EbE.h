#ifndef ALIANAPI0EBE_H
#define ALIANAPI0EBE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
//
// Class for the analysis of high pT pi0 event by event
// Pi0/Eta identified by one of the following:
//  -Invariant mass of 2 cluster in calorimeter
//  -Shower shape analysis in calorimeter
//  -Invariant mass of one cluster in calorimeter and one photon reconstructed in TPC (in near future)
//
//-- Author: Gustavo Conesa (INFN-LNF)  &  Raphaelle Ichou (SUBATECH)
//_________________________________________________________________________


// --- ROOT system ---
class TList ;
class TObjString;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaPi0EbE : public AliAnaCaloTrackCorrBaseClass {

 public: 
  AliAnaPi0EbE() ; // default ctor
  virtual ~AliAnaPi0EbE() { ; } //virtual dtor
	  
  TObjString *   GetAnalysisCuts();
  
  TList      *   GetCreateOutputObjects();
  
  Int_t          GetMCIndex(Int_t aodTag);
  
  void           Init();
  
  void           InitParameters();

  void           MakeAnalysisFillAOD()  ;
   
  void           MakeAnalysisFillHistograms() ; 
  
  void           Print(const Option_t * opt) const;
  
  // Main
  
  void           FillPileUpHistograms(Float_t pt, Float_t time, AliVCluster * c) ;
  
  void           FillRejectedClusterHistograms(TLorentzVector mom, Int_t mctag, Int_t nMaxima);
  
  void           FillSelectedClusterHistograms(AliVCluster* cluster, Float_t pt,
                                               Int_t nLocMax,        Int_t tag,
                                               Float_t asy = 0);
    
  void           FillWeightHistograms(AliVCluster *clus);
    
  void           HasPairSameMCMother(AliAODPWG4Particle * photon1, 
                                     AliAODPWG4Particle * photon2, 
                                     Int_t & label, Int_t & tag);
  
  void           MakeInvMassInCalorimeter() ;
  
  void           MakeInvMassInCalorimeterAndCTS() ;
  
  void           MakeShowerShapeIdentification() ;
          
  //Setters Getters
  
  //Analysis types
  enum anaTypes  {kIMCalo, kSSCalo, kIMCaloTracks};  
  anaTypes       GetAnalysisType()                     const { return fAnaType                 ; }
  void           SetAnalysisType(anaTypes ana)               { fAnaType = ana                  ; }
  
  TString        GetInputAODGammaConvName()            const { return fInputAODGammaConvName   ; }
  void           SetInputAODGammaConvName(TString name)      { fInputAODGammaConvName = name   ; }	
  
  //Only for pi0 SS identification case
  void           SetCalorimeter(TString & det)               { fCalorimeter = det              ; }
  
  void           SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) {
                  fMinDist = m1; fMinDist2 = m2; fMinDist3 = m3                                ; }
  
  void           SetNLMCut(Int_t min, Int_t max)             { fNLMCutMin = min; 
                                                               fNLMCutMax = max                ; }
  Int_t          GetNLMCutMin()                        const { return fNLMCutMin               ; }
  Int_t          GetNLMCutMax()                        const { return fNLMCutMax               ; }	
  
  void           SetNLMMinEnergy(Int_t i, Float_t min)       { if (i < 3 && i >=0 ) fNLMECutMin[i]  = min   ; }
  Float_t        GetNLMMinEnergy(Int_t i) const              { if( i < 3 && i >=0 ) return fNLMECutMin[i]   ;  else return 0 ; }

  void           SetTimeCut(Double_t min, Double_t max)      { fTimeCutMin = min;
                                                               fTimeCutMax = max               ; }
  Double_t       GetTimeCutMin()                       const { return fTimeCutMin              ; }
  Double_t       GetTimeCutMax()                       const { return fTimeCutMax              ; }
 
  Bool_t         IsTrackMatchRejectionOn()             const { return fRejectTrackMatch        ; }
  void           SwitchOnTrackMatchRejection()               { fRejectTrackMatch      = kTRUE  ; }
  void           SwitchOffTrackMatchRejection()              { fRejectTrackMatch      = kFALSE ; }
  
  void           SwitchOnFillPileUpHistograms()              { fFillPileUpHistograms  = kTRUE  ; }
  void           SwitchOffFillPileUpHistograms()             { fFillPileUpHistograms  = kFALSE ; }    
    
  void           SwitchOnFillWeightHistograms()              { fFillWeightHistograms  = kTRUE  ; }
  void           SwitchOffFillWeightHistograms()             { fFillWeightHistograms  = kFALSE ; }  
  
  void           SwitchOnTMHistoFill()                       { fFillTMHisto           = kTRUE  ; }
  void           SwitchOffTMHistoFill()                      { fFillTMHisto           = kFALSE ; }

  void           SwitchOnSelectedClusterHistoFill()          { fFillSelectClHisto     = kTRUE  ; }
  void           SwitchOffSelectedClusterHistoFill()         { fFillSelectClHisto     = kFALSE ; }
  
  void           SwitchOnOnlySimpleSSHistoFill()             { fFillOnlySimpleSSHisto = kTRUE  ; }
  void           SwitchOffOnlySimpleHistoFill()              { fFillOnlySimpleSSHisto = kFALSE ; }

  void           SwitchOnFillEMCALBCHistograms()             { fFillEMCALBCHistograms = kTRUE  ; }
  void           SwitchOffFillEMCALBCHistograms()            { fFillEMCALBCHistograms = kFALSE ; }
  
  void           SwitchOnSplitClusterDistToBad()             { fCheckSplitDistToBad   = kTRUE  ; }
  void           SwitchOffSplitClusterDistToBad()            { fCheckSplitDistToBad   = kFALSE ; }

  void           SetNumberOfSuperModules(Int_t nSM)          { fNSuperModules         = nSM    ; }

  
  //For histograms
  enum mcTypes   { kmcPhoton = 0, kmcConversion = 1, kmcPi0    = 2,  
                   kmcEta    = 3, kmcElectron   = 4, kmcHadron = 5 };

 private:
  
  anaTypes       fAnaType;                 // Select analysis type
  
  //Only for pi0 SS identification case, kSSCalo
  TString        fCalorimeter ;            // Calorimeter where the gamma is searched;
  Float_t        fMinDist ;                // Minimal distance to bad channel to accept cluster
  Float_t        fMinDist2;                // Cuts on Minimal distance to study acceptance evaluation
  Float_t        fMinDist3;                // One more cut on distance used for acceptance-efficiency study
  Int_t          fNLMCutMin  ;             // Remove clusters/cells with number of local maxima smaller than this value
  Int_t          fNLMCutMax  ;             // Remove clusters/cells with number of local maxima larger than this value
  Float_t        fNLMECutMin[3] ;          // Minimum energy of the cluster, depending on nlm.
  Double_t       fTimeCutMin  ;            // Remove clusters/cells with time smaller than this value, in ns
  Double_t       fTimeCutMax  ;            // Remove clusters/cells with time larger than this value, in ns
  Bool_t         fRejectTrackMatch ;       // Remove clusters which have an associated TPC track

  Bool_t         fFillPileUpHistograms;    // Fill pile-up related histograms
  Bool_t         fFillWeightHistograms ;   // Fill weigth histograms
  Bool_t         fFillTMHisto;             // Fill track matching plots
  Bool_t         fFillSelectClHisto;       // Fill selected cluster histograms
  Bool_t         fFillOnlySimpleSSHisto;   // Fill selected cluster histograms, selected SS histograms
  Bool_t         fFillEMCALBCHistograms;   // Fill eta-phi BC dependent histograms

  //Only for combination of calorimeter and conversion photons, kIMCaloTracks
  TString        fInputAODGammaConvName;   //  Name of AOD branch with conversion photons

  Bool_t         fCheckSplitDistToBad;     // Check the distance to bad channel and to EMCal borders of split clusters
  
  Int_t          fNSuperModules;           // Number of supermodules
  
  //Histograms
  
  TH1F         * fhPt  ;                   //! Number of identified  pi0/eta vs pT
  TH1F         * fhE   ;                   //! Number of identified  pi0/eta vs E
  TH2F         * fhPtEta  ;                //! Pt vs eta of identified  pi0/eta
  TH2F         * fhPtPhi  ;                //! Pt vs phi of identified  pi0/eta
  TH2F         * fhEtaPhi  ;               //! eta vs phi of identified  pi0/eta
  TH2F         * fhEtaPhiEMCALBC0  ;       //! Pseudorapidity vs Phi of clusters 
  TH2F         * fhEtaPhiEMCALBC1  ;       //! Pseudorapidity vs Phi of clusters 
  TH2F         * fhEtaPhiEMCALBCN  ;       //! Pseudorapidity vs Phi of clusters 

  TH2F         * fhEtaPhiTriggerEMCALBC[11]  ;    //! Pseudorapidity vs Phi of pi0 for E > 2
  TH2F         * fhTimeTriggerEMCALBC  [11]  ;    //! Time distribution of pi0, when trigger is in a given BC
  TH2F         * fhTimeTriggerEMCALBCPileUpSPD[11] ; //! Time distribution of pi0, when trigger is in a given BC, tagged as pile-up SPD
  TH2F         * fhEtaPhiTriggerEMCALBCUM[11]  ;  //! Pseudorapidity vs Phi of pi0 for E > 2, not matched to trigger
  TH2F         * fhTimeTriggerEMCALBCUM[11]  ;    //! Time distribution of pi0, when trigger is in a given BC, not matched to trigger

  TH2F         * fhTimeTriggerEMCALBC0UMReMatchOpenTime   ; //! Time distribution of pi0s in event, when trigger is not found, rematched open time trigger
  TH2F         * fhTimeTriggerEMCALBC0UMReMatchCheckNeigh ; //! Time distribution of pi0s in event, when trigger is not found, rematched with neigbour patchs
  TH2F         * fhTimeTriggerEMCALBC0UMReMatchBoth       ; //! Time distribution of pi0s in event, when trigger is not found, rematched open both

  TH2F         * fhPtCentrality ;          //! centrality  vs pi0/eta pT
  TH2F         * fhPtEventPlane ;          //! event plane vs pi0/eta pT
  
  TH1F         * fhPtReject  ;             //! Number of rejected as  pi0/eta vs pT
  TH1F         * fhEReject   ;             //! Number of rejected as  pi0/eta vs E
  TH2F         * fhPtEtaReject  ;          //! pT vs eta of rejected as  pi0/eta
  TH2F         * fhPtPhiReject  ;          //! pT vs phi of rejected as  pi0/eta
  TH2F         * fhEtaPhiReject  ;         //! eta vs phi of rejected as  pi0/eta 
  
  TH2F         * fhMass  ;                 //! pair mass vs E, for all pairs
  TH2F         * fhMassPt  ;               //! pair mass vs pT, for all pairs
  TH2F         * fhMassSplitPt  ;          //! pair mass vs pT (split), for all pairs
  TH2F         * fhSelectedMass  ;         //! pair mass vs E, for selected pairs
  TH2F         * fhSelectedMassPt  ;       //! pair mass vs pT, for selected pairs
  TH2F         * fhSelectedMassSplitPt  ;  //! pair mass vs pT (split), for selected pairs
  
  TH2F         * fhMassPtLocMax[3] ;       //! pair mass vs pT, for all pairs
  TH2F         * fhSelectedMassPtLocMax[3] ;//! pair mass vs pT, for selected pairs
  TH2F         * fhMCSelectedMassPtLocMax[6][3] ;//! pair mass vs pT, for selected pairs, vs originating particle

  TH2F         * fhMassNoOverlap  ;                 //! pair mass vs E, for all pairs, no overlap
  TH2F         * fhMassPtNoOverlap  ;               //! pair mass vs pT, for all pairs, no overlap
  TH2F         * fhMassSplitPtNoOverlap  ;          //! pair mass vs pT (split), for all pairs, no overlap
  TH2F         * fhSelectedMassNoOverlap  ;         //! pair mass vs E, for selected pairs, no overlap
  TH2F         * fhSelectedMassPtNoOverlap  ;       //! pair mass vs pT, for selected pairs, no overlap
  TH2F         * fhSelectedMassSplitPtNoOverlap  ;  //! pair mass vs pT (split), for selected pairs, no overlap

  TH2F         * fhMCPi0PtRecoPtPrim;      //! pt reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaPtRecoPtPrim;      //! pt reco vs pt prim for eta mother
  TH2F         * fhMCPi0PtRecoPtPrimNoOverlap; //! pt reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaPtRecoPtPrimNoOverlap; //! pt reco vs pt prim for eta mother

  TH2F         * fhMCPi0SplitPtRecoPtPrim;      //! pt split reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSplitPtRecoPtPrim;      //! pt split reco vs pt prim for eta mother
  TH2F         * fhMCPi0SplitPtRecoPtPrimNoOverlap; //! pt split reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSplitPtRecoPtPrimNoOverlap; //! pt split reco vs pt prim for eta mother

  TH2F         * fhMCPi0SelectedPtRecoPtPrim;      //! pt reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSelectedPtRecoPtPrim;      //! pt reco vs pt prim for eta mother
  TH2F         * fhMCPi0SelectedPtRecoPtPrimNoOverlap; //! pt reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSelectedPtRecoPtPrimNoOverlap; //! pt reco vs pt prim for eta mother
  
  TH2F         * fhMCPi0SelectedSplitPtRecoPtPrim;      //! pt split reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSelectedSplitPtRecoPtPrim;      //! pt split reco vs pt prim for eta mother
  TH2F         * fhMCPi0SelectedSplitPtRecoPtPrimNoOverlap; //! pt split reco vs pt prim for pi0 mother
  TH2F         * fhMCEtaSelectedSplitPtRecoPtPrimNoOverlap; //! pt split reco vs pt prim for eta mother
  
  TH2F         * fhAsymmetry ;             //! cluster pT vs asymmetry of 2 splitted clusters
  TH2F         * fhSelectedAsymmetry  ;    //! cluster pT vs asymmetry of 2 splitted clusters, for selected pairs
  TH1F         * fhSplitE  ;               //! split sub-cluster pair energy sum
  TH1F         * fhSplitPt  ;              //! split sub-cluster pair pT sum
  TH2F         * fhSplitPtEta  ;           //! split sub-cluster pair pT sum vs eta
  TH2F         * fhSplitPtPhi  ;           //! split sub-cluster pair pT sum vs phi
  TH2F         * fhNLocMaxSplitPt  ;       //! split sub-cluster pair pT sum, as a function of n maxima
  
  TH1F         * fhPtDecay  ;              //! Number of identified  pi0/eta decay photons vs pT
  TH1F         * fhEDecay   ;              //! Number of identified  pi0/eta decay photons vs E
  
  TH2F         * fhPtDispersion ;           //! pT vs disp of selected cluster
  TH2F         * fhPtLambda0 ;              //! pT vs lambda0 of selected cluster 
  TH2F         * fhPtLambda1 ;              //! pT vs lambda1 of selected cluster 
  TH2F         * fhPtLambda0NoTRD ;         //! pT vs lambda0 of selected cluster, not behind TRD 
  TH2F         * fhPtLambda0FracMaxCellCut ;//! pT vs lambda0 of selected cluster, fraction of cluster energy in max cell cut 
  TH2F         * fhPtFracMaxCell ;          //! pT vs frac max cell of selected cluster 
  TH2F         * fhPtFracMaxCellNoTRD ;     //! pT vs frac max cell of selected cluster, not behind TRD  
  TH2F         * fhPtNCells;                //! pT vs N cells in selected cluster
  TH2F         * fhPtTime;                  //! pT vs Time of selected cluster 
  TH2F         * fhEPairDiffTime;           //! E pair vs Pair of clusters time difference vs E
  
  TH2F         * fhPtDispEta ;              //! shower dispersion in eta direction
  TH2F         * fhPtDispPhi ;              //! shower dispersion in phi direction
  TH2F         * fhLambda0DispEta[7] ;      //! shower shape correlation l0 vs disp eta
  TH2F         * fhLambda0DispPhi[7] ;      //! shower shape correlation l0 vs disp phi
  TH2F         * fhPtSumEta ;               //! shower dispersion in eta direction
  TH2F         * fhPtSumPhi ;               //! shower dispersion in phi direction
  TH2F         * fhPtSumEtaPhi ;            //! shower dispersion in eta and phi direction
  TH2F         * fhPtDispEtaPhiDiff ;       //! shower dispersion eta - phi
  TH2F         * fhPtSphericity ;           //! shower sphericity in eta vs phi
  TH2F         * fhDispEtaDispPhi[7] ;      //! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F         * fhAsymmetryLambda0[7] ;    //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhAsymmetryDispEta[7] ;    //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhAsymmetryDispPhi[7] ;    //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins

  //MC histograms
  
  TH2F         * fhMCPtLambda0[6] ;              //! pT vs lambda0 of pi0 pairs but really from MC particle
  TH2F         * fhMCPtLambda1[6] ;              //! pT vs lambda1 of pi0 pairs but really from MC particle
  TH2F         * fhMCPtDispersion[6] ;           //! pT vs dispersion of pi0 pairs but really from MC particle
  TH2F         * fhMCPtLambda0NoTRD[6] ;         //! pT vs lambda0 of pi0 pairs but really from MC particle, not behind TRD
  TH2F         * fhMCPtLambda0FracMaxCellCut[6] ;//! pT vs lambda0 of pi0 pairs but really from MC particle, fraction of cluster energy in max cell cut
  TH2F         * fhMCPtFracMaxCell[6] ;       //! pT vs fraction of max cell
  
  TH2F         * fhMCPtDispEta[6] ;           //! shower dispersion in eta direction
  TH2F         * fhMCPtDispPhi[6] ;           //! shower dispersion in phi direction
  TH2F         * fhMCLambda0DispEta[7][6] ;   //! shower shape correlation l0 vs disp eta
  TH2F         * fhMCLambda0DispPhi[7][6] ;   //! shower shape correlation l0 vs disp phi
  TH2F         * fhMCPtSumEtaPhi[6] ;         //! shower dispersion in eta vs phi direction
  TH2F         * fhMCPtDispEtaPhiDiff[6] ;    //! shower dispersion in eta -phi direction
  TH2F         * fhMCPtSphericity[6] ;        //! shower sphericity, eta vs phi
  TH2F         * fhMCDispEtaDispPhi[7][6] ;   //! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F         * fhMCPtAsymmetry[6] ;         //! E asymmetry of 2 splitted clusters vs cluster pT
  TH2F         * fhMCAsymmetryLambda0[7][6] ; //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhMCAsymmetryDispEta[7][6] ; //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhMCAsymmetryDispPhi[7][6] ; //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  
  TH1F         * fhMCE[6];                    //! Number of identified as pi0 vs E coming from X
  TH1F         * fhMCPt[6];                   //! Number of identified as pi0 vs Pt coming from X
  TH2F         * fhMCPtPhi[6];                //! pt vs phi of identified as pi0, coming from X
  TH2F         * fhMCPtEta[6];                //! pt vs eta of identified as pi0, coming from X
  TH1F         * fhMCEReject[6];              //! Number of rejected as pi0 vs E coming from X
  TH1F         * fhMCPtReject[6];             //! Number of rejected as pi0 vs Pt coming from X

  TH1F         * fhMCSplitE[6];               //! Number of identified as pi0 vs sum E  split coming from X
  TH1F         * fhMCSplitPt[6];              //! Number of identified as pi0 vs sum Pt split coming from X
  TH2F         * fhMCSplitPtPhi[6];           //! pt vs phi of identified as pi0, coming from X
  TH2F         * fhMCSplitPtEta[6];           //! pt vs eta of identified as pi0, coming from X
  TH2F         * fhMCNLocMaxSplitPt[6];       //! Number of identified as pi0 vs sum Pt split coming from X, for different NLM
  
  TH2F         * fhMCMassPt[6];               //! pair pT vs Mass coming from X
  TH2F         * fhMCMassSplitPt[6];          //! pair pT (split) vs Mass coming from X
  TH2F         * fhMCSelectedMassPt[6];       //! selected pair pT vs Mass coming from X
  TH2F         * fhMCSelectedMassSplitPt[6];  //! selected pair pT (split) vs Mass coming from X

  TH2F         * fhMCMassPtNoOverlap[6];               //! pair pT vs Mass coming from X, no random particles overlap
  TH2F         * fhMCMassSplitPtNoOverlap[6];          //! pair pT (split) vs Mass coming from X, no random particles overlap
  TH2F         * fhMCSelectedMassPtNoOverlap[6];       //! selected pair pT vs Mass coming from X, no random particles overlap
  TH2F         * fhMCSelectedMassSplitPtNoOverlap[6];  //! selected pair pT (split) vs Mass coming from X, no random particles overlap
  
  TH2F         * fhMCPtCentrality[6] ;        //! centrality  vs pi0/eta pT  coming from X
  
  TH2F         * fhMCPi0PtGenRecoFraction;    //! SS id, clusters id as pi0 (eta), coming from 2 photon, pi0 primary, pt vs E prim pi0 / E reco
  TH2F         * fhMCEtaPtGenRecoFraction;    //! SS id, clusters id as pi0 (eta), coming from 2 photon, eta primary, pt vs E prim eta / E reco  
  TH1F         * fhMCPi0DecayPt;              //! SS id, clusters id as pi0 (eta), coming from 1 photon, pi0 decay primary, pt
  TH2F         * fhMCPi0DecayPtFraction;      //! SS id, clusters id as pi0 (eta), coming from 1 photon, pi0 decay primary, pt vs pt decay / pt mother
  TH1F         * fhMCEtaDecayPt;              //! SS id, clusters id as pi0 (eta), coming from 1 photon, eta decay primary, pt
  TH2F         * fhMCEtaDecayPtFraction;      //! SS id, clusters id as pi0 (eta), coming from 1 photon, eta decay primary, pt vs pt decay / pt mother  
  TH1F         * fhMCOtherDecayPt;            //! SS id, clusters id as pi0 (eta), coming from 1 photon, other decay primary, pt

  TH2F         * fhMassPairMCPi0;             //! pair mass, origin is same pi0
  TH2F         * fhMassPairMCEta;             //! pair mass, origin is same eta
  TH2F         * fhAnglePairMCPi0;            //! pair opening angle, origin is same pi0
  TH2F         * fhAnglePairMCEta;            //! pair opening angle, origin is same eta
  
  // Weight studies
  
  TH2F         * fhECellClusterRatio;      //! e cell / e cluster vs e cluster for selected photons
  TH2F         * fhECellClusterLogRatio;   //! log (e cell / e cluster)  vs e cluster for selected photons
  TH2F         * fhEMaxCellClusterRatio;   //! e max cell / e cluster vs e cluster for selected photons
  TH2F         * fhEMaxCellClusterLogRatio;//! log (e max cell / e cluster) vs e cluster for selected photons
  TH2F         * fhLambda0ForW0[14];       //! L0 for 7 defined w0= 3, 3.5 ... 6 for selected photons
  //TH2F         * fhLambda1ForW0[7];        //! L1 for 7 defined w0= 3, 3.5 ... 6 for selected photons  
  
  // Track Matching
  TH2F         * fhTrackMatchedDEta     ;  //! Eta distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDPhi     ;  //! Phi distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDEtaDPhi ;  //! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV
  TH2F         * fhTrackMatchedDEtaPos  ;  //! Eta distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDPhiPos  ;  //! Phi distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDEtaDPhiPos ; //! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV
  TH2F         * fhTrackMatchedDEtaNeg  ;  //! Eta distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDPhiNeg  ;  //! Phi distance between track and cluster vs cluster E
  TH2F         * fhTrackMatchedDEtaDPhiNeg ; //! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV
  
  TH2F         * fhTrackMatchedMCParticlePt;   //! Trace origin of matched particle, energy
  TH2F         * fhTrackMatchedMCParticleDEta; //! Trace origin of matched particle, eta residual
  TH2F         * fhTrackMatchedMCParticleDPhi; //! Trace origin of matched particle, phi residual
  TH2F         * fhdEdx  ;                 //! matched track dEdx vs cluster E
  TH2F         * fhEOverP;                 //! matched track E cluster over P track vs cluster E
  TH2F         * fhEOverPNoTRD;                 //! matched track E cluster over P track vs cluster E, not behind TRD 

  // Local maxima
  TH2F         * fhNLocMaxPt;               //! number of maxima in selected clusters
  TH2F         * fhNLocMaxPtSM[22] ;        //! Pt of identified clusters, vs NLM
  TH2F         * fhMCNLocMaxPt[6];          //! number of maxima in selected clusters, vs originating particle
  TH2F         * fhPtLambda0LocMax[3] ;     //! pT vs lambda0 of selected cluster, 1,2,>2 local maxima in cluster
  TH2F         * fhMCPtLambda0LocMax[6][3] ;//! pT vs lambda0 of selected cluster, 1,2,>2 local maxima in cluster, vs originating particle
  TH2F         * fhPtLambda1LocMax[3] ;     //! pT vs lambda1 of selected cluster, 1,2,>2 local maxima in cluster
  TH2F         * fhPtDispersionLocMax[3] ;  //! pT vs lambda1 of selected cluster, 1,2,>2 local maxima in cluster 
  TH2F         * fhPtDispEtaLocMax[3] ;     //! pT vs eta dispersion of selected cluster, 1,2,>2 local maxima in cluster 
  TH2F         * fhPtDispPhiLocMax[3] ;     //! pT vs phi dispersion of selected cluster, 1,2,>2 local maxima in cluster 
  TH2F         * fhPtSumEtaPhiLocMax[3] ;   //! pT vs dispersion in eta and phi direction
  TH2F         * fhPtDispEtaPhiDiffLocMax[3] ; //! pT vs dispersion eta - phi
  TH2F         * fhPtSphericityLocMax[3] ;  //! pT vs sphericity in eta vs phi  
  TH2F         * fhPtAsymmetryLocMax[3] ;   //! E asymmetry of 2 splitted clusters vs cluster E for different NLM

  TH2F         * fhMassPairLocMax[8];      //! pair mass, origin is same pi0, combine clusters depending on number of maxima
  
  TH2F         * fhNLocMaxPtReject;              //! number of maxima in selected clusters
  TH2F         * fhMCNLocMaxPtReject[6];         //! number of maxima in selected clusters
  
  // Pile-up
  TH1F         * fhPtPileUp[7];                   //! pT distribution of selected pi0/eta
  TH2F         * fhPtCellTimePileUp[7];           //! pT vs Time inside cluster, before any selection, not max cell
  TH2F         * fhPtTimeDiffPileUp[7];           //! pT vs Time difference inside cluster, before any selection
  TH2F         * fhTimePtNoCut;                   //! time of cluster vs pT, no cut
  TH2F         * fhTimePtSPD;                     //! time of cluster vs pT, IsSPDPileUp
  TH2F         * fhTimePtSPDMulti;                //! time of cluster vs pT, IsSPDPileUpMulti
  TH2F         * fhTimeNPileUpVertSPD;            //! time of cluster vs n pile-up vertices from SPD
  TH2F         * fhTimeNPileUpVertTrack;          //! time of cluster vs n pile-up vertices from Tracks
  TH2F         * fhTimeNPileUpVertContributors;   //! time of cluster vs n pile-up vertex from SPD contributors
  TH2F         * fhTimePileUpMainVertexZDistance; //! time of cluster vs difference of z main vertex and pile-up vertex 
  TH2F         * fhTimePileUpMainVertexZDiamond;  //! time of cluster vs difference of z diamond and pile-up vertex 
  
  TH2F         * fhPtNPileUpSPDVtx;	              //! cluster pt vs number of spd pile-up vertices
  TH2F         * fhPtNPileUpTrkVtx;               //! cluster pt vs number of track pile-up vertices
  TH2F         * fhPtNPileUpSPDVtxTimeCut;	      //! cluster pt vs number of spd pile-up vertices, time cut +-25 ns
  TH2F         * fhPtNPileUpTrkVtxTimeCut;        //! cluster pt vs number of track pile-up vertices, time cut +- 25 ns 		
  TH2F         * fhPtNPileUpSPDVtxTimeCut2;	      //! cluster pt vs number of spd pile-up vertices, time cut +-75 ns
  TH2F         * fhPtNPileUpTrkVtxTimeCut2;       //! cluster pt vs number of track pile-up vertices, time cut +- 75 ns
  
  AliAnaPi0EbE(              const AliAnaPi0EbE & pi0ebe) ; // cpy ctor
  AliAnaPi0EbE & operator = (const AliAnaPi0EbE & pi0ebe) ; // cpy assignment
  
  ClassDef(AliAnaPi0EbE,35)
} ;


#endif //ALIANAPI0EBE_H



