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
  
  Int_t          GetMCIndex(const Int_t aodTag);
  
  void           Init();
  
  void           InitParameters();

  void           MakeAnalysisFillAOD()  ;
   
  void           MakeAnalysisFillHistograms() ; 
  
  void           Print(const Option_t * opt) const;
  
  // Main
  
  void           FillPileUpHistograms(const Float_t energy, const Float_t time) ;
  
  void           FillRejectedClusterHistograms(const TLorentzVector mom, const Int_t mctag);
  
  void           FillSelectedClusterHistograms(AliVCluster* cluster, 
                                               const Int_t nLocMax,
                                               const Int_t tag,
                                               const Float_t asy = 0);
    
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

  
  //Only for combination of calorimeter and conversion photons, kIMCaloTracks
  TString        fInputAODGammaConvName;   //  Name of AOD branch with conversion photons
  
  //Histograms
  
  TH1F         * fhPt  ;                   //! Number of identified  pi0/eta vs pT
  TH1F         * fhE   ;                   //! Number of identified  pi0/eta vs E
  TH2F         * fhEEta  ;                 //! E vs eta of identified  pi0/eta 
  TH2F         * fhEPhi  ;                 //! E vs phi of identified  pi0/eta 
  TH2F         * fhEtaPhi  ;               //! eta vs phi of identified  pi0/eta 

  TH2F         * fhPtCentrality ;          //! centrality  vs pi0/eta pT
  TH2F         * fhPtEventPlane ;          //! event plane vs pi0/eta pT
  
  TH1F         * fhPtReject  ;             //! Number of rejected as  pi0/eta vs pT
  TH1F         * fhEReject   ;             //! Number of rejected as  pi0/eta vs E
  TH2F         * fhEEtaReject  ;           //! E vs eta of rejected as  pi0/eta 
  TH2F         * fhEPhiReject  ;           //! E vs phi of rejected as  pi0/eta 
  TH2F         * fhEtaPhiReject  ;         //! eta vs phi of rejected as  pi0/eta 
  
  TH2F         * fhMass  ;                 //! pair mass vs E, for all pairs
  TH2F         * fhAsymmetry ;             //! cluster E vs asymmetry of 2 splitted clusters 
  TH2F         * fhSelectedMass  ;         //! pair mass vs E, for selected pairs
  TH2F         * fhSelectedAsymmetry  ;    //! cluster E vs asymmetry of 2 splitted clusters, for selected pairs
  TH1F         * fhSplitE  ;       //! split sub-cluster pair energy sum
  TH1F         * fhSplitPt  ;      //! split sub-cluster pair pT sum
  
  TH1F         * fhPtDecay  ;              //! Number of identified  pi0/eta decay photons vs pT
  TH1F         * fhEDecay   ;              //! Number of identified  pi0/eta decay photons vs E
  
  TH2F         * fhEDispersion ;           //! E vs disp of selected cluster
  TH2F         * fhELambda0 ;              //! E vs lambda0 of selected cluster 
  TH2F         * fhELambda1 ;              //! E vs lambda1 of selected cluster 
  TH2F         * fhELambda0NoTRD ;         //! E vs lambda0 of selected cluster, not behind TRD 
  TH2F         * fhELambda0FracMaxCellCut ;//! E vs lambda0 of selected cluster, fraction of cluster energy in max cell cut 
  TH2F         * fhEFracMaxCell ;          //! E vs frac max cell of selected cluster 
  TH2F         * fhEFracMaxCellNoTRD ;     //! E vs frac max cell of selected cluster, not behind TRD  
  TH2F         * fhENCells;                //! E vs N cells in selected cluster
  TH2F         * fhETime;                  //! E vs Time of selected cluster 
  TH2F         * fhEPairDiffTime;          //! E vs Pair of clusters time difference vs E
  
  TH2F         * fhDispEtaE ;              //! shower dispersion in eta direction
  TH2F         * fhDispPhiE ;              //! shower dispersion in phi direction
  TH2F         * fhLambda0DispEta[7] ;     //! shower shape correlation l0 vs disp eta
  TH2F         * fhLambda0DispPhi[7] ;     //! shower shape correlation l0 vs disp phi
  TH2F         * fhSumEtaE ;               //! shower dispersion in eta direction
  TH2F         * fhSumPhiE ;               //! shower dispersion in phi direction
  TH2F         * fhSumEtaPhiE ;            //! shower dispersion in eta and phi direction
  TH2F         * fhDispEtaPhiDiffE ;       //! shower dispersion eta - phi
  TH2F         * fhSphericityE ;           //! shower sphericity in eta vs phi
  TH2F         * fhDispEtaDispPhi[7] ;     //! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F         * fhAsymmetryLambda0[7] ;   //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhAsymmetryDispEta[7] ;   //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhAsymmetryDispPhi[7] ;   //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins

  //MC histograms
  
  TH2F         * fhEMCLambda0[6] ;            //! E vs lambda0 of pi0 pairs but really from MC particle
  TH2F         * fhEMCLambda1[6] ;            //! E vs lambda1 of pi0 pairs but really from MC particle
  TH2F         * fhEMCDispersion[6] ;         //! E vs dispersion of pi0 pairs but really from MC particle
  TH2F         * fhEMCLambda0NoTRD[6] ;         //! E vs lambda0 of pi0 pairs but really from MC particle, not behind TRD
  TH2F         * fhEMCLambda0FracMaxCellCut[6] ;//! E vs lambda0 of pi0 pairs but really from MC particle, fraction of cluster energy in max cell cut
  TH2F         * fhEMCFracMaxCell[6] ;        //! E vs fraction of max cell 
  
  TH2F         * fhMCEDispEta[6] ;            //! shower dispersion in eta direction
  TH2F         * fhMCEDispPhi[6] ;            //! shower dispersion in phi direction
  TH2F         * fhMCLambda0DispEta[7][6] ;   //! shower shape correlation l0 vs disp eta
  TH2F         * fhMCLambda0DispPhi[7][6] ;   //! shower shape correlation l0 vs disp phi
  TH2F         * fhMCESumEtaPhi[6] ;          //! shower dispersion in eta vs phi direction
  TH2F         * fhMCEDispEtaPhiDiff[6] ;     //! shower dispersion in eta -phi direction
  TH2F         * fhMCESphericity[6] ;         //! shower sphericity, eta vs phi
  TH2F         * fhMCDispEtaDispPhi[7][6] ;   //! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F         * fhMCEAsymmetry[6] ;          //! E asymmetry of 2 splitted clusters vs cluster E
  TH2F         * fhMCAsymmetryLambda0[7][6] ; //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhMCAsymmetryDispEta[7][6] ; //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhMCAsymmetryDispPhi[7][6] ; //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  
  TH1F         * fhMCE[6];                    //! Number of identified as pi0 vs E coming from X
  TH1F         * fhMCPt[6];                   //! Number of identified as pi0 vs Pt coming from X
  TH2F         * fhMCPhi[6];                  //! Phi of identified as pi0, coming from X
  TH2F         * fhMCEta[6];                  //! eta of identified as pi0, coming from X
  TH1F         * fhMCEReject[6];              //! Number of rejected as pi0 vs E coming from X
  TH1F         * fhMCPtReject[6];             //! Number of rejected as pi0 vs Pt coming from X

  TH1F         * fhMCSplitE[6];               //! Number of identified as pi0 vs sum E  split coming from X
  TH1F         * fhMCSplitPt[6];              //! Number of identified as pi0 vs sum Pt split coming from X

  
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
  TH2F         * fhTrackMatchedMCParticleE;    //! Trace origin of matched particle, energy
  TH2F         * fhTrackMatchedMCParticleDEta; //! Trace origin of matched particle, eta residual
  TH2F         * fhTrackMatchedMCParticleDPhi; //! Trace origin of matched particle, phi residual
  TH2F         * fhdEdx  ;                 //! matched track dEdx vs cluster E
  TH2F         * fhEOverP;                 //! matched track E cluster over P track vs cluster E
  TH2F         * fhEOverPNoTRD;                 //! matched track E cluster over P track vs cluster E, not behind TRD 

  // Local maxima
  TH2F         * fhNLocMax;                //! number of maxima in selected clusters
  TH2F         * fhELambda0LocMax[3] ;     //! E vs lambda0 of selected cluster, 1,2,>2 local maxima in cluster 
  TH2F         * fhELambda1LocMax[3] ;     //! E vs lambda1 of selected cluster, 1,2,>2 local maxima in cluster 
  TH2F         * fhEDispersionLocMax[3] ;  //! E vs lambda1 of selected cluster, 1,2,>2 local maxima in cluster 
  TH2F         * fhEDispEtaLocMax[3] ;     //! E vs eta dispersion of selected cluster, 1,2,>2 local maxima in cluster 
  TH2F         * fhEDispPhiLocMax[3] ;     //! E vs phi dispersion of selected cluster, 1,2,>2 local maxima in cluster 
  TH2F         * fhESumEtaPhiLocMax[3] ;   //! E vs dispersion in eta and phi direction
  TH2F         * fhEDispEtaPhiDiffLocMax[3] ; //! E vs dispersion eta - phi
  TH2F         * fhESphericityLocMax[3] ;  //! E vs sphericity in eta vs phi  
  TH2F         * fhEAsymmetryLocMax[3] ;   //! E asymmetry of 2 splitted clusters vs cluster E for different NLM

  TH2F         * fhMassPairLocMax[8];      //! pair mass, origin is same pi0, combine clusters depending on number of maxima

  // Pile-up
  TH1F         * fhPtPi0PileUp[7];                //! pT distribution of selected pi0/eta
  TH2F         * fhTimeENoCut;                    //! time of cluster vs E, no cut
  TH2F         * fhTimeESPD;                      //! time of cluster vs E, IsSPDPileUp
  TH2F         * fhTimeESPDMulti;                 //! time of cluster vs E, IsSPDPileUpMulti
  TH2F         * fhTimeNPileUpVertSPD;            //! time of cluster vs n pile-up vertices from SPD
  TH2F         * fhTimeNPileUpVertTrack;          //! time of cluster vs n pile-up vertices from Tracks
  TH2F         * fhTimeNPileUpVertContributors;   //! time of cluster vs n pile-up vertex from SPD contributors
  TH2F         * fhTimePileUpMainVertexZDistance; //! time of cluster vs difference of z main vertex and pile-up vertex 
  TH2F         * fhTimePileUpMainVertexZDiamond;  //! time of cluster vs difference of z diamond and pile-up vertex 
  
  AliAnaPi0EbE(              const AliAnaPi0EbE & pi0ebe) ; // cpy ctor
  AliAnaPi0EbE & operator = (const AliAnaPi0EbE & pi0ebe) ; // cpy assignment
  
  ClassDef(AliAnaPi0EbE,25)
} ;


#endif //ALIANAPI0EBE_H



