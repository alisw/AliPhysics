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
  
  void           Init();
  
  void           InitParameters();

  void           MakeAnalysisFillAOD()  ;
   
  void           MakeAnalysisFillHistograms() ; 
  
  void           Print(const Option_t * opt) const;
  
  // Main
  
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
  anaTypes       GetAnalysisType()                     const { return fAnaType                ; }
  void           SetAnalysisType(anaTypes ana)               { fAnaType = ana                 ; }
  
  TString        GetInputAODGammaConvName()            const { return fInputAODGammaConvName  ; }
  void           SetInputAODGammaConvName(TString name)      { fInputAODGammaConvName = name  ; }	
  
  //Only for pi0 SS identification case
  void           SetCalorimeter(TString & det)               { fCalorimeter = det             ; }
  
  void           SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) {
                  fMinDist = m1; fMinDist2 = m2; fMinDist3 = m3                               ; }
  
  void           SetTimeCut(Double_t min, Double_t max)      { fTimeCutMin = min; 
                                                               fTimeCutMax = max              ; }
  Double_t       GetTimeCutMin()                       const { return fTimeCutMin             ; }
  Double_t       GetTimeCutMax()                       const { return fTimeCutMax             ; }	

  void           SwitchOnFillWeightHistograms()              { fFillWeightHistograms = kTRUE  ; }
  void           SwitchOffFillWeightHistograms()             { fFillWeightHistograms = kFALSE ; }  
  
  void           SwitchOnTMHistoFill()                       { fFillTMHisto          = kTRUE  ; }
  void           SwitchOffTMHistoFill()                      { fFillTMHisto          = kFALSE ; }

  void           SwitchOnSelectedClusterHistoFill()          { fFillSelectClHisto    = kTRUE  ; }
  void           SwitchOffSelectedClusterHistoFill()         { fFillSelectClHisto    = kFALSE ; }
  
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
  Double_t       fTimeCutMin  ;            // Remove clusters/cells with time smaller than this value, in ns
  Double_t       fTimeCutMax  ;            // Remove clusters/cells with time larger than this value, in ns
  
  Bool_t         fFillWeightHistograms ;   // Fill weigth histograms
  Bool_t         fFillTMHisto;             // Fill track matching plots
  Bool_t         fFillSelectClHisto;       // Fill selected cluster histograms

  //Only for combination of calorimeter and conversion photons, kIMCaloTracks
  TString        fInputAODGammaConvName;   //  Name of AOD branch with conversion photons
  
  //Histograms
  
  TH1F         * fhPt  ;                   //! Number of identified  pi0/eta vs pT
  TH1F         * fhE   ;                   //! Number of identified  pi0/eta vs E
  TH2F         * fhEEta  ;                 //! E vs eta of identified  pi0/eta 
  TH2F         * fhEPhi  ;                 //! E vs phi of identified  pi0/eta 
  TH2F         * fhEtaPhi  ;               //! eta vs phi of identified  pi0/eta 

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
  TH2F         * fhLambda0DispEta[5] ;     //! shower shape correlation l0 vs disp eta
  TH2F         * fhLambda0DispPhi[5] ;     //! shower shape correlation l0 vs disp phi
  TH2F         * fhSumEtaE ;               //! shower dispersion in eta direction
  TH2F         * fhSumPhiE ;               //! shower dispersion in phi direction
  TH2F         * fhSumEtaPhiE ;            //! shower dispersion in eta and phi direction
  TH2F         * fhDispEtaPhiDiffE ;       //! shower dispersion eta - phi
  TH2F         * fhSphericityE ;           //! shower sphericity in eta vs phi
  TH2F         * fhDispEtaDispPhi[5] ;     //! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F         * fhAsymmetryE ;            //! E asymmetry of 2 splitted clusters vs cluster E
  TH2F         * fhAsymmetryLambda0[5] ;   //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhAsymmetryDispEta[5] ;   //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhAsymmetryDispPhi[5] ;   //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins

  //MC histograms
  
  TH2F         * fhEMCLambda0[6] ;         //! E vs lambda0 of pi0 pairs but really from MC particle
  TH2F         * fhEMCLambda1[6] ;         //! E vs lambda1 of pi0 pairs but really from MC particle
  TH2F         * fhEMCDispersion[6] ;      //! E vs dispersion of pi0 pairs but really from MC particle
  TH2F         * fhEMCLambda0NoTRD[6] ;         //! E vs lambda0 of pi0 pairs but really from MC particle, not behind TRD
  TH2F         * fhEMCLambda0FracMaxCellCut[6] ;//! E vs lambda0 of pi0 pairs but really from MC particle, fraction of cluster energy in max cell cut
  TH2F         * fhEMCFracMaxCell[6] ;     //! E vs fraction of max cell 
  
  TH2F         * fhMCEDispEta[6] ;         //! shower dispersion in eta direction
  TH2F         * fhMCEDispPhi[6] ;         //! shower dispersion in phi direction
  TH2F         * fhMCLambda0DispEta[5][6] ;//! shower shape correlation l0 vs disp eta
  TH2F         * fhMCLambda0DispPhi[5][6] ;//! shower shape correlation l0 vs disp phi
  TH2F         * fhMCESumEtaPhi[6] ;       //! shower dispersion in eta vs phi direction
  TH2F         * fhMCEDispEtaPhiDiff[6] ;  //! shower dispersion in eta -phi direction
  TH2F         * fhMCESphericity[6] ;      //! shower sphericity, eta vs phi
  TH2F         * fhMCDispEtaDispPhi[5][6] ;//! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F         * fhMCEAsymmetry[6] ;          //! E asymmetry of 2 splitted clusters vs cluster E
  TH2F         * fhMCAsymmetryLambda0[5][6] ; //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhMCAsymmetryDispEta[5][6] ; //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  TH2F         * fhMCAsymmetryDispPhi[5][6] ; //! E asymmetry of 2 splitted clusters vs lam0 for 5 E bins
  
  TH1F         * fhPtMCNo;                 //! Number of identified pi0, not coming from pi0/eta
  TH2F         * fhPhiMCNo;                //! Phi of identified pi0, not coming from pi0/eta
  TH2F         * fhEtaMCNo;                //! eta of identified  pi0, not coming from pi0/eta
  TH1F         * fhPtMC;                   //! Number of identified pi0, coming from pi0/eta
  TH2F         * fhPhiMC;                  //! Phi of identified pi0, coming from pi0/eta
  TH2F         * fhEtaMC;                  //! eta of identified pi0, coming from pi0/eta

  TH2F         * fhMassPairMCPi0;          //! pair mass, origin is same pi0
  TH2F         * fhMassPairMCEta;          //! pair mass, origin is same eta
  TH2F         * fhAnglePairMCPi0;         //! pair opening angle, origin is same pi0
  TH2F         * fhAnglePairMCEta;         //! pair opening angle, origin is same eta
  
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
  TH2F         * fhTrackMatchedMCParticle; //! Trace origin of matched particle
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

  AliAnaPi0EbE(              const AliAnaPi0EbE & pi0ebe) ; // cpy ctor
  AliAnaPi0EbE & operator = (const AliAnaPi0EbE & pi0ebe) ; // cpy assignment
  
  ClassDef(AliAnaPi0EbE,17)
} ;


#endif //ALIANAPI0EBE_H



