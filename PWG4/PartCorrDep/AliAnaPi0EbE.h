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
  
  void           FillSelectedClusterHistograms(AliVCluster* cluster, const Int_t tag);
    
  void           FillWeightHistograms(AliVCluster *clus);
    
  void           MakeInvMassInCalorimeter() ;
  
  void           MakeInvMassInCalorimeterAndCTS() ;
  
  void           MakeShowerShapeIdentification() ;
  
  void           RecalibrateCellAmplitude(Float_t  & amp,  const Int_t absId);
      
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

  void           SwitchOnFillWeightHistograms()              { fFillWeightHistograms = kTRUE  ; }
  void           SwitchOffFillWeightHistograms()             { fFillWeightHistograms = kFALSE ; }  
  
  void           SwitchOnTMHistoFill()                       { fFillTMHisto          = kTRUE  ; }
  void           SwitchOffTMHistoFill()                      { fFillTMHisto          = kFALSE ; }

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
  
  Bool_t         fFillWeightHistograms ;   // Fill weigth histograms
  Bool_t         fFillTMHisto;             // Fill track matching plots

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

  //MC histograms
  
  TH2F         * fhEMCLambda0[6] ;         //! E vs lambda0 of pi0 pairs but really from MC particle
  TH2F         * fhEMCLambda1[6] ;         //! E vs lambda1 of pi0 pairs but really from MC particle
  TH2F         * fhEMCDispersion[6] ;      //! E vs dispersion of pi0 pairs but really from MC particle
  TH2F         * fhEMCLambda0NoTRD[6] ;         //! E vs lambda0 of pi0 pairs but really from MC particle, not behind TRD
  TH2F         * fhEMCLambda0FracMaxCellCut[6] ;//! E vs lambda0 of pi0 pairs but really from MC particle, fraction of cluster energy in max cell cut
  TH2F         * fhEMCFracMaxCell[6] ;     //! E vs fraction of max cell 
  
  TH1F         * fhPtMCNo;                 //! Number of identified pi0, not coming from pi0/eta
  TH2F         * fhPhiMCNo;                //! Phi of identified pi0, not coming from pi0/eta
  TH2F         * fhEtaMCNo;                //! eta of identified  pi0, not coming from pi0/eta
  TH1F         * fhPtMC;                   //! Number of identified pi0, coming from pi0/eta
  TH2F         * fhPhiMC;                  //! Phi of identified pi0, coming from pi0/eta
  TH2F         * fhEtaMC;                  //! eta of identified pi0, coming from pi0/eta
  
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
  
  
  AliAnaPi0EbE(              const AliAnaPi0EbE & g) ; // cpy ctor
  AliAnaPi0EbE & operator = (const AliAnaPi0EbE & g) ; // cpy assignment
  
  ClassDef(AliAnaPi0EbE,11)
} ;


#endif //ALIANAPI0EBE_H



