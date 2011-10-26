#ifndef ALIANAPI0EBE_H
#define ALIANAPI0EBE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnaPi0EbE.h 27413 2008-07-18 13:28:12Z gconesab $ */

//_________________________________________________________________________
//
// Class for the analysis of high pT pi0 event by event
// Pi0 identified by one of the following:
//  -Invariant mass of 2 cluster in calorimeter
//  -Shower shape analysis in calorimeter
//  -Invariant mass of one cluster in calorimeter and one photon reconstructed in TPC (in near future)
//
//-- Author: Gustavo Conesa (INFN-LNF)  &  Raphaelle Ichou (SUBATECH)
//_________________________________________________________________________


// --- ROOT system ---
class TH3F ; 
class TList ;
class TObjString;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"

class AliAnaPi0EbE : public AliAnaPartCorrBaseClass {

 public: 
  AliAnaPi0EbE() ; // default ctor
  virtual ~AliAnaPi0EbE() { ; } //virtual dtor
 private:
  AliAnaPi0EbE(const AliAnaPi0EbE & g) ; // cpy ctor
  AliAnaPi0EbE & operator = (const AliAnaPi0EbE & g) ;//cpy assignment

 public:  
	
  //General
  
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
    
  void           SwitchOnFillWeightHistograms()              { fFillWeightHistograms = kTRUE  ; }
  void           SwitchOffFillWeightHistograms()             { fFillWeightHistograms = kFALSE ; }  
  
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

  //For histograms
  enum mcTypes   { mcPhoton = 0, mcConversion = 1, mcPi0    = 2,  
                   mcEta    = 3, mcElectron   = 4, mcHadron = 5 };

 private:
  
  anaTypes       fAnaType; //Select analysis type
  
  //Only for pi0 SS identification case, kSSCalo
  TString        fCalorimeter ;            // Calorimeter where the gamma is searched;
  Float_t        fMinDist ;                // Minimal distance to bad channel to accept cluster
  Float_t        fMinDist2;                // Cuts on Minimal distance to study acceptance evaluation
  Float_t        fMinDist3;                // One more cut on distance used for acceptance-efficiency study
  
  Bool_t         fFillWeightHistograms ;   // Fill weigth histograms
  
  //Only for combination of calorimeter and conversion photons, kIMCaloTracks
  TString        fInputAODGammaConvName;   //  Name of AOD branch with conversion photons
  
  //Histograms
  
  TH1F         * fhPtPi0  ;                //! Number of identified  pi0 vs pT
  TH1F         * fhEPi0   ;                //! Number of identified  pi0 vs E
  TH3F         * fhEEtaPhiPi0  ;           //! E vs eta phi of identified  pi0 
  
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
  
  TH1F         * fhPtMCNoPi0;              //! Number of identified pi0, not coming from pi0
  TH2F         * fhPhiMCNoPi0;             //! Phi of identified pi0, not coming from pi0
  TH2F         * fhEtaMCNoPi0;             //! eta of identified  pi0, not coming from pi0
  TH1F         * fhPtMCPi0;                //! Number of identified pi0, coming from pi0
  TH2F         * fhPhiMCPi0;               //! Phi of identified pi0, coming from pi0
  TH2F         * fhEtaMCPi0;               //! eta of identified pi0, coming from pi0
  
  // Weight studies
  
  TH2F         * fhECellClusterRatio;      //! e cell / e cluster vs e cluster for selected photons
  TH2F         * fhECellClusterLogRatio;   //! log (e cell / e cluster)  vs e cluster for selected photons
  TH2F         * fhEMaxCellClusterRatio;   //! e max cell / e cluster vs e cluster for selected photons
  TH2F         * fhEMaxCellClusterLogRatio;//! log (e max cell / e cluster) vs e cluster for selected photons
  TH2F         * fhLambda0ForW0[14];        //! L0 for 7 defined w0= 3, 3.5 ... 6 for selected photons
  //TH2F         * fhLambda1ForW0[7];        //! L1 for 7 defined w0= 3, 3.5 ... 6 for selected photons  
  
  ClassDef(AliAnaPi0EbE,10)
} ;


#endif //ALIANAPI0EBE_H



