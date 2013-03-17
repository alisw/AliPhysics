#ifndef ALIANAINSIDECLUSTERINVARIANTMASS_H
#define ALIANAINSIDECLUSTERINVARIANTMASS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
//
// Split clusters with some criteria and calculate invariant mass
// to identify them as pi0 or conversion
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble)  
//_________________________________________________________________________


// --- ROOT system ---
class TList ;
class TObjString;
class TLorentzVector;

// --- ANALYSIS system ---
class AliAODCaloCluster;

#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaInsideClusterInvariantMass : public AliAnaCaloTrackCorrBaseClass {

 public: 
  
  AliAnaInsideClusterInvariantMass() ; // default ctor
  virtual ~AliAnaInsideClusterInvariantMass() { ; } //virtual dtor
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
    
  void         Init();
  
  void         InitParameters();
     
  void         MakeAnalysisFillHistograms() ; 
      
  void         Print(const Option_t * opt) const;

  void         SetCalorimeter(TString & det)             { fCalorimeter = det ; }
    
  void         SetM02Cut(Float_t min=0, Float_t max=10)  { fM02MinCut   = min ; fM02MaxCut  = max ; }

  void         SetMinNCells(Int_t cut)                   { fMinNCells   = cut ; }

  void         SetMinBadChannelDistance(Float_t cut)     { fMinBadDist  = cut ; }

  void         SwitchOnFillAngleHistograms()             { fFillAngleHisto      = kTRUE  ; }
  void         SwitchOffFillAngleHistograms()            { fFillAngleHisto      = kFALSE ; }
  
  void         SwitchOnFillExtraSSHistograms()           { fFillSSExtraHisto    = kTRUE  ; }
  void         SwitchOffFillExtraSSHistograms()          { fFillSSExtraHisto    = kFALSE ; }

  void         SwitchOnFillTMResidualHistograms()        { fFillTMResidualHisto = kTRUE  ; }
  void         SwitchOffFillTMResidualHistograms()       { fFillTMResidualHisto = kFALSE ; }
  
  void         SwitchOnMCFractionHistograms()            { fFillMCFractionHisto = kTRUE  ; }
  void         SwitchOffMCFractionHistograms()           { fFillMCFractionHisto = kFALSE ; }

  
  //For histograms
  enum mcTypes { kmcPhoton = 1, kmcConversion = 2, kmcPi0    = 3,  
                 kmcEta    = 4, kmcElectron   = 5, kmcHadron = 6, kmcPi0Conv = 7 };

 private:
  
  TString      fCalorimeter ;          // Calorimeter where the gamma is searched
  Float_t      fM02MaxCut   ;          // Study clusters with l0 smaller than cut
  Float_t      fM02MinCut   ;          // Study clusters with l0 larger than cut
  Int_t        fMinNCells   ;          // Study clusters with ncells larger than cut
  Float_t      fMinBadDist  ;          // Minimal distance to bad channel to accept cluster
  
  Bool_t       fFillAngleHisto;        // Fill splitted clusters angle histograms
  Bool_t       fFillTMResidualHisto ;  // Fill track matching histos, residuals
  Bool_t       fFillSSExtraHisto ;     // Fill shower shape extra histos
  Bool_t       fFillMCFractionHisto ;  // Fill MC energy fraction histos

  //Histograms
  
  TH2F       * fhMassNLocMax1[8][2]  ;                  //! Mass of 2 highest energy cells when 1 local max vs E, 1-6 for different MC particle types 
  TH2F       * fhMassNLocMax2[8][2]  ;                  //! Mass of 2 cells local maxima vs E,  1-6 for different MC particle types
  TH2F       * fhMassNLocMaxN[8][2]  ;                  //! Mass of >2 cells local maxima vs E, 1-6 for different MC particle types

  TH2F       * fhAsymNLocMax1[8][2]  ;                  //! Asymmetry of 2 highest energy cells when 1 local max vs E, 1-6 for different MC particle types 
  TH2F       * fhAsymNLocMax2[8][2]  ;                  //! Asymmetry of 2 cells local maxima vs E,  1-6 for different MC particle types
  TH2F       * fhAsymNLocMaxN[8][2]  ;                  //! Asymmetry of >2 cells local maxima vs E, 1-6 for different MC particle types
  
  TH2F       * fhSplitEFractionvsAsyNLocMax1[2] ;       //! sum of splitted cluster energy / cluster energy for N Local Maxima = 1 vs |A|
  TH2F       * fhSplitEFractionvsAsyNLocMax2[2] ;       //! sum of splitted cluster energy / cluster energy for N Local Maxima = 2 vs |A|
  TH2F       * fhSplitEFractionvsAsyNLocMaxN[2] ;       //! sum of splitted cluster energy / cluster energy for N Local Maxima > 2 vs |A|  
  
  TH2F       * fhMassM02CutNLocMax1  ;                  //! M02(E) selection, not matched, Mass of split clusters, NLM = 1
  TH2F       * fhMassM02CutNLocMax2  ;                  //! M02(E) selection, not matched, Mass of split clusters, NLM = 1
  TH2F       * fhMassM02CutNLocMaxN  ;                  //! M02(E) selection, not matched, Mass of split clusters, NLM > 2

  TH2F       * fhAsymM02CutNLocMax1  ;                  //! M02(E) selection, not matched, energy asymmetry of split clusters, NLM = 1
  TH2F       * fhAsymM02CutNLocMax2  ;                  //! M02(E) selection, not matched, energy asymmetry of split clusters, NLM = 2
  TH2F       * fhAsymM02CutNLocMaxN  ;                  //! M02(E) selection, not matched, energy asymmetry of split clusters, NLM > 2
  
  TH2F       * fhMassSplitECutNLocMax1 ;                //! 85% of split energy, not matched, Mass of split clusters, NLM = 1
  TH2F       * fhMassSplitECutNLocMax2 ;                //! 85% of split energy, not matched, Mass of split clusters, NLM = 1
  TH2F       * fhMassSplitECutNLocMaxN ;                //! 85% of split energy, not matched, Mass of split clusters, NLM > 2    
    
  TH2F       * fhMassM02NLocMax1[8][2]  ;               //! Mass of splitted clusters when 1  local max vs M02, for E > 8 GeV, 1-6 for different MC particle types
  TH2F       * fhMassM02NLocMax2[8][2]  ;               //! Mass of splitted clusters when 2  local max vs M02, for E > 8 GeV, 1-6 for different MC particle types
  TH2F       * fhMassM02NLocMaxN[8][2]  ;               //! Mass of splitted clusters when >2 local max vs M02, for E > 8 GeV, 1-6 for different MC particle types
  
  TH2F       * fhMassM02NLocMax1Ebin[4] ;               //! Mass of splitted clusters when 1  local max vs M02, 4 E bins, neutral clusters
  TH2F       * fhMassM02NLocMax2Ebin[4] ;               //! Mass of splitted clusters when 2  local max vs M02, 4 E bins, neutral clusters
  TH2F       * fhMassM02NLocMaxNEbin[4] ;               //! Mass of splitted clusters when >2 local max vs M02, 4 E bins, neutral clusters

  TH2F       * fhMassAsyNLocMax1Ebin[4] ;               //! Mass of Mass of splitted clusters when 1  local max vs asymmetry, 4 E bins, neutral clusters
  TH2F       * fhMassAsyNLocMax2Ebin[4] ;               //! Mass of Mass of splitted clusters when 2  local max vs asymmetry, 4 E bins, neutral clusters
  TH2F       * fhMassAsyNLocMaxNEbin[4] ;               //! Mass of Mass of splitted clusters when >2 local max vs asymmetry, 4 E bins, neutral clusters

  TH2F       * fhAsyMCGenRecoNLocMax1EbinPi0[4] ;       //! Generated vs reconstructed asymmetry of splitted clusters from pi0 when 1  local max, 4 E bins, neutral clusters
  TH2F       * fhAsyMCGenRecoNLocMax2EbinPi0[4] ;       //! Generated vs reconstructed asymmetry of splitted clusters from pi0 when 2  local max, 4 E bins, neutral clusters
  TH2F       * fhAsyMCGenRecoNLocMaxNEbinPi0[4] ;       //! Generated vs reconstructed asymmetry of splitted clusters from pi0 when >2 local max, 4 E bins, neutral clusters
  
  TH2F       * fhMassDispEtaNLocMax1[8][2]  ;           //! Mass of 2 highest energy cells when 1 local max, vs M02, for E > 8 GeV, 1-6 for different MC particle types 
  TH2F       * fhMassDispEtaNLocMax2[8][2]  ;           //! Mass of 2 cells local maxima, vs M02, for E > 8 GeV,  1-6 for different MC particle types
  TH2F       * fhMassDispEtaNLocMaxN[8][2]  ;           //! Mass of >2 cells local maxima, vs M02, for E > 8 GeV, 1-6 for different MC particle types  
  
  TH2F       * fhMassDispEtaNLocMax1Ebin[4] ;           //! Mass of 2 highest energy cells when 1 local max, vs M02, 4 E bins, neutral clusters 
  TH2F       * fhMassDispEtaNLocMax2Ebin[4] ;           //! Mass of 2 cells local maxima, vs M02, 4 E bins, neutral clusters
  TH2F       * fhMassDispEtaNLocMaxNEbin[4] ;           //! Mass of >2 cells local maxima, vs M02, 4 E bins, neutral clusters  
  
  TH2F       * fhMassDispPhiNLocMax1[8][2]  ;           //! Mass of 2 highest energy cells when 1 local max, vs M02, for E > 8 GeV, 1-6 for different MC particle types 
  TH2F       * fhMassDispPhiNLocMax2[8][2]  ;           //! Mass of 2 cells local maxima, vs M02, for E > 8 GeV,  1-6 for different MC particle types
  TH2F       * fhMassDispPhiNLocMaxN[8][2]  ;           //! Mass of >2 cells local maxima, vs M02, for E > 8 GeV, 1-6 for different MC particle types  
  
  TH2F       * fhMassDispPhiNLocMax1Ebin[4] ;           //! Mass of 2 highest energy cells when 1 local max, vs M02, 4 E bins, neutral clusters 
  TH2F       * fhMassDispPhiNLocMax2Ebin[4] ;           //! Mass of 2 cells local maxima, vs M02, 4 E bins, neutral clusters
  TH2F       * fhMassDispPhiNLocMaxNEbin[4] ;           //! Mass of >2 cells local maxima, vs M02, 4 E bins, neutral clusters  
  
  TH2F       * fhMassDispAsyNLocMax1[8][2]  ;           //! Mass of 2 highest energy cells when 1 local max, vs M02, for E > 8 GeV, 1-6 for different MC particle types 
  TH2F       * fhMassDispAsyNLocMax2[8][2]  ;           //! Mass of 2 cells local maxima, vs M02, for E > 8 GeV,  1-6 for different MC particle types
  TH2F       * fhMassDispAsyNLocMaxN[8][2]  ;           //! Mass of >2 cells local maxima, vs M02, for E > 8 GeV, 1-6 for different MC particle types  
  
  TH2F       * fhMassDispAsyNLocMax1Ebin[4] ;           //! Mass of 2 highest energy cells when 1 local max, vs M02, 4 E bins, neutral clusters 
  TH2F       * fhMassDispAsyNLocMax2Ebin[4] ;           //! Mass of 2 cells local maxima, vs M02, 4 E bins, neutral clusters
  TH2F       * fhMassDispAsyNLocMaxNEbin[4] ;           //! Mass of >2 cells local maxima, vs M02, 4 E bins, neutral clusters  
  
  TH2F       * fhNLocMax      [8][2] ;                  //! Number of maxima in cluster vs E, 1-6 for different MC particle types
  TH2F       * fhNLocMaxM02Cut[8][2] ;                  //! Number of maxima in cluster vs E, 1-6 for different MC particle types, after SS cut

  TH2F       * fhM02NLocMax1  [8][2] ;                  //! M02 vs E for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhM02NLocMax2  [8][2] ;                  //! M02 vs E for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhM02NLocMaxN  [8][2] ;                  //! M02 vs E for N max in cluster > 2, 1-6 for different MC particle types
  
  TH2F       * fhMCAsymM02NLocMax1MCPi0Ebin[4] ;        //! M02 vs decay asymmetry for N max in cluster = 1, for 4 energy bins
  TH2F       * fhMCAsymM02NLocMax2MCPi0Ebin[4] ;        //! M02 vs decay asymmetry for N max in cluster = 2, for 4 energy bins
  TH2F       * fhMCAsymM02NLocMaxNMCPi0Ebin[4] ;        //! M02 vs decay asymmetry for N max in cluster > 2, for 4 energy bins
  
  TH2F       * fhMCGenFracNLocMax1[8][2] ;              //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenFracNLocMax2[8][2] ;              //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenFracNLocMaxN[8][2] ;              //! E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types  
 
  TH2F       * fhMCGenFracAfterCutsNLocMax1MCPi0 ;      //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, MCPi0 after M02 and asymmetry cut
  TH2F       * fhMCGenFracAfterCutsNLocMax2MCPi0 ;      //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, MCPi0, after M02 and asymmetry cut
  TH2F       * fhMCGenFracAfterCutsNLocMaxNMCPi0 ;      //! E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, MCPi0, after M02 and asymmetry cut

  TH2F       * fhMCGenSplitEFracNLocMax1[8][2] ;        //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenSplitEFracNLocMax2[8][2] ;        //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenSplitEFracNLocMaxN[8][2] ;        //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types  

  TH2F       * fhMCGenSplitEFracAfterCutsNLocMax1MCPi0; //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenSplitEFracAfterCutsNLocMax2MCPi0; //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenSplitEFracAfterCutsNLocMaxNMCPi0; //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types
  
  TH2F       * fhMCGenEFracvsSplitEFracNLocMax1[8][2] ; //! E generated particle / E reconstructed vs E1+E2 reconstructed / E reconstructed for N max in cluster = 1, MC pi0
  TH2F       * fhMCGenEFracvsSplitEFracNLocMax2[8][2] ; //! E generated particle / E reconstructed vs E1+E2 reconstructed / E reconstructed for N max in cluster = 2, MC pi0
  TH2F       * fhMCGenEFracvsSplitEFracNLocMaxN[8][2] ; //! E generated particle / E reconstructed vs E1+E2 reconstructed / E reconstructed for N max in cluster > 2, MC pi0
  
  TH2F       * fhMCGenEvsSplitENLocMax1[8][2] ;         //! E generated particle vs E1+E2 for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenEvsSplitENLocMax2[8][2] ;         //! E generated particle vs E1+E2 for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenEvsSplitENLocMaxN[8][2] ;         //! E generated particle vs E1+E2 for N max in cluster > 2, 1-6 for different MC particle types  
  
  TH2F       * fhMCGenFracNLocMaxEbin[8][4] ;           //! NLM vs E generated particle / E reconstructed vs E reconstructed 1-6 for different MC particle types, not matched to track
  TH2F       * fhMCGenFracNLocMaxEbinMatched[8][4] ;    //! NLM vs E generated particle / E reconstructed vs E reconstructed 1-6 for different MC particle types, matched to track
  
  TH2F       * fhM02MCGenFracNLocMax1Ebin[8][4] ;       //! M02 vs E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types, not track matched
  TH2F       * fhM02MCGenFracNLocMax2Ebin[8][4] ;       //! M02 vs E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types, not track matched
  TH2F       * fhM02MCGenFracNLocMaxNEbin[8][4] ;       //! M02 vs E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types, not track matched  
  
  TH2F       * fhMassMCGenFracNLocMax1Ebin[8][4] ;      //! Mass vs E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types, not track matched
  TH2F       * fhMassMCGenFracNLocMax2Ebin[8][4] ;      //! Mass vs E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types, not track matched
  TH2F       * fhMassMCGenFracNLocMaxNEbin[8][4] ;      //! Mass vs E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types, not track matched  
  
  TH2F       * fhNCellNLocMax1[8][2] ;                  //! n cells in cluster vs E for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhNCellNLocMax2[8][2] ;                  //! n cells in cluster vs E for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhNCellNLocMaxN[8][2] ;                  //! n cells in cluster vs E for N max in cluster > 2, 1-6 for different MC particle types
  
  TH2F       * fhM02Pi0NLocMax1[8][2] ;                 //! M02 for Mass around pi0, N Local Maxima = 1
  TH2F       * fhM02EtaNLocMax1[8][2] ;                 //! M02 for Mass around eta, N Local Maxima = 1
  TH2F       * fhM02ConNLocMax1[8][2] ;                 //! M02 for Mass around close to 0, N Local Maxima = 1
  
  TH2F       * fhM02Pi0NLocMax2[8][2] ;                 //! M02 for Mass around pi0, N Local Maxima = 2
  TH2F       * fhM02EtaNLocMax2[8][2] ;                 //! M02 for Mass around eta, N Local Maxima = 2
  TH2F       * fhM02ConNLocMax2[8][2] ;                 //! M02 for Mass around close to 0, N Local Maxima = 2
  
  TH2F       * fhM02Pi0NLocMaxN[8][2] ;                 //! M02 for Mass around pi0, N Local Maxima > 2
  TH2F       * fhM02EtaNLocMaxN[8][2] ;                 //! M02 for Mass around eta, N Local Maxima > 2
  TH2F       * fhM02ConNLocMaxN[8][2] ;                 //! M02 for Mass around close to 0, N Local Maxima > 2

  TH2F       * fhMassPi0NLocMax1[8][2] ;                //! Mass for selected pi0, N Local Maxima = 1
  TH2F       * fhMassEtaNLocMax1[8][2] ;                //! Mass for selected around eta, N Local Maxima = 1
  TH2F       * fhMassConNLocMax1[8][2] ;                //! Mass for selected around close to 0, N Local Maxima = 1
  
  TH2F       * fhMassPi0NLocMax2[8][2] ;                //! Mass for selected around pi0, N Local Maxima = 2
  TH2F       * fhMassEtaNLocMax2[8][2] ;                //! Mass for selected around eta, N Local Maxima = 2
  TH2F       * fhMassConNLocMax2[8][2] ;                //! Mass for selected around close to 0, N Local Maxima = 2
  
  TH2F       * fhMassPi0NLocMaxN[8][2] ;                //! Mass for selected around pi0, N Local Maxima > 2
  TH2F       * fhMassEtaNLocMaxN[8][2] ;                //! Mass for selected around eta, N Local Maxima > 2
  TH2F       * fhMassConNLocMaxN[8][2] ;                //! Mass for selected around close to 0, N Local Maxima > 2
  
  TH2F       * fhMassAfterCutsNLocMax1[8][2] ;          //! Mass after M02, asymmetry cuts for pi0, N Local Maxima = 1
  TH2F       * fhMassAfterCutsNLocMax2[8][2] ;          //! Mass after M02, asymmetry cuts for pi0, N Local Maxima = 2
  TH2F       * fhMassAfterCutsNLocMaxN[8][2] ;          //! Mass after M02, asymmetry cuts for pi0, N Local Maxima > 2
  
  TH2F       * fhAsyPi0NLocMax1[8][2] ;                 //! Asy for Mass around pi0, N Local Maxima = 1
  TH2F       * fhAsyEtaNLocMax1[8][2] ;                 //! Asy for Mass around eta, N Local Maxima = 1
  TH2F       * fhAsyConNLocMax1[8][2] ;                 //! Asy for Mass around close to 0, N Local Maxima = 1
  
  TH2F       * fhAsyPi0NLocMax2[8][2] ;                 //! Asy for Mass around pi0, N Local Maxima = 2
  TH2F       * fhAsyEtaNLocMax2[8][2] ;                 //! Asy for Mass around eta, N Local Maxima = 2
  TH2F       * fhAsyConNLocMax2[8][2] ;                 //! Asy for Mass around close to 0, N Local Maxima = 2
  
  TH2F       * fhAsyPi0NLocMaxN[8][2] ;                 //! Asy for Mass around pi0, N Local Maxima > 2
  TH2F       * fhAsyEtaNLocMaxN[8][2] ;                 //! Asy for Mass around eta, N Local Maxima > 2
  TH2F       * fhAsyConNLocMaxN[8][2] ;                 //! Asy for Mass around close to 0, N Local Maxima > 2
  
  TH2F       * fhSplitEFractionNLocMax1[8][2] ;         //! sum of splitted cluster energy / cluster energy for N Local Maxima = 1
  TH2F       * fhSplitEFractionNLocMax2[8][2] ;         //! sum of splitted cluster energy / cluster energy for N Local Maxima = 2
  TH2F       * fhSplitEFractionNLocMaxN[8][2] ;         //! sum of splitted cluster energy / cluster energy for N Local Maxima > 2

  TH2F       * fhSplitEFractionAfterCutsNLocMax1[8][2] ; //! sum of splitted cluster energy / cluster energy for N Local Maxima = 1, after M02 and asymmetry cut
  TH2F       * fhSplitEFractionAfterCutsNLocMax2[8][2] ; //! sum of splitted cluster energy / cluster energy for N Local Maxima = 2, after M02 and asymmetry cut
  TH2F       * fhSplitEFractionAfterCutsNLocMaxN[8][2] ; //! sum of splitted cluster energy / cluster energy for N Local Maxima > 2, after M02 and asymmetry cut
  
  TH2F       * fhMassSplitEFractionNLocMax1Ebin[8][4] ; //! Mass vs sum of splitted cluster energy / cluster energy for N max in cluster = 1, 1-6 for different MC particle types, not track matched
  TH2F       * fhMassSplitEFractionNLocMax2Ebin[8][4] ; //! Mass vs sum of splitted cluster energy / cluster energy for N max in cluster = 2, 1-6 for different MC particle types, not track matched
  TH2F       * fhMassSplitEFractionNLocMaxNEbin[8][4] ; //! Mass vs sum of splitted cluster energy / cluster energy for N max in cluster > 2, 1-6 for different MC particle types, not track matched  
    
  TH2F       * fhAnglePairNLocMax1[2] ;                 //! pair opening angle vs E
  TH2F       * fhAnglePairNLocMax2[2] ;                 //! pair opening angle vs E
  TH2F       * fhAnglePairNLocMaxN[2] ;                 //! pair opening angle vs E

  TH2F       * fhAnglePairMassNLocMax1[2] ;             //! pair opening angle vs Mass for E > 7 GeV
  TH2F       * fhAnglePairMassNLocMax2[2] ;             //! pair opening angle vs Mass for E > 7 GeV
  TH2F       * fhAnglePairMassNLocMaxN[2] ;             //! pair opening angle vs Mass for E > 7 GeV
  
  TH2F       * fhTrackMatchedDEtaNLocMax1[8] ;          //! Eta distance between track and cluster vs cluster E, 1 local maximum
  TH2F       * fhTrackMatchedDPhiNLocMax1[8] ;          //! Phi distance between track and cluster vs cluster E, 1 local maximum
  TH2F       * fhTrackMatchedDEtaNLocMax2[8] ;          //! Eta distance between track and cluster vs cluster E, 2 local maximum
  TH2F       * fhTrackMatchedDPhiNLocMax2[8] ;          //! Phi distance between track and cluster vs cluster E, 2 local maximum
  TH2F       * fhTrackMatchedDEtaNLocMaxN[8] ;          //! Eta distance between track and cluster vs cluster E, more than 2 local maximum
  TH2F       * fhTrackMatchedDPhiNLocMaxN[8] ;          //! Phi distance between track and cluster vs cluster E, more than 2 local maximum
  
  TH2F       * fhCentralityPi0NLocMax1[8][2] ;          //! Centrality for selected pi0, N Local Maxima = 1
  TH2F       * fhCentralityEtaNLocMax1[8][2] ;          //! Centrality for selected eta, N Local Maxima = 1
  TH2F       * fhCentralityPi0NLocMax2[8][2] ;          //! Centrality for selected pi0, N Local Maxima = 2
  TH2F       * fhCentralityEtaNLocMax2[8][2] ;          //! Centrality for selected eta, N Local Maxima = 2
  TH2F       * fhCentralityPi0NLocMaxN[8][2] ;          //! Centrality for selected pi0, N Local Maxima > 2
  TH2F       * fhCentralityEtaNLocMaxN[8][2] ;          //! Centrality for selected eta, N Local Maxima > 2

  TH2F       * fhEventPlanePi0NLocMax1 ;                //! Event plane for selected pi0, N Local Maxima = 1
  TH2F       * fhEventPlaneEtaNLocMax1 ;                //! Event plane for selected eta, N Local Maxima = 1
  TH2F       * fhEventPlanePi0NLocMax2 ;                //! Event plane for selected pi0, N Local Maxima = 2
  TH2F       * fhEventPlaneEtaNLocMax2 ;                //! Event plane for selected eta, N Local Maxima = 2
  TH2F       * fhEventPlanePi0NLocMaxN ;                //! Event plane for selected pi0, N Local Maxima > 2
  TH2F       * fhEventPlaneEtaNLocMaxN ;                //! Event plane for selected eta, N Local Maxima > 2

  
  AliAnaInsideClusterInvariantMass(              const AliAnaInsideClusterInvariantMass & split) ; // cpy ctor
  AliAnaInsideClusterInvariantMass & operator = (const AliAnaInsideClusterInvariantMass & split) ; // cpy assignment
  
  ClassDef(AliAnaInsideClusterInvariantMass,19)
  
} ;

#endif //ALIANAINSIDECLUSTERINVARIANTMASS_H



