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
                 kmcEta    = 4, kmcElectron   = 5, kmcHadron = 6 };

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
  
  TH2F       * fhMassNLocMax1[7][2]  ; //! Mass of 2 highest energy cells when 1 local max vs E, 1-6 for different MC particle types 
  TH2F       * fhMassNLocMax2[7][2]  ; //! Mass of 2 cells local maxima vs E,  1-6 for different MC particle types
  TH2F       * fhMassNLocMaxN[7][2]  ; //! Mass of >2 cells local maxima vs E, 1-6 for different MC particle types

  TH2F       * fhAsymNLocMax1[2]     ; //! Asymmetry of 2 highest energy cells when 1 local max vs E, 1-6 for different MC particle types 
  TH2F       * fhAsymNLocMax2[2]     ; //! Asymmetry of 2 cells local maxima vs E,  1-6 for different MC particle types
  TH2F       * fhAsymNLocMaxN[2]     ; //! Asymmetry of >2 cells local maxima vs E, 1-6 for different MC particle types
  
  TH2F       * fhMassM02CutNLocMax1  ; //! M02(E) selection, not matched, Mass of 2 highest energy cells when 1 local max vs E, 1-6 for different MC particle types 
  TH2F       * fhMassM02CutNLocMax2  ; //! M02(E) selection, not matched, Mass of 2 cells local maxima vs E,  1-6 for different MC particle types
  TH2F       * fhMassM02CutNLocMaxN  ; //! M02(E) selection, not matched, Mass of >2 cells local maxima vs E, 1-6 for different MC particle types  

  TH2F       * fhMassSplitECutNLocMax1 ; //! 85% of split energy, not matched, Mass of 2 highest energy cells when 1 local max vs E, 1-6 for different MC particle types 
  TH2F       * fhMassSplitECutNLocMax2 ; //! 85% of split energy, not matched, Mass of 2 cells local maxima vs E,  1-6 for different MC particle types
  TH2F       * fhMassSplitECutNLocMaxN ; //! 85% of split energy, not matched, Mass of >2 cells local maxima vs E, 1-6 for different MC particle types    
  
  TH2F       * fhMassM02NLocMax1[7][2]  ; //! Mass of 2 highest energy cells when 1 local max, vs M02, for E > 7 GeV, 1-6 for different MC particle types 
  TH2F       * fhMassM02NLocMax2[7][2]  ; //! Mass of 2 cells local maxima, vs M02, for E > 7 GeV,  1-6 for different MC particle types
  TH2F       * fhMassM02NLocMaxN[7][2]  ; //! Mass of >2 cells local maxima, vs M02, for E > 7 GeV, 1-6 for different MC particle types  
  
  TH2F       * fhMassM02NLocMax1Ebin[4] ; //! Mass of 2 highest energy cells when 1 local max, vs M02, 4 E bins, neutral clusters 
  TH2F       * fhMassM02NLocMax2Ebin[4] ; //! Mass of 2 cells local maxima, vs M02, 4 E bins, neutral clusters
  TH2F       * fhMassM02NLocMaxNEbin[4] ; //! Mass of >2 cells local maxima, vs M02, 4 E bins, neutral clusters  
  
  TH2F       * fhMassDispEtaNLocMax1[7][2]  ; //! Mass of 2 highest energy cells when 1 local max, vs M02, for E > 7 GeV, 1-6 for different MC particle types 
  TH2F       * fhMassDispEtaNLocMax2[7][2]  ; //! Mass of 2 cells local maxima, vs M02, for E > 7 GeV,  1-6 for different MC particle types
  TH2F       * fhMassDispEtaNLocMaxN[7][2]  ; //! Mass of >2 cells local maxima, vs M02, for E > 7 GeV, 1-6 for different MC particle types  
  
  TH2F       * fhMassDispEtaNLocMax1Ebin[4] ; //! Mass of 2 highest energy cells when 1 local max, vs M02, 4 E bins, neutral clusters 
  TH2F       * fhMassDispEtaNLocMax2Ebin[4] ; //! Mass of 2 cells local maxima, vs M02, 4 E bins, neutral clusters
  TH2F       * fhMassDispEtaNLocMaxNEbin[4] ; //! Mass of >2 cells local maxima, vs M02, 4 E bins, neutral clusters  
  
  TH2F       * fhMassDispPhiNLocMax1[7][2]  ; //! Mass of 2 highest energy cells when 1 local max, vs M02, for E > 7 GeV, 1-6 for different MC particle types 
  TH2F       * fhMassDispPhiNLocMax2[7][2]  ; //! Mass of 2 cells local maxima, vs M02, for E > 7 GeV,  1-6 for different MC particle types
  TH2F       * fhMassDispPhiNLocMaxN[7][2]  ; //! Mass of >2 cells local maxima, vs M02, for E > 7 GeV, 1-6 for different MC particle types  
  
  TH2F       * fhMassDispPhiNLocMax1Ebin[4] ; //! Mass of 2 highest energy cells when 1 local max, vs M02, 4 E bins, neutral clusters 
  TH2F       * fhMassDispPhiNLocMax2Ebin[4] ; //! Mass of 2 cells local maxima, vs M02, 4 E bins, neutral clusters
  TH2F       * fhMassDispPhiNLocMaxNEbin[4] ; //! Mass of >2 cells local maxima, vs M02, 4 E bins, neutral clusters  
  
  TH2F       * fhMassDispAsyNLocMax1[7][2]  ; //! Mass of 2 highest energy cells when 1 local max, vs M02, for E > 7 GeV, 1-6 for different MC particle types 
  TH2F       * fhMassDispAsyNLocMax2[7][2]  ; //! Mass of 2 cells local maxima, vs M02, for E > 7 GeV,  1-6 for different MC particle types
  TH2F       * fhMassDispAsyNLocMaxN[7][2]  ; //! Mass of >2 cells local maxima, vs M02, for E > 7 GeV, 1-6 for different MC particle types  
  
  TH2F       * fhMassDispAsyNLocMax1Ebin[4] ; //! Mass of 2 highest energy cells when 1 local max, vs M02, 4 E bins, neutral clusters 
  TH2F       * fhMassDispAsyNLocMax2Ebin[4] ; //! Mass of 2 cells local maxima, vs M02, 4 E bins, neutral clusters
  TH2F       * fhMassDispAsyNLocMaxNEbin[4] ; //! Mass of >2 cells local maxima, vs M02, 4 E bins, neutral clusters  
  
  TH2F       * fhNLocMax      [7][2] ; //! Number of maxima in cluster vs E, 1-6 for different MC particle types
  TH2F       * fhNLocMaxM02Cut[7][2] ; //! Number of maxima in cluster vs E, 1-6 for different MC particle types, after SS cut

  TH2F       * fhM02NLocMax1  [7][2] ; //! M02 vs E for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhM02NLocMax2  [7][2] ; //! M02 vs E for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhM02NLocMaxN  [7][2] ; //! M02 vs E for N max in cluster > 2, 1-6 for different MC particle types
  
  TH2F       * fhMCAsymM02NLocMax1MCPi0Ebin[4] ; //! M02 vs decay asymmetry for N max in cluster = 1, for 4 energy bins
  TH2F       * fhMCAsymM02NLocMax2MCPi0Ebin[4] ; //! M02 vs decay asymmetry for N max in cluster = 2, for 4 energy bins
  TH2F       * fhMCAsymM02NLocMaxNMCPi0Ebin[4] ; //! M02 vs decay asymmetry for N max in cluster > 2, for 4 energy bins
  
  TH2F       * fhMCGenFracNLocMax1[7][2] ; //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenFracNLocMax2[7][2] ; //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenFracNLocMaxN[7][2] ; //! E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types  
  
  TH2F       * fhMCGenSplitEFracNLocMax1[7][2] ; //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenSplitEFracNLocMax2[7][2] ; //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenSplitEFracNLocMaxN[7][2] ; //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types  
  
  TH2F       * fhMCGenEFracvsSplitEFracNLocMax1[7][2] ; //! E generated particle / E reconstructed vs E1+E2 reconstructed / E reconstructed for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenEFracvsSplitEFracNLocMax2[7][2] ; //! E generated particle / E reconstructed vs E1+E2 reconstructed / E reconstructed for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenEFracvsSplitEFracNLocMaxN[7][2] ; //! E generated particle / E reconstructed vs E1+E2 reconstructed / E reconstructed for N max in cluster > 2, 1-6 for different MC particle types  
  
  TH2F       * fhMCGenEvsSplitENLocMax1[7][2] ; //! E generated particle vs E1+E2 for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenEvsSplitENLocMax2[7][2] ; //! E generated particle vs E1+E2 for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenEvsSplitENLocMaxN[7][2] ; //! E generated particle vs E1+E2 for N max in cluster > 2, 1-6 for different MC particle types  
  
  TH2F       * fhMCGenFracNLocMaxEbin[7][4] ;    //! NLM vs E generated particle / E reconstructed vs E reconstructed 1-6 for different MC particle types, not matched to track
  TH2F       * fhMCGenFracNLocMaxEbinMatched[7][4] ; //! NLM vs E generated particle / E reconstructed vs E reconstructed 1-6 for different MC particle types, matched to track
  
  TH2F       * fhM02MCGenFracNLocMax1Ebin[7][4] ; //! M02 vs E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types, not track matched
  TH2F       * fhM02MCGenFracNLocMax2Ebin[7][4] ; //! M02 vs E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types, not track matched
  TH2F       * fhM02MCGenFracNLocMaxNEbin[7][4] ; //! M02 vs E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types, not track matched  
  
  TH2F       * fhMassMCGenFracNLocMax1Ebin[7][4] ; //! Mass vs E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types, not track matched
  TH2F       * fhMassMCGenFracNLocMax2Ebin[7][4] ; //! Mass vs E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types, not track matched
  TH2F       * fhMassMCGenFracNLocMaxNEbin[7][4] ; //! Mass vs E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types, not track matched  
  
  TH2F       * fhNCellNLocMax1[7][2] ; //! n cells in cluster vs E for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhNCellNLocMax2[7][2] ; //! n cells in cluster vs E for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhNCellNLocMaxN[7][2] ; //! n cells in cluster vs E for N max in cluster > 2, 1-6 for different MC particle types
  
  TH2F       * fhM02Pi0LocMax1[7][2] ; //! M02 for Mass around pi0, N Local Maxima = 1
  TH2F       * fhM02EtaLocMax1[7][2] ; //! M02 for Mass around eta, N Local Maxima = 1
  TH2F       * fhM02ConLocMax1[7][2] ; //! M02 for Mass around close to 0, N Local Maxima = 1

  TH2F       * fhM02Pi0LocMax2[7][2] ; //! M02 for Mass around pi0, N Local Maxima = 2
  TH2F       * fhM02EtaLocMax2[7][2] ; //! M02 for Mass around eta, N Local Maxima = 2
  TH2F       * fhM02ConLocMax2[7][2] ; //! M02 for Mass around close to 0, N Local Maxima = 2
  
  TH2F       * fhM02Pi0LocMaxN[7][2] ; //! M02 for Mass around pi0, N Local Maxima > 2
  TH2F       * fhM02EtaLocMaxN[7][2] ; //! M02 for Mass around eta, N Local Maxima > 2
  TH2F       * fhM02ConLocMaxN[7][2] ; //! M02 for Mass around close to 0, N Local Maxima > 2

  TH2F       * fhMassPi0LocMax1[7][2] ; //! Mass for selected pi0, N Local Maxima = 1
  TH2F       * fhMassEtaLocMax1[7][2] ; //! Mass for selected around eta, N Local Maxima = 1
  TH2F       * fhMassConLocMax1[7][2] ; //! Mass for selected around close to 0, N Local Maxima = 1
  
  TH2F       * fhMassPi0LocMax2[7][2] ; //! Mass for selected around pi0, N Local Maxima = 2
  TH2F       * fhMassEtaLocMax2[7][2] ; //! Mass for selected around eta, N Local Maxima = 2
  TH2F       * fhMassConLocMax2[7][2] ; //! Mass for selected around close to 0, N Local Maxima = 2
  
  TH2F       * fhMassPi0LocMaxN[7][2] ; //! Mass for selected around pi0, N Local Maxima > 2
  TH2F       * fhMassEtaLocMaxN[7][2] ; //! Mass for selected around eta, N Local Maxima > 2
  TH2F       * fhMassConLocMaxN[7][2] ; //! Mass for selected around close to 0, N Local Maxima > 2
  
  TH2F       * fhSplitEFractionNLocMax1[7][2] ; //! sum of splitted cluster energy / cluster energy for N Local Maxima = 1
  TH2F       * fhSplitEFractionNLocMax2[7][2] ; //! sum of splitted cluster energy / cluster energy for N Local Maxima = 2
  TH2F       * fhSplitEFractionNLocMaxN[7][2] ; //! sum of splitted cluster energy / cluster energy for N Local Maxima > 2
  
  TH2F       * fhMassSplitEFractionNLocMax1Ebin[7][4] ; //! Mass vs sum of splitted cluster energy / cluster energy for N max in cluster = 1, 1-6 for different MC particle types, not track matched
  TH2F       * fhMassSplitEFractionNLocMax2Ebin[7][4] ; //! Mass vs sum of splitted cluster energy / cluster energy for N max in cluster = 2, 1-6 for different MC particle types, not track matched
  TH2F       * fhMassSplitEFractionNLocMaxNEbin[7][4] ; //! Mass vs sum of splitted cluster energy / cluster energy for N max in cluster > 2, 1-6 for different MC particle types, not track matched  
    
  TH2F       * fhAnglePairLocMax1[2] ; //! pair opening angle vs E
  TH2F       * fhAnglePairLocMax2[2] ; //! pair opening angle vs E
  TH2F       * fhAnglePairLocMaxN[2] ; //! pair opening angle vs E
  
  TH2F       * fhAnglePairMassLocMax1[2] ; //! pair opening angle vs Mass for E > 7 GeV
  TH2F       * fhAnglePairMassLocMax2[2] ; //! pair opening angle vs Mass for E > 7 GeV
  TH2F       * fhAnglePairMassLocMaxN[2] ; //! pair opening angle vs Mass for E > 7 GeV
  
  TH2F       * fhTrackMatchedDEtaLocMax1[7] ; //! Eta distance between track and cluster vs cluster E, 1 local maximum
  TH2F       * fhTrackMatchedDPhiLocMax1[7] ; //! Phi distance between track and cluster vs cluster E, 1 local maximum
  TH2F       * fhTrackMatchedDEtaLocMax2[7] ; //! Eta distance between track and cluster vs cluster E, 2 local maximum
  TH2F       * fhTrackMatchedDPhiLocMax2[7] ; //! Phi distance between track and cluster vs cluster E, 2 local maximum
  TH2F       * fhTrackMatchedDEtaLocMaxN[7] ; //! Eta distance between track and cluster vs cluster E, more than 2 local maximum
  TH2F       * fhTrackMatchedDPhiLocMaxN[7] ; //! Phi distance between track and cluster vs cluster E, more than 2 local maximum
  
  AliAnaInsideClusterInvariantMass(              const AliAnaInsideClusterInvariantMass & split) ; // cpy ctor
  AliAnaInsideClusterInvariantMass & operator = (const AliAnaInsideClusterInvariantMass & split) ; // cpy assignment
  
  ClassDef(AliAnaInsideClusterInvariantMass,16)
  
} ;

#endif //ALIANAINSIDECLUSTERINVARIANTMASS_H



