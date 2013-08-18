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
  
  void         CheckLocalMaximaMCOrigin(AliVCluster* cluster, const Int_t mcindex, const Int_t noverlaps,
                                        const Float_t e1,     const Float_t e2,    const Float_t mass);
                                        //, Float_t m02, TLorentzVector l1, TLorentzVector l2);
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  void         GetMCIndex(AliVCluster * cluster, Int_t & mcindex, Int_t & mcTag);
  
  void         GetMCPrimaryKine(AliVCluster* cluster, const Int_t mcindex, const Int_t mcTag, const Bool_t matched,
                                Float_t & eprim, Float_t & asymGen, Int_t & noverlaps );
  
  void         FillAngleHistograms(const Int_t nMax, const Bool_t matched,
                                   const Float_t en, const Float_t angle, const Float_t mass);
  
  void         FillEBinHistograms(const Int_t ebin, const Int_t nMax, const Int_t mcindex, const Float_t splitFrac,
                                  const Float_t mass, const Float_t asym, const Float_t l0);
  
  void         FillMCHistograms(const Float_t en,        const Float_t e1  , const Float_t e2,
                                const Int_t ebin,        const Int_t mcindex,const Int_t noverlaps,
                                const Float_t l0,        const Float_t mass,
                                const Int_t nMax,        const Bool_t  matched,
                                const Float_t splitFrac, const Float_t asym,
                                const Float_t eprim,     const Float_t asymGen);
  
  void         FillMCOverlapHistograms(const Float_t en,      const Float_t enprim,
                                       const Int_t   nc,      const Float_t mass,    const Float_t l0,
                                       const Float_t asym,    const Float_t splitFrac,
                                       const Int_t   nlm,     const Int_t ebin,   const Bool_t matched,
                                       const Int_t   mcindex, const Int_t noverlaps);
  
  void         FillSSWeightHistograms(AliVCluster *cluster, const Int_t nlm, const Int_t absId1, const Int_t absId2);
  
  void         FillSSExtraHistograms(AliVCluster *cluster, const Int_t nMax,
                                     const Bool_t matched, const Int_t mcindex,
                                     const Float_t mass  , const Int_t ebin);

  void         FillNCellHistograms(const Int_t   ncells,  const Float_t energy, const Int_t nMax,
                                   const Bool_t  matched, const Int_t mcindex,
                                   const Float_t mass   , const Float_t l0);
  
  void         FillTrackMatchingHistograms(AliVCluster * cluster,const Int_t nMax, const Int_t mcindex);
  
  void         FillHistograms1(const Float_t en,     const Float_t e1,     const Float_t e2,
                               const Int_t nMax,     const Float_t mass,   const Float_t l0,
                               const Float_t eta,    const Float_t phi,
                               const Bool_t matched, const Int_t mcindex);

  
  void         FillHistograms2(const Float_t en,     const Float_t eprim,
                               const Float_t e1,     const Float_t e2,      const Int_t nMax,  
                               const Float_t mass,   const Float_t l0,
                               const Bool_t matched, const Int_t mcindex);
  
  void         FillIdPi0Histograms(const Float_t en,     const Float_t e1,  const Float_t e2,
                                   const Int_t nc,       const Int_t nMax,  const Float_t t12diff,
                                   const Float_t mass,   const Float_t l0,
                                   const Float_t eta,    const Float_t phi,
                                   const Bool_t matched, const Int_t mcindex);
  
  void         FillIdEtaHistograms(const Float_t en,     const Float_t e1,  const Float_t e2,
                                   const Int_t nc,       const Int_t nMax,  const Float_t t12diff,
                                   const Float_t mass,   const Float_t l0,
                                   const Float_t eta,    const Float_t phi,
                                   const Bool_t matched, const Int_t mcindex);
  
  void         FillIdConvHistograms(const Float_t en,    const Int_t nMax, const Float_t asym,
                                    const Float_t mass,   const Float_t l0,
                                    const Bool_t matched, const Int_t mcindex);
  
  void         Init();
  
  void         InitParameters();
  
  void         MakeAnalysisFillHistograms() ;
  
  void         Print(const Option_t * opt) const;

  void         SetCalorimeter(TString & det)             { fCalorimeter = det ; }
  
  void         SetMinNCells(Int_t cut)                   { fMinNCells   = cut ; }

  void         SetMinBadChannelDistance(Float_t cut)     { fMinBadDist  = cut ; }

  void         SetWCorrectionParameter(Float_t p = 0.07) { fWSimu       = p   ; }
  
  void         SwitchOnFillAngleHistograms()             { fFillAngleHisto      = kTRUE  ; }
  void         SwitchOffFillAngleHistograms()            { fFillAngleHisto      = kFALSE ; }
  
  void         SwitchOnFillExtraSSHistograms()           { fFillSSExtraHisto    = kTRUE  ; }
  void         SwitchOffFillExtraSSHistograms()          { fFillSSExtraHisto    = kFALSE ; }
  
  void         SwitchOnFillHighMultHistograms()          { fFillHighMultHisto   = kTRUE  ; }
  void         SwitchOffFillHighMultHistograms()         { fFillHighMultHisto   = kFALSE ; }
  
  void         SwitchOnFillIdConvHistograms()            { fFillIdConvHisto     = kTRUE  ; }
  void         SwitchOffFillIdConvHistograms()           { fFillIdConvHisto     = kFALSE ; }

  void         SwitchOnFillIdEtaHistograms()             { fFillIdEtaHisto      = kTRUE  ; }
  void         SwitchOffFillIdEtaHistograms()            { fFillIdEtaHisto      = kFALSE ; }
  
  void         SwitchOnFillTMHistograms()                { fFillTMHisto         = kTRUE  ; }
  void         SwitchOffFillTMHistograms()               { fFillTMHisto         = kFALSE ; }
  
  void         SwitchOnFillTMResidualHistograms()        { fFillTMResidualHisto = kTRUE  ; }
  void         SwitchOffFillTMResidualHistograms()       { fFillTMResidualHisto = kFALSE ; }
  
  void         SwitchOnFillMCPrimaryHistograms()         { fFillMCHisto         = kTRUE  ; }
  void         SwitchOffFillMCPrimaryHistograms()        { fFillMCHisto         = kFALSE ; }

  void         SwitchOnFillSSWeightHistograms()          { fFillSSWeightHisto   = kTRUE  ; }
  void         SwitchOffFillSSWeightHistograms()         { fFillSSWeightHisto   = kFALSE ; }

  void         SwitchOnFillEbinHistograms()              { fFillEbinHisto       = kTRUE  ; }
  void         SwitchOffFillEbinHistograms()             { fFillEbinHisto       = kFALSE ; }
  
  void         SwitchOnFillMCOverlapHistograms()         { fFillMCOverlapHisto  = kTRUE  ; }
  void         SwitchOffFillMCOverlapHistograms()        { fFillMCOverlapHisto  = kFALSE ; }

  void         SwitchOnFillNCellHistograms()             { fFillNCellHisto      = kTRUE  ; }
  void         SwitchOffFillNCellHistograms()            { fFillNCellHisto      = kFALSE ; }
  
  void         SwitchOnSplitClusterDistToBad()           { fCheckSplitDistToBad = kTRUE  ; }
  void         SwitchOffSplitClusterDistToBad()          { fCheckSplitDistToBad = kFALSE ; }
  
  void         SetNWeightForShowerShape(Int_t n)           { fSSWeightN = n ; }
  void         SetWeightForShowerShape(Int_t i, Float_t v) { if (i < 10) fSSWeight[i] = v ; }

  void         SetNECellCutForShowerShape(Int_t n)           { fSSECellCutN = n ; }
  void         SetECellCutForShowerShape(Int_t i, Float_t v) { if (i < 10) fSSECellCut[i] = v ; }

 
  void         RecalculateClusterShowerShapeParametersWithCellCut(const AliEMCALGeometry * geom, AliVCaloCells* cells, AliVCluster * cluster,
                                                   Float_t & l0,   Float_t & l1,
                                                   Float_t & disp, Float_t & dEta, Float_t & dPhi,
                                                   Float_t & sEta, Float_t & sPhi, Float_t & sEtaPhi,Float_t eCellMin = 0.);

  
  //For histograms
  enum mcTypes { kmcPhoton = 1, kmcConversion = 2, kmcPi0    = 3,  
                 kmcEta    = 4, kmcElectron   = 5, kmcHadron = 6, kmcPi0Conv = 7 };

 private:
  
  TString      fCalorimeter ;          // Calorimeter where the gamma is searched
  Int_t        fMinNCells   ;          // Study clusters with ncells larger than cut
  Float_t      fMinBadDist  ;          // Minimal distance to bad channel to accept cluster
  Float_t      fHistoECut   ;          // Fixed E cut for some histograms
  Bool_t       fCheckSplitDistToBad;   // Check the distance to bad channel and to EMCal borders of split clusters
  
  Bool_t       fFillAngleHisto;        // Fill splitted clusters angle histograms
  Bool_t       fFillTMHisto ;          // Fill track matching histos,
  Bool_t       fFillTMResidualHisto ;  // Fill track matching histos, residuals
  Bool_t       fFillSSExtraHisto ;     // Fill shower shape extra histos
  Bool_t       fFillMCHisto ;          // Fill MC energy fraction histos
  Bool_t       fFillSSWeightHisto ;    // Fill weigth histograms
  Bool_t       fFillEbinHisto ;        // Fill E bin histograms
  Bool_t       fFillMCOverlapHisto ;   // Fill MC particles overlap histograms
  Bool_t       fFillNCellHisto ;       // Fill n cells in cluster dependent histograms
  Bool_t       fFillIdConvHisto ;      // Fill histograms for clusters identified as conversion
  Bool_t       fFillIdEtaHisto ;       // Fill histograms for clusters identified as Eta
  Bool_t       fFillHighMultHisto;     // Fill centrality/event plane histograms
  
  Float_t      fSSWeight[10];          // List of weights to test
  Int_t        fSSWeightN;             // Total number of weights to test
  
  Float_t      fSSECellCut[10];        // List of cell min energy cuts to test
  Int_t        fSSECellCutN;           // Total number of cell min energy cuts to test
  
  Float_t      fWSimu;                 // Slope of the linear correction factor for the shower
                                       // shape weight in simulation, about 0.07
  
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
  
  TH2F       * fhAsyMCGenRecoDiffMCPi0[3];              //! reconstructed-generated asymmetry of splitted clusters vs E from pi0, for 3 NLM cases
  TH2F       * fhAsyMCGenRecoDiffMCPi0Conv[3];          //! reconstructed-generated asymmetry of splitted clusters vs E from converted pi0, for 3 NLM cases
  
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

  TH2F       * fhSplitClusterENLocMax[8][2] ;           //! Number of maxima in cluster vs E of splitted clusters, 1-6 for different MC particle types
  TH2F       * fhSplitClusterEPi0NLocMax[8][2] ;        //! Number of maxima in cluster vs E of splitted clusters when cluster id as pi0, 1-6 for different MC particle types
  
  TH2F       * fhM02NLocMax1  [8][2] ;                  //! M02 vs E for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhM02NLocMax2  [8][2] ;                  //! M02 vs E for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhM02NLocMaxN  [8][2] ;                  //! M02 vs E for N max in cluster > 2, 1-6 for different MC particle types
  
  TH2F       * fhMCAsymM02NLocMax1MCPi0Ebin[4] ;        //! M02 vs decay asymmetry for N max in cluster = 1, for 4 energy bins
  TH2F       * fhMCAsymM02NLocMax2MCPi0Ebin[4] ;        //! M02 vs decay asymmetry for N max in cluster = 2, for 4 energy bins
  TH2F       * fhMCAsymM02NLocMaxNMCPi0Ebin[4] ;        //! M02 vs decay asymmetry for N max in cluster > 2, for 4 energy bins
  
  TH2F       * fhMCGenFracNLocMax1[8][2] ;              //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenFracNLocMax2[8][2] ;              //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenFracNLocMaxN[8][2] ;              //! E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types  

  TH2F       * fhMCGenFracNLocMax1NoOverlap[8][2] ;     //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types, no overlap found
  TH2F       * fhMCGenFracNLocMax2NoOverlap[8][2] ;     //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types, no overlap found
  TH2F       * fhMCGenFracNLocMaxNNoOverlap[8][2] ;     //! E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types, no overlap found
  
  TH2F       * fhMCGenFracAfterCutsNLocMax1MCPi0 ;      //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 1, MCPi0 after M02 and asymmetry cut
  TH2F       * fhMCGenFracAfterCutsNLocMax2MCPi0 ;      //! E generated particle / E reconstructed vs E reconstructed for N max in cluster = 2, MCPi0, after M02 and asymmetry cut
  TH2F       * fhMCGenFracAfterCutsNLocMaxNMCPi0 ;      //! E generated particle / E reconstructed vs E reconstructed for N max in cluster > 2, MCPi0, after M02 and asymmetry cut

  TH2F       * fhMCGenSplitEFracNLocMax1[8][2] ;        //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types
  TH2F       * fhMCGenSplitEFracNLocMax2[8][2] ;        //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types
  TH2F       * fhMCGenSplitEFracNLocMaxN[8][2] ;        //! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types  

  TH2F       * fhMCGenSplitEFracNLocMax1NoOverlap[8][2];//! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 1, 1-6 for different MC particle types, no overlap
  TH2F       * fhMCGenSplitEFracNLocMax2NoOverlap[8][2];//! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster = 2, 1-6 for different MC particle types, no overlap
  TH2F       * fhMCGenSplitEFracNLocMaxNNoOverlap[8][2];//! E generated particle / E1+E2 reconstructed vs E reconstructed for N max in cluster > 2, 1-6 for different MC particle types, no overlap
  
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
  
  TH2F       * fhNCellMassEHighNLocMax1MCPi0 ;          //! n cells in cluster vs mass for high energy clusters,  for N max in cluster = 1, for MC pi0
  TH2F       * fhNCellM02EHighNLocMax1MCPi0  ;          //! n cells in cluster vs m02  for high energy clusters,  for N max in cluster = 1, for MC pi0
  TH2F       * fhNCellMassELowNLocMax1MCPi0  ;          //! n cells in cluster vs mass for low  energy clusters,  for N max in cluster = 1, for MC pi0
  TH2F       * fhNCellM02ELowNLocMax1MCPi0   ;          //! n cells in cluster vs m02  for low  energy clusters,  for N max in cluster = 1, for MC pi0

  TH2F       * fhNCellMassEHighNLocMax2MCPi0 ;          //! n cells in cluster vs mass for high energy clusters,  for N max in cluster = 2, for MC pi0
  TH2F       * fhNCellM02EHighNLocMax2MCPi0  ;          //! n cells in cluster vs m02  for high energy clusters,  for N max in cluster = 2, for MC pi0
  TH2F       * fhNCellMassELowNLocMax2MCPi0  ;          //! n cells in cluster vs mass for low  energy clusters,  for N max in cluster = 2, for MC pi0
  TH2F       * fhNCellM02ELowNLocMax2MCPi0   ;          //! n cells in cluster vs m02  for low  energy clusters,  for N max in cluster = 2, for MC pi0
  
  TH2F       * fhNCellMassEHighNLocMaxNMCPi0 ;          //! n cells in cluster vs mass for high energy clusters,  for N max in cluster > 2, for MC pi0
  TH2F       * fhNCellM02EHighNLocMaxNMCPi0  ;          //! n cells in cluster vs m02  for high energy clusters,  for N max in cluster > 2, for MC pi0
  TH2F       * fhNCellMassELowNLocMaxNMCPi0  ;          //! n cells in cluster vs mass for low  energy clusters,  for N max in cluster > 2, for MC pi0
  TH2F       * fhNCellM02ELowNLocMaxNMCPi0   ;          //! n cells in cluster vs m02  for low  energy clusters,  for N max in cluster > 2, for MC pi0
  
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

  TH2F       * fhNCellPi0NLocMax1[8][2] ;               //! n cells for selected around pi0, N Local Maxima = 1
  TH2F       * fhNCellEtaNLocMax1[8][2] ;               //! n cells for selected around eta, N Local Maxima = 1
  TH2F       * fhNCellPi0NLocMax2[8][2] ;               //! n cells for selected around pi0, N Local Maxima = 2
  TH2F       * fhNCellEtaNLocMax2[8][2] ;               //! n cells for selected around eta, N Local Maxima = 2
  TH2F       * fhNCellPi0NLocMaxN[8][2] ;               //! n cells for selected around pi0, N Local Maxima > 2
  TH2F       * fhNCellEtaNLocMaxN[8][2] ;               //! n cells for selected around eta, N Local Maxima > 2
  
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

  TH2F       * fhTrackMatchedDEtaNLocMax1Pos[8] ;       //! Eta distance between track and cluster vs cluster E, 1 local maximum
  TH2F       * fhTrackMatchedDPhiNLocMax1Pos[8] ;       //! Phi distance between track and cluster vs cluster E, 1 local maximum
  TH2F       * fhTrackMatchedDEtaNLocMax2Pos[8] ;       //! Eta distance between track and cluster vs cluster E, 2 local maximum
  TH2F       * fhTrackMatchedDPhiNLocMax2Pos[8] ;       //! Phi distance between track and cluster vs cluster E, 2 local maximum
  TH2F       * fhTrackMatchedDEtaNLocMaxNPos[8] ;       //! Eta distance between track and cluster vs cluster E, more than 2 local maximum
  TH2F       * fhTrackMatchedDPhiNLocMaxNPos[8] ;       //! Phi distance between track and cluster vs cluster E, more than 2 local maximum

  TH2F       * fhTrackMatchedDEtaNLocMax1Neg[8] ;       //! Eta distance between track and cluster vs cluster E, 1 local maximum
  TH2F       * fhTrackMatchedDPhiNLocMax1Neg[8] ;       //! Phi distance between track and cluster vs cluster E, 1 local maximum
  TH2F       * fhTrackMatchedDEtaNLocMax2Neg[8] ;       //! Eta distance between track and cluster vs cluster E, 2 local maximum
  TH2F       * fhTrackMatchedDPhiNLocMax2Neg[8] ;       //! Phi distance between track and cluster vs cluster E, 2 local maximum
  TH2F       * fhTrackMatchedDEtaNLocMaxNNeg[8] ;       //! Eta distance between track and cluster vs cluster E, more than 2 local maximum
  TH2F       * fhTrackMatchedDPhiNLocMaxNNeg[8] ;       //! Phi distance between track and cluster vs cluster E, more than 2 local maximum

  TH2F       * fhCentralityPi0NLocMax1 ;                //! Centrality for selected pi0, N Local Maxima = 1
  TH2F       * fhCentralityEtaNLocMax1 ;                //! Centrality for selected eta, N Local Maxima = 1
  TH2F       * fhCentralityPi0NLocMax2 ;                //! Centrality for selected pi0, N Local Maxima = 2
  TH2F       * fhCentralityEtaNLocMax2 ;                //! Centrality for selected eta, N Local Maxima = 2
  TH2F       * fhCentralityPi0NLocMaxN ;                //! Centrality for selected pi0, N Local Maxima > 2
  TH2F       * fhCentralityEtaNLocMaxN ;                //! Centrality for selected eta, N Local Maxima > 2

  TH2F       * fhEventPlanePi0NLocMax1 ;                //! Event plane for selected pi0, N Local Maxima = 1
  TH2F       * fhEventPlaneEtaNLocMax1 ;                //! Event plane for selected eta, N Local Maxima = 1
  TH2F       * fhEventPlanePi0NLocMax2 ;                //! Event plane for selected pi0, N Local Maxima = 2
  TH2F       * fhEventPlaneEtaNLocMax2 ;                //! Event plane for selected eta, N Local Maxima = 2
  TH2F       * fhEventPlanePi0NLocMaxN ;                //! Event plane for selected pi0, N Local Maxima > 2
  TH2F       * fhEventPlaneEtaNLocMaxN ;                //! Event plane for selected eta, N Local Maxima > 2

  TH2F       * fhClusterEtaPhiNLocMax1 ;                //! Eta vs Phi of clusters with N Local Maxima = 1, E > 8 GeV
  TH2F       * fhClusterEtaPhiNLocMax2 ;                //! Eta vs Phi of clusters with N Local Maxima = 2, E > 8 GeV
  TH2F       * fhClusterEtaPhiNLocMaxN ;                //! Eta vs Phi of clusters with N Local Maxima > 2, E > 8 GeV
  TH2F       * fhPi0EtaPhiNLocMax1 ;                    //! Eta vs Phi of pi0's with N Local Maxima = 1, E > 8 GeV
  TH2F       * fhPi0EtaPhiNLocMax2 ;                    //! Eta vs Phi of pi0's with N Local Maxima = 2, E > 8 GeV
  TH2F       * fhPi0EtaPhiNLocMaxN ;                    //! Eta vs Phi of pi0's with N Local Maxima > N, E > 8 GeV
  TH2F       * fhEtaEtaPhiNLocMax1 ;                    //! Eta vs Phi of eta's with N Local Maxima = 1, E > 8 GeV
  TH2F       * fhEtaEtaPhiNLocMax2 ;                    //! Eta vs Phi of eta's with N Local Maxima = 2, E > 8 GeV
  TH2F       * fhEtaEtaPhiNLocMaxN ;                    //! Eta vs Phi of eta's with N Local Maxima > N, E > 8 GeV

  TH2F       * fhPi0CellE[3] ;                          //! pi0's energy vs cluster cell energy with NLM = 1, = 2, > 2 
  TH2F       * fhPi0CellEFrac[3] ;                      //! pi0's energy vs cluster cell energy fraction with NLM = 1, = 2, > 2 
  TH2F       * fhPi0CellLogEFrac[3] ;                   //! pi0's energy vs cluster log cell energy fraction with NLM = 1, = 2, > 2
  TH2F       * fhPi0CellEMaxEMax2Frac   [3];            //! pi0's energy vs fraction of 2 main maxima energy with NLM = 1, = 2, > 2
  TH2F       * fhPi0CellEMaxClusterFrac [3];            //! pi0's energy vs energy fraction of main   LM and cluster energy with NLM = 1, = 2, > 2
  TH2F       * fhPi0CellEMax2ClusterFrac[3];            //! pi0's energy vs energy fraction of second LM and cluster energy with NLM = 1, = 2, > 2
  TH2F       * fhPi0CellEMaxFrac [3];                   //! pi0's energy vs energy fraction of main LM and cluster cell energy with NLM = 1, = 2, > 2
  TH2F       * fhPi0CellEMax2Frac [3];                  //! pi0's energy vs energy fraction of second LM and cluster cell energy with NLM = 1, = 2, > 2
  
  TH2F       * fhM02WeightPi0[3][10] ;                  //! M02 for selected pi0 with different weight, with NLM = 1, = 2, > 2
  TH2F       * fhM02ECellCutPi0[3][10] ;                //! M02 for selected pi0 with different cut on cell energy, with NLM = 1, = 2, > 2

  TH2F       * fhPi0EPairDiffTimeNLM1;                  //! E vs Pair of clusters time difference vs E, for selected pi0, NLM=1
  TH2F       * fhPi0EPairDiffTimeNLM2;                  //! E vs Pair of clusters time difference vs E, for selected pi0, NLM=2
  TH2F       * fhPi0EPairDiffTimeNLMN;                  //! E vs Pair of clusters time difference vs E, for selected pi0, NLM>2
  TH2F       * fhEtaEPairDiffTimeNLM1;                  //! E vs Pair of clusters time difference vs E, for selected eta, NLM=1
  TH2F       * fhEtaEPairDiffTimeNLM2;                  //! E vs Pair of clusters time difference vs E, for selected eta, NLM=2
  TH2F       * fhEtaEPairDiffTimeNLMN;                  //! E vs Pair of clusters time difference vs E, for selected eta, NLM>2

  TH2F       * fhMCEM02Overlap0[3][8];                  //! E vs M02 for different MC origin, no other MC particles contributes, neutral cluster
  TH2F       * fhMCEM02Overlap1[3][8];                  //! E vs M02 for different MC origin, 1  other MC particles contributes, neutral cluster
  TH2F       * fhMCEM02OverlapN[3][8];                  //! E vs M02 for different MC origin, N  other MC particles contributes, neutral cluster
  TH2F       * fhMCEM02Overlap0Match[3][8];             //! E vs M02 for different MC origin, no other MC particles contributes, charged cluster
  TH2F       * fhMCEM02Overlap1Match[3][8];             //! E vs M02 for different MC origin, 1  other MC particles contributes, charged cluster
  TH2F       * fhMCEM02OverlapNMatch[3][8];             //! E vs M02 for different MC origin, N  other MC particles contributes, charged cluster
  
  TH2F       * fhMCEMassOverlap0[3][8];                 //! E vs Mass for different MC origin, no other MC particles contributes, neutral cluster
  TH2F       * fhMCEMassOverlap1[3][8];                 //! E vs Mass for different MC origin, 1  other MC particles contributes, neutral cluster
  TH2F       * fhMCEMassOverlapN[3][8];                 //! E vs Mass for different MC origin, N  other MC particles contributes, neutral cluster
  TH2F       * fhMCEMassOverlap0Match[3][8];            //! E vs Mass for different MC origin, no other MC particles contributes, charged cluster
  TH2F       * fhMCEMassOverlap1Match[3][8];            //! E vs Mass for different MC origin, 1  other MC particles contributes, charged cluster
  TH2F       * fhMCEMassOverlapNMatch[3][8];            //! E vs Mass for different MC origin, N  other MC particles contributes, charged cluster

  TH2F       * fhMCESplitEFracOverlap0[3][8];           //! E vs sum of splitted cluster energy / cluster energy for different MC origin, no other MC particles contributes, neutral cluster
  TH2F       * fhMCESplitEFracOverlap1[3][8];           //! E vs sum of splitted cluster energy / cluster energy for different MC origin, 1  other MC particles contributes, neutral cluster
  TH2F       * fhMCESplitEFracOverlapN[3][8];           //! E vs sum of splitted cluster energy / cluster energy for different MC origin, N  other MC particles contributes, neutral cluster
  TH2F       * fhMCESplitEFracOverlap0Match[3][8];      //! E vs sum of splitted cluster energy / cluster energy for different MC origin, no other MC particles contributes, charged cluster
  TH2F       * fhMCESplitEFracOverlap1Match[3][8];      //! E vs sum of splitted cluster energy / cluster energy for different MC origin, 1  other MC particles contributes, charged cluster
  TH2F       * fhMCESplitEFracOverlapNMatch[3][8];      //! E vs sum of splitted cluster energy / cluster energy for different MC origin, N  other MC particles contributes, charged cluster

  TH2F       * fhMCEAsymOverlap0[3][8];                 //! E vs sum of splitted cluster energy asymmetry for different MC origin, no other MC particles contributes, neutral cluster
  TH2F       * fhMCEAsymOverlap1[3][8];                 //! E vs sum of splitted cluster energy asymmetry for different MC origin, 1  other MC particles contributes, neutral cluster
  TH2F       * fhMCEAsymOverlapN[3][8];                 //! E vs sum of splitted cluster energy asymmetry for different MC origin, N  other MC particles contributes, neutral cluster
  TH2F       * fhMCEAsymOverlap0Match[3][8];            //! E vs sum of splitted cluster energy asymmetry for different MC origin, no other MC particles contributes, charged cluster
  TH2F       * fhMCEAsymOverlap1Match[3][8];            //! E vs sum of splitted cluster energy asymmetry for different MC origin, 1  other MC particles contributes, charged cluster
  TH2F       * fhMCEAsymOverlapNMatch[3][8];            //! E vs sum of splitted cluster energy asymmetry for different MC origin, N  other MC particles contributes, charged cluster

  TH2F       * fhMCENCellOverlap0[3][8];                //! E vs sum of splitted cluster energy asymmetry for different MC origin, no other MC particles contributes, neutral cluster
  TH2F       * fhMCENCellOverlap1[3][8];                //! E vs sum of splitted cluster energy asymmetry for different MC origin, 1  other MC particles contributes, neutral cluster
  TH2F       * fhMCENCellOverlapN[3][8];                //! E vs sum of splitted cluster energy asymmetry for different MC origin, N  other MC particles contributes, neutral cluster
  TH2F       * fhMCENCellOverlap0Match[3][8];           //! E vs sum of splitted cluster energy asymmetry for different MC origin, no other MC particles contributes, charged cluster
  TH2F       * fhMCENCellOverlap1Match[3][8];           //! E vs sum of splitted cluster energy asymmetry for different MC origin, 1  other MC particles contributes, charged cluster
  TH2F       * fhMCENCellOverlapNMatch[3][8];           //! E vs sum of splitted cluster energy asymmetry for different MC origin, N  other MC particles contributes, charged cluster
  
  TH2F       * fhMCEEpriOverlap0[3][8];                 //! E reco vs primary for different MC origin, no other MC particles contributes, neutral cluster
  TH2F       * fhMCEEpriOverlap1[3][8];                 //! E reco vs primary for different MC origin, 1  other MC particles contributes, neutral cluster
  TH2F       * fhMCEEpriOverlapN[3][8];                 //! E reco vs primary for different MC origin, N  other MC particles contributes, neutral cluster
  TH2F       * fhMCEEpriOverlap0Match[3][8];            //! E reco vs primary for different MC origin, no other MC particles contributes, charged cluster
  TH2F       * fhMCEEpriOverlap1Match[3][8];            //! E reco vs primary for different MC origin, 1  other MC particles contributes, charged cluster
  TH2F       * fhMCEEpriOverlapNMatch[3][8];            //! E reco vs primary for different MC origin, N  other MC particles contributes, charged cluster
  
  TH2F       * fhMCEEpriOverlap0IdPi0[3][8];            //! E reco vs primary for different MC origin, no other MC particles contributes, neutral cluster, neutral clusters id as pi0
  TH2F       * fhMCEEpriOverlap1IdPi0[3][8];            //! E reco vs primary for different MC origin, 1  other MC particles contributes, neutral cluster, neutral clusters id as pi0
  TH2F       * fhMCEEpriOverlapNIdPi0[3][8];            //! E reco vs primary for different MC origin, 1  other MC particles contributes, neutral cluster, neutral clusters is as pi0
  
  TH2F       * fhMCPi0MassM02Overlap0[3][4];            //! MC Pi0 M02 vs Mass for different MC origin, no other MC particles contributes, neutral cluster, 4 E bins
  TH2F       * fhMCPi0MassM02Overlap1[3][4];            //! MC Pi0 M02 vs Mass for different MC origin, 1  other MC particles contributes, neutral cluster, 4 E bins
  TH2F       * fhMCPi0MassM02OverlapN[3][4];            //! MC Pi0 M02 vs Mass for different MC origin, N  other MC particles contributes, neutral cluster, 4 E bins
  TH2F       * fhMCPi0MassM02Overlap0Match[3][4];       //! MC Pi0 M02 vs Mass for different MC origin, no other MC particles contributes, charged cluster, 4 E bins
  TH2F       * fhMCPi0MassM02Overlap1Match[3][4];       //! MC Pi0 M02 vs Mass for different MC origin, 1  other MC particles contributes, charged cluster, 4 E bins
  TH2F       * fhMCPi0MassM02OverlapNMatch[3][4];       //! MC Pi0 M02 vs Mass for different MC origin, N  other MC particles contributes, charged cluster, 4 E bins
  
  TH2F       * fhMCENOverlaps[3][8];                    //! E vs number of Overlaps in MC, neutral cluster
  TH2F       * fhMCENOverlapsMatch[3][8];               //! E vs number of Overlaps in MC, charged cluster
  
  TH2F       * fhMCPi0HighNLMPair;                      //! E vs NLM when cluster originated in pi0 merging and highest energy local maxima correspond to 2 photons
  TH2F       * fhMCPi0LowNLMPair;                       //! E vs NLM when cluster originated in pi0 merging and a pair of local maxima except highest energy correspond to 2 photons
  TH2F       * fhMCPi0AnyNLMPair;                       //! E vs NLM when cluster originated in pi0 merging and a both highest energy pairs and other pairs correspond to 2 photons
  TH2F       * fhMCPi0NoneNLMPair;                      //! E vs NLM when cluster originated in pi0 merging and a both no NLM corresponds to the photons
  // No match between highest energy local maxima and highest energy MC particle 
  TH2F       * fhMCPi0HighNLMPairNoMCMatch;             //! E vs NLM when cluster originated in pi0 merging and highest energy local maxima correspond to 2 photons
  TH2F       * fhMCPi0LowNLMPairNoMCMatch;              //! E vs NLM when cluster originated in pi0 merging and a pair of local maxima except highest energy correspond to 2 photons
  TH2F       * fhMCPi0AnyNLMPairNoMCMatch;              //! E vs NLM when cluster originated in pi0 merging and a both highest energy pairs and other pairs correspond to 2 photons
  TH2F       * fhMCPi0NoneNLMPairNoMCMatch;             //! E vs NLM when cluster originated in pi0 merging and a both no NLM corresponds to the photons

  TH2F       * fhMCPi0HighNLMPairOverlap;              //! E vs NLM when cluster originated in pi0 merging and highest energy local maxima correspond to 2 photons, overlap
  TH2F       * fhMCPi0LowNLMPairOverlap;               //! E vs NLM when cluster originated in pi0 merging and a pair of local maxima except highest energy correspond to 2 photons, overlap
  TH2F       * fhMCPi0AnyNLMPairOverlap;               //! E vs NLM when cluster originated in pi0 merging and a both highest energy pairs and other pairs correspond to 2 photons, overlap
  TH2F       * fhMCPi0NoneNLMPairOverlap;              //! E vs NLM when cluster originated in pi0 merging and a both no NLM corresponds to the photons, overlap
  // No match between highest energy local maxima and highest energy MC particle
  TH2F       * fhMCPi0HighNLMPairNoMCMatchOverlap;     //! E vs NLM when cluster originated in pi0 merging and highest energy local maxima correspond to 2 photons, overlap
  TH2F       * fhMCPi0LowNLMPairNoMCMatchOverlap;      //! E vs NLM when cluster originated in pi0 merging and a pair of local maxima except highest energy correspond to 2 photons, overlap
  TH2F       * fhMCPi0AnyNLMPairNoMCMatchOverlap;      //! E vs NLM when cluster originated in pi0 merging and a both highest energy pairs and other pairs correspond to 2 photons, overlap
  TH2F       * fhMCPi0NoneNLMPairNoMCMatchOverlap;     //! E vs NLM when cluster originated in pi0 merging and a both no NLM corresponds to the photons, overlap
  
  TH2F       * fhMCPi0DecayPhotonHitHighLM;             //! E vs NLM when cluster originated in pi0 merging and MC photon decay hit the cell local maxima
  TH2F       * fhMCPi0DecayPhotonAdjHighLM;             //! E vs NLM when cluster originated in pi0 merging and MC photon decay hit the adjacent cell local maxima
  TH2F       * fhMCPi0DecayPhotonHitOtherLM;            //! E vs NLM when cluster originated in pi0 merging and MC photon decay hit the cell local maximas, not high
  TH2F       * fhMCPi0DecayPhotonAdjOtherLM;            //! E vs NLM when cluster originated in pi0 merging and MC photon decay do not hit the adjacent cell local maximas, not high
  TH2F       * fhMCPi0DecayPhotonAdjacent;              //! E vs NLM when cluster originated in pi0 merging and MC photon decay hit adjacen cells, not 2 LM
  TH2F       * fhMCPi0DecayPhotonHitNoLM;               //! E vs NLM when cluster originated in pi0 merging and MC photon decay do not hit the cell local maximas
  
  TH2F       * fhMCPi0DecayPhotonHitHighLMOverlap;      //! E vs NLM when cluster originated in pi0 merging and MC photon decay hit the cell local maxima, overlap
  TH2F       * fhMCPi0DecayPhotonAdjHighLMOverlap;      //! E vs NLM when cluster originated in pi0 merging and MC photon decay hit the adjacent cell local maxima, overlap
  TH2F       * fhMCPi0DecayPhotonHitOtherLMOverlap;     //! E vs NLM when cluster originated in pi0 merging and MC photon decay hit the cell local maximas, not high, overlap
  TH2F       * fhMCPi0DecayPhotonAdjOtherLMOverlap;     //! E vs NLM when cluster originated in pi0 merging and MC photon decay do not hit the adjacent cell local maximas, not high, overlap
  TH2F       * fhMCPi0DecayPhotonAdjacentOverlap;       //! E vs NLM when cluster originated in pi0 merging and MC photon decay hit adjacen cells, not 2 LM, overlap
  TH2F       * fhMCPi0DecayPhotonHitNoLMOverlap;        //! E vs NLM when cluster originated in pi0 merging and MC photon decay do not hit the cell local maximas, overlap

  TH2F       * fhMCPi0DecayPhotonHitHighLMDiffELM1[3];             //! E vs Ephoton-Esplit cluster when cluster originated in pi0 merging and MC photon decay hit the cell local maxima
  TH2F       * fhMCPi0DecayPhotonAdjHighLMDiffELM1[3];             //! E vs Ephoton-Esplit cluster when cluster originated in pi0 merging and MC photon decay hit the adjacent cell local maxima
  TH2F       * fhMCPi0DecayPhotonHitOtherLMDiffELM1[3];            //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay hit the cell local maximas, not high
  TH2F       * fhMCPi0DecayPhotonAdjOtherLMDiffELM1[3];            //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay do not hit the adjacent cell local maximas, not high
 
  TH2F       * fhMCPi0DecayPhotonHitHighLMOverlapDiffELM1[3];      //! E vs Ephoton-Esplit cluster when cluster originated in pi0 merging and MC photon decay hit the cell local maxima
  TH2F       * fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM1[3];      //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay hit the adjacent cell local maxima, overlap
  TH2F       * fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM1[3];     //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay hit the cell local maximas, not high, overlap
  TH2F       * fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM1[3];     //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay do not hit the adjacent cell local maximas, not high, overlap
 
  TH2F       * fhMCPi0DecayPhotonHitHighLMDiffELM2[3];             //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay hit the cell local maxima
  TH2F       * fhMCPi0DecayPhotonAdjHighLMDiffELM2[3];             //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay hit the adjacent cell local maxima
  TH2F       * fhMCPi0DecayPhotonHitOtherLMDiffELM2[3];            //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay hit the cell local maximas, not high
  TH2F       * fhMCPi0DecayPhotonAdjOtherLMDiffELM2[3];            //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay do not hit the adjacent cell local maximas, not high
  
  TH2F       * fhMCPi0DecayPhotonHitHighLMOverlapDiffELM2[3];      //! E vs Ephoton-Esplit cluster when cluster originated in pi0 merging and MC photon decay hit the cell local maxima
  TH2F       * fhMCPi0DecayPhotonAdjHighLMOverlapDiffELM2[3];      //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay hit the adjacent cell local maxima, overlap
  TH2F       * fhMCPi0DecayPhotonHitOtherLMOverlapDiffELM2[3];     //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay hit the cell local maximas, not high, overlap
  TH2F       * fhMCPi0DecayPhotonAdjOtherLMOverlapDiffELM2[3];     //! E vs Ephoton-Esplit when cluster originated in pi0 merging and MC photon decay do not hit the adjacent cell local maximas, not high, overlap
  
  TH2F       * fhMCPi0DecayPhotonHitHighLMMass[3];                 //! E vs Mass when cluster originated in pi0 merging and MC photon decay hit the cell local maxima
  TH2F       * fhMCPi0DecayPhotonAdjHighLMMass[3];                 //! E vs Mass when cluster originated in pi0 merging and MC photon decay hit the adjacent cell local maxima
  TH2F       * fhMCPi0DecayPhotonHitOtherLMMass[3];                //! E vs Mass when cluster originated in pi0 merging and MC photon decay hit the cell local maximas, not high
  TH2F       * fhMCPi0DecayPhotonAdjOtherLMMass[3];                //! E vs Mass when cluster originated in pi0 merging and MC photon decay do not hit the adjacent cell local maximas, not high
  TH2F       * fhMCPi0DecayPhotonAdjacentMass[3];                  //! E vs Mass when cluster originated in pi0 merging and MC photon decay hit adjacen cells, not 2 LM
  TH2F       * fhMCPi0DecayPhotonHitNoLMMass[3];                   //! E vs Mass when cluster originated in pi0 merging and MC photon decay do not hit the cell local maximas
               
  TH2F       * fhMCPi0DecayPhotonHitHighLMOverlapMass[3];          //! E vs Mass when cluster originated in pi0 merging and MC photon decay hit the cell local maxima, overlap
  TH2F       * fhMCPi0DecayPhotonAdjHighLMOverlapMass[3];          //! E vs Mass when cluster originated in pi0 merging and MC photon decay hit the adjacent cell local maxima, overlap
  TH2F       * fhMCPi0DecayPhotonHitOtherLMOverlapMass[3];         //! E vs Mass when cluster originated in pi0 merging and MC photon decay hit the cell local maximas, not high, overlap
  TH2F       * fhMCPi0DecayPhotonAdjOtherLMOverlapMass[3];         //! E vs Mass when cluster originated in pi0 merging and MC photon decay do not hit the adjacent cell local maximas, not high, overlap
  TH2F       * fhMCPi0DecayPhotonAdjacentOverlapMass[3];           //! E vs Mass when cluster originated in pi0 merging and MC photon decay hit adjacen cells, not 2 LM, overlap
  TH2F       * fhMCPi0DecayPhotonHitNoLMOverlapMass[3];            //! E vs Mass when cluster originated in pi0 merging and MC photon decay do not hit the cell local maximas, overlap

  
  TH2F       * fhMCEOverlapType;                        //! what particles overlap with pi0, neutral clusters
  TH2F       * fhMCEOverlapTypeMatch;                   //! what particles overlap with pi0, charged clusters
  
  TH2F       * fhMassBadDistClose[3];                   //! split mass of clusters with second LM close to bad channel
  TH2F       * fhM02BadDistClose[3];                    //! m02 of clusters with second LM close to bad channel
  TH2F       * fhMassOnBorder[3];                       //! split mass of clusters with second LM on EMCAL border
  TH2F       * fhM02OnBorder[3];                        //! m02 of clusters with second LM close to EMCAL border

  
  AliAnaInsideClusterInvariantMass(              const AliAnaInsideClusterInvariantMass & split) ; // cpy ctor
  AliAnaInsideClusterInvariantMass & operator = (const AliAnaInsideClusterInvariantMass & split) ; // cpy assignment
  
  ClassDef(AliAnaInsideClusterInvariantMass,24)
  
} ;

#endif //ALIANAINSIDECLUSTERINVARIANTMASS_H



