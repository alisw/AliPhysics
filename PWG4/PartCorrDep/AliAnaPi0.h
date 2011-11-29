#ifndef ALIANAPI0_H
#define ALIANAPI0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Class to fill two-photon invariant mass histograms 
// to be used to extract pi0 raw yield.
// Input is produced by AliAnaPhoton (or any other analysis producing output AliAODPWG4Particles), 
// it will do nothing if executed alone
//
//-- Author: Dmitri Peressounko (RRC "KI")
//-- Adapted to PartCorr frame by Lamia Benhabib (SUBATECH)
//-- and Gustavo Conesa (INFN-Frascati)

//Root
class TList;
class TH3F ;
class TH2F ;
class TObjString;

//Analysis
#include "AliAnaPartCorrBaseClass.h"
class AliAODEvent ;
class AliESDEvent ;
class AliAODPWG4Particle ;

class AliAnaPi0 : public AliAnaPartCorrBaseClass {
  
 public:   
  AliAnaPi0() ; // default ctor
  virtual ~AliAnaPi0() ;//virtual dtor
  
  //-------------------------------
  // General analysis frame methods
  //-------------------------------

  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects(); 
  
  void         Print(const Option_t * opt) const;
  
  void         MakeAnalysisFillHistograms();
  
  void         InitParameters();

  //Calorimeter options
  TString      GetCalorimeter()         const   { return fCalorimeter           ; }
  void         SetCalorimeter(TString & det)    { fCalorimeter         = det    ; }
  void         SetNumberOfModules(Int_t nmod)   { fNModules            = nmod   ; }
  
  //-------------------------------
  // EVENT Bin Methods
  //-------------------------------

  Int_t        GetEventIndex(AliAODPWG4Particle * part, Double_t * vert)  ;

  void         CountAndGetAverages(Int_t &nClus,Int_t &nCell, Float_t &eClusTot,Float_t &eCellTot, Float_t &eDenClus,Float_t &eDenCell) ;
  
  //Switchs for event multiplicity bin option, by default, centrality
  
  void         SwitchOnTrackMultBins()          { fUseTrackMultBins    = kTRUE  ; }
  void         SwitchOffTrackMultBins()         { fUseTrackMultBins    = kFALSE ; }
  
  void         SwitchOnPhotonMultBins()         { fUsePhotonMultBins   = kTRUE  ; }
  void         SwitchOffPhotonMultBins()        { fUsePhotonMultBins   = kFALSE ; }
  
  void         SwitchOnClusterEBins()           { fUseAverClusterEBins = kTRUE  ; }
  void         SwitchOffClusterEBins()          { fUseAverClusterEBins = kFALSE ; }
  
  void         SwitchOnCellEBins()              { fUseAverCellEBins    = kTRUE  ; }
  void         SwitchOffCellEBins()             { fUseAverCellEBins    = kFALSE ; }

  void         SwitchOnClusterEDenBins()        { fUseAverClusterEDenBins = kTRUE  ; }
  void         SwitchOffClusterEDenBins()       { fUseAverClusterEDenBins = kFALSE ; }

  //-------------------------------
	//Opening angle pair selection
  //-------------------------------
  void         SwitchOnAngleSelection()         { fUseAngleCut         = kTRUE  ; }
  void         SwitchOffAngleSelection()        { fUseAngleCut         = kFALSE ; }
  
  void         SwitchOnAngleEDepSelection()     { fUseAngleEDepCut     = kTRUE  ; }
  void         SwitchOffAngleEDepSelection()    { fUseAngleEDepCut     = kFALSE ; }
    
  void         SetAngleCut(Float_t a)           { fAngleCut            = a      ; }
  void         SetAngleMaxCut(Float_t a)        { fAngleMaxCut         = a      ; }

  //-------------------------------
  // Use mixing code of this class
  //-------------------------------
  void         SwitchOnOwnMix()                 { fDoOwnMix            = kTRUE  ; }
  void         SwitchOffOwnMix()                { fDoOwnMix            = kFALSE ; }

  //------------------------------------------
  //Do analysis only with clusters in same SM or different combinations of SM
  //------------------------------------------
  void         SwitchOnSameSM()                 { fSameSM              = kTRUE  ; }
  void         SwitchOffSameSM()                { fSameSM              = kFALSE ; }
  
  void         SwitchOnSMCombinations()         { fFillSMCombinations  = kTRUE  ; }
  void         SwitchOffSMCombinations()        { fFillSMCombinations  = kFALSE ; }
  
  //-------------------------------
  //Histogram filling options off by default
  //-------------------------------
  void         SwitchOnInvPtWeight()            { fMakeInvPtPlots      = kTRUE  ; }
  void         SwitchOffInvPtWeight()           { fMakeInvPtPlots      = kFALSE ; }
  
  void         SwitchOnFillBadDistHisto()       { fFillBadDistHisto    = kTRUE  ; }
  void         SwitchOffFillBadDistHisto()      { fFillBadDistHisto    = kFALSE ; }
  
  //-------------------------------------------
  //Cuts for multiple analysis, off by default
  //-------------------------------------------
  void         SwitchOnMultipleCutAnalysis()    { fMultiCutAna         = kTRUE  ; }
  void         SwitchOffMultipleCutAnalysis()   { fMultiCutAna         = kFALSE ; }

  void         SetNPtCuts   (Int_t s)           { if(s <= 10)fNPtCuts    = s    ; }
  void         SetNAsymCuts (Int_t s)           { if(s <= 10)fNAsymCuts  = s    ; }
  void         SetNNCellCuts(Int_t s)           { if(s <= 10)fNCellNCuts = s    ; }
  void         SetNPIDBits  (Int_t s)           { if(s <= 10)fNPIDBits   = s    ; }
  
  void         SetPtCutsAt  (Int_t p,Float_t v) { if(p < 10)fPtCuts[p]   = v    ; }
  void         SetAsymCutsAt(Int_t p,Float_t v) { if(p < 10)fAsymCuts[p] = v    ; }
  void         SetNCellCutsAt(Int_t p,Int_t v)  { if(p < 10)fCellNCuts[p]= v    ; }
  void         SetPIDBitsAt  (Int_t p,Int_t v)  { if(p < 10)fPIDBits[p]  = v    ; }
  
  //MC analysis related methods
    
  void         SwitchOnConversionChecker()      { fCheckConversion     = kTRUE  ; }
  void         SwitchOffConversionChecker()     { fCheckConversion     = kFALSE ; }  
  
  void         SwitchOnMultipleCutAnalysisInSimulation()  { fMultiCutAnaSim = kTRUE  ; }
  void         SwitchOffMultipleCutAnalysisInSimulation() { fMultiCutAnaSim = kFALSE ; }
  
  void         FillAcceptanceHistograms();
  void         FillMCVersusRecDataHistograms(const Int_t    index1,  const Int_t    index2,
                                             const Float_t  pt1,     const Float_t  pt2, 
                                             const Int_t    ncells1, const Int_t    ncells2, 
                                             const Double_t mass,    const Double_t pt,  const Double_t asym,    
                                             const Double_t deta,    const Double_t dphi);
  
  private:

  Bool_t   fDoOwnMix;                  // Do combinatorial background not the one provided by the frame
  TList ** fEventsList ;               //![GetNCentrBin()*GetNZvertBin()*GetNRPBin()] Containers for photons in stored events

  TString  fCalorimeter ;              // Select Calorimeter for IM
  Int_t    fNModules ;                 // Number of EMCAL/PHOS modules, set as many histogras as modules 
  
  Bool_t   fUseAngleCut ;              // Select pairs depending on their opening angle
  Bool_t   fUseAngleEDepCut ;          // Select pairs depending on their opening angle
  Float_t  fAngleCut ;                 // Select pairs with opening angle larger than a threshold
  Float_t  fAngleMaxCut ;              // Select pairs with opening angle smaller than a threshold
  
  //Multiple cuts analysis
  Bool_t   fMultiCutAna;               // Do analysis with several or fixed cut
  Bool_t   fMultiCutAnaSim;            // Do analysis with several or fixed cut, in the simulation related part
  Int_t    fNPtCuts;                   // Number of pt cuts
  Float_t  fPtCuts[10];                // Array with different pt cuts
  Int_t    fNAsymCuts;                 // Number of assymmetry cuts
  Float_t  fAsymCuts[10];              // Array with different assymetry cuts
  Int_t    fNCellNCuts;                // Number of cuts with number of cells in cluster
  Int_t    fCellNCuts[10];             // Array with different cell number cluster cuts
  Int_t    fNPIDBits ;		             // Number of possible PID bit combinations
  Int_t    fPIDBits[10];               // Array with different PID bits
  
  //Switchs of different analysis options
  Bool_t   fMakeInvPtPlots;            // D plots with inverse pt weight
  Bool_t   fSameSM;                    // Select only pairs in same SM;
  Bool_t   fFillSMCombinations;        // Fill histograms with different cluster pairs in SM combinations
  Bool_t   fCheckConversion;           // Fill histograms with tagged photons as conversion
  Bool_t   fUseTrackMultBins;          // Use track multiplicity and not centrality bins
  Bool_t   fUsePhotonMultBins;         // Use photon multiplicity and not centrality bins
  Bool_t   fUseAverCellEBins;          // Use cell average energy and not centrality bins
  Bool_t   fUseAverClusterEBins;       // Use cluster average energy and not centrality bins
  Bool_t   fUseAverClusterEDenBins;    // Use cluster average energy density and not centrality bins
  Bool_t   fFillBadDistHisto;          // Do plots for different distances to bad channels
  
  //Histograms
  
  //Event characterization
  TH1F *   fhAverTotECluster;          //! Average number of clusters in SM
  TH1F *   fhAverTotECell;             //! Average number of cells    in SM
  TH2F *   fhAverTotECellvsCluster;    //! Average number of cells    in SM
  TH1F *   fhEDensityCluster;          //! Deposited energy in event per cluster
  TH1F *   fhEDensityCell;             //! Deposited energy in event per cell vs cluster
  TH2F *   fhEDensityCellvsCluster;    //! Deposited energy in event per cell vs cluster

  TH2F **  fhReMod ;                  //![fNModules]   REAL  two-photon invariant mass distribution for different calorimeter modules.
  TH2F **  fhReSameSideEMCALMod ;     //![fNModules-2] REAL  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhReSameSectorEMCALMod ;   //![fNModules/2] REAL  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhReDiffPHOSMod ;          //![fNModules]   REAL  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhMiMod ;                  //![fNModules]   MIXED two-photon invariant mass distribution for different calorimeter modules.
  TH2F **  fhMiSameSideEMCALMod ;     //![fNModules-2] REAL  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhMiSameSectorEMCALMod ;   //![fNModules/2] REAL  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhMiDiffPHOSMod ;          //![fNModules-1] REAL  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  
  // Pairs with at least one cluster tagged as conversion
  TH2F *   fhReConv ;                 //! REAL  two-photon invariant mass distribution one of the pair was 2 clusters with small mass 
  TH2F *   fhMiConv ;                 //! MIXED two-photon invariant mass distribution one of the pair was 2 clusters with small mass
  TH2F *   fhReConv2 ;                //! REAL  two-photon invariant mass distribution both pair photons recombined from 2 clusters with small mass 
  TH2F *   fhMiConv2 ;                //! MIXED two-photon invariant mass distribution both pair photons recombined from 2 clusters with small mass

  TH2F **  fhRe1 ;                    //![GetNCentrBin()*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry 
  TH2F **  fhMi1 ;                    //![GetNCentrBin()*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry
  TH2F **  fhRe2 ;                    //![GetNCentrBin()*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry 
  TH2F **  fhMi2 ;                    //![GetNCentrBin()*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry
  TH2F **  fhRe3 ;                    //![GetNCentrBin()*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry 
  TH2F **  fhMi3 ;                    //![GetNCentrBin()*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry

  //Histograms weighted by inverse pT
  TH2F **  fhReInvPt1 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  TH2F **  fhMiInvPt1 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  TH2F **  fhReInvPt2 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT 
  TH2F **  fhMiInvPt2 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  TH2F **  fhReInvPt3 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts] REAL  two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  TH2F **  fhMiInvPt3 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts] MIXED two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT
  
  //Multiple cuts: Assymmetry, pt, n cells, PID
  TH2F **  fhRePtNCellAsymCuts ;       //![fNPtCuts*fNAsymCuts*fNCellNCuts*] REAL two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry
  TH2F **  fhMiPtNCellAsymCuts ;       //![fNPtCuts*fNAsymCuts*fNCellNCuts] Mixed two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry
  TH2F **  fhRePtNCellAsymCutsSM[12] ; //![fNPtCuts*fNAsymCuts*fNCellNCutsfNModules] REAL two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry for each module
 
  TH2F **  fhRePIDBits ;               //![fNPIDBits]  REAL two-photon invariant mass distribution for different PID bits
  TH3F **  fhRePtMult ;                //![fNAsymCuts] REAL two-photon invariant mass distribution for different track multiplicity and assymetry cuts
  TH2F *   fhReSS[3] ;                 //! Combine clusters with 3 different cuts on shower shape
  
  // Asymmetry vs pt, in pi0/eta regions
  TH2F *   fhRePtAsym    ;             //! REAL two-photon pt vs asymmetry
  TH2F *   fhRePtAsymPi0 ;             //! REAL two-photon pt vs asymmetry, close to pi0 mass
  TH2F *   fhRePtAsymEta ;             //! REAL two-photon pt vs asymmetry, close to eta mass
  
  //Centrality, Event plane bins
  TH3F *   fhEvents;                   //! Number of events per centrality, RP, zbin
  TH1F *   fhCentrality;               //! Histogram with centrality bins with at least one pare
  TH1F *   fhCentralityNoPair;         //! Histogram with centrality bins with no pair

  TH1F *   fhEventPlaneAngle;          //! Histogram with Event plane angle
  TH2F *   fhEventPlaneResolution;     //! Histogram with Event plane resolution vs centrality
  
  // Pair opening angle
  TH2F *   fhRealOpeningAngle ;        //! Opening angle of pair versus pair energy
  TH2F *   fhRealCosOpeningAngle ;     //! Cosinus of opening angle of pair version pair energy
  TH2F *   fhMixedOpeningAngle ;       //! Opening angle of pair versus pair energy
  TH2F *   fhMixedCosOpeningAngle ;    //! Cosinus of opening angle of pair version pair energy
  
  //MC analysis histograms
  //Pi0 Acceptance
  TH1F *   fhPrimPi0Pt ;               //! Spectrum of Primary 
  TH1F *   fhPrimPi0AccPt ;            //! Spectrum of primary with accepted daughters 
  TH2F *   fhPrimPi0Y ;                //! Rapidity distribution of primary particles  vs pT
  TH2F *   fhPrimPi0AccY ;             //! Rapidity distribution of primary with accepted daughters  vs pT
  TH2F *   fhPrimPi0Phi ;              //! Azimutal distribution of primary particles  vs pT
  TH2F *   fhPrimPi0AccPhi;            //! Azimutal distribution of primary with accepted daughters  vs pT
  TH2F *   fhPrimPi0OpeningAngle ;     //! Opening angle of pair versus pair energy, primaries
  TH2F *   fhPrimPi0CosOpeningAngle ;  //! Cosinus of opening angle of pair version pair energy, primaries
  //Eta acceptance
  TH1F *   fhPrimEtaPt ;               //! Spectrum of Primary 
  TH1F *   fhPrimEtaAccPt ;            //! Spectrum of primary with accepted daughters 
  TH2F *   fhPrimEtaY ;                //! Rapidity distribution of primary particles vs pT
  TH2F *   fhPrimEtaAccY ;             //! Rapidity distribution of primary with accepted daughters  vs pT
  TH2F *   fhPrimEtaPhi ;              //! Azimutal distribution of primary particles  vs pT
  TH2F *   fhPrimEtaAccPhi;            //! Azimutal distribution of primary with accepted daughters	 vs pT
  
  // Primaries origin
  TH2F *   fhPrimPi0PtOrigin ;         //! Spectrum of generated pi0 vs mother
  TH2F *   fhPrimEtaPtOrigin ;         //! Spectrum of generated eta vs mother
  
  //Pair origin
  //Array of histograms ordered as follows: 0-Photon, 1-electron, 2-pi0, 3-eta, 4-a-proton, 5-a-neutron, 6-stable particles, 
  // 7-other decays, 8-string, 9-final parton, 10-initial parton, intermediate, 11-colliding proton, 12-unrelated
  TH2F *   fhMCOrgMass[13];            //! Mass vs pt of real pairs, check common origin of pair
  TH2F *   fhMCOrgAsym[13];            //! Asymmetry vs pt of real pairs, check common origin of pair
  TH2F *   fhMCOrgDeltaEta[13];        //! Delta Eta vs pt of real pairs, check common origin of pair
  TH2F *   fhMCOrgDeltaPhi[13];        //! Delta Phi vs pt of real pairs, check common origin of pair
  
  //Multiple cuts in simulation, origin pi0 or eta
  TH2F **  fhMCPi0MassPtRec;           //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real pi0 pairs, reconstructed mass vs reconstructed pt of original pair  
  TH2F **  fhMCPi0MassPtTrue;          //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real pi0 pairs, reconstructed mass vs generated pt of original pair  
  TH2F **  fhMCPi0PtTruePtRec;         //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real pi0 pairs, reconstructed pt vs generated pt of pair
  TH2F **  fhMCEtaMassPtRec;           //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real eta pairs, reconstructed mass vs reconstructed pt of original pair  
  TH2F **  fhMCEtaMassPtTrue;          //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real eta pairs, reconstructed mass vs generated pt of original pair  
  TH2F **  fhMCEtaPtTruePtRec;         //![fNPtCuts*fNAsymCuts*fNCellNCuts] Real eta pairs, reconstructed pt vs generated pt of pair

  TH2F *   fhMCPi0PtOrigin ;           //! Mass of reoconstructed pi0 pairs  in calorimeter vs mother
  TH2F *   fhMCEtaPtOrigin ;           //! Mass of reoconstructed pi0 pairs  in calorimeter vs mother

  TH2F *   fhReMCFromConversion ;      //! Invariant mass of 2 clusters originated in conversions
  TH2F *   fhReMCFromNotConversion ;   //! Invariant mass of 2 clusters not originated in conversions
  TH2F *   fhReMCFromMixConversion ;   //! Invariant mass of 2 clusters one from conversion and the other not

  AliAnaPi0(const AliAnaPi0 & g) ; // cpy ctor
  AliAnaPi0 & operator = (const AliAnaPi0 & api0) ;//cpy assignment
  
  ClassDef(AliAnaPi0,21)
} ;


#endif //ALIANAPI0_H



