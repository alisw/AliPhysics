#ifndef ALIANAPI0_H
#define ALIANAPI0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaPi0
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Selected photon clusters invariant mass analysis.
///
/// Class to collect two-photon invariant mass distributions for
/// extracting raw pi0 or eta yields.
/// Input is produced by AliAnaPhoton (or any other analysis producing output AliAODPWG4Particles),
/// it will do nothing if executed alone.
///
/// Original author: Dmitri Peressounko (RRC "KI").
/// Adapted to CaloTrackCorr frame by Lamia Benhabib (SUBATECH).
///
/// More information on the code can be found in this
/// [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this
/// [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaPi0).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// Root
class TList;
class TH3F ;
class TH2F ;
class TObjString;

// Analysis
#include "AliAnaCaloTrackCorrBaseClass.h"
class AliAODEvent ;
class AliESDEvent ;
class AliAODPWG4Particle ;

class AliAnaPi0 : public AliAnaCaloTrackCorrBaseClass {
  
 public:   
  
  AliAnaPi0() ;
    
  virtual ~AliAnaPi0() ;
  
  //-------------------------------
  // General analysis frame methods
  //-------------------------------

  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects(); 
  
  void         Print(const Option_t * opt) const;
  
  void         MakeAnalysisFillHistograms();
  
  void         InitParameters();
  
  //-------------------------------
  // EVENT Bin Methods
  //-------------------------------

  Int_t        GetEventIndex(AliAODPWG4Particle * part, Double_t * vert)  ;  

  //-------------------------------
  // Opening angle pair selection
  //-------------------------------
    
  void         SwitchOnAngleSelection()         { fUseAngleCut         = kTRUE  ; }
  void         SwitchOffAngleSelection()        { fUseAngleCut         = kFALSE ; }
  
  void         SwitchOnAngleEDepSelection()     { fUseAngleEDepCut     = kTRUE  ; }
  void         SwitchOffAngleEDepSelection()    { fUseAngleEDepCut     = kFALSE ; }
    
  void         SetAngleCut(Float_t a)           { fAngleCut            = a      ; }
  void         SetAngleMaxCut(Float_t a)        { fAngleMaxCut         = a      ; }

  void         SwitchOnFillAngleHisto()         { fFillAngleHisto      = kTRUE  ; }
  void         SwitchOffFillAngleHisto()        { fFillAngleHisto      = kFALSE ; }
  
  //------------------------------------------
  // Do analysis only with clusters in same SM or different combinations of SM
  //------------------------------------------
    
  void         SwitchOnSameSM()                 { fSameSM              = kTRUE  ; }
  void         SwitchOffSameSM()                { fSameSM              = kFALSE ; }
  
  void         SwitchOnSMCombinations()         { fFillSMCombinations  = kTRUE  ; }
  void         SwitchOffSMCombinations()        { fFillSMCombinations  = kFALSE ; }
  
  //-------------------------------
  // Histogram filling options off by default
  //-------------------------------
    
  void         SwitchOnInvPtWeight()            { fMakeInvPtPlots      = kTRUE  ; }
  void         SwitchOffInvPtWeight()           { fMakeInvPtPlots      = kFALSE ; }
  
  void         SwitchOnFillBadDistHisto()       { fFillBadDistHisto    = kTRUE  ; }
  void         SwitchOffFillBadDistHisto()      { fFillBadDistHisto    = kFALSE ; }

  void         SwitchOnFillOnlyMCAcceptanceHisto()  { fFillOnlyMCAcceptanceHisto    = kTRUE  ; }
  void         SwitchOffFillOnlyMCAcceptanceHisto() { fFillOnlyMCAcceptanceHisto    = kFALSE ; }
  
  //-------------------------------------------
  // Cuts for multiple analysis, off by default
  //-------------------------------------------
    
  void         SwitchOnMultipleCutAnalysis()    { fMultiCutAna         = kTRUE  ; }
  void         SwitchOffMultipleCutAnalysis()   { fMultiCutAna         = kFALSE ; }
  
  void         SwitchOnMultipleCutAcceptAnalysis()  { fMultiCutAnaAcc  = kTRUE  ; }
  void         SwitchOffMultipleCutAcceptAnalysis() { fMultiCutAnaAcc  = kFALSE ; }

  void         SetNPtCuts   (Int_t s)           { if(s <= 10)fNPtCuts    = s    ; }
  void         SetNAsymCuts (Int_t s)           { if(s <= 10)fNAsymCuts  = s    ; }
  void         SetNNCellCuts(Int_t s)           { if(s <= 10)fNCellNCuts = s    ; }
  void         SetNPIDBits  (Int_t s)           { if(s <= 10)fNPIDBits   = s    ; }
  void         SetNAngleCutBins(Int_t s)        { if(s <= 10)fNAngleCutBins = s ; }
  
  void         SetPtCutsAt   (Int_t p,Float_t v){ if(p < 11)fPtCuts[p]   = v    ; }
  void         SetPtCutsMaxAt(Int_t p,Float_t v){ if(p < 11)fPtCutsMax[p]= v    ; }
  void         SetAsymCutsAt (Int_t p,Float_t v){ if(p < 10)fAsymCuts[p] = v    ; }
  void         SetNCellCutsAt(Int_t p,Int_t v)  { if(p < 10)fCellNCuts[p]= v    ; }
  void         SetPIDBitsAt  (Int_t p,Int_t v)  { if(p < 10)fPIDBits[p]  = v    ; }
  void         SetAngleCutBin(Int_t p,Float_t v){ if(p < 10)fAngleCutBinsArray[p]= v; }
  
  void         SwitchOnFillSSCombinations()     { fFillSSCombinations  = kTRUE  ; }
  void         SwitchOffFillSSCombinations()    { fFillSSCombinations  = kFALSE ; }
  
  void         SwitchOnFillAsymmetryHisto()     { fFillAsymmetryHisto  = kTRUE  ; }
  void         SwitchOffFillAsymmetryHisto()    { fFillAsymmetryHisto  = kFALSE ; }

  void         SwitchOnFillOriginHisto()        { fFillOriginHisto     = kTRUE  ; }
  void         SwitchOffFillOriginHisto()       { fFillOriginHisto     = kFALSE ; }

  void         SwitchOnFillArmenterosThetaStarHisto()  { fFillArmenterosThetaStar = kTRUE  ; }
  void         SwitchOffFillArmenterosThetaStarHisto() { fFillArmenterosThetaStar = kFALSE ; }

  void         SwitchOnFillSecondaryCellTimeSel()      { fFillSecondaryCellTiming = kTRUE  ; }
  void         SwitchOffFillSecondaryCellTimeSel()     { fFillSecondaryCellTiming = kFALSE ; }
  
  void         SwitchOnFillOpAngleCutHisto()    { fFillOpAngleCutHisto = kTRUE  ; }
  void         SwitchOffFillOpAngleCutHisto()   { fFillOpAngleCutHisto = kFALSE ; }
  
  //-------------------------------------------
  // Pair 2 different inputs
  //-------------------------------------------

  void         SwitchOnPairWithOtherDetector() { fPairWithOtherDetector = kTRUE ; }    
  void         SwitchOffPairWithOtherDetector(){ fPairWithOtherDetector = kFALSE; } 
  void         SetOtherDetectorInputName(TString name)
  { fOtherDetectorInputName = name ;   if(name != "") SwitchOnPairWithOtherDetector() ; }

  // MC analysis related methods
    
  void         SwitchOnConversionChecker()      { fCheckConversion     = kTRUE  ; }
  void         SwitchOffConversionChecker()     { fCheckConversion     = kFALSE ; }  
  
  void         SwitchOnMultipleCutAnalysisInSimulation()  { fMultiCutAnaSim = kTRUE  ; }
  void         SwitchOffMultipleCutAnalysisInSimulation() { fMultiCutAnaSim = kFALSE ; }
  
  void         SwitchOnCheckAcceptanceInSector() { fCheckAccInSector   = kTRUE  ; }
  void         SwitchOffCheckAcceptanceInSector(){ fCheckAccInSector   = kFALSE ; }
  
  void         FillAcceptanceHistograms();
    
  void         FillMCVersusRecDataHistograms(Int_t ancLabel ,  Int_t    ancPDG, 
                                             Int_t ancStatus,  Double_t weightPt,
                                             Int_t    iclus1,  Int_t    iclus2,
                                             Int_t    mctag1,  Int_t    mctag2,
                                             Float_t  pt1,     Float_t  pt2,
                                             Int_t    ncells1, Int_t    ncells2,
                                             Double_t mass,    Double_t pt,     Double_t asym,
                                             Double_t deta,    Double_t dphi,   Double_t angle);
  
  void         SetPi0MassSelectionWindow(Float_t min, Float_t max) { fPi0MassWindow[0] = min ; fPi0MassWindow[1] = max ; }
  void         SetEtaMassSelectionWindow(Float_t min, Float_t max) { fEtaMassWindow[0] = min ; fEtaMassWindow[1] = max ; }
  
  void         FillArmenterosThetaStar(Int_t pdg);

  private:

  /// Containers for photons in stored events
  TList ** fEventsList ;               //![GetNCentrBin()*GetNZvertBin()*GetNRPBin()]
  
  Bool_t   fUseAngleCut ;              ///<  Select pairs depending on their opening angle
  Bool_t   fUseAngleEDepCut ;          ///<  Select pairs depending on their opening angle
  Float_t  fAngleCut ;                 ///<  Select pairs with opening angle larger than a threshold
  Float_t  fAngleMaxCut ;              ///<  Select pairs with opening angle smaller than a threshold
  
  Float_t  fPi0MassWindow[2];          ///<  Pi0 mass selection window
  Float_t  fEtaMassWindow[2];          ///<  Eta mass selection window
  
  // Multiple cuts analysis
  Bool_t   fMultiCutAna;               ///<  Do analysis with several or fixed cut
  Bool_t   fMultiCutAnaSim;            ///<  Do analysis with several or fixed cut, in the simulation related part
  Bool_t   fMultiCutAnaAcc;            ///<  Do analysis with several or fixed cut, acceptance plots (eta-phi, col-row)
  Int_t    fNPtCuts;                   ///<  Number of pt cuts
  Float_t  fPtCuts[11];                ///<  Array with different pt cuts, minimum
  Float_t  fPtCutsMax[11];             ///<  Array with different pt cuts, maximum
  Int_t    fNAsymCuts;                 ///<  Number of assymmetry cuts
  Float_t  fAsymCuts[10];              ///<  Array with different assymetry cuts
  Int_t    fNCellNCuts;                ///<  Number of cuts with number of cells in cluster
  Int_t    fCellNCuts[10];             ///<  Array with different cell number cluster cuts
  Int_t    fNPIDBits ;                 ///<  Number of possible PID bit combinations
  Int_t    fPIDBits[10];               ///<  Array with different PID bits

  Int_t    fNAngleCutBins ;            ///<  Number of angle cuts bins
  Float_t  fAngleCutBinsArray[11];     ///<  Array with angle cut bins
  
  // Switchs of different analysis options
  Bool_t   fMakeInvPtPlots;            ///<  Do plots with inverse pt weight
  Bool_t   fSameSM;                    ///<  Select only pairs in same SM;
  Bool_t   fFillSMCombinations;        ///<  Fill histograms with different cluster pairs in SM combinations
  Bool_t   fCheckConversion;           ///<  Fill histograms with tagged photons as conversion
  Bool_t   fFillBadDistHisto;          ///<  Do plots for different distances to bad channels
  Bool_t   fFillSSCombinations;        ///<  Do invariant mass for different combination of shower shape clusters
  Bool_t   fFillAngleHisto;            ///<  Fill histograms with pair opening angle
  Bool_t   fFillAsymmetryHisto;        ///<  Fill histograms with asymmetry vs pt
  Bool_t   fFillOriginHisto;           ///<  Fill histograms depending on their origin
  Bool_t   fFillArmenterosThetaStar;   ///<  Fill armenteros histograms
  Bool_t   fFillOnlyMCAcceptanceHisto; ///<  Do analysis only of MC kinematics input
  Bool_t   fFillSecondaryCellTiming;   ///<  Fill histograms depending on timing of secondary cells in clusters
  Bool_t   fFillOpAngleCutHisto;       ///<  Fill histograms depending on opening angle of pair
  
  Bool_t   fCheckAccInSector;          ///<  Check that the decay pi0 falls in the same SM or sector
  
  Bool_t   fPairWithOtherDetector;     ///<  Pair (DCal and PHOS) or (PCM and (PHOS or DCAL or EMCAL))
  TString  fOtherDetectorInputName;    ///<  String with name of extra detector data
  
  TLorentzVector fPhotonMom1;          //!<! Photon cluster momentum, temporary array
  TLorentzVector fPhotonMom1Boost;     //!<! Photon cluster momentum, temporary array
  TLorentzVector fPhotonMom2;          //!<! Photon cluster momentum, temporary array
  TLorentzVector fMCPrimMesonMom;      //!<! Pi0/Eta MC primary momentum, temporary array
  TVector3       fMCProdVertex;        //!<! Pi0/Eta MC Production vertex, temporary array
    
  // ----------
  // Histograms
  // ----------
        
  /// REAL two-photon invariant mass distribution for different calorimeter modules.
  TH2F **  fhReMod ;                   //![fNModules]
    
  /// REAL two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhReSameSideEMCALMod ;      //![fNModules-2]
    
  /// REAL two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhReSameSectorEMCALMod ;    //![fNModules/2]
    
  /// REAL  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhReDiffPHOSMod ;           //![fNModules]
    
  /// REAL two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhReSameSectorDCALPHOSMod ;    //![6]
  
  /// REAL  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhReDiffSectorDCALPHOSMod ;    //![8]
  
  /// MIXED two-photon invariant mass distribution for different calorimeter modules.
  TH2F **  fhMiMod ;                   //![fNModules]
    
  /// MIXED  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhMiSameSideEMCALMod ;      //![fNModules-2]
    
  /// MIXED two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhMiSameSectorEMCALMod ;    //![fNModules/2]
    
  /// MIXED two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhMiDiffPHOSMod ;           //![fNModules-1]
  
  /// MIXED two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhMiSameSectorDCALPHOSMod ;    //![6]
  
  /// MIXED  two-photon invariant mass distribution for different clusters in different calorimeter modules.
  TH2F **  fhMiDiffSectorDCALPHOSMod ;    //![8]

  // Pairs with at least one cluster tagged as conversion
    
  TH2F *   fhReConv ;                  //!<! REAL  two-photon invariant mass distribution one of the pair was 2 clusters with small mass 
  TH2F *   fhMiConv ;                  //!<! MIXED two-photon invariant mass distribution one of the pair was 2 clusters with small mass
  TH2F *   fhReConv2 ;                 //!<! REAL  two-photon invariant mass distribution both pair photons recombined from 2 clusters with small mass 
  TH2F *   fhMiConv2 ;                 //!<! MIXED two-photon invariant mass distribution both pair photons recombined from 2 clusters with small mass

  /// REAL two-photon invariant mass distribution for different centralities and Asymmetry.
  TH2F **  fhRe1 ;                     //![GetNCentrBin()*fNPIDBits*fNAsymCuts]
    
  /// MIXED two-photon invariant mass distribution for different centralities and Asymmetry.
  TH2F **  fhMi1 ;                     //![GetNCentrBin()*fNPIDBits*fNAsymCuts]
    
   /// REAL  two-photon invariant mass distribution for different centralities and Asymmetry.
   /// Apply strict cut on distance to bad channel than in fhRe1.
  TH2F **  fhRe2 ;                     //![GetNCentrBin()*fNPIDBits*fNAsymCuts]
    
  /// MIXED two-photon invariant mass distribution for different centralities and Asymmetry.
  /// Apply strict cut on distance to bad channel than in fhMi1.
  TH2F **  fhMi2 ;                     //![GetNCentrBin()*fNPIDBits*fNAsymCuts]
    
  /// REAL two-photon invariant mass distribution for different centralities and Asymmetry.
  /// Apply strict cut on distance to bad channel than in fhRe2.
  TH2F **  fhRe3 ;                     //![GetNCentrBin()*fNPIDBits*fNAsymCuts]

  /// MIXED two-photon invariant mass distribution for different centralities and Asymmetry.
  /// Apply strict cut on distance to bad channel than in fhMi2.
  TH2F **  fhMi3 ;                     //![GetNCentrBin()*fNPIDBits*fNAsymCuts]

  // Histograms weighted by inverse pT
    
  /// REAL two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT.
  TH2F **  fhReInvPt1 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts]
    
  /// MIXED two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT.
  TH2F **  fhMiInvPt1 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts]
    
  /// REAL two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT.
  /// Apply strict cut on distance to bad channel than in fhRe1InvPt1.
  TH2F **  fhReInvPt2 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts]
  
  /// MIXED two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT.
  /// Apply strict cut on distance to bad channel than in fhMi1InvPt1.
  TH2F **  fhMiInvPt2 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts]

  /// REAL two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT.
  /// Apply strict cut on distance to bad channel than in fhRe1InvPt2.
  TH2F **  fhReInvPt3 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts]
    
  /// MIXED two-photon invariant mass distribution for different centralities and Asymmetry, inverse pT.
  /// Apply strict cut on distance to bad channel than in fhRe1InvPt2.
  TH2F **  fhMiInvPt3 ;                //![GetNCentrBin()*fNPIDBits*fNAsymCuts]
  
  // Multiple cuts: Assymmetry, pt, n cells, PID
    
  /// REAL two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry
  TH2F **  fhRePtNCellAsymCuts ;       //![fNPtCuts*fNAsymCuts*fNCellNCuts]
    
  /// Mixed two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry.
  TH2F **  fhMiPtNCellAsymCuts ;       //![fNPtCuts*fNAsymCuts*fNCellNCuts]
    
  /// REAL two-photon invariant mass distribution for different pt cut, n cell cuts and assymetry for each module.
  TH2F **  fhRePtNCellAsymCutsSM[20] ; //![fNPtCuts*fNAsymCuts*fNCellNCuts]

  /// REAL two-photon opening angle distribution for different pt cut, n cell cuts and assymetry
  TH2F **  fhRePtNCellAsymCutsOpAngle ;       //![fNPtCuts*fNAsymCuts*fNCellNCuts]
  
  /// Mixed two-photon opening angle distribution for different pt cut, n cell cuts and assymetry.
  TH2F **  fhMiPtNCellAsymCutsOpAngle ;       //![fNPtCuts*fNAsymCuts*fNCellNCuts]
  
  /// REAL two-photon opening angle distribution for different pt cut, n cell cuts and assymetry for each module.
  TH2F **  fhRePtNCellAsymCutsSMOpAngle[20] ; //![fNPtCuts*fNAsymCuts*fNCellNCuts]

  TH2F *   fhReSS[3] ;                 //!<! Combine clusters with 3 different cuts on shower shape
    
  // Asymmetry vs pt, in pi0/eta regions
    
  TH2F *   fhRePtAsym    ;             //!<! REAL two-photon pt vs asymmetry
  TH2F *   fhRePtAsymPi0 ;             //!<! REAL two-photon pt vs asymmetry, close to pi0 mass
  TH2F *   fhRePtAsymEta ;             //!<! REAL two-photon pt vs asymmetry, close to eta mass
  TH2F *   fhMiPtAsym    ;             //!<! Mix two-photon pt vs asymmetry
  TH2F *   fhMiPtAsymPi0 ;             //!<! Mix two-photon pt vs asymmetry, close to pi0 mass
  TH2F *   fhMiPtAsymEta ;             //!<! Mix two-photon pt vs asymmetry, close to eta mass
  
  // Centrality, Event plane bins
    
  TH1I *   fhEventBin;                 //!<! Number of real  pairs in a particular bin (cen,vz,rp)
  TH1I *   fhEventMixBin;              //!<! Number of mixed pairs in a particular bin (cen,vz,rp)
  TH1F *   fhCentrality;               //!<! Histogram with centrality bins with at least one pare
  TH1F *   fhCentralityNoPair;         //!<! Histogram with centrality bins with no pair

  TH2F *   fhEventPlaneResolution;     //!<! Histogram with Event plane resolution vs centrality
  
  // Pair opening angle
    
  TH2F *   fhRealOpeningAngle ;        //!<! Opening angle of pair versus pair energy
  TH2F *   fhRealCosOpeningAngle ;     //!<! Cosinus of opening angle of pair version pair energy
  TH2F *   fhMixedOpeningAngle ;       //!<! Opening angle of pair versus pair energy
  TH2F *   fhMixedCosOpeningAngle ;    //!<! Cosinus of opening angle of pair version pair energy
  
  TH2F *   fhRealOpeningAnglePerSM [20]; //!<! Opening angle of pair versus pair energy, per SM
  TH2F *   fhMixedOpeningAnglePerSM[20]; //!<! Opening angle of pair versus pair energy, per SM
  
  // MC analysis histograms
  // Pi0 Acceptance
    
  TH1F *   fhPrimPi0E ;                //!<! Spectrum of Primary
  TH1F *   fhPrimPi0Pt ;               //!<! Spectrum of Primary
  TH1F *   fhPrimPi0PtInCalo ;         //!<! Spectrum of Primary, meson in calo acceptance
  TH1F *   fhPrimPi0AccE ;             //!<! Spectrum of primary with accepted daughters
  TH1F *   fhPrimPi0AccPt ;            //!<! Spectrum of primary with accepted daughters
  TH1F *   fhPrimPi0AccPtPhotonCuts ;  //!<! Spectrum of primary with accepted daughters, photon pt or angle cuts
  TH1F *   fhPrimPi0AccPtOpAngCuts[10];//!<! Spectrum of primary with accepted daughters, different opening angles
  TH2F *   fhPrimPi0Y ;                //!<! Rapidity distribution of primary particles  vs pT
  TH2F *   fhPrimPi0AccY ;             //!<! Rapidity distribution of primary with accepted daughters  vs pT
  TH2F *   fhPrimPi0Yeta ;             //!<! PseudoRapidity distribution of primary particles  vs pT
  TH2F *   fhPrimPi0YetaYcut ;         //!<! PseudoRapidity distribution of primary particles  vs pT, Y<1
  TH2F *   fhPrimPi0AccYeta ;          //!<! PseudoRapidity distribution of primary with accepted daughters  vs pT
  TH2F *   fhPrimPi0Phi ;              //!<! Azimutal distribution of primary particles  vs pT
  TH2F *   fhPrimPi0AccPhi;            //!<! Azimutal distribution of primary with accepted daughters  vs pT
  TH2F *   fhPrimPi0OpeningAngle ;     //!<! Opening angle of pair versus pair energy, primaries
  TH2F *   fhPrimPi0OpeningAnglePhotonCuts ; //!<! Opening angle of pair versus pair energy, primaries, photon pt cuts
  TH2F *   fhPrimPi0OpeningAngleAsym ; //!<! Opening angle of pair versus pair E asymmetry, pi0 primaries
  TH2F *   fhPrimPi0CosOpeningAngle ;  //!<! Cosinus of opening angle of pair version pair energy, pi0 primaries
  TH2F *   fhPrimPi0PtCentrality ;     //!<! primary pi0 reconstructed centrality  vs pT
  TH2F *   fhPrimPi0PtEventPlane ;     //!<! primary pi0 reconstructed event plane vs pT
  TH2F *   fhPrimPi0AccPtCentrality ;  //!<! primary pi0 with accepted daughters reconstructed centrality  vs pT
  TH2F *   fhPrimPi0AccPtEventPlane ;  //!<! primary pi0 with accepted daughters reconstructed event plane vs pT

  // Eta acceptance
    
  TH1F *   fhPrimEtaE ;                //!<! Spectrum of Primary
  TH1F *   fhPrimEtaPt ;               //!<! Spectrum of Primary
  TH1F *   fhPrimEtaPtInCalo ;         //!<! Spectrum of Primary, meson in calo acceptance
  TH1F *   fhPrimEtaAccE ;             //!<! Spectrum of primary with accepted daughters
  TH1F *   fhPrimEtaAccPt ;            //!<! Spectrum of primary with accepted daughters
  TH1F *   fhPrimEtaAccPtPhotonCuts ;  //!<! Spectrum of primary with accepted daughters, photon pt or angle cuts
  TH1F *   fhPrimEtaAccPtOpAngCuts[10];//!<! Spectrum of primary with accepted daughters, different opening angles
  TH2F *   fhPrimEtaY ;                //!<! Rapidity distribution of primary particles vs pT
  TH2F *   fhPrimEtaAccY ;             //!<! Rapidity distribution of primary with accepted daughters  vs pT
  TH2F *   fhPrimEtaYeta ;             //!<! PseudoRapidity distribution of primary particles vs pT
  TH2F *   fhPrimEtaYetaYcut ;         //!<! PseudoRapidity distribution of primary particles vs pT, Y<1
  TH2F *   fhPrimEtaAccYeta ;          //!<! PseudoRapidity distribution of primary with accepted daughters  vs pT
  TH2F *   fhPrimEtaPhi ;              //!<! Azimutal distribution of primary particles  vs pT
  TH2F *   fhPrimEtaAccPhi;            //!<! Azimutal distribution of primary with accepted daughters	 vs pT
  TH2F *   fhPrimEtaOpeningAngle ;     //!<! Opening angle of pair versus pair energy, eta primaries
  TH2F *   fhPrimEtaOpeningAnglePhotonCuts ; //!<! Opening angle of pair versus pair energy, eta primaries, photon pT cuts
  TH2F *   fhPrimEtaOpeningAngleAsym ; //!<! Opening angle of pair versus pair E asymmetry, eta primaries
  TH2F *   fhPrimEtaCosOpeningAngle ;  //!<! Cosinus of opening angle of pair version pair energy, eta primaries
  TH2F *   fhPrimEtaPtCentrality ;     //!<! primary eta reconstructed centrality  vs pT
  TH2F *   fhPrimEtaPtEventPlane ;     //!<! primary eta reconstructed event plane vs pT
  TH2F *   fhPrimEtaAccPtCentrality ;  //!<! primary eta with accepted daughters reconstructed centrality  vs pT
  TH2F *   fhPrimEtaAccPtEventPlane ;  //!<! primary eta with accepted daughters reconstructed event plane vs pT
  
  // Primaries origin
    
  TH2F *   fhPrimPi0PtOrigin ;         //!<! Spectrum of generated pi0 vs mother
  TH2F *   fhPrimEtaPtOrigin ;         //!<! Spectrum of generated eta vs mother
  TH2F *   fhPrimNotResonancePi0PtOrigin ; //!<! Spectrum of generated pi0 vs mother
  TH2F *   fhPrimPi0PtStatus ;         //!<! Spectrum of generated pi0 vs pi0 status

  // Per Generator in Cocktail
  TH1F *   fhPrimPi0PtPerGenerator[10] ;               //!<! Spectrum of primary with |y| < 1
  TH1F *   fhPrimPi0PtInCaloPerGenerator[10] ;         //!<! Spectrum of primary with pi0 in calo
  TH1F *   fhPrimPi0AccPtPerGenerator[10] ;            //!<! Spectrum of primary with accepted daughters
  TH1F *   fhPrimPi0AccPtPhotonCutsPerGenerator[10] ;  //!<! Spectrum of primary with accepted daughters, photon pt or angle cuts
  TH2F *   fhPrimPi0PhiPerGenerator[10] ;              //!<! Azimutal distribution of primary particles  vs pT
  TH2F *   fhPrimPi0YPerGenerator[10] ;                //!<! Rapidity distribution of primary particles  vs pT

  TH1F *   fhPrimEtaPtPerGenerator[10] ;               //!<! Spectrum of primary with |y| < 1
  TH1F *   fhPrimEtaPtInCaloPerGenerator[10] ;         //!<! Spectrum of primary with eta in calo
  TH1F *   fhPrimEtaAccPtPerGenerator[10] ;            //!<! Spectrum of primary with accepted daughters
  TH1F *   fhPrimEtaAccPtPhotonCutsPerGenerator[10] ;  //!<! Spectrum of primary with accepted daughters, photon pt or angle cuts
  TH2F *   fhPrimEtaPhiPerGenerator[10] ;              //!<! Azimutal distribution of primary particles  vs pT
  TH2F *   fhPrimEtaYPerGenerator[10] ;                //!<! Rapidity distribution of primary particles  vs pT
  
  // Pair origin
  // Array of histograms ordered as follows: 0-Photon, 1-electron, 2-pi0, 3-eta, 4-a-proton, 5-a-neutron, 6-stable particles,
  // 7-other decays, 8-string, 9-final parton, 10-initial parton, intermediate, 11-colliding proton, 12-unrelated
    
  TH2F *   fhMCOrgMass[13];            //!<! Mass vs pt of real pairs, check common origin of pair
  TH2F *   fhMCOrgAsym[13];            //!<! Asymmetry vs pt of real pairs, check common origin of pair
  TH2F *   fhMCOrgDeltaEta[13];        //!<! Delta Eta vs pt of real pairs, check common origin of pair
  TH2F *   fhMCOrgDeltaPhi[13];        //!<! Delta Phi vs pt of real pairs, check common origin of pair
  
  // Multiple cuts in simulation, origin pi0 or eta
    
  /// Real pi0 pairs, reconstructed mass vs reconstructed pt of original pair.
  TH2F **  fhMCPi0MassPtRec;           //![fNPtCuts*fNAsymCuts*fNCellNCuts]
  
  /// Real pi0 pairs, reconstructed mass vs generated pt of original pair.
  TH2F **  fhMCPi0MassPtTrue;          //![fNPtCuts*fNAsymCuts*fNCellNCuts]
    
  /// Real pi0 pairs, reconstructed pt vs generated pt of pair.
  TH2F **  fhMCPi0PtTruePtRec;         //![fNPtCuts*fNAsymCuts*fNCellNCuts]
  
  /// Real pi0 pairs, reconstructed pt vs generated pt of pair, apply cut on pi0 mass
  TH2F **  fhMCPi0PtTruePtRecMassCut;  //![fNPtCuts*fNAsymCuts*fNCellNCuts]
  
  /// Real eta pairs, reconstructed mass vs reconstructed pt of original pair.
  TH2F **  fhMCEtaMassPtRec;           //![fNPtCuts*fNAsymCuts*fNCellNCuts]
  
  /// Real eta pairs, reconstructed mass vs generated pt of original pair.
  TH2F **  fhMCEtaMassPtTrue;          //![fNPtCuts*fNAsymCuts*fNCellNCuts]
    
  /// Real eta pairs, reconstructed pt vs generated pt of pair.
  TH2F **  fhMCEtaPtTruePtRec;         //![fNPtCuts*fNAsymCuts*fNCellNCuts]

  /// Real eta pairs, reconstructed pt vs generated pt of pair, apply cut on eta mass
  TH2F **  fhMCEtaPtTruePtRecMassCut;  //![fNPtCuts*fNAsymCuts*fNCellNCuts]

  TH2F *   fhMCPi0PerCentrality;       //!<! Real pi0 pairs, reco pT vs centrality 
  TH2F *   fhMCPi0PerCentralityMassCut;//!<! Real pi0 pairs, reco pT vs centrality, mass cut around pi0 
  TH2F *   fhMCEtaPerCentrality;       //!<! Real eta pairs, reco pT vs centrality  
  TH2F *   fhMCEtaPerCentralityMassCut;//!<! Real eta pairs, reco pT vs centrality, mass cut around eta 
  
  TH2F *   fhMCPi0PtTruePtRecRat;      //!<! Real pi0 pairs, reco pT vs pT ratio reco / generated 
  TH2F *   fhMCPi0PtTruePtRecDif;      //!<! Real pi0 pairs, reco pT vs pT difference generated - reco 
  TH2F *   fhMCPi0PtRecOpenAngle;      //!<! Real pi0 pairs, reco pT vs reco opening angle 
  
  TH2F *   fhMCEtaPtTruePtRecRat;      //!<! Real pi0 pairs, reco pT vs pT ratio reco / generated 
  TH2F *   fhMCEtaPtTruePtRecDif;      //!<! Real pi0 pairs, reco pT vs pT difference generated - reco 
  TH2F *   fhMCEtaPtRecOpenAngle;      //!<! Real pi0 pairs, reco pT vs reco opening angle 

  TH2F *   fhMCPi0PtTruePtRecRatMassCut; //!<! Real pi0 pairs, reco pT vs pT ratio reco / generated, inside a mass window 
  TH2F *   fhMCPi0PtTruePtRecDifMassCut; //!<! Real pi0 pairs, reco pT vs pT difference generated - reco, inside a mass window 
  TH2F *   fhMCPi0PtRecOpenAngleMassCut; //!<! Real pi0 pairs, reco pT vs reco opening angle, inside a mass window 
  
  TH2F *   fhMCEtaPtTruePtRecRatMassCut; //!<! Real pi0 pairs, reco pT vs pT ratio reco / generated, inside a mass window 
  TH2F *   fhMCEtaPtTruePtRecDifMassCut; //!<! Real pi0 pairs, reco pT vs pT difference generated - reco, inside a mass window 
  TH2F *   fhMCEtaPtRecOpenAngleMassCut; //!<! Real pi0 pairs, reco pT vs reco opening angle, inside a mass window 
  
  TH2F *   fhMCPi0PtOrigin ;           //!<! Mass of reconstructed pi0 pairs in calorimeter vs mother origin ID.
  TH2F *   fhMCEtaPtOrigin ;           //!<! Mass of reconstructed eta pairs in calorimeter vs mother origin ID.
  TH2F *   fhMCNotResonancePi0PtOrigin;//!<! Mass of reconstructed pi0 pairs in calorimeter vs mother origin ID, pi0 status 1.
  TH2F *   fhMCPi0PtStatus ;           //!<! Mass of reconstructed pi0 pairs in calorimeter vs mother status.

  TH2F *   fhMCPi0ProdVertex;          //!<! Spectrum of selected pi0 vs production vertex
  TH2F *   fhMCEtaProdVertex;          //!<! Spectrum of selected eta vs production vertex
  TH2F *   fhPrimPi0ProdVertex;        //!<! Spectrum of primary pi0 vs production vertex
  TH2F *   fhPrimEtaProdVertex;        //!<! Spectrum of primary eta vs production vertex

  TH2F *   fhMCPi0Radius[10] ;         //!<! reconstructed Pi0 production vertex vs pT for mother origin
  TH2F *   fhMCEtaRadius[6] ;          //!<! reconstructed Eta production vertex vs pT for mother origin

  
  TH2F *   fhReMCFromConversion ;      //!<! Invariant mass of 2 clusters originated in conversions
  TH2F *   fhReMCFromNotConversion ;   //!<! Invariant mass of 2 clusters not originated in conversions
  TH2F *   fhReMCFromMixConversion ;   //!<! Invariant mass of 2 clusters one from conversion and the other not

  TH2F *   fhArmPrimPi0[4];            //!<! Armenteros plots for primary pi0 in 6 energy bins
  TH2F *   fhArmPrimEta[4];            //!<! Armenteros plots for primary eta in 6 energy bins
  TH2F *   fhCosThStarPrimPi0;         //!<! cos(theta*) plots vs E for primary pi0, same as asymmetry ...
  TH2F *   fhCosThStarPrimEta;         //!<! cos(theta*) plots vs E for primary eta, same as asymmetry ...
  
  TH2F *   fhEPairDiffTime;            //!<! E pair vs Pair of clusters time difference vs E
  
  // Select clusters depending on cell time content
  TH2F * fhReSecondaryCellInTimeWindow; //!<! Combine clusters when all significant cells in cluster have t < 50 ns, same event
  TH2F * fhMiSecondaryCellInTimeWindow; //!<! Combine clusters when all significant cells in cluster have t < 50 ns, different events
  TH2F * fhReSecondaryCellOutTimeWindow;//!<! Combine clusters when at least one significant cells in cluster has t > 50 ns, same event
  TH2F * fhMiSecondaryCellOutTimeWindow;//!<! Combine clusters when at least one significant cells in cluster has t > 50 ns, different events
  
  // Cluster pair studies depending on opening angle
  TH2F *  fhReOpAngleBinMinClusterEtaPhi       [10] ; //!<! Eta-Phi location of lowest energy cluster in pair, depending on opening angle cut
  TH2F *  fhReOpAngleBinMaxClusterEtaPhi       [10] ; //!<! Eta-Phi location of highest energy cluster in pair, depending on opening angle cut
  TH2F *  fhReOpAngleBinMinClusterColRow       [10] ; //!<! Column and row location of main cell of lowest energy cluster in pair, depending on opening angle cut
  TH2F *  fhReOpAngleBinMaxClusterColRow       [10] ; //!<! Column and row location of main cell of highest energy cluster in pair, depending on opening angle cut
  TH2F *  fhReOpAngleBinMinClusterEPerSM       [10] ; //!<! energy of lowest energy cluster in pair, depending on opening angle cut, y axis is SM number
  TH2F *  fhReOpAngleBinMaxClusterEPerSM       [10] ; //!<! energy of highest energy cluster in pair, depending on opening angle cut, y axis is SM number    
  TH2F *  fhReOpAngleBinMinClusterTimePerSM    [10] ; //!<! time of lowest energy cluster in pair, depending on opening angle cut, y axis is SM number
  TH2F *  fhReOpAngleBinMaxClusterTimePerSM    [10] ; //!<! time of highest energy cluster in pair, depending on opening angle cut, y axis is SM number    
  TH2F *  fhReOpAngleBinMinClusterNCellPerSM   [10] ; //!<! N cells of lowest energy cluster in pair, depending on opening angle cut, y axis is SM number
  TH2F *  fhReOpAngleBinMaxClusterNCellPerSM   [10] ; //!<! N cells of highest energy cluster in pair, depending on opening angle cut, y axis is SM number  
  TH2F *  fhReOpAngleBinPairClusterRatioPerSM  [10] ; //!<! lowest/highest energy cluster in pair, depending on opening angle cut, y axis is SM number  
  TH2F *  fhReOpAngleBinPairClusterMass        [10] ; //!<! cluster pair mass vs pT, depending on opening angle cut  
  TH2F *  fhReOpAngleBinPairClusterMassPerSM   [10] ; //!<! cluster pair mass, depending on opening angle cut, y axis is SM number  
//  TH2F *  fhReOpAngleBinPairClusterAbsIdMaxCell[10] ; //!<! Cluster cell with maximum energy in one selected photon vs the other 

  TH2F *  fhReOpAngleBinPairClusterMassMCTruePi0[10] ; //!<! cluster pair mass vs pT, depending on opening angle cut, true pi0 decay pairs from MC 
  TH2F *  fhReOpAngleBinPairClusterMassMCTrueEta[10] ; //!<! cluster pair mass vs pT, depending on opening angle cut, true eta decay pairs from MC
  
  TH2F *  fhMiOpAngleBinMinClusterEtaPhi       [10] ; //!<! Eta-Phi location of lowest energy cluster in pair, depending on opening angle cut, mixed event
  TH2F *  fhMiOpAngleBinMaxClusterEtaPhi       [10] ; //!<! Eta-Phi location of highest energy cluster in pair, depending on opening angle cut, mixed event
//TH2F *  fhMiOpAngleBinMinClusterColRow       [10] ; //!<! Column and row location of main cell of lowest energy cluster in pair, depending on opening angle cut, mixed event
//TH2F *  fhMiOpAngleBinMaxClusterColRow       [10] ; //!<! Column and row location of main cell of highest energy cluster in pair, depending on opening angle cut, mixed event
  TH2F *  fhMiOpAngleBinMinClusterEPerSM       [10] ; //!<! energy of lowest energy cluster in pair, depending on opening angle cut, y axis is SM number, mixed event
  TH2F *  fhMiOpAngleBinMaxClusterEPerSM       [10] ; //!<! energy of highest energy cluster in pair, depending on opening angle cut, y axis is SM number, mixed event    
  TH2F *  fhMiOpAngleBinMinClusterTimePerSM    [10] ; //!<! time of lowest energy cluster in pair, depending on opening angle cut, y axis is SM number, mixed event
  TH2F *  fhMiOpAngleBinMaxClusterTimePerSM    [10] ; //!<! time of highest energy cluster in pair, depending on opening angle cut, y axis is SM number, mixed event  
  TH2F *  fhMiOpAngleBinMinClusterNCellPerSM   [10] ; //!<! N cells of lowest energy cluster in pair, depending on opening angle cut, y axis is SM number, mixed event
  TH2F *  fhMiOpAngleBinMaxClusterNCellPerSM   [10] ; //!<! N cells of highest energy cluster in pair, depending on opening angle cut, y axis is SM number, mixed event  
  TH2F *  fhMiOpAngleBinPairClusterRatioPerSM  [10] ; //!<! lowest/highest energy cluster in pair, depending on opening angle cut, y axis is SM number, mixed event  
  TH2F *  fhMiOpAngleBinPairClusterMass        [10] ; //!<! cluster pair mass vs pT, depending on opening angle cut, mixed event  
  TH2F *  fhMiOpAngleBinPairClusterMassPerSM   [10] ; //!<! cluster pair mass, depending on opening angle cut, y axis is SM number, mixed event  
//TH2F *  fhMiOpAngleBinPairClusterAbsIdMaxCell[10] ; //!<! Cluster cell with maximum energy in one selected photon vs the other, mixed event 
  
  TH2F *  fhPtBinClusterEtaPhi                 [10] ; //!<! Eta-Phi location of cluster in different energy bins.
  TH2F *  fhPtBinClusterColRow                 [10] ; //!<! Column and row location of cluster in different energy bins.

  TH2F *  fhPairGeneratorsBkgMass               [10][10]; //!<! Mass for a pair of clusters depending bkg type
  TH2F *  fhPairGeneratorsBkgMassMCPi0          [10][10]; //!<! Mass for a pair of clusters with depending bkg type, pi0 true pairs
  TH2F *  fhPairGeneratorsBkgCentMCPi0          [10][10]; //!<! Centrality for a pair of clusters with depending bkg type, pi0 true pairs
  TH2F *  fhPairGeneratorsBkgCentMCPi0MassCut   [10][10]; //!<! Centrality for a pair of clusters with depending bkg type, pi0 true pairs, mass cut
  TH2F *  fhPairGeneratorsBkgEPrimRecoRatioMCPi0[10][10]; //!<! pT reco / pT primary for a pair of clusters with depending bkg type, pi0 true pairs
  TH2F *  fhPairGeneratorsBkgEPrimRecoDiffMCPi0 [10][10]; //!<! pT reco - pT primary for a pair of clusters with depending bkg type, pi0 true pairs  
  TH2F *  fhPairGeneratorsBkgMassMCEta          [10][10]; //!<! Mass for a pair of clusters with depending bkg type, eta true pairs
  TH2F *  fhPairGeneratorsBkgCentMCEta          [10][10]; //!<! Centrality for a pair of clusters with depending bkg type, eta true pairs
  TH2F *  fhPairGeneratorsBkgCentMCEtaMassCut   [10][10]; //!<! Centrality for a pair of clusters with depending bkg type, eta true pairs, mass cut
  TH2F *  fhPairGeneratorsBkgEPrimRecoRatioMCEta[10][10]; //!<! pT reco / pT primary for a pair of clusters with depending bkg type, eta true pairs
  TH2F *  fhPairGeneratorsBkgEPrimRecoDiffMCEta [10][10]; //!<! pT reco - pT primary for a pair of clusters with depending bkg type, eta true pairs
  
  TH2F *  fhPairGeneratorsBkgEPrimRecoRatioMCPi0MassCut[10][10]; //!<! pT reco / pT primary for a pair of clusters with depending bkg type, pi0 true pairs, pi0 mass window
  TH2F *  fhPairGeneratorsBkgEPrimRecoDiffMCPi0MassCut [10][10]; //!<! pT reco - pT primary for a pair of clusters with depending bkg type, pi0 true pairs, pi0 mass window
  TH2F *  fhPairGeneratorsBkgEPrimRecoRatioMCEtaMassCut[10][10]; //!<! pT reco / pT primary for a pair of clusters with depending bkg type, eta true pairs, eta mass window
  TH2F *  fhPairGeneratorsBkgEPrimRecoDiffMCEtaMassCut [10][10]; //!<! pT reco - pT primary for a pair of clusters with depending bkg type, eta true pairs, eta mass window
  
  /// Copy constructor not implemented.
  AliAnaPi0(              const AliAnaPi0 & api0) ;
   
  /// Assignment operator not implemented.
  AliAnaPi0 & operator = (const AliAnaPi0 & api0) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaPi0,35) ;
  /// \endcond
  
} ;

#endif //ALIANAPI0_H



