#ifndef ALIANALYSISTASKEMCALTRIGGERQA_H
#define ALIANALYSISTASKEMCALTRIGGERQA_H

//------------------------------------------------------------------------
/// \class AliAnalysisTaskEMCALTriggerQA
/// \ingroup EMCALPerformance 
/// \brief Fill histograms with basic QA information for EMCAL offline trigger.
///
///  Fill histograms with basic QA information for EMCAL offline trigger.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
/// \author Nicolas Arbor <Nicolas.Arbor@cern.ch>, LPSC-IN2P3-CNRS
/// \author Rachid Guernane <Rachid.Guernane@cern.ch>, LPSC-IN2P3-CNRS
///
//------------------------------------------------------------------------//

//--- Root ---
class TList;
class TH1F;
class TH2I;
class TH2F;
class AliEMCALGeometry;
class TProfile2D;

//--- AliRoot ---
class AliEMCALRecoUtils;
#include "AliEMCALGeoParams.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALTriggerQA : public AliAnalysisTaskSE 
{
public:

  AliAnalysisTaskEMCALTriggerQA();                   // default constructor
  
  AliAnalysisTaskEMCALTriggerQA(const char *name);   // named constructor
  
  virtual ~AliAnalysisTaskEMCALTriggerQA() { ; }     // destructor
  
  void   ClusterAnalysis();
  
  void   FillCellMaps();
  
  void   FillTriggerPatchMaps(TString triggerclasses);
  
  void   FillClusterHistograms(Int_t triggerNumber, Bool_t maxCluster,
                               Float_t e,Float_t eta,Float_t phi,
                               Float_t ietamax,Float_t iphimax, Int_t sm,
                               Float_t centrality, Float_t v0AC);
  
  void   FillCorrelationHistograms();
  
  void   FillEventCounterHistogram();
  
  void   FillL1GammaPatchHistograms();
  
  void   FillL1JetPatchHistograms();
  
  void   FillMapHistograms();
  
  void   FillV0Histograms();
  
  void   Init() ;

  void   InitHistogramArrays() ;

  void   InitCellPatchMaps();
  
  void   LocalInit()                     { Init()                       ; }

  void   UserCreateOutputObjects();    
  
  void   UserExec(Option_t *option);   
  
  AliEMCALRecoUtils* GetRecoUtils()      { if(!fRecoUtils) fRecoUtils = new AliEMCALRecoUtils ;
                                           return fRecoUtils            ; }
  
  void   SetEtaPhiEnMin(Float_t en)      { fEtaPhiEnMin       = en      ; }

  void   SetTriggerEventBit(TString list) ;
  
  // OADB and geometry settings
  
  void   InitGeometry();
  
  void   SetGeometryName(TString name)   { fGeoName           = name    ; }
    
  void   SetEventTriggerL1Bit(Int_t ega, Int_t eje)
                                         { fBitEGA   = ega ; fBitEJE = eje; }

  void   AccessOADB() ;
  
  void   SwitchOnEMCALOADB()             { fAccessOADB        = kTRUE   ; }
  void   SwitchOffEMCALOADB()            { fAccessOADB        = kFALSE  ; }

  void   SwitchOnMCData()                { fMCData            = kTRUE   ; }
  void   SwitchOffMCData()               { fMCData            = kFALSE  ; }

  void   SwitchOnV0SignalHistograms()    { fFillV0SigHisto    = kTRUE   ; }
  void   SwitchOffV0SignalHistograms()   { fFillV0SigHisto    = kFALSE  ; }

  void   SwitchOnCentralityHistograms()  { fFillCenHisto      = kTRUE   ; }
  void   SwitchOffCentralityHistograms() { fFillCenHisto      = kFALSE  ; }
  
  void   SwitchOnAliCentrality ()        { fUseAliCentrality  = kTRUE  ; }
  void   SwitchOffAliCentrality()        { fUseAliCentrality  = kFALSE ; }
  
  void   SwitchOnClusterAcceptanceHistograms()  { fFillClusAcceptHisto = kTRUE   ; }
  void   SwitchOffClusterAcceptanceHistograms() { fFillClusAcceptHisto = kFALSE  ; }
  
  void   SetSuperModulesRange(Int_t min, Int_t max) { fFirstSM = min ; if(max < 20) fLastSM = max ; else fLastSM = 19 ; }
  
  void   SetCentralityEstimator(TString est) { fCentEstimator  = est    ; }
  
  void   SetOADBFilePath(TString path)       { fOADBFilePath   = path   ; }
  
  // Histogram setters

  void   SetTRUTotalSignalHistogramsRange(Int_t nbins,  Float_t max) { fNBinsTRUSignal   = nbins; fMaxTRUSignal   = max ; }
  void   SetSTUTotalSignalHistogramsRange(Int_t nbins,  Float_t max) { fNBinsSTUSignal   = nbins; fMaxSTUSignal   = max ; }
  void   SetV0TotalSignalHistogramsRange (Int_t nbins,  Float_t max) { fNBinsV0Signal    = nbins; fMaxV0Signal    = max ; }
  void   SetSTUFEERatioHistogramsRange   (Int_t nbins,  Float_t max) { fNBinsSTUFEERatio = nbins; fMaxSTUFEERatio = max ; }
  void   SetSTUTRURatioHistogramsRange   (Int_t nbins,  Float_t max) { fNBinsSTUTRURatio = nbins; fMaxSTUFEERatio = max ; }
  void   SetClusterEHistogramsRange      (Int_t nbins,  Float_t max) { fNBinsClusterE    = nbins; fMaxClusterE    = max ; }
  void   SetClusterEtaHistogramsRange    (Int_t nbins,  Float_t max) { fNBinsClusterEta  = nbins; fMaxClusterEta  = max ; }
  void   SetClusterPhiHistogramsRange    (Int_t nbins,  Float_t max, Float_t min)
  { fNBinsClusterPhi  = nbins; fMaxClusterPhi  = max ;  fMinClusterPhi = min ; }
  
private:
    
  TList            *fOutputList;      //!<! Output list
  
  AliEMCALRecoUtils *fRecoUtils;      ///<  RecoUtils

  Bool_t            fGeoSet  ;        ///<  Geometry already set
  AliEMCALGeometry *fGeometry;        ///<  Access to EMCAL geometry utils
  TString           fGeoName;         ///<  Name of geometry used
  
  Bool_t            fOADBSet ;        ///<  AODB parameters already set
  Bool_t            fAccessOADB ;     ///<  Get calibration from OADB for EMCAL
  TString           fOADBFilePath ;   ///<  Default path $ALICE_PHYSICS/OADB/EMCAL, if needed change
  
  Int_t             fBitEGA;          ///<  EGA trigger bit
  Int_t             fBitEJE;          ///<  EJE trigger bit
  
  Float_t           fEtaPhiEnMin;     ///<  Min energy for Eta/Phi histograms
  
  Int_t             fSTUTotal;        ///<  Sum of STU time sums
  Float_t           fTRUTotal;        ///<  Sum of TRU amplitudes
  Float_t           fV0Trigger;       ///<  V0 signal from trigger
  Float_t           fV0A;             ///<  V0 A signal
  Float_t           fV0C;             ///<  V0 C signal
  
  Bool_t            fFillV0SigHisto;  ///<  V0 signal creation and fill
  Bool_t            fFillCenHisto;    ///<  Centrality histograms creation and fill
  Bool_t            fFillClusAcceptHisto; ///<  Fill eta/phi distributions
  Bool_t            fMCData;          ///<   Simulation On/Off
  Int_t             fFirstSM;         ///<  Fill SM histograms for SM >= fFirstSM
  Int_t             fLastSM;          ///<  Fill SM histograms for SM <= fLastSM
  
  TString           fCentEstimator;   ///< Centrality estimator string: V0M, TKL, FMD, ZEMvsZDC, ...
  Bool_t            fUseAliCentrality;///< Use the centrality estimator from AliCentrality or AliMultSelection

  // Event by event trigger recognition bit
  Bool_t            fEventMB   ;      ///<  Bit for MB events
  Bool_t            fEventL0   ;      ///<  Bit for L0 events
  Bool_t            fEventL0D  ;      ///<  Bit for L0 events, DCal
  Bool_t            fEventL1G  ;      ///<  Bit for L1 Gamma 1 events
  Bool_t            fEventL1GD  ;     ///<  Bit for L1 Gamma 1 events, DCal
  Bool_t            fEventL1G2 ;      ///<  Bit for L1 Gamma 2 events
  Bool_t            fEventL1G2D ;     ///<  Bit for L1 Gamma 2 events, DCal
  Bool_t            fEventL1J  ;      ///<  Bit for L1 Jet 1 events
  Bool_t            fEventL1JD  ;     ///<  Bit for L1 Jet 1 events, DCal
  Bool_t            fEventL1J2 ;      ///<  Bit for L1 JEt 2 events
  Bool_t            fEventL1J2D ;     ///<  Bit for L1 JEt 2 events, DCal

  Bool_t            fEventCen  ;      ///<  Bit for Central events
  Bool_t            fEventSem  ;      ///<  Bit for Semi Central events
  
  TLorentzVector    fMomentum ;       ///< Cluster kinematics, temporal, avoid recreation per event.
    
  // Histograms
  
  TH1F             *fhNEvents;                //!<! Number of selected events
  TH2F             *fhFORAmp;                 //!<! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column
  TH2F             *fhFORAmpL1G;              //!<! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Gamma trigger event
  TH2F             *fhFORAmpL1G2;             //!<! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Gamma2 trigger event
  TH2F             *fhFORAmpL1J;              //!<! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Jet trigger event
  TH2F             *fhFORAmpL1J2;             //!<! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Jet2 trigger event
  TH2F             *fhL0Amp;                  //!<! FALTRO signal per Row and Column for FOR involves L0 patch
  TH2F             *fhL0AmpL1G;               //!<! FALTRO signal per Row and Column for FOR involves L0 patch, with L1G trigger event
  TH2F             *fhL0AmpL1J;               //!<! FALTRO signal per Row and Column for FOR involves L0 patch, with L1J trigger event
  TH2F             *fhL1Amp;                  //!<! STU signal per Row and Column for FOR involves L0 patch
  TH2F             *fhL1GAmp;                 //!<! STU signal per Row and Column for FOR position of L1 Gamma patch (top-left)
  TH2F             *fhL1G2Amp;                //!<! STU signal per Row and Column for FOR position of L1 Gamma2 patch (top-left)
  TH2F             *fhL1JAmp;                 //!<! STU signal per Row and Column for FOR position of L1 Jet patch (top-left)
  TH2F             *fhL1J2Amp;                //!<! STU signal per Row and Column for FOR position of L1 Jet2 patch (top-left)
  TH2F             *fhL1FOREnergy;            //!<! STU signal per Row and Column for FOR position vs FOR energy
  
  TH2F             *fhL0Patch;                //!<! FOR with L0 patch associated
  TH2F             *fhL1GPatch;               //!<! FOR with L1 Gamma patch associated
  TH2F             *fhL1G2Patch;              //!<! FOR with L1 Gamma patch associated
  TH2F             *fhL1GPatchNotFake;        //!<! FOR with L1 Gamma patch associated but no energy in the related cells
  TH2F             *fhL1GPatchFake;           //!<! FOR with L1 Gamma patch associated
  TH2F             *fhL1GPatchNotAllFake;     //!<! FOR with at least 1 L1 Gamma patch associated that has energy in the related celles : not a fake event
  TH2F             *fhL1GPatchAllFake;        //!<! FOR without any L1 Gamma patch associated with energy in the related cells: fake patch
  TH2F             *fhL1GPatchNotAllFakeMax;  //!<! FOR with at least one L1 Gamma patch associated with energy in the related cell, maximal energy patch : not fake events
  TH2F             *fhL1GPatchAllFakeMax;     //!<! FOR without any L1 Gamma patch associated with energy in the related cell, maximal energy patch : fake events
  TH1F             *fhL1GPatchNotAllFakeMaxE; //!<! Energy distrib of FOR for non fake events, patch of maximal energy
  TH1F             *fhL1GPatchAllFakeMaxE;    //!<! Energy distrib FOR for fake events, patch of maximal energy
  TH1F             *fhL1GPatchNotAllFakeE;	  //!<! Energy distrib of FOR for non fake events, all patch energy
  TH1F             *fhL1GPatchAllFakeE;       //!<! Energy distrib of FOR forfake events, all patch energy
  TH1F             *fhL1GPatchFakeE;          //!<! Energy distrib of FOR for fake events, all patch energy
  TH1F             *fhL1GPatchNotFakeE;       //!<! Energy distrib of FOR for non fake events, all patch energy
  TH2F             *fhNPatchFake;             //!<! number of fake patchs per event vs. if all were fakes or not
  TH2F             *fhNPatchNotFake;          //!<! number of non fake patchs per events vs. if all were fakes or not
  
  TH2F             *fhL1JPatch;               //!<! FOR with L1 Jet patch associated
  TH2F             *fhL1J2Patch;              //!<! FOR with L1 Jet patch associated
  TH2F             *fhFEESTU;                 //!<! Correlation FEE vs STU
  TH2F             *fhTRUSTU;                 //!<! Correlation TRU vs STU
  TH2I             *fhV0STU;                  //!<! Total signal STU vs V0C+V0S
  
  TH2F             *fhGPMaxVV0TT;             //!<! V0 signal vs maximum gamma L1 patch
  TH2F             *fhJPMaxVV0TT;             //!<! V0 signal vs maximum jet L1 patch
  TProfile2D       *fhFORMeanAmp;             //!<! Mean FastOR(FEE) signal per Row and Column
  TProfile2D       *fhL0MeanAmp;              //!<! Mean FastOR(TRU) signal per Row and Column
  TProfile2D       *fhL1MeanAmp;              //!<! Mean FastOR(STU) signal per Row and Column
  TH2F             *fhL1GPatchMax;            //!<! FOR of max. amplitude patch with L1 Gamma patch associated
  TH2F             *fhL1G2PatchMax;           //!<! FOR of max. amplitude patch with L1 Gamma patch associated
  TH2F             *fhL1JPatchMax;            //!<! FOR of max. amplitude patch with L1 Jet patch associated
  TH2F             *fhL1J2PatchMax;           //!<! FOR of max. amplitude patch with L1 Jet patch associated
  
  ///< Cluster vs trigger histograms enum with trigger types.
  enum triggerType{ kMBTrig                = 0,  
                    kL0Trig                = 1,  kL0TrigD                = 2,
                    kL1GammaTrig           = 3,  kL1GammaTrigD           = 4, 
                    kL1GammaTrig2          = 5,  kL1GammaTrig2D          = 6,
                    kL1JetTrig             = 7,  kL1JetTrigD             = 8,
                    kL1JetTrig2            = 9,  kL1JetTrig2D            = 10,
                    kL1GammaOnlyTrig       = 11, kL1GammaOnlyTrigD       = 12,  
                    kL1JetOnlyTrig         = 13, kL1JetOnlyTrigD         = 14,
                    kL1Gamma2OnlyGammaTrig = 15, kL1Gamma2OnlyGammaTrigD = 16,  
                    kL1Jet2OnlyJetTrig     = 17, kL1Jet2OnlyJetTrigD     = 18,
                    kL0TrigPureEMC         = 19, kL0TrigPureDMC          = 20,
                    kL1GTrigPureEMC        = 21, kL1GTrigPureDMC         = 22,
                    kL1JTrigPureEMC        = 23, kL1JTrigPureDMC         = 24,
                    kCentralTrig           = 25, kSemiCentralTrig        = 26 };
  
  TH1F             *fhClusMBPure[3];                                 //!<! Clusters E distribution for pure MB trigger
  TH1F             *fhClusMaxMBPure[3];                              //!<! Maximum E Cluster per event distribution for pure MB trigger
  
  static const int  fgkTriggerCombi = 27;                            ///<  Total number of trigger combinations defined above
  
  TH1F             *fhClus   [fgkTriggerCombi];                      //!<! Clusters E distribution for a trigger
  TH1F             *fhClusMax[fgkTriggerCombi];                      //!<! Maximum E Cluster per event distribution for MB trigger

  TH2F             *fhClusCen   [fgkTriggerCombi];                   //!<! Clusters Centrality vs E distribution for a trigger
  TH2F             *fhClusCenMax[fgkTriggerCombi];                   //!<! Maximum E Cluster  vs Centrality per event distribution for a trigger
  
  TH2F             *fhClusV0   [fgkTriggerCombi];                    //!<! Clusters V0 vs E distribution for a trigger
  TH2F             *fhClusV0Max[fgkTriggerCombi];                    //!<! Maximum E Cluster  vs Centrality per event distribution for a trigger

  TH1F             *fhClusSM   [fgkTriggerCombi][20];                //!<! Clusters E distribution for a trigger, per SM
  TH2F             *fhClusCenSM[fgkTriggerCombi][20];                //!<! Clusters Centrality vs E distribution for a trigger, per SM
  TH2F             *fhClusV0SM [fgkTriggerCombi][20];                //!<! Clusters V0 vs E distribution for a trigger, per SM
  
  TH2F             *fhClusEta   [fgkTriggerCombi];                   //!<! Clusters eta vs E distribution for a trigger
  TH2F             *fhClusEtaMax[fgkTriggerCombi];                   //!<! Maximum E Cluster  vs Eta per event distribution for a trigger

  TH2F             *fhClusPhi   [fgkTriggerCombi];                   //!<! Clusters Phi vs E distribution for a trigger
  TH2F             *fhClusPhiMax[fgkTriggerCombi];                   //!<! Maximum E Cluster  vs Phi per event distribution for a trigger

  TH2F             *fhClusEtaPhiHigh      [fgkTriggerCombi];         //!<! Clusters eta vs phi distribution for a trigger, energy above fEtaPhiEnMin GeV
  TH2F             *fhClusEtaPhiHighCluMax[fgkTriggerCombi];         //!<! Maximum E Cluster, Phi vs Eta per event distribution for a trigger, energy above fEtaPhiEnMin GeV
 
  TH2F             *fhClusEtaPhiHighCellMax      [fgkTriggerCombi];  //!<! Clusters maximum energy cell index eta vs phi distribution for a trigger, energy above fEtaPhiEnMin GeV
  TH2F             *fhClusEtaPhiHighCellMaxCluMax[fgkTriggerCombi];  //!<! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for MB trigger, energy above fEtaPhiEnMin GeV
 
  TH2F             *fhClusEtaPhiLow      [fgkTriggerCombi];          //!<! Clusters eta vs phi distribution for a trigger, energy below fEtaPhiEnMin GeV
  TH2F             *fhClusEtaPhiLowCluMax[fgkTriggerCombi];          //!<! Maximum E Cluster, Phi vs Eta per event distribution for MB trigger, energy below fEtaPhiEnMin GeV
 
  TH2F             *fhClusEtaPhiLowCellMax      [fgkTriggerCombi];   //!<! Clusters maximum energy cell index eta vs phi distribution for a trigger, energy below fEtaPhiEnMin GeV
  TH2F             *fhClusEtaPhiLowCellMaxCluMax[fgkTriggerCombi];   //!<! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for MB trigger, energy below fEtaPhiEnMin GeV

  TH1F             *fhV0[fgkTriggerCombi];                           //!<! V0 distribution for a triggered event
  TH1F             *fhCentrality[fgkTriggerCombi];                   //!<! Centrality distribution for a triggered event

  // Histograms bins
  
  Int_t             fNBinsSTUSignal   ;     // Number of bins for STU total signal histograms
  Float_t           fMaxSTUSignal     ;     // Maximum value for TRU total signal histograms
  Int_t             fNBinsTRUSignal   ;     // Number of bins for TRU total signal histograms
  Float_t           fMaxTRUSignal     ;     // Maximum value for TRU total signal histograms
  Int_t             fNBinsV0Signal    ;     // Number of bins for V0 total signal histograms
  Float_t           fMaxV0Signal      ;     // Maximum value for V0 total signal histograms
  Int_t             fNBinsSTUFEERatio ;     // Number of bins for STU/FEE ratios histograms
  Float_t           fMaxSTUFEERatio   ;     // Maximum value for STU/FEE ratios histograms
  Int_t             fNBinsSTUTRURatio ;     // Number of bins for STU/TRU ratios histograms
  Float_t           fMaxSTUTRURatio   ;     // Maximum value for STU/TRU ratios histograms
  Int_t             fNBinsClusterE    ;     // Number of bins for E cluster histograms
  Float_t           fMaxClusterE      ;     // Maximum value for E cluster histograms
  Int_t             fNBinsClusterPhi  ;     // Number of bins for Phi cluster histograms
  Float_t           fMaxClusterPhi    ;     // Maximum value for Phi cluster histograms
  Float_t           fMinClusterPhi    ;     // Maximum value for Phi cluster histograms
  Int_t             fNBinsClusterEta  ;     // Number of bins for Eta cluster histograms
  Float_t           fMaxClusterEta    ;     // Maximum value for Eta cluster histograms
    
  //static const int  fgkFALTRORows = AliEMCALGeoParams::fgkEMCALRows*(AliEMCALGeoParams::fgkEMCALModules-7)/2;   // total number
    
  /// Total number of fake altro rows in EMCAL, temporary, not considers DCal yet (ALTRO channels in one SM times 5 SM divided by 2 per FALTRO)
  static const int  fgkFALTRORows = 104; // 60 //AliEMCALGeoParams::fgkEMCALSTURows-4;
  
  /// Total number of fake altro collumns in EMCAL,  (ALTRO channels in one SM times 2 SM divided by 2 per FALTRO)
  static const int  fgkFALTROCols = AliEMCALGeoParams::fgkEMCALSTUCols;
  
  // cell, patch maps
  Double_t fMapCell     [fgkFALTRORows][fgkFALTROCols]; // Cell map
  Double_t fMapCellL1G  [fgkFALTRORows][fgkFALTROCols]; // Cell map for L1G
  Double_t fMapCellL1G2 [fgkFALTRORows][fgkFALTROCols]; // Cell map for L1G2
  Double_t fMapCellL1J  [fgkFALTRORows][fgkFALTROCols]; // Cell map for L1J
  Double_t fMapCellL1J2 [fgkFALTRORows][fgkFALTROCols]; // Cell map for L1J2
  Double_t fMapTrigL0   [fgkFALTRORows][fgkFALTROCols]; // Patch map for L0
  Double_t fMapTrigL1   [fgkFALTRORows][fgkFALTROCols]; // Patch map for L1
  Double_t fMapTrigL0L1G[fgkFALTRORows][fgkFALTROCols]; // Patch map for L0L1G
  Double_t fMapTrigL0L1J[fgkFALTRORows][fgkFALTROCols]; // Patch map for L0L1J
  Double_t fMapTrigL1G  [fgkFALTRORows][fgkFALTROCols]; // Patch map for L1G
  Double_t fMapTrigL1G2 [fgkFALTRORows][fgkFALTROCols]; // Patch map for L1G2
  Double_t fMapTrigL1J  [fgkFALTRORows][fgkFALTROCols]; // Patch map for L1J
  Double_t fMapTrigL1J2 [fgkFALTRORows][fgkFALTROCols]; // Patch map for L1J2

  /// Copy constructor not implemented.
  AliAnalysisTaskEMCALTriggerQA           (const AliAnalysisTaskEMCALTriggerQA&) ;
  
  /// Assignment operator not implemented.
  AliAnalysisTaskEMCALTriggerQA& operator=(const AliAnalysisTaskEMCALTriggerQA&) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALTriggerQA, 18) ;
  /// \endcond

};

#endif 
