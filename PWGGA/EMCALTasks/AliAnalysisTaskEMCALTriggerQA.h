#ifndef ALIANALYSISTASKEMCALTRIGGERQA_H
#define ALIANALYSISTASKEMCALTRIGGERQA_H

// $Id$

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
  
  void   FillClusterHistograms(Int_t triggerNumber, Bool_t maxCluster,
                               Float_t e,Float_t eta,Float_t phi,
                               Float_t ietamax,Float_t iphimax,
                               Float_t centrality, Float_t v0AC);
  
  void   Init() ;

  void   InitHistogramArrays() ;

  void   LocalInit()                     { Init()                       ; }

  void   UserCreateOutputObjects();    
  
  void   UserExec(Option_t *option);   
  
  AliEMCALRecoUtils* GetRecoUtils()      { if(!fRecoUtils) fRecoUtils = new AliEMCALRecoUtils ;
                                           return fRecoUtils            ; }
  
  void   SetEtaPhiEnMin(Float_t en)      { fEtaPhiEnMin       = en      ; }

  // OADB and geometry settings
  
  void   InitGeometry();
  
  void   SetGeometryName(TString name)   { fGeoName           = name    ; }   
  
  void   AccessOADB() ;
  
  void   SwitchOnEMCALOADB()             { fAccessOADB        = kTRUE   ; }
  void   SwitchOffEMCALOADB()            { fAccessOADB        = kFALSE  ; }
  
  void   SetOADBFilePath(TString path)   { fOADBFilePath      = path    ; }
  
  //Histogram setters

  void   SetTRUTotalSignalHistogramsRange(Int_t nbins,  Float_t max) { fNBinsTRUSignal   = nbins; fMaxTRUSignal   = max ; }
  void   SetSTUTotalSignalHistogramsRange(Int_t nbins,  Float_t max) { fNBinsSTUSignal   = nbins; fMaxSTUSignal   = max ; }
  void   SetV0TotalSignalHistogramsRange (Int_t nbins,  Float_t max) { fNBinsV0Signal    = nbins; fMaxV0Signal    = max ; }
  void   SetSTUFEERatioHistogramsRange   (Int_t nbins,  Float_t max) { fNBinsSTUFEERatio = nbins; fMaxSTUFEERatio = max ; }
  void   SetSTUTRURatioHistogramsRange   (Int_t nbins,  Float_t max) { fNBinsSTUTRURatio = nbins; fMaxSTUFEERatio = max ; }
  void   SetClusterEHistogramsRange      (Int_t nbins,  Float_t max) { fNBinsClusterE    = nbins; fMaxClusterE    = max ; }
    
private:
  TList            *fOutputList;      //! Output list
  
  AliEMCALRecoUtils *fRecoUtils;      //  RecoUtils

  Bool_t            fGeoSet  ;        //  Geometry already set
  AliEMCALGeometry *fGeometry;        //  Access to EMCAL geometry utils
  TString           fGeoName;         //  Name of geometry used
  
  Bool_t            fOADBSet ;        //  AODB parameters already set
  Bool_t            fAccessOADB ;     //  Get calibration from OADB for EMCAL
  TString           fOADBFilePath ;   //  Default path $ALICE_ROOT/OADB/EMCAL, if needed change
  
  Int_t             fBitEGA;          //  fBitEGA
  Int_t             fBitEJE;          //  fBitEJE
  
  Float_t           fEtaPhiEnMin;     //  Min energy for Eta/Phi histograms   
  
  TH1F             *fhNEvents;        //! Number of selected events
  TH2F             *fhFORAmp;         //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column
  TH2F             *fhFORAmpL1G;      //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Gamma trigger event
  TH2F             *fhFORAmpL1G2;     //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Gamma2 trigger event
  TH2F             *fhFORAmpL1J;      //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Jet trigger event
  TH2F             *fhFORAmpL1J2;     //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Jet2 trigger event
  TH2F             *fhL0Amp;          //! FALTRO signal per Row and Column for FOR involves L0 patch
  TH2F             *fhL0AmpL1G;       //! FALTRO signal per Row and Column for FOR involves L0 patch, with L1G trigger event
  TH2F             *fhL0AmpL1J;       //! FALTRO signal per Row and Column for FOR involves L0 patch, with L1J trigger event
  TH2F             *fhL1Amp;          //! STU signal per Row and Column for FOR involves L0 patch
  TH2F             *fhL1GAmp;         //! STU signal per Row and Column for FOR position of L1 Gamma patch (top-left)
  TH2F             *fhL1G2Amp;        //! STU signal per Row and Column for FOR position of L1 Gamma2 patch (top-left)
  TH2F             *fhL1JAmp;         //! STU signal per Row and Column for FOR position of L1 Jet patch (top-left)
  TH2F             *fhL1J2Amp;        //! STU signal per Row and Column for FOR position of L1 Jet2 patch (top-left)
  TH2F             *fhL1FOREnergy;    //! STU signal per Row and Column for FOR position vs FOR energy
  
  TH2F             *fhL0Patch;                //! FOR with L0 patch associated
  TH2F             *fhL1GPatch;               //! FOR with L1 Gamma patch associated
  TH2F             *fhL1G2Patch;              //! FOR with L1 Gamma patch associated
  TH2F             *fhL1GPatchNotFake;        //! FOR with L1 Gamma patch associated but no energy in the related cells
  TH2F             *fhL1GPatchFake;           //! FOR with L1 Gamma patch associated
  TH2F             *fhL1GPatchNotAllFake;     //! FOR with at least 1 L1 Gamma patch associated that has energy in the related celles : not a fake event
  TH2F             *fhL1GPatchAllFake;        //! FOR without any L1 Gamma patch associated with energy in the related cells: fake patch
  TH2F             *fhL1GPatchNotAllFakeMax;  //! FOR with at least one L1 Gamma patch associated with energy in the related cell, maximal energy patch : not fake events
  TH2F             *fhL1GPatchAllFakeMax;     //! FOR without any L1 Gamma patch associated with energy in the related cell, maximal energy patch : fake events
  TH1F             *fhL1GPatchNotAllFakeMaxE; //! Energy distrib of FOR for non fake events, patch of maximal energy
  TH1F             *fhL1GPatchAllFakeMaxE;    //! Energy distrib FOR for fake events, patch of maximal energy
  TH1F             *fhL1GPatchNotAllFakeE;	  //! Energy distrib of FOR for non fake events, all patch energy
  TH1F             *fhL1GPatchAllFakeE;       //! Energy distrib of FOR forfake events, all patch energy
  TH1F             *fhL1GPatchFakeE;          //! Energy distrib of FOR for fake events, all patch energy
  TH1F             *fhL1GPatchNotFakeE;       //! Energy distrib of FOR for non fake events, all patch energy
  TH2F             *fhnpatchFake;             //! number of fake patchs per event vs. if all were fakes or not
  TH2F             *fhnpatchNotFake;          //! number of non fake patchs per events vs. if all were fakes or not
  
  TH2F             *fhL1JPatch;       //! FOR with L1 Jet patch associated
  TH2F             *fhL1J2Patch;      //! FOR with L1 Jet patch associated
  TH2F             *fhFEESTU;         //! Correlation FEE vs STU
  TH2F             *fhTRUSTU;         //! Correlation TRU vs STU
  TH2I             *fhV0STU;          //! Total signal STU vs V0C+V0S
  
  TH2F             *fhGPMaxVV0TT;     //! V0 signal vs maximum gamma L1 patch
  TH2F             *fhJPMaxVV0TT;     //! V0 signal vs maximum jet L1 patch
  TProfile2D       *fhFORMeanAmp;     //! Mean FastOR(FEE) signal per Row and Column
  TProfile2D       *fhL0MeanAmp;      //! Mean FastOR(TRU) signal per Row and Column
  TProfile2D       *fhL1MeanAmp;      //! Mean FastOR(STU) signal per Row and Column
  TH1F             *fhV0[10];         //! V0 distribution for a triggered event
  TH2F             *fhL1GPatchMax;    //! FOR of max. amplitude patch with L1 Gamma patch associated
  TH2F             *fhL1G2PatchMax;   //! FOR of max. amplitude patch with L1 Gamma patch associated
  TH2F             *fhL1JPatchMax;    //! FOR of max. amplitude patch with L1 Jet patch associated  
  TH2F             *fhL1J2PatchMax;   //! FOR of max. amplitude patch with L1 Jet patch associated
  
  // Cluster vs trigger histograms
  enum triggerType{kMBTrig = 0, kL0Trig = 1, kL1GammaTrig = 2, kL1GammaTrig2 = 3, kL1JetTrig = 4, kL1JetTrig2 = 5, kL1GammaOnlyTrig = 6, kL1JetOnlyTrig = 7, kCentralTrig = 8, kSemiCentralTrig = 9 };
  
  TH1F             *fhClusMBPure[3];       //! Clusters E distribution for pure MB trigger
  TH1F             *fhClusMaxMBPure[3];    //! Maximum E Cluster per event distribution for pure MB trigger
  
  TH1F             *fhClus[10];            //! Clusters E distribution for a trigger
  TH1F             *fhClusMax[10];         //! Maximum E Cluster per event distribution for MB trigger

  TH2F             *fhClusCen[10];         //! Clusters Centrality vs E distribution for a trigger
  TH2F             *fhClusCenMax[10];      //! Maximum E Cluster  vs Centrality per event distribution for a trigger
  
  TH2F             *fhClusV0[10];          //! Clusters V0 vs E distribution for a trigger
  TH2F             *fhClusV0Max[10];       //! Maximum E Cluster  vs Centrality per event distribution for a trigger

  TH2F             *fhClusEta[10];         //! Clusters eta vs E distribution for a trigger
  TH2F             *fhClusEtaMax[10];      //! Maximum E Cluster  vs Eta per event distribution for a trigger

  TH2F             *fhClusPhi[10];         //! Clusters Phi vs E distribution for a trigger
  TH2F             *fhClusPhiMax[10];      //! Maximum E Cluster  vs Phi per event distribution for a trigger

  TH2F             *fhClusEtaPhiHigh[10];               //! Clusters eta vs phi distribution for a trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCluMax[10];         //! Maximum E Cluster, Phi vs Eta per event distribution for a trigger, energy above 10 GeV
 
  TH2F             *fhClusEtaPhiHighCellMax[10];        //! Clusters maximum energy cell index eta vs phi distribution for MB trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxCluMax[10];  //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for MB trigger, energy above 10 GeV
 
  TH2F             *fhClusEtaPhiLow[10];                //! Clusters eta vs phi distribution for MB trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCluMax[10];          //! Maximum E Cluster, Phi vs Eta per event distribution for MB trigger, energy below 10 GeV
 
  TH2F             *fhClusEtaPhiLowCellMax[10];         //! Clusters maximum energy cell index eta vs phi distribution for MB trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxCluMax[10];   //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for MB trigger, energy below 10 GeV

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

  //Constants needed by the class: EMCAL 
  static const int  fgkFALTRORows = AliEMCALGeoParams::fgkEMCALRows*(AliEMCALGeoParams::fgkEMCALModules-7)/2;   // total number 
  // of fake altro rows    in EMCAL
  // (ALTRO channels in one SM times 5 SM divided by 2 per FALTRO)
  
  static const int  fgkFALTROCols = AliEMCALGeoParams::fgkEMCALCols; // total number of fake altro columns in EMCAL 
  // (ALTRO channels in one SM times 2 SM divided by 2 per FALTRO)
  
  
  AliAnalysisTaskEMCALTriggerQA           (const AliAnalysisTaskEMCALTriggerQA&); // not implemented
  
  AliAnalysisTaskEMCALTriggerQA& operator=(const AliAnalysisTaskEMCALTriggerQA&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALTriggerQA, 11);   
};

#endif 
