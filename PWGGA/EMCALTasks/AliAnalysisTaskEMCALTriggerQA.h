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
  
  void   Init() ;
  
  void   LocalInit()                     { Init()                       ; }

  void   UserCreateOutputObjects();    
  
  void   UserExec(Option_t *option);   
  
  AliEMCALRecoUtils* GetRecoUtils()      { if(!fRecoUtils) fRecoUtils = new AliEMCALRecoUtils ;
                                           return fRecoUtils            ; }
  
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
  
  TH1F             *fhNEvents;        //! Number of selected events
  TH2F             *fhFORAmp;         //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column
  TH2F             *fhFORAmpL1G;      //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Gamma trigger event
  TH2F             *fhFORAmpL1J;      //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Jet trigger event
  TH2F             *fhL0Amp;          //! FALTRO signal per Row and Column for FOR involves L0 patch
  TH2F             *fhL0AmpL1G;       //! FALTRO signal per Row and Column for FOR involves L0 patch, with L1G trigger event
  TH2F             *fhL0AmpL1J;       //! FALTRO signal per Row and Column for FOR involves L0 patch, with L1J trigger event
  TH2F             *fhL1Amp;          //! STU signal per Row and Column for FOR involves L0 patch
  TH2F             *fhL1GAmp;         //! STU signal per Row and Column for FOR position of L1 Gamma patch (top-left)
  TH2F             *fhL1JAmp;         //! STU signal per Row and Column for FOR position of L1 Jet patch (top-left)
  TH2F             *fhL0Patch;        //! FOR with L0 patch associated
  TH2F             *fhL1GPatch;       //! FOR with L1 Gamma patch associated
  TH2F             *fhL1JPatch;       //! FOR with L1 Jet patch associated
  TH2F             *fhFEESTU;         //! Correlation FEE vs STU
  TH2F             *fhTRUSTU;         //! Correlation TRU vs STU
  TH2I             *fhV0STU;          //! Total signal STU vs V0C+V0S
  
  TH2F             *fhGPMaxVV0TT;     //! V0 signal vs maximum gamma L1 patch
  TH2F             *fhJPMaxVV0TT;     //! V0 signal vs maximum jet L1 patch
  TProfile2D       *fhFORMeanAmp;     //! Mean FastOR(FEE) signal per Row and Column
  TProfile2D       *fhL0MeanAmp;      //! Mean FastOR(TRU) signal per Row and Column
  TProfile2D       *fhL1MeanAmp;      //! Mean FastOR(STU) signal per Row and Column
  TH1F             *fhV0MB;           //! V0 distribution for MB triggered event
  TH1F             *fhV0L1G;          //! V0 distribution for L1G triggered event
  TH1F             *fhV0L1J;          //! V0 distribution for L1J triggered event
  TH2F             *fhL1GPatchMax;    //! FOR of max. amplitude patch with L1 Gamma patch associated
  TH2F             *fhL1JPatchMax;    //! FOR of max. amplitude patch with L1 Jet patch associated  
  
  // Cluster vs trigger histograms
  
  TH1F             *fhClusMB;         //! Clusters distribution for MB trigger
  TH1F             *fhClusMBPure;     //! Clusters distribution for MB trigger
  TH1F             *fhClusL0;         //! Clusters distribution for L0 trigger	
  TH1F             *fhClusL1G;        //! Clusters distribution for L1G trigger
  TH1F             *fhClusL1J;        //! Clusters distribution for L1J trigger
  TH1F             *fhClusL1GOnly;    //! Clusters distribution for L1G trigger and not L1J
  TH1F             *fhClusL1JOnly;    //! Clusters distribution for L1J trigger and not L1G
  TH1F             *fhClusMaxMB;      //! Maximum E Cluster per event distribution for MB trigger
  TH1F             *fhClusMaxMBPure;  //! Maximum E Cluster per event distribution for MB trigger
  TH1F             *fhClusMaxL0;      //! Maximum E Cluster per event distribution for L0 trigger	
  TH1F             *fhClusMaxL1G;     //! Maximum E Cluster per event distribution for L1G trigger
  TH1F             *fhClusMaxL1J;     //! Maximum E Cluster per event distribution for L1J trigger
  TH1F             *fhClusMaxL1GOnly; //! Maximum E Cluster per event distribution for L1G trigger and not L1J
  TH1F             *fhClusMaxL1JOnly; //! Maximum E Cluster per event distribution for L1J trigger and not L1G

  TH2F             *fhClusCenMB;            //! Clusters Centrality vs E distribution for MB trigger
  TH2F             *fhClusCenL0;            //! Clusters Centrality vs E distribution for L0 trigger	
  TH2F             *fhClusCenL1G;           //! Clusters Centrality vs E distribution for L1G trigger
  TH2F             *fhClusCenL1J;           //! Clusters Centrality vs E distribution for L1J trigger
  TH2F             *fhClusCenL1GOnly;       //! Clusters Centrality vs E distribution for L1G trigger and not L1J
  TH2F             *fhClusCenL1JOnly;       //! Clusters Centrality vs E distribution for L1J trigger and not L1G
  TH2F             *fhClusCenMaxMB;         //! Maximum E Cluster  vs Centrality per event distribution for MB trigger
  TH2F             *fhClusCenMaxL0;         //! Maximum E Cluster  vs Centrality  per event distribution for L0 trigger	
  TH2F             *fhClusCenMaxL1G;        //! Maximum E Cluster  vs Centrality  per event distribution for L1G trigger
  TH2F             *fhClusCenMaxL1J;        //! Maximum E Cluster  vs Centrality  per event distribution for L1J trigger
  TH2F             *fhClusCenMaxL1GOnly;    //! Maximum E Cluster  vs Centrality  per event distribution for L1G trigger and not L1J
  TH2F             *fhClusCenMaxL1JOnly;    //! Maximum E Cluster  vs Centrality  per event distribution for L1J trigger and not L1G  
  
  TH2F             *fhClusV0MB;             //! Clusters Centrality vs E distribution for MB trigger
  TH2F             *fhClusV0L0;             //! Clusters Centrality vs E distribution for L0 trigger	
  TH2F             *fhClusV0L1G;            //! Clusters Centrality vs E distribution for L1G trigger
  TH2F             *fhClusV0L1J;            //! Clusters Centrality vs E distribution for L1J trigger
  TH2F             *fhClusV0L1GOnly;        //! Clusters Centrality vs E distribution for L1G trigger and not L1J
  TH2F             *fhClusV0L1JOnly;        //! Clusters Centrality vs E distribution for L1J trigger and not L1G
  TH2F             *fhClusV0MaxMB;          //! Maximum E Cluster  vs Centrality per event distribution for MB trigger
  TH2F             *fhClusV0MaxL0;          //! Maximum E Cluster  vs Centrality  per event distribution for L0 trigger	
  TH2F             *fhClusV0MaxL1G;         //! Maximum E Cluster  vs Centrality  per event distribution for L1G trigger
  TH2F             *fhClusV0MaxL1J;         //! Maximum E Cluster  vs Centrality  per event distribution for L1J trigger
  TH2F             *fhClusV0MaxL1GOnly;     //! Maximum E Cluster  vs Centrality  per event distribution for L1G trigger and not L1J
  TH2F             *fhClusV0MaxL1JOnly;     //! Maximum E Cluster  vs Centrality  per event distribution for L1J trigger and not L1G 
  
  TH2F             *fhClusEtaMB;            //! Clusters eta vs E distribution for MB trigger
  TH2F             *fhClusEtaL0;            //! Clusters eta vs E distribution for L0 trigger	
  TH2F             *fhClusEtaL1G;           //! Clusters eta vs E distribution for L1G trigger
  TH2F             *fhClusEtaL1J;           //! Clusters eta vs E distribution for L1J trigger
  TH2F             *fhClusEtaL1GOnly;       //! Clusters eta vs E distribution for L1G trigger and not L1J
  TH2F             *fhClusEtaL1JOnly;       //! Clusters eta vs E distribution for L1J trigger and not L1G
  TH2F             *fhClusEtaMaxMB;         //! Maximum E Cluster  vs Eta per event distribution for MB trigger
  TH2F             *fhClusEtaMaxL0;         //! Maximum E Cluster  vs Eta  per event distribution for L0 trigger	
  TH2F             *fhClusEtaMaxL1G;        //! Maximum E Cluster  vs Eta  per event distribution for L1G trigger
  TH2F             *fhClusEtaMaxL1J;        //! Maximum E Cluster  vs Eta  per event distribution for L1J trigger
  TH2F             *fhClusEtaMaxL1GOnly;    //! Maximum E Cluster  vs Eta  per event distribution for L1G trigger and not L1J
  TH2F             *fhClusEtaMaxL1JOnly;    //! Maximum E Cluster  vs Eta  per event distribution for L1J trigger and not L1G

  TH2F             *fhClusPhiMB;            //! Clusters Phi vs E distribution for MB trigger
  TH2F             *fhClusPhiL0;            //! Clusters Phi vs E distribution for L0 trigger	
  TH2F             *fhClusPhiL1G;           //! Clusters Phi vs E distribution for L1G trigger
  TH2F             *fhClusPhiL1J;           //! Clusters Phi vs E distribution for L1J trigger
  TH2F             *fhClusPhiL1GOnly;       //! Clusters Phi vs E distribution for L1G trigger and not L1J
  TH2F             *fhClusPhiL1JOnly;       //! Clusters Phi vs E distribution for L1J trigger and not L1G
  TH2F             *fhClusPhiMaxMB;         //! Maximum E Cluster  vs Phi per event distribution for MB trigger
  TH2F             *fhClusPhiMaxL0;         //! Maximum E Cluster  vs Phi  per event distribution for L0 trigger	
  TH2F             *fhClusPhiMaxL1G;        //! Maximum E Cluster  vs Phi  per event distribution for L1G trigger
  TH2F             *fhClusPhiMaxL1J;        //! Maximum E Cluster  vs Phi  per event distribution for L1J trigger
  TH2F             *fhClusPhiMaxL1GOnly;    //! Maximum E Cluster  vs Phi  per event distribution for L1G trigger and not L1J
  TH2F             *fhClusPhiMaxL1JOnly;    //! Maximum E Cluster  vs Phi  per event distribution for L1J trigger and not L1G

  TH2F             *fhClusEtaPhiHighMB;            //! Clusters eta vs phi distribution for MB trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighL0;            //! Clusters eta vs phi distribution for L0 trigger, energy above 10 GeV	
  TH2F             *fhClusEtaPhiHighL1G;           //! Clusters eta vs phi distribution for L1G trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighL1J;           //! Clusters eta vs phi distribution for L1J trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighL1GOnly;       //! Clusters eta vs phi distribution for L1G trigger and not L1J, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighL1JOnly;       //! Clusters eta vs phi distribution for L1J trigger and not L1G, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCluMaxMB;      //! Maximum E Cluster, Phi vs Eta per event distribution for MB trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCluMaxL0;      //! Maximum E Cluster, Phi vs Eta per event distribution for L0 trigger, energy above 10 GeV	
  TH2F             *fhClusEtaPhiHighCluMaxL1G;     //! Maximum E Cluster, Phi vs Eta per event distribution for L1G trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCluMaxL1J;     //! Maximum E Cluster, Phi vs Eta per event distribution for L1J trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCluMaxL1GOnly; //! Maximum E Cluster, Phi vs Eta per event distribution for L1G trigger and not L1J, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCluMaxL1JOnly; //! Maximum E Cluster, Phi vs Eta per event distribution for L1J trigger and not L1G, energy above 10 GeV
  
  TH2F             *fhClusEtaPhiHighCellMaxMB;            //! Clusters maximum energy cell index eta vs phi distribution for MB trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxL0;            //! Clusters maximum energy cell index eta vs phi distribution for L0 trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxL1G;           //! Clusters maximum energy cell index eta vs phi distribution for L1G trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxL1J;           //! Clusters maximum energy cell index eta vs phi distribution for L1J trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxL1GOnly;       //! Clusters maximum energy cell index eta vs phi distribution for L1G trigger and not L1J, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxL1JOnly;       //! Clusters maximum energy cell index eta vs phi distribution for L1J trigger and not L1G, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxCluMaxMB;      //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for MB trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxCluMaxL0;      //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L0 trigger, energy above 10 GeV	
  TH2F             *fhClusEtaPhiHighCellMaxCluMaxL1G;     //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L1G trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxCluMaxL1J;     //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L1J trigger, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxCluMaxL1GOnly; //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L1G trigger and not L1J, energy above 10 GeV
  TH2F             *fhClusEtaPhiHighCellMaxCluMaxL1JOnly; //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L1J trigger and not L1G, energy above 10 GeV
  
  TH2F             *fhClusEtaPhiLowMB;                    //! Clusters eta vs phi distribution for MB trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowL0;                    //! Clusters eta vs phi distribution for L0 trigger, energy below 10 GeV	
  TH2F             *fhClusEtaPhiLowL1G;                   //! Clusters eta vs phi distribution for L1G trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowL1J;                   //! Clusters eta vs phi distribution for L1J trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowL1GOnly;               //! Clusters eta vs phi distribution for L1G trigger and not L1J, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowL1JOnly;               //! Clusters eta vs phi distribution for L1J trigger and not L1G, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCluMaxMB;              //! Maximum E Cluster, Phi vs Eta per event distribution for MB trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCluMaxL0;              //! Maximum E Cluster, Phi vs Eta per event distribution for L0 trigger, energy below 10 GeV	
  TH2F             *fhClusEtaPhiLowCluMaxL1G;             //! Maximum E Cluster, Phi vs Eta per event distribution for L1G trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCluMaxL1J;             //! Maximum E Cluster, Phi vs Eta per event distribution for L1J trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCluMaxL1GOnly;         //! Maximum E Cluster, Phi vs Eta per event distribution for L1G trigger and not L1J, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCluMaxL1JOnly;         //! Maximum E Cluster, Phi vs Eta per event distribution for L1J trigger and not L1G, energy below 10 GeV
  
  TH2F             *fhClusEtaPhiLowCellMaxMB;             //! Clusters maximum energy cell index eta vs phi distribution for MB trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxL0;             //! Clusters maximum energy cell index eta vs phi distribution for L0 trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxL1G;            //! Clusters maximum energy cell index eta vs phi distribution for L1G trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxL1J;            //! Clusters maximum energy cell index eta vs phi distribution for L1J trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxL1GOnly;        //! Clusters maximum energy cell index eta vs phi distribution for L1G trigger and not L1J, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxL1JOnly;        //! Clusters maximum energy cell index eta vs phi distribution for L1J trigger and not L1G, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxCluMaxMB;       //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for MB trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxCluMaxL0;       //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L0 trigger, energy below 10 GeV	
  TH2F             *fhClusEtaPhiLowCellMaxCluMaxL1G;      //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L1G trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxCluMaxL1J;      //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L1J trigger, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxCluMaxL1GOnly;  //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L1G trigger and not L1J, energy below 10 GeV
  TH2F             *fhClusEtaPhiLowCellMaxCluMaxL1JOnly;  //! Maximum E Cluster, maximum energy cell index Phi vs Eta per event distribution for L1J trigger and not L1G, energy below 10 GeV
  
  
  
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
  
  ClassDef(AliAnalysisTaskEMCALTriggerQA, 10);   
};

#endif 
