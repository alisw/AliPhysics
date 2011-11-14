#ifndef ALIANALYSISTASKEMCALTRIGGERQA_H
#define ALIANALYSISTASKEMCALTRIGGERQA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------// 
//  Fill histograms with basic QA information for EMCAL offline trigger    //
//  Author: Nicolas Arbor (LPSC-Grenoble), Rachid Guernane (LPSC-Grenoble) //
//          Gustavo Conesa Balbastre  (LPSC-Grenoble)                      //
//                                                                         //
//-------------------------------------------------------------------------//

//--- Root ---
class TList;
class TH1F;
class TH2F;
class AliEMCALGeometry;

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
  
  
  void   UserCreateOutputObjects();    // you should create your output objects in that function if possible
  
  void   UserExec(Option_t *option);   // function called for each event
  
  void   SetGeometryName(TString name)  { fGeoName = name ; } 
  
  void   Terminate(Option_t *option);
  
  AliEMCALRecoUtils* GetRecoUtils()     { return fRecoUtils ; }
  
  //Histogram setters
  
  void   SetTRUTotalSignalHistogramsRange(Int_t nbins,  Float_t max) { fNBinsTRUSignal   = nbins; fMaxTRUSignal   = max ; }
  void   SetSTUTotalSignalHistogramsRange(Int_t nbins,  Float_t max) { fNBinsSTUSignal   = nbins; fMaxSTUSignal   = max ; }
  void   SetV0TotalSignalHistogramsRange (Int_t nbins,  Float_t max) { fNBinsV0Signal    = nbins; fMaxV0Signal    = max ; }
  void   SetSTUFEERatioHistogramsRange   (Int_t nbins,  Float_t max) { fNBinsSTUFEERatio = nbins; fMaxSTUFEERatio = max ; }
  void   SetSTUTRURatioHistogramsRange   (Int_t nbins,  Float_t max) { fNBinsSTUTRURatio = nbins; fMaxSTUFEERatio = max ; }
  void   SetClusterEHistogramsRange      (Int_t nbins,  Float_t max) { fNBinsClusterE    = nbins; fMaxClusterE    = max ; }

  
private:
  TList            *fOutputList;     //! Output list
  
  AliEMCALGeometry *fGeometry;       //  Access to EMCAL geometry utils
  TString           fGeoName;        //  Name of geometry used
  
  AliEMCALRecoUtils *fRecoUtils;     //  RecoUtils
  
  TH1F             *fhNEvents;       //! Number of selected events
  TH2F             *fhFORAmp;        //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column
  TH2F             *fhFORAmpL1G;     //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Gamma trigger event
  TH2F             *fhFORAmpL1J;     //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column, with L1 Jet trigger event
  TH2F             *fhL0Amp;         //! FALTRO signal per Row and Column for FOR involves L0 patch
  TH2F             *fhL0AmpL1G;      //! FALTRO signal per Row and Column for FOR involves L0 patch, with L1G trigger event
  TH2F             *fhL0AmpL1J;      //! FALTRO signal per Row and Column for FOR involves L0 patch, with L1J trigger event
  TH2F             *fhL1Amp;         //! STU signal per Row and Column for FOR involves L0 patch
  TH2F             *fhL1GAmp;        //! STU signal per Row and Column for FOR position of L1 Gamma patch (top-left)
  TH2F             *fhL1JAmp;        //! STU signal per Row and Column for FOR position of L1 Jet patch (top-left)
  TH2F             *fhL0Patch;       //! FOR with L0 patch associated
  TH2F             *fhL1GPatch;      //! FOR with L1 Gamma patch associated
  TH2F             *fhL1JPatch;      //! FOR with L1 Jet patch associated
  TH2F             *fhFEESTU;        //! Correlation FEE vs STU
  TH2F             *fhTRUSTU;        //! Correlation TRU vs STU
  TH2I             *fhV0STU;         //! Total signal STU vs V0C+V0S
  TH2I             *fhFullTRUSTU;    //! Total signal STU vs TRU
  TH2I             *fhSTUChecks;     //! Checks STU/TRU link
  TH1F             *fhClusMB;        //! Clusters distribution for MB trigger
  TH1F             *fhClusL0;        //! Clusters distribution for L0 trigger	
  TH1F             *fhClusL1G;       //! Clusters distribution for L1G trigger
  TH1F             *fhClusL1J;       //! Clusters distribution for L1J trigger
  TH1F             *fhClusL1GOnly;   //! Clusters distribution for L1G trigger and not L1J
  TH1F             *fhClusL1JOnly;   //! Clusters distribution for L1J trigger and not L1G
  TH1F             *fhClusMaxMB;      //! Maximum E Cluster per event distribution for MB trigger
  TH1F             *fhClusMaxL0;      //! Maximum E Cluster per event distribution for L0 trigger	
  TH1F             *fhClusMaxL1G;     //! Maximum E Cluster per event distribution for L1G trigger
  TH1F             *fhClusMaxL1J;     //! Maximum E Cluster per event distribution for L1J trigger
  TH1F             *fhClusMaxL1GOnly; //! Maximum E Cluster per event distribution for L1G trigger and not L1J
  TH1F             *fhClusMaxL1JOnly; //! Maximum E Cluster per event distribution for L1J trigger and not L1G
  TH2F             *fhGPMaxVV0TT;     //! V0 signal vs maximum gamma L1 patch
  TH2F             *fhJPMaxVV0TT;     //! V0 signal vs maximum jet L1 patch
  
  // Histograms bins
  
  Int_t             fNBinsSTUSignal   ; // Number of bins for STU total signal histograms
  Float_t           fMaxSTUSignal     ; // Maximum value for TRU total signal histograms
  Int_t             fNBinsTRUSignal   ; // Number of bins for TRU total signal histograms
  Float_t           fMaxTRUSignal     ; // Maximum value for TRU total signal histograms
  Int_t             fNBinsV0Signal    ; // Number of bins for V0 total signal histograms
  Float_t           fMaxV0Signal      ; // Maximum value for V0 total signal histograms
  Int_t             fNBinsSTUFEERatio ; // Number of bins for STU/FEE ratios histograms
  Float_t           fMaxSTUFEERatio   ; // Maximum value for STU/FEE ratios histograms
  Int_t             fNBinsSTUTRURatio ; // Number of bins for STU/TRU ratios histograms
  Float_t           fMaxSTUTRURatio   ; // Maximum value for STU/TRU ratios histograms
  Int_t             fNBinsClusterE    ; // Number of bins for E cluster histograms
  Float_t           fMaxClusterE      ; // Maximum value for E cluster histograms
  
  //Constants needed by the class: EMCAL 
  static const int  fgkFALTRORows = AliEMCALGeoParams::fgkEMCALRows*(AliEMCALGeoParams::fgkEMCALModules-7)/2; 
  // total number of fake altro rows    in EMCAL
  // (ALTRO channels in one SM times 5 SM divided by 2 per FALTRO)
  
  static const int  fgkFALTROCols = AliEMCALGeoParams::fgkEMCALCols;                                          
  // total number of fake altro columns in EMCAL 
  // (ALTRO channels in one SM times 2 SM divided by 2 per FALTRO)
  
  
  AliAnalysisTaskEMCALTriggerQA(const AliAnalysisTaskEMCALTriggerQA&);            //not implemented
  
  AliAnalysisTaskEMCALTriggerQA& operator=(const AliAnalysisTaskEMCALTriggerQA&); //not implemented
  
  ClassDef(AliAnalysisTaskEMCALTriggerQA, 6);   
};

#endif 
