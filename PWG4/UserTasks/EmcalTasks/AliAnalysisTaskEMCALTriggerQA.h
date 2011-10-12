#ifndef ALIANALYSISTASKEMCALTRIGGERQA_H
#define ALIANALYSISTASKEMCALTRIGGERQA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------// 
//  Fill histograms with basic QA information for EMCAL offline trigger    //
//  Author: Nicola Arbor (LPSC-Grenoble)                                   //
//          Gustavo Conesa Balbastre  (LPSC-Grenoble)                      //
//                                                                         //
//-------------------------------------------------------------------------//

//--- Root ---
class TList;
class TH1F;
class TH2F;
class AliEMCALGeometry;

//--- AliRoot ---
#include "AliEMCALGeoParams.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALTriggerQA : public AliAnalysisTaskSE 
{
public:
  AliAnalysisTaskEMCALTriggerQA();                   // default constructor
  
  AliAnalysisTaskEMCALTriggerQA(const char *name);   // named constructor
  
  virtual ~AliAnalysisTaskEMCALTriggerQA() {;}       // destructor
  
  
  void   UserCreateOutputObjects();    // you should create your output objects in that function if possible
  
  void   UserExec(Option_t *option);   // function called for each event
  
  void   SetGeometryName(TString name)  { fGeoName = name ; } 
  
  void   Terminate(Option_t *option);
  
  //Histogram setters
  
  void   SetTRUTotalSignalHistogramsRange(Int_t nbins,  Float_t max) { fNBinsTRUSignal   = nbins; fMaxTRUSignal   = max ; }
  void   SetSTUTotalSignalHistogramsRange(Int_t nbins,  Float_t max) { fNBinsSTUSignal   = nbins; fMaxSTUSignal   = max ; }
  void   SetV0TotalSignalHistogramsRange (Int_t nbins,  Float_t max) { fNBinsV0Signal    = nbins; fMaxV0Signal    = max ; }
  void   SetSTUFEERatioHistogramsRange   (Int_t nbins,  Float_t max) { fNBinsSTUFEERatio = nbins; fMaxSTUFEERatio = max ; }
  void   SetSTUTRURatioHistogramsRange   (Int_t nbins,  Float_t max) { fNBinsSTUTRURatio = nbins; fMaxSTUFEERatio = max ; }

  
private:
  TList            *fOutputList;  //! Output list
  
  AliEMCALGeometry *fGeometry;       //  Access to EMCAL geometry utils
  TString           fGeoName;        //  Name of geometry used
  
  TH1F             *fhNEvents;       //! Number of selected events
  TH2F             *fhFORAmp;        //! FEE cells deposited energy, grouped like FastOR 2x2 per Row and Column
  TH2F             *fhL0Amp;         //! FALTRO signal per Row and Column for FOR involves L0 patch
  TH2F             *fhL1GAmp;         //! STU signal per Row and Column for FOR involves in L1 Gamma patch
  TH2F             *fhL1JAmp;         //! STU signal per Row and Column for FOR involves in L1 Jet patch
  TH2F             *fhL0Patch;       //! FOR with L0 patch associated
  TH2F             *fhL1GPatch;      //! FOR with L1 Gamma patch associated
  TH2F             *fhL1JPatch;      //! FOR with L1 Jet patch associated
  TH2F             *fhFEESTU;        //! Correlation FEE vs STU
  TH2F             *fhTRUSTU;        //! Correlation TRU vs STU
  TH2I             *fhV0STU;         //! Total signal STU vs V0C+V0S
  TH2I             *fhFullTRUSTU;    //! Total signal STU vs TRU
  TH2I             *fhSTUChecks;     //! Checks STU/TRU link

  
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
  
  //Constants needed by the class: EMCAL 
  static const int  fgkFALTRORows = AliEMCALGeoParams::fgkEMCALRows*(AliEMCALGeoParams::fgkEMCALModules-7)/2; 
  // total number of fake altro rows    in EMCAL
  // (ALTRO channels in one SM times 5 SM divided by 2 per FALTRO)
  
  static const int  fgkFALTROCols = AliEMCALGeoParams::fgkEMCALCols;                                          
  // total number of fake altro columns in EMCAL 
  // (ALTRO channels in one SM times 2 SM divided by 2 per FALTRO)
  
  
  AliAnalysisTaskEMCALTriggerQA(const AliAnalysisTaskEMCALTriggerQA&);            //not implemented
  
  AliAnalysisTaskEMCALTriggerQA& operator=(const AliAnalysisTaskEMCALTriggerQA&); //not implemented
  
  ClassDef(AliAnalysisTaskEMCALTriggerQA, 1);   
};

#endif 
