/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/*************************************
 *   AliFlowCommonHistResults:       *
 *   class to organize the common    *
 *   histograms for Flow Analysis    * 
 *                                   * 
 * authors: Naomi van der Kolk       *
 *           (kolk@nikhef.nl)        *  
 *          Raimond Snellings        *
 *           (snelling@nikhef.nl)    * 
 *          Ante Bilandzic           *
 *           (anteb@nikhef.nl)       * 
 * **********************************/
 
#ifndef ALIFLOWCOMMONHISTRESULTS_H
#define ALIFLOWCOMMONHISTRESULTS_H

class TH1F;
class TH1D;
class TCollection;
class TList;
class TBrowser;

class AliFlowCommonHistResults : public TNamed {

 public:
  AliFlowCommonHistResults();                                                                              // default constructor
  AliFlowCommonHistResults(const char *name, const char *title = "AliFlowCommonHist", Int_t harmonic = 2); // constructor
  virtual ~AliFlowCommonHistResults();                                                                     // destructor

  Bool_t IsFolder() const {return kTRUE;};
  void Browse(TBrowser *b); 

  // Make fill methods for reference flow here: (TBI: rename eventually to FillReferenceFlow)
  Bool_t FillIntegratedFlow(Double_t aV, Double_t anError); // fill fHistIntFlow
  Bool_t FillChi(Double_t aChi);                            // fill fHistChi

  // Make fill methods for differential and integrated flow here:
  Bool_t FillIntegratedFlowRP(Double_t aV, Double_t anError);                   // fill fHistIntFlowRP
  Bool_t FillDifferentialFlowPtRP(Int_t aBin, Double_t av, Double_t anError);   // fill fHistDiffFlowPtRP
  Bool_t FillDifferentialFlowEtaRP(Int_t aBin, Double_t av, Double_t anError);  // fill fHistDiffFlowEtaRP 
  Bool_t FillIntegratedFlowPOI(Double_t aV, Double_t anError);                  // fill fHistIntFlowPOI
  Bool_t FillDifferentialFlowPtPOI(Int_t aBin, Double_t av, Double_t anError);  // fill fHistDiffFlowPtPOI
  Bool_t FillDifferentialFlowEtaPOI(Int_t aBin, Double_t av, Double_t anError); // fill fHistDiffFlowEtaPOI   
             
  // Getters:
  // Reference flow:
  TH1D* GetHistChi(){return fHistChi;};
  TH1D* GetHistIntFlow(){return fHistIntFlow;};    
  // RP = Reference Particles:  
  TH1D* GetHistIntFlowRP(){return fHistIntFlowRP;}; 
  TH1D* GetHistDiffFlowPtRP(){return fHistDiffFlowPtRP;}; 
  TH1D* GetHistDiffFlowEtaRP(){return fHistDiffFlowEtaRP;}; 
  // POI = Particles Of Interest:
  TH1D* GetHistIntFlowPOI(){return fHistIntFlowPOI;};
  TH1D* GetHistDiffFlowPtPOI(){return fHistDiffFlowPtPOI;}; 
  TH1D* GetHistDiffFlowEtaPOI(){return fHistDiffFlowEtaPOI;}; 
  
  TList* GetHistList(){return fHistList;};  

  virtual Double_t Merge(TCollection *aList); // merge function
  void Print(Option_t* option = "") const;    // method to print stats

 private:
  AliFlowCommonHistResults(const AliFlowCommonHistResults& aSetOfResultHists);            // copy constructor
  AliFlowCommonHistResults& operator=(const AliFlowCommonHistResults& aSetOfResultHists); // assignment operator
  // Reference flow:
  TH1D* fHistIntFlow; // reference flow
  TH1D* fHistChi;     // resolution
  // RP = Reference Particles:  
  TH1D* fHistIntFlowRP;     // integrated flow of RPs
  TH1D* fHistDiffFlowPtRP;  // differential flow (Pt) of RPs
  TH1D* fHistDiffFlowEtaRP; // differential flow (Eta) of RPs
  // POI = Particles Of Interest:
  TH1D* fHistIntFlowPOI;     // integrated flow of POIs
  TH1D* fHistDiffFlowPtPOI;  // differential flow (Pt) of POIs
  TH1D* fHistDiffFlowEtaPOI; // differential flow (Eta) of POIs    
  
  TList* fHistList; // list to hold all histograms

  ClassDef(AliFlowCommonHistResults,2) // macro for rootcint
};
 
#endif

