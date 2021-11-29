/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************** 
* TBI document eventually             * 
**************************************/ 

#ifndef ALIANALYSISTASKEBECUMULANTS_H
#define ALIANALYSISTASKEBECUMULANTS_H

#include <AliAnalysisTaskSE.h>
#include <AliAODTrack.h>
#include <AliAODEvent.h>
#include <AliVEvent.h>
#include <TSystem.h>
#include <TH1F.h>

//================================================================================================================

class AliAnalysisTaskEbECumulants : public AliAnalysisTaskSE{
 public:
  
  AliAnalysisTaskEbECumulants();
  AliAnalysisTaskEbECumulants(const char *name);
  virtual ~AliAnalysisTaskEbECumulants(); 
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  
  // 0.) Methods called in the constructor:
  virtual void InitializeArrays();
 
  // 1.) Methods called in UserCreateOutputObjects():
  virtual void BookAndNestAllLists();
  virtual void BookControlHistograms();
  virtual void BookFinalResultsHistograms();

  // 2.) Methods called in UserExec(Option_t *):
  virtual void RandomIndices(AliVEvent *ave);
  virtual void ResetEventByEventQuantities();

  // 3.) Setters and getters:
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 
  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  void SetPtBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbinsPt = nbins;
   this->fMinBinPt = min;
   this->fMaxBinPt = max;
  };

  // Utility:
  void Red(const char* text);
  void Green(const char* text);
  void Yellow(const char* text);
  void Blue(const char* text);
  TObject* GetObjectFromList(TList *list, Char_t *objectName); // see .cxx
  Int_t NumberOfNonEmptyLines(const char *externalFile);  

 private:
  AliAnalysisTaskEbECumulants(const AliAnalysisTaskEbECumulants& aatmpf);
  AliAnalysisTaskEbECumulants& operator=(const AliAnalysisTaskEbECumulants& aatmpf);
  
  // 0.) Base lists:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)
  Bool_t fUseFisherYates; // use SetUseFisherYates(kTRUE); in the steering macro to randomize particle indices
  TArrayI *fRandomIndices; // array to store random indices obtained from Fisher-Yates algorithm 
  
  // 1.) Control histograms:  
  TList *fControlHistogramsList; // list to hold all control histograms
  TH1F *fPtHist;                 // atrack->Pt()
  Int_t fNbinsPt;                // number of bins
  Float_t fMinBinPt;             // min bin
  Float_t fMaxBinPt;             // min bin
  
  // 2.) Final results:
  TList *fFinalResultsList; // list to hold all histograms with final results

  // Increase this counter in each new version:
  ClassDef(AliAnalysisTaskEbECumulants,2);

};

//================================================================================================================

#endif




