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

// Global variables:
const Int_t gEventHistogramsEbE = 7; // total number of event histograms => TBI use the enum trick instead! 

// Enums:
enum eEventHistogramsEbE { eNumberOfEvents, eTotalMultiplicity, eSelectedParticles, eCentrality, eVertex_x, eVertex_y, eVertex_z };
enum eBeforeAfterEbE { eBefore = 0, eAfter = 1 };
enum eRecSimEbE { eRec = 0, eSim = 1 };
enum eMinMaxEbE { eMin = 0, eMax = 1 };
enum eDefaultColorsEbE { COLOREbE = kBlack, FILLCOLOREbE = kGray };

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
  virtual void BookEventHistograms();
  virtual void BookFinalResultsHistograms();

  // 2.) Methods called in UserExec(Option_t *):
  virtual void RandomIndices(AliVEvent *ave);
  virtual void ResetEventByEventQuantities();

  // Utility:
  void Red(const char* text);
  void Green(const char* text);
  void Yellow(const char* text);
  void Blue(const char* text);
  TObject* GetObjectFromList(TList *list, Char_t *objectName); // see .cxx
  Int_t NumberOfNonEmptyLines(const char *externalFile);  

  // *) Setters for event histograms:
  void SetVerbose(Bool_t v) {this->fVerbose = v;};
  void SetBookEventHistograms(const char* type, Bool_t bookOrNot)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("NumberOfEvents")){var = eNumberOfEvents;} 
   else if (TString(type).EqualTo("TotalMultiplicity")){var = eTotalMultiplicity;}
   else if (TString(type).EqualTo("SelectedParticles")){var = eSelectedParticles;}
   else if (TString(type).EqualTo("Centrality")){var = eCentrality;}
   else if (TString(type).EqualTo("Vertex_x")){var = eVertex_x;}
   else if (TString(type).EqualTo("Vertex_y")){var = eVertex_y;}
   else if (TString(type).EqualTo("Vertex_z")){var = eVertex_z;}
   else{exit(1);}
   this->fBookEventHistograms[var] = bookOrNot;
  }
  void SetEventHistogramsBins(const char* type, const Double_t nbins, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("NumberOfEvents")){var = eNumberOfEvents;} 
   else if (TString(type).EqualTo("TotalMultiplicity")){var = eTotalMultiplicity;}
   else if (TString(type).EqualTo("SelectedParticles")){var = eSelectedParticles;}
   else if (TString(type).EqualTo("Centrality")){var = eCentrality;}
   else if (TString(type).EqualTo("Vertex_x")){var = eVertex_x;}
   else if (TString(type).EqualTo("Vertex_y")){var = eVertex_y;}
   else if (TString(type).EqualTo("Vertex_z")){var = eVertex_z;}
   else{exit(1);}
   this->fEventHistogramsBins[var][0] = nbins;
   this->fEventHistogramsBins[var][1] = min;
   this->fEventHistogramsBins[var][2] = max;
  }
  void SetEventCuts(const char* type, const Double_t min, const Double_t max)
  {
   Int_t var = -44;
   if(TString(type).EqualTo("NumberOfEvents")){var = eNumberOfEvents;} 
   else if (TString(type).EqualTo("TotalMultiplicity")){var = eTotalMultiplicity;}
   else if (TString(type).EqualTo("SelectedParticles")){var = eSelectedParticles;}
   else if (TString(type).EqualTo("Centrality")){var = eCentrality;}
   else if (TString(type).EqualTo("Vertex_x")){var = eVertex_x;}
   else if (TString(type).EqualTo("Vertex_y")){var = eVertex_y;}
   else if (TString(type).EqualTo("Vertex_z")){var = eVertex_z;}
   else{exit(1);}
   this->fEventCuts[var][0] = min;
   this->fEventCuts[var][1] = max;
  }

 private:
  AliAnalysisTaskEbECumulants(const AliAnalysisTaskEbECumulants& aatmpf);
  AliAnalysisTaskEbECumulants& operator=(const AliAnalysisTaskEbECumulants& aatmpf);
  
  // *) Base lists:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)
  Bool_t fUseFisherYates; // use SetUseFisherYates(kTRUE); in the steering macro to randomize particle indices
  TArrayI *fRandomIndices; // array to store random indices obtained from Fisher-Yates algorithm 
  Bool_t fVerbose; // print all additional info like Green(__PRETTY_FUNCTION__); etc.

  // *) Event histograms and cuts: 
  TList *fEventHistogramsList; // base list to hold all event histograms
  TProfile *fEventHistogramsPro; // keeps all flags for event histograms
  TH1D *fEventHistograms[gEventHistogramsEbE][2][2]; //! [ type - see enum eEventHistograms ][reco,sim][before,after event cuts]
  Bool_t fBookEventHistograms[gEventHistogramsEbE]; // book or not this histogram, see SetBookEventHistograms
  Double_t fEventHistogramsBins[gEventHistogramsEbE][3]; // [nBins,min,max]
  Double_t fEventCuts[gEventHistogramsEbE][2]; // [min,max]
  
  // 2.) Final results:
  TList *fFinalResultsList; // list to hold all histograms with final results

  // Increase this counter in each new version:
  ClassDef(AliAnalysisTaskEbECumulants,3);

};

//================================================================================================================

#endif




