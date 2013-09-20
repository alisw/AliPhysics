/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

 /************************************ 
 * flow analysis with multi-particle *
 *           correlations            * 
 *                                   * 
 * author: Ante Bilandzic            * 
 *        (abilandzic@gmail.com)     *
 ************************************/ 

#ifndef ALIFLOWANALYSISWITHMULTIPARTICLECORRELATIONS_H
#define ALIFLOWANALYSISWITHMULTIPARTICLECORRELATIONS_H

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "Riostream.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"

class AliFlowAnalysisWithMultiparticleCorrelations{
 public:
  AliFlowAnalysisWithMultiparticleCorrelations();
  virtual ~AliFlowAnalysisWithMultiparticleCorrelations(); 
  // Member functions are grouped as:
  // 0.) Methods called in the constructor;
  // 1.) Method Init() and methods called in it (!);
  // 2.) Method Make() and methods called in it;
  // 3.) Method Finish() and methods called in it;
  // 4.) Method GetOutputHistograms() and methods called in it;
  // 5.) Setters and getters;
  // 6.) The rest.

  // 0.) Methods called in the constructor:
  virtual void InitializeArraysForControlHistograms(); 

  // 1.) Method Init() and methods called in it (!):
  virtual void Init();
   virtual void CrossCheckSettings();
   virtual void BookAndNestAllLists(); 
   virtual void BookEverythingForControlHistograms();

  // 2.) Method Make() and methods called in it:
  virtual void Make(AliFlowEventSimple *anEvent);
   virtual void FillControlHistograms(AliFlowEventSimple *anEvent);

  // 3.) Method Finish() and methods called in it:
  virtual void Finish();

  // 4.) Method GetOutputHistograms() and methods called in it: 
  virtual void GetOutputHistograms(TList *histList);

  // 5.) Setters and getters:
  //  5.0.) Base:
  void SetHistList(TList* const hlist) {this->fHistList = hlist;} 
  TList* GetHistList() const {return this->fHistList;} 
  //  5.1.) Control histograms:  
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 
  void SetControlHistogramsFlagsPro(TProfile* const chfp) {this->fControlHistogramsFlagsPro = chfp;};
  TProfile* GetControlHistogramsFlagsPro() const {return this->fControlHistogramsFlagsPro;}; 

  // 6.) The rest:
  virtual void WriteHistograms(TString outputFileName);
  virtual void WriteHistograms(TDirectoryFile *outputFileName);
    
 private:
  AliFlowAnalysisWithMultiparticleCorrelations(const AliFlowAnalysisWithMultiparticleCorrelations& afawQc);
  AliFlowAnalysisWithMultiparticleCorrelations& operator=(const AliFlowAnalysisWithMultiparticleCorrelations& afawQc); 
  // Data members are grouped as:
  // 0.) Base;
  // 1.) Control histograms. 

  // 0.) Base:
  TList* fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  // 1.) Control histograms:  
  TList *fControlHistogramsList;        // list to hold all control histograms
  TProfile *fControlHistogramsFlagsPro; // profile to hold all flags for control histograms
  TH1D *fKinematicsHist[2][3];          // [RP,POI][phi,pt,eta] distributions
  TH1D *fMultDistributionsHist[3];      // multiplicity distribution [RP,POI,reference multiplicity]
  TH2D *fMultCorrelationsHist[3];       // [RP vs. POI, RP vs. refMult, POI vs. refMult]  

  ClassDef(AliFlowAnalysisWithMultiparticleCorrelations,0);

};

//================================================================================================================

#endif





