/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************** 
* template class for student projects * 
**************************************/ 

#ifndef ALIANALYSISTASKSTUDENTSCM_H		// {Needed with line 12 and the last one to avoid errors from multiple inclusions of the same class (see RC3)}
#define ALIANALYSISTASKSTUDENTSCM_H

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TComplex.h"
#include "TH2F.h"

//================================================================================================================

class AliAnalysisTaskStudentsCM : public AliAnalysisTaskSE{		// {AliAnalysisTaskSE is mandatory for the compiler to know how execute the tasks}
 public:

  // {The six following lines are mandatory}  
  AliAnalysisTaskStudentsCM();				// {Constructors of the class}
  AliAnalysisTaskStudentsCM(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskStudentsCM();		// {Destructor of the class}
  virtual void UserCreateOutputObjects();	// {The three following lines are needed with these exact names to let know to the compiler where the data members are defined, where the tasks to do on the events are and what to do once the run on the events is over}
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  
  // 0.) Methods called in the constructor:
  virtual void InitializeArrays();			// {Every data members are be initialise in the constructor without any problem, except Arrays that need a special function}
 
  // 1.) Methods called in UserCreateOutputObjects():
  virtual void BookAndNestAllLists();		// {The three functions are useful to get a nice organisation in folders and subfolders in the final output root file}
  virtual void BookControlHistograms();
  virtual void BookTestHistograms();		// {TEST: addition of a new subfolder in the output root with personal histograms}
  virtual void BookFinalResultsHistograms();

  // 2.) Methods called in UserExec(Option_t *):
  // ...

  // 3.) Methods called in Terminate():
  // ...

  // 4.) Setters and getters:				// {Used to access the data members and initialise them for example}
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;}
  void SetTestHistogramsList(TList* const thl) {this->fTestHistogramsList = thl;};	// {TEST: setter and getter of the new subfolder}
  TList* GetTestHistogramsList() const {return this->fTestHistogramsList;}
  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  void SetBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbins = nbins;
   this->fMinBin = min;
   this->fMaxBin = max;
  };

 private:
  AliAnalysisTaskStudentsCM(const AliAnalysisTaskStudentsCM& aatmpf);
  AliAnalysisTaskStudentsCM& operator=(const AliAnalysisTaskStudentsCM& aatmpf);
  
  // 0.) Base lists:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  // 1.) Control histograms:  
  TList *fControlHistogramsList; // list to hold all control histograms
  TList *fTestHistogramsList;	 // {TEST: list to hold the personal histograms for playing around ^^}
  TH1F *fPtHistNoCut;			 // Histogram containing the transverse momenta before the application of cuts
  TH1F *fPtHist;                 // atrack->Pt() {Histogram containing the transverse momenta}
  Int_t fNbins;                  // number of bins
  Float_t fMinBin;               // min bin
  Float_t fMaxBin;               // max bin 
  TH1F *fPhiHist;                // atrack->Phi()
  TH1F *fEtaHist;                // atrack->Eta()
  TH1F *fEnergyHist;			 // atrack->E() {TEST, histogram containing the energy}
  TH1F *fEnergyHistNoCut;	     // Histogram containing the energy before the application of cuts
  TH2F *fPtPhiHist;				 // {TEST: 2D histogram containing the azimuthal angle phi as a function of the transverse momentum}
  TH1F *fMassSquareHist;		 // TEST: mass^2 of the particles

  // 2.) Final results:
  TList *fFinalResultsList; // list to hold all histograms with final results
  TH2F *fEtaMassSquareHist;	// TEST: m^2 as a function of eta, offline made

  ClassDef(AliAnalysisTaskStudentsCM,2);	// {Needed to be increase each time a new version that changes the structure of the output file is sent for the daily tag}

};

//================================================================================================================

#endif











