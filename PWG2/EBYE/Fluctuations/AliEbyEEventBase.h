#ifndef ALIEBYEEVENTBASE_H
#define ALIEBYEEVENTBASE_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*-------------------------------------------------------------------------
 *                     AliEbyEEventBase Base Class  
 *       This class deals with Selectoin of Events fot MF and CF Analysis
 *                 origin: Satyajit Jena <sjena@cern.ch>
 * 
 *------------------------------------------------------------------------*/

#include "TObject.h"
#include "TString.h"
class TH1F;
class TH2F;
class TList;

#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliPID.h"
#include "AliESDCentrality.h"
#include "AliCentralitySelectionTask.h"
class AliESDEvent;
class AliESDtrack;
class AliESDVertex;
class AliAODTrack;
class AliMCParticle;

class AliEbyEEventBase : public TObject {
 public:
  enum TriggerMode { kMB1 = 0, kMB2, kSPDFASTOR };
  enum AnalysisMode { kInvalid = -1, kTPC = 0, kITS, kTPCnITS, kForward, kGlobal };
  enum CentralityType { kHardFlat = 0, kHardAcc, kAll, kFlat, kAcc};
  //____________________________________________________________________________//
  AliEbyEEventBase();
  virtual ~AliEbyEEventBase();
  //____________________________________________________________________________//
  void SetAnalysisLevel(const char* type) {fAnalysisLevel = type; }
  const char *GetAnalysisLevel() {return fAnalysisLevel.Data();  }
  void SetAnalysisMode(AnalysisMode analysismode) {fAnalysisMode = analysismode;}
  AnalysisMode GetAnalysisMode() const {return fAnalysisMode;}
  void SetDebugMode() {fDebugMode = kTRUE;}
  //____________________________________________________________________________//
  /* void OfflineTriggerInit() {
    kUseOfflineTrigger = kTRUE;
    fPhySel = new AliPhysicsSelection();
    fPhySel->AddBackgroundIdentification(new AliBackgroundSelection());
    fPhySel->SetAnalyzeMC(fAnalysisMC);
  }

  Bool_t IsOfflineTriggerUsed() {return kUseOfflineTrigger;}
  AliPhysicsSelection *GetPhysicsSelectionObject() {return fPhySel;}*/
 //____________________________________________________________________________//

  void SetPhaseSpace(Int_t nBinsX, Double_t gXmin, Double_t gXmax,
		     Int_t nBinsY, Double_t gYmin, Double_t gYmax) {
    fNBinsX = nBinsX; fMinX   = gXmin; fMaxX   = gXmax;
    fNBinsY = nBinsY; fMinY   = gYmin; fMaxY   = gYmax;
  }

  Int_t    GetNBinsX() const {return fNBinsX;}
  Int_t    GetNBinsY() const {return fNBinsY;}
  Double_t GetMinX()   const {return fMinX;}
  Double_t GetMinY()   const {return fMinY;}
  Double_t GetMaxX()   const {return fMaxX;}
  Double_t GetMaxY()   const {return fMaxY;}

  // Bool_t IsInPhaseSpace(AliESDtrack *track);

  //____________________________________________________________________________//

  const  AliESDVertex *GetVertex(AliESDEvent *esd,
				 AnalysisMode mode,
				 Double_t gVx = 100.,
				 Double_t gVy = 100.,
				 Double_t gVz = 100.);
  
  void SetAcceptedVertexDiamond(Double_t gVx, Double_t gVy, Double_t gVz) {
    fVxMax = gVx; 
    fVyMax = gVy; 
    fVzMax = gVz;
  }
  Double_t GetVxMax() const {return fVxMax;}
  Double_t GetVyMax() const {return fVyMax;}
  Double_t GetVzMax() const {return fVzMax;}
  //____________________________________________________________________________//
  
  TList *GetQA() {return fListQA;}
 //____________________________________________________________________________//
  
  void SetCentralityType(CentralityType centralitytype) {fCentralityType = centralitytype;}
  CentralityType GetCentralityType() const {return fCentralityType;}
  void SetCentralityBin(Int_t cBin) { fCentralityBin = cBin;}
  Int_t GetCentrality(AliESDEvent *esd) const;
  void SetCentralityEstimator(const char *estimator) { fCentralityEstimator = estimator;}
  const char *GetCentralityEstimator() {return fCentralityEstimator.Data();}
  const Int_t GetCentralityBin() {return fCentralityBin;}
  void SetCentralityInputFiles(const char * file1, const char * file2) { fFile1 = file1; fFile2 = file2; 
    fCentralityPointer->SetPercentileFile(fFile1);
    fCentralityPointer->SetPercentileFile2(fFile2);
  }
  AliCentralitySelectionTask *GetCentralityObject() {return fCentralityPointer;}
  Int_t FindCentralityESD(Double_t mult) const;
  Int_t FindCentralityMC(Double_t mult) const;
 //____________________________________________________________________________//
 private:

  TString      fAnalysisLevel; //! ESD or AOD or MC
  AnalysisMode fAnalysisMode;  //! Global, TPC, TPCnITS kForward, kGlobal
  Bool_t fDebugMode; 
  AliPhysicsSelection *fPhySel; //! Physics Seletection 
 

  TList *fListQA;   

  Int_t    fNBinsX;  //! number of bins in y or eta
  Double_t fMinX;    //! max value of y or eta
  Double_t fMaxX;    //! max value of y or eta
  Int_t    fNBinsY;  //! number of bins in pT
  Double_t fMinY;    //! Min value of pT
  Double_t fMaxY;    //! Max value of pT

  Double_t fVxMax;  //!vertex x 
  Double_t fVyMax;  //!vertex x 
  Double_t fVzMax;  //!vertex x 
  
  CentralityType fCentralityType; //! All | EquiDivision | Accumulated | 
  Int_t fCentralityBin; //! Centrality Bin
  TString fCentralityEstimator; //! Centrality Estimator
  AliCentralitySelectionTask *fCentralityPointer;
  TString fFile1; //! file used by centrality task. Set here for bookkeeping
  TString fFile2; //! file used by centrality task. Set here for bookkeeping
  //____________________________________________________________________________//
  AliEbyEEventBase(const AliEbyEEventBase&); // Not implemented
  AliEbyEEventBase& operator=(const AliEbyEEventBase&); // Not implemented


  ClassDef(AliEbyEEventBase,1);
};

#endif
