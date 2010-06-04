#ifndef ALIPROTONANALYSIS_H
#define ALIPROTONANALYSIS_H

/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                       Class AliProtonAnalysis
//   This is the class for the baryon (proton) analysis
//
//    Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"
#include "TH1I.h"

class TF1;
class TH2D;
class TH1F;
class TList;

//#include "AliPID.h"
class AliPID;
#include "AliCFContainer.h"
class AliCFDataGrid;
class AliAODEvent;
class AliAODtrack;
class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;
class AliStack;
class AliESDVertex;
class AliProtonAnalysisBase;

class AliProtonAnalysis : public TObject {
 public:
  enum {
    kStepIdentified      = 0,
    kStepIsPrimary       = 1,
    kStepSurvived        = 2,
    kStepInPhaseSpace    = 3,
    kNSteps = 4
  };
  AliProtonAnalysis();
  AliProtonAnalysis(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
		    Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  AliProtonAnalysis(Int_t nbinsY, Double_t *gY,
		    Int_t nbinsPt, Double_t *gPt);
  virtual ~AliProtonAnalysis();

  void SetBaseAnalysis(AliProtonAnalysisBase * const baseAnalysis) {
    fProtonAnalysisBase = baseAnalysis;}
  AliProtonAnalysisBase *GetProtonAnalysisBaseObject() const {
    return fProtonAnalysisBase;}

  void InitAnalysisHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
			      Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  void InitAnalysisHistograms(Int_t nbinsY, Double_t *gY,
			      Int_t nbinsPt, Double_t *gPt);
  Bool_t ReadFromFile(const char* filename);
  void Analyze(AliESDEvent *fESD, 
	       const AliESDVertex *vertex);
  void Analyze(AliAODEvent *fAOD);
  void Analyze(AliStack *stack, Bool_t iInclusive);
  
  //QA for real data
  void InitQA();
  void FillQA(AliESDEvent *esd,
	      const AliESDVertex *vertex, 
	      AliESDtrack* track);
  TList *GetQAList() {return fGlobalQAList;}

  AliCFContainer *GetProtonContainer() const {return fProtonContainer;}
  AliCFContainer *GetAntiProtonContainer() const {return fAntiProtonContainer;}

  TH2D *GetProtonYPtHistogram() const {return fHistYPtProtons;}
  TH2D *GetAntiProtonYPtHistogram() const {return fHistYPtAntiProtons;}
  TH2D *GetProtonYPtCorrectedHistogram() const {
    return fHistYPtProtonsCorrected;}
  TH2D *GetAntiProtonCorrectedYPtHistogram() const {
    return fHistYPtAntiProtonsCorrected;}
  TH1D *GetProtonYHistogram();
  TH1D *GetAntiProtonYHistogram();
  TH1D *GetProtonPtHistogram();
  TH1D *GetAntiProtonPtHistogram();
  TH1D *GetProtonCorrectedYHistogram();
  TH1D *GetAntiProtonCorrectedYHistogram();
  TH1D *GetProtonCorrectedPtHistogram();
  TH1D *GetAntiProtonCorrectedPtHistogram();

  TH1F *GetEventStatistics() {return fHistEventStats;}

  TList *GetYRatioHistogramsInPtBins();
  TH1D *GetYRatioHistogram();
  TH1D *GetYRatioCorrectedHistogram();
  TH1D *GetPtRatioHistogram();
  TH1D *GetPtRatioCorrectedHistogram();
  TH1D *GetYAsymmetryHistogram();
  TH1D *GetPtAsymmetryHistogram();
  TH2D *GetProtonsAbsorptionMaps() {return fHistEfficiencyYPtProtons;}
  TH2D *GetAntiProtonsAbsorptionMaps() {return fHistEfficiencyYPtAntiProtons;}

  TH1I *GetEventHistogram() const {return fHistEvents;}

  Int_t   GetNumberOfAnalyzedEvents() {return (Int_t)fHistEventStats->GetBinContent(5);}
  Bool_t  PrintMean(TH1 *hist, Double_t edge);
  Bool_t  PrintYields(TH1 *hist, Double_t edge); 

  //interface to the correction framework
  void Correct();
  void Correct(Int_t step);
  Bool_t ReadCorrectionContainer(const char* filename);

  void SetCorrectionMapForFeedDown(const char* filename);
  TH2D *GetCorrectionMapForFeedDownProtons() {
    return fHistYPtCorrectionForFeedDownProtons;}
  TH2D *GetCorrectionMapForFeedDownAntiProtons() {
    return fHistYPtCorrectionForFeedDownAntiProtons;}

  void SetCorrectionMapForCuts(const char* filename);
  TH2D *GetCorrectionMapForCutsProtons() {
    return fHistYPtCorrectionForCutsProtons;}
  TH2D *GetCorrectionMapForCutsAntiProtons() {
    return fHistYPtCorrectionForCutsAntiProtons;}

  void SetCorrectionMapForTracking(const char* filename);
  TH2D *GetCorrectionMapForTrackingProtons() {
    return fHistYPtCorrectionForTrackingProtons;}
  TH2D *GetCorrectionMapForTrackingAntiProtons() {
    return fHistYPtCorrectionForTrackingAntiProtons;}


  void SetCorrectionMapForSecondaries(const char* filename);
  TH2D *GetCorrectionMapForSecondaries() {
    return fHistYPtCorrectionForSecondaries;}

  void SetCorrectionMapForCrossSection(const char* filename);
  TH2D *GetProtonsCorrectionMapForCrossSection() {
    return fHistCorrectionForCrossSectionYPtProtons;}
  TH2D *GetAntiProtonsCorrectionMapForCrossSection() {
    return fHistCorrectionForCrossSectionYPtAntiProtons;}

  //iStep=0->MC - iStep=1->Acceptance - iStep=2->Reconstruction - iStep=3->PID
  TH1D  *GetUncorrectedProtonYHistogram(Int_t iStep) {return fProtonContainer->ShowProjection(0, iStep);}
  TH1D  *GetUncorrectedProtonPtHistogram(Int_t iStep) {return fProtonContainer->ShowProjection(1, iStep);}
  TH1D  *GetUncorrectedAntiProtonYHistogram(Int_t iStep) {return fAntiProtonContainer->ShowProjection(0, iStep);}
  TH1D  *GetUncorrectedAntiProtonPtHistogram(Int_t iStep) {return fAntiProtonContainer->ShowProjection(1, iStep);}

 private:
  AliProtonAnalysis(const AliProtonAnalysis&); // Not implemented
  AliProtonAnalysis& operator=(const AliProtonAnalysis&); // Not implemented
  
  AliProtonAnalysisBase *fProtonAnalysisBase;//base analysis object

  Int_t fNBinsY; //number of bins in y or eta
  Double_t fMinY, fMaxY; //min & max value of y or eta
  Int_t fNBinsPt;  //number of bins in pT
  Double_t fMinPt, fMaxPt; //min & max value of pT
  
  //Analysis containers
  AliCFContainer *fProtonContainer; //container for protons
  AliCFContainer *fAntiProtonContainer; //container for antiprotons
  TH1I *fHistEvents; //event counter
  TH2D *fHistYPtProtons; //Y-Pt of Protons
  TH2D *fHistYPtAntiProtons; // Y-Pt of Antiprotons
  TH2D *fHistYPtProtonsCorrected; //Y-Pt of Protons (corrected)
  TH2D *fHistYPtAntiProtonsCorrected; // Y-Pt of Antiprotons (corrected)
  TH1F *fHistEventStats;//Event statistics
  TList *fYRatioInPtBinsList;//TList of the eta dependent ratios for each pT bin

  //Corrections
  TList *fEffGridListProtons; //list for the efficiency grid - protons 
  TList *fCorrectionListProtons2D; //list for the 2d corrections 
  TList *fEfficiencyListProtons1D; //list for the 1d efficiencies
  TList *fCorrectionListProtons1D; //list for the 1d corrections 
  TList *fEffGridListAntiProtons; //list for the efficiency grid - antiprotons 
  TList *fCorrectionListAntiProtons2D; //list for the 2d corrections 
  TList *fEfficiencyListAntiProtons1D; //list for the 1d efficiencies
  TList *fCorrectionListAntiProtons1D; //list for the 1d corrections 
  AliCFDataGrid *fCorrectProtons; //corrected data grid for protons
  AliCFDataGrid *fCorrectAntiProtons; //corrected data grid for antiprotons
  TH2D *fHistEfficiencyYPtProtons;//efficiency 2D for the corrections - protons
  TH2D *fHistEfficiencyYPtAntiProtons;//efficiency 2D for the corrections - antiprotons
  TH2D *fHistCorrectionForCrossSectionYPtProtons;//correction for the proper cross-section - 2D protons
  TH2D *fHistCorrectionForCrossSectionYPtAntiProtons;//correction for the proper cross-section - 2D antiprotons
  Bool_t fHistCorrectionForCrossSectionFlag;//correct for cross-section
  TH2D *fHistYPtCorrectionForCutsProtons;//correction factors for the cut efficiency (protons)
  TH2D *fHistYPtCorrectionForCutsAntiProtons;//correction factors for the cut efficiency (antiprotons)
  Bool_t fCorrectForCutsFlag;//correct for the cut efficiency
  TH2D *fHistYPtCorrectionForTrackingProtons;//correction factors for the tracking efficiency (protons)
  TH2D *fHistYPtCorrectionForTrackingAntiProtons;//correction factors for the tracking efficiency (antiprotons)
  Bool_t fCorrectForTrackingFlag;//correct for the tracking efficiency
  TH2D *fHistYPtCorrectionForFeedDownProtons;//correction factors for the feed-down contamination (protons)
  TH2D *fHistYPtCorrectionForFeedDownAntiProtons;//correction factors for the feed-down contamination (antiprotons)
  Bool_t fCorrectForFeedDownFlag;//correct for cut efficiency
  TH2D *fHistYPtCorrectionForSecondaries;//correction factors for the background protons
  Bool_t fCorrectForSecondariesFlag;//correct for secondaries

  //QA lists
  TList *fGlobalQAList; //global list
  TList *fQA2DList; //QA 2D list
  TList *fQAProtonsAcceptedList; //accepted protons
  TList *fQAProtonsRejectedList; //rejected protons
  TList *fQAAntiProtonsAcceptedList; //accepted antiprotons
  TList *fQAAntiProtonsRejectedList; //rejected antiprotons
  Bool_t fInitQAFlag;//Init flag

  ClassDef(AliProtonAnalysis,1);
};

#endif
