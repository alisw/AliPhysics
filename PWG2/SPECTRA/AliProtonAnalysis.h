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

#include "AliPID.h"
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
  AliProtonAnalysis();
  AliProtonAnalysis(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
		    Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  virtual ~AliProtonAnalysis();

  void SetBaseAnalysis(AliProtonAnalysisBase *baseAnalysis) {
    fProtonAnalysisBase = baseAnalysis;}
  AliProtonAnalysisBase *GetProtonAnalysisBaseObject() {
    return fProtonAnalysisBase;}

  void InitAnalysisHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
			      Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  Bool_t ReadFromFile(const char* filename);
  void Analyze(AliESDEvent *fESD, 
	       const AliESDVertex *vertex);
  void Analyze(AliAODEvent *fAOD);
  void Analyze(AliStack *stack, Bool_t iInclusive);
  
  AliCFContainer *GetProtonContainer() const {return fProtonContainer;}
  AliCFContainer *GetAntiProtonContainer() const {return fAntiProtonContainer;}

  TH2D *GetProtonYPtHistogram() const {return fHistYPtProtons;}
  TH2D *GetAntiProtonYPtHistogram() const {return fHistYPtAntiProtons;}
  TH1D *GetProtonYHistogram();
  TH1D *GetAntiProtonYHistogram();
  TH1D *GetProtonPtHistogram();
  TH1D *GetAntiProtonPtHistogram();
  TH1D *GetProtonCorrectedYHistogram();
  TH1D *GetAntiProtonCorrectedYHistogram();
  TH1D *GetProtonCorrectedPtHistogram();
  TH1D *GetAntiProtonCorrectedPtHistogram();
  
  TH1D *GetYRatioHistogram();
  TH1D *GetYRatioCorrectedHistogram(TH2D *gCorrectionMapProtons,
				    TH2D *gCorrectionMapAntiProtons);
  TH1D *GetPtRatioHistogram();
  TH1D *GetPtRatioCorrectedHistogram(TH2D *gCorrectionMapProtons,
				     TH2D *gCorrectionMapAntiProtons);
  TH1D *GetYAsymmetryHistogram();
  TH1D *GetPtAsymmetryHistogram();

  TH1I *GetEventHistogram() const {return fHistEvents;}

  Int_t   GetNumberOfAnalyzedEvents() {return (Int_t)fHistEvents->GetEntries();} 
  Bool_t  PrintMean(TH1 *hist, Double_t edge);
  Bool_t  PrintYields(TH1 *hist, Double_t edge); 

  //interface to the correction framework
  void Correct(Int_t step);
  Bool_t ReadCorrectionContainer(const char* filename);
  TList *GetCorrectionListProtons2D() const {return fCorrectionListProtons2D;} 
  TList *GetEfficiencyListProtons1D() const {return fEfficiencyListProtons1D;} 
  TList *GetCorrectionListProtons1D() const {return fCorrectionListProtons1D;} 
  TList *GetCorrectionListAntiProtons2D() const {return fCorrectionListAntiProtons2D;} 
  TList *GetEfficiencyListAntiProtons1D() const {return fEfficiencyListAntiProtons1D;} 
  TList *GetCorrectionListAntiProtons1D() const {return fCorrectionListAntiProtons1D;} 
  
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

  ClassDef(AliProtonAnalysis,1);
};

#endif
