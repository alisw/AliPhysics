#ifndef ALIANALYSISTASKADVVQA_H
#define ALIANALYSISTASKADVVQA_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskADVVQA class
//            This task is for QAing the AD data from ESD/AOD
//		with VZERO veto events
//              Origin: Michal Broz
//-----------------------------------------------------------------

class TString;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TF1;
class AliADCalibData;
class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"

class AliESDEvent;


class AliAnalysisTaskADVVQA : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskADVVQA();
  AliAnalysisTaskADVVQA(const char *name);
 ~AliAnalysisTaskADVVQA();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void InitHistos(TList* list,const char* TriggerName);
  void FillHistos(TList* list);

  TH1F* CreateHist1D(const char* name, const char* title, Int_t nBins, Double_t xMin, Double_t xMax,
		      const char* xLabel = NULL, const char* yLabel = NULL);
  TH2F* CreateHist2D(const char* name, const char* title, Int_t nBinsX, Double_t xMin, Double_t xMax,
		      Int_t nBinsY, Double_t yMin, Double_t yMax,
		      const char* xLabel = NULL, const char* yLabel = NULL, const char* zLabel = NULL);
  TH3F* CreateHist3D(const char* name, const char* title, Int_t nBinsX, Double_t xMin, Double_t xMax,
		      Int_t nBinsY, Double_t yMin, Double_t yMax,
		      Int_t nBinsZ, Double_t zMin, Double_t zMax,
		      const char* xLabel = NULL, const char* yLabel = NULL, const char* zLabel = NULL);

private:

  TList		*fListHist;                       //! List of histograms
  TList		*fList_NoVBA_NoVBC;
  TList		*fList_VBA_NoVBC;
  TList		*fList_NoVBA_VBC;
  TList		*fList_VBA_VBC;
  TList		*fList_TVX;
  TList		*fList_UBA_UBC;
  
  Int_t        fRun;
  Int_t        fOldRun;
  void	       SetCalibData();
  AliADCalibData* fCalibData;      // calibration data
  AliAnalysisUtils* fAnalysisUtils;
   
  AliAnalysisTaskADVVQA(const AliAnalysisTaskADVVQA&);            // not implemented
  AliAnalysisTaskADVVQA& operator=(const AliAnalysisTaskADVVQA&); // not implemented
  
  ClassDef(AliAnalysisTaskADVVQA, 5);
};

#endif
