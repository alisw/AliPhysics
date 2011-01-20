#ifndef ALIEBYECHARGEFLUCTUATIONANALYSIS_H
#define ALIEBYECHARGEFLUCTUATIONANALYSIS_H

/*------------------------------------------------------------------------
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * 
 *-----------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 *
 *             AliEbyEChargeFluctuationAnalysis base class
 *           This class deals with the Charge   fluctuation 
 *                Origin: Satyajit Jena <sjena@cern.ch>
 *
 ------------------------------------------------------------------------*/

#include "TObject.h"
#include "TH1I.h"

class TF1;
class TH2D;
class TH1F;
class TList;
class AliAODEvent;
class AliAODtrack;
class AliESDEvent;
class AliESDtrack;
class AliStack;
class AliESDVertex;

class AliEbyEEventBase;


class AliEbyEChargeFluctuationAnalysis : public TObject {
 public:
  
  AliEbyEChargeFluctuationAnalysis();
  virtual ~AliEbyEChargeFluctuationAnalysis();
  
  void SetBaseAnalysis(AliEbyEEventBase * const baseAnalysis) {
    fEbyEBase = baseAnalysis;}
  AliEbyEEventBase *GetEbyEEventBaseObject() const {
    return fEbyEBase;}
  
  void InitHistos();
  
  void Analyze(AliESDEvent *esd);
  void Analyze(AliAODEvent *aod);
  void Analyze(AliStack *stack);
  
  void Calculate(TObjArray *gTrackArray, Int_t cent);
  TList *GetListCFQA() const {return fListCFQA;}
  TList *GetListMeasureCF() const {return fListMeasureCF;}
    
  void ReadFromFile();
  
 private:
  
  TList *fListCFQA;           //! Global list
  TList *fListMeasureCF;      //! List Of Mesures for Multiplisity Fluctuation
  AliEbyEEventBase *fEbyEBase;//! EbyE Events base
  
  TH1F *fhNpAverage;   
  TH1F *fhNnAverage;   
  TH1F *fhNCAverage;   
  TH1F *fhNetAverage;  
  TH1F *fhNnSNpAverage;

  
  AliEbyEChargeFluctuationAnalysis(const AliEbyEChargeFluctuationAnalysis&);            //! Not implemented
  AliEbyEChargeFluctuationAnalysis& operator=(const AliEbyEChargeFluctuationAnalysis&); //! Not implemented
  
  ClassDef(AliEbyEChargeFluctuationAnalysis,1);
};

#endif
