#ifndef ALIEBYEMULTIPLICITYFLUCTUATIONANALYSIS_H
#define ALIEBYEMULTIPLICITYFLUCTUATIONANALYSIS_H

/*------------------------------------------------------------------------
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * 
 *-----------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 *
 *             AliEbyEMultiplicityFluctuationAnalysis base class
 *           This class deals with the Multiplicity   fluctuation 
 *                Origin: Satyajit Jena <sjena@cern.ch>
 *
 ------------------------------------------------------------------------*/

#include "TObject.h"
#include "TH1I.h"

class TF1;
class TH2D;
class TH1F;
class TList;
class TTree;
class AliAODEvent;
class AliAODtrack;
class AliESDEvent;
class AliESDtrack;
class AliStack;
class AliESDVertex;
class AliESDtrackCuts;
class AliEbyEEventBase;


class AliEbyEMultiplicityFluctuationAnalysis : public TObject {
 public:
  
  AliEbyEMultiplicityFluctuationAnalysis();
  virtual ~AliEbyEMultiplicityFluctuationAnalysis();

  void SetBaseAnalysis(AliEbyEEventBase * const baseAnalysis) {
    fEbyEBase = baseAnalysis;}


  AliEbyEEventBase *GetEbyEEventBaseObject() const {
    return fEbyEBase;}

  void InitHistos();
  
  void Analyze(AliESDEvent *esd);
  void Analyze(AliAODEvent *aod);
  void Analyze(AliStack *stack);
  
  void Calculate(TObjArray *gTrackArray, Int_t cent);
  TList *GetListMFQA() const {return fListMFQA;}
  TList *GetListMeasureMF() const {return fListMeasureMF;}
  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
    fESDtrackCuts = trackCuts;}
  
  void SetTreeMode(Bool_t istree) { fIsTreeMode = istree;}
  
  TTree *GetFluctuationTree() const {return fFluctuationTree;}

 private:
  TList *fListMFQA;           //! Global list
  TList *fListMeasureMF;      //! List Of Mesures for Multiplisity Fluctuation
  AliEbyEEventBase *fEbyEBase;//! EbyE Events base
  AliESDtrackCuts *fESDtrackCuts;  
  Bool_t fIsTreeMode;
  Int_t iNpossitive;
  Int_t iNnegative;
  Int_t iCentrality;
  TTree *fFluctuationTree;

AliEbyEMultiplicityFluctuationAnalysis(const AliEbyEMultiplicityFluctuationAnalysis&);            //! Not implemented
AliEbyEMultiplicityFluctuationAnalysis& operator=(const AliEbyEMultiplicityFluctuationAnalysis&); //! Not implemented

ClassDef(AliEbyEMultiplicityFluctuationAnalysis,1);
};

#endif
