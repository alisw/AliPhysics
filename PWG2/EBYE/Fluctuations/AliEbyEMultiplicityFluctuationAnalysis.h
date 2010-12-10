#ifndef ALIEBYEMULTIPLICITYFLUCTUATIONANALYSIS_H
#define ALIEBYEMULTIPLICITYFLUCTUATIONANALYSIS_H


/*------------------------------------------------------------------------
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * 
 *-----------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 *
 *           AliEbyEMultiplicityFluctuationAnalysis base class
 *        This class deals with the multiplicity  fluctuation 
 *              Origin: Satyajit Jena <sjena@cern.ch>
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

class AliEbyEEventBase;

class AliEbyEMultiplicityFluctuationAnalysis : public TObject {
 public:
  
  AliEbyEMultiplicityFluctuationAnalysis();
  virtual ~AliEbyEMultiplicityFluctuationAnalysis();

  void SetBaseAnalysis(AliEbyEEventBase * const baseAnalysis) { fEbyEBase = baseAnalysis;}
  AliEbyEEventBase *GetEbyEEventBaseObject() const { return fEbyEBase;}

  void InitHistos();
  
  void Analyze(AliESDEvent *esd);
  void Analyze(AliAODEvent *aod);
  void Analyze(AliStack *stack);
  
  TList *GetListMFQA() const {return fListMFQA;}
  TList *GetListMeasureMF() {return fListMeasureMF;}
  
 private:

  TList *fListMFQA;       //! Global list
  TList *fListMeasureMF;  //! List Of Mesures for Multiplisity Fluctuation
  
  AliEbyEEventBase *fEbyEBase;//EbyE Events base
  
  //____________________________________________________________________________//
  AliEbyEMultiplicityFluctuationAnalysis(const AliEbyEMultiplicityFluctuationAnalysis&);            //! Not implemented
  AliEbyEMultiplicityFluctuationAnalysis& operator=(const AliEbyEMultiplicityFluctuationAnalysis&); //! Not implemented
  
  ClassDef(AliEbyEMultiplicityFluctuationAnalysis,1);
};

#endif
