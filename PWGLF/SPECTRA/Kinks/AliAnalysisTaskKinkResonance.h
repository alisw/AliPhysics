#ifndef ALIANALYSISTASKKINKRESONANCE_H
#define ALIANALYSISTASKKINKRESONANCE_H

/*  See cxx source for full Copyright notice */

//------------------------------------------------------------------------------
//                   class AliAnalysisTaskKinkResonance
//         This task is an example of an analysis task
//        for analysing resonances having one kaon kink
//Author: Paraskevi Ganoti, University of Athens (pganoti@phys.uoa.gr)
//------------------------------------------------------------------------------

class TList;
class AliResonanceKink;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskKinkResonance : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskKinkResonance(const char *dname = "AliAnalysisTaskKinkResonance");
  virtual ~AliAnalysisTaskKinkResonance() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetAnalysisKinkObject(AliResonanceKink * const kinkResonance) {
    fKinkResonance=kinkResonance;}
  
 private:
  TList             *fList; // List 
  AliResonanceKink  *fKinkResonance; // resonance object configured externaly
  
  AliAnalysisTaskKinkResonance(const AliAnalysisTaskKinkResonance&); // not implemented
  AliAnalysisTaskKinkResonance& operator=(const AliAnalysisTaskKinkResonance&); // not implemented

  ClassDef(AliAnalysisTaskKinkResonance, 1); // example of analysis
};

#endif
