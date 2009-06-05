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
class AliESDEvent;
class AliMCEvent;
class AliResonanceKink;
class TH1D;

class AliAnalysisTaskKinkResonance : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskKinkResonance();
  AliAnalysisTaskKinkResonance(const char *name);
  virtual ~AliAnalysisTaskKinkResonance() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetAnalysisKinkObject(AliResonanceKink * const kinkResonance) {
    fKinkResonance=kinkResonance;}
  
 private:
  AliESDEvent       *fESD;    //! ESD object
  AliMCEvent        *fmcEventH;
  TList             *fList; //! List 
  AliResonanceKink  *fKinkResonance;
  
  AliAnalysisTaskKinkResonance(const AliAnalysisTaskKinkResonance&); // not implemented
  AliAnalysisTaskKinkResonance& operator=(const AliAnalysisTaskKinkResonance&); // not implemented

  ClassDef(AliAnalysisTaskKinkResonance, 1); // example of analysis
};

#endif
