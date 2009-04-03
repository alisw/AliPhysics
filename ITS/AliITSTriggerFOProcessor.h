#ifndef ALIITS_TRIGGERFOPROCESSOR_H
#define ALIITS_TRIGGERFOPROCESSOR_H

/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class takes care of simulating the output from the pixel   //
// trigger system.                                                 //
// The fast-or signals are given as input and they are processed   //
// to see which algorithm will give a signal to the central        //
// trigger system. To avoid doing the same operations many times,  //
// there is a method called PreprocessFOSignals, which should      //
// always be called for each event before the processing of        //
// each individual algorithm is done.                              //
//                                                                 //
// As soon as a new algorithm has been implemented online, a       //
// corresponding method should be implemented here. Which method   //
// should be used for a given trigger label is taken care of in    //
// ProcessFOSignalsLabel method.                                   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliITSFOSignalsSPD.h"

class AliITSTriggerConditions;

class AliITSTriggerFOProcessor {

 public:
  AliITSTriggerFOProcessor();
  AliITSTriggerFOProcessor(AliITSTriggerConditions* ocdbCond);
  AliITSTriggerFOProcessor(const AliITSTriggerFOProcessor& handle);
  virtual ~AliITSTriggerFOProcessor();
  AliITSTriggerFOProcessor& operator=(const AliITSTriggerFOProcessor& handle);

  virtual void          SetTriggerConditions(AliITSTriggerConditions* ocdbCond);
  virtual Bool_t        TriggerConditionsSet() {return fTriggerCond!=NULL;}

  virtual UInt_t        GetNumOutputs() const;
  virtual const Char_t* GetOutputLabel(Short_t index) const;


  virtual void          PreprocessFOSignals(AliITSFOSignalsSPD* signals); // NB! Call this before calling the below methods!!!

  virtual Bool_t        ProcessFOSignalsIndex(Short_t index, AliITSFOSignalsSPD* signals);
  virtual Bool_t        ProcessFOSignalsLabel(const Char_t* label, AliITSFOSignalsSPD* signals);
  
  virtual Bool_t        ProcessFOSignalsTHRTotal(Short_t index, AliITSFOSignalsSPD* /*signals*/);
  virtual Bool_t        ProcessFOSignalsTHRTotalAndTHRInnerAndTHROuter(Short_t index, AliITSFOSignalsSPD* /*signals*/);
  virtual Bool_t        ProcessFOSignalsTHRInnerAndTHROuter(Short_t index, AliITSFOSignalsSPD* /*signals*/);
  virtual Bool_t        ProcessFOSignalsInnerGTOuterPlusOffset(Short_t index, AliITSFOSignalsSPD* /*signals*/);
  virtual Bool_t        ProcessFOSignalsOuterGTInnerPlusOffset(Short_t index, AliITSFOSignalsSPD* /*signals*/);

  
 protected:
  Int_t fFOInner;
  Int_t fFOOuter;
  AliITSTriggerConditions *fTriggerCond;

};

#endif
