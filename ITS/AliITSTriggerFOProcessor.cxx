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

#include "AliITSTriggerFOProcessor.h"
#include "AliITSTriggerConditions.h"
#include <TError.h>

/* $Id$ */

AliITSTriggerFOProcessor::AliITSTriggerFOProcessor() :
  fFOInner(0), fFOOuter(0), fTriggerCond(NULL)
{
  // default constructor
}
//______________________________________________________________________
AliITSTriggerFOProcessor::AliITSTriggerFOProcessor(AliITSTriggerConditions* ocdbCond) :
  fFOInner(0), fFOOuter(0), fTriggerCond(ocdbCond)
{
  // constructor
}
//______________________________________________________________________
AliITSTriggerFOProcessor::AliITSTriggerFOProcessor(const AliITSTriggerFOProcessor& handle): 
  fFOInner(handle.fFOInner), fFOOuter(handle.fFOOuter), fTriggerCond(handle.fTriggerCond)
{
  // copy constructor
}
//______________________________________________________________________
AliITSTriggerFOProcessor::~AliITSTriggerFOProcessor() {
  // destructor
}
//______________________________________________________________________
AliITSTriggerFOProcessor& AliITSTriggerFOProcessor::operator=(const AliITSTriggerFOProcessor& handle) {
  // assignment operator
  if (this!=&handle) {
    fFOInner = handle.fFOInner;
    fFOOuter = handle.fFOOuter;
    fTriggerCond = handle.fTriggerCond;
  }
  return *this;
}
//______________________________________________________________________
void AliITSTriggerFOProcessor::SetTriggerConditions(AliITSTriggerConditions* const ocdbCond) {
  // Method to give pointer to the OCDB conditions entry
  fTriggerCond = ocdbCond;
}
//______________________________________________________________________
UInt_t AliITSTriggerFOProcessor::GetNumOutputs() const {
  // return number of outputs (algorithms) in use
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::GetNumOutputs","Trigger conditions entry not yet given.");
    return 0;
  }
  return fTriggerCond->GetNumAlgo();
}
//______________________________________________________________________
const Char_t* AliITSTriggerFOProcessor::GetOutputLabel(Short_t index) const {
  // return the label of the index'th algorithm in use
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::GetOutputLabel","Trigger conditions entry not yet given.");
    return "";
  }
  return fTriggerCond->GetAlgoLabelI(index);
}
//______________________________________________________________________
void AliITSTriggerFOProcessor::PreprocessFOSignals(AliITSFOSignalsSPD* const signals) {
  // Pre-process the fast-or signals to retrieve some data needed by most algorithms
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::PreprocessFOSignals","Trigger conditions entry not yet given.");
    return;
  }
  fFOInner = 0;
  fFOOuter = 0;

  Int_t eq   = -1;
  Int_t hs   = -1;
  Int_t chip = -1;
  while (signals->GetNextSignal(eq,hs,chip)) {
    if (fTriggerCond->IsChipActive( (UInt_t)eq, (UInt_t)hs, (UInt_t)chip) ) {
      if (hs<=1) fFOInner++;
      else       fFOOuter++;
    }
  }
}
//______________________________________________________________________
Bool_t AliITSTriggerFOProcessor::ProcessFOSignalsIndex(Short_t index, AliITSFOSignalsSPD* signals) {
  // NB: For every event - Always call PreprocessFOSignals before calling this method
  // Process index'th algorithm - returns true if output active
  return ProcessFOSignalsLabel(GetOutputLabel(index), signals);
}
//______________________________________________________________________
Bool_t AliITSTriggerFOProcessor::ProcessFOSignalsLabel(const Char_t* label, AliITSFOSignalsSPD* signals) {
  // NB: For every event - Always call PreprocessFOSignals before calling this method
  // Process algorithm with label ... - returns true if output active 
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsLabel","Trigger conditions entry not yet given.");
    return kFALSE;
  }

  Short_t index = fTriggerCond->GetAlgoIndexL(label);
  if (index<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsLabel", "No conditions for label '%s'.",label);
    return kFALSE;
  }
  
  if      (strcmp(label, "0SMB") == 0) return ProcessFOSignalsTHRTotalAndTHRInnerAndTHROuter(index, signals);
  else if (strcmp(label, "0SH1") == 0) return ProcessFOSignalsTHRInnerAndTHROuter(index, signals);
  else if (strcmp(label, "0SH2") == 0) return ProcessFOSignalsTHRInnerAndTHROuter(index, signals);
  else if (strcmp(label, "0SH3") == 0) return ProcessFOSignalsTHRInnerAndTHROuter(index, signals);
  else if (strcmp(label, "0SH4") == 0) return ProcessFOSignalsTHRInnerAndTHROuter(index, signals);
  else if (strcmp(label, "0SPF") == 0) return ProcessFOSignalsTHRTotalAndTHRInnerAndTHROuter(index, signals);
  else if (strcmp(label, "0SX1") == 0) return ProcessFOSignalsInnerGTOuterPlusOffset(index, signals);
  else if (strcmp(label, "0SX2") == 0) return ProcessFOSignalsOuterGTInnerPlusOffset(index, signals);
  else if (strcmp(label, "0SBK") == 0) return ProcessFOSignalsTHRTotal(index, signals);
  else if (strcmp(label, "0SCO") == 0) return ProcessFOSignalsCosmic(index, signals);

  else {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsLabel", "Algorithm not yet implemented for label '%s'.",label);
    return kFALSE;
  }
}
//______________________________________________________________________
Bool_t AliITSTriggerFOProcessor::ProcessFOSignalsTHRTotalAndTHRInnerAndTHROuter(Short_t index, AliITSFOSignalsSPD* /*signals*/) {
  // NB: For every event - Always call PreprocessFOSignals before calling this method
  // Process algorithm ' I+O > thr && I > thr && O > thr ' (index is needed to get the correct parameters from the ocdb object)
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsTHRTotalAndTHRInnerAndTHROuter","Trigger conditions entry not yet given.");
    return kFALSE;
  }
  
  // Get parameter values:
  Int_t thIO = fTriggerCond->GetAlgoParamValueIN(index, "total_threshold");
  if (thIO<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsTHRTotalAndTHRInnerAndTHROuter","Parameter 'total_threshold' not defined");
    return kFALSE;
  }
  Int_t thI = fTriggerCond->GetAlgoParamValueIN(index, "inner_threshold");
  if (thI<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsTHRTotalAndTHRInnerAndTHROuter","Parameter 'inner_threshold' not defined");
    return kFALSE;
  }
  Int_t thO = fTriggerCond->GetAlgoParamValueIN(index, "outer_threshold");
  if (thO<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsTHRTotalAndTHRInnerAndTHROuter","Parameter 'outer_threshold' not defined");
    return kFALSE;
  }

  // Evaluate:
  return (fFOInner+fFOOuter >= thIO && fFOInner >= thI && fFOOuter >= thO);
}
//______________________________________________________________________
Bool_t AliITSTriggerFOProcessor::ProcessFOSignalsTHRInnerAndTHROuter(Short_t index, AliITSFOSignalsSPD* /*signals*/) {
  // NB: For every event - Always call PreprocessFOSignals before calling this method
  // Process algorithm ' I > thr && O > thr ' (index is needed to get the correct parameters from the ocdb object)
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsTHRInnerAndTHROuter","Trigger conditions entry not yet given.");
    return kFALSE;
  }
  
  // Get parameter values:
  Int_t thI = fTriggerCond->GetAlgoParamValueIN(index, "inner_threshold");
  if (thI<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsTHRInnerAndTHROuter","Parameter 'inner_threshold' not defined");
    return kFALSE;
  }
  Int_t thO = fTriggerCond->GetAlgoParamValueIN(index, "outer_threshold");
  if (thO<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsTHRInnerAndTHROuter","Parameter 'outer_threshold' not defined");
    return kFALSE;
  }

  // Evaluate:
  return (fFOInner >= thI && fFOOuter >= thO);
}
//______________________________________________________________________
Bool_t AliITSTriggerFOProcessor::ProcessFOSignalsTHRTotal(Short_t index, AliITSFOSignalsSPD* /*signals*/) {
  // NB: For every event - Always call PreprocessFOSignals before calling this method
  // Process algorithm ' I+O > thr' (index is needed to get the correct parameters from the ocdb object)
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsTHRTotal","Trigger conditions entry not yet given.");
    return kFALSE;
  }
  
  // Get parameter values:
  Int_t thIO = fTriggerCond->GetAlgoParamValueIN(index, "background_threshold_both");
  if (thIO<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsTHRTotal","Parameter 'total_threshold' not defined");
    return kFALSE;
  }

  // Evaluate:
  return (fFOInner + fFOOuter >= thIO);
}
//______________________________________________________________________
Bool_t AliITSTriggerFOProcessor::ProcessFOSignalsInnerGTOuterPlusOffset(Short_t index, AliITSFOSignalsSPD* /*signals*/) {
  // NB: For every event - Always call PreprocessFOSignals before calling this method
  // Process algorithm ' I > O+offset ' (index is needed to get the correct parameters from the ocdb object)
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsGTOuterPlusOffset","Trigger conditions entry not yet given.");
    return kFALSE;
  }
  
  // Get parameter values:
  Int_t offset = fTriggerCond->GetAlgoParamValueIN(index, "background_offset_inner");
  if (offset<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsInnerGTOuterPlusOffset","Parameter 'offset' not defined");
    return kFALSE;
  }

  // Evaluate:
  return (fFOInner >=  fFOOuter + offset);
}
//______________________________________________________________________
Bool_t AliITSTriggerFOProcessor::ProcessFOSignalsOuterGTInnerPlusOffset(Short_t index, AliITSFOSignalsSPD* /*signals*/) {
  // NB: For every event - Always call PreprocessFOSignals before calling this method
  // Process algorithm ' O > I+offset ' (index is needed to get the correct parameters from the ocdb object)
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsGTOuterPlusOffset","Trigger conditions entry not yet given.");
    return kFALSE;
  }
  
  // Get parameter values:
  Int_t offset = fTriggerCond->GetAlgoParamValueIN(index, "background_offset_outer");
  if (offset<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsOuterGTInnerPlusOffset","Parameter 'offset' not defined");
    return kFALSE;
  }

  // Evaluate:
  return (fFOOuter >=  fFOInner + offset);
}
//______________________________________________________________________
Bool_t AliITSTriggerFOProcessor::ProcessFOSignalsCosmic(Short_t index, AliITSFOSignalsSPD* const signals) {
  // NB: For every event - Always call PreprocessFOSignals before calling this method
  // Process algorithm 'cosmic' (index is needed to get the correct parameters from the ocdb object)
  if (fTriggerCond==NULL) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsCosmic","Trigger conditions entry not yet given.");
    return kFALSE;
  }

  // Get parameter values:
  Int_t cosmicParam = fTriggerCond->GetAlgoParamValueIN(index, "cosmic_mode");
  if (cosmicParam<0) {
    Error("AliITSTriggerFOProcessor::ProcessFOSignalsCosmic","Parameter 'cosmic_mode' not defined");
    return kFALSE;
  }

  // Evaluate:

  UShort_t topOuter = 0;
  UShort_t topInner = 0;
  UShort_t bottomOuter = 0;
  UShort_t bottomInner = 0;

  Int_t eq   = -1;
  Int_t hs   = -1;
  Int_t chip = -1;
  while (signals->GetNextSignal(eq,hs,chip)) {
    if (fTriggerCond->IsChipActive( (UInt_t)eq, (UInt_t)hs, (UInt_t)chip) ) {
      if (hs<=1) {
        if (eq%10 < 5) topInner++;
        else        bottomInner++;
      }
      else {
        if (eq%10 < 5) topOuter++;
        else        bottomOuter++;
      }
    }
  }

  // top outer & bottom outer
  if (cosmicParam == 0) return (topOuter>0 && bottomOuter>0);
  // inner & outer
  if (cosmicParam == 1) return (fFOInner>0 && fFOOuter>0);
  // double layer ( >=2 of top inner, top outer, bottom inner, bottom outer )
  if (cosmicParam == 2) {
    UShort_t nHalfLayers = 0;
    if (topOuter>0)    nHalfLayers++;
    if (topInner>0)    nHalfLayers++;
    if (bottomOuter>0) nHalfLayers++;
    if (bottomInner>0) nHalfLayers++;
    return (nHalfLayers>=2);
  }
  // top outer & top inner & bottom outer & bottom inner
  if (cosmicParam == 3) return (topOuter>0 && topInner>0 && bottomOuter>0 && bottomInner>0);
  // top outer & bottom outer & inner
  if (cosmicParam == 4) return (topOuter>0 && bottomOuter>0 && fFOInner>0);
  // global or
  if (cosmicParam == 5) return (fFOOuter>0 || fFOInner>0);

  Error("AliITSTriggerFOProcessor::ProcessFOSignalsCosmic","'cosmic_algorithm_parameter' = %d not defined",cosmicParam);
  return kFALSE;
}

