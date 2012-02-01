#ifndef ALIFMDANALYSISTASKBFCORRELATION_H
#define ALIFMDANALYSISTASKBFCORRELATION_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TArrayI.h"
#include "TH1I.h"
#include "TH2.h"
#include "AliMCEvent.h"
#include "AliFMDFloatMap.h"
#include "TCanvas.h"

/**
 * Task to do the backward/forward correlation analysis
 * 
 * Input:
 *   List of histograms from AliFMDAnaysisTaskBackground
 *
 * Output: 
 *   List of histograms of ...
 *
 * Used correction objects:
 *   
 * 
 * @ingroup FMD_ana
 */
class AliFMDAnalysisTaskBFCorrelation : public AliAnalysisTask
{
public:
  /** 
   * Constructor 
   * 
   * 
   */
    AliFMDAnalysisTaskBFCorrelation();
  /** 
   * Constructor
   * 
   * @param name  Name of task 
   * @param SE    Whether we're run from SE task
   */
    AliFMDAnalysisTaskBFCorrelation(const char* name, Bool_t SE = kTRUE);
  /** 
   * Destructor
   */
  virtual ~AliFMDAnalysisTaskBFCorrelation() {;}
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from
   */
  AliFMDAnalysisTaskBFCorrelation(const AliFMDAnalysisTaskBFCorrelation& o) 
    : AliAnalysisTask(),
      fDebug(o.fDebug),
      fOutputList(0),
      fInputList(0),
      fInternalList(0),
      fVertexString(o.fVertexString),
      fStandalone(o.fStandalone),
      fEvent(0),
      fnBinsX(0),
      fXmin(0),
      fXmax(0),
      fnBinsY(0),
      fYmin(0),
      fYmax(0)
  {}
  /** 
   * Assignment operator 
   * 
   * 
   * @return Reference to this. 
   */  
  AliFMDAnalysisTaskBFCorrelation& 
  operator=(const AliFMDAnalysisTaskBFCorrelation&) { return *this; }
  /**
   * @{
   * @name Implementation of interface methods
   */
  virtual void ConnectInputData(Option_t *option = "");
  virtual void CreateOutputObjects();
  virtual void Init() {}
  virtual void LocalInit() {Init();}
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void SetDebugLevel(Int_t level) {fDebug = level;}
  /** 
   * @}
   */
  /** 
   * Set the input list of histograms
   * 
   * @param inputList Input list 
   */
  void SetInputList(TList* inputList) {fInputList = inputList;}
  /** 
   * Set the input vertex 
   * 
   * @param vtxString String 
   */
  void SetInputVertex(TObjString* vtxString) {fVertexString = vtxString;}
  /** 
   * Set the output list
   * 
   * @param outputList Output list 
   */
  void SetOutputList(TList* outputList) {fOutputList = outputList;}
  /** 
   * Project the data, and mirror it.
   * 
   * @param sType 
   */
  void ProjectAndMirror(TString sType);
  /** 
   * Calculate values 
   * 
   * @param sType 
   */
  void CalculateValues(TString sType);
  //  void ProjectAndMirror(TString type);
  //  void CalculateParameters(TString type);
  /** 
   * The multiplicity versus eta 
   * 
   * @param type 
   */
  void MultiplicityVsEta(TString type);
  /** 
   * Create the response matrix 
   * 
   */
  void CreateResponseMatrix();
  /** 
   * Process a primary hit
   * 
   */
  void ProcessPrimary();
  /** 
   * Get the list out out objects 
   * 
   * 
   * @return 
   */  
  TList* GetOutputList() {return fOutputList;}
private:
  
  Int_t         fDebug;        //  Debug flag
  TList*        fOutputList;   // output list
  TList*        fInputList;    // Input list
  TList*        fInternalList; // Internal list
  TObjString*   fVertexString; // Vertex string
  Bool_t        fStandalone;   // Running standalone?
 
  Int_t   fEvent;              // Event number
  Int_t   fnBinsX;             // Number of bins 
  Float_t fXmin;               // Minimum 
  Float_t fXmax;               // Maximum
  Int_t   fnBinsY;             // Number of bins
  Float_t fYmin;               // Minumum
  Float_t fYmax;               // Maximum

  ClassDef(AliFMDAnalysisTaskBFCorrelation, 0); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
