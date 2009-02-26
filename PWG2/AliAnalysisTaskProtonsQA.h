#ifndef AliAnalysisTaskProtonsQA_cxx
#define AliAnalysisTaskProtonsQA_cxx

// Analysis task creating a the 2d y-p_t spectrum of p and antip
// Author: Panos Cristakoglou
class TList;
class AliESDEvent;
class AliMCEvent;
class AliProtonQAAnalysis;

#include "AliAnalysisTask.h"

class AliAnalysisTaskProtonsQA : public AliAnalysisTask {
 public:
  enum TriggerMode { kMB1 = 0, kMB2, kSPDFASTOR }; 
  enum AnalysisMode { kInvalid = -1, kTPC = 0, kHybrid, kGlobal };
  
  AliAnalysisTaskProtonsQA();
  AliAnalysisTaskProtonsQA(const char *name);
  virtual ~AliAnalysisTaskProtonsQA() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetTriggerMode(TriggerMode triggermode) {fTriggerMode = triggermode;}
  void SetAnalysisMode(AnalysisMode analysismode) {fProtonAnalysisMode = analysismode;}
  void SetAcceptedVertexDiamond(Double_t gVx, Double_t gVy, Double_t gVz) {
    fVxMax = gVx; fVyMax = gVy; fVzMax = gVz;}

  static Bool_t IsEventTriggered(const AliESDEvent *esd, 
				 TriggerMode trigger = kMB2);
  static const  AliESDVertex *GetVertex(AliESDEvent *esd, 
					AnalysisMode mode = kTPC,
					Double_t gVx = 100.,
					Double_t gVy = 100.,
					Double_t gVz = 100.);
 
 private:
  AliESDEvent *fESD;    //ESD object
  AliMCEvent  *fMC;     //MC object

  TList  *fList0; //TList output object
  TList  *fList1; //TList output object
  TList  *fList2; //TList output object
  TList  *fList3; //TList output object
  TList  *fList4; //TList output object
  TList  *fList5; //TList output object
  TList  *fList6; //TList output object
  TList  *fList7; //TList output object
  
  AliProtonQAAnalysis *fProtonQAAnalysis; //analysis object
 
  TriggerMode   fTriggerMode; //Trigger mode
  AnalysisMode  fProtonAnalysisMode; //Analysis mode
  Double_t      fVxMax, fVyMax, fVzMax; //vertex diamond constrain

  AliAnalysisTaskProtonsQA(const AliAnalysisTaskProtonsQA&); // not implemented
  AliAnalysisTaskProtonsQA& operator=(const AliAnalysisTaskProtonsQA&); // not implemented
  
  ClassDef(AliAnalysisTaskProtonsQA, 1); // example of analysis
};

#endif
